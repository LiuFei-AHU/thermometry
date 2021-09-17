function img = recon_inccaipi_image(acqp,imgIters)

%|function kspace_hybrid_thermo
%|
%| Inputs:
%|  acqp    Acquisition parameters structure containing (required):
%|              data        [Nk,Nc]       Nc complex k-space data vectors
%|              k           [Nk,Nd]       Nd k-space sample vectors (cycles/cm)
%|                       OR [Nkx,Nky,Nkz] logical Cartesian k-space sampling mask
%|              fov         1             Field of view (cm) (Non-Cartesian only)
%|              mask        [Nx,Ny,Nz]    Binary mask over the FOV (Non-Cartesian only)
%|              L           [Nx*Ny*Nz*Nc,Nl] Multibaseline image library
%|  imgIters    # CG iterations to use for image reconstruction
%|
%| Copyright 2014-05-19, William A Grissom, Pooja Gaur, Vanderbilt University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('imgIters','var')
    imgIters = 20;
end

disp('Performing inccaipi image reconstruction for g-factor measurement.');

Nc = size(acqp.data,2); % Number of rx coils
nSlices = size(acqp.kWts,2); % Number of simultaneous slices

if islogical(acqp.k)
    disp('k-space is logical array; using Gmri_cart.');
    acqp.fov = 0;
    % pre-fftshift the k-data if Cartesian, so we don't have to do it
    % iteratively
    if min(size(acqp.k)) > 1
        acqp.mask = true(size(acqp.k));
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShift = zeros(size(acqp.k));
            dataShift(acqp.k) = acqp.data(:,ii);
            dataShift = fftshift(dataShift);
            acqp.data(:,ii) = dataShift(kShift);
            clear dataShift
        end
        acqp.kWts = fftshift(acqp.kWts);
        acqp.k = kShift; % Now that the data is shifted, we need to always use shifted k mask
        clear kShift
    end
    if size(acqp.k,1) == 1 % we have undersampling only in column dim
        acqp.mask = true(length(acqp.k),length(acqp.k),nSlices); % assume square k-space and images for Cartesian
        % fft the data to image domain in row dim and
        % pre-fftshift the k-space data in column dim
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShiftHybrid = zeros(length(acqp.k)); % assume square k-space + images
            dataShiftHybrid(:,acqp.k) = reshape(acqp.data(:,ii),[length(acqp.k) sum(acqp.k)]);
            dataShiftHybrid = fftshift(ifft(fftshift(dataShiftHybrid),[],1),1)*length(acqp.k); % ifft in row dim
            acqp.data(:,ii) = col(dataShiftHybrid(:,kShift));
            clear dataShiftHybrid;
        end
        for jj = 1:nSlices
            acqp.kWts = col(fftshift(reshape(acqp.kWts(:,jj),[length(acqp.k) sum(acqp.k)])));
        end
        % fftshift the mask
        acqp.k = kShift;
        clear kShift
    end
    if size(acqp.k,2) == 1 % we have undersampling only in row dim
        acqp.mask = true(length(acqp.k),length(acqp.k),nSlices); % assume square k-space and images for Cartesian
        % fft the data to image domain in column dim and
        % pre-fftshift the k-space data in row dim
        kShift = fftshift(acqp.k);
        for ii = 1:Nc
            dataShiftHybrid = zeros(length(acqp.k)); % assume square k-space + images
            dataShiftHybrid(acqp.k,:) = reshape(acqp.data(:,ii),[sum(acqp.k) length(acqp.k)]);
            dataShiftHybrid = fftshift(ifft(fftshift(dataShiftHybrid,2),[],2),2)*length(acqp.k); % ifft in column dim
            dataShiftHybrid = fftshift(dataShiftHybrid,1); % pre-fftshift in row dim
            acqp.data(:,ii) = col(dataShiftHybrid(kShift,:));
            clear dataShiftHybrid;
        end
        for jj = 1:nSlices
            acqp.kWts(:,jj) = col(fftshift(reshape(acqp.kWts(:,jj),[sum(acqp.k) length(acqp.k)])));
        end
        % fftshift the mask
        acqp.k = kShift;
        clear kShift
    end
end

Ns = zeros(nSlices,1);
for ii = 1:nSlices
    Ns(ii) = sum(sum(acqp.mask(:,:,ii),1),2); % Number of spatial locations
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Build objects
%%%%%%%%%%%%%%%%%%%%%%%%%


% build system matrix
G = buildG(1,acqp.k,acqp.fov,acqp.mask,acqp.kWts);
% apply the baselines as sensitivities - assume one baseline
% G should be (Nc*Nk) x (Nx*Ny*Nslice)
L = reshape(acqp.L,...
    [sqrt(length(acqp.L(:,1,ii))/Nc) sqrt(length(acqp.L(:,1,ii))/Nc) Nc size(acqp.L,2) size(acqp.L,3)]);
%L = squeeze(L);
G = block_fatrix(G,'type','row');
Gall = {};
for ii = 1:Nc
    % get sense maps for this coil, for all slices
    tmp = permute(squeeze(L(:,:,ii,:)),[3 1 2]);
    tmp = tmp(:,:).';tmp = tmp(:);
    Gall{ii} = block_fatrix({G,spdiag(tmp)},'type','mult');
end
Gall = block_fatrix(Gall,'type','col');

% do a CG reconstruction of the data
xS = qpwls_pcg(zeros(numel(acqp.mask),1),Gall,1,acqp.data,0,0,1,imgIters);

img = reshape(xS(:,end),size(acqp.mask));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Build the system matrices
%
function G = buildG(Nc,k,fov,mask,kWts)

if ~isempty(k) % image-domain

    nSlices = size(kWts,2);
    Gsl = {};
    for jj = 1:nSlices
        if ~islogical(k) % non-cartesian

            % build system matrix
            if size(k,3) == 1 % 1 shot
                Gsl = Gmri(k,mask,'fov',fov,'basis',{'dirac'});
            else % multishot
                nshot = size(k,3);
                for ii = 1:nshot % build a system matrix for each shot
                    Gsub{ii} = Gmri(k(:,:,ii),mask,'fov',fov,'basis',{'dirac'});
                end
                Gsl = block_fatrix(Gsub,'type','col');
            end

        else % cartesian

            kfftShift = false; % switch to do second fftshift, in frequency domain
            % If false, must make sure k-space data is not centered before starting
            % algorithm
            Gsl{jj} = Gmri_cart(k,[],kfftShift,[],kWts(:,jj));

        end

        if Nc > 1 % multiple coils; replicate the nufft's into a block-diag matrix
            tmp = {};
            for ii = 1:Nc
                tmp{ii} = Gsl{jj};
            end
            Gsl{jj} = block_fatrix(tmp,'type','diag');
        end
    end
    %if nSlices == 1
    %    G = Gsl{1}; % convert back to a single fatrix
    %else
        G = Gsl; % keep it as a cell array
    %end

else
    G = 1; % image domain - no transform
end

