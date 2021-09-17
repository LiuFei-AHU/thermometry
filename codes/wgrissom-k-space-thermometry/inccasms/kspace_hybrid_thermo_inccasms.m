function [theta,A,c,f,Ac,algp,delta] = kspace_hybrid_thermo_inccaipi_cd(acqp,thetainit,algp)

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
%|  thetainit   [Nx,Ny,Nz]  Initial temperature map (real, negative). Real part is
%|                          temperature-induced phase at TE. (optional)
%|  algp    Algorithm parameters structure containing (structure and each entry are optional):
%|              order       1             Polynomial order (default = 0)
%|              lam         [1 2]         l1 penalty weights for real and imaginary parts of m
%|                                        (default = 10^-6)
%|              beta        1             Roughness penalty weight for real
%|                                        and imaginary parts of m (default = 0)
%|                                        for second stage of algorithm. (radians; default = 0.01)
%|              dofigs      1             Display intermediate figures (default = 0)
%|              thiters     1             Number of CG iterations per theta update (default = 10)
%|              citers      1             Number of CG iterations per c update (default = 5)
%|              masknz      [Nx,Ny,Nz]    Mask of non-zero heating
%|                                        locations. This will cause
%|                                        the algorithm to skip the l1-regularized
%|                                        stage and go straight to the masked/unregularized stage
%|                                        (default = [])
%|
%| Outputs:
%|  theta       [Nx,Ny,Nz]    Complex temperature map
%|  A           [Nx*Ny*Nz,Np] Polynomial matrix (may be masked)
%|  c           [Np,Nc]       Polynomial coeffs
%|  f           [Nx,Ny,Nz,Nc] Baseline estimate
%|  Ac          [Nx,Ny,Nc]    Polynomial phase estimate (embedded into original mask)
%|  algp        struct        Final algorithm parameters structure
%|
%| Copyright 2019, William A Grissom, Kristin Quah, Vanderbilt University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('algp','var')
    algp = struct();
end
if ~isfield(algp,'order')
    algp.order = 0; % zeroth-order only (phase drift)
end
if ~isfield(algp,'lam')
    algp.lam = 10^-6; % very small value
end
if ~isfield(algp,'beta')
    algp.beta = -1; % turn off roughness penalty if beta not supplied
end
if ~isfield(algp,'dofigs')
    algp.dofigs = 0;
end
if ~isfield(algp,'thiters')
    algp.thiters = 10; % theta iterations
end
if ~isfield(algp,'citers')
    algp.citers = 5; % c iterations
end
if ~isfield(algp,'FISTAiters')
    algp.FISTAiters = 50; % FISTA iterations
end
if ~isfield(algp,'masknz')
    algp.masknz = [];
else
    if ~isempty(algp.masknz)
        algp.masknz = algp.masknz(acqp.mask);
    end
end
if ~isfield(algp,'maskThresh')
    algp.maskThresh = 0;
end
if ~isfield(algp,'modeltest')
    algp.modeltest = 0;
end
if ~isfield(algp,'sumMask')
    algp.sumMask = false;
end
if ~isfield(algp,'stopFrac')
    algp.stopFrac = 0.0001;
end
if ~isfield(algp,'useGPU')
    algp.useGPU = false;
end
if ~isfield(algp,'useWST')
    algp.useWST = false;
end
if ~isfield(acqp,'dcf')
    acqp.dcf = 1; % optional density compensation function
end
if ~isfield(acqp,'thetaSign')
    acqp.thetaSign = -1; % assume negative phase for positive heating
end
if ~isfield(acqp,'LipConst')
    acqp.LipConst = 0.1; % approximate lipschitz constant of \grad f for complex diff FISTA
end

disp('Performing k-space hybrid thermometry.');

Nc = size(acqp.data,2); % Number of rx coils
nSlices = size(acqp.kWts,2); % Number of simultaneous slices

% calculation data normalization to get on scale with penalties
% we do it here bc the median is sensitive to the Cartesian processing below
%dnorm = median(sqrt(acqp.dcf(:)).*abs(acqp.data(abs(acqp.data) > 0.0001*max(abs(acqp.data(:))))))*sqrt(length(acqp.data(:)));
algp.lam = algp.lam / nSlices;
dnorm = median(sqrt(acqp.dcf(:)).*abs(acqp.data(acqp.data ~= 0)))*sqrt(length(acqp.data(:)));
% In calculating data median, I only consider data within 0.01% of the data max in the above
% bc of a problem I had with the simulated data where I would sometimes
% have identically zero data and sometimes not for very similar settings.
% this shouldn't be necesssary for real data.

if islogical(acqp.k)
    disp('k-space is logical array; using Gmri_cart.');
    acqp.fov = 0;
    % pre-fftshift the k-data if Cartesian, so we don't have to do it
    % iteratively
    if min(size(acqp.k)) > 1
        acqp.mask = true(size(acqp.k));
        if ~exist('thetainit','var')
            thetainit = zeros(size(acqp.k));
        end
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
        if ~exist('thetainit','var')
            thetainit = zeros(length(acqp.k),length(acqp.k),nSlices);
        end
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
        if ~exist('thetainit','var')
            thetainit = zeros(length(acqp.k),length(acqp.k),nSlices);
        end
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
else
    disp('k-space is double array; using Gmri');
    if ~exist('thetainit','var')
        thetainit = zeros(sum(acqp.mask(:)),1);
    end
end

Ns = zeros(nSlices,1);
for ii = 1:nSlices
    Ns(ii) = sum(sum(acqp.mask(:,:,ii),1),2); % Number of spatial locations
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Build objects
%%%%%%%%%%%%%%%%%%%%%%%%%

% build polynomial matrix
A = buildA(acqp.mask,algp.order);

% build system matrix
if ~isfield(acqp,'extraWts')
    G = buildG(Nc,acqp.k,acqp.fov,acqp.mask,acqp.kWts);
else
    G = buildG(Nc,acqp.k,acqp.fov,acqp.mask,acqp.kWts.*acqp.extraWts);
end

% build penalty object
if algp.beta > 0
    R = {};
    for ii = 1:nSlices
        R{ii} = Robject(acqp.mask(:,:,ii),'order',2,'beta',algp.beta,'type_denom','matlab');
    end
else
    R = {};
end

if algp.useWST
    W = Wavelet('Daubechies',4,3); % create Wavelet operator
else
    W = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% get initial f,c,theta
%%%%%%%%%%%%%%%%%%%%%%%%%

% get initial baseline estimate f
c = zeros(size(A{1},2),nSlices,class(A{1}));
Ac = {};
for ii = 1:nSlices
    Ac{ii} = A{ii}*c(:,ii);
end
% get initial theta
theta = {};
for ii = 1:nSlices
    tmp = thetainit(:,:,ii);
    theta{ii} = tmp(acqp.mask(:,:,ii));
    % force negativity
    theta{ii}(acqp.thetaSign*theta{ii} <= 0) = 0;
end
if algp.modeltest
    for ii = 1:nSlices
        theta{ii} = zeros(Ns(ii),1);
    end
end
delta = theta;
f = f_update(acqp.data,Ac,theta,acqp.L,G);

% Normalize data to get on scale with penalties
if dnorm ~= 0
  acqp.data = acqp.data / dnorm;
  f = f / dnorm;
  acqp.L = acqp.L / dnorm;
else
  disp ['Warning: normalization = 0, so not applied. This can ' ...
        'happen when the object has been masked. lam ' ...
        'may need tweaking.'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1-Penalized Complex difference Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(algp.masknz) % only run if we don't already have a heating mask

    itr = 0;
    deltaDiffNorm = Inf;
    deltaNorm = 0;
    fprintf('L1-penalized complex difference iteration %d, delta difference = %f%%\n',itr,deltaDiffNorm/deltaNorm*100);
    while deltaDiffNorm > algp.stopFrac*deltaNorm

        % update baseline image estimates
        % subtract delta signals from data
        deltaSigs = 0;
        for ii = 1:nSlices
            deltaSigs = deltaSigs + G{ii}*delta{ii};
        end
        f = f_update(acqp.data-deltaSigs,Ac,theta,acqp.L,G);

        % update poly coeffs
        c = c_update(acqp.data-deltaSigs,A,c,theta,f,G,algp);
        Ac = {};
        for ii = 1:nSlices
            Ac{ii} = A{ii}*c(:,ii);
        end

        % update temp shift
        deltaOld = delta;
        delta = delta_update(acqp.data,Ac,delta,f,G,algp,algp.lam,acqp.LipConst,W,acqp.mask);

        if algp.dofigs;figure(201);
            maxt = 0;
            for ii = 1:nSlices
                maxt = max(maxt,max(abs(delta{ii})));
            end
            if maxt == 0; maxt = 1; end
            for ii = 1:nSlices
                tmp = zeros(size(acqp.mask(:,:,ii)));tmp(acqp.mask(:,:,ii)) = abs(delta{ii});
                subplot(nSlices*100 + 20 + (ii-1)*2 + 1); imagesc(tmp,[0 maxt]); axis image; title 'Estimated |complex difference|';colorbar
                tmp = zeros(size(acqp.mask(:,:,ii)));tmp(acqp.mask(:,:,ii)) = abs(delta{ii}) > 0;
                subplot(nSlices*100 + 20 + (ii-1)*2 + 2); imagesc(tmp); axis image; title 'Significant |complex difference|';colorbar
            end
            drawnow;
        end

        itr = itr + 1;
        deltaNorm = 0;
        deltaDiffNorm = 0;
        for ii = 1:nSlices
            deltaNorm = deltaNorm + sum(abs(delta{ii}).^2);
            deltaDiffNorm = deltaDiffNorm + sum(abs(delta{ii}-deltaOld{ii}).^2);
        end
        fprintf('L1-penalized complex difference iteration %d, delta difference = %f%%\n',itr,deltaDiffNorm/deltaNorm*100);

    end

    % get a mask of potential temperature shifts.
    algp.masknz = {};
    for ii = 1:nSlices
        % outside brain mask, subtract delta signal from the data
        acqp.data = acqp.data - G{ii}*(col(1-acqp.brainMask(:,:,ii)).*delta{ii});
        % inside brain mask, absorb magnitude changes into baseline images, and absorb
        % phase shift into theta
        theta{ii} = angle(col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii})))+delta{ii}).*conj(col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii}))))));
        theta{ii} = theta{ii}.*col(acqp.brainMask(:,:,ii));
        theta{ii} = theta{ii}.*(acqp.thetaSign*theta{ii} > 0); % enforce positivity
        brainMag = abs(col(acqp.brainMask(:,:,ii)).*col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii})))+delta{ii}));
        for jj = 1:size(acqp.L,2)
            acqp.L(:,jj,ii) = col(1-acqp.brainMask(:,:,ii)).*acqp.L(:,jj,ii) + ...
                brainMag.*exp(1i*angle(acqp.L(:,jj,ii)));
        end
        algp.masknz{ii} = abs(theta{ii}) > abs(algp.maskThresh);
    end
    % update f now that we changed L
    f = f_update(acqp.data,Ac,theta,acqp.L,G);

    if algp.modeltest % get the mask based on significant image signal
        for ii = 1:nSlices
            tmp = sqrt(sum(abs(f(:,:,ii)).^2,2));
            algp.masknz{ii} = (tmp > 0.05*max(tmp)) & col(acqp.brainMask(:,:,ii));
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masked Component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(algp.masknz)

    % run theta_update for nonzero pixels, with no sparsity regularization.
    % we do this because we know that sparsity regularization will
    % attenuate the map somewhat, so we need to relax that effect.
    costOld = Inf;
    cost = cost_eval(acqp.data,G,f,Ac,theta,R);
    itr = 0;
    fprintf('Masked phase difference iteration %d, cost = %f\n',itr,cost);
    while costOld-cost >= algp.stopFrac*costOld

        if ~algp.modeltest
            % update image estimate
            f = f_update(acqp.data,Ac,theta,acqp.L,G);

            % update poly coeffs
            c = c_update(acqp.data,A,c,theta,f,G,algp);
            Ac = {};
            for ii = 1:nSlices
                Ac{ii} = A{ii}*c(:,ii);
            end
        end

        % update temp shift
        theta = theta_update(acqp.data,Ac,theta,f,G,algp.masknz,algp,R,algp.sumMask,algp.modeltest,acqp.thetaSign);

        if algp.dofigs;figure(201);
            maxt = 0;
            for ii = 1:nSlices
                maxt = max(maxt,max(-real(theta{ii}(:))));
            end
            if maxt == 0; maxt = 1; end
            for ii = 1:nSlices
              tmp = zeros(size(acqp.mask(:,:,ii)));tmp(acqp.mask(:,:,ii)) = acqp.thetaSign*real(theta{ii});
              subplot(nSlices*100 + 20 + (ii-1)*2 + 1); imagesc(tmp,[0 maxt]); axis image; title 'Estimated phase';colorbar
              tmp = zeros(size(acqp.mask(:,:,ii)));tmp(acqp.mask(:,:,ii)) = acqp.thetaSign*real(theta{ii}) > 0;
              subplot(nSlices*100 + 20 + (ii-1)*2 + 2); imagesc(tmp); axis image; title 'Significant phase';colorbar
            end
            drawnow;
        end

        % calculate cost with updated parameters
        costOld = cost;
        cost = cost_eval(acqp.data,G,f,Ac,theta,R);

        itr = itr + 1;
        fprintf('Masked phase difference iteration %d, cost = %f\n',itr,cost);

    end
end

% embed final results into full image matrices
Act = zeros(size(acqp.mask));
ft = zeros([size(acqp.mask,1) size(acqp.mask,2) Nc nSlices]);
thetat = zeros(size(acqp.mask));
for ii = 1:nSlices
    tmp = zeros(size(acqp.mask(:,:,ii))); tmp(acqp.mask(:,:,ii)) = gather(theta{ii});thetat(:,:,ii) = tmp;
    tmp = zeros([size(acqp.mask(:,:,ii)) Nc]);
    tmp(repmat(acqp.mask(:,:,1),[1 1 Nc])) = gather(f(:,:,ii)); ft(:,:,:,ii) = tmp;
    tmp = zeros(size(acqp.mask(:,:,1))); tmp(acqp.mask(:,:,1)) = gather(Ac{ii});Act(:,:,ii) = tmp;
end
Ac = Act;
f = ft;
theta = thetat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Build the polynomial matrix
%
function Asl = buildA(mask,order)

nSlices = size(mask,3);

[yc,xc] = meshgrid(linspace(-1/2,1/2,size(mask,2)), ...
    linspace(-1/2,1/2,size(mask,1)));
yc = yc(:);
xc = xc(:);
A = [];
for yp = 0:order
    for xp = 0:(order-yp)
        A = [A (xc.^xp).*(yc.^yp)];
    end
end

Asl = {};
for ii = 1:nSlices
    Asl{ii} = A(col(mask(:,:,ii)),:);
end

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
            G{1} = block_fatrix({diag_sp(kWts(:)),Gsl},'type','mult');
            Gsl = G;

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


%
% Evaluate cost
%
function cost = cost_eval(data,G,f,Ac,theta,R)

nSlices = length(Ac);

err = data(:);
for ii = 1:nSlices
    err = err - G{ii}*col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii}+theta{ii}))));
end
cost = 1/2*real(err'*err);

if exist('R','var')
    if ~isempty(R)
        for ii = 1:nSlices
            cost = cost + R{ii}.penal(R{ii},real(theta{ii})) + R{ii}.penal(R{ii},imag(theta{ii}));
        end
    end
end

%
% Update complex image differences
%
function delta = delta_update(data,Ac,delta,f,G,algp,lam,LipConst,W,mask)

% FISTA with constant stepsize

nSlices = length(G);

% subtract off the drift-corrected baseline signals
b = data(:);
for ii = 1:nSlices
    b = b - G{ii}*col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii}))));
end
%keyboard
y = {};
for ii = 1:nSlices
    y{2,ii} = delta{1,ii}; % y(1) = x(0)
end
t(2) = 1; % stepsize; t(1) = 1

% take some ISTA steps
for ii = 1:algp.FISTAiters

    % update delta based on y; xk = pl(yk), eq 4.1
    res = -b;
    for jj = 1:nSlices
        res = res + G{jj}*y{ii+1,jj};
    end
    % project residual back to image domain and soft-threshold
    %L = 100000;
    for jj = 1:nSlices
        %delta{ii+1,jj} = y{ii+1,jj} - t(ii+1)*(G{jj}'*res);
        % Eqtn 1.4
        delta{ii+1,jj} = y{ii+1,jj} - 2/LipConst*(G{jj}'*res);
        if ~isempty(W)
            % embed delta's in mask, take WT, soft threshold, and mask
            % again
            tmp = zeros(size(mask(:,:,jj)));
            tmp(col(mask(:,:,jj))) = delta{ii+1,jj};
            tmp = W*tmp;
            tmp = exp(1i*angle(tmp)).*max(abs(tmp) - lam/LipConst, 0);
            tmp = W'*tmp;
            delta{ii+1,jj} = tmp(col(mask(:,:,jj)));
        else
            delta{ii+1,jj} = exp(1i*angle(delta{ii+1,jj})).*max(abs(delta{ii+1,jj})-lam/LipConst,0);
        end
    end

    %figure(301);im(reshape(delta{ii+1,1},[128 128]));drawnow;

    t(ii+2) = (1 + sqrt(1+4*t(ii+1)^2))/2; % eq 4.2

    for jj = 1:nSlices % eq. 4.3
        y{ii+2,jj} = delta{ii+1,jj} + (t(ii+1)-1)/t(ii+2)*(delta{ii+1,jj} - delta{ii,jj});
    end

end
tmp = {};
for jj = 1:nSlices
    tmp{1,jj} = delta{end,jj};
end
delta = tmp;



%
% Update heat phase shift vector theta
%
function theta = theta_update(data,Ac,theta,f,G,masknz,algp,R,sumMask,modeltest,thetaSign)

% Polak-Ribiere PCG algorithm from JA Fessler's book, chapter 2, 11.7.13
g = [];
thresh = pi/1000;
nSlices = length(Ac);
for nn = 1:algp.thiters
    gOld = g;
    g = gradcalc_theta(data,Ac,theta,f,G,R);
    for ii = 1:length(g)
        g{ii} = masknz{ii}.*g{ii};
    end
    if exist('sumMask','var')
        % update phase equally over hot spot during masked iterations,
        % since we expect the l1 penalty to just shrink everything by the
        % same amount
        if sumMask == true
            gt = 0;
            for ii = 1:nSlices
                gt = gt + sum(g{ii}(:));
            end
            for ii = 1:nSlices
                g{ii} = gt*masknz{ii};
            end
        end
    end
    if nn == 1
        dir = {};
        for ii = 1:nSlices
            dir{ii} = -g{ii};
        end
    else
        gVec = [];
        gOldVec = [];
        for ii = 1:nSlices
            gVec = [gVec; g{ii}(:)];
            gOldVec = [gOldVec; gOld{ii}(:)];
        end
        gamma = max(0,real(gVec'*(gVec-gOldVec))/real(gOldVec'*gOldVec));
        for ii = 1:nSlices
            dir{ii} = -g{ii} + gamma*dir{ii};
        end
    end
    % dir(dir > 0 & -theta < thresh) = 0; % Fessler 11.11.1
    % WAG 10-18-2016: when running modeltest, also set lam(1) = 0; lam(2) = 10^6
    if ~modeltest
        for ii = 1:nSlices
            dir{ii}(thetaSign*real(dir{ii}) < 0 & thetaSign*real(theta{ii}(:,end)) < thresh) = 1i*imag(dir{ii}(thetaSign*real(dir{ii}) < 0 & thetaSign*real(theta{ii}(:,end)) < thresh));
            dir{ii}(imag(dir{ii}) < 0 & imag(theta{ii}(:,end)) < thresh) = real(dir{ii}(imag(dir{ii}) < 0 & imag(theta{ii}(:,end)) < thresh));
        end
    end

    maxDir = 0;
    for ii = 1:nSlices
        if max(abs(dir{ii})) > maxDir
            maxDir = max(abs(dir{ii}));
        end
    end
    [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,100*min(1,pi/2/maxDir),algp,g);

    if ~modeltest
        z = {};
        constraintViolated = false;
        for ii = 1:nSlices
            z{ii} = theta{ii} + t*dir{ii};
            constraintViolated = constraintViolated || (any(thetaSign*real(z{ii}) < 0) || any(imag(z{ii}) < 0));
        end
      if constraintViolated
          %dir = z.*(z < 0) - theta;
          dir = {};
          for ii = 1:nSlices
            dir{ii} = real(z{ii}).*(thetaSign*real(z{ii}) > 0) - real(theta{ii}(:,end)) + 1i*(imag(z{ii}).*(imag(z{ii}) > 0) - imag(theta{ii}(:,end)));
          end
          [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,1,algp,g);
      end
    end
    if breakOut == true;break;end
    for ii = 1:nSlices
        theta{ii} = theta{ii} + t*dir{ii};
    end
end

%
% Calculate gradient of cost wrt theta
%
function g = gradcalc_theta(data,Ac,theta,f,G,R)

% data fidelity derivatives
Nc = size(f,2);
nSlices = length(Ac);
img = {};
res = data(:);
for ii = 1:nSlices
    img{ii} = col(bsxfun(@times,f(:,:,ii),exp(1i*(Ac{ii}+theta{ii}))));
    res = res - G{ii}*img{ii};
end
g = {};
for ii = 1:nSlices
    g{ii} = real(1i*sum(reshape(conj(img{ii}).*(G{ii}'*res),[length(theta{ii}) Nc]),2));

    if ~isempty(R) % roughness penalty derivatives
        g{ii} = g{ii} + R{ii}.cgrad(R{ii},real(theta{ii}));
    end
end

%
% Calculate theta step size
%
function [t,breakOut] = stepcalc_theta(dir,data,Ac,theta,f,G,R,tmax,algp,thetagrad)

% use boyd's backtracking line search, which usually requires fewer cost evaluations

nSlices = length(Ac);

% calculate current cost
cost = cost_eval(data,G,f,Ac,theta,R);

% line search to get step
costt = cost;
a = 0.5; b = 0.5; t = tmax/b;
thetaGradDir = 0;
for ii = 1:nSlices
    thetaGradDir = thetaGradDir + thetagrad{ii}'*dir{ii};
end
while (costt > cost + a*t*real(thetaGradDir)) && t > 10^-6

    % reduce t
    t = b*t;

    % get test point
    thetat = {};
    for ii = 1:nSlices
        thetat{ii} = theta{ii} + t*dir{ii};
    end

    % calculate cost of test point
    costt = cost_eval(data,G,f,Ac,thetat,R);

end

if t == tmax/b % loop was never entered; return zero step
    t = 0;
end
if cost-costt >= algp.stopFrac*cost
    breakOut = false;
else
    breakOut = true;
end


%
% Update polynomial coefficient vector c
%
function c = c_update(data,A,c,theta,f,G,algp)

g = []; % gradient
for nn = 1:algp.citers
    gold = g;
    g = col(gradcalc_c(data,A,c,theta,f,G));
    if nn == 1
        dir = -g;
    else
        gamma = max(0,real(g'*(g-gold))/real(gold'*gold));
        dir = -g + gamma*dir;
    end
    alpha = stepcalc_c(dir,data,A,c,theta,f,G,min(1,pi/2/max(abs(dir(:)))),g);
    c = c + alpha*reshape(dir,size(c));
end


%
% Calculate gradient of cost wrt c
%
function g = gradcalc_c(data,A,c,theta,f,G)

nSlices = length(A);
Nc = size(f,2);
res = data(:);
for ii = 1:nSlices
    img{ii} = col(bsxfun(@times,f(:,:,ii),exp(1i*(A{ii}*c(:,ii)+theta{ii}))));
    res = res - G{ii}*img{ii};
end
g = zeros(size(c));
for ii = 1:nSlices
    g(:,ii) = A{ii}'*real(sum(reshape(1i*conj(img{ii}).*(G{ii}'*res),...
        [size(A{ii},1) Nc]),2));
end


%
% Calculate step size for c
%
function t = stepcalc_c(dir,data,A,c,theta,f,G,tmax,cgrad)

% use boyd's backtracking line search, which usually requires fewer cost evaluations

nSlices = length(A);
Ac = {};
for ii = 1:nSlices
    Ac{ii} = A{ii}*c(:,ii);
end

% calculate current cost
cost = cost_eval(data,G,f,Ac,theta,[]);

dir = reshape(dir,size(c));

% line search to get step
costt = cost;
a = 0.5; b = 0.5; t = tmax/b;
while (costt > cost + a*t*real(cgrad'*dir(:))) && t > 10^-6

  % reduce t
  t = b*t;

  % get test point
  ct = c + t*dir;

  % calculate cost of test point
  for ii = 1:nSlices
      Act{ii} = A{ii}*ct(:,ii);
  end
  costt = cost_eval(data,G,f,Act,theta,[]);

end

if t == tmax/b % loop was never entered; return zero step
  t = 0;
end

%
% Update baseline image estimates
%
function f = f_update(data,Ac,theta,L,G)

Nc = size(data,2); % # coils
nSlices = length(Ac); % Ac is a cell array of coefficients for each slice

if size(L,2) > 1 % if more than one library image

    % project library images to k-space - sum over slices so we do a joint estimation
    Lk = zeros(size(G{1},1),size(L,2));
    for jj = 1:nSlices % loop over slices
      for ii = 1:size(L,2) % loop over library images
        Lk(:,ii) = Lk(:,ii) + (G{jj}*(L(:,ii,jj).*repmat(exp(1i*(Ac{jj}+theta{jj})),[Nc 1])));
      end
    end
    LtL = real(Lk'*Lk);

    % set up constraint that they must sum to 1
    Ceq = ones(1,size(L,2));
    beq = 1;

    % set up cost
    c = -real(data(:)'*Lk);

    % solve
    options = optimset('Display','off');options.MaxIter = 100000;options.Algorithm = 'interior-point-convex';
    wts = quadprog(double(LtL),double(c),[],[],Ceq,beq,zeros(size(L,2),1),ones(size(L,2),1),[],options);

    % get f
    for jj = 1:nSlices
      f(:,:,jj) = reshape(L(:,:,jj) * wts,[size(L,1)/Nc Nc]);
    end

else

    % only one baseline, so weight vector = 1;
    f = reshape(L,[length(L)/Nc Nc nSlices]);

end


%
% Helper function
%
function out = col(in)

out = in(:);
