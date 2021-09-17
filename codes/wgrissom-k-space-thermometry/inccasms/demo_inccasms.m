%| Demonstration of k-space temperature reconstruction algorithm from
%| two-slice SMS acquisition with incoherent CAIPI blipping
%|
%| Copyright 2018, Kristin Quah, Megan Poorman, and William A Grissom, Vanderbilt University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data and get the baselines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('util');
addpath('../util');
addpath('../');
simulatedData = true; % REDFLAG false = untested!
multiCoil = false;
signPatternType = 'grad'; % switch between 'convCAIPI', 'randBlipCAIPI', 'RF' and 'chaotic'
nSlices = 3; % how many slices to simultaneously image
nLib = 1; % how many library entries to use, for testing multibaseline
inPlaneAcc = 1; % in-plane acceleration factor
noiseStd = 1/50; % standard deviation of noise added to heating data (~1/SNR) 1/50 is a good value
modelTest = false; % mode to measure temperature precision with no actual heating
nNoise = 100; % how many noise realizations to run if modelTest == true
showFigures = true; % display figures
algorithm = 'complexDifference'; % 'complexDifference' or 'phaseDifference'

if multiCoil
    % these sense maps are normalized to an in-brain RMS value of 1
    load 32ch7TRxMaps % 3 mm slices, 2 mm in-plane
    RxB1m3dxyz = permute(RxB1m3dxyz,[2 1 3 4]); % want undersampling A/P
    mask = ssq(RxB1m3dxyz,4) > 0;
    NcAll = size(RxB1m3dxyz,4); % # rx coils
    sz = size(RxB1m3dxyz);
    Nc = 8; % # of compressed channels
    % compress the Rx maps down to the desired # channels
    if Nc < NcAll
        RxB1m3dxyz = permute(RxB1m3dxyz,[4 1 2 3]);
        RxB1m3dxyz = RxB1m3dxyz(:,:).';
        [u,s,v] = svd(RxB1m3dxyz,'econ');
        s = diag(s);s(Nc+1:end) = 0;s = diag(s);
        RxB1m3dxyz = u*s*v';
        RxB1m3dxyz = reshape(RxB1m3dxyz(:,1:Nc),[sz(1:3) Nc]);
    end
else
    Nc = 1; % # rx coils
end

if ~simulatedData
    load khtdemo_data_cart;
    samePhase = true;
    maxtind = 9; % 9th dynamic has highest temp
    dacc = permute(data(:,:,:,maxtind),[3 1 2]); % extract sampled lines for this dynamic
    dacc = dacc(:,:).';
    [Nx,Ny,Nc,Nt] = size(data); % # x,y-locs, coils, dynamics
    
    % recon the baselines & duplicate to other slices
    for ii = 1:Nc
        L(:,:,ii) = fftshift(ifft2(fftshift(sqz(data(:,:,ii,1)))))*Nx*Ny;
    end
    L = repmat(L,[1 1 1 nSlices]);
    
    % Duplicate data & baselines to other slices
    if samePhase
        dacc = repmat(dacc,[1 1 nSlices]);
        hotSpot = repmat(hotSpot,[1 1 nSlices]);
    else
        % copy baseline data to other slices
        dtmp = permute(data(:,:,:,1),[3 1 2]); % extract sampled lines for this dynamic
        dtmp = dtmp(:,:).';
        for ii = 2:nSlices
            dacc(:,:,ii) = dtmp;
        end
        hotSpot = cat(3,hotSpot,zeros([size(hotSpot) nSlices-1]));
    end
else
    
    N = 128;Nx = N;Ny = N;
    yPos = zeros(nSlices,1); % hot spot positions
    
    % define underlying image
    if ~multiCoil
        img = abs(phantom(N,[1 0.72 0.95 0 0 0]));
        img = repmat(img,[1 1 nSlices]);
    else
        % extract flat brain mask as "image" and sensitivity maps
        % for middle three contiguous slices
        img = mask(:,:,23:25);
        sensMaps = RxB1m3dxyz(:,:,23:25,:);
    end
    [x,y] = meshgrid(-N/2:N/2-1);
    sigma = 8*ones(nSlices,1); % stdev of hot spots
    hotSpot = zeros(N,N,nSlices);
    for ii = 1:nSlices
        hotSpot(:,:,ii) = double(~modelTest)*(-3*pi/4)*exp(-(x.^2 + (y-yPos(ii)).^2)./sigma(ii)^2); % 18 degree heat
        hotSpot(:,:,ii) = hotSpot(:,:,ii).*(img(:,:,ii) > 0); % to mask out temp outside volume for error calc
    end
    
    % get hot images as ulderying images times the hot spot phase
    imgHot = img.*exp(1i*hotSpot);
    if multiCoil
        imgHot = repmat(imgHot,[1 1 1 Nc]).*sensMaps; % Nx,Ny,nSlices,Nc
    end
    
    % get k-space data
    dacc = zeros(Nx*Ny,Nc,nSlices);
    for ii = 1:nSlices
        for jj = 1:Nc
            dacc(:,jj,ii) = col(fftshift(fft2(fftshift(imgHot(:,:,ii,jj)))));
        end
    end
    
    % build baseline library
    if ~multiCoil
        L = permute(img,[1 2 4 3])*Nx*Ny; % x/y/#LibImgs/#slices - kspace_hybrid_thermo wants (Nspace*Ncoils) x Nlib x Nslices
        if nLib > 1 % synthesize some library images at different positions, to test multibaseline
            for ii = 2:nLib
                L = cat(3,L,circshift(L,(ii-1)*round(Nx/nLib),2));
            end
        end
        L = permute(L,[3 4 1 2]); % permute to # library images, # slices, space
        L = permute(L(:,:,:),[3 1 2]); % collapse spatial dimensions and permute to space, # library images, # slices
    else
        L = permute(repmat(img,[1 1 1 Nc]).*sensMaps,[1 2 4 5 3])*Nx*Ny; % x/y/#LibImgs/#slices - kspace_hybrid_thermo wants (Nspace*Ncoils) x Nlib x Nslices
        if nLib > 1
            for ii = 2:nLib
                L = cat(4,L,circshift(L,(ii-1)*round(Nx/nLib),2));
            end
        end
        L = permute(L,[4 5 1 2 3]); % permute to # library images, # slices, space
        L = permute(L(:,:,:),[3 1 2]); % collapse spatial & coil dimensions and permute to space, # library images, # slices
    end
    
end

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the 'accelerated' data and sampled k-space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = false(Nx,1);
k(1:inPlaneAcc:end) = true;
dacc = dacc(repmat(k,[Nx 1]),:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the phase shifts & apply to one slice's data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch signPatternType
    case 'chaotic'
        % get the chaotic sign pattern to modulate data with in k-space
        x = 0.5; % initial x
        r = 3.99; % growth rate
        for ii = 1:(nSlices-1)*Nx-1
            x(ii+1) = r*x(ii)*(1-x(ii));
        end
        if nSlices > 1
            kWts = [(-1).^(reshape(x(:),[Nx nSlices-1]) > 0.5) ones(Nx,1)];
        else
            kWts = ones(Nx,1);
        end
        % randomly permute
        %permutation = randperm(N);
        %kWts = kWts(permutation,:);
    case 'convCAIPI'
        % generate phase ramps that shift the slices by fov/sliceIndex
        kWts = exp(1i*2*pi/N*(-N/2:N/2-1)'*(N/nSlices)*(0:nSlices-1));
    case 'grad'
        switch nSlices
            case 1
                encodes = 1;
                kWts = ones(N,1);
            case 2
                encodes = [1 1; 1 1i];% 1 -1]; % works better than [1 -1], i guess since you have high signal for both
                kWts = repmat(encodes,[floor(N/size(encodes,1)) 1]);
            case 3
                encodes = [1 1 1; 1 1i -1];
                kWts = repmat(encodes,[floor(N/size(encodes,1)) 1]);
            case 4
                encodes = [1 1 1 1; 1 (1+1i)/sqrt(2) 1i (-1+1i)/sqrt(2); 1 (-1-1i)/sqrt(2) 1i (1-1i)/sqrt(2)];% 1 (1+1i*sqrt(3))/2 (-1+1i*sqrt(3))/2 -1];
                kWts = repmat(encodes,[floor(N/size(encodes,1)) 1]);
        end
        kWts(end+1:N,:) = encodes(1:N-floor(N/size(encodes,1))*size(encodes,1),:);
        load util/permutation2
        kWts = kWts(permutation,:);
        
    case 'RF'
        switch nSlices
            case 1
                kWts = ones(N,1);
            case 2
                % 2 unique encodings
                encodes = [1 1i; 1 -1i];
                kWts = repmat(encodes,[floor(N/2) 1]);
            case 3
                % 4 unique encodings
                encodes = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1];
                kWts = repmat(encodes,[floor(N/4) 1]);
        end
        kWts(end+1:N,:) = encodes(1:N-floor(N/size(encodes,1))*size(encodes,1),:);
        load util/permutation2
        kWts = kWts(permutation,:);
        
end
% duplicate to all freq enc points
kWts = kron(ones(Nx,1),kWts); % duplicate to all frequency encoded points
kWts = kWts(repmat(k,[Nx 1]),:); % remove unsampled lines
daccCaipi = dacc;%.*repmat(kWts,[1 Nc]);
for ii = 1:nSlices
    daccCaipi(:,:,ii) = bsxfun(@times,dacc(:,:,ii),kWts(:,ii));%.*repmat(kWts,[1 Nc]);
end

if showFigures
    % show slice-collapsed image magnitude and phase
    figure;
    tmp = sum(reshape(daccCaipi,[sum(k) Ny Nc nSlices]),4);
    imagesc(abs(ft2(tmp(:,:,1))));
    axis image,colormap gray,axis off
    figure;
    imagesc(-angle(ft2(tmp(:,:,1))),[-pi/2 pi/2]);
    axis image;axis off;colorbar
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4x-accelerated k-space recon with FFT's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetainit = zeros(Nx,Ny,nSlices); % initialize temp phase shift map with zeros
acqp.k = k(:,1); % k-space sampling mask - passing a vector tells ksh to
% preemptively do FFT in fully-sampled dimension, and run in hybrid space
acqp.kWts = kWts; %(repmat(k,[Nx 1]),1:nSlices); % in case we don't use all the weights for testing
acqp.L = L; % baseline 'library'
acqp.brainMask = true(Nx,Ny,nSlices);
algp.nISTAiters = 50; % for complex difference
algp.order = 0; % polynomial order
if ~modelTest
    switch algorithm
        case 'phaseDifference'
            algp.lam = [10^-4 -1]; % sparsity regularization parameter
        case 'complexDifference'
            algp.lam = 10^-2; % different for complex difference due to different units
    end
    algp.beta = -1; % roughness regularization parameter: -1 turns it off
else
    algp.lam = [0 -1]; % sparsity regularization parameter
    algp.beta = []; % roughness regularization parameter
end
algp.modeltest = modelTest; % measure temperature errors without heating
algp.dofigs = showFigures;
algp.stopFrac = 10^-5;
tic

if ~modelTest
    % get noisy data with hot spot
    acqp.data = sum(daccCaipi,3) + ...
        noiseStd*N/sqrt(2)*(randn(size(sum(daccCaipi,3))) + 1i*randn(size(sum(daccCaipi,3)))); % accelerated data
    % recon temp maps
    switch algorithm
        case 'phaseDifference'
            [thetakacc,~,~,~,Ac] = kspace_hybrid_thermo_inccasms_phsonly(acqp,thetainit,algp);
        case 'complexDifference'
            [thetakacc,~,~,~,Ac] = kspace_hybrid_thermo_inccasms(acqp,thetainit,algp);
    end
else
    cgIters = 10;
    for ii = 1:nNoise
        % get noisy data
        acqp.data = sum(daccCaipi,3) + ...
            noiseStd*N/sqrt(2)*(randn(size(sum(daccCaipi,3))) + 1i*randn(size(sum(daccCaipi,3)))); % accelerated data
        % recon temp maps
        switch algorithm
            case 'phaseDifference'
                [thetakacc(:,:,:,ii),~,~,~,Ac] = kspace_hybrid_thermo_inccasms_phsonly(acqp,thetainit,algp);
            case 'complexDifference'
                [thetakacc(:,:,:,ii),~,~,~,Ac] = kspace_hybrid_thermo_inccasms(acqp,thetainit,algp);
        end
        % calculate image g-factor (conventional approach)
        imgs(:,:,:,ii) = recon_inccaipi_image(acqp,cgIters);
    end
    % calculate g-factor map
    gMap = std(real(imgs),[],4)./(noiseStd/sqrt(2));
end

tempkacc = ct*real(thetakacc);
toc
rmsError = norm(ct*(thetakacc(:)-repmat(hotSpot(:),[size(thetakacc,4) 1])))/N
maxError = max(abs(ct*(thetakacc(:)-repmat(hotSpot(:),[size(thetakacc,4) 1]))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showFigures && ~modelTest
    figure;
    if modelTest
        maxDisplayTemp = 1;minDisplayTemp = -1;
    else
        maxDisplayTemp = max(ct*hotSpot(:));minDisplayTemp = 0;
    end
    for ii = 1:nSlices
        subplot(3,nSlices,ii); imagesc(ct*hotSpot(:,:,ii),[minDisplayTemp maxDisplayTemp]); axis image;axis off
        h = colorbar; ylabel(h,'degrees C');
        title(sprintf('Truth, slice %d',ii));
        subplot(3,nSlices,nSlices+ii); imagesc(tempkacc(:,:,ii),[minDisplayTemp maxDisplayTemp]); axis image;axis off
        h = colorbar; ylabel(h,'degrees C');
        title(sprintf('Incoherent CAIPI, slice %d',ii));
        subplot(3,nSlices,2*nSlices+ii); imagesc(abs(ct*hotSpot(:,:,ii)-tempkacc(:,:,ii)),[minDisplayTemp 1]); axis image;axis off
        h = colorbar; ylabel(h,'degrees C');
        title(sprintf('Error, slice %d',ii));
    end
end
