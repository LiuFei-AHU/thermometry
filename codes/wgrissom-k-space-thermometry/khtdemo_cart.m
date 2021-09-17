%| Demonstration of k-space temperature reconstruction algorithm of
%| 4x-undersampled Cartesian/2DFT data
%|
%| Copyright 2015, William A Grissom, Pooja Gaur, Vanderbilt University

addpath('util');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data and get the baselines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load khtdemo_data_cart2;

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
[Nx,Ny,Nc,Nt] = size(data); % # x,y-locs, coils, dynamics
maxtind = 6; % 9th dynamic has highest temp

% recon the baselines
for ii = 1:Nc
    L(:,:,ii) = fftshift(ifft2(fftshift(sqz(data(:,:,ii,1)))))*Nx*Ny;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decimate the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inds = 1:4:Nx; % 'acquired' k-space phase sampling locations
dacc = permute(data(inds,:,:,maxtind),[3 1 2]); % extract sampled lines for this dynamic
dacc = dacc(:,:).';
k = false(Nx,Nx);
k(inds,:) = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4x-accelerated k-space recon with FFT's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetainit = zeros(Nx,Ny); % initialize temp phase shift map with zeros
acqp.data = dacc; % accelerated data
acqp.k = k(:,1); % k-space sampling mask - passing a vector tells ksh to
% premptively do FFT in fully-sampled dimension, and run in hybrid space
acqp.L = L(:); % baseline 'library'
algp.order = 1; % polynomial order
algp.lam = [10^-2 -1]; % sparsity regularization parameter
algp.beta = 0;%2^-11; % roughness regularization parameter
algp.useGPU = false; % cast all variables to the GPU (Cartesian Only)
algp.stopFrac = 0.001;
tic
[thetakacc,~,~,~,Ac] = kspace_hybrid_thermo(acqp,thetainit,algp);
% kspace_out = [thetakacc,Ac];
% save('kspace_out.mat','kspace_out');
save('thetakacc.mat','thetakacc');
tempkacc = ct*real(thetakacc);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fully-sampled recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recon the images
for ii = 1:Nc
    imgfull(:,:,ii) = fftshift(ifft2(fftshift(sqz(data(:,:,ii,maxtind)))));
end
% calculate an image-magnitude-weighted mean phase shift
tmp = angle(imgfull.*conj(L.*repmat(exp(1i*Ac),[1 1 Nc])));
tempfull = ct*sum(tmp.*abs(L),3)./sum(abs(L),3);
tempfull(tempfull < 0 | isnan(tempfull)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accelerated recon without k-space hybrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recon the images
G = Gmri_cart(k);
for ii = 1:Nc
    imgacc(:,:,ii) = reshape(G'*dacc(:,ii),[Nx Ny]);
end
% calculate an image-magnitude-weighted mean phase shift
tmp = angle(imgacc.*conj(L.*repmat(exp(1i*Ac),[1 1 Nc])));
tempacc = ct*sum(tmp.*abs(L),3)./sum(abs(L),3);
tempacc(tempacc < 0 | isnan(tempacc)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('tempfull.mat','tempfull')
save('tempacc.mat','tempacc')
save('tempkacc.mat','tempkacc')
figure; 
subplot(1,3,1); imagesc(tempfull,[0 18]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, full sampling');
subplot(1,3,2); imagesc(tempacc,[0 18]); axis image 
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, 4x-acc: FFT recon');
subplot(1,3,3); imagesc(tempkacc,[0 18]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, 4x-acc: k-space recon');

