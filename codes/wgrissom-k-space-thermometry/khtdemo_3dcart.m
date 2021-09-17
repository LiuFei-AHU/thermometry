%| Demonstration of k-space temperature reconstruction algorithm of
%| 4x-undersampled Cartesian/2DFT data
%|
%| Copyright 2015, William A Grissom, Pooja Gaur, Vanderbilt University

addpath('util');

sxy = 0.025; % Gaussian in-plane width
sz = 0.1; % Gaussian through-plane width
ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
phscale = -0.9*pi;           % phase map scaling (negative for heating)
Nxy = 64; Nz = 32; % image dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create circular phantom object with a Gaussian hot spot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rdisc1 = 0.7/2; % disc1 radius is fraction of FOV/2
rdisc2 = 0.9/2; % disc2 radius is fraction of FOV/2
hdisc = 0.9/2; % extent of object in z

xy = -1/2:1/Nxy:1/2-1/Nxy;
z = -1/2:1/Nz:1/2-1/Nz;
[X,Y,Z] = meshgrid(xy,xy,z);

% define intensity 1 inside rdisc1 and 0.5 between rdisc1 and rdisc2
mag = 0.5*double(X.^2+Y.^2 <= rdisc1^2 & abs(Z) <= hdisc) + ...
    0.5*double(X.^2+Y.^2 <= rdisc2^2 & abs(Z) <= hdisc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create complex images with heat-induced phase 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
puls = exp( -(X.^2+Y.^2)/(2*sxy^2) -(Z.^2)/(2*sz^2)); % Gaussian-shaped pulse (focal heating)
puls = puls/max(puls(:));           % normalize pulse 

L = mag; % library image
img = mag.*exp(1i*phscale*puls); % dynamic image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data and decimate it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = fftshift(fftn(fftshift(img)))/Nxy^2/Nz; % fully-sampled data
indsky = unique([1:8:Nxy Nxy/2-2:Nxy/2+2]); % 'acquired' k-space phase sampling locations in y
indskz = unique([1:8:Nz Nz/2-2:Nz/2+2]); % 'acquired' k-space phase sampling locations in z
dacc = data(:,indsky,indskz); % extract sampled lines for this dynamic
k = false(Nxy,Nxy,Nz);
k(:,indsky,indskz) = true;
fprintf('Overall acceleration: %0.1fx\n',Nxy^2*Nz/sum(k(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k-space recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetainit = zeros(Nxy,Nxy,Nz); % initialize temp phase shift map with zeros
acqp.data = dacc(:); % accelerated data
acqp.k = k; % k-space sampling mask
acqp.L = L(:); % baseline 'library'
algp.order = 0; % polynomial order
algp.lam = [10^-3 -1]; % sparsity regularization parameters for real (phase)
% and imag (attenuation) parts of temperature map. -1: force imag part to
% zero
algp.beta = 0;%2^-12; % roughness regularization parameter
algp.useGPU = false;
tic
[thetakacc,~,~,~,Ac] = kspace_hybrid_thermo(acqp,thetainit,algp);
tempkacc = ct*real(thetakacc);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerated recon without k-space hybrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recon the images
G = Gmri_cart(k);
imgacc = reshape(G'*dacc(:),[Nxy Nxy Nz]);
% calculate temperature map
tempacc = ct*angle(imgacc.*conj(L.*exp(1i*Ac)));
tempacc(tempacc < 0 | isnan(tempacc)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
% axial view
subplot(2,4,1); imagesc(ct*phscale*puls(:,:,Nz/2+1),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('True temp map');
subplot(2,4,2); imagesc(tempacc(:,:,Nz/2+1),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('Temp map, 4x-acc: FFT recon');
subplot(2,4,3); imagesc(tempkacc(:,:,Nz/2+1),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('Temp map, 4x-acc: k-Space recon');
subplot(2,4,4); imagesc(ct*phscale*puls(:,:,Nz/2+1)-tempkacc(:,:,Nz/2+1),[-1 1]); axis image
h = colorbar; ylabel(h,'degrees C'); title('k-Space recon error');
% sagittal view
subplot(2,4,5); imagesc(ct*phscale*squeeze(puls(:,Nxy/2+1,:)),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('True temp map');
subplot(2,4,6); imagesc(squeeze(tempacc(:,Nxy/2+1,:)),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('Temp map, 4x-acc: FFT recon');
subplot(2,4,7); imagesc(squeeze(tempkacc(:,Nxy/2+1,:)),[0 ct*phscale]); axis image
h = colorbar; ylabel(h,'degrees C'); title('Temp map, 4x-acc: k-Space recon');
subplot(2,4,8); imagesc(squeeze(ct*phscale*puls(:,Nxy/2+1,:)-tempkacc(:,Nxy/2+1,:)),[-1 1]); axis image
h = colorbar; ylabel(h,'degrees C'); title('k-Space recon error');




