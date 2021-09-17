%| Demonstration of k-space temperature reconstruction algorithm
%| 
%| data fields in garad structure:
%| 
%|   (specific to GA radial scans:)
%|   L: baseline library image (columized)
%|   dacc: data with accelerated sampling [samples x Nc]
%|   dfull: data with full sampling [samples x Nc]
%|   fov: field of view for reconstruction
%|   kacc: k-space locations of data with accelerated sampling [samples x number of dimensions]
%|   kfull: k-space locations of data with full sampling [samples x number of dimensions]
%|   mask: mask of field of view
%|   nshotacc: number of lines for accelerated sampling
%|   nshotfull: number of lines for full sampling
%|
%| Copyright 2015, William A Grissom, Pooja Gaur, Vanderbilt University

addpath('util');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data and get the baselines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load khtdemo_data_radial.mat

ct = -7.7871; % degrees C/radian
Nc = size(garad.dacc,2); % Number of rx coils

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerated k-space recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetainit = 0*garad.mask; % initialize with zeros
acqp.data = garad.dacc; % accelerated data
acqp.k = garad.kacc; % undersampled k-space trajectory
acqp.L = garad.L; % baseline 'library'
acqp.mask = garad.mask; % object support mask
acqp.fov = garad.fov; 
acqp.dcf = repmat(abs(acqp.k(:,1)+1i*acqp.k(:,2)),[1 size(acqp.data,2)]); 
% ^ approximate density compensation function
% Using a DCF helps reduce blooming around the edges of a hot spot mask;
% may also reduce the total # of iterations required in some cases
% (unverified).
algp.order = 0; % polynomial order
algp.lam = [0.003 -1]; % sparsity regularization parameter - first entry is sparsity for temp phase, second entry is sparsity for attenuation
algp.beta = 2^-12; % roughness regularization parameter
tic
[thetakacc,~,~,f,Ac] = kspace_hybrid_thermo(acqp,thetainit,algp);
tempkacc = ct*real(thetakacc);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fully-sampled recon using CG and the (f, Ac) estimated by k-space method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = Gmri(garad.kfull,garad.mask,'fov',garad.fov,'basis',{'dirac'});
for jj = 1:Nc
    [xS,info] = qpwls_pcg(zeros(sum(garad.mask(:)),1),G,1,garad.dfull(:,jj),0,0,1,25,garad.mask);
    imgfull(:,:,jj) = embed(xS(:,end),garad.mask);
end
tmp = angle(imgfull.*conj(f.*repmat(exp(1i*Ac),[1 1 Nc])));
tempfull = ct*sum(tmp.*abs(f),3)./sum(abs(f),3);
tempfull(tempfull < 0 | isnan(tempfull)) = 0;

% reshape matrix of k-space locations for plotting
tmpkacc = reshape(garad.kacc,[size(garad.kacc,1)/garad.nshotacc garad.nshotacc 2]);
tmpkfull = reshape(garad.kfull,[size(garad.kfull,1)/garad.nshotfull garad.nshotfull 2]);
d = garad.fov*(-1/2:1/(size(garad.mask,1)-1):1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
subplot(2,2,1); plot(squeeze(complexify(tmpkfull(:,:,1)+1i*tmpkfull(:,:,2))),'k'); 
axis([-2.5 2.5 -2.5 2.5]); xlabel('k_x [1/cm]'); ylabel('k_y [1/cm]'); axis square;
title(sprintf('GA radial kspace, full: %d lines',garad.nshotfull));
subplot(2,2,2); imagesc(d,d,tempfull,[0 18]); axis image; colormap jet; h = colorbar;
ylabel(h,'degrees C'); xlabel('x [cm]'); ylabel('y [cm]');
title('Temperature map, full: CG reconstruction');
subplot(2,2,3); plot(squeeze(complexify(tmpkacc(:,:,1)+1i*tmpkacc(:,:,2))),'k'); 
axis([-2.5 2.5 -2.5 2.5]); xlabel('k_x [1/cm]'); ylabel('k_y [1/cm]'); axis square;
title(sprintf('GA radial kspace, acc: %d lines',garad.nshotacc));
subplot(2,2,4); imagesc(d,d,tempkacc,[0 18]); axis image; colormap jet; h = colorbar;
ylabel(h,'degrees C'); xlabel('x [cm]'); ylabel('y [cm]');
title('Temperature map, acc: k-space reconstruction');