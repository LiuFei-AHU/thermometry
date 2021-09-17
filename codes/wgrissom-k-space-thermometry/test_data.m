addpath('util');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data and get the baselines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load khtdemo_data_cart;

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
[Nx,Ny,Nc,Nt] = size(data); % # x,y-locs, coils, dynamics
maxtind = 9; % 9th dynamic has highest temp

inds = 1:4:Nx; % 'acquired' k-space phase sampling locations
dacc = permute(data(inds,:,:,maxtind),[3 1 2]); % extract sampled lines for this dynamic
dacc = dacc(:,:).';

figure; 
subplot(1,2,1); 
imshow(real(data(:,:,:,9)));
% imshow(real(data(:,:,1,9)));
% data(:,:,1,9)
% fftshift(data(:,:,1,9))
% imshow(real(fftshift(ifft2(fftshift(data(:,:,:,9))))))
imwrite((fftshift(ifft2(fftshift(data(:,:,:,9))))),'../n9.jpg')
n9 = imread('../n9.jpg');
n9 = (fft2(n9));
subplot(1,2,2);
imshow(real(n9(:,:,1)));
% n9(:,:,1)




