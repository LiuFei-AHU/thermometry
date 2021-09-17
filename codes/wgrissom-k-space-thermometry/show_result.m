load tempacc.mat;
load tempfull.mat;
load tempkacc.mat;
load thetakacc.mat

% theta 的内容 
size(thetakacc)
max(max(thetakacc))

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
bt = 26.0;
mint = bt;
maxt = 18+bt;

tempfull = tempfull+bt;
tempacc = tempacc + bt;
tempkacc = tempkacc + bt;

% 显示图片
figure; 
subplot(1,3,1); imagesc(tempfull,[mint maxt]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, full sampling');
subplot(1,3,2); imagesc(tempacc,[mint maxt]); axis image 
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, 4x-acc: FFT recon');
subplot(1,3,3); imagesc(tempkacc,[mint maxt]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map, 4x-acc: k-space recon');


% 打开数据内容（调试用）
% % load khtdemo_data_cart.mat
% % imshow(data(:,:,:,1));
% % data(:,:,:,1)
% 
% % load khtdemo_data_cart2;
% % imshow(data(:,:,:,1))