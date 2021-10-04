% 展示Kspace算法计算的结果,这些数据都来自算法计算过程的输出（在原始代码中添加了这些输出部分）
% 代码是上传到linux服务器上运行的，结果被从服务器上下载到本地解析并显示

% 这是加载运行结果的代码，注意修改路径
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempacc.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempfull.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempkacc.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\thetakacc.mat

% theta 的内容 
disp(size(thetakacc));
disp(['max value in theta:',num2str(max(max(thetakacc)))]);

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
% 设置基线温度，这样温度图显示的就是叠加基线温度的数据
bt = 0; % 例如：26.0; 
mint = bt;
maxt = 20+bt;


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