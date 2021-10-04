% չʾKspace�㷨����Ľ��,��Щ���ݶ������㷨������̵��������ԭʼ�������������Щ������֣�
% �������ϴ���linux�����������еģ�������ӷ����������ص����ؽ�������ʾ

% ���Ǽ������н���Ĵ��룬ע���޸�·��
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempacc.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempfull.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\tempkacc.mat;
load D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\kspace-algrithm\heat\thetakacc.mat

% theta ������ 
disp(size(thetakacc));
disp(['max value in theta:',num2str(max(max(thetakacc)))]);

ct = -7.7871; % degrees C/radian (phase->temp conversion factor)
% ���û����¶ȣ������¶�ͼ��ʾ�ľ��ǵ��ӻ����¶ȵ�����
bt = 0; % ���磺26.0; 
mint = bt;
maxt = 20+bt;


tempfull = tempfull+bt;
tempacc = tempacc + bt;
tempkacc = tempkacc + bt;

% ��ʾͼƬ
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


% ���������ݣ������ã�
% % load khtdemo_data_cart.mat
% % imshow(data(:,:,:,1));
% % data(:,:,:,1)
% 
% % load khtdemo_data_cart2;
% % imshow(data(:,:,:,1))