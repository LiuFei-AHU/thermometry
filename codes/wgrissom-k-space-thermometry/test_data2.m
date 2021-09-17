load khtdemo_data_cart2;
load mask.mat

[Nx,Ny,Nc,Nt] = size(data);

ct = -7.7871;

heat = sqz(data(:,:,:,2));
baseline = sqz(data(:,:,:,1));

% heat = fftshift(ifft2(fftshift(heat)));
% baseline = fftshift(ifft2(fftshift(baseline)));

% 傅里叶逆变换
heat=ifft2(heat);
heat=ifftshift(heat);           
% heat=fliplr(heat);
baseline=ifft2(baseline);
baseline=ifftshift(baseline);           
% baseline=fliplr(baseline);
% baseline = baseline *256*256;

% 正确计算相位差方法 
tmp = angle(heat.*conj(baseline));
tempfull = ct*sum(tmp.*abs(baseline),3)./sum(abs(baseline),3);
tempfull(tempfull < 0 | isnan(tempfull)) = 0;

% 错误的计算相位差的方法
% diff = ((heat) - (baseline));
% diff = heat*conj(baseline);
% diff = diff / abs(diff);
% diff = atan(real(data(:,:,:,2))*real(data(:,:,:,1))/imag(data(:,:,:,2))*imag(data(:,:,:,1)));
% diff = atan((real(data(:,:,:,1))*imag(data(:,:,:,2)) - imag(data(:,:,:,1))*real(data(:,:,:,2)))/(real(data(:,:,:,1))*real(data(:,:,:,2)) + imag(data(:,:,:,1))*imag(data(:,:,:,2))))/(42.58*-0.01*20);
% tempfull = ct*real(sqz(((diff))));
% 
heat = ct*real(sqz(heat));
baseline = ct*real(sqz(baseline));

figure; 
subplot(1,3,1); imagesc(tempfull,[0 18]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map from diff');
subplot(1,3,2); imagesc(heat,[0 18]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map from diff');
subplot(1,3,3); imagesc(baseline,[0 18]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('Temperature map from diff');


% 显示图片
% mag=abs(heat);
% maxi=max(max(mag));
% ff = figure;
% imagesc(mag);
% axis off;
% axis image;
% caxis([0,maxi]);  
% colormap(gray);
% saveas(gcf,'ttt.jpg');

% 计算距离差异
% drift
% dd = zeros(256, 256);
% for i=1:256
%     for j=1:256
%         a = diff(i,j);
%         b = heat(i,j);
%         c = baseline(i,j);
%         dis = sqrt((real(b)-real(c))^2 + (imag(b)-imag(c))^2);
% %         if abs(imag(b)) > 0
% %             drift = abs(imag(a))/abs(imag(b));
% %         else
% %             drift = abs(imag(a));
% %         end
%         dd(i, j) = dis;
%     end
% end
% dd
% max(max(dd))
% dd(120:136,120:136)
% figure; 
% subplot(1,1,1); 
% imshow(dd/max(max(dd))*255)
