load khtdemo_data_cart3;

heat = sqz(data(:,:,:,2));
baseline = sqz(data(:,:,:,1));

% 傅里叶逆变换
heat=ifft2(heat);
heat=ifftshift(heat);           
% heat=fliplr(heat);
baseline=ifft2(baseline);
baseline=ifftshift(baseline);           
% baseline=fliplr(baseline);
% baseline = baseline *256*256;

% 
% phase_a = atan(imag(heat)./real(heat));
% phase_b = atan(imag(baseline)./real(baseline));
% phase_a = acos(real(heat./abs(heat)));
% phase_b = acos(real(baseline./abs(baseline)));
% phase_diff = phase_a-phase_b;
phase_diff = angle(heat.*conj(baseline));
disp(max(max(phase_diff)));

phase_diff = phase_diff * -7.781;
% phase_diff

% max(max(phase_diff)) % 2.7010
% [x y] = find(phase_diff==max(max(phase_diff)));

phase_diff(phase_diff < 0 | isnan(phase_diff)) = 0;

figure; 
subplot(1,1,1); imagesc(phase_diff,[0 max(max(phase_diff))]); axis image
h = colorbar; ylabel(h,'degrees C'); 
% imshow(phase_diff)

% exp(0.5i)
% cos(0.5)
% acos(0.8776)
% abs( exp(0.5i)*2)

% 证明先单独求两个复数的相位角度做差的结果和直接求A*CONJ(b)的相位角度是一样的！！！！
a = exp(0.5i);
b = exp(0.3i);
ac = acos(real(a/abs(a)));
bc = acos(real(b/abs(b)));
cc = angle(a*conj(b));
ac-bc
cc