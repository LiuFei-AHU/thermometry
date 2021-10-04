load khtdemo_data_cart2;
load mask.mat

[Nx,Ny,Nc,Nt] = size(data);

ct = -7.7871;
% ct = -1;

for i = 2:2
    heat = sqz(data(:,:,:,i));
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
%     
%     heat = pp(heat);
%     baseline = pp(baseline);
    % 正确计算相位差方法 
    tmp = angle(heat.*conj(baseline));
    % tempfull = ct*sum(tmp.*abs(baseline),3)./sum(abs(baseline),3);
    tmp = pp1(tmp);
    tempfull = ct*tmp;  % 这样也可以
    tempfull(tempfull < 0 | isnan(tempfull)) = 0;
    
    % 错误的计算相位差的方法
    % diff = ((heat) - (baseline));
    % diff = heat*conj(baseline);
    % diff = diff / abs(diff);
    % diff = atan(real(data(:,:,:,2))*real(data(:,:,:,1))/imag(data(:,:,:,2))*imag(data(:,:,:,1)));
    % diff = atan((real(data(:,:,:,1))*imag(data(:,:,:,2)) - imag(data(:,:,:,1))*real(data(:,:,:,2)))/(real(data(:,:,:,1))*real(data(:,:,:,2)) + imag(data(:,:,:,1))*imag(data(:,:,:,2))))/(42.58*-0.01*20);
    % tempfull = ct*real(sqz(((diff))));
    % 
    heat = abs(sqz(heat));
    baseline = abs(sqz(baseline));

    figure; 
    subplot(1,3,1); imagesc(tempfull,[0 max(max(tempfull))]); axis image
    h = colorbar; ylabel(h,'degrees C'); 
    title('Temperature map from diff');
    subplot(1,3,2); imagesc(heat,[0 max(max(heat))]); axis image
    h = colorbar; ylabel(h,'degrees C'); 
    title('Temperature map from diff');
    subplot(1,3,3); imagesc(baseline,[0 max(max(baseline))]); axis image
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

    tt = tempfull(80:178,86:176);
    % s = 0;
    % c = 0;
    % for i=1:10
    %    for j=1:10
    %       if tt(i,j)>0
    %           s = s+tt(i,j);
    %           c = c+1;
    %       end
    %    end
    % end
    % s/c
    max(max(tt))
    % median(tt)

    % 1.0144 1.0112 0.9558 0.9770 0.9150 1.0787 1.0233
    % 
end


% tmp = reshape(tmp,[256*256,1]);
%     xw = tmp;
%     xu = tmp;
%     for i = 2:length(xw)
%             diff = xw(i) - xw(i-1);
%             if diff > pi/2
%                 xu(i:end) = xu(i:end) - pi;
%             elseif diff < -pi/2
%                 xu(i:end) = xu(i:end) + pi;
%             end
%     end
%     tmp = reshape(tmp,[256,256]);


function val = pp(f_data)
    cnt = 0;
    mag = abs(f_data);
    tmp = angle(f_data);
    xw = reshape(tmp,[256*256,1]);
    
    while cnt < length(xw)
        cnt = 0;
        for i = 1:length(xw)
            tm = xw(i);
            if xw(i) > pi/2
                xw(i) = xw(i) - pi;
            elseif xw(i) < -pi/2
               xw(i) = xw(i) + pi;
            end
            if tm == xw(i)
               cnt = cnt+1; 
            end
        end
        disp(cnt);
    end
    xw = reshape(xw,[256,256]);
    val = mag.* exp(1i*xw);
end

function val = pp1(f_data)
    cnt = 0;
    tmp = f_data;
    xw = reshape(tmp,[256*256,1]);
    
    while cnt < length(xw)
        cnt = 0;
        for i = 1:length(xw)
            tm = xw(i);
            if xw(i) > pi/2
                xw(i) = xw(i) - pi;
            elseif xw(i) < -pi/2
               xw(i) = xw(i) + pi;
            end
            if tm == xw(i)
               cnt = cnt+1; 
            end
        end
        disp(cnt);
    end
    val = reshape(xw,[256,256]);
end