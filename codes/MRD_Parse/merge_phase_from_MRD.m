%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 融合MRD中8通道的相位信息，构造测温算法的所需的数据格式
% V1. 2021.09.10  CB LF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 8个通道的相位数据累加
data = zeros( 256,256,1,11);

root_dir = 'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\';

for d = 4:14
    pathname = [root_dir,num2str(d),'\kspace_data_from_mrd\'];
    maindir = dir(pathname);
    one_data = zeros( 256,256,1);
    for i = 1 : length( maindir )
        if( isequal( maindir(i).name, '.' ) || isequal( maindir(i).name, '..')|| maindir(i).isdir)
            continue;
        end
        disp([pathname, maindir(i).name])
        fname = [pathname, maindir(i).name];
        load(fname);
        % one_data(:,:,i) = KspaceData(:,:,1); % 取第一个slice
        % 按通道累加
        if (i ==1)
            one_data = KspaceData(:,:,1);
        else
            one_data = one_data + KspaceData(:,:,1);
        end
    end
      data(:,:,:,d-3) = one_data;
end
data(:,:,:,1) = data(:,:,:,11);   % 将最后一组加热图作为基线可以去除假体上下位置的噪声
save('khtdemo_data_cart2.mat','data') % 保存数据
disp('Done!')

% data(:,:,:,1)
% 按权重累加
%     one_data2 = zeros( 256,256,4);
%     one_data3 = zeros( 256,256,2);
%     for i = 1:2:8
%         one_data2(:,:,mod(i,2)+1) = one_data(:,:,i)*0.5+one_data(:,:,i+1)*0.5;
%     end
%     one_data3(:,:,1) = one_data2(:,:,1)*0.5+one_data2(:,:,2)*0.5;
%     one_data3(:,:,2) = one_data2(:,:,3)*0.5+one_data2(:,:,4)*0.5;
%     data(:,:,:,d-3) = one_data3(:,:,1)*0.5+one_data3(:,:,2)*0.5;

% 从保存的kspace数据中还原图片
% load 'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\4\kspace_data_from_mrd/Temp1.mat'
% 
% [a,b,c] = size(KspaceData);  % 256 256 4
% 
% data = zeros( a,b,1,2);
% data(:,:,:,1) = KspaceData(:,:,1);
% load 'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\5\kspace_data_from_mrd/Temp1.mat'
% data(:,:,:,2) = KspaceData(:,:,1);
% save('khtdemo_data_cart2.mat','data')

% no_slices = c
% for f=1:no_slices
%         imTemp=fftshift(ifft2(KspaceData(:,:,f)));  % Inverse descrete Fourier transform
% %         imTemp=ifftshift(imTemp);           % Pixel value exchange Left and Right
% %         imTemp=fliplr(imTemp);             % Rotation 旋转
%         imTemp = imTemp*256;
% %         mag = abs(imTemp(:,:,:));
%         data = zeros( a,b,3);
%         data(:,:,1) = abs(imTemp(:,:,:));
%         data(:,:,2) = abs(imTemp(:,:,:));
%         data(:,:,3) = abs(imTemp(:,:,:));
%         
%         figure;
%         imagesc(data);
%         axis off;
%         axis image;
%         caxis([0,3]);
%         colormap(gray);
% end

