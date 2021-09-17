function x = mrd_3D(filename,pathname)
    % addpath(pathname);
    reordering = 'seq'; %input('ordering(seq/cen): ');
    
    % 用于保存提取的图片
    curr = [pathname,'\mrd_pared_images\']; %input('current direction: '); 
    if ~exist(curr)
        mkdir(curr);
    end
    % 用于保存提取的kspace数据
    kspace_filename = [pathname,'\kspace_data_from_mrd\'];
    if ~exist(kspace_filename)
        mkdir(kspace_filename)
    end
    
    % disp([pathname,filename]);
    f=fullfile(pathname,filename);  %
    fid = fopen(f,'r');             % Define the file id
    val = fread(fid,4,'int32');     %val为MRD文件前四个数据 
    xdim = val(1); %128
    ydim = val(2); %128
    zdim = val(3); %1
    dim4 = val(4); %3                                                                                                                                   
    fseek(fid,18,'bof');                   %从文件开始(Beginning of file)偏移18个字节的位置
    datatype=fread(fid,1, 'uint16');       %19
    datatype = dec2hex(datatype);          %10 to 16  '19' to '13'  :13
    fseek(fid,48,'bof');                   %从文件开始偏移48个字节的位置
    fseek(fid,152,'bof');                  %从文件开始偏移152个字节的位置
    val = fread(fid,2, 'int32');
    dim5 = val(1);                         %1
    dim6 = val(2);                         %1
    fseek(fid,256,'bof');                     %从文件开始偏移256个字节的位置
    no_samples = xdim;                        %NOS-X:12
    no_views = ydim;                          %NOV-Y:128
    no_views_2=zdim;                          %NOV2-Z:1
    no_slices = dim4;                         %number of slices:3
    no_echoes = dim5;                         %回波:1
    no_expts = dim6;                          %输出:1

    if size(datatype,2)>1                       %datatype'13'的列数2，大于1
        onlydatatype = datatype(2);             %onlydatatype ='3'
        iscomplex = 2;                          %是否为复数：2为复数
    else
        onlydatatype = datatype(1); 
        iscomplex = 1;
    end
    switch onlydatatype
        case '0' 
            dataformat = 'uchar';
        case '1' 
            dataformat = 'schar';
        case '2' 
            dataformat = 'short';
        case '3' 
            dataformat = 'int16';
        case '4' 
            dataformat = 'int32';
        case '5' 
            dataformat = 'float32';
        case '6' 
            dataformat = 'double';
    end

    num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex; %*datasize;
    [m_total, count] = fread(fid,num2read,dataformat); % reading all the data at once
                                                       % m_total is data and count=98304
    if (count~=num2read)
        msgbox('We have a problem...');         %判断数据点是否对应
    end

    if iscomplex == 2
        a=1:num2read/2;             %num2read=98304;a:<1:49152 double>
        m_real = m_total(2*a-1);    %m_real(1.3.5.7.9...)
        m_imag = m_total(2*a);      %m_imag(2.4.6.8...)
        m_C = m_real+m_imag*1i;
    else

        m_C = m_total;
    end
    n=0;
    if strcmp(reordering, 'seq')
        ord=1:no_views;                    %ord:<1:128>
    elseif strcmp(reordering, 'cen')
        for g=1:no_views/2
            ord(2*g-1)=no_views/2+g;
            ord(2*g)=no_views/2-g+1;
        end
    end                              %循环确定ord

    m_C_1 = zeros(no_expts,no_echoes,no_slices,no_views_2,no_views,no_samples); % added by liufei 20210831.
    for a=1:no_expts                    %1
        for b=1:no_echoes               %1
            for c=1:no_slices           %1
                for e=1:no_views       %128      
                    for d=1:no_views_2  %1
                        m_C_1(a,b,c,d,ord(e),:) = m_C(1+n:no_samples+n); % sequential ordering                    
                        n=n+no_samples;                  %m_C_1
                    end
                end
            end
        end
    end
    m_C= squeeze(m_C_1);                   %squeeze函数，删除单一维
    %m_C(1,:,:)=m_C_2

    im=m_C;
    fclose(fid);                            %关闭指针
    KspaceData = zeros(fix(no_views),fix(no_samples),fix(no_slices));
    for t=1:no_slices
        temp = im(t,:,:);
    %     temp = reshape(temp,1*no_samples,no_views);
        KspaceData(:,:,t) = temp;
        temp=squeeze(temp);
    end

    % 将MRD提取出来的kspace数据保存到文件中
    save([kspace_filename,'/',filename(1:end-4),'.mat'],'KspaceData')

    h_slices = no_slices/2;  %% marked by liufei 20210831.what is the purpose of this code?  

    SpatialComplexP=zeros(fix(no_views),fix(no_samples),fix(h_slices));
    %SpatialComplexN=zeros(fix(no_views),fix(no_samples),fix(h_slices));   %创建三维零矩阵
    % added by liufei 20210831. 注释上面的代码，按实际切片保存，避免混淆
    SpatialComplexN=zeros(fix(no_views),fix(no_samples),fix(no_slices));   
    inter = 1;
    %inter = 0;
    comp = 0;
    comn = 0;
    maxi=0;
    I=zeros(no_views,no_samples,no_slices);

    for f=1:no_slices
        imTemp=ifft2(KspaceData(:,:,f));  % Inverse descrete Fourier transform
        imTemp=ifftshift(imTemp);           % Pixel value exchange Left and Right
        imTemp=fliplr(imTemp);             % Rotation 旋转

        if (inter == 0)
            comp = comp+1;
            SpatialComplexP(:,:,comp)=imTemp;      % define Spatial data for Positive current
            inter =1;
        else
            comn = comn+1;
            SpatialComplexN(:,:,comn)=imTemp;      % define Spatial data for Negative current
            %inter = 0;    
            mag = abs(SpatialComplexN(:,:,comn));
            %save(strcat(curr,'\MR\00',num2str(comn),'.png'),'mag','-ascii');%保存

            maxj=max(max(mag)); %% maxi用于控制下面伪彩色图像的值范围的最大值,
            if(maxi>=maxj)
            else
                maxi=maxj;
                [x, ~]=find(mag==maxi);
            end
            I(:,:,f)=mag;
        end
    end

    for f=1:no_slices
        mag=I(:,:,f);
        ff = figure;
        set(ff,'visible','off');
        [M,N1]=size(mag);
        set(ff,'position',[0 0 N1 M]);
        imagesc(mag);
        axis off;
        axis image;
        caxis([0,maxi]);  
        colormap(gray); % marked by liufei .否则是3通道灰度图
        % saveas(ff,[curr,filename,num2str(f),'.jpg']);
        f1 = getframe(ff);
        imwrite(f1.cdata,[curr,filename,num2str(f),'.jpg']);
       % title(['Magnitude Image', ...'Slice: ',num2str(comn)])
    end

    disp('finish')
end




