function x = mrd_3D(filename,pathname)
    % addpath(pathname);
    reordering = 'seq'; %input('ordering(seq/cen): ');
    
    % ���ڱ�����ȡ��ͼƬ
    curr = [pathname,'\mrd_pared_images\']; %input('current direction: '); 
    if ~exist(curr)
        mkdir(curr);
    end
    % ���ڱ�����ȡ��kspace����
    kspace_filename = [pathname,'\kspace_data_from_mrd\'];
    if ~exist(kspace_filename)
        mkdir(kspace_filename)
    end
    
    % disp([pathname,filename]);
    f=fullfile(pathname,filename);  %
    fid = fopen(f,'r');             % Define the file id
    val = fread(fid,4,'int32');     %valΪMRD�ļ�ǰ�ĸ����� 
    xdim = val(1); %128
    ydim = val(2); %128
    zdim = val(3); %1
    dim4 = val(4); %3                                                                                                                                   
    fseek(fid,18,'bof');                   %���ļ���ʼ(Beginning of file)ƫ��18���ֽڵ�λ��
    datatype=fread(fid,1, 'uint16');       %19
    datatype = dec2hex(datatype);          %10 to 16  '19' to '13'  :13
    fseek(fid,48,'bof');                   %���ļ���ʼƫ��48���ֽڵ�λ��
    fseek(fid,152,'bof');                  %���ļ���ʼƫ��152���ֽڵ�λ��
    val = fread(fid,2, 'int32');
    dim5 = val(1);                         %1
    dim6 = val(2);                         %1
    fseek(fid,256,'bof');                     %���ļ���ʼƫ��256���ֽڵ�λ��
    no_samples = xdim;                        %NOS-X:12
    no_views = ydim;                          %NOV-Y:128
    no_views_2=zdim;                          %NOV2-Z:1
    no_slices = dim4;                         %number of slices:3
    no_echoes = dim5;                         %�ز�:1
    no_expts = dim6;                          %���:1

    if size(datatype,2)>1                       %datatype'13'������2������1
        onlydatatype = datatype(2);             %onlydatatype ='3'
        iscomplex = 2;                          %�Ƿ�Ϊ������2Ϊ����
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
        msgbox('We have a problem...');         %�ж����ݵ��Ƿ��Ӧ
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
    end                              %ѭ��ȷ��ord

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
    m_C= squeeze(m_C_1);                   %squeeze������ɾ����һά
    %m_C(1,:,:)=m_C_2

    im=m_C;
    fclose(fid);                            %�ر�ָ��
    KspaceData = zeros(fix(no_views),fix(no_samples),fix(no_slices));
    for t=1:no_slices
        temp = im(t,:,:);
    %     temp = reshape(temp,1*no_samples,no_views);
        KspaceData(:,:,t) = temp;
        temp=squeeze(temp);
    end

    % ��MRD��ȡ������kspace���ݱ��浽�ļ���
    save([kspace_filename,'/',filename(1:end-4),'.mat'],'KspaceData')

    h_slices = no_slices/2;  %% marked by liufei 20210831.what is the purpose of this code?  

    SpatialComplexP=zeros(fix(no_views),fix(no_samples),fix(h_slices));
    %SpatialComplexN=zeros(fix(no_views),fix(no_samples),fix(h_slices));   %������ά�����
    % added by liufei 20210831. ע������Ĵ��룬��ʵ����Ƭ���棬�������
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
        imTemp=fliplr(imTemp);             % Rotation ��ת

        if (inter == 0)
            comp = comp+1;
            SpatialComplexP(:,:,comp)=imTemp;      % define Spatial data for Positive current
            inter =1;
        else
            comn = comn+1;
            SpatialComplexN(:,:,comn)=imTemp;      % define Spatial data for Negative current
            %inter = 0;    
            mag = abs(SpatialComplexN(:,:,comn));
            %save(strcat(curr,'\MR\00',num2str(comn),'.png'),'mag','-ascii');%����

            maxj=max(max(mag)); %% maxi���ڿ�������α��ɫͼ���ֵ��Χ�����ֵ,
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
        colormap(gray); % marked by liufei .������3ͨ���Ҷ�ͼ
        % saveas(ff,[curr,filename,num2str(f),'.jpg']);
        f1 = getframe(ff);
        imwrite(f1.cdata,[curr,filename,num2str(f),'.jpg']);
       % title(['Magnitude Image', ...'Slice: ',num2str(comn)])
    end

    disp('finish')
end




