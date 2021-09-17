

clear all;
close all;clc;

%filename = 'ge_brain1.MRD'; %input('Enter the filename: ');
[filename,pathname] = uigetfile('*.MRD');
reordering = 'seq'; %input('ordering(seq/cen): ');
curr = 'Images\'; %input('current direction: ');

Tc = 20; %input('current the Tc value(ms): ');
% filename = 'ge.MRD';
% reordering = 'seq';
% curr = 'H';
% Tc = 20;
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
scaling=fread(fid,1, 'float32');       %�߶�:32767=��maybe��128*128*2-1
bitsperpixel=fread(fid,1, 'uchar');    %ÿ�����ص�λ��:0
fseek(fid,152,'bof');                  %���ļ���ʼƫ��152���ֽڵ�λ��
val = fread(fid,2, 'int32');
dim5 = val(1);                         %1
dim6 = val(2);                         %1
fseek(fid,256,'bof');                     %���ļ���ʼƫ��256���ֽڵ�λ��
text=fread(fid,256);                      %text:256*1<double>��0..0��
no_samples = xdim;                        %NOS-X:12
no_views = ydim;                          %NOV-Y:128
no_views_2=zdim;                          %NOV2-Z:1
no_slices = dim4;                         %number of slices:3
no_echoes = dim5;                         %�ز�:1
no_expts = dim6;                          %���:1

dim = [no_expts,no_echoes,no_slices,no_views_2,no_views,no_samples];

if size(datatype,2)>1                       %datatype'13'������2������1
    onlydatatype = datatype(2);             %onlydatatype ='3'
    iscomplex = 2;                          %�Ƿ�Ϊ������2Ϊ����
else
    onlydatatype = datatype(1); 
    iscomplex = 1;
end
switch onlydatatype
    case '0' 
        dataformat = 'uchar';   datasize = 1; % size in bytes
    case '1' 
        dataformat = 'schar';   datasize = 1; % size in bytes
    case '2' 
        dataformat = 'short';   datasize = 2; % size in bytes
    case '3' 
        dataformat = 'int16';   datasize = 2; % size in bytes ��
    case '4' 
        dataformat = 'int32';   datasize = 4; % size in bytes
    case '5' 
        dataformat = 'float32'; datasize = 4; % size in bytes
    case '6' 
        dataformat = 'double';  datasize = 8; % size in bytes
end

num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex; %*datasize;
[m_total, count] = fread(fid,num2read,dataformat); % reading all the data at once
                                                   % m_total is data and count=98304
if (count~=num2read)
    h = msgbox('We have a problem...');         %�ж����ݵ��Ƿ��Ӧ
end

if iscomplex == 2
    a=1:num2read/2;             %num2read=98304;a:<1:49152 double>
    m_real = m_total(2*a-1);    %m_real(1.3.5.7.9...)
    m_imag = m_total(2*a);      %m_imag(2.4.6.8...)
    m_C = m_real+m_imag*i;
else
    
    m_C = m_total;
end
n=0;
if reordering == 'seq'
    ord=1:no_views;                    %ord:<1:128>
elseif reordering == 'cen'
    for g=1:no_views/2
        ord(2*g-1)=no_views/2+g;
        ord(2*g)=no_views/2-g+1;
    end
end                              %ѭ��ȷ��ord

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
par=[];                                 % ����

for t=1:no_slices;
    temp = im(t,:,:);
%     temp = reshape(temp,1*no_samples,no_views);
    KspaceData(:,:,t) = temp;
    temp=squeeze(temp);
    
    for x=1:128;
        for y=1:128;
            a=temp(x,y);
            a=abs(a);
            temp(x,y)=a;
        end
    end
    KspaceData_show(:,:,t)=temp;
    
end




%mkdir(curr,'Images')
%mkdir(curr,'before_unwrap_Bz')
%mkdir(curr,'after_unwrap_Bz')               %��curr��H�д����ļ���


h_slices = no_slices/2;    


SpatialComplexP=zeros(fix(no_views),fix(no_samples),fix(h_slices));
SpatialComplexN=zeros(fix(no_views),fix(no_samples),fix(h_slices));   %������ά�����

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
            
            maxj=max(max(mag));
            if(maxi>=maxj)
            else
                maxi=maxj;
                [x y]=find(mag==maxi);
            end
            I(:,:,f)=mag;
            %figure;imagesc(mag);
            %caxis([0,2]);
            %colorbar;
            %colormap(gray);
            %title(['Magnitude Image',...'Slice: ',num2str(comn)]);

        end
end


%SIGNAL-NOISE RATIO
sig=0;
noi=0;
for f=1:no_slices

            mag=I(:,:,f);
           
           % SNR=sig/noi;
        
            figure;imagesc(mag);
            axis off;
            axis image;
            caxis([0,maxi]);
            %colorbar;
            colormap(gray);
            %title(['SNR=',...
            %    num2str(SNR)]);
           
            saveas(gcf,[curr,filename,num2str(f),'.jpg']);
           
           % title(['Magnitude Image', ...'Slice: ',num2str(comn)])
end





%{
phaseTemp = SpatialComplexP./SpatialComplexN;                                                     %����
for jj = 1:h_slices
    [beforePhaseUnwrap,Bdata] = unwrap_phase(phaseTemp(:,:,jj),no_views,Tc,curr);
    before_unwrap(:,:,jj) = beforePhaseUnwrap;
    after_unwrap(:,:,jj) = Bdata; 
    
    %% save phase data
    if (curr == 'H')
        save(strcat(curr,'\before_unwrap_Bz\00',num2str(jj),'.bz1'),'beforePhaseUnwrap','-ascii');
    else
        save(strcat(curr,'\before_unwrap_Bz\00',num2str(jj),'.bz2'),'beforePhaseUnwrap','-ascii');
    end
    
        if (curr == 'H')
        save(strcat(curr,'\after_unwrap_Bz\00',num2str(jj),'.bz1'),'Bdata','-ascii');
    else
        save(strcat(curr,'\after_unwrap_Bz\00',num2str(jj),'.bz2'),'Bdata','-ascii');
    end
end
%}

%{
aft_unw = [];

for jj = 1:h_slices

kk = jj;
%    
%    figure;imagesc(abs(posComp(:,:,jj)));colorbar;colormap(gray)
%     title(['\bf Magnitude Image with positive current injection',...
%         'Slice: ',num2str(kk)]); 
%     figure;imagesc(abs(negComp(:,:,jj)));colorbar;colormap(gray)
%     title(['\bf Magnitude Image with negative current injection',...
%         'Slice: ',num2str(kk)]); 
    figure;imagesc(before_unwrap(:,:,jj),1E-8*[-1.5 1.5]);colorbar;colormap(gray)
    title(['\bf Bz data before unwrapping',...
        'Slice: ',num2str(kk)]);
    figure;imagesc(after_unwrap(:,:,jj),1E-8*[-1.5 1.5]);colorbar;colormap(gray)
    title(['\bf Bz data after unwrapping',...
        'Slice: ',num2str(kk)]);
    
end
%}




