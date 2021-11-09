%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 融合MRD中8通道的相位信息，构造测温算法的所需的数据格式
% V1. 2021.09.10  CB LF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 8个通道的相位数据累加
data = zeros( 256,256,1,11);

root_dir = 'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\';
subdir = dir(root_dir);
% subdir = sort_nat({subdir.name});
for s = 1:length(subdir)
    if( isequal( subdir(s).name, '.' ) || isequal( subdir(s).name, '..')|| ~subdir(s).isdir)
            continue;
    end
    pathname = [root_dir,subdir(s).name,'\kspace_data_from_mrd\'];
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
      data(:,:,:,s+1) = one_data;ksapce_data1014_h3.mat
end
data(:,:,:,1) = data(:,:,:,10);   % 将最后一组加热图作为基线可以去除假体上下位置的噪声
save('khtdemo_data_cart6.mat','data') % 保存数据
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
%sort_nat具体内容
function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
    mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end
