%%%
%����MRD���ݽű���������Ŀ¼�б���ȫ��MRD,���ý��������������ͼƬ��kspace����
%�ֱ𱣴�����ӦĿ¼�µ����Ŀ¼��
%%%
close all;clc;

% [filename,pathname] = uigetfile('*.MRD');
root_dir = 'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\';
for d = 4:14
    pathname = [root_dir,num2str(d),'\'];
    maindir = dir( pathname  );
    % dir_name = split(pathname,'\');
    % disp(dir_name(length(dir_name)-1));
    for i = 1 : length( maindir )
        % If maindir(i) is not a dir, skip
        if( isequal( maindir(i).name, '.' ) || isequal( maindir(i).name, '..')|| maindir(i).isdir)
            continue;
        end
        disp([pathname, maindir(i).name])
        x = mrd_3D(maindir(i).name, pathname);
    end
end
disp('Done!')