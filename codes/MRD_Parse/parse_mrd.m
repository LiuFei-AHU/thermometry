%%%
%解析MRD数据脚本，从数据目录中遍历全部MRD,调用解析函数处理，获得图片和kspace数据
%分别保存至对应目录下的相关目录中
%%%
close all;clc;

% [filename,pathname] = uigetfile('*.MRD');
root_dir = 'D:\personal\thermometry\codes\MRD_Parse\MRD\20211029\n7\';
subdir = dir(root_dir);
for s = 1 : length( subdir )
    if( isequal( subdir( s ).name, '.' )||...
        isequal( subdir( s ).name, '..')||...
        ~subdir( s ).isdir)
        continue;
    end
    pathname = [root_dir,subdir(s).name,'\'];
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