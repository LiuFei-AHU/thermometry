img4 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\4\image_merged\4-MRD1.jpg');
img5 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\5\image_merged\5-MRD1.jpg');
img6 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\6\image_merged\6-MRD1.jpg');
img7 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\7\image_merged\7-MRD1.jpg');
img8 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\8\image_merged\8-MRD1.jpg');
img9 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\9\image_merged\9-MRD1.jpg');
img10 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\10\image_merged\10-MRD1.jpg');
img11 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\11\image_merged\11-MRD1.jpg');
img12 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\12\image_merged\12-MRD1.jpg');
img13 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\13\image_merged\13-MRD1.jpg');
img14 = imread('D:\personal\thermometry\codes\MRD_Parse\MRD\20210903\14\image_merged\14-MRD1.jpg');

[x,y,z] = size(img4);
data = zeros( y,x,3,10);
% img = fftshift(fft2(fftshift(img14)));
% size(img)
% img4

data(:,:,:,1) = fftshift(fftshift(fft2((img4))));
data(:,:,:,2) = fftshift(fftshift(fft2((img5))));
data(:,:,:,3) = fftshift(fftshift(fft2((img6))));
data(:,:,:,4) = fftshift(fftshift(fft2((img7))));
% data(:,:,:,5) = fftshift(fftshift(fft2((img8))));
% data(:,:,:,6) = fftshift(fftshift(fft2((img9))));
% data(:,:,:,7) = fftshift(fftshift(fft2((img10))));
% data(:,:,:,8) = fftshift(fftshift(fft2((img11))));
% data(:,:,:,9) = fftshift(fftshift(fft2((img12))));
% data(:,:,:,10) =fftshift(fftshift(fft2((img13))));

figure; 
subplot(1,2,1); 
imshow(fftshift(data(:,:,:,1)));
% data(:,:,:,1)

save('khtdemo_data_cart2.mat','data')