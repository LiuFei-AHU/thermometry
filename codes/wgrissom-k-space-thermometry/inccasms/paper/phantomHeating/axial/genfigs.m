saveFigs = false;

tr = 0.0274; % TR
Nt = 21;
t = (0:Nt-1)*128*tr; % seconds

hsmask = false(128);
hsmask(80:81,68:69) = true;
hsmask = hsmask(54:113,43:102);

% figure plan:
% Vertical, left: show individual slice images
% Vertical, middle: plots in three slices
% Temp maps at 3 time points: 15s, 30s (peak heat), 45s

% load the one slice data
oneSlice = load('results_proc_data_1sl');
oneSlice.tempSMS_cd = oneSlice.tempSMS_cd(54:113,43:102,:,:);

% load the two slice data
twoSlice_I15 = load('results_proc_data_2sl');
twoSlice_I15.tempSMS_cd = twoSlice_I15.tempSMS_cd(54:113,43:102,:,:);

% load the three slice data
threeSlice_perm2_2enc = load('results_proc_data_3sl');
threeSlice_perm2_2enc.tempSMS_cd = circshift(threeSlice_perm2_2enc.tempSMS_cd,[0 4 0 0]);
threeSlice_perm2_2enc.tempSMS_cd = threeSlice_perm2_2enc.tempSMS_cd(54:113,43:102,:,:);

N = 60;

objMask = sum(squeeze(oneSlice.tempSMS_cd) ~= 0,3) ~= 0 | ...
    sum(squeeze(twoSlice_I15.tempSMS_cd(:,:,1,:)) ~= 0,3) ~= 0 | ...
    sum(squeeze(twoSlice_I15.tempSMS_cd(:,:,2,:)) ~= 0,3) ~= 0 | ...
    sum(squeeze(threeSlice_perm2_2enc.tempSMS_cd(:,:,2,:)) ~= 0,3) ~= 0 | ...
    sum(squeeze(threeSlice_perm2_2enc.tempSMS_cd(:,:,3,:)) ~= 0,3) ~= 0 | ...
    sum(squeeze(threeSlice_perm2_2enc.tempSMS_cd(:,:,1,:)) ~= 0,3) ~= 0;
load objMask
objMask = objMask(54:113,43:102);


figure

% tick through two-slice time points and compare in plot and map
Nt = size(oneSlice.tempSMS_cd,4);
%N = 128;
for ii = 9
    img = zeros(N,2*N);
    img(1:N,1:N) = oneSlice.tempSMS_cd(:,:,1,ii);
    img(N+1:2*N,1:2*N) = reshape(twoSlice_I15.tempSMS_cd(:,:,:,ii),[N 2*N]);
    %img(2*N+1:3*N,1:2*N) = reshape(twoSlice_S15.tempSMS_cd(:,:,:,ii),[N 2*N]);
    imagesc(img);axis image;colorbar;drawnow;
    %pause;
end

% tick through three-slice time points and compare in plot and map
for ii = 9
    img = zeros(N,3*N);
    img(1:N,1:N) = oneSlice.tempSMS_cd(:,:,1,ii);
    img(N+1:2*N,1:3*N) = reshape(threeSlice_perm2_2enc.tempSMS_cd(:,:,:,ii),[N 3*N]);
    %img(2*N+1:3*N,1:3*N) = reshape(threeSlice_perm2_3enc.tempSMS_cd(:,:,:,ii),[N 3*N]);
    %img(3*N+1:4*N,1:3*N) = reshape(threeSlice_perm1_2enc.tempSMS_cd(:,:,:,ii),[N 3*N]);
    %img(4*N+1:5*N,1:3*N) = reshape(threeSlice_perm1_3enc.tempSMS_cd(:,:,:,ii),[N 3*N]);
    imagesc(img);axis image;colorbar;drawnow;
    %pause;
end

[x,y] = meshgrid(linspace(-2,2,3*N),linspace(-1,1,3*N));

img = zeros(3*N,3*N);
ti = 9;
img(N+1:2*N,1:N) = oneSlice.tempSMS_cd(:,:,1,ti)';
img(N+1:3*N,N+1:2*N) = reshape(twoSlice_I15.tempSMS_cd(:,:,:,ti),[N 2*N])';
img(1:3*N,2*N+1:3*N) = reshape(threeSlice_perm2_2enc.tempSMS_cd(:,:,:,ti),[N 3*N])';
figure;
imagesc(x(:),y(:),flipud(img),[0 12]);axis image;colorbar;
axis off
drawnow;

load imgSepForDisplay
imgSepForDisplay = circshift(imgSepForDisplay,[0 4 0 0]);
imgSepForDisplay = reshape(imgSepForDisplay(54:113,43:102,:),[N 3*N]).';
imgSepForDisplay = repmat(abs(imgSepForDisplay)./8,[1 3 3]);
imgSepForDisplay(imgSepForDisplay > 1) = 1;

threeColor = zeros([3*N 3*N 3]);
map = hot(64);
img(img > 12) = 12;
for kk = 1:3 % loop over color channels
  threeColor(:,:,kk) = reshape(interp1(linspace(0,12,64),map(:,kk),img(:)),[3*N 3*N]);
  tmp = threeColor(:,:,kk);
  tmp2 = imgSepForDisplay(:,:,kk);
  tmp(img < 0.1) = tmp2(img < 0.1);
  threeColor(:,:,kk) = tmp;
end
%load brainMask_3sl_2enc_perm1_128
%brainMask = circshift(brainMask,[0 4 0 0]);
%brainMask = reshape(brainMask(54:113,43:102,:),[N 3*N]).';
%threeColor = threeColor.*repmat(brainMask,[1 3]);
figure;imagesc(x(:),y(:),flipud(threeColor));
axis image;axis off
drawnow
if saveFigs
    mysavefig(gcf,'maxTemps');
end

figure;colormap hot;colorbar
drawnow
if saveFigs
    mysavefig(gcf,'colorbar');
end

figure;
F(21) = struct('cdata',[],'colormap',[]);
for ii = 1:21

    img = zeros(3*N,3*N);
    ti = ii;
    img(N+1:2*N,1:N) = oneSlice.tempSMS_cd(:,:,1,ti)';
    img(N+1:3*N,N+1:2*N) = reshape(twoSlice_I15.tempSMS_cd(:,:,:,ti),[N 2*N])';
    img(1:3*N,2*N+1:3*N) = reshape(threeSlice_perm2_2enc.tempSMS_cd(:,:,:,ti),[N 3*N])';

    load imgSepForDisplay
    imgSepForDisplay = circshift(imgSepForDisplay,[0 4 0 0]);
    imgSepForDisplay = reshape(imgSepForDisplay(54:113,43:102,:),[N 3*N]).';
    imgSepForDisplay = repmat(abs(imgSepForDisplay)./8,[1 3 3]);
    imgSepForDisplay(imgSepForDisplay > 1) = 1;

    % zero out two other slices for 1-slice scan
    imgSepForDisplay(1:60,1:60,:) = 0;
    imgSepForDisplay(121:end,1:60,:) = 0;

    % zero out third slice for 2-slice scan
    imgSepForDisplay(1:60,61:120,:) = 0;

    threeColor = zeros([3*N 3*N 3]);
    map = hot(64);
    img(img > 12) = 12;
    for kk = 1:3 % loop over color channels
        threeColor(:,:,kk) = reshape(interp1(linspace(0,12,64),map(:,kk),img(:)),[3*N 3*N]);
        tmp = threeColor(:,:,kk);
        tmp2 = imgSepForDisplay(:,:,kk);
        tmp(img < 0.1) = tmp2(img < 0.1);
        threeColor(:,:,kk) = tmp;
    end
    %load brainMask_3sl_2enc_perm1_128
    %brainMask = circshift(brainMask,[0 4 0 0]);
    %brainMask = reshape(brainMask(54:113,43:102,:),[N 3*N]).';
    %threeColor = threeColor.*repmat(brainMask,[1 3]);
    clf;imagesc(x(:),y(:),flipud(threeColor));
    axis image;axis off
    colormap(hot);
    fontColor = [0.75 0.75 0.75];
    h = colorbar;
    set(h,'Limits',[0 12],'Ticks',[0 4 8 12],'TickLabels',{'0','4','8','12'},'FontSize',14,'Color',fontColor);
    caxis([0 12]);
    h = text(2.5,0,'\Delta ^{\circ}C','FontSize',20,'Color',fontColor);
    h = text(-4/3,-1.15,'1 Slice','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(0,-1.15,'2 Slices','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(4/3,-1.15,'3 Slices','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-2.4,-2/3,'Slice 1','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-2.4,0,'Slice 2','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-2.4,2/3,'Slice 3','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);

    h = text(0,1.15,['Time = ' num2str(t(ii),3) 's'],'HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);

    axis off
    set(gcf,'color',[0 0 0],'Position',[1000 918 640 300]);

    drawnow,pause
    F(ii) = getframe(gcf);
    %F(ii).cdata = F(ii).cdata(1:400,61:470,:);




end
v = VideoWriter('phantomAxialTemp.mp4','MPEG-4');
v.FrameRate = 5;
v.Quality = 95;
open(v)
writeVideo(v,F)
close(v)

% show all # slices at time point 20
%img = zeros(N,3*N);
%img(1:N,N+1:2*N) = oneSlice.tempSMS_cd(:,:,1,20);
%img(N+1:2*N,N+1:3*N) = reshape(twoSlice_I15.tempSMS_cd(:,:,:,20),[N 2*N]);
%img(2*N+1:3*N,1:3*N) = reshape(threeSlice_perm2_2enc.tempSMS_cd(:,:,:,20),[N 3*N]);

maxTemp_1sl = [];
maxTemp_2sl = [];
maxTemp_3sl = [];
meanTemp_1sl = [];
meanTemp_2sl = [];
meanTemp_3sl = [];
for ii = 1:Nt
    tmp = oneSlice.tempSMS_cd(:,:,1,ii);
    maxTemp_1sl(ii) = max(tmp(hsmask(:)));
    meanTemp_1sl(ii) = mean(tmp(hsmask(:)));
    tmp = twoSlice_I15.tempSMS_cd(:,:,1,ii);
    maxTemp_2sl(ii) = max(tmp(hsmask(:)));
    meanTemp_2sl(ii) = mean(tmp(hsmask(:)));
    tmp = threeSlice_perm2_2enc.tempSMS_cd(:,:,2,ii);
    maxTemp_3sl(ii) = max(tmp(hsmask(:)));
    meanTemp_3sl(ii) = mean(tmp(hsmask(:)));
end
figure;
subplot(312)
plot(t,maxTemp_1sl);
hold on
plot(t,maxTemp_2sl);
plot(t,maxTemp_3sl);
axis([0 70 0 13]);
% subplot(322)
% plot(meanTemp_1sl);
% hold on
% plot(meanTemp_2sl);
% plot(meanTemp_3sl);

maxTemp_2sl_2nd = [];
maxTemp_3sl_2nd = [];
meanTemp_2sl_2nd = [];
meanTemp_3sl_2nd = [];
for ii = 1:Nt
    tmp = twoSlice_I15.tempSMS_cd(:,:,2,ii);
    maxTemp_2sl_2nd(ii) = max(tmp(hsmask(:)));
    meanTemp_2sl_2nd(ii) = mean(tmp(hsmask(:)));
    tmp = threeSlice_perm2_2enc.tempSMS_cd(:,:,3,ii);
    maxTemp_3sl_2nd(ii) = max(tmp(hsmask(:)));
    meanTemp_3sl_2nd(ii) = mean(tmp(hsmask(:)));
end
subplot(311)
hold on
plot(t,maxTemp_2sl_2nd);
plot(t,maxTemp_3sl_2nd);
axis([0 70 0 13]);
% subplot(324)
% hold on
% plot(meanTemp_2sl_2nd);
% plot(meanTemp_3sl_2nd);

maxTemp_3sl_3rd = [];
meanTemp_3sl_3rd = [];
for ii = 1:Nt
    tmp = threeSlice_perm2_2enc.tempSMS_cd(:,:,1,ii);
    maxTemp_3sl_3rd(ii) = max(tmp(hsmask(:)));
    meanTemp_3sl_3rd(ii) = mean(tmp(hsmask(:)));
end
subplot(313)
hold on
plot(t,maxTemp_3sl_3rd);
axis([0 70 0 13]);
% subplot(324)
% hold on
% plot(meanTemp_3sl_3rd);

drawnow
if saveFigs
    mysavefig(gcf,'tempPlots');
end

load imgAlias_2sl
imgAlias_2sl = imgAlias;
%imgAlias_2sl = abs(imgAlias_2sl)./4; imgAlias_2sl(imgAlias_2sl > 1) = 1;
%imgAlias_2sl = repmat(imgAlias_2sl,[1 1 3]);
%imgAlias_2sl = flip(imgAlias_2sl,1);
load imgAlias_3sl
imgAlias_3sl = imgAlias;
figure;imagesc(flipud(abs([imgAlias_2sl; imgAlias_3sl])'),[0 5]);
axis image;
colormap gray;axis off
drawnow
if saveFigs
    mysavefig(gcf,'aliasedImages');
end
