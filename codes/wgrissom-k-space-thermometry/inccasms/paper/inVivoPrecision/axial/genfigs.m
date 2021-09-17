% display all the tempStd images, stacked
N = 256;
N2 = 184;
dispImg = zeros(3*N,3*N2);
dispImg_cd = dispImg;
mask = dispImg_cd;
rInd = 1;

% one slice
load results_proc_data_1sl;
tempSMS_cd = flip(tempSMS_cd,1);
tempSMS_cd = tempSMS_cd(:,39:222,:,:);

tempSMS_cd_1sl = tempSMS_cd;

mask(N+1:2*N,rInd:rInd+N2-1) = std(tempSMS_cd,[],4) ~= 0;
dispImg_cd(N+1:2*N,rInd:rInd+N2-1) = std(tempSMS_cd,[],4);
printf('1 slice std: %f',std(tempSMS_cd(tempSMS_cd ~= 0)));
printf('1 slice mean: %f',mean(tempSMS_cd(tempSMS_cd ~= 0)));
rInd = rInd+N2;

load results_proc_data_2sl;
tempSMS_cd = flip(tempSMS_cd,1);
tempSMS_cd = tempSMS_cd(:,39:222,:,:);

tempSMS_cd_2sl = tempSMS_cd;

stdTemp = [std(tempSMS_cd(:,:,1,:),[],4); std(tempSMS_cd(:,:,2,:),[],4)];
mask(1:2*N,rInd:rInd+N2-1) = stdTemp ~= 0;
dispImg_cd(1:2*N,rInd:rInd+N2-1) = stdTemp; %reshape(std(tempSMS_cd,[],4),[N 2*N]);
printf('2 slice S15 std: %f',std(tempSMS_cd(tempSMS_cd ~= 0)));
printf('2 slice S15 mean: %f',mean(tempSMS_cd(tempSMS_cd ~= 0)));
tmp = tempSMS_cd(:,:,1,:);
printf('2 slice S15 std, slice 1: %f',std(tmp(tmp ~= 0)));
tmp = tempSMS_cd(:,:,2,:);
printf('2 slice S15 std, slice 2: %f',std(tmp(tmp ~= 0)));
rInd = rInd+N2;

load results_proc_data_3sl;
tempSMS_cd = flip(tempSMS_cd,1);
tempSMS_cd = tempSMS_cd(:,39:222,:,:);

tempSMS_cd_3sl = tempSMS_cd;

%dispImg(rInd:rInd+N-1,:) = reshape(std(tempSMS,[],4),[N 3*N]);
stdTemp = [std(tempSMS_cd(:,:,1,:),[],4); std(tempSMS_cd(:,:,2,:),[],4); std(tempSMS_cd(:,:,3,:),[],4)];
mask(:,rInd:rInd+N2-1) = stdTemp ~= 0; %reshape(std(tempSMS_cd,[],4),[N 3*N]);
dispImg_cd(:,rInd:rInd+N2-1) = stdTemp; %reshape(std(tempSMS_cd,[],4),[N 3*N]);
printf('3 slice 2enc perm2 std: %f',std(tempSMS_cd(tempSMS_cd ~= 0)));
printf('3 slice 2enc perm2 mean: %f',mean(tempSMS_cd(tempSMS_cd ~= 0)));
tmp = tempSMS_cd(:,:,1,:);
printf('3 slice 2enc perm2 std, slice 1: %f',std(tmp(tmp ~= 0)));
tmp = tempSMS_cd(:,:,2,:);
printf('3 slice 2enc perm2 std, slice 2: %f',std(tmp(tmp ~= 0)));
tmp = tempSMS_cd(:,:,3,:);
printf('3 slice 2enc perm2 std, slice 3: %f',std(tmp(tmp ~= 0)));
rInd = rInd+N2;


figure
imagesc(dispImg_cd,[0 2]);axis image;colorbar
axis off;

threeColor = zeros([size(mask) 3]);
map = parula(64); % get a colormap
% interpolate temp errors onto colormap
dispImg_cd(dispImg_cd > 2) = 2;
for kk = 1:3 % loop over color channels
  threeColor(:,:,kk) = reshape(interp1(linspace(0,2,64),map(:,kk),dispImg_cd(:)),size(mask));
end
threeColor = threeColor.*repmat(mask,[1 1 3]);
figure;imagesc(threeColor);
axis image;axis off

load imgSep_3sl
imgSep = imgSep(:,39:222,:);
imgSep = abs(imgSep)./0.5; imgSep(imgSep > 1) = 1;
imgSep = [flipud(imgSep(:,:,1));flipud(imgSep(:,:,2));flipud(imgSep(:,:,3))];
imgSep = repmat(imgSep,[1 1 3]);

finalImg = cat(2,imgSep,threeColor);

load imgAlias_2sl
imgAlias_2sl = imgAlias(:,39:222,:);
imgAlias_2sl = abs(imgAlias_2sl)./1; imgAlias_2sl(imgAlias_2sl > 1) = 1;
imgAlias_2sl = repmat(imgAlias_2sl,[1 1 3]);
imgAlias_2sl = flip(imgAlias_2sl,1);

load imgAlias_3sl
imgAlias_3sl = imgAlias(:,39:222,:);
imgAlias_3sl = abs(imgAlias_3sl)./1; imgAlias_3sl(imgAlias_3sl > 1) = 1;
imgAlias_3sl = repmat(imgAlias_3sl,[1 1 3]);
imgAlias_3sl = flip(imgAlias_3sl,1);

imgAlias = cat(2,imgAlias_2sl,imgAlias_3sl);
imgAlias = cat(2,zeros(N,2*N2,3),imgAlias);

finalImg = cat(1,finalImg,imgAlias);

figure;imagesc(finalImg);axis image;
axis off;

figure;
F(38) = struct('cdata',[],'colormap',[]);
for ii = 1:38
    dispImg = zeros(3*N,3*N2);
    dispImg_cd = dispImg;

    rInd = 1;

    % one slice
    dispImg_cd(N+1:2*N,rInd:rInd+N2-1) = tempSMS_cd_1sl(:,:,:,ii);
    rInd = rInd+N2;

    % two slices
    imgTemp = [tempSMS_cd_2sl(:,:,1,ii); tempSMS_cd_2sl(:,:,2,ii)];
    dispImg_cd(1:2*N,rInd:rInd+N2-1) = imgTemp;
    rInd = rInd+N2;

    % three slices
    imgTemp = [tempSMS_cd_3sl(:,:,1,ii); tempSMS_cd_3sl(:,:,2,ii); tempSMS_cd_3sl(:,:,3,ii)];
    dispImg_cd(:,rInd:rInd+N2-1) = imgTemp;

    threeColor = zeros([size(mask) 3]);
    map = parula(64); % get a colormap
    % interpolate temp errors onto colormap
    dispImg_cd(dispImg_cd > 2) = 2;dispImg_cd(dispImg_cd < -2) = -2;
    for kk = 1:3 % loop over color channels
        threeColor(:,:,kk) = reshape(interp1(linspace(-2,2,64),map(:,kk),dispImg_cd(:)),size(mask));
    end
    threeColor = threeColor.*repmat(mask,[1 1 3]);
    clf;imagesc(threeColor);
    axis image;axis off;
    fontColor = [0.75 0.75 0.75];
    h = text(92,-25,'1 Slice','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(276,-25,'2 Slices','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(460,-25,'3 Slices','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-75,128,'Slice 1','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-75,384,'Slice 2','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = text(-75,640,'Slice 3','HorizontalAlignment','Center','FontSize',20,'FontWeight','bold','Color',fontColor);
    h = colorbar;
    set(h,'Limits',[-2 2],'Ticks',[-2 -1 0 1 2],'TickLabels',{'-2','-1','0','1','2'},'FontSize',14,'Color',fontColor);
    caxis([-2 2]);
    h = text(660,384,'^{\circ}C','FontSize',20,'Color',fontColor);
    axis off
    set(gcf,'color',[0 0 0]);
    drawnow
    F(ii) = getframe(gcf);
    F(ii).cdata = F(ii).cdata(1:400,61:470,:);

end

v = VideoWriter('inVivoAxialTempErrors.mp4','MPEG-4');
v.FrameRate = 7;
v.Quality = 95;
open(v)
writeVideo(v,F)
close(v)
