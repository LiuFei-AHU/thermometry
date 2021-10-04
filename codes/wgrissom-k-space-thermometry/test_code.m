% load khtdemo_data_cart2;
% [Nx,Ny,Nc,Nt] = size(data);
% inds = 1:4:Nx;
% k = false(Nx,Nx);
% k(inds,:) = true;
% k = k(:,1);
% order = 1;
% mask = true(length(k));
% 
% [yc,xc] = meshgrid(linspace(-1/2,1/2,size(mask,2)), ...
%         linspace(-1/2,1/2,size(mask,1)));
% yc = yc(:);
% xc = xc(:);
% A = [];
% for yp = 0:order
%     for xp = 0:(order-yp)
%         A = [A (xc.^xp).*(yc.^yp)];
%     end
% end
% A = A(mask(:),:);
% A
%  


load theta-2s;

th = theta;
th2 = tth;

% th = reshape(th,[256, 256, 2]);
% 
% tth = cos(th(:,:,1)) + 1i*sin(th(:,:,2));
% tempfull = angle(tth.')*-7.781;
% tempfull(tempfull < 0 | isnan(tempfull)) = 0;


%%%%%
load theta-100;
th2_1 = tth;
% 
% th = theta;
% 
% th = reshape(th,[256, 256, 2]);
% 
% tth = cos(th(:,:,1)) + 1i*sin(th(:,:,2));
% tempfull1 = angle(tth.')*-7.781;
% tempfull1(tempfull1 < 0 | isnan(tempfull1)) = 0;

load theta_mask;
th2_2 = tth;

figure; 
subplot(1,1,1); imagesc(th2,[0 max(max(th2))]); axis image
h = colorbar; ylabel(h,'degrees C'); 
title('iterate: 500');
% subplot(1,3,2); imagesc(th2_1,[0 max(max(th2_1))]); axis image
% h = colorbar; ylabel(h,'degrees C'); 
% title('iterate: 50000');
% subplot(1,3,3); imagesc(th2_2,[0 max(max(th2_2))]); axis image
% h = colorbar; ylabel(h,'degrees C'); 
% title('heating area');

% saveas(gcf,'./temperature_map5.jpg')
