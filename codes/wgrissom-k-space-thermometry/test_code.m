load khtdemo_data_cart2;
[Nx,Ny,Nc,Nt] = size(data);
inds = 1:4:Nx;
k = false(Nx,Nx);
k(inds,:) = true;
k = k(:,1);
order = 1;
mask = true(length(k));

[yc,xc] = meshgrid(linspace(-1/2,1/2,size(mask,2)), ...
        linspace(-1/2,1/2,size(mask,1)));
yc = yc(:);
xc = xc(:);
A = [];
for yp = 0:order
    for xp = 0:(order-yp)
        A = [A (xc.^xp).*(yc.^yp)];
    end
end
A = A(mask(:),:);
A
 