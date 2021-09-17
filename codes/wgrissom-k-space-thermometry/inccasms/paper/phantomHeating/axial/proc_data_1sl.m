% water bath
N = 256;
nSlices = 1;
nDummies = 3;
pFileName = 'P95232.7';
load([pFileName '.mat'])
load brainMask_1sl
Nt = size(data1,3);
data1 = data1.*repmat(kron(ones(64,1),[1 -1]'),[1 N Nt]);
data1 = conj(data1);

% crop it in freq enc dimension
data1 = fftshift(ifft(data1,[],2),2);
data1(:,1:109,:) = 0;
data1(:,160:end,:) = 0;
data1 = fft(fftshift(data1,2),[],2);

% downsample in freq enc dimension
data1 = data1(:,1:2:end,:) + data1(:,2:2:end,:);
N = 128;

imgSep = zeros(N,N,nSlices,nDummies);
for ii = 1:nDummies
    dataSep = data1(:,:,ii);
    imgSep(:,:,1,ii) = ift2(dataSep);
end

% remove the ones involving the first image
imgSep = imgSep(:,:,:,2:end);

% get kWts
kWts = ones(N*N,1);

% zet up precision calculation
L = permute(imgSep*N*N,[1 2 4 3]);
ct = -1/(3*2*pi*42.58*0.01*0.0139); % phase->temperature

parfor ii = 10:size(data1,3)

    printf('working on dynamic %d',ii);

    dacc = col(data1(:,:,ii));
    acqp = struct('data',dacc);
    acqp.k = true(N,1);
    acqp.kWts = kWts;
    acqp.L = L;
    acqp.L = permute(acqp.L,[3 4 1 2]); % permute to # library images, # slices, space
    acqp.L = permute(acqp.L(:,:,:),[3 1 2]); % collapse spatial dimensions and permute to space, # library images, # slices
    acqp.brainMask = brainMask;

    thetainit = zeros(N,N,nSlices); % initialize temp phase shift map with zeros

    algp = struct('dofigs',1);
    algp.lam = 10^-2;

    [thetakacc,~,~,~,Ac,~,delta] = kspace_hybrid_thermo_inccasms(acqp,thetainit,algp);
    tempSMS_cd(:,:,:,ii) = ct*real(thetakacc);

end
tempSMS_cd = tempSMS_cd(:,:,:,10:end);

save(['results_' mfilename],'tempSMS_cd');
