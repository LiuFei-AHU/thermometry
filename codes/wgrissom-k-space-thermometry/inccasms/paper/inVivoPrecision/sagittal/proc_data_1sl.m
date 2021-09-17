% water bath
N = 256;
nSlices = 1;
nDummies = 4;
nEncodes = 2;
load P11776.7.mat
load brainMask_1sl
Nt = size(data1,3);
data1 = data1.*repmat(kron(ones(64,1),[1 -1]'),[1 N Nt]);
data1 = conj(data1);

imgSep = zeros(N,N,nSlices,nDummies);
for ii = 1:nDummies
    dataSep = data1(:,:,ii);
    datazp = zeros(N,N,size(dataSep,3));
    datazp((N-128)/2+1:(N-128)/2+128,:,:) = dataSep;
    dataSep = datazp; clear datazp
    imgSep(:,:,1,ii) = ift2(dataSep);
end

% zero pad the data to square
datazp = zeros(N,N,size(data1,3));
datazp((N-128)/2+1:(N-128)/2+128,:,:) = data1;
data1 = datazp; clear datazp

% get kWts
kWts = ones(N*N,1);

% zet up precision calculation
L = permute(imgSep*N*N,[1 2 4 3]);
ct = -1/(3*2*pi*42.58*0.01*0.0139);

parfor ii = nDummies+1:size(data1,3)

    printf('working on dynamic %d',ii);

    dacc = col(data1(:,:,ii));
    acqp = struct('data',dacc);
    acqp.k = true(N,1);
    acqp.kWts = kWts;
    acqp.L = L;
    acqp.L = permute(acqp.L,[3 4 1 2]); % permute to # library images, # slices, space
    acqp.L = permute(acqp.L(:,:,:),[3 1 2]); % collapse spatial dimensions and permute to space, # library images, # slices
    acqp.brainMask = brainMask;

    algp = struct('dofigs',1);
    algp.order = 2; % polynomial order
    algp.modeltest = true;
    algp.lam = 5*10^-4;

    thetainit = zeros(N,N,nSlices); % initialize temp phase shift map with zeros

    [thetakacc,~,~,~,Ac] = kspace_hybrid_thermo_inccasms(acqp,thetainit,algp);
    tempSMS_cd(:,:,:,ii) = ct*real(thetakacc);

end
tempSMS_cd = tempSMS_cd(:,:,:,nDummies+1:end);

save(['results_' mfilename],'tempSMS_cd');
