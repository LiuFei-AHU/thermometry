% water bath
N = 256;
nSlices = 2;
nDummies = 4;
nEncodes = 2;
load P06144.7.mat
Nt = size(data1,3);
data1 = data1.*repmat(kron(ones(64,1),[1 -1]'),[1 N Nt]);
data1 = conj(data1);

sliceMid = 0; % mm
sliceGap = 30; % mm
phs = pi/2;
A = [1 1;
    exp(1i*phs/2+1i*(sliceMid)*phs/sliceGap) ...
    exp(-1i*phs/2+1i*(sliceMid)*phs/sliceGap)]; % encoding functions

load permutation2
pulse_index_baseline = zeros(nSlices,128); % which pulses are used for each line of each baseline acquisition
pulse_index_baseline(1,:) = kron(ones(1,64),[1 2]);
pulse_index_baseline(2,:) = kron(ones(1,64),[2 1]);
pulse_index = mod(permutation-1,2)+1;

imgSep = zeros(N,N,nSlices,nSlices*(nDummies+1));
for ii = 1:nDummies
    %jj = ii; % one baseline per dummy; no mixing
    for jj = 1:nDummies
        dataSep = zeros(128,nSlices,N);
        for ll = 1:128
            pulses = pulse_index_baseline(:,ll);
            % what was the pulse order for this
            At = A(pulses,:); % get the active A matrix for this phase encode
            % collect the data into a vector
            datat = squeeze(data1(ll,:,[(ii-1)*nSlices+1 (jj-1)*nSlices+2])).'; % put the third index first
            % invert A to get the separated slice data
            dataSep(ll,:,:) = At\datat;
        end
        dataSep = permute(dataSep,[1 3 2]);
        datazp = zeros(N,N,size(dataSep,3));
        datazp((N-128)/2+1:(N-128)/2+128,:,:) = dataSep;
        dataSep = datazp; clear datazp


        for ll = 1:nSlices
            imgSep(:,:,ll,(ii-1)*nSlices+jj) = ift2(dataSep(:,:,ll));
        end
    end
end

imgBase = zeros(N,N,2);
for ii = 2:3
    dataBase = zeros(N,N);
    dataBase((N-128)/2+1:(N-128)/2+128,:) = data1(:,:,ii);
    imgBase(:,:,ii-1) = ift2(dataBase);
end
save(['baseImgs_' mfilename],'imgSep','imgBase');

addpath ../../Nifti_Analyze/

brainMask = zeros(N,N,nSlices);
for ii = 1:nSlices
    imgwr = make_nii(abs(imgSep(:,:,ii,end)));
    save_nii(imgwr,'imgBet');
    unix('bet imgBet imgBetMask -m');
    unix('gunzip -f imgBetMask_mask.nii.gz');
    maskStruct = load_untouch_nii('imgBetMask_mask.nii');
    brainMask(:,:,ii) = logical(maskStruct.img);
end
brainMask = brainMask & abs(imgSep(:,:,:,end)) < 0.4;
brainMask = brainMask & abs(imgSep(:,:,:,end)) > 0.1;

% zero pad the data to square
datazp = zeros(N,N,size(data1,3));
datazp((N-128)/2+1:(N-128)/2+128,:,:) = data1;
data1 = datazp; clear datazp

% get kWts
kWts = A(pulse_index,:); % active pulses during dynamics are same as first baseline
kWts = [zeros((N-128)/2,nSlices);kWts;zeros((N-128)/2,nSlices)];
kWts = kron(ones(N,1),kWts);

% zet up precision calculation
L = permute(imgSep*N*N,[1 2 4 3]);
ct = -1/(3*2*pi*42.58*0.01*0.0139);
tempSMS_cd = zeros(N,N,2,size(data1,3)-nSlices*nDummies);
parfor ii = nSlices*nDummies+1:size(data1,3)

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
tempSMS_cd = tempSMS_cd(:,:,:,nSlices*nDummies+1:end);

save(['results_' mfilename],'tempSMS_cd');
