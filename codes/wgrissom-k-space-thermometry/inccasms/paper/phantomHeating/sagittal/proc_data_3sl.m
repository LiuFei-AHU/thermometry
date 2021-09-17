% water bath
N = 256;
nSlices = 3;
nDummies = 3;
nEncodes = 2;
pFileName = 'P17408.7';
load([pFileName '.mat'])
Nt = size(data1,3);
data1 = data1.*repmat(kron(ones(64,1),[1 -1]'),[1 N Nt]);
data1 = conj(data1);

datazp = zeros(N,N);
datazp((N-128)/2+1:(N-128)/2+128,:) = data1(:,:,20);
dataAlias = datazp; clear datazp
imgAlias = circshift(ift2(dataAlias),[-30 0]);
save imgAlias_3sl imgAlias

% crop it in freq enc dimension
data1 = fftshift(ifft(data1,[],2),2);
data1(:,1:100,:) = 0;
data1(:,140:end,:) = 0;
data1 = fft(fftshift(data1,2),[],2);

% downsample in freq enc dimension
data1 = data1(:,1:2:end,:) + data1(:,2:2:end,:);
N = 128;

phs = pi;
A = [1 1 1;exp(1i*phs/2) 1 exp(-1i*phs/2); exp(1i*phs) 1 exp(-1i*phs)]; % encoding functions

load permutation2
pulse_index_baseline = zeros(nSlices,129); % which pulses are used for each line of each baseline acquisition
pulse_index_baseline(1,:) = kron(ones(1,43),[1 2 3]);
pulse_index_baseline(2,:) = kron(ones(1,43),[2 3 1]);
pulse_index_baseline(3,:) = kron(ones(1,43),[3 1 2]);
pulse_index = mod(permutation-1,nEncodes)+1;

imgSep = zeros(N,N,nSlices,nDummies^3);
for ii = 1:nDummies
    for jj = 1:nDummies
        for kk = 1:nDummies
            dataSep = zeros(128,nSlices,N);
            for ll = 1:128
                pulses = pulse_index_baseline(:,ll);
                % what was the pulse order for this
                At = A(pulses,:); % get the active A matrix for this phase encode
                % collect the data into a vector
                datat = squeeze(data1(ll,:,[ii nDummies+jj 2*nDummies+kk])).'; % put the third index first
                % invert A to get the separated slice data
                dataSep(ll,:,:) = At\datat;
            end
            dataSep = permute(dataSep,[1 3 2]);

            for ll = 1:nSlices
                imgSep(:,:,ll,(ii-1)*nDummies^2+(jj-1)*nDummies+kk) = ift2(dataSep(:,:,ll));
            end
        end
    end
end

% remove the ones involving the first image
imgSep = imgSep(:,:,:,10:end);

imgSepForDisplay = imgSep(:,:,:,end);
save imgSepForDisplay imgSepForDisplay

% get kWts
kWts = A(pulse_index,:); % active pulses during dynamics are same as first baseline
kWts = kron(ones(N,1),kWts);

% zet up precision calculation
L = permute(imgSep*N*N,[1 2 4 3]);
ct = -1/(3*2*pi*42.58*0.01*0.0139);

load brainMask_3sl

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

    algp = struct('dofigs',0);
    algp.lam = 10^-2;

    [thetakacc,~,~,~,Ac] = kspace_hybrid_thermo_inccasms(acqp,thetainit,algp);
    tempSMS_cd(:,:,:,ii) = ct*real(thetakacc);

end
tempSMS_cd = tempSMS_cd(:,:,:,10:end);

save(['results_' mfilename],'tempSMS_cd');
