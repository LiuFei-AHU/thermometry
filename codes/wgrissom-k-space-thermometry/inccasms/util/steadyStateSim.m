% idea: starting with all M along z, simulate excitation and relaxations of a set of isochromats 
% distributed at frequencies 0-1/TR Hz
T1 = 200; % T1 in units of TR - e.g., if tap water has T1 = 5s, and TR = 25 ms
T2 = 0.5*T1;


nTRs = 1024; % number of TRs to simulate through

flip = 30; % flip angle % Zur fig 5 appears to have used about 60

nIso = 200; % # isochromats

frqs = linspace(0,1,nIso); % cycles per TR

Mz = ones(nIso,1);
Mxy = zeros(nIso,1);


alpha = cosd(flip/2);
beta = sind(flip/2);

phsInc = [0 180 90 117];% Zur Fig 5 values

MxySave = zeros(nTRs,length(phsInc));

flippingPattern = 'everyOther' % 'everyOther', 'random', 'noFlipping'
switch flippingPattern
    case 'everyOther'
        flipped = repmat([1;0],[nTRs/2 1]);
    case 'random'
        flipped = rand(128,nTRs/128) > 0.5;
        flipped = flipped(:); 
    case 'noFlipping'
        flipped = zeros(nTRs,1);
end
flippedPhs = pi;

for jj = 1:length(phsInc)
    
    %phsInc = 117*pi/180;
    Mz = ones(nIso,1);
    Mxy = zeros(nIso,1);
    totalPhs = 0;
    for ii = 1:nTRs
        
        % increment RF phase and apply RF
        totalPhs = (ii-1)*phsInc(jj) + totalPhs;
        %totalPhs = rand(1)*360;
        
        betaPhs = beta * exp(1i*(totalPhs*pi/180 + flippedPhs*flipped(mod(ii-1,length(flipped))+1)));
        Mxyt = Mxy * alpha^2 - conj(Mxy)*betaPhs^2 + 2*alpha*betaPhs*Mz;
        Mzt = -alpha*conj(betaPhs)*Mxy - alpha*betaPhs*conj(Mxy) + (abs(alpha)^2 - abs(betaPhs)^2)*Mz;
        Mxy = Mxyt;
        Mz = Mzt;
        
        % store total Mxy
        MxySave(ii,jj) = sum(Mxy)*exp(-1i*(totalPhs*pi/180 + flippedPhs*flipped(mod(ii-1,length(flipped))+1)));
        
        % apply relaxation
        Mxy = Mxy * exp(-1/T2);
        Mz = (1-exp(-1/T1)) + Mz * exp(-1/T1);
        
        % apply off-resonance phase shifts
        Mxy = Mxy.*exp(1i*2*pi*frqs');
        
        %clf;subplot(121);plot(real(Mxy));hold on;plot(imag(Mxy));plot(abs(Mxy));subplot(122);plot(Mz);pause
        
    end
    
end

figure;plot(real(MxySave)/nIso);
axis([1 1024 -0.1 0.4]);

%figure;plot(real(MxySave(end,:))/nIso);

legend('0','180','90','117');

