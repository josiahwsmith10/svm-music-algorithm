function Psvm_music = svm_musicDoA(R,nSig,antLoc,lambda)

% Decare Uk and rk

nCons = 10; % J - Number of Constraints
scanAng_deg = linspace(-90,90,256);
idxStep = 30;

Psvm_music = zeros(length(scanAng_deg),1);

[Evec,Eval] = eig(R,'vector');
[~,iEvec] = sort(Eval);
Evec = Evec(:,flip(iEvec));
Vn = Evec(:,1+nSig:end);

for iAng = 1:length(scanAng_deg)  % k - Frequency index
    
    uk_correct = exp(1j*2*pi*antLoc/lambda*sin(pi/180*scanAng_deg(iAng)));
    
    %     idxStep = 72;
    idx = [(iAng-idxStep*nCons/2):idxStep:(iAng-idxStep),(iAng+idxStep):idxStep:(iAng-idxStep+idxStep*nCons/2)];
    while min(idx) < 1 || max(idx) > length(scanAng_deg)
        idx(idx < 1) = idx(idx < 1) + length(scanAng_deg);
        idx(idx > length(scanAng_deg)) = idx(idx > length(scanAng_deg)) - length(scanAng_deg);
    end
    if ~isempty(find(idx == iAng, 1))
        warning("Repeated index error!!")
        return;
    end
    
    uk_wrong = exp(1j*2*pi*antLoc/lambda.*sin(pi/180*scanAng_deg(idx)));
    
    Uk = cat(2,uk_correct,uk_wrong);
    
    rk = zeros(nCons,1);
    rk(1) = 1;
    
    % Declare H and f
    
    gamma = 0.01;
    
    H = -(Uk'*(Vn*Vn')*Uk + gamma*eye(nCons));
    
    H = (H+H')/2;
    
    f = -rk;
    
    OPTIONS = optimset('LargeScale','off',...
        'diffmaxchange',1e-4,'Diagnostics','off','MaxIter',5e5,'Display','off');
    
    Aeq = real(Uk'*(Vn*Vn')*Uk);
    phik = quadprog(real(H),f,[],[],Aeq,rk,[],[],[],OPTIONS);
    
    
    wk = (Vn*Vn')*Uk*phik;
    
    Psvm_music(iAng) = mag2db(abs(wk'*R*wk));
end

Psvm_music = Psvm_music - max(Psvm_music);