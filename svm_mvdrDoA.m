%% SVM-MVDR DoA Algorithm

function Psvm_mvdr = svm_mvdrDoA(R,antLoc,lambda)

nCons = 10; % J - Number of Constraints
scanAng_deg = linspace(-90,90,256);

Rreg = R + 1e-6*eye(size(R,1));
invR = pinv(Rreg);
idxStep = 30;

Psvm_mvdr = zeros(length(scanAng_deg),1);

for iAng = 1:length(scanAng_deg)  % k - Frequency index
    
    uk_correct = exp(1j*2*pi*antLoc/lambda*sin(pi/180*scanAng_deg(iAng)));
    
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
    
    H = -(Uk'*invR*Uk + gamma*eye(nCons));
    
    H = -(H+H')/2;
    
    f = rk;
    
    OPTIONS = optimset('LargeScale','off',...
        'diffmaxchange',1e-4,'Diagnostics','off','MaxIter',5e10,'Display','off');
    
    Aeq = Uk'*invR*Uk;
    phik = quadprog(real(H),f,[],[],real(Aeq),rk,[],[],[],OPTIONS);
    
    wk = invR*Uk*phik;
    
    Psvm_mvdr(iAng) = mag2db(abs(wk'*R*wk));
end

Psvm_mvdr = Psvm_mvdr - max(Psvm_mvdr);