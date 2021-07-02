%% MUSIC DoA algorithm

function Pmusic = musicDoA(Rx,nSig,antLoc,lambda)
[Evec,Eval] = eig(Rx,'vector');
[~,iEvec] = sort(Eval);
Evec = Evec(:,flip(iEvec));
scanAng = linspace(-90,90,256);
Pmusic = zeros(size(scanAng));

Vn = Evec(:,1+nSig:end);
for iAng = 1:length(scanAng)
    wk = exp(1j*2*pi*antLoc/lambda*sin(pi/180*scanAng(iAng)));
    Pmusic(iAng) = 1 / ( wk'*(Vn*Vn')*wk );
end

Pmusic = mag2db(abs(Pmusic));
Pmusic = Pmusic - max(Pmusic);