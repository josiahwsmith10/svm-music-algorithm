%% MVDR DoA Algorithm

function Pmvdr = mvdrDoA(Rx,d,lambda)
scanAng = linspace(-90,90,256);
Pmvdr = zeros(size(scanAng));

for iAng = 1:length(scanAng)
    ek = exp(1j*2*pi*d/lambda*sin(pi/180*scanAng(iAng)));
    Pmvdr(iAng) = 1/(ek'*pinv(Rx)*ek);
end

Pmvdr = mag2db(abs(Pmvdr));
Pmvdr = Pmvdr - max(Pmvdr);