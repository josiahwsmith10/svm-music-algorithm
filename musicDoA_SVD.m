%% My music DoA SVD algorithm (Worse than EVD Method)

function [doas,Pmusic] = musicDoA_SVD(Rx,nSig,d,lambda)
[U,~,V] = svd(sqrtm(Rx));
scanAng = linspace(-90,90,256);
Pmusic = zeros(size(scanAng));

Un = U(:,1+nSig:end);
Vn = V(1+nSig:end,:);
for iAng = 1:length(scanAng)
    a = exp(1j*2*pi*d/lambda*sin(pi/180*scanAng(iAng)));
    Pmusic(iAng) = 1 / ( a'*Un*Vn*a );
end
figure;plot(scanAng,abs(Pmusic))
title("Music DoA Pseudo-Spectrum")
xlabel("Angle (degrees)")
xlim([-90,90])
[~,idoas] = maxk(abs(Pmusic),nSig);
doas = sort(scanAng(idoas));
end