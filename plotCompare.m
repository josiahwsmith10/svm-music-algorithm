function plotCompare(Pmvdr,Pmvdrsvm,Pmusic,Pmusicsvm,mindB)

figure;
subplot(121)
plot(linspace(-90,90,256),Pmvdr(:)-max(Pmvdr),'--k');
hold on
plot(linspace(-90,90,256),Pmvdrsvm(:)-max(Pmvdrsvm),'k')
ylabel("dB")
xlabel("Angle (degrees)")
ylim([mindB 0])
xlim([-90 90])
legend("MVDR","SVM-MVDR","Location","southwest")
title("MVDR vs. SVM-MVDR")
subplot(122)
plot(linspace(-90,90,256),Pmusic(:)-max(Pmusic),'--k')
hold on
plot(linspace(-90,90,256),Pmusicsvm(:)-max(Pmusicsvm),'k')
ylabel("dB")
xlabel("Angle (degrees)")
ylim([mindB 0])
xlim([-90 90])
legend("MUSIC","SVM-MUSIC","Location","southwest")
title("MUSIC vs. SVM-MUSIC")