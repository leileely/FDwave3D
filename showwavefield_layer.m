% show wavefield of the layered model
% reproducing Figure 3(a)-(c) of the related manuscript

%% synthetic wavefields
load wavefield_layer3_iso;
waveISO=reshape(wavefield(21:271,21:321,46),251,301);
load wavefield_layer3_vti;
waveVTI=reshape(wavefield(21:271,21:321,46),251,301);
load wavefield_layer3_hti;
waveHTI=reshape(wavefield(21:271,21:321,46),251,301);

dx=10;
dt=0.0005;

figure
h1=subplot(131);
imagesc(waveISO/max(max(waveISO)));
colormap(gray);
caxis([-0.02 0.02]);
%xlabel('X (m)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(a)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1:50:301],'xticklabel',[],'fontsize',16,'fontname','Arial')


h2=subplot(132);
imagesc(waveVTI/max(max(waveVTI)));
colormap(gray);
caxis([-0.02 0.02]);
%xlabel('X (m)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(b)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1:50:301],'xticklabel',[],'fontsize',16,'fontname','Arial')


h3=subplot(133);
imagesc(waveHTI/max(max(waveHTI)));
colormap(gray);
caxis([-0.02 0.02]);
xlabel('X (m)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(c)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1:50:301],'xticklabel',[0:50:300]*dx,'fontsize',16,'fontname','Arial')



set(h1,'pos',[0.2 0.7 0.7 0.28]);
set(h2,'pos',[0.2 0.38 0.7 0.28]);
set(h3,'pos',[0.2 0.06 0.7 0.28]);

set(gcf,'pos',[100 100 550 1400])

