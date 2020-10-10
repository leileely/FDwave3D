% show wavefield of overthrust model
% reproducing Figure 5(a)-(c) of the related manuscript

%% synthetic wavefields
load wavefield_overthrust_iso;
waveISO=reshape(wavefield(:,:,2),94,401);
load wavefield_overthrust_vti;
waveVTI=reshape(wavefield(:,:,2),94,401);
load wavefield_overthrust_hti;
waveHTI=reshape(wavefield(:,:,2),94,401);

dx=10;
dt=0.0005;

figure
h1=subplot(131);
imagesc(waveISO/max(max(waveISO)));
colormap(gray);
caxis([-0.5 0.6]);
xlabel('X (km)','fontsize',16);
ylabel('Z (km)','fontsize',16);
title('(a)','fontweight','normal');
set(gca,'xtick',[1:100:401],'xticklabel',[0:100:400]*dx/1000,'ytick',[1:25:94],'yticklabel',[0:25:93]*dx/1000,'fontsize',16,'fontname','Arial')


h2=subplot(132);
imagesc(waveVTI/max(max(waveVTI)));
colormap(gray);
caxis([-0.5 0.6]);
xlabel('X (km)','fontsize',16);
ylabel('Z (km)','fontsize',16);
title('(b)','fontweight','normal');
set(gca,'xtick',[1:100:401],'xticklabel',[0:100:400]*dx/1000,'ytick',[1:25:94],'yticklabel',[0:25:93]*dx/1000,'fontsize',16,'fontname','Arial')


h3=subplot(133);
imagesc(waveHTI/max(max(waveHTI)));
colormap(gray);
caxis([-0.5 0.6]);
xlabel('X (km)','fontsize',16);
ylabel('Z (km)','fontsize',16);
title('(c)','fontweight','normal');
set(gca,'xtick',[1:100:401],'xticklabel',[0:100:400]*dx/1000,'ytick',[1:25:94],'yticklabel',[0:25:93]*dx/1000,'fontsize',16,'fontname','Arial')



set(h1,'pos',[0.1 0.74 0.8 0.2]);
set(h2,'pos',[0.1 0.42 0.8 0.2]);
set(h3,'pos',[0.1 0.1 0.8 0.2]);

set(gcf,'pos',[100 100 1000 900])

