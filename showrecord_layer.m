% show original seismograms of the layered model
% reproducing Figure 3(d)-(f) of the related manuscript

%% synthetic seismograms
load layer3_iso_dc1;
recISO=Sy;
load layer3_vti_dc1;
recVTI=Sy;
load layer3_hti_dc1;
recHTI=Sy;

dx=10;
dt=0.0005;

figure
h1=subplot(131);
imagesc(recISO'/max(max(recISO')));
colormap(gray);
caxis([-0.08 0.1]);
%xlabel('Time (s)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(a)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1,400:400:2400],'xticklabel',[],'fontsize',16,'fontname','Arial')


h2=subplot(132);
imagesc(recVTI'/max(max(recVTI')));
colormap(gray);
caxis([-0.08 0.1]);
%xlabel('Time (s)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(b)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1,400:400:2400],'xticklabel',[],'fontsize',16,'fontname','Arial')


h3=subplot(133);
imagesc(recHTI'/max(max(recHTI')));
colormap(gray);
caxis([-0.08 0.1]);
xlabel('Time (s)','fontsize',16);
ylabel('Z (m)','fontsize',16);
title('(c)','fontweight','normal');
set(gca,'ytick',[1:50:251],'yticklabel',[0:50:250]*dx,'xtick',[1,400:400:2400],'xticklabel',[0,400:400:2400]*dt,'fontsize',16,'fontname','Arial')



set(h1,'pos',[0.2 0.7 0.68 0.28]);
set(h2,'pos',[0.2 0.38 0.68 0.28]);
set(h3,'pos',[0.2 0.06 0.68 0.28]);

set(gcf,'pos',[100 100 500 1400])

