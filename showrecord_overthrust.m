% show original seismograms of the overthrust model
% reproducing Figure 5(d)-(f) of the related manuscript

%% synthetic seismograms
load overthrust_iso_dc1;
recISO=Sy;
load overthrust_vti_dc1;
recVTI=Sy;
load overthrust_hti_dc1;
recHTI=Sy;

dx=10;
dt=0.0005;

figure
h1=subplot(131);
imagesc(recISO/max(max(recISO)));
colormap(gray);
caxis([-0.08 0.1]);
xlabel('Receiver No.','fontsize',16);
ylabel('Time (s)','fontsize',16);
title('(d)','fontweight','normal');
set(gca,'xtick',[1:50:201],'ytick',[1:400:2000,2000],'yticklabel',[0:400:2000,2000]*dt,'fontsize',16,'fontname','Arial')


h2=subplot(132);
imagesc(recVTI/max(max(recVTI)));
colormap(gray);
caxis([-0.08 0.1]);
xlabel('Receiver No.','fontsize',16);
ylabel('Time (s)','fontsize',16);
title('(e)','fontweight','normal');
set(gca,'xtick',[1:50:201],'ytick',[1:400:2000,2000],'yticklabel',[0:400:2000,2000]*dt,'fontsize',16,'fontname','Arial')


h3=subplot(133);
imagesc(recHTI/max(max(recHTI)));
colormap(gray);
caxis([-0.08 0.1]);
xlabel('Receiver No.','fontsize',16);
ylabel('Time (s)','fontsize',16);
title('(f)','fontweight','normal');
set(gca,'xtick',[1:50:201],'ytick',[1:400:2000,2000],'yticklabel',[0:400:2000,2000]*dt,'fontsize',16,'fontname','Arial')



set(h1,'pos',[0.06 0.15 0.25 0.8]);
set(h2,'pos',[0.4 0.15 0.25 0.8]);
set(h3,'pos',[0.74 0.15 0.25 0.8]);

set(gcf,'pos',[100 100 1400 500])

