% show overthrust model
% reproducing Figure 4 of the related manuscript

load overthrust_3d_vp;
data=permute(data,[1,3,2]);
dx=5;

%% in 2D slices
figure;

h2=subplot(121);
imagesc(reshape(data(:,:,401),187,801))
hold on;
plot(401,94,'rp','markerface','r','markersize',10);
hold on;
plot(20:40:800,5,'yv','markerface','y','markersize',10);
hold on;
rectangle('pos',[300 50 200 100],'linewidth',2);
xlabel('X (km)','fontsize',16);
ylabel('Z (km)','fontsize',16);
title('(b)','fontweight','normal');
set(gca,'xtick',[1:200:801],'xticklabel',[0:200:800]*dx/1000,'ytick',[1:50:187],'yticklabel',[0:50:186]*dx/1000,'fontsize',16,'fontname','Arial')
colormap(jet);

h3=subplot(122);
imagesc(reshape(data(:,401,:),187,801))
hold on;
plot(401,94,'rp','markerface','r','markersize',10);
hold on;
plot(401,5,'yv','markerface','y','markersize',10);
hold on;
rectangle('pos',[300 50 200 100],'linewidth',2);
xlabel('Y (km)','fontsize',16);
ylabel('Z (km)','fontsize',16);
title('(c)','fontweight','normal');
set(gca,'xtick',[1:200:801],'xticklabel',[0:200:800]*dx/1000,'ytick',[1:50:187],'yticklabel',[0:50:186]*dx/1000,'fontsize',16,'fontname','Arial')
colormap(jet);

set(h2,'pos',[0.15 0.65 0.8 0.3]);
set(h3,'pos',[0.15 0.15 0.8 0.3]);

set(gcf,'pos',[100 100 600 500])

%% in 3D slices
figure
dx=10;
data=permute(data,[3,2,1]);
xslice=[1:1:1];
yslice=[1:1:1];
zslice=[1:1:1];
[x,y,z] = meshgrid(dx*(0:400),dx*(0:400),dx*(0:93));
h1=slice(x,y,z,data(1:2:end,1:2:end,1:2:end),xslice,yslice,zslice);
shading interp;
set(gca,'fontsize',16,'xdir','reverse','ydir','reverse','zdir','reverse');

xlabel('X (km)','fontsize',16);
ylabel('Y (km)','fontsize',16);
zlabel('Z (km)','fontsize',16);
title('(a)','fontweight','normal');
set(gca,'xtick',[1:100:401]*dx,'xticklabel',[0:100:400]*dx/1000,'ytick',[1:100:401]*dx,'yticklabel',...
    [0:100:400]*dx/1000,'ztick',[1:25:94]*dx,'zticklabel',[0:25:93]*dx/1000,'fontsize',16,'fontname','Arial')
axis([0 4000 0 4000 0 930])
view(120,18)
colormap(jet);

set(gcf,'pos',[100 100 600 500])
ch = colorbar('horiz');%
set(get(ch,'title'),'string','km/s','fontsize',16);
set(ch,'ticks',(3000:1000:6000),'fontsize',14,'ticklabels',(3:6));

