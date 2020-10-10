function model_3Dplot(varargin)       %  wave_type,  dh,dv,nh,nv,   vpm,rhom,vsm,qpm,qsm)
% This function plots 3D model
for i=1:2:length(varargin)
    switch varargin{i}
        case 'wave_type';       wave_type=varargin{i+1};
        case 'wfp';             wfp=varargin{i+1};
    end
end

if ~exist('wfp','var');              wfp=pwd;                   end 
  
str=strcat(wfp,[filesep,'model']);
load(str);

figure()
xslice=round(nx/2)*dx;
yslice=round(ny/2)*dy;
zslice=round(nz/2)*dz;
[x,y,z] = meshgrid(dx*(0:nx-1),dy*(0:ny-1),dz*(0:nz-1));

if strcmp(wave_type,'acoustic1') % acoustic/scaler wave
    
    vpm=permute(vpm,[2,3,1]);
    slicen(x,y,z,vpm,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)'); title ('Vp (m/s)')
    axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
    colormap (jet); colorbar;set(gca,'fontsize',14);
    
    %     subplot(1,2,1); imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    
elseif strcmp(wave_type,'acoustic2') % acoustic wave
    subplot(2,2,1)
    vpm=permute(vpm,[2,3,1]);
    slicen(x,y,z,vpm,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)');title ('Vp (m/s)')
    axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
    colormap (jet); colorbar;set(gca,'fontsize',14);
    
    subplot(2,2,2)
    rhom=permute(rhom,[2,3,1]);
    slicen(x,y,z,rhom,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)');title ('/rho (kg/m^3)')
    axis([0 dx*(nx-1) 0 dx*(ny-1) 0 dx*(nz-1)]);
    colormap (jet); colorbar;set(gca,'fontsize',14);
    %     subplot(2,2,1); imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    %     subplot(2,2,2); imagesc(dh*(0:nh-1),dv*(0:nv-1),rhom);axis ij; colormap(flipud(gray)); colorbar;title('Rho');set(gca,'XAxisLocation','top')
    
elseif strcmp(wave_type,'elastic')  % elastic wave
    subplot(2,2,1)
    vpm=permute(vpm,[2,3,1]);
    slicen(x,y,z,vpm,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)');title ('Vp (m/s)')
    axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)]);
    colormap (jet); colorbar;set(gca,'fontsize',14);
    
    subplot(2,2,2)
    vsm=permute(vsm,[2,3,1]);
    slicen(x,y,z,vsm,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)');title ('Vs (m/s)')
    axis([0 dx*(nx-1) 0 dx*(ny-1) 0 dx*(nz-1)]);
    colormap (jet); colorbar;set(gca,'fontsize',14);
    
    subplot(2,2,3)
    rhom=permute(rhom,[2,3,1]);
    slicen(x,y,z,rhom,xslice,yslice,zslice);
    shading interp;
    set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
    xlabel('X (m)');ylabel('Y (m)');
    zlabel('Z (m)');title ('/rho (kg/m^3)')
    axis([0 dx*(nx-1) 0 dx*(ny-1) 0 dx*(nz-1)]);
    colormap (jet); colorbar;set(gca,'fontsize',14);
    %subplot(2,2,1);    imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    % subplot(2,2,2);   imagesc(dh*(0:nh-1),dv*(0:nv-1),vsm);axis ij; colormap(flipud(gray)); colorbar;title('Vs');set(gca,'XAxisLocation','top')
    % subplot(2,2,3);   imagesc(dh*(0:nh-1),dv*(0:nv-1),rhom);axis ij; colormap(flipud(gray)); colorbar;title('Rho');set(gca,'XAxisLocation','top')
    
elseif strcmp(wave_type,'viscoelastic') % viscoelastic wave
    plotmat2(2,3,1,dx,dz,nx,nz,vpm2,'V_p')
    plotmat2(2,3,2,dx,dz,nx,nz,vsm2,'V_s')
    plotmat2(2,3,3,dx,dz,nx,nz,rhom2,'Density')
    plotmat2(2,3,4,dx,dz,nx,nz,qpm2,'Q_p')
    plotmat2(2,3,5,dx,dz,nx,nz,qsm2,'Q_s')
    %     subplot(2,3,1); imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    %     subplot(2,3,2); imagesc(dh*(0:nh-1),dv*(0:nv-1),vsm);axis ij; colormap(flipud(gray)); colorbar;title('Vs');set(gca,'XAxisLocation','top')
    %     subplot(2,3,3); imagesc(dh*(0:nh-1),dv*(0:nv-1),rhom);axis ij; colormap(flipud(gray)); colorbar;title('Rho');set(gca,'XAxisLocation','top')
    %     subplot(2,3,4); imagesc(dh*(0:nh-1),dv*(0:nv-1),qpm);axis ij; colormap(flipud(gray)); colorbar;title('Qp');set(gca,'XAxisLocation','top')
    %     subplot(2,3,5); imagesc(dh*(0:nh-1),dv*(0:nv-1),qsm);axis ij; colormap(flipud(gray)); colorbar;title('Qs');set(gca,'XAxisLocation','top')
    
else
    warning('    Wrong model selected')
end


