function model_2Dplot(varargin)       %  wave_type,  dh,dv,nh,nv,   vpm,rhom,vsm,qpm,qsm)
% This function plots 2D model
for i=1:2:length(varargin)
    switch varargin{i}
        case 'wave_type';       wave_type=varargin{i+1};
        case 'wfp';             wfp=varargin{i+1};
    end
end

if ~exist('wfp','var');              wfp=pwd;                   end;    
  
str=strcat(wfp,[filesep,'model']);
load(str);
vpm2=reshape(vpm(:,:,round(ny/2)),nz,nx);
vsm2=reshape(vsm(:,:,round(ny/2)),nz,nx);
rhom2=reshape(rhom(:,:,round(ny/2)),nz,nx);

figure()

if strcmp(wave_type,'acoustic1') % acoustic/scaler wave
    plotmat2(1,2,1,dx,dz,nx,nz,vpm2,'V_p')
    %     subplot(1,2,1); imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    
elseif strcmp(wave_type,'acoustic2') % acoustic wave
    plotmat2(2,2,1,dx,dz,nx,nz,vpm2,'V_p')
    plotmat2(2,2,2,dx,dz,nx,nz,rhom2,'Density')
    %     subplot(2,2,1); imagesc(dh*(0:nh-1),dv*(0:nv-1),vpm);axis ij; colormap(flipud(gray)); colorbar;title('Vp');set(gca,'XAxisLocation','top')
    %     subplot(2,2,2); imagesc(dh*(0:nh-1),dv*(0:nv-1),rhom);axis ij; colormap(flipud(gray)); colorbar;title('Rho');set(gca,'XAxisLocation','top')
    
elseif strcmp(wave_type,'elastic')  % elastic wave
    plotmat2(2,2,1,dx,dz,nx,nz,vpm2,'V_p')
    plotmat2(2,2,2,dx,dz,nx,nz,vsm2,'V_s')
    plotmat2(2,2,3,dx,dz,nx,nz,rhom2,'\rho')
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


