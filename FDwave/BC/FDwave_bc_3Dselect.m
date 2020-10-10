function FDwave_bc_3Dselect(varargin)     %BCname,BCtype,nAB  )
% This function generates the necessary data for different BC
% Complete Syntax:
%       select_bc('WFP',path,'BCNAME',value,'BCTYPE',value,'NAB',value,'PlotON',option)
% Description of parameters:
%       WFP    :  Path to working directory
%       BCNAME = 'ABC1'   for Englist & Clayton EnglistClayton (acoustic wave eqn)
%                'ABC2'   for Raynold EnglistClayton (acoustic wave eqn)
%                'ABL'    for Damping EnglistClayton, after Cerjan
%                'PML'    for Perfectly Matched Layer, after Collino
%       BCTYPE = 'topFS'    for free surface at top
%              = 'topABC'   for absorbing type boundary at surface
%       NAB    =  no of layers to be used (for damping BC only)
%       PlotON = 'y',  To plot the boundaries
% Example:
%       select_bc('WFP',pwd,'BCNAME','ABL','BCTYPE','topABC','NAB',40)


for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';        wfp=varargin{i+1};
        case 'bcname';     BCname=varargin{i+1};
        case 'bctype';     BCtype=varargin{i+1};
        case 'nab';        nAB=varargin{i+1};
        case 'verbose';    verbose=varargin{i+1};
        case 'ploton';     plotON=varargin{i+1};
    end
end

str0='        ';
str1=str0;          str2=str0;      str3=str0;
if ~exist('wfp','var');           wfp=pwd;              end;
if ~exist('BCname','var');        BCname='ABL';         str1=[str0,'(Default)'];    end
if ~exist('nAB','var');           nAB= 40;              str2=[str0,'(Default)'];    end
if ~exist('BCtype','var');        BCtype='topABC';      str3=[str0,'(Default)'];    end
if ~exist('verbose','var');       verbose='n'  ;       end
if ~exist('plotON','var');        plotON='n';          end
if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
    disp('    FUNC: Select the boundary types');
end

if strcmp(BCname,'ABC1')||strcmp(BCname,'ABC2');           nAB = 2;         end


if strcmp(BCname,'ABL')
    if strcmp(verbose,'y')
        disp([str0,'BC Method    :  Damping Layer (Cerjan)']);
    end
    BC = bc_3Ddamp(wfp, BCtype,nAB );  % options: topFS, allABC

elseif strcmp(BCname,'PML')
    if strcmp(verbose,'y')
        disp([str0,'BC Method    :  Perfectly Matched Layer (Collino)']);
    end
    BC = bc_3Dpml(wfp, BCtype,nAB ); 
    
elseif strcmp(BCname,'ABC1')
    disp([str0,'BC Method    :  Absorbing BC (EnglistClayton)']);
    load('model','nx','ny','nz');
    BC = ones(nz,nx,ny);
    tmp=ceil(.01*nz);
    if strcmp(BCtype,'topFS')
        BC(1:nz-tmp,1+tmp:nx-tmp,1+tmp:ny-tmp)= 0;
    elseif strcmp(BCtype,'topABC')
        BC(1+tmp:nz-tmp,1+tmp:nx-tmp,1+tmp:ny-tmp)= 0;
    end
    
elseif strcmp(BCname,'ABC2')
    disp([str0,'BC Method    :  TransparantBC(Cerjan)']);
    load('model','nx','ny','nz');
    BC = ones(nz,nx,ny);
    tmp=ceil(.01*nz);
    if strcmp(BCtype,'topFS')
        BC(1:nz-tmp,1+tmp:nx-tmp,1+tmp:ny-tmp)= 0;
    elseif strcmp(BCtype,'topABC')
        BC(1+tmp:nz-tmp,1+tmp:nx-tmp,1+tmp:ny-tmp)= 0;
    end
end

str=strcat(wfp,[filesep,'BC']);
save(str,'BC','nAB','BCtype','BCname');

if strcmp(verbose,'y')
    disp([str0,'No of Layers :  ', num2str(nAB),str1]);
    disp([str0,'BC type      :  ', BCtype,str2]);
    disp([str0,'BC name      :  ',BCname,str3])
    disp(['       Boundaries data is saved in "',str,'"'])
end


if strcmp(plotON,'y2d')
    
    load([wfp,[filesep,'model']],'dx','dz')
    if strcmp(BCname,'PML')
        figure(); [nz,nx,ny]=size(BC.pmlx); plotmat2(1,1,1,dx,dz,nx,nz,BC.pmlx(:,:,round(ny/2))+BC.pmly(:,:,round(ny/2))+BC.pmlz(:,:,round(ny/2)),'Boundaries');
        % reshape(BC.pmly(:,round(nx/2),:),nz,ny)+reshape(BC.pmlz(:,round(nx/2),:),nz,ny)+reshape(BC.pmlx(:,round(nx/2),:),nz,ny)
    else
        figure(); [nz,nx,ny]=size(BC.abl); plotmat2(1,1,1,dx,dz,nx,nz,BC.abl(:,:,round(ny/2)),'Boundaries')
    end
    axis image;  colormap(flipud(jet))
    
    else if strcmp(plotON,'y3d')
            load([wfp,[filesep,'model']],'nx','ny','nz','dx','dy','dz')
            xslice=round(nx/2)*dx;
            yslice=round(ny/2)*dy;
            zslice=round(nz/2)*dz;
            [x,y,z] = meshgrid(dx*(0:nx-1),dy*(0:ny-1),dz*(0:nz-1));
            if strcmp(BCname,'PML')
                figure(); 
                bcpml=permute(BC.pmlx+BC.pmly+BC.pmlz,[2,3,1]);
                slicen(x,y,z,bcpml,xslice,yslice,zslice);
                shading interp;
                set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
                xlabel('X (m)');ylabel('Y (m)');
                zlabel('Z (m)');title ('Boundaries')
                colormap(flipud(jet))
                axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
                colorbar;set(gca,'fontsize',14);
            else
                figure(); 
                bcabl=permute(BC.abl,[2,3,1]);
                slicen(x,y,z,bcabl,xslice,yslice,zslice);
                shading interp;
                set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
                xlabel('X (m)');ylabel('Y (m)');
                zlabel('Z (m)');title ('Boundaries')
                colormap(flipud(jet))
                axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
                colorbar;set(gca,'fontsize',14);
            end
        end
            
            
end

