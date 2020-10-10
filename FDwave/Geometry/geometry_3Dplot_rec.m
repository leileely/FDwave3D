function geometry_3Dplot_rec(varargin)%figno, hop, rnx_vec, rny_vec, rnz_vec,plots)
% This function plots the recievers over the with/without model.
% Complete syntax:
%      geometry_plot_rec('WFP',path,'FIGNO',value,'HOP',value,'RNX_VEC',vector,'RNY_VEC',vector,'RNZ_VEC',vector)
%   Description of parameters:
%           WFP      :  Path to working directory
%           FIGNO    :  Figure no in which to be plotted.
%                       If blank then plotted in new figure.
%           HOP      :  No of object (receivers) to be skipped
%           RNX_VEC  :  A vector containing X location of all sources
%           RNY_VEC  :  A vector containing Y location of all sources
%           RNZ_VEC  :  A vector containing Z location of all sources
%   Example:
%           nvpair_geometry_plot_rec('WFP',pwd,'figno',3,'hop',5)

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'figno';       figno=varargin{i+1};
        case 'wfp';         wfp=varargin{i+1};
        case 'hop';         hop=varargin{i+1};
        case 'rnx_vec';     rnx_vec=varargin{i+1};
        case 'rny_vec';     rny_vec=varargin{i+1};
        case 'rnz_vec';     rnz_vec=varargin{i+1};
        case 'bcplot';      BCPlot=varargin{i+1};
        case 'modelplot';   ModelPlot=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

if ~exist('figno','var');        hfig=figure();
else hfig=figure(figno);
end

if ~exist('wfp','var');          wfp=pwd;                           end
if ~exist('hop','var');          hop=1;                             end
if ~exist('rnx_vec','var');      load([wfp,[filesep,'geometry_rec']],'rnx_vec');    end
if ~exist('rny_vec','var');      load([wfp,[filesep,'geometry_rec']],'rny_vec');    end
if ~exist('rnz_vec','var');      load([wfp,[filesep,'geometry_rec']],'rnz_vec');    end
if ~exist('BCPlot','var');       BCPlot='n';                           end
if ~exist('ModelPlot','var');    ModelPlot='y';                           end

load([wfp,[filesep,'model']],'vpm','dx','dy','dz','nx','ny','nz')
load([wfp,[filesep,'BC']],'BC','BCname','nAB')

if xor(strcmp(BCPlot,'y'),strcmp(ModelPlot,'y'));
    i1=1; i2=1;  i3=1;
elseif strcmp(BCPlot,'y')&&strcmp(ModelPlot,'y');
    if nx>nz;   i1=2; i2=1;  i3=1;    end
    if nx<=nz;  i1=1; i2=2;  i3=1;      end
end

load([wfp,[filesep,'model']],'nx','ny','nz','dx','dy','dz')
xslice=dx*(nx-nAB);%round(nx/2);
yslice=dy*nAB;%round(ny/2);
zslice=dz*round(nz/2);
[x,y,z] = meshgrid(dx*(0:nx-1),dy*(0:ny-1),dz*(0:nz-1));


    if strcmp(BCPlot,'y')
        h=subplot(i1,i2,i3); 
        if strcmp(BCname,'PML')
        bcpml=permute(BC.pmlx+BC.pmly+BC.pmlz,[2,3,1]);
        slicen(x,y,z,bcpml,xslice,yslice,zslice);
        else
            bcabl=permute(BC.abl,[2,3,1]);
            slicen(x,y,z,bcabl,xslice,yslice,zslice); 
        end
            
        shading interp;
        set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
        xlabel('X (m)');ylabel('Y (m)');
        zlabel('Z (m)');title ('Boundaries')
        axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
        colorbar;set(gca,'fontsize',14);
        hold on;
        plot3(dx*rnx_vec(1:hop:end),dy*rny_vec(1:hop:end),dz*rnz_vec(1:hop:end),'bv' ); colormap(flipud(gray))
        title('Receivers with Boundaries Only')
        
%         xtk = dx*get(h,'XTick');        set(h,'XTickLabel',xtk)
%         ytk = dy*get(h,'YTick');        set(h,'YTickLabel',ytk)
%         ztk = dz*get(h,'ZTick');        set(h,'ZTickLabel',ztk)
%         
%         ax = gca;                       set(ax,'UserData',[dx,dy,dz]);
%         h=gcf;                          set(h,'ResizeFcn','axTickUpdate');
    end
    
    if strcmp(ModelPlot,'y')
        vpm=vpm/max(vpm(:));
        h=subplot(i1,i2,i3+1);  
        if strcmp(BCname,'PML')
        bcpml=permute(BC.pmlx+BC.pmly+BC.pmlz+vpm,[2,3,1]);
        slicen(x,y,z,bcpml,xslice,yslice,zslice);
        else
            bcabl=permute(BC.abl+vpm,[2,3,1]);
            slicen(x,y,z,bcabl,xslice,yslice,zslice); 
        end
        shading interp;
        set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
        xlabel('X (m)');ylabel('Y (m)');
        zlabel('Z (m)');title ('Boundaries')
        axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
        colorbar;set(gca,'fontsize',14);
        hold on;
        plot3(dx*rnx_vec(1:hop:end),dy*rny_vec(1:hop:end),dz*rnz_vec(1:hop:end),'bv' ); colormap(flipud(gray))
        title('Receivers plotted on Boundaries & Model')
                
%         xtk = dx*get(h,'XTick');        set(h,'XTickLabel',xtk)
%         ytk = dy*get(h,'YTick');        set(h,'YTickLabel',ytk)
%         ztk = dz*get(h,'ZTick');        set(h,'ZTickLabel',ztk)
%         
%         ax = gca;                       set(ax,'UserData',[dx,dy,dz]);
%         h=gcf;                          set(h,'ResizeFcn','axTickUpdate');
    end
    
end

