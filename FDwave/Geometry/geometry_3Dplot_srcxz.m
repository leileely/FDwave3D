function geometry_3Dplot_srcxz( varargin)%figno, hop, snx_vec, sny_vec, snz_vec )
% This function plots the source over the with/without model.
%   Complete syntax:
%      geometry_plot_src('figno',value,'hop',value,'snh_vec',vector,'snv_vec',vector)
%   Description of parameters:
%      FIGNO    :  Figure no in which to be plotted.
%                  If blank then plotted in new figure.
%      HOP      :  No of object (sources) to be skipped
%      snx_vec  :  A vector containing X location of all sources
%      sny_vec  :  A vector containing Y location of all sources
%      snz_vec  :  A vector containing Z location of all sources
%      BCPLOT   :  To plot the source with boundries
%      MODELPLOT:  To plot the source with boundaries and model
%   Example:
%       nvpair_geometry_plot_src('figno',3,'hop',5)

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'figno';       figno=varargin{i+1};
        case 'wfp';         wfp=varargin{i+1};
        case 'hop';         hop=varargin{i+1};
        case 'snx_vec';     snx_vec=varargin{i+1};
        case 'sny_vec';     sny_vec=varargin{i+1};
        case 'snz_vec';     snz_vec=varargin{i+1};
        case 'bcplot';      BCPlot=varargin{i+1};
        case 'modelplot';   ModelPlot=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end
if ~exist('figno','var');        hfig=figure();
else hfig=figure(figno);
end

if ~exist('wfp','var');        wfp=pwd;   end;
if ~exist('hop','var');          hop=1;                           end
if ~exist('snx_vec','var');      load([wfp,[filesep,'geometry_src']],'snx_vec');     end
if ~exist('sny_vec','var');      load([wfp,[filesep,'geometry_src']],'sny_vec');     end
if ~exist('snz_vec','var');      load([wfp,[filesep,'geometry_src']],'snz_vec');     end
if ~exist('BCPlot','var');       BCPlot='n';                           end
if ~exist('ModelPlot','var');    ModelPlot='y';                           end

load([wfp,[filesep,'model']],'vpm','nx','ny','nz','dx','dy','dz')
load([wfp,[filesep,'BC']],'BC','BCname','nAB')

if xor(strcmp(BCPlot,'y'),strcmp(ModelPlot,'y'));
   i1=1; i2=1;  i3=1;
elseif strcmp(BCPlot,'y')&&strcmp(ModelPlot,'y');
    if nx>nz;   i1=2; i2=1;  i3=1;    end
    if nx<=nz;  i1=1; i2=2;  i3=1;      end
end


if verLessThan('matlab','8.4')
    if strcmp(BCPlot,'y')
        
        h=subplot(i1,i2,i3);  
        if strcmp(BCname,'PML')
        imagesc(0:nx-1,  0:nz-1,BC.pmlx(:,:,round(ny/2))+BC.pmly(:,:,round(ny/2))+BC.pmlz(:,:,round(ny/2)));hold on;
        else
            imagesc(0:nx-1,  0:nz-1,BC.abl(:,:,round(ny/2)));hold on;
        end
               plot(snx_vec(1:hop:end),snz_vec(1:hop:end),'r*' ); colormap('gray')
               title('Source with Boundaries Only')
               xlabel('X (m)');    ylabel('Z (m)')
        
        xtk = dz*get(h,'XTick');        set(h,'XTickLabel',xtk)
        ytk = dz*get(h,'YTick');        set(h,'YTickLabel',ytk)
        ax = gca;                       set(ax,'UserData',[dx,dz]);
        h=gcf;                          set(h,'ResizeFcn','axTickUpdate');
    end
    
    if strcmp(ModelPlot,'y')
        vpm=vpm/max(vpm(:));
        if strcmp(BCname,'PML')
            test= BC.pmlx+BC.pmly+BC.pmlz+vpm;
        else
            test= BC.abl+vpm;
        end
        h=subplot(i1,i2,i3+1);   imagesc(0:nx-1,  0:nz-1,  test(:,:,round(ny/2))); hold on;  colormap('gray')
        plot(snx_vec(1:hop:end),  snz_vec(1:hop:end),'r*' ); colormap('gray')
        title('Source plotted on Boundaries & Model')
        xlabel('X (m)');    ylabel('Z (m)')
        
        xtk = dx*get(h,'XTick');        set(h,'XTickLabel',xtk)
        ytk = dz*get(h,'YTick');        set(h,'YTickLabel',ytk)
        ax = gca;                       set(ax,'UserData',[dx,dz]);
        h=gcf;                          set(h,'ResizeFcn','axTickUpdate');
    end
else
    if strcmp(BCPlot,'y')
        h=subplot(i1,i2,i3);  
        if strcmp(BCname,'PML')
        imagesc(0:nx-1,  0:nz-1,BC.pmlx(:,:,round(ny/2))+BC.pmly(:,:,round(ny/2))+BC.pmlz(:,:,round(ny/2)));hold on;
        else
            imagesc(0:nx-1,  0:nz-1,BC.abl(:,:,round(ny/2)));hold on;
        end
        plot(snx_vec(1:hop:end),snz_vec(1:hop:end),'r*' ); colormap('gray')
        title('Receiver with Boundaries Only')
        xlabel('X (m)');    ylabel('Z (m)')
        h.XTickLabel = dh*h.XTick;
        h.YTickLabel = dh*h.YTick;
    end
    
    if strcmp(ModelPlot,'y')
        vpm=vpm/max(vpm(:));
        test= BC+vpm;
        h=subplot(i1,i2,i3+1);   imagesc(0:nx-1,  0:nz-1,  test(:,:,round(ny/2))); hold on;  colormap('gray')
        plot(snx_vec(1:hop:end),  snz_vec(1:hop:end),'r*' ); colormap('gray')
        title('Receivers plotted on Boundaries & Model')
        xlabel('X (m)');    ylabel('Z (m)')
        
        h.XTickLabel = dh*h.XTick;
        h.YTickLabel = dh*h.YTick;
    end
end

