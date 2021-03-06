function geometry_3Dplot_recxz(varargin)%figno, hop, rnx_vec, rny_vec, rnz_vec,plots)
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

if verLessThan('matlab','8.4')
    if strcmp(BCPlot,'y')
        h=subplot(i1,i2,i3);  
        if strcmp(BCname,'PML')
        imagesc(0:nx-1,  0:nz-1,BC.pmlx(:,:,round(ny/2))+BC.pmly(:,:,round(ny/2))+BC.pmlz(:,:,round(ny/2)));hold on;
        else
            imagesc(0:nx-1,  0:nz-1,BC.abl(:,:,round(ny/2)));hold on;
        end
        plot(rnx_vec(1:hop:end),rnz_vec(1:hop:end),'bv' ); colormap('gray')
        title('Receiver with Boundaries Only')
        xlabel('X (m)');    ylabel('Z (m)')
        
        xtk = dx*get(h,'XTick');        set(h,'XTickLabel',xtk)
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
        h=subplot(i1,i2,i3+1);   imagesc(0:nx-1,  0:nz-1,test(:,:,round(ny/2))); hold on;  colormap('gray')
        plot(rnx_vec(1:hop:end), rnz_vec(1:hop:end),'bv' );
        title('Receivers plotted on Boundaries & Model')
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
        plot(rnx_vec(1:hop:end),rnz_vec(1:hop:end),'bv' ); colormap('gray')
        title('Receiver with Boundaries Only')
        xlabel('X (m)');    ylabel('Z (m)')
        h.XTickLabel = dx*h.XTick;
        h.YTickLabel = dz*h.YTick;
    end
    
    if strcmp(ModelPlot,'y')
        vpm=vpm/max(vpm(:));
        if strcmp(BCname,'PML')
            test= BC.pmlx+BC.pmly+BC.pmlz+vpm;
        else
            test= BC.abl+vpm;
        end
        h=subplot(i1,i2,i3+1);   imagesc(0:nx-1,  0:nz-1,test(:,:,round(ny/2))); hold on;  colormap('gray')
        plot(rnx_vec(1:hop:end), rnz_vec(1:hop:end),'bv' );
        title('Receivers plotted on Boundaries & Model')
        xlabel('X (m)');    ylabel('Z (m)')
        h.XTickLabel = dx*h.XTick;
        h.YTickLabel = dz*h.YTick;
    end
end

end

