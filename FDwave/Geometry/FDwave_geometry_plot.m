function FDwave_geometry_plot( varargin)%figno, hop_s, hop_r)
% This function plots the source and reciever with/without model.
% Complete Syntax:
%          geometry_plot('FIGNO',value,'HOP_S',value,'HOP_R',value)
% Description of parameters:
%                FIGNO    :  Plot in figure number given by.
%                HOP_S    :  No of receiver to be skipped
%                HOP_R    :  No of receiver to be skipped
% Note: 
% 1) Don't place source very close to surface.
% 2) Give all the position after adding/subtracting the ABC layer thickness according to side where they are added.
% 3) All parameters are mandatory.
%       In case any parameter is not provided the program will try to use the values stored in input folder. 
%       If it cannot find any stored value then it shows error.
% Example:
%       geometry_plot('FIGNO',3,'HOP_S',1,'HOP_R',5);
%       geometry_plot()

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'figno';      figno=varargin{i+1};
        case 'wfp';        wfp=varargin{i+1};            
        case 'hop_s';      hop_s=varargin{i+1};
        case 'hop_r';      hop_r=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

if exist('figno','var');         figure(figno);         
elseif ~exist('figno','var');         figure();         
end

if ~exist('wfp','var');         wfp=pwd;        end;
if ~exist('hop_s','var');       hop_s=1;        end
if ~exist('hop_r','var');       hop_r=4;        end


load([wfp,[filesep,'model']],'vpm','nh','nv','dh','dv')
load([wfp,[filesep,'geometry_rec']],'rnh_vec','rnv_vec');
load([wfp,[filesep,'geometry_src']],'snh_vec','snv_vec');
load([wfp,[filesep,'BC']],'BC')

if nh<=nv
    i1=1;  i2=2;
else 
    i1=2; i2=1;
end

h=subplot(i1,i2,1);   imagesc(0:nh-1,0:nv-1,(BC));     
hold on;  colormap('gray')
plot(rnh_vec(1:hop_r:end),rnv_vec(1:hop_r:end),'bv' ); 
plot(snh_vec(1:hop_s:end),snv_vec(1:hop_s:end),'r*' )
title('Source and Reciever with Boundaries Only');
axis ij;


if verLessThan('matlab','8.4')
   xtk= dh*get(h,'XTick');      set(h,'XTickLabel',xtk);
   ytk= dv*get(h,'YTick');      set(h,'YTickLabel',ytk); 
   ax = gca;                    set(ax,'UserData',[dh,dv]);
   h=gcf;                       set(h,'ResizeFcn','axTickUpdate');
else
   h.XTickLabel=dh*h.XTick;            
   h.YTickLabel=dv*h.YTick;
   ax = gca;
   ax.UserData.dh = dh;         
   ax.UserData.dv = dv;
   h=gcf;     
   h.SizeChangedFcn = 'axTickUpdate';
end

vpm=vpm/max(max(vpm));
test= BC+vpm;
h=subplot(i1,i2,2);   imagesc(0:nh-1,  0:nv-1,test); hold on;  colormap('gray')
plot(rnh_vec(1:hop_r:end),rnv_vec(1:hop_r:end),'bv' ); 
plot(snh_vec(1:hop_s:end),snv_vec(1:hop_s:end),'r*' )
title('Recievers plotted on Boundaries & Model')
xlabel('X (m)');    ylabel('Z(m)')
axis ij;

if verLessThan('matlab','8.4')
   xtk= dh*get(h,'XTick');      set(h,'XTickLabel',xtk);
   ytk= dv*get(h,'YTick');      set(h,'YTickLabel',ytk); 
   ax = gca;                    set(ax,'UserData',[dh,dv]);
   h=gcf;                       set(h,'ResizeFcn','axTickUpdate');
else
   h.XTickLabel=dh*h.XTick;            
   h.YTickLabel=dv*h.YTick;
   ax = gca;
   ax.UserData.dh = dh;         
   ax.UserData.dv = dv;
   h=gcf;     
   h.SizeChangedFcn = 'axTickUpdate';
end


end

