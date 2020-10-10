function FDwave_calculation_3Dimage_plot(varargin)%(ipstr,figNO, shotno)
% This function plots the synthetic seimogram with scaling and clipping the amplitudes.
% Complete syntax: 
%   calculation_plot('Wave_Type',option,'ShotGather',matrix,'FigNo',value,...
%        'ShotNo',value,'clip',value,'Scale')
% Description of Parameters
%       ShotGather  :  The synthetic seimogram (Optional)
%       Wave_Type   :  'Acoustic', 'Elastic', etc.
%       FigNo       :  Figure No on which you want to plot
%       ShotNo      :  Source No, if there are more than one (e.g. for an array)
%       Clip        :  To remove very high amplitudes
%       Scale       :  To enhance deeper/later reflections.
% Example
%    calculation_plot('wave_type','elastic','shotno',5)
%    calculation_plot('wave_type','elastic','shotno',5, 'clip',.95,'scale',6)

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                        wfp=varargin{i+1};
        case 'ssfilename';            FileName=varargin{i+1};
        case 'wave_type';           wave_type = lower(varargin{i+1});
        case 'geometry_type';    geometry_type= lower(varargin{i+1});
        case 'dx';                         dx= lower(varargin{i+1});
        case 'dt';                          dt= lower(varargin{i+1});
        case 'figno';                     figNO = varargin{i+1};
        case 'clip';                        clip = varargin{i+1};
        case 'scale';                      scale = varargin{i+1};
        otherwise
           error('%s is not a valid argument name',varargin{i});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('FileName','var');                error('Please enter the file name');         end
if ~exist('geometry_type','var');      geometry_type='surface';     end
if ~exist('dx','var');      load(FileName,'dx');     end
if ~exist('dt','var');      load(FileName,'dt');     end

if ~exist('figNO','var');       h=figure();             end
if  exist('figNO','var');       h=figure(figNO);        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(FileName,'Sz','dN_SS');             
set(0,'CurrentFigure',h);
load([wfp,[filesep,'geometry_rec']],'geometry_rec_type');
[N,rec_n]=size(Sz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('clip','var')
    if isnumeric(clip)
        Sz_max=clip*max(max(Sz));
        [r,c]=find(Sz>Sz_max);
        Sz(r,c)=Sz(r,c)*clip;
    else
        error('Error: Parameter "clip" should be a number. ')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('scale','var')
    if isnumeric(scale)
        scaling= linspace(1,scale,N);
        for j=1:rec_n
            Sz(:,j)=Sz(:,j).*(2.^scaling');
        end
    else
       error('Parameter "scale" is not a number. ') 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if  strcmp(geometry_rec_type,'downhole')
%     
% elseif strcmp(geometry_rec_type,'surface')
%     imagesc(1:rec_n,(0:(N))*dN_SS*dt,Sz);
% end
imagesc(1:rec_n,(0:(N-1))*dN_SS*dt,Sz);

colormap(flipud(gray)); 
set(gca,'xaxislocation','top');
set(gca,'FontSize',14)
xlabel('Receiver number','FontSize',14)
ylabel('T (s)','FontSize',14)
title('Synthetic seismogram','FontSize',16)

%grid on;

end