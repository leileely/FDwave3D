function  FDwave_calculation_3Dwiggle_plot( varargin)%ipstr,    xF,xL,xD,    tF,tL,tD,  figNo, shotno )
% This function plots seismograms in wiggle style.
%   This function is under testing
%   Detailed explanation goes here

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                       wfp=varargin{i+1};
        case 'ssfilename';           FileName=varargin{i+1};
        case 'wave_type';           wave_type = lower(varargin{i+1});
        case 'geometry_type';    geometry_type= lower(varargin{i+1});
        case 'dx';                        dx= lower(varargin{i+1});
        case 'dt';                         dt= lower(varargin{i+1});
        case 'figno';                    figNO = varargin{i+1};
        case 'rnf';                          rnF = varargin{i+1};
        case 'rnd';                         rnD = varargin{i+1};
        case 'rnl';                          rnL = varargin{i+1};
        case 'tf';                          tF = varargin{i+1};
        case 'tl';                          tL = varargin{i+1};
        case 'td';                          tD = varargin{i+1};
        case 'clip';                      clip = varargin{i+1};
        case 'scale';                   scale = varargin{i+1};
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

load(FileName,'Sz','dN_SS'); 
if ~exist('rnF','var');      rnF=1;       end
if ~exist('rnL','var');      rnL=size(Sz,2);       end
if ~exist('rnD','var');      rnD=1;       end

if ~exist('tF','var');      tF=1;       end
if ~exist('tL','var');      tL= size(Sz,1) ;       end
if ~exist('tD','var');      tD=1;       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,rec_n]=size(Sz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('clip','var')
    if isnumeric(clip)
        SS_max=clip*max(max(Sz));
        [r,c]=find(Sz>SS_max);
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
% xvecn =  rnF:rnD:rnL;
% 
% tvec = (tF:tD:tL) + dt;
% tvecn = round(tvec/dt);

wiggle(Sz);

set(gca,'FontSize',14)
xlabel('Receiver number','FontSize',14)
ylabel('T (s)','FontSize',14)
set(gca,'ytick',1:round(N/5):N,'yticklabel',(0:round(N/5):N-1)*dt);
title('Synthetic seismogram','FontSize',16)
end

