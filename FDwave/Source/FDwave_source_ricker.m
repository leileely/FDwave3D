function [src,t] = FDwave_source_ricker(varargin)   %T,dt,f0,t0,src_scale)

% This function generates the ricker wavelet signature.
% Complete Syntax:
%       source_ricker('WFP',path,'T',value,'DT',value,'F0',value,'T0',value,'SRC_SCALE',value,'PlotON','y');
% Description of parameters:
%       WFP         :  Path to working directory
%       Ts          : Total time duration for the source wavelet.
%       DT          : Time step
%       F0          : Central frequency of source
%       T0          : Zero offset time aka Lag time (optional)
%       SRC_SCALE   : Scaling of amplitude (optional)
%       PlotON      : 'y'/'n' 
% Example:
%       source_ricker('WFP',pwd,'T',2,'Ts',0.2,'DT',.0004,'F0',10)



for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';        wfp=varargin{i+1};
        case 't';       T=varargin{i+1};
        case 'dt';      dt=varargin{i+1};
        case 'f0';      f0=varargin{i+1};
        case 't0';      t0=varargin{i+1};
        case 'src_scale'; src_scale=varargin{i+1};
        case 'ploton';   plotON=varargin{i+1};
        case 'verbose';         verbose=varargin{i+1};
        otherwise;      error('%s is not a valid argument name',varargin{i});
            
    end
end

if ~exist('wfp','var');     wfp=pwd;                                                    end;
if ~exist('T','var');       error('Provide total time duration, T, for simulation');    end 
if ~exist('dt','var');      error('Provide time step, dt, for simulation');             end
if  ~exist('f0','var');     error('Provide source frequency, f0, for simulation ');     end

str0='        ';   
str1=str0;          str2=str0;      

if  ~exist('t0','var')||strcmp(num2str(t0),'-9999');    t0 = 1/(sqrt(2)*pi*f0);    str1=[str0,'(Default)'];    end
if ~exist('src_scale','var')||strcmp(num2str(src_scale),'-9999');   src_scale=1;          str2=[str0,'(Default)'];   end
if  ~exist('plotON','var');     plotON='n';     end
if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
    disp('    FUNC: Source parameters')
    
    disp([str0,'Source type                :    Ricker Wavelet'])
    disp([str0,'Source Frequency,       f0 =  ',num2str(f0)] )
    disp([str0,'Time step size,         dt =  ',num2str(dt)] )
    disp([str0,'Total time duration,     T =  ',num2str(T)] )
    disp([str0,'Time shift,             t0 =  ',num2str(t0),str1] )
    disp([str0,'Source amplitude scaling   =  ',num2str(src_scale),str2] )
end
N = round(T/dt)+1;    % No of total time steps


t = dt*(0:N-1);
tau = t-t0;     % Create the wavelet and shift in time

src = src_scale*(1 -2* tau.*tau * f0^2 * pi^2).*exp(-tau.^2 * pi^2 * f0^2);

% used for plotting the wavelet. td is time difference between side lobes of ricker.
tdn =round(t0/dt) + 2*round(sqrt(6)/(pi*f0)/dt);  
Ns=tdn;
src_name='Ricker';
str = strcat(wfp,[filesep,'source']);
save(str,'src_name','T','dt','f0','t0','src_scale','src','t','N','Ns','tdn');

if strcmp(verbose,'y')
    disp(['        Source saved in "',str])
end
if strcmp(plotON,'y')
    source_plot('wfp',wfp)
end

