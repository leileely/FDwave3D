function FDwave_analyse_3Delastic(varargin)%dt,f0,dx,dy,dz,vpm,vsm)

%       This function is used to check if parameters are satisfying CFL and dispersion condition in 3D scenario.
%       In case the conditions does not satisfy, this program will show a warning and continue exectuion.
% Complete Syntax:
%       analyse_elastic('WFP',path,'DT',value, 'F0',value, 'DX',value, 'DY',value, 'DZ',value, 'VP',value,'VS',value)
% Description of parameters:
%       WFP : Path to working directory
%       DT  : Time step
%       F0  : Central frequency of source
%       DX  : Grid spacing along x
%       DY  : Grid spacing along y
%       DZ  : Grid spacing along z
%       VP  : Velocity of P wave
%       VS  : Velocity of S wave
%
%       Note: All parameters are mandatory.
%       In case any parameter is not provided the program will try to use the values present in input folder.
%       If it cannot find any stored value then it shows error.
% Example:
%       analyse_3Delastic('WFP',pwd,'DT',.0002,'F0',20,'DX',5,'DY',5,'DZ',5,'VP',2300,'VS',2100)
%       analyse_3Delastic()



for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                 wfp=varargin{i+1};
        case 'dt';                  dt=varargin{i+1};
        case 'f0';                  f0=varargin{i+1};
        case 'dx';                  dx=varargin{i+1};
        case 'dy';                  dy=varargin{i+1};
        case 'dz';                  dz=varargin{i+1};
        case 'vp';                  vpm=varargin{i+1};
        case 'vs';                  vsm=varargin{i+1};            
        case '2n';                  M=varargin{i+1};
        case 'verbose';             verbose=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';
str1=str0;          str2=str0;      str3=str0;      str4=str0;      str5=str0;      str6=str0;      str7=str0;    

if ~exist('wfp','var');                                    wfp=pwd;                end;
if ~exist('dt','var')||strcmp(num2str(dt),'-9999');        load([wfp,[filesep,'source']],'dt');         str1=[str0,'(Default, stored)'];         end;
if ~exist('f0','var')||strcmp(num2str(f0),'-9999');        load([wfp,[filesep,'source']],'f0');         str2=[str0,'(Default, stored)'];         end;
if ~exist('dx','var')||strcmp(num2str(dh),'-9999');        load([wfp,[filesep,'model']],'dx');          str3=[str0,'(Default, stored)'];         end;
if ~exist('dy','var')||strcmp(num2str(dh),'-9999');        load([wfp,[filesep,'model']],'dy');          str4=[str0,'(Default, stored)'];         end;
if ~exist('dz','var')||strcmp(num2str(dv),'-9999');        load([wfp,[filesep,'model']],'dz');          str5=[str0,'(Default, stored)'];         end;
if ~exist('vpm','var')||strcmp(num2str(vpm),'-9999');      load([wfp,[filesep,'model']],'vpm');         str6=[str0,'(Default, stored)'];         end;
if ~exist('vsm','var')||strcmp(num2str(vsm),'-9999');      load([wfp,[filesep,'model']],'vsm');         str7=[str0,'(Default, stored)'];         end;
if ~exist('M','var')    M=4 ; end
if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
    disp('    FUNC: Analysis for Elastic wave')
    disp([str0,'Time step size,         dt =  ',num2str(dt),str1] )
    disp([str0,'Source Frequency,       f0 =  ',num2str(f0),str2] )
    disp([str0,'Grid spacing Taken,     dx =  ',num2str(dx),str3] )
    disp([str0,'Grid spacing Taken,     dy =  ',num2str(dy),str4] )
    disp([str0,'Grid spacing Taken,    dz =  ',num2str(dz),str5] )
    disp([str0,'Vp model supplied',str6] )
    disp([str0,'Vs model supplied',str7] )
    disp([str0,str0,'Analysis Results:'])
end
%%%%%%%%%%%%%%%%%%%%%%stability condition%%%%%%%
% see eq.(13)of the related manuscript
NN=M/2;
fdc=DiffCoef(NN,'s');% coefficients of FD operator
p=min(min(dx,dy),dz); 
%check=0.606*p/max(max(vpm));%2D(~ 1/sqrt(2))
check=(1/(sum(abs(fdc))*sqrt(3)))*p/max(vpm(:));%3D(~ 1/sqrt(3))fourth-order=0.495

if strcmp(verbose,'y')
    if(dt<=check)
        fprintf('        Given time step: ok;        "dt" Present value( %f) < Required( %f). \n',dt,check);
    else
        error('        Revise time step...!!!      "dt" Present value( %f) > Required( %f). \n',dt,check);
        %     temp= input('');
    end
end
%%%%%%%%%%%%%%%%%%%dispersion condition%%%%%%%
% see eq.(14)of the related manuscript
p=max(max(dx,dy),dz);
vsm_min=min(vsm(:));
nwave=8-(M-4)/2;% number of grid points per minimum wavelength (refer to Bohlen et al.,2015)
if vsm_min~=0
    check = vsm_min/nwave/f0;
else
    if strcmp(verbose,'y');        disp('        Since, min(Vs)=0 so using the min(Vp)') ; end
    vpm_min=min(vpm(:));
    check = vpm_min/nwave/f0;
end
%%%%%%%%%%%%%%%%%%%
if strcmp(verbose,'y')
    if (p<=check )
        fprintf('        Given grid spacing:  ok;      "dx" Present value( %f) < Required( <%f). \n',p,check);
    else
        warning('        Revise grid spacing...!!!!      "dx" Present value( %f) > Required( <%f). \n',p,check);
        %     temp= input('');
    end
    
    fprintf('        CFL and Grid size requirement analysis done.  \n \n');
end
end