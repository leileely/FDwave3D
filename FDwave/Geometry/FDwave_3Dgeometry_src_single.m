function FDwave_3Dgeometry_src_single(varargin)  %sv,shn,dh,dv,BCtype,nAB )
% This function places a single source at the given location
% Check is applied so that no reciever is placed in absorbing boundary.
% Complete Syntax:
%     geometry_src_Single('WFP',path,'SX',value,'SY',value, 'SZ',value,'DX',value,'DY',value,'DZ',value,
%          'NX',value,'NY',value,'NZ',value,'BCTYPE',value, 'NAB',value)
% Description of parameters:
%         WFP   :  Path to working directory
%         SX    :  Node Location of source (x-position)
%         SY    :  Node Location of source (y-position)
%         SZ    :  Node Location of source (z-position)
%         DX    :  Grid size in x direction
%         DY    :  Grid size in y direction
%         DZ    :  Grid size in z direction
%         NX    :  No of nodes in model along X axis
%         NY    :  No of nodes in model along Y axis
%         NZ    :  No of nodes in model along Z axis
%         BCTYPE:  Type of boundary used
%         NAB   :  Number of grid nodes used for boundaries.
%         PlotON:  'y'      to plot the source positions
% Note: 
% 1) Do not place source very close to surface.
% 2) Give all the position after adding/subtracting the ABC layer thickness 
%    according to side where they are added. 
%    e.g. if you want to place source at some depth then grid location for
%         topABC is then, Sz = nAB + dept/dz 
% 3) All parameters are mandatory.
%       In case any parameter is not provided the program will try to use 
%       the values stored in input folder. 
%       If it cannot find any stored value then it shows error.
% Example:
%       geometry_src_Single('WFP',pwd,'Sx',2000,'Sy',1000,'Sz',300,'dx',5,'dy',5,'dz',5,'BCtype','topABC','nAB',50,'PlotON','Y');
%       geometry_src_Single('Sx',2000,'Sy',1000,'Sz',300)
%       geometry_src_Single()


for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                 wfp=varargin{i+1};
        case 'sx';                  snx=varargin{i+1};
        case 'sy';                  sny=varargin{i+1};
        case 'sz';                  snz=varargin{i+1};
        case 'dx';                  dx=varargin{i+1};
        case 'dy';                  dy=varargin{i+1};
        case 'dz';                  dz=varargin{i+1};
        case 'nx';                  nx=varargin{i+1};
        case 'ny';                  ny=varargin{i+1};
        case 'nz';                  nz=varargin{i+1};
        case 'bctype';              BCtype=varargin{i+1};
        case 'nab';                 nAB=varargin{i+1};
        case 'ploton';              plotON=varargin{i+1};
        case 'verbose';             verbose=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';   str1=str0;      str2=str0;      str3=str0;      str4=str0;      str5=str0;      str6=str0;      str7=str0;      str8=str0;         str9=str0;      str10=str0;      str11=str0; 

if ~exist('wfp','var');        wfp=pwd;   end;
if ~exist('dx','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dx');            str1=[str0,'(Default, stored) '];          end;
if ~exist('dy','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dy');            str2=[str0,'(Default, stored) '];          end;
if ~exist('dz','var')||strcmp(num2str(dv),'-9999');         load([wfp,[filesep,'model']],'dz');            str3=[str0,'(Default, stored) '];          end;
if ~exist('nx','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'nx');            str4=[str0,'(Default, stored) '];          end;
if ~exist('ny','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'ny');            str5=[str0,'(Default, stored) '];          end;
if ~exist('nz','var')||strcmp(num2str(nv),'-9999');         load([wfp,[filesep,'model']],'nz');            str6=[str0,'(Default, stored) '];          end;
if ~exist('BCtype','var')||strcmp(BCtype,'-9999');          load([wfp,[filesep,'BC']],'BCtype');           str7=[str0,'(Default, stored) '];          end;
if ~exist('nAB','var')||strcmp(num2str(nAB),'-9999');       load([wfp,[filesep,'BC']],'nAB');              str8=[str0,'(Default, stored) '];          end;
if ~exist('verbose','var');        verbose='n';              end;

if ~exist('snx','var')||strcmp(num2str(snx),'-9999');         snx = round(.5*nx);             str1=[str0,'(Default) '];        end;
if ~exist('sny','var')||strcmp(num2str(sny),'-9999');         sny = round(.5*ny);             str2=[str0,'(Default) '];        end;
if ~exist('snz','var')||strcmp(num2str(snz),'-9999');         snz= 6;                       str3=[str0,'(Default) '];        end;
if ~exist('plotON','var');   plotON='n'; end

% if strcmp(BCtype,'topABC');    
%  snv= snv+nAB;
% end
if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
disp('    FUNC: Source geometry (Single at surface)')
disp([str0,'Source location along x(m)       :  ', num2str(snx*dx),str1])
disp([str0,'Source location along y(m)       :  ', num2str(sny*dy),str2])
disp([str0,'Source location along z(m)       :  ', num2str(snz*dz),str3])
disp([str0,'Source node location along x     :  ', num2str(snx),str1])
disp([str0,'Source node location along y     :  ', num2str(sny),str2])
disp([str0,'Source node location along z     :  ', num2str(snz),str3])
disp([str0,'Note: the program take care of absorbing boundary nodes (if present).'])
if strcmp(verbose,'y')
    disp([str0,'Used parameters are'])
    disp([str0,str0,'Model dx     :  ' , num2str(dx),str4])
    disp([str0,str0,'Model dy     :  ' , num2str(dy),str5])
    disp([str0,str0,'Model dz     :  ' , num2str(dz),str6])
    disp([str0,str0,'Model nx     :  ' , num2str(nx),str7])
    disp([str0,str0,'Model ny     :  ' , num2str(ny),str8])
    disp([str0,str0,'Model nz     :  ' , num2str(nz),str9])
    disp([str0,str0,'BCtype       :  ' , BCtype,str10])
    disp([str0,str0,'nAB          :  ' , num2str(nAB),str11])
end

end 

snx_vec= snx;   % source location in x-direction
sny_vec= sny;   % source location in y-direction
snz_vec= snz;   % source location in z-direction

geometry_src=sub2ind([nz,nx,ny],snz_vec,snx_vec,sny_vec);
src_n=length(geometry_src);
geometry_src_type = 'single';

str=strcat(wfp,[filesep,'geometry_src']);
save(str,'geometry_src','snz_vec','snx_vec','sny_vec','src_n','geometry_src_type')

if strcmp(verbose,'y')
    disp(['        Source geometry saved in "',str,'"'])
end

if strcmp(plotON,'y2d')
     geometry_3Dplot_srcxz('BCPlot','y');
else if strcmp(plotON,'y3d')
        geometry_3Dplot_src('BCPlot','y');
    end
end
end

