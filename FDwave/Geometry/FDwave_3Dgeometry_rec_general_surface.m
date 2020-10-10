function FDwave_3Dgeometry_rec_general_surface(varargin) %
% This function places recievers along the surface at a more general format (e.g., scattering distributed).
% Check is applied so that no reciever is placed in absorbing boundary.
% Complete Syntax:
%        3Dgeometry_rec_general('WFP',path,'RNX_VEC',value, 'RNY_VEC',value,...
%        'RNZ_VEC',value, 'DX',value,'DY',value, 'DZ',value, 'BCTYPE',value, 'NAB',value,'PlotON',option)
% Description of parameters:
%         WFP   :  Path to working directory
%         RNX_VEC :  Location vector of Receivers along x in terms of nodes
%         RNY_VEC :  Location vector of Receivers along y in terms of nodes
%         RNZ_VEC :  Location vector of Receivers along z in terms of nodes
%         DX    :  Grid size in x direction
%         DY    :  Grid size in y direction
%         DZ    :  Grid size in z direction
%         BCTYPE:  Type of boundary used
%         NAB   :  Number of grid nodes used for boundaries.
%         PlotON:  'y'/'n'
% Note: 
% 1) Give all the position after adding/subtracting the ABC layer thickness according to side where they are added.
% 2) All parameters are mandatory.
%       In case any parameter is not provided the program will try to use the values stored in input folder. 
%       If it cannot find any stored value then it shows error.
% Example:
%       3Dgeometry_rec_general



for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                    wfp=varargin{i+1};
        case 'rnx_vec';                rnx_vec=varargin{i+1};
        case 'rny_vec';                rny_vec=varargin{i+1};
        case 'rnz_vec';                rnz_vec=varargin{i+1};
        case 'dx';                     dx=varargin{i+1};
        case 'dy';                     dy=varargin{i+1};
        case 'dz';                     dz=varargin{i+1};
        case 'nx';                     nx=varargin{i+1};
        case 'ny';                     ny=varargin{i+1};
        case 'nz';                     nz=varargin{i+1};
        case 'bctype';                 BCtype=varargin{i+1};
        case 'nab';                    nAB=varargin{i+1};
        case 'ploton';                 plotON=varargin{i+1};
        case 'verbose';                verbose=varargin{i+1};
        otherwise;  error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';    str1=str0;      str2=str0;      str3=str0;      str4=str0;      
str5=str0;      str6=str0;      str7=str0;      str8=str0;     
if ~exist('wfp','var');                                     wfp=pwd;                       end;
if ~exist('plotON','var');                                  plotON='y';                       end;
if ~exist('dx','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dx');            str1=[str0,'(Default, stored) '];          end;
if ~exist('dy','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dy');            str2=[str0,'(Default, stored) '];          end;
if ~exist('dz','var')||strcmp(num2str(dv),'-9999');         load([wfp,[filesep,'model']],'dz');            str3=[str0,'(Default, stored) '];          end;
if ~exist('nx','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'nx');            str4=[str0,'(Default, stored) '];          end;
if ~exist('ny','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'ny');            str5=[str0,'(Default, stored) '];          end;
if ~exist('nz','var')||strcmp(num2str(nv),'-9999');         load([wfp,[filesep,'model']],'nz');            str6=[str0,'(Default, stored) '];          end;
if ~exist('BCtype','var')||strcmp(BCtype,'-9999');          load([wfp,[filesep,'BC']],'BCtype');           str7=[str0,'(Default, stored) '];          end;
if ~exist('nAB','var')||strcmp(num2str(nAB),'-9999');       load([wfp,[filesep,'BC']],'nAB');              str8=[str0,'(Default, stored) '];            end;


if ~exist('rnx_vec','var'); error('rnx_vec is missing!');        end;
if ~exist('rny_vec','var'); error('rny_vec is missing!');        end;
if ~exist('rnz_vec','var'); error('rnz_vec is missing!');        end;
if ~exist('verbose','var');        verbose='n';              end;

fprintf('\n')

if ~exist('verbose','var');    verbose='y';     end

geometry_rec=sub2ind([nz,nx,ny],rnz_vec,rnx_vec,rny_vec);
rec_n=length(geometry_rec);

geometry_rec_type = 'surface';


str=strcat(wfp,[filesep,'geometry_rec']);
save(str,'geometry_rec','rnz_vec','rnx_vec','rny_vec', 'rec_n','geometry_rec_type');

if strcmp(verbose,'y')
    disp([str0,'Receiver geometry saved in "',str,'"'])
end

if strcmp(plotON,'y2d')
    geometry_3Dplot_rexz('BCPlot','y');
else if strcmp(plotON,'y3d')
        geometry_3Dplot_rec('BCPlot','y');
    end
end
end

