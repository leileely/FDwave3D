function FDwave_3Dgeometry_rec_st_line_surf(varargin) %Dn, HFn,HLn,HDn)
% This function places recievers along the surface at a given depth.
% Check is applied so that no reciever is placed in absorbing boundary.
% Complete Syntax:
%        geometry_rec_st_line_surf('WFP',path,'DEPTH',value, 'FIRST',value, 'LAST',value,...
%        'DIFF',value, 'DX',value,'DY',value, 'DZ',value, 'BCTYPE',value, 'NAB',value,'PlotON',option)
% Description of parameters:
%         WFP   :  Path to working directory
%         DEPTH :  Depth in terms of nodes
%         FIRST :  Location of first Receiver along x in terms of nodes
%         LAST  :  Location of last Receiver  along x in terms of nodes
%         DIFF  :  Difference in consecutive Receiver  along x in terms of nodes
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
%       geometry_rec_st_line_surf('WFP',pwd,'DEPTH',300,'FIRST',260,'LAST',4400,'DIFF',20,'DX',5,'DY',5,'DZ',5,'BCTYPE','topABC','NAB',50);
%       geometry_rec_st_line_surf('DEPTH',300,'FIRST',260,'LAST',4400,'DIFF',20)
%       geometry_rec_st_line_surf()



for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                    wfp=varargin{i+1};
        case 'depth';                  Dn=varargin{i+1};
        case 'first';                  HFn=varargin{i+1};
        case 'last';                   HLn=varargin{i+1};
        case 'diff';                   HDn=varargin{i+1};
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
str5=str0;      str6=str0;      str7=str0;      str8=str0;      str9=str0;      str10=str0;            str11=str0;      str12=str0;  

if ~exist('wfp','var');                                     wfp=pwd;                       end;
if ~exist('plotON','var');                                  plotON='y';                       end;
if ~exist('dx','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dx');            str1=[str0,'(Default, stored) '];          end;
if ~exist('dy','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dy');            str2=[str0,'(Default, stored) '];          end;
if ~exist('dz','var')||strcmp(num2str(dv),'-9999');         load([wfp,[filesep,'model']],'dz');            str3=[str0,'(Default, stored) '];          end;
if ~exist('nx','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'nx');            str4=[str0,'(Default, stored) '];          end;
if ~exist('ny','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'ny');            str5=[str0,'(Default, stored) '];          end;
if ~exist('nz','var')||strcmp(num2str(nv),'-9999');         load([wfp,[filesep,'model']],'nz');            str6=[str0,'(Default, stored) '];          end;
if ~exist('BCtype','var')||strcmp(BCtype,'-9999');          load([wfp,[filesep,'BC']],'BCtype');           str7=[str0,'(Default, stored) '];          end;
if ~exist('nAB','var')||strcmp(num2str(nAB),'-9999');       load([wfp,[filesep,'BC']],'nAB');              str8=[str0,'(Default, stored) '];          end;

if ~exist('Dn','var')||strcmp(num2str(Dn),'-9999');
    Dn=10;
    if  strcmp(BCtype,'topABC');     Dn=Dn+nAB;    end
    str9=[str0,'(Default) '];
end;

if ~exist('HFn','var')||strcmp(num2str(HFn),'-9999');        HFn=1+nAB;           str10=[str0,'(Default) '];        end;
if ~exist('HLn','var')||strcmp(num2str(HLn),'-9999');        HLn=nh-nAB;          str11=[str0,'(Default) '];          end;
if ~exist('HDn','var')||strcmp(num2str(HDn),'-9999');        HDn=1;           str12=[str0,'(Default) '];          end;
if ~exist('verbose','var');        verbose='n';              end;

fprintf('\n')
if strcmp('BCtype','topFS')
    disp([str0,'Depth of all receiver (m)   :  ' , num2str(Dn*dv),str1])
elseif strcmp('BCtype','topABC')
    disp([str0,'Depth of all receiver (m)   :  ' , num2str((Dn-nAB)*dv),str1])
end

if ~exist('verbose','var');    verbose='y';     end

if strcmp(verbose,'y')
    disp('    FUNC: Reciever geometry (Along surface)')
    disp([str0,'First receiver location (m) :  ', num2str(HFn*dx),str2])
    disp([str0,'Last receiver location (m)  :  ' ,num2str(HLn*dx),str3])
    disp([str0,'Distance between two consecutive Receiver :  ' , num2str(HDn*dx),str4])
    
    fprintf('\n')
    disp([str0,'Position of receiver on grid'])
    disp([str0,'Depth of all receiver (in terms of node)   :  ' , num2str(Dn),str1])
    disp([str0,'First receiver location (in terms of node) :  ', num2str(HFn),str2])
    disp([str0,'Last receiver location (in terms of node)  :  ' ,num2str(HLn),str3])
    disp([str0,'Distance between two consecutive Receiver (in terms of node):  ' , num2str(HDn),str4])
    fprintf('\n')
    disp('        Note: Node position in grid(N) & distances(X) related by (N= nAB + X/h), due to ABC layer.');
    fprintf('\n')
    if strcmp(verbose,'y')
        disp([str0,'Used parameters are: '])
        disp([str0,str0,'dx     :  ' , num2str(dx),str5])
        disp([str0,str0,'dy     :  ' , num2str(dy),str6])
        disp([str0,str0,'dz     :  ' , num2str(dz),str7])
        disp([str0,str0,'nx     :  ' , num2str(nx),str8])
        disp([str0,str0,'ny     :  ' , num2str(ny),str9])
        disp([str0,str0,'nz     :  ' , num2str(nz),str10])
        disp([str0,str0,'BCtype :  ' , BCtype,str11])
        disp([str0,str0,'nAB    :  ' , num2str(nAB),str12])
    end
end
rnx_vec = HFn:HDn:HLn;
rny_vec = round(ny/2)*ones(size(rnx_vec));
rnz_vec= Dn*ones(size(rnx_vec));

geometry_rec=sub2ind([nz,nx,ny],rnz_vec,rnx_vec,rny_vec);
rec_n=length(geometry_rec);

geometry_rec_type = 'surface';


str=strcat(wfp,[filesep,'geometry_rec']);
save(str,'geometry_rec','rnz_vec','rnx_vec','rny_vec', 'rec_n','geometry_rec_type');

if strcmp(verbose,'y')
    disp([str0,'Receiver geometry saved in "',str,'"'])
end

if strcmp(plotON,'y2d')
    geometry_3Dplot_recxz('BCPlot','y');
else if strcmp(plotON,'y3d')
        geometry_3Dplot_rec('BCPlot','y');
    end
end
end

