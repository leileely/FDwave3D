function FDwave_3Dgeometry_rec_downhole(varargin)%H,DF,DL,DD)
% This function places recievers vertically along the depth at some given spacing.
% The configration is similar to the VSP configration of receivers.
%   Complete Syntax:
%       FDwave_3Dgeometry_rec_downhole('XLOC',value,'YLOC',value,'FIRST',value,'LAST',value,'DIFF',value...
%                'DX',value,'DY',value,'DZ',value,'BCTYPE',option,'NAB',value,'PlotON',option)
% Description of parameters:
%                XLOC   :  Horizontal location in terms of nodes (X)
%                YLOC   :  Horizontal location in terms of nodes (Y)
%                FIRST  :  Location of first Receiver along x in terms of nodes
%                LAST   :  Location of last Receiver  along x in terms of nodes
%                DIFF   :  Difference in consecutive Receiver along x in terms of nodes
%                DX     :  Grid size in x direction
%                DY    :  Grid size in y direction
%                DZ     :  Grid size in y direction
%                BCTYPE :  Type of boundary used
%                NAB    :  Number of grid nodes used for boundaries.
%                PlotON :  'y'/'n'
% Note: 
% 1) Give all the position after adding/subtracting the ABC layer thickness according to side where they are added.
%    Node position in grid(N) & distances(X) related by (N= nAB + X/h), due to ABC layer.'
% 2) In case of free surface two topmost node layers are ignored (imaging technique for Free surface)
%    Hence a node vertical position is given by  N= (X/h)-2;
%
% 3) All parameters are mandatory.
%       In case any parameter is not provided the program will try to use the values stored in input folder. 
%       If it cannot find any stored value then it shows error.


disp('    FUNC: Receiver geometry Along Depth (or VSP type)')

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wfp';                    wfp=varargin{i+1};
        case 'xloc';                   Xn=varargin{i+1};            
        case 'yloc';                   Yn=varargin{i+1};
        case 'first';                  DFn=varargin{i+1};
        case 'last';                   DLn=varargin{i+1};
        case 'diff';                   DDn=varargin{i+1};
        case 'dx';                     dx=varargin{i+1};
        case 'dy';                     dy=varargin{i+1};
        case 'dz';                     dz=varargin{i+1};
        case 'nx';                     nh=varargin{i+1};
        case 'ny';                     ny=varargin{i+1};
        case 'nz';                     nv=varargin{i+1};
        case 'bctype';                 BCtype=varargin{i+1};
        case 'nab';                    nAB=varargin{i+1};
        case 'ploton';                 plotON=varargin{i+1};
        case 'verbose';                verbose=varargin{i+1};
        otherwise;  error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';    str1=str0;      str2=str0;      str3=str0;      str4=str0;      
str5=str0;          str6=str0;      str7=str0;      str8=str0;      str9=str0;      str10=str0;       str11=str0;      str12=str0;         str13=str0; 

if ~exist('wfp','var');                                     wfp=pwd;                       end;
if ~exist('dx','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dx');            str1=[str0,'(Default, stored) '];          end;
if ~exist('dy','var')||strcmp(num2str(dh),'-9999');         load([wfp,[filesep,'model']],'dy');            str2=[str0,'(Default, stored) '];          end;
if ~exist('dz','var')||strcmp(num2str(dv),'-9999');         load([wfp,[filesep,'model']],'dz');            str3=[str0,'(Default, stored) '];          end;
if ~exist('nx','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'nx');            str4=[str0,'(Default, stored) '];          end;
if ~exist('ny','var')||strcmp(num2str(nh),'-9999');         load([wfp,[filesep,'model']],'ny');            str5=[str0,'(Default, stored) '];          end;
if ~exist('nz','var')||strcmp(num2str(nv),'-9999');         load([wfp,[filesep,'model']],'nz');            str6=[str0,'(Default, stored) '];          end;
if ~exist('BCtype','var')||strcmp(BCtype,'-9999');          load([wfp,[filesep,'BC']],'BCtype');           str7=[str0,'(Default, stored) '];          end;
if ~exist('nAB','var')||strcmp(num2str(nAB),'-9999');       load([wfp,[filesep,'BC']],'nAB');              str8=[str0,'(Default, stored) '];          end;
if ~exist('Xn','var')||strcmp(num2str(Xn),'-9999');         Xn=round(nx/2);                str9=[str0,'(Default) '];                  end;
if ~exist('Yn','var')||strcmp(num2str(Yn),'-9999');         Yn=round(ny/2);                str10=[str0,'(Default) '];                  end;
if ~exist('DFn','var')||strcmp(num2str(DFn),'-9999');        
    if strcmp(BCtype,'topFS');              DFn = 3;
    elseif strcmp(BCtype, 'topABC');        DFn = 1+nAB;
    end
    str11=[str0,'(Default) '];        
end;
if ~exist('DLn','var')||strcmp(num2str(DLn),'-9999');        DLn=nz-nAB;          str12=[str0,'(Default) '];          end;
if ~exist('DDn','var')||strcmp(num2str(DDn),'-9999');        DDn=1;               str13=[str0,'(Default) '];          end;
if ~exist('verbose','var');                                  verbose='n';              end;
if ~exist('plotON','var');   plotON='y';  end;

% fprintf('\n')
    disp([str0,'Horizontal position of receiver (m)        :  ', num2str(Xn*dx),str2])    
    disp([str0,'Horizontal position of receiver (m)        :  ', num2str(Yn*dy),str3])
if strcmp(BCtype,'topFS')
    disp([str0,'First receiver depth (m)                   :  ', num2str((DFn-2)*dx),str3])
    disp([str0,'Last receiver depth (m)                    :  ' ,num2str((DLn-2)*dx),str4])
    disp([str0,'Distance between two consecutive receivers :  ' , num2str((DDn)*dx),str5])
elseif strcmp(BCtype,'topABC')
    disp([str0,'First receiver depth (m)                   :  ', num2str((DFn-nAB)*dx),str3])
    disp([str0,'Last receiver depth (m)                    :  ' ,num2str((DLn-nAB)*dx),str4])
    disp([str0,'Distance between two consecutive receivers :  ' , num2str((DDn)*dx),str5])
end

fprintf('\n')

disp([str0,'Position of receivers on grid'])
disp([str0,'Horizontal position of all receivers (in terms of node)     :  ', num2str(Xn),str1])
disp([str0,'Horizontal position of all receivers (in terms of node)     :  ', num2str(Yn),str2])
disp([str0,'First receiver location (in terms of node)                  :  ', num2str(DFn),str3])
disp([str0,'Last receiver location (in terms of node)                   :  ', num2str(DLn),str4])
disp([str0,'Distance between two consecutive receiver (in terms of node):  ' , num2str(DDn),str5])
fprintf('\n')

% fprintf('\n')

if strcmp(verbose,'y')
    disp([str0,'Used parameters are: '])
    disp([str0,str0,'dx     :  ' , num2str(dx),str6])
    disp([str0,str0,'dy     :  ' , num2str(dy),str7])
    disp([str0,str0,'dz     :  ' , num2str(dz),str8])
    disp([str0,str0,'nx     :  ' , num2str(nx),str9])
    disp([str0,str0,'ny     :  ' , num2str(ny),str10])
    disp([str0,str0,'nz     :  ' , num2str(nz),str11])
    disp([str0,str0,'BCtype :  ' , BCtype,str12])
    disp([str0,str0,'nAB    :  ' , num2str(nAB),str13])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if strcmp(BCtype,'topFS')
%     rnv_vec= 1:nv-nAB;         % all recievers location in vertical direction for VSP
% elseif strcmp(BCtype,'topABC')
%     rnv_vec= 1+nAB:nv-nAB;         % all recievers location in vertical direction for VSP
% end
rnz_vec= DFn:DDn:DLn;
rnx_vec= Xn*ones(size(rnz_vec));    % all recievers location in horizontal direction for VSP (constant)
rny_vec= Yn*ones(size(rnz_vec));    % all recievers location in horizontal direction for VSP (constant)



geometry_rec=sub2ind([nz,nx,ny],rnz_vec,rnx_vec,rny_vec);
rec_n=length(geometry_rec);

geometry_rec_type = 'downhole';

str=strcat(wfp,[filesep,'geometry_rec']);
save(str,'geometry_rec','rnz_vec','rnx_vec','rny_vec','rec_n','geometry_rec_type')

if strcmp(plotON,'y2d')
     geometry_3Dplot_recxz('BCPlot','y');
else if strcmp(plotON,'y3d')
        geometry_3Dplot_rec('BCPlot','y');
    end
end
end

