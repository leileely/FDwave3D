function FDwave_model_3Dload(varargin)  % wave_type,   dx,dy,dz,HV_ratio,tl,  vpv,vsv,rhov,qpv,qsv)

% This function can load a 3D model 
% Complete Syntax:
% FDwave_model_3Dload('WFP',wf_path,'WAVE_TYPE','Elastic','DX',10,'DY',10,'DZ',10,...
%     'VFilename',vfile,'PlotON','y3d','verbose','y')%
%       Description of parameters:
%       WFP          :  Path to working folder
%       WAVE_TYPE    :  'acoustic1', 'acoustic2', 'elastic', 'viscoelastic'
%       DX, DY, DZ       :  Grid size in horizontal and vertical direction in meters
%       VP, VS       :  Velocity of P and S wave in form of vector
%       RHO          :  Density of medium in form of vector
%       PlotON       :  'y2d'/'y3d'/'n' for plotting
% Note:
%       acoustic1 - requires VP
%       acoustic2 - requires VP, RHO
%       Elastic   - requires VP, VS, RHO
%       viscoelastic - require VP, VS, RHO, QP, QS
%       The vector is in the form of e.g. [600,700,800,1000]
% Example:

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wave_type';       wave_type=varargin{i+1};
        case 'wfp';             wfp=varargin{i+1};
        case 'dx';              dx=varargin{i+1};
        case 'dy';              dy=varargin{i+1};
        case 'dz';              dz=varargin{i+1};            
        case 'nab';             nab=varargin{i+1};
        case 'vp';              vp=varargin{i+1};
        case 'vs';              vs=varargin{i+1};
        case 'rho';             rho=varargin{i+1};
        case 'c';               C=varargin{i+1};
        case 'verbose';         verbose=varargin{i+1};
        case 'ploton';          plotON=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';
str1=str0;      str2=str0;      str3=str0;      str4=str0;

if ~exist('wfp','var');             wfp=pwd;                  end

if ~exist('wave_type','var');       wave_type = 'elastic';      str1=[str0,'(Default)'];
else        wave_type=lower(wave_type);    end

if ~exist('dx','var');              dx= 4;                      str2=[str0,'(Default)'];  end
if ~exist('dy','var');              dy= 4;                      str3=[str0,'(Default)'];  end
if ~exist('dz','var');              dz= 4;                      str4=[str0,'(Default)'];  end

if ~exist('verbose','var');         verbose='y';                end
if ~exist('plotON','var');          plotON='n';                 end

if strcmp(verbose,'y')
    disp('    FUNC: Model Building')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([str0,'Model Selected: loading model'])
    disp([str0,str0,'Provided parameters']);
    disp([str0,'Type of wave         : ',wave_type,str1] )
    disp([str0,'Grid spacing along x : ',num2str(dx),str2])
    disp([str0,'Grid spacing along y : ',num2str(dy),str3])
    disp([str0,'Grid spacing along z : ',num2str(dz),str4])
    
end

[nz,nx,ny]=size(vp);


vpm=vp;
vsm=vs;
rhom=rho;
C11=C{1};
C12=C{2};
C13=C{3};
C22=C{4};
C23=C{5};
C33=C{6};
C44=C{7};
C55=C{8};
C66=C{9};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str=strcat(wfp,[filesep,'model']);
m_name='load_model';

if strcmp(wave_type,'acoustic1')    
    save(str,'dx','dy','dz','nx','ny','nz','vpm','wave_type','m_name');
    
elseif strcmp(wave_type,'acoustic2')
    save(str,'dx','dy','dz','nx','ny','nz','vpm','rhom','wave_type','m_name');
    
elseif strcmp(wave_type,'elastic')
    save(str,'dx','dy','dz','nx','ny','nz','vpm','vsm','rhom','C11','C12','C13','C22','C23','C33','C44','C55','C66','wave_type','m_name');
    
elseif strcmp(wave_type,'viscoelastic')
    save(str,'dx','dy','dz','nx','ny','nz','vpm','vsm','rhom','C11','C12','C13','C22','C23','C33','C44','C55','C66','qpm','qsm','wave_type','m_name');
else
    warning('      Wrong name for model entered, No model saved')
end

if strcmp(verbose,'y')
    if strcmp(wave_type,'acoustic1')||strcmp(wave_type,'acoustic2')||strcmp(wave_type,'elastic')||strcmp(wave_type,'viscoelastic')
        disp(['        dx =  ',num2str(dx)] )
        disp(['        dy =  ',num2str(dy)] )
        disp(['        dz =  ',num2str(dz)] )
        disp(['        nx =  ',num2str(nx)] )
        disp(['        ny =  ',num2str(ny)] )
        disp(['        nz =  ',num2str(nz)] )
        disp(['        Model saved in',str])
        
        if strcmp(plotON,'y2d')
            model_2Dplot('wfp',wfp);
        else if strcmp(plotON,'y3d')
            model_3Dplot('wfp',wfp);
            end
        end
    end
end
