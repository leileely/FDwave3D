function FDwave_model_n_3Dlayers_TI(varargin)  % wave_type,   dx,dy,dz,HV_ratio,tl,  vpv,vsv,rhov,qpv,qsv)

% This function can create a 3D model consisting of 'n' horizontal layers in transversely isotropic case.
% Complete Syntax:
%     model_N_layers('WFP',path,'WAVE_TYPE',options,'DX',value,'DY',value,'DZ',value,'THICKNESS',value,...
%         'HV_RATIO',value,'VP',value,'VS',value,'RHO',value,'QP',value,'QS',value,'PlotON',option )
%
%       Description of parameters:
%       WFP          :  Path to working folder
%       WAVE_TYPE    :  'acoustic1', 'acoustic2', 'elastic', 'viscoelastic'
%       DX, DY, DZ       :  Grid size in horizontal and vertical direction in meters
%       THICKNESS    :  Thickness of each layer in form of vector
%       XZ_RATIO     :  Total size of the model in respect to
%       VP, VS       :  Velocity of P and S wave in form of vector
%       RHO          :  Density of medium in form of vector
%       C             :  Elastic tensor
%       QP, QS       :  Attenutation of P and S value in form of vector
%       PlotON       :  'y'/'n' for plotting
% Note:
%       acoustic1 - requires VP
%       acoustic2 - requires VP, RHO
%       Elastic   - requires VP, VS, RHO
%       viscoelastic - require VP, VS, RHO, QP, QS
%       The vector is in the form of e.g. [600,700,800,1000]
% Example:
%     model_n_3Dlayers('WFP',pwd,'WAVE_TYPE','Elastic','DX',4,'DY',4,'DZ',4,'THICKNESS',[80,50,60,40,90],...
%          'HV_RATIO',1,'VP',2000:250:3000,'VS',1700:250:2700,'RHO',2300:100:2700)

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'wave_type';       wave_type=varargin{i+1};
        case 'wfp';             wfp=varargin{i+1};
        case 'dx';              dx=varargin{i+1};
        case 'dy';              dy=varargin{i+1};
        case 'dz';              dz=varargin{i+1};
        case 'thickness';       tl=varargin{i+1};
        case 'xz_ratio';        XZ_ratio=varargin{i+1};
        case 'nab';             nab=varargin{i+1};
        case 'vp';              vpv=varargin{i+1};
        case 'vs';              vsv=varargin{i+1};
        case 'rho';             rhov=varargin{i+1};
        case 'c';               C=varargin{i+1};
        case 'qp';              qpv=varargin{i+1};
        case 'qs';              qsv=varargin{i+1};
        case 'verbose';         verbose=varargin{i+1};
        case 'ploton';          plotON=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end

str0='        ';
str1=str0;      str2=str0;      str3=str0;      str4=str0;
str5=str0;      str6=str0;      str7=str0;      str8=str0;
str9=str0;      str10=str0;      str11=str0;   str12=str0;

if ~exist('wfp','var');             wfp=pwd;                  end;

if ~exist('wave_type','var');       wave_type = 'elastic';      str1=[str0,'(Default)'];
else        wave_type=lower(wave_type);    end

if ~exist('dx','var');              dx= 4;                      str2=[str0,'(Default)'];  end
if ~exist('dy','var');              dy= 4;                      str3=[str0,'(Default)'];  end
if ~exist('dz','var');              dz= 4;                      str4=[str0,'(Default)'];  end
if ~exist('XZ_ratio','var');        XZ_ratio=1;                 str5=[str0,'(Default)'];  end

if strcmp(wave_type,'acoustic1')&&(~exist('tl','var')||~exist('vpv','var'))
    tl = [100,46,30,43,100];    str6=[str0,'(Default)'];
    vpv = 2000:250:3000;        str7=[str0,'(Default)'];
elseif strcmp(wave_type,'acoustic2')&&(~exist('tl','var')||~exist('vpv','var')||~exist('rhov','var'))
    tl = [100,46,30,43,100];    str6=[str0,'(Default)'];
    vpv = 2000:250:3000;        str7=[str0,'(Default)'];
    rhov = 2300:100:2700;       str9=[str0,'(Default)'];
elseif strcmp(wave_type,'elastic')&&(~exist('tl','var')||~exist('vpv','var')||~exist('vsv','var')||~exist('rhov','var')||~exist('C','var'))
    tl = [100,46,30,43,100];    str6=[str0,'(Default)'];
    vpv = 2000:250:3000;        str7=[str0,'(Default)'];
    vsv = 1700:250:2700;        str8=[str0,'(Default)'];
    rhov = 2300:100:2700;       str9=[str0,'(Default)'];
    C = [26.39 12.71 6.11 26.39 6.11 15.6 4.38 4.38 6.84]*10^9;       str10=[str0,'(Default)'];
elseif strcmp(wave_type,'viscoelastic')&&(~exist('tl','var')||~exist('vpv','var')||~exist('vsv','var')||~exist('rhov','var')||~exist('qpv','var')||~exist('qsv','var')||~exist('C','var'))
    tl = [100,46,30,43,100];    str6=[str0,'(Default)'];
    vpv = 2000:250:3000;        str7=[str0,'(Default)'];
    vsv = 1700:250:2700;        str8=[str0,'(Default)'];
    rhov = 2300:100:2700;       str9=[str0,'(Default)'];    
    C = [26.39 12.71 6.11 26.39 6.11 15.6 4.38 4.38 6.84]*10^9;       str10=[str0,'(Default)'];
    qpv = 80:5:100;             str11=[str0,'(Default)'];
    qsv = 60:5:80 ;             str12=[str0,'(Default)'];
end

if ~exist('verbose','var');         verbose='y';                end
if ~exist('plotON','var');          plotON='n';                 end

if strcmp(verbose,'y')
    disp('    FUNC: Model Building')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([str0,'Model Selected: N-Layers'])
    disp([str0,str0,'Provided parameters']);
    disp([str0,'Type of wave         : ',wave_type,str1] )
    disp([str0,'Grid spacing along x : ',num2str(dx),str2] )
    disp([str0,'Grid spacing along y : ',num2str(dy),str3])
    disp([str0,'Grid spacing along z : ',num2str(dz),str4])
    disp([str0,'X/Z dimension Ratio  : ',num2str(XZ_ratio),str5])
    
    if strcmp(wave_type,'acoustic1')
        disp([str0,'Thickness of layer(s): ',num2str(tl),str6])
        disp([str0,'Vp of layer(s): ',num2str(vpv),str7])
    elseif strcmp(wave_type,'acoustic2')
        disp([str0,'Thickness of layer(s): ',num2str(tl),str6])
        disp([str0,'Vp of layer(s): ',num2str(vpv),str7])
        disp([str0,'Density of layer(s)  : ',num2str(rhov),str9])
    elseif strcmp(wave_type,'elastic')
        disp([str0,'Thickness of layer(s): ',num2str(tl),str6])
        disp([str0,'Vp of layer(s)       : ',num2str(vpv),str7])
        disp([str0,'Vs of layer(s)       : ',num2str(vsv),str8])
        disp([str0,'Density of layer(s)  : ',num2str(rhov),str9])
        disp([str0,'Elastic tensor of layer(s)  : ',num2str(reshape(C,1,size(C,1)*size(C,2))),str10])
    elseif strcmp(wave_type,'viscoelastic')
        disp([str0,'Thickness of layer(s): ',num2str(tl),str6])
        disp([str0,'Vp of layer(s)       : ',num2str(vpv),str7])
        disp([str0,'Vs of layer(s)       : ',num2str(vsv),str8])
        disp([str0,'Density of layer(s)  : ',num2str(rhov),str9])        
        disp([str0,'Elastic tensor of layer(s)  : ',num2str(reshape(C,1,size(C,1)*size(C,2))),str10])
        disp([str0,'Qp of layer(s)       : ',num2str(qpv),str11])
        disp([str0,'Qs of layer(s)       : ',num2str(qsv),str12])
    end
end

tln = round(tl/dz);
nz = sum(tln)+1;
nx=round(XZ_ratio*nz);     % No of points or nodes along x AND depth,z (Including absorbing boundary nodes)
ny=round(XZ_ratio*nz);     % No of points or nodes along y AND depth,z (Including absorbing boundary nodes)

tl(1)=tl(1)+nab*dx;
tl(end)=tl(end)+nab*dx;
tln = round(tl/dz);
nz = sum(tln)+1;
nx=nx+2*nab;
ny=ny+2*nab;

zn = [1,cumsum(tln)+1];
vpm=zeros(nz,nx,ny);   vsm=zeros(nz,nx,ny);   qpm=zeros(nz,nx,ny);   qsm=zeros(nz,nx,ny);   rhom=zeros(nz,nx,ny);
C11=zeros(nz,nx,ny);    C12=zeros(nz,nx,ny);   C13=zeros(nz,nx,ny);   
C22=zeros(nz,nx,ny);    C23=zeros(nz,nx,ny);   C33=zeros(nz,nx,ny);  
C44=zeros(nz,nx,ny);    C55=zeros(nz,nx,ny);   C66=zeros(nz,nx,ny);
%
% for i=1:length(vpv)
%     vpm(zn(i):zn(i+1),:) = vpv(i);
%     vsm(zn(i):zn(i+1),:) = vsv(i);
%     qpm(zn(i):zn(i+1),:) = qpv(i);
%     qsm(zn(i):zn(i+1),:) = qsv(i);
%     rhom(zn(i):zn(i+1),:)= rhov(i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str=strcat(wfp,[filesep,'model']);
m_name='N_3Dlayers';

if strcmp(wave_type,'acoustic1')
    for i=1:length(vpv);        vpm(zn(i):zn(i+1),:,:) = vpv(i);          end
    save(str,'dx','dy','dz','nx','ny','nz','vpm','wave_type','m_name');
    
elseif strcmp(wave_type,'acoustic2')
    for i=1:length(vpv)
        vpm(zn(i):zn(i+1),:,:) = vpv(i);
        rhom(zn(i):zn(i+1),:,:)= rhov(i);
    end
    save(str,'dx','dy','dz','nx','ny','nz','vpm','rhom','wave_type','m_name');
    
elseif strcmp(wave_type,'elastic')
    for i=1:length(vpv)
        vpm(zn(i):zn(i+1),:,:) = vpv(i);
        vsm(zn(i):zn(i+1),:,:) = vsv(i);
        rhom(zn(i):zn(i+1),:,:)= rhov(i);
        C11(zn(i):zn(i+1),:,:)= C(i,1);
        C12(zn(i):zn(i+1),:,:)= C(i,2);
        C13(zn(i):zn(i+1),:,:)= C(i,3);
        C22(zn(i):zn(i+1),:,:)= C(i,4);
        C23(zn(i):zn(i+1),:,:)= C(i,5);
        C33(zn(i):zn(i+1),:,:)= C(i,6);
        C44(zn(i):zn(i+1),:,:)= C(i,7);
        C55(zn(i):zn(i+1),:,:)= C(i,8);
        C66(zn(i):zn(i+1),:,:)= C(i,9);
    end
    save(str,'dx','dy','dz','nx','ny','nz','vpm','vsm','rhom','C11','C12','C13','C22','C23','C33','C44','C55','C66','wave_type','m_name');
    
elseif strcmp(wave_type,'viscoelastic')
    for i=1:length(vpv)
        vpm(zn(i):zn(i+1),:,:) = vpv(i);
        vsm(zn(i):zn(i+1),:,:) = vsv(i);
        qpm(zn(i):zn(i+1),:,:) = qpv(i);
        qsm(zn(i):zn(i+1),:,:) = qsv(i);
        rhom(zn(i):zn(i+1),:,:)= rhov(i);        
        C11(zn(i):zn(i+1),:,:)= C(i,1);
        C12(zn(i):zn(i+1),:,:)= C(i,2);
        C13(zn(i):zn(i+1),:,:)= C(i,3);
        C22(zn(i):zn(i+1),:,:)= C(i,4);
        C23(zn(i):zn(i+1),:,:)= C(i,5);
        C33(zn(i):zn(i+1),:,:)= C(i,6);
        C44(zn(i):zn(i+1),:,:)= C(i,7);
        C55(zn(i):zn(i+1),:,:)= C(i,8);
        C66(zn(i):zn(i+1),:,:)= C(i,9);
    end
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

