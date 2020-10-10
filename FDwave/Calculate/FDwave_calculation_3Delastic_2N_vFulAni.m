function [Sx,Sy,Sz,wavefield,simtime]=FDwave_calculation_3Delastic_2N_vFulAni(varargin) %src_i,dN_W,dN_SS,plotON)

% The solution of wave equation is done with VECTORIZED FD operator.
% Description of parameters:
%        SRC_I    :  Source No for which the simulation will be carried out.
%        MT     :  Source mechanism in moment tensor
%        DN_W     :  No time steps of wavefield to skip
%        DN_SS    :  No time steps of synthetic seismogram to skip
%        PlotON   :  Show wave propagation while simulation
%        DN_P     :  No of time steps to skip plotting.

% Note:
% 1) It is advised to keep dN_W to very large value if you don't want to save entire wavefield
% 2) Plotting itself takes much time so it is advisable to skip few steps in case of very fine steps.
% 3) All other parameters are taken from the input folder.
% 4) The grid arrangement used in calculation is given following
%
%%%%%%%%%%%%%%%%%%%%%% Grid arrangement%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        txx,tzz__________ vx ___________ txx,tzz------>
%                lbd,mu  |                        bh                            |
%                             |                         |                             |
%                             |                         |                             |
%                             |                         |                             |
%                        vz  |____________txz                           |
%                       bv  |                       muvh                        |
%                             |                                                        |
%                             |                                                        |
%                             |                                                        |
%                 txx,tzz  |____________________________|
%                          |
%                         \|/
%

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'src_i';       src_i=varargin{i+1};
        case 'mt';         MT=varargin{i+1};
        case 'wfp';         wfp=varargin{i+1};
        case 'dn_w';        dN_W=varargin{i+1};
        case 'dn_ss';       dN_SS=varargin{i+1};
        case 'dn_p';        dN_P=varargin{i+1};
        case '2n';          M=varargin{i+1};
        case 'ploton';      plotON=varargin{i+1};
        case 'verbose';     verbose=varargin{i+1};
        case 'filename';    FileName=varargin{i+1};                  %WITH PATH
        case 'opname';      opname=varargin{i+1};
        otherwise
            error('%s is not a valid argument name',varargin{i});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Check IP parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('src_i','var')
    src_i=1   ;          
end

if ~exist('MT','var')
    MT=[1 0 0; 0 1 0; 0 0 1];             
end

if ~exist('wfp','var')
    wfp=pwd ;            
end

if ~exist('dN_W','var')
    load([wfp,[filesep,'source']],'N')
    dN_W=N+1;
end

if ~exist('dN_SS','var')
    dN_SS=1;             
end

if ~exist('M','var')
    M=4 ;            
end

if ~exist('dN_P','var')         
    dN_P=1;              
end

if ~exist('plotON','var')
    plotON='y' ;        
end

if ~exist('verbose','var')
    verbose='y';         
end

if ~exist('opname','var')
    opname='';    
else
    opname=[opname,'_'];
end

if strcmp(verbose,'y')
    disp('    FUNC: FD calculation begins (Elastic)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([wfp,[filesep,'model']],'dx','dy','dz','nx','ny','nz');

load([wfp,[filesep,'derived_param']],'mu','muxyz','lam','lamu','bx','by','bz');

load([wfp,[filesep,'source']],'dt','T','N','Ns','src');

load([wfp,[filesep,'BC']],'BC','BCname','BCtype');

load([wfp,[filesep,'geometry_src']],'geometry_src','snz_vec','snx_vec','sny_vec');

load([wfp,[filesep,'geometry_rec']],'geometry_rec','rec_n');


[vx,   vy,   vz,   txx,   tyy,   tzz,   txy,   txz,   tyz  ] = deal(zeros(nz,nx,ny));
[vx_x, vy_x, vz_x, txx_x, tyy_x, tzz_x, txy_x, txz_x, tyz_x] = deal(zeros(nz,nx,ny));
[vx_y, vy_y, vz_y, txx_y, tyy_y, tzz_y, txy_y, txz_y, tyz_y] = deal(zeros(nz,nx,ny));
[vx_z, vy_z, vz_z, txx_z, tyy_z, tzz_z, txy_z, txz_z, tyz_z] = deal(zeros(nz,nx,ny));

%% The 21 independent elastic parameters for full anisotropic models
[C11,   C12,   C13,   C14,   C15,   C16  ] = deal(zeros(nz,nx,ny));
[C22,   C23,   C24,   C25,   C26  ] = deal(zeros(nz,nx,ny));
[C33,   C34,   C35,   C36  ] = deal(zeros(nz,nx,ny));
[C44,   C45,   C46  ] = deal(zeros(nz,nx,ny));
[C55,   C56] = deal(zeros(nz,nx,ny));
[C66] = deal(zeros(nz,nx,ny));
%% loading elastic parameters for different types of anisotropic models
% here we use TI models as an example
load([wfp,[filesep,'model']]);
load([wfp,[filesep,'derived_param']]);

nframes=length(1:dN_W:N);      
wavefield=zeros(nz,nx,nframes);%ny,
nsteps=length(1:dN_SS:N);        
Sx = zeros(nsteps,rec_n);
Sy = zeros(nsteps,rec_n);
Sz = zeros(nsteps,rec_n);

xslice=snx_vec*dx;
yslice=sny_vec*dy;
zslice=snz_vec*dz;
[x,y,z] = meshgrid(dx*(0:nx-1),dy*(0:ny-1),dz*(0:nz-1));

%%%%%%%%%%%%%% FD operator of M=2N order
NN=M/2;
fdc=DiffCoef(NN,'s');% coefficients of FD operator
switch M
    case 4
        Dxfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj+1,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj+2,kk)-a(ii,jj-1,kk)) );
        Dyfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk+1)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj,kk+2)-a(ii,jj,kk-1)) );
        Dzfm=@(a,ii, jj, kk)( fdc(1)*(a(ii+1,jj,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii+2,jj,kk)-a(ii-1,jj,kk)) );
        Dxbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj-1,kk))  +  fdc(2)*(a(ii,jj+1,kk)-a(ii,jj-2,kk)) );
        Dybm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj,kk-1))  +  fdc(2)*(a(ii,jj,kk+1)-a(ii,jj,kk-2)) );
        Dzbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii-1,jj,kk))  +  fdc(2)*(a(ii+1,jj,kk)-a(ii-2,jj,kk)) );
            
    case 6
        Dxfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj+1,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj+2,kk)-a(ii,jj-1,kk)) ...
            +  fdc(3)*(a(ii,jj+3,kk)-a(ii,jj-2,kk)));
        Dyfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk+1)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj,kk+2)-a(ii,jj,kk-1)) ...
            +  fdc(3)*(a(ii,jj,kk+3)-a(ii,jj,kk-2)));
        Dzfm=@(a,ii, jj, kk)( fdc(1)*(a(ii+1,jj,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii+2,jj,kk)-a(ii-1,jj,kk)) ...
            +  fdc(3)*(a(ii+3,jj,kk)-a(ii-2,jj,kk)));
        Dxbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj-1,kk))  +  fdc(2)*(a(ii,jj+1,kk)-a(ii,jj-2,kk)) ...
            +  fdc(3)*(a(ii,jj+2,kk)-a(ii,jj-3,kk)));
        Dybm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj,kk-1))  +  fdc(2)*(a(ii,jj,kk+1)-a(ii,jj,kk-2)) ...
            +  fdc(3)*(a(ii,jj,kk+2)-a(ii,jj,kk-3)));
        Dzbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii-1,jj,kk))  +  fdc(2)*(a(ii+1,jj,kk)-a(ii-2,jj,kk)) ...
            +  fdc(3)*(a(ii+2,jj,kk)-a(ii-3,jj,kk)));
    
    case 8
        Dxfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj+1,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj+2,kk)-a(ii,jj-1,kk)) ...
            +  fdc(3)*(a(ii,jj+3,kk)-a(ii,jj-2,kk))+  fdc(4)*(a(ii,jj+4,kk)-a(ii,jj-3,kk)));
        Dyfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk+1)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj,kk+2)-a(ii,jj,kk-1)) ...
            +  fdc(3)*(a(ii,jj,kk+3)-a(ii,jj,kk-2))+  fdc(4)*(a(ii,jj,kk+4)-a(ii,jj,kk-3)));
        Dzfm=@(a,ii, jj, kk)( fdc(1)*(a(ii+1,jj,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii+2,jj,kk)-a(ii-1,jj,kk)) ...
            +  fdc(3)*(a(ii+3,jj,kk)-a(ii-2,jj,kk))+  fdc(4)*(a(ii+4,jj,kk)-a(ii-3,jj,kk)));
        Dxbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj-1,kk))  +  fdc(2)*(a(ii,jj+1,kk)-a(ii,jj-2,kk)) ...
            +  fdc(3)*(a(ii,jj+2,kk)-a(ii,jj-3,kk))+  fdc(4)*(a(ii,jj+3,kk)-a(ii,jj-4,kk)));
        Dybm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj,kk-1))  +  fdc(2)*(a(ii,jj,kk+1)-a(ii,jj,kk-2)) ...
            +  fdc(3)*(a(ii,jj,kk+2)-a(ii,jj,kk-3))+  fdc(4)*(a(ii,jj,kk+3)-a(ii,jj,kk-4)));
        Dzbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii-1,jj,kk))  +  fdc(2)*(a(ii+1,jj,kk)-a(ii-2,jj,kk)) ...
            +  fdc(3)*(a(ii+2,jj,kk)-a(ii-3,jj,kk))+  fdc(4)*(a(ii+3,jj,kk)-a(ii-4,jj,kk)));
    
    case 10
        Dxfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj+1,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj+2,kk)-a(ii,jj-1,kk)) ...
            +  fdc(3)*(a(ii,jj+3,kk)-a(ii,jj-2,kk))+  fdc(4)*(a(ii,jj+4,kk)-a(ii,jj-3,kk))+  fdc(5)*(a(ii,jj+5,kk)-a(ii,jj-4,kk))) ;
        Dyfm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk+1)-a(ii,jj,kk))  +  fdc(2)*(a(ii,jj,kk+2)-a(ii,jj,kk-1)) ...
            +  fdc(3)*(a(ii,jj,kk+3)-a(ii,jj,kk-2))+  fdc(4)*(a(ii,jj,kk+4)-a(ii,jj,kk-3))+  fdc(5)*(a(ii,jj,kk+5)-a(ii,jj,kk-4))) ;
        Dzfm=@(a,ii, jj, kk)( fdc(1)*(a(ii+1,jj,kk)-a(ii,jj,kk))  +  fdc(2)*(a(ii+2,jj,kk)-a(ii-1,jj,kk)) ...
            +  fdc(3)*(a(ii+3,jj,kk)-a(ii-2,jj,kk))+  fdc(4)*(a(ii+4,jj,kk)-a(ii-3,jj,kk))+  fdc(5)*(a(ii+5,jj,kk)-a(ii-4,jj,kk))) ;
        Dxbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj-1,kk))  +  fdc(2)*(a(ii,jj+1,kk)-a(ii,jj-2,kk)) ...
            +  fdc(3)*(a(ii,jj+2,kk)-a(ii,jj-3,kk))+  fdc(4)*(a(ii,jj+3,kk)-a(ii,jj-4,kk))+  fdc(5)*(a(ii,jj+4,kk)-a(ii,jj-5,kk))) ;
        Dybm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii,jj,kk-1))  +  fdc(2)*(a(ii,jj,kk+1)-a(ii,jj,kk-2)) ...
            +  fdc(3)*(a(ii,jj,kk+2)-a(ii,jj,kk-3))+  fdc(4)*(a(ii,jj,kk+3)-a(ii,jj,kk-4))+  fdc(5)*(a(ii,jj,kk+4)-a(ii,jj,kk-5))) ;
        Dzbm=@(a,ii, jj, kk)( fdc(1)*(a(ii,jj,kk)-a(ii-1,jj,kk))  +  fdc(2)*(a(ii+1,jj,kk)-a(ii-2,jj,kk)) ...
            +  fdc(3)*(a(ii+2,jj,kk)-a(ii-3,jj,kk))+  fdc(4)*(a(ii+3,jj,kk)-a(ii-4,jj,kk))+  fdc(5)*(a(ii+4,jj,kk)-a(ii-5,jj,kk))) ;
        
end


h=NN+1;
ii = (NN+1):(nz-NN);
jj = (NN+1):(nx-NN);
kk = (NN+1):(ny-NN);
    
Pmlxd = (1+0.5*dt*BC.pmlx(ii, jj, kk)); 
Pmlyd = (1+0.5*dt*BC.pmly(ii, jj, kk)); 
Pmlzd = (1+0.5*dt*BC.pmlz(ii, jj, kk)); 
Pmlxn = (1-0.5*dt*BC.pmlx(ii, jj, kk));
Pmlyn = (1-0.5*dt*BC.pmly(ii, jj, kk));
Pmlzn = (1-0.5*dt*BC.pmlz(ii, jj, kk));

srcnz= snz_vec(src_i);
srcnx= snx_vec(src_i);
srcny= sny_vec(src_i);

tic
h_wt = waitbar(0,' completed...','Name','FDwave_3D');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Simulation starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(plotON,'y')
    figure;
end
    
for t=1:N -3              
    
    %% PML boundary
    if strcmp(BCname,'PML')  
    %% source types
    % The source implementation details can be found in the section of 
    % Moment tensor source implementation  eqs. (10)(11) in the related manuscript 
    if t<=Ns
        % four nodes
        %{
    txx_x(srcnz,srcnx,srcny)=txx_x(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);
    txx_y(srcnz,srcnx,srcny)=txx_y(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);  
    txx_z(srcnz,srcnx,srcny)=txx_z(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);  
    
    tyy_x(srcnz,srcnx,srcny)=tyy_x(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);  
    tyy_y(srcnz,srcnx,srcny)=tyy_y(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);  
    tyy_z(srcnz,srcnx,srcny)=tyy_z(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);   
    
    tzz_x(srcnz,srcnx,srcny)=tzz_x(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    tzz_y(srcnz,srcnx,srcny)=tzz_y(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    tzz_z(srcnz,srcnx,srcny)=tzz_z(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    
    txy_x(srcnz,srcnx,srcny)=txy_x(srcnz,srcnx,srcny)+(-MT(1,2)/8)*src(t);
    txy_x(srcnz,srcnx-1,srcny)=txy_x(srcnz,srcnx-1,srcny)+(-MT(1,2)/8)*src(t);
    txy_x(srcnz,srcnx,srcny-1)=txy_x(srcnz,srcnx,srcny-1)+(-MT(1,2)/8)*src(t);
    txy_x(srcnz,srcnx-1,srcny-1)=txy_x(srcnz,srcnx-1,srcny-1)+(-MT(1,2)/8)*src(t);
    txy_y(srcnz,srcnx,srcny)=txy_y(srcnz,srcnx,srcny)+(-MT(1,2)/8)*src(t);
    txy_y(srcnz,srcnx-1,srcny)=txy_y(srcnz,srcnx-1,srcny)+(-MT(1,2)/8)*src(t);
    txy_y(srcnz,srcnx,srcny-1)=txy_y(srcnz,srcnx,srcny-1)+(-MT(1,2)/8)*src(t);
    txy_y(srcnz,srcnx-1,srcny-1)=txy_y(srcnz,srcnx-1,srcny-1)+(-MT(1,2)/8)*src(t);
    
    txz_x(srcnz,srcnx,srcny)=txz_x(srcnz,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz_x(srcnz-1,srcnx,srcny)=txz_x(srcnz-1,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz_x(srcnz,srcnx-1,srcny)=txz_x(srcnz,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    txz_x(srcnz-1,srcnx-1,srcny)=txz_x(srcnz-1,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    txz_z(srcnz,srcnx,srcny)=txz_z(srcnz,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz_z(srcnz-1,srcnx,srcny)=txz_z(srcnz-1,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz_z(srcnz,srcnx-1,srcny)=txz_z(srcnz,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    txz_z(srcnz-1,srcnx-1,srcny)=txz_z(srcnz-1,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    
    tyz_y(srcnz,srcnx,srcny)=tyz_y(srcnz,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz_y(srcnz-1,srcnx,srcny)=tyz_y(srcnz-1,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz_y(srcnz,srcnx,srcny-1)=tyz_y(srcnz,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
    tyz_y(srcnz-1,srcnx,srcny-1)=tyz_y(srcnz-1,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
    tyz_z(srcnz,srcnx,srcny)=tyz_z(srcnz,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz_z(srcnz-1,srcnx,srcny)=tyz_z(srcnz-1,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz_z(srcnz,srcnx,srcny-1)=tyz_z(srcnz,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
    tyz_z(srcnz-1,srcnx,srcny-1)=tyz_z(srcnz-1,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
        
        
    txx=txx_x+txx_y+txx_z;
    tyy=tyy_x+tyy_y+tyy_z;
    tzz=tzz_x+tzz_y+tzz_z;
    txy=txy_x+txy_y;
    txz=txz_x+txz_z;
    tyz=tyz_y+tyz_z;  
    %}
    
        % one node
    txx_x(srcnz,srcnx,srcny)=txx_x(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);
    txx_y(srcnz,srcnx,srcny)=txx_y(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);  
    txx_z(srcnz,srcnx,srcny)=txx_z(srcnz,srcnx,srcny)+(-MT(1,1)/3)*src(t);  
    
    tyy_x(srcnz,srcnx,srcny)=tyy_x(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);  
    tyy_y(srcnz,srcnx,srcny)=tyy_y(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);  
    tyy_z(srcnz,srcnx,srcny)=tyy_z(srcnz,srcnx,srcny)+(-MT(2,2)/3)*src(t);   
    
    tzz_x(srcnz,srcnx,srcny)=tzz_x(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    tzz_y(srcnz,srcnx,srcny)=tzz_y(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    tzz_z(srcnz,srcnx,srcny)=tzz_z(srcnz,srcnx,srcny)+(-MT(3,3)/3)*src(t);
    
    txy_x(srcnz,srcnx,srcny)=txy_x(srcnz,srcnx,srcny)+(-MT(1,2)/3)*src(t);
    txy_y(srcnz,srcnx,srcny)=txy_y(srcnz,srcnx,srcny)+(-MT(1,2)/3)*src(t);
    txy_z(srcnz,srcnx,srcny)=txy_z(srcnz,srcnx,srcny)+(-MT(1,2)/3)*src(t);
    
    txz_x(srcnz,srcnx,srcny)=txz_x(srcnz,srcnx,srcny)+(-MT(1,3)/3)*src(t);
    txz_y(srcnz,srcnx,srcny)=txz_y(srcnz,srcnx,srcny)+(-MT(1,3)/3)*src(t);
    txz_z(srcnz,srcnx,srcny)=txz_z(srcnz,srcnx,srcny)+(-MT(1,3)/3)*src(t);
    
    tyz_x(srcnz,srcnx,srcny)=tyz_x(srcnz,srcnx,srcny)+(-MT(2,3)/3)*src(t);
    tyz_y(srcnz,srcnx,srcny)=tyz_y(srcnz,srcnx,srcny)+(-MT(2,3)/3)*src(t);
    tyz_z(srcnz,srcnx,srcny)=tyz_z(srcnz,srcnx,srcny)+(-MT(2,3)/3)*src(t);
        
        
    txx=txx_x+txx_y+txx_z;
    tyy=tyy_x+tyy_y+tyy_z;
    tzz=tzz_x+tzz_y+tzz_z;
    txy=txy_x+txy_y+txy_z;
    txz=txz_x+txz_y+txz_z;
    tyz=tyz_x+tyz_y+tyz_z;   
       
    end
    
   
    %%%%%%%%%%%%%%%%%% velocity component %%%%%%%%%%%%%%%%%
    % corresponds to eq. (3) in the manuscript
    %vx(3:nz-2, 3:nx-2, 3:ny-2) = vx(3:nz-2, 3:nx-2, 3:ny-2) +  dt*bx(3:nz-2, 3:nx-2, 3:ny-2).*( Dxfm(txx,nx,ny,nz)/dx + Dzbm(txz,nx,ny,nz)/dz  + Dybm(txy,nx,ny,nz)/dy);
    vx_x(ii, jj, kk) = (Pmlxn.*vx_x(ii, jj, kk) +  dt*bx(ii, jj, kk).* Dxfm(txx,ii, jj, kk)/dx)./Pmlxd;
    vx_y(ii, jj, kk) = (Pmlyn.*vx_y(ii, jj, kk) +  dt*bx(ii, jj, kk).* Dybm(txy,ii, jj, kk)/dy)./Pmlyd;
    vx_z(ii, jj, kk) = (Pmlzn.*vx_z(ii, jj, kk) +  dt*bx(ii, jj, kk).* Dzbm(txz,ii, jj, kk)/dz)./Pmlzd;
    
    %vy(3:nz-2, 3:nx-2, 3:ny-2) = vy(3:nz-2, 3:nx-2, 3:ny-2) +  dt*by(3:nz-2, 3:nx-2, 3:ny-2).*( Dxbm(txy,nx,ny,nz)/dx + Dzbm(tyz,nx,ny,nz)/dz  + Dyfm(tyy,nx,ny,nz)/dy);
    vy_x(ii, jj, kk) = (Pmlxn.*vy_x(ii, jj, kk) +  dt*by(ii, jj, kk).*Dxbm(txy,ii, jj, kk)/dx)./Pmlxd;
    vy_y(ii, jj, kk) = (Pmlyn.*vy_y(ii, jj, kk) +  dt*by(ii, jj, kk).*Dyfm(tyy,ii, jj, kk)/dy)./Pmlyd;
    vy_z(ii, jj, kk) = (Pmlzn.*vy_z(ii, jj, kk) +  dt*by(ii, jj, kk).*Dzbm(tyz,ii, jj, kk)/dz)./Pmlzd;
    
    %vz(3:nz-2, 3:nx-2, 3:ny-2) = vz(3:nz-2, 3:nx-2, 3:ny-2) +  dt*bz(3:nz-2, 3:nx-2, 3:ny-2).*( Dxbm(txz,nx,ny,nz)/dx + Dzfm(tzz,nx,ny,nz)/dz  + Dybm(tyz,nx,ny,nz)/dy);
    vz_x(ii, jj, kk) =  (Pmlxn.*vz_x(ii, jj, kk) +  dt*bz(ii, jj, kk).* Dxbm(txz,ii, jj, kk)/dx)./Pmlxd;
    vz_y(ii, jj, kk) =  (Pmlyn.*vz_y(ii, jj, kk) +  dt*bz(ii, jj, kk).* Dybm(tyz,ii, jj, kk)/dy)./Pmlyd;
    vz_z(ii, jj, kk) =  (Pmlzn.*vz_z(ii, jj, kk) +  dt*bz(ii, jj, kk).* Dzfm(tzz,ii, jj, kk)/dz)./Pmlzd;
    
    
    vx=vx_x+vx_y+vx_z;
    vy=vy_x+vy_y+vy_z;
    vz=vz_x+vz_y+vz_z;
    
    % topFS with the assumption of weak anisotropy near the surface
    if strcmp(BCtype,'topFS')
            vz(h-1,NN+1:nx-NN, NN+1:ny-NN) = vz(h,NN+1:nx-NN, NN+1:ny-NN) +... 
                (lam(h,NN+1:nx-NN, NN+1:ny-NN)./lamu(h,NN+1:nx-NN, NN+1:ny-NN)).*(...
                vx(h,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,(NN+1:nx-NN)-1, NN+1:ny-NN)...
                +vx(h,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, (NN+1:ny-NN)-1)   ); 
            %vx
            %vy
             vx(h-1,NN+1:nx-NN, NN+1:ny-NN) = vx(h,NN+1:nx-NN, NN+1:ny-NN) ... 
                +vz(h-1,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h-1,NN+1:nx-NN, NN+1:ny-NN)...
                +vz(h  ,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h  ,NN+1:nx-NN, NN+1:ny-NN)...
                +vx(h+1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, NN+1:ny-NN)  ;
            
            vy(h-1,NN+1:nx-NN, NN+1:ny-NN) = vy(h,NN+1:nx-NN, NN+1:ny-NN) ... 
                +vz(h-1,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h-1,NN+1:nx-NN, NN+1:ny-NN)...
                +vz(h  ,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h  ,NN+1:nx-NN, NN+1:ny-NN)...
                +vx(h+1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, NN+1:ny-NN)  ;
                
            vz(h-2,NN+1:nx-NN, NN+1:ny-NN) = vz(h-1,NN+1:nx-NN, NN+1:ny-NN) +... 
                (lam(h,NN+1:nx-NN, NN+1:ny-NN)./lamu(h,NN+1:nx-NN, NN+1:ny-NN)).*(...
                 vx(h-1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h-1,(NN+1:nx-NN)-1, NN+1:ny-NN)...
                +vy(h-1,NN+1:nx-NN, NN+1:ny-NN)  -  vy(h-1,NN+1:nx-NN, (NN+1:ny-NN)-1)   ); 
           
    end
    
    
    
    %%%%%%%%%%%%%%%%%% stress component %%%%%%%%%%%%%%%%%
    % corresponds to eq. (4) in the manuscript

    %txx(3:nz-2, 3:nx-2, 3:ny-2) = txx(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*( lamu(3:nz-2, 3:nx-2, 3:ny-2).*dxvx + lam(3:nz-2, 3:nx-2, 3:ny-2).*(dyvy + dzvz) );
    txx_x(ii, jj, kk) =  (Pmlxn.*txx_x(ii, jj, kk)  +  dt* (C11(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C16(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C15(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)./Pmlxd;
    txx_y(ii, jj, kk) =  (Pmlyn.*txx_y(ii, jj, kk)  +  dt* (C16(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C12(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C14(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)./Pmlyd;
    txx_z(ii, jj, kk) =  (Pmlzn.*txx_z(ii, jj, kk)  +  dt* (C15(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C14(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C13(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    %tyy(3:nz-2, 3:nx-2, 3:ny-2) = tyy(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*( lam(3:nz-2, 3:nx-2, 3:ny-2).*(dxvx + dzvz) + lamu(3:nz-2, 3:nx-2, 3:ny-2).*dyvy );
    tyy_x(ii, jj, kk) = (Pmlxn.*tyy_x(ii, jj, kk)  +  dt* (C12(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C26(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C25(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)./Pmlxd;
    tyy_y(ii, jj, kk) = (Pmlyn.*tyy_y(ii, jj, kk)  +  dt* (C26(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C22(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C24(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)./Pmlyd;
    tyy_z(ii, jj, kk) = (Pmlzn.*tyy_z(ii, jj, kk)  +  dt* (C25(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C24(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C23(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    %tzz(3:nz-2, 3:nx-2, 3:ny-2) = tzz(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*( lam(3:nz-2, 3:nx-2, 3:ny-2).*(dxvx + dyvy)+ lamu(3:nz-2, 3:nx-2, 3:ny-2).*dzvz );
    tzz_x(ii, jj, kk) = (Pmlxn.*tzz_x(ii, jj, kk)  +  dt* (C13(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C36(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C35(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)./Pmlxd;
    tzz_y(ii, jj, kk) = (Pmlyn.*tzz_y(ii, jj, kk)  +  dt* (C36(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C23(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C34(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)./Pmlyd;
    tzz_z(ii, jj, kk) = (Pmlzn.*tzz_z(ii, jj, kk)  +  dt* (C35(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C34(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C33(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    %txy(3:nz-2, 3:nx-2, 3:ny-2) = txy(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*muxyz(3:nz-2, 3:nx-2, 3:ny-2).*(  Dyfm(vx,nx,ny,nz)/dy + Dxfm(vy,nx,ny,nz)/dx   );
    txy_x(ii, jj, kk) = (Pmlxn.*txy_x(ii, jj, kk)  +  dt*(C16(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C66(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C56(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)./Pmlxd;
    txy_y(ii, jj, kk) = (Pmlyn.*txy_y(ii, jj, kk)  +  dt*(C66(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C26(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C46(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)./Pmlyd;
    txy_z(ii, jj, kk) = (Pmlzn.*txy_z(ii, jj, kk)  +  dt*(C56(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C46(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C36(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    %txz(3:nz-2, 3:nx-2, 3:ny-2) = txz(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*muxyz(3:nz-2, 3:nx-2, 3:ny-2).*(  Dzfm(vx,nx,ny,nz)/dz + Dxfm(vz,nx,ny,nz)/dx   );
    txz_x(ii, jj, kk) = (Pmlxn.*txz_x(ii, jj, kk)  +  dt*(C15(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C56(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C55(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)./Pmlxd;
    txz_y(ii, jj, kk) = (Pmlyn.*txz_y(ii, jj, kk)  +  dt*(C56(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C25(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C45(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)./Pmlyd;
    txz_z(ii, jj, kk) = (Pmlzn.*txz_z(ii, jj, kk)  +  dt*(C55(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C45(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C35(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    %tyz(3:nz-2, 3:nx-2, 3:ny-2) = tyz(3:nz-2, 3:nx-2, 3:ny-2)  +  dt*muxyz(3:nz-2, 3:nx-2, 3:ny-2).*(  Dzfm(vy,nx,ny,nz)/dz + Dyfm(vz,nx,ny,nz)/dy   );
    tyz_x(ii, jj, kk) = (Pmlxn.*tyz_x(ii, jj, kk)  +  dt*(C14(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C46(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C45(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)./Pmlxd;
    tyz_y(ii, jj, kk) = (Pmlyn.*tyz_y(ii, jj, kk)  +  dt*(C46(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C24(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C44(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)./Pmlyd;
    tyz_z(ii, jj, kk) = (Pmlzn.*tyz_z(ii, jj, kk)  +  dt*(C45(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C44(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C34(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz)./Pmlzd;
    
    txx=txx_x+txx_y+txx_z;
    tyy=tyy_x+tyy_y+tyy_z;
    tzz=tzz_x+tzz_y+tzz_z;
    txy=txy_x+txy_y+txy_z;
    txz=txz_x+txz_y+txz_z;
    tyz=tyz_x+tyz_y+tyz_z;      
    
   % topFS with the assumption of weak anisotropy near the surface
    if strcmpi(BCtype,'topFS')
            tzz(h,:,:)=0;
            tzz(h-1,:,:)=-tzz(h+1,:,:);
            txz(h-1,:,:)=-txz(h,:,:);
            txz(h-2,:,:)=-txz(h+1,:,:);
            tyz(h-1,:,:)=-tyz(h,:,:);
            tyz(h-2,:,:)=-tyz(h+1,:,:);
            %thh(h,3:nx-2)= thh(h,3:nx-2) + dt*FScoeff(h, 3:nx-2).*Dxbv(vh(h,:),nx);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% ABL boundary
    else
        
    % source types
    % The source implementation details can be found in the section of 
    % Moment tensor source implementation  eqs. (10)(11) in the related manuscript 
    if t<=Ns
        % four nodes
        %{
    txx(srcnz,srcnx,srcny)=txx(srcnz,srcnx,srcny)+(-MT(1,1))*src(t);
    tyy(srcnz,srcnx,srcny)=tyy(srcnz,srcnx,srcny)+(-MT(2,2))*src(t);   
    tzz(srcnz,srcnx,srcny)=tzz(srcnz,srcnx,srcny)+(-MT(3,3))*src(t);
    
    txy(srcnz,srcnx,srcny)=txy(srcnz,srcnx,srcny)+(-MT(1,2)/8)*src(t);
    txy(srcnz,srcnx-1,srcny)=txy(srcnz,srcnx-1,srcny)+(-MT(1,2)/8)*src(t);
    txy(srcnz,srcnx,srcny-1)=txy(srcnz,srcnx,srcny-1)+(-MT(1,2)/8)*src(t);
    txy(srcnz,srcnx-1,srcny-1)=txy(srcnz,srcnx-1,srcny-1)+(-MT(1,2)/8)*src(t);
   
    txz(srcnz,srcnx,srcny)=txz(srcnz,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz(srcnz-1,srcnx,srcny)=txz(srcnz-1,srcnx,srcny)+(-MT(1,3)/8)*src(t);
    txz(srcnz,srcnx-1,srcny)=txz(srcnz,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    txz(srcnz-1,srcnx-1,srcny)=txz(srcnz-1,srcnx-1,srcny)+(-MT(1,3)/8)*src(t);
    
    tyz(srcnz,srcnx,srcny)=tyz(srcnz,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz(srcnz-1,srcnx,srcny)=tyz(srcnz-1,srcnx,srcny)+(-MT(2,3)/8)*src(t);
    tyz(srcnz,srcnx,srcny-1)=tyz(srcnz,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
    tyz(srcnz-1,srcnx,srcny-1)=tyz(srcnz-1,srcnx,srcny-1)+(-MT(2,3)/8)*src(t);
     
    %}
    
        % one node
    txx(srcnz,srcnx,srcny)=txx(srcnz,srcnx,srcny)+(-MT(1,1))*src(t);
    tyy(srcnz,srcnx,srcny)=tyy(srcnz,srcnx,srcny)+(-MT(2,2))*src(t);   
    tzz(srcnz,srcnx,srcny)=tzz(srcnz,srcnx,srcny)+(-MT(3,3))*src(t);
    
    txy(srcnz,srcnx,srcny)=txy(srcnz,srcnx,srcny)+(-MT(1,2))*src(t);
    txz(srcnz,srcnx,srcny)=txz(srcnz,srcnx,srcny)+(-MT(1,3))*src(t);
    tyz(srcnz,srcnx,srcny)=tyz(srcnz,srcnx,srcny)+(-MT(2,3))*src(t);
        
        
       
    end
    
    %%%%%%%%%%%%%%%%%% velocity component %%%%%%%%%%%%%%%%%
    % corresponds to eq. (3) in the manuscript
    vx(ii, jj, kk) = vx(ii, jj, kk) +  dt*bx(ii, jj, kk).*( Dxfm(txx,ii, jj, kk)/dx + Dzbm(txz,ii, jj, kk)/dz  + Dybm(txy,ii, jj, kk)/dy);
    vy(ii, jj, kk) = vy(ii, jj, kk) +  dt*by(ii, jj, kk).*( Dxbm(txy,ii, jj, kk)/dx + Dzbm(tyz,ii, jj, kk)/dz  + Dyfm(tyy,ii, jj, kk)/dy);
    vz(ii, jj, kk) = vz(ii, jj, kk) +  dt*bz(ii, jj, kk).*( Dxbm(txz,ii, jj, kk)/dx + Dzfm(tzz,ii, jj, kk)/dz  + Dybm(tyz,ii, jj, kk)/dy);
    
    % topFS with the assumption of weak anisotropy near the surface
    if strcmp(BCtype,'topFS')
              vz(h-1,NN+1:nx-NN, NN+1:ny-NN) = vz(h,NN+1:nx-NN, NN+1:ny-NN) +... 
                (lam(h,NN+1:nx-NN, NN+1:ny-NN)./lamu(h,NN+1:nx-NN, NN+1:ny-NN)).*(...
                vx(h,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,(NN+1:nx-NN)-1, NN+1:ny-NN)...
                +vx(h,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, (NN+1:ny-NN)-1)   ); 
            %vx
            %vy
             vx(h-1,NN+1:nx-NN, NN+1:ny-NN) = vx(h,NN+1:nx-NN, NN+1:ny-NN) ... 
                +vz(h-1,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h-1,NN+1:nx-NN, NN+1:ny-NN)...
                +vz(h  ,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h  ,NN+1:nx-NN, NN+1:ny-NN)...
                +vx(h+1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, NN+1:ny-NN)  ;
            
            vy(h-1,NN+1:nx-NN, NN+1:ny-NN) = vy(h,NN+1:nx-NN, NN+1:ny-NN) ... 
                +vz(h-1,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h-1,NN+1:nx-NN, NN+1:ny-NN)...
                +vz(h  ,(NN+1:nx-NN)+1, NN+1:ny-NN)  -  vz(h  ,NN+1:nx-NN, NN+1:ny-NN)...
                +vx(h+1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h,NN+1:nx-NN, NN+1:ny-NN)  ;
                
            vz(h-2,NN+1:nx-NN, NN+1:ny-NN) = vz(h-1,NN+1:nx-NN, NN+1:ny-NN) +... 
                (lam(h,NN+1:nx-NN, NN+1:ny-NN)./lamu(h,NN+1:nx-NN, NN+1:ny-NN)).*(...
                 vx(h-1,NN+1:nx-NN, NN+1:ny-NN)  -  vx(h-1,(NN+1:nx-NN)-1, NN+1:ny-NN)...
                +vy(h-1,NN+1:nx-NN, NN+1:ny-NN)  -  vy(h-1,NN+1:nx-NN, (NN+1:ny-NN)-1)   ); 
           
    end
    
    
    vx=BC.abl.*vx;   		  
    vy=BC.abl.*vy;   
    vz=BC.abl.*vz; 
    
    %%%%%%%%%%%%%%%%%% stress component %%%%%%%%%%%%%%%%%
    % corresponds to eq. (4) in the manuscript
    
    txx(ii, jj, kk) = txx(ii, jj, kk)  +  dt*( ((C11(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C16(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C15(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)+...
        ((C16(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C12(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C14(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)+...
        ((C15(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C14(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C13(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz) );
    tyy(ii, jj, kk) = tyy(ii, jj, kk)  +  dt*( ((C12(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C26(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C25(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)+...
        ((C26(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C22(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C24(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)+...
        ((C25(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C24(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C23(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz) );
    tzz(ii, jj, kk) = tzz(ii, jj, kk)  +  dt*( ((C13(ii, jj, kk).*Dxbm(vx,ii, jj, kk)+C36(ii, jj, kk).*Dxbm(vy,ii, jj, kk)+C35(ii, jj, kk).*Dxbm(vz,ii, jj, kk))/dx)+...
        ((C36(ii, jj, kk).*Dybm(vx,ii, jj, kk)+C23(ii, jj, kk).*Dybm(vy,ii, jj, kk)+C34(ii, jj, kk).*Dybm(vz,ii, jj, kk))/dy)+...
        ((C35(ii, jj, kk).*Dzbm(vx,ii, jj, kk)+C34(ii, jj, kk).*Dzbm(vy,ii, jj, kk)+C33(ii, jj, kk).*Dzbm(vz,ii, jj, kk))/dz) );
    txy(ii, jj, kk) = txy(ii, jj, kk)  +  dt*( ((C16(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C66(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C56(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)+...
        ((C66(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C26(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C46(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)+...
        ((C56(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C46(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C36(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz) );
    txz(ii, jj, kk) = txz(ii, jj, kk)  +  dt*( ((C15(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C56(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C55(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)+...
        ((C56(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C25(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C45(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)+...
        ((C55(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C45(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C35(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz) );
    tyz(ii, jj, kk) = tyz(ii, jj, kk)  +  dt*( ((C14(ii, jj, kk).*Dxfm(vx,ii, jj, kk)+C46(ii, jj, kk).*Dxfm(vy,ii, jj, kk)+C45(ii, jj, kk).*Dxfm(vz,ii, jj, kk))/dx)+...
        ((C46(ii, jj, kk).*Dyfm(vx,ii, jj, kk)+C24(ii, jj, kk).*Dyfm(vy,ii, jj, kk)+C44(ii, jj, kk).*Dyfm(vz,ii, jj, kk))/dy)+...
        ((C45(ii, jj, kk).*Dzfm(vx,ii, jj, kk)+C44(ii, jj, kk).*Dzfm(vy,ii, jj, kk)+C34(ii, jj, kk).*Dzfm(vz,ii, jj, kk))/dz) );
    
    % topFS with the assumption of weak anisotropy near the surface 
    if strcmpi(BCtype,'topFS')
            tzz(h,:,:)=0;
            tzz(h-1,:,:)=-tzz(h+1,:,:);
            txz(h-1,:,:)=-txz(h,:,:);
            txz(h-2,:,:)=-txz(h+1,:,:);
            tyz(h-1,:,:)=-tyz(h,:,:);
            tyz(h-2,:,:)=-tyz(h+1,:,:);
            %thh(h,3:nx-2)= thh(h,3:nx-2) + dt*FScoeff(h, 3:nx-2).*Dxbv(vh(h,:),nx);
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    txx=BC.abl.*txx;        
    tyy=BC.abl.*tyy;             
    tzz=BC.abl.*tzz;  
    txy=BC.abl.*txy;        
    txz=BC.abl.*txz;             
    tyz=BC.abl.*tyz;
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(plotON,'y')&&(mod(t-1,dN_P)==0)
        
        subplot(2,1,1);     
		slicen(x,y,z,permute(vz,[3,2,1]),xslice,yslice,zslice);
        %imagesc((0:nx-1)*dx,(0:nz-1)*dz,vz(:,:,round(ny/2)));
        shading interp;
        set(gca,'fontsize',14,'ydir','reverse','zdir','reverse');
        xlabel('X (m)');ylabel('Y (m)');
        zlabel('Z (m)');
        colormap(flipud(jet))
        axis([0 dx*(nx-1) 0 dy*(ny-1) 0 dz*(nz-1)])
		title(['Vertical velocity at time =',num2str(dt*(t-1),'%3.6f'),'s'])
        %caxis([cmin cmax]);
        %set(gca,'XAxisLocation','Top');
        
        subplot(2,1,2);     
		imagesc(1:rec_n,(1:(N-1))*dt,Sz);
        set(gca,'fontsize',14);
        xlabel('Receiver Number');  
		ylabel('T (s)');
        title(['Synthetic seismogram at time =',num2str(dt*(t-1),'%3.6f'),'s'])
        set(gca,'XAxisLocation','Top')
        %caxis([cmin cmax]);
        %colormap(flipud(gray));
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mod(t-1,dN_W)==0)
        dN_Wi = round(t/dN_W);
        wavefield(:,:,dN_Wi+1)= vy(:,:,round(ny/2));
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mod(t-1,dN_SS)==0)
        dN_SSi = round(t/dN_SS);
        Sx(dN_SSi+1,:)= vx(geometry_rec);
        Sy(dN_SSi+1,:)= vy(geometry_rec);
        Sz(dN_SSi+1,:)= vz(geometry_rec);
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perc=round(t/N*100);
    waitbar(perc/100,h_wt,sprintf('%d%% completed...',perc));
    
end

close(h_wt)
simtime=toc;
str=strcat(wfp,[filesep,'SS_'],opname,num2str(src_i));
save(str,'Sx','Sy','Sz','dN_SS','dN_W','dt','N','M','dx','dy','dz','nx','ny','nz','-v7.3');

str=strcat(wfp,[filesep,'wavefield_'],opname,num2str(src_i));
if dN_W<N
    save(str,'wavefield','dN_W','dt','N','M','dx','dy','dz','nx','ny','nz','-v7.3')
end






