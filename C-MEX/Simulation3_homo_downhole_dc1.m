% This code is for efficiency comparison between Matlab and C++ codes


clc; close all; clear all;

%% Grid setting
G=GridClass;
% [m]       [m]         [m]         [s]
G.x0=0;     G.y0=0;     G.z0=0;     G.t0 = 0.00; % initial point
G.nx=121;   G.ny=121;   G.nz=121;   G.nt = 801; % grid size
G.dx=2.5;   G.dy=2.5;   G.dz=2.5;   G.dt = 0.0003; % grid step (meter)
G.gridInfo;
G.setGrid;
Gold = oldGrid(G);

%% Creat velocity model
Vp = 3000*ones(G.nx,G.ny,G.nz);
Vs = Vp/1.67;
eps = 0*ones(G.nx,G.ny,G.nz);%.334
del = 0*ones(G.nx,G.ny,G.nz);%.575
gam = 0*ones(G.nx,G.ny,G.nz);%.75
rho = 2500*ones(G.nx,G.ny,G.nz);

Model.C33 = Vp.^2.*rho;
Model.C44 = Vs.^2.*rho;
Model.C55 = Model.C44;
Model.C11 = (1+2*eps).*Model.C33;
Model.C22 = Model.C11;
Model.C13 =-Model.C44 + sqrt((Model.C33-Model.C44).*(Model.C33.*(1+2*del) - Model.C44));
Model.C23 = Model.C13;
Model.C66 = (1+2*gam).*Model.C44;
Model.C12 = Model.C11 - 2*Model.C66;

Model.bx = (1./rho(1:end-1,:,:) + 1./rho(2:end,:,:))/2;
Model.by = (1./rho(:,1:end-1,:) + 1./rho(:,2:end,:))/2;
Model.bz = (1./rho(:,:,1:end-1) + 1./rho(:,:,2:end))/2;
Model.mu = Model.C44;
Model.lam = Model.C33 - 2*Model.C44;
Model.lamu = Model.C33;

%% Source and receiver parameters
% source in moment tensor style
Source.MT=sqrt(1/2)*[0 1 0;1 0 0;0 0 0];%DC1
% source in fault plane parameters
%{
% input angles are in degrees
strike =30;
dip = 70; 
rake = 90; 
gamma = 0; %tensile angle
sigma = 0.25; %Poisson¡¯s ratio
Source.MT=FPgenMT(strike, dip, rake, gamma, sigma)
%}
%Source.FT = get_ricker(G.tt,0,60);
Source.FT = get_ASOFI_signal(G,5,0,60,1); 

Acq.sx = 61;
Acq.sy = 61;
Acq.sz = 101;
Acq.rx = 21*ones(1,81);
Acq.ry = 21*ones(1,81);
Acq.rz = 21:101;
Acq.nrece = length(Acq.rx);
Acq.nshot = length(Acq.sx);
%% BC condition
BC.N =20;
BC.type = 'topABC';

ppml = -log(1e-6)*3*max(Vp(:))/(2*BC.N^3);
[BCx, BCy, BCz] = deal(zeros(G.nx,G.ny,G.nz));
for k=1:BC.N
BCx([k G.nx-k+1],:,:) = (BC.N-k+1)^2*ppml/G.dx;
BCy(:,[k G.ny-k+1],:) = (BC.N-k+1)^2*ppml/G.dy;
if strcmp(BC.type,'topFS')
BCz(:,:,G.nz-k+1) = (BC.N-k+1)^2*ppml/G.dz;
elseif strcmp(BC.type,'topABC')
BCz(:,:,[k G.nz-k+1]) = (BC.N-k+1)^2*ppml/G.dz;
end
end

BC.pmlx = BCx;
BC.pmly = BCy;
BC.pmlz = BCz;

%% Simulation
%mex FD_TC_cpp.cpp;
tic;
[SX,SY,SZ] = FD_TC_cpp(G,Model,Acq,BC,Source);
disp('C++ code finished. ')
toc
prot=toc;

%% Post processing

fig = figure('Position', [1 1 1200 900]);


subplot(1,3,1)
imagesc(Acq.rz,G.tt,SX);
title('C++ code, X component')
xlabel('Depth, m')
ylabel('Time, s')
colorbar; 

subplot(1,3,2)
imagesc(Acq.rz,G.tt,SY);
title('C++ code, Y component')
xlabel('Depth, m')
ylabel('Time, s')
colorbar; 

subplot(1,3,3)
imagesc(Acq.rz,G.tt,SZ);
title('C++ code, Z component')
xlabel('Depth, m')
ylabel('Time, s')
colorbar; 


colormap jet

%save('homo_iso_dc1_C2n','SX','SY','SZ','G','Model','Acq','Source','BC','prot','-v7.3');



