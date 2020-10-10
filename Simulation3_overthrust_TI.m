% This code is for generating the synthetic waveforms for overthrust model
% This part will setup all the FDwave program
clc; close all; clear all;
code_path=['.' filesep, 'FDwave'];       % Path of FD code files 
addpath(genpath(code_path));        % Add the code folder to the current command space
wf_path=pwd;                % where you want to store your data 


%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocessing  %%%%%%%%%%%%%%%%%%%%%
% create/ modify model
load overthrust_3d_vp;
vp=data(1:2:end,1:2:end,1:2:end);
clear data;
vp=permute(vp,[1,3,2]);
vs=vp./1.7;
rho=1.741*(vp/1000).^(0.25)*1000;

for testno=1:3

%% ISO model  
if testno==1    
% TI parameters (for isotropic media, the three parameters are equal to zero)
eps=zeros(size(vp));%[0,0,0];%.334
gamma=zeros(size(vp));%[0,0,0];%.575
delta=zeros(size(vp));%[0,0,0];%.75


for i=1:size(vp,1)
    for j=1:size(vp,2)
        for k=1:size(vp,3)
            thomsen_parameters=[vp(i,j,k),vs(i,j,k),eps(i,j,k),gamma(i,j,k),delta(i,j,k)];
            cc = thomsen_to_c(thomsen_parameters);
            cc11(i,j,k)=cc(1)*rho(i,j,k);
            cc33(i,j,k)=cc(2)*rho(i,j,k);
            cc44(i,j,k)=cc(3)*rho(i,j,k);
            cc66(i,j,k)=cc(4)*rho(i,j,k);
            cc13(i,j,k)=cc(5)*rho(i,j,k);
        end
    end
end
% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
CC={cc11 cc12 cc13 cc22 cc23 cc33 cc44 cc55 cc66};    % VTI & HTI
record_out='overthrust_iso_dc1';
wavefield_out='wavefield_overthrust_iso';

%% VTI model
else if testno==2
% TI parameters (for isotropic media, the three parameters are equal to zero)
eps=zeros(size(vp));%[0,0,0];%.334
gamma=zeros(size(vp));%[0,0,0];%.575
delta=zeros(size(vp));%[0,0,0];%.75

% anisotropy surrounding the source area (for isotropic media, the three parameters are equal to zero)
eps(26:75,151:250,151:250)=.334;
gamma(26:75,151:250,151:250)=.575;
delta(26:75,151:250,151:250)=.75;

for i=1:size(vp,1)
    for j=1:size(vp,2)
        for k=1:size(vp,3)
            thomsen_parameters=[vp(i,j,k),vs(i,j,k),eps(i,j,k),gamma(i,j,k),delta(i,j,k)];
            cc = thomsen_to_c(thomsen_parameters);
            cc11(i,j,k)=cc(1)*rho(i,j,k);
            cc33(i,j,k)=cc(2)*rho(i,j,k);
            cc44(i,j,k)=cc(3)*rho(i,j,k);
            cc66(i,j,k)=cc(4)*rho(i,j,k);
            cc13(i,j,k)=cc(5)*rho(i,j,k);
        end
    end
end
% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
CC={cc11 cc12 cc13 cc22 cc23 cc33 cc44 cc55 cc66};    % VTI & HTI
record_out='overthrust_vti_dc1';
wavefield_out='wavefield_overthrust_vti';

%% HTI model
else if testno==3 
% TI parameters (for isotropic media, the three parameters are equal to zero)
eps=zeros(size(vp));%[0,0,0];%.334
gamma=zeros(size(vp));%[0,0,0];%.575
delta=zeros(size(vp));%[0,0,0];%.75

% anisotropy surrounding the source area (for isotropic media, the three parameters are equal to zero)
eps(26:75,151:250,151:250)=.334;
gamma(26:75,151:250,151:250)=.575;
delta(26:75,151:250,151:250)=.75;

for i=1:size(vp,1)
    for j=1:size(vp,2)
        for k=1:size(vp,3)
            thomsen_parameters=[vp(i,j,k),vs(i,j,k),eps(i,j,k),gamma(i,j,k),delta(i,j,k)];
            cc = thomsen_to_c(thomsen_parameters);
            cc11(i,j,k)=cc(1)*rho(i,j,k);
            cc33(i,j,k)=cc(2)*rho(i,j,k);
            cc44(i,j,k)=cc(3)*rho(i,j,k);
            cc66(i,j,k)=cc(4)*rho(i,j,k);
            cc13(i,j,k)=cc(5)*rho(i,j,k);
        end
    end
end
% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
%CC={cc11 cc12 cc13 cc22 cc23 cc33 cc44 cc55 cc66};    % VTI & HTI
% Option II: define the HTI model by rotating the VTI anticlockwize (Y) pi/2
CC={cc33 cc13 cc13 cc11 cc12 cc11 cc66 cc55 cc55};   % HTI rotated by VTI 
record_out='overthrust_hti_dc1';
wavefield_out='wavefield_overthrust_hti';
    end
    end
end

nab=20;             % Number of nodes for absorbing boundaries
M=10;               % the order of FD operator M=2N

FDwave_model_3Dload('WFP',wf_path,'WAVE_TYPE','Elastic','DX',10,'DY',10,'DZ',10,... 
         'NAB',nab,'Vp',vp,'VS',vs,'RHO',rho,'C',CC,'PlotON','y3d','verbose','y')

% source wavelet
FDwave_source_ricker('WFP',wf_path,'T',1,'DT',.0005,'F0',30,'T0',0.03,'PlotON','n','verbose','y');

% analyze/check the parameters
FDwave_analyse_3Delastic('WFP',wf_path,'2N',M,'verbose','y')

% derived model
FDwave_model_derived_3Delastic_2N('WFP',wf_path,'verbose','y')

% boundary conditions
FDwave_bc_3Dselect('WFP',wf_path,'BCNAME','PML','BCTYPE','topABC','NAB',nab,'PlotON','y3d','verbose','y') ;

% source and receiver geometry
load([wf_path,[filesep,'model']],'dx','dy','dz','nx','ny','nz');
FDwave_3Dgeometry_src_single('WFP',wf_path,'SX',round(nx/2),'SY',round(ny/2),'SZ',round(nz/2),'PlotON','y3d','verbose','y');  

FDwave_3Dgeometry_rec_st_line_surf('WFP',wf_path,'DEPTH',nab+1,'FIRST',nab+1,'LAST',nx-nab,'DIFF',2,'PlotON','y3d','verbose','y');  
%FDwave_3Dgeometry_rec_downhole('XLOC',round(nx/2)+20,'YLOC',round(ny/2)+20,'FIRST',1+nab,'LAST',nz-nab,'DIFF',1,'PlotON','y3d')

% source implementation
MT0=sqrt(1/2)*[0 1 0;1 0 0;0 0 0];%DC1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the time stepping calculations of wavefield
tic;

FDwave_calculation_3Delastic_2N_vTI('src_i',1,'MT',MT0,'wfp',pwd,'dN_W',1000,'DN_SS',1,'DN_P',10,'2N',M,'plotON','y','verbose','n');

prot=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-processing  %%%%%%%%%%%%%%%%%%%%%%%
strrec='SS_1.mat';          % this can be changed by user if required
strwav='wavefield_1.mat'; 
% plotting the seismogram in image form
%FDwave_calculation_3Dimage_plot('wfp',wf_path,'SSFileName',strrec)

% plotting the seismogram in wiggle traces form
%FDwave_calculation_3Dwiggle_plot('wfp',wf_path,'SSFileName',strrec,'scale',3)

load (strrec);
load (strwav);
save(record_out,'Sx','Sy','Sz','dN_SS','dt','N','M','dx','dy','dz','nx','ny','nz','prot','-v7.3');
save(wavefield_out,'wavefield');

end

