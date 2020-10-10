% This code is for generating the synthetic waveforms for 3 layers (anisotropic) case
% This part will setup all the FDwave program
clc; close all; clear all;
code_path=['.' filesep, 'FDwave'];      	% Path of FD code files 
addpath(genpath(code_path));      	% Add the code folder to the current command space
wf_path=pwd;               	% where you want to store your data 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocessing  %%%%%%%%%%%%%%%%%%%%%%
% create/ modify model
scale=1; % scale size of the small model
T=[750,1000,750]/scale;
vp=[3724,4640,5854];%m/s
vs=[1944,2583,3251];
rho=[2450,2490,2680];%kg/m3


for testno=1:3

%% ISO model  
if testno==1    
% TI parameters (for isotropic media, the three parameters are equal to zero)
eps=[0,0,0];%
gamma=[0,0,0];%
delta=[0,0,0];%
for i=1:length(vp)
        thomsen_parameters=[vp(i),vs(i),eps(i),gamma(i),delta(i)];
        cc = thomsen_to_c(thomsen_parameters);
        cc11(i)=cc(1)*rho(i);
        cc33(i)=cc(2)*rho(i);
        cc44(i)=cc(3)*rho(i);
        cc66(i)=cc(4)*rho(i);
        cc13(i)=cc(5)*rho(i);
end
% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
CC=[cc11' cc12' cc13' cc22' cc23' cc33' cc44' cc55' cc66'];     % VTI & HTI, given by Thomsen parameters
record_out='layer3_iso_dc1';
wavefield_out='wavefield_layer3_iso';

%% VTI model
else if testno==2
% TI parameters
eps=[0,0.334,0];%
gamma=[0,0.575,0];%
delta=[0,0.75,0];%
for i=1:length(vp)
        thomsen_parameters=[vp(i),vs(i),eps(i),gamma(i),delta(i)];
        cc = thomsen_to_c(thomsen_parameters);
        cc11(i)=cc(1)*rho(i);
        cc33(i)=cc(2)*rho(i);
        cc44(i)=cc(3)*rho(i);
        cc66(i)=cc(4)*rho(i);
        cc13(i)=cc(5)*rho(i);
end

% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
CC=[cc11' cc12' cc13' cc22' cc23' cc33' cc44' cc55' cc66'];     % VTI & HTI, given by Thomsen parameters
record_out='layer3_vti_dc1';
wavefield_out='wavefield_layer3_vti';

%% HTI model
else if testno==3 
% TI parameters
eps=[0,0.334,0];%
gamma=[0,0.575,0];%
delta=[0,0.75,0];%
for i=1:length(vp)
        thomsen_parameters=[vp(i),vs(i),eps(i),gamma(i),delta(i)];
        cc = thomsen_to_c(thomsen_parameters);
        cc11(i)=cc(1)*rho(i);
        cc33(i)=cc(2)*rho(i);
        cc44(i)=cc(3)*rho(i);
        cc66(i)=cc(4)*rho(i);
        cc13(i)=cc(5)*rho(i);
end

% Option I: define the TI model directly
cc22=cc11;cc55=cc44;cc23=cc13;cc12=cc11-2*cc66;    %% VTI
%cc22=cc33;cc12=cc13;cc23=cc22-2*cc44;cc55=cc66;   %% HTI
%CC=[cc11' cc12' cc13' cc22' cc23' cc33' cc44' cc55' cc66'];     % VTI & HTI, given by Thomsen parameters

% Option II: define the HTI model by rotating the VTI anticlockwize (Y) pi/2
CC=[cc33' cc13' cc13' cc11' cc12' cc11' cc66' cc55' cc55'];    % HTI rotated by VTI anticlockwize (Y) pi/2
record_out='layer3_hti_dc1';
wavefield_out='wavefield_layer3_hti';
    end
    end
end


nab=20;
M=10;            % the order of FD operator M=2N

FDwave_model_n_3Dlayers_TI('WFP',wf_path,'WAVE_TYPE','Elastic','DX',10,'DY',10,...
	'DZ',10,'THICKNESS',T,'XZ_RATIO',1.2,'NAB',nab,'Vp',vp,'VS',vs,'RHO',rho,...
	'C',CC,'PlotON','y3d','verbose','y')

% source wavelet
FDwave_source_ricker('WFP',wf_path,'T',1.2,'DT',.0005,'F0',30,'T0',0.03,'PlotON','y','verbose','y');

% analyze/check the parameters
FDwave_analyse_3Delastic('WFP',wf_path,'2N',M,'verbose','y')

% derived model
FDwave_model_derived_3Delastic_2N('WFP',wf_path,'verbose','y')

% boundary conditions
FDwave_bc_3Dselect('WFP',wf_path,'BCNAME','PML','BCTYPE','topABC','NAB',nab,'PlotON','y3d','verbose','y') ;

% source and receiver geometry
load([wf_path,[filesep,'model']],'dx','dy','dz','nx','ny','nz');
FDwave_3Dgeometry_src_single('WFP',wf_path,'SX',round(nx/2),'SY',round(ny/2),'SZ',round(nz/2),'PlotON','y3d','verbose','y');  

%FDwave_3Dgeometry_rec_st_line_surf('WFP',wf_path,'DEPTH',nab+1,'FIRST',nab+1,'LAST',nx-nab,'DIFF',5,'PlotON','y3d','verbose','y');  
FDwave_3Dgeometry_rec_downhole('XLOC',round(nx/2)+nab,'YLOC',round(ny/2)+nab,'FIRST',1+nab,'LAST',nz-nab,'DIFF',1,'PlotON','y3d')


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
