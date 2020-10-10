% Script for generating the synthetic waveforms for isotropic and homogeneous model
% This part will setup all the FDwave program
clc; close all; clear all;
code_path=['.' filesep, 'FDwave'];      	% Path of FD code files 
addpath(genpath(code_path));        	% Add the code folder to the current command space
wf_path=pwd;               	% where you want to store your data 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocessing  %%%%%%%%%%%%%%%%%%%%%%%
% create/modify model
T=200;
vp=3000;
vs=vp/1.67;
rho=2500;

eps=[0];	% TI parameters (for isotropic media, equal to zero)
gamma=[0];	% TI parameters (for isotropic media, equal to zero)
delta=[0]; 	% TI parameters (for isotropic media, equal to zero)

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

% Option II: define the HTI model by rotating the VTI anticlockwize (Y) pi/2
%CC=[cc33' cc13' cc13' cc11' cc12' cc11' cc66' cc55' cc55'];    % HTI rotated by VTI anticlockwize (Y) pi/2

nab=20;				% Number of nodes for absorbing boundaries
M=10;				% the order of FD operator M=2N

FDwave_model_n_3Dlayers_TI('WFP',wf_path,'WAVE_TYPE','Elastic','DX',2.5,'DY',2.5,...
	'DZ',2.5,'THICKNESS',T,'XZ_RATIO',1,'NAB',nab,'Vp',vp,'VS',vs,'RHO',rho,...
	'C',CC,'PlotON','y3d','verbose','y')

% source wavelet
FDwave_source_ricker('WFP',wf_path,'T',0.24,'DT',.0003,'F0',60,'T0',0.03,'PlotON','y','verbose','y');

% analyze/check the parameters
FDwave_analyse_3Delastic('WFP',wf_path,'2N',M,'verbose','y')

% derived model
FDwave_model_derived_3Delastic_2N('WFP',wf_path,'verbose','y')

% boundary conditions
FDwave_bc_3Dselect('WFP',wf_path,'BCNAME','PML','BCTYPE','topABC','NAB',nab,'PlotON','y3d','verbose','y') ;

% source and receiver geometry
load([wf_path,[filesep,'model']],'dx','dy','dz','nx','ny','nz');
FDwave_3Dgeometry_src_single('WFP',wf_path,'SX',round(nx/2),'SY',round(ny/2),'SZ',nz-nab,'PlotON','y3d','verbose','y');  

FDwave_3Dgeometry_rec_downhole('XLOC',1+nab,'YLOC',1+nab,'FIRST',1+nab,'LAST',nz-nab,'DIFF',1,'PlotON','y3d')

for testno=1:2
% source implementation
if testno==1
    MT0=sqrt(1/2)*[0 1 0;1 0 0;0 0 0];%DC1
    outputname='homo_iso_dc1_2n';
else
    MT0=sqrt(1/2)*[0 0 0;0 0 -1;0 -1 0];%DC2
    outputname='homo_iso_dc2_2n';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the time stepping calculations of wavefield
tic;
FDwave_calculation_3Delastic_2N_vTI('src_i',1,'MT',MT0,'wfp',pwd,'dN_W',1000,'DN_SS',1,'DN_P',10,'2N',M,'plotON','y','verbose','n');

prot=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-processing  %%%%%%%%%%%%%%%%%%%%%%%

str='SS_1.mat';          % this can be changed by user if required

%plotting the seismogram in image form
FDwave_calculation_3Dimage_plot('wfp',wf_path,'SSFileName',str)

%plotting the seismogram in wiggle traces form
FDwave_calculation_3Dwiggle_plot('wfp',wf_path,'SSFileName',str,'scale',3)

load (str);
save(outputname,'Sx','Sy','Sz','dN_SS','dt','N','M','dx','dy','dz','nx','ny','nz','prot','-v7.3');

end

