%% simultion for a homogeneous (and VTI) model with general anisotropic parameters

clc; close all; clear all;
code_path=['.' filesep, 'FDwave'];      	% Path of FD code files 
addpath(genpath(code_path));      	% Add the code folder to the current command space
wf_path=pwd;               	% where you want to store your data 

T=1000;
vp=4640;
vs=2583;
rho=2490;


%% The 21 independent elastic parameters for full anisotropic models
[cc11,   cc12,   cc13,   cc14,   cc15,   cc16  ] = deal(zeros(size(vp)));
[cc22,   cc23,   cc24,   cc25,   cc26  ] = deal(zeros(size(vp)));
[cc33,   cc34,   cc35,   cc36  ] = deal(zeros(size(vp)));
[cc44,   cc45,   cc46  ] = deal(zeros(size(vp)));
[cc55,   cc56] = deal(zeros(size(vp)));
[cc66] = deal(zeros(size(vp)));



% TI parameters
eps=[0.334];%
gamma=[0.575];%
delta=[0.75];%
for i=1:length(vp)
        thomsen_parameters=[vp(i),vs(i),eps(i),gamma(i),delta(i)];
        cc = thomsen_to_c(thomsen_parameters);
        cc11(i)=cc(1)*rho(i);
        cc33(i)=cc(2)*rho(i);
        cc44(i)=cc(3)*rho(i);
        cc66(i)=cc(4)*rho(i);
        cc13(i)=cc(5)*rho(i);
        
        
        cc22(i)=cc11(i);% VTI
        cc55(i)=cc44(i);
        cc23(i)=cc13(i);
        cc12(i)=cc11(i)-2*cc66(i);
        
%         cc22(i)=cc33(i); % HTI
%         cc12(i)=cc13(i);
%         cc23(i)=cc22(i)-2*cc44(i);
%         cc55(i)=cc66(i);
        
        Cv=[cc11(i),  cc12(i),  cc13(i),  cc14(i),  cc15(i), cc16(i);...
           cc12(i),  cc22(i),  cc23(i),  cc24(i),  cc25(i), cc26(i);...
           cc13(i),  cc23(i),  cc33(i),  cc34(i),  cc35(i), cc36(i);...
           cc14(i),  cc24(i),  cc34(i),  cc44(i),  cc45(i), cc46(i);...
           cc15(i),  cc25(i),  cc35(i),  cc45(i),  cc55(i), cc56(i);...
           cc16(i),  cc26(i),  cc36(i),  cc46(i),  cc56(i), cc66(i)];
       
       CC{i}=Cv;
        
              
end


nab=20;
M=10;            % the order of FD operator M=2N

FDwave_model_n_3Dlayers_FulAni('WFP',wf_path,'WAVE_TYPE','Elastic','DX',10,'DY',10,...
	'DZ',10,'THICKNESS',T,'XZ_RATIO',1,'NAB',nab,'Vp',vp,'VS',vs,'RHO',rho,...
	'C',CC,'PlotON','y3d','verbose','y')

% source wavelet
FDwave_source_ricker('WFP',wf_path,'T',0.35,'DT',.0005,'F0',30,'T0',0.03,'PlotON','y','verbose','y');

% analyze/check the parameters
FDwave_analyse_3Delastic('WFP',wf_path,'2N',M,'verbose','y')

% derived model
FDwave_model_derived_3Delastic_2N('WFP',wf_path,'verbose','y')

% boundary conditions
FDwave_bc_3Dselect('WFP',wf_path,'BCNAME','PML','BCTYPE','topABC','NAB',nab,'PlotON','y3d','verbose','y') ;

% source and receiver geometry
load([wf_path,[filesep,'model']],'dx','dy','dz','nx','ny','nz');
FDwave_3Dgeometry_src_single('WFP',wf_path,'SX',round(nx/2),'SY',round(ny/2),'SZ',round(nz/2),'PlotON','y3d','verbose','y');  

FDwave_3Dgeometry_rec_st_line_surf('WFP',wf_path,'DEPTH',nab+1,'FIRST',nab+1,'LAST',nx-nab,'DIFF',1,'PlotON','y3d','verbose','y');  

% source implementation
MT0=sqrt(1/3)*[1 0 0;0 1 0;0 0 1];%ISO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the time stepping calculations of wavefield
FDwave_calculation_3Delastic_2N_vFulAni('src_i',1,'MT',MT0,'wfp',pwd,'dN_W',1000,'DN_SS',1,'DN_P',10,'2N',M,'plotON','y','verbose','n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-processing  %%%%%%%%%%%%%%%%%%%%%%%
strrec='SS_1.mat';          % this can be changed by user if required

%plotting the seismogram in image form
load (strrec);
save('homo_GenAni','Sx','Sy','Sz','-v7.3');




