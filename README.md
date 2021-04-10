# FDwave3D


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FDWAVE3D IS DEVELOPED BASED ON A 2D VERSION PACKAGE OF VECTORIZED FD OPERATOR 
THE PACKAGE AND THE TESTS ARE CONDUCTED WITH MATLAB 2016b UNDER BOTH LINUX and WINDOWS SYSTEMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CHARACTERSTICS OF THIS 3D FINITE DIFFERENCE MODELLING CODE PACKAGE ARE

	1) In- time domain 
	2) For- anisotropic elastic media and moment tensor sources
	3) Over- staggered grid 
	4) Uses- vectorized finite-difference operator
	5) With- Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CONTENT IN PACKAGE

	manual.pdf:
	Instruction and documentation of the code package.

	FDwave: 
	This directory contains all the programs & functions related to seismic modeling.
	
	Matlab_vs_CPP:
	The codes of the test for the efficiency comparison between Matlab and C++ codes.
 
	The scripts used to reproduce the records and figures in the related paper are also contained.

	SEG/EAGE overthrust model 
	https://wiki.seg.org/wiki/SEG/EAGE_Salt_and_Overthrust_Models
	or download via the link below (about 140 M):
	https://drive.google.com/file/d/1vTemFS0poXAUMhHfea-nVGRSvuLPNa5K/view?usp=sharing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


HOW TO RUN THE PROGRAM
	
	More detailed instructions are presented in the manual file
	The figures can be found in the related manuscript

	e.g.,
	Figure 2: 
	Simulation3_homo_downhole.m: 		for the 3D homogeneous and isotropic model
	compare_homo.m

	Figure 3:
	Simulation3_layer3_TI.m: 			for the layered anisotropic model
	showwavefield_layer.m:			show original wavefields of the layered model (Figures 3a-c)
	showrecord_layer.m:			show original seismograms of the layered model (Figures 3d-f)

	Figure 4: show_overthrust.m

	Figure 5:
	Simulation3_overthrust_TI.m: 			for the anisotropic overthrust model 
	showwavefield_overthrust.m:			show original wavefields of the overthrust model (Figures 5a-c)
	showrecord_overthrust.m:			show original seismograms of the overthrust model (Figures 5d-f)

	Figure 6: compare_cost.m


	All the results can be directly downloaded via the link below (about 150 MB) for a quick validation.
	https://drive.google.com/file/d/1mCKpKfma-oWOW9pfu3bvknVzmX1hsEFk/view?usp=sharing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GETTING HELP

	Details about a function in the FDwave folder can be obtained by typing "help fun_name" in the command window.
	More details of the package can be found in above references and the manual file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


REFERENCE

  Malkoti, A., Vedanti, N., Tiwari, R.K.: An algorithm for fast elastic wave simulation using a vectorized finite difference operator. Computers & Geosciences. 116, 23-31 (2018). https:
//doi.org/10.1016/j.cageo.2018.04.002.

  Li, L., Tan, J., Zhang, D., Malkoti, A., Abakumov, I., Xie, Y.: FDwave3D: A MATLAB solver for the 3D anisotropic wave equation using the finite-difference method. Computational Geosciences. 2021. Accepted Manuscript.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CONTACT

  Lei Li: leileely@126.com
  Ajay Malkoti: ajmalkoti@gmail.com


