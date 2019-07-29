# 5GNR-HetNet_Model
we (5G NR HetNet Project) present the following; 
1.	A simulation setup created for
  o	Hexagonal grid with wraparound
  o	HetNet of sub-6 (4GHz) + above-6 (28GHz) layers; Where small cells are UEs of macro layer (for backhauling)
  o	Use of 5GNR channel models 3GPP_38.901_UMa_NLOS and 3GPP_38.901_UMi_LOS (for sub-6 and above-6 respectively) per 3GPP TR 38.901
2.	Aligned system parameters of sub-6 and mmW layer (frame structure, freq, BW, etc) to 5GNR design
3.	Changed duplexing from TDD to FDD in sub-6 layer
4.	Implemented CSI acquisition for sub-6 layer (FDD duplexing)
5.	Implemented modifications for CSI acquisition for ’mobile’ and 'backhaul subscribers' (i.e. mmWave small cells) of sub-6 tier 
6.	Studied (KPIs/plots/calculations) on DL and UL capacity of 2-layer system with MMIMO backhauling:
  o	User throughputs for macro UEs ; CDFs of DL and UL UE throughputs, also separately for all UEs, mobile UEs only and “backhaul” UEs only
  o	User throughputs for small cell UEs; CDFs of DL and UL UE throughputs, also scaled according to calculated backhaul availability
  o	Generated CDF plots for cell throughput; macro layer, small cell layer and small cell layer – scaled to account for backhaul
  o	Implemented calculation of various statistics for user throughput; Median and average DL/UL user throughput, 5th and 95th percentile DL/UL user throughput
  o	Implemented calculation of statistics for cell throughput; Median and average DL/UL cell throughput
  o	Calculation for capacity loss due to wireless backhaul is calculated, based on median or average cell throughputs.
  o	Calculation for area throughput (DL and UL, 2-layer-system)

Running the Code: 
Here; each script has been used to reproduce a simulation-generated figure in the plots(s), paper publications, other sources wherever this repository shall be referenced in relation to 5G NR HetNet project. 
The naming convention of the figures (plots) may differ due to the varying requirements of each publication, however the purpose remains consistent.
The functions are used by the scripts to carry out certain tasks, such as initiating a simulation setup, generating channel correlation matrices, generating channel realizations, computing channel estimates, computing SEs, throughput (CDFs), etc.
Note that some of the functions listed, use [CVX](http://cvxr.com/cvx/) and [QuaDRiGa](http://quadriga-channel-model.de), which shall need to be installed separately.
## Software and Hardware Requirements
 The code was written to be used in Matlab and has been tested using Matlab 2018b on Academic License . Some of the scripts and functions might also work in Octave, but there is no guarantee of compatibility.
 Since the running HetNet project considers a setup with 7-19 cells, 100/128 antennas per BS tier, and 15/3 UEs per MacroCell BS or SmallCell, some of the simulations require a lot of RAM to store the channel correlation matrices and channel realizations. Moreover, this code has been tested successfully on an "HP ENVY 23-d120d TouchSmart All-in-One Desktop PC" with Windows 8.1 (64bit), with 8 GB DDR3 RAM and a 3.10 GHz Intel Core i5 processor, which should be viewed as a minimum requirement for using this code. Some of the simulations could take up to thirteen(13) hours of computation time to run, therefore we recommend that you first set nbrOfSetups = 1 to check how much time it takes for each realization of random UE location and shadow fading.

Acknowledgement 
This work was motivated by similar works published by Prof. Emil Bjornson and thus we have used code attained from his textbook (as basis) titled Massive MIMO Networks and which is available at https://github.com/emilbjornson/massivemimobook/tree/master/Code
The previous works presented a Case -study – Chapter 7 – Practical deployments) which covered the following; 
1.  A simulation setup for 
  a.  square grid with wraparound
  b.	sub-6 (2GHz) frequency range
  c.	3GPP_3D_UMi_NLOS channel model per 3GPP TR 36.873 (3D channel model for LTE)
  d.	TDD duplexing
  e.	MMIMO with 100/200 antennas on Tx side

2.	A Study on DL and UL system capacity for 
  a.	MR/RZF precoding, 
  b.	LS channel estimation
  c.	max product SINR and max-min fairness power allocation
