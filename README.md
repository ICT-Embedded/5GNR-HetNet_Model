## Running the Code: 
Here; each script has been used to reproduce a simulation-generated figure in the plots(s), paper publications, other sources wherever this repository shall be referenced in relation to 5G NR HetNet project. 
The naming convention of the figures (plots) may differ due to the varying requirements of each publication, however the purpose remains consistent.
The functions are used by the scripts to carry out certain tasks, such as initiating a simulation setup, generating channel correlation matrices, generating channel realizations, computing channel estimates, computing SEs, throughput (CDFs), etc.
Note that some of the functions listed, use [CVX](http://cvxr.com/cvx/) and [QuaDRiGa](http://quadriga-channel-model.de), which shall need to be installed separately.
## Software and Hardware Requirements
 The code was written to be used in Matlab and has been tested using Matlab 2018b on Academic License . Some of the scripts and functions might also work in Octave, but there is no guarantee of compatibility.
 Since the running HetNet project considers a setup with 7-19 cells, 100/128 antennas per BS tier, and 15/3 UEs per MacroCell BS or SmallCell, some of the simulations require a lot of RAM to store the channel correlation matrices and channel realizations. Moreover, this code has been tested successfully on an "HP ENVY 23-d120d TouchSmart All-in-One Desktop PC" with Windows 8.1 (64bit), with 8 GB DDR3 RAM and a 3.10 GHz Intel Core i5 processor, which should be viewed as a minimum requirement for using this code. Some of the simulations could take up to thirteen(13) hours of computation time to run, therefore we recommend that you first set nbrOfSetups = 1 to check how much time it takes for each realization of random UE location and shadow fading.

