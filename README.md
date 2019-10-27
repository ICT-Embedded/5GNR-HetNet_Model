## Description and usage: 
This Matlab 5GNR heterogeneous model simulation includes analysis for sub-6 layer of macro cells and above-6(mmW) layer of small cells with CSI acquisition design.
The launcher script for the simulation is HetNet_5GNR_simulation_launcher.m
Simulation parameters are described and defined in HetNet_5GNR_simulation_launcher.m and functionNetworkSetup_Quadriga.m

## Configuration of Matlab's search path:
Should include folders added by cvx_setup script (see documentation for installation of CVX)
Include paths to QuaDRiGa's "quadriga_src" and "quadriga_src\config" folders
Should include path to HybridPrecodingExample folder

## License:
This code is licensed under the GPLv3 license. Refer to folder name 'COPYING'.

## Software and Hardware Requirements
The code was written to be used in Matlab and has been tested using Matlab 2018b on Academic License . Some of the scripts and functions might also work in Octave, but there is no guarantee of compatibility.

Since the running HetNet project considers a setup with 7-19 cells, 100/128 antennas per BS tier, and 15/3 UEs per MacroCell BS or SmallCell, some of the simulations require a lot of RAM to store the channel correlation matrices and channel realizations. 

Moreover, this code has been tested successfully on an "HP ENVY 23-d120d TouchSmart All-in-One Desktop PC" with Windows 8.1 (64bit), with 8 GB DDR3 RAM and a 3.10 GHz Intel Core i5 processor, which should be viewed as a minimum requirement for using this code. 

Some of the simulations could take up to thirteen(13) hours of computation time to run, therefore we recommend that you first set nbrOfSetups = 1 to check how much time it takes for each realization of random UE location and shadow fading.

## Notes 
Script is based on the reference code from the following monograph: Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/2000000093. source: https://github.com/emilbjornson/massivemimobook

Script requires additional software packages to be used, which need to be downloaded and installed separately. These packages are developed independently and are delivered with separate licenses.

The channels are generated using QuaDRiGa from the Fraunhofer Heinrich Hertz Institute (http://www.quadriga-channel-model.de). This script has been tested with QuaDRiGa version 2.0.0-664.

Downlink channel matrix quantization is performed using CVX optimization from CVX Research, Inc. (http://cvxr.com/cvx/). This script has been tested with CVX 2.1, using the solver Mosek, version 8.0.0.60.
