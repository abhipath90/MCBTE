# MCBTE
## Introduction
 MCBTE solves linearized Boltzmann transport equation for phonons in three-dimensional domain using variance-reduces Monte Carlo simulation. The code is written in MATLAB's distributed computing framework and scale almost linearly with problem size. The open source code is suitable for making small changes to the code to suit the problem that is being solved. The package contains following directories
* 'src' contains the source-code which is portable and tested on MATLAB/R2019b for serial and share-memory parallel computations. For distributed computing the code is tested on MATLAB/R2017b.
* 'example_input_files' contains few examples elucidating capabilities of the code. These are described in more details below.
* 'utils' contains useful utilities that are not part of the source codes but can be useful in preparation of input data or analysis of raw data.

## Features
* Both transient and steady-state simulations can be performed.
* Material data can be used that is derived from empirical relations or obtained from DFT calculation
* Space- and time-resolved output of heat flux and temperature deviation from equilibrium for transient simulations
* Frequency-resolved contribution to heat flux and temperature deviation from equilibrium for steady-state simulations


## Supplementary examples
To facilitate a smooth adaptability to the code, input files for a few example problems and sample instructions to run the code in serial, shared-memory parallel and distributed-memory parallel modes are inculed in 'example_input_files' directory. A brief discription of the example problems is as follows

* **1D_ballistic_transient:** Simulates ballistic heat transport in 1D rod with impulsivly applied temperature at the ends. Modified material data to simulate ballistic transport and other input files are provided for simulating 3000nm long rod. The equilibrium temperature is 300K and ends are kept at 303K and 297K respectively. The example also serves for validation against analytical expression.

* **1D_transient:** A quasi-ballistic heat transport in 1D rod with impulsively applied temperatures on both ends. This type of heat transport is characteristic of nano-structures with domain size comparable to phonon mean free path. All input data is provided for simulating a 100nm long rod with equilibirum temperature 300K. The ends are kept at 303K and 297K respectively.

* **nanomesh:** Thermal conductivity of a nanomesh strucutre is calculated using steady-state simulation. The example input files are provided to simulate one unit-cell of a nanomesh made of periodic arrangements of through-the-thickness square holes in Si film. The simulation calculates heat flux and termperature deviation for equilibrium in steady state under expernally applied thermal gradient. A small MATLAB code snippet is also included in README.txt file to calculate thermal conductivity from the output data of heat flux. The equilibrium temperature is chosen as 300K but it can easily be modified in Sim_param.txt

* **nanomesh_accumu_2THz:** This example elucidats how easily this open-source code can be edited make it suit unique requirements of a particular analysis. In this example same nanomesh problem described above is simulated but with an additional assumption that long wavelength phonons can follow the modified dispersion due to nano-structuring and reflect specularly from the nanomesh boundries. Only one line change in BTE_solution_3D.m is made to account for this change. The input files are essentially same as '**nanomesh**' example. The MATLAB code snippet provide for this problem calculates thermal conductivity accumulation with phonon frequency.

* **thin_film_100nm_ref** Thermal conductivity of 100nm thin film of Si at room temperature(300K) is calculated. All input files are provided and a MATLAB code snippet is given in README.txt that calculates thermal conductivity from heat flux output.

* **thin_film_100nm_DFT** This example elucidates the integration of the material properties calculated using DFT with the code. Thermal conductivity of 100nm Si thin film is calculated at room temperature(300K). The values of frequency, group velocity and relaxation times are calculated using DFT simulations and heat capacity is calculated using a utility 'modal_heat_capacity.mm' supplied in 'utils' directory. A MATLAB code snippet is provided that calculates thermal conductivity from the heat flux data.

## Utilities
'utils' directory contains 'modal_heat_capacity.m' file that is used to calculate heat capacity of each phonon frequency bin using Bose-Einstein statistics and high temperature limit of heat capacity =3kB, where 3kB is boltzmann constant.