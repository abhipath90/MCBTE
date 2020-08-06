# MCBTE
## Introduction
 MCBTE solves the linearized Boltzmann transport equation for phonons in three-dimension using a variance-reduced Monte Carlo solution approach. The code is written in MATLAB's distributed computing framework and scales linearly with size. The source code can be modified to adapt to new computational studies. The package contains following directories -
* 'src' contains the source code which is portable and tested on MATLAB/R2019b for serial and share-memory parallel computations. For distributed computing, the code is tested on MATLAB/R2017b.
* 'example_input_files' contains few illustrating examples highlighting the capabilities of the code. See below for more details.
* 'utils' contains useful utilities that are not part of the source code but are useful for the preparation of input data or analysis of raw data.

## Features
* Simulation of 3D geometry. Modeling of 1D and 2D geometry using appropriate boundary conditions, such as periodic boundaries.
* Simulation of nanostructures using periodic boundary conditions under a constant thermal gradient.
* Integrated with both first-principles density functional theory calculations (such as obtained from almaBTE, https://almabte.bitbucket.io/) and empirical relations for the input of phonon frequency, group velocity, and mean free path.
* Space- and time-resolved output of heat flux and temperature for transient simulations.
* Phonon frequency-resolved output of heat flux and temperature for steady-state simulations.
* Parallelized using MATLAB's distributed computing framework providing near-linear scaling with size.
* Benchmarked against available analytical/computational/theoretical results in the literature.


## Illustrative examples
We provide input files in 'example_input_files' directory for selected example problems and instructions to run the code in serial, shared-memory parallel, and distributed-memory parallel modes. A brief description of the example problems is as follows.

* **1D_ballistic_transient:** Simulates ballistic heat transport in a 1D rod of length 3000 nm for a constant temperature gradient. Here phonon mean free path is much larger than the domain size.

* **1D_transient:** Simulates quasi-ballistic heat transport in a 1D rod of length 100 nm for a constant temperature gradient. Here phonon mean free path is comparable to the domain size.

* **nanomesh:** Calculates the thermal conductivity of a periodic nanomesh structure under steady-state conditions. In input files, the equilibrium temperature is 300K but can be modified in Sim_param.txt to obtain thermal conductivity as a function of temperature.

* **nanomesh_accumu_2THz:** This example illustrates how easily the source code can be edited to suit the unique requirements of a particular analysis. Here, we solve for the same nanomesh structure described above with an additional constraint that the long-wavelength phonons (< 2THz) scatter differently from the nanomesh boundaries in comparison to > 2THz phonons. We make a one-line change in the BTE_solution_3D.m to account for this change.

* **thin_film_100nm_ref:** Calculates the thermal conductivity of a 100nm thin film of Si at 300K.

* **thin_film_100nm_DFT:** This example illustrates the integration of code with first-principles density functional theory (DFT) simulations. The thermal conductivity of 100 nm Si thin film is calculated at 300K. 

## Utilities
'utils' directory contains 'modal_heat_capacity.m' file, which calculates the heat capacity for a given phonon frequency, density, and atomic mass. This is useful when using phonon frequency, group velocity, and relaxation times calculated from first-principles DFT calculations.
