# MCBTE
## Introduction
 MCBTE solves the linearized Boltzmann transport equation for phonons in three-dimension using a variance-reduced Monte Carlo solution approach. For a detailed description of the formulation and its implementation, please refere to "Pathak *et al.*, "**MCBTE: A variance-reduced Monte Carlo solution of the linearized Boltzmann transport equation for phonons**", http://arxiv.org/abs/2104.03573, https://doi.org/10.1016/j.cpc.2021.108003". The code is written in MATLAB's distributed computing framework and scales linearly with size. The source code can be modified to adapt to new computational studies. The package contains following directories -
* 'src' contains the source code which is portable and tested on MATLAB/R2019b for serial and share-memory parallel computations. For distributed computing, the code is tested on MATLAB/R2017b.
* 'example_input_files' contains few illustrating examples highlighting the capabilities of the code. See below for more details.
* 'utils' contains useful utilities that are not part of the source code but are useful for the preparation of input data or analysis of raw data.
## How to cite
https://doi.org/10.1016/j.cpc.2021.108003

Pathak, A., Pawnday, A., Roy, A. P., Aref, A. J., Dargush, G. F., & Bansal, D. (2021). MCBTE: A variance-reduced Monte Carlo solution of the linearized Boltzmann transport equation for phonons. Computer Physics Communications, 265, 108003.
````
@article{pathak2021mcbte,
  title={MCBTE: A variance-reduced Monte Carlo solution of the linearized Boltzmann transport equation for phonons},
  author={Pathak, Abhishek and Pawnday, Avinash and Roy, Aditya Prasad and Aref, Amjad J and Dargush, Gary F and Bansal, Dipanshu},
  journal={Computer Physics Communications},
  volume={265},
  pages={108003},
  year={2021},
  publisher={Elsevier}
  }
  ````
## Features
* Simulation of 3D geometry. Modeling of 1D and 2D geometry using appropriate boundary conditions, such as periodic boundaries.
* Simulation of nanostructures using periodic boundary conditions under a constant thermal gradient.
* Integrated with both first-principles density functional theory calculations (such as obtained from almaBTE, https://almabte.bitbucket.io/) and empirical relations for the input of phonon frequency, group velocity, and mean free path.
* Space- and time-resolved output of heat flux and temperature for transient simulations.
* Phonon frequency-resolved output of heat flux and temperature for steady-state simulations.
* Parallelized using MATLAB's distributed computing framework providing near-linear scaling with size.
* Benchmarked against available analytical/computational/theoretical results in the literature.

## Instructions to run the program
The program can be executed either on a single node with multiple processors sharing the same memory using 'Single_node_multiple_proc.m' or on multiple nodes using 'Distributed_computing.m'. Both of the files are available in the 'example_input_files' directory. The MATLAB package requires access to 'Parallel Computing Toolbox' for execution and can be run from GUI or command-line. An open-source alternative, an Octave implementation, is also provided in 'Octave_implementation' directory. The Octave program can either be run from GUI or on command-line using *octave --persist BTE_solution_3D.m*.

## Illustrative examples
We provide input files in 'example_input_files' directory for selected example problems and instructions to run the MATLAB code in serial, shared-memory parallel, and distributed-memory parallel modes. The same example problems are provided for Octave implementation separately inside 'Octave_implementation/Octave_example_files'. A brief description of the example problems is as follows. A more detailed README.txt is provided in each folder.

* **1D_ballistic_transient:** Simulates ballistic heat transport in a 1D rod of length 3000 nm for a constant temperature gradient. Here phonon mean free path is much larger than the domain size.

* **1D_transient:** Simulates quasi-ballistic heat transport in a 1D rod of length 100 nm for a constant temperature gradient. Here phonon mean free path is comparable to the domain size.

* **nanomesh:** Calculates the thermal conductivity of a periodic nanomesh structure under steady-state conditions. In input files, the equilibrium temperature is 300 K but can be modified in Sim_param.txt to obtain thermal conductivity as a function of temperature.

* **nanomesh_accumu_2THz:** This example illustrates how easily the source code can be edited to suit the unique requirements of a particular analysis. Here, we solve for the same nanomesh structure described above with an additional constraint that the long-wavelength phonons (< 2 THz) scatter differently from the nanomesh boundaries in comparison to > 2 THz phonons. We make a one-line change in the BTE_solution_3D.m to account for this change.

* **thin_film_100nm_ref:** Calculates the thermal conductivity of a 100 nm thin film of Si at 300 K.

* **thin_film_100am_DFT:** This example illustrates the integration of code with first-principles density functional theory (DFT) simulations. The thermal conductivity of 100 nm Si thin film is calculated at 300 K. 

## Utilities
'utils' directory contains 'modal_heat_capacity.m' file, which calculates the heat capacity for a given phonon frequency, density, and atomic mass. This is useful when using phonon frequency, group velocity, and relaxation times calculated from first-principles DFT calculations.

## Developers
Abhishek Pathak apathak@buffalo.edu, abhishekpathak9066@gmail.com \
Avinash Pawnday \
Dipanshu Bansal dipanshu@iitb.ac.in
