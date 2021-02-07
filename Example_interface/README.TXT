This example simulates the steady-state heat transfer in a 100 nm silicon bilayer at equilibrium temperature 300 K. Both the layers are 50 nm thick and interface between them is modelled as a high-pass filter.
The high pass filter lets phonons pass through the interface in the high frequency part of spectrum. In this case we have chosen frequency cutoff at 0.6 of maximum frequency for filter.
We choose a 3D cube of length 100 nm and apply periodic conditions in the x and y directions. Walls at z = 0 is istothermal with T = 303 K and wall at z = 100 nm is isothermal with T = 297 K.
The interface is located at z = 50 nm. 

INPUT DESCRIPTION

Boundary_prop.txt defines the properties of all six outer boundaries. The periodic boundaries are specified with a translational vector and isothermal boundaries are specified by defining their temperature.

In_bnd.txt is empty since there are no internal boundaries

mat_data.txt contains material properties of Si at 300 K -- Phonon frequency (rad/s), the density of states (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), relaxation time (s), and polarization in the same order.

Measure_regions.txt defines the domain with a degree of refinement = 3, discretizing the domain into 8^3 =512 equal spatial cells for output of the heat flux and temperature.

Measure_times.txt is empty since this is a steady-state simulation.

Out_bnd.txt contains the extents of the domain in x, y, and z axes.

Sim_param.txt lists the simulation parameters. We take N = 1000000 computational particles, number of maximum scattering event = 30, the volume of domain = 10e-22m^3 and the equilibrium temperature = 300 K.

Thermal_gradient.txt is empty and the source of thermal energy are two isothermal boundaries.

OUTPUT DESCRIPTION

detector_location.txt lists the location of all spatial cells in the order x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in the detector_location.txt and each column corresponding to each frequency listed in the mat_data.txt. The output is the contribution of each phonon mode to the flux at a given spatial cell. Total flux at any spatial cell is calculated by summing up the contribution of all phonons at that location.

Temperature output T*.txt has the same format as Q*.txt but contains the deviation of temperature from the equilibrium (300K) for each spatial cell at each phonon frequency. Net deviation from equilibrium temperature is calculated by summing up the contribution of all phonons in the spatial cell.

Output directory contains all these output files along with a plot of temperature profile along z axis in the domain. The temperature profile shows the temperature jump expected at the interface location.

*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

This directory contains all the required files neede to run this example. Most of these files are identical to corresponding files located in 'src' directory. A couple of files are changed to incorporate interface scattering.

Case 1: Running on a single processor
Run the file 'BTE_solution_3D.m' from GUI or command-line.

Case 2: Running on multiple processors on a single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory.

