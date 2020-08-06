This example simulates the steady-state heat transfer in a periodic nanomesh structure to calculate its thermal conductivity at 300K. In this particular case, phonons of frequency less than 2 THz are specularly reflected from the internal boundaries of the pore. BTE_solution_3D.m is modified to reflect this change (only a single line is modified). We choose a 3D unit-cell of length 34 nm, 34 nm, and 22 nm along x, y, and z directions, respectively. We apply periodic boundary conditions in x and y directions and diffusely reflective boundary conditions on boundaries at z = 0 and z = 22 nm.
The thermal gradient of 0.1/34 = 2.94e6 K/m is applied along the y-direction. Heat conduction is along the y-direction. The square pores of 11 nm length in the x-y plane are drilled through the thickness, i.e.,  z-direction. The boundaries of pore are treated as diffusively reflective.

INPUT DESCRIPTION

Boundary_prop.txt defines the properties of all 10 boundaries - 6 outer boundaries followed by 4 internal boundaries from the pore. The periodic boundaries are specified with a translational vector and diffusive boundaries are specified by taking the degree of specularity = 0.

In_bnd.txt defines the internal boundaries (boundaries of a pore) by specifying their start and end x and y coordinates followed by the normal pointing inward in the simulation domain.

mat_data.txt contains material properties of Si at 300 K -- Phonon frequency (rad/s), the density of state (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), relaxation time (s), and polarization in the same order. Impurity scattering relaxation times are specified separately in an additional column in the end.

Measure_regions.txt contains 8 regions created by extending the planes defined by the square pore boundaries. Each region is assigned a degree of refinement = 1; hence is subdivided into 8 equal spatial cells. The total number of spatial cells is 8*8 =64 at which the flux and temperature data is calculated.

Measure_times.txt is empty since this is a steady-state simulation.

Out_bnd.txt contains the extents of the domain in x, y, and z axes.

Sim_param.txt lists the simulation parameters. We take N = 1000000 computational particles, number of maximum scattering event = 20, the volume of domain = 2.277e-23 m^3 and the equilibrium temperature = 300 K.

Thermal_gradient.txt specifies the pair of boundaries (1,3) and the thermal gradient vector parallel to the vector connecting these boundaries.

OUTPUT DESCRIPTION

detector_location.txt lists the location of all spatial cells in the order x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in the detector_location.txt and each column corresponding to each frequency listed in the mat_data.txt. The output is the contribution of each phonon mode to the flux at a given spatial cell. Total flux at any spatial cell is calculated by summing up the contribution of all phonons at that location.

Temperature output T*.txt has the same format as Q*.txt but contains the deviation of temperature from the equilibrium (300K) for each spatial cell at each phonon frequency. Net deviation from equilibrium temperature is calculated by summing up the contribution of all phonons in the spatial cell.

THERMAL CONDUCTIVITY CALCULATION

Thermal conductivity of the periodic nanomesh structure at 300K is calculated using additional two lines of code as follows:
Qy = mean(sum(load('Qy300.txt'),2));
kappa = abs(Qy/2.94e6);

*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

Case 1: Running on a single processor
Include the 'src' directory in the PATH variable or copy all the files from 'src' to the directory containing input files.
Run the file 'BTE_solution_3D.m' from GUI or command-line.

Case 2: Running on multiple processors on a single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory.
