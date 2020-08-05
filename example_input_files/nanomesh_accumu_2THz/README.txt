This example simulates steady state heat trasnfer in a nanomesh to calculate its thermal condictivity  accumulation at 300K. In this particular case phonons with less that 2THz frequency are assumed to reflect specularly from the internal boundaries of the hole. BTE_solution_3D.m is modified to reflect that change.
We choose a 3D unit-cell lengths 34nm, 34nm and 22nm in x, y and z respectively. we apply periodic conditions in x and y directions and diffusive reflection on walls at z =0 and z =22nm.
Heat conduction happens along y direction under the influence of applied thermal gradient of 0.1K/34nm =2.94e6 K/m in y direction. The square holes of sides 11nm in x and y direction and through the thickness in z direction is at the center of unit-cell. The boundaries of internal holes are treated as diffusive mirror.

INPUT DESCRIPTION

Boundary_prop.txt defines properties of all the 10 boundaries: 6 outer boundaries followed by 4 internal boundaried from hole. The periodic boudaries are specified with translational vector and Diffusive boundaries are described with their degree of specularity =0.

In_bnd.txt contains definition of internal boundaries by specifying their starting and ending (x,y) coordinates followed by the normal pointing inward to the simulation domain.

mat_data.txt contains material properties of Si at room temperature. It has Frequency (rad/s), Density of state (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), realaxation time (s) and polarization the same order. Impurity scattering relaxation times are specified separatly in an additional column in the end.

Measure_regions.txt contains 8 regions created by extending the planes defined by internal hole boundaries. Each region is assigned a degree of refinement 1, and thus will be subdivided into 8 equal spatial cells. Total number of spatial cells will be 8*8 =64 at which flux and temperature data will be reported

Measure_times.txt is empty as this is a steady state simulation.

Out_bnd.txt contains the extents of the domain in all coordianted axes.

Sim_param.txt contains N=1000000, number of maximum scattering event = 20, volume of the domain = 2.277e-23 m^3 and equilibrium temperature = 300K.

Thermal_gradient.txt contains the pair of boundaries (1,3) and the thermal gradient vector along those boundaries.



OUTPUT DESCRIPTION

detector_location.txt will contain the location of all the spatial cells in the format x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in detector_location.txt and each colomun corresponding to each frequency defined in mat_data.txt. The results are reported as contribution of each phonon mode to the flux at a given spatial cell. Total flux at any spatical cell is calculated by summing up the contribution of all phonons at that location.

Temperature output T*.txt will have the same format as Q*.txt but it will contain the deviation of temperature from equilibrium (300K) for each spatial cell for each phonon frequency. Net deviation from equilibrium is calculated by summing up contribution of all phonons in that spatial cell.


THERMAL CONDUCTIVITY CALCULATION

Thermal conductivity accumulation of the nanomesh at 300K is calculated in matlab as
Qy = mean(load('Qy300.txt'),2);
kappa = abs(Qy/2.94e6);



*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

Case 1: Running on single processor
Include the 'src' directory in the PATH variable or copy all the file from 'src' to directory containing input files.
Run the file 'BTE_solution_3D.m' from GUI or commandline.

Case 2: Running on multiple processors on single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory
