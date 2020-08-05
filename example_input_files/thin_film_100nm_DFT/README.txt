This example simulates steady state heat trasnfer in a 100nm thin film to calculate its thermal condictivity at 300K.
We choose a 3D cube of side length 100nm, apply periodic conditions in x and y directions and diffusive reflection on walls at z =0 and z =100nm. Heat conduction happens along y direction under the influence of applied thermal gradient of 1e7 K/m in y direction.

INPUT DESCRIPTION

Boundary_prop.txt defines properties of all the 6 outer boundaries. The periodic boudaries are specified with translational vector and Diffusive boundaries are described with their degree of specularity =0.

In_bnd.txt is empty since there are no internal boundaries

mat_data.txt contains material properties of Si at room temperature calculated using Density Functional Theory (DFT) simulation. It contains Frequency (rad/s), group velocity (m/s), realaxation time (s) and heat capacity (J/m^3K) in the same order.

Measure_regions.txt contains the whole domain as one region and a degree of refinement 3. This will discretize the domain into 8^3 =512 equal spatial cells for reporting final values of heat flux and temperature.

Measure_times.txt is empty as this is a steady state simulation.

Out_bnd.txt contains the extents of the domain in all coordianted axes.

Sim_param.txt contains N=1000000, number of maximum scattering event = 20, volume of the domain = 10e-22m^3 and equilibrium temperature = 300K.

Thermal_gradient.txt contains the pair of boundaries (1,3) and the thermal gradient vector along those boundaries.


OUTPUT DESCRIPTION

detector_location.txt will contain the location of all the spatial cells in the format x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in detector_location.txt and each colomun corresponding to each frequency defined in mat_data.txt. The results are reported as contribution of each phonon mode to the flux at a given spatial cell. Total flux at any spatical cell is calculated by summing up the contribution of all phonons at that location.

Temperature output T*.txt will have the same format as Q*.txt but it will contain the deviation of temperature from equilibrium (300K) for each spatial cell for each phonon frequency. Net deviation from equilibrium is calculated by summing up contribution of all phonons in that spatial cell.


THERMAL CONDUCTIVITY CALCULATION

Thermal conductivity of the nanomesh at 300K is calculated in matlab as
Qy = mean(sum(load('Qy300.txt'),2));
kappa = abs(Qy/1e7);


*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

Case 1: Running on single processor
Include the 'src' directory in the PATH variable or copy all the file from 'src' to directory containing input files.
Run the file 'BTE_solution_3D.m' from GUI or commandline.

Case 2: Running on multiple processors on single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory
