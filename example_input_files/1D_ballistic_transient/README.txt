This example simulates 1D ballistic heat transfer at equilibrium temperature 300K.
We choose a 3D cube of side length 3000nm and apply periodic conditions in x and y directions.
Heat conduction happens in z direction by impulsively applying wall temperatures at z=0 and z =3000nm to 303K and 297K respectively.

INPUT DESCRIPTION

Boundary_prop.txt defines properties of all the six boundaries. The periodic boudaries are specified with translational vector and prescribed boundaries are described with their temperature.

In_bnd.txt is empty as there are no internal boundaries in the domain

mat_data.txt contains material properties Frequency (rad/s), Density of state (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), realaxation time (s) and polarization the same order.
A very large relaxation time (1s) for all the frequency is used to simulate ballistic transport. Debye approximation of dispersion relation is assumed with group velocity of 12360m/s.

Measure_regions.txt contains the whole domain as one region and a degree of refinement 4. This will discretize the domain into 8^4 = 4096 equal spatial cells for reporting final values of heat flux and temperature.

Measure_times.txt contains time stamps from t=0 to t=13ns at 5ps interval. The output data will be generated at these time steps.

Out_bnd.txt contains the extents of the domain in all coordianted axes.

Sim_param.txt contains N=8000000, number of maximum scattering event = 10 (not important here), volume of the domain = 27e-18 m^3 and equilibrium temperature = 300K.

Thermal_gradient.txt is empty as there is no externally applied thermal gradient in the problem.



OUTPUT DESCRIPTION

detector_location.txt will contain the location of all the spatial cells in the format x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in detector_location.txt and each colomun corresponding to each entry in Measurement_times.txt

Temperature output T*.txt will have the same format as Q*.txt but it will contain the deviation of temperature from equilibrium (300K) for each spatial cell at each time stamp.



*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

Case 1: Running on single processor
Include the 'src' directory in the PATH variable or copy all the file from 'src' to directory containing input files.
Run the file 'BTE_solution_3D.m' from GUI or commandline.

Case 2: Running on multiple processors on single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory
