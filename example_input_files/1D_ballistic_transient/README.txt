This example simulates 1D ballistic heat transfer at equilibrium temperature 300K.
We choose a 3D cube of side length 3000 nm and apply periodic conditions in x and y directions.
Heat conduction is in z-direction. Temperature at z = 0 and z = 3000 nm is set to 303K and 297K, respectively.

INPUT DESCRIPTION

Boundary_prop.txt defines properties of all six boundaries. The periodic boundaries are specified with a translational vector, and prescribed/isothermal boundaries are specified by temperature.

In_bnd.txt is empty as there are no internal boundaries in the domain

mat_data.txt contains material properties - phonon frequency (rad/s), density of state (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), realaxation time (s) and polarization in the same order.
A very large relaxation time (1 sec) for all the frequency is used to simulate the ballistic transport. We take group velocity to be 12360 m/s.

Measure_regions.txt defines the domain with a degree of refinement = 4, discretizing the domain into 8^4 = 4096 equal spatial cells for output of the heat flux and temperature.

Measure_times.txt defines the timestamps from t = 0 to t = 13 ns at 5 ps interval. The output data will at these time steps.

Out_bnd.txt contains the extents of the domain in x, y, and z axes.

Sim_param.txt lists the simulation parameters. We take N = 8000000 computational particles, number of maximum scattering event = 10 (not important here), the volume of domain = 27e-18 m^3, and the equilibrium temperature = 300 K.

Thermal_gradient.txt is empty as there are no periodic boundaries.


OUTPUT DESCRIPTION

detector_location.txt lists the location of all spatial cells in the order x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in the detector_location.txt and each column corresponding to each entry in the Measurement_times.txt

Temperature output T*.txt has the same format as Q*.txt but contains the deviation of temperature from the equilibrium (300 K) for each spatial cell at each timestamp.



*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

Case 1: Running on a single processor
Include the 'src' directory in the PATH variable or copy all the files from 'src' to the directory containing input files.
Run the file 'BTE_solution_3D.m' from GUI or command-line.

Case 2: Running on multiple processors on a single computer (Shared memory)
Use example instructions in Single_node_multiple_proc.m provided in the directory.

Case 3: Running on multiple processors on multiple nodes (Distributed memory)
Use example instructions in Distributed_computing.m provided in the directory.
