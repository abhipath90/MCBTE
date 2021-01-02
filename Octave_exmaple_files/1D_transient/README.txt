This example simulates 1D quasi-ballistic heat transfer at equilibrium temperature 300 K.
We choose a 3D cube of side length 100 nm and apply periodic conditions in x and y directions.
Heat conduction is in the z-direction. Temperature at z = 0 and z = 100 nm is set to 303 K and 297 K, respectively.

INPUT DESCRIPTION

Boundary_prop.txt defines properties of all six boundaries. The periodic boundaries are specified with a translational vector, and prescribed/isothermal boundaries are specified by temperature.

In_bnd.txt is empty as there are no internal boundaries in the domain

mat_data.txt contains material properties of Si at 300 K - phonon frequency (rad/s), density of states (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), realaxation time (s) and polarization in the same order.

Measure_regions.txt defines the domain with a degree of refinement = 4, discretizing the domain into 8^4 = 4096 equal spatial cells for output of the heat flux and temperature.

Measure_times.txt defines the timestamps from t = 0 to t = 13 ns at 5 ps interval. The output data will at these time steps.

Out_bnd.txt contains the extents of the domain in x, y, and z axes.

Sim_param.txt lists the simulation parameters. We take N = 100000 computational particles, number of maximum scattering event = 10 (not important here), the volume of domain = 1e-21 m^3, and the equilibrium temperature = 300K.

Thermal_gradient.txt is empty as there are no periodic boundaries.


OUTPUT DESCRIPTION

detector_location.txt lists the location of all spatial cells in the order x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in the detector_location.txt and each column corresponding to each entry in the Measurement_times.txt

Temperature output T*.txt has the same format as Q*.txt but contains the deviation of temperature from the equilibrium (300 K) for each spatial cell at each timestamp.


*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM
This ocatve script is written to run on single processor as serial program. This is only intended for testing the code with an more easily avialable software than MATLAB. Number of simulation particles is chosen arbitarily low number for testing the error free run of the code. For representative numbers used in our calculations, please refer to the example files provided with MATLAB source code. To run this example problem use the following in command-line

octave --persist BTE_solution_3D.m

Alternatively, This programm can be run in Octve GUI.
