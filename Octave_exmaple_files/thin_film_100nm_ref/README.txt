This example simulates the steady-state heat transfer in a 100 nm silicon thin-film and calculates its thermal conductivity at 300 K.
We choose a 3D cube of length 100 nm and apply periodic conditions in the x and y directions. Walls at z = 0 and z = 100 nm are taken to be diffusively reflective. The thermal gradient of 1e7 K/m is applied in the y-direction. Heat conduction along the y-direction is probed. 

INPUT DESCRIPTION

Boundary_prop.txt defines the properties of all six outer boundaries. The periodic boundaries are specified with a translational vector and diffusive boundaries are specified by taking the degree of specularity = 0.

In_bnd.txt is empty since there are no internal boundaries

mat_data.txt contains material properties of Si at 300 K -- Phonon frequency (rad/s), the density of states (s/rad m^3), group velocity (m/s), frequency bin size (rad/s), relaxation time (s), and polarization in the same order.

Measure_regions.txt defines the domain with a degree of refinement = 3, discretizing the domain into 8^3 =512 equal spatial cells for output of the heat flux and temperature.

Measure_times.txt is empty since this is a steady-state simulation.

Out_bnd.txt contains the extents of the domain in x, y, and z axes.

Sim_param.txt lists the simulation parameters. We take N = 50000 computational particles, number of maximum scattering event = 20, the volume of domain = 10e-22m^3 and the equilibrium temperature = 300 K.

Thermal_gradient.txt specifies the pair of boundaries (1,3) and the thermal gradient vector parallel to the vector connecting these boundaries.


OUTPUT DESCRIPTION

detector_location.txt lists the location of all spatial cells in the order x_min, x_max, y_min, y_max, z_min, z_max.

Flux output Q*.txt will have each row corresponding to each spatial cell in the detector_location.txt and each column corresponding to each frequency listed in the mat_data.txt. The output is the contribution of each phonon mode to the flux at a given spatial cell. Total flux at any spatial cell is calculated by summing up the contribution of all phonons at that location.

Temperature output T*.txt has the same format as Q*.txt but contains the deviation of temperature from the equilibrium (300K) for each spatial cell at each phonon frequency. Net deviation from equilibrium temperature is calculated by summing up the contribution of all phonons in the spatial cell.


THERMAL CONDUCTIVITY CALCULATION

Thermal conductivity of the silicon thin-film at 300K is calculated using additional two lines of code as follows:
Qy = mean(sum(load('Qy300.txt'),2));
kappa = abs(Qy/1e7);


*******************************************************************************************************
INSTRUCTION TO RUN THE PROGRAM

This ocatve script is written to run on single processor as serial program. This is only intended for testing the code with an more easily avialable software than MATLAB. Number of simulation particles is chosen arbitarily low number for testing the error free run of the code. For representative numbers used in our calculations, please refer to the example files provided with MATLAB source code. To run this example problem use the following in command-line

octave --persist BTE_solution_3D.m

Alternatively, This programm can be run in Octve GUI.
