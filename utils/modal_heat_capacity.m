%{
    This utility is to calculate the modal heat capacity at each phonon 
    frequency in the material data. Input data contains phonon
    frequency (rad/s), group velocity (m/s), and relaxation time (s)
    listed in 3 columns. In the output mat_data.txt file, a fourth column is
    added that lists the heat capacity at the equilibrium temperature

%}
%% Input to the script
file_name = 'Input_mat.txt';
rho = 2330; % kg/m^3 Density of the material 
at_mass = 28.0855; % atomic mass g/mol
Teq = 300; % K temperature at which heat capacity is desired

%% Constants

hbar=1.054517e-34; % J s = m^2 kg s-1
boltz=1.38065e-23; % m2 kg s-2 K-1 = J k-1
avo = 6.02214e23; %mol-1 Avogadro number
new_file = 'mat_data.txt';

%% Calculating modal heat capacity
material = load(file_name);
F = material(:,1);
V = material(:,2);
tau = material(:,3);
de_dT_DFT = (hbar*F/Teq).^2/boltz.*exp(hbar*F/boltz/Teq)./(exp(hbar*F/boltz/Teq)-1).^2; %derivative of Bose-Einstein

% factor 3*3/61954 comes from matching high temperature limit to 3kB.
C_data = (3*3/61954)*avo*1000*rho*de_dT_DFT/at_mass;
K = sum(C_data.*V.*V.*tau)/3;

material = [material(:,1:3) C_data];
[row,col]=size(material);
fid = fopen(new_file,'w');
for ii=1:row
    for jj=1:col
        fprintf(fid,'%d ',material(ii,jj));
    end
    fprintf(fid,'\n');
end
