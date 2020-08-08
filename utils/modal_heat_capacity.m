%{
    This utility is to calculate modal heat capacity for each phonon
    frequency in the material data. input data contains phonon
    frequencies (rad/s), group velocity (m/s) and realaxation time (s)
    arranged in 3 coulumns. In the output mat_data.txt fourth column is
    added that calculates heat capacity for the material at temperature Teq

%}
%% Input to the script
file_name = 'Input_data.txt';
rho = 2330; % kg/m^3 Density of the material 
at_mass = 28.0855; % atomic mass g/mol
Teq = 300; % K temperature at which heat capacity is desired

%% Constants

hbar=1.054517e-34; % J s = m^2 kg s-1
Kb=1.38065e-23; % m2 kg s-2 K-1 = J k-1
avo = 6.02214e23; %mol-1 Avogadro number
e = 1.6*10^(-19); % charge of electron
c = 3*10^8; % m/sec
h = 6.626068*10^(-34); % m^2kg/sec
T_high = 5000; % K
new_file = 'mat_data.txt';


%% Calculating correction factor from high temperature limit
material = load(file_name);
F = material(:,1);
V = material(:,2);
tau = material(:,3);

F_mev = 1e-12*4.1357*material(:,1)/(2*pi()); % convering freuency to meV units
Heat_capacity = 0;
for ii=1:length(F_mev)
    X_hw = e*10^(-3)*F_mev(ii)/Kb/T_high;
    exp_Xhw = exp(X_hw);
    dE = Kb*X_hw^2*exp_Xhw/(exp_Xhw - 1)^2;
    fact = 1; % assuming 1, will calculate from high temp limit
    Heat_capacity = Heat_capacity + 3*fact*dE;    
end

fact = 3*(3*Kb)/Heat_capacity;

%% Calculating modal heat capacity

de_dT_DFT = (hbar*F/Teq).^2/Kb.*exp(hbar*F/Kb/Teq)./(exp(hbar*F/Kb/Teq)-1).^2; %derivative of Bose-Einstein

% factor 3*3/61954 comes from matching high temperature limit to 3kB.
C_data = fact*avo*1000*rho*de_dT_DFT/at_mass;
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