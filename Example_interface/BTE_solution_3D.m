function [time]=BTE_solution_3D(~)
    %{
      This method solves the Boltzmann Transport Equation (BTE) for phonons in three-dimensional (3D) geometries.
      For running the simulation using distributed computing a dummy input (neglected) and a dummy output (simulation time, again not used) are assigned.
      Material properties, geometric description, and simulation parameters are supplied using designated text files.
      By supplying Measure_times.txt, a transient simulation is performed otherwise a steady-state simulation is performed.
      Output files contain three components of heat flux and deviation of temperature from equilibrium/linearization value.
    %}
%% Reading all the input files to pass it to all the workers

% Defining file names
File_mat = 'mat_data.txt'; % contains material data
File_out_bnd = 'Out_bnd.txt'; % The outer extents of the domain 
File_in_bnd = 'In_bnd.txt'; % Definition of any internal boundaries
File_bnd_prop = 'Boundary_prop.txt'; % Properties of all the boundaries
File_measure_reg = 'Measure_regions.txt'; % Regions where output is calculated
File_thermal_grad = 'Thermal_gradient.txt'; % Thermal gradient source
File_measure_time = 'Measure_times.txt'; % Measurement times for transient simulation
File_param = 'Sim_param.txt'; % Other simulation paramters

interface = 50e-9; % z coordinate of interface

% Reading all the files at the head processor and distributing it to all the workers
if labindex==1
    %% Reading file if they exist
    Mat_data = ReadFile(File_mat);
    Out_bnd = ReadFile(File_out_bnd);
    In_bnd = ReadFile(File_in_bnd);
    Bnd_prop = ReadFile(File_bnd_prop);
    Measure_reg = ReadFile(File_measure_reg);
    Thermal_grad = ReadFile(File_thermal_grad);
    Measure_time = ReadFile(File_measure_time);
    Param = ReadFile(File_param);
    
    data_coll = {Mat_data;Out_bnd;In_bnd;Bnd_prop;Measure_reg;Thermal_grad;Measure_time;Param};
    data_sent = labBroadcast(1,data_coll);
else
    data_sent = labBroadcast(1);
end % if labindex==1
    
% unzipping the cell array
Mat_data = data_sent{1};
Out_bnd = data_sent{2};
In_bnd = data_sent{3};
Bnd_prop = data_sent{4};
Measure_reg = data_sent{5};
Thermal_grad = data_sent{6};
Measure_time = data_sent{7};
Param = data_sent{8};

% Defining debug files, can be commented out
DebugFile = fopen(['Monte_carlo' num2str(labindex) '.out'],'w');
fprintf(DebugFile, ['Hello! I am ' num2str(labindex) ' of ' num2str(numlabs) '\n']);


labBarrier;

% Reading equilibrium/linearization temperature
Teq = Param(4);

%% Importing material data
if labindex ==1
    if(isempty(Mat_data) || strcmp(Mat_data, File_mat))
        fprintf(DebugFile,'No material data is provided. Please check the input \n');
        return;
    end
    
    [dataSi,Impurity]=Load_Mat_data(Teq,Mat_data); % Loads material properties
                                        % after making suitable
                                        % adjustments for the
                                        % current temperature.
    dataSi_get = labBroadcast(1,dataSi); % Sending data to all the workers
else
    dataSi_get = labBroadcast(1); % recieving data on all the workers
end % if labindex ==1

if labindex==1
    Impurity_get = labBroadcast(1,Impurity);
else
    Impurity_get = labBroadcast(1);
end

dataSi = dataSi_get;
Impurity = Impurity_get;
labBarrier; % making sure that all workers get the data before proceeding

F = dataSi(:,1); % Frequencies radians/sec
V = dataSi(:,2); % velocities
tau = dataSi(:,3); % Relaxation times
C_data = dataSi(:,4); % Modal heat capacities
tau_inv = 1./tau; % relaxation rates
tau_im = dataSi(:,5); % impurity scattering rates
Nmodes = length(F);

% sanity check for the material data. Theoretical values of thermal conductivity
ktest = sum(C_data.*V.*V.*tau)/3;
fprintf(DebugFile,'The material thermal conductivity is %f \n',ktest);

%% Defining geometry domain
if(isempty(Out_bnd) || strcmp(Out_bnd,File_out_bnd))
    fprintf(DebugFile, 'Outer boundaries are not specified. Please check the input.\n');
    return;
end % if(isempty(Out_bnd) || strcmp(Out_bnd,

% The outer boundaries are fixed in terms of their orientation. Assigning their normals.
Normals =  [0,1,0;
            -1,0,0;
            0,-1,0;
            1,0,0;
            0,0,1;
            0,0,-1];

% Reading internal boundaries if they exist
if(isempty(In_bnd) || strcmp(In_bnd,File_in_bnd))
    fprintf(DebugFile, 'Inner boundaries are not specified. Homogeneous domain\n');
else
    Normals = [Normals; In_bnd(:,5:7)]; % adding to normals
    In_bnd = In_bnd(:,1:4); % keeping only x and y coordinates
end % if(isempty(In_bnd) || strcmp(In_bnd,

% Reading boundary properties
if(isempty(Bnd_prop) || strcmp(Bnd_prop,File_bnd_prop))
    fprintf(DebugFile, 'Boundary properties are not specified.\n');
    return;
else
    Prop_plane = Bnd_prop(:,2:5);
end % if(isempty(Bnd_prop) || strcmp,

% Reading reagions where output is requested
if(isempty(Measure_reg) || strcmp(Measure_reg,File_measure_reg))
    fprintf(DebugFile, 'No measurement regions are specified.\n');
    return;
else
    Region_data = Measure_reg;
end % if(isempty(Measure_reg) || strcmp(Measure_reg,

[num_bnd,~] = size(Prop_plane);

% reading volume of the domain from the file for now.
if(isempty(Param) || strcmp(Param,File_param))
    fprintf(DebugFile, 'Error in Simulation paramter file Sim_param.txt.\n');
    return;
else
    N = Param(1);
    max_scat = Param(2);
    Vol_dom = Param(3);
end % if(isempty(vol) || strcmp(vol,

%% Source terms
% Reading source terms other that isothermal boundaries. For now only thermal gradient is implemented

if(isempty(Thermal_grad) || strcmp(Thermal_grad,File_thermal_grad))
    fprintf(DebugFile, 'No additional source is specified.\n');
    num_sources = 0;
    Source_data = [];
else
    Source_data = Thermal_grad;
    [num_sources,~] = size(Source_data);
    
end % if(isempty(Thermal_grad) || strcmp(Thermal_grad,

% Reading measurement times if exist and deciding on 'Steady state' vs 'Transient' simulation
if(isempty(Measure_time) || strcmp(Measure_time,File_measure_time))
    Type = 'Steady';
     fprintf(DebugFile, 'Either no measurement times are specified or there is an error reading that file.\n Proceeding as Steady state simulation.\n');
else
    Type = 'Transient';
    time_points = Measure_time;
    [noTime,~] = size(time_points);
    fprintf(DebugFile, 'Measurement times are found. Proceeding with Transient simulation.\n');
end % if(isempty(Measure_time) || strcmp(Measure_time,


    
%****************************************************************************************************
% All data is read. Setting up sumulation

%% Deviational energy calculation

% Cumulative distributions calculation
cumul_base = zeros(Nmodes,1);
cumul_coll = zeros(Nmodes,1);
cumul_vel  = zeros(Nmodes,1);
cumul_base(1) = C_data(1);
cumul_coll(1) = C_data(1)*tau_inv(1);
cumul_vel(1) = C_data(1)*V(1);

for i=2:Nmodes
    cumul_base(i) = cumul_base(i-1)+C_data(i);
    cumul_coll(i) = cumul_coll(i-1)+C_data(i)*tau_inv(i);
    cumul_vel(i) = cumul_vel(i-1) + C_data(i)*V(i);
end
C = cumul_base(Nmodes); % Heat capacity at Teq

%Y=0, X=max, Y=max, X=0, Z=0, Z=max
% finding area of all the outer boundaries
Area = [Out_bnd(1)*Out_bnd(3);
        Out_bnd(2)*Out_bnd(3);
        Out_bnd(1)*Out_bnd(3);
        Out_bnd(2)*Out_bnd(3);
        Out_bnd(1)*Out_bnd(2);
        Out_bnd(1)*Out_bnd(2)];

% Deviational energy of prescribed boundaries
Pres_index = find(Prop_plane(:,1)==1);
if isempty(Pres_index)
    num_pres=0;
else
    [num_pres,~] = size(Pres_index);
    DevE_pres = zeros(num_pres,1);
     for ii=1:num_pres
         ID = Pres_index(ii);
         if(strcmpi(Type,'Steady'))
             DevE_pres(ii) = Area(ID)*abs(Teq-Prop_plane(ID,2))*cumul_vel(Nmodes)/4;
         else
             DevE_pres(ii) = time_points(end,1)*Area(ID)*abs(Teq-Prop_plane(ID,2))*cumul_vel(Nmodes)/4;
         end % if(strcmpi(Type,
         
     end % for ii=1:num_pres
     
end % if isempty(Pres)


% Deviational energy of body force
if num_sources~=0
    DevE_body = zeros(num_sources,1);
    for ii=1:num_sources
        grad_x = Source_data(ii,3);
        grad_y = Source_data(ii,4);
        grad_z = Source_data(ii,5);
        if(strcmpi(Type,'Steady'))
            DevE_body(ii) = sqrt(grad_x^2 + grad_y^2 + grad_z^2)*Vol_dom*cumul_vel(Nmodes)/2;
        else
            DevE_body(ii) = time_points(end,1)*sqrt(grad_x^2 + grad_y^2 + grad_z^2)*Vol_dom*cumul_vel(Nmodes)/2;
        end % if(strcmpi(Type,
        
    end % for ii=1:num_sources
    
end % if num_sources~=0

num_dev = num_pres + num_sources;

if (num_dev<1)
    % No sources found need to abort
    fprintf(DebugFile,'No sources are found!! exiting the simulation\n');
    return ;
else
    if num_pres==0
        DevE = DevE_body;
    elseif num_sources == 0
        DevE = DevE_pres;
    else
        DevE = [DevE_pres; DevE_body];
    end % if num_pres==0    
       
end % if (num_dev<1)

Cum_devE = cumsum(DevE);
Eeff = Cum_devE(end)/N;

% normalizing Deviational energy
Cum_devE = Cum_devE./Cum_devE(end);

prob = 0;
total_scat=0;
imp_scat =0;

% File to report real time progress at each worker
prog_ID = fopen(['Progress' num2str(labindex) '.txt'],'w');
labBarrier;

%% setting up the measurement variables 
Detector = Produce_detect(Region_data(:,1:6),Region_data(:,7));
% stores centroids of the Detectors, useful in transient cases
Centroids = [(Detector(:,1)+Detector(:,2))/2, (Detector(:,3)+Detector(:,4))/2, (Detector(:,5)+Detector(:,6))/2 ];

[noDetect, ~] = size(Detector);
if(strcmpi(Type,'Steady'))
    T = zeros(noDetect,Nmodes); % Temperature solutions
    Qx = zeros(noDetect,Nmodes); % Heat flux x component
    Qy = zeros(noDetect,Nmodes); % Heat flux y component
    Qz = zeros(noDetect,Nmodes); % Heat flux z component
elseif(strcmpi(Type,'Transient'))    
    T = zeros(noDetect,noTime); % Temperature solutions
    Qx = zeros(noDetect,noTime); % Heat flux x component
    Qy = zeros(noDetect,noTime); % Heat flux y component
    Qz = zeros(noDetect,noTime); % Heat flux z component
end % if(strcmpi(Type,

 tic;
for ii=labindex:numlabs:N % loop over N particles
    
    try % for debugging any edge case that is left out
        
        if(mod(ii*1000,N)==0)
            fprintf(prog_ID,'Current %f  \n',ii*100/N);
            disp(num2str(ii*100/N));
        end
        Scatt_count =0;
        
        % choosing source to emit particle from
        R = rand();
        for jj=1:num_dev
            if R<Cum_devE(jj)
                break;
            end % if R<Cum_devE(jj)
        end % for jj=1:num_sources
        
        if (jj<=num_pres) % emit from prescribed bnd
            mode = select_mode(cumul_vel,Nmodes); % calculating index of frequency
            speed = V(mode);
            [psign,Vx,Vy,Vz,x0,y0,z0] = emit_from_bnd(jj,Pres_index,Prop_plane(Pres_index,:),Out_bnd,Normals(Pres_index,:),speed,Teq);
            
        else % emit from body  source
            jj = jj-num_pres;
            mode = select_mode(cumul_vel,Nmodes); % calculating index of frequency
            speed = V(mode);
            [psign,Vx,Vy,Vz,x0,y0,z0] = emit_from_body(jj,Source_data,Detector,speed);
            
        end % if (jj<=num_pres) % emit from prescribed bnd
        
        finished = false; % current particle is active as long as this is false
        
        mode_change = true;
        imp_draw = true;

        if(strcmpi(Type,'Transient'))
            
            % selecting time for emission
            t0 = rand()*time_points(end);
            
            % Finding time coordinate for the start
            it=1;
            while (time_points(it)< t0)
                it = it + 1;
            end
            
        elseif(strcmpi(Type,'Steady'))
            t0 = 0; % Initialization always at zero pseudo-time
        end % if(strcmpi(Type,
        
        scat_type = 0;
        hit_bnd = 0;
        while ~finished
                        
            if(strcmpi(Impurity,'yes')) % impurity scattering is included
                if(mode_change)
                    Delt1 = -tau(mode)*log(1-rand());
                    Delt2 = -tau_im(mode)*log(1-rand()); % impurity scattering time
                elseif(imp_draw)
                    Delt2 = -tau_im(mode)*log(1-rand()); % impurity scattering time
                end % if(mode_change)
                
                total_scat = total_scat +1 ;
                
                if(Delt1<Delt2)
                    Delt=Delt1;
                    scat_type = 1;  % 1 for 3 phonon, 2 for isothermal, 3 for adiabatic 4 periodic 5 impurity
                else
                    Delt = Delt2;
                    scat_type = 5;
                end

                
            else % impurity scattering is not explicitly included set Delt2 very high

                if(mode_change)
                    Delt1 = -tau(mode)*log(1-rand());
                    Delt2 = 1; % impurity scattering time
                elseif(imp_draw)
                    Delt2 = 1; % impurity scattering time
                end % if(mode_change)
                
                total_scat = total_scat +1 ;
                
                if(Delt1<Delt2)
                    Delt=Delt1;
                    scat_type = 1;  % 1 for 3 phonon, 2 for isothermal, 3 for adiabatic 4 periodic 5 impurity
                else
                    Delt = Delt2;
                    scat_type = 5;
                end
                
            end % if(strcmpi(Impurity,
                       
            t1 = t0 + Delt;
            x1 = x0 + Vx*Delt;
            y1 = y0 + Vy*Delt;
            z1 = z0 + Vz*Delt;

            % Searching for encounter with external
            % boundaries. This redundent step is required because
            % in run time there are some edge cases thown by method
            % "Interact_out"
            
            in=true;
            if(x1>Out_bnd(1) || x1<0 || y1>Out_bnd(2) || y1<0 || z1>Out_bnd(3) || z1<0)
                in=false;
            end
            
            
            % Executed only if determined that particle is going outside
            if(~in)
                
                try % Catching exceptions for debugging purpose
                    [xout,yout,zout,frac_in,hit_bnd_out,scat_type_out] = Interact_out(x0,y0,z0,x1,y1,z1,Out_bnd,Prop_plane(1:6,1),Normals(1:6,:));
                catch e
                    fprintf(DebugFile,'Error with Check_intersect external, terminating particle\n');
                    fprintf(DebugFile,'particle id is = %d\n',ii);
                    fprintf(DebugFile,'The error information is as follows \n');
                    fprintf(DebugFile,'The identifier was: \n %s \n', e.identifier);
                    fprintf(DebugFile,'The error message was : \n %s \n', e.message);
                    finished = true;
                    prob=prob+1;
                    continue;
                    
                end % try % Catching exceptions for debugging purpose
                
                if(scat_type_out==0)
                    fprintf(DebugFile,'Scattering type is = %d\n',scat_type);
                end
                
                x1=xout;
                y1=yout;
                z1=zout;
                t1 = t0 + Delt*frac_in;
                hit_bnd = hit_bnd_out; % Outer boundaries are
                                       % stored first in the
                                       % sequence
                scat_type = scat_type_out;
%                 fprintf(DebugFile,'Scattering type is = %d\n',scat_type);
                
            end % if(~in)
            
            if(num_bnd>6)
                % Checking for interactions with internal boundaries
                % Currently only works for one hole in the unit
                % cell that which has its walls aligned with the
                % coordinate planes.
                
                try
                   
                    % finding the coordinates of the extents of
                    % internal boundaries
                    Xpt = unique(In_bnd(:,1));
                    Ypt = unique(In_bnd(:,2));
                    Zpt = [0;Out_bnd(3)];
                    
                    [isContr, len1] = findContr3D(x0,y0,z0,x1,y1,z1,Xpt,Ypt,Zpt);
                    
                    if(isContr && len1~=0) % particle crosses internal boundary
                        [frac_in, hit_bnd_in, scat_type_in,isReject] = Interact_in(x0,y0,x1,y1,Xpt,Ypt,Prop_plane(7:end,1),Normals(7:end,1:2));                      
                                                
% $$$                         % another piece for debugging
% $$$                         if(isempty(frac_in))
% $$$                             disp('empty variable');
% $$$                         end % if(isempty(frac_in))
                        
                        if(~isReject)
                            xin = x0 + frac_in*(x1-x0);
                            yin = y0 + frac_in*(y1-y0);
                            zin = z0 + frac_in*(z1-z0);
                            
                            x1=xin;
                            y1=yin;
                            z1=zin;
                            hit_bnd = hit_bnd_in + 6; % 6 outer
                                                      % boundaies
                                                      % are stored
                                                      % first
                            scat_type = scat_type_in;
                            t1 = t0 + frac_in*(t1-t0); % t1-t0 ~=Delt
                                                       % here. it could
                                                       % have been
                                                       % updated by
                                                       % outer boundary
                        end % if(~isReject)
                        
                    end % if(isContr && len1~=0) 
                    
                

                catch e
                    fprintf(DebugFile,'Error with internal bnd, terminating particle\n');
                    fprintf(DebugFile,'particle id is = %d\n',ii);
                    fprintf(DebugFile,'The error information is as follows \n');
                    fprintf(DebugFile,'The identifier was: \n %s \n', e.identifier);
                    fprintf(DebugFile,'The error message was : \n %s \n', e.message);
                    finished = true;
                    prob=prob+1;
                    continue;
                end % try
                
                
            end % if(num_bnd>6)
            
            % Test for interface
            if((interface-z0)*(interface-z1)<0)
                % particle crosses interface
                scat_type = 6;
                frac = (interface-z0)/(z1-z0);
                x1 = x0 + frac*(x1-x0);
                y1 = y0 + frac*(y1-y0);
                z1 = interface;
                t1 = t0 + frac*(t1-t0);
            end
            
            % Sampling step depends on whether it is a steady
            % state simulation or transient one
            
            if (strcmpi(Type,'Steady'))
                % It is a Steady-state simulation
                
                % In most cases when we perform steady state
                % simulation, we also require modal
                % decomposition. It is easy to turn that off by
                % editing the following code 
                
                for jj=1:noDetect
                    Xpt=[Detector(jj,1);Detector(jj,2)];
                    Ypt=[Detector(jj,3);Detector(jj,4)];
                    Zpt=[Detector(jj,5);Detector(jj,6)];
                    [isContr, len] = findContr3D(x0,y0,z0,x1,y1,z1,Xpt,Ypt,Zpt);
                    
% $$$                     if(isnan(len)) % For debugging
% $$$                         fprintf(DebugFile,'I have got a bad contribution\n');
% $$$                     end % if(isnan(len)) % For debugging
                    
                    if isContr==1
                        Dx = abs(Detector(jj,1)-Detector(jj,2));
                        Dy = abs(Detector(jj,3)-Detector(jj,4));
                        Dz = abs(Detector(jj,5)-Detector(jj,6));
                        
                        T(jj,mode) = T(jj,mode) + psign*Eeff*len/C/(Dx*Dy*Dz)/V(mode); % temperature
                        Qx(jj,mode) = Qx(jj,mode) + psign*Eeff*len*Vx/(Dx*Dy*Dz)/V(mode); % heat flux
                        Qy(jj,mode) = Qy(jj,mode) + psign*Eeff*len*Vy/V(mode)/(Dx*Dy*Dz); % heat flux
                        Qz(jj,mode) = Qz(jj,mode) + psign*Eeff*len*Vz/V(mode)/(Dx*Dy*Dz); % heat flux
                    end % if isContr==1
                    
                end % for jj=1:noDetect
                
                
                
            elseif(strcmpi(Type,'Transient'))
                % It is Transient simulation
                
                while( it<(length(time_points)+1) && time_points(it)<=t1 )
                    
                    X_c = x0 + Vx*(time_points(it) - t0);
                    Y_c = y0 + Vy*(time_points(it) - t0);
                    Z_c = z0 + Vz*(time_points(it) - t0);
                    
                    % finding the cell to which the particle belong
                    index_pos = dsearchn(Centroids,[X_c,Y_c,Z_c]);
                    
                    Dx = abs(Detector(index_pos,1)-Detector(index_pos,2));
                    Dy = abs(Detector(index_pos,3)-Detector(index_pos,4));
                    Dz = abs(Detector(index_pos,5)-Detector(index_pos,6));
                    
                    % Sampling the observables
                    T(index_pos,it) = T(index_pos, it) + psign*Eeff/C/(Dx*Dy*Dz);
                    Qx(index_pos,it) = Qx(index_pos,it) + psign*Eeff*Vx/(Dx*Dy*Dz);
                    Qy(index_pos,it) = Qy(index_pos,it) + psign*Eeff*Vy/(Dx*Dy*Dz);
                    Qz(index_pos,it) = Qz(index_pos,it) + psign*Eeff*Vz/(Dx*Dy*Dz);
                    
                    it = it +1;
                    
                end % while(time_points(it)<=t1)
                
                
            end % if(strcmpi(Type,
            

            % Scattering implementation
            
            switch scat_type
                
                case 1 % Three phonon process
                    
                    % select post-collision mode
                    mode = select_mode(cumul_coll,Nmodes);
                    
                    % Time to next 3 phonon scattering will change
                    % Time to next impurity scattering will also change
                    mode_change = true;
                    imp_draw = false;
                    
                    % Sampling new random direction of travel
                    Theta = 2*pi()*rand(); % azimuthal angle
                    cos_phi = 2*rand()-1; % Polar angle
                    
                    Vx = V(mode)*cos_phi;
                    Vy = V(mode)*sqrt(1-cos_phi^2)*cos(Theta);
                    Vz = V(mode)*sqrt(1-cos_phi^2)*sin(Theta);
                    x0=x1;
                    y0=y1;
                    z0=z1;
                    t0=t1;
                    Scatt_count = Scatt_count +1;
                    
                case 2 % isothermal boundary
                    finished = true;
                    
                case 3 % adiabatic boundary is hit
                    mode_change = false;
                    imp_draw = false;
                    
                    % Updating time available for next collision
                    Delt1 = Delt1 - (t1-t0);
                    Delt2 = Delt2 - (t1-t0);
                    
                    spec = Prop_plane(hit_bnd,2);
                    nx = Normals(hit_bnd,1);
                    ny = Normals(hit_bnd,2);
                    nz = Normals(hit_bnd,3);
                    
                    if (rand()<spec ) % specular reflection
                        
                        V_new = reflection([Vx;Vy;Vz],[nx;ny;nz]);
                        
                        Vx = V_new(1); Vy=V_new(2); Vz=V_new(3);
                    else % diffusive reflection
                        
                        % Finding random directions
                        Theta = 2*pi()*rand(); % azimuthal angle
                        
                        U = rand(); % for sampling polar angle
                        
                        cos_phi = sqrt(U);
                        sin_phi = sqrt(1-U);
                        
                        % mode does not change
                        VX = V(mode)*cos_phi;
                        VY = V(mode)*sin_phi*cos(Theta);
                        VZ = V(mode)*sin_phi*sin(Theta);
                        
                        % Orienting so that emitted properly.
                        Rot = Orient([1;0;0],[nx;ny;nz]);
                        V_new = Rot*[VX;VY;VZ];
                        Vx=V_new(1); Vy=V_new(2); Vz=V_new(3);
                        
                    end % if (rand()<spec )
                    
                    x0=x1;
                    y0=y1;
                    z0=z1;
                    t0=t1;
                    
                    % Scattering count is not increased as this
                    % does not result in relaxation
                    
                    
                case 4 % Periodic boundary
                    
                    % this implementation assumes that Teq is same
                    % at both the boudaries in pair
                    mode_change = false;
                    imp_draw = false;
                    % Updating time available for next collision
                    Delt1 = Delt1 - (t1-t0);
                    Delt2 = Delt2 - (t1-t0);
                    
                    % everything will be same only the translation vector will
                    % be added to final position
                    x0 = x1 + Prop_plane(hit_bnd,2);
                    y0 = y1 + Prop_plane(hit_bnd,3);
                    z0=z1 + Prop_plane(hit_bnd,4);
                    t0 = t1;
                    
                    % Again this scattering does not result in relaxation
                    
                case 5 % Impurity scattering
                    
                    % this only randomizes the travelling direction
                    % without changing the mode
                    mode_change = false;
                    imp_draw = true;
                    % Updating time available for next collision
                    Delt1 = Delt1 - (t1-t0);
                    % New impurity scattering time will be drawn
                    
                    Theta = 2*pi()*rand();
                    cos_phi = 2*rand()-1;
                    
                    Vx = V(mode)*cos_phi;
                    Vy = V(mode)*sqrt(1-cos_phi^2)*cos(Theta);
                    Vz = V(mode)*sqrt(1-cos_phi^2)*sin(Theta);
                    x0=x1;
                    y0=y1;
                    z0=z1;
                    t0=t1;
                    
                    imp_scat = imp_scat +1;
                    % this scattering also does not result in relaxation
                    
                case 6  % interface scattering
                    % apply interface scattering changes
                    mode_change = false;
                    imp_draw = false;
                    % Updating time available for next collision
                    Delt1 = Delt1 - (t1-t0);
                    Delt2 = Delt2 - (t1-t0);
                    
                    % transmission/reflection, specular/diffusive
                    [Vx,Vy,Vz]=interface_implement(Vx,Vy,Vz,F(mode));
                    
                    x0=x1;
                    y0=y1;
                    z0=z1;
                    t0=t1;
                otherwise % added for debugging
                    fprintf(DebugFile,'No scattering event is found. \n ');
            end % switch scat_type


            % Other termination conditions

            if(strcmpi(Type,'Transient'))
                if(t0 >= time_points(end))
                    finished = true;
                end % if(t0 >= time_points(end))
                
            elseif (strcmpi(Type,'Steady'))
                if(Scatt_count > max_scat)
                    finished = true;
                end % if(Scatt_count > max_scat)
                
            end % if(strcmpi(Type,
            
        end % while ~finished
        
    catch e
        fprintf(DebugFile,'Some problem occured in runtime, particel id= %d \n',ii);
        fprintf(DebugFile,'The error information is as follows \n');
        fprintf(DebugFile,'The identifier was: \n %s \n', e.identifier);
        fprintf(DebugFile,'The error message was : \n %s \n', e.message);
        fprintf(DebugFile,getReport(e));
        fprintf(DebugFile,'\n');
        finished = true;
        prob=prob+1;

        continue;

    end % try % for debugging any edge case that is left out

end % for ii=labindex:numlabs:N % loop over N particles

time = toc;

labBarrier;

fclose(prog_ID);

fprintf(DebugFile,'The calculation is done\n');

%% Collecting all data and cosolidating results

try % For debugging issues in collecting the data and writing it to the files

    if labindex == 1
        for ii=2:numlabs
            GetT = labReceive(ii,1);
            GetQx = labReceive(ii,2);
            GetQy = labReceive(ii,3);
            GetQz = labReceive(ii,4);
            Getprob = labReceive(ii,5);
            
            T = T + GetT;
            Qx = Qx + GetQx;
            Qy = Qy + GetQy;
            Qz = Qz + GetQz;
            prob = prob + Getprob;
            
        end % for ii=2:numlabs
        
    else
        labSend(T,1,1);
        labSend(Qx,1,2);
        labSend(Qy,1,3);
        labSend(Qz,1,4);
        labSend(prob,1,5);
        
    end % if labindex == 1
    
    % Writing all results to txt files
    
    if labindex == 1
        
        Tempfile = fopen(['T' num2str(Teq) '.txt'], 'w');
        Qxfile = fopen(['Qx' num2str(Teq) '.txt'], 'w');
        Qyfile = fopen(['Qy' num2str(Teq) '.txt'], 'w');
        Qzfile = fopen(['Qz' num2str(Teq) '.txt'], 'w');
        
        if(strcmpi(Type,'Transient'))
            data_points = length(time_points);
            add_T = Teq;
        elseif (strcmpi(Type, 'Steady'))
            data_points = Nmodes;
            add_T = 0;
        end % if(strcmpi(Type,
        
        for jj=1:noDetect
            for kk=1:data_points
                fprintf(Tempfile, '%15.5e', T(jj,kk) + add_T);
                fprintf(Qxfile, '%15.5e', Qx(jj,kk));
                fprintf(Qyfile, '%15.5e', Qy(jj,kk));
                fprintf(Qzfile, '%15.5e', Qz(jj,kk));
                
            end % for kk=1:data_points
            
            fprintf(Tempfile, '\n');
            fprintf(Qxfile, '\n');
            fprintf(Qyfile, '\n');
            fprintf(Qzfile, '\n');
            
        end % for jj=1:noDetect
        
        fclose(Tempfile);
        fclose(Qxfile);
        fclose(Qyfile);
        fclose(Qzfile);
        
        % Writing down detector locations 
        detectfile = fopen('detector_location.txt','w');
        
        for jj=1:noDetect
            fprintf(detectfile, '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n', Detector(jj,:));
        end % for jj=1:noDetect
        
        fclose(detectfile);
        
        fid = fopen('problem_count.txt','w');
        fprintf(fid, '%d\n', prob);
        fclose(fid);
        
    end % if labindex == 1
        
        
    
catch e
     fprintf(DebugFile,'Some problem occured in runtime, terminating the run and saving results till now. \n');
    fprintf(DebugFile,'The error information is as follows \n');
    fprintf(DebugFile,'The identifier was: \n %s \n', e.identifier);
    fprintf(DebugFile,'The error message was : \n %s \n', e.message);
    fprintf(DebugFile,getReport(e));
    fprintf(DebugFile,'\n');
    
end % try % For debugging issues in collecting the data and writing it to the files
    
fclose(DebugFile);



end % function [time]=BTE_solution_3D(varargin)


