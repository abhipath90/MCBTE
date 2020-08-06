function [psign,Vx,Vy,Vz,x0,y0,z0] = emit_from_bnd(jj,Pres_index,Props,Out_bnd_data,Normals,speed,Teq)

    %{
      This method assigns initial properties to a computational particle when it is emitted from one of the prescribed boundaries.
      'jj' tells information about which boundary the particle will originate from      
    %}
    
    
% jj is the relative index of origination boundary among Prescribed boundaries.
% Pres_index stores indices of all the prescribed boundaries
    
    origin = Pres_index(jj);
    R1 = rand();
    R2 = rand();
    
    % calculating initial location of particle
    switch origin
      
      case 1
         % its xz plane
        x0 = R1*Out_bnd_data(1);
        y0 = 0;
        z0 = R2*Out_bnd_data(3);
        
      case 2
        % its yz plane
        x0 = Out_bnd_data(1);
        y0 = R1*Out_bnd_data(2);
        z0 = R2*Out_bnd_data(3);
        
      case 3
        % its xz plane
        x0 = R1*Out_bnd_data(1);
        y0 = Out_bnd_data(2);
        z0 = R2*Out_bnd_data(3);
        
      case 4
        % its yz plane
        x0 = 0;
        y0 = R1*Out_bnd_data(2);
        z0 = R2*Out_bnd_data(3);
        
      case 5
        % its z=0 plane
        x0 = R1*Out_bnd_data(1);
        y0 = R2*Out_bnd_data(2);
        z0 = 0;
        
      case 6
        % its z=max_z plane
        x0 = R1*Out_bnd_data(1);
        y0 = R2*Out_bnd_data(2);
        z0 = Out_bnd_data(3);
        
    end % switch jj
    
    Tbc = Props(jj,2);    
    
    % Calculating initial travelling directions 
    Theta = 2*pi()*rand();
    U = rand();
    cos_phi = sqrt(U);
    sin_phi = sqrt(1-U);
    
    VX = speed*cos_phi;
    VY = speed*sin_phi*cos(Theta);
    VZ = speed*sin_phi*sin(Theta);
    
    % Orienting so that emitted properly.
    Rot = Orient([1;0;0],(Normals(jj,:))');
    V_new = Rot*[VX;VY;VZ];
    Vx=V_new(1); Vy=V_new(2); Vz=V_new(3);
    
    psign = sign(Tbc-Teq);
    
end % function             [psign,
