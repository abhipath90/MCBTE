function [psign,Vx,Vy,Vz,x0,y0,z0] = emit_from_body(jj,Source_data,Detector,speed)

    %{
      This method assigns initial properties to the particle when it originates from body source (Thermal gradient)
            
    %}
    
    [noDetect,~] = size(Detector);
    psign = sign(rand()-0.5); % particle can have positive or negative sign equally likely

    % finding random location in the domain
    cell = randi(noDetect);
    x0 = Detector(cell,1) + rand()*(Detector(cell,2)-Detector(cell,1));
    y0 = Detector(cell,3) + rand()*(Detector(cell,4)-Detector(cell,3));
    z0 = Detector(cell,5) + rand()*(Detector(cell,6)-Detector(cell,5));

    % selecting initial travelling direction
    Theta = 2*pi()*rand();
    U = rand();
    cos_phi = sqrt(U);
    sin_phi = sqrt(1-U);
    
    % - sign is because body force is always in opposite direction
    % to the thermal gradient
    VX = -psign*speed*cos_phi;
    VY = -psign*speed*sin_phi*cos(Theta);
    VZ = -psign*speed*sin_phi*sin(Theta);
    
    grad_x = Source_data(jj,3); 
    grad_y = Source_data(jj,4);
    grad_z = Source_data(jj,5);
    % Orienting so that gradient biases are taken into account.
    Rot = Orient([1;0;0],[grad_x;grad_y;grad_z]);
    V_new = Rot*[VX;VY;VZ];
    Vx=V_new(1); Vy=V_new(2); Vz=V_new(3);

    
end % function
