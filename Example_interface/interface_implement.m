function [Vx,Vy,Vz]=interface_implement(Vx,Vy,Vz,F)

V = sqrt(Vx^2+Vy^2+Vz^2);
% interface normal in the direction of travel
nx = 0; ny = 0;
nz = sign(Vz);
%% define interface conditions here
spec_inter = 1; % degree of specularity of interface

%% defining transimission probability
% for high pass filter
  omega_cutoff = 0.6*(2*pi()*12e12); % 0.6 of max frequency of 12 THz
  if(F>omega_cutoff)
      T = 1;
  else
      T = 0;
  end

% % incident angle based filter
% % since the interface is perpendicular to z direction
% T = abs(Vz/V);

trans = rand();
spec = rand();

%%% ****** Block based on what Hao et al. did in 2009 paper
% $$$ 
% $$$     Theta = 2*pi()*rand();
% $$$     cos_phi = 2*rand()-1;
% $$$     % Randomizing the direction of travel
% $$$     Vx = V(mode)*cos_phi;
% $$$     Vy = V(mode)*sqrt(1-cos_phi^2)*cos(Theta);
% $$$     % since interface is in z=a, a being variable
% $$$     % special attention is needed for the sign of Vz
% $$$     Vz_new = V(mode)*sqrt(1-cos_phi^2)*sin(Theta);
% $$$ if(trans<T) % transmission    
% $$$     Vz = sign(Vz)*abs(Vz_new);
% $$$ else % reflection
% $$$     Vz = -sign(Vz)*abs(Vz_new);
% $$$ end

    

%%% ******* Block based on our earlier understanding

% $$$ 
  if(trans<T) % transmission
      if(spec<spec_inter) % specular
          % No change in the velocity
      else
          % randomize velocity in the direction of travel       
          Theta = 2*pi()*rand(); % azimuthal angle
          U = rand(); % for sampling polar angle
          cos_phi = sqrt(U);
          sin_phi = sqrt(1-U);
          
          % mode does not change
          VX = V*cos_phi;
          VY = V*sin_phi*cos(Theta);
          VZ = V*sin_phi*sin(Theta);
          
          % Orienting so that emitted properly.
          Rot = Orient([1;0;0],[nx;ny;nz]);
          V_new = Rot*[VX;VY;VZ];
          Vx=V_new(1); Vy=V_new(2); Vz=V_new(3);
      end
  else
      % reflection
      nz = -nz; % will reflect in opposite direction than travel
      if(spec<spec_inter) % specular
          V_new = reflection([Vx;Vy;Vz],[nx;ny;nz]);
          Vx = V_new(1); Vy=V_new(2); Vz=V_new(3);
      else
          % diffusive
          % randomize velocity in the direction of travel       
          Theta = 2*pi()*rand(); % azimuthal angle
          U = rand(); % for sampling polar angle
          cos_phi = sqrt(U);
          sin_phi = sqrt(1-U);
          
          % mode does not change
          VX = V*cos_phi;
          VY = V*sin_phi*cos(Theta);
          VZ = V*sin_phi*sin(Theta);
          
          % Orienting so that emitted properly.
          Rot = Orient([1;0;0],[nx;ny;nz]);
          V_new = Rot*[VX;VY;VZ];
          Vx=V_new(1); Vy=V_new(2); Vz=V_new(3);
      end
  end