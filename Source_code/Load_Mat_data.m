function [Data_out, Impurity]=Load_Mat_data(Teq,Mat_data)
    %{
      The source of material properties that we have contains following two possible format
      1) Either 6 columns with Omega,DOS,Velocity,dOmega,tau_inv,tau,polarization format
      2) Or 4 columns with Omega, Velocity, tau, specific heat
      If impurity scattering times are specified, they will be listed in an additional column in the end.
      Both are calculated at 300K temperature so need to adjust for desired temperature
      
    %}
    
            % define reduced Planck constant and Boltzmann constant
            hbar=1.054517e-34; % J s = m^2 kg s-1
            boltz=1.38065e-23; % m2 kg s-2 K-1
                               %    Mat_data = load('dataSi.txt');
            Impurity = 'No';
    [~,col] = size(Mat_data);

    if(col==6 || col==7)
        type=1; % The data needs to be converted
        if(col==7)
            Impurity = 'Yes';
        end
    elseif(col==4 || col==5)
        type=2; % The data only needs adjustment for temperature
        if(col==5)
            Impurity = 'Yes';
        end
    end % if(col==4)
        
    switch type
      case 1
        F = Mat_data(:,1); % frequencies radian/sec
        SD = Mat_data(:,2);  % Density of states
        V = Mat_data(:,3); % velocities
        Dom = Mat_data(:,4); % Delta frequencies
        tau = Mat_data(:,5); % relaxation times
        de_dT = (hbar*F/Teq).^2/boltz.*exp(hbar*F/boltz/Teq)./(exp(hbar*F/boltz/Teq)-1).^2; %derivative of Bose-Einstein
        C_data = SD.*Dom.*de_dT;
      case 2
        F = Mat_data(:,1);
        V = Mat_data(:,2);
        tau = Mat_data(:,3);
        C_data = Mat_data(:,4);        
        % Although an adjustment needs to be applied to relaxation
        % time as well but for now it has not been done. For more
        % sophisticated data available at different temperature
        % this adjustement won't even be required.
        
      otherwise
        fprintf(2,"Unrecognized material data input file format. Please see the text file again \n");
    end % switch type
        

    if(strcmpi(Impurity,'Yes'))
        tau_im = Mat_data(:,col);%1./((2e-44)*F.^4); % impurity scattering time scale
    elseif(strcmpi(Impurity,'No'))
        tau_im = (-1)*ones(length(F),1); 
        % dummy value of -1 to keep the material property array
        % dimensions consistent.
    else
        fprintf(2,"Not able to determine if impurity scattering to be added or not.\n Specify the variable correctly");
    end % if(~strcmpi(Impurity,
    
    Data_out = [F,V,tau,C_data,tau_im];
    end % function [dataSi]=Load_Mat_data(Teq)
    
    
    %{
            % Calculating heat capacity
        % in high temperature limit C=3*kB
        Ref_temp = 5000; % 5000K as high temperature limit
        X_hw_ref = hbar*F/(boltz*REf_temp);
        e_X_ref = exp(X_hw);
        C_data_ref = fact*boltz*(X_hw.^2).*e_X./(e_X-1).^2;
        %The factor required to represent finite
        % volume represented by mesh in q space while calculating
        % DFT data
        fact =  3*boltz/sum(C_data_ref);
        %Calculating adjusted data for modal heat capacity
        X_hw = hbar*F/(boltz*Teq);
        e_X = exp(X_hw);
        C_data = fact*boltz*(X_hw.^2).*e_X./(e_X-1).^2;
    %}