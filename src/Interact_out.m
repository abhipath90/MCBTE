function [xout,yout,zout,frac_in,hit_bnd_out,scat_type_out] = Interact_out(x0,y0,z0,x1,y1,z1,Out_bnd_data,plane_type,Normals)
X_high = Out_bnd_data(1);
Y_high = Out_bnd_data(2);
Z_high = Out_bnd_data(3);

xout = x1; yout=y1; zout=z1; frac_in=1; hit_bnd_out=0; scat_type_out=0;
flight_time = zeros(length(plane_type),2);
flight_time(:,1) = 1000;
flight_time(:,2) = (1:length(plane_type))';


% The order of planes is y=0;x=W;y=W;x=0;z=0;z=thick;
if(y1~=y0)
    flight_time(1,1) = (0-y0)/(y1-y0);
    flight_time(3,1) = (Y_high-y0)/(y1-y0);
end
if(x1~=x0)
    flight_time(4,1) = (0-x0)/(x1-x0);
    flight_time(2,1) = (X_high-x0)/(x1-x0);
end

if(z1~=z0)
    flight_time(5,1) = (0-z0)/(z1-z0);
    flight_time(6,1) = (Z_high-z0)/(z1-z0);
end


Flight_vec = [x1-x0,y1-y0,z1-z0]'; % default vector is column vector
% Repeating trajectory vector to calculate dot product with each normal
Flight_mat = repmat(Flight_vec,[1,length(plane_type)]);

% it can only interact with planes with which it has negative dot product.
flight_time = flight_time(dot(Flight_mat,Normals')<0,:);


flight_time = flight_time(flight_time(:,1)>0,:);
flight_time = flight_time(flight_time(:,1)<1,:);
if(~isempty(flight_time))
    [frac_in,index]=min(flight_time(:,1));
    hit_bnd_out = flight_time(index,2);
    scat_type_out = plane_type(hit_bnd_out,1) + 1;
    xout = x0 + frac_in*(x1-x0);
    yout = y0 + frac_in*(y1-y0);
    zout = z0 + frac_in*(z1-z0);
end