function [frac_in, hit_bnd_in, scat_type_in,isReject] = Interact_in(x0,y0,x1,y1,Xpt,Ypt,Prop_plane_in,Normals)

    %{
      This method checks for interaction of particle trajectory with internal boundaries
      It returns 'frac_in', the fraction of trajectory that is inside the domain after curtailment due to boundary interaction
    %}
isReject = false;
time = [[1;2;3;4],100*ones(4,1)];
tx = (Xpt-x0)./(x1-x0);
ty = (Ypt-y0)./(y1-y0);

% time of flight to boundaries
time(1,2) = ty(1);
time(2,2) = tx(2);
time(3,2) = ty(2);
time(4,2) = tx(1);

time = time(time(:,2)>0,:);
time = time(time(:,2)<1,:);

[frac_in, index] = min(time(:,2));
hit_bnd_in = time(index,1);
scat_type_in = Prop_plane_in(hit_bnd_in) + 1;

% There are few edge cases when the code takes the particle just scattered
% from the boundary and scatters it again thinking it is coming to boundary
% while it is actually going away. The following flags that
Flight_vec = [x1-x0,y1-y0]'; % default vector is column vector
if(dot(Flight_vec,Normals(hit_bnd_in,:)')>0)
    isReject = true;
end
