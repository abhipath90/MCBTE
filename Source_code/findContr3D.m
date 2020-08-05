function [isContr, len] = findContr3D(x0,y0,z0,x1,y1,z1,Xpt,Ypt,Zpt)
isContr=0;
len=0;
% X_pt = [X-Dx/2; X+Dx/2];
% Y_pt = [Y-Dy/2; Y+Dy/2];
% Z_pt = [Z-Dz/2; Z+Dz/2];

tx = (Xpt-x0)./(x1-x0);
ty = (Ypt-y0)./(y1-y0);
tz = (Zpt-z0)./(z1-z0);

maxes = [max(tx);max(ty);max(tz)];
mins = [min(tx); min(ty); min(tz)];
minofmax = min(maxes);
maxofmin = max(mins);
if(minofmax<maxofmin)
    % Line will not pass through the cell
    isContr=0;
    len=0;
    return;
else
    % The line will pass through the cell
    % now checking if it will pass within the flight time
    if(minofmax<0 || maxofmin>1)
        % particle will not intersect in flight time
        isContr = 0;
        len=0;
        return;
    elseif(minofmax>1)
        % set it to 1 as its ending inside the cell
        minofmax = 1;
    elseif(maxofmin<0)
        % set it to 0 as its starting from the cell
        maxofmin=0;
    end
end
isContr=1;
delta_t = minofmax - maxofmin;
seg_len = sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2);
len = delta_t*seg_len;