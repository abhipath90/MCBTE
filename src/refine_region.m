function [new_regions, num] = refine_region(region,order)

    %{
      This method subdivides spatial cells into smaller ones based on the specified order.
      It returns new_regions and num of new regions created.
      
    %}
    
% Adding this part to deal with fully defined regions specified with order
% 0
if(order==0)
    new_regions = region;
    num = 1;
    
else
    % calculate subregions
    X0=region(1,1); X1=region(1,2); DX = X1-X0;
    Y0=region(1,3); Y1=region(1,4); DY = Y1-Y0;
    Z0=region(1,5); Z1=region(1,6); DZ = Z1-Z0;
    
    subreg = [X0,X0+DX/2,Y0,Y0+DY/2,Z0,Z0+DZ/2;
        X0+DX/2,X0+DX,Y0,Y0+DY/2,Z0,Z0+DZ/2;
        X0,X0+DX/2,Y0+DY/2,Y0+DY,Z0,Z0+DZ/2;
        X0+DX/2,X0+DX,Y0+DY/2,Y0+DY,Z0,Z0+DZ/2;
        X0,X0+DX/2,Y0,Y0+DY/2,Z0+DZ/2,Z0+DZ;
        X0+DX/2,X0+DX,Y0,Y0+DY/2,Z0+DZ/2,Z0+DZ;
        X0,X0+DX/2,Y0+DY/2,Y0+DY,Z0+DZ/2,Z0+DZ;
        X0+DX/2,X0+DX,Y0+DY/2,Y0+DY,Z0+DZ/2,Z0+DZ];
    
    
    if(order==1)
        % return refined regions with num=8
        new_regions = subreg;
        num=8;
    else % recursive call
        % call the function with order-1
        % collect the regions and num
        num = 8^order;
        checksum=0;
        new_regions = zeros(num,6);
        for ii=1:8
            [getregion, getnum] = refine_region(subreg(ii,:),order-1);
            new_regions((checksum+1): (checksum + getnum),:) = getregion;
            checksum = checksum + getnum;
        end
        
    end
end

end

