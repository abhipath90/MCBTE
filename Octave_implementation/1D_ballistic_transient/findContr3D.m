function [isContr, len] = findContr3D(x0,y0,z0,x1,y1,z1,Xpt,Ypt,Zpt)

    %{
      This method calculates if the particle trajectory intersects with 
      the current spatial cell, and further calculates the length of 
      the trajectory inside the spatial cell
      
    %}

isContr=0;
len=0;
len_sec = sqrt((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2);                                                                    
% direction cosins                                                                                                    
lx=(x1-x0)/len_sec; ly=(y1-y0)/len_sec; lz=(z1-z0)/len_sec;                                                           
                                                                                                                      
% extents of the box                                                                                                  
Xmax = max(Xpt); Xmin = min(Xpt);                                                                                     
Ymax = max(Ypt); Ymin = min(Ypt);                                                                                     
Zmax = max(Zpt); Zmin = min(Zpt);                                                                                     
                                                                                                                      
if (lx>=0)                                                                                                            
    txmin = (Xmin - x0)/lx;                                                                                           
    txmax = (Xmax - x0)/lx;                                                                                           
else                                                                                                                  
    txmin = (Xmax - x0)/lx;                                                                                           
    txmax = (Xmin - x0)/lx;                                                                                           
end                                                                                                                   
                                                                                                                      
tmin = txmin;                                                                                                         
tmax = txmax;                                                                                                         
                                                                                                                      
if (ly>=0)                                                                                                            
    tymin = (Ymin - y0)/ly;                                                                                           
    tymax = (Ymax - y0)/ly;                                                                                           
else                                                                                                                  
    tymin = (Ymax - y0)/ly;                                                                                           
    tymax = (Ymin -y0)/ly;                                                                                            
end                                                                                                                   
                                                                                                                      
if ( (tmin>tymax) || (tymin > tmax))                                                                                  
    isContr =0;                                                                                                       
    len = 0;                                                                                                          
    return;                                                                                                           
end                                                                                                                   
                                                                                                                      
if ( tymin > tmin)                                                                                                    
    tmin = tymin;                                                                                                     
end                                                                                                                   
                                                                                                                      
if (tymax < tmax)                                                                                                     
    tmax = tymax;                                                                                                     
end                                                                                                                   
                                                                                                                      
if (lz>=0)                                                                                                            
    tzmin = (Zmin - z0)/lz;                                                                                           
    tzmax = (Zmax -z0)/lz;                                                                                            
else                                                                                                                  
    tzmin = (Zmax - z0)/lz;                                                                                           
    tzmax = (Zmin - z0)/lz;                                                                                           
end                                                                                                                   
                                                                                                                      
if ( (tmin > tzmax) || (tzmin > tmax) )                                                                               
    isContr = 0;                                                                                                      
    len = 0;                                                                                                          
    return;                                                                                                           
end                                                                                                                   
                                                                                                                      
if (tzmin > tmin)                                                                                                     
    tmin = tzmin;                                                                                                     
end                                                                                                                   
                                                                                                                      
if (tzmax < tmax)                                                                                                     
    tmax = tzmax;                                                                                                     
end                      

% The detector is in the path of the ray. The actual intersection depends on the values of tmin and tmax              
                                                                                                                      
% Calculating length of the intersection                                                                              
% Box does not lie inside the segment                                                                                 
if (tmin > len_sec || tmax <0)                                                                                        
    isContr =0;                                                                                                       
    len = 0;                                                                                                          
    return;                                                                                                           
end                                                                                                                   
% Box emcompasses the whole segment                                                                                   
if ( (tmin<=0) && (tmax >= len_sec))                                                                                  
    isContr = 1;                                                                                                      
    len = len_sec;                                                                                                    
    return;                                                                                                           
end                                                                                                                   
                                                                                                                      
% Currently tmin<len_sec, tmax >0                                                                                     
% Segment cuts box twice                                                                                              
if ( (tmin>0) && (tmax<len_sec))                                                                                      
    isContr = 1;                                                                                                      
    len = tmax - tmin ;                                                                                               
    return;                                                                                                           
elseif ((tmin>0) && (tmax>=len_sec))                                                                                  
    % segment enters but doesn't go out i.e. ends inside                                                              
    isContr=1;                                                                                                        
    len = len_sec - tmin;                                                                                             
    return;                                                                                                           
elseif ((tmin<=0) && (tmax<len_sec))                                                                                  
    % segment starts inside and goes out                                                                              
    isContr = 1;                                                                                                      
    len = tmax;                                                                                                       
    return;                                                                                                           
else                                                                                                                  
    % throw error                                                                                                     
    isContr = [];                                                                                                     
    len = [];                                                                                                         
end                                                                                                                   
                  
