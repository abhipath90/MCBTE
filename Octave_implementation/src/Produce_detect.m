function Detectors = Produce_detect(regions,refinement)
% Input definition
% The simulation domain contains cuboidal regions with lower and upper bound on x, y, and z
% coordinates
% The refinement index indicates the desired level of refinement given by 2^{3n}. Each level of
% refinement produces 8 cells for each region.

[row,~] = size(regions);
if(row~=length(refinement))
     % Throw exception
end

No_detect = sum(8.^refinement);
Detectors = zeros(No_detect,6);
count =1;
for ii=1:row
    [Add,num] = refine_region(regions(ii,:),refinement(ii,1));
    Detectors(count:(count+num-1),:) = Add;
    count = count + num;
end
