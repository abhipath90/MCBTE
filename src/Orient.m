function R = orient(ref, dest)

% this code generates roation matrix to rotate a vector ref to vector dest
% This is used many times to orient the velocity vector appropriately when
% emmited or diffusively reflected from a plane
% WORKS FOR ANY DIMENSION

% This procedure only works when ref~=-(dest);

if(~all(ref==-dest))
    u = ref/norm(ref);
    v = dest/norm(dest);
    N = length(u);
    S = reflection(eye(N),v+u); % SU=-v; Sv=-u;
    R = reflection(S,v);        % v=Su;
else
    R = -eye(length(ref));
end