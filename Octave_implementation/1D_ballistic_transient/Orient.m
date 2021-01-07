function R = orient(ref, dest)

% This code generates the rotation matrix to rotate a vector ref to vector dest
% This is used many times to orient the velocity vector appropriately when
% emitted or diffusively reflected from a plane
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
