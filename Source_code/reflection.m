function out_vec = reflection(in_vec,n)

% reflection from a plane that has normal n. (any dimension will do)
% householder transformation
out_vec = in_vec - 2*n*(n'*in_vec)/(n'*n);
