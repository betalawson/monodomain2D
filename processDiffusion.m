function V = processDiffusion(V,A,b, dt, scale)
% This function simply solves the linear system defined by the input matrix
% and vector, thus processing a diffusive update using implicit
% timestepping.
% 
% The inputs are the matrix A and vector b such that
%
% (dV/dt)_diffusive = A V - b,
% 
% so that the implicit timestepping solution is obtained by solving
%
%  ( I - dt * A ) * V = V_old - dt * b

% Solve the system to update voltage values
%V = lsqminnorm( eye(size(A)) - dt * scale * A, V - dt * b );
V = (speye(size(A)) - dt * scale * A) \ (V - dt * b);


end

