function  V = takeTimestep(V_old, J, A_new, A_old, A_J, exact)
% This function takes as inputs the current 'state' (voltage and current
% vectors) and matrices that define the numerical representation of the
% problem such that:
%
%   A_new * V_new =   A_old * V_old  -  A_J * J .
%
% The user can also specify if the resulting linear system is to be solved
% exactly or by use of the conjugate gradient method.
% This setup assumes a fixed timestep, which is then 'baked in' to the
% forms of A_new, A_old and A_J and thus does not appear explicitly here.
% 
% The matrices in question are setup to only process the active nodes, so
% the inputs V_old and J should be lists of these values at active nodes,
% and the output V will be only at the active nodes.

% Calculate the RHS vector
b_vec = A_old * V_old - A_J * J;

% Solve the linear system
if exact
    V = A_new \ b_vec;
else
    V = conjugateGradient(A_new, b_vec, V_old, 1e-9);
end


end

