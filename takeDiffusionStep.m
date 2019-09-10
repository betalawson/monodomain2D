function  V = takeDiffusionStep(V_old, A_new, A_old, exact)
% This function takes as inputs the current 'state' of voltage values
% and the matrices that define the diffusive portion of the timestep in
% form:
%
%   A_new * V_new =   A_old * V_old 
%
% The matrices in question are setup to only process the active nodes, so
% the inputs V_old and J should be lists of these values at active nodes,
% and the output V will be only at the active nodes.
%
% The timestep is not provided here, because it is 'baked in' to the
% definitions of A_new and A_old (these matrices can be set up using dt/2
% in the case of a Strang splitting)

% Calculate the RHS vector
b_vec = A_old * V_old;

% Solve the linear system
if exact
    V = A_new \ b_vec;
else
    V = conjugateGradient(A_new, b_vec, V_old, 1e-9);
end


end

