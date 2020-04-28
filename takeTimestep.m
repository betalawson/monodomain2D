function  V = takeTimestep(V_old, J, J_old, A_new, A_old, A_J, exact, preconditioning, second_order)
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
%
% The user specifies if second order updates are to be used. These are of
% Adams-Bashforth type.

% If first timestep, just use J_old = J
if isempty(J_old)
    J_old = J;
end

% Calculate the RHS vector
if second_order
    b_vec = A_old * V_old - A_J * (3/2 * J - 1/2 * J_old);   % Adams-Bashforth update for reaction term
else
    b_vec = A_old * V_old - A_J * J;                         % Standard explicit update (with mass matrix)
end

% Solve the linear system
if exact
    V = A_new \ b_vec;
else
    
    % Take a timestep using the biconjugate gradient method, with basic
    % preconditioning if requested
    if preconditioning
        [L,U] = ilu(A_new);
        [V, flag, relres] = bicgstab(A_new, b_vec, 1e-9, 20, L, U, V_old);
    else
        [V, flag, relres] = bicgstab(A_new, b_vec, 1e-9, 20, speye(size(A_new)), speye(size(A_new)), V_old);
    end
    
    if flag ~= 0
        warning('Biconjugate gradient failed to converge to requested tolerance in specified number of maximum iterations. Residual was %g',relres);
    end
    
end


end

