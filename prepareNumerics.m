function [A_new, A_old, A_J] = prepareNumerics(K, M, dt, timestepping, nonlocal_timederiv, nonlocal_reaction)
% This matrix takes an input stiffness and mass matrix as well as the 
% desired timestep (assumed fixed) and creates instead a problem expressed 
% as:
%
%     A_new * V_new = A_old * V_old - A_J * J .
%
% This simply prevents re-calculation of unchanging matrices. A set of
% additional options can be provided to decide which integrations should be
% performed 'non-locally' (involving the mass matrix). The timestepping
% method can also be selected, but note that this only applies to the
% effects of the stiffness matrix (reaction term is always solved using
% fully explicit + Rush-Larsen)

% Non-local integration of the time derivative turn affects the effective
% mass matrix on this term
if nonlocal_timederiv
    M_time = M;
else
    M_time = speye(size(M));
end

% Handling of matrices A_new and A_old depends on the selected timestepping
% method
switch timestepping
    case {'implicit','Implicit'}
        A_new = M_time - dt * K;
        A_old = M_time;
    case {'crank-nicholson', 'Crank-Nicholson', 'Crank-nicholson'}
        A_new = M_time - dt/2 * K;
        A_old = M_time + dt/2 * K;
    otherwise
        error('Timestepping must be ''implicit'' or ''crank-nicholson''');
end 

% Reaction term unaffected by choice of timestepping
if nonlocal_reaction
    A_J = dt * M;
else
    A_J = dt * speye(size(M));
end


end