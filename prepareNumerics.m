function [A_new, A_old, A_J] = prepareNumerics(K, M, dt, second_order)
% This matrix takes an input stiffness and mass matrix as well as the 
% desired timestep (assumed fixed) and creates instead a problem expressed 
% as:
%
%     A_new * V_new = A_old * V_old - A_J * J
%
% that represents the diffusive updates to the system as calculated by
% replacement of the time derivative with a Crank-Nicholson timestepping
% scheme. No similar operation is required for the reaction terms because 
% the mass matrix applies equally to the time derivative and the reaction 
% update. This code is used simply to avoid re-calculation of unchanging 
% matrices.


if second_order  % Crank-Nicholson for the diffusive portion
    A_new = M - dt/2 * K;
    A_old = M + dt/2 * K;

else  % Fully implicit for the diffusive portion
    A_new = M - dt * K;
    A_old = M;
end

% Reactive term always just uses mass matrix for non-local integration
A_J = dt * M;
