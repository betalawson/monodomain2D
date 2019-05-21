function V = processDiffusion(V_old,A,B,b,exact)
% This function simply solves the linear system defined by the input matrix
% and vector, thus processing a diffusive update using implicit
% timestepping.
% 
% The inputs are the matrices A and B, and vector b such that:
%     A V_new = B V_old - b
%
% The user specifies via the 'exact' flag whether to use a direct solution
% method (MATLAB's backslash operator) or a sparse conjugate gradient
% solution technique

% Calculate RHS vector
b_vec = B * V_old - b;

% Solve system using the requested method
if exact
    V = A \ b_vec;
else
    V = conjugateGradient(A, b_vec, V_old, 1e-9);
end