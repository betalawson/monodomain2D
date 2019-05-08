function [V, S] = initialiseTT06endo(N_cells)
% This function initialises V and S to the starting values for the TT04
% model (epicardial), creating vectors of length specified by N_cells

% Initial values
%V = -85.23;
%S = [ 8.604, 136.89, 0.000126, 0.00036, 3.64, 0.00172, 0.7444, 0.7045,
%2.42e-8, 0.999998, 0.00621, 0.4712, 0.0095, 3.373e-5, 0.7888, 0.9755,0.9953, 0.9073];     % CellML version
V = -86.2;
S = [ 7.67, 138.3, 0.00007, 0.00007, 1.3, 0, 0.75, 0.75, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1];     % Online code version

% Build matrix of these values repeated the requested number of times
V = repmat(V, N_cells, 1);
S = repmat(S, N_cells, 1);

end

