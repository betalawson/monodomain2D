function [V, S] = initialiseTT04epi(N_cells)
% This function initialises V and S to the starting values for the TT04
% model (epicardial), creating vectors of length specified by N_cells

% Initial values
V = -86.2;
S = [ 11.6, 138.3, 0.0002, 0.2, 0, 0.75, 0.75, 0, 1, 0, 1, 0, 0, 1, 1, 1];

% Build matrix of these values repeated the requested number of times
V = repmat(V, N_cells, 1);
S = repmat(S, N_cells, 1);

end