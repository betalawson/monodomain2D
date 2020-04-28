function [V, S] = initialiseFHN(N_cells)
% This function initialises V and S to the starting values for the
% FitzHugh-Nagumo model

% Initial values
V = -85;
S = [0];

% Build matrix of these values repeated the requested number of times
V = repmat(V, N_cells, 1);
S = repmat(S, N_cells, 1);

end

