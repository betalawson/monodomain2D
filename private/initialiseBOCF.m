function [V, S] = initialiseBOCF(N_cells)
% This function initialises V and S to the starting values for the Fenton
% Karma model, creating vectors of length specified by N_cells

% Initial values
V = -84;
S = [1, 1, 0];

% Build matrix of these values repeated the requested number of times
V = repmat(V, N_cells, 1);
S = repmat(S, N_cells, 1);

end

