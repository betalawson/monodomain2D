function [V, S] = initialiseCRN(N_cells)
% This function initialises V and S to the starting values for the CRN
% model, creating vectors of length specified by N_cells

% Initial values
V = -81.18;
S = [ 11.17, 139, 1.013e-4, 0.002908, 0.9649, 0.9775, 3.043e-2, 0.9992, 4.966e-3, 0.9986, 3.296e-5, 1.869e-2, 1.367e-4, 0.9996, 0.7755, 2.35e-112, 1, 0.9992, 1.488, 1.488];

% Build matrix of these values repeated the requested number of times
V = repmat(V, N_cells, 1);
S = repmat(S, N_cells, 1);

end