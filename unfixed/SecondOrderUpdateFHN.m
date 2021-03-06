function [I_ion, S, Sinf, invtau, b] = SecondOrderUpdateFHN(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, params)
% This function performs a second order timestep of the specified length 
% for the Fitzhugh-Nagumo modelm as defined on page:
%  http://www.scholarpedia.org/article/Models_of_cardiac_cell
%
% The form of the equations is
%
% du/dt = (u_thresh - u) * (u - 1) * u - v
% dv/dt = epsilon * ( beta * u - gamma * v - delta )
%
% The user may specify additional parameters through the 'params' input
% (currently unused)

% Define voltage scaling
V_rest = -85;
V_max = 15;
V_thresh = -40;

% Define parameters
epsilon = 0.01;
beta = 0.5;
gamma = 1;
delta = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All variables are of gating type in the Fitzhugh-Nagumo model
gating = logical( [1] );
%                  v

% This also means that the rate of change info for non-gating variables is
% blank
b = [];


%%% Read out variables from state variable matrix
v = S(:,1);


%%% Calculate dimensionless membrane potential
u = (V - V_rest) / (V_max - V_rest);
u_thresh = (V_thresh - V_rest) / (V_max - V_rest);


%%% Change in u, scaled back into the corresponding amount of V, is current
I_ion = -( (u_thresh - u) .* (u - 1) .* u - v) * (V_max - V_rest);


%%% Update the state variable (Rush-Larsen holds V fixed, so this is a
%%% linear equation)
Sinf = (beta * u - delta) / gamma;
invtau = gamma * epsilon;

% If dummy information was provided for invtau and Sinf (on first timestep,
% this is not available), just populate them with the current values
if isempty(Sinf_old)
    Sinf_old = Sinf;
end
if isempty(invtau_old)
    invtau_old = invtau;
end

% Create a matrix "A" that stores the approximations of the linear
% coefficients (a in the above) for each variable at the half step
%  a(n+1/2) = 3/2 a(n) - 1/2 a(n-1),    with   a = diag( -1/tau )
% So here, each column of a corresponds to a different variable, and each
% row to another node.
%%% ONLY USED (NONZERO) FOR GATING VARIABLES
A = -3/2 * invtau + 1/2 * invtau_old;
% Ensure no estimated time constants go negative - estimate below zero
% is bad so ensure in this scenario it is calculated using only the current
% estimate (destroys 2nd order to rescue dire situations)
A(A > 0) = -invtau(A > 0);

% Do the same for "B" which is the remainder (here just the constant terms)
%  b(n+1/2) = 3/2 b(n) - 1/2 b(n-1),    with   b = diag( Sinf / tau )
B(:,gating) = 3/2 * (Sinf .* invtau) - 1/2 * (Sinf_old .* invtau_old);

%%% Update gating variables using a second order exponential integration:
%%%  g_new = g_old + dt * ( exp( dt a ) - 1 ) / ( dt a ) * ( a g_old +  b )  - update for gating variables
%%%  S_new = S_old + b                                                       - update for non-gating variables
S(:,gating) = S(:,gating) + dt * ( exp( dt * A ) - 1 ) ./ ( dt * A ) .* ( A .* S(:,gating) + B(:,gating) );

end