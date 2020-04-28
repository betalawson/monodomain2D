function [I_ion, S] = RLUpdateFHN(V, S, dt, I_stim, params)
% This function performs a Rush Larsen timestep of the specified length for
% the Fitzhugh-Nagumo modelm as defined on page:
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
V_max = 40;
V_thresh = -67.5;

% Define parameters
c1 = 1;
epsilon = 0.005;
beta = 0.5;
gamma = 1;
delta = 0;


%%% Read out variables from state variable matrix
v = S(:,1);


%%% Calculate dimensionless membrane potential
u = (V - V_rest) / (V_max - V_rest);
u_thresh = (V_thresh - V_rest) / (V_max - V_rest);


%%% Change in u, scaled back into the corresponding amount of V, is current
I_ion = -( c1 * (u_thresh - u) .* (u - 1) .* u - v) * (V_max - V_rest);


%%% Update the state variable (Rush-Larsen holds V fixed, so this is a
%%% linear equation)
v_inf = (beta * u - delta) / gamma;
tau_v = 1 / (gamma * epsilon);
v = v_inf + (v - v_inf) .* exp( -dt ./ tau_v );

%%% Repack all state variables back into the state matrix
S = [v];


end