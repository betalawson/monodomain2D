function [I_ion, S, Sinf, invtau, b] = SecondOrderUpdateFK(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, params)
% This function performs a Rush Larsen timestep of the specified length for
% the Fenton-Karma (3V) model. Current parameters are Beeler-Reuter 
% parameters (Fenton and Karma, 1998).
%
% The user may specify additional parameters through the 'params' input
% (currently unused)

% Define tissue capacitance (in micro-Farads per square centimetre)
C_m = 1;

% Define voltage scaling
V_rest = -85;
V_max = 15;

% Define reference voltages (all scaled to [0,1])
u_c = 0.13;             % Threshold for activation
u_v = 0.04;            % Threshold for altered recovery
u_c_si = 0.85;          % Slow inward activation half point

% Define time constants (all in ms)
tau_vplus = 3.33;       % Inactivation of fast inward current
tau_vminus1 = 1250;     % Recovery of excitability - controls refractory period
tau_vminus2 = 19.6;     % Recovery of excitability - controls restitution slope
tau_wplus = 870;        % Inactivation of slow inward current
tau_wminus = 41;        % Recovery of slow inward current

% Current time constants (inverse of conductances)
tau_d = C_m/4;          % Fast inward current
tau_r = 33;             % Outward current (during repolarisation)
tau_o = 12.5;            % Outward current (return to rest)
tau_si = 30;            % Slow inward current

% Other parameters
k = 10;                 % Slow inward current activation slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All variables are of gating type in the Fenton-Karma model
gating = logical( [1 1] );
%                  v w

% This also means that the rate of change info for non-gating variables is
% blank
b = [];


%%% Read out variables from state variable matrix
v = S(:,1);
w = S(:,2);


%%% Calculate dimensionless membrane potential
u = (V - V_rest) / (V_max - V_rest);


%%% Steady state values and time constants

% Create a switch variable that informs which cells are above threshold
u_uc = (u >= u_c);
u_uv = (u >= u_v);

% Define inverse time constants and steady state values in one matrix each
tau_vminus = u_uv * tau_vminus1 + ~u_uv * tau_vminus2;
invtau = 1./[ u_uc .* tau_vplus + ~u_uc .* tau_vminus, u_uc .* tau_wplus + ~u_uc .* tau_wminus];
Sinf = [ ~u_uc , ~u_uc ];

%%% Calculate current strengths 
% Fast inward current
J_fi = - u_uc .* v / tau_d .* (1 - u) .* (u - u_c); 
% Slow outward current
J_so = u_uc / tau_r  + ~u_uc .* (u / tau_o);
% Slow inward current
J_si = -w / (2 * tau_si) .* ( 1 + tanh(k*(u-u_c_si)) );

% Calculate ionic flux (includes an un-scaling)
I_ion = (J_fi + J_so + J_si) * (V_max - V_rest);


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