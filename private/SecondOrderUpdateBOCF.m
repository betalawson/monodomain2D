function [I_ion, S, Sinf, invtau, b] = SecondOrderUpdateBOCF(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, params)
% This function performs a Rush Larsen timestep of the specified length for
% the BuenoOrovio-Cherry-Fenton (4V) model. Current parameters correspond
% to the epicardial cell parameters given in the 4V model paper.
%
% The user may specify additional parameters through the 'params' input
% (currently unused)
%
% b_old is not used here because the BOCF model contains no state variables
% that are not gating variables. 

% Define tissue capacitance (in micro-Farads per square centimetre)
C_m = 1;

% Define voltage scaling
V_rest = -84;
V_max = 1.7;

% Define reference voltages (all scaled to [0,1])
u_v = 0.3;              % Threshold for activation
u_vminus = 0.006;       % Threshold for altered recovery
u_w = 0.13;            % Slow inward activation half point
u_wminus = 0.03;
u_o = 0.006;            % Threshold for slow inward activation gate activation
u_s = 0.9087;
u_so = 0.65;
u_u = 1.55;

% Define time constants (all in ms)
tau_vplus = 1.4506;     % Inactivation of fast inward current
tau_vminus1 = 60;       % Recovery of excitability - controls restitution slope
tau_vminus2 = 1150;     % Recovery of excitability - controls refractory period
tau_wplus = 280;        % Inactivation of slow inward current
tau_wminus1 = 60;       % Recovery of slow inward current
tau_wminus2 = 15;       % Recovery of slow inward current
tau_w_inf = 0.07;
tau_o1 = 400;
tau_o2 = 6;
tau_s1 = 2.7342;
tau_s2 = 16;

% Current time constants (inverse of conductances)
tau_fi = C_m * 0.11;    % Fast inward current
tau_si = 1.8875;        % Slow inward current
tau_so1 = 30.0181;      % Slow outward current (above relevant threshold)
tau_so2 = 0.9957;       % Slow outward current (above relevant threshold)

% Other parameters
k_s = 2.0994;           % Slow inward current activation slope
k_wminus = 65;          % Steady state w gate opening slope
k_so = 2.0458;
w_star_inf = 0.94;      % Steady state value for w above u_o


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All variables are of gating type in the BuenoOrovio-Cherry-Fenton model
gating = logical([1, 1, 1]);
%                 v  w  s

% This also means that the rate of change info for non-gating variables is
% blank
b = [];


%%% Read out variables from state variable matrix
v = S(:,1);
w = S(:,2);
s = S(:,3);


%%% Calculate dimensionless membrane potential
u = (V - V_rest) / (V_max - V_rest);


%%% Determine the rate constants and steady state values for the state
%%% variables

% Create a switch variable that informs which cells are above threshold
u_uv = (u >= u_v);
u_uvminus = (u >= u_vminus);
u_uw = (u >= u_w);
u_uo = (u >= u_o);

% Steady state values
v_inf = ~u_uvminus;
w_inf = ( u_uo .* w_star_inf + ~u_uo .* (1 - u/tau_w_inf) ) .* ~u_uw;
s_inf = ( 1 + tanh(k_s * (u - u_s) ) ) / 2;
Sinf = [v_inf, w_inf, s_inf];

% Time constants
tau_vminus = u_uvminus * tau_vminus2 + ~u_uvminus * tau_vminus1;
tau_wminus = tau_wminus1 + (tau_wminus2 - tau_wminus1) * ( 1 + tanh( k_wminus * (u - u_wminus) ) ) / 2;

tau_v = u_uv * tau_vplus + ~u_uv .* tau_vminus;
tau_w = u_uw * tau_wplus + ~u_uw .* tau_wminus;
tau_s = u_uw * tau_s2 + ~u_uw * tau_s1;
tau_o = u_uo * tau_o2 + ~u_uo * tau_o1;
tau_so = tau_so1 + (tau_so2 - tau_so1) * ( 1 + tanh( k_so * (u - u_so) ) ) / 2;
invtau = 1./[ tau_v, tau_w, tau_s ];

%%% Calculate current strengths 

% Fast inward current
J_fi = -u_uv .* v .* (u - u_v) .* (u_u - u) / tau_fi;
% Slow outward current
J_so = u .* ~u_uw ./ tau_o + u_uw ./ tau_so;
% Slow inward current
J_si = -u_uw .* w .* s / tau_si;

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