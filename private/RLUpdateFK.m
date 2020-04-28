function [I_ion, S] = RLUpdateFK(V, S, dt, I_stim, params)
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


%%% Read out variables from state variable matrix
v = S(:,1);
w = S(:,2);


%%% Calculate dimensionless membrane potential
u = (V - V_rest) / (V_max - V_rest);


%%% Calculate current strengths 

% Create a switch variable that informs which cells are above threshold
u_uc = (u >= u_c);
u_uv = (u >= u_v);

% Fast inward current
J_fi = - u_uc .* v / tau_d .* (1 - u) .* (u - u_c);
% Slow outward current
J_so = u_uc / tau_r  + ~u_uc .* (u / tau_o);
% Slow inward current
J_si = -w / (2 * tau_si) .* ( 1 + tanh(k*(u-u_c_si)) );

% Calculate ionic flux (includes an un-scaling)
I_ion = (J_fi + J_so + J_si) * (V_max - V_rest);

%%% Update the state variables using Rush-Larsen update
tau_vminus = u_uv * tau_vminus1 + ~u_uv * tau_vminus2;
v = u_uc .* v * exp(-dt/tau_vplus) + ~u_uc .* ( 1 - (1-v).*exp(-dt./tau_vminus) );
w = u_uc .* w * exp(-dt/tau_wplus) + ~u_uc .* ( 1 - (1-w)*exp(-dt/tau_wminus) );


%%% Repack all state variables back into the state matrix
S = [v, w];


end