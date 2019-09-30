function [I_ion, S, I_Na, I_CaL, I_Kr, I_Ks, I_to, I_K1, I_NaK, I_NaCa] = RLUpdateTT3endo(V, S, dt, I_stim)
% This function performs a Rush Larsen timestep of the specified length for
% the reduced Ten-Tusscher 2006 model for ventricular myocytes. The
% implementation is according to the online source code, except with the
% changes made to reflect the different cell types implemented.


% Define basic constants
R = 8314.472;         % Gas constant (J K^-1 mol^-1 LOL NOT REALLY)
F = 96485.3415;       % Faraday constant (C/mol  BUT CELLML SAYS MILLIMOLE)
T = 310;              % Temperature (K)


% Define basic biophysical properties
Na_o = 140;           % Extracellular Na+ concentration
K_o = 5.4;            % Extracellular K+ concentration
Ca_o = 2;             % Extracellular Ca2+ concentration
Na_i = 7.67;          % Intracellular Na+ concentration (fixed in reduced model)
K_i = 138.3;          % Intracellular K+ concentration (fixed in reduced model)
Ca_i = 0.00007;       % Intracellular Ca2+ concentration (fixed in reduced model)


% Define channel conductances and other current strengths
g_Na = 14.838;        % Maximum conductance of I_Na (nS/pF)
g_to = 0.073;         % Maximum conductance of I_to (nS/pF)
g_Kr = 0.101;         % Maximum conductance of I_Kr (nS/pF)
g_Ks = 0.257;         % Maximum conductance of I_Ks (nS/pF)
g_K1 = 5.405;         % Maximum conductance of I_K1 (nS/pF)
g_CaL = 0.2786;       % Maximum conductance of I_CaL (cm^3 ?F^-1 s^-1)        
k_NaK = 2.724;        % Maximum I_NaK (pA/pF)
k_NaCa = 1000;        % Maximum I_NaCa (pA/pF)

g_pK = 0.0293;        % Maximum conductance of I_pK (nS/pF)
g_pCa = 0.1238;       % Maximum conductance of I_pCa (nS/pF)
g_bNa = 0.00029;      % Maximum conductance of I_bNa (nS/pF)
g_bCa = 0.000592;     % Maximum conductance of I_bCa (nS/pF)


% Define other specific biophysical properties
p_KNa = 0.03;         % Relative permeability to Na+ for I_Ks
gamma = 0.35;         % Voltage dependence parameter for I_NaCa
K_pCa = 0.0005;       % Half-saturation constant for I_pCa - Ca2+ conc. (mM)
K_mK = 1;             % Half-saturation constant for I_NaK - K+ conc. (mM)
K_mNa_NaK = 40;       % Half-saturation constant for I_NaK - Na+ conc. (mM)
K_mCa = 1.38;         % Half-saturation constant for I_NaCa - Ca2+ conc. (mM)
K_mNa_NaCa = 87.5;    % Half-saturation constant for I_NaCa - Na+ conc. (mM)
k_sat = 0.1;          % Saturation factor for I_NaCa
alpha = 2.5;          % Enhancement factor for outward I_NaCa


% Calculate useful basic quantities
RTonF = R * T / F;
FonRT = 1 / RTonF;
expFVonRT = exp(FonRT * V);
K_ofactor = sqrt(K_o / 5.4);


% Read out variables from state variable matrix
m = S(:,1);
h = S(:,2);
j = S(:,3);
s = S(:,4);
xr1 = S(:,5);
xs = S(:,6);
f = S(:,7);
f2 = S(:,8);


% Calculate reversal potentials
E_Na = RTonF * log( Na_o ./ Na_i );
E_K = RTonF * log( K_o ./ K_i );
E_Ks = RTonF * log( (K_o + p_KNa * Na_o) ./ (K_i + p_KNa * Na_i) );
E_Ca = 0.5 * RTonF * log( Ca_o ./ Ca_i );
% Convert these to potential differences
dV_Na = V - E_Na;
dV_K = V - E_K;
dV_Ca = V - E_Ca;


%%% Calculate steady state values and rate constants for the gating variables

% Fast Na+ gates
V_m40 = (V >= -40);
m_inf = 1 ./ ( 1 + exp( (-56.86 - V ) / 9.03 ) ).^2;
tau_m = 1 ./ (1 + exp( -(V+60) / 5 ) ) .* 0.1 .* ( 1 ./ (1 + exp( (V+35)/5 ) ) + 1 ./ ( 1 + exp( (V-50)/200 ) ) );
h_inf = 1 ./ ( 1 + exp( (V + 71.55) / 7.43 ) ).^2;
tau_h = 1./ (  ~V_m40 .* ( 0.057 * exp( -(V+80) / 6.8 ) + 2.7 * exp(0.079 * V) + 310000 * exp(0.3485 * V) ) + V_m40 .* 0.77 ./ (0.13 * (1 + exp(-(V+10.66)/11.1) ) ) );
j_inf = h_inf;
tau_j = 1 ./ ( ~V_m40 .* ( ( -25428 * exp(0.2444 * V) - 0.000006948 * exp(-0.04391 * V) ) .* (V + 37.78) ./ (1 + exp( 0.311 * (V + 79.23) ) )  +  0.02424 * exp(-0.01052 * V) ./ ( 1 + exp(-0.1378 * (V + 40.14) ) ) ) + V_m40 .* 0.6 .* exp(0.057 * V) ./ (1 + exp(-0.1 * (V + 32) ) ) );

% Transient outward K+ gates
r_inf = 1 ./ ( 1 + exp( (20-V) / 6 ) );
s_inf = 1 ./ ( 1 + exp( (V + 28) / 5 ) );
tau_s = 1000 * exp( -(V+67).^2 / 1000 ) + 8;

% Delayed rectifier K+ gates
xr1_inf = 1 ./ ( 1 + exp( -( V + 26) / 7 ) );
tau_xr1 = 450 ./ ( 1 + exp( -(V + 45) / 10 ) ) .* 6 ./ ( 1 + exp( (V+30) / 11.5 ) );
xr2_inf = 1 ./ ( 1 + exp( (V + 88) / 24 ) );
xs_inf = 1 ./ ( 1 + exp( -(V+5)/14 ) );
tau_xs = 1400 ./ sqrt( 1 + exp( (5 - V)/6 ) ) ./ ( 1 + exp( (V-35)/15 ) ) + 80;

% Inward rectifier K+ gate (assumed always at steady state)
alpha_K1 = 0.1 ./ ( 1 + exp( 0.06 * (dV_K - 200 ) ) );
beta_K1 = ( 3 * exp( 0.0002 * (dV_K + 100 ) ) + exp( 0.1 * (dV_K - 10) ) ) ./ ( 1 + exp( -0.5 * dV_K ) );
xK1_inf = alpha_K1 ./ ( alpha_K1 + beta_K1 );

% L-type Ca2+ gates
d_inf = 1 ./ ( 1 + exp( -(V+8)/7.5 ) );       % Re-using calculations to evaluated exp( -(5 + V)/7.5 )
f_inf = 1 ./ ( 1 + exp( (V+20)/7 ) );
tau_f = 1102.5 * exp( -( (V + 27)/15).^2 ) + 200 ./ ( 1 + exp( (13 - V) / 10 ) ) + 180 ./ ( 1 + exp( (V+30)/10 ) ) + 20;
f2_inf = 0.67 ./ (1 + exp( (V + 35)/7 ) ) + 0.33;
tau_f2 = 600 * exp( -(V+27).^2 / 170 ) + 7.75 ./ ( 1 + exp( (25-V)/10 ) ) + 16 ./ ( 1 + exp( (V+30)/10 ) );



%%% Calculate strengths of all currents for the current (V,S) state (post gating updates)

%%% Na+ currents
% Fast inward Na+ current
I_Na = g_Na * m .* m .* m .* h .* j .* dV_Na;
% Background Na+ current
I_bNa = g_bNa * dV_Na;

%%% K+ currents
% Transient outward K+ current
I_to = g_to * r_inf .* s .* dV_K;
% Delayed rectifier K+ currents (rapid and slow)
I_Kr = g_Kr * K_ofactor * xr1 .* xr2_inf .* dV_K;
I_Ks = g_Ks .* xs .* xs .* (V - E_Ks);
% Inward rectifier K+ current
I_K1 = g_K1 * K_ofactor * xK1_inf .* dV_K;
% K+ pump
I_pK = g_pK * dV_K ./ ( 1 + exp((25 - V)/5.98) );

%%% Ca2+ currents
% L-type Ca2+ current
I_CaL = g_CaL * d_inf .* f .* f2 .* (V - 60);
% Ca2+ pump
I_pCa = g_pCa * Ca_i ./ ( K_pCa + Ca_i );
% Background Ca2+ current
I_bCa = g_bCa * dV_Ca;

%%% Ion transfer pumps/exchangers
% Na+/K+ pump
I_NaK = k_NaK * K_o * Na_i ./ ( ( K_o + K_mK ) * (Na_i + K_mNa_NaK) .* (1 + 0.1245 * expFVonRT.^(-0.1) + 0.0353 ./ expFVonRT ) );
% Na+/Ca2+ exchanger
expFVonRTgamma_minus_one = expFVonRT .^ (gamma-1);
I_NaCa = k_NaCa * ( expFVonRTgamma_minus_one .* expFVonRT .* Na_i .* Na_i .* Na_i * Ca_o - expFVonRTgamma_minus_one * Na_o^3 .* Ca_i * alpha ) ./ ( ( K_mNa_NaCa^3 + Na_o^3 ) .* ( K_mCa + Ca_o ) .* (1 + k_sat * expFVonRTgamma_minus_one) );


%%% Total ion transfer in/out of cell membrane
I_ion = I_Na + I_bNa + I_CaL + I_pCa + I_bCa + I_to + I_Kr + I_Ks + I_K1 + I_pK + I_NaK + I_NaCa;


%%% Update gating variables using their steady state values and rate
%%% constants (exponential integration as per Rush Larsen)
m = m_inf + (m - m_inf) .* exp( -dt ./ tau_m );
h = h_inf + (h - h_inf) .* exp( -dt ./ tau_h );
j = j_inf + (j - j_inf) .* exp( -dt ./ tau_j );
s = s_inf + (s - s_inf) .* exp( -dt ./ tau_s );
xr1 = xr1_inf + (xr1 - xr1_inf) .* exp( -dt ./ tau_xr1 );
xs = xs_inf + (xs - xs_inf) .* exp( -dt ./ tau_xs );
f = f_inf + (f - f_inf) .* exp( -dt ./ tau_f );
f2 = f2_inf + (f2 - f2_inf) .* exp( -dt ./ tau_f2 );


%%% Repack all state variables back into the state matrix
S = [m, h, j, s, xr1, xs, f, f2];

end

