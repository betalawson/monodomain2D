function [I_ion, S, I_Na, I_CaL, I_Kr, I_Ks, I_to, I_K1, I_NaK, I_NaCa, I_up, I_rel] = RLUpdateTT04endo(V, S, dt, I_stim)
% This function performs a Rush Larsen timestep of the specified length for
% the Ten-Tusscher 2004 model for ventricular myocytes. The "axial current"
% has not been included for now.


% Define basic constants
R = 8314.472;         % Gas constant (J K^-1 mol^-1 LOL NOT REALLY)
F = 96485.3415;       % Faraday constant (C/mol  BUT CELLML SAYS MILLIMOLE)
T = 310;              % Temperature (K)


% Define basic biophysical properties
Vol_c = 0.016404;     % Volume of cytoplasm (?m^3)
Vol_sr = 0.001094;    % Volume of sarcoplasmic reticulum (?m^3)
Na_o = 140;           % Extracellular Na+ concentration
K_o = 5.4;            % Extracellular K+ concentration
Ca_o = 2;             % Extracellular Ca2+ concentration
Cm = 0.185;           % Cell capacitance (microfarads)


% Define channel conductances and other current strengths
g_Na = 14.838;        % Maximum conductance of I_Na (nS/pF)
g_to = 0.073;         % Maximum conductance of I_to (nS/pF)
g_Kr = 0.096;         % Maximum conductance of I_Kr (nS/pF)
g_Ks = 0.245;         % Maximum conductance of I_Ks (nS/pF)
g_K1 = 5.405;         % Maximum conductance of I_K1 (nS/pF)
g_CaL = 0.000175;     % Maximum conductance of I_CaL (cm^3 ?F^-1 s^-1)        
k_NaK = 1.362;        % Maximum I_NaK (pA/pF)
k_NaCa = 1000;        % Maximum I_NaCa (pA/pF)
k_maxup = 0.000425;   % Maximum value for I_up (mM/ms)
I_rel_scale = 1;      % Scaling factor that may be applied to I_rel (for variability studies)

g_pK = 0.0146;        % Maximum conductance of I_pK (nS/pF)
g_pCa = 0.825;        % Maximum conductance of I_pCa (nS/pF)
g_bNa = 0.00029;      % Maximum conductance of I_bNa (nS/pF)
g_bCa = 0.000592;     % Maximum conductance of I_bCa (nS/pF)
k_leak = 0.00008;     % Rate of leak current (ms^-1)


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
K_up = 0.00025;       % Half-saturation constant for I_up (mM)
a_rel = 0.016464;     % Maximal Ca_SR dependent I_rel (mM s^-1)
b_rel = 0.25;         % Ca_sr half-saturation constant for I_rel (mM)
c_rel = 0.008232;     % Maximal Ca_sr independent I_rel (mM s^-1)
buf_c = 0.15;         % Cytoplasmic Ca2+ buffer concentration (mM)
K_bufc = 0.001;       % Half saturation constant for cytoplasmic buffer - Ca2+ conc. (mM)
buf_sr = 10;          % Sarcoplasmic reticulum Ca2+ buffer concentration (mM)
K_bufsr = 0.3;        % Half saturation constant for sarcoplasmic reticulum buffer - Ca2+ conc. (mM)


% Calculate useful basic quantities
RTonF = R * T / F;
FonRT = 1 / RTonF;
expFVonRT = exp(FonRT * V);
K_ofactor = sqrt(K_o / 5.4);


% Read out variables from state variable matrix
Na_i = S(:,1);
K_i = S(:,2);
Ca_i = S(:,3);
Ca_sr = S(:,4);
m = S(:,5);
h = S(:,6);
j = S(:,7);
r = S(:,8);
s = S(:,9);
xr1 = S(:,10);
xr2 = S(:,11);
xs = S(:,12);
d = S(:,13);
f = S(:,14);
fCa = S(:,15);
g = S(:,16);


% Calculate reversal potentials
E_Na = RTonF * log( Na_o ./ Na_i );
E_K = RTonF * log( K_o ./ K_i );
E_Ks = RTonF * log( (K_o + p_KNa * Na_o) ./ (K_i + p_KNa * Na_i) );
E_Ca = 0.5 * RTonF * log( Ca_o ./ Ca_i );
% Convert these to potential differences
dV_Na = V - E_Na;
dV_K = V - E_K;
dV_Ca = V - E_Ca;


% Create switch variables that select whether or not variables are above
% certain threshold
V_m40 = (V >= -40);
V_m60 = (V >= -60);
Ca_35em5 = (Ca_i >= 0.00035);

%%% Calculate steady state values and rate constants for the gating variables

% Fast Na+ gates
m_inf = 1 ./ ( 1 + exp( (-56.86 - V ) / 9.03 ) ).^2;
tau_m = 1 ./ (1 + exp( -(V+60) / 5 ) ) .* 0.1 .* ( 1 ./ (1 + exp( (V+35)/5 ) ) + 1 ./ ( 1 + exp( (V-50)/200 ) ) );
h_inf = 1 ./ ( 1 + exp( (V + 71.55) / 7.43 ) ).^2;
tau_h = 1./ (  ~V_m40 .* ( 0.057 * exp( -(V+80) / 6.8 ) + 2.7 * exp(0.079 * V) + 310000 * exp(0.3485 * V) ) + V_m40 .* 0.77 ./ (0.13 * (1 + exp(-(V+10.66)/11.1) ) ) );
j_inf = h_inf;
tau_j = 1 ./ ( ~V_m40 .* ( ( -25428 * exp(0.2444 * V) - 0.000006948 * exp(-0.04391 * V) ) .* (V + 37.78) ./ (1 + exp( 0.311 * (V + 79.23) ) )  +  0.02424 * exp(-0.01052 * V) ./ ( 1 + exp(-0.1378 * (V + 40.14) ) ) ) + V_m40 .* 0.6 .* exp(0.057 * V) ./ (1 + exp(-0.1 * (V + 32) ) ) );

% Transient outward K+ gates
r_inf = 1 ./ ( 1 + exp( (20-V) / 6 ) );
tau_r = 9.5 * exp( -( (V+40).^2 / 1800 ) ) + 0.8;
s_inf = 1 ./ ( 1 + exp( (V + 28) / 5 ) );
tau_s = 1000 * exp( -(V+67).^2 / 1000 ) + 8;

% Delayed rectifier K+ gates
xr1_inf = 1 ./ ( 1 + exp( -( V + 26) / 7 ) );
tau_xr1 = 450 ./ ( 1 + exp( -(V + 45) / 10 ) ) .* 6 ./ ( 1 + exp( (V+30) / 11.5 ) );
xr2_inf = 1 ./ ( 1 + exp( (V + 88) / 24 ) );
tau_xr2 = 3 ./ ( 1 + exp( -(V+60) / 20 ) ) .* 1.12 ./ ( 1 + exp( (V-60)/20 ) );
xs_inf = 1 ./ ( 1 + exp( -(V+5)/14 ) );
tau_xs = 1100 ./ sqrt( 1 + exp( -(V+10)/6 ) ) ./ ( 1 + exp( (V-60)/20 ) );

% Inward rectifier K+ gate (assumed always at steady state)
alpha_K1 = 0.1 ./ ( 1 + exp( 0.06 * (dV_K - 200 ) ) );
beta_K1 = ( 3 * exp( 0.0002 * (dV_K + 100 ) ) + exp( 0.1 * (dV_K - 10) ) ) ./ ( 1 + exp( -0.5 * dV_K ) );
xK1_inf = alpha_K1 ./ ( alpha_K1 + beta_K1 );

% L-type Ca2+ gates
d_inf = 1 ./ ( 1 + exp( -(V+5)/7.5 ) );       % Re-using calculations to evaluated exp( -(5 + V)/7.5 )
tau_d = ( 0.25 + 1.4 ./ (1 + exp( -(35 + V)/13 ) ) ) .* ( 1.4 ./ (1 + exp( (V+5)/5) ) ) + 1 ./ ( 1 + exp( (50 - V)/20 ) );
f_inf = 1 ./ ( 1 + exp( (V+20)/7 ) );
tau_f = 1125 * exp( -( (V + 27).^2 )/240 ) + 165 ./ ( 1 + exp( (25 - V) / 10 ) ) + 80;
fCa_inf = ( 1 ./ ( 1 + ( Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i ) / (0.000325^8) ) +  0.1 ./ ( 1 + exp( (Ca_i - 0.0005) / 0.0001) )  +  0.2 ./ ( 1 + exp( (Ca_i - 0.00075) / 0.0008 ) ) + 0.23 ) / 1.46;
tau_fCa = 2;

% Calcium release gate
g_inf = 1 ./ ( 1 + Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i / (0.00035^6) .* ( Ca_35em5 .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i .* Ca_i / (0.00035^10) + ~Ca_35em5 ) );
tau_g = 2;


%%% Update gating variables using their steady state values and rate
%%% constants (exponential integration as per Rush Larsen)
m = m_inf + (m - m_inf) .* exp( -dt ./ tau_m );
h = h_inf + (h - h_inf) .* exp( -dt ./ tau_h );
j = j_inf + (j - j_inf) .* exp( -dt ./ tau_j );
r = r_inf + (r - r_inf) .* exp( -dt ./ tau_r );
s = s_inf + (s - s_inf) .* exp( -dt ./ tau_s );
xr1 = xr1_inf + (xr1 - xr1_inf) .* exp( -dt ./ tau_xr1 );
xr2 = xr2_inf + (xr2 - xr2_inf) .* exp( -dt ./ tau_xr2 );
xs = xs_inf + (xs - xs_inf) .* exp( -dt ./ tau_xs );
d = d_inf + (d - d_inf) .* exp( -dt ./ tau_d );
f = f_inf + (f - f_inf) .* exp( -dt ./ tau_f );
fCa_switch = ( fCa >= fCa_inf | ~V_m60 );
fCa = fCa_switch .* ( fCa_inf + (fCa - fCa_inf) .* exp( -dt ./ tau_fCa ) ) + ~fCa_switch .* fCa;    % Change turned off if fCa_inf > fCa and V > -60
g_switch = ( g >= g_inf | ~V_m60 );
g = g_switch .* ( g_inf + (g - g_inf) .* exp( -dt ./ tau_g ) ) + ~g_switch .* g;                    % Change turned off if g_inf > g and V > -60


%%% Calculate strengths of all currents for the current (V,S) state (post gating updates)

%%% Na+ currents
% Fast inward Na+ current
I_Na = g_Na * m .* m .* m .* h .* j .* dV_Na;
% Background Na+ current
I_bNa = g_bNa * dV_Na;

%%% K+ currents
% Transient outward K+ current
I_to = g_to * r .* s .* dV_K;
% Delayed rectifier K+ currents (rapid and slow)
I_Kr = g_Kr * K_ofactor * xr1 .* xr2 .* dV_K;
I_Ks = g_Ks .* xs .* xs .* (V - E_Ks);
% Inward rectifier K+ current
I_K1 = g_K1 * K_ofactor * xK1_inf .* dV_K;
% K+ pump
I_pK = g_pK * dV_K ./ ( 1 + exp((25 - V)/5.98) );

%%% Ca2+ currents
% L-type Ca2+ current
I_CaL = g_CaL * d .* f .* fCa .* 4 * FonRT * F .* V .* ( Ca_i .* expFVonRT .* expFVonRT - 0.341 * Ca_o) ./ ( expFVonRT .* expFVonRT - 1);
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

%%% Intracellular currents
% Ca2+ leak current
I_leak = k_leak .* (Ca_sr - Ca_i);
% Ca2+ uptake
I_up = k_maxup ./ ( 1 + K_up.^2 ./ ( Ca_i .* Ca_i ) );
% Ca2+ release
I_rel = I_rel_scale * ( a_rel .* Ca_sr .* Ca_sr ./ (  b_rel.*b_rel + Ca_sr .* Ca_sr ) + c_rel ) .* d .* g;


%%% Total ion transfer in/out of cell membrane
I_ion = I_Na + I_bNa + I_CaL + I_pCa + I_bCa + I_to + I_Kr + I_Ks + I_K1 + I_pK + I_NaK + I_NaCa;



%%% Updates to ion concentrations in terms of currents

% Intracellular Na+
Na_i = Na_i + dt * ( - ( I_Na + I_bNa + 3 * I_NaK + 3 * I_NaCa ) / (Vol_c * F) * Cm );

% Intracellular K+
K_i = K_i + dt * ( - ( I_to + I_Kr + I_Ks - 2 * I_NaK + I_pK + I_stim ) / (Vol_c * F) * Cm );

% Intracellular Ca2+
Ca_i = Ca_i + dt ./ ( 1 + buf_c * K_bufc ./ ( ( Ca_i + K_bufc ).^2 ) ) .* ( - (I_CaL + I_bCa + I_pCa - 2 * I_NaCa ) / ( 2 * Vol_c * F ) * Cm + I_leak - I_up + I_rel );

% SR Ca2+
Ca_sr = Ca_sr + dt ./ ( 1 + buf_sr * K_bufsr ./ ( Ca_sr + K_bufsr ).^2 ) .* Vol_c / Vol_sr .* ( -I_leak + I_up - I_rel );


%%% Repack all state variables back into the state matrix
S = [Na_i, K_i, Ca_i, Ca_sr, m, h, j, r, s, xr1, xr2, xs, d, f, fCa, g];

end

