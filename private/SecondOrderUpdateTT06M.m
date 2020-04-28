function [I_ion, S, Sinf, invtau, b, I_Na, I_CaL, I_Kr, I_Ks, I_to, I_K1, I_NaK, I_NaCa] = SecondOrderUpdateTT06M(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, params)
% This function performs a Rush Larsen timestep of the specified length for
% the Ten-Tusscher 2006 model for ventricular myocytes. The "axial current"
% has not been included for now. The implementation of the f2 gate is the
% one from the CellML version, which differs somewhat from the published
% version.

% Define which state variables are gating variables
gating = logical([ 0,   0,   0,    0,     0,    1, 1, 1, 1, 1, 1,   1,  1,  1, 1, 1,  1,  1]);
%                 Na_i K_i  Ca_i  Ca_ss  Ca_sr  m  h  j  r  s  xr1 xr2  xs  d  f  f2 fCa  Rb

% Define basic constants
R = 8314.472;         % Gas constant (J K^-1 mol^-1 LOL NOT REALLY)
F = 96485.3415;       % Faraday constant (C/mol  BUT CELLML SAYS MILLIMOLE)
T = 310;              % Temperature (K)


% Define basic biophysical properties
Vol_c = 0.016404;     % Volume of cytoplasm (?m^3)
Vol_sr = 0.001094;    % Volume of sarcoplasmic reticulum (?m^3)
Vol_ss = 0.00005468;  % Volume of subspace (?m^3)
Na_o = 140;           % Extracellular Na+ concentration
K_o = 5.4;            % Extracellular K+ concentration
Ca_o = 2;             % Extracellular Ca2+ concentration
Cm = 0.185;           % Cell capacitance (microfarads)


% Define channel conductances and other current strengths
g_Na = 14.838;        % Maximum conductance of I_Na (nS/pF)
g_to = 0.294;         % Maximum conductance of I_to (nS/pF)
g_Kr = 0.153;         % Maximum conductance of I_Kr (nS/pF)
g_Ks = 0.098;         % Maximum conductance of I_Ks (nS/pF)
g_K1 = 5.405;         % Maximum conductance of I_K1 (nS/pF)
g_CaL = 0.00003980;   % Maximum conductance of I_CaL (cm^3 ?F^-1 s^-1)        
k_NaK = 2.724;        % Maximum I_NaK (pA/pF)
k_NaCa = 1000;        % Maximum I_NaCa (pA/pF)
k_maxup = 0.006375;   % Maximum value for I_up (mM/ms)
k_rel = 0.102;        % Maximum value for I_rel (mM/ms)

g_pK = 0.0146;        % Maximum conductance of I_pK (nS/pF)
g_pCa = 0.1238;       % Maximum conductance of I_pCa (nS/pF)
g_bNa = 0.00029;      % Maximum conductance of I_bNa (nS/pF)
g_bCa = 0.000592;     % Maximum conductance of I_bCa (nS/pF)
k_leak = 0.00036;     % Rate of leak current (ms^-1)
k_xfer = 0.0038;      % Rate of transfer current (ms^-1)


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
buf_c = 0.2;          % Cytoplasmic Ca2+ buffer concentration (mM)
K_bufc = 0.001;       % Half saturation constant for cytoplasmic buffer - Ca2+ conc. (mM)
buf_ss = 0.4;         % Subspace Ca2+ buffer concentration (mM)
K_bufss = 0.00025;    % Half saturation constant for subspace buffer - Ca2+ conc. (mM)
buf_sr = 10;          % Sarcoplasmic reticulum Ca2+ buffer concentration (mM)
K_bufsr = 0.3;        % Half saturation constant for sarcoplasmic reticulum buffer - Ca2+ conc. (mM)
k_1base = 0.15;       % CICR-related constant 1
k_2base = 0.045;      % CICR-related constant 2
k_3 = 0.06;           % CICR-related constant 3
k_4 = 0.005;          % CICR-related constant 4
EC = 1.5;             % Half saturation constant for k_casr - Ca2+ SR conc. (mM)
k_casrmax = 2.5;      % Maximum value of k_sr
k_casrmin = 1;        % Minimum value of k_sr


% Calculate useful basic quantities
RTonF = R * T / F;
FonRT = 1 / RTonF;
expFVonRT = exp(FonRT * V);
K_ofactor = sqrt(K_o / 5.4);


% Read out variables from state variable matrix
Na_i = S(:,1);
K_i = S(:,2);
Ca_i = S(:,3);
Ca_ss = S(:,4);
Ca_sr = S(:,5);
m = S(:,6);
h = S(:,7);
j = S(:,8);
r = S(:,9);
s = S(:,10);
xr1 = S(:,11);
xr2 = S(:,12);
xs = S(:,13);
d = S(:,14);
f = S(:,15);
f2 = S(:,16);
fCass = S(:,17);
Rb = S(:,18);


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
tau_r = 9.5 * exp( -( (V+40).^2 / 1800 ) ) + 0.8;
s_inf = 1 ./ ( 1 + exp( (V + 20) / 5 ) );
tau_s = 85 * exp( -( (V+45).^2 / 320) ) + 5 ./ ( 1 + exp( (V - 20) / 5 ) ) + 3;

% Delayed rectifier K+ gates
xr1_inf = 1 ./ ( 1 + exp( -( V + 26) / 7 ) );
tau_xr1 = 450 ./ ( 1 + exp( -(V + 45) / 10 ) ) .* 6 ./ ( 1 + exp( (V+30) / 11.5 ) );
xr2_inf = 1 ./ ( 1 + exp( (V + 88) / 24 ) );
tau_xr2 = 3 ./ ( 1 + exp( -(V+60) / 20 ) ) .* 1.12 ./ ( 1 + exp( (V-60)/20 ) );
xs_inf = 1 ./ ( 1 + exp( -(V+5)/14 ) );
tau_xs = 1400 ./ sqrt( 1 + exp( (5 - V)/6 ) ) ./ ( 1 + exp( (V-35)/15 ) ) + 80;

% Inward rectifier K+ gate (assumed always at steady state)
alpha_K1 = 0.1 ./ ( 1 + exp( 0.06 * (dV_K - 200 ) ) );
beta_K1 = ( 3 * exp( 0.0002 * (dV_K + 100 ) ) + exp( 0.1 * (dV_K - 10) ) ) ./ ( 1 + exp( -0.5 * dV_K ) );
xK1_inf = alpha_K1 ./ ( alpha_K1 + beta_K1 );

% L-type Ca2+ gates
d_inf = 1 ./ ( 1 + exp( -(V+8)/7.5 ) );       % Re-using calculations to evaluated exp( -(5 + V)/7.5 )
tau_d = ( 0.25 + 1.4 ./ (1 + exp( -(35 + V)/13 ) ) ) .* ( 1.4 ./ (1 + exp( (V+5)/5) ) ) + 1 ./ ( 1 + exp( (50 - V)/20 ) );
f_inf = 1 ./ ( 1 + exp( (V+20)/7 ) );
tau_f = 1102.5 * exp( -( (V + 27)/15).^2 ) + 200 ./ ( 1 + exp( (13 - V) / 10 ) ) + 180 ./ ( 1 + exp( (V+30)/10 ) ) + 20;
f2_inf = 0.67 ./ (1 + exp( (V + 35)/7 ) ) + 0.33;
%tau_f2 = 600 * exp( -(V+25).^2 / 170 ) + 31 ./ ( 1 + exp( (25-V)/10 ) ) + 16 ./ ( 1 + exp( (V+30)/10 ) );
tau_f2 = 562 * exp( -(V+27).^2 / 240 ) + 31 ./ ( 1 + exp( (25-V)/10 ) ) + 80 ./ ( 1 + exp( (V+30)/10 ) );
fCass_inf = 0.6 ./ ( 1 + (Ca_ss/0.05).^2 ) + 0.4;
tau_fCass = 2 + 80 ./ ( 1 + (Ca_ss/0.05).^2 );

%%% Calculate some state-dependent values for Ca2+ handling
k_casr = k_casrmax - ( k_casrmax - k_casrmin) ./ ( 1 + (EC ./ Ca_sr ).^2 );
k_1 = k_1base ./ k_casr;
k_2 = k_2base * k_casr;
Rb_inf = k_4 ./ (k_4 + k_2 .* Ca_ss); 
tau_Rb = 1 ./ (k_4 + k_2 .* Ca_ss);


% Store all of the inverse time constants and steady state values in one 
% big matrix each
invtau = 1./[tau_m, tau_h, tau_j, tau_r, tau_s, tau_xr1, tau_xr2, tau_xs, tau_d, tau_f, tau_f2, tau_fCass, tau_Rb];
Sinf = [m_inf, h_inf, j_inf, r_inf, s_inf, xr1_inf, xr2_inf, xs_inf, d_inf, f_inf, f2_inf, fCass_inf, Rb_inf];


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
exp2FVm15onRT = expFVonRT .* expFVonRT .* exp(-30 * FonRT);
I_CaL = g_CaL * d .* f .* f2 .* fCass .* 4 * FonRT * F .* (V - 15) .* ( 0.25 * Ca_ss .* exp2FVm15onRT - Ca_o) ./ ( exp2FVm15onRT - 1);
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
I_rel = k_rel * (k_1 .* Ca_ss .* Ca_ss .* Rb) ./ (k_3 + k_1 .* Ca_ss .* Ca_ss) .* (Ca_sr - Ca_ss);
% Ca2+ transfer from subspace to cytoplasm
I_xfer = k_xfer .* (Ca_ss - Ca_i);


%%% Total ion transfer in/out of cell membrane
I_ion = I_Na + I_bNa + I_CaL + I_pCa + I_bCa + I_to + I_Kr + I_Ks + I_K1 + I_pK + I_NaK + I_NaCa;



%%% Updates to ion concentrations in terms of currents

% Intracellular Na+
b(:,1) = ( - ( I_Na + I_bNa + 3 * I_NaK + 3 * I_NaCa ) / (Vol_c * F) * Cm );

% Intracellular K+
b(:,2) = K_i + dt * ( - ( I_to + I_Kr + I_Ks - 2 * I_NaK + I_pK + I_stim ) / (Vol_c * F) * Cm );

% Intracellular Ca2+
b(:,3) = Ca_i + dt ./ ( 1 + buf_c * K_bufc ./ ( ( Ca_i + K_bufc ).^2 ) ) .* ( - (I_bCa + I_pCa - 2 * I_NaCa ) / ( 2 * Vol_c * F ) * Cm + Vol_sr / Vol_c * (I_leak - I_up) + I_xfer );

% Subspace Ca2+
b(:,4) = Ca_ss + dt ./ ( 1 + buf_ss * K_bufss ./ ( ( Ca_ss + K_bufss).^2 ) ) .* ( -I_CaL / (2 *  F) * Cm + Vol_sr * I_rel - Vol_c * I_xfer ) / Vol_ss;

% SR Ca2+
b(:,5) = Ca_sr + dt ./ ( 1 + buf_sr * K_bufsr ./ ( Ca_sr + K_bufsr ).^2 ) .* ( -I_leak + I_up - I_rel );

% If dummy information was provided for invtau and Sinf (on first timestep,
% this is not available), just populate them with the current values
if isempty(Sinf_old)
    Sinf_old = Sinf;
end
if isempty(invtau_old)
    invtau_old = invtau;
end
if isempty(b_old)
    b_old = b;
end


% Create a matrix "A" that stores the approximations of the linear
% coefficients (a in the above) for each variable at the half step
%  a(n+1/2) = 3/2 a(n) - 1/2 a(n-1),    with   a = diag( -1/tau )
% So here, each column of a corresponds to a different variable, and each
% row to another node.
%%% ONLY USED (NONZERO) FOR GATING VARIABLES
A = -3/2 * invtau + 1/2 * invtau_old;
% Ensure no estimated time constants go negative - estimate below zero
% corresponds to a rapidly decreasing value for invtau and so it is set to
% a small positive value
A(A > 0) = -1e-6;

% Do the same for "B" which is the remainder (here just the constant terms)
%  b(n+1/2) = 3/2 b(n) - 1/2 b(n-1),    with   b = diag( Sinf / tau )
B(:,gating) = 3/2 * (Sinf .* invtau) - 1/2 * (Sinf_old .* invtau_old);
B(:,~gating) = 3/2 * b - 1/2 * b_old;

%%% Update gating variables using a second order exponential integration:
%%%  g_new = g_old + dt * ( exp( dt a ) - 1 ) / ( dt a ) * ( a g_old +  b )  - update for gating variables
%%%  S_new = S_old + b                                                       - update for non-gating variables
S(:,gating) = S(:,gating) + dt * ( exp( dt * A ) - 1 ) ./ ( dt * A ) .* ( A .* S(:,gating) + B(:,gating) );
S(:,~gating) = S(:,~gating) + dt * B(:,~gating);

end