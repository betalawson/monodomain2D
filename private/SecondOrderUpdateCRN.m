function [I_ion, S, Sinf, invtau, b, I_Na, I_CaL, I_Kr, I_Ks, I_to, I_K1, I_NaK, I_NaCa] = SecondOrderUpdateCRN(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, params)
% This function performs a Rush Larsen timestep of the specified length for
% the Courtemanche-Ramirez-Nattel 1996 model for atrial myocytes.

% Define which variables are gating variables
gating = logical( [  0,    0,   0,   1, 1, 1, 1,  1,  1,  1,  1,  1,  1, 1,  1,  1, 1, 1,   0,     0] );
%                   Na_i  K_i  Ca_i  m  h  j  oa  oi  ua  ui  xr  xs  d  f  fCa  u  v  w  Ca_rel  Ca_up

% Define basic constants
R = 8.3143;            % Gas constant
T = 309.65;            % Temperature (K)
F = 96.4867;           % Faraday constant

% Cell properties
Cm = 100;                       % Membrane capacitance
Vol_tot = 20100;				% Cell volume
Vol_i = 0.68 * Vol_tot;			% Cell volume for intracellular ions
Vol_rel = 0.0048 * Vol_tot;		% Cell volume for 'rel' component
Vol_up = 0.0552 * Vol_tot;		% Cell volume for 'up' component

% Extracellular ion concentrations
Na_o = 140;         	% Extracellular Na concentration
K_o = 5.4;				% Extracellular K concentration	
Ca_o = 1.8;				% Extracellular Ca2+ concentration

% Effect of temperature
KQ10 = 3;				% Q10 value (should depend on temperature)

% Current strengths
g_Na = 7.8;
g_bNa = 0.0006744375;
g_to = 0.1652;
g_Kur_scale = 1;
g_Kr = 0.029411765;
g_Ks = 0.12941176;
g_K1 = 0.09;
g_CaL = 0.12375;
g_bCa = 0.001131;
g_NaK = 0.59933874;
g_NaCa = 1600;
I_CaP_max = 0.275;
I_up_max = 0.005;
J_rel = 30;
Ca_up_max = 15;

% Other constants

% Na+/K+ pump
Km_Nai = 10;
Km_Ko = 1.5;
% Na+/Ca2+ exchanger
gamma = 0.35;
Km_Nao = 87.5;
Km_Cao = 1.38;
k_sat = 0.1;
% Ca2+ handling/buffering
K_up = 0.00092;
CMDN_max = 0.05;
Km_CMDN = 0.00238;
CSQN_max = 10;
Km_CSQN = 0.8;
TRPN_max = 0.07;
Km_TRPN = 0.0005;


% Read out variables from state variable matrix
Na_i = S(:,1);
K_i = S(:,2);
Ca_i = S(:,3);
m = S(:,4);
h = S(:,5);
j = S(:,6);
oa = S(:,7);
oi = S(:,8);
ua = S(:,9);
ui = S(:,10);
xr = S(:,11);
xs = S(:,12);
d = S(:,13);
f = S(:,14);
fCa = S(:,15);
u = S(:,16);
v = S(:,17);
w = S(:,18);
Ca_rel = S(:,19);
Ca_up = S(:,20);


% Calculate the Nernst Potential differentials
RTonF = R * T / F;
VminusENa = V - RTonF * log( Na_o ./ Na_i );
VminusEK = V - RTonF * log( K_o ./ K_i );
VminusECa = V - RTonF / 2 * log( Ca_o ./ Ca_i );
	
% Calculate some useful constants to speed things up a little
expFVonRT = exp(V / RTonF);
expFVonRTtoGammaMinusOne = expFVonRT.^(gamma-1);
sigma = (1/7) * ( exp( Na_o / 67.3 ) - 1);


% Create switch variables that select whether or not variables are above
% certain threshold
V_m40 = (V >= -40);


%%% Calculate the values of the derived currents

%%% Na+ currents
% Fast Na+ current
I_Na = Cm * g_Na * m .* m .* m .* h .* j .* VminusENa;
% Background Na+ current
I_bNa = Cm * g_bNa * VminusENa;

%%% K+ currents
% Transient outward K+ current
I_to = Cm * g_to * oa .* oa .* oa .* oi .* VminusEK;
% Ultrarapid rectifier K+ current
g_Kur = g_Kur_scale * (0.005 + 0.05 ./ ( 1 + exp( -(V - 15)/13) ) );
I_Kur = Cm * g_Kur .* ua .* ua .* ua .* ui .* VminusEK;
% Delayed rectifier K+ current - rapid component
I_Kr = Cm * g_Kr * xr .* VminusEK ./ (1 + exp( (V+15)/22.4) );
% Delayed rectifier K+ current - slow component
I_Ks = Cm * g_Ks * xs .* xs .* VminusEK;
% Inward rectifier K+ current
I_K1 = Cm * g_K1 * VminusEK ./ (1 + exp(0.07 * (V+80) ) );

%%% Ca2+ currents
% L-type Ca2+ current
I_CaL = Cm * g_CaL * d .* f .* fCa .* (V - 65);
% Background Ca2+ current
I_bCa = Cm * g_bCa * VminusECa;

%%% Pump and exchange currents
% Na+/K+ pump
I_NaK = Cm * g_NaK * ( K_o / (K_o + Km_Ko) ) ./ ( (1 + (Km_Nai ./ Na_i).^1.5) .* (1 + 0.1245 * expFVonRT.^-0.1 + 0.0365 * sigma ./ expFVonRT) );
% Na+/Ca2+ exchanger
I_NaCa = Cm * g_NaCa * (expFVonRT .* expFVonRTtoGammaMinusOne .* Na_i .* Na_i .* Na_i * Ca_o - expFVonRTtoGammaMinusOne * Na_o^3 .* Ca_i ) ./ ( (Km_Nao^3 + Na_o^3) * (Km_Cao + Ca_o) * (1 + k_sat * expFVonRTtoGammaMinusOne) );
% Ca2+ pump
I_CaP = Cm * I_CaP_max * Ca_i ./ ( 0.0005 + Ca_i);

%%% Internal currents
I_tr = ( Ca_up - Ca_rel ) / 180;
I_rel = J_rel * u .* u .* v .* w .* ( Ca_rel - Ca_i );
I_up = I_up_max ./ ( 1 + K_up ./ Ca_i );
I_up_leak = I_up_max * Ca_up / Ca_up_max;


% Total membrane current
I_ion = (I_Na + I_bNa + I_to + I_Kur + I_Kr + I_Ks + I_K1 + I_CaL + I_bCa + I_NaK + I_NaCa + I_CaP) / Cm;


%%% Calculate steady state values and rate constants for the gating variables
%%% This is performed after current definition, because some values depend
%%% on currents


%%% Fast Na+ gates
% m gate - activation
Vshift = V + 47.13;
alpha_m = 0.32 * Vshift ./ ( 1 - exp( - 0.1 * Vshift ) );
alpha_m( abs(Vshift) < 1e-9 ) = 3.2;                           % Remove divide by zero using limiting value
beta_m = 0.08 * exp(-V/11);
invtau_m = alpha_m + beta_m;
m_inf = alpha_m ./ invtau_m;
% h gate - initial inactivation
alpha_h = (~V_m40) .* ( 0.135 * exp( -(V+80) / 6.8 ) );
beta_h = (~V_m40) .* ( 3.56 * exp(0.079 * V) + 310000 * exp( 0.35 * V ) ) + (V_m40) ./ ( 0.13 * ( 1 + exp( -(V+10.66)/11.1 ) ) );
invtau_h = alpha_h + beta_h;
h_inf = alpha_h ./ invtau_h;
% j gate - slowly recovered inactivation
alpha_j = (~V_m40) .* ( -127140 * exp( 0.2444 * V ) - 0.00003474 * exp( -0.04391 * V ) ) .* (V + 37.78) ./ ( 1 + exp( 0.311 * (V + 79.23) ) );
beta_j = (~V_m40) .* ( 0.1212 * exp( -0.01052 * V) ./ ( 1 + exp( -0.1378 * (V + 40.14) ) ) ) + (V_m40) .* (0.3 * exp( -2.535e-7 * V ) ./ ( 1 + exp( -0.1 * (V + 32) ) ) );
invtau_j = alpha_j + beta_j;
j_inf = alpha_j ./ invtau_j;

%%% Transient outward K+ gates
% oa gate - activation
oa_inf = 1 ./ ( 1 + exp( -(V+20.47) / 17.54) );
invtau_oa = KQ10 * ( 0.65 ./ ( exp( -(V+10) / 8.5 ) + exp( -(V-30) / 59 ) ) + 0.65 ./ ( 2.5 + exp( (V+82) / 17 ) ) );
% oi gate - inactivation
oi_inf = 1 ./ ( 1 + exp( (V+43.1) / 5.3 ) );
invtau_oi = KQ10 * ( 1 ./ ( 18.53 + exp( (V+113.7) / 10.95 ) ) + 1 ./ ( 35.56 + exp( -(V+1.26) / 7.44 ) ) );

%%% Ultrarapid rectifier K+ gates
% ua gate - activation
ua_inf = 1 ./ ( 1 + exp( -(V+30.3) / 9.6 ) );
invtau_ua = invtau_oa;
% ui gate - inactivation
ui_inf = 1 ./ ( 1 + exp( (V-99.45) / 27.48 ) );
invtau_ui = KQ10 * (  1 ./ (21 + exp( -(V-185) / 28 ) ) + exp( (V-158) / 16 ) );

%%% Delayed rectifier fast component gate - xr
Vshift = V + 14.1;
xr_inf = 1 ./ ( 1 + exp( -Vshift / 6.5 ) );
alpha_xr = 0.0003 * Vshift ./ (1 - exp( -Vshift / 5) );
alpha_xr( abs(Vshift) < 1e-9 ) = 0.0015;                   % Remove divide by zero using limiting value
Vshift = V - 3.3328;
beta_xr = 7.3898e-5 * Vshift ./ ( exp(Vshift / 5.1237) - 1);
beta_xr( abs(Vshift) < 1e-9 ) = 0.000378361;               % Remove divide by zero using limiting value
invtau_xr = alpha_xr + beta_xr;

%%% Delayed rectifier slow component gate - xs
Vshift = V - 19.9;
xs_inf = 1 ./ sqrt( 1 + exp(-Vshift / 12.7) );
invtau_xs = 2 * ( 4e-5 * Vshift ./ ( 1 - exp( -Vshift / 17) ) + 3.5e-5 * Vshift ./ (exp(Vshift / 9) - 1 ) );
invtau_xs( abs(Vshift) < 1e-9 ) = 2 * (0.00068 + 0.000315);    % Remove divide by zero using limiting value

%%% L-type Ca2+ current gates
% d gate - activation
Vshift = V + 10;
expval = exp(-Vshift / 6.24);
d_inf = 1 ./ ( 1 + exp( -Vshift / 8 ) );
invtau_d = 0.035 * Vshift .* (1 + expval) ./ (1 - expval);
invtau_d( abs(Vshift) < 1e-9 ) = 0.035 * 2 * 6.24;            % Remove divide by zero using limiting value
% f gate - voltage-based inactivation
f_inf = 1 ./ ( 1 + exp( (V+28) / 6.9 ) );
invtau_f = (1/9) * ( 0.0197 * exp( -(0.0337^2 * (V + 10).^2 ) ) + 0.02 );
% fCa gate - [Ca2+]_i inactivation
fCa_inf = 1 ./ ( 1 + Ca_i / 0.00035);
invtau_fCa = 0.5 * ones(size(V));

%%% SR Ca2+ release gates
% u gate
Fn = 1e-12 * (Vol_rel * I_rel - (0.25 * I_CaL - 0.1 * I_NaCa) / F);
u_inf = 1 ./ ( 1 + exp( -(Fn - 3.4175e-13) / 1.367e-15 ) );
invtau_u = 1/8 * ones(size(V));
% v gate
v_inf = 1 - 1 ./ ( 1 + exp( -(Fn - 6.835e-14) / 1.367e-15 ) );
invtau_v = 1 ./ ( 1.91 + 2.09 * u_inf );
% w gate
Vshift = V - 7.9;
expval = exp(-Vshift / 5);
w_inf = 1 - 1 ./ ( 1 + exp( -(V-40) / 17 ) );
invtau_w = Vshift .* ( 1 + 0.3 * expval ) ./ ( 6 * (1 - expval ) );
invtau_w( abs(Vshift) < 1e-9 ) = 13 / 12;

% Rates of change for non-gating variables (all ion concentrations in this case)
b(:,1) = -( I_Na + I_bNa + 3*I_NaK + 3*I_NaCa ) / (F * Vol_i);
b(:,2) = -( I_to + I_Kur + I_Kr + I_Ks + I_K1 - 2*I_NaK ) / (F * Vol_i);
b(:,3) = ( -( I_CaL + I_bCa + I_CaP - 2*I_NaCa ) / (2 * Vol_i * F) + ( (I_up_leak - I_up) * Vol_up + I_rel * Vol_rel ) / Vol_i ) ./ ( 1 + TRPN_max * Km_TRPN ./ ( Ca_i + Km_TRPN ).^2 + CMDN_max * Km_CMDN ./ ( Ca_i + Km_CMDN ).^2 );
b(:,4) = ( I_tr - I_rel ) ./ ( 1 + CSQN_max * Km_CSQN ./ ( Ca_rel + Km_CSQN).^2 );
b(:,5) = I_up - I_up_leak - I_tr * Vol_rel / Vol_up;

% Gather rate constants for vectorisation purposes
Sinf = [ m_inf, h_inf, j_inf, oa_inf, oi_inf, ua_inf, ui_inf, xr_inf, xs_inf, d_inf, f_inf, fCa_inf, u_inf, v_inf, w_inf ];
invtau = [ invtau_m, invtau_h, invtau_j, invtau_oa, invtau_oi, invtau_ua, invtau_ui, invtau_xr, invtau_xs, invtau_d, invtau_f, invtau_fCa, invtau_u, invtau_v, invtau_w ];

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