function [I_ion, S, I_Na, I_CaL, I_Kr, I_Ks, I_to, I_K1, I_NaK, I_NaCa, I_up, I_rel] = RLUpdateTT06M(V, S, dt, I_stim)







%%% Variable current strengths
g_to=inputs(1);       			% 0.1652
g_Kur_scale=inputs(2);     		% 1
g_Ca_L=inputs(3);     			% 0.12375
g_Kr=inputs(4);      			% 0.029411765
g_Ks=inputs(5);      			% 0.12941176
g_K1=inputs(6);      			% 0.09
g_Na=inputs(7);     			% 7.8
g_NaK=inputs(8);    			% 0.59933874
g_NaCa=inputs(9);    			% 1600
i_up_max_scale = inputs(10);	% 1
J_rel = inputs(11);				% 30

stim_start = inputs(12);
stim_BCL = inputs(13);
stim_dur = inputs(14);
stim_ampl = inputs(15);
    
    R = 8.3143;            % Gas constant
    T = 309.65;            % Temperature (K)
    F = 96.4867;           % Faraday constant
    Cm = 100;              % Membrane capacitance
    
    Na_o = 149.42;  		% Extracellular Na concentration
	K_o = 4.5;				% Extracellular K concentration	
    KQ10 = 3;				% Q10 value
    k_NaK_Na = 10;			% Saturation constant for Na in NaK pump
	k_NaK_K = 1.5;			% Saturation constant for K in NaK pump
    g_B_Na = 6.744375e-4;	% Background Na+ current
	g_B_Ca = 1.131e-3;		% Background Ca2+ current
	Ca_o = 1.8;				% Extracellular Ca2+ concentration
    k_NaCa_Na = 87.5;		% Saturation constant for Na in NaCa pump
	k_NaCa_Ca = 1.38;		% Saturation constant for Ca in NaCa pump
	k_sat = 0.1;			% Other saturation constant in NaCa pump
	gamma = 0.35;			% 
	i_CaP_max = 0.275;
	i_up_max = 0.005;
	K_up = 0.00092;
	Ca_up_max = 15;
    CMDN_max = 0.05;
	TRPN_max = 0.07;
	CSQN_max = 10;
    k_CMDN = 0.00238;			% Saturation constant for calmodulin
	k_TRPN = 0.0005;			% Saturation constant for troponin
	k_CSQN = 0.8;				% Saturation constant for calcium sequestrin
    Vol_tot = 20100;				% Cell volume
	Vol_i = 0.68 * Vol_tot;			% Cell volume for intracellular ions
	Vol_rel = 0.0048 * Vol_tot;		% Cell volume for 'rel' component
	Vol_up = 0.0552 * Vol_tot;		% Cell volume for 'up' component

	% Calculate the Nernst Potential differentials
	RTonF = R * T / F;
	VminusENa = y(1) - RTonF * log( Na_o / y(2) );
	VminusEK = y(1) - RTonF * log( K_o / y(6) );
	VminusECa = y(1) - RTonF / 2 * log( Ca_o / y(13) );
	
	% Calculate some useful constants to speed things up a little
	expFVonRT = exp(y(1) / RTonF);
	sigma = (1/7) * ( exp( Na_o / 67.3 ) - 1);
	
	% Calculate the values of the derived currents
	i_Na = Cm * g_Na * y(3)^3 * y(4) * y(5) * VminusENa;
	i_K1 = Cm * g_K1 * VminusEK / (1 + exp(0.07*(y(1)+80)));
	i_to = Cm * g_to * y(7)^3 * y(8) * VminusEK;
	g_Kur = g_Kur_scale * (0.005 + 0.05 / (1 + exp(-(y(1) - 15)/13)) );
	i_Kur = Cm * g_Kur * y(9)^3 * y(10) * VminusEK;
	i_Kr = Cm * g_Kr * y(11) * VminusEK / (1 + exp((y(1)+15)/22.4));
	i_Ks = Cm * g_Ks * y(12)^2 * VminusEK;
	i_B_Na = Cm * g_B_Na * VminusENa;
	i_B_Ca = Cm * g_B_Ca * VminusECa;
	i_NaK = Cm * g_NaK * K_o / ( (1 + (k_NaK_Na / y(2))^1.5) * (K_o + k_NaK_K) * (1 + 0.1245 * expFVonRT^-0.1 + 0.0365 * sigma / expFVonRT) );
	i_NaCa = Cm * g_NaCa * (expFVonRT^gamma * y(2)^3 * Ca_o - expFVonRT^(gamma-1) * Na_o^3 * y(13) ) / ( (k_NaCa_Na^3 + Na_o^3) * (k_NaCa_Ca + Ca_o) * (1 + k_sat * expFVonRT^(gamma-1)) );
	i_CaP = Cm * i_CaP_max * y(13) / ( 0.0005 + y(13) );
	i_Ca_L = Cm * g_Ca_L * y(14) * y(15) * y(16) * (y(1) - 65);
	i_rel = J_rel * y(18)^2 * y(19) * y(20) * (y(17) - y(13));
	i_tr = ( y(21) - y(17) ) / 180;
	i_up = i_up_max_scale * i_up_max / (1 + K_up / y(13));
	i_up_leak = i_up_max * y(21) / Ca_up_max;

	%cur_time = (t - stim_start) - floor( (t - stim_start) / stim_BCL ) * stim_BCL;
	cur_time = t - stim_start - floor( (t - stim_start) / stim_BCL) * stim_BCL;

	%if cur_time <= stim_dur && cur_time >= 0
	if cur_time <= stim_dur && t >= stim_start
		i_stim = stim_ampl;
	else
		i_stim = 0;
	end;
	
	% Build the ODE RHS
	dy = zeros(21,1);
	
	% Rate of change for V
	dy(1) = -(i_Na + i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_NaCa + i_CaP + i_Ca_L + i_stim) / Cm;
	
	% Rate of change for Na_i
	dy(2) = (-3 * i_NaK - 3 * i_NaCa - i_B_Na - i_Na) / (F * Vol_i);
	
	% Rate of change for m (very cleaned up)
	vshift = y(1) + 47.13;
	if abs(vshift) < 1e-10
		alpha_m = 3.2;
	else
		alpha_m = 0.32 * vshift / (1 - exp(-0.1 * vshift));
	end;
	dy(3) = alpha_m  - y(3) * (alpha_m + 0.08 * exp(-y(1)/11));
		
	% Rate of change for h
	if y(1) < -40
		alpha_h = 0.135 * exp( -(y(1) + 80) / 6.8 );
		beta_h = 3.56 * exp(0.079 * y(1)) + 310000 * exp(0.35 * y(1));	
	else
		alpha_h = 0;
		beta_h = 1 / (0.13 * (1 + exp( -(y(1) + 10.66)/11.1)));
	end;
	dy(4) = (alpha_h - y(4) * (alpha_h + beta_h));
		
	% Rate of change for j	
	if y(1) < -40
		alpha_j = (-127140 * exp(0.2444 * y(1)) - 3.474e-5 * exp(-0.04391 * y(1))) * (y(1) + 37.78) / (1 + exp(0.311 * (y(1) + 79.23)));
		beta_j = 0.1212 * exp(-0.01052 * y(1)) / (1 + exp(-0.1378 * (y(1) + 40.14) ) );	
	else
		alpha_j = 0;
		beta_j = 0.3 * exp(-2.535e-7 * y(1)) / (1 + exp(-0.1 * (y(1) + 32) ) );
	end;
	dy(5) = alpha_j - y(5) * (alpha_j + beta_j);

	% Rate of change for K_i
	dy(6) = (2 * i_NaK - i_K1 - i_to - i_Kur - i_Kr - i_Ks) / (F * Vol_i);

	% Rate of change for oa
	dy(7) = (1 / (1 + exp(-(y(1)+20.47)/17.54) ) - y(7)) * KQ10 * ( 0.65 / ( exp(-(y(1)+10)/8.5) + exp(-(y(1)-30)/59) ) + 0.65 / (2.5 + exp((y(1)+82)/17) ) );
	
	% Rate of change for oi
	dy(8) = ( 1 / (1 + exp((y(1) + 43.1)/5.3)) - y(8) ) * KQ10 * ( 1 / (18.53 + exp((y(1)+113.7)/10.95) )  +  1 / (35.56 + exp(-(y(1)+1.26)/7.44)) );
	
	% Rate of change for ua
	dy(9) = ( 1 / (1 + exp(-(y(1) + 30.3)/9.6)) - y(9) ) * KQ10 * ( 0.65 / ( exp(-(y(1)+10)/8.5) + exp(-(y(1)-30)/59) )  +  0.65 / (2.5 + exp((y(1)+82)/17)) );
	
	% Rate of change for ui
	dy(10) = ( 1 / (1 + exp((y(1) - 99.45)/27.48)) - y(10) ) * KQ10 * ( 1 / ( 21 + exp(-(y(1)-185)/28) )  +  exp((y(1)-158)/16) );

	% Rate of change for xr
	vshift = y(1) - 3.3328;
	if abs(vshift) < 1e-10
		beta_xr = 0.000378361;
	else
		beta_xr = 7.3898e-5 * vshift / (exp(vshift / 5.1237)-1);
	end;
	vshift = y(1) + 14.1;
	if abs(vshift) < 1e-10
		alpha_xr = 0.0015;
	else
		alpha_xr = 0.0003 * vshift / (1 - exp(-vshift / 5));
	end;
	dy(11) = ( 1 / (1 + exp(-vshift/6.5)) - y(11) ) * (alpha_xr + beta_xr);
	
	% Rate of change for xs
	vshift = y(1) - 19.9;
	if abs(vshift) < 1e-10
		alpha_xs = 0.00068;
		beta_xs = 0.000315;
	else
		alpha_xs = 4e-5 * vshift / (1 - exp(-vshift / 17));
		beta_xs = 3.5e-5 * vshift / (exp(vshift / 9)-1);
	end;
	dy(12) = ( (1 + exp(-vshift/12.7))^-0.5 - y(12) ) * 2 * (alpha_xs + beta_xs);
	
	% Rate of change for Ca_i
	B1 = ( (2*i_NaCa - i_CaP - i_Ca_L - i_B_Ca) / (2 * F) + (Vol_up * (i_up_leak - i_up) + i_rel * Vol_rel) ) / Vol_i;
	B2 = 1 + TRPN_max * k_TRPN / (y(13) + k_TRPN)^2 + CMDN_max * k_CMDN / (y(13) + k_CMDN)^2;
	dy(13) = B1 / B2;
	
	% Rate of change for d
	vshift = y(1) + 10;
	expval = exp(-vshift/6.24);
	if abs(vshift) < 1e-10
		tau_d_inv = (1 + expval)/4.579;
	else
		tau_d_inv = 0.035 * vshift * (1 + expval) / (1 - expval);
	end;
	dy(14) = ( 1 / (1 + exp(-vshift/8)) - y(14)) * tau_d_inv;
	
	% Rate of change for f
	dy(15) = ( 1 / (1 + exp((y(1)+28)/6.9)) - y(15)) * (0.0197 * exp(-(0.0337 * (y(1)+10))^2) +0.02) / 9;
	
	% Rate of change for f_Ca
	dy(16) = 0.5 * ( 1 / (1 + y(13)/0.00035) - y(16) );
	
	% Rate of change for Ca_rel
	dy(17) = (i_tr - i_rel) / (1 + CSQN_max * k_CSQN / (y(17) + k_CSQN)^2 );
	
	% Rate of change for u
	Fn = 1000 * 1e-15 * (Vol_rel * i_rel - (0.25 * i_Ca_L - 0.1 * i_NaCa) / F);
	invexpval = 1 / (1 + exp(-(Fn-3.4175e-13)/1.367e-15) );
	dy(18) = ( invexpval - y(18) ) / 8;
	
	% Rate of change for v
	dy(19) = ( 1 - 1 / (1 + exp(-(Fn-6.835e-14)/1.367e-15) ) - y(19) ) / (1.91 + 2.09 * invexpval);
	
	% Rate of change for w
	vshift = y(1) - 7.9;
	if abs(vshift) < 1e-10
		tau_w_inv = 1.083333;
	else
		expval = exp(-vshift/5);
		tau_w_inv = vshift * (1 + 0.3 * expval) / (6 * (1 - expval));
	end;
	dy(20) = ( ( 1 - 1 / (1 + exp(-(y(1)-40)/17))) - y(20) ) * tau_w_inv;
	
	% Rate of change for Ca_up
	dy(21) = i_up - i_up_leak - i_tr * Vol_rel / Vol_up;