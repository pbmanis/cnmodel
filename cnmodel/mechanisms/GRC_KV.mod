TITLE Cerebellum Granule Cell Model

COMMENT
        KDr channel
	Gutfreund parametrization
   
	Author: A. Fontana
	Last revised: 12.12.98
ENDCOMMENT

NEURON { 
	SUFFIX GRC_KV 
	USEION k READ ek WRITE ik 
	RANGE gbar, ik, gk, alpha_n, beta_n 
	RANGE Aalpha_n, Kalpha_n, V0alpha_n
	RANGE Abeta_n, Kbeta_n, V0beta_n
	RANGE n_inf, tau_n 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	:Kbeta_n = -0.0125 (/mV)
	
	Aalpha_n = -0.01 (/ms-mV)
	Kalpha_n = -10 (mV)
	V0alpha_n = -25 (mV)
	Abeta_n = 0.125 (/ms)
	
	Kbeta_n = -80 (mV)
	V0beta_n = -35 (mV)
	v (mV)  
	gbar= 0.003 (mho/cm2) : 0.0015
	ek = -84.69 (mV) 
	celsius = 30 (degC) 
} 

STATE { 
	n 
} 

ASSIGNED { 
	ik (mA/cm2) 
	n_inf 
	tau_n (ms) 
	gk (mho/cm2) 
	alpha_n (/ms) 
	beta_n (/ms) 
} 
 
INITIAL { 
	rate(v) 
	n = n_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	gk = gbar*n*n*n*n 
	ik = gk*(v - ek) 
	alpha_n = alp_n(v) 
	beta_n = bet_n(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	n' =(n_inf - n)/tau_n 
} 
 
FUNCTION alp_n(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-6.3(degC))/10(degC)) 
	alp_n = Q10*Aalpha_n*linoid(v-V0alpha_n, Kalpha_n)
} 
 
FUNCTION bet_n(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-6.3(degC))/10(degC)) 
	bet_n = Q10*Abeta_n*exp((v-V0beta_n)/Kbeta_n) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_n, b_n 
	TABLE n_inf, tau_n 
	DEPEND Aalpha_n, Kalpha_n, V0alpha_n, 
               Abeta_n, Kbeta_n, V0beta_n, celsius FROM -100 TO 30 WITH 13000 
	a_n = alp_n(v)  
	b_n = bet_n(v) 
	tau_n = 1/(a_n + b_n) 
	n_inf = a_n/(a_n + b_n) 
} 

FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
                linoid = x/(exp(x/y) - 1)
        }
}
