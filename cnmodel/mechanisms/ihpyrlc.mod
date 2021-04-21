TITLE ihpyrlc.mod
 
COMMENT

Ih current
based on mathematical model from Leao et al (2012)

ENDCOMMENT
 
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

 
? interface
NEURON {
    THREADSAFE
    SUFFIX ihpyrlc
    NONSPECIFIC_CURRENT i
    RANGE gh, gbar, eh, i   : Ih channels
    RANGE kh_m_tau, kh_n_tau  : time constants for Ih channels
    GLOBAL kh_m_inf, kh_n_inf, aih 
}

 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    dt (ms)
    gbar = 0.00054 (mho/cm2) <0,1e9>
}


STATE {
    khm khn
}
 
ASSIGNED {
    celsius (degC)
    gh (mho/cm2)
    eh (mV)
    ih (mA/cm2)
    kh_m_inf 
    kh_n_inf
    kh_m_tau (ms)
    kh_n_tau (ms)
    aih
    i (ma/cm2)

}


LOCAL mexp, hexp, nexp, kh_m_exp, kh_n_exp, mexp_p, nexp_ir
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp

    aih = 0.5*khm+0.5*khn
    gh = gbar*aih
    i = gh*(v - eh)
}
? currents

UNITSOFF 
 

INITIAL {
    rates(v)

    khm = kh_m_inf
    khn = kh_n_inf
}


? states
DERIVATIVE states {  
    rates(v)

    khm' = (kh_m_inf - khm) / kh_m_tau
    khn' = (kh_n_inf - khn) / kh_n_tau

}
 
LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
    LOCAL  alpha, beta, sum
    TABLE kh_m_inf, kh_n_inf, kh_m_tau, kh_n_tau DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
    q10 = 3^((celsius - 32)/10)

    :"kh" adaptation of Destexhe hyp-activated cation current by Patrick Kanold
    kh_m_inf = kh_m(v) 
    kh_n_inf = kh_n(v)
    kh_m_tau = kh_mt(v)
    kh_n_tau = kh_nt(v) 

}

FUNCTION kh_m(x) { : h activation
    kh_m = 1/(1+exp((x+87)/8.9))
}

FUNCTION kh_n(x) { : h activation
    kh_n = 1/(1+exp((x+87)/8.9)) : same as kh_m
}

FUNCTION kh_mt(v) { : h fast time constant
    kh_mt = 100 + exp((v+183.6)/30.48)
}

FUNCTION kh_nt(v) { : h slow time constant
    kh_nt = 700 + exp((v+188.6)/11.2)/(1+exp((v+105)/5.5))
}

UNITSON