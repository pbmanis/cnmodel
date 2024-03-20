TITLE Cerebellum Granule Cell Model

COMMENT
        Gaba A leakage
   
	Author: A. Fontana
	Last revised: 18.2.99
ENDCOMMENT
 
NEURON { 
	THREADSAFE
	SUFFIX GRCLKGG 
	NONSPECIFIC_CURRENT i
	RANGE e_gaba, gbar, i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gbar = 3e-5 (mho/cm2)  : 2.17e-5 
	celsius = 30 (degC)
	e_gaba = -65 (mV)
} 

ASSIGNED { 
	i (mA/cm2) 
}

INITIAL {

}

BREAKPOINT { 
	i = gbar*(v - e_gaba) 
} 
