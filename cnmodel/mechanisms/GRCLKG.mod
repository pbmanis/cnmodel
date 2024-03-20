TITLE Cerebellum Granule Cell Model

COMMENT
        Leakage
   
	Author: A. Fontana
	Last revised: 18.12.98
ENDCOMMENT
 
NEURON { 
	THREADSAFE
	SUFFIX GRCLKG
	NONSPECIFIC_CURRENT i
	RANGE e_l, gbar, i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gbar = 5.68e-5 (mho/cm2)
	celsius = 30 (degC)
	e_l =  -58.0 : -16.5 (mV) : resting at -70 mV	:-11 resting at -68 mV
	: -58 : to make it 
} 

ASSIGNED { 
	i (mA/cm2) 
}

INITIAL {

}

BREAKPOINT { 
	i = gbar*(v - e_l) 
} 
