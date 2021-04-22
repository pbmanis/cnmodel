TITLE kir.mod 
 
COMMENT

Inwardly rectifying potassium (Kir) currents
based on mathematical model from Leao et al (2012)
 Ceballos CC, Li S, Roque AC, Tzounopoulos T, Leão RM (2016) Ih Equalizes Membrane Input Resistance in a Heterogeneous Population of Fusiform Neurons in the Dorsal Cochlear Nucleus. Front Cell Neurosci 10:249 [PubMed]
ModelDB Accession: 206252
Reduced from fusiform.mod from Leao et al. 2012, 

ENDCOMMENT
 
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

 
NEURON {
    THREADSAFE
    SUFFIX kir
    USEION kir READ ekir WRITE ikir VALENCE 1
	RANGE gkir, gbar, ikir 	: Kir channels 
	GLOBAL ntau, ninf	: Kir channels 
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    dt (ms)
    gbar = 0.0005 (mho/cm2)    <0,1e9>
    q10tau = 1.0
    q10g = 1.0
}


STATE {
    n
}
 
ASSIGNED {
    celsius (degC)
    ikir (mA/cm2)
	ekir (mV)
    gkir (mho/cm2)
	ninf
    ntau (ms)
    qg()
    q10 ()
}


BREAKPOINT {
    SOLVE states METHOD cnexp

    gkir = qg*gbar*n
	ikir = gkir*(v - ekir)
}

INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
	q10 = q10tau^((celsius - 32)/10 (degC))
    ntau = 0.5 (ms)
	rates(v)
	n = ninf
}

UNITSOFF 

DERIVATIVE states {  
    rates(v)
    n' = (ninf - n) / ntau
}


PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	:"n" potassium activation system
    ninf = 1/(1+exp((v+85.48)/12.0)) 
    ntau = ntau/q10
}


UNITSON