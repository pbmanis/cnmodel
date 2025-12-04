:  Vector stream of events
:  From NEURON source: nrn/examples/nrniv/netcon/vecevent.mod
:  Modified to be compatible with Neuron V9.0 pbm 2025-02-12
    
NEURON {
	ARTIFICIAL_CELL VecStim
}

ASSIGNED {
	index
	etime (ms)
	space
}

INITIAL {
	index = 0
	element()
	if (index > 0) {
		net_send(etime - t, 1)
	}
}

NET_RECEIVE (w) {
	if (flag == 1) {
		net_event(t)
		element()
		if (index > 0) {
			net_send(etime - t, 1)
		}
	}
}

VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
#endif
ENDVERBATIM

PROCEDURE element() {
VERBATIM	
  { int i, size; double* px;
	i = (int)index;
	if (i >= 0) {
		IvocVect* vv = *((IvocVect**)(&space));
		if (vv) {
			size = vector_capacity(vv);
			px = vector_vec(vv);
			if (i < size) {
				etime = px[i];
				index += 1.;
			}else{
				index = -1.;
			}
		}else{
			index = -1.;
		}
	}
  }
ENDVERBATIM
}

PROCEDURE play() {
VERBATIM
	IvocVect** vv = (IvocVect**)(&space);
	*vv = (IvocVect*)0;
	if (ifarg(1)) {
		*vv = vector_arg(1);
	}
ENDVERBATIM
}
        

