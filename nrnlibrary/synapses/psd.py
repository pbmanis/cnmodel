
class PSD(object):
    """
    Base class for postsynaptic density mechanisms, possibly including cleft. 
    May accept either NetCon or pointer inputs from a Terminal, and directly
    modifies the membrane potential and/or ion concentrations of the 
    postsynaptic cell.    
    """
    def __init__(self, pre_sec, post_sec, terminal):
        pass
    
