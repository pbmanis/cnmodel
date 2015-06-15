from neuron import h
import neuron as nrn
import numpy as np
import scipy.optimize

from .cell import Cell
from .. import synapses
from ..util import nstomho

__all__ = ['Generic', 'GenericCell']

class Generic(Cell):

    @classmethod
    def create(cls, model='generic', **kwds):
        if model == 'generic':
            return GenericCell(**kwds)
        else:
            raise ValueError ('Cell type %s is unknown' % model)

    def make_psd(self, terminal, **kwds):
        from .. import cells
        
        pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma
        
        if isinstance(pre_cell, cells.SGC):
            return synapses.GluPSD(post_sec, terminal,
                                   ampa_gmax=1700.,
                                   nmda_ampa_ratio = 0.36,
                                   )
        elif isinstance(pre_cell, cells.DStellate):
            return synapses.GlyPSD(post_sec, terminal,
                                   psdType='glyslow',
                                   )
        else:
            raise TypeError("Cannot make PSD for %s => %s" % 
                            (pre_cell.__class__.__name__, 
                             self.__class__.__name__))


class GenericCell(Generic):
    """
    Generic cell model.
    Does not instantiate any particular conductances.
    Mainly used as a framework for detailed models, allowing integration with nrnlibrary.
    Pass the soma section so that it is used here (otherwise, one will be made for you)
    """

    def __init__(self, soma=None, nach='na', ttx=False, debug=False, species='mouse', type=None):
        """
        initialize the generic cell, using default parameters
        Modifications to the cell can be made by calling methods below.
        """
        super(GenericCell, self).__init__()
        if type == None:
            type = 'Generic'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Generic'}
        self.i_test_range=(-0.5, 0.5, 0.05)
        self.spike_threshold = -40
        self.vrange = [-85., -50.]  # set a default vrange for searching for rmp

        if soma is None:
            soma = h.Section(name="Generic_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1

        self.mechanisms = self.get_mechs(soma)
        # for mech in self.mechanisms:
        #     soma.insert(mech)
        # soma.ena = self.e_na
        # soma.ek = self.e_k
        # print 'soma dir: ', dir(soma().leak)
        # print soma().leak.e_leak
        # soma().leak.e_leak = self.e_leak
        self.add_section(soma, 'soma')
#        self.species_scaling(silent=True, species=species, type=type)  # set the default type II cell parameters
        self.get_mechs(soma)
        self.cell_initialize(vrange=self.vrange)
        if debug:
            print "<< generic: Generic cell model created >>"
        #print 'Cell created: ', self.status

    def species_scaling(self, species='guineapig', type='II', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        For a generic cell, nothing specific is changed
        """
        #print '\nSpecies scaling: %s   %s' % (species, type)

        self.status['species'] = species
        self.status['type'] = type

    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        adjust the sodium channel conductance
        :param soma: a soma object whose sodium channel complement will have it's 
        conductances adjusted depending on the channel type
        :return nothing:
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)

    def add_axon(self):
        pass

    def add_pumps(self):
        """
        Insert mechanisms for potassium ion, sodium ion, and a
        sodium-potassium pump at the soma.
        """
        pass

    def add_dendrites(self, debug=False):
        """
        Add a simple dendrite to the bushy cell.
        """
        pass
