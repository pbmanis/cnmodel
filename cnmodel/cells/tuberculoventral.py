from neuron import h
import numpy as np
#import neuron as nrn

from .cell import Cell
#from .. import synapses
from ..util import nstomho
from ..util import Params
from .. import data

__all__ = ['Tuberculoventral'] 


class Tuberculoventral(Cell):
    
    type = 'Tuberculoventral'

    @classmethod
    def create(cls, modelType='TVmouse', **kwds):
        if modelType in ['TVmouse', 'I']:
            return Tuberculoventral(**kwds)
        else:
            raise ValueError ('Tuberculoventral type %s is unknown', modelType)

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is to try to pass the default unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dict of options. Two are currently handled:
        postsize : expect a list consisting of [sectionno, location (float)]
        AMPAScale : float to scale the ampa currents
        
        """
        if 'postsite' in kwds:  # use a defined location instead of the default (soma(0.5)
            postsite = kwds['postsite']
            loc = postsite[1]  # where on the section?
            uname = 'sections[%d]' % postsite[0]  # make a name to look up the neuron section object
            post_sec = self.hr.get_section(uname)  # Tell us where to put the synapse.
        else:  # place synapses on the soma at the middle of the section.
            loc = 0.5
            post_sec = self.soma
        
        if psd_type == 'simple':
            return self.make_exp2_psd(post_sec, terminal, loc=loc)
        elif psd_type == 'multisite':
            if terminal.cell.type == 'sgc':
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                if 'AMPAScale' in kwds:
                    AMPA_gmax = AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    NMDA_gmax = NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, AMPAR_gmax, NMDAR_gmax, loc=loc)
            elif terminal.cell.type == 'dstellate':
                # WBI input determines tone/noise response for TV cells (Nelken and Young)
                return self.make_gly_psd(post_sec, terminal, type='glyfast', loc=loc)
            elif terminal.cell.type == 'tuberculoventral': 
                # Evidence from Kuo et al that TV cells have synapses with each other
                return self.make_gly_psd(post_sec, terminal, type='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.type, self.type))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class Tuberculoventral(Tuberculoventral):
    """
    Tuberculoventral Neuron (DCN) base model
    Adapted from T-stellate model, using target parameters from Kuo et al. J. Neurophys. 2012
    """
    def __init__(self, morphology=None, decorator=None, nach=None, ttx=False,
                species='mouse', modelType=None, debug=False):
        """
        Initialize a DCN Tuberculoventral cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
        Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
        Changing "species" to mouse or cat (scales conductances)
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
        
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels aer inserted into the first soma section, and the
            rest of the structure is "bare".
    
        nach : string (default: 'na')
            nach selects the type of sodium channel that will be used in the model. A channel mechanims
            by that name must exist. 
    
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
    
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.
        
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point models.
        
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
        
        Returns
        -------
        Nothing
        """
        super(Tuberculoventral, self).__init__()
        if modelType == None:
            modelType = 'TVmouse'
        if nach == None:
            nach = 'nacncoop'
            
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Tuberculoventral',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}

        self.i_test_range = {'pulse': (-0.4, 0.6, 0.02)}
        self.vrange = [-80., -60.]  # set a default vrange for searching for rmp
        
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            print "<< Tuberculoventral model: Creating point cell >>"
            soma = h.Section(name="Tuberculoventral_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            print "<< Tuberculoventral model: Creating structured cell >>"
            self.set_morphology(morphology_file=morphology)

        if decorator is None:   # basic model, only on the soma ("point")
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True, species=species, modelType=modelType)  # adjust the default parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        
        self.get_mechs(self.soma)
#        self.cell_initialize(vrange=self.vrange)
        if debug:
                print "<< Tuberculoventral cell model created >>"

    def species_scaling(self, species='guineapig', modelType='TVmouse', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I-c')
            definition of model type from RM03 models, type I-c or type I-t
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        soma = self.soma
        print 'species: ', species, 'modelType: ', modelType
        if species == 'mouse' and modelType == 'TVmouse':
            """#From Kuo 150 Mohm, 10 msec tau
            Firing at 600 pA about 400 Hz
            The values for the conductances were determined from
            a brute_force exploration of the parameter space
            for Na and K channel conductances, getting 380 Hz at 600 pA at 34C
            Input resistance and Vm is ok, time constant is short
                *** Rin:       143  tau:       2.9   v:  -68.4
            Attempts to get longer time constant - cannot keep rate up.
            """
            # Adapted from TStellate model type I-c'
            print 'decorate as TVmouse'
            self.vrange=[-80., -55.]
            self.set_soma_size_from_Cm(35.0)
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.adjust_na_chans(soma, gbar=8000.)
            #soma.ek = -77.0 
            soma().kht.gbar = nstomho(2000.0, self.somaarea)
            soma().ka.gbar = nstomho(65.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.5, self.somaarea)  # 1.25
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(4.25, self.somaarea)  # 5.5
            soma().leak.erev = -70.0
            self.axonsf = 0.5
        else:
            raise ValueError('Species %s or species-type %s is not recognized for Tuberculoventralcells' % (species, type))

        self.status['species'] = species
        self.status['modelType'] = modelType
        self.check_temperature()


    def channel_manager(self, modelType='RM03'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, it's private
        _biophys function) to decorate the cell membrane.
        
        Parameters
        ----------
        modelType : string (default: 'RM03')
            A string that defines the type of the model. Currently, 3 types are implemented:
            RM03: Rothman and Manis, 2003 somatic densities for guinea pig
            XM13: Xie and Manis, 2013, somatic densities for mouse
            XM13PasDend: XM13, but with only passive dendrites, no channels.
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        
        This routine defines the following variables for the class:
            
            - conductances (gBar)
            - a channelMap (dictonary of channel densities in defined anatomical compartments)
            - a current injection range for IV's (when testing)
            - a distance map, which defines how selected conductances in selected compartments
                will change with distance. This includes both linear and exponential gradients,
                the minimum conductance at the end of the gradient, and the space constant or
                slope for the gradient.
        
        """
        if modelType == 'TVmouse':
            print 'decorate as tvmouse'
#            totcap = 95.0E-12  # Tuberculoventral cell (type I), based on stellate, adjusted for Kuo et al. TV firing
            self.set_soma_size_from_Section(self.soma)
            totcap = self.totcap
            refarea = self.somaarea # totcap / self.c_m  # see above for units
            self.gBar = Params(nabar=1520.0E-9/refarea, # these values may not corrrespond to the values for point model.
                               khtbar=160.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               kabar=65.0/refarea,
                               ihbar=1.25E-9/refarea,
                               leakbar=5.5E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
                         'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                             'ihvcn': 0., 'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar,
                         'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.3, 0.6, 10)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }
        else:
            raise ValueError('model type %s is not implemented' % modelType)
        
    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        Adjust the sodium channel conductance, depending on the type of conductance
        
        Parameters
        ----------
        soma : NEURON section object (required)
            This identifies the soma object whose sodium channel complement will have it's
            conductances adjusted depending on the sodium channel type
        gbar : float (default: 1000.)
            The "maximal" conductance to be set in the model.
        debug : boolean (default: False)
            A flag the prints out messages to confirm the operations applied.
            
        Returns
        -------
        Nothing
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)
        nach = self.status['na']
        if nach == 'nacncoop':
            print('TV model using nacncoop')
            soma().nacncoop.gbar = gnabar
            soma().nacncoop.KJ = 2000.
            soma().nacncoop.p = 0.25
            soma.ena = self.e_na
            if debug:
                print 'nacncoop gbar: ', soma().nacncoop.gbar
        elif nach == 'jsrna':
            soma().jsrna.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'jsrna gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "Tuberculoventral using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        elif  nach == 'nacn':
            soma().nacn.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().nacn.gbar
        else:
            raise ValueError("Tuberculoventral setting Na channels: channel %s not known" % nach)


