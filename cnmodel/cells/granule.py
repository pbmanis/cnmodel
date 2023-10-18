from neuron import h
from ..util import nstomho
import numpy as np
from .cell import Cell
from .. import data
from .. import synapses

__all__ = ['Granule', 'GranuleDefault']

class Granule(Cell):

    celltype = 'granule'
    scaled = False
    
    @classmethod
    def create(cls, model='GRC', **kwds):
        if model == 'GRC':
            return GranuleDefault(**kwds)
        else:
            raise ValueError ('Granule cell model %s is unknown', model)

    def __init__(self):
        Cell.__init__(self)
        self.spike_source = None  # used by DummyDStellate to connect VecStim to terminal

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
        self.pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma

        if psd_type == 'simple':
            if terminal.cell.celltype in ['sgc', 'dstellate', 'tuberculoventral', 'cartwheel']:
                weight = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='weight')
                tau1 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau1')
                tau2 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau2')
                erev = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='erev')
                return self.make_exp2_psd(post_sec, terminal, weight=weight, loc=loc,
                        tau1=tau1, tau2=tau2, erev=erev)
            else:
                raise TypeError("Cannot make simple PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))

        else:
            raise ValueError("Unsupported psd type %s for cartwheel cell (inputs not implemented yet)" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        if term_type == 'simple':
            return synapses.SimpleTerminal(self.pre_sec, post_cell, 
                                            **kwds)
        elif term_type == 'multisite':
            if post_cell.celltype in ['tuberculoventral', 'pyramidal']:
                nzones = data.get('cartwheel_synapse', species=self.species,
                        post_type=post_cell.celltype, field='n_rsites')
                delay = data.get('cartwheel_synapse', species=self.species,
                        post_type=post_cell.celltype, field='delay')
            else:
                raise NotImplementedError("No knowledge as to how to connect cartwheel cell to cell type %s" %
                                        type(post_cell))
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, spike_source=self.spike_source,
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)

class GranuleDefault(Granule):
    """
    Cochlear Nucleus Granule cell model.
    
    """
    def __init__(self, morphology=None, decorator=None, nach=None,
                ttx=False, temperature=None,
                species='mouse', modelType=None, modelName=None, 
                debug=False):
        """        
        Create a granule cell model, based on a the granule cell model frim
        Diwakar et al., J. Physiology.
        There are no variations available for this model.
        
        Parameters
        ----------
        morphology : string (default: None)
            Name of a .hoc file representing the morphology. This file is used to constructe
            an electrotonic (cable) model. 
            If None (default), then a "point" (really, single cylinder) model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
            by that name must exist. The default is naRsg, a resurgent sodium channel model.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            This flag duplicates the effects of tetrodotoxin in the model. Currently, the flag is not implemented.
        
        species: string (default 'rat')
            species defines the pattern of ion channel densities that will be inserted, according to 
            prior measurements in various species. Note that
            if a decorator function is specified, this argument is ignored as the decorator will
            specify the channel density.
            
        modelType: string (default: None)
            modelType specifies the subtype of the cell model that will be used.
            modelType is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            Only type "I" is recognized for the cartwheel cell model.

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
    
        debug: boolean (default: False)
            When True, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        """
        super(GranuleDefault, self).__init__()
        self.i_test_range={'pulse': (-0.02, 0.02, 0.002)}  # note that this might get reset with decorator according to channels
                                                    # The default values are set in the species_scaling routine
        if species == 'mouse':
            if modelType == None or modelType == 'GRC':
                modelName = 'granule'
                modelType = 'GRC'
                dataset = 'GRC_channels'
                temp = 34.

            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for {self.celltype:s} cells")

        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': self.celltype,
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.debug = debug
        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            # v_potassium = -80       # potassium reversal potential
            # v_sodium = 50           # somaodium reversal potential

            self.mechanisms = ['GRC_NA', 'GRC_LKG1', 'GRC_LKG2', 'GRC_KIR', 'GRC_KA', 'GRC_KM', 'GRC_KV',
                               'GRC_KCA', 'GRC_CA'] # , 'GRC_CALC']
            for mech in self.mechanisms:
                self.soma.insert(mech)
            # self.soma.insert('cadiff')
            self.soma.insert('GRC_CALC')
            self.species_scaling()  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        
        if debug:
            print( "<< Granule: Diwakar et al. model created >>")

    def get_cellpars(self, dataset, species='mouse', modelType='GRC'):
        pars = self.get_initial_pars(dataset, species, modelType)

        for g in ['soma_GRC_NA_gbar', 'soma_GRC_KV_gbar', 'soma_GRC_KM_gbar', 'soma_GRC_KA_gbar',
                  'soma_GRC_KCA_gbar',
                  'soma_GRC_KIR_gbar', 'soma_GRC_CA_gbar',
                #   'soma_GRC_LGK1_gl', #'soma_GRC_LGK2_gl',
                  'soma_e_k', 'soma_e_na', 'soma_e_ca', 'soma_e_leak',
                  'soma_Dia',
                  ]:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        return pars
        
    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        This scaling should be used ONLY for point models, as no other compartments
        are scaled.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I')
            definition of model type from RM03 models, type II or type II-I
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        
        Note
        ----
            For the cartwheel cell model, there is only a single scaling recognized. 
        """        
        assert self.scaled is False  # block double scaling!
        self.scaled = True
        
        if self.status['species'] != 'mouse':
            raise ValueError ('Granule,  species: only "mouse" is recognized')
        if self.status['modelType'] != 'GRC':
            raise ValueError ('Granule modelType: only "GRC" is recognized, got %s', modelType)
        self._valid_temperatures = (34.,)
        if self.status['temperature'] is None:
            self.set_temperature(34.)

        # self.i_test_range = {'pulse': (-0.05, 0.05, 0.005)}
       # self.spike_threshold = 0
        self.vrange = [-75., -52.]  # set a default vrange for searching for rmp
        
        # self.set_soma_size_from_Cm(self.pars.cap)
        print(dir(self.pars))
        self.set_soma_size_from_Diam(self.pars.soma_Dia)

        self.soma().GRC_NA.gbar = self.g_convert(self.pars.soma_GRC_NA_gbar, self.pars.units, self.somaarea)
        self.soma().GRC_KV.gbar = self.g_convert(self.pars.soma_GRC_KV_gbar, self.pars.units, self.somaarea)
        self.soma().GRC_KA.gbar = self.g_convert(self.pars.soma_GRC_KA_gbar, self.pars.units, self.somaarea)
        self.soma().GRC_KIR.gbar = self.g_convert(self.pars.soma_GRC_KIR_gbar, self.pars.units, self.somaarea)
        self.soma().GRC_KCA.gbar = self.g_convert(self.pars.soma_GRC_KCA_gbar, self.pars.units, self.somaarea)
        self.soma().GRC_CA.gbar = self.g_convert(self.pars.soma_GRC_CA_gbar, self.pars.units, self.somaarea)
        # self.soma().GRC_LGK1.gl = self.g_convert(self.pars.soma_GRC_LGK1_gl, self.pars.units, self.somaarea)
      
        self.soma().ena = self.pars.soma_e_na # 50
        self.soma().ek = self.pars.soma_e_k # -80
        self.soma().eca = self.pars.soma_e_ca # 50
        # self.soma().el = self.pars.soma_e_leak
        
        self.check_temperature()
        if not silent:
            print( 'set cell as: ', self.status['species'])
       # print 'set up'
        
    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point.
        Overrides i_currents in cells.py, because this model uses conductances
        that are not specified in the default cell mode.
        
        Parameters
        ----------
        V : float, mV (no default)
            Voltage at which the current for each conductance is computed.
        
        Returns
        -------
        I : float, nA
             The sum of the currents at steady-state for all of the conductances.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.celsius = self.status['temperature']
        h.finitialize()
        self.ix = {}

        if 'GRC_NA' in self.mechanisms:
            self.ix['GRC_NA'] = self.soma().GRC_NA.gna*(V - self.soma().ena)
        if 'GRC_KV' in self.mechanisms:
             self.ix['GRC_KV'] = self.soma().GRC_KV.gk*(V - self.soma().ek)
        if 'GRC_KA' in self.mechanisms:
             self.ix['GRC_KA'] = self.soma().GRC_KA.gk*(V - self.soma().ek)
        if 'GRC_KM' in self.mechanisms:
             self.ix['GRC_KM'] = self.soma().GRC_KM.gk*(V - self.soma().ek)
        if 'GRC_KCA' in self.mechanisms:
             self.ix['GRC_KCA'] = self.soma().GRC_KCA.gk*(V - self.soma().ek)
        if 'GRC_KIR' in self.mechanisms:
             self.ix['GRC_KIR'] = self.soma().GRC_KIR.gk*(V - self.soma().ek)
        if 'GRC_CA' in self.mechanisms:
             self.ix['GRC_CA'] = self.soma().GRC_CA.gca*(V - self.soma().eca)
        if 'GRC_LKG1' in self.mechanisms:  # resistive leak
             self.ix['GRC_LKG1'] = self.soma().GRC_LKG1.gl*(V - self.soma().GRC_LKG1.el)
        if 'GRC_LKG2' in self.mechanisms:  # resting GABA conductance
             self.ix['GRC_LKG2'] = self.soma().GRC_LKG2.ggaba*(V - self.soma().GRC_LKG2.egaba)
        # if 'GRC_CALC' in self.mechanisms:
        #      self.ix['GRC_CALC'] = self.soma().GRC_CALC.gca*(V - self.soma().GRC_CALC.eca)
        # leak
        # if 'lkpkj' in self.mechanisms:
        #     self.ix['lkpkj'] = self.soma().lkpkj.gbar*(V - self.soma().GRC.ek)
        return np.sum([self.ix[i] for i in self.ix])

    def ghk(self, v, ci, co, z):
        """
        GHK flux equation, used to calculate current density through calcium channels
        rather than standard Nernst equation.
        
        Parameters
        ----------
        v : float, mV
            voltage for GHK calculation
        ci : float, mM
            internal ion concentration
        co : float, mM
            external ion concentraion
        z : float, no units
            valence
        
        Returns
        -------
        flux : A/m^2

        """
        F = 9.6485e4  # (coul)
        R = 8.3145 # (joule/degC)
        T = h.celsius + 273.19  # Kelvin
        E = (1e-3) * v  # convert mV to V
        Ci = ci + (self.soma().cap.monovalPerm) * (self.soma().cap.monovalConc)  #       : Monovalent permeability
        if (np.fabs(1-np.exp(-z*(F*E)/(R*T))) < 1e-6):  #denominator is small -> Taylor series
            ghk = (1e-6) * z * F * (Ci-co*np.exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
        else:
            ghk = (1e-6) * z**2.*(E*F**2.)/(R*T)*(Ci-co*np.exp(-z*(F*E)/(R*T)))/(1-np.exp(-z*(F*E)/(R*T)))
        return ghk
