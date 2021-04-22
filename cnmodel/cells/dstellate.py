from __future__ import print_function
from neuron import h
from .cell import Cell
from ..util import Params
from .. import synapses
from .. import data

__all__ = ['DStellate', 'DStellateRothman', 'DStellateEager']


class DStellate(Cell):
    
    celltype = 'dstellate'
    scaled = False
    
    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return DStellateRothman(**kwds)
        elif model == 'Eager':
            return DStellateEager(**kwds)
        elif model == 'dummy':
            return DummyDStellate(**kwds)
        else:
            raise ValueError ('DStellate type %s is unknown', type)

    def __init__(self):
        Cell.__init__(self)
        self.spike_source = None  # used by DummyDStellate to connect VecStim to terminal

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is designed to pass the unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dictionary of options. 
            Two are currently handled:
            postsize : expect a list consisting of [sectionno, location (float)]
            AMPAScale : float to scale the ampa currents
        
        """
        if 'postsite' in kwds:  # use a defined location instead of the default (soma(0.5)
            postsite = kwds['postsite']
            loc = postsite[1]  # where on the section?
            uname = 'sections[%d]' % postsite[0]  # make a name to look up the neuron section object
            post_sec = self.hr.get_section(uname)  # Tell us where to put the synapse.
        else:
            loc = 0.5
            post_sec = self.soma
        
        if psd_type == 'simple':
            if terminal.cell.celltype in ['sgc', 'dstellate', 'tuberculoventral']:
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

        
        elif psd_type == 'multisite':
            if terminal.cell.celltype == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                # old values:
                # AMPA_gmax = 0.22479596944138733*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                # NMDA_gmax = 0.12281291946623739*1e3
                if 'AMPAScale' in kwds:
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)

            elif terminal.cell.celltype == 'dstellate':
                # Get GLY kinetic constants from database 
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            elif terminal.cell.celltype == 'tuberculoventral':
                # Get GLY kinetic constants from database 
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                                (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        if term_type == 'simple':
            return synapses.SimpleTerminal(self.soma, post_cell, spike_source=self.spike_source,
                     **kwds)
        elif term_type == 'multisite':
            if post_cell.celltype in ['dstellate', 'tuberculoventral', 'pyramidal', 'bushy', 'tstellate']:
                nzones = data.get('dstellate_synapse', species=self.species,
                        post_type=post_cell.celltype, field='n_rsites')
                delay = data.get('dstellate_synapse', species=self.species,
                        post_type=post_cell.celltype, field='delay')
            else:
                raise NotImplementedError("No knowledge as to how to connect D stellate cell to cell type %s" %
                                        type(post_cell))
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, spike_source=self.spike_source,
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)


class DStellateRothman(DStellate):
    """
    VCN D-stellate model:
    as a type I-II from Rothman and Manis, 2003
    """
    def __init__(self, morphology=None, decorator=None,  nach=None, ttx=False,
                species='guineapig', modelType=None, modelName=None, debug=False):
        """
        initialize a radial stellate (D-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I-II cell.
        Modifications to the cell can be made by calling methods below. These include:
                
            * changing the sodium channel
            * Changing "species" to mouse or cat (scales conductances)
            * Shifting model type
                
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
            by that name must exist. A value of None will set the channel to a default for the model (nacn).
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            This flag duplicates the effects of tetrodotoxin in the model. Currently, the flag is not implemented.
        
        species: string (default 'guineapig')
            species defines the pattern of ion channel densities that will be inserted, according to 
            prior measurements in various species. Note that
            if a decorator function is specified, this argument is ignored as the decorator will
            specify the channel density.
            
        modelType: string (default: None)
            modelType specifies the subtype of the cell model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
        debug: boolean (default: False)
            When True, there will be multiple printouts of progress and parameters.
            
            
        Returns
        -------
            Nothing
        """
        
        super(DStellateRothman, self).__init__()
        if modelType == None:
            modelType = 'I-II'
        if species == 'guineapig':
            modelName = 'RM03'
            dataset = 'RM03_channels'
            temp = 22.
        elif species == 'mouse':
            temp = 34.
            if modelName is None:
                modelName = 'XM13'
                dataset = 'XM13_channels'
            elif modelName  == 'XM13nacncoop':
                dataset = 'XM13_channels_nacncoop'
            else:
                raise ValueError(f"ModelName {modelName:s} not recognized for mouse {self.celltype:s} cells")

        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        self.debug = debug

        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName,
                       'ttx': ttx, 'name': 'DStellate',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.spike_threshold = -40.  # matches threshold in released CNModel (set in base cell class)
        
        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype
        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', self.pars.natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ena = self.e_na
            self.soma.ek = self.e_k
            self.soma().leak.erev = self.e_leak
            self.c_m = 0.9
            self.set_soma_size_from_Cm(12.0)
            self.species_scaling(silent=True)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)

        if self.debug:
            print('selfcelltype: ', self.celltype.title())
            print('morpholog: ', morphology)
            print (f"<< {self.celltype.title():s} model: Creating cell with morphology from {str(morphology):s} >>" )

    def get_cellpars(self, dataset, species='guineapig', modelType='I-II'):
        """
        Retrieve parameters for the specifed model type and species from the data tables
        
        dataset : str (no default)
            name of the data table to use
        
        species : str (default: 'guineapig')
            Species table to use
        
        modelType : str (default: 'I-II')
            Model type to get parameters from the table.
        
        """
        pars = self.get_initial_pars(dataset, species, modelType)

        if self.status['modelName'] == 'RM03':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ka_gbar', 'ih_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'mGBC':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        else:
            raise ValueError(f"get_cellpars: Model name {self.status['modelName']} is not yet implemented for cell type {self.celltype.title():s}")
        
        if self.debug:
            pars.show()
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        
        Parameters
        ----------

        silent : boolean (default: True)
            Flag for printing debugging information.
            
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True
        
        soma = self.soma
        if self.status['species'] == 'mouse':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            self.i_test_range = {'pulse': [(-0.3, 0.3, 0.03), (-0.05, 0., 0.005)]}  # set range for ic command test
            self.vrange = [-75., -55.]
            self._valid_temperatures = (34., )
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.c_m = 0.9
            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=1.0)
            soma().kht.gbar = self.g_convert(self.pars.kht_gbar, self.pars.units, self.somaarea)
            soma().klt.gbar = self.g_convert(self.pars.klt_gbar, self.pars.units, self.somaarea)
            soma().ihvcn.gbar = self.g_convert(self.pars.ihvcn_gbar, self.pars.units, self.somaarea)
            soma().ihvcn.eh = self.pars.ih_eh # Rodrigues and Oertel, 2006
            soma().leak.gbar = self.g_convert(self.pars.leak_gbar, self.pars.units, self.somaarea)
            soma().leak.erev = self.pars.leak_erev
            self.e_k = self.pars.e_k
            self.e_na = self.pars.e_na
            soma.ena = self.e_na
            soma.ek = self.e_k
           # soma().leak.erev = pars.leak_erev
            self.axonsf = 0.5
            
        elif self.status['species'] == 'guineapig':  # values from R&M 2003, Type II-I
            self.i_test_range = {'pulse': [(-0.3, 0.3, 0.03), (-0.05, 0., 0.005)]}  # set range for ic command test
            self.vrange = [-75., -55.]
            self.c_m = 0.9
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                self.i_test_range={'pulse': (-0.3, 0.3, 0.03)}
            self.vrange = [-75., -55.]
            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = self.g_convert(self.pars.kht_gbar, self.pars.units, self.somaarea)
            soma().klt.gbar = self.g_convert(self.pars.klt_gbar, self.pars.units, self.somaarea)
            #soma().ka.gbar = self.g_convert(self.pars.ka_gbar, self.somaarea)
            soma().ihvcn.gbar = self.g_convert(self.pars.ih_gbar, self.pars.units, self.somaarea)
            soma().leak.gbar = self.g_convert(self.pars.leak_gbar, self.pars.units, self.somaarea)
            soma().leak.erev = self.pars.leak_erev
            self.axonsf = 0.5

        else:
            raise ValueError(f"Species {self.status['species']} or type {self.status['modelType']}  is not recognized for {self.celltype:s} cells")
        self.check_temperature()

    def add_axon(self):
        """
        Add a default axon from the generic cell class to the bushy cell (see cell class).
        """
        Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)

    def add_dendrites(self):
        """
        Add simple unbranched dendrites to basic Rothman Type I-II model.
        The dendrites have some kht and ih current
        """
        cs = False  # not implemented outside here - internal Cesium.
        nDend = range(4) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 300 # length of the dendrite (not tapered)
            dendrites[i].diam = 1.25 # dendrite diameter
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('kht')
            if cs is False:
                dendrites[i]().kht.gbar = 0.005 # a little Ht
            else:
                dendrites[i]().kht.gbar = 0.0
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.0001
            dendrites[i].insert('ihvcn') # some H current
            dendrites[i]().ihvcn.gbar = 0.# 0.001
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


class DummyDStellate(DStellate):
    """ DStellate class with no cell body; this cell only replays a predetermined
    spike train. Useful for testing, or replacing spike trains to determine
    the importance of spike structures within a network.
    """
    def __init__(self, cf=None, species='mouse'):
        """
        Parameters
        ----------
        cf : float (default: None)
            Required: the characteristic frequency for the DStellate
            Really just for reference.

        """

        DStellate.__init__(self)
        self.vecstim = h.VecStim()
        
        # this causes the terminal to receive events from the VecStim:
        self.spike_source = self.vecstim
        
        # just an empty section for holding the terminal
        self.add_section(h.Section(), self.somaname)
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': None, 'species': species, 'modelType': 'Dummy', 'modelName': 'DummyDStellate',
                       'ttx': None, 'name': 'DummyDStellate',
                       'morphology': None, 'decorator': None, 'temperature': None}
        print("<< DStellate: Dummy DStellate Cell created >>")
        

    def set_spiketrain(self, times):
        """ Set the times of spikes (in seconds) to be replayed by the cell.
        """
        self._spiketrain = times
        self._stvec = h.Vector(times)
        self.vecstim.play(self._stvec)


class DStellateEager(DStellate):
    """
    This is a model of the VCN D-Stellate cells as proposed by
    Eager, M.A., Grayden, D.B., Burkitt, A.N., and Meffin, H.,
    "A neural circuit model of the ventral cochlear nucleus",
    Internet:
    http://citeseerx.ist.pus.edu/viewdoc/download?doi=10.1.79.9620.pdf&rep
    =rep&type=pdf
    also cited as:
    Proceedings of the 10th Australian International Conference on
    Speech Science and Technology, pp. 539-544, 2004.
    It is based on the Rothman and Manis (2003c) model,
    with small modifications.
    Their model includes dendrites and an axon, which are added in this version
    """
    def __init__(self, nach='na', ttx=False, species='guineapig', modelType='I-II', debug=False):
        """
        Initialize the VCN D-stellate model of Eager et al. Some model parameters may be modified.
        
        Parameters
        ----------
        nach : string (default: 'na')
            Set the sodium channel model. Choices are 'na', 'nav11', 'jsrna'
        
        ttx : boolean (default: False)
            ttx sets the sodium channel conductance to 0
        
        species : string (default: 'guineapig')
            species to use for conductance scaling
        
        modelType : string (default: 'I-II')
            RM03 model type to use for conductances. 
        
        debug : boolean (default: False)
            Flag to use to enable print statements for debugging purposes.
        
        """
        super(DStellateEager, self).__init__()

        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'DStellateEager'}
        self.i_test_range=(-0.25, 0.25, 0.025)  # set range for ic command test

        soma = h.Section(name="DStellateEager_Soma_%x" % id(self)) # one compartment

        soma.nseg = 1

        if nach in ['nacn', 'na']:
            soma.insert('na')
        elif nach == 'nav11':
            soma.insert('nav11')
        elif nach == 'jsrna':
            soma.insert('jsrna')
        else:
            raise ValueError('Sodium channel %s in type 1 cell not known' % nach)
        self.debug = debug
        soma.insert("kht")
        soma.insert('klt')
        soma.insert('ihvcn')
        soma.insert('leak')
        soma.ek = self.e_k
        soma().leak.erev = self.e_leak
        self.mechanisms = ['kht', 'klt', 'ihvcn', 'leak', nach]
        self.add_section(soma, self.somaname)
        self.species_scaling(silent=False, species=species, modelType=modelType)  # set the default type II cell parameters
        self.add_axon()  # must follow species scaling so that area parameters are available
        self.add_dendrites()   # similar for dendrites
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(soma)

        if self.debug:
                print("<< D-stellateEager: Eager DStellate Type I-II cell model created >>")

    def species_scaling(self, species='guineapig', modelType='I-II', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            A string specifying the species used for scaling. Recognized values are
            'mouse', 'guineapig', and 'cat' (cat is just a larger version of the guineapig)
        
        modelType : string (default: 'I-II')
            A string specifying the version of the model to use. 
            Current choices are 'I-II' (others need to be implemented)
        
        silent : boolean (default: True)
            Flag for printing debugging information.
            
        """
        soma = self.soma
        if species == 'mouse' and modelType == 'I-II':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma, gbar=800.)
            soma().kht.gbar = self.g_convert(150.0, 'nS', self.somaarea)
            soma().klt.gbar = self.g_convert(20.0, 'nS', self.somaarea)
            soma().ihvcn.gbar = self.g_convert(2.0, 'nS', self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = self.g_convert(2.0, 'nS', self.somaarea)
            self.axonsf = 0.5
        elif species == 'guineapig' and modelType == 'I-II':  # values from R&M 2003, Type II-I
            self.set_soma_size_from_Diam(25.0)
            self.adjust_na_chans(soma, gbar=1000.*0.75)
            soma().kht.gbar = 0.02  # self.g_convert(150.0, 'nS', self.somaarea)
            soma().klt.gbar = 0.005  #self.g_convert(20.0, 'nS', self.somaarea)
            soma().ihvcn.gbar = 0.0002  #self.g_convert(2.0, 'nS', self.somaarea)
            soma().leak.gbar = 0.0005  # self.g_convert(2.0,'nS',  self.somaarea)
            self.axonsf = 1.0
        elif species == 'cat' and modelType == 'I=II':  # a cat is a big guinea pig Type I
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = self.g_convert(150.0, 'nS', self.somaarea)
            soma().klt.gbar = self.g_convert(20.0, 'nS', self.somaarea)
            soma().ihvcn.gbar = self.g_convert(2.0, 'nS', self.somaarea)
            soma().leak.gbar = self.g_convert(2.0, 'nS', self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species %s or species-type %s is not recognized for D-StellateEager cells' %  (species, type))
        self.status['species'] = species
        self.status['type'] = modelType
        self.cell_initialize(showinfo=True)
        if not silent:
            print(' set cell as: ', species)
            print(' with Vm rest = %6.3f' % self.vm0)

    def adjust_na_chans(self, soma, gbar=1000.):
        """
        adjust the sodium channel conductance
        
        Parameters
        ----------
        soma : soma object (no default)
            soma object whose sodium channel complement will have it's
            conductances adjusted depending on the channel type
        
        gbar : float (default: 1000.)
            The conductance to be set for the sodium channel
        
        Returns
        -------
            Nothing
        """
        
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = self.g_convert(gbar, self.somaarea)
        nach = self.status['na']
        if nach == 'jsrna':
            soma().jsrna.gbar = gnabar
            soma.ena = self.e_na
            if self.debug:
                print('using jsrna with gbar: ', soma().jsrna.gbar)
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar * 0.5
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if self.debug:
                print("using inva11 with gbar:", soma().na.gbar)
            print('nav11 gbar: ', soma().nav11.gbar)
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if self.debug:
                print('using na with gbar: ', soma().na.gbar)
        elif nach == 'nach':
            soma().nach.gbar = gnabar
            soma.ena = self.e_na
            if self.debug:
                print(('uwing nacn with gbar: ', soma().nacn.gbar))
        else:
            raise ValueError("DstellateEager setting Na channels: channel %s not known" % nach)
        #print soma().na.gbar

    def add_axon(self):
        """
        Adds an axon to the Eager. et al model
        Cell.add_axon(self, nodes=1, c_m=self.c_m, R_a=self.R_a, axonsf=self.axonsf, dia=3.0, len=70, seg=2)
        The Eager et al model just uses one cable, 70 microns long and 3 microns in dameter.
        
        Parameters
        ----------
            None
        
        Returns
        -------
            Nothing
        """
        
        naxons = 1
        axon = []
        for i in range(naxons):
            axon.append(h.Section(cell=self.soma))
        for i in range(naxons):
            axon[i].connect(self.soma)
            axon[i].L = 70
            axon[i].diam = 3.0
            axon[i].Ra = 500
            axon[i].cm = 0.9
            axon[i].nseg = 2
            axon[i].insert("kht")
            axon[i].insert('klt')
            axon[i].insert('ihvcn')
            axon[i].insert('leak')
            axon[i].insert('na')
            axon[i].ek = self.e_k
            axon[i].ena = self.e_na
            axon[i]().leak.erev = self.e_leak
            axon[i]().na.gbar = 0.5
            axon[i]().klt.gbar = 0.005
            axon[i]().kht.gbar = 0.02
            axon[i]().ihvcn.gbar = 0.0002
            axon[i]().leak.gbar = 0.0005
        self.status['axon'] = True
        self.add_section(axon, 'axon')

    def add_dendrites(self):
        """
        Adds dendrites to the Eager model. The Eager model uses simple passive dendrites.
        
        Parameters
        ----------
            None
        
        Returns
        -------
            Nothing
        """
        
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 1100 # length of the dendrite (not tapered)
            dendrites[i].diam = 3.5 # dendrite diameter
            dendrites[i].nseg = 5 # # segments in dendrites
            dendrites[i].Ra = 1500 # ohm.cm
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00025
            dendrites[i]().leak.erev = self.e_leak
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


