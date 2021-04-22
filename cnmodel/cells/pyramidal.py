from __future__ import print_function
from neuron import h

import numpy as np
from .cell import Cell
from ..util import Params
from .. import data

__all__ = ['Pyramidal', 'PyramidalKanold', 'PyramidalCeballos']

class Pyramidal(Cell):

    celltype = 'pyramidal'
    scaled = False
    
    @classmethod
    def create(cls, model='POK', **kwds):
        if model == 'POK': 
            return PyramidalKanold(**kwds)
        if model == 'Ceballos':
            return PyramidalCeballos(**kwds)
        else:
            raise ValueError ("Pyramidal model '%s' is unknown", model)

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
        else:
            loc = 0.5
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
                if 'AMPAScale' in kwds:
                    self.AMPA_gmax = self.AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDA_gmax = self.NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.celltype == 'dstellate':  # WBI input -Voigt, Nelken, Young
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            elif terminal.cell.celltype == 'tuberculoventral':  # TV cells talk to each other-Kuo et al.
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

class PyramidalKanold(Pyramidal, Cell):
    """
    DCN pyramidal cell
    Kanold and Manis, 1999, 2001, 2005
    """
    def __init__(self,  morphology=None, decorator=None, nach=None,
                 ttx=False, species='rat', modelType=None, modelName=None,
                 debug=False, temperature=None):
        """
        initialize a pyramidal cell, based on the Kanold-Manis (2001) pyramidal cell model.
        Modifications to the cell can be made by calling methods below. These include
        converting to a model with modified size and conductances (experimental).
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanim
            by that name must exist. None implies the default channel, 'napyr'.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'rat')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored (overridden by decorator).

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point models.
            
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        
        """
        super(PyramidalKanold, self).__init__()
        if modelType == None or modelType == 'I':
            modelName = 'POK'
            modelType = 'pyramidal'
            dataset = 'POK_channels'
            temp = 32.  # as defined in the ms. and in the hoc files on modelDB

        else:
            raise ValueError(f"Species {species:s} and modeltype {modelType:s} not recognized for {self.celltype:s} cells")

        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName, 'ttx': ttx, 'name': 'Pyramidal',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None,
                   }
        self.debug=debug
        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp

        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'leak']
            for mech in self.mechanisms:
                try:
                    self.soma.insert(mech)
                except ValueError:
                    print('WARNING: Mechanism %s not found' % mech)
            self.soma().kif.kif_ivh = -89.6
            self.species_scaling(silent=True)  # set the default type I-c  cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if debug:
            print("<< PYR: POK Pyramidal Cell created >>")


    def get_cellpars(self, dataset, species='guineapig', modelType='II'):
        pars = self.get_initial_pars(dataset, species, modelType)

        for g in ['soma_napyr_gbar', 'soma_kdpyr_gbar', 'soma_kif_gbar', 'soma_kis_gbar',
                  'soma_ihpyr_gbar', 'soma_leak_gbar',
                  'soma_e_h','soma_leak_erev', 'soma_e_k', 'soma_e_na']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        if self.debug:
            pars.show()
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be 'rat'
        
        modelType: string (default: 'I')
            definition of model type from Kanold and Manis, 2001
            choices are 'I' or 'POK' (canonical model) or
            'II', a modified model with more physiological surface area and KCNQ channels
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True
            
        soma = self.soma
        if self.status['species'] in ['rat', 'mouse']:
            print('for species: ', self.status['species'])
            if self.status['modelType'] not in ['pyramidal']:  # canonical K&M2001 model cell
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            
            self._valid_temperatures = (32.,)
            if self.status['temperature'] is None:
              self.set_temperature(32.)
            self.i_test_range = {'pulse': (-0.3, 0.401, 0.02)}
            self.vrange = [-75., -60.]
            self.set_soma_size_from_Cm(self.pars.cap)
            soma().napyr.gbar = self.g_convert(self.pars.soma_napyr_gbar, self.pars.units, self.somaarea)
            # soma().nap.gbar = self.g_convert(self.pars.soma_nap_gbar, self.somaarea) # does not exist in canonical model
            soma().kdpyr.gbar = self.g_convert(self.pars.soma_kdpyr_gbar, self.pars.units, self.somaarea)
            # soma().kcnq.gbar = self.g_convert(self.pars.soma_kcnq_gbar, self.somaarea) # does not exist in canonical model.
            # soma().kpksk.gbar = self.g_convert(self.pars.soma_kpksk_gbar, self.somaarea) # does not exist in canonical model.
            # soma().kir.gbar = self.g_convert(self.pars.soma_kir_gbar, self.somaarea)
            soma().kif.gbar = self.g_convert(self.pars.soma_kif_gbar, self.pars.units, self.somaarea)
            soma().kis.gbar = self.g_convert(self.pars.soma_kis_gbar, self.pars.units, self.somaarea)
            soma().ihpyr.gbar = self.g_convert(self.pars.soma_ihpyr_gbar, self.pars.units, self.somaarea)
#            soma().ihpyr_adj.q10 = 3.0  # no temp scaling to sta
            soma().leak.gbar = self.g_convert(self.pars.soma_leak_gbar, self.pars.units, self.somaarea)
            
            # Reversal potentials
            soma().leak.erev = self.pars.soma_leak_erev
            soma().ena = self.pars.soma_e_na
            soma().ek = self.pars.soma_e_k
            soma().ihpyr.eh = self.pars.soma_e_h

        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for T-stellate cells")

        self.check_temperature()
        if not silent:
            print(f"Set cell as: {self.status['species']:s}, {self.status['modelType']:s}")
            print(self.status)
            for m in self.mechanisms:
                print('%s.gbar = %f' % (m, eval('soma().%s.gbar' % m)))

    def add_dendrites(self):
        """
        Add simple unbranched dendrite.
        The dendrites have some kd, kif and ih current
        """
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 250 # length of the dendrite (not tapered)
            dendrites[i].diam = 1
            dendrites[i].cm = self.c_m
            #h('dendrites[i].diam(0:1) = 2:1') # dendrite diameter, with tapering
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('napyr')
            dendrites[i]().napyr.gbar = 0.00
            dendrites[i].insert('kdpyr')
            dendrites[i]().kdpyr.gbar = 0.002 # a little Ht
            dendrites[i].insert('kif')
            dendrites[i]().kif.gbar = 0.0001 # a little Ht
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00001
            dendrites[i].insert('ihpyr') # some H current
            dendrites[i]().ihvcn.gbar =  0. # 0.00002
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


class PyramidalCeballos(Pyramidal, Cell):
    """
    DCN pyramidal cell
    Ceballos et al., Front Cell Neurosci. 2016
    (based on POK cell, above)
    """
    def __init__(self,  morphology=None, decorator=None, nach=None,
                 ttx=False, species='mouse', modelType=None, modelName=None,
                 debug=False, temperature=None):
        """
        initialize a pyramidal cell, based on the Kanold-Manis (2001) pyramidal cell model.
        Modifications to the cell can be made by calling methods below. These include
        converting to a model with modified size and conductances (experimental).
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanim
            by that name must exist. None implies the default channel, 'napyr'.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'rat')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored (overridden by decorator).

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point models.
            
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        
        """
        super(PyramidalCeballos, self).__init__()

        dataset = 'Ceballos_channels'
        temp = 34.
        print(modelType, modelName)
        if modelType in ['quiet', 'active']:
            modelName = 'Ceballos'
        else:
            raise ValueError(f"Species {species:s} and modeltype {modelType:s} not recognized for {self.celltype:s} cells")
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName, 'ttx': ttx, 'name': 'Pyramidal',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None,
                   }
        self.debug=debug
        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp
        self.set_cm(1.0)  # original Ceballos model uses this value
        
        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['napyr', 'nappyr', 'kdpyr', 'kif', 'kis', 'kcnq', 'kir', 'ihpyrlc', 'leak']
            for mech in self.mechanisms:
                try:
                    self.soma.insert(mech)
                except:
                    raise ValueError('WARNING: Mechanism %s not found' % mech)
            self.soma().kif.kif_ivh = -89.6
            self.species_scaling(silent=True)  # set the default type I-c  cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if debug:
            print("<< PYR: Ceballos Pyramidal Cell created >>")


    def get_cellpars(self, dataset, species='mouse', modelType='I'):
        pars = self.get_initial_pars(dataset, species, modelType)

        for g in ['soma_napyr_gbar', 'soma_nappyr_gbar', 
                  'soma_kdpyr_gbar', 'soma_kif_gbar', 'soma_kis_gbar',
                  'soma_kcnq_gbar', 'soma_kir_gbar', 'soma_ihpyrlc_gbar', 
                  'soma_leak_gbar',
                  'soma_e_h','soma_leak_erev', 'soma_e_k', 'soma_e_na']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        if self.debug:
            pars.show()
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be 'rat'
        
        modelType: string (default: 'I')
            definition of model type from Kanold and Manis, 2001
            choices are 'I' or 'POK' (canonical model) or
            'II', a modified model with more physiological surface area and KCNQ channels
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True

        soma = self.soma
        if self.status['species'] in ['rat', 'mouse']:
            if self.status['modelType'] not in ['quiet', 'active']:  # canonical K&M2001 model cell
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
              self.set_temperature(34.)
            self.i_test_range = {'pulse': (-0.3, 0.401, 0.02)}
            self.vrange = [-75., -56.]  # tricky as there is a negative slope above about -58 mV - bistability
            self.set_soma_size_from_Cm(self.pars.cap)
            soma().napyr.gbar = self.g_convert(self.pars.soma_napyr_gbar, self.pars.units, self.somaarea)
            soma().nappyr.gbar = self.g_convert(self.pars.soma_nappyr_gbar, self.pars.units, self.somaarea) # does not exist in canonical model
            soma().kdpyr.gbar = self.g_convert(self.pars.soma_kdpyr_gbar, self.pars.units, self.somaarea)
            soma().kif.gbar = self.g_convert(self.pars.soma_kif_gbar, self.pars.units, self.somaarea)
            soma().kis.gbar = self.g_convert(self.pars.soma_kis_gbar, self.pars.units, self.somaarea)
            soma().kcnq.gbar = self.g_convert(self.pars.soma_kcnq_gbar, self.pars.units, self.somaarea) # does not exist in canonical model.
            # soma().kpksk.gbar = self.g_convert(self.pars.soma_kpksk_gbar, self.somaarea) # does not exist in canonical model.
            soma().kir.gbar = self.g_convert(self.pars.soma_kir_gbar,self.pars.units,  self.somaarea)
            soma().ihpyrlc.gbar = self.g_convert(self.pars.soma_ihpyrlc_gbar, self.pars.units, self.somaarea)
#            soma().ihpyr_adj.q10 = 3.0  # no temp scaling to sta
            soma().leak.gbar = self.g_convert(self.pars.soma_leak_gbar,self.pars.units,  self.somaarea)
          
            # Reversal Potentials for various ions
            soma().leak.erev = self.pars.soma_leak_erev
            soma().ena = self.pars.soma_e_na
            soma().ek = self.pars.soma_e_k
            soma().ihpyrlc.eh = self.pars.soma_e_h

            print("pyr gbar napyr: ", self.pars.soma_napyr_gbar)
            print("somaarea: ", self.somaarea)
            print("pyr gbar", soma().napyr.gbar)            

        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for T-stellate cells")

        self.check_temperature()
        if not silent:
            print(f"Set cell as: {self.status['species']:s}, {self.status['modelType']:s}")
            print(self.status)
            for m in self.mechanisms:
                print('%s.gbar = %f' % (m, eval('soma().%s.gbar' % m)))

