from neuron import h
# from ..util.pynrnutilities import nstomho
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
        # Diwakar HOC files use 'Soma' (capital S) as the
        # SectionList name; base class defaults somaname to 'soma'.  Override here
        # so both point and stick paths find the soma in all_sections['Soma'].
        self.somaname = 'Soma'
        self.i_test_range={'pulse': (-0.02, 0.02, 0.005)}  # note that this might get reset with decorator according to channels
                                                    # The default values are set in the species_scaling routine
        if species == 'mouse' and morphology is not None:
            if modelType == None or modelType == 'GRC':
                modelName = 'granule'
                modelType = 'GRC'
                dataset = 'GRC_channels'
                temp = 34.

            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for {self.celltype:s} cells with morphology {morphology:s}")

        elif species == 'mouse' and morphology is None:
            if modelType == None or modelType == 'GRC':
                modelName = 'granule'
                modelType = 'GRC'
                dataset = 'GRC_channels_singlecompartment'
                temp = 34.

            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for {self.celltype:s} cells with morphology {morphology:s}")
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        print("Decorator: ", decorator, "Morphology: ", morphology)
        # the table name is 'GRC_channels' or 'GRC_channels_singlecompartment', so modelName must be 'GRC' here.
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': 'GRC',
                       'ttx': ttx, 'name': self.celltype,
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
            # species_scaling sets _valid_temperatures, vrange,
            # and calls set_temperature, but is NOT called in the decorator path.
            # Initialize here so check_temperature() at end of channel_manager works.
            self._valid_temperatures = (34.,)
            self.set_temperature(34.)
            self.vrange = [-90., -30.]
            self.decorate()
            # GRC_CALC is a Ca2+ buffer required for Ca2+
            # dynamics; it has no gbar so the decorator cannot insert it.
            for sec in self.all_sections[self.somaname]:
                sec.insert('GRC_CALC')
            # reversal potentials are never applied in the
            # decorator path (species_scaling is not called). Set ena and ek from
            # the GRC_channels table. GRC_LKG1.el default is -16.5 mV (MOD file);
            # must be set to soma_e_leak (-75.0 mV).
            e_na = self.pars.soma_e_na    # 87.39 mV
            e_k = self.pars.soma_e_k      # -84.69 mV
            e_leak = self.pars.soma_e_leak  # -75.0 mV
            for sec_list in self.all_sections.values():
                for sec in sec_list:
                    try:
                        sec.ena = e_na
                    except (AttributeError, RuntimeError):
                        pass
                    try:
                        sec.ek = e_k
                    except (AttributeError, RuntimeError):
                        pass
                    try:
                        sec.GRC_LKG1.el = e_leak
                    except (AttributeError, RuntimeError):
                        pass
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)

        if debug:
            print("<< Granule: Diwakar et al. model created >>")
            print("    Using ionchannel table: ", dataset)
            if self.status['decorator'] is not None:
                self._print_conductance_table()

    def get_cellpars(self, dataset, species='mouse', modelType='GRC'):
        pars = self.get_initial_pars(dataset, species, modelType)

        for g in ['soma_GRC_NA_gbar', 'soma_GRC_KV_gbar', 'soma_GRC_KM_gbar', 'soma_GRC_KA_gbar',
                  'soma_GRC_KCA_gbar',
                  'soma_GRC_KIR_gbar', 'soma_GRC_CA_gbar',
                  'soma_GRC_LKG1_gl',  'soma_GRC_LKG2_ggaba',
                  'soma_e_k', 'soma_e_na', 'soma_e_ca', 'soma_e_leak',
                  'soma_Dia', 'soma_GRC_LKG1_erev', 'soma_GRC_LKG2_egaba',
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
        self.vrange = [-90., -30.] 
        
        # self.set_soma_size_from_Cm(self.pars.cap)
        # print("granule cell pars: ", dir(self.pars))
        self.set_soma_size_from_Diam(self.pars.soma_Dia)

        print("G units: ", self.pars.units)
        print("Soma area: ", self.somaarea, "cm2")


        if self.pars.units == 'nS':
            self.soma().GRC_NA.gbar = self.g_convert(self.pars.soma_GRC_NA_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_KV.gbar = self.g_convert(self.pars.soma_GRC_KV_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_KA.gbar = self.g_convert(self.pars.soma_GRC_KA_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_KIR.gbar = self.g_convert(self.pars.soma_GRC_KIR_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_KCA.gbar = self.g_convert(self.pars.soma_GRC_KCA_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_CA.gbar = self.g_convert(self.pars.soma_GRC_CA_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_KM.gbar = self.g_convert(self.pars.soma_GRC_KM_gbar, self.pars.units, self.somaarea)
            self.soma().GRC_LKG1.gl = self.g_convert(self.pars.soma_GRC_LKG1_gl, self.pars.units, self.somaarea) 
            self.soma().GRC_LKG2.ggaba = self.g_convert(self.pars.soma_GRC_LKG2_ggaba, self.pars.units, self.somaarea)  # Claude fixed 2026-07-10: GRC_LKG2 param is ggaba not gl
        elif self.pars.units == 'mmho/cm2':
            self.soma().GRC_NA.gbar = self.pars.soma_GRC_NA_gbar
            self.soma().GRC_KV.gbar = self.pars.soma_GRC_KV_gbar
            self.soma().GRC_KA.gbar = self.pars.soma_GRC_KA_gbar
            self.soma().GRC_KIR.gbar = self.pars.soma_GRC_KIR_gbar
            self.soma().GRC_KCA.gbar = self.pars.soma_GRC_KCA_gbar
            self.soma().GRC_CA.gbar = self.pars.soma_GRC_CA_gbar
            self.soma().GRC_KM.gbar = self.pars.soma_GRC_KM_gbar
            self.soma().GRC_LKG1.gl = self.pars.soma_GRC_LKG1_gl
            self.soma().GRC_LKG2.ggaba = self.pars.soma_GRC_LKG2_ggaba  # Claude fixed 2026-07-10: GRC_LKG2 param is ggaba not gl
        else:
            raise ValueError ('Granule,  species: only "nS" or "mmho/cm2" are recognized for units')
      
        self.soma().ena = self.pars.soma_e_na # 50
        self.soma().ek = self.pars.soma_e_k # -80
        self.soma().eca = self.pars.soma_e_ca # 50
        self.soma().GRC_LKG1.el = self.pars.soma_GRC_LKG1_erev # -16.5 mV
        self.soma().GRC_LKG2.egaba = self.pars.soma_GRC_LKG2_egaba # -65 mV
        
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
        h.finitialize(V)
        self.ix = {}

        if self.status['morphology'] is not None:
            # Claude 2026-07-09: morphological model — sum i_membrane_ (mA/cm²)
            # over all segments to get total ionic current in nA across all compartments.
            # use_fast_imem enables seg.i_membrane_ without requiring extracellular.
            # self.mechanisms is a dict {section_name: set} after decorate(); populate
            # self.ix with those keys so find_i0's completeness check passes.
            self.ix = {name: 0.0 for name in self.mechanisms}
            h.cvode.use_fast_imem(1)
            h.fcurrent()
            result = sum(seg.i_membrane_ * seg.area() * 1e-2
                         for sec in h.allsec() for seg in sec)
            h.cvode.use_fast_imem(0)
            return result

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
        #  decorator path inserts 'leak' (gbar=0); must appear in ix
        # or find_i0's completeness check raises ValueError.
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar * (V - self.soma().leak.erev)
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

    def _print_conductance_table(self):
        """
        Print two side-by-side diagnostic tables after decoration.

        Table 1 — channelMap: what channel_manager computed (decorator's intent).
        Table 2 — psection:   what NEURON actually has in each section.

        Comparing them identifies whether a mismatch is in channel_manager
        (wrong channelMap values) or in _biophys (decorator failed to set gbar).

        Values are in mS/cm2 with nS in parentheses.  Area uses NEURON's own
        seg.area() so tapered sections are correct.
        """
        cond_var_map = {'GRC_LKG1': 'gl', 'GRC_LKG2': 'ggaba'}

        # One representative section per compartment type (first = most proximal)
        comp_names = [n for n, secs in self.all_sections.items() if secs]
        rep_secs   = {n: self.all_sections[n][0] for n in comp_names}

        # --- Build mechanism list from channelMap (conductance parameters only) ---
        cond_mechs = []   # ordered list of mechanism names
        cond_param  = {}  # mech → channelMap parameter key (e.g. 'GRC_NA_gbar')
        has_map = hasattr(self, 'channelMap') and bool(self.channelMap)
        if has_map:
            all_params = set()
            for cdata in self.channelMap.values():
                all_params.update(cdata.keys())
            for param in sorted(all_params):
                if '_gbar' in param or '_gl' in param or '_ggaba' in param:  # Claude fixed 2026-07-10: include _ggaba for GRC_LKG2
                    mech = param.rsplit('_', 1)[0]
                    if mech not in cond_mechs:
                        cond_mechs.append(mech)
                        cond_param[mech] = param

        # Fall back to psection if no channelMap
        if not cond_mechs:
            for sec in rep_secs.values():
                for m in sec.psection().get('density_mechs', {}):
                    if m not in cond_mechs:
                        cond_mechs.append(m)

        if not cond_mechs:
            print("  *** No mechanisms found after decoration! ***")
            return

        # --- Formatting helpers ---
        lw = 20    # compartment name column
        aw =  8    # area column
        cw = 22    # per-mechanism column

        def _fmt(g_S, area_cm2):
            if g_S is None:
                return '---'
            return f"{g_S*1e3:.4f}({g_S*area_cm2*1e9:.4f})"

        def _header(title):
            print(f"\n  {title}")
            print(f"  {'Compartment':<{lw}} {'Area(um2)':>{aw}}  " +
                  ''.join(f"{m:>{cw}}" for m in cond_mechs))
            print(f"  {'':>{lw}} {'':>{aw}}  " +
                  ''.join(f"{'mS/cm2(nS)':>{cw}}" for _ in cond_mechs))
            print("  " + "-" * (lw + aw + 2 + cw * len(cond_mechs)))

        def _row(comp_name, area_um2, area_cm2, get_g):
            row = f"  {comp_name:<{lw}} {area_um2:>{aw}.2f}  "
            for mech in cond_mechs:
                row += f"{_fmt(get_g(mech), area_cm2):>{cw}}"
            print(row)

        # ─── Table 1: channelMap (intended by channel_manager) ────────────────
        if has_map:
            _header("=== [1] channelMap — intended conductances ===")
            for cn in comp_names:
                sec      = rep_secs[cn]
                area_um2 = sum(seg.area() for seg in sec)
                area_cm2 = area_um2 * 1e-8
                cdata    = self.channelMap.get(cn, {})
                _row(cn, area_um2, area_cm2,
                     lambda mech, cd=cdata: cd.get(cond_param.get(mech)))
            print()
        else:
            print("\n  (channelMap not available — point model path)")

        # ─── Table 2: psection (actual NEURON values) ─────────────────────────
        _header("=== [2] psection  — actual NEURON conductances ===")
        any_nonzero = False
        for cn in comp_names:
            sec      = rep_secs[cn]
            area_um2 = sum(seg.area() for seg in sec)
            area_cm2 = area_um2 * 1e-8
            dm       = sec.psection().get('density_mechs', {})
            def get_actual(mech, dm=dm):
                cv = cond_var_map.get(mech, 'gbar')
                if mech not in dm:
                    return None
                g = dm[mech].get(cv, [0.0])[0]
                return g
            _row(cn, area_um2, area_cm2, get_actual)
            for mech in cond_mechs:
                g = get_actual(mech)
                if g is not None and g > 0:
                    any_nonzero = True
        if not any_nonzero:
            print("  *** ALL conductances are zero — decoration failed to set gbar! ***")
        print()
