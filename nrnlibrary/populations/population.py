import scipy.stats

"""
Todo: 

* Need to work out most of the API details here, probably by starting with a
  specific use case and working backward

* Distributions of cell properties

* Mechanisms for automatically recording from neurons
    - record all Vm for all real neurons
    - record spike trains
    - record per-synapse currents
    



"""

class Population(object):
    """
    A Population represents a group of cell all having the same type. 

    Populations provide methods for:
    
    * Adding cells to the population with characteristic distributions.
    * Connecting the cells in one population to the cells in another.
    * Automatically adding cells to satisfy connectivity requirements when 
      connecting populations together.
    
    Populations have a concept of a "natural" underlying distribution of
    neurons, and behave as if all neurons in this distribution already exist
    in the model. However, initially all neurons are virtual, and are only 
    instantiated to become a part of the running model if the neuron provides 
    synaptic input to another non-virtual neuron, or if the user explicitly 
    requests a recording of the neuron.
    
    """
    def __init__(self, species, size, fields):
        self._species = species
        self._connections = []
        
        # numpy record array with information about each cell in the 
        # population
        fields = [('cell', object), ('inputs_resolved', bool)] + fields
        self._cells = np.zeros(size, dtype=fields)

    @property
    def cells(self):
        """ The array of cells in this population. 
        
        For all populations, this array has a 'cell' field that is either 0
        (for virtual cells) or a Cell instance (for real cells). 
        
        Extra fields may be added by each Population subclass.
        """
        return self._cells.copy()
    
    @property
    def species(self):
        return self._species
    
    def unresolved_cells(self):
        """ Return indexes of all real cells whose inputs have not been 
        resolved.
        """
        real = self._cells['cell'] != 0
        unresolved = self._cells['input_resolved'] == False
        return np.argwhere(real & unresolved)[0]

    def connect(self, *pops):
        """ Connect this population to any number of other populations. 
        
        A connection is unidirectional; calling ``pop1.connect(pop2)`` can only
        result in projections from pop1 to pop2.
        
        Note that the connection is purely symbolic at first; no cells are 
        actually connected by synapses at this time.
        """
        self._connections.extend(pops)

    @property
    def connections(self):
        return self._connections[:]

    def resolve_inputs(self):
        """ For each _real_ cell in the population, select a set of 
        presynaptic partners from each connected population and generate a 
        synapse from each.
        
        Although it is allowed to call ``resolve_inputs`` multiple times for
        a single population, each individual cell will only resolve its inputs
        once. Therefore, it is recommended to create and connect all 
        populations before making any calls to ``resolve_inputs``.
        
        """
        raise NotImplementedError()
    
    def select(self, size, create=False, **kwds):
        """ Return a list of indexes for cells matching the selection criteria.
        
        The *size* argument specifies the number of cells to return.
        
        If *create* is True, then any selected cells that are virtual will be
        instantiated.
        
        Each keyword argument must be the name of a field in self.cells. Values
        may be either a number, in which case the cell with the closest match 
        is returned, or a distribution (see scipy.stats), in which case random
        values will be selected from the distribution.
        """
        if len(kwds) == 0:
            raise TypeError("Must specify at least one selection criteria")
        if len(kwds) > 1:
            raise NotImplementedError("Multiple selection criteria not yet 
                supported.")
        
        field, values = list(kwds.items())[0]
        if isinstance(values, scipy.stats.distributions.rv_frozen):
            values = values.rvs(size=size)
        elif np.isscalar(values):
            values = [values]
        
        cells = []
        mask = np.ones(self._cells.shape, dtype=bool)
        for val in values:
            err = np.abs(self._cells[field] - val)
            cell = np.argmin(err[mask])
            mask[cell] = False
            cells.append(cell)
            
        return cells

    def get_cell(self, i):
        """ Return the cell at index i. If the cell is virtual, then it will 
        be instantiated first.
        """
        if self._cells[i]['cell'] == 0:
            self.create_cells([i])
        return self._cells[i]['cell']
        
    def create_cells(self, cells):
        """ Instantiate each cell in *cells*, which is a list of indexes into
        self.cells.
        
        Subclasses must implement this method.
        """
        raise NotImplementedError()
