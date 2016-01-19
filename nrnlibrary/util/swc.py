# -*- coding: utf8 -*-
import numpy as np
from neuron import h


class SWC(object):
    """Encapsulates a morphology tree as defined by the SWC standard.
    
    Parameters
    ----------
    filename : str or None
        The name of an swc file to load
    types : dict or None
        A dictionary mapping {type_id: type_name} that describes the type IDs
        in the swc data (second column).
    data : ndarray or None
        Optionally, a data array may be provided instead of an swc file. This
        is used internally.
    """
    def __init__(self, filename=None, types=None, data=None):
        self._dtype = [
            ('id', int), 
            ('type', int), 
            ('x', float), 
            ('y', float), 
            ('z', float), 
            ('r', float), 
            ('parent', int)
        ]
        
        self._id_lookup = None
        self._sections = None
        self._children = None
        self.filename = None
        self.log = []

        self.sectypes = {
            # SWC-defined node types
            0: 'undefined',
            1: 'soma',
            2: 'axon',
            3: 'basal_dendrite',
            4: 'apical_dendrite',
            
            # Custom node types
            7: 'main_dendrite', 
            8: 'secondary_dendrite',
            9: 'axon_hillock',
            10: 'initial_segment',
            11: 'axon_node',
            12: 'internode',
        }
        if types is not None:
            self.sectypes.update(types)
        
        if data is not None:
            self.data = data
        elif filename is not None:
            self.load(filename)
        else:
            raise TypeError("Must initialize with filename or data array.")
        
        self.sort()
        
    def load(self, filename):
        self.filename = filename
        self.data = np.loadtxt(filename, dtype=self._dtype)
        self.log.append("Loaded %d nodes from %s" % (len(self.data), filename))
        
    def copy(self):
        s = SWC(data=self.data.copy(), types=self.sectypes)
        s.filename = self.filename
        s.log = self.log[:]
        return s

    @property
    def lookup(self):
        """Return a dict that maps *id* to *index* in the data array.
        """
        if self._id_lookup is None:
            self._id_lookup = dict([(rec['id'], i) for i, rec in enumerate(self.data)])
            #self._id_lookup = {}
            #for i, rec in enumerate(self.data):
                #self._id_lookup[rec['id']] = i
        return self._id_lookup

    def children(self, id):
        """Return a list of all children of the node *id*.
        """
        if self._children is None:
            self._children = {}
            for rec in self.data:
                ch = self._children.setdefault(rec['parent'], [])
                ch.append(rec['id'])
        return self._children.get(id, [])

    def __getitem__(self, id):
        """Return record for node *id*.
        """
        return self.data[self.lookup[id]]

    def reparent(self, id):
        """Rearrange tree to make *id* the new root parent.
        """
        d = self.data
        
        # bail out if this is already the root
        if self[id]['parent'] == -1:
            return
        
        self.log.append("Reparented tree using %d as new root node" % id)
        parent = -1
        while id != -1:
            oldparent = self[id]['parent']
            self[id]['parent'] = parent
            parent = id
            id = oldparent
            
        self._children = None
        self.sort()
        
    @property
    def sections(self):
        """Return lists of IDs grouped by topological section.
        
        The first item in each list connects to the last item in a previous
        list.
        """
        if self._sections is None:
            sections = []
            sec = []
            
            # find all nodes with nore than 1 child
            branchpts = set()
            endpoints = set(self.data['id'])
            endpoints.add(-1)
            seen = set()
            for r in self.data:
                p = r['parent']
                if p in seen:
                    branchpts.add(p)
                else:
                    seen.add(p)
                    endpoints.remove(p)
            
            # build lists of unbranched node chains
            for r in self.data:
                sec.append(r['id'])
                if r['id'] in branchpts or r['id'] in endpoints:
                    sections.append(sec)
                    sec = []
            
            self._sections = sections
            
        return self._sections
        
    def connect(self, parent_id, swc):
        """Combine this tree with another by attaching the root of *swc* as a 
        child of *parent_id*.
        """
        data = swc.data.copy()
        shift = self.data['id'].max() + 1 - data['id'].min()
        data['id'] += shift
        rootmask = data['parent'] == -1
        data['parent'] += shift
        data['parent'][rootmask] = parent_id
        
        self.data = np.concatenate([self.data, data])
        self._children = None
        self.sort()
        if len(swc.log) == 0:
            self.log.append('Connected new branch (with no log data)')
        else:
            self.log.append('Connected new branch with this log:')
            for line in swc.log:
                self.log.append('|   ' + line)
        
    def set_type(self, typ):
        self.data['type'] = typ
        self.log.append('Set type ID for all nodes to %d' % typ)
        
    def write_hoc(self, filename, types=None):
        """Write data to a HOC file.
        
        Each node type is written to a separate section list.
        """
        hoc = []
        
        if len(self.log) == 0:
            hoc.append('// HOC generated from SWC (with no log data)')
        else:
            hoc.append('// HOC generated from SWC:')
            for line in self.log:
                hoc.append('// ' + line)
        
        sectypes = self.sectypes.copy()
        for t in np.unique(self.data['type']):
            if t not in sectypes:
                sectypes[t] = 'type_%d' % t

        # create section lists
        for t in sectypes.values():
            hoc.extend(['objref %s' % t,
                        '%s = new SectionList()' % t])
        hoc.append('')
            
        # create sections
        sects = self.sections
        hoc.append('create sections[%d]' % len(sects))
        sec_ids = {}
        for i, sec in enumerate(sects):
            # remember hoc index for this section
            endpt = self[sec[-1]]['id']
            sec_id = len(sec_ids)
            sec_ids[endpt] = sec_id
            
            # add section to list
            hoc.append('access sections[%d]' % sec_id)
            typ = self[sec[0]]['type']
            hoc.append('%s.append()' % sectypes[typ])
            
            # connect section to parent
            p = self[sec[0]]['parent']
            if p != -1:
                hoc.append('connect sections[%d](0), sections[%d](1)' % (sec_id, sec_ids[p]))

            # set up geometry for this section
            hoc.append('sections[%d] {' % sec_id)
            hoc.append('  pt3dclear()')
            for i, seg in enumerate(sec):
                rec = self[seg]
                hoc.append('  pt3dadd(%f, %f, %f, %f)' % (rec['x'], rec['y'], rec['z'], rec['r']*2))
                p = self[sec[0]]['parent']
                if i == 0 and p != -1:
                    # Don't let NEURON force our sections to appear attached
                    prec = self[p]
                    dx = rec['x'] + prec['x']
                    dy = rec['y'] + prec['y']
                    dz = rec['z'] + prec['z']
                    hoc.append('  pt3dstyle(1, %f, %f, %f)' % (dx, dy, dz))
            hoc.append('}')
            
            hoc.append('')
        
        open(filename, 'w').write('\n'.join(hoc))

    def write_swc(self, filename):
        """Write SWC data to *filename*.
        """
        with open(filename, 'w') as fh:
            for line in self.log:
                fh.write('# ' + line + '\n')
            
            fh.write('\n# Section types:\n')
            for k,v in self.sectypes.items():
                fh.write("# %d: %s\n" % (k, v))
            
            fh.write("\n# id  type  x   y   z   r   parent\n")
            for rec in self.data:
                fh.write('%d %d %f %f %f %f %d\n' % tuple(rec))

    def create_sections(self):
        """Create NEURON sections for this tree.
        """
        sectypes = self.sectypes.copy()
        for t in np.unique(self.data['type']):
            if t not in sectypes:
                sectypes[t] = 'type_%d' % t

        # Store sections by type and by ID as they are created.
        sec_by_type = dict([(t, []) for t in sectypes.values()])
        sec_by_id = {}
            
        # create sections
        for i, sec in enumerate(self.sections):
            # remember hoc index for this section
            endpt = self[sec[-1]]['id']
            
            # add section to set
            typ_id = self[sec[0]]['type']
            typ = sectypes[typ_id]
            name = '%s_%d' % (typ, len(sec_by_type[typ]))
            section = h.Section(name=name)
            sec_by_type[typ].append(section)
            sec_by_id[endpt] = section
            
            # connect section to parent
            pid = self[sec[0]]['parent']
            if pid != -1:
                psec = sec_by_id[pid]
                section.connect(psec(1.0), 0.0)

            # set up geometry for this section
            h.pt3dclear(sec=section)
            for i, seg in enumerate(sec):
                rec = self[seg]
                h.pt3dadd(rec['x'], rec['y'], rec['z'], rec['r']*2, sec=section)
                #p = self[sec[0]]['parent']
                #if i == 0 and p != -1:
                    ## Don't let NEURON force our sections to appear attached
                    #prec = self[p]
                    #dx = rec['x'] + prec['x']
                    #dy = rec['y'] + prec['y']
                    #dz = rec['z'] + prec['z']
                    #hoc.append('  pt3dstyle(1, %f, %f, %f)' % (dx, dy, dz))
                    
        return sec_by_id, sec_by_type

    @property
    def root(self):
        """ID of the root node of the tree.
        """
        ind = np.argwhere(self.data['parent'] == -1)[0, 0]
        return self.data[ind]['id']

    def sort(self):
        """Sort the tree in topological order.
        """
        order = self.branch(self.root)
        lt = self.lookup
        indexes = np.array([lt[i] for i in order], dtype=int)
        self.data = self.data[indexes]
        
        self._id_lookup = None
        self._sections = None
        
    def path(self, node):
        path = [node]
        while True:
            node = self[node]['parent']
            if node < 0:
                return path
            path.append(node)

    def scale(self, x, y, z, r):
        self.data['x'] *= x
        self.data['y'] *= y
        self.data['z'] *= z
        self.data['r'] *= r
        self.log.append('Applied scale factors: x=%f y=%f z=%f r=%f' % (x, y, z, r))
        
    def translate(self, x, y, z):
        self.data['x'] += x
        self.data['y'] += y
        self.data['z'] += z
        self.log.append('Applied translation: x=%f y=%f z=%f' % (x, y, z))
        
    def branch(self, id):
        """Return a list of IDs in the branch beginning at *id*.
        """
        branch = [id]
        for ch in self.children(id):
            branch.extend(self.branch(ch))
        return branch
    
    def topology(self):
        """Print the tree topology.
        """
        path = []
        indent = ''
        secparents = [self[s[0]]['parent'] for s in self.sections]
        
        for i, sec in enumerate(self.sections):
            p = secparents[i]
            if p != -1:
                ind = path.index(p)
                path = path[:ind+1]
                indent = indent[:(ind+1) * 3]
            path.append(self[sec[-1]]['id'])

            # look ahead to see whether subsequent sections are children
            if p in secparents[i+1:]:
                this_indent = indent[:-2] + u"├─ "
                indent =      indent[:-2] + u"│  │  "
            else:
                this_indent = indent[:-2] + u"└─ "
                indent =      indent[:-2] + u"   │  "
                
                
            typ = self.sectypes[self[sec[0]]['type']]
            if len(sec) > 10:
                secstr = "%s,...%s" % (str(tuple(sec[:3]))[:-1], str(tuple(sec[-3:]))[1:])
            else:
                secstr = str(tuple(sec))
            print "%ssections[%d] type=%s parent=%d %s" % (this_indent, i, typ, p, secstr)



if __name__ == '__main__':
    soma = SWC('data/cellbody.swc', types={1:'soma', 2:'axon', 3:'dendrite'})
    soma.set_type(1)
    soma.data['r'] *= 0.5  # this data was recorded as diameter
    axon = SWC('data/axonnonscaled.swc')
    axon.set_type(2)
    dend = SWC('data/dendnonscaled.swc')
    dend.set_type(7)
    dend.reparent(755)
    
    cell = soma.copy()
    cell.connect(57, axon)
    cell.connect(39, dend)
    # correct for voxel size
    cell.scale(0.11, 0.11, 0.06, 0.11)
    # correct for shrinkage
    s = 1.0 / 0.75
    cell.scale(s, s, s, s)
    cell.translate(-70, -90, -60)
    cell.write_hoc('test.hoc')
    cell.write_swc('test.swc')
    #soma.topology()
    
    #dend.data[748]['type'] = 2
    #dend.scale(0.11, 0.11, 0.06, 0.11)
    #dend.data['x'] -= 45
    #dend.data['y'] -= 60
    #dend.data['z'] -= 60
    #dend.write_hoc('test.hoc')
    

    #s = SWC('data/test.swc')
    #s.topology()
