
#!/usr/bin/python

"""
Test the basic membrane physiology of cell types.

Basic Usage:   python test_cells.py celltype species [--cc | --vc]

This script generates a cell of the specified type and species, then tests the
cell with a series of current/voltage pulses to produce I/V, F/I, and spike
latency analyses.

"""

import argparse
import os, sys
import numpy as np
from pathlib import Path
from neuron import h
import pyqtgraph as pg

import cnmodel.cells as cells
import cnmodel


from cnmodel.protocols import IVCurve, VCCurve



debugFlag = True
ax = None
h.celsius = 22
cclamp = False

cellinfo = {'types': ['bushy', 'bushycoop', 'tstellate', 'tstellatenav11', 'dstellate', 'dstellateeager', 'sgc',
                      'cartwheel', 'pyramidal', 'pyramidalceballos', 'octopus', 'tuberculoventral', 'mso',
                      'granule'],
            'morphology': ['point', 'waxon', 'stick'],
            'nav': ['std', 'jsrna', 'nav11', 'nacncoop'],
            'species': ['guineapig', 'cat', 'rat', 'mouse'],
            'pulse': ['step', 'pulse']}

# Format for ivranges is list of tuples. This allows finer increments in selected ranges, such as close to rest
ccivrange = {'mouse':
                {'bushy': {'pulse': [(-1, 1.2, 0.05)]},
                 'bushycoop': {'pulse': [(-0.5, 0.7, 0.02)]},
                 'tstellate': {'pulse': [(-1., 1.01, 0.05), (-0.015, 0, 0.005)]},
                 'tstellatenav11': {'pulse': [(-1, 1., 0.1)]},
                 'dstellate': {'pulse': [(-0.3, 0.301, 0.015)]},
                 'octopus': {'pulse': [(-1.0, 1.0, 0.05)]},
                 'sgc': {'pulse': [(-0.3, 0.6, 0.02)]},
                 'cartwheel': {'pulse': [(-0.5, 0.5, 0.05)]},
                 'pyramidal': {'pulse': [(-0.3, 0.3, 0.025), (-0.040, 0.025, 0.005)]},
                 'pyramidalceballos': {'pulse': [(-0.09, 0.00, 0.09), (0, 0.008, 0.008)]}, #, 'prepulse': [(-0.25, -0.25, 0.25)]},
                 'tuberculoventral': {'pulse': [(-0.35, 1.0, 0.05), (-0.040, 0.01, 0.005)]},
                 'granule': {'pulse': [(-0.02, 0.02, 0.002)]}
                },

            'guineapig':
                {'bushy': {'pulse': [(-1, 1.2, 0.05)]},
                'tstellate': {'pulse': [(-0.15, 0.15, 0.01)]},
                'dstellate': {'pulse': [(-0.25, 0.25, 0.025)]},
                'dstellateeager': {'pulse': [(-0.6, 1.0, 0.025)]},
                'octopus': {'pulse': [(-2., 6., 0.2)]},
                'sgc': {'pulse': [(-0.3, 0.6, 0.02)]},
                'mso': {'pulse': [(-1, 1.2, 0.05)]},
                },
            'rat':
                {'pyramidal': {'pulse': [(-0.3, 0.3, 0.025), (-0.040, 0.025, 0.005)]}, # 'prepulse': [(-0.25, -0.25, 0.25)]},
                }
            }

# scales holds some default scaling to use in the cciv plots
# argument is {cellname: (xmin, xmax, IVymin, IVymax, FIspikemax,
# offset(for spikes), crossing (for IV) )}
## the "offset" refers to setting the axes back a bit
scale = {'bushy': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
        'bushycoop': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                'crossing', [0, -60]),
        'tstellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'tstellatenav11': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'tstellatedend': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'dstellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'dstellateeager': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'sgc:': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'cartwheel': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'pyramidal': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'pyramidalceballos': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'tuberculoventral': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'granule': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'octopus': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
            'crossing', [0, -60]),
        'mso': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                        'crossing', [0, -60]),
        }

class Tests():
    """
    Class to select cells for tests
    """
    def __init__(self):
        pass

    def selectCell(self, args):
        """
        Parameters
        ----------
        args : argparse args from command line
        
        Returns
        -------
        cell
            Instantiated cell of the selected celltype
        """
        h.celsius = float(args.temp)
        #
        # Spiral Ganglion cell tests
        #
        if args.celltype == 'sgc':  # morphology is always "point" for SGCs
            cell = cells.SGC.create(debug=debugFlag, species=args.species,
                nach=args.nav, ttx=args.ttx, modelType=args.type)
        #
        # Bushy tests
        #
        elif args.celltype == 'bushy' and args.morphology == 'point':
            cell = cells.Bushy.create(model='RM03', species=args.species, modelType=args.type,
                ttx=args.ttx, nach=args.nav, debug=debugFlag)
#            cell.soma().klt.gbar = 0.0003

        elif args.celltype == 'bushy' and args.morphology == 'waxon':
            cell = cells.Bushy.create(model='RM03', species=args.species, modelType=args.type,
                nach=args.nav, ttx=args.ttx, debug=debugFlag)
            cell.add_axon()

        elif args.celltype == 'bushy' and args.morphology == 'stick':
            cell = cells.Bushy.create(model='RM03', species=args.species, modelType=args.type, 
                morphology='cnmodel/morphology/bushy_stick.hoc', decorator=True,
                nach=args.nav, ttx=args.ttx, debug=debugFlag)
            h.topology()

        elif args.celltype == 'bushycoop' and args.morphology == 'point':
            cell = cells.Bushy.create(model='RM03', species=args.species, modelType=args.type,
                ttx=args.ttx, nach=args.nav, debug=debugFlag)

        #
        # T-stellate tests
        #
        elif args.celltype == 'tstellate' and args.morphology == 'point':
             cell = cells.TStellate.create(model='RM03', species=args.species, modelType=args.type,
                 nach=args.nav, ttx=args.ttx, debug=debugFlag)

        elif args.celltype == 'tstellate' and args.morphology == 'stick':
             cell = cells.TStellate.create(model='RM03', species=args.species, modelType=args.type,
                 nach=args.nav, ttx=args.ttx, debug=debugFlag, 
                 morphology='cnmodel/morphology/tstellate_stick.hoc', decorator=True)

        elif args.celltype == 'tstellatenav11' and args.morphology == 'point':  # note this uses a different model...
            print('test_cells: Stellate NAV11')
            cell = cells.TStellateNav11.create(model='Nav11', species=args.species, modelType=None,
                ttx=args.ttx, debug=debugFlag)

        elif args.celltype == 'tstellatenav11' and args.morphology == 'stick':  # note this uses a different model...
            cell = cells.TStellateNav11.create(model='Nav11', species=args.species, modelType=None, 
                morphology='cnmodel/morphology/tstellate_stick.hoc', decorator=True, 
                ttx=args.ttx, debug=debugFlag,)
            h.topology()

        #
        # Octopus cell tests
        #
        elif args.celltype == 'octopus' and args.morphology == 'point':
            cell = cells.Octopus.create(species=args.species, modelType='RM03', # args.type,
                nach=args.nav, ttx=args.ttx, debug=debugFlag)

        elif args.celltype == 'octopus' and args.morphology == 'stick':  # Go to spencer et al. model
            cell = cells.Octopus.create(modelType='Spencer', species=args.species,
                morphology='cnmodel/morphology/octopus_spencer_stick.hoc', decorator=True,
                nach=args.nav, ttx=args.ttx, debug=debugFlag)
            h.topology()

        #
        # D-stellate tests
        #
        elif args.celltype == 'dstellate':
            cell = cells.DStellate.create(debug=debugFlag, species=args.species, ttx=args.ttx, modelType=args.type)

        elif args.celltype == 'dstellateeager':
            cell = cells.DStellateEager.create(debug=debugFlag, ttx=args.ttx, modelType=args.type)

        #
        # DCN pyramidal cell tests
        #
        elif args.celltype == 'pyramidal':
            cell = cells.Pyramidal.create(model='POK', modelType=args.type, 
                ttx=args.ttx, debug=debugFlag)

        elif args.celltype == 'pyramidalceballos':
            cell = cells.Pyramidal.create(modelType=args.type, model='Ceballos',
                ttx=args.ttx, debug=debugFlag)
            cell.vm0 = -65.0 # cell may not have a stable state, so force rmp

        #
        # DCN tuberculoventral cell tests
        #
        elif args.celltype == 'tuberculoventral' and args.morphology == 'point':
            cell = cells.Tuberculoventral.create(species='mouse', modelType='TVmouse',
                     ttx=args.ttx, nach=args.nav, debug=debugFlag)

        elif args.celltype == 'tuberculoventral' and args.morphology == 'stick':
            cell = cells.Tuberculoventral.create(species='mouse', modelType='TVmouse', 
                    morphology='cnmodel/morphology/tv_stick.hoc', decorator=True,
                    ttx=args.ttx, debug=debugFlag)
            h.topology()

       #
        # DCN granule cell tests
        #
        elif args.celltype == 'granule' and args.morphology == 'point':
            cell = cells.Granule.create(model="GRC", species='mouse', modelType='GRC',
                     ttx=args.ttx, nach=args.nav, debug=debugFlag)

        elif args.celltype == 'granule' and args.morphology == 'stick':
            grc_stick = Path('cnmodel/morphology/granule_stick_simple.hoc')
            print("find morphology: ", grc_stick.is_file())
            cell = cells.Granule.create(species=args.species, 
                modelName="GRC_channels_compartments", modelType='GRC', 
                morphology=str(grc_stick), decorator=True,
                ttx=args.ttx, debug=debugFlag)
            h.topology()
        # elif args.celltype == 'granule' and args.morphology == 'stick':
        #     grc = Path('cnmodel/morphology/granule_stick_simple.hoc')
        #     print("find morphology: ", grc.is_file())
        #     cell = cells.Granule.create(species='mouse', modelType='GRC_channels_compartment', 
        #             morphology='cnmodel/morphology/granule_stick_simple.hoc', decorator=True,
        #             ttx=args.ttx, debug=debugFlag)
        #     h.topology()

        #
        # DCN cartwheel cell tests
        #
        elif args.celltype == 'cartwheel':
            cell = cells.Cartwheel.create(modelType=args.type, ttx=args.ttx, debug=debugFlag)

        #
        # MSO principal neuron tests
        #
        elif args.celltype == 'mso' and args.morphology == 'point':
            print("MSO Creation")
            cell = cells.MSO.create(model='MSO-principal', species=args.species, modelType=args.type,
                ttx=args.ttx, nach=args.nav, debug=debugFlag)    

        else:
            raise ValueError ("Cell Type %s and configurations nav=%s or config=%s are not available" %
                 (args.celltype, args.nav, args.morphology))

        print(cell.__doc__)
        self.cell = cell
        

    def run_test(self, sites, ptype, args):
        """
        Run either vc or cc test, and plot the result
        
        Parameters
        ----------
        args : argparse args from command line
        
        """
        self.cell.set_temperature(float(args.temp))
        print("Cell status dictionary:\n", self.cell.status)
        print("Resting potential: ", self.cell.vm0)
        durations = eval(args.durations)
        if self.cell.vm0 is None:
            V0 = self.cell.find_i0(showinfo=True, vrange=[-85, -55])
        else:
            V0 = self.cell.vm0
#        self.cell.cell_initialize()
        print('Currents at nominal Vrest= %.2f I = 0: I = %g ' % (V0, self.cell.i_currents(V=V0)))
        self.cell.print_mechs(self.cell.soma)
        instant = self.cell.compute_rmrintau(auto_initialize=False, vrange=None)
        print('    From Inst: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}'.format(instant['Rin'], instant['tau'], instant['v']))
        if args.cc is True:
            # define the current clamp electrode and default settings
            self.iv = IVCurve()
            self.iv.run(ccivrange[args.species][args.celltype], self.cell, durs=durations, 
                   sites=sites, reppulse=ptype, temp=float(args.temp))
            ret = self.iv.input_resistance_tau()
            if not np.isnan(ret['slope']):
                print('    From IV: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}'.format(ret['slope'], ret['tau'], ret['intercept']))
            self.iv.show(cell=self.cell)

        elif args.rmp is True:
            print('temperature: ', self.cell.status['temperature'])
            self.iv = IVCurve()
            self.iv.run({'pulse': (0, 0, 1)}, self.cell, durs=durations,
                   sites=sites, reppulse=ptype, temp=float(args.temp))
            self.iv.show(cell=self.cell, rmponly=True)

        elif args.vc is True:
            # define the voltage clamp electrode and default settings
            self.vc = VCCurve()
            self.vc.run((-120, 40, 10), self.cell)
            self.vc.show(cell=self.cell)

        else:
            raise ValueError("Nothing to run. Specify one of --cc, --vc, --rmp.")

def main():
    parser = argparse.ArgumentParser(description=('test_cells.py:',
    ' Biophysical representations of neurons (mostly auditory), test file'))
    parser.add_argument('celltype', action='store')
    parser.add_argument('species', action='store', default='guineapig')
    parser.add_argument('--type', action='store', default=None)
    parser.add_argument('--temp', action='store', default=22.0,
                        help=("Temp DegC (22 default)"))
    parser.add_argument('-m', action="store", dest="morphology",
        default='point', help=("Set morphology: %s " %
            [morph for morph in cellinfo['morphology']]))
    parser.add_argument('--nav', action="store", dest="nav", default=None,
        help=("Choose sodium channel: %s " % [ch for ch in cellinfo['nav']]))
    parser.add_argument('--ttx', action="store_true", dest="ttx", default=False,
        help=("Use TTX (no sodium current"))
    parser.add_argument('-p', action="store", dest="pulsetype", default="step",
        help=("Set CCIV pulse to step or repeated pulse"))
    parser.add_argument('-d', '--durations', action="store", dest="durations",
        default='[10, 100, 50]',
        help=("Set pulse durations in format '[10, 100, 20]' (as a string)"))
    clampgroup = parser.add_mutually_exclusive_group()
    clampgroup.add_argument('--vc', action='store_true',
        help="Run in voltage clamp mode")
    clampgroup.add_argument('--cc', action='store_true',
        help="Run in current clamp mode")
    clampgroup.add_argument('--rmp', action='store_true',
        help="Run to get RMP in current clamp mode")

    
    args = parser.parse_args()
    try:
        eval(args.durations)
    except:
        raise ValueError("Durations values could not be parsed\nFor example, use: '[10,100,10]' in quotes")
        
    if args.celltype not in cellinfo['types']:
        print('cell: %s is not in our list of cell types' % (args.celltype))
        print('celltypes: ', cellinfo['types'])
        sys.exit(1)

    path = os.path.dirname(cnmodel.__file__)
#    h.load_file("stdrun.hoc")
#    h.load_file(os.path.join(path, "custom_init.hoc")) # replace init with one that gets closer to steady state

    # print 'Species: ', args.species
    # print 'Morphology configuration: ', args.morphology
    sites = None
    if args.pulsetype == 'step':
        ptype = None
    else:
        ptype = 'pulses'
    if args.morphology in cellinfo['morphology']:
        print('Morphological configuration %s is ok' % args.morphology)

    t = Tests()
    t.selectCell(args)
    app = pg.mkQApp()
    t.run_test(sites, ptype, args)
    print("run complete")
    app.exec()
    # if sys.flags.interactive == 0:
    #     pg.QtWidgets.QApplication.exec() 

if __name__ == '__main__':
    main()
