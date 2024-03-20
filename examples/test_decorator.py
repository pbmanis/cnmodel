#!/usr/bin/python
"""
Basic test of cnmodel decorator on simple cell
Uses LC_bushy.hoc and XM13 mechanisms.

"""

__author__ = "pbmanis"


import sys
import numpy as np
import timeit

import pyqtgraph as pg
import cnmodel.cells as cells
import cnmodel.decorator as Decorator
from cnmodel.util import pyqtgraphPlotHelpers as PH
from cnmodel.protocols import IVCurve


class test_decorator:
    def __init__(self, filename, celltype="bushy"):
        #
        self.filename = filename
        self.iv = IVCurve()  # use standard IVCurve here...
        self.temperature = 37
        self.initdelay = 20.0
        self.cell_type = celltype

    def run(self):
        if self.cell_type == "bushy":
            self.post_cell = cells.Bushy.create(
            morphology=self.filename,
            decorator=Decorator,
            species="mouse",
            modelName="XM13_nacncoop_compartments",
            modelType="II",
        )
        elif self.cell_type == "granule":
            self.post_cell = cells.Granule.create(
                morphology=self.filename,
                decorator=Decorator,
                species="mouse",
                modelName="granule",
                modelType="GRC",
            )
            for s in self.post_cell.all_sections:
                if s in ['axon', 'ais', 'hillock', 'soma', 'primarydendrites', 'preclaw', 'dendriticclaw']:
                    for i, p in enumerate(self.post_cell.all_sections[s]):
                        for dm in p.psection()['density_mechs']:
                            for m in p.psection()['density_mechs'][dm]:
                                if m == 'gbar':
                                    print(f"     {s:>16}  {i:2d} {dm:<8s} {m:<12s} {1e3*p.psection()['density_mechs'][dm][m][0]:>10.4f} mmho/cm2")



        self.post_cell.set_temperature(float(self.temperature))
        self.post_cell.set_nseg(
            freq=2000.0
        )  # necessary to ensure appropriate spatial discreetization
        self.iv.reset()
        irange = self.post_cell.i_test_range
        if self.cell_type == "bushy":
            irange = {"pulse": (-1.0, 1.2, 0.2)}  # for Figure 5 of paper
        elif self.cell_type == "granule":
            irange = {"pulse": (-0.02, 0.02, 0.002)}
        self.durs = [self.initdelay, 100.0, 50.0] # (self.initdelay + 20.0, 100.0, 50.0)
        ptype = None  # 'pulses'
        sites = None
        # self.iv.run(
        #     irange,
        #     self.post_cell,
        #     durs=self.durs,
        #     temp=float(self.temperature),
        #     initdelay=self.initdelay,
        # )
        self.iv.run(irange, self.post_cell, durs=self.durs, 
                   sites=sites, reppulse=ptype, temp=float(self.temperature))
        ret = self.iv.input_resistance_tau()
        if not np.isnan(ret['slope']):
            print('    From IV: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}'.format(ret['slope'], ret['tau'], ret['intercept']))
        self.iv.show(cell=self.post_cell)
        return 
        

    def plot(self):
        pg.setConfigOption("background", "w")  # set background to white
        pg.setConfigOption("foreground", "k")
        self.app = pg.mkQApp()

        wintitle = "test_decorator"
        self.view = pg.GraphicsView()
        self.layout = pg.GraphicsLayout()  # (border=(100,100,100))
        self.view.setCentralItem(self.layout)
        self.view.resize(800, 600)
        self.view.show()
        self.plots = {}
        nr1 = 6
        nc = 10
        for i in range(1, nr1):
            self.plots["p%d" % i] = None
        for i in range(1, nr1 + 1):
            self.layout.addLayout(row=i, col=nc)
        for i in range(1, nc + 1):
            self.layout.addLayout(row=nr1 + 2, col=i)

        self.plots["p1"] = self.layout.addPlot(
            row=1,
            col=1,
            rowspan=6,
            colspan=9,
            labels={"left": "V (mV)", "bottom": "Time (ms)"},
        )
        self.plots["p2"] = self.layout.addPlot(
            row=7,
            col=1,
            rowspan=1,
            colspan=9,
            labels={"left": "I (nA)", "bottom": "Time (ms)"},
        )

        for k in range(len(self.iv.voltage_traces)):
            self.plots["p1"].plot(
                self.iv.time_values,
                np.array(self.iv.voltage_traces[k]),
                pen=pg.mkPen("k", width=0.75),
            )
            self.plots["p2"].plot(
                self.iv.time_values,
                np.array(self.iv.current_traces[k]),
                pen=pg.mkPen("k", width=0.75),
            )
        self.plots["p1"].setRange(
            xRange=(0.0, np.sum(self.durs) - self.initdelay), yRange=(-160.0, 40.0)
        )
        self.plots["p2"].setRange(
            xRange=(0.0, np.sum(self.durs) - self.initdelay), yRange=(-1, 1)
        )
        PH.noaxes(self.plots["p1"])
        PH.calbar(
            self.plots["p1"],
            calbar=[125.0, -120.0, 10.0, 20.0],
            unitNames={"x": "ms", "y": "mV"},
        )
        PH.noaxes(self.plots["p2"])
        PH.calbar(
            self.plots["p2"],
            calbar=[125, 0.1, 0.0, 0.5],
            unitNames={"x": "ms", "y": "nA"},
        )

        text = "{0:2d}\u00b0C {1:.2f}-{2:.2f} nA".format(
            int(self.temperature),
            np.min(self.iv.current_cmd),
            np.max(self.iv.current_cmd),
        )
        ti = pg.TextItem(text, anchor=(1, 0))
        ti.setFont(pg.QtGui.QFont("Arial", 9))
        ti.setPos(120.0, -120.0)
        self.plots["p1"].addItem(ti)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'lc_bushy':
            decorator = test_decorator("examples/LC_bushy.hoc", celltype="bushy")
        elif sys.argv[1] == 'granule':
            decorator = test_decorator("cnmodel/morphology/grc_stick_simple.hoc", celltype="granule")
        else:
            print(f"Unknown cell type: {sys.argv[1]}")
            sys.exit(0)
    print(decorator)
    # else:
    #     fn = ('/Users/pbmanis/Desktop/Python/VCNModel/VCN_Cells/VCN_c{0:02d}/Morphology/VCN_c{0:02d}.hoc'.format(int(sys.argv[1])))
    #     fig5 = F5(fn)
    start_time = timeit.default_timer()
    decorator.run()
    elapsed = timeit.default_timer() - start_time
    print("Elapsed time for simulation: %f" % (elapsed))
    # decorator.plot()
    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()
