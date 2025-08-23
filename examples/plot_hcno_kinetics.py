#!/usr/bin/python
"""
test and plot hcno kinetics

"""
import sys
import numpy as np

# import matplotlib.pyplot as plt
import pyqtgraph as pg

# Parameters#
gbar = 0.0005  # (mho/cm2)

vhalf1 = -50  # (mV)        # v 1/2 for forward
vhalf2 = -84  # (mV)        # v 1/2 for backward
gm1 = 0.3  ## (mV)           # slope for forward
gm2 = 0.6  #    # (mV)      # slope for backward
zeta1 = 3  #   (/ms)
#    zeta2   = 3 #   (/ms)
a01 = 0.008  # (/ms)
a02 = 0.0029  # (/ms)
frac = 0.0
c0 = 273.16  # (degC)
thinf = -66  # (mV)        # inact inf slope
qinf = 7  # (mV)        # inact inf slope
q10tau = 4.5  # from Magee (1998)
# v       # (mV)
q10g = 4.5  # Cao and Oertel
celsius = 35.0

F = 9.648e4
R = 8.314

ct = 1e-3 * zeta1 * F / (R * (c0 + celsius))

q10 = q10tau ** (
    (celsius - 33.0) / 10.0
)  # (degC)) : if you don't like room temp, it can be changed!


def rates(v: np.ndarray):  # (mV)) {
    tau1 = bet1(v) / (q10 * a01 * (1.0 + alp1(v)))
    tau2 = bet2(v) / (q10 * a02 * (1.0 + alp2(v)))
    hinf = 1.0 / (1.0 + np.exp((v - thinf) / qinf))
    return tau1, tau2, hinf


def alp1(v: np.ndarray):  # (mV)) {
    alp1 = np.exp((v - vhalf1) * ct)
    return alp1


def bet1(v: np.ndarray):  # (mV)) {
    bet1 = np.exp(gm1 * (v - vhalf1) * ct)
    return bet1


def alp2(v: np.ndarray):  # (mV)) {
    alp2 = np.exp((v - vhalf2) * ct)
    return alp2


def bet2(v:np.ndarray):  # (mV)) {
    bet2 = np.exp(gm2 * (v - vhalf2) * ct)
    return bet2


def plots():

    v = np.linspace(-120.0, 0.0, 120)
    t1, t2, hi = rates(v)

    win = pg.GraphicsLayoutWidget()
    win.setBackground("w")
    p1 = win.addPlot(title="tau\sub1, tau\sub2", row=0, col=0, labels={"bottom": "V (mV)", "left": "tau1"})
    p1.plot(v, t1, pen=pg.mkPen("k", width=0.75))
    p1.plot(v, t2, pen=pg.mkPen("r", width=0.75))
    p3 = win.addPlot(title="h", row=2, col=0, labels={"bottom": "V (mV)", "left": "h"})
    p3.plot(v, hi, pen=pg.mkPen("k", width=0.75))
    return win

# prevent sphinx from running this script at first build
if __name__ == "__main__":
    win = plots()
    win.show()
    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()
