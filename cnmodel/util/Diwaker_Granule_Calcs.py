"""Convert nS to mmho/cm2 given an area."""

import numpy as np
import math

pi = math.pi

# values taken from the ionchannels.py file, in mmho/cm2
soma_GRC_NA_gbar   =   27.122 # [0]  
hillock_GRC_NA_gbar =  238.0 # [0] 
soma_GRC_KV_gbar  =    3.0 # [1]   
soma_GRC_KA_gbar     = 3.201 # [1]   
soma_GRC_KM_gbar      = 0.250 # [2]   
soma_GRC_KIR_gbar     = 0.900 # [1] 
soma_GRC_KCA_gbar     = 3.000 # [3]   
soma_GRC_CA_gbar      = 0.460 # [4]   
soma_GRC_CALC_gbar    = 0.0095 # [1]
soma_GRC_LKG1_gl      = 0.107 # [1]
soma_GRC_LKG2_ggaba      = 0.040 # [1]

gs = {
    "GRC_NA_gbar": soma_GRC_NA_gbar,
    "GRC_NA_hillock_gbar": hillock_GRC_NA_gbar,
    "GRC_KV_gbar": soma_GRC_KV_gbar,
    "GRC_KA_gbar": soma_GRC_KA_gbar,
    "GRC_KM_gbar": soma_GRC_KM_gbar,
    "GRC_KIR_gbar": soma_GRC_KIR_gbar,
    "GRC_KCA_gbar": soma_GRC_KCA_gbar,
    "GRC_CA_gbar": soma_GRC_CA_gbar,
    "GRC_CALC_gbar": soma_GRC_CALC_gbar,
    "GRC_LKG1_gl": soma_GRC_LKG1_gl,
    "GRC_LKG2_ggaba": soma_GRC_LKG2_ggaba,
}


# routine to convert conductances from nS as given elsewhere
#   to mho/cm2 as required by NEURON 1/28/99 P. Manis
# units: nano siemens, soma area in um^2
#
def nstomho(ns, somaarea, refarea=None):
    if refarea == None:
        return 1e-9 * float(ns) / float(somaarea)
    else:
        return 1e9 * float(ns) / float(refarea)


def mho2ns(mho, somaarea):
    return float(mho) * somaarea / 1e-9


def spherearea_um_to_cm2(dia):
    """
    given diameter in microns, return sphere area in cm2
    """
    r = dia/2 * 1e-4  # convert to cm
    return 4 * pi * r**2


def spherearea_um_to_um2(dia):
    """
    given diameter in microns, return sphere area in um2
    """
    r = dia/2
    return 4 * pi * r**2

soma_area = spherearea_um_to_cm2(5.8)
print(f"soma area: {soma_area:.9f} cm2")
print(f"Soma cap: {1e-6 * soma_area:.4e} uF")
soma_area_um2 = spherearea_um_to_um2(5.8)
print(f"soma area: {soma_area_um2:.9f} um2")
# for ng in gs:
#     print(f"var: {ng}, value: {gs[ng]:6.3f} nS = {nstomho(gs[ng]*1e3, soma_area):.9f} mmho/cm2")


# values and calculations from the Diwaker model, as given in Parametri.hoc
soma_L = 5.8
soma_diam = 5.8
dend_1_L = 5
dend_1_diam = 0.75
dend_4_L = 2.5
dend_4_diam = 0.75
axon_L = 20
axon_diam = 0.3
n_axon = 30
n_dend = 4

# Compartment area estimation
SomaArea=soma_L * soma_diam*pi
Dend12Area=dend_1_L * dend_1_diam * pi
Dend34Area=dend_4_L * dend_4_diam * pi
SomascArea=pi*9.76*9.76

#Scale factors for compartaments
RappSomaDend12=SomascArea/(4*Dend12Area)
RappSomaDend34=SomascArea/(4*Dend34Area)
RappSomaNew=SomascArea/SomaArea
RappSomahill=SomascArea/(3.75*pi) 
RappAH = 3.75/(axon_L*axon_diam)
RappAxon = (9.76*9.76)/(n_axon*axon_L*axon_diam*pi)

gamma = 0.5

gbar_axon = soma_GRC_NA_gbar

axon_gNa = soma_GRC_NA_gbar * (1-gamma) * RappAxon-0.00232
hillock_gNa = hillock_GRC_NA_gbar * (gamma) * RappSomahill-0.00232

axon_gK = soma_GRC_KV_gbar * (1-gamma) * RappAxon-0.00232
hillock_gK = soma_GRC_KV_gbar * (gamma) * RappSomahill-0.00232

print(f"\naxon gNa: {axon_gNa:.6f} mmho/cm2 [vs {soma_GRC_NA_gbar} in the table]. Factor = axon_gNa/soma_GRC_NA_gbar = {axon_gNa/soma_GRC_NA_gbar:.3f}")
print(f"hillock gNa: {hillock_gNa:.6f} mmho/cm2 [vs {hillock_GRC_NA_gbar} in the table]. Factor = hillock_gNa/hillock_GRC_NA_gbar = {hillock_gNa/hillock_GRC_NA_gbar:.3f}"
      )
print("Hillock NA factor relative to reference axon: ", hillock_gNa/axon_gNa, " Reference to base soma: ", hillock_gNa/soma_GRC_NA_gbar)

print(f"\naxon gK: {axon_gK:.6f} mmho/cm2 [vs {soma_GRC_KV_gbar} in the table]. Factor = axon_gK/soma_GRC_KV_gbar = {axon_gK/soma_GRC_KV_gbar:.3f}")
print(f"hillock gK: {hillock_gK:.6f} mmho/cm2 [vs {soma_GRC_KV_gbar} in the table]. Factor = hillock_gK/soma_GRC_KV_gbar = {hillock_gK/soma_GRC_KV_gbar:.3f}"
      )
print("Hillock K factor relative to reference axon: ", hillock_gK/axon_gK, " Reference to base soma: ", hillock_gK/soma_GRC_KV_gbar)

# Leakage
soma_gl = soma_GRC_LKG1_gl * RappSomaNew* (2/3.)  # why 2/3? 
print(f"\nSoma gl: {soma_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG1_gl} in the table]. Factor = soma_gl/soma_GRC_LKG1_gl = {soma_gl/soma_GRC_LKG1_gl:.3f}")
axon_gl = soma_GRC_LKG1_gl * RappAxon * (1/30.)
print(f"\nAxon gl: {axon_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG1_gl} in the table]. Factor = axon_gl/soma_GRC_LKG1_gl = {axon_gl/soma_GRC_LKG1_gl:.3f}")
print("Axon gl factor relative to reference soma: ", axon_gl/soma_GRC_LKG1_gl, " Reference to base soma: ", axon_gl/soma_GRC_LKG1_gl)

hillock_gl = soma_GRC_LKG1_gl * RappSomahill * (1/15.)
print(f"\nHillock gl: {hillock_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG1_gl} in the table]. Factor = hillock_gl/soma_GRC_LKG1_gl = {hillock_gl/soma_GRC_LKG1_gl:.3f}")
print("Hillock gl factor relative to reference soma: ", hillock_gl/soma_GRC_LKG1_gl, " Reference to base soma: ", hillock_gl/soma_GRC_LKG1_gl)

dend_proximal_gl = soma_GRC_LKG1_gl * RappSomaDend12 * (1/16.)
print(f"\nDendrite proximal gl: {dend_proximal_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG1_gl} in the table]. Factor = dend_proximal_gl/soma_GRC_LKG1_gl = {dend_proximal_gl/soma_GRC_LKG1_gl:.3f}")
print(f"Dendrite proximal gl factor relative to reference soma: {dend_proximal_gl/soma_GRC_LKG1_gl:.3f} Reference to base soma: {dend_proximal_gl/soma_GRC_LKG1_gl:.3f}")

dend_distal_gl = soma_GRC_LKG1_gl * RappSomaDend34 * (1/16.)
print(f"\nDendrite distal LKG1: {dend_distal_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG1_gl} in the table]. Factor = dend_distal_gl/soma_GRC_LKG1_gl = {dend_distal_gl/soma_GRC_LKG1_gl:.3f}")
print(f"Dendrite distal LKG1 factor relative to reference soma: {dend_distal_gl/soma_GRC_LKG1_gl:.3f} Reference to base soma: {dend_distal_gl/soma_GRC_LKG1_gl:.3f}")


# GABA leak

dend_proximal_gl = soma_ggaba * (1/n_dend) * (RappSomaDend12)
print(f"\nDendrite proximal LKG2: {dend_proximal_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG2_ggaba} in the table]. Factor = dend_proximal_gl/soma_GRC_LKG2_ggaba = {dend_proximal_gl/soma_GRC_LKG2_ggaba:.3f}")
print(f"Dendrite proximal LKG2 factor relative to reference soma: {dend_proximal_gl/soma_GRC_LKG2_ggaba:.3f} Reference to base soma: {dend_proximal_gl/soma_GRC_LKG2_ggaba:.3f}")

dend_distal_gl = soma_GRC_LKG2_ggaba * (1/n_dend) * RappSomaDend34
print(f"\nDendrite distal LKG2: {dend_distal_gl:.6f} mmho/cm2 [vs {soma_GRC_LKG2_ggaba} in the table]. Factor = dend_distal_gl/soma_GRC_LKG2_ggaba = {dend_distal_gl/soma_GRC_LKG2_ggaba:.3f}")
print(f"Dendrite distal LKG2 factor relative to reference soma: {dend_distal_gl/soma_GRC_LKG2_ggaba:.3f} Reference to base soma: {dend_distal_gl/soma_GRC_LKG2_ggaba:.3f}")