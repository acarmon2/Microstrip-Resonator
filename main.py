#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################################################################################
# main function to calculate the generators
import numpy as np
import time
from MicrostripTech import *

###################################################################################################
# Main variables of the transmission line (SI units)
W = 1.230e-3                                                                # Width of the conductor
H = 0.508e-3                                                                # Height of the substrate
T = 17.5e-6                                                                 # Metallization profile
tand = 0.0035                                                               # Loss tangent for ROS4350B
er = 3.48                                                                   # Dielectric constant
Delta = 0.4e-6                                                              # Rugosity profile of Copper electrodeposited for RO4350B
sigma = 5.983e7                                                             # Conductivity of the Copper
c = 299792548.0                                                             # Light velocity in vacuum
mu0 = 4*np.pi*1e-7                                                          # Magnetic permeability of vacuum
e0 = 8.8541878176e-12                                                       # Electric permitivity of vacuum
f0 = 2.45e9                                                                 # Frequency of interest in X Band
L0 = c/f0                                                                   # Wavelength in the vacuum

# Main characteristics of LT
ef = Effective(er, W, H)                                                    # Effective dielectric constant
Gap = 0.25e-3                                                               # Gap of the resonator and feed lines
Z0 = Impedance(ef, W, H)                                                    # Characteristic impedance
Ceq = CapGap(er, G, W, H)                                                   # Both capacitances required to find the extra length of the tranmission line
C1 = Ceq[0]                                                                 # Capacitance between end line and ground
C2 = Ceq[1]                                                                 # Capacitance between couple of end lines
Lr = (c/(2*f0*sqrt(ef))) - (2*(c*Z0*(C1+C2)/np.raiz(ef)))                   # Length "real" of the resonator


Z0 = Impedance(ef, W, H)                                                    # Characteristic impedance
Cl = Capacitance(er, W, H)                                                  # Capacitance per length (F/m)
Il = Inductance(W, H)                                                       # Inductance per length (H/m)
Rl = Resistance(W, H, T, f0, mu0, sigma)                                    # Frequency depend line resistance (Resistance per length)(Ohm/m)                                                                 % Total resistance of the resonator section
alphaRT = ConductorLoss(ef, W, H, T, Delta, f0, Z0, mu0, sigma)             # Total conductor losses (Conductor + rugosity sections) per length (dB/m)
Qc = (27.3*sqrt(ef))/(alphaRT*L0)
alphaDT = DielectricLoss(L0, er, ef, tand)                                  # Dielectric total losses per length unit (dB/m)
Qd = (27.3*sqrt(ef))/(alphaDT*L0)
Qu = (((1/Qc)+(1/Qd))**(-1))                                                # Limit of material properties

