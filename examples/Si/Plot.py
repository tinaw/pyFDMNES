# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:14:41 2013

@author: weiget
"""
import os
import sys
import fdmnes
import matplotlib.pyplot as plt
#import numpy as np

sim = fdmnes.fdmnes("Si.cif")

#sim.Range = np.array([-50, 0.2, 50., 0.5, 100., 2.0, 200])
sim.radius = 2.0
sim.Rpotmax = 8.50

sim.Quadrupole = False
sim.Density = True
sim.Green = True
#sim.Self-absorption = True


sim.P.Atom["Ti1"] = [3, 3,2,0., 4,0,2., 4,1,2.]
sim.P.Atom["O1"] = [2, 2,0,2., 2,1,4.]

#sim.SCF = True
#sim.N_self = 5

#sim.convolution = True
sim.P.cartesian = False
sim.verbose = True
#sim.extract = True

sim.FileOut("Si_inp.txt", overwrite=True)

sim.do_convolution("Si_conv.txt")

sim.FDMNESfile()

sim.retrieve()

data = sim.get_XANES()
plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.legend()
plt.plot (data[:,0], data[:,1], label = "with Convolution")
plt.show()

data = sim.get_XANES(conv = False)
plt.plot (data[:,0], data[:,1], label = "without Convolution")
plt.legend()
plt.show()