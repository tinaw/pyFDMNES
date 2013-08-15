# -*- coding: utf-8 -*-
"""
Created on Thu Aug 08 15:44:09 2013

@author: weiget
"""
import os
import sys
sys.path.append(os.path.abspath(os.path.pardir))

import pyFDMNES
import matplotlib.pyplot as plt

sim = pyFDMNES.pyFDMNES("TiO2.cif")

#sim.Range ='-50 0.2 50. 0.5 100. 2.0 200'
sim.radius = 2.0
sim.Rpotmax = 8.50

sim.Quadrupole = False
sim.Density = True
sim.Green = True
#sim.Self-absorption = True

sim.Atom["Ti1"] = [3, 3,2,0., 4,0,2., 4,1,2.]
sim.Atom["O1"] = [2, 2,0,2., 2,1,4.]


#sim.SCF = True
#sim.N_self = 5

sim.cartesian =False
sim.verbose = True
#sim.extract = True

sim.Filout("tio2_dafs_py_inp.txt")

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

