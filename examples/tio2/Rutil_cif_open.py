# -*- coding: utf-8 -*-
"""
Created on Thu Aug 08 15:44:09 2013

@author: weiget
"""
import os
import sys
sys.path.append("../../")

import pyFDMNES
import matplotlib.pyplot as plt
#import numpy as np

sim = pyFDMNES.pyFDMNES("TiO2.cif")

#sim.Range = np.array([-50, 0.2, 50., 0.5, 100., 2.0, 200])
sim.P.Radius = 2.0
sim.Rpotmax = 8.50

sim.Quadrupole = False
sim.Density = True
sim.Green = True
#sim.Self-absorption = True

#sim.Atom["Ti1"] = [3, 3,2,0., 4,0,2., 4,1,2.]
#sim.Atom["O1"] = [2, 2,0,2., 2,1,4.]

#sim.SCF = True
#sim.N_self = 5

sim.convolution = False
sim.cartesian =False
#sim.verbose = True
#sim.extract = True


sim.FileOut("tio2_dafs_py_inp.txt", overwrite=True)

#sim.do_convolution("tio2_dafs_py_inp_conv.txt")

sim.FDMNESfile()

sim.retrieve()

if sim.convolution:
    data = sim.get_XANES()
    plt.plot (data[:,0], data[:,1], label = "with Convolution")

data = sim.get_XANES(conv = False)
plt.plot (data[:,0], data[:,1], label = "without Convolution")

plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
plt.legend()
plt.show()



sim.do_convolution("tio2_dafs_py_inp_conv.txt")

sim.FDMNESfile()

sim.retrieve()


data = sim.get_XANES()
plt.plot (data[:,0], data[:,1], label = "with Convolution")

data = sim.get_XANES(conv = False)
plt.plot (data[:,0], data[:,1], label = "without Convolution")

plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
plt.legend()
plt.show()
