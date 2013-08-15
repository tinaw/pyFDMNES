# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 18:22:28 2013

@author: weiget
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.pardir))
import pyFDMNES


sim = pyFDMNES.pyFDMNES("tio2_dafs_py_inp.txt")

sim.extract = True

sim.Reflex = np.array([
                      [0,0,0, 1,1,  45.],
                      [0,0,0, 1,2,  45.],
                      [0,0,1, 1,1,  45.],
                      [0,0,1, 1,2,  45.],
                      [1,1,1, 1,1,  45.],
                      [1,1,1, 1,2,  45.]
                      ])

sim.Filout("tio2_dafs_py_inp2.txt")

sim.FDMNESfile()

sim.retrieve()

dafs = sim.get_dafs([0,0,0], 1,2,45.)
plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.legend()
plt.figure(0)
plt.plot (dafs[:,0], dafs[:,(sim.column[0])], label = "%s"%sim.index)
plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.legend()
plt.figure(1)
plt.plot(dafs[:,0], dafs[:,(sim.column[1])], label = "%s"%sim.index)
plt.show()

"""dafs = sim.get_dafs([0,0,0], conv = False)
y = self.
plt.plot (dafs[:,0], dafs[:,y], label = "")
plt.legend()
plt.show()"""

