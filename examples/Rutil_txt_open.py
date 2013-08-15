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

dafs = sim.get_dafs((0,0,2), 1,2,45.,conv=False)
"""
plt.title("DAFS  %s"%sim.index)
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.plot (dafs[:,0], dafs[:,(sim.column_real)], label = sim.index_real)
plt.plot (dafs[:,0], dafs[:,(sim.column_im)], label = sim.index_im)
plt.legend(loc = 3)
plt.show()


plt.title("DAFS %s"%sim.index)
plt.xlabel("Energy")
plt.ylabel("Intensity")
y = sim.index
plt.plot (dafs[:,0], dafs[:,(sim.column)], label = sim.index)
plt.legend(loc = 2)
plt.show()
"""

