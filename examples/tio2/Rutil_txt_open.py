# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 18:22:28 2013

@author: weiget
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import fdmnes

sim = fdmnes.fdmnes("tio2_dafs_py_inp.txt")

sim.WriteInputFile("tio2_daf_2_py_inp.txt", overwrite=True)
sim.Run(wait=True)

sim.Status()
sim.DoConvolution(overwrite=True)
sim.Status()

dafs = sim.get_dafs((0,0,2), 1,1,45., conv=True)


sim.WriteInputFile("tio2_dafs_3_py_inp.txt", overwrite=True)
sim.Run()
sim.Status()
sim.DoConvolution(overwrite=True)

plt.title("RXS %s" %sim.index)
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.plot (dafs[:,0], dafs[:,(sim.column)], label = sim.index)
plt.legend(loc = 2)
plt.show()
