# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 18:22:28 2013

@author: weiget
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

import pyFDMNES

sim = pyFDMNES.pyFDMNES("tio2_dafs_py_inp.txt")



sim.Fileout("tio2_daf_2_py_inp.txt")

sim.FDMNESfile()

sim.retrieve()

sim.do_convolution()

sim.FDMNESfile()

sim.retrieve()

dafs = sim.get_dafs((0,0,2), 1,1,45., conv=True)
"""
plt.title("RXS  %s"%sim.index)
plt.xlabel("Energy")
plt.ylabel("Amplitude")
plt.plot (dafs[:,0], dafs[:,(sim.column_real)], label = sim.index_real)
plt.plot (dafs[:,0], dafs[:,(sim.column_im)], label = sim.index_im)
plt.legend(loc = 1)
plt.show()

plt.title("RXS %s"%sim.index)
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.plot (dafs[:,0], dafs[:,(sim.column)], label = sim.index)
plt.legend(loc = 2)
plt.show()

"""
sim.Filout("tio2_dafs_3_py_inp.txt")
sim.FDMNESfile()
sim.retrieve()
sim.do_convolution()
sim.FDMNESfile()
sim.retrieve()

dafs = sim.get_dafs((0,0,2), 1,1,45., conv=True)

plt.title("RXS %s" %sim.index)
plt.xlabel("Energy")
plt.ylabel("Intensity")
plt.plot (dafs[:,0], dafs[:,(sim.column)], label = sim.index)
plt.legend(loc = 2)
plt.show()

"""
plt.title("RXS %s"%sim.index)
plt.xlabel("Energy")
plt.ylabel("Amplitude")
plt.plot (dafs[:,0], dafs[:,(sim.column_real)], label = sim.index_real)
plt.plot (dafs[:,0], dafs[:,(sim.column_im)], label = sim.index_im)
plt.legend(loc = 2)
#plt.legend(loc = 2)
plt.show()"""