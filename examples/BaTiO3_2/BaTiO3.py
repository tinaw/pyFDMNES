# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:03:25 2013

@author: weiget
"""
import os
import sys
#sys.path.append("../../")

import numpy as np
import fdmnes
import matplotlib.pyplot as plt

sim = fdmnes.fdmnes("BaTiO3.cif", resonant="Ti")

sim.P.radius = 2.0
sim.P.Range = -10., 0.1, 10, 0.2, 30, 1, 50
sim.P.Rpotmax = 8.50
sim.P.Green = True

zpos = np.linspace(0.3, 0.5, 11) # array of z-positions

plt.title("Convoluted BaTiO$_3$ XANES for different Ti z position")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
    
for z in zpos:
    sim.positions["Ti1"][2] = z # set z-coordinate of titanium atom
    sim.FileOut("BaTiO3_py_2_inp.txt", overwrite=True)
    
    sim.FDMNESfile()
    
    sim.retrieve()
    sim.do_convolution()
    sim.FDMNESfile()
    data = sim.get_XANES(conv=True) # fetch convoluted xanes spectrum
    
    thislabel = " z = %1.2f"%z
    plt.plot(data[:,0], data[:,1], label = thislabel)
    

plt.legend(loc = 2)
plt.show()
