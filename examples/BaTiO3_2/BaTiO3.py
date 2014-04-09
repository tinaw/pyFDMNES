# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:03:25 2013

@author: weiget
"""
import os
import fdmnes
import matplotlib.pyplot as plt
import numpy as np


DIR = "output"
if not os.path.isdir(DIR):
    os.mkdir(DIR)

sim = fdmnes.fdmnes("BaTiO3.cif", resonant="Ti")

os.chdir(DIR) # move to output folder

sim.P.radius = 3.5
sim.P.Range = -10., 0.1, 10, 0.2, 30, 1, 50
#sim.P.Rpotmax = 8.50
sim.P.Green = True

zpos = np.linspace(0.4, 0.5, 6) # array of z-positions

plt.title("Convoluted BaTiO$_3$ XANES for different Ti z position")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
    
for z in zpos:
    sim.positions["Ti1"][2] = z # set z-coordinate of titanium atom
    sim.WriteInputFile("BaTiO3_py_%.2f_inp.txt"%z, overwrite=True)
    
    sim.Run(wait=True)
    
    assert sim.Status()
    
    sim.DoConvolution(overwrite=True)
    sim.Run(wait=True)
    
    data = sim.get_XANES(conv=True) # fetch convoluted xanes spectrum
    
    thislabel = " z = %1.2f"%z
    plt.plot(data[:,0], data[:,1], label = thislabel)
    

plt.legend(loc = 2)
plt.show()
