# -*- coding: utf-8 -*-
"""
Created on Thu Aug 08 15:44:09 2013

@author: weiget
"""
import fdmnes
import os
import matplotlib.pyplot as plt

DIR = "output"
if not os.path.isdir(DIR):
    os.mkdir(DIR)


sim = fdmnes.fdmnes("BaTiO3.cif")

os.chdir(DIR) # move to output folder

sim.P.Range = (-15, 0.5, 50)
sim.P.radius = 4.0
sim.P.Rpotmax = 8.50
sim.P.Green = True
sim.P.cartesian =False

sim.WriteInputFile("BaTiO3_py_inp.txt", overwrite=True)

sim.Run(wait=True)
sim.Status()
data = sim.get_XANES()
plt.plot (data[:,0], data[:,1], label="Green")

sim.LoadInputFile("BaTiO3_py_inp.txt")
sim.P.Green = False
sim.WriteInputFile("BaTiO3_py_2_inp.txt", overwrite=True)

sim.Run(wait=True)
sim.Status()
data = sim.get_XANES()
plt.plot (data[:,0], data[:,1], label="No Green")

sim.DoConvolution(overwrite=True)

sim.Status()
data = sim.get_XANES(conv=True)
plt.plot (data[:,0], data[:,1], label="No Green convoluted")

plt.title("XANES with Convolution")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
plt.legend(loc = 1)
plt.show()
