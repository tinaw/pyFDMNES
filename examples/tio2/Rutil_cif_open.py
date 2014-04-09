# -*- coding: utf-8 -*-
"""
Created on Thu Aug 08 15:44:09 2013

@author: weiget
"""
import os
import fdmnes
import matplotlib.pyplot as plt


DIR = "output"
if not os.path.isdir(DIR):
    os.mkdir(DIR)


sim = fdmnes.fdmnes("TiO2.cif", resonant="Ti")

os.chdir(DIR) # move to output folder

#sim.Range = -50, 0.2, 50., 0.5, 100., 2.0, 200
sim.P.Radius = 2.0
sim.P.Rpotmax = 8.50

sim.P.Quadrupole = False
sim.P.Density = True
sim.P.Green = True
#sim.P.Self-absorption = True

#sim.P.SCF = True
#sim.P.N_self = 5

sim.P.cartesian =False


sim.WriteInputFile("tio2_dafs_py_inp.txt", overwrite=True)

sim.Run(wait=True)
sim.Status()

#sim.DoConvolution()

data = sim.get_XANES(conv = False)
plt.plot(data[:,0], data[:,1], label = "without Convolution")

sim.DoConvolution("tio2_dafs_py_inp_conv.txt", overwrite=True)
sim.Status()

data = sim.get_XANES(conv=True)
plt.plot (data[:,0], data[:,1], label = "with Convolution")
data = sim.get_XANES(conv = False)
plt.plot (data[:,0], data[:,1], label = "without Convolution")

plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
plt.legend()
plt.show()
