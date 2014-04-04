# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:03:25 2013

@author: weiget
"""

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
#import numpy as np

sim = pyFDMNES.pyFDMNES("BaTiO3.cif")

#sim.Range = np.array([-50, 0.2, 50., 0.5, 100., 2.0, 200])
sim.radius = 2.0
sim.Rpotmax = 8.50

#sim.Absorber = np.array([1, 3])
#sim.Quadrupole = False
#sim.Density = True
sim.Green = True

sim.convolution = True
sim.cartesian =False

sim.checkvars()
sim.Filout("BaTiO3_py_inp.txt")
sim.FDMNESfile()
sim.retrieve()

#sim.Green = True
sim.Filein("BaTiO3_py_inp.txt")
sim.Filout("BaTiO3_py_inp.txt")

sim.FDMNESfile()

sim.retrieve()
"""
data = sim.get_XANES()

plt.title("XANES with Convolution")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
x = "translation about"
plt.plot (data[:,0], data[:,1], label = "with Convolution")
plt.legend()
plt.show()
"""
i = 0
while i<1:
    num, pos = sim.pos[1]
    x,y,z = pos  
    position = [x+i, y+i, z+i]
    pos = (num, position)
  #  i += 0.5
    sim.pos.pop(1)
    sim.pos.insert(1, pos)
 #   sim.pos.pop(1)

    sim.Filout("BaTiO3_py_inp.txt")

#sim.do_convolution("tio2_dafs_py_inp_conv.txt")

    sim.FDMNESfile()

    sim.retrieve()

    data = sim.get_XANES()

    plt.title("XANES with Convolution")
    plt.xlabel("Energy")
    plt.ylabel("Absorption Cross Section")
    x = " transposition about %1.1f"%i
    plt.plot (data[:,0], data[:,1], label = x)
    plt.legend(loc = 1)
    #plt.legend(bbox_to_anchor=(1, 1), loc=2,label = transpoition)
    plt.show()
    
    i += 0.1