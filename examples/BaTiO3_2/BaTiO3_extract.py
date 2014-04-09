# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:03:25 2013

@author: weiget
"""
import fdmnes
import pylab as pl
import os

DIR = "output"
if not os.path.isdir(DIR):
    os.mkdir(DIR)

os.chdir(DIR) # move to output folder

sim = fdmnes.fdmnes("BaTiO3_py_0.50_inp.txt")

sim.P.Extract = "BaTiO3_py_0.50_out_bav.txt"

sim.P.RXS = [(0,0,1,1,1,45.)]
sim.P.RXS.append((0,0,1,1,1,0.))

sim.WriteInputFile("text.txt", overwrite=True)

sim.Run(wait=True)

convpath = sim.DoConvolution(overwrite=True)

xanes   = sim.get_XANES(conv=False)
xanes_c = sim.get_XANES(conv=True)

pl.plot(xanes[:,0], xanes[:,1], label = "before convolution")
pl.plot(*xanes_c.T, label = "after convolution")


pl.show()
