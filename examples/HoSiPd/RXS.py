# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:44:32 2013

@author: weiget
"""


import numpy as np
import matplotlib.pyplot as pl
import os
import sys
import fdmnes

DIR = "output/"
if not os.path.isdir(DIR):
    os.mkdir(DIR)

os.chdir(DIR) # change to output directory

sim = fdmnes.fdmnes("A-L3_3_inp.txt")


dafs220 = sim.get_DAFS((2,2,0), "sigma", "sigma", 0., verbose=True)
pl.plot(dafs220.Energy, dafs220[1])

dafs001 = sim.get_DAFS((0,0,1), "sigma", "sigma", 0., verbose=True)
pl.plot(dafs001.Energy, dafs001[1])


pl.show()
