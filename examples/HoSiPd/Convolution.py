# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:44:32 2013

@author: weiget
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import fdmnes

DIR = "output/"
if not os.path.isdir(DIR):
    os.mkdir(DIR)

sim = fdmnes.fdmnes("A-L3_inp.txt")


os.chdir(DIR) # change to output directory

sim.P.Range = (-50., 2., -30., 0.4, 30., 1., 100.)
#sim.P.Range = np.array([7900., 100.,10000.])
sim.P.Radius = 2.0
#sim.P.Self_absorption = False
sim.P.convolution = False
sim.P.Memory_save = True
sim.P.Rpotmax = 15.4
sim.P.Edge = "L3"

sim.WriteInputFile("A-L3_2_inp.txt", overwrite=True)

sim.Run(wait=True)

sim.Status(verbose=True)

sim.LoadInputFile("A-L3_2_inp.txt")

#sim.convolution = True
sim.P.RXS = [(2,2,0, 1,1,0),
             (2,2,0, 1,2,0),
             (0,0,1, 1,1,0),
             (0,0,1, 1,2,0),
             (6,0,1, 1,1,0),
             (6,0,1, 1,2,0)]

sim.P.Extract = sim.bavfile
sim.WriteInputFile("A-L3_3_inp.txt", overwrite=True)
sim.Run(wait=True)

sim.DoConvolution(overwrite=True)
sim.Status()
