#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Fr 11. Apr 12:41:19 CEST 2014
# Computer: haso227r 
# System: Linux 3.2.0-60-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import os
import fdmnes # better before numpy!!

import numpy as np
import time

sim = fdmnes.fdmnes("Anatas.cif", resonant="Ti")
DIR = "output"


sim.P.Range = (-20, 0.1, 0, 0.25, 10, 1, 50)

sim.P.Green = True
sim.P.Quadrupole = True
sim.P.Convolution = True
sim.P.Radius = 3.5
sim.P.SCF = False


# 3  0  0.5    2.000 
# 3  1  0.5    2.000 
# 3  1  1.5    4.000 
# 4  0  0.5    2.000 
# 3  2  1.5    0.800 
# 3  2  2.5    1.200

#sim.P.Check_all = True

sim.P.R_self = 4.
if not os.path.isdir(DIR):
    os.mkdir(DIR)

for x in np.linspace(0,4.,21):
    sim.P.Atom_conf = [
        (1, 1, 3, 4,0,x/2., 3,2,x/2., 4,1,0),
        (1, 2, 1, 2,1,6-x/2.)
    ]
    fout = os.path.join(DIR, "Ti_atomc_x%.2f_inp.txt"%x)
    sim.WriteInputFile(fout, overwrite=True)
    sim.Run(verbose=False, wait=False)
    while True:
        NumRunning = sum([sim.Status(i)==False for i in range(len(sim.proc))])
        print("%i procs running"%NumRunning)
        if NumRunning < 2:
            break
        time.sleep(5)

#sim.P.pop("Atom_conf")
#
#fout = os.path.join(DIR, "Ti_inp.txt")
#sim.WriteInputFile(fout, overwrite=True)
#sim.Run(verbose=True, wait=True)

