# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:44:32 2013

@author: weiget
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.abspath(os.path.pardir))

import pyFDMNES

sim = pyFDMNES.pyFDMNES("A-L3_inp.txt")

sim.Range = np.array([-50., 2., -30., 0.4, 30., 1., 100.])
#sim.Range = np.array([7900., 100.,10000.])
sim.Radius = 2.0
#sim.Self_absorption = False
#sim.convolution = True
sim.Memory_save = True
sim.Rpotmax = 15.4
#sim.Absorber = ' '

sim.Filout("A-L3_2_inp.txt")

sim.FDMNESfile()

sim.retrieve()

sim.Filein("A-L3_2_inp.txt")

sim.extract = True
#sim.convolution = True
sim.Reflex = np.array([
                      [2,2,0, 1,1],
                      [2,2,0, 1,2],
                      [0,0,1, 1,1],
                      [0,0,1, 1,2],
                      [6,0,1, 1,1],
                      [6,0,1, 1,2]
                      ])

sim.Filout("A-L3_3_inp.txt")

#sim.Cal = 'a.txt'
#sim.Cal["a.txt"] = [1.0]#, 0.0]
#sim.Cal["b.txt"] = [3.0]#, 0.0]
#sim.Conv_out = 'A-L3-v23.txt'
#sim.Scan = ['c.txt', 'd.txt']
#sim.Scan_conv = "A-L3-v23_scan_conv.txt"
sim.Exp["exp.dat "] = [3, 2]

sim.check_conv = True
sim.Gen_shift = '8050  8090 40'
sim.Efermi = -5
sim.Ecent = 30
sim.Elarg = 30
sim.Gamma_max = 15 
#sim.Memory_save = True
sim.Fprime = True
sim.Fprime_atom = True

sim.do_convolution()#"A-L3_2_inp_conv.txt")

sim.FDMNESfile()

sim.retrieve()

#sim.do_convolution("A-L3_inp_conv_2_.txt")
"""
data = sim.get_XANES()
plt.title("XANES")
plt.xlabel("Energy")
plt.ylabel("Absorption Cross Section")
plt.legend()
plt.plot (data[:,0], data[:,1], label = "simulation")
plt.show()

exp_data = np.loadtxt("exp_data.dat",skiprows =1)
plt.legend()
plt.plot(exp_data[:,0], exp_data[:,8], label = "experimental data")
plt.show()"""