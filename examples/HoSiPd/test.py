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

os.chdir(DIR) # change to output directory

sim = fdmnes.fdmnes("A-L3_3_inp.txt")


