# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 18:22:28 2013

@author: weiget
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import fdmnes

sim = fdmnes.fdmnes("tio2_dafs_2_py_inp.txt")



sim.WriteInputFile("tio2_dafs_4_py_inp.txt", overwrite=True)

