#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mi 4. Mar 19:40:23 CET 2015
# Computer: haso227r 
# System: Linux 3.13.0-45-generic on x86_64
#
# Copyright (c) 2015 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import os
import pylab as pl
#import evaluationtools as et
import argparse
#import ConfigParser


parser = argparse.ArgumentParser(description=
   "Vizualization of shells and charge from FDMNES *bav* files as a "
   "function of distance from the resonant atom.")
#parser.add_argument("-s", "--sort", choices=["time", "name"], default=None,
#                    help="sort input files")
parser.add_argument("bavfile", nargs=1, default=[], 
                    help="bav file to parse")
#parser.add_argument("scanno", type=int, default=None,
#                    help="Numbers of scans to plot", nargs="*")
#parser.add_argument("-c", "--col", type=str, default=None,
#                    help="Name of columns to plot for each scan")
#parser.add_argument("-o", "--outfile", type=str, 
#                    default="fdmnes_overview.csv", help="output file")
args = parser.parse_args()

fname = os.path.abspath(args.bavfile[0])

if not os.path.isfile(fname):
    raise ValueError("File not found: %s"%fname)

Z, x, y, z, dist, numato = [], [], [], [], [], []
charges = []
ionicradii = []

with open(fname) as bf:
    bfiter = iter(bf)
    for line in bfiter:
        if "atom positions in order, in the internal r2 bases" in line.lower():
            break
        else:
            continue
    header = bfiter.next()
    for line in bfiter:
        if not line.strip():
            break
        if not line.strip()[0].isdigit():
            break
        Z.append(int(line[:4]))
        x.append(float(line[4:15]))
        y.append(float(line[15:26]))
        z.append(float(line[26:37]))

        dist.append(float(line[41:53]))
        numato.append(int(line[67:71]))
    
    for line in bfiter:
        if "-- Potrmt --" in line:
            bfiter.next()
            header = bfiter.next()
            line = bfiter.next()
            if len(line)>5 and line[4:5]=="*":
                break
    
    charges.append(float(line[13:22]))
    ionicradii.append(float(line[31:40]))
    for line in bfiter:
        if not line.strip():
            break
        if not line.strip()[0].isdigit():
            break
        charges.append(float(line[13:22]))
        ionicradii.append(float(line[31:40]))

#print charges
#print charges, numato
charge = [charges[i] for i in numato]

data = pl.array([numato, Z, x, y, z, dist, charge])
#for dat in data.T:
#    print dat
pl.hist(data[5], bins = 200)
pl.xlabel("distance / $\\AA$")
#pl.ylabel("atom count")

pl.plot(data[5], data[6], "xk", label="charge")
pl.plot(data[5], data[6].cumsum(), ".-r", label="accumulated charge")
pl.plot(data[5], pl.ones(len(data[0])).cumsum(), "g", label="no. atoms")
pl.legend(loc=0)
pl.show()
