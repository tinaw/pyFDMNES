#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Mi 11. Jun 15:26:16 CEST 2014
# Computer: haso227r 
# System: Linux 3.2.0-64-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import os
import argparse
#import pylab as pl
#import evaluationtools as et
import fdmnes
import collections

sep = "; "

parser = argparse.ArgumentParser(description=
   "Overview creator for fdmnes input files")
parser.add_argument("-s", "--sort", choices=["time", "name"], default=None,
                    help="sort input files")
parser.add_argument("file", nargs="*", default=[], 
                    help="list of files")
parser.add_argument("-o", "--outfile", type=str, 
                    default="fdmnes_overview.csv", help="output file")
args = parser.parse_args()


flist = filter(os.path.isfile, args.file)

if args.sort=="time":
    flist = sorted(flist, key=os.path.getmtime)
elif args.sort=="name":
    flist = sorted(flist)
    

default = dict()
diff = []

for fname in flist:
    sim = fdmnes.fdmnes(fname)
    diff.append(collections.defaultdict(str))
    for key in sim.P:
        if key not in default or default[key]!=sim.P[key]:
            default[key] = sim.P[key]
            diff[-1][key] = str(sim.P[key])


with open(args.outfile, "w") as fh:
    fh.write(sep)
    fh.write(sep.join([key for key in default]))
    fh.write(os.linesep)
    for i, fname in enumerate(flist):
        fh.write(fname + sep)
        fh.write(sep.join([diff[i][key] for key in default]))
        fh.write(os.linesep)

