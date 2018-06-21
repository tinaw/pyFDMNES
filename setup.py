#!/usr/bin/env python
import sys
import os
import distutils.core

try:
    import ConfigParser as configparser
except ImportError:
    import configparser


conf = configparser.ConfigParser()
conf.read("setup.cfg")
if not conf.has_option("global", "fdmnes_path"):
    raise ValueError(
        "No entry for ``fdmnes_path'' found in setup.cfg")
else:
    fdmnes_path = conf.get("global", "fdmnes_path")
fdmnes_path = os.path.expanduser(fdmnes_path)
fdmnes_path = os.path.abspath(fdmnes_path)



if not os.path.isfile(fdmnes_path):
    print("File not found: %s"%fdmnes_path)
    print("Please edit file ``setup.cfg''")
else:
    confsave = configparser.RawConfigParser()
    confsave.add_section('global')
    confsave.set('global', 'fdmnes_path', fdmnes_path)
    with open(os.path.join("fdmnes", "config.ini"), 'w') as configfile:
        confsave.write(configfile)
    if len(sys.argv)<2:
        print("see install.txt for installation instructions.")
    instpackage = ["fdmnes"]
    #try:
    #    import CifFile
    #except:
    #    instpackage.append("CifFile")
    instpackage.append("CifFile")
    distutils.core.setup( name = "fdmnes", 
       version = "0.2",
       packages = instpackage,
       package_data={'fdmnes':['config.ini']},
       author = "Carsten Richter",
       author_email = "carsten.richter@desy.de",
       description = """
        Python interface to the X-Ray Spectroscopy Simulation software FDMNES
        """,
       long_description = """
        Python interface to the X-Ray Spectroscopy Simulation software FDMNES
        """
     )

