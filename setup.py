#!/usr/bin/env python
import sys
import os
import distutils.core
import ConfigParser



conf = ConfigParser.ConfigParser()
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
    confsave = ConfigParser.RawConfigParser()
    confsave.add_section('global')
    confsave.set('global', 'fdmnes_path', fdmnes_path)
    with open(os.path.join("fdmnes", "config.ini"), 'wb') as configfile:
        confsave.write(configfile)
    if len(sys.argv)<2:
        print("see install.txt for installation instructions.")
    
    distutils.core.setup( name = "fdmnes", 
       version = "0.2",
       packages = ["fdmnes", "CifFile"],
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

