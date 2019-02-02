import os
import subprocess
import sys
from setuptools import setup
from setuptools import Command
from setuptools.command.build_py import build_py
from setuptools import find_packages


#########################
## config.ini patching ## 
#########################

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

def fdmnes_path():
    conf = configparser.ConfigParser()
    conf.read("setup.cfg")
    if not conf.has_option("global", "fdmnes_path"):
        raise ValueError(
            "No entry for ``fdmnes_path'' found in setup.cfg")
    else:
        fdmnes_path = conf.get("global", "fdmnes_path")
    if not os.path.isfile(fdmnes_path):
        raise IOError("File not found: {}\Please edit file ``setup.cfg''".format(fdmnes_path))
    return fdmnes_path

def create_ini(fdmnes_path):
    confsave = configparser.RawConfigParser()
    confsave.add_section('global')
    confsave.set('global', 'fdmnes_path', fdmnes_path)
    with open(os.path.join("fdmnes", "resources", "config.ini"), 'w') as configfile:
        confsave.write(configfile)

#####################
## "build" command ## 
#####################

cmdclass = {}

class BuildWithConfig(build_py):
    description = 'build with generated resources'
    
    def run(self):
        create_ini(fdmnes_path())
        build_py.run(self)

cmdclass['build_py'] = BuildWithConfig


###################
## Package setup ## 
###################

setup( name = "fdmnes", 
       version = "0.2",
       packages = find_packages(),
       package_data={'fdmnes.resources': ['spacegroup.txt', 'config.ini']},
       author = "Carsten Richter",
       author_email = "carsten.richter@desy.de",
       description = """
        Python interface to the X-Ray Spectroscopy Simulation software FDMNES
        """,
       long_description = """
        Python interface to the X-Ray Spectroscopy Simulation software FDMNES
        """,
        install_requires=[
                      'numpy',
                      'PyCifRW'
                     ],
       cmdclass=cmdclass,
       test_suite="fdmnes.tests.test_all.test_suite"
      )

