pyFDMNES
========

Python interface to fdmnes simulations. Simplifies and automates data input/output.
- load cif files
- Read/Modify/Write FDMNES input files
- easily perform repeated `Colvolution` on existing DOS calculation
- flexibly extract `DAFS`/`XANES` curves from existing calculation

A new FDMNES simulation is represented by the `fdmnes` object which  can be created by loading an existing (manually created) fdmnes input file, a .cif file of by defining space-group or metric and subsequently fill the unit cell with atoms.

After loading an existing input file, one can for example change some parameters in the `fdmnes.P` Parameter instance and save/create new inputfiles with WriteInputFile. One can add the saved file to the job queue which can be executed using sim.Run().

Furthermore, the `fdmnes.P` Parameters Object similar to a python dict but does some consistency check of all given parameters, to a certain extent. Basically it checks if the given parameter has the proper shape (2d, 1d, 0d, or bool, etc..). How the definitions are stored, can be seen in the settings file:
https://github.com/tinaw/pyFDMNES/blob/master/fdmnes/settings.py The idea is also to identify parameters which are crucial for the computation of the matrix elements and others which are used for post-processing (e.g. Convolution). This would allow to let pyFDMNES decide when a new computation of the matrix elements or electron density is necessary.


Certainly, FDMNES evolved and not all parameters are correctly represented in the settings file. Therefore it is not complete yet and updates are needed.


# Requirements:
- python 2.7 or 3
- numpy
- PyCifRW
- fdmnes (version > March 28th 2014) from e.g. https://fdmnes.neel.cnrs.fr/


# Installation
Change dir to the downloaded repository (where the setup.py is found) and modify the `setup.cfg`.

## using pip (recommended)
    pip install . [--user]

## using setup script
    python setup.py install [--user]

# Usage
See `examples` folder.

You may set the `FDMNES` environment variable to specify the FDMNES executable command,
for instance when using MPI:
```
FDMNES='mpirun fdmnes'
FDMNES='/path/to/fdmnes_linux_serial'
FDMNES='/path/to/fdmnes_linux_parallel'
```
The default is to use the executable defined in the `setup.cfg` file of the pyFDMNES package.

You may as well set the `FDMNES_SPACEGROUP` environment variable to specify the `spacegroup.txt` file which defines the space-group operators. The default is to use the file provided by the pyFDMNES package.

