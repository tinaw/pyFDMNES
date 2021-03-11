#!/usr/bin/env python
# -*- coding: utf-8 -*-F
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Author: Tina Weigel <tina.weigel@student.tu-freiberg.de>
# Created on Thu Jul 25 17:03:33 2013
# System: Linux 3.2.0-60-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter, Tina Weigel All rights reserved.
#----------------------------------------------------------------------

import os
os.environ["OPENBLAS_MAIN_FREE"] = '1'
import numpy as np
#import BytesIO
import subprocess
import collections
import itertools
import time

from . import settings
from . import elements
from .resources import resource_filename
from .parse import keyword_exists, parse_input_file, parse_bavfile, param2str, mkfloat


if os.sys.version_info[0] > 2:
    long = int

try:
    import ConfigParser as configparser
except ImportError:
    import configparser


conffile = resource_filename("config.ini")
conf = configparser.ConfigParser()
conf.read(conffile)


class ConsistencyError(Exception):
    def __init__(self, value="Inconsistent parameter values found.", 
                       errmsg="", identifier=None):
        self.value = value
        self.errmsg = errmsg
    def __str__(self):
        return self.value + os.linesep + "Message: %s" %self.errmsg


def get_polarization(pol):
    options = ["sigma", "pi", "right", "left"]
    try:
        pol = float(pol)
        pol = (5, pol)
    except:
        if pol in options:
            pol = options.index(pol) + 1
        else:
            pol = None
    return pol



def get_energy(Range):
    if len(Range)%2 != 1:
        raise ConsistencyError("Invalid input for Range."
            "tuple of odd length expected.")
    energy = []
    for i in range(0, len(Range)-2, 2):
        energy.append(np.arange(Range[i], Range[i+2], Range[i+1]))
    energy.append(Range[-1])
    return np.hstack(energy)


class Parameters(dict):
    """
        Modified dictionary class that supplies content as attributes
        and takes standard values from settings.Defaults.
    """
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError('expected dict')

    def __setitem__(self, key, value):
        keyset = keyword_exists(key)
        if keyset == None:
            raise ValueError(
                "'%s' not found in list of valid FDMNES parameters."%key
            )
        else:
            GroupName = settings.GroupMembers[keyset]
            DefVal = getattr(settings.Defaults, GroupName)[keyset]
            Type = type(DefVal)
            try:
                if Type is tuple:
                    if isinstance(value, tuple):
                        pass
                    else:
                        value = tuple((value,))
                else:
                    value = Type(value)
            except TypeError:
                raise TypeError(
                    'Expected %s for parameter %s'%(Type.__name__, keyset))
            # deep check:
            if Type is tuple:
                DefLen = len(DefVal)
                if DefLen > 0 and len(value)!=DefLen:
                    raise ValueError(
                        "Length of %s is expected to be %i"%(keyset, DefLen))
                for i in range(DefLen):
                    try:
                        value = type(DefVal)(value)
                    except:
                        raise TypeError(
                            "Element %i of %s has to be of type %s"\
                            %(i, keyset, type(DefVal)))
            dict.__setitem__(self, keyset, value)


    def __getitem__(self, key):
        if key in self:
            return dict.__getitem__(self, key)
        elif key in settings.GroupMembers:
            Group = getattr(settings.Defaults, settings.GroupMembers[key])
            return Group[key]
        else:
            raise ValueError(
                "'%s' not found in list of valid FDMNES parameters."%key)
    __setattr__ = __setitem__
    __getattr__ = __getitem__

    def __dir__(self):
        return dir(dict) + list(self.keys())



class CurrentCalculation(object):
    def __init__(self):
        self.infile = None
        self.outfile = None
        self.bavfile = None
        self.convfile = None
        self.rxsfile = None
        self.jobID = None

    @property
    def infile(self):
        return self._infile

    @infile.setter
    def infile(self, path):
        self.outfile = None
        self.convfile = None
        self.bavfile = None
        self.rxsfile = None
        self.jobID = None
        self._infile = path

    @property
    def outfile(self):
        if self._outfile is None and not self._infile is None:
            self._outfile = parse_input_file(self._infile, outpath_only=True)
        return self._outfile

    @outfile.setter
    def outfile(self, path):
        self._outfile = path







class Job(object):
    _ids = itertools.count(0)
    def __init__(self, path, command, stdout=None):
        self.id = _id = next(self._ids)
        self.path = path
        self.command = command
        self.process = subprocess.Popen(command, stdout=stdout, shell=True)
        self.finished = False
        self._successful = False

    @property
    def is_running(self):
        if self.finished:
            return False
        if self.process.poll()==None:
            return True
        else:
            self.finished = True
            return False


    @property
    def successful(self):
        if not self.finished:
            return None
        if self._successful:
            return True

        ParamIn = parse_input_file(self.path)
        P = ParamIn["Param"]
        path_out = ParamIn["path_out"]
        conv = ParamIn["conv"]
        bavfile = ParamIn["bavfile"]
        if os.path.exists(bavfile):
            self.bavfile = bavfile
            self.bavinfo = parse_bavfile(bavfile)
            if self.bavinfo["success"]:
                self._successful = True
                return True

        elif conv and "Calculation" in P and os.path.isfile(path_out):
            return True

        return False


    def wait(self):
        self.process.wait()
        self.finished = True

    def __repr__(self):
        result = "Job #%i"%self.id
        if self.finished:
            result +=" (finished)"
        result += ": %s - %s"%(self.process.__repr__(), self.path)
        return result


class Jobs(dict):
    def add(self, job):
        self[job.id] = job

    def get_running(self):
        return (self[job] for job in self if self[job].is_running)

    def get_finished(self):
        return (self[job] for job in self if not self[job].is_running)

    def get_successful(self):
        return (self[job] for job in self.get_finished() if self[job].successful)





class fdmnes(object):
    """ 
        Python interface to the x-ray spectroscopy simulation software FDMNES.
        Easy handling of FDMNES data input and output.
    """

    def __init__(self, structure, resonant="", verbose=False, 
                       fdmnes_path=None):
        """
            Calculation of spectra for x-ray absorption spectroscopy with
            FDMNES for given parameters and structures in the following steps:
                - compile input file for FDMNES 
                - process and control the FDMNES calculation
                - fit and plot the calculate datas

            Input parameters:
                structure : either
                                - metric (cell dimensions) of the crystal
                            or
                                - path to .cif-file
                            or
                                - detailed name of space group 
                                  (including setting)
        """
        self.positions = collections.OrderedDict() # symmetric unit
        self.Defaults = settings.Defaults
        self.elements = {}
        self.verbose = verbose
        self.Z = {}
        self.occupancy = {}
        self.Crystal = True
        # self.convolution = True
        self.P = Parameters()
        self.jobs = Jobs()
        self.current = CurrentCalculation() # contains paths belonging to current simulation
        self.NumCPU = 1


        ### FIND FDMNES ####################################################
        if fdmnes_path==None:
            if not conf.has_option("global", "fdmnes_path"):
                raise ValueError(
                    "No entry for ``fdmnes_path'' found in config file:%s%s"\
                    %(os.linesep, conffile))
            else:
                fdmnes_path = conf.get("global", "fdmnes_path")

        fdmnes_path = os.path.realpath(fdmnes_path)
        print("Using FDMNES at %s"%fdmnes_path)
        self.fdmnes_dir = os.path.dirname(fdmnes_path)
        self.fdmnes_exe = fdmnes_path
        fdmnes_bin = os.path.basename(fdmnes_path)

        for fname in [fdmnes_bin]:
            fpath = os.path.join(self.fdmnes_dir, fname)
            if not os.path.isfile(fpath):
                raise ValueError(
                    """
                        File %s not found in %s

                        Have you entered a valid path in %s?
                        It must point on the fdmnes executable.
                    """%(fname, self.fdmnes_dir, conffile))

        fpath = os.path.join(self.fdmnes_dir, "spacegroup.txt")
        if not os.path.isfile(fpath):
            fpath = resource_filename("spacegroup.txt")

        with open(fpath, "r") as fh:
            sgcont = filter(lambda s: s.startswith("*"), fh.readlines())
        sgcont = map(str.strip, sgcont)
        sgcont = map(lambda s: s.strip("*").replace("=",""), sgcont)
        sgcont = map(lambda s: s.split(), sgcont)
        self.spacegroups = sgcont
        self._spacegroups_flat = list(itertools.chain.from_iterable(sgcont))

        ### LOAD STRUCTURE #################################################
        try:
            if structure in self._spacegroups_flat:
                self.sg = structure
                self.a, self.b, self.c, \
                self.alpha, self.beta, self.gamma = 1, 1, 1, 90, 90, 90
                self.cif = False
            elif hasattr(structure, "__iter__") and len(structure)==6:
                structure = map(float, structure)
                self.a, self.b, self.c, \
                self.alpha, self.beta, self.gamma = structure
                self.cif = False
            elif os.path.isfile(structure):
                 end = os.path.splitext(structure)[1].lower()
                 if end == ".cif":
                     self.load_cif(structure, resonant)
                     self.cif = True
                 elif end == ".txt":
                     self.LoadInputFile(structure)
            else:
                raise ValueError("Input not understood")
        except Exception as emsg:
            print("Invalid input for structure (file not found): %s%s%s"
                %(str(structure), os.linesep, emsg))
            raise


    
    def add_atom(self, label, position, resonant=False, occupancy=1.):
        """
            Method to give parameters for the FDMNES calculation.

            Inputs:
            -------
                label : string
                    The label of the atom.
                    It has to be unique and to start with the symbold of the
                    chemical element.
                position : iterable of length 3
                    Postion of the atoms in the structure.
                resonant : bool
                    Specifies whether this atom is taking part in resonant
                    scattering / absorption.
                occupancy : float
                    Determines the probability of this site to be occupied
                    by the atom.
        """
        if type(label) is not str:
            raise TypeError("Invalid label. Need string.")
        if len(position) is not 3:
            raise TypeError("Enter 3D position object!")

        position = tuple(position)
        #label = label.replace("_", "")
        SymList = list(elements.Z)

        labeltest = label.capitalize()

        if labeltest in SymList:
            # the label consists only of the element symbol
            num = list(self.elements.values()).count(labeltest) + 1
            label += str(num) # produce a unique label
            self.elements[label] = labeltest

        elif labeltest[:2] in SymList:
            self.elements[label] = labeltest[:2]

        elif labeltest[:1] in SymList:
            self.elements[label] = labeltest[:1]

        else:
            raise ValueError("Atom label shall start with the symbol of the"
                             " chemical element. "
                             "Chemical element not found in %s"%label)

        if label in self.positions:
            raise ValueError("Atom `%s` already exists in structure. "
                             "Remove atom using .remove_atom to overwrite")

        self.Z[label] = elements.Z[self.elements[label]]
        self.positions[label] = position
        self.occupancy[label] = occupancy
        if resonant:
            if "Absorber" in self.P:
                self.P.Absorber += (len(self.positions),)
            else:
                self.P.Absorber = (len(self.positions),)

    def remove_atom(label):
        out = []
        out.append(self.positions.pop(label))
        out.append(self.elements.pop(label))
        out.append(self.Z.pop(label))
        out.append(self.occupancy.pop(label))
        return out



    def load_cif(self, path, resonant=""):
        """ Method to loads a structure from a CIF file.
            Read-out parameters are space group number, spac egroup name, 
            space group origin, metric and atom positions.

            Input:
            ------
            path : string
                The path to the CIF file.
            resonant: string
                Symbol of the resonant atom.

            Informations about CIF file:
            -----------------------------
                Hall SR, Allen FH, Brown ID (1991).
                "The Crystallographic Information File (CIF): a new standard 
                archive file for crystallography".
                Acta Crystallographica A47 (6): 655-685
        """
        if not os.path.isfile(path):
            raise IOError("File not found!")
        import CifFile
        try:
            cf = CifFile.ReadCif(path)
            cb = cf.first_block()
        except Exception as e:
            print("File doesn't seem to be a valid .cif file: %s"%path)
            print(e)
            raise

        self.Crystal = True
        # Reset Structure:
        self.positions.clear()
        self.elements.clear()
        self.P.Atom = []
        self.P.Absorber = ()
        self.P.Z_absorber = 0

        if cb.has_key("_symmetry_int_tables_number"):
            sg_num = int(cb["_symmetry_int_tables_number"])
            self.sg_num = str(sg_num)
        #_space_group_name_H-M_alt
        if cb.has_key("_symmetry_space_group_name_h-m"):
            self.sg_name = cb["_symmetry_space_group_name_h-m"]
            self.sg_name = "".join(self.sg_name.split())
            i = self.sg_name.find(":")
            if i>=0 and hasattr(self, "sg_num"):
                self.sg_num += self.sg_name[i:]

        if not hasattr(self, "sg_num") and not hasattr(self, "sg_name"):
            raise IOError("No space group found in .cif file: %s"%path)

        if hasattr(self, "sg_num") and \
            self._check_sg(self.sg_num, raiseError=False):
            self.sg = self.sg_num
        elif hasattr(self, "sg_name") and \
            self._check_sg(self.sg_name, raiseError=True):
            self.sg = self.sg_name
        else:
            raise ValueError("No valid space group given in .cif file.")

        self.a = mkfloat(cb["_cell_length_a"])
        self.b = mkfloat(cb["_cell_length_b"])
        self.c = mkfloat(cb["_cell_length_c"])
        self.alpha = mkfloat(cb["_cell_angle_alpha"])
        self.beta  = mkfloat(cb["_cell_angle_beta"])
        self.gamma = mkfloat(cb["_cell_angle_gamma"])

        if isinstance(resonant, str):
            resonant = [resonant]
        elif hasattr(resonant, "__iter__"):
            pass

        for line in cb.GetLoop("_atom_site_label"):
            label = str(line._atom_site_label)
            if hasattr(line, "_atom_site_type_symbol"):
                symbol = str(line._atom_site_type_symbol)
                symbol = filter(str.isalpha, symbol)
            else:
                symbol = filter(str.isalpha, label)
            symbol = "".join(symbol)
            px = mkfloat(line._atom_site_fract_x)
            py = mkfloat(line._atom_site_fract_y)
            pz = mkfloat(line._atom_site_fract_z)
            if hasattr(line, "_atom_site_occupancy"):
                occ = mkfloat(line._atom_site_occupancy)
            else:
                occ = 1.
            position = (px, py, pz)
            self.add_atom(label, position, (label in resonant), occ)
            if symbol in resonant and symbol!=label:
                self.P.Z_absorber = elements.Z[symbol]

    def check_parameters(self, keyw, Group=None):
        """
            Checks configured FDMNES parameters for consistency.
        """
        if not keyw in self.P:
            return True
        else:
            value = self.P[keyw]
        if keyw=="Atom" and len(value):
            numspecies = len(set(self.elements.values()))
            numconf = len(value)
            if not numconf == numspecies:
                raise ConsistencyError(errmsg="Number of given electron "
                    "configurations (%i) does not match number of different "
                    "species in structure (%i)"%(numconf,numspecies))
            if "Atom_conf" in self.P and len(self.P.Atom_conf):
                raise ConsistencyError(errmsg="Electronic configuration"
                    "given twice: Parameters Atom and Atom_conf")
        return True

    def _check_sg(self, sg, raiseError=True):
        sg = str(sg)
        if sg in self._spacegroups_flat:
            return True
        else:
            sgstr = map(lambda s: ", ".join(s), self.spacegroups)
            foundgroup = list(filter(lambda s: sg in s, sgstr))
            message = os.linesep.join(foundgroup)
            if len(foundgroup):
                message = "Did you mean one of the following groups?:%s%s"\
                          %(os.linesep, message)
            errmsg="""
            The given Space group ``%s'' was not found in the database.

            Please enter the full number of the space group that you desire.
            See the International Tables or the attribute ``spacegroups''
            for more details.

            %s
            """%(sg, message)
            if raiseError:
                raise ConsistencyError(errmsg = errmsg)
            else:
                print(errmsg)
                return False

    
    def _skip_group(self, Group, TestGroup="all"):
        ConvGroups = ["Convolution", "Experiment", "Parameter"]
        if TestGroup=="all":
            TestGroup = ["SCF", "Convolution", "Extract", "RXS"]
        if Group in TestGroup and not bool(self.P[Group]):
            return True
        if self.P.Convolution and len(self.P.Calculation) and Group not in ConvGroups:
            return True
        return False

    def write_structure(self):
        """
            Returns the structure as it is understood by FDMNES as a list
            of lines.
        """
        output = []
        if hasattr(self, "sg"):
            self._check_sg(self.sg)
            output.append("Spgroup")
            output.append(self.sg)
            output.append("")

        if any(map(lambda x: x<1, self.occupancy.values())):
            suffix = "_t"
        else:
            suffix = ""
        if self.Crystal == True:
            output.append("Crystal"+suffix)
        else:
            output.append("Molecule"+suffix)

        cell = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        output.append("  %g %g %g %g %g %g" %cell)

        self.check_parameters("Atom", settings.Defaults.Basic)
        for label in self.positions:
            if hasattr(self.P, "Atom") and len(self.P.Atom):
                atomnum = [at[0] for at in self.P.Atom].index(self.Z[label])+1
            else:
                atomnum = self.Z[label]
            pos = tuple(self.positions[label])
            line = (atomnum,) + pos
            if suffix != "":
                occ = self.occupancy[label]
                line += (occ,label)
                output.append("%i %.10g %.10g %.10g %.10g !%s"%line)
            else:
                line += (label,)
                output.append("%i %.10g %.10g %.10g !%s"%line)

        return output

    def WriteInputFile(self, path, overwrite=False, update=True):
        """ 
            Method writes an input file for the FDMNES calculation.
            The calculation parameters, for example Range and Radius etc., 
            are taken from the ``P'' object. 

            Input:
            ------
            path : string
                Path to of the input file.
            overwrite : bool
                Has to be True to overwrite existing input file.
        """
        if os.path.isfile(path) and not overwrite:
            raise IOError("File %s already exists. "
                          "Use overwrite=True to replace it."%path)

        basepath, ext = os.path.splitext(path)
        dirname, basename = os.path.split(basepath)

        conv = "Convolution" in self.P and self.P.Convolution
        convonly = bool(conv and len(self.P.Calculation))

        basename = basename.replace("_conv", "")
        suffix = "_out_conv" if convonly else "_out"
        if "_inp" in basename:
            base_out = basename.replace("_inp", suffix)
        else:
            base_out = basename + suffix

        path_out = os.path.join(dirname, base_out)

        output = ["Filout"] if not convonly else ["Conv_out"]
        output.append(path_out + ".txt"*convonly)
        output.append("")

        #output.append("Folder_dat") # not needed anymore
        #output.append(self.fdmnes_dir)
        #output.append("")

        if self.P.Extract:
            bavfile = self.P.Extract
        else:
            bavfile = path_out + "_bav.txt"

        if not convonly:
            self.P.Calculation = []
            output.extend(self.write_structure())

        for GroupName in settings.Defaults._fields:
            Group = getattr(settings.Defaults, GroupName)
            if self._skip_group(GroupName):
                continue
            for ik, keyw in enumerate(Group):
                if keyw in self.P and self.P[keyw]!=Group[keyw]:
                    if GroupName=="Parameter":
                        output.append("")
                        output.append("Parameter")
                    self.check_parameters(keyw, Group)
                    print("-> %s"%keyw)
                    value = param2str(self.P[keyw])
                    output.append(keyw)
                    if value:
                        output.append(value)

        # print empty line before each keyword
        for keyw in self.P:
            if keyw in output and not keyw.startswith("Par_"):
                ind = output.index(keyw)
                output.insert(ind, "")
        output.append("")
        output.append("! Wrote file at %s"%time.ctime())
        output.append("End")
        if convonly:
            self.current.convfile = path
        elif update:
            self.current.infile = path
            self.current.bavfile = bavfile
            self.current.outfile = path_out
        
        with open(path, "wb") as f:
            f.write(os.linesep.join(output).encode())


    
    def Run(self, path=None, wait=True, delay=5, logpath=None, verbose=False,
                  writeonly=False, command=None):
        """
            Method to write the ``fdmfile.txt'' and, subsequently, to start
            the FDMNES simulation. The simulation will be perfomed for the
            specified input file path.  If it is None, the current one will be
            run.

            Inputs:
                path : string
                    Path for fdmnes input file to be run
                wait : bool
                    If True, waits for the process to finish.
                delay : float
                    Time in seconds to wait after sending Job to background.
                    This has proven to be important if starting several jobs
                    in a row. It is ignored if `wait==True`.
                logpath : string
                    Path for log file. It is derived from input file path
                    by default.
                verbose : bool

                writeonly : book
                    If True, only writes the `fdmfile.txt`
                command : string
                    path to fdmnes executable. By default this
                    is given by self.fdmnes_exe
        """
        if path is None:
            path = self.current.infile
        if not os.path.isfile(path):
            raise ValueError("File not found:"%path)

        
        fdmfile = "fdmfile.txt"

        output = []
        output.append("1")
        output.append("")
        output.append(path)

        print("Processing: %s"%path)
        with open(fdmfile, "wb") as f:
            f.write(os.linesep.join(output).encode())

        if writeonly:
            return path

        if logpath==None and not verbose:
            basename = os.path.basename(path)
            basename, ext = os.path.splitext(basename)
            outdir = os.path.dirname(path)
            if not outdir:
                outdir = os.curdir
            flist = os.listdir(outdir)
            flist = filter(lambda s: s.startswith(basename), flist)
            flist = filter(lambda s: s.endswith(".log"), flist)
            flist =  [os.path.splitext(F)[0] for F in flist]
            lognum = [os.path.splitext(F)[1].strip(os.extsep) for F in flist]
            lognum = sorted(map(int, lognum))
            lognum = (lognum[-1]+1) if len(lognum) else 0
            logpath = os.path.join(outdir, basename + ".%i.log"%lognum)

        #errfile = open(os.path.join(self.workdir, "fdmnes_error.txt", "w"))
        if verbose:
            stdout = None
        else:
            logfile = open(logpath, "wb")
            stdout = logfile
            print("Writing output to: %s"%logpath)
        try:
            if command==None:
                command = self.fdmnes_exe
            job = Job(path, command, stdout)
            self.jobs.add(job)
            if path == self.current.infile:
                self.current.jobID = job.id
            if wait:
                job.wait()
                print("FDMNES simulation finished.")
            else:
                print("FDMNES process #%i started in background."%job.id)
                print("See the ``Status'' method for details.")
                time.sleep(delay)
        except Exception as e:
            print("An error occured when running FDMNES command:")
            print(command)
            print("Message: %s"%e)
        finally:
            if not verbose:
                logfile.close()
        return job

    def wait(self):
        for job in self.jobs.get_running():
            self.Status(job.id)
            job.wait()
        return True
            

    def Status(self, jobID=None, full_output=False, verbose=True):
        """ 
            Method to check the status of running simulations.

            Returns:
                False, if a process is running
                True, if no process is running  
                2, if the specified job is finished (last by default)
        """
        message = []
        njobs = len(self.jobs)
        if jobID==None and njobs:
            jobID = max(self.jobs)
        if not njobs:
            result = True
            message.append("No process was started.")
        elif self.jobs[jobID].is_running:
            result = False
            message.append("Simulation is running for job #%i: %s"\
                                        %(jobID, self.jobs[jobID]))
        else:
            result = True
            job = self.jobs[jobID]
            message.append("Current Job: %s"%job.__repr__())
            if job.successful:
                if hasattr(job, "bavfile"):
                    message.append("Job bav-file: %s"%job.bavfile)
                result = 2
            else:
                message.append("Job aborted.")
        message.append("")
        if full_output:
            return result, message
        else:
            if verbose:
                print(os.linesep.join(message))
            return result

        
    def LoadInputFile(self, fpath):
        """
            Method to read an existing input file.

            Input:
            ------
                fpath: string
                    Path to the input file.
        """
        self.P.clear()
        self.positions.clear()
        self.elements.clear()
        self.cif = False

        ParamIn = parse_input_file(fpath)
        self.ParamIn = ParamIn

        
        path_out = ParamIn["path_out"]

        if ParamIn["sg"]!=None:
            self.sg = ParamIn["sg"]
        structure = ParamIn["structure"]
        ParamIn = ParamIn["Param"]

        if "Extract" in ParamIn and ParamIn["Extract"]:
            bavfile = ParamIn["Extract"][0]
        else:
            bavfile = path_out + "_bav.txt"

        if structure and list(structure)[0].startswith("Crystal"):
            self.Crystal = True

        if all([k in ParamIn for k in ["Atom", "Atom_conf"]]):
            raise ConsistencyError(
                "Only one of the keywords ``Atom'' and "\
                "``Atom_conf'' may be given in the input file!")

        for keyw in ParamIn:
            if not len(ParamIn[keyw]):
                self.P[keyw] = True
                continue
            Values = ParamIn[keyw]
            GroupName = settings.GroupMembers[keyw]
            Group = getattr(self.Defaults, GroupName)
            Type = type(Group[keyw])
            for j in range(len(Values)):
                Value = Values[j]
                Value = Value.split() if Type in [list, tuple] else [Value]
                for i in range(len(Value)):
                    try:
                        Value[i] = int(Value[i])
                        continue
                    except:
                        pass
                    try:
                        Value[i] = float(Value[i])
                        continue
                    except:
                        pass
                if len(Value) > 1:
                    Value = tuple(Value)
                else:
                    Value = Value[0]
                Values[j] = Value
            if Type is not list:
                assert len(Values) == 1, "Many lines of input given for " \
                    "Parameter %s. Only one line of input accepted! " \
                    "Or unsupported keyword given: %s"%(keyw, Values[1])
                Values = Values[0]
            self.P[keyw] = Values

        if structure and not list(structure)[0].endswith("_p"):
            structure = list(structure.values())[0]
            cell = map(float, structure.pop(0).split())
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma = cell
            if "Atom" in self.P and len(self.P.Atom):
                Z = lambda i: self.P.Atom[i-1][0]
            else:
                Z = lambda i: i

            for Atom in structure:
                Atom = Atom.split()
                symbol = elements.symbols[Z(int(Atom[0]))]
                position = tuple(map(float, Atom[1:4]))
                occ = float(Atom[4]) if (len(Atom)>4) else 1
                self.add_atom(symbol, position, occupancy = occ)

        self.current.infile = os.path.relpath(fpath)
        self.current.bavfile = bavfile


    def DoConvolution(self, convpath=None, overwrite=False, writeonly=False):
        """
            Method to write a special convolution file. It is also possible to
            define the convolution parameters in the WriteInputFile method.

            Input:
            ------
                path: string
                    Name of the convolution file.
        """

        path_out = self.current.outfile
        job = self.jobs[self.current.jobID]


        self.P.Convolution = True
        if not "Calculation" in self.P:
            self.P.Calculation = []
        if not self.P.Calculation and job.successful:
            num_absorb = job.bavinfo["num_absorber"]
            if num_absorb > 1:
                calclist = ["_%.i"%(i+1) for i in range(num_absorb)]
            else:
                calclist = [""]
            calclist = map(lambda s: path_out + s + ".txt", calclist)
            self.P.Calculation.extend(calclist)

        foundCalc = list(map(os.path.isfile, self.P.Calculation))
        if not foundCalc or not foundCalc[0] or not any(foundCalc):      
            raise ValueError(
                "``Calculation'' files for Convolution not found or "\
                "Calculation not finished. "\
                "Check P.Calculation list or Run() Simulation in advance .")

        if convpath == None:
            curdir = os.path.dirname(self.current.infile)
            basename = os.path.basename(self.current.infile)
            if not "_conv" in basename:
                basename = "_conv".join(os.path.splitext(basename))
            if not "_inp" in basename:
                basename = basename.replace("_conv", "_inp_conv")
            convpath = os.path.join(curdir, basename)

        self.current.convfile = convpath
        self.WriteInputFile(convpath, overwrite)
        if not writeonly:
            self.Run(convpath, wait=True)
        self.P.Convolution = False
        return convpath

    
    def get_XANES(self, conv=False, fpath=None):
        """
            Method to read-out the calculated XANES datas of the output file 
            created by FDMNES.

            Input:
            -----
                conv: bool (default: False)
                    If the calculation after convolution is desired.

                fpath: string (optional)
                    A specific FDMNES output can be given. By default, if
                    choosing unconvoluted data, all absorbing atoms will be
                    returned.

        """
        conv = bool(conv)

        path_out = self.current.outfile
        if path_out is None:
            raise ValueError("No calculation loaded.")

        if os.path.isfile(path_out + "_tddft.txt") and self.P.TDDFT:
            path_out += "_tddft"

        if fpath==None:
            fpath = path_out + ".txt"
        if conv and not fpath.endswith("_conv.txt"):
            if os.path.isfile(path_out + "_conv.txt"):
                fpath = path_out + "_conv.txt"
            else:
                print("Doing Convolution...")
                fpath = self.DoConvolution(overwrite=True)
                #fpath = path_out + ".txt"
        # Can be the case when writeonly == True:
        elif not conv and fpath.endswith("_conv.txt"):
            if len(self.P.Calculation):
                fpath = self.P.Calculation
            else:
                fpath = fpath.replace("_conv", "")

        if not isinstance(fpath, list):
            fpath = [fpath]

        output = []
        for fp in fpath:
            if os.path.isfile(fp):
                print("Using data from %s"%fp)
                with open(fp, "r") as fh:
                    content = fh.readlines()
                content = map(str.strip, content)
                findEne = list(map(lambda s: s.startswith("Energy"), content))
                ind = (findEne.index(True)+1) if any(findEne) else 0
                data = np.loadtxt(fp, skiprows = ind)
            else:
                raise IOError(
                    "FDMNES output file not found:%s%s"%(os.linesep, fp))
            if not len(output):
                output.append(data[:,0])
            output.append(data[:,1])

        return np.vstack((output)).T

        
    def get_DAFS(self, miller, pol_in, pol_out, azimuth, conv=True,
                 verbose=False):
        """

            Method to read-out the calculated DAFS datas. It's possible to get
            the DAFS data's for given paramters with and without convolusion.
            If the given parameters doesn't exist, it is possible to add these
            parameter to the input file.

            Input:
            ------
                miller: 3-tuple (int)
                    reflection indexes
                pol_in and pol_out: 
                    orientation of the polarization of incident and scattered
                    photons. Possible values are:
                        "sigma" - perpendicular to plane of scattering
                        "pi"    - parallel to plane of scattering
                        "right" - circular right polarized
                        "left"  - circular left polarized
                        float   - any float value giving the angle of 
                                  polarization with respect to the scattering
                                  plane normal in degree.
                                  (A value of 90 corresponds to ``pi'', 
                                               0 corresponds to ``sigma'')
                azimuth: float
                    azimuthal angle
                conv: bool
                    pick convoluted data?
        """

        assert all(map(lambda x: isinstance(x, (int, long)), miller)) \
            and len(miller)==3,\
            "Wrong input for ``miller''- (reflection-) indices. "\
            "A 3-tuple of type int is required."
        miller = tuple(miller)

        conv = bool(conv)

        self.P.Extract = self.current.bavfile
        bavinfo = parse_bavfile(self.current.bavfile)
        wasconv = self.P.Convolution
        self.P.Convolution = conv


        pol_in = get_polarization(pol_in)
        if pol_in == None:
            raise ValueError("Invalid input for pol_in")
        pol_out = get_polarization(pol_out)
        if pol_out == None:
            raise ValueError("Invalid input for pol_out")

        angles = [0,90] # what to do with left and right?
        if isinstance(pol_in, tuple) and not isinstance(pol_out, tuple):
            pol_out = (pol_out, angles[pol_out-1])
        if not isinstance(pol_in, tuple) and isinstance(pol_out, tuple):
            pol_in = (pol_in, angles[pol_in-1])


        azimuth = float(azimuth)

        self.P.RXS = [miller + (pol_in, pol_out, azimuth)]
        path = "_rxs".join(os.path.splitext(self.current.infile))
        self.WriteInputFile(path, overwrite=True, update=False)

        self.Run(path, wait=True, verbose=verbose)
        self.Extract = ""
        self.P.Convolution = wasconv



        if conv:
            fpath = self.current.outfile + "_rxs_conv.txt"
        else:
            num_absorb = bavinfo["num_absorber"]
            if num_absorb > 1:
                calclist = ["_%.i"%(i+1) for i in range(num_absorb)]
            else:
                calclist = [""]
            fpath = ["%s_rxs%s.txt"%(self.current.outfile, s)  for s in calclist]

        self.current.rxsfile = fpath

        with open(fpath, "r") as fh:
            content = fh.readlines()
        content = list(map(str.strip, content))
        findEne = list(map(lambda s: s.startswith("Energy"), content))
        ind = (findEne.index(True)+1) if any(findEne) else 0
        header = content[ind-1] if ind>0 else ""
        header = header.replace("(","").replace(")", "")
        header = header.split()
        #header = map(lambda s: "F"+s if s[0].isdigit() else s, header)
        header.pop(1) # XANES
        #DAFSdata = collections.namedtuple("DAFSdata", header)
        data = np.loadtxt(fpath, skiprows = ind)
        data = list(data.T)
        data.pop(1) # XANES
        DAFS = collections.OrderedDict(zip(header, data))
        #DAFS = DAFSdata(*data)

        return DAFS
