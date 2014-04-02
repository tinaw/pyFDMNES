#!/usr/bin/env python
# -*- coding: utf-8 -*-F
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Author: Tina Weigel <tina.weigel@student.tu-freiberg.de>
# Created on Thu Jul 25 17:03:33 2013
# System: Linux 3.2.0-60-generic on x86_64
#
# Copyright (c) 2014 Carsten Richter, Tina Weigel  All rights reserved.
#----------------------------------------------------------------------

import numpy as np
import StringIO
import os
import elements
import subprocess
import string
import collections
import settings
import ConfigParser


conffile = os.path.join(os.path.dirname(__file__), "config.ini")
conf = ConfigParser.ConfigParser()
conf.read(conffile)



def array2str(array, precision=4):
    array = np.array(array)
    if array.ndim<2:
            array = array[:,np.newaxis].T
    if array.dtype == int:
        zahl = "%i "
    elif array.dtype == float:
        zahl = "%." + str(int(precision)) + "f " 
    save_value = StringIO.StringIO()
    np.savetxt(save_value, array, fmt = zahl)   
    values = save_value.getvalue()
    return values


def param2str(param):
    if isinstance(param, bool):
        return ""
    if isinstance(param, str):
        return param
    if isinstance(param, int):
        return "%i"%param
    if isinstance(param, float):
        return "%g"%param
    if isinstance(param, list):
        assert all(map(lambda x: not isinstance(x, list), param))
        return os.linesep.join(map(param2str, param))
    if isinstance(param, tuple):
        assert all(map(lambda x: isinstance(x, (bool,str,int,float,tuple)),
                       param))
        return " ".join(map(param2str, param))

def mkfloat(string):
    i = string.find("(")
    if i>=0:
        string = string[:i]
    return float(string)


def keyword_exists(key):
    key = key.lower()
    keyset = None
    for Group in settings.Defaults:
        keys = Group.keys()
        keysl = map(str.lower, keys)
        if key in keysl:
            i = keysl.index(key)
            keyset = keys[i]
            break
    return keyset


class ConsistencyError(Exception):
    def __init__(self, value="Inconsistent parameter values found.", 
                       errmsg="", identifier=None):
        self.value = value
        self.errmsg = errmsg
    def __str__(self):
        return self.value + os.linesep + "Message: %s" %self.errmsg




class Parameters(dict):
    """
        Modified dictionary class that supplies content as attributes.
    """
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError, 'expected dict'

    def __setitem__(self, key, value):
        keyset = keyword_exists(key)
        if keyset == None:
            raise ValueError(
                "'%s' not found in list of valid FDMNES parameters."%key
            )
        else:
            dict.__setitem__(self, keyset, value)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)
    __setattr__ = __setitem__
    __getattr__ = __getitem__
    
    def __dir__(self):
        return dir(dict) + self.keys()



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
                - process and controll the FDMNES calculation
                - fit and plot the calculate datas

            Input parameters:
                structure : either
                                - metric (cell dimensions) of the crystal
                            or
                                - path to .cif-file
        """
        self.positions = collections.OrderedDict() # symmetric unit
        self.Defaults = settings.Defaults
        self.Groups = self.Defaults._fields
        self.GroupMembers = dict()
        for Group in self.Groups:
            for Member in getattr(self.Defaults, Group):
                self.GroupMembers[Member] = Group
        self.elements = {}
        self.verbose = verbose
        self.resonant = []
        self.Z = {}
        self.occupancy = {}
        self.Crystal = True
        # self.convolution = True
        self.Exp = {}
        self.P = Parameters()
        self.Jobs = []
        
        ### FIND FDMNES ####################################################
        if fdmnes_path==None:
            if not conf.has_option("global", "fdmnes_path"):
                raise ValueError(
                    "No entry for ``fdmnes_path'' found in config file:%s%s"\
                    %(os.linesep, conffile))
            else:
                fdmnes_path = conf.get("global", "fdmnes_path")
                fdmnes_path = os.path.expanduser(fdmnes_path)
                fdmnes_path = os.path.abspath(fdmnes_path)
        
        self.fdmnes_dir = os.path.dirname(fdmnes_path)
        self.fdmnes_exe = fdmnes_path
        fdmnes_bin = os.path.basename(fdmnes_path)
        
        for fname in [fdmnes_bin, "xsect.dat", "spacegroup.txt"]:
            fpath = os.path.join(self.fdmnes_dir, fname)
            if not os.path.isfile(fpath):
                raise ValueError(
                    """
                        File %s not found in %s
                        
                        Have you entered a valid path in %s?
                        It must point on the fdmnes executable.
                        Further, the files xsect.dat and spacegroup.txt are
                        required in the same folder.
                    """%(fname, self.fdmnes_dir, conffile))
        with open(fpath, "r") as fh:
            sgcont = filter(lambda s: s.startswith("*"), fh.readlines())
            sgcont = map(str.strip, sgcont)
            sgcont = map(lambda s: s.strip("*"), sgcont)
            sgcont = map(lambda s: s.split()[:3], sgcont)
            self.spacegroups = np.vstack(sgcont).astype(str)
        
        ### LOAD STRUCTURE #################################################
        if structure in self.spacegroups:
            self.sg = structure
            self.a, self.b, self.c, \
            self.alpha, self.beta, self.gamma = 1, 1, 1, 90, 90, 90
            self.cif = False
        elif isinstance(structure, (tuple, list, np.ndarray)):
            structure = np.ravel(structure).astype(float)
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
            raise ValueError(
                "Invalid input for structure (file not found): %s"\
                %str(structure))
        
        
    
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
        
        position = np.array(position)
        #label = label.replace("_", "")
        labeltest = label[0].upper()
        if len(label) > 1:
            labeltest += label[1].lower()
        if labeltest in elements.Z.keys():
            self.elements[label] = labeltest
        elif labeltest[:1] in elements.Z.keys():
            self.elements[label] = labeltest[:1]
        else:
            raise ValueError("Atom label shall start with the symbol of the"
                             " chemical element. "
                             "Chemical element not found in %s"%label)
        self.Z[label] = elements.Z[self.elements[label]]
        self.positions[label] = position
        self.occupancy[label] = occupancy
        if resonant:
            if self.P.has_key("Absorber"):
                self.P.Absorber += (len(self.positions),)
            else:
                self.P.Absorber = (len(self.positions),)
        
        
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
            print e
            return
        
        self.Crystal = True
        # Reset Structure:
        self.positions.clear()
        self.P.Atom = []
        self.P.Absorber = ()
        self.P.Z_absorber = 0
        
        if cb.has_key("_symmetry_int_tables_number"):
            sg_num = int(cb["_symmetry_int_tables_number"])
            self.sg_num = str(sg_num)
        
        if cb.has_key("_symmetry_space_group_name_h-m"):
            self.sg_name = cb["_symmetry_space_group_name_h-m"]
            self.sg_name = "".join(self.sg_name.split())
            i = self.sg_name.find(":")
            if i>=0:
                self.sg_num += self.sg_name[i:]
                self.sg_name += self.sg_name[i:]
        
        if not hasattr(self, "sg_num") and not hasattr(self, "sg_name"):
            raise IOError("No space group found in .cif file:%s%s"\
                          %os.linesep, path)
        
        if hasattr(self, "sg_num") and \
            self.check_sg(self, self.sg_num, raiseError=False):
            self.sg = self.sg_num
        elif hasattr(self, "sg_name") and \
             self.check_sg(self, self.sg_name, raiseError=True):
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
            self.resonant = [resonant]
        elif hasattr(resonant, "__iter__"):
            self.resonant.extend(resonant)
        
        for line in cb.GetLoop("_atom_site_label"):
            symbol = line._atom_site_type_symbol
            label = line._atom_site_label
            px = mkfloat(line._atom_site_fract_x)
            py = mkfloat(line._atom_site_fract_y)
            pz = mkfloat(line._atom_site_fract_z)
            occ = mkfloat(line._atom_site_occupancy)
            position = (px, py, pz)
            self.add_atom(label, position, (label in self.resonant), occ)
            if symbol == resonant:
                self.P.Z_absorber = elements.Z[symbol]
        
    
    def check_parameters(self, keyw, Group=None):
        """
            Checks configured FDMNES parameters for consistency.
        """
        if not self.P.has_key(keyw):
            return True
        else:
            value = self.P[keyw]
        if Group!=None:
            deftype = type(Group[keyw])
            assert isinstance(self.P[keyw], deftype),\
                "Wrong type: Parameter %s has to be of type %s"\
                %(keyw, deftype)
        if keyw=="Atom" and len(value):
            numspecies = len(np.unique(self.elements.values()))
            numconf = len(value)
            if not numconf == numspecies:
                raise ConsistencyError(errmsg="Number of given electron "
                    "configurations (%i) does not match number of different "
                    "species in structure (%i)"%(numconf,numspecies))
            if self.P.has_key("Atom_conf") and len(self.P.Atom_conf):
                raise ConsistencyError(errmsg="Electronic configuration"
                    "given twice: Parameters Atom and Atom_conf")
        return True
    
    def _check_sg(self, sg, raiseError=True):
        sg = str(sg)
        if sg in self.spacegroups:
            return True
        else:
            message = ""
            for group in self.spacegroups:
                groupstr = ", ".join(group)
                if sg in groupstr:
                    message += groupstr + os.linesep
            if message:
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
                print errmsg
                return False
    
    def _skip_group(self, Group, TestGroup="all"):
        if TestGroup=="all":
            TestGroup = ["SCF", "Convolution", "Extract", "RXS"]
        for Name in TestGroup:
            if Group.has_key(Name):
                if not hasattr(self.P, Name) or not self.P[Name]:
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
        for label in self.positions.iterkeys():
            if hasattr(self.P, "Atom") and len(self.P.Atom):
                atomnum = self.P.Atom.index(self.Z[label]) + 1
            else:
                atomnum = self.Z[label]
            pos = tuple(self.positions[label])
            line = (atomnum,) + pos
            if suffix != "":
                occ = self.occupancy[label]
                line += (occ,label)
                output.append("%i %.10g %.10g %.10g %.g !%s"%line)
            else:
                line += (label,)
                output.append("%i %.10g %.10g %.10g !%s"%line)
        
        return output
    
    def WriteInputFile(self, path, overwrite=False):
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
        
        
        
        self.path = path
        basepath = os.path.splitext(path)[0]
        
        if "_inp" in path:
            self.path_out = basepath.replace("_inp","_out")
        else:
            self.path_out = basepath + "_out"
        output = ["Filout"]
        output.append(self.path_out)
        output.append("")
        
        output.append("Folder_dat")
        output.append(self.fdmnes_dir)
        output.append("")
        
        self.bavfile = self.path_out + "_bav.txt"
        
        output.extend(self.write_structure())
        
        for Group in settings.Defaults:
            if self._skip_group(Group):
                continue
            for keyw in Group.iterkeys():
                if self.P.has_key(keyw) and self.P[keyw]!=Group[keyw]:
                    self.check_parameters(keyw, Group)
                    print("-> %s"%keyw)
                    value = param2str(self.P[keyw])
                    output.append(keyw)
                    if value:
                        output.append(value)
        
        # print empty line before each keyword
        for keyw in self.P.viewkeys():
            if keyw in output:
                ind = output.index(keyw)
                output.insert(ind, "")
        with open(path, "w") as f:
            f.writelines(os.linesep.join(output))
    
    def AddToJobs(self):
        """
            Add the current job (the current input file) to the list of jobs.
        """
        assert hasattr(self, "path"), \
            "Attribute `path` has not been defined. "\
            "Run method WriteInputFile() first."
        self.Jobs.append(self.path)
        
    
    def Run(self, jobs="last", wait=False):
        """
            Method to write the ``fdmfile.txt'' and, subsequently, to start
            the FDMNES simulation.
            By default, the last job (input file) will be run. By giving 
            the additional argument ``jobs'' one can define a list of jobs
            (= paths to input files) that will be processed.
        """
        if jobs=="last":
            assert hasattr(self, "path"), \
                "Attribute `path` has not been defined. "\
                "Run method self.WriteInputFile() first."
            jobs = [self.path]
        if len(jobs)==0:
            print("Warning: No jobs to process!")
            return None
        
        self.Jobs = jobs
        
        fdmfile = "fdmfile.txt"
        
        output = []
        Run_File = 1
        output.append("%i"%len(jobs))
        output.append("")
        output.extend(jobs)
        
        with open(fdmfile, "w") as f:
            f.writelines(os.linesep.join(output))
        
        
        logpath = "fdmnes.log"
        logfile = open(logpath, "w")
        #errfile = open(os.path.join(self.workdir, "fdmnes_error.txt", "w"))
        if self.verbose:
            stdout = None
        else:
            stdout = logfile
            print("Writing output to:%s%s"%(os.linesep, logpath))
        try:
            self.proc = subprocess.Popen(self.fdmnes_exe, stdout=stdout)
                                                      #stderr = errfile)
            if wait:
                self.proc.wait()
                print("FDMNES simulation finished.")
            else:
                print("FDMNES process started in background.")
                print("See the ``Status'' method for details.")
        except Exception as e:
            print("An error occured when running FDMNES executable at")
            print(self.fdmnes_exe)
            print("Message: %s"%e)
        finally:
            logfile.close()
        
    def Status(self, full_output=False):
        """ 
            Method to check the status of running simulations.
            
            Returns:
                True, if a process is running
                False, if no process is running
        """
        message = []
        if not hasattr(self, "proc"):
            result = False
            message.append("No process was started.")
        elif self.proc.poll() == None:
            result = True
            message.append("FDMNES simulation is running.")
            k = len(self.Jobs)
            for i in range(k):
                path = self.Jobs[i]
                basepath = os.path.splitext(path)[0]
                if "_inp" in path:
                    path_out = basepath.replace("_inp","_out")
                else:
                    path_out = basepath + "_out"
                bavfile = path_out + "_bav.txt"
                if os.path.exists(self.bavfile):
                    with open(self.bavfile, "r") as bf:
                        bf.seek(-60, 2)
                        tail = map(str.strip, bf.readlines())
                    if tail[-1]=='Have a beautiful day !':
                        message.append("Job %i/%i finished: %s"%(i,k,path))
                    else:
                        message.append("Job %i/%i running: %s"%(i,k,path))
                else:
                    message.append("Job %i/%i waiting: %s"%(i,k,path))
        else:
            result = False
            message.append("All %i jobs finished!"%len(self.Jobs))
        if full_output:
            return result, message
        else:
            print os.linesep.join(message)
            return result
        
        
    def LoadInputFile(self, fpath):
        """
            Method to read an existing input file.

            Input:
            ------
                fpath: string
                    Path to the input file.
        """
        with open(fpath, "r") as fh:
            content = fh.readlines()
        
        content = filter(lambda s: not s.startswith("!"), content)
        content = map(lambda s: s.split("!")[0], content)
        content = map(str.strip, content)
        content = filter(lambda s: bool(s.strip()), content)
        
        
        NewParam = dict()
        keyw = None
        while content:
            line = content.pop(0)
            lline = line.lower()
            if lline == "filout":
                self.path_out = content.pop(0)
                self.bavfile = self.path_out + "_bav.txt"
            elif lline == "spgroup":
                self.sg = content.pop(0)
            elif lline.startswith("crystal") or lline.startswith("molecule"):
                keyw = lline.capitalize()
                NewParam[keyw] = []
            elif lline=="end":
                break
            else:
                nextkeyw = keyword_exists(line)
                if nextkeyw:
                    keyw = nextkeyw
                    NewParam[keyw] = []
                elif keyw!=None:
                    NewParam[keyw].append(line)
        
        structure_keyw = ["Crystal", "Crystal_t", "Crystal_p",
                          "Molecule", "Molecule_t"]
        found_structure = map(NewParam.has_key, structure_keyw)
        if sum(found_structure) > 1:
            raise ConsistencyError(
                "Several structure definitions found in input-file. "\
                "Structure can be specified in input-file by the "\
                "keywords Crystal or Molecule, not both!")
        elif sum(found_structure) == 1:
            ind = found_structure.index(True)
            structure = NewParam.pop(structure_keyw[ind])
            self.cif = False
            self.crystal = bool(ind < 3)
        elif not any(found_structure):
            raise IOError("No structure has been defined in input-file")
            
        self.LoadedParam = NewParam

        if all(map(NewParam.has_key, ["Atom", "Atom_conf"])):
            raise ConsistencyError(
                "Only one of the keywords ``Atom'' and "\
                "``Atom_conf'' may be given in the input file!")
            
        for keyw in NewParam.iterkeys():
            if not len(NewParam[keyw]):
                self.P[keyw] = True
                continue
            Values = NewParam[keyw]
            GroupName = self.GroupMembers[keyw]
            Group = getattr(self.Defaults, GroupName)
            Type = type(Group[keyw])
            for j in range(len(Values)):
                Value = Values[j]
                Value = Value.split() if Type in [list, tuple] else [Value]
                for i in range(len(Value)):
                    try:
                        Value[i] = int(Value[i])
                    except:
                        pass
                    try:
                        Value[i] = float(Value[i])
                    except:
                        pass
                if len(Value) > 1:
                    Value = tuple(Value)
                else:
                    Value = Value[0]
                Values[j] = Value
            if Type is not list:
                assert len(Values) == 1, "Many lines of input given for "\
                    "Parameter %s. Only one line of input accepted!"%keyw
                Values = Type(Values[0])
            self.P[keyw] = Values
            
        if not ind == 2:
            cell = map(float, structure.pop(0).split())
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma = cell
            if self.P.has_key("Atom"):
                Z = lambda i: self.P.Atom[i][0]
            else:
                Z = lambda i: i
            
            for Atom in structure:
                Atom = Atom.split()
                symbol = elements.symbols[Z(int(Atom[0]))]
                num = self.elements.values().count(symbol) + 1
                label = symbol + str(num)
                position = tuple(map(float, Atom[1:4]))
                occ = float(Atom[4]) if (len(Atom)>4) else 1
                self.add_atom(label, position, occupancy = occ)
            
            
            
        
    def get_XANES(self, conv = True):
        """
            Method to read-out the calculated XANES datas of the output file 
            created by FDMNES.
            
            Input:
            -----
                conv: bool
                    If the calculation with or without convolution is wanted.
            
        """
        if conv == True:
            fname = self.path_out + "_conv.txt" 
            skiprows = 1
        else:
            fname = self.path_out + ".txt"
            skiprows = 4
        
        data = np.loadtxt(fname, skiprows = skiprows)
        return data
        
        
    def get_dafs(self, Reflex, pol_in, pol_out, azimuth, conv = True):
        """
            Method to read-out the calculated DAFS datas. It's possible to get 
            the DAFS data's for given paramters with and without convolusion.
            If the given parameters doesn't exist, it is possible to add these 
            parameter to the input file.
            
            Input:
            ------
                Reflex: string
                    reflection indexes
                pol_in and pol_out: int
                    orientation of the polarization
                azimuth: int
                    azimuthal angle
                conv: bool
                    convoluted data
        """
        
        self.Reflex = self.Reflex.round(0)
    
        if conv == True:
            fname = os.path.abspath(self.path_out) + "_conv.txt"  
            skiprows = 1  
            a = 0

        else:
            fname = os.path.abspath(self.path_out) + ".txt"
            skiprows = 4
            a = 3
                   
        fobject = open(fname,"r") 
        lines = fobject.readlines()
        fobject.close()
        headline = lines[a]
        headline_keywords = headline.split()
            
        Reflex_a = np.array(Reflex)
        if np.ndim(self.Reflex)>1:
            columns = (self.Reflex[:,0:3] == Reflex).all(1)
        
            if pol_in != None:
                col_pol_in = (self.Reflex[:,3] == pol_in)
                columns *= col_pol_in
            if pol_out != None:
                col_pol_out = (self.Reflex[:,4] == pol_out)
                columns *= col_pol_out
            if azimuth != None:
                col_azimuth = (self.Reflex[:,5].round(0) == np.round(azimuth))
                columns *= col_azimuth
                
        else:
            columns = (self.Reflex[0:3] == Reflex).all(0)
            
            col_pol_in = (self.Reflex[3] == pol_in)
            columns *= col_pol_in
            
            col_pol_out = (self.Reflex[4] == pol_out)
            columns *= col_pol_out
            
            col_azimuth = (self.Reflex[5].round(0) == np.round(azimuth))
            columns *= col_azimuth
            
        columns_a = np.where(columns)[0]
        
        if len(columns_a) == 0:
            Reflex_new = np.append(Reflex_a, np.array([pol_out, pol_in, azimuth]))
            self.Reflex = Reflex_new
            return self.Reflex

        else:
            if conv == True:
                self.column = columns_a+2
                self.index = headline_keywords[self.column]
            else:
                column = (2*columns_a+2)
                index = headline_keywords[column]
                index_a = index.replace("r","I")
                self.index = index_a+" without convolution"
                
                self.column_real = column
                self.index_real = headline_keywords[self.column_real]
                
                self.column_im = column + 1
                self.index_im = headline_keywords[self.column_im]

            
        dafs = np.loadtxt(fname, skiprows = skiprows)
        return dafs
        
    def do_convolution(self, path = None):
        """
            Method to write a special convolution file. It is also possible to 
            define the convolution parameters in the WriteInputFile method.
            
            Input:
            ------
                path: string
                    Name of the convolution file.
        """
        if path == None:
            path = "_conv.".join(self.path.split("."))
            
        try: f = open(path, "w")
        except IOError:
            print "Error: No new file open"
        
        else: 
            print "Written content in the Conv_file successfully"
            
            
            if hasattr(self, "Cal") and isinstance(self.Cal, dict):
                key = self.Cal.keys()
                for element in self.Cal:
                    f.write("Calculation\n%s\n" %element)
                    value = self.Cal[element]
                    val = array2str(value, precision=1)
                   # value = self.Cal_val.values()   
                    f.write("%s\n" %val)
            elif hasattr(self, "Cal"): 
                f.write("Calculation\n%s\n" %self.Cal)
            else:
                pass
            #print("Using previous calculation in %s"%cal_path) 
            
            conv_out = self.path_out + "_conv.txt"
            f.write("\nConv_out\n%s\n" %conv_out)
                
            if hasattr(self, "Exp") and len(self.Exp)>0:
                f.write("\nExperiment\n")
                key = self.Exp.keys()
                for element in self.Exp:
                    exp_path = os.path.abspath(element) 
                    f.write(" %s\n" %exp_path)
                    value = self.Exp[element]
                    val = array2str(value, precision=1)
                   # value = self.Cal_val.values()   
                    f.write(" %s" %val)
                    
            f.write("\nConvolution \n\n")

            for key in conv_flags:
                if hasattr(self,key):
                    value =  getattr(self,key)
                    if isinstance(value, bool):
                        f.write("%s\n\n" %key)
                    else:
                        f.write("%s\n %s\n\n" %(key, getattr(self, key))) 
                        
            if hasattr(self, "Thompson") and len(self.Thompson)>0:
                self.Thompson = array2str(self.Thompson)
                f.write("\nThompson\n % \n" %self.Thompson)
                
            if hasattr(self, "Table") and len(self.Table)>0:
                self.Table = array2str(self.Table)  
                f.write("\nTable\n %s\n" %self.Table)

            f.write("\nEnd")
                
        finally: 
                f.close()

