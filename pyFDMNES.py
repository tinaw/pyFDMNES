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


class ConsistencyError(Exception):
    def __init__(self, value="Inconsistent parameter values found.", 
                       errmsg="", identifier=None):
        self.value = value
        self.errmsg = errmsg
    def __str__(self):
        return self.value + os.linesep + "Message: %s" %self.errmsg


class Parameters(dict):
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError, 'expected dict'

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)
    __setattr__ = __setitem__
    __getattr__ = __getitem__



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
        self.elements = {}
        self.verbose = verbose
        self.resonant = []
        self.Z = {}
        self.occupancy = {}
        self.Crystal = True
        self.Polarise = []
        self.Atom = {}
        self.extract = False
        self.Absorber = ()
        # self.convolution = True
        self.Exp = {}
        self.P = Parameters()
        self._WD = os.getcwd()
        
        
        ### LOAD STRUCTURE #################################################
        if str(structure).isdigit():
            self.sg_num = int(structure)
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
             elif end == ".txt":
                 self.Filein(structure)
        
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
    
    def add_atom(self, label, position, resonant=False):
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
            self.sg_num = mkfloat(cb["_symmetry_int_tables_number"])
        
        if cb.has_key("_symmetry_space_group_name_h-m"):
            self.sg_name = cb["_symmetry_space_group_name_h-m"]
            self.sg_name = self.sg_name.replace(" ", "")
            i = self.sg_name.find(":")
            if i>=0:
                self.cscode = self.sg_name[i:]
        
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
            position = (px, py, pz)
            self.add_atom(label, position, (label in self.resonant))
            if symbol == resonant:
                self.P.Z_absorber = elements.Z[symbol]
        
    
    def check_parameters(self, keyw):
        """
            Checks configured FDMNES parameters for consistency.
        """
        if not self.P.has_key(keyw):
            return True
        if keyw=="Atom" and len(self.P.Atom):
            numspecies = len(np.unique(self.elements.values()))
            numconf = len(self.P.Atom)
            if not numconf == numspecies:
                raise ConsistencyError(errmsg="Number of given electron "
                    "configurations (%i) does not match number of different "
                    "species in structure (%i)"%(numconf,numspecies))
            if self.P.has_key("Atom_conf") and len(self.P.Atom_conf):
                raise ConsistencyError(errmsg="Electronic configuration"
                    "given twice: Parameters Atom and Atom_conf")
        return True

    def _skip_group(self, Group, TestGroup="all"):
        if TestGroup=="all":
            TestGroup = ["SCF", "Convolution", "Extract", "Reflection"]
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
        if hasattr(self,"sg_num"):    
            sg = "Spgroup%s%i"%(os.linesep, self.sg_num)
        if hasattr(self,"cscode"):
            sg += os.linesep + self.cscode
        
        output.append(sg)
        output.append(os.linesep)
        
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
        
        self.check_parameters("Atom")
        for label in self.positions.iterkeys():
            if len(self.P.Atom):
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
    
    def FileOut(self, path, overwrite=False):
        """ 
            Method writes an input file for the FDMNES calculation.
            The calculation parameters, for example Range and Radius etc., 
            are taken from the fdmnes.P NamedTuple. 
         
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
        output.append(os.linesep)
        
        output.append("Folder_dat")
        output.append(self.fdmnes_dir)
        output.append(os.linesep)
        
        self.bavfile = self.path_out + "_bav.txt"
        
        output.extend(self.write_structure())
        
        for Group in settings.Defaults:
            if self._skip_group(Group):
                continue
            for keyw in Group.iterkeys():
                if self.P.has_key(keyw) and self.P[keyw]!=Group[keyw]:
                    deftype = type(Group[keyw])
                    assert isinstance(self.P[keyw], deftype),\
                        "Wrong type: Parameter %s has to be of type %s"\
                            %(keyw, deftype)
                    self.check_parameters(keyw) # t.b.w.
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
        
    def FDMNESfile(self):
        """
            Method to start the FDMNES calculation of the created input file.
        """
        
        assert hasattr(self, "path"), \
            "Attribute `path` has not been defined. "\
            "Try fdmnes.FileOut() first."
        
        fdmfile = "fdmfile.txt"
        
        output = []
        if hasattr(self,"conv_path"):
            Run_File = 2 
            output.append("%i"%Run_File)
            output.append("")
            output.append("%s"%self.path)
            output.append("%s"%self.conv_path)
        else:
            Run_File = 1
            output.append("%i"%Run_File)
            output.append("")
            output.append("%s"%self.path)
        
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
            self.proc.wait()
            print("FDMNES simulation finished")
        except Exception as e:
            print("An error occured when running FDMNES executable at")
            print(self.fdmnes_exe)
            print("Message: %s"%e)
        finally:
            logfile.close()
        
    def retrieve(self):
        """ Method to check the excisting of the created BAV file.
        """
        if (not hasattr(self, "proc") or self.proc.poll() != None) \
            and os.path.exists(self.bavfile):
            bav_open  = open(self.bavfile)
            line = bav_open.readlines()
            if ("Have a beautiful day !") in line[-1]:
                print("Process finished successfully")
        else: print("Process wasn't sucessfully")
            
            
    def Filein(self, fname):
        """
            Method to readout exist input files. It is possible to get 
            informations about the Radius, Range, Atom configuration and Crystal
            structure and also different calculation parameters.
            

            Input:
            ------
                fname: string
                    Name of the input file.
        """
        fobject = open(fname, "r")
        content = fobject.read()
        fobject.close()
        
        content_red = map(string.strip,content.splitlines())
        content_red = filter(lambda x: x is not "", content_red)
        
        keywords = ["SCF", "SCF_exc", "SCF_mag_free","Memory_save",
                    "Density"]
        keywords_float = ["Radius", "Efermie", "Estart","Rpotmax",
                          "N_self","P_self","R_self", "Delta_E_conv","Hubbard"]
        keywords_string = ["Absorber","Edge","Range"]
                
        self.path_out = os.path.splitext(fname)[0].replace("_inp","_out")
        self.bavfile = self.path_out + "_bav.txt"

       # if "Range" in content_red:
           #     self.Range = content_red[content_red.index("Range")+1]
                
        for line in content_red:
            if line in flags + Multi_Exp_flags + keywords:
                 setattr(self, line, True)
            
            if line in keywords_float:
                setattr(self, line, float(content_red[content_red.index(line)+1]))
                
            if line in keywords_string:
                setattr(self, line, str(content_red[content_red.index(line)+1]))
                
            if line == "Spgroup":
                line = content_red[content_red.index(line)+1]
                if ":" in line:
                    line.split(':')
                    self.sg_num = int(line[0])
                    self.cscode = ":" + line[1]
                else:
                    self.sg_num = int(line)
                
            if line == "Atom":
                Atom = []
                linenum = content_red.index(line)
                atom_num_list = []
                while True:
                    try:
                        linenum += 1
                        atomline = content_red[linenum]
                        atomline = atomline.split()
                        atomline = map(float, atomline)
                        atom_num = int(atomline[0])
                        atom_num_list.append(atom_num)
                        anzahl = atom_num_list.count(atom_num)
                        symbol = elements.symbols[atom_num] + str(anzahl)
                        atm_conf = atomline[1:]
                        Atom.append((symbol,atm_conf))
                    except ValueError:
                        break
                self.Atom = dict(Atom)
               # while :
                # self.atm_conf = content_red[content_red.index(line)+1]
                
            if line in ["Crystal", "Molecule"]:
                if line=="Crystal":
                    self.crystal=True
                linenum = content_red.index(line)+1
                cellline = content_red[linenum]
                cellline = cellline.split()
                cellline = map(float, cellline)
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma = cellline
                positions = []
                while True:
                    try:
                        linenum += 1
                        atomline = content_red[linenum]
                        atomline = atomline.split()
                        atomline = map(float, atomline)
                        num = int(atomline[0])
                        position = atomline[1:4]
                        #position = array2str(position)
                        positions.append((num, position))
                    except ValueError:
                        break
                #self.pos = array2str(positions)
        if "Atom" in content_red :
            pass
        else:
            self.pos = positions
            self.atm_num = num
        
        self.positions = {}
        self.Z = {}
        for atom in positions:
            num = atom[0]
            position = atom[1:]
            symbols = []
            if "Atom" in content_red:
                symbol = Atom[num-1][0]
                self.Z[symbol] = elements.Z[symbol[:-1]]
            else:
                symbol = elements.symbols[num]
                self.Z[symbol] = elements.Z[symbol]
                symbols.append(symbol)
                anzahl = symbols.count(symbol)
                symbol += str(anzahl)

            self.positions[symbol] = position
      #  self.Atom = dict(Atom)
             #   self.pos = positions
             #   return self.pos
            
       #         self.Crysatal = line
        #        self.cell = content_red[content_red.index(line)+1]
         #       while 
                
                                  
        #return content_red

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
            define the convolution parameters in the FileOut method.
            
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

