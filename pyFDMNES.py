# -*- coding: utf-8 -*-F
"""
Created on Thu Jul 25 17:03:33 2013

@author: weiget
"""

import numpy as np
import StringIO
import os
import elements
import subprocess
import string
import collections
from settings import Defaults, ParamTypes




EXEfile = os.path.abspath('/space/crichter/backup/Programme/FDMNES/fdmnes_2013_12_12/fdmnes_linux64')


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



class pyFDMNES(object):
    """ 
        Python programm  for FDMNES application.
    """
    
    def __init__(self, structure, resonant="", verbose=False, 
                       fdmnes_exe=EXEfile):
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
        self.convolution = True
        self.Cal = {}
        self.Exp = {}
        self.Scan = []
        self.Scan_conv = []
        self.P = Parameters()
        self._WD = os.getcwd()
        
        if str(structure).isdigit():
            int(structure)
            self.a, self.b, self.c, \
            self.alpha, self.beta, self.gamma = structure
            self.cif = False
        elif os.path.isfile(structure):
             end = os.path.splitext(structure)[1].lower()
             if end == ".cif":
                 self.load_cif(structure, resonant)
             elif end == ".txt":
                 self.Filein(structure)
        
        self.fdmnes_dir = os.path.dirname(fdmnes_exe)
        self.fdmnes_exe = fdmnes_exe
    
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
                    Postion of the atoms in the structure
        """
        if type(label) is not str:
            raise TypeError("Invalid label. Need string.")
        if len(position) is not 3:
            raise TypeError("Enter 3D position object!")
        
        position = tuple(position)
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
        try:
            cf = CifFile.ReadCif(path)
            cb = cf.first_block()
        except Exception as e:
            print("File doesn't seem to be a valid .cif file: %s"%path)
            print e
        
        import CifFile
        self.Crystal = True
        # Reset Structure:
        self.positions.clear()
        self.P.Atom = []
        self.P.Absorber = ()
        self.P.Z_absorber = 0
        
        if cb.has_key("_symmetry_int_tables_number")
            self.sg_num = mkfloat(cb["_symmetry_int_tables_number"])
        
        if cb.has_key("_symmetry_space_group_name_h-m")
            self.sg_name = mkfloat(cb["_symmetry_space_group_name_h-m"])
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
        
    
    def check_parameters(self):
        if self.P.has_key("Atom"):
            numspecies = len(np.unique(self.elements.values()))
            numconf = len(self.P.Atom)
            assert numconf == numspecies,
             "Number of given electron configurations (%i) does not "%numconf\
             "match number of different species in structure (%i)"%numspecies
            

    def FileOut(self, path, overwrite=False):
        """ Method writes an input file for the FDMNES calculation.
            Further calculation parameters, for example Range and Radius
            have to given of the user. It is also possible to get parameter
            from an exist input file and add to new parameters. 
            For this use the keyword "extract". 
         
            Input:
            ------
            path : string
                Path to of the input file.
            overwrite : bool
                Has to be True to overwrite existing input file.
        """
        if os.path.isfile(path) and not overwrite:
            raise IOError("File %s already exists. "%path
                          "Use overwrite=True to replace it.")
        
        self.path = path
        
        output = []
        
        basepath = os.path.splitext(path)[0]
        if "_inp" in path:
            self.path_out = basepath.replace("_inp","_out")
        else:
            self.path_out = basepath + "_out"
            self.path_out = os.path.relpath(path_out, self.fdmnes_dir)
        output.append("Filout")
        output.append(self.path_out)
        
        self.bavfile = self.path_out + "_bav.txt"
        
        
        if hasattr(self,"sg_num"):    
            sg = "Spgroup%s%i"%(os.linesep, self.sg_num)
        if hasattr(self,"cscode"):
            sg += os.linesep + self.cscode
        else:
            output.append(sg)
        
        
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
        
        for label in self.positions.iterkeys():
            if len(self.P.Atom) and not len(self.P.Atom_conf):
                atomnum = self.P.Atom.index(self.Z[label]) + 1
            else:
                atomnum = self.Z[label]
            pos = self.positions[label]
            line = (atomnum,) + pos
            if suffix != "":
                occ = self.occupancy[label]
                line += (occ,label)
                output.append("%i %.10g %.10g %.10g %.g !%s"%line)
            else:
                line += (label,)
                output.append("%i %.10g %.10g %.10g !%s"%line)

        for Group in Defaults:
            if Group.has_key("SCF") and not Group["SCF"]:
                continue
            if Group.has_key("Convolution") and not Group["Convolution"]:
                continue
            if Group.has_key("Extract") and not Group["Extract"]:
                continue
            if Group.has_key("Reflection") and not len(Group["Reflection"])>0:
                continue
            for keyw in Group.iterkeys():
                if self.P.has_key(keyw) and self.P[keyw]!=Group[keyw]:
                    deftype = type(Group[keyw])
                    assert isinstance(self.P[keyw], deftype),
                        "Wrong type: Parameter %s has to be of type %s"\
                            %(keyw, deftype)
                    check_parameters(keyw, self.P[keyw]) # t.b.w.
                    output.append(keyw)
                    output.append(param2str(self.P[keyw]))

        with open(path, "w") as f:
            f.writelines(output)
        ########################################################################
        
    def FDMNESfile(self):
        """ Method to start the FDMNES calculation of the createt input file.
            During the calculation it isn't possible to work. If the calculation 
            finished, "FDMNES simulation finished" is printing.
        """
        
        assert hasattr(self, "path"), 
            "Attribute `path` has not been defined. Try self.FileOut() first"
        
        File = os.path.split(EXEfile)
        fdmfile =os.path.join(File[0],"fdmfile.txt")
        
        
        Input = os.path.relpath(self.path, os.path.dirname(EXEfile))
        #Input = os.path.relpath(self.path, os.getcwd())
        
        try: f = open(fdmfile,"w")
        except IOError:
            print "Error: fdmfile not open"
        
        else: 
            print "Written fdmfile was successfully"
            
            if hasattr(self,"conv_path"):
                Run_File = 2 
                Input_conv = os.path.relpath(self.conv_path, os.path.dirname(EXEfile))
                f.write("\n%i\n\n%s\n%s\n" %(Run_File,Input,Input_conv))
            else:
                Run_File = 1
                f.write("\n%i\n\n%s\n" %(Run_File,Input))
        
        finally: 
            f.close()
        
        os.chdir(File[0])
        logfile = open(os.path.join(self.workdir, "fdmnfile.txt"), "w")
        #errfile = open(os.path.join(self.workdir, "fdmnes_error.txt", "w"))
        if self.verbose:
            stdout = None
        else:
            stdout = logfile
        try:
            self.proc = subprocess.Popen(EXEfile, stdout=stdout) #stderr = errfile)
            self.proc.wait()
            print("FDMNES simulation finished")
        except:
            pass
        finally:
            logfile.close()
        os.chdir(self.workdir)
        
    def retrieve(self):
        """ Method to check the excisting of the created BAV file.
        """
        if (not hasattr(self, "proc") or self.proc.poll() != None) and os.path.exists(self.bav):
            bav_open  = open(self.bav)
            line = bav_open.readlines()
            if ("Have a beautiful day !") in line[-1]:
                print("Process finished successfully")
        else: print("Process wasn't sucessfully")
            
            
    def Filein (self, fname):
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
                
        self.new_name = os.path.splitext(fname)[0].replace("_inp","_out")
        self.bav = self.new_name + "_bav.txt"

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
            fname = self.new_name + "_conv.txt" 
            skiprows = 1
        else:
            fname = self.new_name + ".txt"
            skiprows = 4
        
        data = np.loadtxt(fname, skiprows = skiprows)
        return data
        
        
    def get_dafs (self, Reflex, pol_in, pol_out, azimuth, conv = True):
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
                    convolution information
        """
        
        self.Reflex = self.Reflex.round(0)
    
        if conv == True:
            fname = os.path.abspath(self.new_name) + "_conv.txt"  
            skiprows = 1  
            a = 0

        else:
            fname = os.path.abspath(self.new_name) + ".txt"
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
        
    def do_convolution(self, path):
        """
            Method to write a special convolution file. It is also possible to 
            define the convolution parameters in the FileOut method.
            
            Input:
            ------
                path: string
                    Name of the convolution file.
        """
        
        self.conv_path = os.path.abspath(path)
        EXE = os.path.dirname(EXEfile)
        rel_conv_name = os.path.relpath(self.conv_path,EXE)
        
        
        try: f = open(path, "w")
        except IOError:
            print "Error: No new file open"
        
        else: 
            print "Written content in the Conv_file successfully"
            
            
            if hasattr(self, "Cal") and isinstance(self.Cal, dict):
                key = self.Cal.keys()
                f.write("Calculation\n")
                for element in self.Cal:
                    cal_path = os.path.abspath(element) 
                    Calculation = os.path.relpath(cal_path,EXE)
                    f.write("%s\n" %Calculation)
                    value = self.Cal[element]
                    val = array2str(value, precision=1)
                   # value = self.Cal_val.values()   
                    f.write("%s\n" %val)
            else: 
                cal_path = os.path.abspath(self.Cal) 
                Calculation = os.path.relpath(cal_path,EXE)
                f.write("Calculation\n%s\n" %Calculation)
             
            if hasattr(self,"Scan") and len(self.Reflex)>0 and len(self.Reflex[0,:])<6:
                for element in self.Scan:
                    scan_path = os.path.abspath(element) 
                    Scan = os.path.relpath(scan_path,EXE)
                    f.write("\nScan\n%s\n" %Scan)
              
            if hasattr(self,"Scan_conv"):
                for element in self.Scan_conv:
                    scan_conv_path = os.path.abspath(element) 
                    Scan_conv = os.path.relpath(scan_conv_path,EXE)
                    f.write("\nScan_conv\n%s\n" %Scan_conv)
        
            if hasattr(self,"Conv_out"):
                name = "\\" + self.Conv_out
                Conv_out = os.path.dirname(rel_conv_name) + name
                f.write("\nConv_out\n%s\n" %Conv_out)
                
            if hasattr(self, "Exp") and len(self.Exp)>0:
                f.write("Experiment\n")
                key = self.Exp.keys()
                for element in self.Exp:
                    exp_path = os.path.abspath(element) 
                    Exp = os.path.relpath(exp_path,EXE)
                    f.write(" %s\n" %Exp)
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
        
    def checkvars(self):
        """
            Method to check the right attribut of the calculation parameters.
        """
        
        vars_float = ["N_self","P_self","R_self", "Delta_E_conv","Radius",
                      "Estart","Efermi"]
        vars_int = ["Run_Fil"]
        vars_bool = ["Quadrupole","Octupole","Dimag", "E1E2","E1E3","E2E2","E3E3",
                   "E1M1","M1M1","No_E1E1","No_E2E2","No_E1E2","No_E1E3","Green",
                   "Magnetism", "Density","Self-absorption","cartesian","nodipole",
                   "SCF_exc", "SCF_mag_free"]
        vars_string = ["Range","Absorber"]
        
        for varname in vars_float + vars_int + vars_bool + vars_string:
            if hasattr (self,varname):
                if varname in vars_float:
                    setattr(self, varname, float(getattr(self,varname)))

                if varname in vars_int:
                    setattr(self, varname, int(getattr(self,varname)))
                    
                if varname in vars_bool:
                    setattr(self, varname, float(getattr(self,varname)))
                
             #   if varname in vars_string:
                #    setattr(self, varname, float(getattr(self,varname)))
       
       
