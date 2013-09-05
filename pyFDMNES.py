# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 17:03:33 2013

@author: weiget
"""

import numpy as np
import StringIO
import logging
import os
import elements
import subprocess
import string


flags = ["Green","Magnetism", "Density","Self_absorption","cartesian",
         "nodipole","Memory_save"] 

string_flag = ["Hubbard","Edge"]         
         
conv_flags = ["Gamma_fix", "Fprime", "Fprime_atom", "Estart","Efermi",
              "check_conv","Gen_shift", "S0_2", "Selec_core", "Photoemission",
              "Forbidden","Gaussian", "Seah","Gamma_var", "Gamma_max", "Elarg",
              "Ecent", "Gamma_hole","Dec", "Memory_save"]

SCF_flags = ["N_self","P_self","R_self", "Delta_E_conv", "SCF_exc",
             "SCF_mag_free"]
        
Multi_Exp_flags = ["Quadrupole","Octupole","Dimag", "E1E2","E1E3","E2E2","E3E3",
                   "E1M1","M1M1","No_E1E1","No_E2E2","No_E1E2","No_E1E3"]

EXEfile = os.path.abspath('\\\\win.desy.de\\home\\weiget\\My Documents\\fdmnes_2013_05_27\\fdmnes\\fdmnes.exe')

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


class pyFDMNES(object):
    """ Python programm  for FDMNES application.
    """
    
    def __init__(self, structure, resonant="", verbose=False):
        """
            Calculation of spectra for x-ray absorption spectroscopy with FDMNES for
            given parameters and structures in the following steps:
                - compile input file for FDMNES 
                - process and controll the FDMNES calculation
                - fit and plot the calculate datas

            Input parameters:
                structure : either
                                - metric (cell dimensions) of the crystal
                            or
                                - path to .cif-file
        """
        self.positions = {} # symmetric unit
        self.elements = {}
        self.verbose = verbose
        self.resonant = []
        self.atom_num = {}
        self.occupancy = {}
        self.Crystal = True
        self.Radius = 2.5
        self.Range = np.array([-10., 0.1, -2, 0.2, 0., 0.5, 20., 1., 40.])
        self.Polarise = []
        self.Reflex = []
        self.Atom = {}
        self.dafs = []
        self.Rpotmax = 0
        self.extract = False
        self.Absorber = ()
       # self.convolution = True
        self.Exp = {}
      #  self.Scan = []
       # self.Scan_conv = []
        
        if str(structure).isdigit():
            int(structure)
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma = structure
            self.cif = False
        elif os.path.isfile(structure):
             end = os.path.splitext(structure)[1].lower()
             if end == ".cif":
                 self.load_cif(structure, resonant)
             elif end == ".txt":
                 self.Filein(structure)
        
    
    def add_atom(self, label, position):
        """
            Method to give parameters for the FDMNES calculation.
            
            Inputs:
            -------
                label : string
                    The label of the atom.
                    It has to be unique and to start with the symbold of the
                    chemical element.
                position : iterable of length 3, string
                    Postion of the atoms in the structure
        """
        if type(label) is not str: raise TypeError("Invalid label. Need string.")
        if len(position) is not 3: raise TypeError("Enter 3D position object!")
        
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
            raise ValueError("Atom label shall start with the symbol of the chemical element" 
                             "Chemical element not found in %s"%label)                    
    
        self.atom_num[label] = elements.Z[self.elements[label]]
        self.positions[label] = position
        
    def load_cif(self, fname, resonant=""):
        """ Method to loads a structure from a CIF file.
            Read-out parameters are space group number, spac egroup name, 
            space group origin, metric and atom positions.
         
            Input:
            ------
            fname: string
                The name of the CIF file.
            resonant: string
                Symbol of the resonant atom.

            Informations about CIF file:
            -----------------------------
                Hall SR, Allen FH, Brown ID (1991).
                "The Crystallographic Information File (CIF): a new standard 
                archive file for crystallography".
                Acta Crystallographica A47 (6): 655-685
        """
        fobject = open(fname,"r")
        lines = fobject.readlines()
        fobject.close()
        lines.reverse()
        num_atom = 0
        while lines:
            Line = lines.pop()
            Line = Line.replace("\t", " ")
            line = Line.lower()
            if line.startswith("_symmetry_int_tables_number"):
                self.sg_num = int(Line.split()[1])
            elif line.startswith("_symmetry_space_group_name_h-m"):
                self.sg_name = Line.split("'")[1].replace(" ","")
                if ":" in self.sg_name:
                    self.cscode = ":" + self.sg_name.split(':')[1]
            elif line.startswith("_cell_length_a"):
                self.a = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_cell_length_b"):
                self.b = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_cell_length_c"):
                self.c = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_cell_angle_alpha"):
                self.alpha = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_cell_angle_beta"):
                self.beta = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_cell_angle_gamma"):
                self.gamma = float(Line.split()[1].partition("(")[0]) #)
            elif line.startswith("_atom_site"):
                if line.startswith("_atom_site_label"): col_label = num_atom
                elif line.startswith("_atom_site_type_symbol"): col_symbol = num_atom
                elif line.startswith("_atom_site_fract_x"): col_x = num_atom
                elif line.startswith("_atom_site_fract_y"): col_y = num_atom
                elif line.startswith("_atom_site_fract_z"): col_z = num_atom
                num_atom+=1
            elif num_atom>0 and len(Line.split())==num_atom:
                atomline = Line.split()
                label = atomline[col_label]
                symbol = atomline[col_symbol]
                if symbol[:2].isalpha(): 
                   symbol = symbol[:2]
                else: symbol = symbol[:1]
                if symbol==resonant:
                    self.resonant.append(label)
                px = float(atomline[col_x].partition("(")[0]) #)
                py = float(atomline[col_y].partition("(")[0]) #)
                pz = float(atomline[col_z].partition("(")[0]) #)
                position = (px, py, pz)       
                if logging.DEBUG: 
                    self.add_atom(label, position)


    def Filout(self, path):
        """ Method write a input file for the FDMNES calculation.
            Edit the given structure parameter or parameter of the .cif file 
            to the new inpout file.
            More useful calculation parameters, for example Range and Radius
            have to given of the user. It is also possible to get parameter from 
            an exist input file and add to new parameters. For this use the 
            keyword "extract". 
         
            Input:
            ------
            path: string
                Name of the input file.
        """
        
        self.path = path
        
       # if os.path.exists(path):
            # print ("path exist!")
            # return
        
        try: f = open(path, "w")
        except IOError:
            print "Error: No new file open"
        
        else: 
            print "Written content in the file successfully"
            
            if "_inp" in path:
                path_part = os.path.splitext(path)[0]
                self.new_name = path_part.replace("_inp","_out")
                self.new_path = os.path.relpath(self.new_name, os.path.dirname(EXEfile))
                f.write("Filout \n %s \n\n" %self.new_path)
            else:
                path_basename = os.path.splitext(path)[0]
                self.new_name = path_basename + "_out"
                self.new_path = os.path.relpath(self.new_name, os.path.dirname(EXEfile))
                f.write("Filout \n %s \n\n" %self.new_path)
            
            if self.extract == False:
                self.bav = self.new_name + "_bav.txt"

        #if hasattr(self, "Range"):
            
            if hasattr (self,"Range") and not self.extract:
                if isinstance(self.Range, str):
                    f.write("Range \n %s \n" % self.Range)
                else:
                    self.Range = array2str(self.Range)
                    f.write("Range \n %s \n" % self.Range)
        
        #if self.radius == isinstance (float):
            if isinstance(self.Radius, float) and not self.extract:
                f.write("Radius \n %f \n" % self.Radius)
            
            if self.Rpotmax != 0 and not self.extract:
                f.write("\nRpotmax \n %s\n" %self.Rpotmax)
        
            for key in flags:
                if hasattr(self,key) and getattr(self,key) ==True:
                        f.write("\n%s \n" %key)
              
            for key in Multi_Exp_flags :
                 if hasattr(self,key) and getattr(self,key)==True:
                     f.write("\n%s \n" %key)
            
            for key in string_flag:
                if hasattr(self, key): 
                    value = getattr(self, key)
                    f.write("\n%s \n %s \n" %(key, value)) 
                    
                #if getattr(self,key) == True:
            if hasattr(self, "SCF") and self.SCF and not self.extract:
                f.write("\nSCF \n")
                for key in SCF_flags:
                    if hasattr(self, key):
                        value = getattr(self, key)
                        if isinstance(value, bool):
                            f.write("%s \n" %key)
                        else:
                            f.write("%s \n %f \n" %(key, getattr(self, key))) 
                
            if hasattr(self,"sg_num"):    
                f.write("\nSpgroup \n %i" %self.sg_num)
            if hasattr(self,"cscode"):
                f.write("%s\n\n" %self.cscode)
            else:f.write("\n")
                
            if hasattr(self,"Absorber") and len(self.Absorber)>0:
                if not isinstance(self.Absorber, str):
                    self.Absorber = array2str(self.Absorber)
                f.write("Absorber\n %s \n" %self.Absorber)

            if hasattr(self, "Atom") and len(self.Atom)>0:
                f.write("\nAtom\n")
                #if isinstance(self.Atom, dict):
                for label in self.Atom.keys():
                    f.write(" %i" %self.atom_num[label])
                    atm_conf = array2str(self.Atom[label])
                    f.write(" %s" %atm_conf)
          
            if self.Crystal == True: f.write("\nCrystal \n")
            else: f.write("Molecule \n")
            
            cell = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
            f.write(" %f %f %f %f %f %f\n" %cell)
                                                               
            for label in self.resonant:
                if len(self.Atom)>0:
                    num = self.Atom.keys().index(label) + 1
                    f.write(" %i " %num)
                    atm_pos = array2str(self.positions[label])
                    f.write("%s" %atm_pos)
                else:
                    if hasattr(self,"pos"):
                        if label[:2].isalpha(): 
                            symbol = label[:2]
                        else: symbol = label[:1]
                        if self.atm_num == elements.Z[symbol] and len(self.Atom) == 0:
                            for pos in self.pos:
                                num, position = pos
                                x, y, z = position
                                f.write("%s %s %s %s\n" %(num, x, y, z))
                    else:    
                        f.write(" %i " %self.atom_num[label])
                        atm_pos = array2str(self.positions[label])
                        f.write("%s" %atm_pos)
               
            for label in self.positions.iterkeys():
                if label in self.resonant:
                    pass
                else:
                    if len(self.Atom)>0:
                        num = self.Atom.keys().index(label) + 1
                        f.write(" %i " %num)
                        atm_pos = array2str(self.positions[label])
                        f.write("%s" %atm_pos)
                    else:
                      if hasattr(self,"pos"):
                          if label[:2].isalpha(): 
                               symbol = label[:2]
                          else: symbol = label[:1]
                          if self.atm_num == elements.Z[symbol] and len(self.Atom) == 0:
                            for pos in self.pos:
                                num, position = pos
                                x,y,z = position
                                f.write("%s %s %s %s\n" %(num, x, y, z))
                      else: 
                            f.write(" %i " %self.atom_num[label])
                            atm_pos = array2str(self.positions[label])
                            f.write("%s" %atm_pos)
              
    
            if self.extract == True and os.path.exists(self.bav):
                bav_file = os.path.abspath(self.bav)
                f.write("\nextract\n %s\n" %bav_file)
                
                if hasattr(self,"ext_Absorber") and len(self.ext_Absorber)>0:
                    f.write("Absorber\n%s \n\n" %self.ext_Absorber)
                    f.write("Extractpos\n%s \n\n"%self.Absorber)
                    
                if hasattr(self,"Rotsup") and len(self.Rotsup)>0:
                   f.write("Rotsup\n%s \n\n" %self.Rotsup)
                
                if hasattr(self,"Extractsym") and len(self.Extractsym)>0:
                   f.write("Extractsym\n%s \n\n" %self.Extractsym)
                 
                if hasattr(self,"Reflex") and len(self.Reflex)>0:
                    Reflex_val = array2str(self.Reflex, precision=0)
                    f.write("\nRXS\n%s\n" %Reflex_val)
                
                if hasattr(self,"Polarise") and len(self.Polarise)>0:
                    Polarise_val = array2str(self.Polarise)
                    f.write("\nPolarise\n%s\n" %Polarise_val)
                    
                if hasattr(self,"dafs") and len(self.dafs)>0:
                    dafs_val = array2str(self.dafs)
                    f.write("\nAtom\n%s \n" %dafs_val)
                    
     #       if self.convolution == True:
     #           f.write("\nConvolution\n")
     #           
     #       if hasattr(self,"Efermi"):
     #           f.write("Efermi \n %f\n\n" %self.Efermi)
     #           
     #       if hasattr(self,"Estart"):
     #           f.write("Estart \n %f\n\n" %self.Estart)
            
            f.write("\nEnd")
    
       # if self.R_self != self.Radius: oder 3.5???
            # f.write("...")
       #if self.N_self != 30:
            #f.write("...")

        finally:
                f.close()
            
    def FDMNESfile (self):
        """ Method to start the FDMNES calculation of the createt input file.
            During the calculation it isn't possible to work. If the calculation 
            finished, "FDMNES simulation finished" is printing.
        """
        
        assert hasattr(self, "path"), "Attribute `path` has not been defined. Try self.Filout() first"
        self.workdir = os.getcwd()
        
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
        
    def retrieve (self):
        """ Method to check the excisting of the created BAV file.
        """
        if (not hasattr(self, "proc") or self.proc.poll() != None) and os.path.exists(self.bav):
            bav_open  = open(self.bav)
            line = bav_open.readlines()
            if ("Have a beautiful day !") in line[-1]:
                print("Process was sucessfully")
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
        self.atom_num = {}
        for atom in positions:
            num = atom[0]
            position = atom[1:]
            symbols = []
            if "Atom" in content_red:
                symbol = Atom[num-1][0]
                self.atom_num[symbol] = elements.Z[symbol[:-1]]
            else:
                symbol = elements.symbols[num]
                self.atom_num[symbol] = elements.Z[symbol]
                symbols.append(symbol)
                anzahl = symbols.count(symbol)
                symbol += str(anzahl)

            self.positions[symbol] = position
            self.atom_num
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
        
    def do_convolution(self, path = None):
        """
            Method to write a special convolution file. It is also possible to 
            define the convolution parameters in the Filout method.
            
            Input:
            ------
                path: string
                    Name of the convolution file.
        """
        if path == None:
            path = "_conv.".join(self.path.split("."))
            
        self.path = os.path.abspath(path)
        EXE = os.path.dirname(EXEfile)
        rel_conv_name = os.path.relpath(self.path,EXE)
        
        try: f = open(path, "w")
        except IOError:
            print "Error: No new file open"
        
        else: 
            print "Written content in the Conv_file successfully"
            
            
            if hasattr(self, "Cal") and isinstance(self.Cal, dict):
                key = self.Cal.keys()
                for element in self.Cal:
                    cal_path = os.path.abspath(element) 
                    Calculation = os.path.relpath(cal_path,EXE)
                    f.write("Calculation\n%s\n" %Calculation)
                    value = self.Cal[element]
                    val = array2str(value, precision=1)
                   # value = self.Cal_val.values()   
                    f.write("%s\n" %val)
            elif hasattr(self, "Cal"): 
                cal_path = os.path.abspath(self.Cal) 
                Calculation = os.path.relpath(cal_path,EXE)
                f.write("Calculation\n%s\n" %Calculation)
            else:
                cal_path = os.path.relpath(self.new_name + ".txt", EXE)
                f.write("Calculation\n%s\n" %cal_path)
            #print("Using previous calculation in %s"%cal_path) 
            
            if hasattr(self,"Scan") and len(self.Reflex)>0 and len(self.Reflex[0,:])<6:
                for element in self.Scan:
                    scan_path = os.path.abspath(element) 
                    Scan = os.path.relpath(scan_path,EXE)  
                    f.write("\nScan\n%s\n" %Scan)
            elif hasattr(self,"Scan"):
                Scan = os.path.relpath(self.new_name + "_scan.txt", EXE)     
                f.write("\nScan\n%s\n" %Scan)
              
            if hasattr(self,"Scan_conv") and path!= None:
                for element in self.Scan_conv:
                    scan_conv_path = os.path.abspath(element) 
                    Scan_conv = os.path.relpath(scan_conv_path,EXE)
                    f.write("\nScan_conv\n%s\n" %Scan_conv)
            elif hasattr(self,"Scan_conv"):
                Scan_conv = os.path.relpath(self.new_name + "_scan_conv.txt", EXE)
                f.write("\nScan_conv\n%s\n" %Scan_conv)
        
            if hasattr(self,"Conv_out"):
                name = "\\" + self.Conv_out
            else:
                name = self.new_name + "_conv.txt"
            Conv_out = os.path.relpath(name, EXE)
            f.write("\nConv_out\n%s\n" %Conv_out)
                
            if hasattr(self, "Exp") and len(self.Exp)>0:
                f.write("\nExperiment\n")
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
                    setattr(self, varname, bool(getattr(self,varname)))
                
             #   if varname in vars_string:
                #    setattr(self, varname, float(getattr(self,varname)))
       
       