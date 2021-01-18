from collections import namedtuple

nan = float("nan")

Groups = ["Basic", 
          "SCF", 
          "Multipole",
          "RXS",
          "Spin",
          "Convolution", 
          "Extract", 
          "Experiment",
          "Parameter",
          "Technical",
          "Fit"]


ParamTypes = namedtuple("ParamTypes", Groups)

Defaults = ParamTypes(*[{} for i in range(len(Groups))])

Defaults.SCF['Delta_E_conv'] = 1.
Defaults.SCF['N_self'] = 30
Defaults.SCF['P_self'] = 0.1
Defaults.SCF['SCF'] = False
Defaults.SCF['SCF_exc'] = False
Defaults.SCF['SCF_mag_free'] = False

Defaults.Basic["Absorber"] = ()
Defaults.Basic["Atom"] = []
Defaults.Basic["Atom_conf"] = []
Defaults.Basic["Cartesian"] = False
Defaults.Basic["Crystal_p"] = ""
Defaults.Basic["Density"] = False
Defaults.Basic["Density_all"] = False
Defaults.Basic["Edge"] = "K"
Defaults.Basic["Energpho"] = False
Defaults.Basic["Green"] = False
Defaults.Basic["Hubbard"] = ()
Defaults.Basic["Hubbard_z"] = [] # list of 2-tuples
Defaults.Basic["Magnetism"] = False
Defaults.Basic["Memory_save"] = False
Defaults.Basic["Polarize"] = []
Defaults.Basic['R_self'] = 0. # not only for SCF but also Fermi energy
Defaults.Basic["Radius"] = 5.67
Defaults.Basic["Range"] = () # eV with respect to edge
Defaults.Basic["Relativism"] = False
Defaults.Basic["RPA"] = False
Defaults.Basic["Rpotmax"] = 0.
Defaults.Basic["Screening"] = []
Defaults.Basic["TDDFT"] = False
Defaults.Basic["Z_absorber"] = 0
Defaults.Basic["Full_atom"] = False
#Defaults.Basic["fprim"] = False

Defaults.Multipole['Dipmag'] = False
Defaults.Multipole['E1E2'] =  False
Defaults.Multipole['E1E3'] =  False
Defaults.Multipole['E1M1'] =  False
Defaults.Multipole['E2E2'] =  False
Defaults.Multipole['E3E3'] =  False
Defaults.Multipole['M1M1'] =  False
Defaults.Multipole['No_E1E1'] =  False
Defaults.Multipole['No_E1E2'] =  False
Defaults.Multipole['No_E1E3'] =  False
Defaults.Multipole['No_E2E2'] =  False
Defaults.Multipole['Octupole'] =  False
Defaults.Multipole['Quadrupole'] =  False

Defaults.RXS['Azimuth'] = ()
Defaults.RXS['Circular'] = False
Defaults.RXS['Dead_layer'] = 0.
Defaults.RXS['Double_cor'] = False
Defaults.RXS['Full_self_abs'] = False
Defaults.RXS['RXS'] = []
Defaults.RXS['Reflection'] = ()
Defaults.RXS['Self_abs'] = False
Defaults.RXS['Step_azim'] = 2.
Defaults.RXS['Zero_azim'] = (0.,0.,0.) # 3-tuple float
Defaults.RXS['Selec_input'] = ()
Defaults.RXS['Selec_out'] = ()
Defaults.RXS['Energy'] = ()


Defaults.Spin["Axe_spin"] = (0.,0.,1.) # 3-tuple float
Defaults.Spin["Ang_spin"] = (0.,0.,0.) # 3-tuple float
Defaults.Spin["Magnetism"] = False
Defaults.Spin["Nonrelat"] = False
Defaults.Spin["Spinorbite"] = False


Defaults.Technical["Check"] = ()
Defaults.Technical["check_all"] = False #equivalent to icheck = 3 for all the subroutines
Defaults.Technical["no_check"] = False #equivalent to icheck = 0 for all the subroutines
Defaults.Technical["check_pot"] = False  #equivalent to icheck = 3 for all the subroutines concerning the potential
Defaults.Technical["check_mat"] = False  #equivalent to icheck = 3 for the matrices MSM or FDM
Defaults.Technical["check_sph"] = False  #equivalent to icheck = 3 for the subroutine ``sphere''
Defaults.Technical["check_coabs"] = False #equivalent to icheck = 3 for the subroutine ``coabs''

Defaults.Technical["Sym"] = "" # point symmetry of the resonant atom
Defaults.Technical["Nonexc"] = False # Make calculation for non excited absorbing atom for K, L1, M1 edges
Defaults.Technical["Excited"] = False # Make calculation for excited absorbing atom for all other edges
Defaults.Technical["Rmt"] = nan # muffin tin radius, around 0.65
Defaults.Technical["Iord"] = 4 # order of taylor expansion for fdm
Defaults.Technical["Adimp"] = 0.25 # inter point distance for fdm

Defaults.Technical["Lmax"] = -1 # number of spherical harmotics, negative values mean additional


Defaults.Convolution["Calculation"] = []
Defaults.Convolution["Check_conv"] = False
Defaults.Convolution["Convolution"] = False
Defaults.Convolution["Dec"] = False
Defaults.Convolution["Ecent"] = 30. 
Defaults.Convolution["Efermi"] = -5. 
Defaults.Convolution["Elarg"] = 30. 
Defaults.Convolution["Estart"] = nan 
Defaults.Convolution["Forbidden"] = False 
Defaults.Convolution["Fprime"] = False 
Defaults.Convolution["Fprime_atom"] = False 
#Defaults.Convolution["Gamma_fix"] = False # obsolete, use Gamma_var instead
Defaults.Convolution["Gamma_hole"] = -1. 
Defaults.Convolution["Gamma_max"] = 15. 
Defaults.Convolution["Gamma_var"] = False 
Defaults.Convolution["Gaussian"] = 0. 
Defaults.Convolution["Nocut"] = False
Defaults.Convolution["Photoemission"] = False
Defaults.Convolution["S0_2"] = 1.
Defaults.Convolution["Scan"] = [] 
Defaults.Convolution["Scan_conv"] = ""
Defaults.Convolution["Seah"] = (0., 0.)
Defaults.Convolution["Selec_core"] = -1
Defaults.Convolution["Thomson"] = nan
Defaults.Convolution["Xan_atom"] = False

Defaults.Parameter["Par_ecent"] = (0.,0.,0)
Defaults.Parameter["Par_elarg"] = (0.,0.,0)
Defaults.Parameter["Par_efermi"] = (0.,0.,0)
Defaults.Parameter["Par_gamma_hole"] = (0.,0.,0)
Defaults.Parameter["Par_gamma_max"] = (0.,0.,0)
Defaults.Parameter["Par_gauss"] = (0.,0.,0)
Defaults.Parameter["Par_shift"] = (0.,0.,0)
Defaults.Parameter["Par_aseah"] = (0.,0.,0)
Defaults.Parameter["Par_abs_u_iso"] = (0.,0.,0)


Defaults.Extract["Extract"] = ""
Defaults.Extract["Extractpos"] = ()
Defaults.Extract["Extractsym"] = []
Defaults.Extract["Rotsup"] = []

Defaults.Experiment["Gen_shift"] = (0.,0.,0)
Defaults.Experiment["Experiment"] = ""
Defaults.Experiment["Emin"] = nan
Defaults.Experiment["Emax"] = nan
Defaults.Experiment["Kev"] = False



# mostly french dictionary
synonyms = dict()
synonyms["dafs"] = "RXS"
synonyms["ecrantage"] = "Screening"
synonyms["gamme"] = "Range"
synonyms["rayon"] = "Radius"
synonyms["seuil"] = "Edge"
synonyms["self_absorption"] = "Self_abs"
synonyms["gamma_fix"] = "Gamma_var"

GroupMembers = dict()
for Group in Groups:
    for Member in getattr(Defaults, Group):
        GroupMembers[Member] = Group



