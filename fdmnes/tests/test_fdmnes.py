import os
import time
import unittest
from testfixtures import TempDirectory
import urllib
from ..pyFDMNES import fdmnes

class test_fdmnes(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
        self.url = "http://rruff.geo.arizona.edu/AMS/CIF_text_files/13750_cif.txt"
        self.element = "Ti"
        self.edge = "K"
        self.fermilevel = 4.9664
        self.MM = 47.867
        self.ciffile = os.path.join(self.dir.path,"rutile.cif")
        self.inputfile = os.path.join(self.dir.path,"rutile.txt")

    def tearDown(self):
        self.dir.cleanup()

    def test_simulate(self):
        
        urllib.urlretrieve(self.url, self.ciffile)
        sim = fdmnes(self.ciffile, resonant=self.element)

        # Settings
        sim.P.Energpho = False # relative energies as output
        sim.P.Radius = 3.5 # Radius of the cluster for calculation
        sim.P.Rpotmax = sim.P.Radius + 5 # Radius of the cluster for potential calculation
        sim.P.Quadrupole = True # multipole approximation
        sim.P.Green = False # MS instead of FDM (faster but less accurate)
        sim.P.TDDFT = True # multi electron correction
        sim.P.Convolution = True # 
        sim.P.Density = False # save density of states
        sim.P.Edge = self.edge
        sim.P.Range = (-10,2,100) # eV

        sim.WriteInputFile(self.inputfile, overwrite=True)

        sim.Run(wait=False)
        while True:
            NumRunning = sum(sim.Status(i)==False for i in range(len(sim.proc)))
            if NumRunning == 0:
                break
            time.sleep(5)

        data = sim.get_XANES(conv = True)

        energy = data[:,0]/1000. + self.fermilevel # ev -> keV
        mu = data[:,1]*6.022140857e5/self.MM # Absorption cross section (Mbarn/atom) -> mass absorption coefficient (cm^2/g)

        #import matplotlib.pyplot as plt
        #plt.plot(energy,mu)
        #plt.show()

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fdmnes("test_simulate"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

