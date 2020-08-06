import os
import shutil
import casac
from tasks import *
from taskinit import *
from __main__ import *
import unittest
import subprocess

import time
import string
import importlib
import inspect
import shutil
import logging
import hashlib
import threading
import imp
from contextlib import contextmanager
from contextlib import closing
from functools import wraps
import unittest
from numpy import count_nonzero



def compVarColTables(referencetab, testtab, varcol, tolerance=0.):
    '''Compare a variable column of two tables.
       referencetab  --> a reference table
       testtab       --> a table to verify
       varcol        --> the name of a variable column (str)
       Returns True or False.
    '''
    
    retval = True
    tb2 = casac.table()

    tb.open(referencetab)
    cnames = tb.colnames()

    tb2.open(testtab)
    col = varcol
    if tb.isvarcol(col) and tb2.isvarcol(col):
        try:
            # First check
            if tb.nrows() != tb2.nrows():
                print('Length of '+ str(referencetab) +' differ from '+ str(testtab)+','+ str(tb.nrows())+ '!=' + str(tb2.nrows()))
                retval = False
            else:
                for therow in range(tb.nrows()):
            
                    rdata = tb.getcell(col,therow)
                    tdata = tb2.getcell(col,therow)
                    if not rdata.all()==tdata.all():
                        if (tolerance>0.):
                            differs=False
                            for j in range(0,len(rdata)):

                                if ((isinstance(rdata[j],float)) or (isinstance(rdata[j],int))):
                                    if (abs(rdata[j]-tdata[j]) > tolerance*abs(rdata[j]+tdata[j])):
                                        differs = True
                                elif (isinstance(rdata[j],list)) or (isinstance(rdata[j],np.ndarray)):
                                    for k in range(0,len(rdata[j])):
                                        if (abs(rdata[j][k]-tdata[j][k]) > tolerance*abs(rdata[j][k]+tdata[j][k])):
                                            differs = True
                                if differs:
                                    print('ERROR: Column ' + str(col) + ' of '  + str(referencetab) +  ' and ' + str(testtab)+  ' do not agree within tolerance '+ str(tolerance))
                                    break
                        else:
                            print('ERROR: Column ' +str(col)+ ' of ' +str(referencetab)+ ' and ' +str(testtab) + ' do not agree.')
                            print('ERROR: First row to differ is row=' + str(therow))
                            retval = False
                            break
        finally:
            tb.close()
            tb2.close()
    
    else:
        print('Column: ' +str(col) + 'are not varcolumns.')
        retval = False

    if retval:
        print('Column ' + str(col) + ' of '  + str(referencetab) +  ' and ' + str(testtab) + ' agree')
        
    return retval
 
def assert_file(file):
    return os.access(file, os.F_OK)

def openTable(tableName):
    try:
        import casac
        from casac import casac
        tb = casac.table()
        tb.open(str(tableName))
        tb.close()
        return True
    except:
        return False

def msHandler(file):
    table_instance = tbtool()
    table_instance.open(file)
    return table_instance

def suite():
    return [Test010_NGC612,Test020_NGC612,Test021_NGC612]


version = '4.7'
class Test010_NGC612(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "ngc612.uv",os.getcwd()+'/ngc612.uv')

        if os.path.isfile("/usr/bin/wget"):
            os.system("/usr/bin/wget -O ngc612region.crtf https://casaguides.nrao.edu/images/4/48/Ngc612region.txt")
        elif os.path.isfile("/usr/bin/curl"):
            os.system("/usr/bin/curl https://casaguides.nrao.edu/images/4/48/Ngc612region.txt -o ngc612region.crtf")
        else:
            print("Could not download ngc612region.crtf!")

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/ATCA_Continuum_Polarization_Tutorial_NGC612-CASA%s'"%(version))
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/ATCA_Continuum_Polarization_Tutorial_NGC612-CASA%s'"%(version))
    
    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide: NGC612'''

        exec(compile(open('ATCAContinuumPolarizationTutorialNGC612-CASA' + version +'.py', "rb").read(), 'ATCAContinuumPolarizationTutorialNGC612-CASA' + version +'.py', 'exec'))
                
        return True

class Test020_NGC612(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass

    def test_8_ngc612_ms(self):
        '''Test 8: Check ngc612.ms'''
        tableName = 'ngc612.ms'
        self.assertTrue(openTable(tableName))

    def test_9_ngc612h_P(self):
        '''Test 9: Check ngc612h.P'''
        tableName = 'ngc612h.P'
        self.assertTrue(openTable(tableName))

    def test_10_ngc612h_X(self):
        '''Test 10: Check ngc612h.X'''
        tableName = 'ngc612h.X'
        self.assertTrue(openTable(tableName))

    def test_11_ngc612h_flux(self):
        '''Test 11: Check ngc612h.flux'''
        tableName = 'ngc612h.flux'
        self.assertTrue(openTable(tableName))

    def test_12_ngc612h_flux_pbcoverage(self):
        '''Test 12: Check ngc612h.flux.pbcoverage'''
        tableName = 'ngc612h.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_13_ngc612h_image(self):
        '''Test 13: Check ngc612h.image'''
        tableName = 'ngc612h.image'
        self.assertTrue(openTable(tableName))

    def test_14_ngc612h_mask(self):
        '''Test 14: Check ngc612h.mask'''
        tableName = 'ngc612h.mask'
        self.assertTrue(openTable(tableName))

    def test_15_ngc612h_model(self):
        '''Test 15: Check ngc612h.model'''
        tableName = 'ngc612h.model'
        self.assertTrue(openTable(tableName))

    def test_16_ngc612h_psf(self):
        '''Test 16: Check ngc612h.psf'''
        tableName = 'ngc612h.psf'
        self.assertTrue(openTable(tableName))

    def test_17_ngc612h_residual(self):
        '''Test 17: Check ngc612h.residual'''
        tableName = 'ngc612h.residual'
        self.assertTrue(openTable(tableName))
    """
    def test_18_ngc612h_flagversions(self):
        '''Test 17: Check ngc612.ms.flagversions'''
        tableName = 'ngc612.ms.flagversions'
        self.assertTrue(openTable(tableName))
    """

class Test021_NGC612(unittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        os.unlink(os.getcwd()+'/ngc612.uv')
        rmtables("ngc612h*")
        rmtables("ngc612.ms")
        os.system("rm -rf *.last")
        os.system("rm -rf ngc612region.crtf")
        os.system("rm -rf *.flagversions")

    def test_19_Total_polarized_flux_density(self):
        '''Total polarized flux density'''
        mystats1=imstat(imagename='ngc612h.P',mask='ngc612h.P>0.002')
        expected = 1.166
        val = float(mystats1['flux'])
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, 1.164)
        self.assertLessEqual(val, 1.168)
        #assert val >= 1.164 and val <= 1.168, "Total polarized flux density Error Should Equal: 1.166"

    def test_20_Pol_angles_in_western_lobe(self):
        '''First Pol. angle in western lobe'''
        mystats2=imstat(imagename='ngc612h.X',region='box [ [ 325pix , 260pix] , [345pix, 275pix ] ]')
        expected = ((-1.0) * (61.7))
        val = mystats2['mean']
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, ((-1.0)*(63.7)))
        self.assertLessEqual(val, ((-1.0)*(59.7)))
        #assert val >= ((-1.0)*(61.9)) and val <= ((-1.0)*(61.5)) , "Pol. angles in western lobe Error Should Equal: -61.7 Degrees"

    def test_21_Pol_angles_in_western_lobe(self):
        '''Second Pol. angle in western lobe'''
        mystats3=imstat(imagename='ngc612h.X',region='box [ [ 365pix , 250pix] , [385pix,  265pix ] ]')
        expected = 22.8
        val = mystats3['mean'][0]
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, 21.2)
        self.assertLessEqual(val, 25.2)
        #assert val >= 23.0 and val <= 23.4 , "Pol. angles in western lobe Error Should Equal: 23.2 Degrees"

if __name__ == '__main__':
    unittest.main()
