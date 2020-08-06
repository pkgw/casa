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
    return [Test010_NGC612_ADVANCED,Test020_NGC612_ADVANCED,Test021_NGC612_ADVANCED]

version = '4.7'
class Test010_NGC612_ADVANCED(unittest.TestCase):
    def setUp(self):

        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "datafiles/2012-10-25_0707.C2728",os.getcwd()+'/2012-10-25_0707.C2728')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "datafiles/2012-10-25_0903.C2728",os.getcwd()+'/2012-10-25_0903.C2728')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "datafiles/2012-10-25_1304.C2728",os.getcwd()+'/2012-10-25_1304.C2728')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "datafiles/2012-10-25_1705.C2728",os.getcwd()+'/2012-10-25_1705.C2728')


        if os.path.isfile("/usr/bin/wget"):
            os.system("/usr/bin/wget -O ngc612region.crtf https://casaguides.nrao.edu/images/4/48/Ngc612region.txt")
        elif os.path.isfile("/usr/bin/curl"):
            os.system("/usr/bin/curl https://casaguides.nrao.edu/images/4/48/Ngc612region.txt -o ngc612region.crtf")
        else:
            print("Could not download ngc612region.crtf!")

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/CASA_Guides:ATCA_Advanced_Continuum_Polarization_Tutorial_NGC612-CASA%s'"%(version))
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/CASA_Guides:ATCA_Advanced_Continuum_Polarization_Tutorial_NGC612-CASA%s'"%(version))
    
    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide: NGC612 ADVANCED'''

        exec(compile(open('CASAGuidesATCAAdvancedContinuumPolarizationTutorialNGC612-CASA'+ version +'.py', "rb").read(), 'CASAGuidesATCAAdvancedContinuumPolarizationTutorialNGC612-CASA'+ version +'.py', 'exec'))
                
        return True

class Test020_NGC612_ADVANCED(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass

    def test_6_cal_B0(self):
        '''Test 6: Check cal.B0'''
        tableName = 'cal.B0'
        self.assertTrue(openTable(tableName))

    def test_7_cal_B1(self):
        '''Test 7: Check cal.B1'''
        tableName = 'cal.B1'
        self.assertTrue(openTable(tableName))

    def test_8_cal_D0(self):
        '''Test 8: Check cal.D0'''
        tableName = 'cal.D0'
        self.assertTrue(openTable(tableName))

    def test_9_cal_D1(self):
        '''Test 9: Check cal.D1'''
        tableName = 'cal.D1'
        self.assertTrue(openTable(tableName))

    def test_10_cal_F0(self):
        '''Test 10: Check cal.F0'''
        tableName = 'cal.F0'
        self.assertTrue(openTable(tableName))

    def test_11_cal_G0(self):
        '''Test 11: Check cal.G0'''
        tableName = 'cal.G0'
        self.assertTrue(openTable(tableName))

    def test_12_cal_G1(self):
        '''Test 12: Check cal.G1'''
        tableName = 'cal.G1'
        self.assertTrue(openTable(tableName))

    def test_13_cal_G2(self):
        '''Test 13: Check cal.G2'''
        tableName = 'cal.G2'
        self.assertTrue(openTable(tableName))

    def test_18_ngc612_ms(self):
        '''Test 18: Check ngc612.ms'''
        tableName = 'ngc612.ms'
        self.assertTrue(openTable(tableName))

    def test_19_ngc612_ms_0(self):
        '''Test 19: Check ngc612.ms.0'''
        tableName = 'ngc612.ms.0'
        self.assertTrue(openTable(tableName))

    def test_20_ngc612h_P(self):
        '''Test 20: Check ngc612h.P'''
        tableName = 'ngc612h.P'
        self.assertTrue(openTable(tableName))

    def test_21_ngc612h_Q(self):
        '''Test 21: Check ngc612h.Q'''
        tableName = 'ngc612h.Q'
        self.assertTrue(openTable(tableName))

    def test_22_ngc612h_U(self):
        '''Test 22: Check ngc612h.U'''
        tableName = 'ngc612h.U'
        self.assertTrue(openTable(tableName))

    def test_23_ngc612h_X(self):
        '''Test 23: Check ngc612h.X'''
        tableName = 'ngc612h.X'
        self.assertTrue(openTable(tableName))

    def test_24_ngc612h_flux(self):
        '''Test 24: Check ngc612h.flux'''
        tableName = 'ngc612h.flux'
        self.assertTrue(openTable(tableName))

    def test_25_ngc612h_flux_pbcoverage(self):
        '''Test 25: Check ngc612h.flux.pbcoverage'''
        tableName = 'ngc612h.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_26_ngc612h_image(self):
        '''Test 26: Check ngc612h.image'''
        tableName = 'ngc612h.image'
        self.assertTrue(openTable(tableName))

    def test_27_ngc612h_mask(self):
        '''Test 27: Check ngc612h.mask'''
        tableName = 'ngc612h.mask'
        self.assertTrue(openTable(tableName))

    def test_28_ngc612h_model(self):
        '''Test 28: Check ngc612h.model'''
        tableName = 'ngc612h.model'
        self.assertTrue(openTable(tableName))

    def test_29_ngc612h_psf(self):
        '''Test 29: Check ngc612h.psf'''
        tableName = 'ngc612h.psf'
        self.assertTrue(openTable(tableName))

    def test_30_ngc612h_residual(self):
        '''Test 30: Check ngc612h.residual'''
        tableName = 'ngc612h.residual'
        self.assertTrue(openTable(tableName))

class Test021_NGC612_ADVANCED(unittest.TestCase):
    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):
        os.unlink(os.getcwd()+'/2012-10-25_0707.C2728')
        os.unlink(os.getcwd()+'/2012-10-25_0903.C2728')
        os.unlink(os.getcwd()+'/2012-10-25_1304.C2728')
        os.unlink(os.getcwd()+'/2012-10-25_1705.C2728')
        rmtables("ngc612h*")
        rmtables("cal*")
        rmtables("ngc612*")
        os.system("rm -rf *.last")
        os.system("rm -rf ngc612region.crtf")
        os.system("rm -rf *.flagversions")


    def test_19_Total_polarized_flux_density(self):
        '''Total polarized flux density'''
        mystats1=imstat(imagename='ngc612h.P',mask='ngc612h.P>0.002')
        expected = 1.183
        val = mystats1['flux']
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, 1.181)
        self.assertLessEqual(val, 1.185)
        #assert val >= 1.164 and val <= 1.168, "Total polarized flux density Error Should Equal: 1.183"

    def test_20_Pol_angles_in_western_lobe(self):
        '''First Pol. angle in western lobe'''
        mystats2=imstat(imagename='ngc612h.X',region='box [ [ 325pix , 260pix] , [345pix, 275pix ] ]')
        expected = ((-1.0) * (60.5))
        val = mystats2['mean']
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, ((-1.0)*(62.5)))
        self.assertLessEqual(val, ((-1.0)*(58.5)))
        #assert val >= ((-1.0)*(60.7)) and val <= ((-1.0)*(60.3)) , "Pol. angles in western lobe Error Should Equal: -60.5 Degrees"

    def test_21_Pol_angles_in_western_lobe(self):
        '''Second Pol. angles in western lobe'''
        mystats3=imstat(imagename='ngc612h.X',region='box [ [ 365pix , 250pix] , [385pix,  265pix ] ]')
        expected = 23.8
        val = mystats3['mean']
        print("Expected: %s , Actual: %s"%(expected,val))
        self.assertGreaterEqual(val, 21.8)
        self.assertLessEqual(val, 25.8)
        #assert val >= 23.0 and val <= 23.6 , "Pol. angles in western lobe Error Should Equal: 23.8 Degrees"

if __name__ == '__main__':
    unittest.main()
