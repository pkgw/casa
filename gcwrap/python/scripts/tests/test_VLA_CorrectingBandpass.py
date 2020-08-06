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
    import casac
    from casac import casac
    import numpy as np
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
    return [Test010_CorrectingBandpassCalib,Test020_CorrectingBandpassCalib,Test021_CorrectingBandpassCalib]

class Test010_CorrectingBandpassCalib(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        # Untar data. Dataset needs write permissions so copy a full set to working dir
        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/G192-BP.ms.tar.gz")
        tar.extractall()
        tar.close()


        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Bandpass_Slope'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Bandpass_Slope'")

        time.sleep(5) # Allow extract time to download script
        lines = open('VLACASABandpassSlope.py')
        file = open("newfile.txt", "w")
        for line in lines:
            if line.startswith("pl."):
                continue
            if "%cpaste" in line: 
                continue
            file.write(line)
        file.close()
        os.remove('VLACASABandpassSlope.py')
        os.rename("newfile.txt",'VLACASABandpassSlope.py')

    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide:  Topical Guide Correcting Bandpass VLA Data'''

        exec(compile(open('VLACASABandpassSlope.py', "rb").read(), 'VLACASABandpassSlope.py', 'exec'))
                
        return True

class Test020_CorrectingBandpassCalib(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass


    def test_1_3C84_fluxinfo(self):
        '''Test 1: Check 3C84.fluxinfo'''
        self.assertTrue(assert_file('3C84.fluxinfo'))

    def test_2_G192_BP_ms(self):
        '''Test 2: Check G192-BP.ms'''
        tableName = 'G192-BP.ms'
        self.assertTrue(openTable(tableName))

    def test_3_calG192_B0(self):
        '''Test 3: Check calG192.B0'''
        tableName = 'calG192.B0'
        self.assertTrue(openTable(tableName))

    def test_4_calG192_B0_b(self):
        '''Test 4: Check calG192.B0.b'''
        tableName = 'calG192.B0.b'
        self.assertTrue(openTable(tableName))

    def test_5_calG192_F1(self):
        '''Test 5: Check calG192.F1'''
        tableName = 'calG192.F1'
        self.assertTrue(openTable(tableName))

    def test_6_calG192_G0(self):
        '''Test 6: Check calG192.G0'''
        tableName = 'calG192.G0'
        self.assertTrue(openTable(tableName))

    def test_7_calG192_G0_b(self):
        '''Test 7: Check calG192.G0.b'''
        tableName = 'calG192.G0.b'
        self.assertTrue(openTable(tableName))

    def test_8_calG192_G1(self):
        '''Test 8: Check calG192.G1'''
        tableName = 'calG192.G1'
        self.assertTrue(openTable(tableName))

    def test_9_calG192_G1p(self):
        '''Test 9: Check calG192.G1p'''
        tableName = 'calG192.G1p'
        self.assertTrue(openTable(tableName))

    def test_10_calG192_K0(self):
        '''Test 10: Check calG192.K0'''
        tableName = 'calG192.K0'
        self.assertTrue(openTable(tableName))

    def test_11_calG192_K0_b(self):
        '''Test 11: Check calG192.K0.b'''
        tableName = 'calG192.K0.b'
        self.assertTrue(openTable(tableName))

    """
    def test_12_plotG192_plotcal_G0p1_png(self):
        '''Test 12: Check plotG192_plotcal_G0p1.png'''
        self.assertTrue(assert_file('plotG192_plotcal_G0p1.png'))

    def test_14_plotG192_plotcal_G0p2_png(self):
        '''Test 14: Check plotG192_plotcal_G0p2.png'''
        self.assertTrue(assert_file('plotG192_plotcal_G0p2.png'))
    """
####################################################################################################

class Test021_CorrectingBandpassCalib(unittest.TestCase):

    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):

        rmtables("calG192*")
        os.system("rm -rf 3C84.fluxinfo")
        rmtables("G192-BP.ms*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.flagversions")


    def test_16_Flux_density_3c84_J0319_413(self):
        '''Flux density for 3c84-J0319+413'''
        flux1 = fluxscale(vis='G192-BP.ms', caltable='calG192.G1',  fluxtable='calG192.F1', reference='0', transfer='1', listfile='3C84.fluxinfo', fitorder=1)
        fluxDensity = float(flux1['1']['fitFluxd'])
        expected = 29.0285 
        print("Expected: %s , Actual: %s"%(expected,fluxDensity))
        assert 28.9976336 <= fluxDensity <= 29.0593664, "Error in Flux density for 3c84-J0319+413 with 28.9976336 <= %s <= 29.0593664"%(fluxDensity)

    def test_17_Spidx_3c84_J0319_413(self):
        '''Spidx for 3c84-J0319+413'''
        flux1 = fluxscale(vis='G192-BP.ms', caltable='calG192.G1',  fluxtable='calG192.F1', reference='0', transfer='1', listfile='3C84.fluxinfo', fitorder=1)
        Spidx = float(flux1['1']['spidx'][1])
        expected = (-1.0) * 0.538791
        print("Expected: %s , Actual: %s"%(expected,Spidx))
        assert  ((-1.0)* 0.54762051) <= Spidx <= ((-1.0)* 0.52996149), "Error in Flux Spidx for 3c84-J0319+413 -0.54762051 <= %s <= -0.52996149" %(Spidx)

class Test022_CorrectingBandpassCalib(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass

    def test_15_compare_plotG192_plotcal_G0p2_png(self):
        '''Test 15: Image Comparison plotG192_plotcal_G0p2.png'''
        imageName = 'plotG192_plotcal_G0p2.png'
        self.assertTrue(compareImages('plotG192_plotcal_G0p2.png'))

####################################################################################################

