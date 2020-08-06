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

indir = os.environ.get('CASAPATH').split()[0] + '/../../'
"""
def compareTables(dataset):
    refmeaset = indir + "pybotWorkspace/prime_FirstLookatSelfCalibration/" + str(dataset)
    measet = indir + "pybotWorkspace/FirstLookatSelfCalibration/" + str(dataset)
    tb = casac.table()
    tb.open(str(dataset))
    cnames = tb.colnames()
    tb.close()
    cnames = ['FLAG', 'WEIGHT', 'SIGMA','DATA']
    #cnames = ['WEIGHT_SPECTRUM']
    for i in range(0,len(cnames)):
        boolean = compVarColTables(refmeaset, measet, str(cnames[i]), tolerance=0.001)
        clearstat()
        if boolean == False:
            return False
        gc.collect()
    return True
"""


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
    return [Test010_FirstLookatSelfCalibration,Test020_FirstLookatSelfCalibration]
    #return [Test010_FirstLookatSelfCalibration,Test020_FirstLookatSelfCalibration,Test021_FirstLookatSelfCalibration] # Test021_FirstLookatSelfCalibration Requires a method to compare images

class Test010_FirstLookatSelfCalibration(unittest.TestCase):
    def setUp(self):

        dirsToRemove = ['sis14_twhya_selfcal.ms','sis14_twhya_selfcal_2.ms','sis14_twhya_selfcal_3.ms']
        for dataset in dirsToRemove:
            if os.path.isdir(os.getcwd()+'/%s'%(dataset)):
                try: os.unlink(os.getcwd()+'/%s'%(dataset))
                except: shutil.rmtree(os.getcwd()+'/%s'%(dataset))

        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/ss_alma_data_v1p2.tar.gz")
        tar.extractall()
        tar.close()

        shutil.move(os.getcwd()+"/working_data/sis14_twhya_selfcal.ms",os.getcwd()+'/sis14_twhya_selfcal.ms')
        shutil.move(os.getcwd()+"/working_data/sis14_twhya_calibrated_flagged.ms",os.getcwd()+'/sis14_twhya_calibrated_flagged.ms')

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Self_Calibration'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Self_Calibration'")
    
        time.sleep(5) # Allow extract time to download script

        lines = open('FirstLookatSelfCalibration.py')
        file = open("newfile.txt", "w")
        for line in lines:

            if "rm -rf sis14_twhya_calibrated_flagged.ms" in line:
                continue

            pattern = r'''niter\ *=\ *(5000)'''
            if re.search(pattern,line):
                line = re.sub( pattern, 'niter=250', line )

            file.write(line)
        file.close()
        os.remove('FirstLookatSelfCalibration.py')
        os.rename("newfile.txt",'FirstLookatSelfCalibration.py')

        time.sleep(5)

    def tearDown(self):
        pass

    def test_00_runGuide(self):
        '''Run Casa Guide: First Look at Self-Calibration'''

        exec(compile(open('FirstLookatSelfCalibration.py', "rb").read(), 'FirstLookatSelfCalibration.py', 'exec'))
                
        return True

class Test020_FirstLookatSelfCalibration(unittest.TestCase):

    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):
        rmtables("first_image*")
        rmtables("second_image*")
        rmtables("third_image*")
        rmtables("fourth_image*")

        rmtables("sis14_twhya*")
        rmtables("*.cal")
        os.system("rm -rf *.last")

        os.system("rm -rf lessons")
        os.system("rm -rf working_data")
        os.system("rm -rf *.flagversions")

    def test_1_amp_cal(self):
        '''Test 1: Check amp.cal'''
        tableName = 'amp.cal'
        self.assertTrue(openTable(tableName))

    def test_2_first_image_pb(self):
        '''Test 2: Check first_image.pb'''
        tableName = 'first_image.pb'
        self.assertTrue(openTable(tableName))

    def test_3_first_image_image(self):
        '''Test 3: Check first_image.image'''
        tableName = 'first_image.image'
        self.assertTrue(openTable(tableName))

    def test_4_first_image_model(self):
        '''Test 4: Check first_image.model'''
        tableName = 'first_image.model'
        self.assertTrue(openTable(tableName))

    def test_5_first_image_psf(self):
        '''Test 5: Check first_image.psf'''
        tableName = 'first_image.psf'
        self.assertTrue(openTable(tableName))

    def test_6_first_image_residual(self):
        '''Test 6: Check first_image.residual'''
        tableName = 'first_image.residual'
        self.assertTrue(openTable(tableName))

    def test_7_first_image_mask(self):
        '''Test 7: Check first_image.mask'''
        tableName = 'first_image.mask'
        self.assertTrue(openTable(tableName))

    def test_8_first_image_sumwt(self):
        '''Test 8: Check first_image.sumwt'''
        tableName = 'first_image.sumwt'
        self.assertTrue(openTable(tableName))

    def test_9_fourth_image_pb(self):
        '''Test 9: Check fourth_image.pb'''
        tableName = 'fourth_image.pb'
        self.assertTrue(openTable(tableName))

    def test_10_fourth_image_image(self):
        '''Test 10: Check fourth_image.image'''
        tableName = 'fourth_image.image'
        self.assertTrue(openTable(tableName))

    def test_11_fourth_image_model(self):
        '''Test 11: Check fourth_image.model'''
        tableName = 'fourth_image.model'
        self.assertTrue(openTable(tableName))

    def test_12_fourth_image_psf(self):
        '''Test 12: Check fourth_image.psf'''
        tableName = 'fourth_image.psf'
        self.assertTrue(openTable(tableName))

    def test_13_fourth_image_residual(self):
        '''Test 13: Check fourth_image.residual'''
        tableName = 'fourth_image.residual'
        self.assertTrue(openTable(tableName))

    def test_14_fourth_image_mask(self):
        '''Test 14: Check fourth_image.mask'''
        tableName = 'fourth_image.mask'
        self.assertTrue(openTable(tableName))

    def test_15_fourth_image_sumwt(self):
        '''Test 8: Check fourth_image.sumwt'''
        tableName = 'fourth_image.sumwt'
        self.assertTrue(openTable(tableName))

    def test_16_phase_cal(self):
        '''Test 12: Check phase.cal'''
        tableName = 'phase.cal'
        self.assertTrue(openTable(tableName))

    def test_17_phase_2_cal(self):
        '''Test 13: Check phase_2.cal'''
        tableName = 'phase_2.cal'
        self.assertTrue(openTable(tableName))

    def test_18_second_image_pb(self):
        '''Test 18: Check second_image.pb'''
        tableName = 'second_image.pb'
        self.assertTrue(openTable(tableName))

    def test_19_second_image_image(self):
        '''Test 19: Check second_image.image'''
        tableName = 'second_image.image'
        self.assertTrue(openTable(tableName))

    def test_20_second_image_model(self):
        '''Test 20: Check second_image.model'''
        tableName = 'second_image.model'
        self.assertTrue(openTable(tableName))

    def test_21_second_image_psf(self):
        '''Test 21: Check second_image.psf'''
        tableName = 'second_image.psf'
        self.assertTrue(openTable(tableName))

    def test_22_second_image_residual(self):
        '''Test 22: Check second_image.residual'''
        tableName = 'second_image.residual'
        self.assertTrue(openTable(tableName))

    def test_23_second_image_mask(self):
        '''Test 23: Check second_image.mask'''
        tableName = 'second_image.mask'
        self.assertTrue(openTable(tableName))

    def test_24_second_image_sumwt(self):
        '''Test 24: Check second_image.sumwt'''
        tableName = 'second_image.sumwt'
        self.assertTrue(openTable(tableName))

    def test_25_sis14_twhya_selfcal_ms(self):
        '''Test 25: Check sis14_twhya_selfcal.ms'''
        tableName = 'sis14_twhya_selfcal.ms'
        self.assertTrue(openTable(tableName))

    def test_26_sis14_twhya_selfcal_2_ms(self):
        '''Test 26: Check sis14_twhya_selfcal_2.ms'''
        tableName = 'sis14_twhya_selfcal_2.ms'
        self.assertTrue(openTable(tableName))

    def test_27_sis14_twhya_selfcal_3_ms(self):
        '''Test 27: Check sis14_twhya_selfcal_3.ms'''
        tableName = 'sis14_twhya_selfcal_3.ms'
        self.assertTrue(openTable(tableName))

    def test_28_third_image_pb(self):
        '''Test 28: Check third_image.pb'''
        tableName = 'third_image.pb'
        self.assertTrue(openTable(tableName))

    def test_29_third_image_image(self):
        '''Test 29: Check third_image.image'''
        tableName = 'third_image.image'
        self.assertTrue(openTable(tableName))

    def test_30_third_image_model(self):
        '''Test 30: Check third_image.model'''
        tableName = 'third_image.model'
        self.assertTrue(openTable(tableName))

    def test_31_third_image_psf(self):
        '''Test 31: Check third_image.psf'''
        tableName = 'third_image.psf'
        self.assertTrue(openTable(tableName))

    def test_32_third_image_residual(self):
        '''Test 32: Check third_image.residual'''
        tableName = 'third_image.residual'
        self.assertTrue(openTable(tableName))

    def test_33_second_image_mask(self):
        '''Test 33: Check third_image.mask'''
        tableName = 'third_image.mask'
        self.assertTrue(openTable(tableName))

    def test_34_third_image_sumwt(self):
        '''Test 34: Check third_image.sumwt'''
        tableName = 'third_image.sumwt'
        self.assertTrue(openTable(tableName))

    # Turn these tests on if Plotting is on.
    #def test_35_sis14_selfcal_phase_scan_png(self):
    #    '''Test 35: Check sis14_selfcal_phase_scan.png'''
    #    self.assertTrue(assert_file('sis14_selfcal_phase_scan.png'))

    #def test_37_sis14_selfcal_phase_scan_2_png(self):
    #    '''Test 37: Check sis14_selfcal_phase_scan_2.png'''
    #    self.assertTrue(assert_file('sis14_selfcal_phase_scan_2.png'))

####################################################################################################
# Needs to be Fixed
class Test021_FirstLookatSelfCalibration(unittest.TestCase):
    def test_36_compare_sis14_selfcal_phase_scan_png(self):
        '''Test 36: Image Comparison sis14_selfcal_phase_scan.png'''
        imageName = 'sis14_selfcal_phase_scan.png'
        self.assertTrue(compareImages('sis14_selfcal_phase_scan.png'))


    def test_38_compare_sis14_selfcal_phase_scan_2_png(self):
        '''Test 38: Image Comparison sis14_selfcal_phase_scan_2.png'''
        imageName = 'sis14_selfcal_phase_scan_2.png'
        self.assertTrue(compareImages('sis14_selfcal_phase_scan_2.png'))

if __name__ == '__main__':
    unittest.main()
