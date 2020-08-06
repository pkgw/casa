import os
import shutil
import casac
from tasks import *
from taskinit import *
from __main__ import *
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

import tarfile

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
    return [Test010_FirstLookatImaging,Test020_FirstLookatImaging]

class Test010_FirstLookatImaging(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/twhya_smoothed.ms'):
            shutil.rmtree(os.getcwd()+'/twhya_smoothed.ms')  

        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"
    
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_calibrated.ms",os.getcwd()+'/sis14_twhya_calibrated.ms')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_uncalibrated.ms",os.getcwd()+'/sis14_twhya_uncalibrated.ms')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_calibrated_flagged.ms",os.getcwd()+'/sis14_twhya_calibrated_flagged.ms')

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Imaging'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Imaging'")
        
        file = open("my_script.py", "w");file.write("print 'Execfile executed'");file.close()
        time.sleep(5) # Allow extract time to download script

        lines = open('FirstLookatImaging.py')
        file = open("newfile.txt", "w")
        for line in lines:

            if "rm -rf sis14_twhya_calibrated.ms" in line:
                continue

            if "rm -rf sis14_twhya_uncalibrated.ms" in line:
                continue

            if "rm -rf sis14_twhya_calibrated_flagged.ms" in line:
                continue

            pattern = r'''niter\ *=\ *(5000)'''
            if re.search(pattern,line):
                line = re.sub( pattern, 'niter=250', line )

            file.write(line)
        file.close()
        os.remove('FirstLookatImaging.py')
        os.rename("newfile.txt",'FirstLookatImaging.py')

        time.sleep(5)
    def tearDown(self):

        pass
    
    def test_00_runGuide(self):
        '''Run Casa Guide: First Look at Imaging'''


        exec(compile(open('FirstLookatImaging.py', "rb").read(), 'FirstLookatImaging.py', 'exec'))
                
        return True

class Test020_FirstLookatImaging(unittest.TestCase):

    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):
        os.unlink(os.getcwd()+'/sis14_twhya_calibrated.ms')
        os.unlink(os.getcwd()+'/sis14_twhya_uncalibrated.ms')
        os.unlink(os.getcwd()+'/sis14_twhya_calibrated_flagged.ms')

        rmtables("amp_cal*")
        rmtables("phase_cal*")
        rmtables("twhya*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.flagversions")
        os.system("rm -rf my_script.py")

    def test_1_amp_cal_bigpix_image(self):
        tableName = 'amp_cal_bigpix.image'
        self.assertTrue(openTable(tableName))

    def test_2_amp_cal_bigpix_mask(self):
        tableName = 'amp_cal_bigpix.mask'
        self.assertTrue(openTable(tableName))

    def test_3_amp_cal_bigpix_model(self):
        tableName = 'amp_cal_bigpix.model'
        self.assertTrue(openTable(tableName))

    def test_4_amp_cal_bigpix_bp(self):
        tableName = 'amp_cal_bigpix.pb'
        self.assertTrue(openTable(tableName))

    def test_5_amp_cal_bigpix_psf(self):
        tableName = 'amp_cal_bigpix.psf'
        self.assertTrue(openTable(tableName))

    def test_6_amp_cal_bigpix_residual(self):
        tableName = 'amp_cal_bigpix.residual'
        self.assertTrue(openTable(tableName))

    def test_7_amp_cal_bigpix_sumwt(self):
        tableName = 'amp_cal_bigpix.sumwt'
        self.assertTrue(openTable(tableName))

    def test_8_amp_cal_robust_image(self):
        tableName = 'amp_cal_robust.image'
        self.assertTrue(openTable(tableName))

    def test_9_amp_cal_robust_mask(self):
        tableName = 'amp_cal_robust.mask'
        self.assertTrue(openTable(tableName))

    def test_10_amp_cal_robust_model(self):
        tableName = 'amp_cal_robust.model'
        self.assertTrue(openTable(tableName))

    def test_11_amp_cal_robust_bp(self):
        tableName = 'amp_cal_robust.pb'
        self.assertTrue(openTable(tableName))

    def test_12_amp_cal_robust_psf(self):
        tableName = 'amp_cal_robust.psf'
        self.assertTrue(openTable(tableName))

    def test_13_amp_cal_robust_residual(self):
        tableName = 'amp_cal_robust.residual'
        self.assertTrue(openTable(tableName))

    def test_14_amp_cal_robust_sumwt(self):
        tableName = 'amp_cal_robust.sumwt'
        self.assertTrue(openTable(tableName))

    def test_15_phase_cal_image(self):
        tableName = 'phase_cal.image'
        self.assertTrue(openTable(tableName))

    def test_16_phase_cal_mask(self):
        tableName = 'phase_cal.mask'
        self.assertTrue(openTable(tableName))

    def test_17_phase_cal_model(self):
        tableName = 'phase_cal.model'
        self.assertTrue(openTable(tableName))

    def test_18_phase_cal_bp(self):
        tableName = 'phase_cal.pb'
        self.assertTrue(openTable(tableName))

    def test_19_phase_cal_psf(self):
        tableName = 'phase_cal.psf'
        self.assertTrue(openTable(tableName))

    def test_20_phase_cal_residual(self):
        tableName = 'phase_cal.residual'
        self.assertTrue(openTable(tableName))

    def test_21_phase_cal_sumwt(self):
        tableName = 'phase_cal.sumwt'
        self.assertTrue(openTable(tableName))

    def test_22_phase_cal_robust_image(self):
        tableName = 'phase_cal_robust.image'
        self.assertTrue(openTable(tableName))

    def test_23_phase_cal_robust_mask(self):
        tableName = 'phase_cal_robust.mask'
        self.assertTrue(openTable(tableName))

    def test_24_phase_cal_robust_model(self):
        tableName = 'phase_cal_robust.model'
        self.assertTrue(openTable(tableName))

    def test_25_phase_cal_robust_bp(self):
        tableName = 'phase_cal_robust.pb'
        self.assertTrue(openTable(tableName))

    def test_26_phase_cal_robust_psf(self):
        tableName = 'phase_cal_robust.psf'
        self.assertTrue(openTable(tableName))

    def test_27_phase_cal_robust_residual(self):
        tableName = 'phase_cal_robust.residual'
        self.assertTrue(openTable(tableName))

    def test_28_phase_cal_robust_sumwt(self):
        tableName = 'phase_cal_robust.sumwt'
        self.assertTrue(openTable(tableName))

    def test_29_phase_cal_uncalibrated_image(self):
        tableName = 'phase_cal_uncalibrated.image'
        self.assertTrue(openTable(tableName))

    def test_30_phase_cal_uncalibrated_mask(self):
        tableName = 'phase_cal_uncalibrated.mask'
        self.assertTrue(openTable(tableName))

    def test_31_phase_cal_uncalibrated_model(self):
        tableName = 'phase_cal_uncalibrated.model'
        self.assertTrue(openTable(tableName))

    def test_32_phase_cal_uncalibrated_bp(self):
        tableName = 'phase_cal_uncalibrated.pb'
        self.assertTrue(openTable(tableName))

    def test_33_phase_cal_uncalibrated_psf(self):
        tableName = 'phase_cal_uncalibrated.psf'
        self.assertTrue(openTable(tableName))

    def test_34_phase_cal_uncalibrated_residual(self):
        tableName = 'phase_cal_uncalibrated.residual'
        self.assertTrue(openTable(tableName))

    def test_35_phase_cal_uncalibrated_sumwt(self):
        tableName = 'phase_cal_uncalibrated.sumwt'
        self.assertTrue(openTable(tableName))

    def test_36_phase_cal_unflagged_image(self):
        tableName = 'phase_cal_unflagged.image'
        self.assertTrue(openTable(tableName))

    def test_37_phase_cal_unflagged_mask(self):
        tableName = 'phase_cal_unflagged.mask'
        self.assertTrue(openTable(tableName))

    def test_38_phase_cal_unflagged_model(self):
        tableName = 'phase_cal_unflagged.model'
        self.assertTrue(openTable(tableName))

    def test_39_phase_cal_unflagged_bp(self):
        tableName = 'phase_cal_unflagged.pb'
        self.assertTrue(openTable(tableName))

    def test_40_phase_cal_unflagged_psf(self):
        tableName = 'phase_cal_unflagged.psf'
        self.assertTrue(openTable(tableName))

    def test_41_phase_cal_unflagged_residual(self):
        tableName = 'phase_cal_unflagged.residual'
        self.assertTrue(openTable(tableName))

    def test_42_phase_cal_unflagged_sumwt(self):
        tableName = 'phase_cal_unflagged.sumwt'
        self.assertTrue(openTable(tableName))

    def test_43_twhya_cont_auto_image(self):
        tableName = 'twhya_cont_auto.image'
        self.assertTrue(openTable(tableName))

    def test_44_twhya_cont_auto_mask(self):
        tableName = 'twhya_cont_auto.mask'
        self.assertTrue(openTable(tableName))

    def test_45_twhya_cont_auto_model(self):
        tableName = 'twhya_cont_auto.model'
        self.assertTrue(openTable(tableName))

    def test_46_twhya_cont_auto_bp(self):
        tableName = 'twhya_cont_auto.pb'
        self.assertTrue(openTable(tableName))

    def test_47_twhya_cont_auto_psf(self):
        tableName = 'twhya_cont_auto.psf'
        self.assertTrue(openTable(tableName))

    def test_48_twhya_cont_auto_residual(self):
        tableName = 'twhya_cont_auto.residual'
        self.assertTrue(openTable(tableName))

    def test_49_twhya_cont_auto_sumwt(self):
        tableName = 'twhya_cont_auto.sumwt'
        self.assertTrue(openTable(tableName))

    def test_50_twhya_cont_image(self):
        tableName = 'twhya_cont.image'
        self.assertTrue(openTable(tableName))

    def test_51_twhya_cont_mask(self):
        tableName = 'twhya_cont.mask'
        self.assertTrue(openTable(tableName))

    def test_52_twhya_cont_model(self):
        tableName = 'twhya_cont.model'
        self.assertTrue(openTable(tableName))

    def test_53_twhya_cont_bp(self):
        tableName = 'twhya_cont.pb'
        self.assertTrue(openTable(tableName))

    def test_54_twhya_cont_psf(self):
        tableName = 'twhya_cont.psf'
        self.assertTrue(openTable(tableName))

    def test_55_twhya_cont_residual(self):
        tableName = 'twhya_cont.residual'
        self.assertTrue(openTable(tableName))

    def test_56_twhya_cont_sumwt(self):
        tableName = 'twhya_cont.sumwt'
        self.assertTrue(openTable(tableName))

    def test_57_twhya_smoothed_ms(self):
        '''Test 51: Check twhya_smoothed.ms'''
        tableName = 'twhya_smoothed.ms'
        self.assertTrue(openTable(tableName))


if __name__ == '__main__':
    unittest.main()
