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
    return [Test010_ImagingVLAData,Test020_ImagingVLAData]


class Test010_ImagingVLAData(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        # Untar data. Dataset needs write permissions so copy a full set to working dir
        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/SNR_G55_10s.calib.tar.gz")
        tar.extractall()
        tar.close()
        shutil.copyfile(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "outliers.txt",os.getcwd()+'/outliers.txt')


        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Imaging'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Imaging'")

        time.sleep(5) # Allow extract time to download script


    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide:  Topical Guide Imaging VLA Data'''

        exec(compile(open('VLACASAImaging.py', "rb").read(), 'VLACASAImaging.py', 'exec'))
                
        return True

class Test020_ImagingVLAData(unittest.TestCase):

    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):
        rmtables("Outlier1.MS*")
        rmtables("Outlier2.MS*")
        rmtables("SNR.MS.MTMFS*")
        rmtables("SNR_G55_10s*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.flagversions")


    def test_1_Outlier1_MS_MTMFS_flux(self):
        '''Test 1: Check Outlier1.MS.MTMFS.flux'''
        tableName = 'Outlier1.MS.MTMFS.flux'
        self.assertTrue(openTable(tableName))

    def test_2_Outlier1_MS_MTMFS_image(self):
        '''Test 2: Check Outlier1.MS.MTMFS.image'''
        tableName = 'Outlier1.MS.MTMFS.image'
        self.assertTrue(openTable(tableName))

    def test_3_Outlier1_MS_MTMFS_model(self):
        '''Test 3: Check Outlier1.MS.MTMFS.model'''
        tableName = 'Outlier1.MS.MTMFS.model'
        self.assertTrue(openTable(tableName))

    def test_4_Outlier1_MS_MTMFS_psf(self):
        '''Test 4: Check Outlier1.MS.MTMFS.psf'''
        tableName = 'Outlier1.MS.MTMFS.psf'
        self.assertTrue(openTable(tableName))

    def test_5_Outlier1_MS_MTMFS_residual(self):
        '''Test 5: Check Outlier1.MS.MTMFS.residual'''
        tableName = 'Outlier1.MS.MTMFS.residual'
        self.assertTrue(openTable(tableName))

    def test_6_Outlier2_MS_MTMFS_flux(self):
        '''Test 6: Check Outlier2.MS.MTMFS.flux'''
        tableName = 'Outlier2.MS.MTMFS.flux'
        self.assertTrue(openTable(tableName))

    def test_7_Outlier2_MS_MTMFS_image(self):
        '''Test 7: Check Outlier2.MS.MTMFS.image'''
        tableName = 'Outlier2.MS.MTMFS.image'
        self.assertTrue(openTable(tableName))

    def test_8_Outlier2_MS_MTMFS_model(self):
        '''Test 8: Check Outlier2.MS.MTMFS.model'''
        tableName = 'Outlier2.MS.MTMFS.model'
        self.assertTrue(openTable(tableName))

    def test_9_Outlier2_MS_MTMFS_psf(self):
        '''Test 9: Check Outlier2.MS.MTMFS.psf'''
        tableName = 'Outlier2.MS.MTMFS.psf'
        self.assertTrue(openTable(tableName))

    def test_10_Outlier2_MS_MTMFS_residual(self):
        '''Test 10: Check Outlier2.MS.MTMFS.residual'''
        tableName = 'Outlier2.MS.MTMFS.residual'
        self.assertTrue(openTable(tableName))

    def test_11_SNR_MS_MTMFS_Main_flux(self):
        '''Test 11: Check SNR.MS.MTMFS-Main.flux'''
        tableName = 'SNR.MS.MTMFS-Main.flux'
        self.assertTrue(openTable(tableName))

    def test_12_SNR_MS_MTMFS_Main_image(self):
        '''Test 12: Check SNR.MS.MTMFS-Main.image'''
        tableName = 'SNR.MS.MTMFS-Main.image'
        self.assertTrue(openTable(tableName))

    def test_13_SNR_MS_MTMFS_Main_model(self):
        '''Test 13: Check SNR.MS.MTMFS-Main.model'''
        tableName = 'SNR.MS.MTMFS-Main.model'
        self.assertTrue(openTable(tableName))

    def test_14_SNR_MS_MTMFS_Main_psf(self):
        '''Test 14: Check SNR.MS.MTMFS-Main.psf'''
        tableName = 'SNR.MS.MTMFS-Main.psf'
        self.assertTrue(openTable(tableName))

    def test_15_SNR_MS_MTMFS_Main_residual(self):
        '''Test 15: Check SNR.MS.MTMFS-Main.residual'''
        tableName = 'SNR.MS.MTMFS-Main.residual'
        self.assertTrue(openTable(tableName))

    def test_16_SNR_G55_10s_MS_pbcorr_image(self):
        '''Test 16: Check SNR_G55_10s.MS.pbcorr.image'''
        tableName = 'SNR_G55_10s.MS.pbcorr.image'
        self.assertTrue(openTable(tableName))

    def test_17_SNR_G55_10s_MultiScale_flux(self):
        '''Test 17: Check SNR_G55_10s.MultiScale.flux'''
        tableName = 'SNR_G55_10s.MultiScale.flux'
        self.assertTrue(openTable(tableName))

    def test_18_SNR_G55_10s_MultiScale_image(self):
        '''Test 18: Check SNR_G55_10s.MultiScale.image'''
        tableName = 'SNR_G55_10s.MultiScale.image'
        self.assertTrue(openTable(tableName))

    def test_19_SNR_G55_10s_MultiScale_model(self):
        '''Test 19: Check SNR_G55_10s.MultiScale.model'''
        tableName = 'SNR_G55_10s.MultiScale.model'
        self.assertTrue(openTable(tableName))

    def test_20_SNR_G55_10s_MultiScale_psf(self):
        '''Test 20: Check SNR_G55_10s.MultiScale.psf'''
        tableName = 'SNR_G55_10s.MultiScale.psf'
        self.assertTrue(openTable(tableName))

    def test_21_SNR_G55_10s_MultiScale_residual(self):
        '''Test 21: Check SNR_G55_10s.MultiScale.residual'''
        tableName = 'SNR_G55_10s.MultiScale.residual'
        self.assertTrue(openTable(tableName))

    def test_22_SNR_G55_10s_Reg_Clean_niter10K_flux(self):
        '''Test 22: Check SNR_G55_10s.Reg.Clean.niter10K.flux'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter10K.flux'
        self.assertTrue(openTable(tableName))

    def test_23_SNR_G55_10s_Reg_Clean_niter10K_image(self):
        '''Test 23: Check SNR_G55_10s.Reg.Clean.niter10K.image'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter10K.image'
        self.assertTrue(openTable(tableName))

    def test_24_SNR_G55_10s_Reg_Clean_niter10K_model(self):
        '''Test 24: Check SNR_G55_10s.Reg.Clean.niter10K.model'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter10K.model'
        self.assertTrue(openTable(tableName))

    def test_25_SNR_G55_10s_Reg_Clean_niter10K_psf(self):
        '''Test 25: Check SNR_G55_10s.Reg.Clean.niter10K.psf'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter10K.psf'
        self.assertTrue(openTable(tableName))

    def test_26_SNR_G55_10s_Reg_Clean_niter10K_residual(self):
        '''Test 26: Check SNR_G55_10s.Reg.Clean.niter10K.residual'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter10K.residual'
        self.assertTrue(openTable(tableName))

    def test_27_SNR_G55_10s_Reg_Clean_niter1K_flux(self):
        '''Test 27: Check SNR_G55_10s.Reg.Clean.niter1K.flux'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter1K.flux'
        self.assertTrue(openTable(tableName))

    def test_28_SNR_G55_10s_Reg_Clean_niter1K_image(self):
        '''Test 28: Check SNR_G55_10s.Reg.Clean.niter1K.image'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter1K.image'
        self.assertTrue(openTable(tableName))

    def test_29_SNR_G55_10s_Reg_Clean_niter1K_model(self):
        '''Test 29: Check SNR_G55_10s.Reg.Clean.niter1K.model'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter1K.model'
        self.assertTrue(openTable(tableName))

    def test_30_SNR_G55_10s_Reg_Clean_niter1K_psf(self):
        '''Test 30: Check SNR_G55_10s.Reg.Clean.niter1K.psf'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter1K.psf'
        self.assertTrue(openTable(tableName))

    def test_31_SNR_G55_10s_Reg_Clean_niter1K_residual(self):
        '''Test 31: Check SNR_G55_10s.Reg.Clean.niter1K.residual'''
        tableName = 'SNR_G55_10s.Reg.Clean.niter1K.residual'
        self.assertTrue(openTable(tableName))

    def test_32_SNR_G55_10s_briggs_flux(self):
        '''Test 32: Check SNR_G55_10s.briggs.flux'''
        tableName = 'SNR_G55_10s.briggs.flux'
        self.assertTrue(openTable(tableName))

    def test_33_SNR_G55_10s_briggs_image(self):
        '''Test 33: Check SNR_G55_10s.briggs.image'''
        tableName = 'SNR_G55_10s.briggs.image'
        self.assertTrue(openTable(tableName))

    def test_34_SNR_G55_10s_briggs_model(self):
        '''Test 34: Check SNR_G55_10s.briggs.model'''
        tableName = 'SNR_G55_10s.briggs.model'
        self.assertTrue(openTable(tableName))

    def test_35_SNR_G55_10s_briggs_psf(self):
        '''Test 35: Check SNR_G55_10s.briggs.psf'''
        tableName = 'SNR_G55_10s.briggs.psf'
        self.assertTrue(openTable(tableName))

    def test_36_SNR_G55_10s_briggs_residual(self):
        '''Test 36: Check SNR_G55_10s.briggs.residual'''
        tableName = 'SNR_G55_10s.briggs.residual'
        self.assertTrue(openTable(tableName))


    def test_38_SNR_G55_10s_dirty_flux(self):
        '''Test 38: Check SNR_G55_10s.dirty.flux'''
        tableName = 'SNR_G55_10s.dirty.flux'
        self.assertTrue(openTable(tableName))

    def test_39_SNR_G55_10s_dirty_image(self):
        '''Test 39: Check SNR_G55_10s.dirty.image'''
        tableName = 'SNR_G55_10s.dirty.image'
        self.assertTrue(openTable(tableName))

    def test_40_SNR_G55_10s_dirty_model(self):
        '''Test 40: Check SNR_G55_10s.dirty.model'''
        tableName = 'SNR_G55_10s.dirty.model'
        self.assertTrue(openTable(tableName))

    def test_41_SNR_G55_10s_dirty_psf(self):
        '''Test 41: Check SNR_G55_10s.dirty.psf'''
        tableName = 'SNR_G55_10s.dirty.psf'
        self.assertTrue(openTable(tableName))

    def test_42_SNR_G55_10s_dirty_residual(self):
        '''Test 42: Check SNR_G55_10s.dirty.residual'''
        tableName = 'SNR_G55_10s.dirty.residual'
        self.assertTrue(openTable(tableName))

    def test_43_SNR_G55_10s_ms_MTMFS_flux(self):
        '''Test 43: Check SNR_G55_10s.ms.MTMFS.flux'''
        tableName = 'SNR_G55_10s.ms.MTMFS.flux'
        self.assertTrue(openTable(tableName))

    def test_44_SNR_G55_10s_ms_MTMFS_image_alpha(self):
        '''Test 44: Check SNR_G55_10s.ms.MTMFS.image.alpha'''
        tableName = 'SNR_G55_10s.ms.MTMFS.image.alpha'
        self.assertTrue(openTable(tableName))

    def test_45_SNR_G55_10s_ms_MTMFS_image_alpha_error(self):
        '''Test 45: Check SNR_G55_10s.ms.MTMFS.image.alpha.error'''
        tableName = 'SNR_G55_10s.ms.MTMFS.image.alpha.error'
        self.assertTrue(openTable(tableName))

    def test_46_SNR_G55_10s_ms_MTMFS_image_tt0(self):
        '''Test 46: Check SNR_G55_10s.ms.MTMFS.image.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.image.tt0'
        self.assertTrue(openTable(tableName))

    def test_47_SNR_G55_10s_ms_MTMFS_image_tt1(self):
        '''Test 47: Check SNR_G55_10s.ms.MTMFS.image.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.image.tt1'
        self.assertTrue(openTable(tableName))

    def test_48_SNR_G55_10s_ms_MTMFS_model_tt0(self):
        '''Test 48: Check SNR_G55_10s.ms.MTMFS.model.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.model.tt0'
        self.assertTrue(openTable(tableName))

    def test_49_SNR_G55_10s_ms_MTMFS_model_tt1(self):
        '''Test 49: Check SNR_G55_10s.ms.MTMFS.model.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.model.tt1'
        self.assertTrue(openTable(tableName))

    def test_50_SNR_G55_10s_ms_MTMFS_psf_tt0(self):
        '''Test 50: Check SNR_G55_10s.ms.MTMFS.psf.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.psf.tt0'
        self.assertTrue(openTable(tableName))

    def test_51_SNR_G55_10s_ms_MTMFS_psf_tt1(self):
        '''Test 51: Check SNR_G55_10s.ms.MTMFS.psf.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.psf.tt1'
        self.assertTrue(openTable(tableName))

    def test_52_SNR_G55_10s_ms_MTMFS_residual_tt0(self):
        '''Test 52: Check SNR_G55_10s.ms.MTMFS.residual.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.residual.tt0'
        self.assertTrue(openTable(tableName))

    def test_53_SNR_G55_10s_ms_MTMFS_residual_tt1(self):
        '''Test 53: Check SNR_G55_10s.ms.MTMFS.residual.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.residual.tt1'
        self.assertTrue(openTable(tableName))

    def test_54_SNR_G55_10s_ms_MTMFS_wProj_flux(self):
        '''Test 54: Check SNR_G55_10s.ms.MTMFS.wProj.flux'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.flux'
        self.assertTrue(openTable(tableName))

    def test_55_SNR_G55_10s_ms_MTMFS_wProj_image_alpha(self):
        '''Test 55: Check SNR_G55_10s.ms.MTMFS.wProj.image.alpha'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.image.alpha'
        self.assertTrue(openTable(tableName))

    def test_56_SNR_G55_10s_ms_MTMFS_wProj_image_alpha_error(self):
        '''Test 56: Check SNR_G55_10s.ms.MTMFS.wProj.image.alpha.error'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.image.alpha.error'
        self.assertTrue(openTable(tableName))

    def test_57_SNR_G55_10s_ms_MTMFS_wProj_image_tt0(self):
        '''Test 57: Check SNR_G55_10s.ms.MTMFS.wProj.image.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.image.tt0'
        self.assertTrue(openTable(tableName))

    def test_58_SNR_G55_10s_ms_MTMFS_wProj_image_tt0_Tb(self):
        '''Test 58: Check SNR_G55_10s.ms.MTMFS.wProj.image.tt0-Tb'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.image.tt0-Tb'
        self.assertTrue(openTable(tableName))

    def test_59_SNR_G55_10s_ms_MTMFS_wProj_image_tt1(self):
        '''Test 59: Check SNR_G55_10s.ms.MTMFS.wProj.image.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.image.tt1'
        self.assertTrue(openTable(tableName))

    def test_60_SNR_G55_10s_ms_MTMFS_wProj_model_tt0(self):
        '''Test 60: Check SNR_G55_10s.ms.MTMFS.wProj.model.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.model.tt0'
        self.assertTrue(openTable(tableName))

    def test_61_SNR_G55_10s_ms_MTMFS_wProj_model_tt1(self):
        '''Test 61: Check SNR_G55_10s.ms.MTMFS.wProj.model.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.model.tt1'
        self.assertTrue(openTable(tableName))

    def test_62_SNR_G55_10s_ms_MTMFS_wProj_pbcor_image_alpha(self):
        '''Test 62: Check SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.alpha'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.alpha'
        self.assertTrue(openTable(tableName))

    def test_63_SNR_G55_10s_ms_MTMFS_wProj_pbcor_image_alpha_error(self):
        '''Test 63: Check SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.alpha.error'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.alpha.error'
        self.assertTrue(openTable(tableName))

    def test_64_SNR_G55_10s_ms_MTMFS_wProj_pbcor_image_tt0(self):
        '''Test 64: Check SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.tt0'
        self.assertTrue(openTable(tableName))

    def test_65_SNR_G55_10s_ms_MTMFS_wProj_pbcor_image_tt1(self):
        '''Test 65: Check SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.pbcor.image.tt1'
        self.assertTrue(openTable(tableName))

    def test_67_SNR_G55_10s_ms_MTMFS_wProj_psf_tt0(self):
        '''Test 67: Check SNR_G55_10s.ms.MTMFS.wProj.psf.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.psf.tt0'
        self.assertTrue(openTable(tableName))

    def test_68_SNR_G55_10s_ms_MTMFS_wProj_psf_tt1(self):
        '''Test 68: Check SNR_G55_10s.ms.MTMFS.wProj.psf.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.psf.tt1'
        self.assertTrue(openTable(tableName))

    def test_69_SNR_G55_10s_ms_MTMFS_wProj_residual_tt0(self):
        '''Test 69: Check SNR_G55_10s.ms.MTMFS.wProj.residual.tt0'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.residual.tt0'
        self.assertTrue(openTable(tableName))

    def test_70_SNR_G55_10s_ms_MTMFS_wProj_residual_tt1(self):
        '''Test 70: Check SNR_G55_10s.ms.MTMFS.wProj.residual.tt1'''
        tableName = 'SNR_G55_10s.ms.MTMFS.wProj.residual.tt1'
        self.assertTrue(openTable(tableName))

    def test_71_SNR_G55_10s_ms_tclean_MTMFS_wProj_alpha(self):
        '''Test 71: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.alpha'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.alpha'
        self.assertTrue(openTable(tableName))

    def test_72_SNR_G55_10s_ms_tclean_MTMFS_wProj_alpha_error(self):
        '''Test 72: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.alpha.error'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.alpha.error'
        self.assertTrue(openTable(tableName))

    def test_73_SNR_G55_10s_ms_tclean_MTMFS_wProj_image_tt0(self):
        '''Test 73: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.image.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.image.tt0'
        self.assertTrue(openTable(tableName))

    def test_74_SNR_G55_10s_ms_tclean_MTMFS_wProj_image_tt1(self):
        '''Test 74: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.image.tt1'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.image.tt1'
        self.assertTrue(openTable(tableName))

    def test_75_SNR_G55_10s_ms_tclean_MTMFS_wProj_mask(self):
        '''Test 75: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.mask'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.mask'
        self.assertTrue(openTable(tableName))

    def test_76_SNR_G55_10s_ms_tclean_MTMFS_wProj_model_tt0(self):
        '''Test 76: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.model.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.model.tt0'
        self.assertTrue(openTable(tableName))

    def test_77_SNR_G55_10s_ms_tclean_MTMFS_wProj_model_tt1(self):
        '''Test 77: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.model.tt1'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.model.tt1'
        self.assertTrue(openTable(tableName))

    def test_78_SNR_G55_10s_ms_tclean_MTMFS_wProj_pb_tt0(self):
        '''Test 78: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.pb.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.pb.tt0'
        self.assertTrue(openTable(tableName))

    def test_79_SNR_G55_10s_ms_tclean_MTMFS_wProj_psf_tt0(self):
        '''Test 79: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt0'
        self.assertTrue(openTable(tableName))

    def test_80_SNR_G55_10s_ms_tclean_MTMFS_wProj_psf_tt1(self):
        '''Test 80: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt1'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt1'
        self.assertTrue(openTable(tableName))

    def test_81_SNR_G55_10s_ms_tclean_MTMFS_wProj_psf_tt2(self):
        '''Test 81: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt2'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.psf.tt2'
        self.assertTrue(openTable(tableName))

    def test_82_SNR_G55_10s_ms_tclean_MTMFS_wProj_residual_tt0(self):
        '''Test 82: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.residual.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.residual.tt0'
        self.assertTrue(openTable(tableName))

    def test_83_SNR_G55_10s_ms_tclean_MTMFS_wProj_residual_tt1(self):
        '''Test 83: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.residual.tt1'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.residual.tt1'
        self.assertTrue(openTable(tableName))

    def test_84_SNR_G55_10s_ms_tclean_MTMFS_wProj_sumwt_tt0(self):
        '''Test 84: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt0'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt0'
        self.assertTrue(openTable(tableName))

    def test_85_SNR_G55_10s_ms_tclean_MTMFS_wProj_sumwt_tt1(self):
        '''Test 85: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt1'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt1'
        self.assertTrue(openTable(tableName))

    def test_86_SNR_G55_10s_ms_tclean_MTMFS_wProj_sumwt_tt2(self):
        '''Test 86: Check SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt2'''
        tableName = 'SNR_G55_10s.ms.tclean.MTMFS.wProj.sumwt.tt2'
        self.assertTrue(openTable(tableName))

    def test_87_SNR_G55_10s_ms_wProj_flux(self):
        '''Test 87: Check SNR_G55_10s.ms.wProj.flux'''
        tableName = 'SNR_G55_10s.ms.wProj.flux'
        self.assertTrue(openTable(tableName))

    def test_88_SNR_G55_10s_ms_wProj_image(self):
        '''Test 88: Check SNR_G55_10s.ms.wProj.image'''
        tableName = 'SNR_G55_10s.ms.wProj.image'
        self.assertTrue(openTable(tableName))

    def test_89_SNR_G55_10s_ms_wProj_model(self):
        '''Test 89: Check SNR_G55_10s.ms.wProj.model'''
        tableName = 'SNR_G55_10s.ms.wProj.model'
        self.assertTrue(openTable(tableName))

    def test_90_SNR_G55_10s_ms_wProj_psf(self):
        '''Test 90: Check SNR_G55_10s.ms.wProj.psf'''
        tableName = 'SNR_G55_10s.ms.wProj.psf'
        self.assertTrue(openTable(tableName))

    def test_91_SNR_G55_10s_ms_wProj_residual(self):
        '''Test 91: Check SNR_G55_10s.ms.wProj.residual'''
        tableName = 'SNR_G55_10s.ms.wProj.residual'
        self.assertTrue(openTable(tableName))

    def test_92_SNR_G55_10s_natural_flux(self):
        '''Test 92: Check SNR_G55_10s.natural.flux'''
        tableName = 'SNR_G55_10s.natural.flux'
        self.assertTrue(openTable(tableName))

    def test_93_SNR_G55_10s_natural_image(self):
        '''Test 93: Check SNR_G55_10s.natural.image'''
        tableName = 'SNR_G55_10s.natural.image'
        self.assertTrue(openTable(tableName))

    def test_94_SNR_G55_10s_natural_model(self):
        '''Test 94: Check SNR_G55_10s.natural.model'''
        tableName = 'SNR_G55_10s.natural.model'
        self.assertTrue(openTable(tableName))

    def test_95_SNR_G55_10s_natural_psf(self):
        '''Test 95: Check SNR_G55_10s.natural.psf'''
        tableName = 'SNR_G55_10s.natural.psf'
        self.assertTrue(openTable(tableName))

    def test_96_SNR_G55_10s_natural_residual(self):
        '''Test 96: Check SNR_G55_10s.natural.residual'''
        tableName = 'SNR_G55_10s.natural.residual'
        self.assertTrue(openTable(tableName))

    def test_97_SNR_G55_10s_uniform_flux(self):
        '''Test 97: Check SNR_G55_10s.uniform.flux'''
        tableName = 'SNR_G55_10s.uniform.flux'
        self.assertTrue(openTable(tableName))

    def test_98_SNR_G55_10s_uniform_image(self):
        '''Test 98: Check SNR_G55_10s.uniform.image'''
        tableName = 'SNR_G55_10s.uniform.image'
        self.assertTrue(openTable(tableName))

    def test_99_SNR_G55_10s_uniform_model(self):
        '''Test 99: Check SNR_G55_10s.uniform.model'''
        tableName = 'SNR_G55_10s.uniform.model'
        self.assertTrue(openTable(tableName))

    def test_100_SNR_G55_10s_uniform_psf(self):
        '''Test 100: Check SNR_G55_10s.uniform.psf'''
        tableName = 'SNR_G55_10s.uniform.psf'
        self.assertTrue(openTable(tableName))

    def test_101_SNR_G55_10s_uniform_residual(self):
        '''Test 101: Check SNR_G55_10s.uniform.residual'''
        tableName = 'SNR_G55_10s.uniform.residual'
        self.assertTrue(openTable(tableName))

####################################################################################################

