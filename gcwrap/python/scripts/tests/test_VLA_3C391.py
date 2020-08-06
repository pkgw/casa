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
    return [Test010_VLAContinuum3C391,Test020_VLAContinuum3C391,Test021_VLAContinuum3C391]


class Test010_VLAContinuum3C391(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/3c391_ctm_mosaic_10s_spw0.ms.tgz")
        tar.extractall()
        tar.close()

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391'")

        time.sleep(5) # Allow extract time to download script

        lines = open('VLAContinuumTutorial3C391.py')
        file = open("newfile.txt", "w")
        for line in lines:

            pattern = r'''niter\ *=\ *(25000)'''
            if re.search(pattern,line):
                line = re.sub( pattern, 'niter=15000', line )

            pattern = r'''niter\ *=\ *(5000)'''
            if re.search(pattern,line):
                line = re.sub( pattern, 'niter=3000', line )

            file.write(line)
        file.close()
        os.remove('VLAContinuumTutorial3C391.py')
        os.rename("newfile.txt",'VLAContinuumTutorial3C391.py')

    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide:  VLA Continuum Tutorial 3C391'''

        exec(compile(open('VLAContinuumTutorial3C391.py', "rb").read(), 'VLAContinuumTutorial3C391.py', 'exec'))
                
        return True

class Test020_VLAContinuum3C391(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass

    def test_1_3c391_ctm_mosaic_10s_spw0_B0(self):
        '''Test 1: Check 3c391_ctm_mosaic_10s_spw0.B0'''
        tableName = '3c391_ctm_mosaic_10s_spw0.B0'
        self.assertTrue(openTable(tableName))

    def test_2_3c391_ctm_mosaic_10s_spw0_D1(self):
        '''Test 2: Check 3c391_ctm_mosaic_10s_spw0.D1'''
        tableName = '3c391_ctm_mosaic_10s_spw0.D1'
        self.assertTrue(openTable(tableName))

    def test_3_3c391_ctm_mosaic_10s_spw0_D2(self):
        '''Test 3: Check 3c391_ctm_mosaic_10s_spw0.D2'''
        tableName = '3c391_ctm_mosaic_10s_spw0.D2'
        self.assertTrue(openTable(tableName))

    def test_4_3c391_ctm_mosaic_10s_spw0_G0(self):
        '''Test 4: Check 3c391_ctm_mosaic_10s_spw0.G0'''
        tableName = '3c391_ctm_mosaic_10s_spw0.G0'
        self.assertTrue(openTable(tableName))

    def test_5_3c391_ctm_mosaic_10s_spw0_G0all(self):
        '''Test 5: Check 3c391_ctm_mosaic_10s_spw0.G0all'''
        tableName = '3c391_ctm_mosaic_10s_spw0.G0all'
        self.assertTrue(openTable(tableName))

    def test_6_3c391_ctm_mosaic_10s_spw0_G1(self):
        '''Test 6: Check 3c391_ctm_mosaic_10s_spw0.G1'''
        tableName = '3c391_ctm_mosaic_10s_spw0.G1'
        self.assertTrue(openTable(tableName))

    def test_7_3c391_ctm_mosaic_10s_spw0_K0(self):
        '''Test 7: Check 3c391_ctm_mosaic_10s_spw0.K0'''
        tableName = '3c391_ctm_mosaic_10s_spw0.K0'
        self.assertTrue(openTable(tableName))

    def test_8_3c391_ctm_mosaic_10s_spw0_Kcross(self):
        '''Test 8: Check 3c391_ctm_mosaic_10s_spw0.Kcross'''
        tableName = '3c391_ctm_mosaic_10s_spw0.Kcross'
        self.assertTrue(openTable(tableName))

    def test_9_3c391_ctm_mosaic_10s_spw0_X1(self):
        '''Test 9: Check 3c391_ctm_mosaic_10s_spw0.X1'''
        tableName = '3c391_ctm_mosaic_10s_spw0.X1'
        self.assertTrue(openTable(tableName))

    def test_10_3c391_ctm_mosaic_10s_spw0_antpos(self):
        '''Test 10: Check 3c391_ctm_mosaic_10s_spw0.antpos'''
        tableName = '3c391_ctm_mosaic_10s_spw0.antpos'
        self.assertTrue(openTable(tableName))

    def test_11_3c391_ctm_mosaic_10s_spw0_fluxscale1(self):
        '''Test 11: Check 3c391_ctm_mosaic_10s_spw0.fluxscale1'''
        tableName = '3c391_ctm_mosaic_10s_spw0.fluxscale1'
        self.assertTrue(openTable(tableName))

    def test_12_3c391_ctm_mosaic_10s_spw0_ms(self):
        '''Test 12: Check 3c391_ctm_mosaic_10s_spw0.ms'''
        tableName = '3c391_ctm_mosaic_10s_spw0.ms'
        self.assertTrue(openTable(tableName))

    def test_13_3c391_ctm_mosaic_spw0_ms(self):
        '''Test 13: Check 3c391_ctm_mosaic_spw0.ms'''
        tableName = '3c391_ctm_mosaic_spw0.ms'
        self.assertTrue(openTable(tableName))

    def test_14_3c391_ctm_mosaic_spw0_selfcal1(self):
        '''Test 14: Check 3c391_ctm_mosaic_spw0.selfcal1'''
        tableName = '3c391_ctm_mosaic_spw0.selfcal1'
        self.assertTrue(openTable(tableName))

    def test_15_3c391_ctm_spw0_F(self):
        '''Test 15: Check 3c391_ctm_spw0.F'''
        tableName = '3c391_ctm_spw0.F'
        self.assertTrue(openTable(tableName))

    def test_16_3c391_ctm_spw0_I(self):
        '''Test 16: Check 3c391_ctm_spw0.I'''
        tableName = '3c391_ctm_spw0.I'
        self.assertTrue(openTable(tableName))

    def test_17_3c391_ctm_spw0_P(self):
        '''Test 17: Check 3c391_ctm_spw0.P'''
        tableName = '3c391_ctm_spw0.P'
        self.assertTrue(openTable(tableName))

    def test_18_3c391_ctm_spw0_P_unbias(self):
        '''Test 18: Check 3c391_ctm_spw0.P_unbias'''
        tableName = '3c391_ctm_spw0.P_unbias'
        self.assertTrue(openTable(tableName))

    def test_19_3c391_ctm_spw0_Q(self):
        '''Test 19: Check 3c391_ctm_spw0.Q'''
        tableName = '3c391_ctm_spw0.Q'
        self.assertTrue(openTable(tableName))

    def test_20_3c391_ctm_spw0_Qflux(self):
        '''Test 20: Check 3c391_ctm_spw0.Qflux'''
        tableName = '3c391_ctm_spw0.Qflux'
        self.assertTrue(openTable(tableName))

    def test_21_3c391_ctm_spw0_U(self):
        '''Test 21: Check 3c391_ctm_spw0.U'''
        tableName = '3c391_ctm_spw0.U'
        self.assertTrue(openTable(tableName))

    def test_22_3c391_ctm_spw0_V(self):
        '''Test 22: Check 3c391_ctm_spw0.V'''
        tableName = '3c391_ctm_spw0.V'
        self.assertTrue(openTable(tableName))

    def test_23_3c391_ctm_spw0_X(self):
        '''Test 23: Check 3c391_ctm_spw0.X'''
        tableName = '3c391_ctm_spw0.X'
        self.assertTrue(openTable(tableName))

    def test_24_3c391_ctm_spw0_pbcorI(self):
        '''Test 24: Check 3c391_ctm_spw0.pbcorI'''
        tableName = '3c391_ctm_spw0.pbcorI'
        self.assertTrue(openTable(tableName))

    def test_25_3c391_ctm_spw0_pbcorP(self):
        '''Test 25: Check 3c391_ctm_spw0.pbcorP'''
        tableName = '3c391_ctm_spw0.pbcorP'
        self.assertTrue(openTable(tableName))

    def test_26_3c391_ctm_spw0_I_flux(self):
        '''Test 26: Check 3c391_ctm_spw0_I.flux'''
        tableName = '3c391_ctm_spw0_I.flux'
        self.assertTrue(openTable(tableName))

    def test_27_3c391_ctm_spw0_I_flux_pbcoverage(self):
        '''Test 27: Check 3c391_ctm_spw0_I.flux.pbcoverage'''
        tableName = '3c391_ctm_spw0_I.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_28_3c391_ctm_spw0_I_image(self):
        '''Test 28: Check 3c391_ctm_spw0_I.image'''
        tableName = '3c391_ctm_spw0_I.image'
        self.assertTrue(openTable(tableName))

    def test_29_3c391_ctm_spw0_I_model(self):
        '''Test 29: Check 3c391_ctm_spw0_I.model'''
        tableName = '3c391_ctm_spw0_I.model'
        self.assertTrue(openTable(tableName))

    def test_30_3c391_ctm_spw0_I_psf(self):
        '''Test 30: Check 3c391_ctm_spw0_I.psf'''
        tableName = '3c391_ctm_spw0_I.psf'
        self.assertTrue(openTable(tableName))

    def test_31_3c391_ctm_spw0_I_residual(self):
        '''Test 31: Check 3c391_ctm_spw0_I.residual'''
        tableName = '3c391_ctm_spw0_I.residual'
        self.assertTrue(openTable(tableName))

    def test_32_3c391_ctm_spw0_IQUV_flux(self):
        '''Test 32: Check 3c391_ctm_spw0_IQUV.flux'''
        tableName = '3c391_ctm_spw0_IQUV.flux'
        self.assertTrue(openTable(tableName))

    def test_33_3c391_ctm_spw0_IQUV_flux_pbcoverage(self):
        '''Test 33: Check 3c391_ctm_spw0_IQUV.flux.pbcoverage'''
        tableName = '3c391_ctm_spw0_IQUV.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_34_3c391_ctm_spw0_IQUV_image(self):
        '''Test 34: Check 3c391_ctm_spw0_IQUV.image'''
        tableName = '3c391_ctm_spw0_IQUV.image'
        self.assertTrue(openTable(tableName))

    def test_35_3c391_ctm_spw0_IQUV_model(self):
        '''Test 35: Check 3c391_ctm_spw0_IQUV.model'''
        tableName = '3c391_ctm_spw0_IQUV.model'
        self.assertTrue(openTable(tableName))

    def test_36_3c391_ctm_spw0_IQUV_pbcorimage(self):
        '''Test 36: Check 3c391_ctm_spw0_IQUV.pbcorimage'''
        tableName = '3c391_ctm_spw0_IQUV.pbcorimage'
        self.assertTrue(openTable(tableName))

    def test_37_3c391_ctm_spw0_IQUV_psf(self):
        '''Test 37: Check 3c391_ctm_spw0_IQUV.psf'''
        tableName = '3c391_ctm_spw0_IQUV.psf'
        self.assertTrue(openTable(tableName))

    def test_38_3c391_ctm_spw0_IQUV_residual(self):
        '''Test 38: Check 3c391_ctm_spw0_IQUV.residual'''
        tableName = '3c391_ctm_spw0_IQUV.residual'
        self.assertTrue(openTable(tableName))

    def test_39_3c391_ctm_spw0_IQUV_selfcal1_flux(self):
        '''Test 39: Check 3c391_ctm_spw0_IQUV_selfcal1.flux'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.flux'
        self.assertTrue(openTable(tableName))

    def test_40_3c391_ctm_spw0_IQUV_selfcal1_flux_pbcoverage(self):
        '''Test 40: Check 3c391_ctm_spw0_IQUV_selfcal1.flux.pbcoverage'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_41_3c391_ctm_spw0_IQUV_selfcal1_image(self):
        '''Test 41: Check 3c391_ctm_spw0_IQUV_selfcal1.image'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.image'
        self.assertTrue(openTable(tableName))

    def test_42_3c391_ctm_spw0_IQUV_selfcal1_model(self):
        '''Test 42: Check 3c391_ctm_spw0_IQUV_selfcal1.model'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.model'
        self.assertTrue(openTable(tableName))

    def test_43_3c391_ctm_spw0_IQUV_selfcal1_psf(self):
        '''Test 43: Check 3c391_ctm_spw0_IQUV_selfcal1.psf'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.psf'
        self.assertTrue(openTable(tableName))

    def test_44_3c391_ctm_spw0_IQUV_selfcal1_residual(self):
        '''Test 44: Check 3c391_ctm_spw0_IQUV_selfcal1.residual'''
        tableName = '3c391_ctm_spw0_IQUV_selfcal1.residual'
        self.assertTrue(openTable(tableName))

    def test_45_3c391_ctm_spw0_ms_I_flux(self):
        '''Test 45: Check 3c391_ctm_spw0_ms_I.flux'''
        tableName = '3c391_ctm_spw0_ms_I.flux'
        self.assertTrue(openTable(tableName))

    def test_46_3c391_ctm_spw0_ms_I_flux_pbcoverage(self):
        '''Test 46: Check 3c391_ctm_spw0_ms_I.flux.pbcoverage'''
        tableName = '3c391_ctm_spw0_ms_I.flux.pbcoverage'
        self.assertTrue(openTable(tableName))

    def test_47_3c391_ctm_spw0_ms_I_image(self):
        '''Test 47: Check 3c391_ctm_spw0_ms_I.image'''
        tableName = '3c391_ctm_spw0_ms_I.image'
        self.assertTrue(openTable(tableName))

    def test_48_3c391_ctm_spw0_ms_I_model(self):
        '''Test 48: Check 3c391_ctm_spw0_ms_I.model'''
        tableName = '3c391_ctm_spw0_ms_I.model'
        self.assertTrue(openTable(tableName))

    def test_49_3c391_ctm_spw0_ms_I_psf(self):
        '''Test 49: Check 3c391_ctm_spw0_ms_I.psf'''
        tableName = '3c391_ctm_spw0_ms_I.psf'
        self.assertTrue(openTable(tableName))

    def test_50_3c391_ctm_spw0_ms_I_residual(self):
        '''Test 50: Check 3c391_ctm_spw0_ms_I.residual'''
        tableName = '3c391_ctm_spw0_ms_I.residual'
        self.assertTrue(openTable(tableName))

####################################################################################################

class Test021_VLAContinuum3C391(unittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):

        rmtables("3c391_ctm_*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.flagversions")


    def test_16_Flux_density_3c391_ctm_spw0_IQUV(self):
        '''Flux density for 3c391_ctm_spw0_IQUV'''
        mystat = imstat(imagename='3c391_ctm_spw0_IQUV.pbcorimage',stokes='')
        peak_flux_density = float(mystat['max'][0])
        expected = 0.15447656810283661 
        print("Expected: %s , Actual: %s"%(expected,peak_flux_density))
        assert 0.146 <= peak_flux_density <= 0.163, "Error in Flux density for 3c391_ctm_spw0 with 0.146 <= %s <= 0.163"%(peak_flux_density)



