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
    return [Test010_VLAIRC10216 ,Test020_VLAIRC10216]


class Test010_VLAIRC10216(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/day2_TDEM0003_10s_norx.tar.gz")
        tar.extractall()
        tar.close()


        #Fitting an average spectrum
        file = open("specfit.crtf", "w")
        file.write("#CRTFv0\n")
        file.write("box[[09:47:59.2, 13.16.24], [09:47:55.8, 13.17.09]]")
        file.close()

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_high_frequency_Spectral_Line_tutorial_-_IRC%2B10216'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_high_frequency_Spectral_Line_tutorial_-_IRC%2B10216'")

        time.sleep(5) # Allow extract time to download script


    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide:  VLA high frequency Spectral Line tutorial - IRC+10216'''

        exec(compile(open('VLAhighfrequencySpectralLinetutorial-IRC_2B10216.py', "rb").read(), 'VLAhighfrequencySpectralLinetutorial-IRC_2B10216.py', 'exec'))
                
        return True

class Test020_VLAIRC10216(unittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        os.system("rm -rf day2_TDEM0003_10s_norx")
        rmtables("IRC10216*")
        rmtables("*.cal")
        rmtables("*.bcal")
        rmtables("*.gcal")
        rmtables("pcal_ch19one_10min")
        os.system("rm -rf *.last")
        os.system("rm -rf specfit.crtf")

    def test_1_IRC10216(self):
        '''Test 1: Check IRC10216'''
        tableName = 'IRC10216'
        self.assertTrue(openTable(tableName))

    def test_2_IRC10216_36GHzcont_flux(self):
        '''Test 2: Check IRC10216.36GHzcont.flux'''
        tableName = 'IRC10216.36GHzcont.flux'
        self.assertTrue(openTable(tableName))

    def test_3_IRC10216_36GHzcont_image(self):
        '''Test 3: Check IRC10216.36GHzcont.image'''
        tableName = 'IRC10216.36GHzcont.image'
        self.assertTrue(openTable(tableName))

    def test_4_IRC10216_36GHzcont_model(self):
        '''Test 4: Check IRC10216.36GHzcont.model'''
        tableName = 'IRC10216.36GHzcont.model'
        self.assertTrue(openTable(tableName))

    def test_5_IRC10216_36GHzcont_psf(self):
        '''Test 5: Check IRC10216.36GHzcont.psf'''
        tableName = 'IRC10216.36GHzcont.psf'
        self.assertTrue(openTable(tableName))

    def test_6_IRC10216_36GHzcont_residual(self):
        '''Test 6: Check IRC10216.36GHzcont.residual'''
        tableName = 'IRC10216.36GHzcont.residual'
        self.assertTrue(openTable(tableName))

    def test_7_IRC10216_cont(self):
        '''Test 7: Check IRC10216.cont'''
        tableName = 'IRC10216.cont'
        self.assertTrue(openTable(tableName))

    def test_8_IRC10216_contsub(self):
        '''Test 8: Check IRC10216.contsub'''
        tableName = 'IRC10216.contsub'
        self.assertTrue(openTable(tableName))

    def test_9_IRC10216_contsub_cveled(self):
        '''Test 9: Check IRC10216.contsub-cveled'''
        tableName = 'IRC10216.contsub-cveled'
        self.assertTrue(openTable(tableName))

    def test_10_IRC10216_HC3N_cube_r0_5_flux(self):
        '''Test 10: Check IRC10216_HC3N.cube_r0.5.flux'''
        tableName = 'IRC10216_HC3N.cube_r0.5.flux'
        self.assertTrue(openTable(tableName))

    def test_11_IRC10216_HC3N_cube_r0_5_image(self):
        '''Test 11: Check IRC10216_HC3N.cube_r0.5.image'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image'
        self.assertTrue(openTable(tableName))

    def test_12_IRC10216_HC3N_cube_r0_5_image_extramoms_maximum(self):
        '''Test 12: Check IRC10216_HC3N.cube_r0.5.image.extramoms.maximum'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image.extramoms.maximum'
        self.assertTrue(openTable(tableName))

    def test_13_IRC10216_HC3N_cube_r0_5_image_extramoms_median(self):
        '''Test 13: Check IRC10216_HC3N.cube_r0.5.image.extramoms.median'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image.extramoms.median'
        self.assertTrue(openTable(tableName))

    def test_14_IRC10216_HC3N_cube_r0_5_image_extramoms_weighted_dispersion_coord(self):
        '''Test 14: Check IRC10216_HC3N.cube_r0.5.image.extramoms.weighted_dispersion_coord'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image.extramoms.weighted_dispersion_coord'
        self.assertTrue(openTable(tableName))

    def test_15_IRC10216_HC3N_cube_r0_5_image_mom0(self):
        '''Test 15: Check IRC10216_HC3N.cube_r0.5.image.mom0'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image.mom0'
        self.assertTrue(openTable(tableName))

    def test_16_IRC10216_HC3N_cube_r0_5_image_mom1(self):
        '''Test 16: Check IRC10216_HC3N.cube_r0.5.image.mom1'''
        tableName = 'IRC10216_HC3N.cube_r0.5.image.mom1'
        self.assertTrue(openTable(tableName))

    def test_17_IRC10216_HC3N_cube_r0_5_model(self):
        '''Test 17: Check IRC10216_HC3N.cube_r0.5.model'''
        tableName = 'IRC10216_HC3N.cube_r0.5.model'
        self.assertTrue(openTable(tableName))

    def test_18_IRC10216_HC3N_cube_r0_5_pselfcal(self):
        '''Test 18: Check IRC10216_HC3N.cube_r0.5.pselfcal'''
        tableName = 'IRC10216_HC3N.cube_r0.5.pselfcal'
        self.assertTrue(openTable(tableName))

    def test_19_IRC10216_HC3N_cube_r0_5_pselfcal_flux(self):
        '''Test 19: Check IRC10216_HC3N.cube_r0.5.pselfcal.flux'''
        tableName = 'IRC10216_HC3N.cube_r0.5.pselfcal.flux'
        self.assertTrue(openTable(tableName))

    def test_20_IRC10216_HC3N_cube_r0_5_pselfcal_mask(self):
        '''Test 20: Check IRC10216_HC3N.cube_r0.5.pselfcal.mask'''
        tableName = 'IRC10216_HC3N.cube_r0.5.pselfcal.mask'
        self.assertTrue(openTable(tableName))

    def test_21_IRC10216_HC3N_cube_r0_5_psf(self):
        '''Test 21: Check IRC10216_HC3N.cube_r0.5.psf'''
        tableName = 'IRC10216_HC3N.cube_r0.5.psf'
        self.assertTrue(openTable(tableName))

    def test_22_IRC10216_HC3N_cube_r0_5_residual(self):
        '''Test 22: Check IRC10216_HC3N.cube_r0.5.residual'''
        tableName = 'IRC10216_HC3N.cube_r0.5.residual'
        self.assertTrue(openTable(tableName))

    def test_23_IRC10216_SiS_cube_r0_5_flux(self):
        '''Test 23: Check IRC10216_SiS.cube_r0.5.flux'''
        tableName = 'IRC10216_SiS.cube_r0.5.flux'
        self.assertTrue(openTable(tableName))

    def test_24_IRC10216_SiS_cube_r0_5_image(self):
        '''Test 24: Check IRC10216_SiS.cube_r0.5.image'''
        tableName = 'IRC10216_SiS.cube_r0.5.image'
        self.assertTrue(openTable(tableName))

    def test_25_IRC10216_SiS_cube_r0_5_image_mom0(self):
        '''Test 25: Check IRC10216_SiS.cube_r0.5.image.mom0'''
        tableName = 'IRC10216_SiS.cube_r0.5.image.mom0'
        self.assertTrue(openTable(tableName))

    def test_26_IRC10216_SiS_cube_r0_5_image_mom1(self):
        '''Test 26: Check IRC10216_SiS.cube_r0.5.image.mom1'''
        tableName = 'IRC10216_SiS.cube_r0.5.image.mom1'
        self.assertTrue(openTable(tableName))

    def test_27_IRC10216_SiS_cube_r0_5_model(self):
        '''Test 27: Check IRC10216_SiS.cube_r0.5.model'''
        tableName = 'IRC10216_SiS.cube_r0.5.model'
        self.assertTrue(openTable(tableName))

    def test_28_IRC10216_SiS_cube_r0_5_pselfcal(self):
        '''Test 28: Check IRC10216_SiS.cube_r0.5.pselfcal'''
        tableName = 'IRC10216_SiS.cube_r0.5.pselfcal'
        self.assertTrue(openTable(tableName))

    def test_29_IRC10216_SiS_cube_r0_5_pselfcal_flux(self):
        '''Test 29: Check IRC10216_SiS.cube_r0.5.pselfcal.flux'''
        tableName = 'IRC10216_SiS.cube_r0.5.pselfcal.flux'
        self.assertTrue(openTable(tableName))

    def test_30_IRC10216_SiS_cube_r0_5_pselfcal_mask(self):
        '''Test 30: Check IRC10216_SiS.cube_r0.5.pselfcal.mask'''
        tableName = 'IRC10216_SiS.cube_r0.5.pselfcal.mask'
        self.assertTrue(openTable(tableName))

    def test_31_IRC10216_SiS_cube_r0_5_psf(self):
        '''Test 31: Check IRC10216_SiS.cube_r0.5.psf'''
        tableName = 'IRC10216_SiS.cube_r0.5.psf'
        self.assertTrue(openTable(tableName))

    def test_32_IRC10216_SiS_cube_r0_5_residual(self):
        '''Test 32: Check IRC10216_SiS.cube_r0.5.residual'''
        tableName = 'IRC10216_SiS.cube_r0.5.residual'
        self.assertTrue(openTable(tableName))

    def test_33_J0954(self):
        '''Test 33: Check J0954'''
        tableName = 'J0954'
        self.assertTrue(openTable(tableName))

    def test_34_amp_gcal(self):
        '''Test 34: Check amp.gcal'''
        tableName = 'amp.gcal'
        self.assertTrue(openTable(tableName))

    def test_35_amp_redo_gcal(self):
        '''Test 35: Check amp_redo.gcal'''
        tableName = 'amp_redo.gcal'
        self.assertTrue(openTable(tableName))

    def test_36_antpos_cal(self):
        '''Test 36: Check antpos.cal'''
        tableName = 'antpos.cal'
        self.assertTrue(openTable(tableName))

    def test_37_bandpass_bcal(self):
        '''Test 37: Check bandpass.bcal'''
        tableName = 'bandpass.bcal'
        self.assertTrue(openTable(tableName))

    def test_38_bandpass_redo_bcal(self):
        '''Test 38: Check bandpass_redo.bcal'''
        tableName = 'bandpass_redo.bcal'
        self.assertTrue(openTable(tableName))

    def test_39_bpphase_gcal(self):
        '''Test 39: Check bpphase.gcal'''
        tableName = 'bpphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_40_bpphase_redo_gcal(self):
        '''Test 40: Check bpphase_redo.gcal'''
        tableName = 'bpphase_redo.gcal'
        self.assertTrue(openTable(tableName))

    def test_41_delays_cal(self):
        '''Test 41: Check delays.cal'''
        tableName = 'delays.cal'
        self.assertTrue(openTable(tableName))

    def test_42_flux_cal(self):
        '''Test 42: Check flux.cal'''
        tableName = 'flux.cal'
        self.assertTrue(openTable(tableName))

    def test_43_flux_redo_cal(self):
        '''Test 43: Check flux_redo.cal'''
        tableName = 'flux_redo.cal'
        self.assertTrue(openTable(tableName))

    def test_44_gaincurve_cal(self):
        '''Test 44: Check gaincurve.cal'''
        tableName = 'gaincurve.cal'
        self.assertTrue(openTable(tableName))

    def test_45_intphase_gcal(self):
        '''Test 45: Check intphase.gcal'''
        tableName = 'intphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_46_intphase_redo_gcal(self):
        '''Test 46: Check intphase_redo.gcal'''
        tableName = 'intphase_redo.gcal'
        self.assertTrue(openTable(tableName))

    def test_47_opacity_cal(self):
        '''Test 47: Check opacity.cal'''
        tableName = 'opacity.cal'
        self.assertTrue(openTable(tableName))

    def test_48_pcal_ch19one_10min(self):
        '''Test 48: Check pcal_ch19one_10min'''
        tableName = 'pcal_ch19one_10min'
        self.assertTrue(openTable(tableName))

    def test_49_scanphase_gcal(self):
        '''Test 49: Check scanphase.gcal'''
        tableName = 'scanphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_50_scanphase_redo_gcal(self):
        '''Test 50: Check scanphase_redo.gcal'''
        tableName = 'scanphase_redo.gcal'
        self.assertTrue(openTable(tableName))

####################################################################################################

    def test_51_day2_TDEM0003_10s_norx_plotweather_png(self):
        '''Test 51: Check day2_TDEM0003_10s_norx.plotweather.png'''
        self.assertTrue(assert_file('day2_TDEM0003_10s_norx.plotweather.png'))

    #def test_52_compare_day2_TDEM0003_10s_norx_plotweather_png(self):
    #    '''Test 52: Image Comparison day2_TDEM0003_10s_norx.plotweather.png'''
    #    imageName = 'day2_TDEM0003_10s_norx.plotweather.png'
    #    self.assertTrue(compareImages('day2_TDEM0003_10s_norx.plotweather.png'))

