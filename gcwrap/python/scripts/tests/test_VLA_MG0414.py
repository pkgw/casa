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
    return [Test010_VLAMG0414,Test020_VLAMG0414]


class Test010_VLAMG0414(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/MG0414_d1_data.ms.tgz")
        tar.extractall()
        tar.close()


        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/MG0414%2B0534_P-band_Spectral_Line_Tutorial'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/MG0414%2B0534_P-band_Spectral_Line_Tutorial'")

        time.sleep(5) # Allow extract time to download script

        lines = open('MG0414_2B0534P-bandSpectralLineTutorial.py')
        file = open("newfile.txt", "w")
        for line in lines:

            if all(["tec_image" in line, "tec_rms_image" in line, "plotname" not in line]) :
                continue

            file.write(line)
        file.close()
        os.remove('MG0414_2B0534P-bandSpectralLineTutorial.py')
        os.rename("newfile.txt",'MG0414_2B0534P-bandSpectralLineTutorial.py')
    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide: MG0414+0534 P-band Spectral Line Tutorial'''

        exec(compile(open('MG0414_2B0534P-bandSpectralLineTutorial.py', "rb").read(), 'MG0414_2B0534P-bandSpectralLineTutorial.py', 'exec'))
                
        return True

class Test020_VLAMG0414(unittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        rmtables("MG0414_d1*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.flagversions")
        os.system("rm -rf *.cal")
        os.system("rm -rf *.gcal")

    def test_1_MG0414_d1_calibrated_ms(self):
        '''Test 1: Check MG0414_d1_calibrated.ms'''
        tableName = 'MG0414_d1_calibrated.ms'
        self.assertTrue(openTable(tableName))

    def test_2_MG0414_d1_cont_R03_image(self):
        '''Test 2: Check MG0414_d1_cont_R03.image'''
        tableName = 'MG0414_d1_cont_R03.image'
        self.assertTrue(openTable(tableName))

    def test_3_MG0414_d1_cont_R03_mask(self):
        '''Test 3: Check MG0414_d1_cont_R03.mask'''
        tableName = 'MG0414_d1_cont_R03.mask'
        self.assertTrue(openTable(tableName))

    def test_4_MG0414_d1_cont_R03_model(self):
        '''Test 4: Check MG0414_d1_cont_R03.model'''
        tableName = 'MG0414_d1_cont_R03.model'
        self.assertTrue(openTable(tableName))

    def test_5_MG0414_d1_cont_R03_pb(self):
        '''Test 5: Check MG0414_d1_cont_R03.pb'''
        tableName = 'MG0414_d1_cont_R03.pb'
        self.assertTrue(openTable(tableName))

    def test_6_MG0414_d1_cont_R03_psf(self):
        '''Test 6: Check MG0414_d1_cont_R03.psf'''
        tableName = 'MG0414_d1_cont_R03.psf'
        self.assertTrue(openTable(tableName))

    def test_7_MG0414_d1_cont_R03_residual(self):
        '''Test 7: Check MG0414_d1_cont_R03.residual'''
        tableName = 'MG0414_d1_cont_R03.residual'
        self.assertTrue(openTable(tableName))

    def test_8_MG0414_d1_cont_R03_sumwt(self):
        '''Test 8: Check MG0414_d1_cont_R03.sumwt'''
        tableName = 'MG0414_d1_cont_R03.sumwt'
        self.assertTrue(openTable(tableName))

    def test_9_MG0414_d1_cont_sc1_R03_image(self):
        '''Test 9: Check MG0414_d1_cont_sc1.R03.image'''
        tableName = 'MG0414_d1_cont_sc1.R03.image'
        self.assertTrue(openTable(tableName))

    def test_10_MG0414_d1_cont_sc1_R03_mask(self):
        '''Test 10: Check MG0414_d1_cont_sc1.R03.mask'''
        tableName = 'MG0414_d1_cont_sc1.R03.mask'
        self.assertTrue(openTable(tableName))

    def test_11_MG0414_d1_cont_sc1_R03_model(self):
        '''Test 11: Check MG0414_d1_cont_sc1.R03.model'''
        tableName = 'MG0414_d1_cont_sc1.R03.model'
        self.assertTrue(openTable(tableName))

    def test_12_MG0414_d1_cont_sc1_R03_pb(self):
        '''Test 12: Check MG0414_d1_cont_sc1.R03.pb'''
        tableName = 'MG0414_d1_cont_sc1.R03.pb'
        self.assertTrue(openTable(tableName))

    def test_13_MG0414_d1_cont_sc1_R03_psf(self):
        '''Test 13: Check MG0414_d1_cont_sc1.R03.psf'''
        tableName = 'MG0414_d1_cont_sc1.R03.psf'
        self.assertTrue(openTable(tableName))

    def test_14_MG0414_d1_cont_sc1_R03_residual(self):
        '''Test 14: Check MG0414_d1_cont_sc1.R03.residual'''
        tableName = 'MG0414_d1_cont_sc1.R03.residual'
        self.assertTrue(openTable(tableName))

    def test_15_MG0414_d1_cont_sc1_R03_sumwt(self):
        '''Test 15: Check MG0414_d1_cont_sc1.R03.sumwt'''
        tableName = 'MG0414_d1_cont_sc1.R03.sumwt'
        self.assertTrue(openTable(tableName))

    def test_16_MG0414_d1_cont_sc2_R03_image(self):
        '''Test 16: Check MG0414_d1_cont_sc2_R03.image'''
        tableName = 'MG0414_d1_cont_sc2_R03.image'
        self.assertTrue(openTable(tableName))

    def test_17_MG0414_d1_cont_sc2_R03_mask(self):
        '''Test 17: Check MG0414_d1_cont_sc2_R03.mask'''
        tableName = 'MG0414_d1_cont_sc2_R03.mask'
        self.assertTrue(openTable(tableName))

    def test_18_MG0414_d1_cont_sc2_R03_model(self):
        '''Test 18: Check MG0414_d1_cont_sc2_R03.model'''
        tableName = 'MG0414_d1_cont_sc2_R03.model'
        self.assertTrue(openTable(tableName))

    def test_19_MG0414_d1_cont_sc2_R03_pb(self):
        '''Test 19: Check MG0414_d1_cont_sc2_R03.pb'''
        tableName = 'MG0414_d1_cont_sc2_R03.pb'
        self.assertTrue(openTable(tableName))

    def test_20_MG0414_d1_cont_sc2_R03_psf(self):
        '''Test 20: Check MG0414_d1_cont_sc2_R03.psf'''
        tableName = 'MG0414_d1_cont_sc2_R03.psf'
        self.assertTrue(openTable(tableName))

    def test_21_MG0414_d1_cont_sc2_R03_residual(self):
        '''Test 21: Check MG0414_d1_cont_sc2_R03.residual'''
        tableName = 'MG0414_d1_cont_sc2_R03.residual'
        self.assertTrue(openTable(tableName))

    def test_22_MG0414_d1_cont_sc2_R03_sumwt(self):
        '''Test 22: Check MG0414_d1_cont_sc2_R03.sumwt'''
        tableName = 'MG0414_d1_cont_sc2_R03.sumwt'
        self.assertTrue(openTable(tableName))

    def test_23_MG0414_d1_cont_sc3ap_R03_image(self):
        '''Test 23: Check MG0414_d1_cont_sc3ap_R03.image'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.image'
        self.assertTrue(openTable(tableName))

    def test_24_MG0414_d1_cont_sc3ap_R03_mask(self):
        '''Test 24: Check MG0414_d1_cont_sc3ap_R03.mask'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.mask'
        self.assertTrue(openTable(tableName))

    def test_25_MG0414_d1_cont_sc3ap_R03_model(self):
        '''Test 25: Check MG0414_d1_cont_sc3ap_R03.model'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.model'
        self.assertTrue(openTable(tableName))

    def test_26_MG0414_d1_cont_sc3ap_R03_pb(self):
        '''Test 26: Check MG0414_d1_cont_sc3ap_R03.pb'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.pb'
        self.assertTrue(openTable(tableName))

    def test_27_MG0414_d1_cont_sc3ap_R03_psf(self):
        '''Test 27: Check MG0414_d1_cont_sc3ap_R03.psf'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.psf'
        self.assertTrue(openTable(tableName))

    def test_28_MG0414_d1_cont_sc3ap_R03_residual(self):
        '''Test 28: Check MG0414_d1_cont_sc3ap_R03.residual'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.residual'
        self.assertTrue(openTable(tableName))

    def test_29_MG0414_d1_cont_sc3ap_R03_sumwt(self):
        '''Test 29: Check MG0414_d1_cont_sc3ap_R03.sumwt'''
        tableName = 'MG0414_d1_cont_sc3ap_R03.sumwt'
        self.assertTrue(openTable(tableName))

    def test_30_MG0414_d1_data_ms_IGS_RMS_TEC_im(self):
        '''Test 30: Check MG0414_d1_data.ms.IGS_RMS_TEC.im'''
        tableName = 'MG0414_d1_data.ms.IGS_RMS_TEC.im'
        self.assertTrue(openTable(tableName))

    def test_31_MG0414_d1_data_ms_IGS_TEC_im(self):
        '''Test 31: Check MG0414_d1_data.ms.IGS_TEC.im'''
        tableName = 'MG0414_d1_data.ms.IGS_TEC.im'
        self.assertTrue(openTable(tableName))

    def test_32_MG0414_d1_line_sc3ap_vel_R03_image(self):
        '''Test 32: Check MG0414_d1_line_sc3ap_vel_R03.image'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.image'
        self.assertTrue(openTable(tableName))

    def test_33_MG0414_d1_line_sc3ap_vel_R03_image_pbcor(self):
        '''Test 33: Check MG0414_d1_line_sc3ap_vel_R03.image.pbcor'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.image.pbcor'
        self.assertTrue(openTable(tableName))

    def test_34_MG0414_d1_line_sc3ap_vel_R03_mask(self):
        '''Test 34: Check MG0414_d1_line_sc3ap_vel_R03.mask'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.mask'
        self.assertTrue(openTable(tableName))

    def test_35_MG0414_d1_line_sc3ap_vel_R03_model(self):
        '''Test 35: Check MG0414_d1_line_sc3ap_vel_R03.model'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.model'
        self.assertTrue(openTable(tableName))

    def test_36_MG0414_d1_line_sc3ap_vel_R03_pb(self):
        '''Test 36: Check MG0414_d1_line_sc3ap_vel_R03.pb'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.pb'
        self.assertTrue(openTable(tableName))

    def test_37_MG0414_d1_line_sc3ap_vel_R03_psf(self):
        '''Test 37: Check MG0414_d1_line_sc3ap_vel_R03.psf'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.psf'
        self.assertTrue(openTable(tableName))

    def test_38_MG0414_d1_line_sc3ap_vel_R03_residual(self):
        '''Test 38: Check MG0414_d1_line_sc3ap_vel_R03.residual'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.residual'
        self.assertTrue(openTable(tableName))

    def test_39_MG0414_d1_line_sc3ap_vel_R03_sumwt(self):
        '''Test 39: Check MG0414_d1_line_sc3ap_vel_R03.sumwt'''
        tableName = 'MG0414_d1_line_sc3ap_vel_R03.sumwt'
        self.assertTrue(openTable(tableName))

    def test_40_MG0414_d1_sc1_ms(self):
        '''Test 40: Check MG0414_d1_sc1.ms'''
        tableName = 'MG0414_d1_sc1.ms'
        self.assertTrue(openTable(tableName))

    def test_41_MG0414_d1_sc2_ms(self):
        '''Test 41: Check MG0414_d1_sc2.ms'''
        tableName = 'MG0414_d1_sc2.ms'
        self.assertTrue(openTable(tableName))

    def test_42_MG0414_d1_sc3ap_ms(self):
        '''Test 42: Check MG0414_d1_sc3ap.ms'''
        tableName = 'MG0414_d1_sc3ap.ms'
        self.assertTrue(openTable(tableName))

    def test_43_amp_gcal(self):
        '''Test 43: Check amp.gcal'''
        tableName = 'amp.gcal'
        self.assertTrue(openTable(tableName))

    def test_44_antpos_cal(self):
        '''Test 44: Check antpos.cal'''
        tableName = 'antpos.cal'
        self.assertTrue(openTable(tableName))

    def test_45_bandpass_cal(self):
        '''Test 45: Check bandpass.cal'''
        tableName = 'bandpass.cal'
        self.assertTrue(openTable(tableName))

    def test_46_bpphase_gcal(self):
        '''Test 46: Check bpphase.gcal'''
        tableName = 'bpphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_47_delays_cal(self):
        '''Test 47: Check delays.cal'''
        tableName = 'delays.cal'
        self.assertTrue(openTable(tableName))

    def test_49_intphase_gcal(self):
        '''Test 49: Check intphase.gcal'''
        tableName = 'intphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_50_rq_cal(self):
        '''Test 50: Check rq.cal'''
        tableName = 'rq.cal'
        self.assertTrue(openTable(tableName))

    def test_51_sc1_gcal(self):
        '''Test 51: Check sc1.gcal'''
        tableName = 'sc1.gcal'
        self.assertTrue(openTable(tableName))

    def test_52_sc2_gcal(self):
        '''Test 52: Check sc2.gcal'''
        tableName = 'sc2.gcal'
        self.assertTrue(openTable(tableName))

    def test_53_sc3ap_gcal(self):
        '''Test 53: Check sc3ap.gcal'''
        tableName = 'sc3ap.gcal'
        self.assertTrue(openTable(tableName))

    def test_54_scanphase_gcal(self):
        '''Test 54: Check scanphase.gcal'''
        tableName = 'scanphase.gcal'
        self.assertTrue(openTable(tableName))

    def test_55_tecim_cal(self):
        '''Test 55: Check tecim.cal'''
        tableName = 'tecim.cal'
        self.assertTrue(openTable(tableName))

####################################################################################################

    def test_56_MG0414_d1_data_ms_IGS_TEC_at_site_png(self):
        '''Test 56: Check MG0414_d1_data.ms.IGS_TEC_at_site.png'''
        self.assertTrue(assert_file('MG0414_d1_data.ms.IGS_TEC_at_site.png'))

    #def test_57_compare_MG0414_d1_data_ms_IGS_TEC_at_site_png(self):
    #    '''Test 57: Image Comparison MG0414_d1_data.ms.IGS_TEC_at_site.png'''
    #    imageName = 'MG0414_d1_data.ms.IGS_TEC_at_site.png'
    #    self.assertTrue(compareImages('MG0414_d1_data.ms.IGS_TEC_at_site.png'))

