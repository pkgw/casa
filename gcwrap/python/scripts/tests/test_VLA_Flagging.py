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
    return [Test010_FlaggingVLAData,Test020_FlaggingVLAData,Test021_FlaggingVLAData]


class Test010_FlaggingVLAData(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"

        # Untar data. Dataset needs write permissions so copy a full set to working dir

        import tarfile
        tar = tarfile.open(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "raw/SNR_G55_10s.tar.gz")
        tar.extractall()
        tar.close()

        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Flagging'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/VLA_CASA_Flagging'")

        time.sleep(5) # Allow extract time to download script

        lines = open('VLACASAFlagging.py')
        file = open("newfile.txt", "w")
        for line in lines:
            pattern = r"""display\ *=\ *['"].*['"]"""
            if re.search( pattern, line ):
                line = re.sub( pattern, "display=''", line )
            if "%cpaste" in line: 
                continue
            file.write(line)
        file.close()
        os.remove('VLACASAFlagging.py')
        os.rename("newfile.txt",'VLACASAFlagging.py')

    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide:  Topical Guide Flagging VLA Data'''

        exec(compile(open('VLACASAFlagging.py', "rb").read(), 'VLACASAFlagging.py', 'exec'))
                
        return True

class Test020_FlaggingVLAData(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass

    def test_1_SNR_G55_10s_hanning_initBP(self):
        '''Test 1: Check SNR_G55_10s-hanning.initBP'''
        tableName = 'SNR_G55_10s-hanning.initBP'
        self.assertTrue(openTable(tableName))

    def test_2_SNR_G55_10s_hanning_initPh(self):
        '''Test 2: Check SNR_G55_10s-hanning.initPh'''
        tableName = 'SNR_G55_10s-hanning.initPh'
        self.assertTrue(openTable(tableName))

    def test_3_SNR_G55_10s_hanning_ms(self):
        '''Test 3.1: Check SNR_G55_10s-hanning.ms'''
        tableName = 'SNR_G55_10s-hanning.ms'
        self.assertTrue(openTable(tableName))

    #TODO: Check if <Dataset.ms>.flagversions only gets created when interactive

    """
    def test_SNR_G55_10s_hanning_ms_flagversions(self):
        '''Test 3.2: Check SNR_G55_10s-hanning.ms.flagversions'''
        tableName = 'SNR_G55_10s-hanning.ms.flagversions'
        self.assertTrue(openTable(tableName))

    def test_4_SNR_G55_10s_ms_flagversions(self):
        '''Test 4: Check SNR_G55_10s.ms.flagversions'''
        tableName = 'SNR_G55_10s.ms.flagversions'
        self.assertTrue(openTable(tableName))
    """
    #TODO: Commented Out while plotting is turned off
    """
    def test_5_SNR_G55_10s_plotants_png(self):
        '''Test 5: Check SNR_G55_10s.plotants.png'''
        self.assertTrue(assert_file('SNR_G55_10s.plotants.png'))

    def test_7_amp_v_freq_afterHanning_png(self):
        '''Test 7: Check amp_v_freq.afterHanning.png'''
        self.assertTrue(assert_file('amp_v_freq.afterHanning.png'))

    def test_9_amp_v_freq_afterRFlag_png(self):
        '''Test 9: Check amp_v_freq.afterRFlag.png'''
        self.assertTrue(assert_file('amp_v_freq.afterRFlag.png'))

    def test_11_amp_v_freq_beforeHanning_png(self):
        '''Test 11: Check amp_v_freq.beforeHanning.png'''
        self.assertTrue(assert_file('amp_v_freq.beforeHanning.png'))

    def test_13_amp_v_freq_beforeRFlag_png(self):
        '''Test 13: Check amp_v_freq.beforeRFlag.png'''
        self.assertTrue(assert_file('amp_v_freq.beforeRFlag.png'))


    def test_15_amp_v_freq_after_tfcrop_Scan190_png(self):
        '''Test 15: Check amp_v_freq_after_tfcrop_Scan190.png'''
        self.assertTrue(assert_file('amp_v_freq_after_tfcrop_Scan190.png'))

    def test_17_amp_v_freq_before_tfcrop_Scan190_png(self):
        '''Test 17: Check amp_v_freq_before_tfcrop_Scan190.png'''
        self.assertTrue(assert_file('amp_v_freq_before_tfcrop_Scan190.png'))

    def test_18_flaggingreason_vs_time_png(self):
        '''Test 18: Check flaggingreason_vs_time.png'''
        self.assertTrue(assert_file('flaggingreason_vs_time.png'))
    """
####################################################################################################
class Test021_FlaggingVLAData(unittest.TestCase):
    def setUp(self):
        pass
    @classmethod
    def tearDownClass(cls):

        rmtables("SNR_G55_10s-hanning*")
        os.system("rm -rf *.last")
        os.system("rm -rf ngc612region.crtf")
        os.system("rm -rf *.flagversions")

    def test_19_Flagging_Summary_G55(self):
        '''Flagging Summary G55.7+3.4'''
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['field']['G55.7+3.4']['flagged'] / flagInfo['field']['G55.7+3.4']['total'])
        print("Actual: %s"%(val))
        assert val >= 40 and val <= 45, "Percent Flagged Error G55.7+3.4"

    def test_20_Flagging_Summary_3C147(self):
        '''Flagging Summary 0542+498=3C147'''
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['field']['0542+498=3C147']['flagged'] / flagInfo['field']['0542+498=3C147']['total'])
        print("Actual: %s"%(val))
        assert val >= 65 and val <= 70, "Percent Flagged Error 0542+498=3C147"

    def test_21_Flagging_Summary_J1925(self):
        '''Flagging Summary J1925+2106'''
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['field']['J1925+2106']['flagged'] / flagInfo['field']['J1925+2106']['total'])
        print("Actual: %s"%(val))
        assert val >= 45 and val <= 50, "Percent Flagged Error J1925+2106"

    def test_22_spw_0_Flagging(self):
        '''Flagging Summary SPW 0'''
        spw = 0
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw'][str(spw)]['total'])
        print("Actual: %s"%(val))
        assert val >= 60 and val <= 65 ,"Spectral Window 0 Flagging Error"

    def test_23_spw_1_Flagging(self):
        '''Flagging Summary SPW 1'''
        spw = 1
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw'][str(spw)]['total'])
        print("Actual: %s"%(val))
        assert val >= 25 and val <= 30 ,"Spectral Window 1 Flagging Error"

    def test_24_spw_2_Flagging(self):
        '''Flagging Summary SPW 2'''
        spw = 2
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw'][str(spw)]['total'])
        print("Actual: %s"%(val))
        assert val >= 62 and val <= 67 ,"Spectral Window 2 Flagging Error"

    def test_25_spw_3_Flagging(self):
        '''Flagging Summary SPW 3'''
        spw = 3
        flagInfo = flagdata(vis='SNR_G55_10s-hanning.ms', mode='summary', action='calculate',  spwchan=True)
        val = 100.0 * (flagInfo['spw'][str(spw)]['flagged'] / flagInfo['spw'][str(spw)]['total'])
        print("Actual: %s"%(val))
        assert val >= 20 and val <= 25 ,"Spectral Window 3 Flagging Error"


####################################################################################################

