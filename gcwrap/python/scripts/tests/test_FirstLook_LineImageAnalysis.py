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

# Needs to Be Updated
"""
def compareTables(dataset):
    refmeaset = indir + "pybotWorkspace/prime_FirstLookatImageAnalysis/" + str(dataset)
    measet = indir + "pybotWorkspace/FirstLookatImageAnalysis/" + str(dataset)
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
    return [Test010_FirstLookatImageAnalysis,Test020_FirstLookatImageAnalysis,Test021_FirstLookatImageAnalysis]

class Test010_FirstLookatImageAnalysis(unittest.TestCase):
    def setUp(self):

        if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
            casaguidedata_path = "/data/casaguidedata/"
        else:
            casaguidedata_path = "/casaguidedata/"


        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_cont.image",os.getcwd()+'/sis14_twhya_cont.image')
        os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_n2hp.image",os.getcwd()+'/sis14_twhya_n2hp.image')
        if os.uname()[0] == 'Darwin':
            os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Image_Analysis'")
        else:
            os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Image_Analysis'")

        time.sleep(5) # Allow extract time to download script

        lines = open('FirstLookatImageAnalysis.py')
        file = open("newfile.txt", "w")
        for line in lines:

            if "rm -rf sis14_twhya_cont.image" in line:
                continue
            if "rm -rf sis14_twhya_n2hp.image" in line:
                continue

            pattern = r'''niter\ *=\ *(5000)'''
            if re.search(pattern,line):
                line = re.sub( pattern, 'niter=250', line )

            file.write(line)
        file.close()
        os.remove('FirstLookatImageAnalysis.py')
        os.rename("newfile.txt",'FirstLookatImageAnalysis.py')

        time.sleep(5)
    def tearDown(self):
        pass
    def test_00_runGuide(self):
        '''Run Casa Guide: First Look at Line Imaging'''


        exec(compile(open('FirstLookatImageAnalysis.py', "rb").read(), 'FirstLookatImageAnalysis.py', 'exec'))
                
        return True

class Test020_FirstLookatImageAnalysis(unittest.TestCase):

    def setUp(self):
        pass
    def tearDown(self):
        pass


    def test_7_sis14_twhya_n2hp_mom0(self):
        '''Test 7: Check sis14_twhya_n2hp.mom0'''
        tableName = 'sis14_twhya_n2hp.mom0'
        self.assertTrue(openTable(tableName))

    def test_8_sis14_twhya_n2hp_mom1(self):
        '''Test 8: Check sis14_twhya_n2hp.mom1'''
        tableName = 'sis14_twhya_n2hp.mom1'
        self.assertTrue(openTable(tableName))


    def test_12_twhya_cont_fits(self):
        '''Test 12: Check twhya_cont.fits'''
        self.assertTrue(assert_file('twhya_cont.fits'))

    def test_13_twhya_n2hp_fits(self):
        '''Test 13: Check twhya_n2hp.fits'''
        self.assertTrue(assert_file('twhya_n2hp.fits'))

####################################################################################################
class Test021_FirstLookatImageAnalysis(unittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        os.unlink(os.getcwd()+'/sis14_twhya_cont.image')
        os.unlink(os.getcwd()+'/sis14_twhya_n2hp.image')
        rmtables("sis14_twhya*")
        os.system("rm -rf *.last")
        os.system("rm -rf *.fits")
        os.system("rm -rf *.flagversions")
    def test_1_Statistics_RMS(self):
        '''Statistics RMS'''
        my_stats = imstat("sis14_twhya_n2hp.image", chans="0~4")
        rms_stat = my_stats['rms'][0]
        expected = 0.02 # 20mJy
        print("Expected: %s , Actual: %s"%(expected,rms_stat))
        assert rms_stat >= 0.018 and rms_stat <= 0.022, "Error in RMS Statistics"

    def test_2_Statistics_FLUX(self):
        '''Statistics Flux'''
        my_stats = imstat("sis14_twhya_cont.image", box="100,100,150,150")
        flux_stat = my_stats['flux'][0]
        expected = 1.5 # 1.5Jy
        print("Expected: %s , Actual: %s"%(expected,flux_stat))
        assert flux_stat >= 1.3 and flux_stat <= 1.7, "Error in Flux Statistics"

if __name__ == '__main__':
    unittest.main()
