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

			file.write(line)
		file.close()
		os.remove('FirstLookatImaging.py')
		os.rename("newfile.txt",'FirstLookatImaging.py')

		time.sleep(15)
	def tearDown(self):

		pass
	
	def test_00_runGuide(self):
		'''Run Casa Guide: First Look at Imaging'''


		exec(compile(open('FirstLookatImaging.py').read(), 'FirstLookatImaging.py', 'exec'))
                
		return True

class Test020_FirstLookatImaging(unittest.TestCase):

	def setUp(self):
		pass
	def tearDown(self):
		pass
	def test_1_primary_robust_flux(self):
		'''Test 1: Check primary_robust.flux'''
		tableName = 'primary_robust.flux'
		self.assertTrue(openTable(tableName))

	def test_2_primary_robust_image(self):
		'''Test 2: Check primary_robust.image'''
		tableName = 'primary_robust.image'
		self.assertTrue(openTable(tableName))

	def test_3_primary_robust_model(self):
		'''Test 3: Check primary_robust.model'''
		tableName = 'primary_robust.model'
		self.assertTrue(openTable(tableName))

	def test_4_primary_robust_psf(self):
		'''Test 4: Check primary_robust.psf'''
		tableName = 'primary_robust.psf'
		self.assertTrue(openTable(tableName))

	def test_5_primary_robust_residual(self):
		'''Test 5: Check primary_robust.residual'''
		tableName = 'primary_robust.residual'
		self.assertTrue(openTable(tableName))

	def test_6_secondary_flux(self):
		'''Test 6: Check secondary.flux'''
		tableName = 'secondary.flux'
		self.assertTrue(openTable(tableName))

	def test_7_secondary_image(self):
		'''Test 7: Check secondary.image'''
		tableName = 'secondary.image'
		self.assertTrue(openTable(tableName))

	def test_8_secondary_model(self):
		'''Test 8: Check secondary.model'''
		tableName = 'secondary.model'
		self.assertTrue(openTable(tableName))

	def test_9_secondary_psf(self):
		'''Test 9: Check secondary.psf'''
		tableName = 'secondary.psf'
		self.assertTrue(openTable(tableName))

	def test_10_secondary_residual(self):
		'''Test 10: Check secondary.residual'''
		tableName = 'secondary.residual'
		self.assertTrue(openTable(tableName))

	def test_11_secondary_bigpix_flux(self):
		'''Test 11: Check secondary_bigpix.flux'''
		tableName = 'secondary_bigpix.flux'
		self.assertTrue(openTable(tableName))

	def test_12_secondary_bigpix_image(self):
		'''Test 12: Check secondary_bigpix.image'''
		tableName = 'secondary_bigpix.image'
		self.assertTrue(openTable(tableName))

	def test_13_secondary_bigpix_model(self):
		'''Test 13: Check secondary_bigpix.model'''
		tableName = 'secondary_bigpix.model'
		self.assertTrue(openTable(tableName))

	def test_14_secondary_bigpix_psf(self):
		'''Test 14: Check secondary_bigpix.psf'''
		tableName = 'secondary_bigpix.psf'
		self.assertTrue(openTable(tableName))

	def test_15_secondary_bigpix_residual(self):
		'''Test 15: Check secondary_bigpix.residual'''
		tableName = 'secondary_bigpix.residual'
		self.assertTrue(openTable(tableName))

	def test_16_secondary_robust_flux(self):
		'''Test 16: Check secondary_robust.flux'''
		tableName = 'secondary_robust.flux'
		self.assertTrue(openTable(tableName))

	def test_17_secondary_robust_image(self):
		'''Test 17: Check secondary_robust.image'''
		tableName = 'secondary_robust.image'
		self.assertTrue(openTable(tableName))

	def test_18_secondary_robust_model(self):
		'''Test 18: Check secondary_robust.model'''
		tableName = 'secondary_robust.model'
		self.assertTrue(openTable(tableName))

	def test_19_secondary_robust_psf(self):
		'''Test 19: Check secondary_robust.psf'''
		tableName = 'secondary_robust.psf'
		self.assertTrue(openTable(tableName))

	def test_20_secondary_robust_residual(self):
		'''Test 20: Check secondary_robust.residual'''
		tableName = 'secondary_robust.residual'
		self.assertTrue(openTable(tableName))

	def test_21_secondary_uncalibrated_flux(self):
		'''Test 21: Check secondary_uncalibrated.flux'''
		tableName = 'secondary_uncalibrated.flux'
		self.assertTrue(openTable(tableName))

	def test_22_secondary_uncalibrated_image(self):
		'''Test 22: Check secondary_uncalibrated.image'''
		tableName = 'secondary_uncalibrated.image'
		self.assertTrue(openTable(tableName))

	def test_23_secondary_uncalibrated_model(self):
		'''Test 23: Check secondary_uncalibrated.model'''
		tableName = 'secondary_uncalibrated.model'
		self.assertTrue(openTable(tableName))

	def test_24_secondary_uncalibrated_psf(self):
		'''Test 24: Check secondary_uncalibrated.psf'''
		tableName = 'secondary_uncalibrated.psf'
		self.assertTrue(openTable(tableName))

	def test_25_secondary_uncalibrated_residual(self):
		'''Test 25: Check secondary_uncalibrated.residual'''
		tableName = 'secondary_uncalibrated.residual'
		self.assertTrue(openTable(tableName))

	def test_26_secondary_unflagged_flux(self):
		'''Test 26: Check secondary_unflagged.flux'''
		tableName = 'secondary_unflagged.flux'
		self.assertTrue(openTable(tableName))

	def test_27_secondary_unflagged_image(self):
		'''Test 27: Check secondary_unflagged.image'''
		tableName = 'secondary_unflagged.image'
		self.assertTrue(openTable(tableName))

	def test_28_secondary_unflagged_model(self):
		'''Test 28: Check secondary_unflagged.model'''
		tableName = 'secondary_unflagged.model'
		self.assertTrue(openTable(tableName))

	def test_29_secondary_unflagged_psf(self):
		'''Test 29: Check secondary_unflagged.psf'''
		tableName = 'secondary_unflagged.psf'
		self.assertTrue(openTable(tableName))

	def test_30_secondary_unflagged_residual(self):
		'''Test 30: Check secondary_unflagged.residual'''
		tableName = 'secondary_unflagged.residual'
		self.assertTrue(openTable(tableName))

	def test_40_twhya_cont_flux(self):
		'''Test 40: Check twhya_cont.flux'''
		tableName = 'twhya_cont.flux'
		self.assertTrue(openTable(tableName))

	def test_41_twhya_cont_image(self):
		'''Test 41: Check twhya_cont.image'''
		tableName = 'twhya_cont.image'
		self.assertTrue(openTable(tableName))

	def test_42_twhya_cont_model(self):
		'''Test 42: Check twhya_cont.model'''
		tableName = 'twhya_cont.model'
		self.assertTrue(openTable(tableName))

	def test_43_twhya_cont_pbcor_image(self):
		'''Test 43: Check twhya_cont.pbcor.image'''
		tableName = 'twhya_cont.pbcor.image'
		self.assertTrue(openTable(tableName))

	def test_44_twhya_cont_psf(self):
		'''Test 44: Check twhya_cont.psf'''
		tableName = 'twhya_cont.psf'
		self.assertTrue(openTable(tableName))

	def test_45_twhya_cont_residual(self):
		'''Test 45: Check twhya_cont.residual'''
		tableName = 'twhya_cont.residual'
		self.assertTrue(openTable(tableName))

	def test_46_twhya_cont_auto_flux(self):
		'''Test 46: Check twhya_cont_auto.flux'''
		tableName = 'twhya_cont_auto.flux'
		self.assertTrue(openTable(tableName))

	def test_47_twhya_cont_auto_image(self):
		'''Test 47: Check twhya_cont_auto.image'''
		tableName = 'twhya_cont_auto.image'
		self.assertTrue(openTable(tableName))

	def test_48_twhya_cont_auto_model(self):
		'''Test 48: Check twhya_cont_auto.model'''
		tableName = 'twhya_cont_auto.model'
		self.assertTrue(openTable(tableName))

	def test_49_twhya_cont_auto_psf(self):
		'''Test 49: Check twhya_cont_auto.psf'''
		tableName = 'twhya_cont_auto.psf'
		self.assertTrue(openTable(tableName))

	def test_50_twhya_cont_auto_residual(self):
		'''Test 50: Check twhya_cont_auto.residual'''
		tableName = 'twhya_cont_auto.residual'
		self.assertTrue(openTable(tableName))

	def test_51_twhya_smoothed_ms(self):
		'''Test 51: Check twhya_smoothed.ms'''
		tableName = 'twhya_smoothed.ms'
		self.assertTrue(openTable(tableName))

