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
	refmeaset = indir + "pybotWorkspace/prime_FirstLookatLineImaging/" + str(dataset)
	measet = indir + "pybotWorkspace/FirstLookatLineImaging/" + str(dataset)
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
                print 'Length of '+ str(referencetab) +' differ from '+ str(testtab)+','+ str(tb.nrows())+ '!=' + str(tb2.nrows())
                retval = False
            else:
                for therow in xrange(tb.nrows()):
            
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
                                    print 'ERROR: Column ' + str(col) + ' of '  + str(referencetab) +  ' and ' + str(testtab)+  ' do not agree within tolerance '+ str(tolerance)
                                    break
                        else:
                            print 'ERROR: Column ' +str(col)+ ' of ' +str(referencetab)+ ' and ' +str(testtab) + ' do not agree.'
                            print 'ERROR: First row to differ is row=' + str(therow)
                            retval = False
                            break
        finally:
            tb.close()
            tb2.close()
    
    else:
        print 'Column: ' +str(col) + 'are not varcolumns.'
        retval = False

    if retval:
	print 'Column ' + str(col) + ' of '  + str(referencetab) +  ' and ' + str(testtab) + ' agree'
        
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
	return [Test010_FirstLookatLineImaging,Test020_FirstLookatLineImaging,Test021_FirstLookatLineImaging]

class Test010_FirstLookatLineImaging(unittest.TestCase):
	def setUp(self):

		if os.path.isdir(os.environ.get('CASAPATH').split()[0] + "/data/casaguidedata"):
			casaguidedata_path = "/data/casaguidedata/"
		else:
			casaguidedata_path = "/casaguidedata/"

		os.symlink(os.environ.get('CASAPATH').split()[0] + casaguidedata_path + "working_data/sis14_twhya_selfcal.ms",os.getcwd()+'/sis14_twhya_selfcal.ms')
		if os.uname()[0] == 'Darwin':
			os.system(os.environ.get('CASAPATH').split()[0] +"/Resources/python/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Line_Imaging'")
		else:
			os.system(os.environ.get('CASAPATH').split()[0] +"/lib/python2.7/extractCASAscript.py -n -p -d 'https://casaguides.nrao.edu/index.php/First_Look_at_Line_Imaging'")
		time.sleep(5) # Allow extract time to download script

		lines = open('FirstLookatLineImaging.py')
		file = open("newfile.txt", "w")
		for line in lines:

			if "rm -rf sis14_twhya_selfcal.ms" in line:
				continue

			if "niter=5000)" in line:
				file.write("niter=250)\n")
				continue
			file.write(line)
		file.close()
		os.remove('FirstLookatLineImaging.py')
		os.rename("newfile.txt",'FirstLookatLineImaging.py')

		time.sleep(15)

	def tearDown(self):
		pass
	def test_00_runGuide(self):
		'''Run Casa Guide: First Look at Line Imaging'''


		execfile('FirstLookatLineImaging.py')
                
		return True




class Test020_FirstLookatLineImaging(unittest.TestCase):

	def setUp(self):
		pass
	def tearDown(self):
		pass


	def test_10_twhya_n2hp_flux(self):
		'''Test 10: Check twhya_n2hp.flux'''
		tableName = 'twhya_n2hp.flux'
		self.assertTrue(openTable(tableName))

	def test_11_twhya_n2hp_image(self):
		'''Test 11: Check twhya_n2hp.image'''
		tableName = 'twhya_n2hp.image'
		self.assertTrue(openTable(tableName))

	def test_12_twhya_n2hp_model(self):
		'''Test 12: Check twhya_n2hp.model'''
		tableName = 'twhya_n2hp.model'
		self.assertTrue(openTable(tableName))

	def test_13_twhya_n2hp_pbcor_image(self):
		'''Test 13: Check twhya_n2hp.pbcor.image'''
		tableName = 'twhya_n2hp.pbcor.image'
		self.assertTrue(openTable(tableName))

	def test_14_twhya_n2hp_psf(self):
		'''Test 14: Check twhya_n2hp.psf'''
		tableName = 'twhya_n2hp.psf'
		self.assertTrue(openTable(tableName))

	def test_15_twhya_n2hp_residual(self):
		'''Test 15: Check twhya_n2hp.residual'''
		tableName = 'twhya_n2hp.residual'
		self.assertTrue(openTable(tableName))

####################################################################################################


class Test021_FirstLookatLineImaging(unittest.TestCase):
	def setUp(self):
		pass
	def tearDown(self):
		pass
