##########################################################################
# test_req_task_clearcal.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_clearcal/about
#
# test_takesMS: Checks that clearcal only takes a valid MS
# test_modeldata: Checks that the MODEL_DATA column is generated with the proper values
# test_corrdata: Checks that the CORRECTED_DATA column is the same as the DATA columns
# test_fieldselect: Checks that the field parameter makes proper selections
# test_spwselect: Checks that the spw param makes the proper selections
# test_selectintent: Checks that the intent parameter make the proper selections
# test_addmodel: Check that a MODEL_DATA column is added if the parameter addmodel is True
# test_addcorr: Check that a CORRECTED_DATA columns is added if there was non before
#
##########################################################################
CASA6 = False
try:
    import casatools
    from casatasks import clearcal, casalog, rmtables
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import os
import numpy as np
import unittest
# Try this instead of os.system
import shutil

# DATA #
if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/nep2-shrunk.ms/')
    workingdir = casatools.ctsys.resolve('nep2-shrunk.ms')
    filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/nep2-shrunk.ms'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/nep2-shrunk.ms'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/listobs/nep2-shrunk.ms'
    workingdir = 'nep2-shrunk.ms'
    filepath = 'testlog.log'
    
clearMS = 'nep2-shrunk.ms'

logpath = casalog.logfile()

class clearcal_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.mkdir('fake.ms')
        if not CASA6:
            default(clearcal)
    
    def setUp(self):
        shutil.copytree(datapath, clearMS)
        os.chmod(clearMS, 493)
        for root, dirs, files in os.walk(clearMS):
            for d in dirs:
                os.chmod(os.path.join(root, d), 493)
            for f in files:
                os.chmod(os.path.join(root, f), 493)
    
    def tearDown(self):
        casalog.setlogfile(logpath)
        rmtables('nep2-shrunk.ms')
        if os.path.exists(filepath):
            os.remove(filepath)
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree('fake.ms')
    
    def test_takesMS(self):
        '''test takeMS: Checks that a MS is accepeted and invalid inputs are refused by clearcal'''
        casalog.setlogfile('testlog.log')
        clearcal(clearMS)
        print((os.path.exists('nep2_shrunk.ms')))
        # Need to check logs for all of these bc it will always return NoneType if it passes or fails
        self.assertFalse('SEVERE' in open('testlog.log').read())
        # In CASA 6 assertion Errors will be raised with improper inputs
        if CASA6:
            with self.assertRaises(AssertionError, msg='An int was accepted as an input'):
                clearcal(1)
            with self.assertRaises(AssertionError, msg='A list was accepted as an input'):
                clearcal([])
            with self.assertRaises(AssertionError, msg='A non-existing ms was accepeted'):
                clearcal('foo.ms')
            with self.assertRaises(RuntimeError, msg='A fake ms was accepeted'):
                clearcal('fake.ms')
        else:
            # Check the log for SEVERE failures, assertFalse cannot be used here.
            clearcal(1)
            self.assertTrue('Argument vis failed to verify' in open('testlog.log').read(), msg='An int was accepeted as an input')
            os.system('rm -rf testlog.log')
            casalog.setlogfile('testlog.log')
            clearcal([])
            self.assertTrue('Argument vis failed to verify' in open('testlog.log').read(), msg='A list was accepted as an input')
            os.system('rm -rf testlog.log')
            casalog.setlogfile('testlog.log')
            clearcal('foo.ms')
            self.assertTrue('Argument vis failed to verify' in open('testlog.log').read(), msg='A list was accepted as an input')
            clearcal('fake.ms')
            self.assertTrue('table.dat does not exist' in open('testlog.log').read(), msg='A list was accepted as an input')
       
    def test_modeldata(self):
        '''test modeldata: Checks that the DATA_MODEL column is set to unity in total intensity and zero in polarization'''
        clearcal(clearMS, addmodel=True)
        # Open table with model column generated
        tb.open(clearMS)
        modelCol = tb.getcol('MODEL_DATA')
        # Create random index to check if the value is what we expect
        randrow1 = np.random.randint(0,2)
        randcol = np.random.randint(0,9)
        randrow2 = np.random.randint(0,(modelCol.shape[2]-1))
        # Select data from random index and check values
        real = modelCol[randrow1,randcol,randrow2].real
        imag = modelCol[randrow1,randcol,randrow2].imag
        self.assertTrue(real == 1 and imag == 0, msg = 'MODEL_DATA values not properly generated')
        # Close the table to remove from cache
        tb.close()
        
    def test_corrdata(self):
        '''test Corr data: Check that the CORRECTED_DATA column is the same as the standard DATA column'''
        clearcal(clearMS)
        # Open the table
        tb.open(clearMS)
        # Check that the CORRECTED_DATA and DATA columns are the same
        cor = tb.getcol('CORRECTED_DATA')
        data = tb.getcol('DATA')
        self.assertTrue(np.all(cor == data), msg='CORRECTED_DATA column not the same as the DATA column')
        # Close the table to remove from cache
        tb.close()
    
    def test_fieldselect(self):
        '''test field select: Check that field select only changes a portion of the data'''
        clearcal(clearMS)
        # Open table and allow modifications
        tb.open(clearMS, nomodify=False)
        corrCol = tb.getcol('CORRECTED_DATA')
        rownum = corrCol.shape[2]
        # Make array of new data. All zeros
        newdata = np.array([[[0 for i in range(rownum)] for k in range(9)] for j in range(2)])
        # Put this new data into the CORRECTED_DATA column
        tb.putcol('CORRECTED_DATA',newdata)
        zeroColOld = tb.getcol('CORRECTED_DATA')[0,0,0]
        endColOld = tb.getcol('CORRECTED_DATA')[0,0,(rownum-1)]
        clearcal(clearMS, field='0')
        zeroColNew = tb.getcol('CORRECTED_DATA')[0,0,0]
        endColNew = tb.getcol('CORRECTED_DATA')[0,0,(rownum-1)]
        # Compare the old col to new in the field that should be changed, then the field that shouldn't
        self.assertFalse(zeroColNew.real == zeroColOld.real, msg='values were not corrected')
        self.assertTrue(endColNew.real == endColOld.real, msg='values were corrected when they should not have been')
        tb.close()
        
        if CASA6:
            with self.assertRaises(AssertionError, msg='An int was accepted as input'):
                clearcal(clearMS, field=2)
        else:
            casalog.setlogfile('testlog.log')
            clearcal(clearMS, field=2)
            self.assertTrue('SEVERE' in open('testlog.log').read(), msg='An int was accepted as input')
            
    def test_spwselect(self):
        '''test spw select: Check that spw select only changes correct spw'''
        clearcal(clearMS)
        # open the table and get number of rows
        tb.open(clearMS, nomodify=False)
        corrCol = tb.getcol('CORRECTED_DATA')
        rownum = corrCol.shape[2]
        # create an empty column and replace CORRECTED_DATA with it
        newdata = np.array([[[0 for i in range(rownum)] for k in range(9)] for j in range(2)])
        tb.putcol('CORRECTED_DATA', newdata)
        # Check that spw select only corrects the selected spw
        colOld = tb.getcol('CORRECTED_DATA')[0,0,0]
        clearcal(clearMS, spw='100')
        colNew1 = tb.getcol('CORRECTED_DATA')[0,0,0]
        self.assertTrue(colNew1.real == colOld.real, msg='Data was corrected even when not selected')
        clearcal(clearMS, spw='0')
        colNew2 = tb.getcol('CORRECTED_DATA')[0,0,0]
        self.assertFalse(colNew2.real == colOld.real, msg='Data was not corrected when selected')
        # Close the table
        tb.close()
        # check invalid inputs
        if CASA6:
            with self.assertRaises(AssertionError, msg='An int was accepted as input'):
                clearcal(clearMS, spw=2)
        else:
            casalog.setlogfile('testlog.log')
            clearcal(clearMS, spw=2)
            self.assertTrue('SEVERE' in open('testlog.log').read(), msg='An int was accepted as input')
            
    def test_selectintent(self):
        '''test select intent: Check that the intent param correcty preforms for only the provided intent'''
        clearcal(clearMS)
        # Open the table and get number of rows
        tb.open(clearMS, nomodify=False)
        corrCol = tb.getcol('CORRECTED_DATA')
        rownum = corrCol.shape[2]
        # Create an empty table and replace CORRECTED_DATA with it
        newdata = np.array([[[0 for i in range(rownum)] for k in range(9)] for j in range(2)])
        tb.putcol('CORRECTED_DATA', newdata)
        # Check that clearcal only corrects the selected portion of the data
        colOld = tb.getcol('CORRECTED_DATA')[0,0,0]
        clearcal(clearMS, intent='*BANDPASS*')
        newCol1 = tb.getcol('CORRECTED_DATA')[0,0,0]
        newCol2 = tb.getcol('CORRECTED_DATA')[0,0,(rownum-1)]
        self.assertFalse(newCol1.real == colOld.real, msg='Selected data was not changed')
        self.assertTrue(newCol2.real == colOld.real, msg='Unselectd data was changed')
        # Close the table
        tb.close()
        # Check invalid inputs
        if CASA6:
            with self.assertRaises(AssertionError, msg='An int was accepted as input'):
                clearcal(clearMS, intent=2)
        else:
            casalog.setlogfile('testlog.log')
            clearcal(clearMS, intent=2)
            self.assertTrue('SEVERE' in open('testlog.log').read(), msg='An int was selected as input')
            
    def test_addmodel(self):
        '''test add model: Check that the addmodel parameter must be set to True for the MODEL_DATA column to be added'''
        tb.open(clearMS)
        ### Is this needed, or is this a test of table tools?
        #with self.assertRaises(RuntimeError):
            #tb.getcol('MODEL_DATA')
        # Set addmodel to True and check if a column is generated
        clearcal(clearMS, addmodel=True)
        self.assertTrue(type(tb.getcol('MODEL_DATA')) == np.ndarray, msg='The MODEL_DATA column was not created properly')
        # Close the table
        tb.close()
        # Test invalid inputs
        if CASA6:
            with self.assertRaises(AssertionError, msg='An int was accepted as input'):
                clearcal(clearMS, addmodel=2)
        else:
            casalog.setlogfile('testlog.log')
            clearcal(clearMS, addmodel=2)
            self.assertTrue('SEVERE' in open('testlog.log').read(), msg='An ine was accepted as input')
            
    def test_addcorr(self):
        '''test add corr: Check that a CORRECTED_DATA column is added if none existed before'''
        tb.open(clearMS)
        ### This might be just testing tabletools again
        #with self.assertRaises(RuntimeError):
            #tb.getcol('CORRECTED_DATA')
        # Check that a CORRECTED_DATA column was generated
        clearcal(clearMS)
        self.assertTrue(type(tb.getcol('CORRECTED_DATA')) == np.ndarray, msg='The CORRECTED_DATA column was not created properly')
        # Close the table
        tb.close()
            
def suite():
    return[clearcal_test]

if __name__ == '__main__':
    unittest.main()
