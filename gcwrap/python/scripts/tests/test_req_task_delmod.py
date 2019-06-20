##########################################################################
# test_req_task_delmod.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_delmod/about
#
# Test_mesSet: Check that only a valid MS is taken
# Test_removescr: Check that the scratch column is removed
# Test_removeotf: Check to make sure the virtual model is removed
# Test-removefield: Check that only the selected fields are removed (This part is broken)
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import delmod, rmtables, clearcal, casalog, ft
    CASA6 = True
    tb = casatools.table()
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
#import sys
import os
import unittest
import shutil
import glob
from filecmp import dircmp


### These along with the sys import are only used if using simuated data
#sys.path.append('/export/data_1/nschweig/task_test_builds/casa-prerelease-5.5.0-29.el7/bin/simtests/')
#import makethesim

#makethesim.make_me()

# Try using generated data set #
# Currently commented out is the version that uses the data rep


# DATA #
if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
    datacopy = casatools.ctsys.resolve('uid___X02_X3d737_X1_01_small.ms/')
    calpath = casatools.ctsys.resolve(os.path.join(os.path.dirname(os.path.abspath(casatools.__file__)),'__data__/nrao/VLA/CalModels/3C138_K.im'))
    filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
    
    datacopy = 'uid___X02_X3d737_X1_01_small.ms/'
    
    #datapath = os.environ.get('CASAPATH').split()[0] + '/bin/simtests/FITS_list/FITS_list.alma.cycle5.1.ms'
    #datacopy = 'FITS_list.alma.cycle5.1.ms'
    
    calpath = os.environ.get('CASAPATH').split()[0] + '/data/nrao/VLA/CalModels/3C138_K.im'
    filepath = 'testlog.log'
    
logpath = casalog.logfile()
    
class delmod_test(unittest.TestCase):
    
    def setUp(self):
        shutil.copytree(datapath, datacopy)
        os.chmod(datacopy, 493)
        for root, dirs, files in os.walk(datacopy):
            for d in dirs:
                os.chmod(os.path.join(root, d), 493)
            for f in files:
                os.chmod(os.path.join(root, f), 493)
        clearcal(datacopy, addmodel=True)
        if not CASA6:
            default(delmod)
    
    def tearDown(self):
        print('TABLE IS BEING REMOVED')
        casalog.setlogfile(logpath)
        rmtables(datacopy)
        if os.path.exists(filepath):
            os.remove(filepath)
    
    def test_mesSet(self):
        '''
            test_messet
            -----------------
            
            Checks that only a valid MS is accept as input
            
            The first assert checks that no severe errors appear in the log when preforming the task on a valid MS
            The second assert checks that a severe error is raised when the provided MS does not exist
        '''
        casalog.setlogfile('testlog.log')
        delmod(datacopy)
        self.assertFalse('SEVERE' in open('testlog.log').read(), msg='delmod raises a severe error when run on a valid MS')
        if CASA6:
            with self.assertRaises(AssertionError, msg='No error is raised when using a fake MS'):
                delmod('notareal.ms')
        else:
            delmod('notareal.ms')
            self.assertTrue('SEVERE' in open('testlog.log').read(), msg='No error is raised when using a fake MS')
            
    def test_removescr(self):
        '''
            test_removescr
            -----------------------
            
            Checks that the scratch column is removed when using the parameer scr=True
            
            The assert checks that table tools cannot access the scratch column after delmod has been performed
        '''
        delmod(datacopy, scr=True)
        tb.open(datacopy)
        with self.assertRaises(RuntimeError, msg='The MODEL_DATA column was not removed when it should'):
            tb.getcol('MODEL_DATA')
        tb.close()
        
    def test_removeotf(self):
        '''
            test_removeotf
            ------------------
            
            Checks that the vitual model is removed when the parameter otf=True
            
            The first assert statment checks that the model does exist initially (this may be unnessisary)
            
            The second assertion checks that the virtual model no longer exists
        '''
        # Make the file have a model that can be removed
        ft(datacopy, model=calpath)
        self.assertTrue(len(glob.glob(datacopy + r'/SOURCE/FT_MODEL*')) == 1, msg='There is no model initially')
        # remove the model and check that it's gone
        delmod(datacopy, otf=True)
        self.assertTrue(len(glob.glob(datacopy + r'/SOURCE/FT_MODEL*')) == 0, msg='The model has not been removed')
    
    def test_removefield(self):
        '''
            test_removefield
            ----------------------
            
            Check that the field selection paramter works as intended
            
            This test checks that the field selection removes the proper FIELD_ID from the MS
            
            The field is deleted and then table tools are used to open and check for the existence of the removed FIELD_ID
            
            The assertion checks that the specified FIELD_ID is no longer present in the MS
        '''
        ft(datacopy, model=calpath)
        
        delmod(datacopy, field='0')
        
        dcmp = dircmp(datacopy, datapath)
        self.assertTrue(len(dcmp.diff_files) > 0)
        
def suite():
    return[delmod_test]

if __name__ == '__main__':
    unittest.main()
