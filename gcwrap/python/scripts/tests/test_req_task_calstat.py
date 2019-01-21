##########################################################################
# test_req_task_calstat.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_calstat/about
#
# Test_logreturn checks to make sure a logfile is generated and populated
# Test_dictreturn checks that the result is a python dict object containing keys specified in the documentation
# Test_takescal checks that a caltable is accepted and non-cal tables are rejected
# Test_axis checks that different axis vaules will provide different information
# Test_axisvals checks that the values for axis provided in the documentatin are accepted as valid values
# Test_datacolumn checks that different datacolumn values provide different information
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import casalog, calstat
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil

if CASA6:
    datapath = casatools.ctsys.resolve('caltables/ggtau.1mm.amp.gcal')
    #filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/caltables/ggtau.1mm.amp.gcal'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/caltables/ggtau.1mm.amp.gcal'
    #filepath = os.environ.get('CASAPATH').split()[0] + '/bin/nosedir/testlog.log'
        
logpath = casalog.logfile()
contained = ['max', 'mean', 'medabsdevmed', 'median', 'min', 'npts', 'quartile', 'rms', 'stddev', 'sum', 'sumsq', 'var']

class calstat_test(unittest.TestCase):
     
    def setUp(self):
        if not CASA6:
            default(calstat)
     
    def tearDown(self):
        if os.path.exists('testlog.log'):
            os.remove('testlog.log')
     
    def test_logreturn(self):
        '''logreturn test: Test that a logfile is written and populated with the expected information'''
        casalog.setlogfile('testlog.log')
        # Check that a logfile is populated and contains the correct infromation
        calstat(datapath, axis='TIME')
        for item in contained:
            self.assertTrue( item in open('testlog.log').read(), msg='Fails to write required information to the log')
        
    def test_dictreturn(self):
        '''dictionary test: Test that calstat makes a python dict with the expected keys'''
        # Check that type of returned object is a python dict
        caldict = calstat(datapath)
        self.assertTrue(isinstance(caldict, dict), msg='calstat does not return a python dict')
        # Check that the dict contains the correct values
        self.assertTrue(all (k in caldict['GAIN'] for k in contained), msg='Dictionalry does not contain all the correct values')
            
    def test_takescal(self):
        '''takes cal test: Test that calstat only takes cal tables'''
        self.assertTrue(calstat(datapath), msg='calstat fails to take a caltable')
        ### Needs work. try and produce some fake cal tables to test with
        # No type checking for CASA 6, but some for CASA 5
        if CASA6:
            with self.assertRaises(RuntimeError, msg='Fails to recognize non-caltable'):
                calstat(casatools.ctsys.resolve('visibilities/'))
        else:
            self.assertFalse(calstat(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/'), msg='Fails to recognize non-caltable')
            self.assertFalse(calstat(1), msg='Fails to refuse int type')
            self.assertFalse(calstat([]), msg='Fails to refuse a list')
            
    def test_axis(self):
        '''axis test: Check that the axis param makes different selections'''
        timeAxis = calstat(datapath, axis='TIME')
        ampAxis = calstat(datapath, axis='amp')
        self.assertFalse(timeAxis == ampAxis, msg='different axis selections should give different values')
        
    def test_axisvals(self):
        '''axis val test: Test that the inputs listed in the documentatin function properly'''
        # Some type checking for CASA 5
        if not CASA6:
            self.assertFalse(calstat(datapath, axis='abc'), msg='Fails to recognize non-existing axis')
            self.assertFalse(calstat(datapath, axis=1), msg='Takes int as input when it should not')
            self.assertFalse(calstat(datapath, axis=[]), msg='Takes list as input when it should not')
        # Check that the documented inputs function
        for item in ['amp', 'amplitude', 'phase', 'real', 'imag', 'imaginary']:
            self.assertTrue(calstat(datapath, axis=item), msg='axis {} is not recognized'.format(item))
    
    def test_dataColumn(self):
        '''datacolumn test: Check that the datacolumn param makes unique selections'''
        # Some type checking for CASA 5
        if not CASA6:
            self.assertFalse(calstat(datapath, axis='amp', datacolumn='abc'), msg='Fails to recognize non-existing datacolumn')
            self.assertFalse(calstat(datapath, axis='amp', datacolumn=1), msg='Takes an int as input when it should not')
            self.assertFalse(calstat(datapath, axis='amp', datacolumn=[]), msg='Takes list as input when it should not')
        # Check that datacolumn selects different data when it should
        timeCol = calstat(datapath, axis='amp', datacolumn='TIME')
        gainCol = calstat(datapath, axis='amp', datacolumn='GAIN')
        self.assertFalse(timeCol == gainCol, msg='different datacolumns should give different results')
        # Added for additional test coverage. These should correspond to CORRECTED_DATA and MODEL_DATA
        # These columns are not present for the current cal table
        #self.assertTrue(calstat(datapath, axis='amp', datacolumn='corrected'))
        #self.assertTrue(calstat(datapath, axis='amp', datacolumn='model'))
        

def suite():
    return[calstat_test]

if __name__ == '__main__':
    unittest.main()
