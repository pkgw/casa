##########################################################################
# test_req_task_listhistory.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_listhistory/about
#
#
##########################################################################
CASA6 = False
try:
    print("Importing CASAtools")
    import casatools
    from casatasks import listhistory, casalog
    CASA6 = True
except ImportError:
    print ("Cannot import CASAtools using taskinit")
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil

logpath = casalog.logfile()

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/Itziar.ms')
    fakepath = casatools.ctsys.resolve('visibilities/')
    #filepath = casatools.ctsys.resolve('testlog.log')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        dataroot = os.environ.get('CASAPATH').split()[0] + '/'
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/Itziar.ms'
        fakepath = dataroot + 'data/casa-data-req/visibilities/'
    else:
        dataroot = os.environ.get('CASAPATH').split()[0] + '/'
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/Itziar.ms'
        fakepath = dataroot + 'casa-data-req/visibilities/'
    #filepath = 'testlog.log'
    
class listhistory_test(unittest.TestCase):
    
    def setUp(self):
        if not CASA6:
            default(listhistory)
        else:
            pass
    
    def tearDown(self):
        casalog.setlogfile(logpath)
        if os.path.exists('testlog,log'):
            os.remove('testlog.log')
    
    def test_takesMS(self):
        '''test takesMS: Check that list history takes a valid MS and refuses incorrect inputs'''
        casalog.setlogfile('testlog.log')
        listhistory(datapath)
        self.assertFalse('SEVERE' in open('testlog.log').read())
        listhistory(fakepath)
        self.assertTrue('SEVERE' in open('testlog.log').read())
        if CASA6:
            with self.assertRaises(AssertionError):
                listhistory('fake')
        else:
            listhistory('fake')
            self.assertTrue('Argument vis failed to verify' in open('testlog.log').read())
        
    
    def test_logfile(self):
        '''test logfile: Checks to see that a log file is written and populated'''
        casalog.setlogfile('testlog.log')
        listhistory(datapath)
        self.assertTrue('History table entries' in open('testlog.log').read())
    
def suite():
    return[listhistory_test]

if __name__ == '__main__':
    unittest.main()
