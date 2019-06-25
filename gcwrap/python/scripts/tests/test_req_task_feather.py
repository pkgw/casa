##########################################################################
# test_req_task_feather.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_feather/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import feather, casalog
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np
from filecmp import dircmp

### DATA ###

if CASA6:
    interpath = casatools.ctsys.resolve('image/orion_tfeather.im/')
    sdpath = casatools.ctsys.resolve('image/orion_tsdmem.image/')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        interpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/orion_tfeather.im/'
        sdpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/orion_tsdmem.image/'
        
    else:
        interpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/orion_tfeather.im/'
        sdpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/orion_tsdmem.image/'

output = 'feathered.im'
output2 = 'feathered2.im'

logpath = casalog.logfile()
logname = 'testlog.log'

class feather_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        if not CASA6:
            default(feather)
            
    def tearDown(self):
        if os.path.exists(output):
            shutil.rmtree(output)
            
        if os.path.exists(output2):
            shutil.rmtree(output2)
            
        if os.path.exists(logname):
            os.remove(logname)
            
        casalog.setlogfile(logpath)
    
    @classmethod
    def tearDownClass(cls):
        pass
    
    def test_combine(self):
        '''
            test_combine
            --------------
            
            Check that interferometric and Single dish images can be combined
        '''
        
        feather(imagename=output, highres=interpath, lowres=sdpath)
        self.assertTrue(os.path.exists(output))
        
    def test_imagename(self):
        '''
            test_imagename
            ----------------
            
            Check that the imagename parameter gives the name of the output image file
        '''
        
        feather(imagename=output, highres=interpath, lowres=sdpath)
        feather(imagename=output2, highres=interpath, lowres=sdpath)
        
        self.assertTrue(os.path.exists(output))
        self.assertTrue(os.path.exists(output2))
        
    def test_highres(self):
        '''
            test_highres
            --------------
            
            Check that the interferometric image is provided with this parameter
            This parameter is nessisary to run the task
        '''
        
        if CASA6:
            with self.assertRaises(AssertionError):
                feather(imagename=output, lowres=sdpath)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output, lowres=sdpath)
            self.assertTrue(('SEVERE' in open(logname).read()))
            
        
    def test_lowres(self):
        '''
            test_lowres
            -------------
            
            Check that the single dish image is provided with this parameter
            This parameter is nessisary to run the task
        '''
        
        if CASA6:
            with self.assertRaises(AssertionError):
                feather(imagename=output, highres=interpath)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output, highres=interpath)
            self.assertTrue('SEVERE' in open(logname).read())
        
        
    def test_sdfactor(self):
        '''
            test_sdfactor
            ---------------
            
            Check that differing sdfactors results in differing image files
        '''
        
        feather(imagename=output, highres=interpath, lowres=sdpath)
        feather(imagename=output2, highres=interpath, lowres=sdpath, sdfactor=0.5)
        
        dcmp = dircmp(output, output2)
        self.assertTrue(len(dcmp.diff_files) > 0)
        
    def test_effdishdiam(self):
        '''
            test_effdishdiam
            ------------------
            
            Check that chaging the effective dish diameter results in differing image files
        '''
        
        feather(imagename=output, highres=interpath, lowres=sdpath)
        feather(imagename=output2, highres=interpath, lowres=sdpath, effdishdiam=1)
        
        dcmp = dircmp(output, output2)
        self.assertTrue(len(dcmp.diff_files) > 0)
        
        if CASA6:
            with self.assertRaises(RuntimeError):
                feather(imagename=output2, highres=interpath, lowres=sdpath, effdishdiam=1000)
        else:
            casalog.setlogfile(logname)
            feather(imagename=output2, highres=interpath, lowres=sdpath, effdishdiam=1000)
            self.assertTrue('SEVERE' in open(logname).read())
        
    def test_lowpassfiltersd(self):
        '''
            test_lowpassfiltersd
            ----------------------
            
            Check that lowpassfiltersd = True results in a different image than the default
        '''
        
        feather(imagename=output, highres=interpath, lowres=sdpath)
        feather(imagename=output2, highres=interpath, lowres=sdpath, lowpassfiltersd=True)
        
        dcmp = dircmp(output, output2)
        self.assertTrue(len(dcmp.diff_files) > 0)
    
    
def suite():
    return[feather_test]

if __name__ == '__main__':
    unittest.main()
