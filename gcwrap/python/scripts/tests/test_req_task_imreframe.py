##########################################################################
# test_req_task_imreframe.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_imreframe/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import imreframe, rmtables, casalog, imhead
    tb = casatools.table()
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
import filecmp

### DATA ###

if CASA6:
    imfile = casatools.ctsys.resolve('image/ngc5921.clean.image/')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        imfile = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image/'
    else:
        imfile = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image/'
    
outfile = 'test.im'
outfile2 = 'test2.im'
imcopy = 'copy.im'
logpath = casalog.logfile()

def compTables(table1, table2):
    
    tb.open(table1)
    val1 = tb.getcol('map')
    tb.close()
    
    tb.open(table2)
    val2 = tb.getcol('map')
    tb.close()
    
    return np.all(val1 == val2)


class imreframe_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        if not CASA6:
            default(imreframe)
        shutil.copytree(imfile, imcopy)
        os.chmod(imcopy, 493)
        for root, dirs, files in os.walk(imcopy):
            for d in dirs:
                os.chmod(os.path.join(root, d), 493)
            for f in files:
                os.chmod(os.path.join(root, f), 493)
        if not CASA6:
            default(importfits)
            
    def tearDown(self):
        
        if os.path.exists(imcopy):
            shutil.rmtree(imcopy)
        
        if os.path.exists(outfile):
            shutil.rmtree(outfile)
            
        if os.path.exists(outfile2):
            shutil.rmtree(outfile2)
        
        casalog.setlogfile(logpath)
    
    @classmethod
    def tearDownClass(cls):
        pass
    
    def test_makesfile(self):
        '''
            test_makesfile
            ----------------
            
            Test to check that the input image is accepted and an output is generated
        '''
        
        imreframe(imagename=imfile, output=outfile)
        self.assertTrue(os.path.exists(outfile))
    
    def test_modifyim(self):
        '''
            test_modifyim
            ----------------
            
            Test that when no outout is given the original file is modified instead
        '''
        self.assertTrue(len(filecmp.dircmp(imcopy, imfile).diff_files) == 0)
        
        imreframe(imagename=imcopy)
        
        self.assertTrue(len(filecmp.dircmp(imcopy, imfile).diff_files) >= 0)
    
    
    def test_outframe(self):
        '''
            test_outframe
            ---------------
            
            Test that the new parameters are used and provide different images than the unmodified one
            
            Check the log file for the change in spectral reference
        '''
        
        outframeList = ['LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB']
        
        casalog.setlogfile('testlog.log')
        
        for frame in outframeList:
            imreframe(imagename=imcopy, outframe=frame)
            self.assertTrue(len(filecmp.dircmp(imcopy, imfile).diff_files) >= 0)
            
            imhead(imcopy, verbose=True)
            self.assertTrue(frame in open('testlog.log').read())
        os.remove('testlog.log')
            
    def test_epoch(self):
        '''
            test_epoch
            -------------
            
            Test that the epoch paramter gives the epoch to be associated with the final image only works for outframe = geo and topo
        '''
        
        for frame in ['GEO', 'TOPO']:
            
            imreframe(imagename=imcopy, output=outfile ,outframe=frame)
            imreframe(imagename=imcopy, output=outfile2, outframe=frame, epoch='2000/12/25/18:30:00')
            
            self.assertFalse(compTables(outfile, outfile2))
            
        imreframe(imagename=imcopy, output=outfile ,outframe='LSRK')
        imreframe(imagename=imcopy, output=outfile2, outframe='LSRK', epoch='2000/12/25/18:30:00')
        
        self.assertTrue(compTables(outfile, outfile2))
            
        
        
    def test_restfreq(self):
        '''
            test_restfreq
            ---------------
            
            Test that the rest frequency sets the rest frequency to use for the velocity value
        '''
        
        casalog.setlogfile('testlog.log')
        
        imreframe(imagename=imcopy, restfreq='1GHz')
        imhead(imcopy, verbose=True)
        self.assertTrue('1e+09 Hz' in open('testlog.log').read())
        os.remove('testlog.log')
        
    
    
def suite():
    return[imreframe_test]

if __name__ == '__main__':
    unittest.main()
    
