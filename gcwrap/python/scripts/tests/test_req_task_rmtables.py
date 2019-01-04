##########################################################################
# test_req_task_listobs.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.0/global-task-list/task_rmtables/about
#
# test_removeMS: Checks that a MS is removed properly
# test_removeCAL: Checks that a CAL table is removed properly
# test_removeIMG: Checks that an Image is removed properly
# test_cacheMS: Check that opening then removing an MS clears it from cache
# test_cacheCAL: Checks that opening then removing a CAL table clears it from cache
# tesst_cacheIMG: Checks that opening then removing an Image clears it from cache
#
##########################################################################
CASA6 = False
try:
    print("Importing CASAtools")
    import casatools
    from casatasks import rmtables
    tb = casatools.table()
    CASA6 = True
except ImportError:
    print ("Cannot import CASAtools using taskinit")
    from __main__ import default
    from tasks import *
    from taskinit import *
import os
import unittest
import shutil

if CASA6:
    mesSet = casatools.ctsys.resolve('visibilities/vla/ngc7538_ut.ms')
    calTab = casatools.ctsys.resolve('caltables/anpos.manual.cal')
    imfile = casatools.ctsys.resolve('image/ngc5921.clean.image')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/'
    mesSet = datapath + 'visibilities/vla/ngc7538_ut.ms'
    calTab = datapath + 'caltables/anpos.manual.cal/'
    imfile = datapath + 'image/ngc5921.clean.image/'
    
mesCopy = 'ngc7538_ut.ms'
calCopy = 'anpos.manual.cal'
imCopy = 'ngc5921.clean.image'

def file_copy(filename, perm):
    os.chmod(filename, perm)
    for root, dirs, files in os.walk(filename):
        for d in dirs:
            os.chmod(os.path.join(root, d), perm)
        for f in files:
            os.chmod(os.path.join(root, f), perm)


class rmtables_test(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(rmtables)

    def test_removeMS(self):
        '''
            test_removeMS
            -----------------
            
            This test checks to make sure that a copy of the mesurement set is removed by rmtables, and that the cache is clear after.
            
            The first assert checks that the path to the copy no longer exists, while the second checks the length of the casa cache
        '''
        
        # There is nothing in the cache before this may not be needed
        
        # self.assertTrue(len(tb.showcache()) == 0, msg='The cache was not empty to start with')
        # Copy the MS over and make sure it is removed
        
        if CASA6:
            shutil.copytree(mesSet, mesCopy)
        else:
            shutil.copytree(mesSet, mesCopy)
        
        file_copy(mesCopy, 493)
        rmtables(mesCopy)
        self.assertFalse(os.path.exists(mesCopy), msg='The MS has not been removed')
        
        # There is nothing in the cache after
        self.assertTrue(len(tb.showcache()) == 0, msg='The cache is not empty')

    def test_removeCAL(self):
        '''
            test_removeCAL
            ---------------------
            
            This test checks to make sure that a copy of a calibration table is removed by rmtables, and that the cache is clear after.
            
            The first assert checks that the path the the copy no longer exists, while the second checks the the length of the casa cache
        '''
        
        # There is nothing in the cache before this may not be needed, just commented out for now
        
        #self.assertTrue(len(tb.showcache()) == 0, msg='The cache was not empty to start with')
        
        # Copy the Cal table and make sure it is removed
        if CASA6:
            shutil.copytree(calTab, calCopy)
        else:
            shutil.copytree(calTab, calCopy)
            
        file_copy(calCopy, 493)
        rmtables(calCopy)
        self.assertFalse(os.path.exists(calCopy), msg='The cal table has not been removed')
        
        # There is nothing in the cache after
        self.assertTrue(len(tb.showcache()) == 0, msg='The cache is not empty')

    def test_removeIMG(self):
        '''
            test_removeIMG
            --------------------
            
            This test checks to make sure that a copy of a casa image is removed by rmtables, and that the cache is clear after.
            
            The first assert checks that the path to the copy no longer exists, while the second checks the length of the casa cache
        '''
        
        # There is nothing in the cache before. This line might not be needed
        
       #self.assertTrue(len(tb.showcache()) == 0, msg='The cache is not empty to start with')
        
        # Copy the image and make sure it is removed
        if CASA6:
            shutil.copytree(imfile, imCopy)
        else:
            shutil.copytree(imfile, imCopy)
            
        file_copy(imCopy, 493)
        rmtables(imCopy)
        self.assertFalse(os.path.exists(imCopy), msg='The image file was not removed')
        
        # There is nothing in the cache after
        self.assertTrue(len(tb.showcache()) == 0, msg='The cache was not cleared')

    def test_cacheMS(self):
        '''
            test_cacheMS
            -----------------
            
            This test checks that rmtables will remove a MS and its cached information.
            
            The first assert checks that the copy of the MS is removed and the second checks that the cache is empty
        '''
        if CASA6:
            shutil.copytree(mesSet, mesCopy)
        else:
            shutil.copytree(mesSet, mesCopy)
            
        file_copy(mesCopy, 493)
        
        # Make sure cache is occuped when open
        tb.open(mesCopy)
        self.assertTrue(len(tb.showcache()) > 0, msg='The cache is not occupied when it should be')
        tb.close()
        rmtables(mesCopy)
        
        # Make sure cache is cleared and file is removed
        self.assertFalse(os.path.exists(mesCopy), msg='The MS was not removed')
        self.assertTrue(len(tb.showcache()) == 0, msg='The cache was not cleared')
    
    def test_cacheCAL(self):
        '''
            test_cacheCAL
            -----------------
            
            This test checks that rmtables will remove a CAL table and its cached information.
            
            The first assert checks that the copy of the CAL table is removed and the second checks that the cache is empty
        '''
        if CASA6:
            shutil.copytree(calTab, calCopy)
        else:
            shutil.copytree(calTab, calCopy)
        
        file_copy(calCopy, 493)
        
        # Make sure cache is occupied when open
        tb.open(calCopy)
        self.assertTrue(len(tb.showcache()) > 0, msg='The cache is not occupied when it should be')
        tb.close()
        
        # Make sure cache is cleared and file is removed
        rmtables(calCopy)
        self.assertFalse(os.path.exists(calCopy), msg='The cal table was not removed')
        self.assertTrue(len(tb.showcache()) == 0, msg='The cache was not cleared')
    
    def test_cacheIMG(self):
        '''
            test_cacheIMG
            -----------------
            
            This test checks that rmtables will remove an image and its cached information.
            
            The first assert checks that the copy of the image is removed and the second checks that the cache is empty
        '''
        if CASA6:
            shutil.copytree(imfile, imCopy)
        else:
            shutil.copytree(imfile,imCopy)
        
        file_copy(imCopy, 493)
        
        # Make sure cache is occupied when open
        tb.open(imCopy)
        self.assertTrue(len(tb.showcache()) > 0, msg='The cache is not occupied when it should be')
        tb.close()
        rmtables(imCopy)
        
        # Make sure cache is cleared and file is removed
        self.assertFalse(os.path.exists(imCopy), msg='The file was not removed')
        self.assertTrue(len(tb.showcache()) == 0, msg='the cache was not cleared')


def suite():
    return[rmtables_test]

if __name__ == '__main__':
    unittest.main()
