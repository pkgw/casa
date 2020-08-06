##########################################################################
# test_req_task_imhead.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.1/global-task-list/task_imhead/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import casalog, imhead, rmtables
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

    #DATA#
if CASA6:
    datapath = casatools.ctsys.resolve('image/ngc5921.clean.image/')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image/'

logfile = casalog.logfile()
datacopy = 'clean.image'
testlog = 'testlog.log'

expectedKeys = ['beammajor', 'beamminor', 'beampa', 'bunit', 'cdelt1', 'cdelt2', 'cdelt3', 'cdelt4', 'crpix1', 'crpix2', 'crpix3', 'crpix4', 'crval1', 'crval2', 'crval3', 'crval4', 'ctype1', 'ctype2', 'ctype3', 'ctype4', 'cunit1', 'cunit2', 'cunit3', 'cunit4', 'datamax', 'datamin', 'date-obs', 'equinox', 'imtype', 'masks', 'maxpixpos', 'maxpos', 'minpixpos', 'minpos', 'object', 'observer', 'projection', 'reffreqtype', 'restfreq', 'shape', 'telescope']

class imhead_test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
 
    def setUp(self):
        if not CASA6:
            default(imhead)
        shutil.copytree(datapath, datacopy)
        os.chmod(datacopy, 493)
        for root, dirs, files in os.walk(datacopy):
            for d in dirs:
                os.chmod(os.path.join(root, d), 493)
            for f in files:
                os.chmod(os.path.join(root, f), 493)
 
    def tearDown(self):
        rmtables(datacopy)
        if os.path.exists(testlog):
            os.remove(testlog)
 
    @classmethod
    def tearDownClass(cls):
        pass
    
    
    def test_listkeys(self):
        '''
            test_listkeys
            ------------------------
            
            This test checks that the list mode displays all the expected keys from the imhead function
            
            A list is generated and the keys of that list are compared to what is expected.
            
            All items in the expected keys must be in the listed keys to pass.
            The length of the expected and recieved keys must also be the same.
        '''
        listed = imhead(datacopy, mode='list')
        self.assertTrue(all([True for item in expectedKeys if item in list(listed.keys())]))
        self.assertTrue(len(expectedKeys) == len(list(listed.keys())))
        
    def test_listkeysHdkeyVal(self):
        '''
            test_listkeysHdkeyVal
            ------------------------------
            
            This test checks to make sure that hdkey and hdvalue have no effect on the output of mode list
            
            A key and value is provided and then the keys are checked to make sure they are the expected values.
            The lists of expected values and recived ones must also be of equal length.
        '''
        listed = imhead(datacopy, mode='list', hdkey='beammajor', hdvalue=1)
        self.assertTrue(all([True for item in expectedKeys if item in list(listed.keys())]))
        self.assertTrue(len(expectedKeys) == len(list(listed.keys())))
    
    def test_history(self):
        '''
            test_history
            ---------------------
            
            This checks to make sure that when the mode is history a log is populated
            
            When the function is ran with mode history a logfile is populated with information
            Check that info is written to the file.
            
            
        '''
        casalog.setlogfile(testlog)
        imhead(datacopy, mode='history')
        self.assertTrue('INFO' in open(testlog).read())
        
    def test_historyHdkeyVal(self):
        '''
            test_historyHdkeyVal
            ------------------------
            
            This checks to make sure this task still functions when provided with Hdkey and value.
            
            These parameters should have no effect, so it should have the identical output to test_history.
        '''
        casalog.setlogfile(testlog)
        imhead(datacopy, mode='history', hdkey='beammajor', hdvalue=1)
        self.assertTrue('INFO' in open(testlog).read())
        
    def test_del(self):
        '''
            test_del
            ---------------
            
            This checks that the result of running the delete mode on all the hd keys has the effect that we expected
            
            Certain values should be removed where others will remain unchanged.
            
            The assertion checks that the image header after removal looks like the expected dictionary
        '''
        
        endDict = {'bunit': '', 'cdelt1': -7.27220521664304e-05, 'cdelt2': 7.27220521664304e-05, 'cdelt3': 1.0, 'cdelt4': 24414.0625, 'crpix1': 128.0, 'crpix2': 128.0, 'crpix3': 0.0, 'crpix4': 0.0, 'crval1': 4.022983925846928, 'crval2': 0.08843001543437938, 'crval3': 1.0, 'crval4': 1412787144.0812755, 'ctype1': 'Right Ascension', 'ctype2': 'Declination', 'ctype3': 'Stokes', 'ctype4': 'Frequency', 'cunit1': 'rad', 'cunit2': 'rad', 'cunit3': '', 'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682, 'date-obs': '1995/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Intensity', 'masks': np.array([], dtype='<U16'), 'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'), 'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz', 'minpixpos': np.array([230,   0,   0,  15], dtype='int32'), 'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz', 'object': '', 'observer': '', 'projection': 'SIN', 'reffreqtype': 'LSRK', 'restfreq': np.array([1420405752.0]), 'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': ''}
        
        for key in expectedKeys:
            imhead(datacopy, mode='del', hdkey=key)
            #print(key, imhead(datacopy, mode='list')[key] == endDict[key])
        #for key in endDict.keys():
            #print(key, np.all(imhead(datacopy, mode='list')[key] == endDict[key]))
            
        #self.assertDictEqual(endDict, imhead(datacopy, mode='list'))
        
        self.assertTrue(len(list(endDict.keys())) == len(list(imhead(datacopy, mode='list').keys())))
        self.assertTrue(np.all([np.all(imhead(datacopy, mode='list')[key]==endDict[key]) for key in list(imhead(datacopy, mode='list').keys())]))
        
    def test_add(self):
        '''
            test_add
            ---------------
            
            This test makes sure that add can add to all the keys specified in the documentation
            
            The endDict is the dictionary we expect out of the list function.
            we first delete all the values we can, and then try to add them back.
            Upon failing to add a section back the dictionaries will become differnt, so my asserting the dictionary is the same as it started assures add worked.
            This relies on del working as well, however if that fails so will the test_del function.
            
            imtype and restfreq cannot be tested, because there is no way to remove the values for them.
            Also adding a beam auto matically sets beamminor and beampa, those values cannot be all chosen using add
        '''
        # There are fields you can add to but can't delete from?
        # How are you supposed to use add when add requires th field to be empty
    
        endDict = {'beammajor': {'unit': 'arcsec', 'value': 51.7}, 'beamminor': {'unit': 'arcsec', 'value': 51.7}, 'beampa': {'unit': 'deg', 'value': 0.0}, 'bunit': 'Jy/beam', 'cdelt1': -7.27220521664304e-05, 'cdelt2': 7.27220521664304e-05, 'cdelt3': 1.0, 'cdelt4': 24414.0625, 'crpix1': 128.0, 'crpix2': 128.0, 'crpix3': 0.0, 'crpix4': 0.0, 'crval1': 4.022983925846928, 'crval2': 0.08843001543437938, 'crval3': 1.0, 'crval4': 1412787144.0812755, 'ctype1': 'Right Ascension', 'ctype2': 'Declination', 'ctype3': 'Stokes', 'ctype4': 'Frequency', 'cunit1': 'rad', 'cunit2': 'rad', 'cunit3': '', 'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682, 'date-obs': '1995/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Intensity', 'masks': np.array([], dtype='<U16'), 'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'), 'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz', 'minpixpos': np.array([230,   0,   0,  15], dtype='int32'), 'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz', 'object': 'N5921_2', 'observer': 'TEST', 'projection': 'SIN', 'reffreqtype': 'LSRK', 'restfreq': np.array([1420405752.0]), 'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': 'VLA'}
        
        for key in expectedKeys:
            imhead(datacopy, mode='del', hdkey=key)
            
        for key in expectedKeys:
            imhead(datacopy, mode='add', hdkey=key, hdvalue=endDict[key])
        
        self.assertTrue(len(list(endDict.keys())) == len(list(imhead(datacopy, mode='list').keys())))
        self.assertTrue(np.all([np.all(imhead(datacopy, mode='list')[key]==endDict[key]) for key in list(imhead(datacopy, mode='list').keys())]))
        
    def test_put(self):
        '''
            test_put
            ----------------
            
            This test checks to see that the put mode will place values into the specified slots.
            
            In this case we remove from each slot if we can. Some will be added anew while other values will be replaced by the end Dict value
            At the end we assert that the final header looks like the dictionary we expected to see.
        '''
        
        endDict = {'beammajor': {'unit': 'arcsec', 'value': 59.2}, 'beamminor': {'unit': 'arcsec', 'value': 42.2}, 'beampa': {'unit': 'deg', 'value': 8.0}, 'bunit': 'Jy/beam', 'cdelt1': float('-8.27220521664304e-05'), 'cdelt2': float('7.17220521664304e-05'), 'cdelt3': float('1.0'), 'cdelt4': float('24413.0625'), 'crpix1': float(127.0), 'crpix2': 127.0, 'crpix3': float(0.0), 'crpix4': float(1.0), 'crval1': 3.02, 'crval2': 0.078, 'crval3': 0.0, 'crval4': 1412787143.0812755, 'ctype1': 'Declination', 'ctype2': 'Right Ascension', 'ctype3': 'Frequency', 'ctype4': 'Stokes', 'cunit1': 'deg', 'cunit2': 'deg', 'cunit3': '', 'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682, 'date-obs': '1996/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Velocity', 'masks': np.array([], dtype='<U16'), 'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'), 'maxpos': '00:12:04.661 +00.04.42.607 ?? 1.41369e+09Hz', 'minpixpos': np.array([230,   0,   0,  15], dtype='int32'), 'minpos': '00:12:02.755 +00.04.08.009 ?? 1.41313e+09Hz', 'object': 'N5921_22', 'observer': 'TESTING', 'projection': 'SIN', 'reffreqtype': 'LSRK', 'restfreq': np.array([1.3]), 'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': 'EVLA'}
        
        InDict = {'beammajor': {'unit': 'arcsec', 'value': 59.2}, 'beamminor': {'unit': 'arcsec', 'value': 42.2}, 'beampa': {'unit': 'deg', 'value': 8.0}, 'bunit': 'Jy/beam', 'cdelt1': '-8.27220521664304e-05deg', 'cdelt2': '7.17220521664304e-05deg', 'cdelt3': '1.0', 'cdelt4': '24413.0625', 'crpix1': float(127.0), 'crpix2': 127.0, 'crpix3': float(0.0), 'crpix4': float(1.0), 'crval1': '3.02deg', 'crval2': '0.078deg', 'crval3': '0.0', 'crval4': '1412787143.0812755', 'ctype1': 'Declination', 'ctype2': 'Right Ascension', 'ctype3': 'Frequency', 'ctype4': 'Stokes', 'cunit1': 'deg', 'cunit2': 'deg', 'cunit3': '', 'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682, 'date-obs': '1996/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Velocity', 'masks': np.array([], dtype='<U16'), 'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'), 'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz', 'minpixpos': np.array([230,   0,   0,  15], dtype='int32'), 'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz', 'object': 'N5921_22', 'observer': 'TESTING', 'projection': 'SIN', 'reffreqtype': 'LSRK', 'restfreq': '1.3', 'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': 'EVLA'}

        for key in expectedKeys:
            imhead(datacopy, mode='put', hdkey=key, hdvalue=InDict[key])
            
        
        self.assertTrue(len(list(endDict.keys())) == len(list(imhead(datacopy, mode='list').keys())))
        self.assertTrue(np.all([np.all(imhead(datacopy, mode='list')[key]==endDict[key]) for key in list(imhead(datacopy, mode='list').keys())]))
        
    def test_get(self):
        '''
            test_get
            --------------
            
                This test makes sure that the get function returns the datatypes(?) at the given keys
                The documentation makes it seem like it should only return the data type, or that only the datatype should be checked
                
                This test creates a list of all the data types that get returns and then compares to a list of the expected datatypes
        '''
        ###NOTE: ask further questions about what exactly this should be testing. Is checking the datatypes good enough?
        
        typeList = [type({}), type({}), type({}), type(''), type({}), type({}), type({}), type({}), type(0.00), type(0.00), type(0.00), type(0.00), type({}), type({}), type(np.ndarray([])), type({}), type(''), type(''), type(''), type(''), type(''), type(''), type(''), type(''), type(0.00), type(0.00), type(''), type(''), type(''), type(np.ndarray([])), type(np.ndarray([])), type(''), type(np.ndarray([])), type(''), type(''), type(''), type(''), type(''), type({}), type(np.ndarray([])), type('')]
        
        getTypes = [type(imhead(datacopy, mode='get', hdkey=x)) for x in expectedKeys]
        
        self.assertTrue(getTypes == typeList)
            
    def test_summary(self):
        '''
            test_summary
            -----------------
        '''
        
        ###NOTE: to test the function of the verbose parameter an image with multiple beams is required
        
        summaryDictnoV = {'axisnames': np.array(['Right Ascension', 'Declination', 'Stokes', 'Frequency'],
      dtype='<U16'), 'axisunits': np.array(['rad', 'rad', '', 'Hz'], dtype='<U16'), 'defaultmask': '', 'hasmask': False, 'imagetype': 'Intensity', 'incr': np.array([-7.27220521664304e-05,  7.27220521664304e-05,  1.00000000e+00,  24414.0625]), 'masks': np.array([], dtype='<U16'), 'messages': np.array([], dtype='<U16'), 'ndim': 4, 'refpix': np.array([128., 128.,   0.,   0.]), 'refval': np.array([4.022983925846928, 0.08843001543437938, 1.00000000e+00, 1412787144.0812755]), 'restoringbeam': {'major': {'unit': 'arcsec', 'value': 51.7}, 'minor': {'unit': 'arcsec', 'value': 47.2}, 'positionangle': {'unit': 'deg', 'value': -171.0}}, 'shape': np.array([256, 256,   1,  46], dtype='int32'), 'tileshape': np.array([64, 64,  1,  8], dtype='int32'), 'unit': 'Jy/beam'}
        
        compared = imhead(datacopy, mode='summary')
        
        #difference = [summaryDictnoV[key] == compared[key] for key in summaryDictnoV.keys()]
        
        #self.assertTrue(str(summaryDictnoV) == str(compared))
        #self.assertDictEqual(summaryDictnoV, imhead(datacopy, mode='summary'))
        self.assertTrue(len(list(compared.keys())) == len(list(summaryDictnoV.keys())))
        self.assertTrue(np.all([np.all(compared[key]==summaryDictnoV[key]) for key in list(compared.keys())]))
        
def suite():
    return[imhead_test]
 
# Main #
if __name__ == '__main__':
    unittest.main()
