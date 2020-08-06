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
# https://casa.nrao.edu/casadocs/casa-5.3.0/global-task-list/task_listobs/about
#
# Each requirement is broken up into 4 versions testing it on a MS, MMS, time-averaged MS, and time-averaged MMS
#
# test_logread checks that the information is output in a form that follows logger guidelines
# test_file checks that the listfile parameter operates as intended by creating a file
# test_input checks that listobs can take a MS, MMS, time-averaged MS, or time-averaged MMS in the input file
# test_logfile checks that the logfile is populated when new listobs commands are run
# test_indivrow checks that one row is displayed per scan
# test_novis checks to be sure that listobs gives no visibility information
# test_scan checks that the scan information exists and the parameter accepts proper inputs
# test_field checks that the field information exists and the parameter accepts proper inputs
# test_spw checks that the spectral window information exists and the parameter accepts proper inputs
# test_corrs checks that the correlation information exists and the parameter accepts proper inputs
# test_ant checks that the antenna information exists and the parameter accepts proper inputs
# test_selectdata checks that the selectdata parameter works as intended and properly handles inputs
# test_uvrange checks that the parameter properly handles inputs
# test_timerange checks that the parameter properly handles inputs
# test_intent checks that intent information exists and inputs are properly handled
# test_array checks for the existence of the array parameter and that inputs are properly handled
# test_observation checks for the existence of the parameter and that inputs are properly handled
# test_verbose checks that verbose outputs are longer than non-verbose
# test_overwrite checks that the overwritten file is the same as the original
# test_CAS_6733 checks for an infinite loop bug
# test_avgInterval checks for the existence of the average int information
# test_listunfl checks that unflagged information is displayed by listobs
#
###########################################################################
CASA6 = False
import sys
import os
import unittest
import hashlib
import subprocess
import shutil
try:
    import casatools
    from casatasks import partition, split, listobs, casalog
    from casatools.platform import bytes2str
    ms = casatools.ms()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *

# If the test is being run in CASA6 use the new method to get the CASA path
if CASA6:
    datapath = casatools.ctsys.resolve('/data/regression/unittest/listobs')

else:
    dataroot = os.environ.get('CASAPATH').split()[0] + '/data/regression/'
    datapath = dataroot + 'unittest/listobs/'

    # Generate the test data

if CASA6:
    mesSet = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms'):
        mesSet = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms'
    else:
        mesSet = datapath + 'uid___X02_X3d737_X1_01_small.ms'

partition(vis=mesSet, outputvis='genmms.mms', createmms=True)
multiMesSet = 'genmms.mms'
split(vis=mesSet, outputvis='gentimeavgms.ms', datacolumn='DATA',timebin='1s')
timeavg_ms = 'gentimeavgms.ms'
split(vis=multiMesSet, outputvis='gentimeavgmms.mms', datacolumn='DATA', timebin='1s')
timeavg_mms = 'gentimeavgmms.mms'

logpath = casalog.logfile()


def _sha1it(filename):
    blocksize = 65536
    hasher = hashlib.sha1()
    with open(filename, 'rb') as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
    return hasher.hexdigest()


class listobs_test_base(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(listobs)

    def tearDown(self):
        # remove files and temp logs
        os.system('rm -rf ' + 'listobs*.txt')
        os.system('rm -rf testlog.log')
        casalog.setlogfile(str(logpath))

    def filefunc(self, dataset, filename):
        listobs(vis=dataset, listfile=filename)
        self.assertTrue(os.path.isfile('listobs.txt'))

    def logreadfunc(self, dataset):
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset)
        casalog.setlogfile(logpath)
        
        if sys.version_info[0] == 3:
            print(('VERSION', ' ', sys.version_info))
            # Check that the file can be read in python session default encoding
            with open('testlog.log','r') as fout:
                list(map(bytes2str, fout.readlines()))
                print((list(map(bytes2str, fout.readlines()))))
                
        else:
            # Check if the file can be decoded as ascii for python 2.7
            with open('testlog.log','r') as log:
                for data in log:
                    try:
                        data.decode('ASCII')
                    except:
                        self.fail()

    def logfilecontain(self, dataset):
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset)
        self.assertTrue('listobs' in open('testlog.log').read(), msg='logfile not populated by listobs command on a MS')

    def novis(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertFalse('Amp' in open('listobs.txt').read())

    def scancheck(self, dataset):
        # Check that listobs runs and that the scan column exists
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Scan' in open('listobs.txt').read(), msg='scan does not exist in output')
        # Check to see if a valid input and invalid input is returns the proper value (pass/fail)
        self.assertTrue(listobs(vis=dataset, scan='1'), msg='Scan fails to select')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, scan=1)
        else:
            self.assertFalse(listobs(vis=dataset, scan=1), msg='Scan incorrectly accepts an int')

        # Make temp log
        casalog.setlogfile('testlog.log')

        # test for certain warnings appearing in the log
        listobs(vis=dataset, scan='1,2')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='A warning is raised for multiple scans')
        listobs(vis=dataset, scan=['1', '2'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='fails to label incorrect data type')
        listobs(vis=dataset, scan='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(), msg='fails to recognize improper string')

    def fieldcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, field=1)
        else:
            self.assertFalse(listobs(vis=dataset, field=1), msg='An int was given to the field (requires a string)')
        self.assertTrue('FldId' in open('listobs.txt').read(), msg='Field Id does not exist in a MS')
        # section that should raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, field='1')
        listobs(vis=dataset, field='0~2')
        listobs(vis=dataset, field='0,2')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for field in a MS')
        # section that should be raising warnings
        listobs(vis=dataset, field='0-2')
        self.assertTrue('No match found for name "0-2"' in open('testlog.log').read(), msg='Failed to identify improper delimiter')
        listobs(vis=dataset, field='abc')
        self.assertTrue('No match found for name "abc"' in open('testlog.log').read(), msg='Failed to identify improper string')
        # return to default log file

    def spwcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, spw=1)
        else:
            self.assertFalse(listobs(vis=dataset, spw=1), msg='An int was given to spw (requires string)')
        self.assertTrue('SpwID' in open('listobs.txt').read(), msg='Spw does not exist in a MS')
        # section that should raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, spw='0')
        listobs(vis=dataset, spw='0,1')
        listobs(vis=dataset, spw='0~1')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for spw in a MS')
        # section that should raise warnings
        listobs(vis=dataset, spw='3')
        self.assertTrue('No match found for 3' in open('testlog.log').read(), msg='fails to recognize out of range values')
        listobs(vis=dataset, spw='0')
        self.assertTrue('-1' in open('testlog.log').read(), msg='Fails to recognize improper delimiter')
        listobs(vis=dataset, spw='abc')
        self.assertTrue('No match found for "abc"' in open('testlog.log').read(), msg='Fails to recognize improper string')
        # return to default log file

    def corrcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Corrs' in open('listobs.txt').read(), msg='Corrs does not exist in a MS')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, correlation=1)
        else:
            self.assertFalse(listobs(vis=dataset, correlation=1), msg='An int was accepted when a string is required')
        # section that should not raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, correlation='XX')
        listobs(vis=dataset, correlation='XX,YY')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='not accepting proper input for correlation in a MS')
        # section that should raise warnings
        listobs(vis=dataset, correlation=['XX', 'YY'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='No warning for using a list was given')
        listobs(vis=dataset, correlation='RR')
        self.assertTrue('named RR' in open('testlog.log').read(), msg='No warning for using a absent correlation')

    def antcheck(self, dataset):
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Antennas' in open('listobs.txt').read(), msg='Antennas section does not exist in MS')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, antenna=0)
        else:
            self.assertFalse(listobs(vis=dataset, antenna=0),
                             msg='Accepts an int as an argument when a string is required')
        # section that should not raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, antenna='0')
        listobs(vis=dataset, antenna='0,DV01')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='Proper inputs raised warnings')
        # section that should raise warnings
        listobs(vis=dataset, antenna='abc')
        self.assertTrue('Antenna Expression: No match found for token(s)' in open('testlog.log').read(), msg='No warning raise for invalid string')
        listobs(vis=dataset, antenna='3')
        self.assertTrue('No match found for the antenna specificion [ID(s): [3]]' in open('testlog.log').read(), msg='No warning for ID out of range')
        # This one is marked as correct by the documentation, but CASA disagrees
        listobs(vis=dataset, antenna=['0,DV01'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='Failed to recognize list as incorrect data type')
        # return to default log file

    def selectcheck(self, dataset):
        # selectdata should fail for all invalid inputs
        self.assertTrue(listobs(vis=dataset, selectdata=False), msg='Passing False to select data fails on a MS')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, selectdata=1)
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, selectdata='str')
        else:
            self.assertFalse(listobs(vis=dataset, selectdata=1), msg='An int is accepted for the select data parameter')
            self.assertFalse(listobs(vis=dataset, selectdata='str'),
                             msg='A string is accepted for the select data parameter')
            self.assertFalse(listobs(vis=dataset, selectdata=[False]),
                             msg='Array data type is accepted when it contains values')
        # I'm really not sure why these happen. especially the empty list
        if not CASA6:
            self.assertTrue(listobs(vis=dataset, selectdata=[]))
            self.assertTrue(listobs(vis=dataset, selectdata=None))

    def uvrangecheck(self, dataset):
        self.assertTrue(listobs(vis=dataset, uvrange='0~100klambda'), msg='fails to read valid input for uvrange in a MS')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, uvrange=0)
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, uvrange=[1, 2])
        else:
            self.assertFalse(listobs(vis=dataset, uvrange=0), msg='Accepts an int when only str should be accepted')
            self.assertFalse(listobs(vis=dataset, uvrange=[1, 2]),
                             msg='accepts an array of ints when only str should be accepted')
        # Use temp log
        casalog.setlogfile('testlog.log')
        # shouldn't raise Warning
        listobs(vis=dataset, uvrange='0~100')
        listobs(vis=dataset, uvrange='0~100klambda')
        listobs(vis=dataset, uvrange='0~50,60~100')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='Warnings are raised for valid inputs')
        # should raise warnings
        listobs(vis=dataset, uvrange=['0~50', '60~100'])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='Fails to raise warning for wrong data type')
        listobs(vis=dataset, uvrange='0-100')
        self.assertTrue('near char. 2 in string "0-100"' in open('testlog.log').read(), msg='Fails to raise warning for wrong delimiter')
        listobs(vis=dataset, uvrange='abc')
        self.assertTrue('near char. 1 in string "abc"' in open('testlog.log').read(), msg='Fails to raise warning for improper string')
        # restore default log path

    def timerangecheck(self, dataset):
        #check valid entry
        self.assertTrue(listobs(vis=dataset, timerange='03:00:00~04:00:00'))
        # create temp log and check inputs that raise no warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, timerange='3:0:0~4:0:0,4:0:0~5:0:0')
        self.assertFalse('WARN' in open('testlog.log').read())
        # check that specific warnings are raised
        listobs(vis=dataset, timerange='abc')
        self.assertTrue('Parse error at or near ' in open('testlog.log').read())
        listobs(vis=dataset, timerange=[])
        self.assertTrue('incorrect data type' in open('testlog.log').read())
        listobs(vis=dataset, timerange='03:00:00-04:00:00')
        self.assertTrue('near char. 9 in string "03:00:00-04:00:00"' in open('testlog.log').read())
        listobs(vis=dataset, timerange='3~4')
        self.assertTrue('MSSelectionNullSelection' in open('testlog.log').read())
        # Check that passing an int fails
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, timerange=4)
        else:
            self.assertFalse(listobs(vis=dataset, timerange=4))

    def intentcheck(self, dataset):
        # Returns true with a valid input and false with an int
        self.assertTrue(listobs(vis=dataset, intent='CALIBRATE_PHASE.ON_SOURCE,CALIBRATE_POINTING.ON_SOURCE,CALIBRATE_WVR.ON_SOURCE'), msg='Fails with valid input')
        if CASA6:
            with self.assertRaises(AssertionError):
                listobs(vis=dataset, intent=1)
        else:
            self.assertFalse(listobs(vis=dataset, intent=1), msg='Accepts int when it should only accept str')
        # Test for the existence of the column scan intent
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('ScanIntent' in open('listobs.txt').read(), msg='There is no ScanIntent information for a MS')
        # These shouldn't raise any warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, intent='CALIBRATE_PHASE.ON_SOURCE')
        listobs(vis=dataset, intent='CALIBRATE_PHASE.ON_SOURCE,CALIBRATE_POINTING.ON_SOURCE,CALIBRATE_WVR.ON_SOURCE')
        self.assertFalse('WARN' in open('testlog.log').read(), msg='There are warnings for inputs that should raise none')
        # These should raise a warning
        listobs(vis=dataset, intent=[])
        self.assertTrue('incorrect data type' in open('testlog.log').read(), msg='Incorrect data type list accepted')
        listobs(vis=dataset, intent='abc')
        self.assertTrue('No match found for "abc"' in open('testlog.log').read(), msg='Invalid string accepted without warning')
        # Set log path back to default

    def arraycheck(self, dataset):
        self.assertTrue(listobs(vis=dataset, array='0'), msg='Listobs fails to recognize valid array in a MS')
        if CASA6:
            with self.assertRaises(AssertionError):
               listobs(vis=dataset, array=0)
        else:
            self.assertFalse(listobs(vis=dataset, array=0), msg='Listobs fails to recognize invalid data type in a MS')
        self.res = listobs(vis=dataset, listfile='listobs.txt')
        # These should raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, array='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(), msg='Listobs fails to recognize invalid string in a MS')
        listobs(vis=dataset, array='10')
        self.assertTrue('The selected table has zero rows' in open('testlog.log').read(), msg='Listobs fails to recognize empty table from a MS')

        self.assertTrue('ArrayID' in open('listobs.txt').read(), msg='There is no Array information for a MS')

    def obscheck(self, dataset):
        self.assertTrue(listobs(vis=dataset, observation=0), msg='Observation fails to accept Int for a MS')
        self.assertTrue(listobs(vis=dataset, observation='0'), msg='Observation fails to accept proper string for a MS')
        # These should raise warnings
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, observation='abc')
        self.assertTrue('Parse error' in open('testlog.log').read(), msg='Listobs fails to identify improper string')
        listobs(vis=dataset, observation='10')
        self.assertTrue(('The selected table has zero rows') in open('testlog.log').read(),
                    msg='Listobs fails to identify an empty table')

        self.res = listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('ObservationID' in open('listobs.txt').read(), msg='There is no Observation information')

    def verbosecheck(self, dataset):
        casalog.setlogfile('testlog.log')
        listobs(vis=dataset, verbose=False)
        nonverb = ['name', 'station']
        self.assertTrue(all(x in open('testlog.log').read() for x in nonverb), msg='non-verbose not showing proper info')
        listobs(vis=dataset, verbose=True)
        items = ['Name', 'Station', r'Diam.', r'Long.', r'Lat.', 'Offset', 'ITRF', 'Scan', 'FieldName', 'SpwIds']
        self.assertTrue(all(x in open('testlog.log').read() for x in items))

    def overwritecheck(self, dataset):
        listfile = "listobs.txt"
        self.assertTrue(listobs(vis=dataset, listfile=listfile))
        # test default value is overwrite=False
        self.assertFalse(listobs(vis=dataset, listfile=listfile))
        self.assertFalse(listobs(vis=dataset, listfile=listfile, overwrite=False))
        expec = _sha1it(listfile)
        self.assertTrue(listobs(vis=dataset, listfile=listfile, overwrite=True))
        got = _sha1it(listfile)
        self.assertTrue(got == expec)

    def avgIntervalcheck(self, dataset):
        listobs(vis=dataset, listfile='listobs.txt')
        self.assertTrue('Average Interval' in open('listobs.txt').read(), msg='There is no average interval column in a MS')

    def unfcheck(self, dataset):
        self.assertTrue(listobs(vis=dataset, listunfl=True))
        listobs(vis=dataset, listfile='listobs1.txt', listunfl=True)
        self.assertTrue('nUnflRows' in open('listobs1.txt').read())

        listobs(vis=dataset, listfile='listobs2.txt', listunfl=False)
        self.assertFalse('nUnflRows' in open('listobs2.txt').read())


class test_listobs(listobs_test_base):

    @classmethod
    def setUpClass(cls):

        pass

    def setUp(self):
        self.res = None
        if not CASA6:
            default(listobs)

    def tearDown(self):
        # remove files and temp logs
        os.system('rm -rf ' + 'listobs*.txt')
        os.system('rm -rf testlog.log')
        casalog.setlogfile(str(logpath))

    @classmethod
    def tearDownClass(cls):
        # remove all the generated data
        os.system('rm -rf genmms.mms')
        os.system('rm -rf gentimeavgms.ms')
        os.system('rm -rf gentimeavgmms.mms')
        os.system('rm -rf genmms.mms.flagversions')

    # Test that Different ms are taken
    def test_wrongInp(self):
        if not CASA6:
            self.assertFalse(listobs(vis='foo.ms'))

    def test_input(self):
        '''Listobs test: Check all if listobs can take ms, mms, time averaged ms, and time averaged mms'''
        # See if listobs can take all these forms of MS
        self.assertTrue(listobs(vis=mesSet), msg='Fails to take MS file')
        self.assertTrue(listobs(vis=multiMesSet), msg='Fails to take MMS file')
        self.assertTrue(listobs(vis=timeavg_ms), msg='Fails to take time-averaged MS file')
        self.assertTrue(listobs(vis=timeavg_mms), msg='Fails to take time-averaged MMS file')

    # Test the list file
    def test_fileMS(self):
        '''Listobs test: Check to see if list file is generated from a MS'''
        self.filefunc(mesSet, 'listobs.txt')

    def test_fileMMS(self):
        '''Listobs test: Check to see if list file is generated from a MMS'''
        self.filefunc(multiMesSet, 'listobs.txt')

    def test_fileTimeAvgMS(self):
        '''Listobs test: Check to see if list file is generated from a time-averaged MS'''
        self.filefunc(timeavg_ms, 'listobs.txt')

    def test_fileTimeAvgMMS(self):
        '''Listobs test: Check to see if list file is generated from a time-averaged MMS'''
        self.filefunc(timeavg_mms, 'listobs.txt')

    # Test the log file encoding
    def test_logreadMS(self):
        '''Listobs test: Check that the log file from a MS is human readable'''
        self.logreadfunc(mesSet)

    def test_logreadMMS(self):
        '''Listobs test: Check that the log file from a MMS is human readable'''
        self.logreadfunc(multiMesSet)

    def test_logreadTimeAvgMS(self):
        '''Listobs test: Check that the log file from a time-averaged MS is human readable'''
        self.logreadfunc(timeavg_ms)

    def test_logreadTimeAvgMMS(self):
        '''Listobs test: Check that the log file from a time-averaged MMS is human readable'''
        self.logreadfunc(timeavg_mms)

    # Test log file written to

    def test_logfileMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.logfilecontain(mesSet)

    def test_logfileMMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.logfilecontain(multiMesSet)

    def test_logfileTimeAvgMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.logfilecontain(timeavg_ms)

    def test_logfileTimeAvgMMS(self):
        '''Listobs test: Check to see if the logger generates a logfile entries'''
        self.logfilecontain(timeavg_mms)

    # Test one row displayed per scan

    def test_indivrow(self):
        '''Listobs test: Check if one row is displayed per scan'''
        ms.open(mesSet)
        self.res = ms.summary(True, listunfl=True)
        ms.close()
        # beginning and ending times should be different
        btime = self.res['scan_1']['0']['BeginTime']
        etime = self.res['scan_1']['0']['EndTime']
        self.assertNotEqual(btime, etime, msg='Beginning and Ending times should not be equal')

    # Test that no visibility information is given

    def test_novisMS(self):
        '''Listobs test: Check if there is any visibility information when using a MS'''
        self.novis(mesSet)

    def test_novisMMS(self):
        '''Listobs test: Check if there is any visibility information when using a MMS'''
        self.novis(multiMesSet)

    def test_novisTimeAvgMS(self):
        '''Listobs test: Check if there is any visibility information when using a time-averaged MS'''
        self.novis(timeavg_ms)

    def test_novisTimeAvgMMS(self):
        '''Listobs test: Check if there is any visibility information when using a time-averaged MMS'''
        self.novis(timeavg_mms)

    # Test the scan parameter

    def test_scanMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a MS'''
        self.scancheck(mesSet)

    def test_scanMMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a MMS'''
        self.scancheck(multiMesSet)

    def test_scanTimeAvgMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a time-averaged MS'''
        self.scancheck(timeavg_ms)

    def test_scanTimeAvgMMS(self):
        '''Listobs test: Check to make sure the scan information exists in all possible inputs in a time-averaged MMS'''
        self.scancheck(timeavg_mms)

    # Test the field parameter

    def test_fieldMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a MS'''
        self.fieldcheck(mesSet)

    def test_fieldMMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a MMS'''
        self.fieldcheck(multiMesSet)

    def test_fieldTimeAvgMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a time-averaged MS'''
        self.fieldcheck(timeavg_ms)

    def test_fieldTimeAvgMMS(self):
        '''Listobs test: Check to make sure the field information exists in all possible inputs in a time-averaged MMS'''
        self.fieldcheck(timeavg_mms)

    # Test the spw parameter

    def test_spwMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a MS'''
        self.spwcheck(mesSet)

    def test_spwMMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a MMS'''
        self.spwcheck(multiMesSet)

    def test_spwTimeAvgMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a time-averaged MS'''
        self.spwcheck(timeavg_ms)

    def test_spwTimeAvgMMS(self):
        '''Listobs test: Check to make sure the spw information exists in all possible inputs in a time-averaged MMS'''
        self.spwcheck(timeavg_mms)

    # Test the Correlation parameter

    def test_corrsMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a MS'''
        self.corrcheck(mesSet)

    def test_corrsMMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a MMS'''
        self.corrcheck(multiMesSet)

    def test_corrsTimeAvgMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a time-averaged MS'''
        self.corrcheck(timeavg_ms)

    def test_corrsTimeAvgMMS(self):
        '''Listobs test: Check to make sure the correlation information exists in all possible inputs in a time-averaged MMS'''
        self.corrcheck(timeavg_mms)

    # Test the antenna parameter

    def test_antMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a MS'''
        self.antcheck(mesSet)

    def test_antMMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a MMS'''
        self.antcheck(multiMesSet)

    def test_antTimeAvgMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a time-averaged MS'''
        self.antcheck(timeavg_ms)

    def test_antTimeAvgMMS(self):
        '''Listobs test: Check to make sure the antenna information exists in all possible inputs in a time-averaged MMS'''
        self.antcheck(timeavg_mms)

    # Test the selectdata parameter

    def test_selectdataMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a MS'''
        self.selectcheck(mesSet)

    def test_selectdataMMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a MMS'''
        self.selectcheck(multiMesSet)

    def test_selectdataTimeAvgMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a time-averaged MS'''
        self.selectcheck(timeavg_ms)

    def test_selectdataTimeAvgMMS(self):
        '''Listobs test: Check to see if the selectdata parameter functions with a time-averaged MMS'''
        self.selectcheck(timeavg_mms)

    # Test the uvrange parameter

    def test_uvrangeMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a MS'''
        self.uvrangecheck(mesSet)

    def test_uvrangeMMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a MMS'''
        self.uvrangecheck(multiMesSet)

    def test_uvrangeTimeAvgMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a time-averaged MS'''
        self.uvrangecheck(timeavg_ms)

    def test_uvrangeTimeAvgMMS(self):
        '''Listobs test: check that the proper inputs work for uvrange on a time-averaged MMS'''
        self.uvrangecheck(timeavg_mms)

    # Test time range parameter

    def test_timerangeMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a MS'''
        self.timerangecheck(mesSet)

    def test_timerangeMMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a MMS'''
        self.timerangecheck(multiMesSet)

    def test_timerangeTimeAvgMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a time-averaged MS'''
        self.timerangecheck(timeavg_ms)

    def test_timerangeTimeAvgMMS(self):
        '''Listobs test: check that the proper inputs work for timerange on a time-averaged MMS'''
        self.timerangecheck(timeavg_mms)

    # Test the intent parameter

    def test_intentMS(self):
        '''Listobs test: Check to see that Intent info exists for a MS and accepts correct inputs'''
        self.intentcheck(mesSet)

    def test_intentMMS(self):
        '''Listobs test: Check to see that Intent info exists for a MMS and accepts correct inputs'''
        self.intentcheck(multiMesSet)

    def test_intentTimeAvgMS(self):
        '''Listobs test: Check to see that Intent info exists for a time-averaged MS and accepts correct inputs'''
        self.intentcheck(timeavg_ms)

    def test_intentTimeAvgMMS(self):
        '''Listobs test: Check to see that Intent info exists for a time-averaged MMS and accepts correct inputs'''
        self.intentcheck(timeavg_mms)

    # Test the array parameter

    def test_arrayMS(self):
        '''Listobs test: Check for the existence of the array parameter in a MS and accepts proper inputs'''
        self.arraycheck(mesSet)

    def test_arrayMMS(self):
        '''Listobs test: Check for the existence of the array parameter in a MMS and accepts proper inputs'''
        self.arraycheck(multiMesSet)

    def test_arrayTimeAvgMS(self):
        '''Listobs test: Check for the existence of the array parameter in a time-averaged MS and accepts proper inputs'''
        self.arraycheck(timeavg_ms)

    def test_arrayTimeAvgMMS(self):
        '''Listobs test: Check for the existence of the array parameter in a time-averaged MMS and accepts proper inputs'''
        self.arraycheck(timeavg_mms)

    # Test the observation parameter

    def test_observationMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a MS and check for proper inputs'''
        self.obscheck(mesSet)

    def test_observationMMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a MMS and check for proper inputs'''
        self.obscheck(multiMesSet)

    def test_observationTimeAvgMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a time-averaged MS and check for proper inputs'''
        self.obscheck(timeavg_ms)

    def test_observationTimeAvgMMS(self):
        '''Listobs test: Check for the existence of the Observation parameter for a time-averaged MMS and check for proper inputs'''
        self.obscheck(timeavg_mms)

    # Test the Verbose parameter

    def test_verboseMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a MS'''
        self.verbosecheck(mesSet)

    def test_verboseMMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a MMS'''
        self.verbosecheck(multiMesSet)

    def test_verboseTimeAvgMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a time-averaged MS'''
        self.verbosecheck(timeavg_ms)

    def test_verboseTimeAvgMMS(self):
        '''Listobs test: Check that a verbose file is larger than a non-verbose one for a time-averaged MMS'''
        self.verbosecheck(timeavg_mms)

    # Test the overwrite param

    def test_overwriteMS(self):
        """Test overwrite parameter - CAS-5203: test for MS"""
        self.overwritecheck(mesSet)

    def test_overwriteMMS(self):
        """Test overwrite parameter - CAS-5203: test for MMS"""
        self.overwritecheck(multiMesSet)

    def test_overwriteTimeAvgMS(self):
        """Test overwrite parameter - CAS-5203: test for time-averaged MS"""
        self.overwritecheck(timeavg_ms)

    def test_overwriteTimeAvgMMS(self):
        """Test overwrite parameter - CAS-5203: test for time-averaged MMS"""
        self.overwritecheck(timeavg_mms)

    # Test the inf loop bug CAS-6733

    def test_CAS_6733(self):
        """Verify listobs runs to completion on data set in CAS-6733. This was an infinite loop bugfix"""
        if CASA6:
            vis = casatools.ctsys.resolve('visibilities/evla/CAS-6733.ms')

        elif os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
            vis = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/evla/CAS-6733.ms'
        else:
            vis = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/evla/CAS-6733.ms'
            
        self.assertTrue(listobs(vis=vis))

    # Test average interval

    def test_avgIntervalMS(self):
        '''Listobs test: Check that the Int (s) column exists for a MS'''
        self.avgIntervalcheck(mesSet)

    def test_avgIntervalMMS(self):
        '''Listobs test: Check that the Int (s) column exists for a MMS'''
        self.avgIntervalcheck(multiMesSet)

    def test_avgIntervalTimeAvgMS(self):
        '''Listobs test: Check that the Int (s) column exists for a time-averaged MS'''
        self.avgIntervalcheck(timeavg_ms)

    def test_avgIntervalTimeAvgMMS(self):
        '''Listobs test: Check that the Int (s) column exists for a time-averaged MMS'''
        self.avgIntervalcheck(timeavg_mms)

    # Test list unflagged parameter

    def test_listunflMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a MS'''
        self.unfcheck(mesSet)

    def test_listunflMMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a MMS'''
        self.unfcheck(multiMesSet)

    def test_listunflTimeAvgMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a time-averaged MS'''
        self.unfcheck(timeavg_ms)

    def test_listunflTimeAvgMMS(self):
        '''Listobs test: Check that the list unflagged column shows up in a time-averaged MMS'''
        self.unfcheck(timeavg_ms)

def suite():
    return[test_listobs]

if __name__ == '__main__':
    unittest.main()
