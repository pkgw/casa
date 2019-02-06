##########################################################################
# test_req_task_visstat.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_visstat/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import visstat
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

### Data ###
if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/other/outlier_ut.ms/')
    mms_data = casatools.ctsys.resolve('visibilities/other/outlier_mms.mms/')
    selectiondata = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
    mms_select = casatools.ctsys.resolve('visibilities/alma/uid_mms.mms')
    singledish = casatools.ctsys.resolve('visibilities/other/analytic_spectra_tsys.ms')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/other/outlier_ut.ms/'
        mms_data = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/other/outlier_mms.mms/'
        selectiondata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
        mms_select = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid_mms.mms'
        singledish = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/other/analytic_spectra_tsys.ms'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/other/outlier_ut.ms/'
        mms_data = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/other/outlier_mms.mms/'
        selectiondata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
        mms_select = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid_mms.mms'
        singledish = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/other/analytic_spectra_tsys.ms'
    
axislist = ['flag', 'antenna1', 'antenna2', 'feed1', 'feed2', 'field_id', 'array_id', 'data_desc_id', 'flag_row', 'interval', 'scan', 'scan_number', 'time', 'weight_spectrum', 'amp', 'amplitude', 'phase', 'real', 'imag', 'imaginary', 'uvrange']
 
keylist = ['firstquartile', 'isMasked', 'isWeighted', 'max', 'maxDatasetIndex', 'maxIndex', 'mean', 'medabsdevmed', 'median', 'min', 'minDatasetIndex', 'minIndex', 'npts', 'rms', 'stddev', 'sum', 'sumOfWeights', 'sumsq', 'thirdquartile', 'variance']

# NOTE: I also need to get and/or make a multimesset to test this data from. I could usr the data from listobs to do this.

nostat = visstat(selectiondata)
nostatmms = visstat(mms_select)

class visstat_test(unittest.TestCase):
     
    def setUp(self):
        if not CASA6:
            default(visstat)
   
    def test_axis(self):
        '''
            test_axis
            ------------------------------
            
            Test the axis parameter values.
            visstat should return a dict and the keys should match the provided key list.
            
            This test iterates over all the possible axis values
        '''
        for axis in axislist:

            axisstat = visstat(datapath, axis=axis)
            axismms = visstat(mms_data, axis=axis)
            
            self.assertTrue( type(axisstat) == type(dict()), msg='output is not a dict for axis {}'.format(axis) )
            self.assertTrue( type(axismms) == type(dict()), msg='output is not a dict for axis {} on mms'.format(axis) )
        
            self.assertTrue(sorted(list(axisstat[list(axisstat.keys())[0]].keys())) == keylist, msg='keys do not match the key list for axis {}'.format(axis))
            self.assertTrue(sorted(list(axismms[list(axismms.keys())[0]].keys())) == keylist, msg='keys do not match the key list for axis {} for mms'.format(axis))
            
    def test_reportingaxes(self):
        '''
            test_reportingaxes
            -----------------------------
            
            Test the reportingaxes parameter.
            The output should be a dict and contain all the expected keys.
            
            Iterate over all the possible values.
        '''
        
        for axes in ['ddid', 'field', 'integration']:
            
            reportstat = visstat(datapath, reportingaxes=axes)
            reportmms = visstat(mms_data, reportingaxes=axes)
            
            self.assertTrue( type(reportstat) == type(dict()) )
            self.assertTrue( type(reportmms) == type(dict()) )
            
            self.assertTrue(sorted(list(reportstat[list(reportstat.keys())[0]].keys())) == keylist)
            self.assertTrue(sorted(list(reportmms[list(reportmms.keys())[0]].keys())) == keylist)
            
    def test_useflags(self):
        '''
            test_useflags
            ----------------------
            
            Check that the useflags parameter produces different results then when useflags = False.
        '''
        
        withflags = visstat(datapath, useflags=True)
        withoutflags = visstat(datapath, useflags = False)
        
        flagsmms = visstat(mms_select, useflags=False)
        noflagmms = visstat(mms_select, useflags=True)
        
        self.assertTrue(withflags != withoutflags)
        self.assertTrue(flagsmms != noflagmms)
        
    def test_datacolumn(self):
        '''
            test_datacolumn
            ----------------------------
            
            Check the data column parameter.
            
            Iterate over possible data column inputs and check that a dictionary is created by visstat
            also check that all the keys that should be present are there.
            
            (This last step may not be nessisary for this test)
        '''
        
        for col in [ 'data', 'corrected', 'model' ]:
            colstat = visstat(datapath, datacolumn=col)
            colmms = visstat(mms_data, datacolumn=col)
            
            self.assertTrue( type(colstat) == type(dict()), msg = 'Fails for column: {}'.format(col) )
            self.assertTrue( type(colmms) == type(dict()), msg = 'Fails for column: {} on mms'.format(col) )
            
            self.assertTrue(sorted(list(colstat[list(colstat.keys())[0]].keys())) == keylist, msg = 'Fails for column: {}'.format(col))
            self.assertTrue(sorted(list(colmms[list(colmms.keys())[0]].keys())) == keylist, msg = 'Fails for column: {} on mms'.format(col))
            
        floatstat = visstat(singledish, datacolumn='float_data')
        self.assertTrue( type(floatstat) == type(dict()) )
            
    def test_spw(self):
        '''
            test_spw
            ---------------
            
            Test the spectral window selection parameter.
            
            Assert that a selection using the spw parameter returns a different result than no selection.
        '''
        
        spwstat = visstat(selectiondata, spw='0')
        
        spwmms = visstat(mms_select, spw='0')
        
        self.assertTrue(spwstat != nostat)
        self.assertTrue(spwmms != nostatmms)
        
    def test_field(self):
        '''
            test_field
            --------------
            
            Test the field selection parameter.
            
            Assert that a selection using the field parameter returns a different result than no selection.
        '''
        
        fieldstat = visstat(selectiondata, field='0')
        
        fieldmms = visstat(mms_select, field='0')
        
        self.assertTrue( fieldstat != nostat )
        self.assertTrue( fieldmms != nostatmms )
        
    def test_selectdata(self):
        '''
            test_selectdata
            -----------------------
            
            Test the selectdata parameter
            
            Assert that the select data parameter prevents other selection fields from having an affect
            
            Assert that with selectdata=False and an active selection produces the same results as the task with no selection
        '''
        
        selectstat = visstat(selectiondata, selectdata=False, antenna='DV01')
        
        selectmms = visstat(mms_select, selectdata=False, antenna='DV01')
        
        self.assertTrue( selectstat == nostat )
        self.assertTrue( selectmms == nostatmms )
        
    def test_antenna(self):
        '''
            test_antenna
            ---------------------------
            
            Test the antenna selection parameter
            
            Assert that selection with this parameter will return a different result than no selection.
        '''
        
        antennastat = visstat(selectiondata, antenna='DV01')
        
        antennamms = visstat(mms_select, antenna='DV01')
        
        self.assertTrue( antennastat != nostat )
        self.assertTrue( antennamms != nostatmms )
        
    def test_uvange(self):
        '''
            test_uvrange
            --------------------
            
            Test the uvrange selection parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        uvRangeStat = visstat(selectiondata, uvrange='0~10')
        
        uvRangemms = visstat(mms_select, uvrange='0~10')
        
        self.assertTrue( uvRangeStat != nostat )
        self.assertTrue( uvRangemms != nostatmms)
        
    def test_timerange(self):
        '''
            test_timerange
            ------------------
            
            Test the timerange selection parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        timerangeStat = visstat(selectiondata, timerange='03:01:30~03:05:00')
        
        timerangemms = visstat(mms_select, timerange='03:01:30~03:05:00')
        
        self.assertTrue( timerangeStat != nostat )
        self.assertTrue( timerangemms != nostatmms )
        
    def test_correlation(self):
        '''
            test_correlation
            -------------------
            
            Test the correlation parameter
            
            Assert that the selection with this parameter will return a different value than no selection.
        '''
        
        corrStat = visstat(selectiondata, correlation='XX')
        
        corrmms = visstat(mms_select, correlation='XX')
        
        self.assertTrue( corrStat != nostat )
        self.assertTrue( corrmms != nostatmms )
        
    def test_scan(self):
        '''
            test_scan
            ------------
        
            Test the scan selection parameter
        
            Assert that the selction with this parameter will return a different result than no selection.
        '''
    
        scanStat = visstat(selectiondata, scan='1')
        
        scanmms = visstat(mms_select, scan='1')
        
        self.assertTrue( scanStat != nostat )
        self.assertTrue( scanmms != nostatmms )
        
    def test_array(self):
        '''
            test_array
            -------------
            
            Test the array selection parameter.
            
            Assert that checking an out of range array returns a NoneType, and valid selections retrun a dictionary
        '''
        if CASA6:
            with self.assertRaises(RuntimeError):
                arrayFail = visstat(selectiondata, array='1')
            with self.assertRaises(RuntimeError):
                arrayFailmms = visstat(mms_select, array='1')
        else:
            arrayFail = visstat(selectiondata, array='1')
            arrayFailmms = visstat(mms_select, array='1')
            
            self.assertTrue( type(arrayFail) == type(None) )
            self.assertTrue( type(arrayFailmms) == type(None) )
        
        arrayPass = visstat(selectiondata, array='0')
        arrayPassmms = visstat(mms_select, array='0')
        
        self.assertTrue( type(arrayPass) == type(dict()) )
        self.assertTrue( type(arrayPassmms) == type(dict()) )
        
    def test_observation(self):
        '''
            test_observation
            -----------------------
            
            Test the observation selection parameter
            
            Assert that checking an out of range observation ID returns a NoneType, and valid selections return a dictionary
        '''
        if CASA6:
            with self.assertRaises(RuntimeError):
                observationFail = visstat(selectiondata, observation=1)
            with self.assertRaises(RuntimeError):
                observationFailmms = visstat(mms_select, observation=1)
        else:
            observationFail = visstat(selectiondata, observation=1)
            observationFailmms = visstat(mms_select, observation=1)
            self.assertTrue( type(observationFail) == type(None) )
            self.assertTrue( type(observationFailmms) == type(None) )
        
        observationPass = visstat(selectiondata, observation=0)
        observationPassmms = visstat(mms_select, observation=0)
    
        self.assertTrue( type(observationPass) == type(dict()) )
        self.assertTrue( type(observationPassmms) == type(dict()) )
        
    def test_timeavg(self):
        '''
            test_timeaverage
            ---------------------
            
            Test the timeaverage parameter
            
            Assert that the dict produced when timeaverage = True is different from the one produced when timeaverage=False
        '''
        
        timeavgTrue = visstat(selectiondata, timeaverage=True)
        timeavgFalse = visstat(selectiondata, timeaverage=False)
        
        timeavgTruemms = visstat(mms_select, timeaverage=True)
        timeavgFalsemms = visstat(mms_select, timeaverage=False)
        
        self.assertTrue( timeavgFalse != timeavgTrue )
        self.assertTrue( timeavgFalsemms != timeavgTruemms)
        
    def test_timebin(self):
        '''
            test_timebin
            -------------------
            
            Test the timebin parameter
            
            Assert that the result when given a bin width for averaging is different than when none is given.
        '''
        
        timebinSelect = visstat(selectiondata, timeaverage=True, timebin='10s')
        
        timebinmms = visstat(selectiondata, timeaverage=True, timebin='10s')
        
        self.assertTrue( timebinSelect != nostat )
        self.assertTrue( timebinmms != nostatmms )
        
    def test_timespan(self):
        '''
            test_timespan
            ------------------
            
            Test the timespan parameter
            
            Assert that all parameter settings give different results than the default output
        '''
        
        scanSelect = visstat(selectiondata, timeaverage=True, timespan='scan')
        stateSelect = visstat(selectiondata, timeaverage=True, timespan='state')
        bothSelect = visstat(selectiondata, timeaverage=True, timespan='scan, state')
        
        scanSelectmms = visstat(mms_select, timeaverage=True, timespan='scan')
        stateSelectmms = visstat(mms_select, timeaverage=True, timespan='state')
        bothSelectmms = visstat(mms_select, timeaverage=True, timespan='scan, state')
        
        for item in [scanSelect, stateSelect, bothSelect]:
            self.assertTrue( item != nostat )
            
        for item in [scanSelectmms, stateSelectmms, bothSelectmms]:
            self.assertTrue( item != nostatmms )
        
    def test_maxuvwdistance(self):
        '''
            test_maxuvwdistance
            -----------------------
            
            Test the maxuvwdistance parameter
            
            Assert that the output is a python dict. Once again this selection seems to not change the values that are returned
        '''
        
        uvwSelect = visstat(selectiondata, timeaverage=True, maxuvwdistance=10.0)
        
        uvwSelectmms = visstat(mms_select, timeaverage=True, maxuvwdistance=10.0)
        
        self.assertTrue( uvwSelect != nostat )
        self.assertTrue( uvwSelectmms != nostatmms )
        
    def test_intent(self):
        '''
            test_intent
            -----------------
            
            Test the intent parameter
            
            Assert that the specified selection creates a different dict than the default values.
        '''
        
        intentSelect = visstat(selectiondata, intent='CALIBRATE_AMPLI.ON_SOURCE')
        
        intentSelectmms = visstat(mms_select, intent='CALIBRATE_AMPLI.ON_SOURCE')
        
        self.assertTrue( intentSelect != nostat )
        self.assertTrue( intentSelectmms != nostatmms )
        
        
def suite():
    return[visstat_test]

if __name__ == '__main__':
    unittest.main()
