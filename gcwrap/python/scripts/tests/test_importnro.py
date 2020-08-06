import os
import sys
import shutil
import inspect
import re
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import sha
import time
import numpy
import itertools

from importnro import importnro

myms = gentools(['ms'])[0]


# Utilities
def get_antenna_position(vis, row):
    (mytb,) = gentools(['tb'])
    antenna_table = os.path.join(vis, 'ANTENNA')
    
    mytb.open(antenna_table)
    try:
        pos = mytb.getcell('POSITION', row)
        poskey = mytb.getcolkeywords('POSITION')
    finally:
        mytb.close()
        
    posref = poskey['MEASINFO']['Ref']
    qpos = [qa.quantity(v,u) for v,u in zip(pos, poskey['QuantumUnits'])]
    mantpos = me.position(rf=posref, v0=qpos[0], v1=qpos[1], v2=qpos[2])
    
    return mantpos

def get_valid_pointing_info(vis):
    (mytb,) = gentools(['tb'])
    pointing_table = os.path.join(vis, 'POINTING')
    
    mytb.open(pointing_table)
    try:
        timekey = mytb.getcolkeywords('TIME')
        dirkey = mytb.getcolkeywords('DIRECTION')
        nrow = mytb.nrows()

        irow = 0
        pdir = mytb.getcell('DIRECTION', irow)
        ptime = mytb.getcell('TIME', irow)

        while numpy.all(pdir == 0.0) and irow < nrow:
            irow += 1
            pdir = mytb.getcell('DIRECTION', irow)
            ptime = mytb.getcell('TIME', irow)
    finally:
        mytb.close()
        
    dirref = dirkey['MEASINFO']['Ref']
    qdir = [qa.quantity(v,u) for v,u in zip(pdir[:,0], dirkey['QuantumUnits'])]
    mpdir = me.direction(rf=dirref, v0=qdir[0], v1=qdir[1])

    timeref = timekey['MEASINFO']['Ref']
    qtime = qa.quantity(ptime, timekey['QuantumUnits'][0])
    mepoch = me.epoch(rf=timeref, v0=qtime)
    
    return mepoch, mpdir
        

class importnro_test(unittest.TestCase):
    """
       test_overwrite -- File existence check
       test_invaliddata -- Invalid data check
       test_normal -- Normal data import
    """
    # Input and output names
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/importnro/'
    infile='orixa.OrionKL.20151209212931.16.Y'
    prefix='importnro_test'
    outfile=prefix+'.ms'

    def setUp(self):
        self.res=None
        if (not os.path.exists(self.infile)):
            shutil.copy(self.datapath+self.infile, self.infile)

        default(importnro)

    def tearDown(self):
        if (os.path.exists(self.infile)):
            os.remove(self.infile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test_overwrite(self):
        """test_overwrite: File existence check"""
        shutil.copy(self.infile, self.outfile)
        with self.assertRaisesRegex(RuntimeError, '.* exists\.$') as cm:
            importnro(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_invaliddata(self):
        """test_invaliddata: Invalid data check"""
        with open(self.infile, 'wb') as f: f.write('AA')
        #os.remove(os.path.join(self.infile, 'table.info'))
        with self.assertRaisesRegex(RuntimeError, '.* is not a valid NOSTAR data\.$') as cm:
            importnro(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_normal(self):
        """test_normal: Normal data import"""
        ret = importnro(infile=self.infile, outputvis=self.outfile, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        try:
            # to check if outfile is valid MS
            myms.open(self.outfile)
            myms.close()
            
        except Exception as e:
            print(e)
            self.fail('outputvis is not a valid ms')
        
        # check weight initialization
        self._check_weights(self.outfile)
        
        # check subtables
        self._check_optional_subtables(self.outfile)
        
        # check SOURCE INTERVAL (CAS-11442)
        self._check_source_interval(self.outfile)
        
        # check WEATHER table
        self._check_weather(self.outfile)
        
    def _check_weights(self, vis):
        _tb = gentools(['tb'])[0]
        take_diff = lambda actual, expected: numpy.abs((actual - expected) / expected)
        tolerance = 1.0e-7
        try:
            _tb.open(os.path.join(vis, 'DATA_DESCRIPTION'))
            spwids = _tb.getcol('SPECTRAL_WINDOW_ID')
            _tb.close()
            
            _tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
            nrow = _tb.nrows()
            g = (numpy.mean(_tb.getcell('EFFECTIVE_BW', irow)) for irow in range(nrow))
            effbws = numpy.fromiter(g, dtype=float)
            _tb.close()
            
            _tb.open(vis)
            nrow = _tb.nrows()
            for irow in range(nrow):
                weight = _tb.getcell('WEIGHT', irow)
                sigma = _tb.getcell('SIGMA', irow)
                interval = _tb.getcell('INTERVAL', irow)
                ddid = _tb.getcell('DATA_DESC_ID', irow)
                spwid = spwids[ddid]
                effbw = effbws[spwid]
                weight_expected = interval * effbw
                sigma_expected = 1.0 / numpy.sqrt(weight_expected)
                #print irow, 'weight', weight, 'sigma', sigma, 'expected', weight_expected, ' ', sigma_expected
                weight_diff = take_diff(weight, weight_expected)
                sigma_diff = take_diff(sigma, sigma_expected)
                #print irow, 'weight_diff', weight_diff, 'sigma_diff', sigma_diff
                self.assertTrue(all(weight_diff < tolerance), msg='Row %s: weight verification failed'%(irow))
                self.assertTrue(all(sigma_diff < tolerance), msg='Row %s: sigma verification failed'%(irow))
            _tb.close()
        finally:
            _tb.close()
            
    def _check_optional_subtables(self, vis):
        """Check if optional subtables are valid"""
        self._check_NRO_ARRAY(vis)
        
    def _check_NRO_ARRAY(self, vis):
        """Check if NRO_ARRAY table is valid"""
        table_name = 'NRO_ARRAY'
        (mytb,) = gentools(['tb'])
        mytb.open(vis)
        try:
            self.assertTrue(table_name in mytb.keywordnames())
        finally:
            mytb.close()
            
        table_path = os.path.join(vis, table_name)
        self.assertTrue(os.path.exists(table_path))
        mytb.open(table_path)
        try:
            cols = set(['BEAM', 'POLARIZATION', 'SPECTRAL_WINDOW', 'ARRAY'])
            self.assertEqual(set(mytb.colnames()), cols)
            beam = mytb.getcol('BEAM')
            pol = mytb.getcol('POLARIZATION')
            spw = mytb.getcol('SPECTRAL_WINDOW')
            arr = mytb.getcol('ARRAY')
            nrow = mytb.nrows()
        finally:
            mytb.close()
            
        arr_expected = numpy.arange(nrow, dtype=int)
        self.assertTrue(numpy.all(arr_expected == arr))
        beam_expected = numpy.empty_like(arr_expected)
        beam_expected[0:4] = 0
        beam_expected[4:8] = 1
        beam_expected[8:12] = 2
        beam_expected[12:16] = 3
        beam_expected[16:] = -1
        self.assertTrue(numpy.all(beam_expected == beam))
        spw_expected = numpy.empty_like(arr_expected)
        spw_expected[0:16:2] = 0
        spw_expected[1:16:2] = 1
        spw_expected[16:] = -1
        self.assertTrue(numpy.all(spw_expected == spw))
        pol_expected = numpy.empty_like(arr_expected)
        pol_expected[0:16:4] = 12
        pol_expected[1:16:4] = 12
        pol_expected[2:16:4] = 9
        pol_expected[3:16:4] = 9
        pol_expected[16:] = -1
        self.assertTrue(numpy.all(pol_expected == pol))
        
    def _check_source_interval(self, vis):
        """Check if SOURCE INTERVAL is consistent with OBSERVATION TIME_RANGE"""
        (mytb,) = gentools(['tb'])
        
        source_table = os.path.join(vis, 'SOURCE')
        observation_table = os.path.join(vis, 'OBSERVATION')
        
        # read OBSERVATION.TIME_RANGE
        mytb.open(observation_table)
        try:
            time_range = mytb.getcell('TIME_RANGE', 0)
        finally:
            mytb.close()
            
        # read SOURCE.TIME and SOURCE.INTERVAL
        mytb.open(source_table)
        try:
            source_time = mytb.getcol('TIME')
            source_interval = mytb.getcol('INTERVAL')
        finally:
            mytb.close()
            
        for t, dt in zip(source_time, source_interval):
            source_time_range = numpy.asarray([t-dt/2, t+dt/2])
            diff = numpy.abs((source_time_range - time_range) / time_range)
            #print 'diff={}'.format(diff)
            self.assertTrue(numpy.all(diff < 1.0e-16))
            
    def test_timestamp(self):
        """test_timestamp: Check if timestamp is properly converted to UTC"""
        ret = importnro(infile=self.infile, outputvis=self.outfile, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        
        # test if telescope elevation has reasonable value
        (myme,) = gentools(['me'])
        # make sure me tool is initialized
        myme.done()
        
        # antenna_position should be a position measure
        antenna_position = get_antenna_position(self.outfile, 0)
        self.assertTrue(myme.ismeasure(antenna_position))
        self.assertTrue(antenna_position['type'] == 'position')
        
        # pointing_time should be a time (epoch) measure
        # pointing_direction should be a direction measure
        # pointing_direction should not be [0,0]
        pointing_time, pointing_direction = get_valid_pointing_info(self.outfile)
        self.assertTrue(myme.ismeasure(pointing_time))
        self.assertTrue(pointing_time['type'] == 'epoch')
        self.assertTrue(myme.ismeasure(pointing_direction))
        self.assertTrue(pointing_direction['type'] == 'direction')
        self.assertFalse(pointing_direction['m0']['value'] == 0.0 and pointing_direction['m1']['value'])
        
        # convert pointing_direction (J2000) to AZELGEO
        # frame configuration
        myme.doframe(pointing_time)
        myme.doframe(antenna_position)
        
        # frame cnversion
        azel = myme.measure(v=pointing_direction, rf='AZELGEO')
        myme.done()
        
        # check if elevation is in range [0deg, 90deg]
        elevation = qa.convert(azel['m1'], 'deg')
        msg = 'Timestamp used for the conversion could be wrong.: calculated elevation={value}{unit}'.format(**elevation)
        self.assertLessEqual(elevation['value'], 90.0, msg='Elevation is above the upper limit (> 90deg). {}'.format(msg))
        self.assertGreaterEqual(elevation['value'], 0.0, msg='Elevation is below the lower limit (< 0deg). {}'.format(msg))

    def _check_weather(self, vis):
        # check PRESSURE
        self._check_weather_column(vis, 'PRESSURE', 'hPa', 400.0, 1100.0)
        
        # check TEMPERATURE
        self._check_weather_column(vis, 'TEMPERATURE', 'K', 243.0, 313.0)
        
    def _check_weather_column(self, vis, colname, unit, value_min, value_max):
        weather_table = os.path.join(vis, 'WEATHER')
        tb.open(weather_table)
        try:
            # column should exist
            self.assertTrue(colname in tb.colnames())
            
            # unit check
            colkeys = tb.getcolkeywords(colname)
            self.assertTrue('QuantumUnits' in colkeys)
            column_unit = colkeys['QuantumUnits'][0]
            print(('{0} unit is {1}'.format(colname, column_unit)))
            self.assertEqual(column_unit, unit)
            
            # value should be in reasonable range
            column_value = tb.getcol(colname)
            self.assertTrue(numpy.all(value_min < column_value))
            self.assertTrue(numpy.all(column_value < value_max))
        finally:
            tb.close()
        

def suite():
    return [importnro_test]
