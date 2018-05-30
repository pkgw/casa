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
    qpos = [qa.quantity(v,u) for v,u in itertools.izip(pos, poskey['QuantumUnits'])]
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
    qdir = [qa.quantity(v,u) for v,u in itertools.izip(pdir[:,0], dirkey['QuantumUnits'])]
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
        with self.assertRaisesRegexp(RuntimeError, '.* exists\.$') as cm:
            importnro(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_invaliddata(self):
        """test_invaliddata: Invalid data check"""
        with open(self.infile, 'wb') as f: f.write('AA')
        #os.remove(os.path.join(self.infile, 'table.info'))
        with self.assertRaisesRegexp(RuntimeError, '.* is not a valid NOSTAR data\.$') as cm:
            importnro(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_normal(self):
        """test_normal: Normal data import"""
        ret = importnro(infile=self.infile, outputvis=self.outfile, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        try:
            # to check if outfile is valid MS
            myms.open(self.outfile)
            myms.close()
            
        except Exception, e:
            print e
            self.fail('outputvis is not a valid ms')
        
        # check weight initialization
        self._check_weights(self.outfile)
        
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
            g = (numpy.mean(_tb.getcell('EFFECTIVE_BW', irow)) for irow in xrange(nrow))
            effbws = numpy.fromiter(g, dtype=float)
            _tb.close()
            
            _tb.open(vis)
            nrow = _tb.nrows()
            for irow in xrange(nrow):
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


def suite():
    return [importnro_test]
