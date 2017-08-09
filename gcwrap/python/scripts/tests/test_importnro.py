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

from importnro import importnro

myms = gentools(['ms'])[0]


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

def suite():
    return [importnro_test]
