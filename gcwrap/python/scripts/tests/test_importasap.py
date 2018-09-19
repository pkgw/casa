import os
import sys
import shutil
import inspect
import re
from __main__ import default
from tasks import *
from taskinit import *
from casac import casac
import unittest
import sha
import time
import numpy

from importasap import importasap

myms = gentools(['ms'])[0]


class importasap_test(unittest.TestCase):
    """
       test_overwrite -- File existence check
       test_invaliddata -- Invalid data check
       test_normal -- Normal data import
       test_flagversions -- Check if existing flagversions file is overwritten
       test_noflagversions -- Do not create flagversions file
    """
    # Input and output names
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/importasap/'
    infile='uid___A002_X85c183_X36f.test.asap'
    prefix='importasap_test'
    outfile=prefix+'.ms'

    def setUp(self):
        self.res=None
        if (not os.path.exists(self.infile)):
            shutil.copytree(self.datapath+self.infile, self.infile)

        default(importasap)

    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test_overwrite(self):
        """test_overwrite: File existence check"""
        shutil.copytree(self.infile, self.outfile)
        with self.assertRaisesRegex(RuntimeError, '.* exists\.$') as cm:
            importasap(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_invaliddata(self):
        """test_invaliddata: Invalid data check"""
        os.remove(os.path.join(self.infile, 'table.info'))
        with self.assertRaisesRegex(RuntimeError, '.* is not a valid Scantable\.$') as cm:
            importasap(infile=self.infile, outputvis=self.outfile, overwrite=False)
    
    def test_normal(self):
        """test_normal: Normal data import"""
        ret = importasap(infile=self.infile, outputvis=self.outfile,
                         flagbackup=True, overwrite=True)
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

        # check flagversions
        self._check_flagversions(self.outfile)
        
        # check pressure unit and its value
        self._check_atm_pressure(self.outfile)

    def test_flagversions(self):
        """test_flagversions -- Check if existing flagversions file is overwritten"""
        # create flagversions file
        ret = importasap(infile=self.infile, outputvis=self.outfile,
                         flagbackup=True, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        self._check_flagversions(self.outfile)

        # save state with different name
        version_name = 'OverwriteTest'
        aflocal = casac.agentflagger()
        aflocal.open(self.outfile)
        aflocal.saveflagversion(versionname=version_name,
                                comment='flag state for testing',
                                merge='save')
        aflocal.done()

        # then, create flagversions file
        ret = importasap(infile=self.infile, outputvis=self.outfile,
                         flagbackup=True, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        self._check_flagversions(self.outfile)

        # verification
        aflocal = casac.agentflagger()
        aflocal.open(self.outfile)
        versions_list = aflocal.getflagversionlist()
        self.assertTrue(all([re.match('^%s : .*$'%(version_name), x) is None for x in versions_list]))
        aflocal.done()

    def test_noflagversions(self):
        """test_noflagversions -- Do not create flagversions file"""
        flagversions = self._flagversions(self.infile)
        
        # create flagversions file
        ret = importasap(infile=self.infile, outputvis=self.outfile,
                         flagbackup=True, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        self._check_flagversions(self.outfile)

        # run importasap
        # existing flagversions file must be deleted
        ret = importasap(infile=self.infile, outputvis=self.outfile,
                         flagbackup=False, overwrite=True)
        self.assertTrue(os.path.exists(self.outfile))
        self.assertFalse(os.path.exists(flagversions))
    
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

    def _check_flagversions(self, vis):
        flagversions = self._flagversions(vis)
        self.assertTrue(os.path.exists(flagversions))

        # keep current flag state
        tb.open(vis)
        nrow = tb.nrows()
        flag_row_org = tb.getcol('FLAG_ROW')
        flag_org = tb.getvarcol('FLAG')
        tb.close()
        
        # flag version named 'Original' should be created
        # its content should match with current flag status
        version_name = 'Original'
        aflocal = casac.agentflagger()
        aflocal.open(vis)
        versions_list = aflocal.getflagversionlist()
        self.assertTrue(any([re.match('^%s : .*$'%(version_name), x) is not None for x in versions_list]))

        # restore flag state
        aflocal.restoreflagversion(versionname=version_name, merge='replace')
        aflocal.done()

        # get restored flag state
        tb.open(vis)
        flag_row = tb.getcol('FLAG_ROW')
        flag = tb.getvarcol('FLAG')
        tb.close()

        # verification
        for irow in range(nrow):
            self.assertEqual(flag_row_org[irow], flag_row[irow], msg='row %s: FLAG_ROW is different'%(irow))

        for (row, val_org) in list(flag_org.items()):
            val = flag[row]
            self.assertTrue(all((val_org == val).flatten()), msg='row %s: FLAG is different'%(row))

    def _flagversions(self, vis):
        return vis.rstrip('/') + '.flagversions'
    
    def _check_atm_pressure(self, vis):
        weather_table = os.path.join(vis, 'WEATHER')
        tb.open(weather_table)
        try:
            # PRESSURE column should exist
            self.assertTrue('PRESSURE' in tb.colnames())
            
            # unit should be hPa
            colkeys = tb.getcolkeywords('PRESSURE')
            self.assertTrue('QuantumUnits' in colkeys)
            pressure_unit = colkeys['QuantumUnits'][0]
            print(('Pressure unit is {0}'.format(pressure_unit)))
            self.assertEqual(pressure_unit, 'hPa')
            
            # value should be in reasonable range
            pressure_min = 400.0
            pressure_max = 1100.0
            pressure_value = tb.getcol('PRESSURE')
            self.assertTrue(numpy.all(pressure_min < pressure_value))
            self.assertTrue(numpy.all(pressure_value < pressure_max))
        finally:
            tb.close()


def suite():
    return [importasap_test]
