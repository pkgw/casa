import shutil
import unittest
import os
import numpy
import math
import sys
import exceptions
import filecmp
import glob
from tasks import nrobeamaverage
from taskinit import mstool, tbtool
from __main__ import default
import testhelper as th
from sdutil import tbmanager, toolmanager, table_selector

# Define the root for the data files
datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/nrobeamaverage/"

def check_eq(val, expval, tol=None):
    """Checks that val matches expval within tol."""
#    print val
    if type(val) == dict:
        for k in val:
            check_eq(val[k], expval[k], tol)
    else:
        try:
            if tol and hasattr(val, '__rsub__'):
                are_eq = abs(val - expval) < tol
            else:
                are_eq = val == expval
            if hasattr(are_eq, 'all'):
                are_eq = are_eq.all()
            if not are_eq:
                raise ValueError, '!='
        except ValueError:
            errmsg = "%r != %r" % (val, expval)
            if (len(errmsg) > 66): # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError, errmsg
        except Exception, e:
            print "Error comparing", val, "to", expval
            raise e

class test_nrobeamaverage(unittest.TestCase):
    def setUp(self):
        self.inputms  = "onon.ms"
        self.outputms = "bave.ms"
        os.system('cp -RL '+ datapath + self.inputms +' '+ self.inputms)
        default(nrobeamaverage)

    def tearDown(self):
        os.system('rm -rf ' + self.inputms)
        os.system('rm -rf ' + self.outputms)

    def test_default(self):
        # timebin='0s' : no time averaging, only rewriting beam IDs
        nrobeamaverage(infile=self.inputms, outfile=self.outputms)
        with tbmanager(self.inputms) as tb:
            indata = tb.getcol('FLOAT_DATA')
        with tbmanager(self.outputms) as tb:
            outdata = tb.getcol('FLOAT_DATA')
        """
        self.assertEqual(len(indata), len(outdata), 'Input and output data have different shape.')
        for i in range(len(indata)): #pol
            for j in range(len(indata[i])): #chan
                for k in range(len(indata[i][j])): #row
                    self.assertEqual(indata[i][j][k], outdata[i][j][k], 'Input and output data unidentical.[%d,%d,%d]' % (i, j, k))
        """

    def atest_stokes_float_data(self):
        nrobeamaverage(infile=self.inputms, outfile=self.outputms, polaverage='stokes', datacolumn='float_data')
        #check data
        with tbmanager(self.inputms) as tb:
            indata = tb.getcell('FLOAT_DATA', 0)
        with tbmanager(self.outputms) as tb:
            outdata = tb.getcell('FLOAT_DATA', 0)
        
        self.assertEqual(len(outdata), 1, 'No averaging over polarization?')
        tol = 1e-5
        for i in range(len(indata[0])):
            mean = 0.5 * (indata[0][i] + indata[1][i])
            check_eq(outdata[0][i], mean, tol)

        #check polarization id (should be 1)
        with tbmanager(self.outputms) as tb:
            outddesc = tb.getcell('DATA_DESC_ID', 0)
        with tbmanager(self.outputms + '/DATA_DESCRIPTION') as tb:
            outpolid = tb.getcol('POLARIZATION_ID')
        with tbmanager(self.outputms + '/POLARIZATION') as tb:
            outpoltype = tb.getcell('CORR_TYPE', outpolid[outddesc])

        self.assertEqual(len(outpoltype), 1, 'Polarization id is inconsistent with data.')
        self.assertEqual(outpoltype[0], 1, 'Has wrong polarization id.')

    def atest_stokes_corrected_data(self):
        nrobeamaverage(infile=self.inputms, outfile=self.outputms, polaverage='stokes', datacolumn='corrected')
        #check data
        with tbmanager(self.inputms) as tb:
            indata = tb.getcell('CORRECTED_DATA', 0)
        with tbmanager(self.outputms) as tb:
            outdata = tb.getcell('DATA', 0)
        
        self.assertEqual(len(outdata), 1, 'No averaging over polarization?')
        tol = 1e-5
        for i in range(len(indata[0])):
            mean = 0.5 * (indata[0][i] + indata[1][i])
            check_eq(outdata[0][i].real, mean.real, tol)
            check_eq(outdata[0][i].imag, mean.imag, tol)

        #check polarization id (should be 1)
        with tbmanager(self.outputms) as tb:
            outddesc = tb.getcell('DATA_DESC_ID', 0)
        with tbmanager(self.outputms + '/DATA_DESCRIPTION') as tb:
            outpolid = tb.getcol('POLARIZATION_ID')
        with tbmanager(self.outputms + '/POLARIZATION') as tb:
            outpoltype = tb.getcell('CORR_TYPE', outpolid[outddesc])

        self.assertEqual(len(outpoltype), 1, 'Polarization id is inconsistent with data.')
        self.assertEqual(outpoltype[0], 1, 'Has wrong polarization id.')


def suite():
    return [test_nrobeamaverage]
