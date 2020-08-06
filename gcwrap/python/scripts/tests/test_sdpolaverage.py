import shutil
import unittest
import os
import numpy
import math
import sys
import exceptions
import filecmp
import glob
from tasks import sdpolaverage
from taskinit import mstool, tbtool, msmdtool, aftool
from __main__ import default
import testhelper as th
from sdutil import tbmanager, toolmanager, table_selector

# Define the root for the data files
datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/tsdfit/"

aflocal = aftool()

def weighToSigma(weight):
    if weight > sys.float_info.min:
        return 1.0/math.sqrt(weight)
    else:
        return -1.0

def sigmaToWeight(sigma):
    if sigma > sys.float_info.min:
        return 1.0/math.pow(sigma,2)
    else:
        return 0.0


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
                raise ValueError('!=')
        except ValueError:
            errmsg = "%r != %r" % (val, expval)
            if (len(errmsg) > 66): # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError(errmsg)
        except Exception as e:
            print("Error comparing", val, "to", expval)
            raise e

class test_sdpolaverage(unittest.TestCase):
    def setUp(self):
        self.inputms  = "analytic_type1.fit.ms"
        self.outputms = "polave.ms"
        datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/tsdfit/"
        os.system('cp -RL '+datapath + self.inputms +' '+ self.inputms)
        default(sdpolaverage)

    def tearDown(self):
        os.system('rm -rf ' + self.inputms)
        os.system('rm -rf ' + self.outputms)

    def test_default(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms, datacolumn='float_data')
        with tbmanager(self.inputms) as tb:
            indata = tb.getcell('FLOAT_DATA', 0)
        with tbmanager(self.outputms) as tb:
            outdata = tb.getcell('FLOAT_DATA', 0)
        
        self.assertEqual(len(indata), len(outdata), 'Input and output data have different shape.')
        for i in range(len(indata)):
            for j in range(len(indata[0])):
                self.assertEqual(indata[i][j], outdata[i][j], 'Input and output data unidentical.')

    def test_stokes_float_data(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms, polaverage='stokes', datacolumn='float_data')
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

    def test_stokes_corrected_data(self):
        sdpolaverage(infile=self.inputms, outfile=self.outputms, polaverage='stokes', datacolumn='corrected')
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
    return [test_sdpolaverage]
