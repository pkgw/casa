import os
import sys
import shutil
import numpy
import math

from __main__ import default
from tasks import *
from taskinit import *
import unittest
#
import listing
from numpy import array

import asap as sd
from sdfitold import sdfitold

try:
    from . import selection_syntax
except:
    import tests.selection_syntax as selection_syntax

class sdfitold_unittest_base:
    """
    Base class for sdfitold unit test.
    """
    # Data path of input/output
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdfit/'

    def setUp(self):
        if hasattr(self, 'inputs'):
            for infile in self.inputs:
                #print 'copying %s...'%(infile)
                if os.path.exists(infile):
                    shutil.rmtree(infile)
                shutil.copytree(os.path.join(self.datapath, infile), infile)
        default(sdfitold)

    def tearDown(self):
        if hasattr(self, 'inputs'):
            for infile in self.inputs:
                #print 'removing %s...'%(infile)
                if os.path.exists(infile):
                    shutil.rmtree(infile)

        if hasattr(self, 'outputs'):
            for outfile in self.outputs:
                #print 'removing %s...'%(outfile)
                os.system('rm -rf %s'%(outfile))

    def _check_file( self, filename, fail=True ):
        """
        Check if filename exists.
        The test fails if fail=True. Otherwise the method returns
        a bool which indicates if exists or not.
        """
        exists = os.path.exists(filename)
        if fail:
            self.assertTrue(exists,"'%s' does not exists." % filename)
        return exists

    def read_result(self, outfile):
        # basic check
        #   - check if self.outfile exists
        #   - check if self.outfile is a regular file
        self.assertTrue(os.path.exists(outfile))
        self.assertTrue(os.path.isfile(outfile))

        result = {}
        with open(outfile, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                s = line.split()
                scanno = int(s[0])
                ifno = int(s[1])
                polno = int(s[2])
                peak = float(s[4])
                center = float(s[5])
                fwhm = float(s[6])
                key = (scanno, ifno, polno)
                #self.assertFalse(result.has_key(key))
                if key in result:
                    result[key].extend([peak, center, fwhm])
                else:
                    result[key] = [peak, center, fwhm]
        return result

    def _verify_outfile(self, expected, outfile=None):
        # dictionary for verification
        # key: tuple consisting of SCANNO, IFNO, and POLNO
        #     (SCANNO,IFNO,POLNO), e.g.,
        #
        #         (0,0,0)
        #
        # value: flattened list of 'peak', 'cent', and 'fwhm'.
        #     If there are multiple rows for one key, the list
        #     is extended, e.g.,
        #
        #         [peak0, cent0, fwhm0, peak1, cent1, fwhm1, ...]
        #
        if outfile is None:
            if hasattr(self, 'outfile'): outfile = self.outfile
            else: self.fail("outfile is not specified.")
        result = self.read_result(self.outfile)
        tol = 1.0e-4
        properties = ['peak', 'cent', 'fwhm']
        nprop = len(properties)
        for key in list(expected.keys()):
            self.assertTrue(key in result,
                            msg='result does\'t have key %s'%(list(key)))
            e = expected[key]
            r = result[key]
            ne = len(e)
            nr = len(r)
            self.assertEqual(nr, ne,
                             msg='%s: number of row mismatch (expected %s, actual %s)'%(list(key), ne/nprop, nr/nprop))
            for irow in range(len(e)):
                diff = abs((r[irow] - e[irow])/e[irow])
                self.assertLessEqual(diff, tol,
                                     msg='%s (%s): result differ (expected %s, actual %s)'%(properties[irow % nprop], irow, e[irow], r[irow]))


class sdfitold_test(sdfitold_unittest_base, unittest.TestCase):
    """
    Unit tests for task sdfitold. No interactive testing.

    The list of tests:
    testGaussian00 --- test fitting a broad Gaussian profile (centre=4000, fwhm=1000, ampl=10)
    testGaussian01 --- test fitting a broad Gaussian profile (centre= 500, fwhm=1000, ampl=10) : on spectral edge
    testGaussian02 --- test fitting a narrow Gaussian profile (centre=4000, fwhm=100, ampl=10)
    testGaussian03 --- test fitting a combination of broad and narrow Gaussian profiles
                       (cen1=3000,fwhm1=1000,ampl1=10,cen2=6000,fwhm2=100,ampl2=10) : separated
    testGaussian04 --- test fitting a combination of broad and narrow Gaussian profiles
                       (cen1=4000,fwhm1=1000,ampl1=10,cen2=4700,fwhm2=100,ampl2=10) : overlapped
    testLorentzian00 --- test fitting a broad Lorentzian profile (centre=4000, fwhm=1000, ampl=10)
    testLorentzian01 --- test fitting a broad Lorentzian profile (centre= 500, fwhm=1000, ampl=10) : on spectral edge
    testLorentzian02 --- test fitting a narrow Lorentzian profile (centre=4000, fwhm=100, ampl=10)
    testLorentzian03 --- test fitting a combination of broad and narrow Lorentzian profiles
                       (cen1=3000,fwhm1=1000,ampl1=10,cen2=6000,fwhm2=100,ampl2=10) : separated
    testLorentzian04 --- test fitting a combination of broad and narrow Lorentzian profiles
                       (cen1=4000,fwhm1=1000,ampl1=10,cen2=4700,fwhm2=100,ampl2=10) : overlapped

    Note: (1) the rms noise is 1.0 for all data.

    created 21/04/2011 by Wataru Kawasaki
    """
    # Input and output names
    infile_gaussian   = 'Artificial_Gaussian.asap'
    infile_lorentzian = 'Artificial_Lorentzian.asap'
    inputs = [infile_gaussian, infile_lorentzian]

    def testGaussian00(self):
        """Test Gaussian00: single broad profile """
        infile = self.infile_gaussian
        scan = '0'
        fitfunc = "gauss"
        fitmode = "list"
        nfit = 1

        res = sdfitold(infile=infile,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0]],
               'fwhm': [[1000.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3997.420166015625, 2.2180848121643066]]],
               'fwhm': [[[1006.1046142578125, 5.2231903076171875]]],
               'nfit': [1],
               'peak': [[[9.9329404830932617, 0.044658195227384567]]]}
        """

        self.checkResult(res, ans)

    def testGaussian01(self):
        """Test Gaussian01: single broad profile on spectral edge"""
        infile = self.infile_gaussian
        scan = '1'
        spw = '0:0~2000'
        fitfunc = "gauss"
        fitmode = "list"
        #maskline = [0,2000]
        nfit = 1

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[500.0]],
               'fwhm': [[1000.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[504.638427734375, 2.7173392772674561]]],
               'fwhm': [[[998.78643798828125, 7.1386871337890625]]],
               'nfit': [1],
               'peak': [[[10.030097961425781, 0.047238241881132126]]]}
        """

        self.checkResult(res, ans)

    def testGaussian02(self):
        """Test Gaussian02: single narrow profile """
        infile = self.infile_gaussian
        scan = '2'
        fitfunc = "gauss"
        fitmode = "list"
        nfit = 1

        res = sdfitold(infile=infile,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0]],
               'fwhm': [[100.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3999.159912109375, 0.68400073051452637]]],
               'fwhm': [[[98.87506103515625, 1.6106985807418823]]],
               'nfit': [1],
               'peak': [[[9.9385099411010742, 0.14021013677120209]]]}
        """

        self.checkResult(res, ans)

    def testGaussian03(self):
        """Test Gaussian03: broad/narrow combination : separated """
        infile = self.infile_gaussian
        scan = '3'
        spw = '0:2000~4000;5500~6500'
        fitfunc = "gauss"
        fitmode = "list"
        #maskline = [[2000,4000],[5500,6500]]
        nfit = [1,1]

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[3000.0, 6000.0]],
               'fwhm': [[1000.0, 100.0]],
               'peak': [[10.0, 10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[2996.004638671875, 2.2644386291503906],
                         [5999.11181640625, 0.70802927017211914]]],
               'fwhm': [[[1001.549560546875, 5.4809303283691406],
                         [99.259437561035156, 1.6672815084457397]]],
               'nfit': [2],
               'peak': [[[9.899937629699707, 0.04574853926897049],
                         [9.9107418060302734, 0.14416992664337158]]]}
        """

        self.checkResult(res, ans)

    def testGaussian04(self):
        """Test Gaussian04: broad/narrow combination : overlapped """
        infile = self.infile_gaussian
        scan = '4'
        spw = '0:3000~4400;4500~5000'
        fitfunc = "gauss"
        fitmode = "list"
        #maskline = [[3000,4400],[4500,5000]]
        nfit = [1,1]

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0, 4700.0]],
               'fwhm': [[1000.0, 100.0]],
               'peak': [[10.0, 10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[4001.522216796875, 2.6332762241363525],
                         [4699.75732421875, 0.6802678108215332]]],
               'fwhm': [[[999.63507080078125, 6.4683256149291992],
                         [97.721427917480469, 1.7482517957687378]]],
               'nfit': [2],
               'peak': [[[9.9929990768432617, 0.04641139879822731],
                         [10.233022689819336, 0.15014420449733734]]]}
        """

        self.checkResult(res, ans)

    def testLorentzian00(self):
        """Test Lorentzian00: single broad profile """
        infile = self.infile_lorentzian
        scan = '0'
        fitfunc = "lorentz"
        fitmode = "list"
        nfit = 1

        res = sdfitold(infile=infile,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0]],
               'fwhm': [[1000.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3997.696044921875, 2.5651662349700928]]],
               'fwhm': [[[1010.3181762695312, 7.2803301811218262]]],
               'nfit': [1],
               'peak': [[[9.9210958480834961, 0.05041566863656044]]]}
        """

        self.checkResult(res, ans)

    def testLorentzian01(self):
        """Test Lorentzian01: single broad profile on spectral edge"""
        infile = self.infile_lorentzian
        scan = '1'
        spw = '0:0~2000'
        fitfunc = "lorentz"
        fitmode = "list"
        #maskline = [0,2000]
        nfit = 1

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[500.0]],
               'fwhm': [[1000.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[500.99105834960938, 2.8661653995513916]]],
               'fwhm': [[[995.85455322265625, 9.5194911956787109]]],
               'nfit': [1],
               'peak': [[[10.041034698486328, 0.053434751927852631]]]}
        """

        self.checkResult(res, ans)

    def testLorentzian02(self):
        """Test Lorentzian02: single narrow profile """
        infile = self.infile_lorentzian
        scan = '2'
        fitfunc = "lorentz"
        fitmode = "list"
        nfit = 1

        res = sdfitold(infile=infile,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0]],
               'fwhm': [[100.0]],
               'peak': [[10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3999.230224609375, 0.79903918504714966]]],
               'fwhm': [[[102.48796081542969, 2.2600326538085938]]],
               'nfit': [1],
               'peak': [[[9.9708395004272461, 0.1554737389087677]]]}
        """

        self.checkResult(res, ans)

    def testLorentzian03(self):
        """Test Lorentzian03: broad/narrow combination : separated """
        infile = self.infile_lorentzian
        scan = '3'
        spw = '0:2000~4000;5500~6500'
        fitfunc = "lorentz"
        fitmode = "list"
        #maskline = [[2000,4000],[5500,6500]]
        nfit = [1,1]

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[3000.0, 6000.0]],
               'fwhm': [[1000.0, 100.0]],
               'peak': [[10.0, 10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3001.23876953125, 2.5231325626373291],
                         [5999.01953125, 0.82874661684036255]]],
               'fwhm': [[[990.19671630859375, 8.1528301239013672],
                         [102.10212707519531, 2.3521277904510498]]],
               'nfit': [2],
               'peak': [[[9.9958734512329102, 0.051685664802789688],
                         [9.6133279800415039, 0.1561257392168045]]]}
        """

        self.checkResult(res, ans)

    def testLorentzian04(self):
        """Test Lorentzian04: broad/narrow combination : overlapped """
        infile = self.infile_lorentzian
        scan = '4'
        spw = '0:3000~4400;4500~5000'
        fitfunc = "lorentz"
        fitmode = "list"
        #maskline = [[3000,4400],[4500,5000]]
        nfit = [1,1]

        res = sdfitold(infile=infile,spw=spw,scan=scan,fitfunc=fitfunc,fitmode=fitmode,nfit=nfit)
        self.assertNotEqual(res, None, msg="The task returned None. Fit failed.")

        ans = {'cent': [[4000.0, 4700.0]],
               'fwhm': [[1000.0, 100.0]],
               'peak': [[10.0, 10.0]]}

        """ the result (RHEL5 64bit)
        ref = {'cent': [[[3995.85693359375, 3.0016641616821289],
                         [4699.53271484375, 0.82658475637435913]]],
               'fwhm': [[[972.22833251953125, 10.149419784545898],
                         [104.71010589599609, 2.7239837646484375]]],
               'nfit': [2],
               'peak': [[[10.013784408569336, 0.053735069930553436],
                         [9.9273672103881836, 0.15813499689102173]]]}
        """

        self.checkResult(res, ans)

    def checkResult(self, result, answer):
        for key in ['cent', 'fwhm', 'peak']:
            for i in range(len(result[key][0])):
                val = result[key][0][i][0]
                err = result[key][0][i][1]
                ans = answer[key][0][i]

                #check if result is consistent with answer in 3-sigma level
                threshold = 3.0

                within_errorrange = (abs(ans - val) <= abs(err * threshold))
                self.assertTrue(within_errorrange)

class sdfitold_test_exceptions(sdfitold_unittest_base, unittest.TestCase):
    """
    test the case when sdfitold throws exception.
    """
    # Input and output names
    infile_gaussian   = 'Artificial_Gaussian.asap'
    inputs = [infile_gaussian]

    def testNoData(self):
        try:
            res = sdfitold(infile=self.infile_gaussian,spw='99')
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            #pos=str(e).find('Invalid spectral window selection. Selection contains no data.')
            pos=str(e).find('No valid spw.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

class sdfitold_selection_syntax(sdfitold_unittest_base, selection_syntax.SelectionSyntaxTest):

    # Data path of input/output
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/singledish/'

    # Input and output names
    infile_convolve = 'sd_analytic_type3-1.asap'
    infile_shift = 'sd_analytic_type4-1.asap'
    infile_duplicate = 'sd_analytic_type5-1.asap'
    prefix = 'sdfitold_selection_syntax'
    inputs = [infile_convolve, infile_shift, infile_duplicate]
    outputs = [prefix+'*']

    # line information
    # | row | line channel | intensity |
    # | 0   | 20           | 5         |
    # | 1   | 40           | 10        |
    # | 2   | 60           | 20        |
    # | 3   | 80           | 30        |
    line_location = [[20, 5.0],
                     [40, 10.0],
                     [60, 20.0],
                     [80, 30.0]]

    # reference values for baseline fit
    # (scan,beam,if,pol): [peak, center, fwhm]
    fit_ref_convolve = {(15,23,0): [[0.9394372701644897, 20.0, 5.0]],
                        (16,25,1): [[1.8788745403289795, 40.0, 5.0]],
                        (16,21,0): [[3.757749080657959, 60.0, 5.0]],
                        (17,23,1): [[5.636623859405518, 80.0, 5.0]]}
    fit_ref_shift = {(15,23,0): [[0.9394372701644897, 50.0, 5.0]],
                     (16,25,1): [[1.8788745403289795, 50.0, 5.0]],
                     (16,21,0): [[3.757749080657959, 50.0, 5.0]],
                     (17,23,1): [[5.636623859405518, 50.0, 5.0]]}
    fit_ref_duplicate = {(15,23,0): [[0.9394372701644897, 50.0, 5.0],
                                      [0.9394372701644897, 75.0, 5.0]],
                         (16,25,1): [[1.8788745403289795, 50.0, 5.0],
                                     [1.8788745403289795, 75.0, 5.0]],
                         (16,21,0): [[3.757749080657959, 50.0, 5.0],
                                     [3.757749080657959, 75.0, 5.0]],
                         (17,23,1): [[5.636623859405518, 50.0, 5.0],
                                     [5.636623859405518, 75.0, 5.0]]}
    fit_ref = {infile_convolve: fit_ref_convolve,
               infile_shift: fit_ref_shift,
               infile_duplicate: fit_ref_duplicate}

    # tolerance
    tol = 1.0e-7

    @property
    def task(self):
        return sdfitold

    @property
    def spw_channel_selection(self):
        return True

    def __test_result(self, infile, result_ret, result_out, rows):
        casalog.post('result=%s'%(result_ret))
        s = sd.scantable(infile, average=False)
        self.assertTrue(infile in self.fit_ref)
        fit_ref = self.fit_ref[infile]

        for irow in range(len(rows)):
            row = rows[irow]
            scanno = s.getscan(row)
            ifno = s.getif(row)
            polno = s.getpol(row)
            key = (scanno, ifno, polno)
            ref = fit_ref[key]

            # check nfit
            nfit = result_ret['nfit'][0]
            self.assertEqual(nfit, len(ref))
            for icomp in range(len(ref)):
                comp = ref[icomp]
                # check peak
                peak = result_ret['peak'][irow][icomp][0]
                diff = abs((peak - comp[0]) / comp[0])
                self.assertLess(diff, self.tol)
                # check center
                center = result_ret['cent'][irow][icomp][0]
                self.assertEqual(center, comp[1])
                # check fwhm
                fwhm = result_ret['fwhm'][irow][icomp][0]
                self.assertEqual(fwhm, comp[2])

        for (k,v) in list(result_out.items()):
            ref = fit_ref[k]

            self.assertEqual(len(v), 3 * len(ref))
            for icomp in range(len(ref)):
                offset = icomp * 3
                _ref = ref[icomp]
                # check peak
                diff = abs((v[offset] - _ref[0]) / _ref[0])
                self.assertLess(diff, self.tol)
                # check center
                self.assertEqual(v[offset+1], _ref[1])
                # check fwhm
                self.assertEqual(v[offset+2], _ref[2])

    def __exec_complex_test(self, infile, params, exprs, rows, regular_test=True):
        num_param = len(params)
        test_name = self._get_test_name(regular_test)
        outfile = '.'.join([self.prefix, test_name])
        #print 'outfile=%s'%(outfile)
        casalog.post('%s: %s'%(test_name, ','.join(['%s = \'%s\''%(params[i],exprs[i]) for i in range(num_param)])))
        nfit = [1,1] if infile == self.infile_duplicate else [1]
        kwargs = {'infile': infile,
                  'nfit': nfit,
                  'fitfunc': 'gauss',
                  'fitmode': 'list',
                  'outfile': outfile,
                  'overwrite': True}

        for i in range(num_param):
            kwargs[params[i]] = exprs[i]

        regular_test = False
        if regular_test:
            result = self.run_task(**kwargs)
        else:
            result = sdfitold(**kwargs)

        # read outfile
        result_out = self.read_result(outfile)

        self.__test_result(infile, result, result_out, rows)

        return outfile

    def __exec_simple_test(self, infile, param, expr, rows, regular_test=True):
        return self.__exec_complex_test(infile, [param], [expr],
                                        rows, regular_test)

    ### field selection syntax test ###
    def test_field_value_default(self):
        """test_field_value_default: Test default value for field"""
        infile = self.infile_shift
        field = ''
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_exact(self):
        """test_field_id_exact: Test field selection by id"""
        infile = self.infile_convolve
        field = '5'
        spw = '23:0~40'
        rows = [0]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_lt(self):
        """test_field_id_lt: Test field selection by id (<N)"""
        infile = self.infile_convolve
        field = '<7'
        spw = '23:0~40,25:20~60'
        rows = [0,1]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_gt(self):
        """test_field_id_gt: Test field selection by id (>N)"""
        infile = self.infile_convolve
        field = '>6'
        spw = '23:60~100,21:40~80'
        rows = [2,3]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_range(self):
        """test_field_id_range: Test field selection by id ('N~M')"""
        infile = self.infile_convolve
        field = '5~6'
        spw = '23:0~40,25:20~60'
        rows = [0,1]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_list(self):
        """test_field_id_list: Test field selection by id ('N,M')"""
        infile = self.infile_convolve
        field = '5,7'
        spw = '23:0~40,21:40~80'
        rows = [0,2]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_id_exprlist(self):
        """test_field_id_exprlist: Test field selection by id ('EXPR0,EXPR1')"""
        infile = self.infile_convolve
        field = '6~7,>7'
        spw = '23:60~100,25:20~60,21:40~80'
        rows = [1,2,3]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_value_exact(self):
        """test_field_value_exact: Test field selection by name"""
        infile = self.infile_convolve
        field = 'M100'
        spw = '23:0~40,25:20~60'
        rows = [0,1]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_value_pattern(self):
        """test_field_value_pattern: Test field selection by pattern match"""
        infile = self.infile_convolve
        field = 'M*'
        spw = '23:0~40,25:20~60,21:40~80'
        rows = [0,1,2]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_value_list(self):
        """test_field_value_list: Test field selection by name list"""
        infile = self.infile_convolve
        field = 'M30,3C273'
        spw = '23:60~100,21:40~80'
        rows = [2,3]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    def test_field_mix_exprlist(self):
        """test_field_mix_list: Test field selection by name and id"""
        infile = self.infile_convolve
        field = '6,M30,3C273'
        spw = '25:20~60,23:60~100,21:40~80'
        rows = [1,2,3]

        self.__exec_complex_test(infile, ['field', 'spw'], [field, spw], rows)

    ### spw selection syntax test ###
    def test_spw_id_default(self):
        """test_spw_id_default: Test default value for spw"""
        infile = self.infile_shift
        spw = ''
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exact(self):
        """test_spw_id_exact: Test spw selection by id ('N')"""
        infile = self.infile_shift
        spw = '21'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_lt(self):
        """test_spw_id_lt: Test spw selection by id ('<N')"""
        infile = self.infile_shift
        spw = '<24'
        rows = [0,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_gt(self):
        """test_spw_id_lt: Test spw selection by id ('>N')"""
        infile = self.infile_shift
        spw = '>22'
        rows = [0,1,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_range(self):
        """test_spw_id_range: Test spw selection by id ('N~M')"""
        infile = self.infile_shift
        spw = '22~25'
        rows = [0,1,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_list(self):
        """test_spw_id_list: Test spw selection by id ('N,M')"""
        infile = self.infile_shift
        spw = '21,23,25'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exprlist(self):
        """test_spw_id_exprlist: Test spw selection by id ('EXP0,EXP1')"""
        infile = self.infile_shift
        spw = '<22,23~25'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_pattern(self):
        """test_spw_id_pattern: Test spw selection by wildcard"""
        infile = self.infile_shift
        spw = '*'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_frequency(self):
        """test_spw_value_frequency: Test spw selection by frequency range ('FREQ0~FREQ1')"""
        infile = self.infile_shift
        spw = '299.4~299.6GHz'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_velocity(self):
        """test_spw_value_velocity: Test spw selection by velocity range ('VEL0~VEL1')"""
        infile = self.infile_shift
        spw = '-100~100km/s'
        rows = [0,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_mix_exprlist(self):
        """test_spw_mix_exprlist: Test spw selection by id and frequency/velocity range"""
        infile = self.infile_shift
        spw = '<22,-100~100km/s'
        rows = [0,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    ### spw (channel) selection syntax test ###
    def test_spw_id_default_channel(self):
        """test_spw_id_default_channel: Test spw selection with channel range (':CH0~CH1')"""
        infile = self.infile_shift
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_default_frequency(self):
        """test_spw_id_default_frequency: Test spw selection with channel range (':FREQ0~FREQ1')"""
        infile = self.infile_convolve
        spw = ':300470~300510MHz'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_default_velocity(self):
        """test_spw_id_default_velocity: Test spw selection with channel range (':VEL0~VEL1')"""
        infile = self.infile_convolve
        spw = ':-509.6~-469.7km/s'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_default_list(self):
        """test_spw_id_default_list: Test spw selection with multiple channel range (':CH0~CH1;CH2~CH3')"""
        infile = self.infile_duplicate
        spw = ':40~60;65~85'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exact_channel(self):
        """test_spw_id_exact_channel: Test spw selection with channel range ('N:CH0~CH1')"""
        infile = self.infile_convolve
        spw = '25:20~60'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exact_frequency(self):
        """test_spw_id_exact_frequency: Test spw selection with channel range ('N:FREQ0~FREQ1')"""
        infile = self.infile_convolve
        spw = '25:3.0047e5~3.0051e5MHz'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exact_velocity(self):
        """test_spw_id_exact_velocity: Test spw selection with channel range ('N:VEL0~VEL1')"""
        infile = self.infile_convolve
        spw = '25:-5.09647e2~-4.69675e2km/s'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_exact_list(self):
        """test_spw_id_exact_list: Test spw selection with channel range ('N:CH0~CH1;CH2~CH3')"""
        infile = self.infile_duplicate
        spw = '23:40~60;65~85'
        rows = [0,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_pattern_channel(self):
        """test_spw_id_pattern_channel: Test spw selection with channel range ('*:CH0~CH1')"""
        infile = self.infile_shift
        spw = '*:30~70'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_pattern_frequency(self):
        """test_spw_id_pattern_frequency: Test spw selection with channel range ('*:FREQ0~FREQ1')"""
        infile = self.infile_convolve
        spw = '*:300470~300510MHz'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_pattern_velocity(self):
        """test_spw_id_pattern_velocity: Test spw selection with channel range ('*:VEL0~VEL1')"""
        infile = self.infile_convolve
        spw = ':-509.6~-469.7km/s'
        rows = [1]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_pattern_list(self):
        """test_spw_id_pattern_list: Test spw selection with channel range ('*:CH0~CH1;CH2~CH3')"""
        infile = self.infile_duplicate
        spw = '*:40~60;65~85'
        rows = [0,1,2,3]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_frequency_channel(self):
        """test_spw_value_frequency_channel: Test spw selection with channel range ('FREQ0~FREQ1:CH0~CH1')"""
        infile = self.infile_convolve
        spw = '299.4~299.6GHz:40~80'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_frequency_frequency(self):
        """test_spw_value_frequency_frequency: Test spw selection with channel range ('FREQ0~FREQ1:FREQ2~FREQ3')"""
        infile = self.infile_convolve
        spw = '299.4~299.6GHz:299.49~299.53GHz'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_frequency_velocity(self):
        """test_spw_value_frequency_velocity: Test spw selection with channel range ('FREQ0~FREQ1:VEL0~VEL1')"""
        infile = self.infile_convolve
        spw = '299.4~299.6GHz:469.67~509.65km/s'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_frequency_list(self):
        """test_spw_value_frequency_list: Test spw selection with channel range ('FREQ0~FREQ1:CH0~CH1;CH2~CH3')"""
        infile = self.infile_duplicate
        spw = '299.4~299.6GHz:40~60;65~85'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_velocity_channel(self):
        """test_spw_value_velocity_channel: Test spw selection with channel range ('VEL0~VEL1:CH0~CH1')"""
        infile = self.infile_convolve
        spw = '400~600km/s:40~80'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_velocity_frequency(self):
        """test_spw_value_velocity_frequency: Test spw selection with channel range ('VEL0~VEL1:FREQ0~FREQ1')"""
        infile = self.infile_convolve
        spw = '400~600km/s:299.49~299.53GHz'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_velocity_velocity(self):
        """test_spw_value_velocity_velocity: Test spw selection with channel range ('VEL0~VEL1:VEL2~VEL3')"""
        infile = self.infile_convolve
        spw = '400~600km/s:469.67~509.65km/s'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_value_velocity_list(self):
        """test_spw_value_velocity_list: Test spw selection with channel range ('VEL0~VEL1:CH0~CH1;CH2~CH3')"""
        infile = self.infile_duplicate
        spw = '400~600km/s:40~60;65~85'
        rows = [2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    def test_spw_id_list_channel(self):
        """test_spw_id_list_channel: Test spw selection with channnel range ('ID0:CH0~CH1,ID1:CH2~CH3')"""
        infile = self.infile_convolve
        spw = '21:40~80,25:20~60'
        rows = [1,2]

        self.__exec_simple_test(infile, 'spw', spw, rows)

    ### scan selection syntax test ###
    def test_scan_id_default(self):
        """test_scan_id_default: Test default value for scan"""
        infile = self.infile_shift
        scan = ''
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_exact(self):
        """test_scan_id_exact: Test scan selection by id ('N')"""
        infile = self.infile_convolve
        scan = '15'
        spw = ':0~40'
        rows = [0]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_lt(self):
        """test_scan_id_lt: Test scan selection by id ('<N')"""
        infile = self.infile_convolve
        scan = '<16'
        spw = ':0~40'
        rows = [0]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_gt(self):
        """test_scan_id_gt: Test scan selection by id ('>N')"""
        infile = self.infile_convolve
        scan = '>16'
        spw = ':60~100'
        rows = [3]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_range(self):
        """test_scan_id_range: Test scan selection by id ('N~M')"""
        infile = self.infile_convolve
        scan = '15~16'
        spw = '23:0~40,25:20~60,21:40~80'
        rows = [0,1,2]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_list(self):
        """test_scan_id_list: Test scan selection by id ('N,M')"""
        infile = self.infile_convolve
        scan = '15,16'
        spw = '21:40~80,23:0~40,25:20~60'
        rows = [0,1,2]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    def test_scan_id_exprlist(self):
        """test_scan_id_exprlist: Test scan selection by id ('EXP0,EXP1')"""
        infile = self.infile_convolve
        scan = '<16,16'
        spw = '21:40~80,23:0~40,25:20~60'
        rows = [0,1,2]

        self.__exec_complex_test(infile, ['scan', 'spw'], [scan, spw], rows)

    ### pol selection syntax test ###
    def test_pol_id_default(self):
        """test_pol_id_default: Test default value for pol"""
        infile = self.infile_shift
        pol = ''
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_exact(self):
        """test_pol_id_exact: Test pol selection by id ('N')"""
        infile = self.infile_convolve
        pol = '0'
        spw = '21:40~80,23:0~40'
        rows = [0,2]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_lt(self):
        """test_pol_id_lt: Test pol selection by id ('<N')"""
        infile = self.infile_convolve
        pol = '<1'
        spw = '21:40~80,23:0~40'
        rows = [0,2]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_gt(self):
        """test_pol_id_gt: Test pol selection by id ('>N')"""
        infile = self.infile_convolve
        pol = '>0'
        spw = '25:20~60,23:60~100'
        rows = [1,3]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_range(self):
        """test_pol_id_range: Test pol selection by id ('N~M')"""
        infile = self.infile_shift
        pol = '0~1'
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_list(self):
        """test_pol_id_list: Test pol selection by id ('N,M')"""
        infile = self.infile_shift
        pol = '0,1'
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

    def test_pol_id_exprlist(self):
        """test_pol_id_exprlist: Test pol selection by id ('EXP0,EXP1')"""
        infile = self.infile_shift
        pol = '>0,<1'
        spw = ':30~70'
        rows = [0,1,2,3]

        self.__exec_complex_test(infile, ['pol', 'spw'], [pol, spw], rows)

class sdfitold_average(sdfitold_unittest_base, unittest.TestCase):
    """
    """
    infile = 'sdfit_average.asap'
    inputs = [infile]
    outfile = 'sdfitold_average.txt'
    outputs = [outfile]
    fitmode = 'list'
    spw = ':0~40'
    nfit = [1]
    tweight = 'tint'
    pweight = 'tsys'

    def _execute_task(self, timeaverage, scanaverage, polaverage):
        kwargs={'infile': self.infile,
                'outfile': self.outfile,
                'spw': self.spw,
                'nfit': self.nfit,
                'fitmode': self.fitmode,
                'timeaverage': timeaverage,
                'polaverage': polaverage}
        if timeaverage is True:
            kwargs['scanaverage'] = scanaverage
            kwargs['tweight'] = self.tweight
        if polaverage is True:
            kwargs['pweight'] = self.pweight

        return sdfitold(**kwargs)

    def _verify(self, expected):
        # dictionary for verification
        # key: tuple consisting of SCANNO, IFNO, and POLNO
        #     (SCANNO,IFNO,POLNO), e.g.,
        #
        #         (0,0,0)
        #
        # value: flattened list of 'peak', 'cent', and 'fwhm'.
        #     If there are multiple rows for one key, the list
        #     is extended, e.g.,
        #
        #         [peak0, cent0, fwhm0, peak1, cent1, fwhm1, ...]
        #
        result = self.read_result(self.outfile)
        tol = 1.0e-4
        properties = ['peak', 'cent', 'fwhm']
        nprop = len(properties)
        for key in list(expected.keys()):
            self.assertTrue(key in result,
                            msg='result does\'t have key %s'%(list(key)))
            e = expected[key]
            r = result[key]
            ne = len(e)
            nr = len(r)
            self.assertEqual(nr, ne,
                             msg='%s: number of row mismatch (expected %s, actual %s)'%(list(key), ne/nprop, nr/nprop))
            for irow in range(len(e)):
                diff = abs((r[irow] - e[irow])/e[irow])
                self.assertLessEqual(diff, tol,
                                     msg='%s (%s): result differ (expected %s, actual %s)'%(properties[irow % nprop], irow, e[irow], r[irow]))

    def testaverageFFF(self):
        """testaverageFFF: test average (timeaverage=False, scanaverage=False, polaverage=False)"""
        result = self._execute_task(timeaverage=False, scanaverage=False, polaverage=False)
        expected = {(0,0,0): [8.0, 20.0, 5.0, 16.0, 20.0, 5.0],
                    (0,0,1): [16.0, 20.0, 5.0, 32.0, 20.0, 5.0],
                    (1,0,0): [64.0, 20.0, 5.0, 128.0, 20.0, 5.0],
                    (1,0,1): [128.0, 20.0, 5.0, 256.0, 20.0, 5.0]}
        self._verify(expected)

    def testaverageTFF(self):
        """testaverageTFF: test average (timeaverage=True, scanaverage=False, polaverage=False)"""
        result = self._execute_task(timeaverage=True, scanaverage=False, polaverage=False)
        expected = {(0,0,0): [54.0, 20.0, 5.0],
                    (0,0,1): [108.0, 20.0, 5.0]}
        self._verify(expected)

    def testaverageTTF(self):
        """testaverageTTF: test average (timeaverage=True, scanaverage=True, polaverage=False)"""
        result = self._execute_task(timeaverage=True, scanaverage=True, polaverage=False)
        expected = {(0,0,0): [12.0, 20.0, 5.0],
                    (0,0,1): [24.0, 20.0, 5.0],
                    (1,0,0): [96.0, 20.0, 5.0],
                    (1,0,1): [192.0, 20.0, 5.0]}
        self._verify(expected)

    def testaverageTTT(self):
        """testaverageTTT: test average (timeaverage=True, scanaverage=True, polaverage=True)"""
        result = self._execute_task(timeaverage=True, scanaverage=True, polaverage=True)
        expected = {(0,0,0): [18.0, 20.0, 5.0],
                    (1,0,0): [144.0, 20.0, 5.0]}
        self._verify(expected)

    def testaverageFFT(self):
        """testaverageFFT: test average (timeaverage=False, scanaverage=False, polaverage=True)"""
        result = self._execute_task(timeaverage=False, scanaverage=False, polaverage=True)
        expected = {(0,0,0): [81.0, 20.0, 5.0]}
        self._verify(expected)


    def testaverageTFT(self):
        """testaverageTFT: test average (timeaverage=True, scanaverage=False, polaverage=True)"""
        result = self._execute_task(timeaverage=True, scanaverage=False, polaverage=True)
        expected = {(0,0,0): [81.0, 20.0, 5.0]}
        self._verify(expected)

class sdfitold_flag(sdfitold_unittest_base, unittest.TestCase):
    """
    Test flag handling in sdfitold
    * channel flag with list
    * channel flag with auto (line finder)
    * row flag with list
    * row flag with auto (line finder)
    * channel flag and timeaverage (scanaverage=True)
    * channel flag and timeaverage (scanaverage=False)
    * channel flag and polaverage
    * row flag and timeaverage (scanaverage=True)
    * row flag and timeaverage (scanaverage=False)
    * row flag and polaverage
    """
    avefile = 'sdfit_average.asap'
    fitfile1 = 'sd_analytic_type4-1.asap'
    fitfile2 = 'sd_analytic_type5-1.asap'
    inputs = [avefile, fitfile1, fitfile2]
    outfile = 'fit_result.txt'
    outputs = [outfile]
    spw = '23'
    fitmode = 'list'
    nfit = [1]
    tweight = 'tint'
    pweight = 'tsys'

    def _run_test(self, refdata, **kwargs):
        # construct task parameters
        if 'infile' not in kwargs:
            self.fail("infile is not specified")
        predefined = ['outfile', 'spw', 'fitmode', 'nfit',
                      'tweight', 'pweight']
        for p in predefined:
            if not hasattr(self, p):
                self.fail("Internal error: %s not defined in class" % p)
            if p not in kwargs: kwargs[p] = getattr(self, p)
        # save flag for comparison
        tbname = kwargs['infile']
        self._check_file(tbname)
        tb.open(tbname)
        flagtra_pre = tb.getcol('FLAGTRA')
        flagrow_pre = tb.getcol('FLAGROW')
        tb.close()
        # invoke task
        retval = sdfitold(**kwargs)
        # verify fit results
        self._verify_outfile(refdata, kwargs['outfile'])
        # make sure FLAGTRA and FLAGROW are not changed
        tb.open(tbname)
        flagtra_post = tb.getcol('FLAGTRA')
        flagrow_post = tb.getcol('FLAGROW')
        tb.close()
        self.assertTrue((flagrow_post==flagrow_pre).all(),
                        "FLAGROW has been changed by task operation")
        self.assertTrue((flagtra_post==flagtra_pre).all(),
                        "FLAGTRA has been changed by task operation")

    def _flag_table(self, infile, row=[], chan=[]):
        """
        flag infile.
        row: a list of row ids to flag
        chan : a list of [start, end] channels to flag
        when chan=[], FLAGROW is set. Otherwise, channel flag is set
        for the row list.
        """
        if len(row) == 0: return
        flagrow = (len(chan)==0)
        self._check_file(infile)
        tb.open(infile, nomodify=False)
        try:
            if flagrow:
                fl = tb.getcol('FLAGROW')
                for idx in row:
                    fl[idx] = 1
                tb.putcol('FLAGROW', fl)
            else:
                if type(chan[0]) in [int, float]:
                    chan = [chan]
                fl = tb.getcol('FLAGTRA').transpose()
                for irow in row:
                    for (schan, echan) in chan:
                        fl[irow][schan:echan+1] = 1
                tb.putcol('FLAGTRA', fl.transpose())
        except:
            raise
        finally:
            tb.flush()
            tb.close()

    def testChanFlagList(self):
        """test channel flag with fitmode='list'"""
        fitmode='list'
        nfit = [1]
        spw = '21'
        infile = self.fitfile2 # spectra with 2 lines
        self._flag_table(infile, row=[2], chan=[40,65]) #flag 1st line
        refdata = {(16,21,0): [3.757749080657959, 75.0, 5.0]}
        # the first line at 50chan will be detected if not flagged
        self._run_test(refdata,infile=infile,spw=spw,fitmode=fitmode,nfit=nfit)

    def testChanFlagAuto(self):
        """test channel flag with fitmode='auto'"""
        fitmode='auto'
        infile = self.fitfile2  # spectra with 2 lines
        box_size = 0.3
        spw = '21'
        self._flag_table(infile, row=[2], chan=[40,65]) #flag 1st line
        refdata = {(16,21,0): [3.757749080657959, 75.0, 5.0]}
        # two lines would be detected if not flagged
        self._run_test(refdata,infile=infile,spw=spw,
                      fitmode=fitmode,nfit=[],box_size=box_size)

    def testRowFlagList(self):
        """test row flag with fitmode='list'"""
        fitmode='list'
        nfit = [1]
        spw = '*:20~80'
        infile = self.fitfile1 # spectra with 1 line
        self._flag_table(infile, row=[0,1,2])
        refdata = {# flagged row will not be saved in outfile
                   (17,23,1): [5.636623859405518, 50.0, 5.0]}
        # the first line at 50chan will be detected if not flagged
        self._run_test(refdata,infile=infile,spw=spw,fitmode=fitmode,nfit=nfit)

    def testRowFlagAuto(self):
        """test channel flag with fitmode='auto'"""
        fitmode='auto'
        infile = self.fitfile1  # spectra with 1 line
        spw = '*:20~80'
        self._flag_table(infile, row=[0, 1, 2])
        refdata = {# flagged row will not be saved in outfile
                   (17,23,1): [5.636623859405518, 50.0, 5.0]}
        # two lines would be detected if not flagged
        self._run_test(refdata,infile=infile,spw=spw,
                      fitmode=fitmode,nfit=[])

    def testChanFlagSave(self):
        """test channel flag with averaging time=T, scan=T, pol=F"""
        spw = '0:0~40'
        pol='0'
        (tave, save, pave) = (True, True, False)
        infile = self.avefile
        self._flag_table(infile, row=[0], chan=[0,40])
        refdata = {(0,0,0): [16.0, 20.0, 5.0],
                   (1,0,0): [96.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,pol=pol,
                       timeaverage=tave,scanaverage=save,polaverage=pave)

    def testRowFlagSave(self):
        """test row flag with averaging time=T, scan=T, pol=F"""
        spw = '0:0~40'
        pol='0'
        (tave, save, pave) = (True, True, False)
        infile = self.avefile
        self._flag_table(infile, row=[0])
        refdata = {(0,0,0): [16.0, 20.0, 5.0],
                   (1,0,0): [96.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,pol=pol,
                       timeaverage=tave,scanaverage=save,polaverage=pave)

    def testChanFlagTave(self):
        """test channel flag with averaging time=T, scan=F, pol=F"""
        spw = '0:0~40'
        pol='0'
        (tave, save, pave) = (True, False, False)
        infile = self.avefile # spectra with 2 lines
        self._flag_table(infile, row=[0,5], chan=[0,40]) #flag 1st line
        refdata = {(0,0,0): [40.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,pol=pol,
                       timeaverage=tave,scanaverage=save,polaverage=pave)

    def testRowFlagTave(self):
        """test row flag with averaging time=T, scan=F, pol=F"""
        spw = '0:0~40'
        pol='0'
        (tave, save, pave) = (True, False, False)
        infile = self.avefile
        self._flag_table(infile, row=[0,5])
        refdata = {(0,0,0): [40.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,pol=pol,
                       timeaverage=tave,scanaverage=save,polaverage=pave)

    def testChanFlagPave(self):
        """test channel flag with averaging time=F, scan=F, pol=T"""
        spw = '0:0~40'
        (tave, save, pave) = (False, False, True)
        infile = self.avefile
        self._flag_table(infile, row=[0,2,5,7], chan=[0,40])
        refdata = {(0,0,0): [60.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,
                       timeaverage=tave,scanaverage=save,polaverage=pave)

    def testRowFlagPave(self):
        """test row flag with averaging time=F, scan=F, pol=T"""
        spw = '0:0~40'
        (tave, save, pave) = (False, False, True)
        infile = self.avefile
        self._flag_table(infile, row=[0,2,5,7])
        refdata = {(0,0,0): [60.0, 20.0, 5.0]}
        self._run_test(refdata,infile=infile,spw=spw,
                       timeaverage=tave,scanaverage=save,polaverage=pave)


def suite():
    return [sdfitold_test, sdfitold_test_exceptions,
            sdfitold_selection_syntax, sdfitold_average,
            sdfitold_flag]
