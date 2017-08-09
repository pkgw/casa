import os
import sys
import shutil
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import numpy
import re
import math

try:
    from . import selection_syntax
except:
    import tests.selection_syntax as selection_syntax

import asap as sd
from sdflagold import sdflagold
from sdutil import tbmanager

def selection_to_list(row):
    def _selection_to_list(sel):
        elements = sel.split(',')
        for elem in elements:
            if elem.isdigit():
                yield int(elem)
            elif re.match('^[0-9]+~[0-9]+$', elem):
                s = [int(e) for e in elem.split('~')]
                for i in range(s[0], s[1]+1):
                    yield i
    l = set(_selection_to_list(row))
    return list(l)    

class sdflagold_unittest_base:
    """
    Base class for sdflagold unit test
    """
    # Data path of input/output
    datapath = os.environ.get('CASAPATH').split()[0] + \
               '/data/regression/unittest/sdflag/'
    sddatapath = os.environ.get('CASAPATH').split()[0] + \
                 '/data/regression/unittest/singledish/'
    taskname = "sdflagold"

    def _checkfile( self, name ):
        isthere=os.path.exists(name)
        self.assertTrue(isthere,
                         msg='output file %s was not created because of the task failure'%(name))

class sdflagold_test(unittest.TestCase):
    """
    Basic unit tests for task sdflagold.

    The list of tests:
    test00   --- test channel flagging/unflagging
    test01   --- test row flagging/unflagging
    test02   --- test clipping flagging/unflagging
    test03   --- test no data after selection (raises an exception)

    ***NOTE*** These tests are for Scantable only. Tests for the other formats
               which ASAP supports, including MS and SDFITS, are to be made later.

    NOTE2: 'Artificial_Flat.asap' has 6 flat spectra with rms of 1.0.

    """
    # Data path of input/output
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdflag/'
    # Input and output names
    infile = 'Artificial_Flat.asap'

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)

        default(sdflagold)

    def tearDown(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)

    def test00(self):
        """Test 0: channel flagging"""
        infile = self.infile
        mode = 'manual'
        spw = "*:1~3;10~15"
        #maskflag = [[1,3],[10,15]]

        #flag
        result = sdflagold(infile=infile, mode=mode, spw=spw)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        scan = sd.scantable(filename=infile, average=False)
        maskdict = scan.parse_spw_selection(spw)
        maskflag = maskdict[list(maskdict.keys())[0]]
        ansmask = scan.create_mask(maskflag, invert=True)
        for i in range (scan.nrow()):
            mask = scan.get_mask(i)
            res = True
            for j in range(len(mask)):
                if mask[j] != ansmask[j]:
                    res = False
                    break
            self.assertTrue(res)
        del scan

        #unflag
        result = sdflagold(infile=infile, mode=mode, spw=spw,
                         unflag=True)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        scan = sd.scantable(filename=infile, average=False)
        for i in range (scan.nrow()):
            mask = scan.get_mask(i)
            res = True
            for j in range(len(mask)):
                if mask[j] != True:
                    res = False
                    break
            self.assertTrue(res)
        del scan

    def test01(self):
        """Test 1: row flagging"""
        infile = self.infile
        mode = 'rowid'
        flagrow = '2,4'

        #flag
        result = sdflagold(infile=infile, mode=mode, row=flagrow)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        scan = sd.scantable(filename=infile, average=False)
        rowlist = scan.parse_idx_selection('row', flagrow)
        for i in range (scan.nrow()):
            resrflag = scan._getflagrow(i)
            ansrflag = False
            for j in rowlist:
                if i == j:
                    ansrflag = True
                    break
            res = (resrflag == ansrflag)
            self.assertTrue(res)
        del scan

        #unflag
        result = sdflagold(infile=infile, mode=mode, row=flagrow,
                         unflag=True)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        scan = sd.scantable(filename=infile, average=False)
        for i in range (scan.nrow()):
            resrflag = scan._getflagrow(i)
            res = (resrflag == False)
            self.assertTrue(res)
        del scan

    def test02(self):
        """Test 2: clipping"""
        infile = self.infile
        mode = 'clip'
        #clip = True
        clipminmax = [-3.5, 3.5] #clip at 3.5-sigma level, i.e., flag channels at which abs(value) exceeds 3.5.

        #flag
        result = sdflagold(infile=infile, mode=mode, clipminmax=clipminmax)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        ansmlist = [
            [[0,301],[303,8191]],
            [[0,3326],[3328,7753],[7755,8191]],
            [[0,6545],[6547,8191]],
            [[0,3996],[3998,7920],[7922,8191]],
            [[0,3274],[3276,8191]],
            [[0,1708],[1710,7710],[7712,7799],[7801,8191]]
            ]

        scan = sd.scantable(filename=infile, average=False)
        for i in range (scan.nrow()):
            mask = scan.get_mask(i)
            ansmask = scan.create_mask(ansmlist[i])
            res = True
            for j in range(len(mask)):
                if mask[j] != ansmask[j]:
                    res = False
                    break
            self.assertTrue(res)
        del scan

        #unflag
        result = sdflagold(infile=infile, mode=mode, clipminmax=clipminmax,
                         unflag=True)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        scan = sd.scantable(filename=infile, average=False)
        for i in range (scan.nrow()):
            mask = scan.get_mask(i)
            res = True
            for j in range(len(mask)):
                if mask[j] != True:
                    res = False
                    break
            self.assertTrue(res)
        del scan

    def test03(self):
        """Test 3: Failure case (No data to flag)"""
        #test for the default parameters (raises an exception)
        infile = self.infile
        spw = '10' # non-existent IF value
        mode = 'manual'
        try:
            result = sdflagold(infile=infile, mode=mode, spw=spw)
        except Exception as e:
            pos = str(e).find('No valid spw.')
            self.assertNotEqual(pos, -1, msg='Unexpected exception was thrown: %s'%(str(e)))

class sdflagold_test_flagged_data(unittest.TestCase):
    """
    Test for proper handling of channel/row flags.

    The list of tests:
    test_channel_flag   --- test channel flagging
    test_channel_unflag --- test channel unflagging
    test_row_flag       --- test row flagging
    test_row_unflag     --- test row unflagging
    test_clip_flag      --- test clipping flagging
    test_clip_unflag    --- test clipping unflagging

    NOTE:
        'Artificial_Flat_Flagged.asap' has 6 flat spectra with rms of 1.0.
        Channels 8000~8010 are flagged in all rows. Row 3 is row flagged.
        Put spurious data (+/-10.0) in channels 1024 and 8001 for clipping.
    """
    # Data path of input/output
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdflag/'
    # Input and output names
    infile = 'Artificial_Flat_Flagged.asap'

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)

        default(sdflagold)

        self._preserve_flagstate()

    def tearDown(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)

    def _preserve_flagstate(self):
        with tbmanager(self.infile) as tb:
            self.flagrow_org = tb.getcol('FLAGROW')
            self.flagtra_org = tb.getcol('FLAGTRA')

    def _get_maskdict(self, spw):
        s = sd.scantable(self.infile, average=False)
        return s.parse_spw_selection(spw)

    def _prepare_for_verify(self):
        # get resulting flags
        with tbmanager(self.infile) as tb:
            flagrow = tb.getcol('FLAGROW')
            flagtra = tb.getcol('FLAGTRA')

        # basic check
        self.assertEqual(len(flagrow), len(self.flagrow_org), msg='FLAGROW length mismatch: result %s expected %s'%(len(flagrow),len(self.flagrow_org)))
        self.assertEqual(flagtra.shape, self.flagtra_org.shape, msg='FLAGTRA shape mismatch: result %s expected %s'%(list(flagtra.shape),list(self.flagtra_org.shape)))

        # check flagged area
        nrow = len(flagrow)

        return nrow, flagrow, flagtra

    def __verify_channelflag(self, masklist, unflag):
        nrow, flagrow, flagtra = self._prepare_for_verify()
        
        # check flagged area
        flag_value = 0 if unflag else 128
        for irow in range(nrow):
            # row flags should not be affected
            self.assertEqual(flagrow[irow], self.flagrow_org[irow], msg='Row %s: FLAGROW differ (result %s expected %s)'%(irow,flagrow[irow],self.flagrow_org[irow]))

            # verify channel flags
            expected = self.flagtra_org[:,irow].copy()
            # flagged rows should be skipped
            if self.flagrow_org[irow] == 0:
                if len(masklist) == 1:
                    _masklist = masklist[0]
                else:
                    _masklist = masklist[irow]
                # flag should be updated
                for m in _masklist:
                    print('row %s: setting value 128 to range %s'%(irow,m))
                    expected[m[0]:m[1]+1] = flag_value                        
            self.assertTrue(all(flagtra[:,irow] == expected), msg='Row %s: FLAGTRA differ'%(irow))
        
    
    def _verify_channelflag(self, spw, unflag=False):
        maskdict = self._get_maskdict(spw)
        # data only contain IFNO 0
        masklist = [maskdict[0]]

        self.__verify_channelflag(masklist, unflag)


    def _verify_clip(self, threshold, unflag):
        # need spectral data for verifying clip
        with tbmanager(self.infile) as tb:
            spectra = tb.getcol('SPECTRA')

        # detect clipped data
        def gen_clipped(spectra, threshold):
            rms = lambda x: math.sqrt(x.std()**2 + x.mean()**2)
            nrow = spectra.shape[1]
            for irow in range(nrow):
                sp = spectra[:,irow]
                _threshold = abs(threshold) * rms(sp)
                yield numpy.where(abs(sp) > _threshold)[0]
        clipped_channels = list(gen_clipped(spectra, threshold))

        # map to masklist
        masklist = [[[y,y] for y in x] for x in clipped_channels]

        self.__verify_channelflag(masklist, unflag)


    def _verify_rowflag(self, row, unflag):
        rowlist = selection_to_list(row)
        
        nrow, flagrow, flagtra = self._prepare_for_verify()

        # check flagged data
        flag_value = 0 if unflag else 1
        for irow in range(nrow):
            # channel flags should not be affected
            self.assertTrue(all(flagtra[:,irow] == self.flagtra_org[:,irow]), msg='Row %s: FLAGTRA differ'%(irow))

            # verify row flags
            if irow in rowlist:
                expected = flag_value
            else:
                expected = self.flagrow_org[irow]
            self.assertEqual(flagrow[irow], expected, msg='Row %s: FLAGROW differ (result %s expected %s)'%(irow,flagrow[irow],expected))
            

    def test_channel_flag(self):
        """test_channel_flag: channel flagging (unflag=False)"""
        infile = self.infile
        mode = 'manual'
        spw = "*:1~3;10~15"
        #maskflag = [[1,3],[10,15]]
        unflag=False

        #flag
        result = sdflagold(infile=infile, mode=mode, spw=spw, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_channelflag(spw, unflag)

    def test_channel_unflag(self):
        """test_channel_unflag: channel flagging (unflag=True)"""
        infile = self.infile
        mode = 'manual'
        spw = "*:1~3;7995~8005"
        #maskflag = [[1,3],[7995,8005]]
        unflag=True

        #flag
        result = sdflagold(infile=infile, mode=mode, spw=spw, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_channelflag(spw, unflag)

 
    def test_row_flag(self):
        """test_row_flag:: row flagging (unflag=False)"""
        infile = self.infile
        mode = 'rowid'
        flagrow = '2,4'
        unflag=False

        #flag
        result = sdflagold(infile=infile, mode=mode, row=flagrow, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_rowflag(flagrow, unflag)

    def test_row_unflag(self):
        """test_row_unflag:: row flagging (unflag=True)"""
        infile = self.infile
        mode = 'rowid'
        flagrow = '2~4'
        unflag=True

        #flag
        result = sdflagold(infile=infile, mode=mode, row=flagrow, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_rowflag(flagrow, unflag)

    def test_clip_flag(self):
        """test_clip_flag: clipping (unflag=False)"""
        infile = self.infile
        mode = 'clip'
        #clip = True
        threshold = 5.0
        clipminmax = [-threshold, threshold] #clip at 5-sigma level, i.e., flag channels at which abs(value) exceeds 5.0.
        unflag = False

        #flag
        result = sdflagold(infile=infile, mode=mode, clipminmax=clipminmax, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_clip(threshold, unflag)

    def test_clip_unflag(self):
        """test_clip_unflag: clipping (unflag=True)"""
        infile = self.infile
        mode = 'clip'
        #clip = True
        threshold = 5.0
        clipminmax = [-threshold, threshold] #clip at 5-sigma level, i.e., flag channels at which abs(value) exceeds 5.0.
        unflag = True

        #flag
        result = sdflagold(infile=infile, mode=mode, clipminmax=clipminmax, unflag=unflag)

        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")

        self._verify_clip(threshold, unflag)

class sdflagold_test_timerange(unittest.TestCase):
    """
    Basic unit tests for task sdflagold.

    The list of tests:
    test01   --- test row flagging with selection by timerange 'T0~T1'
    test02   --- test row flagging with selection by timerange 'T0'
    test03   --- test row flagging with selection by timerange 'T0+dT'
    test04   --- test row flagging with selection by timerange '>T0'
    test05   --- test row flagging with selection by timerange '<T0'
    test06   --- test row flagging with selection by timerange '>T0'
                 and scan number
    test07   --- test row flagging with selection by timerange '>T0'
                 and field name
    test08   --- test row flagging with selection by timerange 'T0~T1'
                 where T0 and T1 are incomplete time strings
    test09   --- test row flagging with selection by timerange '>T0'
                 where T0 is incomplete time string
    test10   --- test row flagging with selection by timerange 'T0~T1'
                 where only T1 is incomplete time strings
    test11   --- test selection by timerange 'T0~T1' is inclusive

    ***NOTE*** These tests are for Scantable only. Tests for the other formats
               which ASAP supports, including MS and SDFITS, are to be made later.

    NOTE2: 'Artificial_Flat.asap' has 6 flat spectra with rms of 1.0.

    """
    # Data path of input/output
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdflag/'
    # Input and output names
    infile = 'Artificial_Flat.asap'

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)

        # edit TIME column for testing
        table = gentools(['tb'])[0]
        table.open(self.infile, nomodify=False)
        time_column = table.getcol('TIME')
        for irow in range(table.nrows()):
            time_column[irow] += irow * 1000.0 / 86400.0
        table.putcol('TIME', time_column)
        table.close()

        default(sdflagold)

    def tearDown(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)

    def verify(self, data, flag_row_expected):
        scan = sd.scantable(self.infile, average=False)
        flag_row = numpy.array([scan._getflagrow(i) for i in range(scan.nrow())], dtype=bool)
        del scan
        self.assertTrue(all(flag_row_expected==flag_row),
                        msg='FLAGROW is different from expected value: %s (expected: %s)'%(flag_row, flag_row_expected))
        
    def test01(self):
        """timerange test01: test row flagging with selection by timerange 'T0~T1'"""
        # first two rows should be flagged
        timerange = '2006/01/19/01:50:00~2006/01/19/02:24:00'
        flag_row_expected = numpy.array([True, True, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)
        
    def test02(self):
        """timerange test02: test row flagging with selection by timerange 'T0'"""
        # only second row should be flagged
        timerange = '2006/01/19/02:16:45.0'
        flag_row_expected = numpy.array([False, True, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test03(self):
        """timerange test03: test row flagging with selection by timerange 'T0+dT'"""
        # second and third rows should be flagged
        timerange = '2006/01/19/02:08:00.0+0:30:00'
        flag_row_expected = numpy.array([False, True, True, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test04(self):
        """timerange test04: test row flagging with selection by timerange '>T0'"""
        # fourth and subsequent rows should be flagged
        timerange = '>2006/01/19/02:42:00'
        flag_row_expected = numpy.array([False, False, False, True, True, True], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test05(self):
        """timerange test05: test row flagging with selection by timerange '<T0'"""
        # only first row should be flagged
        timerange = '<2006/01/19/02:00:00'
        flag_row_expected = numpy.array([True, False, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test06(self):
        """timerange test06: test row flagging with selection by timerange '>T0' and scan number"""
        # only fifth row should be flagged
        timerange = '>2006/01/19/02:42:00'
        flag_row_expected = numpy.array([False, False, False, False, True, False], dtype=bool)
        scanno = '4'
        sdflagold(infile=self.infile, mode='manual', scan=scanno, timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test07(self):
        """timerange test06: test row flagging with selection by timerange '>T0' and field name"""
        # edit infile
        table = gentools(['tb'])[0]
        table.open(self.infile, nomodify=False)
        table.putcell('SRCNAME', 3, 'SDFLAG2_TEST')
        table.close()
        
        # only fourth row should be flagged
        timerange = '>2006/01/19/02:42:00'
        field = 'SDFLAG2*'
        flag_row_expected = numpy.array([False, False, False, True, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', field=field, timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test08(self):
        """timerange test08: test row flagging with selection by timerange 'T0~T1' where T0 and T1 are incomplete time string"""
        # first two rows should be flagged
        timerange = '01:50:00~02:24:00'
        flag_row_expected = numpy.array([True, True, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)
        
    def test09(self):
        """timerange test09: test row flagging with selection by timerange '>T0'"""
        # fourth and subsequent rows should be flagged
        timerange = '>02:42:00'
        flag_row_expected = numpy.array([False, False, False, True, True, True], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test10(self):
        """timerange test10: test row flagging with selection by timerange 'T0~T1' where only T1 is incomplete time string"""
        # first two rows should be flagged
        timerange = '2006/01/19/01:50:00~02:24:00'
        flag_row_expected = numpy.array([True, True, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test11(self):
        """test11: test selection by timerange 'T0~T1' is inclusive"""
        # edit infile
        T0 = '2006/01/19/01:52:05'
        T1 = '2006/01/19/02:08:45'
        qa = qatool()
        mjd0 = qa.totime(T0)['value']
        mjd1 = qa.totime(T1)['value']
        table = gentools(['tb'])[0]
        table.open(self.infile, nomodify=False)
        table.putcell('TIME', 0, mjd0)
        table.putcell('TIME', 1, mjd1)
        table.close()
        
        # first two rows should be flagged
        timerange = '%s~%s'%(T0,T1)
        flag_row_expected = numpy.array([True, True, False, False, False, False], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)

    def test12(self):
        """test12: test with the data that all rows are flagged"""
        # edit infile
        table = gentools(['tb'])[0]
        table.open(self.infile, nomodify=False)
        flagrow = table.getcol('FLAGROW')
        flagrow[:] = 1
        table.putcol('FLAGROW', flagrow)
        table.close()
        
        # all rows are flagged
        T0 = '2006/01/19/01:52:05'
        T1 = '2006/01/19/02:08:45'
        timerange = '%s~%s'%(T0,T1)
        flag_row_expected = numpy.array([True, True, True, True, True, True], dtype=bool)
        sdflagold(infile=self.infile, mode='manual', timerange=timerange)

        # verification
        self.verify(self.infile, flag_row_expected)
        
        
class sdflagold_selection(selection_syntax.SelectionSyntaxTest,
                       sdflagold_unittest_base,unittest.TestCase):
    """
    Test selection syntax in sdflagold.
    Selection parameters to test are:
    field, spw (with channel selection), timerange, scan, pol,
    beam, row
    """
    rawfile = 'sd_analytic_type1-3.bl.asap'
    ntbrow = 4 # number of rows in rawfile
    #single_chan = {'syntax': "40~60", 'chanlist': [(40,60)]}
    #single_freq = {'syntax': "", 'chanlist': [(,)]}
    #single_velo = {'syntax': "", 'chanlist': [(,)]}
    #flagchan = ( (80,90), (60,70), (40,50), (20,30) )
    #selchan = ()
    outfile = ''
    
    @property
    def task(self):
        return sdflagold
    
    @property
    def spw_channel_selection(self):
        return True

    def setUp(self):
        self.res=None
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)

        shutil.copytree(self.sddatapath+self.rawfile, self.rawfile)

        default(sdflagold)
        self.mode = 'manual'
        self.outfile = ''

    def tearDown(self):
        pass
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)

    ####################
    # Additional tests
    ####################
    # Tests of row ID selection
    def test_row_id_default(self):
        """test row selection (row='')"""
        self.mode='rowid'
        row = ''
        try:
            self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_row_id_exact(self):
        """test row selection (row='2')"""
        self.mode='rowid'
        row = '2'
        ref_row = [2]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_row_id_lt(self):
        """test row selection (row='<2')"""
        self.mode='rowid'
        row = '<2'
        ref_row = [0,1]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)

    def test_row_id_gt(self):
        """test row selection (row='>1')"""
        self.mode='rowid'
        row = '>1'
        ref_row = [2,3]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)

    def test_row_id_range(self):
        """test row selection (row='1~2')"""
        self.mode='rowid'
        row = '1~2'
        ref_row = [1,2]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)

    def test_row_id_list(self):
        """test row selection (row='0,2')"""
        self.mode='rowid'
        row = '0,2'
        ref_row = [0,2]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)

    def test_row_id_exprlist(self):
        """test row selection (row='>2,1')"""
        self.mode='rowid'
        row = '>2,1'
        ref_row = [1,3]
        self.res=self.run_task(infile=self.rawfile,row=row,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)

    # Tests of channel flag with other selection
    def test_chanflag_scan(self):
        """test channel flagging with scan selection (scan='', spw='*:20~30;80-90')"""
        scan = '16'
        spw = '*:20~30;80~90'
        ref_row = [1,2]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=sdflagold(infile=self.rawfile,scan=scan,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)
    
    def test_chanflag_pol(self):
        """test channel flagging with pol selection (pol='1', spw='*:20~30;80~90')"""
        pol = '1'
        spw = '*:20~30;80~90'
        ref_row = [1,3]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=sdflagold(infile=self.rawfile, pol=pol,spw=spw, mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)
    
    def test_chanflag_beam(self):
        """test channel flagging with beam selection (beam='13', spw='*:20~30;80~90')"""
        beam = '13'
        spw = '*:20~30;80~90'
        ref_row = [2,3]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=sdflagold(infile=self.rawfile, beam=beam,spw=spw, mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)
    
    def test_chanflag_field(self):
        """test channel flagging with field selection (field='M*', spw='*:20~30;80~90')"""
        field = 'M*'
        spw = '*:20~30;80~90'
        ref_row = [0,1,2]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=sdflagold(infile=self.rawfile, field=field,spw=spw, mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)
    
    def test_chanflag_timerange(self):
        """test channel flagging with timerange selection (timerange='2011/11/11/02:33:03.47', spw='*:20~30;80~90')"""
        timerange = '2011/11/11/02:33:03.47'
        spw = '*:20~30;80~90'
        ref_row = [1,2]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=sdflagold(infile=self.rawfile, timerange=timerange,spw=spw, mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)
    
    ####################
    # scan
    ####################
    def test_scan_id_default(self):
        """test scan selection (scan='') This should raise error."""
        scan = ''
        try:
            self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_scan_id_exact(self):
        """test scan selection (scan='16')"""
        scan = '16'
        ref_row = [1,2]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_scan_id_lt(self):
        """test scan selection (scan='<16')"""
        scan = '<16'
        ref_row = [0]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_scan_id_gt(self):
        """test scan selection (scan='>16')"""
        scan = '>16'
        ref_row = [3]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_scan_id_range(self):
        """test scan selection (scan='16~17')"""
        scan = '16~17'
        ref_row = [1,2,3]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_scan_id_list(self):
        """test scan selection (scan='15,17')"""
        scan = '15,17'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_scan_id_exprlist(self):
        """test scan selection (scan='15,>16')"""
        scan = '15,>16'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,scan=scan,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    ####################
    # beam
    ####################
    def test_beam_id_default(self):
        """test beam selection (beam='') This should raise error."""
        beam = ''
        try:
            self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_beam_id_exact(self):
        """test beam selection (beam='13')"""
        beam = '13'
        ref_row = [2,3]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_beam_id_lt(self):
        """test beam selection (beam='<12')"""
        beam = '<12'
        ref_row = [0]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_beam_id_gt(self):
        """test beam selection (beam='>12')"""
        beam = '>12'
        ref_row = [2,3]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_beam_id_range(self):
        """test beam selection (beam='11~12')"""
        beam = '11~12'
        ref_row = [0,1]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_beam_id_list(self):
        """test beam selection (beam='11,12')"""
        beam = '11,12'
        ref_row = [0,1]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_beam_id_exprlist(self):
        """test beam selection (beam='<12,13')"""
        beam = '<12,13'
        ref_row = [0,2,3]
        self.res=self.run_task(infile=self.rawfile,beam=beam,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    ####################
    # pol
    ####################
    def test_pol_id_default(self):
        """test pol selection (pol='') This should raise error."""
        pol = ''
        try:
            self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_pol_id_exact(self):
        """test pol selection (pol='1')"""
        pol = '1'
        ref_row = [1,3]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_pol_id_lt(self):
        """test pol selection (pol='<1')"""
        pol = '<1'
        ref_row = [0,2]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_pol_id_gt(self):
        """test pol selection (pol='>0')"""
        pol = '>0'
        ref_row = [1,3]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_pol_id_range(self):
        """test pol selection (pol='0~1')"""
        pol = '0~1'
        # it's possible to flag all rows if user explicitly specify IDs
        ref_row = [0,1,2,3]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_pol_id_list(self):
        """test pol selection (pol='0,1')"""
        pol = '0,1'
        # it's possible to flag all rows if user explicitly specify IDs
        ref_row = [0,1,2,3]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_pol_id_exprlist(self):
        """test pol selection (pol='1,<1')"""
        pol = '1,<1'
        # it's possible to flag all rows if user explicitly specify IDs
        ref_row = [0,1,2,3]
        self.res=self.run_task(infile=self.rawfile,pol=pol,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    ####################
    # field
    ####################
    def test_field_value_default(self):
        """test field selection (field='') This should raise error."""
        field = ''
        try:
            self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_field_id_exact(self):
        """test field selection (field='6')"""
        field = '6'
        ref_row = [1]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_id_lt(self):
        """test field selection (field='<6')"""
        field = '<6'
        ref_row = [0]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_id_gt(self):
        """test field selection (field='>7')"""
        field = '>7'
        ref_row = [3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_id_range(self):
        """test field selection (field='6~8')"""
        field = '6~8'
        ref_row = [1,2,3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_id_list(self):
        """test field selection (field='5,8')"""
        field = '5,8'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_id_exprlist(self):
        """test field selection (field='5,>7')"""
        field = '5,>7'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_value_exact(self):
        """test field selection (field='M30')"""
        field = 'M30'
        ref_row = [2]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_value_pattern(self):
        """test field selection (field='M*')"""
        field = 'M*'
        ref_row = [0,1,2]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_value_list(self):
        """test field selection (field='3C273,M30')"""
        field = '3C273,M30'
        ref_row = [2,3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_field_mix_exprlist(self):
        """test field selection (field='<6,3*')"""
        field = '<6,3*'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,field=field,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    ####################
    # spw 
    ####################
    def test_spw_id_default(self):
        """test spw selection (spw='', all channels) This should raise error."""
        spw = ''
        try:
            self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_spw_id_exact(self):
        """test spw selection (spw='23', all channels)"""
        spw = '23'
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_lt(self):
        """test spw selection (spw='<23', all channels)"""
        spw = '<23'
        ref_row = [2]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_gt(self):
        """test spw selection (spw='>23', all channels)"""
        spw = '>23'
        ref_row = [1]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_range(self):
        """test spw selection (spw='23~25', all channels)"""
        spw = '23~25'
        ref_row = [0,1,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_list(self):
        """test spw selection (spw='21,25', all channels)"""
        spw = '21,25'
        ref_row = [1,2]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_exprlist(self):
        """test spw selection (spw='23,>24', all channels)"""
        spw = '23,>24'
        ref_row = [0,1,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_id_pattern(self):
        """test spw selection (spw='*', all channels)  This should raise error."""
        spw = '*'
        # it's possible to flag all rows if user explicitly specify IDs
        ref_row = [0,1,2,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_value_frequency(self):
        """test spw selection (spw='299.52~300.47GHz', all channels)"""
        spw = '299.52~300.47GHz' # IFNO=23 will be selected
        ref_row = [0,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_value_velocity(self):
        """test spw selection (spw='-510.~470.km/s', all channels)"""
        spw = '-510.~470km/s' # IFNO=23,25 will be selected
        ref_row = [0,1,3]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_spw_mix_exprlist(self):
        """test spw selection (spw='25,0~501km/s', all channels)"""
        spw = '25,0~501km/s' # all IFs will be selected
        ref_row = [0,1,2,3]
        # it's possible to flag all rows if user explicitly specify IDs
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    #########################
    # spw with channel range
    #########################
    def test_spw_id_default_channel(self):
        """test spw selection w/ channel selection (spw=':40~60')"""
        spw = ':40~60'
        ref_row = [0,1,2,3]
        ref_chanlist = [ (40,60) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)

    def test_spw_id_default_frequency(self):
        """test spw selection w/ channel selection (spw=':299.490~299.500GHz')"""
        spw = ':299.490~299.500GHz' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_default_velocity(self):
        """test spw selection w/ channel selection (spw=':-519.650~-509.640km/s')"""
        spw = ':-519.650~-509.640km/s' # IFNO=25, channel=60~70 will be selected
        ref_row = [1]
        ref_chanlist = [ (60,70) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_default_list(self):
        """test spw selection w/ channel selection (spw=':20~30;80~90')"""
        spw = ':20~30;80~90'
        ref_row = [0,1,2,3]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_exact_channel(self):
        """test spw selection w/ channel selection (spw='23:40~60')"""
        spw = '23:40~60' # IFNO=23, channel 40~60 will be selected
        ref_row = [0,3]
        ref_chanlist = [ (40,60) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_exact_frequency(self):
        """test spw selection w/ channel selection (spw='21:299.490~299.500GHz')"""
        spw = '21:299.4899~299.5001GHz' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_exact_velocity(self):
        """test spw selection w/ channel selection (spw='25:-519.650~-509.640km/s')"""
        spw = '25:-519.650~-509.640km/s'
        ref_row = [1]
        ref_chanlist = [ (60,70) ] # IFNO=25, channel=60~70 will be selected
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_exact_list(self):
        """test spw selection w/ channel selection (spw='25:20~30;80~90')"""
        spw = '25:20~30;80~90'
        ref_row = [1]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_pattern_channel(self):
        """test spw selection w/ channel selection (spw='*:40~60')"""
        spw = '*:40~60'
        ref_row = [0,1,2,3]
        ref_chanlist = [ (40,60) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_pattern_frequency(self):
        """test spw selection w/ channel selection (spw='*:299.490~299.500GHz')"""
        spw = '*:299.490~299.500GHz' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_pattern_velocity(self):
        """test spw selection w/ channel selection (spw='*:-519.650~-509.640km/s')"""
        spw = '*:-519.650~-509.640km/s' # IFNO=25, channel=60~70 will be selected
        ref_row = [1]
        ref_chanlist = [ (60,70) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_pattern_list(self):
        """test spw selection w/ channel selection (spw='*:20~30;80~90')"""
        spw = '*:20~30;80~90'
        ref_row = [0,1,2,3]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_frequency_channel(self):
        """test spw selection w/ channel selection (spw='299.490~299.510GHz:40~50')"""
        spw = '299.490~299.510GHz:40~50' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_frequency_frequency(self):
        """test spw selection w/ channel selection (spw='299.490~299.510GHz:299.490~299.500GHz')"""
        spw = '299.490~299.510GHz:299.490~299.500GHz' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_frequency_velocity(self):
        """test spw selection w/ channel selection (spw='300.49~300.51GHz:-519.650~-509.640km/s')"""
        spw = '300.49~300.51GHz:-519.650~-509.640km/s'  # IFNO=25, channel=60~70 will be selected
        ref_row = [1]
        ref_chanlist = [ (60,70) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_frequency_list(self):
        """test spw selection w/ channel selection (spw='300.49~300.51GHz:20~30;80~90')"""
        spw = '300.49~300.51GHz:20~30;80~90' # IFNO=25, channel=20~30 and 80~90 will be selected
        ref_row = [1]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_velocity_channel(self):
        """test spw selection w/ channel selection (spw='-30~30km/s:40~60')"""
        spw = '-30~30km/s:40~60' # IFNO=23, channel 40~60 will be selected
        ref_row = [0,3]
        ref_chanlist = [ (40,60) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_velocity_frequency(self):
        """test spw selection w/ channel selection (spw='490~510km/s:299.490~299.500GHz')"""
        spw = '490~510km/s:299.490~299.500GHz' # IFNO=21, channel 40~50 will be selected
        ref_row = [2]
        ref_chanlist = [ (40,50) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_velocity_velocity(self):
        """test spw selection w/ channel selection (spw='-510~-490km/s:-519.650~-509.640km/s')"""
        spw = '-510~-490km/s:-519.650~-509.640km/s'  # IFNO=25, channel=60~70 will be selected
        ref_row = [1]
        ref_chanlist = [ (60,70) ]
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_value_velocity_list(self):
        """test spw selection w/ channel selection (spw='-510~-490km/s:20~30;80~90')"""
        spw = '-510~-490km/s:20~30;80~90' # IFNO=25, channel=20~30 and 80~90 will be selected
        ref_row = [1]
        ref_chanlist = ( (20,30), (80,90) )
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row, chanlist=ref_chanlist)        

    def test_spw_id_list_channel(self):
        """test spw selection w/ channel selection (spw='21:40~50,25:60~70')"""
        spw = '21:40~50,25:60~70'
        ref_dict = {2:[ (40,50) ], 1: [ (60,70) ]}
        self.res=self.run_task(infile=self.rawfile,spw=spw,mode=self.mode,outfile=self.outfile,outform='ASAP')
        # Test
        tb.open(self.rawfile)
        nrow = tb.nrows()
        flagrow = tb.getcol("FLAGROW")
        flagtra = tb.getcol("FLAGTRA").transpose()
        tb.close()
        self.assertEqual(nrow, self.ntbrow, "The number of rows in table changed.")
        print("Testing FLAGROW")
        for fr in flagrow:
            self.assertEqual(fr, 0, "FLAGROW should be all 0")
        print("Testing FLAGROW")
        for irow in range(self.ntbrow):
            if irow in list(ref_dict.keys()): currlist = ref_dict[irow]
            else: currlist = []
            flag = flagtra[irow]
            # convert flagtra to bool list
            bflag = [ True if f!=0 else False for f in flag ]
            nchan = len(flag)
            ref_flag = self._get_bool_array(nchan, masklist=currlist)
            self._exact_compare_array(bflag, ref_flag)


    ####################
    # timerange
    ####################
    def test_timerange_value_default(self):
        """test timerange selection (timerange='') This should raise error."""
        timerange = ''
        try:
            self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
            self.fail('The task must throw exception')
        except Exception as e:
            pos=str(e).find('Trying to flag whole scantable.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test_timerange_value_exact(self):
        """test timerange selection (timerange='2011/11/11/02:33:03.47')"""
        timerange = '2011/11/11/02:33:03.47' # SCANNO=16 will be selected
        ref_row = [1,2]
        self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_timerange_value_lt(self):
        """test timerange selection (timerange='<2011/11/11/02:33:03.47')"""
        timerange = '<2011/11/11/02:33:03.47' # SCANNO=15 will be selected
        ref_row = [0]
        self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_timerange_value_gt(self):
        """test timerange selection (timerange='>2011/11/11/02:33:03.5')"""
        timerange = '>2011/11/11/02:33:03.5' # SCANNO=17 will be selected
        ref_row = [3]
        self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_timerange_value_range(self):
        """test timerange selection (timerange='2011/11/11/02:33:03.47~02:34:03.48')"""
        timerange = '2011/11/11/02:33:03.47~02:34:03.48' # SCANNO=16,17 will be selected
        ref_row = [1,2,3]
        self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        

    def test_timerange_value_interval(self):
        """test timerange selection (timerange='2011/11/11/02:33:03.47+00:01:00.1')"""
        timerange = '2011/11/11/02:33:03.47+00:01:00.1' # SCANNO=16,17 will be selected
        ref_row = [1,2,3]
        self.res=self.run_task(infile=self.rawfile,timerange=timerange,mode=self.mode,outfile=self.outfile,outform='ASAP')
        self._test_flag(self.rawfile, ref_row)        


    ####################
    # Helper functions
    ####################
    def _test_flag(self, name, rowlist, chanlist=None):
        """
        Test FLAGROW and FLAGTRA
        Arguments
            name     : name of scantable to test
            rowlist  : a list of row IDs to be flagged
            chanlist : a list of channel ranges to be flagged
                       If None or [], row flagging is assumed
        """
        if not os.path.exists(name):
            self.fail("Scantable '%s' does not exist." % data)
        if chanlist is None: chanlist = []
        
        tb.open(name)
        flagtra = tb.getcol("FLAGTRA").transpose()
        flagrow = tb.getcol("FLAGROW")
        nrow = tb.nrows()
        tb.close()
        self.assertEqual(nrow, self.ntbrow, "The number of rows in table changed.")
        idx_flagrow = rowlist if len(chanlist) == 0 else []
        # compare FLAGROW
        # create reference flagrow array
        ref_flagrow = [ 0 for idx in range(self.ntbrow) ]
        for irow in idx_flagrow: ref_flagrow[irow] = 1
        print("Testing FLAGROW")
        self._exact_compare_array(flagrow, ref_flagrow)
        # compare FLAGTRA
        print("Testing FLAGTRA")
        for irow in range(self.ntbrow):
            flag = flagtra[irow]
            # convert flagtra to bool list
            bflag = [ True if f!=0 else False for f in flag ]
            nchan = len(flag)
            currlist = chanlist if irow in rowlist else []
            ref_flag = self._get_bool_array(nchan, masklist=currlist)
            self._exact_compare_array(bflag, ref_flag)

    def _exact_compare_array(self, data, reference):
        self.assertEqual(len(data), len(reference), "Number of elements in two arrays are different.")
        for idx in range(len(data)):
            self.assertEqual(data[idx], reference[idx], "Value in idx=%d differs: %s (expected: %s)" % (idx, str(data[idx]), str(reference[idx])) )

    def _get_bool_array(self, nelements, masklist=None):
        """
        Create a bool array of length, nelements.
        Elements of the array are set True if the indices are in range
        of masklist. Otherwise, False.
        """
        self.assertTrue(nelements>0, "Internal error. Length of array should be positive value")
        if masklist is None: masklist=[]
        ret_array = numpy.array([ False for idx in range(nelements) ])
        for irange in range(len(masklist)):
            curr_range = masklist[irange]
            if len(curr_range) < 2:
                self.fail("Internal error. masklist should be a list of 2 elements array.")
            schan = curr_range[0]
            echan = curr_range[1]
            ret_array[schan:echan+1] = True
        
        return ret_array
            
        
def suite():
    return [sdflagold_test, sdflagold_test_flagged_data,
            sdflagold_test_timerange, sdflagold_selection]
