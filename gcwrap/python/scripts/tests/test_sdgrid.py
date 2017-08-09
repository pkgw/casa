import os
import sys
import shutil
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import sha
import time
import numpy
import re
import string
from casa_stack_manip import stack_frame_find

try:
    from . import selection_syntax
except:
    import tests.selection_syntax as selection_syntax

#to rethrow exception
g = stack_frame_find( )
g['__rethrow_casa_exceptions'] = True
from sdgridold import sdgridold
from sdutil import tbmanager
import asap as sd

#
# Unit test of sdgridold task.
# 

###
# Base class for sdimaging unit test
###
class sdgridold_unittest_base(object):
    """
    """
    taskname='sdgridold'
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdgrid/'
    data=None
    tolerance=0.01
    outfile='sdgridold.asap.grid'

    def _checkfile( self, name ):
        isthere=os.path.exists(name)
        self.assertEqual(isthere,True,
                         msg='output file %s was not created because of the task failure'%(name))

    def _getdim( self, sp ):
        dim = 0
        import copy
        a = copy.deepcopy(sp)
        while isinstance(a, list):
            dim = dim + 1
            a = a[0]
        return dim
    
    def _checkshape( self, sp, ref ):
        # check array dimension
        self.assertEqual( self._getdim(sp), self._getdim(ref),
                          msg='array dimension differ' )
        # check number of spectra
        self.assertEqual( len(sp), len(ref),
                          msg='number of spectra differ' )
        # check number of channel
        self.assertEqual( len(sp[0]), len(ref[0]),
                          msg='number of channel differ' )
    def _diff(self, sp, ref):
        diff=abs((sp-ref)/ref)
        idx=numpy.argwhere(numpy.isnan(diff))
        #print idx
        if len(idx) > 0:
            diff[idx]=sp[idx]
        return diff
        
    def getdata(self):
        with tbmanager(self.outfile) as tb:
            self.data = tb.getcol('SPECTRA')

    def check(self,ref,val):
        diff=abs((val-ref)/ref)
        #print 'diff=',diff
        self.assertTrue(diff < self.tolerance,
                        msg='grid result differ: ref %s, val %s'%(ref,val))

    def nonzero(self,ref,index):
        #refpix=ref[1]
        #resultpix=index[1]
        refpix = ref
        resultpix = index
        msglt = 'There are nonzero pixels that should be zero'
        msggt = 'There are zero pixels that should be nonzero'
        self.assertEqual(len(refpix),len(resultpix),
                         msg=(msglt if len(refpix) < len(resultpix) else msggt))
        for i in range(len(refpix)):
            self.assertEqual(refpix[i],resultpix[i],
                             msg='Index doesn\'t match: ref %s, result %s'%(refpix[i],resultpix[i]))

    def generateNonzeroPix(self,npol,npix,width):
        index=[]
        start=(npix-1)/2-(width-1)
        end=(npix-1)/2+(width-1)
        #print 'start=',start,',end=',end
        for i in range(start,end+1):
            tweak=npol if (width>=4 and (i==start or i==end)) else 0
            ifrom=npol*npix*i+npol*start+tweak
            ito=ifrom+npol*2*(width-1)-2*tweak
            index+=list(range(ifrom,ito+npol))
        #print 'index=',index
        #nonzeropix_ref=(numpy.zeros(len(index),int),numpy.array(index))
        nonzeropix_ref=numpy.array(index)
        return nonzeropix_ref

    def addrow(self,val):
        tb.open(self.datapath+'/'+self.rawfile)
        tb.copyrows(self.rawfile,0,-1,1)
        tb.close()
        tb.open(self.rawfile,nomodify=False)
        tb.putcell('SPECTRA',tb.nrows()-1,val)
        tb.flush()
        tb.close()

###
# Test on bad parameter settings
###
class sdgridold_failure_case(sdgridold_unittest_base,unittest.TestCase):
    """
    Test on bad parameter setting
    """
    # Input and output names
    #prefix=sdgridold_unittest_base.taskname+'Test0'
    prefix='sdgridoldTest0'
    badid='99'
    rawfile='testimage1chan.1point.asap'

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test000(self):
        """Test 000: Default parameters"""
        # argument verification error
        res=sdgridold()
        self.assertFalse(res)

    def test001(self):
        """Test001: Invalid SPW"""
        try:
            res=sdgridold(infiles=self.rawfile,spw=self.badid,npix=16,cell='20arcsec',outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            #pos=str(e).find('No corresponding rows for given selection: SPW %s'%(self.badid))
            pos=str(e).find('No valid spw')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test002(self):
        """Test002: Invalid POLNO"""
        try:
            res=sdgridold(infiles=self.rawfile,pol=self.badid,npix=16,cell='20arcsec',outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('Empty pol')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test003(self):
        """Test003: Invalid gridfunction"""
        # argument verification error
        res=sdgridold(infiles=self.rawfile,gridfunction='NONE',npix=16,cell='20arcsec',outfile=self.outfile)
        self.assertFalse(res)

    def test004(self):
        """Test004: Invalid weight type"""
        # argument verification error
        res=sdgridold(infiles=self.rawfile,weight='NONE',npix=16,cell='20arcsec',outfile=self.outfile)
        self.assertFalse(res)

    def test005(self):
        """Test005: Check overwrite option"""
        shutil.copytree(self.rawfile,self.outfile)
        try:
            res=sdgridold(infiles=self.rawfile,npix=16,cell='20arcsec',outfile=self.outfile,overwrite=False)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('Output file \'%s\' exists.'%(self.outfile))
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

# Those two tests are meaningless
#    def test006(self):
#        """Test006: Invalid npix"""
#        res=sdgridold(infiles=self.rawfile,npix=-99,cell='',outfile=self.outfile)
#        self.assertFalse(res)
#
#    def test007(self):
#        """Test007: Invalid unit for cell"""
#        res=sdgridold(infiles=self.rawfile,npix=16,cell='20none',outfile=self.outfile)
#        self.assertFalse(res)

    def test008(self):
        """Test008: Invalid format for center coordinate"""
        try:
            res=sdgridold(infiles=self.rawfile,npix=16,cell='20arcsec',outfile=self.outfile,center='Invalid format')
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('Empty QuantumHolder argument for asQuantumDouble')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

        

###
# Test simple gridding
###
class sdgridold_single_integ(sdgridold_unittest_base,unittest.TestCase):
    """
    Test simple gridding using data containing only one integration.
    """
    # Input and output names
    rawfile='testimage1chan.1point.asap'
    #prefix=sdgridold_unittest_base.taskname+'Test1'
    #outfile=prefix+'.asap'

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test100(self):
        """Test 100: Box kernel"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=2
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 10.0
        pol0=self.data[0,nonzeropix[0]]
        self.check(10.0,pol0)

        # pol1 must be 1.0
        pol1=self.data[0,nonzeropix[1]]
        self.check(1.0,pol1)
        

    def test101(self):
        """Test101: SF kernel"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='SF',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # default width for SF is 3
        width=3
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 10.0 while pol1 must be 1.0
        for i in range(0,len(nonzeropix),npol):
            pol0=self.data[0,nonzeropix[i]]
            self.check(10.0,pol0)
            pol1=self.data[0,nonzeropix[i+1]]
            self.check(1.0,pol1)
        

    def test102(self):
        """Test102: Gaussian kernel"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='GAUSS',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # default width for GAUSS is 4
        width=2
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref=numpy.array([218, 219, 220, 221, 222, 223, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 354, 355, 356, 357, 358, 359])
        #nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)
        
        # pol0 must be 10.0 while pol1 must be 1.0
        for i in range(0,len(nonzeropix),npol):
            pol0=self.data[0,nonzeropix[i]]
            self.check(10.0,pol0)
            pol1=self.data[0,nonzeropix[i+1]]
            self.check(1.0,pol1)

    def test103(self):
        """Test103: Gaussian*Jinc kernel"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='GJINC',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # default width for GAUSS is 4
        width=2
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref=numpy.array([252, 253, 254, 255, 256, 257, 286, 287, 288, 289, 290, 291, 320, 321, 322, 323, 324, 325])
        #nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)
        
        # pol0 must be 10.0 while pol1 must be 1.0
        for i in range(0,len(nonzeropix),npol):
            pol0=self.data[0,nonzeropix[i]]
            self.check(10.0,pol0)
            pol1=self.data[0,nonzeropix[i+1]]
            self.check(1.0,pol1)
        

###
# Test clipminmax
###
class sdgridold_clipping(sdgridold_unittest_base,unittest.TestCase):
    """
    Test clipminmax
    """
    rawfile='testimage1chan.1point.asap'
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

        # modification of input file along with this test
        #   - add row to enable clipping
        #   - all polno set to 0
        self.addrow([0.1])
        tb.open(self.rawfile,nomodify=False)
        polno=tb.getcol('POLNO')
        polno[:]=0
        tb.putcol('POLNO',polno)
        tb.close()

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test200(self):
        """Test 200: test clipping"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False,clipminmax=True)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()
        
        # center is only nonzero pixel
        npol=1
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 1.0
        pol0=self.data[0,nonzeropix[0]]
        self.check(1.0,pol0)
        

###
# Test for flag
###
class sdgridold_flagging(sdgridold_unittest_base,unittest.TestCase):
    """
    Test for flag
    test300. test channel flag
    test301. test row flag
    test302. test row flag (all rows of a pol flagged)
    test303. test combination of channel and row flags 

    Input data has 2 pol x 4channels x 2 rows.
    [How to generate input data]
    nchan = 4
    orgscan = 'testimage1chan.1point.asap'
    outscan = ('testgrid%dchan.1point.asap' % nchan)
    sdcoadd(infiles = [orgscan, orgscan], outfile = outscan)
    tb.open(outscan, nomodify=False)
    subtb = tb.query('POLNO==0') # nrow=2
    sp0 = numpy.array([[10.]*nchan, [1.]*nchan])
    subtb.putcol('SPECTRA', sp0.transpose())
    subtb.putcol('FLAGTRA', sp0.transpose()*0)
    subtb.close()
    subtb = tb.query('POLNO==1') # nrow=2
    sp1 = numpy.array([[15.]*nchan, [7.]*nchan])
    subtb.putcol('SPECTRA', sp1.transpose())
    subtb.putcol('FLAGTRA', sp1.transpose()*0)
    subtb.close()
    tb.close()
    """
    rawfile='testgrid4chan.1point.asap'
    outfile='sdgridold_flagging.asap'
    nchan = 4
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test300(self):
        """Test 300: Test channel flagging"""
        # channel flag pol0 data chan=0~1 @row0 and chan=1~2@row1 
        tb.open(self.rawfile,nomodify=False)
        subtb = tb.query("POLNO==0")
        for irow in range(subtb.nrows()):
            fl=subtb.getcell('FLAGTRA',irow)
            fl[irow:irow+2]=1
            subtb.putcell('FLAGTRA',irow,fl)
        subtb.close()
        tb.close()
        # exec task
        # POL0 should be [1., flagged, 10., 5.5]
        # POL1 should be [11.]*nchan
        refdata = [self._create_masked_array([1., numpy.nan, 10., 5.5]),
                   self._create_masked_array([11.]*self.nchan)]
        self.run_test(refdata)

    def test301(self):
        """Test 301: Test row flagging (one of two rows)"""
        # row flag the first row of pol1 data
        tb.open(self.rawfile,nomodify=False)
        subtb = tb.query("POLNO==1")
        fl=subtb.getcell('FLAGROW',0)
        fl=1
        subtb.putcell('FLAGROW',0,fl)
        subtb.close()
        tb.close()
        # exec task
        # POL0 should be [5.5]*nchan
        # POL1 should be [7.]*nchan
        refdata = [self._create_masked_array([5.5]*self.nchan),
                   self._create_masked_array([7.]*self.nchan)]
        self.run_test(refdata)

    def test302(self):
        """Test 302: Test row flagging all rows"""
        # row flag all rows of pol0 data
        tb.open(self.rawfile,nomodify=False)
        subtb = tb.query("POLNO==0")
        fl=subtb.getcol('FLAGROW')
        fl[:]=1
        subtb.putcol('FLAGROW',fl)
        subtb.close()
        tb.close()

        # exec task
        # POL0 should be flagged
        # POL1 should be [11.]*nchan
        refdata = [self._create_masked_array([numpy.nan]*self.nchan, True),
                   self._create_masked_array([11.]*self.nchan)]
        self.run_test(refdata)

    def test303(self):
        """Test 303: Test combination of row and channel flagging"""
        # flag chans=1~2 in the second row of pol0 data and
        # row flag the first row
        tb.open(self.rawfile,nomodify=False)
        subtb = tb.query("POLNO==0")
        fl=subtb.getcell('FLAGTRA',0)
        fl[1:3] = 1
        fl=subtb.putcell('FLAGTRA',0,fl)
        fl=subtb.getcell('FLAGROW',1)
        fl=1
        subtb.putcell('FLAGROW',1,fl)
        subtb.close()
        tb.close()

        # exec task
        # POL0 should be [10., flagged, flagged, 10.]
        # POL1 should be [11.]*nchan
        refdata = [self._create_masked_array([10., numpy.nan, numpy.nan, 10.]),
                   self._create_masked_array([11.]*self.nchan)]
        self.run_test(refdata)

    ####################
    # Helper functions
    ####################
    def _create_masked_array(self, data, mask=None):
        if mask is None: mask=numpy.isnan(numpy.array(data))
        return numpy.ma.masked_array(data, mask)

    def run_test(self, refdata):
        npix=1
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',weight='UNIFORM',
                   npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self._test_results(self.outfile, refdata)
    
    def _test_results(self, filename, refdata):
        # assert output file is generated
        self._checkfile(filename)
        # test result
        tb.open(filename)
        try:
            nrow = tb.nrows()
            spec = tb.getcol('SPECTRA').transpose()
            cflag = tb.getcol('FLAGTRA').transpose()
            rflag = tb.getcol('FLAGROW')
        except: raise
        finally: tb.close()
        # check row number
        self.assertEqual(nrow,len(refdata),
                         msg='output row number differs: %d (expected: %d)' % (nrow, len(refdata)))
        # spectra and flagtra
        for irow in range(nrow):
            ref = refdata[irow]
            sp = spec[irow]
            # FLAGROW should be set if all channels are flagged.
            if ((cflag[irow]!=0).all()):
                self.assertTrue(rflag[irow]!=0, "All channels are flagged but FLAGROW is not set")
            flg = (rflag[irow]!=0) or (cflag[irow]!=0)
            # check num chan
            self.assertEqual(len(sp),len(ref),
                             msg='nchan in row=%d differs: %d (expected: %d)' % (irow, len(sp), len(ref)))
            # check flag and spectra
            spma = self._create_masked_array(sp, flg)
            if ref.mask.all():
                self.assertTrue(spma.mask.all(), msg='spectram in row=%d should been all flagged' % irow)
            else:
                self.assertTrue((spma==ref).all(),
                                msg='spectrum or flag in row=%d differs: %s (expected: %s)' % (irow, str(spma), str(ref)))


###
# Test various weighting
###
class sdgridold_weighting(sdgridold_unittest_base,unittest.TestCase):
    """
    Test various weighting: UNIFORM, TSYS, TINTSYS
    """
    rawfile='testimage1chan.1point.asap'
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

        # modification of input file along with this test
        #   - all polno set to 0
        tb.open(self.rawfile,nomodify=False)
        polno=tb.getcol('POLNO')
        polno[:]=0
        tb.putcol('POLNO',polno)
        tb.close()

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test400(self):
        """Test 400: test UNIFORM weighting"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False,weight='UNIFORM')
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=1
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 5.5 (={10.0+1.0}/{1.0+1.0}
        pol0=self.data[0,nonzeropix[0]]
        self.check(5.5,pol0)
        
    def test401(self):
        """Test 401: test TINT weighting"""
        # modify INTERVAL
        tb.open(self.rawfile,nomodify=False)
        integ=tb.getcol('INTERVAL')
        integ[0]=0.5
        integ[1]=1.0
        tb.putcol('INTERVAL',integ)
        tb.close()

        # exec task
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False,weight='TINT')
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=1
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 4.0 (={10.0*0.5+1.0*1.0}/{0.5+1.0})
        pol0=self.data[0,nonzeropix[0]]
        self.check(4.0,pol0)
        

    def test402(self):
        """Test402: test TSYS weighting"""
        # modify TSYS
        tb.open(self.rawfile,nomodify=False)
        tsys=tb.getcol('TSYS')
        tsys[:,0]=numpy.sqrt(2.0)
        tsys[:,1]=1.0
        tb.putcol('TSYS',tsys)
        tb.close()

        # exec task
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False,weight='TSYS')
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=1
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 4.0 (={10.0*0.5+1.0*1.0}/{0.5+1.0})
        pol0=self.data[0,nonzeropix[0]]
        self.check(4.0,pol0)
        

    def test403(self):
        """Test403: test TINTSYS weighting"""
        # modify TSYS and INTERVAL
        tb.open(self.rawfile,nomodify=False)
        tsys=tb.getcol('TSYS')
        tsys[:,0]=numpy.sqrt(2.0)
        tsys[:,1]=1.0
        tb.putcol('TSYS',tsys)
        integ=tb.getcol('INTERVAL')
        integ[0]=0.5
        integ[1]=1.0
        tb.putcol('INTERVAL',integ)
        tb.close()

        # exec task
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False,weight='TINTSYS')
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=1
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

        # pol0 must be 2.8 (={10.0*0.5*0.5+1.0*1.0*1.0}/{0.5*0.5+1.0*1.0})
        pol0=self.data[0,nonzeropix[0]]
        self.check(2.8,pol0)
        
###
# Test grid map data
###
class sdgridold_map(sdgridold_unittest_base,unittest.TestCase):
    """
    Test grid map data
    """
    rawfile='testimage1chan.map.asap'

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test500(self):
        """Test BOX gridding for map data"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # center is only nonzero pixel
        npol=2
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        nonzeropix=self.data.nonzero()[1]
        self.nonzero(nonzeropix_ref,nonzeropix)

        pol0=self.data[0,nonzeropix[0]]
        #self.check(0.625,pol0)
        #self.check(0.5,pol0)
        self.check(0.6666666667,pol0)
        
        pol1=self.data[0,nonzeropix[1]]
        #self.check(0.0625,pol1)
        #self.check(0.05,pol1)
        self.check(0.06666666667,pol1)

    def test501(self):
        """Test SF gridding for map data"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='SF',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()

        # default width for SF is 3
        width=3
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)

        # check nonzero values
        refdata=[  1.54954410e-04,   1.54954414e-05,   4.63147834e-03,
                   4.63147851e-04,   9.89488605e-03,   9.89488559e-04,
                   4.63147834e-03,   4.63147851e-04,   1.54954410e-04,
                   1.54954414e-05,   4.63147834e-03,   4.63147851e-04,
                   3.81659232e-02,   3.81659227e-03,   6.86512142e-02,
                   6.86512096e-03,   3.81659232e-02,   3.81659227e-03,
                   4.63147834e-03,   4.63147851e-04,   9.89488605e-03,
                   9.89488559e-04,   6.86512142e-02,   6.86512096e-03,
                   1.19758800e-01,   1.19758807e-02,   6.86512142e-02,
                   6.86512096e-03,   9.89488605e-03,   9.89488559e-04,
                   4.63147834e-03,   4.63147851e-04,   3.81659232e-02,
                   3.81659227e-03,   6.86512142e-02,   6.86512096e-03,
                   3.81659232e-02,   3.81659227e-03,   4.63147834e-03,
                   4.63147851e-04,   1.54954410e-04,   1.54954414e-05,
                   4.63147834e-03,   4.63147851e-04,   9.89488605e-03,
                   9.89488559e-04,   4.63147834e-03,   4.63147851e-04,
                   1.54954410e-04,   1.54954414e-05]
        nonzerodata=numpy.take(self.data,nonzeropix,axis=1).squeeze()
        for i in range(len(nonzerodata)):
            self.check(refdata[i],nonzerodata[i])

    def test502(self):
        """Test GAUSS gridding for map data"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='GAUSS',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()
        
        # default width for GAUSS is 4
        width=3
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref = numpy.array([218, 219, 220, 221, 222, 223, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 354, 355, 356, 357, 358, 359])
        #nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)

        refdata = [1.37290766e-03,   1.37290757e-04,   3.63217224e-03,
         3.63217230e-04,   1.37290766e-03,   1.37290757e-04,
         1.37290766e-03,   1.37290757e-04,   2.71596070e-02,
         2.71596084e-03,   7.29541257e-02,   7.29541294e-03,
         2.71596070e-02,   2.71596084e-03,   1.37290766e-03,
         1.37290757e-04,   3.63217224e-03,   3.63217230e-04,
         7.29541257e-02,   7.29541294e-03,   1.98309869e-01,
         1.98309869e-02,   7.29541257e-02,   7.29541294e-03,
         3.63217224e-03,   3.63217230e-04,   1.37290766e-03,
         1.37290757e-04,   2.71596070e-02,   2.71596084e-03,
         7.29541257e-02,   7.29541294e-03,   2.71596070e-02,
         2.71596084e-03,   1.37290766e-03,   1.37290757e-04,
         1.37290766e-03,   1.37290757e-04,   3.63217224e-03,
         3.63217230e-04,   1.37290766e-03,   1.37290757e-04]
        nonzerodata=numpy.take(self.data,nonzeropix,axis=1).squeeze()
        for i in range(len(nonzerodata)):
            self.check(refdata[i],nonzerodata[i])

    def test503(self):
        """Test GJINC gridding for map data"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='GJINC',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')
        self.getdata()
        
        # default width for GAUSS is 4
        width=3
        npol=2
        nonzeropix=self.data.nonzero()[1]
        nonzeropix_ref = numpy.array([252, 253, 254, 255, 256, 257, 286, 287, 288, 289, 290, 291, 320, 321, 322, 323, 324, 325])
        #nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        self.nonzero(nonzeropix_ref,nonzeropix)

        refdata = [0.0337296, 0.00337296, 0.0818698, 0.00818698, 0.0337296,
                   0.00337296, 0.0818698, 0.00818698, 0.16894495, 0.0168945,
                   0.0818698, 0.00818698, 0.0337296, 0.00337296, 0.0818698,
                   0.00818698, 0.0337296, 0.00337296]
        nonzerodata=numpy.take(self.data,nonzeropix,axis=1).squeeze()
        for i in range(len(nonzerodata)):
            self.check(refdata[i],nonzerodata[i])

###
# Test DEC correction
###
class sdgridold_dec_correction(sdgridold_unittest_base,unittest.TestCase):
    """
    Test DEC correction factor for horizontal (R.A.) auto grid setting.
    """
    rawfile='testimage1chan.1point.asap'
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

        # modification of input file along with this test
        #   - declination set to 60deg
        tb.open(self.rawfile,nomodify=False)
        dir=tb.getcol('DIRECTION')
        dir[1,:]=60.0*numpy.pi/180.0
        tb.putcol('DIRECTION',dir)
        tb.close()

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test600(self):
        """Test 600: Test DEC correction factor"""
        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='20arcsec',outfile=self.outfile,plot=False)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')

        # check horizontal and vertical grid
        tb.open(self.outfile)
        dir0=tb.getcell('DIRECTION',0)
        dir1=tb.getcell('DIRECTION',2)
        dir2=tb.getcell('DIRECTION',npix*2)
        tb.close()
        dx=dir0[0]-dir1[0]
        dy=dir2[1]-dir0[1]
        #print 'dx=',dx,',dy=',dy
        diff=abs((0.5*dx-dy)/dy)
        self.assertTrue(diff<0.01,
                        msg='DEC correction is not correct.')


###
# Test to change center for gridding
###
class sdgridold_grid_center(sdgridold_unittest_base,unittest.TestCase):
    """
    Test to change center for gridding
    """
    rawfile='testimage1chan.1point.asap'
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        if (os.path.exists(self.outfile)):
            shutil.rmtree(self.outfile)

    def test700(self):
        """Test 700: Test to change center for gridding"""
        tb.open(self.rawfile)
        dir=tb.getcell('DIRECTION',0)
        tb.close()

        #shift center 3 pixels upward
        nshift=3
        pix=20.0
        dir[1]+=nshift*pix/3600.0*numpy.pi/180.0

        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='%sarcsec'%(pix),outfile=self.outfile,plot=False,center=dir)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')

        self.getdata()
        
        # center is only nonzero pixel
        npol=2
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        #print nonzeropix_ref
        # shift 3 pixels downward
        nonzeropix_ref[0]-=nshift*npol*npix
        nonzeropix_ref[1]-=nshift*npol*npix
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix_ref
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)

    def test701(self):
        """Test 701: Test to change center for gridding (RA)"""
        tb.open(self.rawfile)
        dir=tb.getcell('DIRECTION',0)
        tb.close()

        #shift center 3 pixels upward
        nshift=3
        pix=20.0
        dir[0]+=nshift*pix/3600.0*numpy.pi/180.0

        npix=17
        res=sdgridold(infiles=self.rawfile,gridfunction='BOX',npix=npix,cell='%sarcsec'%(pix),outfile=self.outfile,plot=False,center=dir)
        self.assertEqual(res,None,
                         msg='Any error occurred during gridding')

        self.getdata()
        
        # center is only nonzero pixel
        npol=2
        width=1
        nonzeropix_ref=self.generateNonzeroPix(npol,npix,width)
        #print nonzeropix_ref
        # shift 3 pixels rightwards
        nonzeropix_ref[0]+=nshift*npol
        nonzeropix_ref[1]+=nshift*npol
        nonzeropix=self.data.nonzero()[1]
        #print nonzeropix_ref
        #print nonzeropix
        self.nonzero(nonzeropix_ref,nonzeropix)
    
class sdgridold_selection(selection_syntax.SelectionSyntaxTest,
                       sdgridold_unittest_base,unittest.TestCase):
    """
    Test selection syntax. Selection parameters to test are:
    spw (no channel selection), scan, pol
    """
    # Input and output names
    rawfile='sd_analytic_type1-3.bl.asap'
    infiles=[rawfile]
    reffile='sd_analytic_type1-3.bl.asap'
    prefix=sdgridold_unittest_base.taskname+'TestSel'
    postfix='.grid.asap'
    outname=prefix+postfix

    gfunc='box'
    npix=1
    
    @property
    def task(self):
        return sdgridold
    
    @property
    def spw_channel_selection(self):
        return False

    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        default(sdgridold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def _compare_OLD( self, name, ref ):
        self._checkfile( name )
        # reference data
        tb.open(ref)
        nrow0=tb.nrows()
        rsp=[]
        for irow in range(nrow0):
            rsp.append(tb.getcell('SPECTRA',irow))
        tb.close()
        # check shape
        tb.open(name)
        nrow=tb.nrows()
        self.assertEqual(nrow,nrow0,msg='number of rows mismatch')
        sp=[]
        for irow in range(nrow):
            sp.append(tb.getcell('SPECTRA',irow))
            self.assertEqual(len(sp[irow]),len(rsp[irow]),
                             msg='SPECTRA: number of channel mismatch in row%s'%(irow)) 
        tb.close()
        # check data
        valuetype=type(sp[0][0])
        #print ''
        for irow in range(nrow):
            #print 'irow=%s'%(irow)
            #print '  rsp=%s'%(rsp[irow])
            #print '   sp=%s'%(sp[irow])
            ret=numpy.allclose(rsp[irow],sp[irow])
            self.assertEqual(ret,True,
                             msg='SPECTRA: data differ in row%s'%(irow))

    def _compare( self, name, ref ):
        self._checkfile( name )
        # reference data
        rsp=[]
        with tbmanager(ref) as tb:
            nrow0=tb.nrows()
            for irow in range(nrow0):
                rsp.append(tb.getcell('SPECTRA',irow))
        # check shape
        sp=[]
        with tbmanager(name) as tb:
            nrow=tb.nrows()
            self.assertEqual(nrow,nrow0,msg='number of rows mismatch')
            for irow in range(nrow):
                sp.append(tb.getcell('SPECTRA',irow))
                self.assertEqual(len(sp[irow]),len(rsp[irow]),
                                 msg='SPECTRA: number of channel mismatch in row%s'%(irow)) 
        
        # check data
        valuetype=type(sp[0][0])
        #print ''
        for irow in range(nrow):
            #print 'irow=%s'%(irow)
            #print '  rsp=%s'%(rsp[irow])
            #print '   sp=%s'%(sp[irow])
            ret=numpy.allclose(rsp[irow],sp[irow])
            self.assertEqual(ret,True,
                             msg='SPECTRA: data differ in row%s'%(irow))

    ####################
    # scan
    ####################
    def test_scan_id_default(self):
        """test scan selection (scan='')"""
        scan=''
        self.res=sdgridold(scan=scan,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15,16,17], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)
    
    def test_scan_id_exact(self):
        """ test scan selection (scan='15')"""
        scan = '15'
        self.res=sdgridold(scan=scan,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)

    def test_scan_id_lt(self):
        """ test scan selection (scan='<16')"""
        scan = '<16'
        self.res=sdgridold(scan=scan,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)

    def test_scan_id_gt(self):
        """ test scan selection (scan='>15')"""
        scan = '>15'
        self.res=sdgridold(scan=scan,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [16,17], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)
    
    def test_scan_id_range(self):
        """ test scan selection (scan='15~16')"""
        scan = '15~16'
        self.res=sdgridold(scan=scan,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15,16], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)
    
    def test_scan_id_list(self):
        """ test scan selection (scan='15,17')"""
        scan = '15,17'
        self.res=sdgridold(scan=scan,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15,17], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)
    
    def test_scan_id_exprlist(self):
        """ test scan selection (scan='<16, 17')"""
        scan = '<16, 17'
        self.res=sdgridold(scan=scan,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'SCANNO': [0], 'IFNO': [23]}
        tbselref = {'SCANNO': [15,17], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel, self.reffile, tbselref)

    ####################
    # pol
    ####################
    def test_pol_id_default(self):
        """test pol selection (pol='')"""
        pol=''
        self.res=sdgridold(pol=pol,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [0,1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)

    def test_pol_id_exact(self):
        """ test pol selection (pol='1')"""
        pol = '1'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    def test_pol_id_lt(self):
        """ test pol selection (pol='<1')"""
        pol = '<1'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [0], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)

    def test_pol_id_gt(self):
        """ test pol selection (pol='>0')"""
        pol = '>0'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    def test_pol_id_range(self):
        """ test pol selection (pol='0~1')"""
        pol = '0~1'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [0,1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    def test_pol_id_list(self):
        """ test pol selection (pol='0,1')"""
        pol = '0,1'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [0,1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    def test_pol_id_exprlist(self):
        """test pol selection (pol='<1,1')"""
        pol = '<1,1'
        self.res=sdgridold(pol=pol,spw='23',infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'POLNO': [0,1], 'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    ####################
    # spw
    #
    # note: implement cases of selecting single spw only, since sd.gridder only works for one spw.
    ####################
    def test_spw_id_default_value(self):
        """ test spw selection (when spw not given, spw is set '-1' and ifno of the first row is adopted for spw)"""
        self.res=sdgridold(infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)

    def test_spw_id_default(self):
        """test spw selection (spw='')"""
        pass

    def test_spw_id_exact(self):
        """ test spw selection (spw='23')"""
        spw='23'
        self.res=sdgridold(spw=spw,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)

    def test_spw_id_lt(self):
        """ test spw selection (spw='<25')"""
        pass

    def test_spw_id_gt(self):
        """ test spw selection (spw='>21')"""
        pass

    def test_spw_id_range(self):
        """ test spw selection (spw='21~24')"""
        pass

    def test_spw_id_list(self):
        """ test spw selection (spw='21,22,23,25')"""
        pass

    def test_spw_id_exprlist(self):
        """ test spw selection (spw='<22,>24')"""
        pass

    def test_spw_id_pattern(self):
        """test spw selection (spw='*')"""
        pass
    
    def test_spw_value_frequency(self):
        """test spw selection (spw='300~400GHz')"""
        spw='300~400GHz'
        self.res=sdgridold(spw=spw,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'IFNO': [25]}
        self._comparecal_with_selection(self.outname, tbsel)
    
    def test_spw_value_velocity(self):
        """test spw selection (spw='-50~50km/s')"""
        spw='-50~50km/s'
        self.res=sdgridold(spw=spw,infiles=self.infiles,outfile=self.outname,gridfunction=self.gfunc,npix=self.npix)
        self.assertEqual(self.res,None, msg='Any error occurred during calibration')
        tbsel = {'IFNO': [23]}
        self._comparecal_with_selection(self.outname, tbsel)

    def test_spw_mix_exprlist(self):
        """test spw selection (spw='150~550km/s,>23')"""
        pass
    
    ####################
    # Helper functions
    ####################
    def _comparecal_with_selection( self, name, tbsel={}, ref=None, tbselref=None):
        if ref is None: ref = self.reffile
        if tbselref is None: tbselref = tbsel
        self._checkfile(name)
        sp=self._getspectra_selected(name, tbsel)
        spref=self._getspectra_selected(ref, tbselref)

        self._checkshape( sp, spref )
        
        for irow in range(len(sp)):
            diff=self._diff(numpy.array(sp[irow]),numpy.array(spref[irow]))
            retval=numpy.all(diff<0.01)
            maxdiff=diff.max()
            self.assertEqual( retval, True,
                             msg='calibrated result is wrong (irow=%s): maxdiff=%s'%(irow,diff.max()) )
        del sp, spref

    def _getspectra_selected( self, name, tbsel={} ):
        """
        Returns an array of spectra in rows selected in table.
        
        name  : the name of scantable
        tbsel : a dictionary of table selection information.
                The key should be column name and the value should be
                a list of column values to select.
        """
        isthere=os.path.exists(name)
        self.assertEqual(isthere,True,
                         msg='file %s does not exist'%(name))

        with tbmanager(name) as tb:
            sp = []
            if len(tbsel) == 0:
                for i in range(tb.nrows()):
                    sp.append(tb.getcell('SPECTRA', i).tolist())
            else:
                command = ''
                for key, val in list(tbsel.items()):
                    if len(command) > 0:
                        command += ' AND '
                    command += ('%s in %s' % (key, str(val)))
                try:
                    newtb = tb.query(command)
                    for i in range(newtb.nrows()):
                        sp.append(newtb.getcell('SPECTRA', i).tolist())
                finally:
                    newtb.close()

        return sp

def suite():
    return [sdgridold_failure_case, sdgridold_single_integ,
            sdgridold_clipping, sdgridold_flagging,
            sdgridold_weighting, sdgridold_map,
            sdgridold_dec_correction, sdgridold_grid_center,
            sdgridold_selection]
