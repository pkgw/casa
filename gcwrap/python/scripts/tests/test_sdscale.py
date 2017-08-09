import itertools
import numpy
import os
import re
import sha
import shutil
import string
import sys
import time

from __main__ import default
from taskinit import *
from tasks import *
import unittest

from sdscaleold import sdscaleold
from sdutil import tbmanager
import asap as sd

#
# Unit test of sdscaleold task.
# 

###
# Base class for sdscaleold unit test
###
class sdscaleold_unittest_base:
    """
    Base class for sdscaleold unit test
    """
    taskname='sdscaleold'
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdscale/'
    
    def _checkfile( self, name ):
        isthere=os.path.exists(name)
        self.assertEqual(isthere,True,
                         msg='output file %s was not created because of the task failure'%(name))

    def _compare( self, name, ref, factor, scaletsys ):
        self._checkfile( name )
        # get original nchan and nrow
        tb.open(ref)
        nrow0=tb.nrows()
        rspchans=[]
        rtsyschans=[]
        for irow in range(nrow0):
            rspchans.append(len(tb.getcell('SPECTRA',irow)))
            rtsyschans.append(len(tb.getcell('TSYS',irow)))
        tb.close()
        # check shape
        tb.open(name)
        nrow=tb.nrows()
        self.assertEqual(nrow,2,
                         msg='number of rows mismatch')
        sp=[]
        tsys=[]
        for irow in range(nrow):
            sp.append(tb.getcell('SPECTRA',irow))
            tsys.append(tb.getcell('TSYS',irow))
            self.assertEqual(len(sp[irow]),rspchans[irow],
                             msg='SPECTRA: number of channel mismatch in row%s'%(irow)) 
            self.assertEqual(len(tsys[irow]),rtsyschans[irow],
                             msg='TSYS: number of channel mismatch in row%s'%(irow))
        tb.close()
        # check data
        valuetype=type(sp[0][0])
        if type(factor) is not list:
            # scalar factor
            for irow in range(nrow):
                arrs = numpy.ones(rspchans[irow],dtype=valuetype)*factor
                ret=numpy.allclose(arrs,sp[irow])
                self.assertEqual(ret,True,
                                 msg='SPECTRA: data differ in row%s'%(irow))
                arrt = numpy.ones(rtsyschans[irow],dtype=valuetype)
                if scaletsys:
                    arrt *= factor
                ret=numpy.allclose(arrt,tsys[irow])
                self.assertEqual(ret,True,
                                 msg='TSYS: data differ in row%s'%(irow))
        elif type(factor[0]) is not list:
            # 1D array factor
            factor = numpy.array(factor)
            for irow in range(nrow):
                ret=numpy.allclose(factor,sp[irow])
                self.assertEqual(ret,True,
                                 msg='SPECTRA: data differ in row%s'%(irow))
                arrt = numpy.ones(rtsyschans[irow],dtype=valuetype)
                if scaletsys:
                    arrt *= factor
                ret=numpy.allclose(arrt,tsys[irow])
                self.assertEqual(ret,True,
                                 msg='TSYS: data differ in row%s'%(irow))
        elif len(factor[0]) == 1:
            # 2D array with shape [nrow,1]
            for irow in range(nrow):
                arrs = numpy.ones(rspchans[irow],dtype=valuetype)*factor[irow]
                ret=numpy.allclose(arrs,sp[irow])
                self.assertEqual(ret,True,
                                 msg='SPECTRA: data differ in row%s'%(irow))
                arrt = numpy.ones(rtsyschans[irow],dtype=valuetype)
                if scaletsys:
                    arrt *= factor[irow]
                ret=numpy.allclose(arrt,tsys[irow])
                self.assertEqual(ret,True,
                                 msg='TSYS: data differ in row%s'%(irow))
        else:
            # 2D array with shape [nrow,nchan]
            for irow in range(nrow):
                arrs = numpy.array(factor[irow])
                ret=numpy.allclose(arrs,sp[irow])
                self.assertEqual(ret,True,
                                 msg='SPECTRA: data differ in row%s'%(irow))
                arrt = numpy.ones(rtsyschans[irow],dtype=valuetype)
                if scaletsys:
                    arrt *= numpy.array(factor[irow])
                ret=numpy.allclose(arrt,tsys[irow])
                self.assertEqual(ret,True,
                                 msg='TSYS: data differ in row%s'%(irow))

                
###
# Test on bad parameter settings
###
class sdscaleold_test0(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on bad parameter setting

    Test data, sdscale1.asap, is artificial data with the following
    status:

       - nrow = 2
       - nchan = 4 for spectral data
       - 1 Tsys value for each spectrum
       - all spectral values are 1.0
       - all Tsys values are 1.0
       
    """
    # Input and output names
    rawfile='sdscale1.asap'
    prefix=sdscaleold_unittest_base.taskname+'Test0'
    outfile=prefix+'.asap'

    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdscaleold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test000(self):
        """Test 000: Default parameters"""
        # argument verification error
        res=sdscaleold()
        self.assertFalse(res)        

    def test001(self):
        """Test 001: Existing outfile with overwrite=False"""
        os.system('cp -r %s %s'%(self.rawfile,self.outfile))
        try:
            res=sdscaleold(infile=self.rawfile,factor=2.0,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('Output file \'%s\' exists.'%(self.outfile))
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test002(self):
        """Test 002: Bad shaped factor"""
        factor = [2.0,3.0]
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=False,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('Vector size must be 1 or be same as number of channel.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))
        
    def test003(self):
        """Test 003: Try to scale non-conform Tsys"""
        factor=[1.0,2.0,3.0,4.0]
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=True,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('SPECTRA and TSYS must conform in shape if you want to apply operation on Tsys.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

###
# Test on scaling 1
###
class sdscaleold_test1(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on actual scaling

    Test data, sdscale0.asap, is artificial data with the following
    status:

       - nrow = 2
       - nchan = 4
       - all spectral values are 1.0
       - all Tsys values are 1.0

    scaletsys=True should be working.
       
    """
    # Input and output names
    rawfile='sdscale0.asap'
    prefix=sdscaleold_unittest_base.taskname+'Test1'
    outfile=prefix+'.asap'
    
    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdscaleold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test100(self):
        """Test 100: scalar factor with Tsys scaling"""
        factor = 2.0
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test101(self):
        """Test 101: scalar factor without Tsys scaling"""
        factor = 2.0
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test102(self):
        """Test 102: 1D array factor with Tsys scaling"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test103(self):
        """Test 103: 1D array factor without Tsys scaling"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test104(self):
        """Test 104: 2D array ([nrow,1]) factor with Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test105(self):
        """Test 105: 2D array ([nrow,1]) factor without Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test106(self):
        """Test 106: 2D array ([nrow,nchan]) factor with Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0,5.0,7.0,9.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test107(self):
        """Test 107: 2D array ([nrow,nchan]) factor without Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0,5.0,7.0,9.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

###
# Test on scaling 2
###
class sdscaleold_test2(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on actual scaling

    Test data, sdscale1.asap, is artificial data with the following
    status:

       - nrow = 2
       - nchan = 4 
       - 1 Tsys value for each spectrum
       - all spectral values are 1.0
       - all Tsys values are 1.0

    scaletsys=True should be working.
       
    """
    # Input and output names
    rawfile='sdscale1.asap'
    prefix=sdscaleold_unittest_base.taskname+'Test2'
    outfile=prefix+'.asap'
    
    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdscaleold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test200(self):
        """Test 200: scalar factor with Tsys scaling"""
        factor = 2.0
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test201(self):
        """Test 201: scalar factor without Tsys scaling"""
        factor = 2.0
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)
        
    def test202(self):
        """Test 202: 1D array factor with Tsys scaling (must fail)"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=True
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('SPECTRA and TSYS must conform in shape if you want to apply operation on Tsys.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test203(self):
        """Test 203: 1D array factor without Tsys scaling"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test204(self):
        """Test 204: 2D array ([nrow,1]) factor with Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test205(self):
        """Test 205: 2D array ([nrow,1]) factor without Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test206(self):
        """Test 206: 2D array ([nrow,nchan]) factor with Tsys scaling (must fail)"""
        factor = [[2.0,4.0,6.0,8.0],[3.0,5.0,7.0,9.0]]
        scaletsys=True
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('SPECTRA and TSYS must conform in shape if you want to apply operation on Tsys.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test207(self):
        """Test 207: 2D array ([nrow,nchan]) factor without Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0,5.0,7.0,9.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

###
# Test on scaling 3
###
class sdscaleold_test3(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on actual scaling

    Test data, sdscale2.asap, is artificial data with the following
    status:

       - nrow = 2
       - nchan = 4 in row0, 1 in row1
       - all spectral values are 1.0
       - all Tsys values are 1.0

    scaletsys=True should be working.
       
    """
    # Input and output names
    rawfile='sdscale2.asap'
    prefix=sdscaleold_unittest_base.taskname+'Test3'
    outfile=prefix+'.asap'
    
    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdscaleold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test300(self):
        """Test 300: scalar factor with Tsys scaling"""
        factor = 2.0
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test301(self):
        """Test 301: scalar factor without Tsys scaling"""
        factor = 2.0
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)
        
    def test302(self):
        """Test 302: 1D array factor with Tsys scaling (must fail)"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=True
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('ArrayColumn::getColumn cannot be done for column SPECTRA; the array shapes vary: Table array conformance error')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test303(self):
        """Test 303: 1D array factor without Tsys scaling (must fail)"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=False
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('All spectra in the input scantable must have the same number of channel for vector operation.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test304(self):
        """Test 304: 2D array ([nrow,1]) factor with Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test305(self):
        """Test 305: 2D array ([nrow,1]) factor without Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test306(self):
        """Test 306: 2D array ([nrow,nchan]) factor with Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test307(self):
        """Test 307: 2D array ([nrow,nchan]) factor without Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

###
# Test on scaling 4
###
class sdscaleold_test4(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on actual scaling

    Test data, sdscale3.asap, is artificial data with the following
    status:

       - nrow = 2
       - nchan = 4 in row0, 1 in row1
       - 1 Tsys value for each spectrum
       - all spectral values are 1.0
       - all Tsys values are 1.0

    scaletsys=True should be working.
       
    """
    # Input and output names
    rawfile='sdscale3.asap'
    prefix=sdscaleold_unittest_base.taskname+'Test4'
    outfile=prefix+'.asap'
    
    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdscaleold)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def test400(self):
        """Test 400: scalar factor with Tsys scaling"""
        factor = 2.0
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test401(self):
        """Test 401: scalar factor without Tsys scaling"""
        factor = 2.0
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test402(self):
        """Test 402: 1D array factor with Tsys scaling (must fail)"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=True
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('ArrayColumn::getColumn cannot be done for column SPECTRA; the array shapes vary: Table array conformance error')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test403(self):
        """Test 403: 1D array factor without Tsys scaling (must fail)"""
        factor = [2.0,3.0,4.0,5.0]
        scaletsys=False
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('All spectra in the input scantable must have the same number of channel for vector operation.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))

    def test404(self):
        """Test 404: 2D array ([nrow,1]) factor with Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test405(self):
        """Test 405: 2D array ([nrow,1]) factor without Tsys scaling"""
        factor = [[2.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

    def test406(self):
        """Test 406: 2D array ([nrow,nchan]) factor with Tsys scaling (must fail)"""
        factor = [[2.0,4.0,6.0,8.0],[3.0]]
        scaletsys=True
        try:
            res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
            self.assertTrue(False,
                            msg='The task must throw exception')
        except Exception as e:
            pos=str(e).find('SPECTRA and TSYS must conform in shape if you want to apply operation on Tsys.')
            self.assertNotEqual(pos,-1,
                                msg='Unexpected exception was thrown: %s'%(str(e)))            

    def test407(self):
        """Test 407: 2D array ([nrow,nchan]) factor without Tsys scaling"""
        factor = [[2.0,4.0,6.0,8.0],[3.0]]
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res,None,
                         msg='Any error occurred during calibration')
        self._compare(self.outfile,self.rawfile,factor,scaletsys)

###
# Test on flag
###
class sdscaleold_testflag(unittest.TestCase,sdscaleold_unittest_base):
    """
    Test on flag information handling

    Test data, sdscale_flagtest.asap, is artificial data with the following
    status:

       - nrow = 4
       - nchan = 100
       - rows #0 and #1 are row-flagged
       - channels 5 to 9, 15 to 19, and 94 to 98 are flagged in rows #1 and #2
       - all spectral values are 1.0
       - all tsys values are 1.0

    if a spectrum is row-flagged, scaling must not be applied on the spectrum.
    as for spectra which are not row-flagged, scaling must be applied to the 
    spectra/tsys values for all channels regardless of channel-flag values. 
    """
    # Input and output names
    rawfile='sdscale_flagtest.asap'
    prefix=sdscaleold_unittest_base.taskname+'TestFlag'
    outfile=prefix+'.asap'
    factor = 2.0

    def setUp(self):
        self.res=None
        if (not os.path.exists(self.rawfile)):
            shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        default(sdscaleold)
        self._getinfo(self.rawfile)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

    def testflag01(self):
        """Testflag: verify proper handling of flag information for scaletsys=True"""
        scaletsys=True
        res=sdscaleold(infile=self.rawfile,factor=self.factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res, None, msg='Any error occurred during calibration')

        self._check_flags_no_change()
        self._check_values(scaletsys)

    def testflag02(self):
        """Testflag: verify proper handling of flag information for scaletsys=False"""
        scaletsys=False
        res=sdscaleold(infile=self.rawfile,factor=self.factor,scaletsys=scaletsys,outfile=self.outfile)
        self.assertEqual(res, None, msg='Any error occurred during calibration')

        self._check_flags_no_change()
        self._check_values(scaletsys)

    def _getinfo(self, infile):
        with tbmanager(infile) as tb:
            self.nchan_orig = len(tb.getcell('FLAGTRA', 0))
            self.rowid_rflag_orig = tb.getcol('FLAGROW')
            cfraw = tb.getcol('FLAGTRA').sum(axis=0)
            self.rowid_cflag_orig = [cfraw[i] > 0 for i in range(len(cfraw))]
            self.cflag_orig = tb.getcell('FLAGTRA', numpy.where(self.rowid_cflag_orig)[0][0])

    def _check_flags_no_change(self):
        """check if no changes applied on flag values"""
        with tbmanager(self.outfile) as tb:
            self.assertTrue(all(tb.getcol('FLAGROW')==self.rowid_rflag_orig))
            for i in range(tb.nrows()):
                chanflag = tb.getcell('FLAGTRA', i)
                chanflag_ref = self.cflag_orig if self.rowid_cflag_orig[i] else numpy.zeros(self.nchan_orig, numpy.int32)
                self.assertTrue(all(chanflag==chanflag_ref))
        
    def _check_values(self, scaletsys):
        """check spectra and tsys values"""
        with tbmanager(self.outfile) as tb:
            for irow, col in itertools.product(range(tb.nrows()), ['SPECTRA', 'TSYS']):
                data = tb.getcell(col, irow)
                data_ref = numpy.ones(self.nchan_orig, numpy.float)
                if self.rowid_rflag_orig[irow] == 0:
                    if scaletsys or col=='SPECTRA':
                        data_ref *= self.factor
                self.assertTrue(all(data==data_ref))


def suite():
    return [sdscaleold_test0,sdscaleold_test1,
            sdscaleold_test2,sdscaleold_test3,sdscaleold_test4,
            sdscaleold_testflag]
