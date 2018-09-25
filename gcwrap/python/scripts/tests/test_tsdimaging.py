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
import math
import re
import string

_ia = iatool( )
_rg = rgtool( )

try:
    from . import selection_syntax
except:
    import tests.selection_syntax as selection_syntax

try:
    from . import testutils
except:
    import tests.testutils as testutils

from tsdimaging import tsdimaging as sdimaging
from sdutil import tbmanager, toolmanager, table_selector
from task_tsdimaging import image_suffix

#
# Unit test of sdimaging task.
# 
def construct_refstat_uniform(fluxval, blc_data, trc_data):
    """
    Return a dictionary of analytic reference statistics of uniform image
    Arguments:
        fluxval  : the uniform flux of the image
        blc_data : blc of un-masked pixel (e.g., [0,0,0,0] for whole image)
        trc_data : trc of un-masked pixel
    Returns:
    a dictionary of statistics, 'min', 'max', 'rms', 'sigma', 'mean',
    'npts', 'sum', and 'sumsq'
    """
    # the number of valid (unmasked) pixels
    ndim = len(blc_data)
    nvalid = 1
    for idim in range(ndim):
        nvalid *= abs(trc_data[idim]-blc_data[idim]+1)
    retstat = {'min': [fluxval], 'max': [fluxval], 'rms': [fluxval],
               'sigma': [0.], 'mean': [fluxval], 'npts': [nvalid],
               'sum': [fluxval*nvalid], 'sumsq': [fluxval**2*nvalid]}
    return retstat

def merge_dict(d1, d2):
    """
    Out of place merge of two dictionaries.
    If both dictionary has the same keys, value of the second
    dictionary is adopted.
    """
    if type(d1) != dict or type(d2) != dict:
        raise ValueError("Internal error. inputs should be dictionaries.")
    d12 = d1.copy()
    d12.update(d2)
    return d12

def get_table_cache():
    (mytb,) = gentools(['tb'])
    cache = mytb.showcache()
    #print 'cache = {}'.format(cache)
    return cache

###
# Base class for sdimaging unit test
###
class sdimaging_standard_paramset(object):
    rawfile='sdimaging.ms'
    phasecenter='J2000 17:18:29 +59.31.23'
    imsize=[75,75]
    cell=['3.0arcmin','3.0arcmin']
    gridfunction='PB'
    minweight0 = 0.
    mode='channel'
    nchan=40
    start=400
    width=10    

class sdimaging_unittest_base(unittest.TestCase, sdimaging_standard_paramset):
    """
    Base class for sdimaging unit test

    Test data is originally FLS3_all_newcal_SP and created
    by the following script:

    asap_init()
    sd.rc('scantable',storage='disk')
    s=sd.scantable('FLS3_all_newcal_SP',average=False,getpt=False)
    s0=s.get_scan('FLS3a')
    s0.save('FLS3a_HI.asap')
    del s,s0
    s=sd.scantable('FLS3a_HI.asap',average=False)
    s.set_fluxunit('K')
    scannos=s.getscannos()
    res=sd.calfs(s,scannos)
    del s
    res.save('FLS3a_calfs','MS2')
    tb.open('FLS3a_calfs')
    tbsel=tb.query('SCAN_NUMBER<100')
    tbc=tbsel.copy('sdimaging.ms',deep=True,valuecopy=True)
    tbsel.close()
    tb.close()
    tbc.close()

    Furthermore, SYSCAL and POINTING tables are downsized.
    
    """
    taskname='sdimaging'
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdimaging/'
    #rawfile='sdimaging.ms'
    postfix='.im'
    ms_nchan = 1024
#     phasecenter='J2000 17:18:29 +59.31.23'
#     imsize=[75,75]
#     cell=['3.0arcmin','3.0arcmin']
#     gridfunction='PB'
#     minweight0 = 0.
    statsinteg={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                'blcf': '17:32:18.690, +57.37.28.536, I, 1.42064e+09Hz',
                'max': numpy.array([ 0.6109162]),
                'maxpos': numpy.array([4, 62,  0,  0], dtype=numpy.int32),
                'maxposf': '17:31:59.439, +60.43.52.421, I, 1.42064e+09Hz',
                'mean': numpy.array([ 0.39524983]),
                'min': numpy.array([ 0.]),
                'minpos': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                'minposf': '17:32:18.690, +57.37.28.536, I, 1.42064e+09Hz',
                'npts': numpy.array([ 5625.]),
                'rms': numpy.array([ 0.43127564]),
                'sigma': numpy.array([ 0.17257331]),
                'sum': numpy.array([ 2223.28028646]),
                'sumsq': numpy.array([ 1046.2425779]),
                'trc': numpy.array([74, 74,  0,  0], dtype=numpy.int32),
                'trcf': '17:03:03.151, +61.19.10.757, I, 1.42064e+09Hz'}
    keys=['max','mean','min','npts','rms','blc','blcf','trc','trcf','sigma','sum','sumsq']

    def run_test_common(self, task_param, refstats, shape, refbeam=None,
                        atol=1.e-8, rtol=1.e-5, compstats=None, ignoremask=True):
        """
        Run sdimaging and test results.
        A list of tests:
        (1) task completes without an error
        (2) ouput image and weight image exist
        (3) image shape
        (4) image statistics
        (5) reference beam of image (optional)
        """
        res=sdimaging(**task_param)
        outprefix = task_param['outfile']
        outfile = outprefix + image_suffix
        # Tests
        self.assertEqual(res,None,
                         msg='Any error occurred during imaging')
        self._checkfile(outfile)
        self._checkfile(outprefix+".weight")
        self._checkframe(outfile)
        self._checkshape(outfile, shape[0], shape[1],shape[2],shape[3])
        self._checkstats(outfile, refstats, compstats=compstats,
                         atol=atol, rtol=rtol, ignoremask=ignoremask)
        if refbeam is not None:
            self._check_beam(outfile, refbeam)

    def _checkfile( self, name ):
        isthere=os.path.exists(name)
        self.assertEqual(isthere,True,
                         msg='output file %s was not created because of the task failure'%(name))
        
    def _checkframe(self, name):
        _ia.open(name)
        csys = _ia.coordsys()
        _ia.close()
        spectral_frames = csys.referencecode('spectral')
        csys.done()
        self.assertEqual(1, len(spectral_frames))
        spectral_frame = spectral_frames[0]
        self.assertEqual('LSRK', spectral_frame)

    def _checkshape(self,name,nx,ny,npol,nchan):
        self._checkfile(name)
        _ia.open(name)
        imshape=_ia.shape()
        _ia.close()
        self.assertEqual(nx,imshape[0],
                    msg='nx does not match')
        self.assertEqual(ny,imshape[1],
                    msg='ny does not match')
        self.assertEqual(npol,imshape[2],
                    msg='npol does not match')
        self.assertEqual(nchan,imshape[3],
                    msg='nchan does not match')

    def _checkstats(self,name, ref, compstats=None, atol=1.e-8, rtol=1.e-5, region=None, ignoremask=False):
        """
        A test function to compare statistics of an image with reference
        values.
        Arguments:
            name  :  name of an image to test statistics
            ref   :  a record (dictionary) of the reference statistic values
            compstats : a list of names of statistis to compare. By default,
                        the list is taken from all keys in ref
            atol  : absolute tolerance (see help in numpy.allclose)
            rtol  : relative tolerance (see help in numpy.allclose)
            region : region of image to calculate statistics. a CASA region
                     record should be specified (see help of , e.g., rg.box).
                     default is whole image.
            ignoremask : when True, mask in image is ignored and statistics
                         are calculated from whole pixels in image. default
                         is False (take image mask into account).
        """
        self._checkfile(name)
        if compstats is None: compstats = list(ref.keys())
        if region is None: region = ""
        _ia.open(name)
        try:
            if ignoremask:
                def_mask = _ia.maskhandler('default')
                _ia.calcmask('T')
            stats=_ia.statistics(region=region, list=True, verbose=True)
            if ignoremask:
                _ia.maskhandler('set',def_mask)
        except: raise
        finally: _ia.close()
        #for key in stats.keys():
        for key in compstats:
            message='statistics \'%s\' does not match: %s (expected: %s)' % ( key, str(stats[key]), str(ref[key]) )
            if type(stats[key])==str:
                self.assertEqual(stats[key],ref[key],
                                 msg=message)
            else:
                #print stats[key]-ref[key]
                ret=numpy.allclose(stats[key],ref[key], atol=atol, rtol=rtol)
                self.assertEqual(ret,True,
                                 msg=message)

    def _checkdirax(self, imagename, center, cell, imsize):
        """ Test image center, cell size and imsize"""
        cell = self._format_dir_list(cell)
        imsize = self._format_dir_list(imsize)
        _ia.open(imagename)
        csys = _ia.coordsys()
        ret = _ia.summary()
        _ia.close()
        ra_idx = csys.findaxisbyname('ra')
        dec_idx = csys.findaxisbyname('dec')
        ra_unit = ret['axisunits'][ra_idx]
        dec_unit = ret['axisunits'][dec_idx]
        # imsize
        self.assertEqual(imsize[0], ret['shape'][ra_idx],\
                         msg="nx = %d (expected: %d)" % \
                         (imsize[0], ret['shape'][ra_idx]))
        self.assertEqual(imsize[1], ret['shape'][dec_idx],\
                         msg="nx = %d (expected: %d)" % \
                         (imsize[1], ret['shape'][dec_idx]))
        # image center
        tol = "1arcsec"
        cen_arr = center.split()
        cen_ref = me.direction(*cen_arr)
        cen_x = (qa.convert(cen_ref['m0'], 'rad')['value'] % (numpy.pi*2))
        cen_y = qa.convert(cen_ref['m1'], 'rad')['value']
        ref_x = qa.convert(qa.quantity(ret['refval'][ra_idx],ra_unit),'rad')['value']
        ref_x = (ref_x % (numpy.pi*2))
        ref_y = qa.convert(qa.quantity(ret['refval'][dec_idx],dec_unit),'rad')['value']
        tol_val = qa.convert(tol, 'rad')['value']
        self.assertTrue(abs(ref_x-cen_x) < tol_val,
                        msg="center_x = %f %s (expected: %f)" % \
                        (ref_x, ra_unit, cen_x))
        self.assertTrue(abs(ref_y-cen_y) < tol_val,
                        msg="center_y = %f %s (expected: %f)" % \
                        (ref_x, ra_unit, cen_x))
        
        # cell (imager seems to set negative incr for dx)
        dx = - qa.convert(cell[0], ra_unit)['value']
        dy = qa.convert(cell[1], dec_unit)['value']
        incx = ret['incr'][ra_idx]
        incy = ret['incr'][dec_idx]
        self.assertAlmostEqual((incx-dx)/dx, 0., places=5, \
                               msg="cellx = %f %s (expected: %f)" % \
                               (incx, ra_unit, dx))
        self.assertAlmostEqual((incy-dy)/dy, 0., places=5, \
                               msg="celly = %f %s (expected: %f)" % \
                               (incy, dec_unit, dy))

    def _format_dir_list(self, inval):
        if type(inval) == str:
            return [inval, inval]
        elif len(inval) == 1:
            return [inval[0], inval[0]]
        return inval[0:2]

    def _check_beam(self, image, ref_beam):
        """Check image beam size"""
        _ia.open(image)
        beam = _ia.restoringbeam()
        _ia.close()
        maj_asec = qa.getvalue(qa.convert(beam['major'], 'arcsec'))
        min_asec = qa.getvalue(qa.convert(beam['minor'], 'arcsec'))
        maj_asec_ref = qa.getvalue(qa.convert(ref_beam['major'], 'arcsec'))
        min_asec_ref = qa.getvalue(qa.convert(ref_beam['minor'], 'arcsec'))
        self.assertAlmostEqual(abs(maj_asec-maj_asec_ref)/max(maj_asec_ref,1.e-12), 0., places=3, msg="major axis = %f arcsec (expected: %f)" % (maj_asec, maj_asec_ref))
        self.assertAlmostEqual(abs(min_asec-min_asec_ref)/max(min_asec_ref,1.e-12), 0., places=3, msg="minor axis = %f arcsec (expected: %f)" % (min_asec, min_asec_ref))
    
    def _check_restfreq(self, imagename, restfreq):
        """ Test image rest frequency"""
        self.assertTrue(qa.compare(restfreq, 'Hz'))
        myunit = qa.getunit(restfreq)
        refval = qa.getvalue(restfreq)[0]
        _ia.open(imagename)
        csys = _ia.coordsys()
        _ia.close()
        testval = qa.getvalue(qa.convert(csys.restfrequency(), myunit))
        csys.done()
        ret=numpy.allclose(testval,refval, atol=1.e-5, rtol=1.e-5)
        self.assertTrue(ret)
        
        
###
# Test on bad parameter settings
###
class sdimaging_test0(sdimaging_unittest_base):
    """
    Test on bad parameter setting
    """
    # Input and output names
    prefix = sdimaging_unittest_base.taskname+'Test0'
    badid = '99'
    outfile = prefix+sdimaging_unittest_base.postfix

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)

        default(sdimaging)
        self.task_param = dict(infiles=self.rawfile,mode='channel',
                               outfile=self.outfile,intent='',
                               cell=self.cell,imsize=self.imsize,
                               phasecenter=self.phasecenter,
                               gridfunction=self.gridfunction,
                               nchan=40,start=400,
                               width=10,
                               minweight=self.minweight0)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )
        
        self.assertEqual(len(get_table_cache()), 0)
        
    def run_exception_case(self, task_param, expected_msg, expected_type=RuntimeError):
        with self.assertRaises(RuntimeError) as cm:
            res=sdimaging(**task_param)
        the_exception = cm.exception
        pos=str(the_exception).find(expected_msg)
        self.assertNotEqual(pos,-1,
                            msg='Unexpected exception was thrown: {0}'.format(str(the_exception)))
        
    def test000(self):
        """Test 000: Default parameters"""
        # argument verification error
        res=sdimaging()
        self.assertFalse(res)

    def test001(self):
        """Test001: Bad mode"""
        # argument verification error
        res=sdimaging(infiles=self.rawfile,mode='badmode',intent='',outfile=self.outfile)
        self.assertFalse(res)

    def test002(self):
        """Test002: Bad field id"""
        self.task_param['field'] = self.badid
        msg = 'Field Expression: Partial or no match for Field ID list [{0}]'.format(self.badid)
        self.run_exception_case(self.task_param, msg)

    def test003(self):
        """Test003: Bad spectral window id"""
        self.task_param['spw'] = self.badid
        msg = 'Spw Expression: No match found for {0}'.format(self.badid)
        self.run_exception_case(self.task_param, msg)

    def test004(self):
        """Test004: Bad antenna id"""
        self.task_param['antenna'] = self.badid
        msg = 'No match found for the antenna specificion'
        self.run_exception_case(self.task_param, msg)
        
    def test005(self):
        """Test005: Bad stokes parameter"""
        self.task_param['stokes'] = 'BAD'
        # argument verification error
        res = sdimaging(**self.task_param)
        self.assertFalse(res)
        
    def test006(self):
        """Test006: Bad gridfunction"""
        # argument verification error
        res=sdimaging(infiles=self.rawfile,gridfunction='BAD',intent='',outfile=self.outfile)
        self.assertFalse(res)

    def test007(self):
        """Test007: Bad scanlist"""
        self.task_param['scan'] = self.badid
        msg = 'Data selection ended with 0 rows'
        self.run_exception_case(self.task_param, msg)

    def test008(self):
        """Test008: Existing outfile with overwrite=False"""
        outfile = self.outfile + image_suffix
        f=open(outfile, 'w')
        print('existing file', file=f)
        f.close()
        self.task_param['overwrite'] = False
        msg = 'Output file \'{0}\' exists.'.format(outfile)
        self.run_exception_case(self.task_param, msg)

    def test009(self):
        """Test009: Bad phasecenter string"""
        self.task_param['phasecenter'] = 'This is bad'
        msg = 'Error in converting \'{0}\' to MDirection.'.format(self.task_param['phasecenter'])
        self.run_exception_case(self.task_param, msg)

    def test010(self):
        """Test010: Bad phasecenter reference (CHANGED: raise an error)"""
        # older sdimaging was so kind that it assumed J2000 when unrecognized direction frame was given
        # in the new tsdimaging raises an error in such case
        false_phasecenter = self.phasecenter.replace('J2000', 'J3000')
        self.task_param['phasecenter'] = false_phasecenter
        msg = 'Invalid Image Parameter set : Error in converting \'{0}\' to MDirection.'.format(false_phasecenter)
        self.run_exception_case(self.task_param, msg)
#         # default for unknown direction frame is J2000 
#         refimage=self.outfile+'2'
#         sdimaging(infiles=self.rawfile,outfile=self.outfile,intent='',cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter.replace('J2000','J3000'),minweight=self.minweight0)
#         sdimaging(infiles=self.rawfile,outfile=refimage,intent='',cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter,minweight=self.minweight0)
#         tb.open(self.outfile)
#         chunk=tb.getcol('map')
#         tb.close()
#         tb.open(refimage)
#         refchunk=tb.getcol('map')
#         tb.close()
#         ret=all(chunk.flatten()==refchunk.flatten())
#         #print ret
#         self.assertTrue(ret)

    def test011(self):
        """Test011: Bad pointingcolumn name"""
        # argument verification error
        res=sdimaging(infiles=self.rawfile,outfile=self.outfile,intent='',cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter,pointingcolumn='non_exist')
        self.assertFalse(res)

    def test012(self):
        """Test012: Bad imsize"""
        self.task_param['imsize'] = [1,0]
        msg = 'Error in building Coordinate System and Image Shape : Internal Error : Image shape is invalid :'
        self.run_exception_case(self.task_param, msg)

    def test013(self):
        """Test013: Bad cell size"""
        self.task_param['cell'] = [0., 0.]
        msg = 'Error in building Coordinate System and Image Shape : wcs wcsset_error: Linear transformation matrix is singular'
        self.run_exception_case(self.task_param, msg)

    def test014(self):
        """Test014: Too fine resolution (smaller than original channel width"""
        specunit = 'GHz'
        start = '%f%s' % (1.4202, specunit)
        width = '%e%s' % (1.0e-10, specunit)
        self.task_param['mode'] = 'frequency'
        self.task_param['start'] = start
        self.task_param['width'] = width
        msg = 'calcChanFreqs failed, check input start and width parameters'
        self.run_exception_case(self.task_param, msg)

    def test015(self):
        """Test015: negative minweight"""
        res=sdimaging(infiles=self.rawfile,outfile=self.outfile,intent='',cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter,minweight=-1.)
        self.assertFalse(res)


###
# Test channel imaging
###
class sdimaging_test1(sdimaging_unittest_base):
    """
    Test channel imaging

       - integrated image
       - full channel image
       - selected channel image
       - BOX and SF imaging (default is PB)
       - two polarization imaging (XX and YY, default is Stokes I)
       - empty phasecenter
       - settting minweight = 0.2
       
    """
    # Input and output names
    prefix=sdimaging_unittest_base.taskname+'Test1'
    outfile=prefix+sdimaging_unittest_base.postfix
#     mode='channel'
#     nchan=40
#     start=400
#     width=10

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        # Common task parameters of the class
        self.task_param = dict(infiles=self.rawfile,mode=self.mode,
                               spw='0',
                               outfile=self.outfile,intent='',
                               cell=self.cell,imsize=self.imsize,
                               phasecenter=self.phasecenter,
                               gridfunction=self.gridfunction,
                               nchan=self.nchan,start=self.start,
                               width=self.width,
                               minweight=self.minweight0)

        default(sdimaging)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)

    def test100(self):
        """Test 100: Integrated image"""
        self.task_param.update(dict(nchan=1,start=0,width=self.ms_nchan))
        outshape = (self.imsize[0],self.imsize[1],1,1)
        self.run_test_common(self.task_param, self.statsinteg, outshape, compstats=self.keys,
                             ignoremask=True)

    def test101(self):
        """Test 101: Full channel image (nchan = -1)"""
        self.task_param.update(dict(nchan=-1,start="",width=""))
        outshape = (self.imsize[0],self.imsize[1],1,self.ms_nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.419395e+09Hz',
                  'max': numpy.array([ 24.77152824]),
                  'maxpos': numpy.array([ 59,  21,   0, 605], dtype=numpy.int32),
                  'maxposf': '17:10:00.642, +58.42.19.808, I, 1.420872e+09Hz',
                  'mean': numpy.array([ 0.39542111]),
                  'min': numpy.array([-1.84636593]),
                  'minpos': numpy.array([  73,    6,    0, 1023], dtype=numpy.int32),
                  'minposf': '17:04:54.966, +57.55.36.907, I, 1.421893e+09Hz',
                  'npts': numpy.array([ 5760000.]),
                  'rms': numpy.array([ 1.01357317]),
                  'sigma': numpy.array([ 0.93325921]),
                  'sum': numpy.array([ 2277625.60731485]),
                  'sumsq': numpy.array([ 5917423.42281288]),
                  'trc': numpy.array([  74,   74,    0, 1023], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.421893e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape, compstats=self.keys,
                             ignoremask=True)

        
    def test102(self):
        """Test 102: Full channel image"""
        tb.open(self.rawfile)
        if 'FLOAT_DATA' in tb.colnames():
            nchan=tb.getcell('FLOAT_DATA').shape[1]
        else:
            nchan=tb.getcell('DATA').shape[1]
        tb.close()
        self.task_param.update(dict(nchan=nchan,start=0,width=1))
        # for testing 
        #self.task_param['gridfunction'] = 'BOX'
        for (k,v) in self.task_param.items():
            casalog.post('test102: {0} = \'{1}\' (type {2})'.format(k,v,type(v)))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.419395e+09Hz',
                  'max': numpy.array([ 24.77152824]),
                  'maxpos': numpy.array([ 59,  21,   0, 605], dtype=numpy.int32),
                  'maxposf': '17:10:00.642, +58.42.19.808, I, 1.420872e+09Hz',
                  'mean': numpy.array([ 0.39542111]),
                  'min': numpy.array([-1.84636593]),
                  'minpos': numpy.array([  73,    6,    0, 1023], dtype=numpy.int32),
                  'minposf': '17:04:54.966, +57.55.36.907, I, 1.421893e+09Hz',
                  'npts': numpy.array([ 5760000.]),
                  'rms': numpy.array([ 1.01357317]),
                  'sigma': numpy.array([ 0.93325921]),
                  'sum': numpy.array([ 2277625.60731485]),
                  'sumsq': numpy.array([ 5917423.42281288]),
                  'trc': numpy.array([  74,   74,    0, 1023], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.421893e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape, compstats=self.keys,
                             ignoremask=True)

    def test103(self):
        """Test 103: Selected channel image"""
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max': numpy.array([ 14.79568005]),
                  'maxpos': numpy.array([57, 20,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:47.496, +58.39.30.813, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.82293006]),
                  'min': numpy.array([-0.08763941]),
                  'minpos': numpy.array([61, 71,  0, 35], dtype=numpy.int32),
                  'minposf': '17:08:30.980, +61.12.02.893, I, 1.42124e+09Hz',
                  'npts': numpy.array([ 225000.]),
                  'rms': numpy.array([ 1.54734671]),
                  'sigma': numpy.array([ 1.31037237]),
                  'sum': numpy.array([ 185159.263672]),
                  'sumsq': numpy.array([ 538713.45272028]),
                  'trc': numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        # beam size from r32523
        ref_beam=dict(major='661.858412arcsec',minor='661.858412arcsec')
        self.run_test_common(self.task_param, refstats, outshape,  refbeam=ref_beam,
                             compstats=self.keys, ignoremask=True)

    def test104(self):
        """Test 104: Box-car gridding"""
        self.task_param.update(dict(gridfunction='BOX'))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max': numpy.array([ 15.64525127]),
                  'maxpos': numpy.array([58, 20,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:24.433, +58.39.25.476, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.66097592]),
                  'min': numpy.array([-0.42533547]),
                  'minpos': numpy.array([69, 62,  0, 38], dtype=numpy.int32),
                  'minposf': '17:05:23.086, +60.44.01.427, I, 1.42131e+09Hz',
                  'npts': numpy.array([ 225000.]),
                  'rms': numpy.array([ 1.38591599]),
                  'sigma': numpy.array([ 1.2181464]),
                  'sum': numpy.array([ 148719.58227018]),
                  'sumsq': numpy.array([ 432171.72687429]),
                  'trc': numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        # beam size from r32523
        ref_beam=dict(major='503.181345arcsec',minor='503.181345arcsec')
        self.run_test_common(self.task_param, refstats, outshape, refbeam=ref_beam,
                             compstats=self.keys, ignoremask=True)

    def test105(self):
        """Test 105: Prolate Spheroidal gridding"""
        self.task_param.update(dict(gridfunction='SF'))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max': numpy.array([ 15.13189793]),
                  'maxpos': numpy.array([58, 21,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:23.737, +58.42.25.413, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.773227]),
                  'min': numpy.array([-0.07284018]),
                  'minpos': numpy.array([ 5, 67,  0, 30], dtype=numpy.int32),
                  'minposf': '17:31:41.090, +60.59.00.556, I, 1.42112e+09Hz',
                  'npts': numpy.array([ 225000.]),
                  'rms': numpy.array([ 1.49926317]),
                  'sigma': numpy.array([ 1.28449107]),
                  'sum': numpy.array([ 173976.07570213]),
                  'sumsq': numpy.array([ 505752.74505987]),
                  'trc': numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        # beam size from analysisUtils.
        #aU.sfBeam(1.42038,diameter=104.9,pixelsize=180.,
        #          xSamplingArcsec=354.16985191848795,
        #          ySamplingArcsec=180.0432853343201,
        #          convsupport=3,obscuration=0.0)
        ref_beam=dict(major='618.853892arcsec',minor='618.853892arcsec')
        self.run_test_common(self.task_param, refstats, outshape,  refbeam=ref_beam,
                             compstats=self.keys, ignoremask=True)

    def test106(self):
        """Test 106: Imaging two polarization separately (XX and YY, not Stokes I)"""
        self.task_param.update(dict(stokes='XXYY',gridfunction='PB'))
        outshape = (self.imsize[0],self.imsize[1],2,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, XX, 1.42038e+09Hz',
                  'max': numpy.array([ 15.057868]),
                  'maxpos': numpy.array([57, 20,  1, 20], dtype=numpy.int32),
                  'maxposf': '17:10:47.496, +58.39.30.813, YY, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.82292841]),
                  'min': numpy.array([-0.41953856]),
                  'minpos': numpy.array([10,  3,  1, 31], dtype=numpy.int32),
                  'minposf': '17:28:37.170, +57.47.49.422, YY, 1.42114e+09Hz',
                  'npts': numpy.array([ 450000.]),
                  'rms': numpy.array([ 1.55436146]),
                  'sigma': numpy.array([ 1.31864787]),
                  'sum': numpy.array([ 370317.78554221]),
                  'sumsq': numpy.array([ 1087217.77687839]),
                  'trc': numpy.array([74, 74,  1, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, YY, 1.42133e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape, compstats=self.keys,
                             ignoremask=True)

    def test107(self):
        """Test 107: Gaussian gridding"""
        self.task_param.update(dict(gridfunction='GAUSS'))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max': numpy.array([ 15.28046036]),
                  'maxpos': numpy.array([58, 21,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:23.737, +58.42.25.413, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.75082603]),
                  'min': numpy.array([-0.14009152]),
                  'minpos': numpy.array([34, 69,  0, 33], dtype=numpy.int32),
                  'minposf': '17:19:43.545, +61.07.22.487, I, 1.42119e+09Hz',
                  'npts': numpy.array([ 225000.]),
                  'rms': numpy.array([ 1.47686982]),
                  'sigma': numpy.array([ 1.2717751]),
                  'sum': numpy.array([ 168935.85698331]),
                  'sumsq': numpy.array([ 490757.49952306]),
                  'trc': numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        # beam size from r32523
        ref_beam=dict(major='510.142597arcsec',minor='510.142597arcsec')
        self.run_test_common(self.task_param, refstats, outshape, refbeam=ref_beam,
                             compstats=self.keys, ignoremask=True)

    def test108(self):
        """Test 108: Gaussian*Jinc gridding"""
        self.task_param.update(dict(gridfunction='GJINC'))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc':numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max':numpy.array([ 15.31498909]),
                  'maxpos':numpy.array([58, 21,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:23.737, +58.42.25.413, I, 1.42087e+09Hz',
                  'mean':numpy.array([ 0.72415226]),
                  'min':numpy.array([-0.16245638]),
                  'minpos':numpy.array([68, 69,  0, 36], dtype=numpy.int32),
                  'minposf': '17:05:39.206, +61.05.09.055, I, 1.42126e+09Hz',
                  'npts':numpy.array([ 225000.]),
                  'rms':numpy.array([ 1.44985926]),
                  'sigma':numpy.array([ 1.25606618]),
                  'sum':numpy.array([ 162934.25891985]),
                  'sumsq':numpy.array([ 472970.63791706]),
                  'trc':numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        # beam size from analysisUtils.
        #aU.gjincBeam(1.42038,diameter=104.9,pixelsize=180.,geometricMean=True,
        #             xSamplingArcsec=354.16985191848795,
        #             ySamplingArcsec=180.0432853343201,
        #             obscuration=0.0,widthMultiplier=1.0)
        ref_beam=dict(major='580.094135arcsec',minor='580.094135arcsec')
        self.run_test_common(self.task_param, refstats, outshape, refbeam=ref_beam,
                             compstats=self.keys, ignoremask=True)

    def test109(self):
        """Test 109: Empty phasecenter (auto-calculation)"""
        self.task_param.update(dict(phasecenter="",gridfunction="BOX"))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:31:48.220, +57.36.09.784, I, 1.42038e+09Hz',
                  'max': numpy.array([ 15.64525127]),
                  'maxpos': numpy.array([57, 20,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:17.816, +58.38.11.961, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.66039867]),
                  'min': numpy.array([-0.42533547]),
                  'minpos': numpy.array([68, 63,  0, 38], dtype=numpy.int32),
                  'minposf': '17:05:16.976, +60.45.51.215, I, 1.42131e+09Hz',
                  'npts': numpy.array([ 225000.]),
                  'rms': numpy.array([ 1.38517249]),
                  'sigma': numpy.array([ 1.21761365]),
                  'sum': numpy.array([ 148589.70138012]),
                  'sumsq': numpy.array([ 431708.13145918]),
                  'trc': numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:02:33.828, +61.17.52.040, I, 1.42133e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)

    def test110(self):
        """Test 110: setting minweight=70."""
        self.task_param.update(dict(gridfunction='GJINC', minweight=70.))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats={'blc':numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.42038e+09Hz',
                  'max':numpy.array([ 15.31498909]),
                  'maxpos':numpy.array([58, 21,  0, 20], dtype=numpy.int32),
                  'maxposf': '17:10:23.737, +58.42.25.413, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.96643395]),
                  'min': numpy.array([-0.01385191]),
                  'minpos': numpy.array([19, 63,  0, 33], dtype=numpy.int32),
                  'minposf': '17:25:51.974, +60.48.38.410, I, 1.42119e+09Hz',
                  'npts': numpy.array([ 143920.]),
                  'rms': numpy.array([ 1.66819704]),
                  'sigma': numpy.array([ 1.35974246]),
                  'sum': numpy.array([ 139089.17359187]),
                  'sumsq': numpy.array([ 400512.27532199]),
                  'trc':numpy.array([74, 74,  0, 39], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42133e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=False)


    def test111(self):
        """imsize in float (ntegrated image)"""
        outshape = (self.imsize[0],self.imsize[1],1,1)
        imsize = [ float(v) for v in self.imsize ]
        self.task_param.update(dict(nchan=1,start=0,width=self.ms_nchan,
                                    imsize=imsize))
        self.run_test_common(self.task_param, self.statsinteg, outshape, compstats=self.keys,
                             ignoremask=True)


    def test112(self):
        """round-up imsize in float (integrated image)"""
        outshape = (self.imsize[0],self.imsize[1],1,1)
        imsize = [ float(v)-0.8 for v in self.imsize ]
        self.task_param.update(dict(nchan=1,start=0,width=self.ms_nchan,
                                    imsize=imsize))
        self.run_test_common(self.task_param, self.statsinteg, outshape, compstats=self.keys,
                             ignoremask=True)


###
# Test frequency imaging
###
class sdimaging_test2(sdimaging_unittest_base):
    """
    Test frequency imaging

       - integrated image
       - selected frequency image
       
    """
    # Input and output names
    prefix=sdimaging_unittest_base.taskname+'Test2'
    outfile=prefix+sdimaging_unittest_base.postfix
    unit='GHz'
    mode = "frequency"

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        # Common task parameters of the class
        self.task_param = dict(infiles=self.rawfile,mode=self.mode,
                               outfile=self.outfile,intent='',
                               cell=self.cell,imsize=self.imsize,
                               phasecenter=self.phasecenter,
                               gridfunction=self.gridfunction,
                               minweight=self.minweight0)

        default(sdimaging)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)

    def test200(self):
        """Test 200: Integrated image"""
        nchan = 1
        ms.open(self.rawfile)
        spwinfo =  ms.getspectralwindowinfo()
        ms.close()
        spwid0 = list(spwinfo.keys())[0]
        start = '%fHz' % (spwinfo[spwid0]['Chan1Freq']+0.5*(spwinfo[spwid0]['TotalWidth']-spwinfo[spwid0]['ChanWidth'])) 
        width = '%fHz' % (spwinfo[spwid0]['TotalWidth'])
        self.task_param.update(dict(nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        self.run_test_common(self.task_param, self.statsinteg, outshape,
                             compstats=self.keys, ignoremask=True)
        
    def test201(self):
        """Test 201: Full channel image (mode='frequency', nchan = -1)"""
        self.task_param.update(dict(nchan = -1, start = '', width = ''))
        # workaround for new imager framework
        # New imager looks SPECTRAL_WINDOW table to get whole frequency range 
        # regardless of whether associating data exist in the MAIN table or not.
        # As a result, resulting image has 2048 channels instead of 1024.
        # To proceed this test, spw is explicitly specified here.
        spw = '0'
        self.task_param.update(dict(spw=spw))
        outshape = (self.imsize[0],self.imsize[1],1,self.ms_nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.419395e+09Hz',
                  'max': numpy.array([ 24.77152824]),
                  'maxpos': numpy.array([ 59,  21,   0, 605], dtype=numpy.int32),
                  'maxposf': '17:10:00.642, +58.42.19.808, I, 1.420872e+09Hz',
                  'mean': numpy.array([ 0.39542111]),
                  'min': numpy.array([-1.84636593]),
                  'minpos': numpy.array([  73,    6,    0, 1023], dtype=numpy.int32),
                  'minposf': '17:04:54.966, +57.55.36.907, I, 1.421893e+09Hz',
                  'npts': numpy.array([ 5760000.]),
                  'rms': numpy.array([ 1.01357317]),
                  'sigma': numpy.array([ 0.93325921]),
                  'sum': numpy.array([ 2277625.60731485]),
                  'sumsq': numpy.array([ 5917423.42281288]),
                  'trc': numpy.array([  74,   74,    0, 1023], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.421893e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)

    def test202(self):
        """Test 202: Selected frequency image"""
        nchan = 100
        start = "%f%s" % (1.4202, self.unit)
        width = "%f%s" % (1.0e-5, self.unit)
        self.task_param.update(dict(nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.4202e+09Hz',
                  'max': numpy.array([ 21.55560875]),
                  'maxpos': numpy.array([59, 21,  0, 67], dtype=numpy.int32),
                  'maxposf': '17:10:00.642, +58.42.19.808, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.80467233]),
                  'min': numpy.array([-0.27736959]),
                  'minpos': numpy.array([58, 71,  0, 10], dtype=numpy.int32),
                  'minposf': '17:09:45.684, +61.12.21.875, I, 1.4203e+09Hz',
                  'npts': numpy.array([ 562500.]),
                  'rms': numpy.array([ 1.56429076]),
                  'sigma': numpy.array([ 1.3414586]),
                  'sum': numpy.array([ 452628.18628213]),
                  'sumsq': numpy.array([ 1376440.6075593]),
                  'trc': numpy.array([74, 74,  0, 99], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42119e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)
        
    def test203(self):
        """Test 203: Selected frequency image with other frequency unit"""
        nchan=100
        loc_unit='MHz'
        start = "%f%s" % (1420.2, loc_unit)
        width = "%f%s" % (0.01, loc_unit)
        self.task_param.update(dict(nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.4202e+09Hz',
                  'max': numpy.array([ 21.55560875]),
                  'maxpos': numpy.array([59, 21,  0, 67], dtype=numpy.int32),
                  'maxposf': '17:10:00.642, +58.42.19.808, I, 1.42087e+09Hz',
                  'mean': numpy.array([ 0.80467233]),
                  'min': numpy.array([-0.27736959]),
                  'minpos': numpy.array([58, 71,  0, 10], dtype=numpy.int32),
                  'minposf': '17:09:45.684, +61.12.21.875, I, 1.4203e+09Hz',
                  'npts': numpy.array([ 562500.]),
                  'rms': numpy.array([ 1.56429076]),
                  'sigma': numpy.array([ 1.3414586]),
                  'sum': numpy.array([ 452628.18628213]),
                  'sumsq': numpy.array([ 1376440.6075593]),
                  'trc': numpy.array([74, 74,  0, 99], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.42119e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)
        
###
# Test velocity imaging
###
class sdimaging_test3(sdimaging_unittest_base):
    """
    Test velocity imaging

       - integrated image
       - selected velocity image
       
    """
    # Input and output names
    prefix=sdimaging_unittest_base.taskname+'Test3'
    outfile=prefix+sdimaging_unittest_base.postfix
    unit='km/s'
    mode = "velocity"

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        # Common task parameters of the class
        self.task_param = dict(infiles=self.rawfile,mode=self.mode,
                               outfile=self.outfile,intent='',
                               cell=self.cell,imsize=self.imsize,
                               phasecenter=self.phasecenter,
                               gridfunction=self.gridfunction,
                               minweight=self.minweight0)

        default(sdimaging)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)

    def test300(self):
        """Test 300: Integrated image"""
        spwid = '0'
        nchan = 1
        restfreq = '1420405800.0Hz'
        ms.open(self.rawfile)
        spwinfo =  ms.getspectralwindowinfo()
        ms.close()
        chan0_freq = spwinfo[spwid]['Chan1Freq']
        bandwidth = spwinfo[spwid]['TotalWidth']
        chanwidth = spwinfo[spwid]['ChanWidth']
        cent_freq = me.frequency(spwinfo[spwid]['Frame'],
                                 qa.quantity(chan0_freq+0.5*(bandwidth-chanwidth),'Hz'))
        cent_vel = me.todoppler('radio', cent_freq, restfreq)
        # band-edge frequencies
        start_freq = me.frequency(spwinfo[spwid]['Frame'],
                                  qa.quantity(chan0_freq-0.5*chanwidth,'Hz'))
        start_vel = me.todoppler('radio', start_freq, restfreq)
        end_freq = me.frequency(spwinfo[spwid]['Frame'],
                                qa.add(start_freq['m0'],
                                       qa.quantity(bandwidth,'Hz')))
        end_vel = me.todoppler('radio', end_freq, restfreq)
        start = qa.tos(cent_vel['m0'])
        width = qa.tos(qa.sub(start_vel['m0'],end_vel['m0']))
        
        self.task_param.update(dict(restfreq=restfreq,spw=spwid,
                                    nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,1)
        self.run_test_common(self.task_param, self.statsinteg, outshape,
                             compstats=self.keys, ignoremask=True)
        
    def test301(self):
        """Test 301: Selected velocity image"""
        nchan=100
        start = "%f%s" % (-200.0, self.unit)
        width = "%f%s" % (2.0, self.unit)
        self.task_param.update(dict(nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.421353e+09Hz',
                  'max': numpy.array([ 21.97223091]),
                  'maxpos': numpy.array([ 4,  5,  0, 50], dtype=numpy.int32),
                  'maxposf': '17:30:54.243, +57.53.03.440, I, 1.42088e+09Hz',
                  'mean': numpy.array([ 0.84673187]),
                  'min': numpy.array([-0.27300295]),
                  'minpos': numpy.array([61, 71,  0, 16], dtype=numpy.int32),
                  'minposf': '17:08:30.980, +61.12.02.893, I, 1.421202e+09Hz',
                  'npts': numpy.array([ 562500.]),
                  'rms': numpy.array([ 1.6305207]),
                  'sigma': numpy.array([ 1.3934297]),
                  'sum': numpy.array([ 476286.67594505]),
                  'sumsq': numpy.array([ 1495461.22406453]),
                  'trc': numpy.array([74, 74,  0, 99], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.420415e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)

    def test302(self):
        """Test 302: Selected velocity image (different rest frequency)"""
        nchan = 100
        start = "%f%s" % (-100.0, self.unit)
        width = "%f%s" % (2.0, self.unit)
        self.task_param.update(dict(restfreq='1.420GHz',nchan=nchan,start=start,width=width))
        outshape = (self.imsize[0],self.imsize[1],1,nchan)
        refstats={'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
                  'blcf': '17:32:18.690, +57.37.28.536, I, 1.420474e+09Hz',
                  'max': numpy.array([ 1.61916351]),
                  'maxpos': numpy.array([ 4, 52,  0, 33], dtype=numpy.int32),
                  'maxposf': '17:31:47.043, +60.13.54.473, I, 1.420161e+09Hz',
                  'mean': numpy.array([ 0.12395606]),
                  'min': numpy.array([-0.41655564]),
                  'minpos': numpy.array([60, 71,  0, 93], dtype=numpy.int32),
                  'minposf': '17:08:55.879, +61.12.09.501, I, 1.419593e+09Hz',
                  'npts': numpy.array([ 562500.]),
                  'rms': numpy.array([ 0.19268371]),
                  'sigma': numpy.array([ 0.14751931]),
                  'sum': numpy.array([ 69725.28195545]),
                  'sumsq': numpy.array([ 20883.94443161]),
                  'trc': numpy.array([74, 74,  0, 99], dtype=numpy.int32),
                  'trcf': '17:03:03.151, +61.19.10.757, I, 1.419536e+09Hz'}
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=True)

###
# Test auto-resolution of spatial gridding parameters
###
class sdimaging_autocoord(sdimaging_unittest_base):
    """
    Test auto-resolution of spatial gridding parameters

       - manual setting
       - all
       - phasecenter
       - cell (get rest freq from table)
       - imsize
    """
    prefix=sdimaging_unittest_base.taskname+'Test4'
    outfile=prefix+sdimaging_unittest_base.postfix
    nchan=1
    start=0
    width=1024
    # auto calculation result of imsize
    cell_auto = "162.545308arcsec"
    imsize_auto = [73, 68]
    phasecenter_auto = "J2000 17:17:59.03 59.30.04.104"
    # manual setup
    imsize = [40, 35]
    cell = ["320arcsec", "350arcsec"]
    phasecenter = "J2000 17:18:05 59.30.05"

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile)
        # Common task parameters of the class
        self.task_param = dict(infiles=self.rawfile,outfile=self.outfile,
                               intent="",nchan=self.nchan,start=self.start,
                               width=self.width,minweight=self.minweight0)

        default(sdimaging)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, task_param, shape, dirax):
        """
        Run sdimaging and test results.
        A list of tests:
        (1) task completes without an error
        (2) ouput image and weight image exist
        (3) image shape
        (4) image direction axis
        """
        res=sdimaging(**task_param)
        #outfile = task_param['outfile']
        outprefix = task_param['outfile']
        outfile = outprefix + image_suffix
        # Tests
        self.assertEqual(res,None,
                         msg='Any error occurred during imaging')
        self._checkfile(outfile)
        self._checkfile(outprefix+".weight")
        self._checkshape(outfile,shape[0],shape[1],shape[2],shape[3])
        self._checkdirax(outfile,dirax[0], dirax[1], dirax[2])

    def test401(self):
        """test 401: Set phasecenter, cell, and imsize manually"""
        self.task_param.update(dict(cell=self.cell,imsize=self.imsize,
                                    phasecenter=self.phasecenter))
        outshape = (self.imsize[0], self.imsize[1], 1, self.nchan)
        dirax = (self.phasecenter, self.cell, self.imsize)
        self.run_test(self.task_param, outshape, dirax)

    def test402(self):
        """test 402: Automatic resolution of phasecenter, cell, and imsize"""
        self.task_param.update(dict(cell="",imsize=[],phasecenter=""))
        outshape = (self.imsize_auto[0],self.imsize_auto[1],1,self.nchan)
        dirax = (self.phasecenter_auto,self.cell_auto,self.imsize_auto)
        self.run_test(self.task_param, outshape, dirax)

    def test403(self):
        """test 403: Resolve phasecenter"""
        self.task_param.update(dict(cell=self.cell,imsize=self.imsize,phasecenter=""))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        dirax = (self.phasecenter_auto,self.cell, self.imsize)
        self.run_test(self.task_param, outshape, dirax)

    def test404(self):
        """test 404: Resolve cell"""
        self.task_param.update(dict(cell="",imsize=self.imsize,
                                    phasecenter=self.phasecenter))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        dirax = (self.phasecenter,self.cell_auto,self.imsize)
        self.run_test(self.task_param, outshape, dirax)

    def test405(self):
        """test 405: Resolve imsize"""
        ref_imsize = [38, 32]
        self.task_param.update(dict(cell=self.cell,imsize=[],phasecenter=self.phasecenter))
        outshape = (ref_imsize[0],ref_imsize[1],1,self.nchan)
        dirax = (self.phasecenter,self.cell,ref_imsize)
        self.run_test(self.task_param, outshape, dirax)

###
# Test data selection
###
class sdimaging_test_selection(selection_syntax.SelectionSyntaxTest,sdimaging_unittest_base):
    """
    Test selection syntax. Selection parameters to test are:
    field, spw (with selection), scan, stokes, and antenna
    """

    prefix = sdimaging_unittest_base.taskname+'TestSel'
    outfile = prefix+sdimaging_unittest_base.postfix
    # input MS names
    miscsel_ms = "selection_misc.ms"
    spwsel_ms = "selection_spw.ms"
    unifreq_ms = "selection_spw_unifreq.ms"
    intentsel_ms = "selection_intent.ms"
    rawfiles = [miscsel_ms, spwsel_ms, unifreq_ms, intentsel_ms]
    # default task parameters
    mode_def = "channel"
    kernel = "BOX"
    # 
    # auto calculation result of imsize
    cell_auto = "6.7275953729549656arcsec"
    imsize_auto = [21, 21]
    phasecenter_auto = "J2000 00:00:00.0 00.00.00.00"
    blc_auto = [0, 0, 0, 0]
    trc_auto = [20, 20, 0, 0]
    blcf_auto = '00:00:04.485, -00.01.07.276, I, 3e+11Hz'
    trcf_auto = '23:59:55.515, +00.01.07.276, I, 3e+11Hz'
    # Reference Statistics
    # blcf and trcf => qa.formxxx(+-qa.mul(cell_auto, 10.), "hms"/"dms", prec=3)
    # --- for "selection_misc.ms"
    unif_flux = 25.
    stat_common = {'blc': blc_auto,'trc': trc_auto,
                   'blcf': blcf_auto, 'trcf': trcf_auto}
    region_all = {'blc': blc_auto, 'trc': trc_auto}
    region_bottom = {'blc': [0, 0, 0, 0], 'trc': [20, 11, 0, 0]}
    region_top = {'blc': [0, 9, 0, 0], 'trc': [20, 20, 0, 0]}
    region_left = {'blc': [0, 0, 0, 0], 'trc': [11, 20, 0, 0]}
    region_right = {'blc': [9, 0, 0, 0], 'trc': [20, 20, 0, 0]}
    region_topleft = {'blc': [0, 9, 0, 0], 'trc': [11, 20, 0, 0]}
    region_topright = {'blc': [9, 9, 0, 0], 'trc': [20, 20, 0, 0]}
    region_bottomleft = {'blc': [0, 0, 0, 0], 'trc': [11, 11, 0, 0]}
    region_bottomright = {'blc': [9, 0, 0, 0], 'trc': [20, 11, 0, 0]}
    # --- for "selection_spw_unifreq.ms" and "selection_spw.ms"
    # flux taken from ms.statistics((column='CORRECTED_DATA', complex_value='amp', spw=idx)['mean']
    spw_flux_unifreq = [3.0008814930915833, 5.0014331340789795, 6.001709461212158]
    spw_flux = [5.001473307609558, 5.982952607795596, 3.011193051868015]  #NOTE spw=1 and 2 has relatively large O(10^-4) dispersion in intensity.
    spw_imsize_auto = [12, 12]
    spw_nchan = 10
    spw_blc_auto = [0, 0, 0, 0]
    spw_trc_auto = [11, 11, 0, 9]
    # blcf and trcf => qa.formxxx(+-qa.mul(cell_auto, 6.), "hms"/"dms", prec=3)
    spw_stat_common = {'blc': spw_blc_auto,'trc': spw_trc_auto}
    spw_region_all = {'blc': [1,1,0,0], 'trc': [11,11,0,9]}
    # select channels 2 - 7
    spw_region_chan1 = {'blc': [1,1,0,2], 'trc': [11,11,0,7]}

    @property
    def task(self):
        return sdimaging
    
    @property
    def spw_channel_selection(self):
        return True

    def setUp(self):
        for name in self.rawfiles:
            if os.path.exists(name):
                shutil.rmtree(name)
            shutil.copytree(self.datapath+name, name)
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile)
        # Common task parameters of the class
        self.task_param = dict(mode=self.mode_def,intent="",
                               gridfunction=self.kernel,outfile=self.outfile,
                               phasecenter=self.phasecenter_auto,
                               cell=self.cell_auto,imsize=self.imsize_auto)

        default(sdimaging)
        os.system( 'rm -rf '+self.prefix+'*' )

    def tearDown(self):
        for name in self.rawfiles:
            if (os.path.exists(name)):
                shutil.rmtree(name)
        os.system( 'rm -rf '+self.prefix+'*' )
        
        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, task_param, refstats, shape,
                 atol=1.e-8, rtol=1.e-5, box=None):
        self.res=self.run_task(**task_param)
        # Tests
        imsize = [shape[0], shape[1]]
        outfile = self.outfile + image_suffix
        self._checkshape(outfile,shape[0], shape[1],shape[2],shape[3])
        self._checkdirax(outfile,self.phasecenter_auto,self.cell_auto,imsize)
        self._checkstats(outfile,refstats,atol=atol,rtol=rtol)
        if box is not None:
            self._checkstats_box(outfile,refstats,box=box,
                                 atol=atol,rtol=rtol)
        

    ####################
    # Additional tests
    ####################
    #N/A Stokes & antenna selection

    ####################
    # scan
    ####################
    def test_scan_id_default(self):
        """test scan selection (scan='')"""
        scan = ''
        region =  self.region_all
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        outfile = self.outfile + image_suffix
        self._checkstats(outfile,refstats,atol=1.e-5)

    def test_scan_id_exact(self):
        """test scan selection (scan='16')"""
        scan = '16'
        region =  self.region_topright
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_scan_id_lt(self):
        """test scan selection (scan='<16')"""
        scan = '<16'
        region =  self.region_left
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_scan_id_gt(self):
        """test scan selection (scan='>16')"""
        scan = '>16'
        region =  self.region_bottomright
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_scan_id_range(self):
        """test scan selection (scan='16~17')"""
        scan = '16~17'
        region =  self.region_right
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_scan_id_list(self):
        """test scan selection (scan='16,17')"""
        scan = '16,17'
        region =  self.region_right
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_scan_id_exprlist(self):
        """test scan selection (scan='16,>16')"""
        scan = '16,>16'
        region =  self.region_right
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,scan=scan))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    ####################
    # intent
    ####################        
    def test_intent_value_default(self):
        """test intent selection (intent='')"""
        intent = ''
        region =  self.region_all
        infile = self.intentsel_ms
        self.task_param.update(dict(infiles=infile,intent=intent))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        outfile = self.outfile + image_suffix
        self._checkstats(outfile,refstats,atol=1.e-5)
        
    def test_intent_value_exact(self):
        """test intent selection (intent='OBSERVE_TARGET.ON_SOURCE')"""
        intent = 'OBSERVE_TARGET.ON_SOURCE'
        region =  self.region_bottomright
        infile = self.intentsel_ms
        self.task_param.update(dict(infiles=infile,intent=intent))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        outfile = self.outfile + image_suffix
        self._checkstats(outfile,refstats,atol=1.e-5)

    def test_intent_value_pattern(self):
        """test intent selection (intent='*CALIBRATE_PHASE*')"""
        intent = '*CALIBRATE_PHASE*'
        region =  self.region_bottomleft
        infile = self.intentsel_ms
        self.task_param.update(dict(infiles=infile,intent=intent))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        outfile = self.outfile + image_suffix
        self._checkstats(outfile,refstats,atol=1.e-5)

    ####################
    # field
    ####################
    def test_field_value_default(self):
        """test field selection (field='')"""
        field = ''
        region =  self.region_all
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_field_id_exact(self):
        """test field selection (field='6')"""
        field = '6'
        region =  self.region_bottomleft
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_id_lt(self):
        """test field selection (field='<7')"""
        field = '<7'
        region =  self.region_bottom
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_id_gt(self):
        """test field selection (field='>6')"""
        field = '>6'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_id_range(self):
        """test field selection (field='7~8')"""
        field = '7~8'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_id_list(self):
        """test field selection (field='5,7')"""
        field = '5,7'
        region =  self.region_right
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_id_exprlist(self):
        """test field selection (field='7,>7')"""
        field = '7,>7'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_value_exact(self):
        """test field selection (field='bottom')"""
        field = 'bottom'
        region =  self.region_bottom
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_value_pattern(self):
        """test field selection (field='top*')"""
        field = 'top*'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_value_list(self):
        """test field selection (field='topright,topleft')"""
        field = 'topright,topleft'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_field_mix_exprlist(self):
        """test field selection (field='topr*,>7')"""
        field = 'topr*,>7'
        region =  self.region_top
        infile = self.miscsel_ms
        self.task_param.update(dict(infiles=infile,field=field))
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, region['blc'], region['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    ####################
    # spw 
    ####################
    def test_spw_id_default(self):
        """test spw selection (spw='')"""
        spw = ''
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        region =  self.spw_region_all
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_exact(self):
        """test spw selection (spw='1')"""
        spw = '1'
        selspw = [1]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_lt(self):
        """test spw selection (spw='<2')"""
        spw = '<2'
        selspw = [0,1]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_gt(self):
        """test spw selection (spw='>0')"""
        spw = '>0'
        selspw = [1,2]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_range(self):
        """test spw selection (spw='1~2')"""
        spw = '1~2'
        selspw = [1,2]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_list(self):
        """test spw selection (spw='0,2')"""
        spw = '0,2'
        selspw = [0,2]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_exprlist(self):
        """test spw selection (spw='0,>1')"""
        spw = '0,>1'
        selspw = [0,2]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    def test_spw_id_pattern(self):
        """test spw selection (spw='*')"""
        spw = '*'
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    @unittest.expectedFailure
    def test_spw_value_frequency(self):
        """test spw selection (spw='299.4~299.6GHz')"""
        spw = '299.4~299.6GHz'
        selspw = [0]
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        region =  self.spw_region_all
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    @unittest.expectedFailure
    def test_spw_value_velocity(self):
        """test spw selection (spw='-550~-450km/s') NOT SUPPORTED YET"""
        self._default_test()

    @unittest.expectedFailure
    def test_spw_mix_exprlist(self):
        """test spw selection (spw='299.99~300.01GHz,0')"""
        spw = '299.99~300.01GHz,0'
        selspw = [0,1]
        region =  self.spw_region_all
        infile = self.unifreq_ms
        flux_list = self.spw_flux_unifreq
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,atol=1.e-5)

    #########################
    # spw with channel range
    #########################
    def test_spw_id_default_channel(self):
        """test spw selection w/ channel selection (spw=':2~7')"""
        spw = ':2~7'   #chan=2-7 in all spws should be selected
        region =  self.spw_region_chan1
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)
    
    def test_spw_id_default_frequency(self):
        """test spw selection w/ channel selection (spw=':300.4749~300.5251GHz')"""
#         spw = ':300.4749~300.5251GHz'   #chan=2-7 in spw=1 should be selected
#         selspw = [1]
        region =  self.spw_region_chan1
#         infile = self.spwsel_ms
#         flux_list = self.__get_flux_value(infile)
        ##### TEMPORARY CHANGING INPUT DATA due to seg fault in sdimaging caused by a bug in ms.msseltoindex() #####
        infile = self.unifreq_ms
        spw = '*:299.9749~300.0251GHz'   #chan=2-7 of spw=1 should be selected
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        # end of temporal change
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    @unittest.expectedFailure
    def test_spw_id_default_velocity(self):
        """test spw selection w/ channel selection (spw='X~Ykm/s') NOT SUPPORTED YET"""
        self._default_test()
        
    def test_spw_id_default_list(self):
        """test spw selection w/ channel selection (spw=':6~7;2~5')"""
        spw = ':6~7;2~5'   #chan=2-7 in all spws should be selected
        region =  self.spw_region_chan1
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)
        
    def test_spw_id_exact_channel(self):
        """test spw selection w/ channel selection (spw='2:2~7')"""
        spw = '2:2~7'   #chan=2-7 of spw=2 should be selected
        selspw = [2]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    def test_spw_id_exact_frequency(self):
        """test spw selection w/ channel selection (spw='1:300.4749~300.5251GHz')"""
        spw = '1:300.4749~300.5251GHz'   #chan=2-7 of spw=1 should be selected
        selspw = [1]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    @unittest.expectedFailure
    def test_spw_id_exact_velocity(self):
        """test spw selection w/ channel selection (spw='0:X~Ykm/s') NOT SUPPORTED YET"""
        self._default_test()
        
    def test_spw_id_exact_list(self):
        """test spw selection w/ channel selection (spw='2:6~7;2~5')"""
        spw = '2:6~7;2~5'   #chan=2-7 of spw=2 should be selected
        selspw = [2]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    def test_spw_id_pattern_channel(self):
        """test spw selection w/ channel selection (spw='*:2~7')"""
        spw = '*:2~7'
        region =  self.spw_region_chan1
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)

    def test_spw_id_pattern_frequency(self):
        """test spw selection w/ channel selection (spw='*:300.4749~300.5251GHz')"""
        #spw = '*:300.4749~300.5251GHz'   #chan=2-7 of spw=1 should be selected
        #selspw = [1]
        region =  self.spw_region_chan1
        #infile = self.spwsel_ms
        #flux_list = self.__get_flux_value(infile)
        ##### TEMPORARY CHANGING INPUT DATA due to seg fault in sdimaging caused by a bug in ms.msseltoindex() #####
        infile = self.unifreq_ms
        spw = '*:299.9749~300.0251GHz'   #chan=2-7 of spw=1 should be selected
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        # end of temporal change
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
 
    @unittest.expectedFailure
    def test_spw_id_pattern_velocity(self):
        """test spw selection w/ channel selection (spw='*:X~Ykm/s') NOT SUPPORTED YET"""
        self._default_test()
        
    def test_spw_id_pattern_list(self):
        """test spw selection w/ channel selection (spw='*:6~7;2~5')"""
        spw = '*:6~7;2~5'
        region =  self.spw_region_chan1
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        selspw = list(range(len(flux_list)))
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)
        
    def test_spw_value_frequency_channel(self):
        """test spw selection w/ channel selection (spw='300.4~300.5GHz:2~7')"""
        spw = '300.4~300.5GHz:2~7'
        selspw = [1]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    def test_spw_value_frequency_frequency(self):
        """test spw selection w/ channel selection (spw='300.4~300.5GHz:300.4749~300.5251GHz')"""
        spw = '300.4~300.5GHz:300.4749~300.5251GHz'   #chan=2-7 of spw=1 should be selected'
        selspw = [1]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-3,rtol=1.e-3)
        
    @unittest.expectedFailure
    def test_spw_value_frequency_velocity(self):
        """test spw selection w/ channel selection (spw='A~BHz:X~Ykm/s') NOT SUPPORTED YET"""
        self._default_test()
        
    @unittest.expectedFailure
    def test_spw_value_frequency_list(self):
        """test spw selection w/ channel selection (spw='299.9~300.1GHz:6~7;2~5')"""
        spw = '299.9~300.1GHz:6~7;2~5'
        selspw = [0]
        region =  self.spw_region_chan1
        infile = self.spwsel_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)
        
    @unittest.expectedFailure
    def test_spw_value_velocity_channel(self):
        """test spw selection w/ channel selection (spw='X~Ykm/s:A~B') NOT SUPPORTED YET"""
        self._default_test()
        
    @unittest.expectedFailure
    def test_spw_value_velocity_frequency(self):
        """test spw selection w/ channel selection (spw='X~Ykm/s:A~BHz') NOT SUPPORTED YET"""
        self._default_test()
        
    @unittest.expectedFailure
    def test_spw_value_velocity_velocity(self):
        """test spw selection w/ channel selection (spw='X~Ykm/s:Z~Wkm/s') NOT SUPPORTED YET"""
        self._default_test()
        
    @unittest.expectedFailure
    def test_spw_value_velocity_list(self):
        """test spw selection w/ channel selection (spw='X~Ykm/s:A~B;C~D') NOT SUPPORTED YET"""
        self._default_test()
        
    def test_spw_id_list_channel(self):
        """test spw selection w/ channel selection (spw='1:2~7,2:2~7')"""
        spw = '0:2~7,2:2~7'
        selspw = [0, 2]
        region =  self.spw_region_chan1
        infile = self.unifreq_ms
        flux_list = self.__get_flux_value(infile)
        self.task_param.update(dict(infiles=infile,spw=spw,imsize=self.spw_imsize_auto))
        flux = sum([flux_list[idx] for idx in selspw])/float(len(selspw))
        refstats = merge_dict(self.spw_stat_common, construct_refstat_uniform(flux, region['blc'], region['trc']) )
        out_shape = (self.spw_imsize_auto[0],self.spw_imsize_auto[1],1,self.spw_nchan)
        # Tests
        self.run_test(self.task_param,refstats,out_shape,box=region,atol=1.e-5)
        
    ####################
    # Helper functions
    ####################
    def _checkstats_box(self,name, ref, compstats=None, atol=1.e-8, rtol=1.e-5, box=None, ignoremask=False):
        """
        A test function to compare statistics of a box region of an image
        with reference values.
        Arguments:
            name  :  name of an image to test statistics
            ref   :  a record (dictionary) of the reference statistic values
            compstats : a list of names of statistis to compare. By default,
                        the list is taken from all keys in ref
            atol  : absolute tolerance (see help in numpy.allclose)
            rtol  : relative tolerance (see help in numpy.allclose)
            box   : a dictionary that specifies a box region of image to
                    calculate statistics. it should be a dictionary with keys,
                    'blc' and 'trc' in pixel unit.
            ignoremask : when True, mask in image is ignored and statistics
                         are calculated from whole pixels in image. default
                         is False (take image mask into account).
        """
        boxreg = _rg.box(**box) if box is not None else None
        refstats = ref.copy()
        refstats.update(box)
        for stats in ['blcf', 'trcf']:
            if stats in refstats: refstats.pop(stats)
        self._checkstats(name,refstats,region=boxreg,
                         compstats=compstats,atol=atol,rtol=rtol,
                         ignoremask=ignoremask)

    def __get_flux_value(self, infile):
        """ returns proper flux list """
        if infile == self.miscsel_ms:
            return self.unif_flux
        elif infile == self.spwsel_ms:
            return self.spw_flux
        elif infile == self.unifreq_ms:
            return self.spw_flux_unifreq
        else: raise Exception("Internal error: invalid input file to get flux value.")

###
# Test to verify if flag information is handled properly
###
class sdimaging_test_flag(sdimaging_unittest_base):
    """
    Test to verify if flag information is handled properly
       
       - If a channel is flagged, the data of the channel must not be
         added to the output CASA image.
       - If all channels of a spectrum are flagged (i.e., row-flagged
         in original Scantable), the whole data of the spectrum must
         not be added to the output CASA image.
       - Flagged channels must not be modified by sdimaging.

       The input data: sdimaging_flagtest.ms
       - this data contains 768 spectra covering 32x24 grid-like points 
         with interval of 10^-5 radian (the position of the bottom-right
         corner is (00:00:00.0, 00.00.00.0) in (RA, Dec)). Amongst the
         spectra, 32x8 spectra corresponding to the middle 1/3 of
         the survey area, all channels are flagged, and for the half of 
         the rest spectra at the smaller side in RA (16x8x2 spectra) are
         flagged at channels 2 to 6 (5 out of 10 channels are flagged).
         The rest (16x8x2 spectra) have no flagged channels (see below).

         The row index, data value (constant along spectral channel), and
         flag status of the data spectra distribute as follows. North is
         up, and east is to the left:

         Row index
         ---------
         
          767, 766, 765, ..., 752, 751, ..., 738, 737, 736,
           ...............................................
          543, 542, 541, ..., 528, 527, ..., 514, 513, 512,
          511, 510, 509, ..., 496, 495, ..., 482, 481, 480,
           ...............................................
          287, 286, 285, ..., 272, 271, ..., 258, 257, 256,
          255, 254, 253, ..., 240, 239, ..., 226, 225, 224,
           ...............................................
           31,  30,  29, ...,  16,  15, ...,   2,   1,   0
         
         Data Value
         ----------
         (A:6.00171, B, 5.00143, C:3.00088)

              -----                     ------
                   A, A, A, ..., A, A, A
                   .....................   8 points
                   A, A, A, ..., A, A, A
              -----                     ------
                   B, B, B, ..., B, B, B
                   .....................   8 points
                   B, B, B, ..., B, B, B
              -----                     ------
                   C, C, C, ..., C, C, C
                   .....................   8 points
                   C, C, C, ..., C, C, C
              -----                     ------

         Flag Status
         -----------
         (A:no flag, B:partly flagged (2-6), C:fully flagged)

     -----|                |               |-----
           A, A, A, ..., A, B, ..., B, B, B
           ................................      8 points
           A, A, A, ..., A, B, ..., B, B, B
     -----                                 ------
           C, C, C, ..., C, C, ..., C, C, C
           ................................      8 points
           C, C, C, ..., C, C, ..., C, C, C
     -----                                 ------
           A, A, A, ..., A, B, ..., B, B, B
           ................................      8 points
           A, A, A, ..., A, B, ..., B, B, B
     -----|                |               |-----
          |   16 points    |   16 points   |
          

    """
    rawfile='sdimaging_flagtest.ms'
    prefix=sdimaging_unittest_base.taskname+'TestFlag'
    outfile=prefix+sdimaging_unittest_base.postfix
    maskfile = outfile + image_suffix + '/mask0'
    weightfile = outfile + '.weight'
    
    gridfunction = "BOX"
    imsize = [32, 24]
    cellarcsec = 2.062648 #= 0.00001*180.0/3.1415926535897932384*3600.0
    cell = [str(cellarcsec)+'arcsec', str(cellarcsec)+'arcsec']
    pcra = cellarcsec*15.0/15.0
    pcdec = cellarcsec*11.8
    phasecenter = "J2000 00:00:0"+str(pcra)+" 00.00."+str(pcdec)

    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile)
        default(sdimaging)
        with tbmanager(self.rawfile) as tb:
            self.nchan = len(tb.getcell('DATA', 0)[0])
            
        # fix timestamp issue
        self.fix_timestamp()

    def tearDown(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)
        
    def fix_timestamp(self):
        # fix duplicated timestamp issue
        # data taken by three spws have essentially same timestamp 
        # but they are intended to be allocated to different pointing 
        # directions. to enable it, timestamps are artificially shifted 
        # by a value significantly larger than integration time.
        with tbmanager(self.rawfile, nomodify=False) as tb:
            nrow = tb.nrows()
            ddid = numpy.unique(tb.getcol('DATA_DESC_ID'))
            nchunk = len(ddid)
            nrow_chunk = nrow / nchunk
            torig = tb.getcol('TIME')
            interval = tb.getcol('INTERVAL')
            max_interval = interval.max()
            tshift = numpy.empty_like(torig)
            for ichunk in range(nchunk):
                ifrom = ichunk * nrow_chunk
                ito = (ichunk+1) * nrow_chunk
                tshift[ifrom:ito] = torig[ifrom:ito] + ichunk * max_interval
            tb.putcol('TIME', tshift)
            
        with tbmanager(os.path.join(self.rawfile, 'POINTING'), nomodify=False) as tb:
            tb.putcol('TIME', tshift) 

    def testFlag01(self):
        """testFlag01: """
        res=sdimaging(infiles=self.rawfile,outfile=self.outfile,intent="",gridfunction=self.gridfunction,cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter,minweight=self.minweight0)
        self.assertEqual(res,None,
                         msg='Any error occurred during imaging')
        outfile = self.outfile + image_suffix
        self._checkshape(outfile,self.imsize[0],self.imsize[1],1,self.nchan)
        self._set_data_ranges()
        self._check_data()
        self._check_mask()
        self._check_weight()

    def testFlag02(self):
        res=sdimaging(infiles=self.rawfile,outfile=self.outfile,intent="",width=10,gridfunction=self.gridfunction,cell=self.cell,imsize=self.imsize,phasecenter=self.phasecenter,minweight=self.minweight0)
        self.assertEqual(res,None,
                         msg='Any error occurred during imaging')
        outfile = self.outfile + image_suffix
        self._checkshape(outfile,self.imsize[0],self.imsize[1],1,1)
        self._set_data_ranges(True)
        self._check_data(True)
        self._check_mask(True)
        self._check_weight(True)
    
    def _set_data_ranges(self, chanmerge=False):
        xn = 2
        xw = self.imsize[0]/xn
        self.x_range = []
        for i in range(xn):
            self.x_range.append([xw*i, xw*(i+1)])
        yn = 3
        yw = self.imsize[1]/yn
        self.y_range = []
        for i in range(yn):
            self.y_range.append([yw*i, yw*(i+1)])
        self.f_range = [[0,1]] if chanmerge else [[0,2],[2,7],[7,10]]

    def _check_data(self, chanmerge=False):
        val = self._get_refvalues(self.rawfile, chanmerge)
        idx = 0
        outfile = self.outfile + image_suffix
        for i in range(len(self.x_range)):
            for j in range(len(self.y_range)):
                for k in range(len(self.f_range)):
                    self._checkvalue(outfile, False, self.x_range[i], self.y_range[j], self.f_range[k],  val[idx], chanmerge)
                    idx += 1

    def _check_mask(self, chanmerge=False):
        val = self._get_refmask(self.maskfile, chanmerge)
        idx = 0
        for i in range(len(self.x_range)):
            for j in range(len(self.y_range)):
                for k in range(len(self.f_range)):
                    self._checkvalue(self.maskfile, True, self.x_range[i], self.y_range[j], self.f_range[k],  val[idx], chanmerge)
                    idx += 1

    def _check_weight(self, chanmerge=False):
        val = self._get_refweight(self.weightfile, chanmerge)
        idx = 0
        for i in range(len(self.x_range)):
            for j in range(len(self.y_range)):
                for k in range(len(self.f_range)):
                    self._checkvalue(self.weightfile, False, self.x_range[i], self.y_range[j], self.f_range[k],  val[idx], chanmerge)
                    idx += 1

    def _get_refmask(self, file, chanmerge=False):
        res = []
        with tbmanager(file) as tb:
            for i in [0, self.imsize[0]/2]:
                for j in [0, self.imsize[1]/3, self.imsize[1]*2/3]:
                    k_range = [0] if chanmerge else [0, 5, 9]
                    for k in k_range:
                        res.append(tb.getcell('PagedArray', 0)[i][j][0][k].real)
        return res
            
    def _get_refweight(self, file, chanmerge=False):
        res = []
        with tbmanager(file) as tb:
            for i in [0, self.imsize[0]/2]:
                for j in [0, self.imsize[1]/3, self.imsize[1]*2/3]:
                    k_range = [0] if chanmerge else [0, 5, 9]
                    for k in k_range:
                        res.append(tb.getcell('map', 0)[i][j][0][k].real)
        return res
            
    def _get_refvalues(self, file, chanmerge=False):
        res = []
        with tbmanager(file) as tb:
            for i in [self.imsize[0]/2, 0]:
                for j in [0, self.imsize[1]/3, self.imsize[1]*2/3]:
                    irow = self.imsize[0]*j+i
                    if chanmerge:
                        if (tb.getcell('FLAG', irow)[0]==True).all():
                            res.append(0.0)
                        else:
                            res.append(tb.getcell('DATA', irow)[0][0].real)
                    else:
                        if (tb.getcell('FLAG', irow)[0]==True).all():
                            for k in range(3): res.append(0.0)
                        else:
                            res.append(tb.getcell('DATA', irow)[0][0].real)
                            if (tb.getcell('FLAG', irow)[0][5]):
                                res.append(0.0)
                            else:
                                res.append(tb.getcell('DATA', irow)[0][5].real)
                            res.append(tb.getcell('DATA', irow)[0][9].real)
        return res

    def _checkvalue(self, file, is_maskfile, x_range, y_range, f_range, ref_value, chanmerge=False):
        tol=1e-5
        colname = 'PagedArray' if is_maskfile else 'map'
        with tbmanager(file) as tb:
            val = tb.getcell(colname, 0)

        for i in range(x_range[0], x_range[1]):
            for j in range(y_range[0], y_range[1]):
                for k in range(f_range[0], f_range[1]):
                    diff_value = abs(val[i][j][0][k]-ref_value)
                    self.assertTrue(diff_value < tol)


class sdimaging_test_polflag(sdimaging_unittest_base):
    """
    Test imaging of an MS one of polarization (XX) is completely flagged.
    """
    prefix = sdimaging_unittest_base.taskname+'TestPol'
    outfile = prefix+sdimaging_unittest_base.postfix
    # input MS names
    infiles = "selection_misc.ms"
    # default task parameters
    mode = "channel"
    kernel = "BOX"
    # 
    # auto calculation result of imsize
    cell_auto = "6.7275953729549656arcsec"
    imsize_auto = [21, 21]
    phasecenter_auto = "J2000 00:00:00.0 00.00.00.00"
    blc_auto = [0, 0, 0, 0]
    trc_auto = [20, 20, 0, 0]
    blcf_auto = '00:00:04.485, -00.01.07.276, I, 3e+11Hz'
    trcf_auto = '23:59:55.515, +00.01.07.276, I, 3e+11Hz'
    # Reference Statistics
    # blcf and trcf => qa.formxxx(+-qa.mul(cell_auto, 10.), "hms"/"dms", prec=3)
    # --- for "selection_misc.ms"
    unif_flux = 32.
    stat_common = {'blc': blc_auto,'trc': trc_auto,
                   'blcf': blcf_auto, 'trcf': trcf_auto}
    region_all = {'blc': blc_auto, 'trc': trc_auto}

    def setUp(self):
        if os.path.exists(self.infiles):
            shutil.rmtree(self.infiles)
        shutil.copytree(self.datapath+self.infiles, self.infiles)
        os.system( 'rm -rf '+self.prefix+'*' )

        # Common task parameters of the class
        self.task_param = dict(infiles=self.infiles,mode=self.mode,
                               gridfunction=self.kernel,
                               outfile=self.outfile,intent="")
        # flag ALL POL='XX'
        flagdata(vis=self.infiles,mode='manual',correlation='XX',action='apply')

        default(sdimaging)

    def tearDown(self):
        if os.path.exists(self.infiles):
            shutil.rmtree(self.infiles)
        # Since the data is flagged by flagdata, flagversions directory 
        # is automatically created. This must be removed
        flagversions = self.infiles + '.flagversions'
        if os.path.exists(flagversions):
            shutil.rmtree(flagversions)
        # By executing flagdata task, flagdata.last is created automatically
        # This must also be removed
        flagdata_last = 'flagdata.last'
        if os.path.exists(flagdata_last):
            os.remove(flagdata_last)
        # Remove test image and its weight image
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, task_param, refstats, shape,
                 atol=1.e-8, rtol=1.e-5, box=None):
        self.res=sdimaging(**task_param)
        # Tests
        imsize = [shape[0], shape[1]]
        outfile = self.outfile + image_suffix
        self._checkshape(outfile,shape[0], shape[1],shape[2],shape[3])
        self._checkdirax(outfile,self.phasecenter_auto,self.cell_auto,imsize)
        self._checkstats(outfile,refstats,atol=atol,rtol=rtol)

    def test_i(self):
        """test stokes='I': image constructed by unflagged YY pol (NB after imager migration: image weights all zero)"""
        self.task_param['stokes'] = 'I'
        # Tests
        # In true Stokes mode, all correlation components are flagged when any of them is flagged
        # Stokes value consistent with older imager based implementation is pseudo-Stokes
        #refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, self.region_all['blc'], self.region_all['trc']) )
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(0.0, self.region_all['blc'], self.region_all['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        
    def test_pseudo_i(self):
        """test pseudo stokes I: image constructed by unflagged YY pol"""
        self.task_param['stokes'] = 'pseudoI'
        # Tests
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, self.region_all['blc'], self.region_all['trc']) )
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
       

    def test_xx(self):
        """test stokes='XX' (flagged): image weights all zero"""
        stokes = 'XX'
        self.task_param['stokes'] = stokes
        # Tests
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(0.0, self.region_all['blc'], self.region_all['trc']) )
        refstats['blcf'] = refstats['blcf'].replace('I', stokes)
        refstats['trcf'] = refstats['trcf'].replace('I', stokes)
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)

    def test_yy(self):
        """test stokes='YY': image constructed by YY pol"""
        stokes = 'YY'
        self.task_param['stokes'] = stokes
        # Tests
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, self.region_all['blc'], self.region_all['trc']) )
        refstats['blcf'] = refstats['blcf'].replace('I', stokes)
        refstats['trcf'] = refstats['trcf'].replace('I', stokes)
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],1,1)
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)

    def test_xxyy(self):
        """test stokes='XXYY': """
        self.task_param['stokes'] = 'XXYY'
        # Tests
        refstats = merge_dict(self.stat_common, construct_refstat_uniform(self.unif_flux, self.region_all['blc'], self.region_all['trc']) )
        refstats['blcf'] = refstats['blcf'].replace('I', 'XX')
        refstats['trcf'] = refstats['trcf'].replace('I', 'YY')
        refstats['trc'][2] = 1 # the image is in 2 polarization
        out_shape = (self.imsize_auto[0],self.imsize_auto[1],2,1)
        self.run_test(self.task_param, refstats, out_shape,atol=1.e-5)
        # statistics of YY only
        refstats['blc'][2] = 1
        for key in ['blcf', 'trcf']: refstats.pop(key)
        box = _rg.box(blc=refstats['blc'],trc=refstats['trc'])
        outfile = self.outfile + image_suffix
        self._checkstats(outfile,refstats,atol=1.e-5,region=box)

class sdimaging_test_mslist(sdimaging_unittest_base):
    """
    Test more than one MSes as inputs
    
    """
    prefix = sdimaging_unittest_base.taskname+'TestListMS'
    outfile = prefix+sdimaging_unittest_base.postfix
    clearup = True
    # input MS names
    org_ms = "selection_misc.ms"
    # imaging parameter
    mode = "channel"
    kernel = "BOX"
    infiles = ["multi-in1", "multi-in2"]
    outfile = prefix+".im"
    # auto calculation result of imsize
    cell = "6.7275953729549656arcsec"
    imsize = [21, 21]
    phasecenter = "J2000 00:00:00.0 00.00.00.00"
    blc = [0, 0, 0, 0]
    trc = [20, 20, 0, 0]
    blcf = '00:00:04.485, -00.01.07.276, I, 3e+11Hz'
    trcf = '23:59:55.515, +00.01.07.276, I, 3e+11Hz'
    # selection
    field_list = ['8,6', '7,5']
    scan_list = ['15', '16,17']
    spw = '0'    
    # Reference Statistics
    # blcf and trcf => qa.formxxx(+-qa.mul(cell_auto, 10.), "hms"/"dms", prec=3)
    # --- for "selection_misc.ms"
    unif_flux = 25.
    refstats = merge_dict({'blc': blc,'trc': trc, 'blcf': blcf, 'trcf': trcf},
                          construct_refstat_uniform(unif_flux, blc, trc) )
    #nvalid = imsize[0]*imsize[1] #21*21
    #{'min': [unif_flux], 'max': [unif_flux], 'rms': [unif_flux],
    # 'sigma': [0.], 'mean': [unif_flux], 'npts': [nvalid],
    # 'sum': [unif_flux*nvalid], 'sumsq': [unif_flux**2*nvalid],
    # 'blc': blc,'trc': trc, 'blcf': blcf, 'trcf': trcf}
    
    def setUp(self):
        if os.path.exists(self.outfile):
            os.system('rm -rf %s*' % self.outfile)
        for name in self.infiles:
            if os.path.exists(name):
                shutil.rmtree(name)
            shutil.copytree(self.datapath+self.org_ms, name)

        default(sdimaging)
        self.default_param = dict(infiles = self.infiles,
                                  outfile = self.outfile,
                                  intent="",
                                  cell = self.cell,
                                  imsize = self.imsize,
                                  phasecenter = self.phasecenter,
                                  mode=self.mode,
                                  gridfunction=self.kernel,
                                  minweight=0.0)

    def tearDown(self):
        if self.clearup:
            outfile = self.outfile + image_suffix
            if os.path.exists(outfile):
                os.system('rm -rf %s*' % self.outfile)
            for name in self.infiles:
                if os.path.exists(name):
                    shutil.rmtree(name)
                    
        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, task_param=None,refstats=None):
        if task_param is None:
            task_param = self.default_param
        if refstats is None:
            refstats = self.refstats
        res=sdimaging(**task_param)
        outfile = self.outfile + image_suffix
        self._checkshape(outfile,self.imsize[0],self.imsize[1],1,1)
        self._checkdirax(outfile,self.phasecenter,self.cell,self.imsize)
        self._checkstats(outfile,refstats,atol=1.e-5)

    ###########################
    # Tests
    ###########################
    def multi_input(self):
        """Test two MSes as input"""
        self.run_test()

    def test_string_selection(self):
        """Test data selection by string (2 MS inputs)"""
        self.default_param['field'] = str(',').join(self.field_list)
        self.default_param['scan'] = str(',').join(self.scan_list)
        self.default_param['spw'] = self.spw
        self.run_test()
        
    def test_1elemlist_selection(self):
        """Test data selection by single element list (2 MS inputs)"""
        self.default_param['field'] = [str(',').join(self.field_list)]
        self.default_param['scan'] = [str(',').join(self.scan_list)]
        self.default_param['spw'] = [self.spw]
        self.run_test()
        
    def test_2elemlist_selection(self):
        """Test data selection by 2 elements list (2 MS inputs)"""
        self.default_param['field'] = self.field_list
        self.default_param['scan'] = self.scan_list
        self.default_param['spw'] = [self.spw, self.spw]
        self.run_test()

###
#
# Test ways to define image rest frequency
#
###
class sdimaging_test_restfreq(sdimaging_unittest_base):
    """
    Unit test for task sdimaging 
    
    Test 3 different ways to define image rest frequency.
    (1) defined by parameter restfreq (test_restfreq_param)
    (2) obtain a value in REST_FREQUENCY column in SOURCE subtable
        (test_restfreq_source)
    (3) define as mean frequency of representative SPW (test_restfreq_mean)
    
    The rest frequency value will affect three numbers in image:
    - the rest frequency of the image
    - the default cell size of the image
    - the beam size of the image
    """
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdimaging/'
    infiles = 'selection_spw.ms'
    outfile = 'sdimaging_restfreq.im'
    param_base = dict(infiles=infiles,outfile=outfile,intent="",
                      outframe='lsrk',stokes='I',spw='1',
                      phasecenter='J2000 00:00:00 00.00.00',
                      restfreq='',overwrite=True)
    unifval = 5.98155

    def setUp(self):
        if os.path.exists(self.infiles):
            shutil.rmtree(self.infiles)
        shutil.copytree(self.datapath+self.infiles, self.infiles)
        default(sdimaging)
        self.param = self.param_base.copy()
        
    def tearDown(self):
        if os.path.exists(self.infiles):
            shutil.rmtree(self.infiles)
        os.system('rm -rf {0}*'.format(self.outfile))

        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, restfreq_ref, beam_ref, cell_ref, stats, **kwargs):
        self.param.update(**kwargs)
        sdimaging(**self.param)
        stats.pop('sumsq')
        outfile = self.outfile + image_suffix
        self._checkstats(outfile, stats, atol=1.e-3, rtol=1.e-3)
        self._check_beam(outfile, beam_ref)
        # check restfreq
        self._check_restfreq(outfile,restfreq_ref)
        # check cell size
        self._checkdirax(outfile, self.param['phasecenter'],
                         cell_ref, self.param['imsize'])

    def test_restfreq_param(self):
        """Rest frequency from restfreq parameter"""
        restfreq='200GHz'
        beam_ref = dict(major='30.276442arcsec',minor='30.276442arcsec')
        cell_ref = '10.091393059432447arcsec'
        stats = construct_refstat_uniform(self.unifval,[0, 0, 0, 0],
                                          [7 , 7 ,  0,  9])
        self.run_test(restfreq, beam_ref, cell_ref, stats,
                      restfreq=restfreq,imsize=[8,8])
        
    def test_restfreq_source(self):
        """Rest Frequency from SOURCE table"""
        restfreq='300GHz'
        beam_ref = dict(major='20.339973arcsec',minor='20.339973arcsec')
        cell_ref = '6.727595372954963arcsec'
        stats = construct_refstat_uniform(self.unifval,[0, 0, 0, 0],
                                          [10, 10,  0,  9])
        self.run_test(restfreq, beam_ref, cell_ref, stats,
                      restfreq='',imsize=[11,11])

    def test_restfreq_mean(self):
        """Rest frequency from mean of SPW frequencies"""
        restfreq='300.5GHz'
        beam_ref = dict(major='20.303418arcsec', minor='20.303418arcsec')
        cell_ref = '6.716401370670513arcsec'
        stats = construct_refstat_uniform(self.unifval,[0, 0, 0, 0],
                                          [10, 10,  0,  9])
        # remove REST_REQUENCY in SOURCE TABLE
        tb.open(self.infiles+'/SOURCE', nomodify=False)
        rf = tb.getcell('REST_FREQUENCY',0)
        rf.resize(0)
        for idx in range(tb.nrows()):
            tb.putcell('REST_FREQUENCY', idx, rf)
            self.assertTrue(len(tb.getcell('REST_FREQUENCY',idx))==0)
        tb.flush()
        tb.close()
        self.run_test(restfreq, beam_ref, cell_ref, stats,
                      restfreq='', imsize=[11,11])

###
#
# Test case for automatic phasecenter calculation
#
###
class sdimaging_test_mapextent(unittest.TestCase):
    """
    Unit test for task sdimaging 
    
    This test case defines automatic calculation of phasecenter and 
    imsize. Basic tests has already been defined in test109 and test402 
    so that some specific tests are defined here:
    
        test_azel_pointing -- Verify phasecenter in J2000 is properly calculated 
                              from AZELGEO pointing direction
        test_data_selection -- Verify phasecenter is properly calculated from 
                               only selected data
        test_ephemeris -- Verify phasecenter for ephemeris source
    """
    datapath=os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdimaging/'
    infiles_ephem = ['Uranus1.cal.Ant0.spw34.ms',
                     'Uranus2.cal.Ant0.spw34.ms']
    infiles_selection = 'selection_misc.ms'
    infiles_azel = 'azelpointing.ms'
    outfile = 'sdimaging_test_mapextent.im'
    
    scan = '16'
    region_topright = {'blc': [9, 9, 0, 0], 'trc': [20, 20, 0, 0]}
    
    param_base = {'mode': 'channel',
                  'start': 0,
                  'nchan': 1,
                  'width': 1, 
                  'cell': '6.7arcsec',
                  'gridfunction': 'BOX',
                  'outfile': outfile,
                  'intent': ""}
    
    def __remove_table(self, f):
        if os.path.exists(f):
            shutil.rmtree(f)
    
    def __copy_table(self, f):
        self.__remove_table(f)
        testutils.copytree_ignore_subversion(self.datapath, f)
        
    def setUp(self):
        default(sdimaging)
        self.param = self.param_base.copy()
        
    def tearDown(self):
        for infile in self.infiles_ephem:
            self.__remove_table(infile)
        self.__remove_table(self.infiles_selection)
        self.__remove_table(self.infiles_azel)
        #self.__remove_table(self.outfile)
        os.system('rm -rf %s*'%(self.outfile))
        
        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, **kwargs):
        self.param.update(**kwargs)
        status = sdimaging(**self.param)
        self.assertIsNone(status, msg='sdimaging failed to execute')
        outfile = self.outfile + image_suffix
        self.assertTrue(os.path.exists(outfile), msg='output image is not created.')
        
    def verify_mapextent(self, npix_ref, blc_ref, trc_ref):
        outfile = self.outfile + image_suffix
        self.assertTrue(os.path.exists(outfile), msg='output image is not created.')
        stats = calc_statistics(outfile)
        map_property = calc_mapproperty(stats)
        npix = map_property['npix']
        blc = map_property['blc']
        trc = map_property['trc']
        extent = map_property['extent']
        #blc_ref = numpy.array([0.0, 0.0])
        #trc_ref = numpy.array(map(str_to_deg, ['23:59:55.515', '+00.01.07.276']))
        if trc_ref[0] > 180.0:
            trc_ref[0] -= 360.0
        if blc_ref[0] > 180.0:
            blc_ref[0] -= 360.0
        #self.verify_mapextent(npix_ref, blc_ref, trc_ref)
        # resulting map contain reference position
        print('npix', npix, 'npix_ref', npix_ref)
        print('blc', blc, 'blc_ref', blc_ref)
        print('trc', trc, 'trc_ref', trc_ref) 
        print('extent', extent)
        # check if map area covers whole pointing data
        # this is done by comparing blc and trc with their references 
        # that are usually computed from actual distribution of 
        # pointing direction (which is calculated by sdsave task) 
        self.assertTrue(all(npix == npix_ref), msg='Unexpected image pixel number')
        self.assertTrue(blc[0] >= blc_ref[0], msg='Unexpected coordinate (blc RA is too narrow)')
        self.assertTrue(blc[1] <= blc_ref[1], msg='Unexpected coordinate (blc DEC is too narrow)')
        self.assertTrue(trc[0] <= trc_ref[0], msg='Unexpected coordinate (trc RA is too narrow)')
        self.assertTrue(trc[1] >= trc_ref[1], msg='Unexpected coordinate (trc DEC is too narrow)')
        # also check if resulting map is not too wide
        # acceptable margin is 5% of the map extent
        margin = 0.05
        self.assertTrue(blc[0] < blc_ref[0] + margin * extent[0], msg='Unexpected coordinate (blc RA is too wide)')
        self.assertTrue(blc[1] > blc_ref[1] - margin * extent[1], msg='Unexpected coordinate (blc DEC is too wide)')
        self.assertTrue(trc[0] > trc_ref[0] - margin * extent[0], msg='Unexpected coordinate (trc RA is too wide)')
        self.assertTrue(trc[1] < trc_ref[1] + margin * extent[1], msg='Unexpected coordinate (trc DEC is too wide)')              
        
    def test_azel_pointing(self):
        """test_azel_pointing: Verify phasecenter in J2000 is properly calculated from AZELGEO pointing direction"""
        self.__copy_table(self.infiles_azel)
        self.run_test(infiles=self.infiles_azel)
        npix_ref = numpy.array([27,37])
        #blc_ref, trc_ref = get_mapextent(self.infiles_azel) #CAS-10301
        blc_ref = numpy.array([-85.2565977,  -13.87524395]) #CAS-10301
        trc_ref = numpy.array([-85.30504227, -13.80972133]) #CAS-10301
        self.verify_mapextent(npix_ref, blc_ref, trc_ref)
    
    def test_data_selection(self):
        self.__copy_table(self.infiles_selection)
        # here imsize is explicitly set to 13 
        # this is because that auto-calculated imsize is 12 (even number) 
        # it is known that phasecenter will not be a map center when 
        # imsize is even number so that expected map coverage is shifted 
        # by an order of 0.5 pixel 
        # this effect causes unexpected failure of the test
        self.run_test(infiles=self.infiles_selection, scan='16', imsize=13)
        npix_ref = numpy.array([13,13])
        #blc_ref, trc_ref = get_mapextent(self.infiles_selection, scan='16') #CAS-10301
        blc_ref = numpy.array([ 0.00202179, -0.00202178]) #CAS-10301
        trc_ref = numpy.array([-0.01819663,  0.01819663]) #CAS-10301
        self.verify_mapextent(npix_ref, blc_ref, trc_ref)

    def test_ephemeris(self):
        for infile in self.infiles_ephem:
            self.__copy_table(infile)
        self.run_test(infiles=self.infiles_ephem, ephemsrcname='Uranus', restfreq='230GHz')
        npix_ref = numpy.array([37,26])
        # set reference value manually since expected map area for 
        # ephemeris object is difficult to calculate 
        blcf_ref = '00:46:43.672 +04.14.51.504'
        trcf_ref = '00:46:27.547 +04.17.39.004'
        blc_ref = numpy.array([qa.quantity(x)['value'] for x in blcf_ref.split()])
        trc_ref = numpy.array([qa.quantity(x)['value'] for x in trcf_ref.split()])
        #blc_ref, trc_ref = get_mapextent_ephemeris(self.infiles_ephem)
        self.verify_mapextent(npix_ref, blc_ref, trc_ref)
    
###
#
# Test case for checking if spline interpolation works for fast scan data
#
###
class sdimaging_test_interp(unittest.TestCase):
    """
    tests:
    test_spline_interp_single_infiles: check if spline interpolation works for single MS
    test_spline_interp_multiple_infiles: check if spline interpolation works for multiple MSs

    data:
    Both 'pointing6.ms' and 'pointing6-2.ms' contain 1000 rows for TP data, while only 10 
    rows given for POINTING data. 
    The pointing data is given as corner points of a hexagon centered at (RA, Dec) = 
    (0h00m00s, 0d00m00s) and with side of 0.001 radian and 0.0008 radian for 'pointing6.ms' 
    and 'pointing6-2.ms', respectively.
    The resulting pattern of weight image should be nearly circular if spline interpolation
    does work, while it should be hexagonal if linear interpolation, the old algorithm, is
    applied.
    Also, 'pointing6-2.ms' has 5 hours lag behind 'pointing6.ms'. 
    """
    datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdimaging/'
    params = dict(antenna = "0",
                  intent  = "*ON_SOURCE*",
                  gridfunction = "SF",
                  convsupport = 6,
                  imsize = [512, 512],
                  cell = "2arcsec",
                  phasecenter = "J2000 0:00:00.0 00.00.00.0",
                  ephemsrcname = '',
                  pointingcolumn = "direction",
                  stokes = 'I')
    infiles = []
    outfiles = [] # have a list of outfiles as multiple task execution may occur in a test

    def __remove_table(self, f):
        if os.path.exists(f):
            shutil.rmtree(f)
    
    def __copy_table(self, f):
        self.__remove_table(f)
        testutils.copytree_ignore_subversion(self.datapath, f)
        
    def setUp(self):
        self.infiles = []
        self.outfiles = []
        default(sdimaging)
        
    def tearDown(self):
        for infile in self.infiles:
            self.__remove_table(infile)
        for outfile in self.outfiles:
            os.system('rm -rf %s*'%(outfile))
        self.assertEqual(len(get_table_cache()), 0)

    def run_task(self, infiles, outfile, **kwargs):
        if isinstance(infiles, str):
            infiles = [ infiles ]
        for i in range(len(infiles)):
            self.infiles.append(infiles[i])
        self.outfiles.append(outfile)

        for infile in infiles:
            self.__copy_table(infile)
        self.params.update(**kwargs)

        status = sdimaging(infiles=infiles, outfile=outfile, **self.params)
        self.assertIsNone(status, msg = 'sdimaging failed to execute')
        self.assertTrue(os.path.exists(outfile+'.image'), msg='output image is not created.')
        
    def check_spline_works(self, outfile, multiple_ms=False):
        weightfile = outfile + '.weight'
        with tbmanager(weightfile) as tb:
            mapdata = tb.getcell('map', 0)
        # for pixels with strong weight value(>14), collect their distance from the image
        # center and then compute the mean and sigma of their distribution.
        dist_answer = [0.0, 0.0]
        dist_answer[0] = 0.001*180.0/numpy.pi*3600.0/float(self.params['cell'][0])
        dist_answer[1] = dist_answer[0]*0.8
        dist_sep = (dist_answer[0] + dist_answer[1])/2.0

        dist_list = [[], []]
        for i in range(self.params['imsize'][0]):
            for j in range(self.params['imsize'][1]):
                if mapdata[i][j][0][0] > 14.0:
                    cenx = float(self.params['imsize'][0])/2.0
                    ceny = float(self.params['imsize'][1])/2.0
                    dx = float(i) - cenx
                    dy = float(j) - ceny
                    dr = numpy.sqrt(dx*dx + dy*dy)
                    idx = 0 if (dist_sep < dr) else 1
                    dist_list[idx].append(dr)
        dist_mean1 = [0.0, 0.0]
        dist_mean2 = [0.0, 0.0]
        dist_sigma = [0.0, 0.0]
        dist_llim = [0.0, 0.0]
        dist_ulim = [0.0, 0.0]
        idx2 = 2 if multiple_ms else 1
        for i in range(idx2):
            for j in range(len(dist_list[i])):
                dist_mean1[i] += dist_list[i][j]
                dist_mean2[i] += dist_list[i][j] * dist_list[i][j]
            dist_mean1[i] = dist_mean1[i] / float(len(dist_list[i]))
            dist_mean2[i] = dist_mean2[i] / float(len(dist_list[i]))
            dist_sigma[i] = numpy.sqrt(dist_mean2[i] - dist_mean1[i] * dist_mean1[i])

            dist_llim[i] = dist_mean1[i] - dist_sigma[i]
            dist_ulim[i] = dist_mean1[i] + dist_sigma[i]

            """
            if spline interpolation is done, the range [dist_llim[0], dist_ulim[0]]
            will be a narrow range (102.683 ~ 103.318) and encloses the answer
            value (103.132), while linear interpolation will result in a wider
            range (94.240 +- 4.281) and depart from the answer value at 2-sigma
            level.
            FYI, [dist_llim[1], dist_ulim[1]] and dist_answer[1] will be 
            (81.609 - 83.151) and (82.506), respectively.
            """
            self.assertTrue(((dist_llim[i] < dist_answer[i]) and (dist_answer[i] < dist_ulim[i])),
                            msg = 'spline interpolation seems not working.')
            #print '['+str(i)+'] --- ' + str(dist_llim[i]) + ' - ' + str(dist_ulim[i])


    def check_images_identical(self, image1, image2, weight_image=False):
        suffix = '.weight' if weight_image else '.image'
        img1 = image1 + suffix
        img2 = image2 + suffix

        with tbmanager(img1) as tb:
            mapdata1 = tb.getcell('map', 0)
        with tbmanager(img2) as tb:
            mapdata2 = tb.getcell('map', 0)

        self.assertTrue(numpy.allclose(mapdata1, mapdata2, rtol=1.0e-5, atol=1.0e-5),
                        msg="%s and %s are not identical" % (img1, img2))

    def test_spline_interp_single_infiles(self):
        """test_spline_interp_single_infiles: Check if spline interpolation works for single fast-scan data."""
        outfile = 'pointing6.out'
        self.run_task(infiles=['pointing6.ms'], outfile=outfile)
        self.check_spline_works(outfile)

    def test_spline_interp_multiple_infiles(self):
        """test_spline_interp_multiple_infiles: Check if spline interpolation works for multiple fast-scan data."""
        outfile12 = "1and2.out"
        self.run_task(infiles=['pointing6.ms', 'pointing6-2.ms'], outfile=outfile12)
        outfile21 = "2and1.out"
        self.run_task(infiles=['pointing6-2.ms', 'pointing6.ms'], outfile=outfile21)

        #check if spline interpolation works
        self.check_spline_works(outfile12, True)
        #check if the results (both image and weight) don't change when infiles has inversed order
        self.check_images_identical(outfile12, outfile21)
        self.check_images_identical(outfile12, outfile21, True)


class sdimaging_test_interp_old(unittest.TestCase):
    """
    The test data 'pointing6.ms' contains 1000 rows for TP data, while only 10 rows given
    for POINTING data. The pointing data is given as corner points of a hexagon centered at
    (RA, Dec) = (0h00m00s, 0d00m00s) and with side of 0.001 radian.
    The resulting pattern of weight image should be nearly circular if spline interpolation
    does work, while it should be hexagonal if linear interpolation, the old algorithm, is
    applied.
    """
    datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdimaging/'
    params = dict(infiles = ['pointing6.ms'],
                  outfile = "pointing6.out",
                  antenna = "0",
                  intent  = "*ON_SOURCE*",
                  gridfunction = "SF",
                  convsupport = 6,
                  imsize = [512, 512],
                  cell = "2arcsec",
                  phasecenter = "J2000 0:00:00.0 00.00.00.0",
                  ephemsrcname = "Venus",
                  pointingcolumn = "direction")
    outfile = params['outfile']

    def __remove_table(self, f):
        if os.path.exists(f):
            shutil.rmtree(f)
    
    def __copy_table(self, f):
        self.__remove_table(f)
        testutils.copytree_ignore_subversion(self.datapath, f)
        
    def setUp(self):
        for infile in self.params['infiles']:
            self.__copy_table(infile)
        default(sdimaging)
        
    def tearDown(self):
        for infile in self.params['infiles']:
            self.__remove_table(infile)
        os.system('rm -rf %s*'%(self.outfile))
        
        self.assertEqual(len(get_table_cache()), 0)

    def run_test(self, **kwargs):
        self.params.update(**kwargs)
        status = sdimaging(**self.params)
        self.assertIsNone(status, msg = 'sdimaging failed to execute')
        outfile = self.outfile + image_suffix
        self.assertTrue(os.path.exists(outfile), msg='output image is not created.')
        
    def test_spline_interp(self):
        """test_spline_interp: Check if spline interpolation works for fast scan data."""
        
        self.run_test()

        weightfile = self.outfile + '.weight'
        with tbmanager(weightfile) as tb:
            mapdata = tb.getcell('map', 0)
        # for pixels with strong weight value(>14), collect their distance from the image
        # center and then compute the mean and sigma of their distribution.
        dist_list = []
        for i in range(self.params['imsize'][0]):
            for j in range(self.params['imsize'][1]):
                if mapdata[i][j][0][0] > 14.0:
                    cenx = float(self.params['imsize'][0])/2.0
                    ceny = float(self.params['imsize'][1])/2.0
                    dx = float(i) - cenx
                    dy = float(j) - ceny
                    dist_list.append(numpy.sqrt(dx*dx + dy*dy))
        dist_mean = 0.0
        dist_mean2 = 0.0
        for i in range(len(dist_list)):
            dist_mean += dist_list[i]
            dist_mean2 += dist_list[i] * dist_list[i]
        dist_mean = dist_mean / float(len(dist_list))
        dist_mean2 = dist_mean2 / float(len(dist_list))
        dist_sigma = numpy.sqrt(dist_mean2 - dist_mean * dist_mean)

        dist_llim = dist_mean - dist_sigma
        dist_ulim = dist_mean + dist_sigma
        dist_answer = 0.001*180.0/numpy.pi*3600.0/float(self.params['cell'][0])

        """
        if spline interpolation is done, the range [dist_llim, dist_ulim]
        will be a narrow range (102.683 ~ 103.318) and encloses the answer
        value (103.132), while linear interpolation will result in a wider
        range (94.240 +- 4.281) and depart from the answer value at 2-sigma
        level.
        """
        self.assertTrue(((dist_llim < dist_answer) and (dist_answer < dist_ulim)),
                        msg = 'spline interpolation seems not working.')
    
    
class sdimaging_test_clipping(sdimaging_unittest_base):
    """
    test_1row: check if clipping is not activated (1 spectrum)
    test_2rows: check if clipping is not activated (2 spectra)
    test_3rows: check if clipping is activated (3 spectra)
    test_multivis: check if clipping properly handles multiple ms inputs
    test_clip: check if clipping is applied to every image pixel separately
    test_clip2: check if clipping is activated on one pixel but is not on others
    test_suprious: check if clipping properly handles suprious data 
    test_multichan: check if clipping handles multi-channel data properly
    """
    data_list = ['clipping_1row.ms', 'clipping_2rows.ms', 'clipping_3rows.ms', 
                 'clipping_3rows_suprious.ms', 'clipping_3rows_2chans.ms']
    outfile = 'sdimaging_test_clipping.im'
    outfile_ref = 'sdimaging_test_clipping.ref.im'
    def setUp(self):
        default(sdimaging)
        
        # clear up test data
        self.__clear_up()
    
    def tearDown(self):
        # remove test data
        self.__clear_up()
        
        self.assertEqual(len(get_table_cache()), 0)

    def __clear_up(self):
        for data in self.data_list:
            if os.path.exists(data):
                shutil.rmtree(data)
        outfile = self.outfile + image_suffix
        if os.path.exists(outfile):
            shutil.rmtree(outfile)
            shutil.rmtree(self.outfile + '.weight')
            shutil.rmtree(self.outfile + '.sumwt')
            shutil.rmtree(self.outfile + '.psf')
        outfile_ref = self.outfile_ref + image_suffix
        if os.path.exists(outfile_ref):
            shutil.rmtree(outfile_ref)
            shutil.rmtree(self.outfile_ref + '.weight')
            shutil.rmtree(self.outfile_ref + '.sumwt')
            shutil.rmtree(self.outfile_ref + '.psf')
    
    def _test_clipping(self, infiles, is_clip_effective=True):
        if isinstance(infiles, str):
            self._test_clipping([infiles], is_clip_effective)
            return
        
        for infile in infiles:
            self.assertTrue(infile in self.data_list)
            self.assertFalse(os.path.exists(infile))
            testutils.copytree_ignore_subversion(self.datapath, infile)
            
        # image with clipping
        outfile = self.outfile
        overwrite = False
        mode = 'channel'
        nchan = -1
        start = 0
        width = 1
        gridfunction = 'BOX'
        imsize = 3
        cell = '1arcmin'
        phasecenter = 'J2000 0h0m0s 0d0m0s'
        sdimaging(infiles=infiles, outfile=outfile, overwrite=overwrite, 
                  mode=mode, nchan=nchan, start=start, width=width, 
                  gridfunction=gridfunction, imsize=imsize, cell=cell,
                  phasecenter=phasecenter, clipminmax=True)
        _outfile = outfile + image_suffix
        self.assertTrue(os.path.exists(_outfile))
        
        if is_clip_effective == True:
            # pre-flag the data to be clipped
            myme, mymsmd, mytb = gentools(['me', 'msmd', 'tb'])
            myqa = qa
            center = myme.direction('J2000', myqa.quantity(0, 'rad'), myqa.quantity(0, 'rad'))
            offset_plus = myqa.convert(myqa.quantity('1arcmin'), 'rad')
            offset_minus = myqa.mul(offset_plus, -1)
            grid = [[[], [], []],
                    [[], [], []],
                    [[], [], []]]
            gridmeta = [[[], [], []],
                        [[], [], []],
                        [[], [], []]]
            ra_list = [offset_plus['value'], 0, offset_minus['value']]
            dec_list = [offset_minus['value'], 0, offset_plus['value']]
            for infile in infiles:
                mymsmd.open(infile)
                try:
                    for irow in range(int(mymsmd.nrows())):
                        pointingdirection = mymsmd.pointingdirection(irow)['antenna1']['pointingdirection']
                        ra = pointingdirection['m0']['value']
                        dec = pointingdirection['m1']['value']
                        min_separation = 1e10
                        min_ra = -1
                        min_dec = -1
                        for ira in range(imsize):
                            for idec in range(imsize):
                                gra = ra_list[ira]
                                gdec = dec_list[idec]
                                separation = math.sqrt(math.pow(ra - gra, 2) + math.pow(dec - gdec, 2))
                                if separation < min_separation:
                                    min_ra = ira
                                    min_dec = idec
                                    min_separation = separation
                        gridmeta[min_ra][min_dec].append((infile, irow))
                finally:
                    mymsmd.close()
            
            print('### gridmeta', gridmeta)
            for ira in range(imsize):
                for idec in range(imsize):
                    meta = gridmeta[ira][idec]
                    for imeta in range(len(meta)):
                        infile, irow = meta[imeta]
                        mytb.open(infile)
                        try:
                            data = mytb.getcell('FLOAT_DATA', irow)[0]
                        finally:
                            mytb.close()
                        grid[ira][idec].append(data)
            
            for ira in range(imsize):
                for idec in range(imsize):
                    data = numpy.asarray(grid[ira][idec], dtype=numpy.float64)
                    if len(data) < 3:
                        continue
                    print('### ira', ira, 'idec', idec, 'data', data)
                    for ichan in range(data.shape[1]):
                        slice = data[:,ichan]
                        argmin = numpy.argmin(slice)
                        argmax = numpy.argmax(slice)
                        print('### ira', ira, 'idec', idec, 'argmin', argmin, 'argmax', argmax)
                        for imeta in (argmin, argmax):
                            infile, irow = gridmeta[ira][idec][imeta]
                            mytb.open(infile, nomodify=False)
                            try:
                                print('### clip', infile, 'row', irow, 'chan', ichan, 'data', mytb.getcell('FLOAT_DATA', irow))
                                #mytb.putcell('FLAG_ROW', irow, True)
                                flag = mytb.getcell('FLAG', irow)
                                print('### flag (before)', flag)
                                flag[0,ichan] = True
                                print('### flag (after)', flag)
                                mytb.putcell('FLAG', irow, flag)
                            finally:
                                mytb.close()
                    
        outfile = self.outfile_ref
        sdimaging(infiles=infiles, outfile=outfile, overwrite=overwrite, 
                  mode=mode, nchan=nchan, start=start, width=width, 
                  gridfunction=gridfunction, imsize=imsize, cell=cell,
                  phasecenter=phasecenter, clipminmax=False)
        _outfile_ref = outfile + image_suffix
        self.assertTrue(os.path.exists(_outfile_ref))
            
        # compare
        myia = gentools(['ia'])[0]
        myia.open(_outfile)
        result = myia.getchunk()
        result_mask = myia.getchunk(getmask=True)
        myia.close()
        
        myia.open(_outfile_ref)
        reference = myia.getchunk()
        reference_mask = myia.getchunk(getmask=True)
        myia.close()
        
        print('### result', result.flatten())
        print('### mask', result_mask.flatten())
        print('### reference', reference.flatten())
        print('### mask', reference_mask.flatten())
        
        self.assertTrue(numpy.all(result_mask == reference_mask))
        
        mresult = result[result_mask]
        mreference = reference[reference_mask]
        self.assertTrue(mresult.shape == mreference.shape)
        #self.assertTrue(numypy.all(result == reference))
        diff = lambda v, r: abs((v - r) / r) if r != 0.0 else abs(v)
        vdiff = numpy.vectorize(diff)
        err = vdiff(mresult, mreference)
        eps = 1.0e-6
        print('err = %s (max %s min %s)'%(err, err.max(), err.min()))
        self.assertTrue(numpy.all(err < eps)) 
    
    def test_1row(self):
        """test_1row: check if clipping is not activated (1 spectrum)"""
        infile = 'clipping_1row.ms'
        self._test_clipping(infile, is_clip_effective=False)
    
    def test_2rows(self):
        """test_2rows: check if clipping is not activated (2 spectra)"""
        infile = 'clipping_2rows.ms'
        self._test_clipping(infile, is_clip_effective=False)
    
    def test_3rows(self):
        """test_3rows: check if clipping is activated (3 spectra)"""
        infile = 'clipping_3rows.ms'
        self._test_clipping(infile, is_clip_effective=True)
    
    def test_multivis(self):
        """test_multivis: check if clipping properly handles multiple ms inputs"""
        infiles = ['clipping_1row.ms', 'clipping_2rows.ms']
        self._test_clipping(infiles, is_clip_effective=True)
        
    def test_clip(self):
        """test_clip: check if clipping is applied to every image pixel separately"""
        infiles = ['clipping_1row.ms', 'clipping_2rows.ms', 'clipping_3rows.ms']
        self._test_clipping(infiles, is_clip_effective=True)
        
    def test_clip2(self):
        """test_clip2: check if clipping is activated on one pixel but is not on others"""
        infiles = ['clipping_1row.ms', 'clipping_3rows.ms']
        self._test_clipping(infiles, is_clip_effective=True)
        
    def test_suprious(self):
        """test_suprious: check if clipping properly handles suprious data"""
        # This test is defined to verify new clipping algorithm
        #
        # Test data contains suprious. It is 10 orders of magnitude larger 
        # than orginary data so that ordinary data will disappear due to 
        # the loss of trailing digits when suprious data is accumulated to grid.
        # (NOTE: grid data is signle-precision)
        #
        # Old algorithm keeps track of minimum and maximum data during accumulation. 
        # However, it accumulates whole data once, then subtract minimum and maximum 
        # from accumulated result. In this procedure, suprious data must be accumulated 
        # to grid. Thus, the result is suffered from the loss of trailing digits.
        #  
        # On the other hand, new algorithm doesn't accumulate mininum and maximum. 
        # If clipping cannot apply (i.e., number of accumulated data is less than 
        # 3), these values are accumulated at the post-accumulation step.
        infile = 'clipping_3rows_suprious.ms'
        self._test_clipping(infile, is_clip_effective=True)
        
    def test_multichan(self):
        """test_multichan: check if clipping handles multi-channel data properly"""
        infile = 'clipping_3rows_2chans.ms'
        self._test_clipping(infile, is_clip_effective=True)
        

class sdimaging_test_projection(sdimaging_unittest_base):
    """
    Test projection

       - test_projection_GSL: unsupported projection type
       - test_projection_SIN: create image with SIN (Slant Orthographic) projection 
       - test_projection_TAN: create image with TAN (Gnomonic) projection
       - test_projection_CAR: create image with CAR (Plate Caree) projection
       - test_projection_SFL: create image with SFL (Sanson-Flamsteed) projection
       
    """
    # Input and output names
    prefix=sdimaging_unittest_base.taskname+'ProjectionTest'
    outfile=prefix+sdimaging_unittest_base.postfix
    mode = 'channel'
    cell = ['3.0arcmin', '3.0arcmin']
    imsize = [75, 75]
    phasecenter = 'J2000 17:18:29 +59.31.23'
    gridfunction = 'PB'
    start = 604
    nchan = 1

    keys=['max','maxpos','maxposf','mean','min','minpos','minposf',
          'npts','rms','blc','blcf','trc','trcf','sigma','sum','sumsq']
    
    def setUp(self):
        if os.path.exists(self.rawfile):
            shutil.rmtree(self.rawfile)
        shutil.copytree(self.datapath+self.rawfile, self.rawfile)
        # Common task parameters of the class
        self.task_param = dict(infiles=self.rawfile,mode=self.mode,
                               outfile=self.outfile,intent='OBSERVE_TARGET_ON_SOURCE',
                               cell=self.cell,imsize=self.imsize,
                               nchan=self.nchan,start=self.start,
                               phasecenter=self.phasecenter,
                               gridfunction=self.gridfunction)

        default(sdimaging)

    def tearDown(self):
        if (os.path.exists(self.rawfile)):
            shutil.rmtree(self.rawfile)
        os.system( 'rm -rf '+self.prefix+'*' )

        self.assertEqual(len(get_table_cache()), 0)
        
    def run_test_common(self, task_param, refstats, shape, refbeam=None,
                        atol=1.e-8, rtol=1.e-5, compstats=None, ignoremask=True,
                        projection='SIN'):
        
        # call super class's run_test_common 
        super(sdimaging_test_projection, self).run_test_common(task_param, refstats, shape, refbeam,
                                                               atol, rtol, compstats, ignoremask)

        # check projection
        outfile = task_param['outfile'] + image_suffix
        _ia.open(outfile)
        try:
            result_projection = _ia.coordsys().projection()['type']
        finally:
            _ia.close()
        self.assertEqual(projection, result_projection)
        
    def test_projection_GSL(self):
        """test_projection_GSL: unsupported projection type"""
        projection = 'GSL'
        spw = '0'
        self.task_param.update(dict(projection=projection, spw=spw))
        res=sdimaging(**self.task_param)
        self.assertFalse(res)
        self.assertFalse(os.path.exists(self.outfile))

    def test_projection_SIN(self):
        """test_projection_SIN: create image with SIN (Slant Orthographic) projection"""
        projection = 'SIN'
        spw = '0'
        self.task_param.update(dict(projection=projection, spw=spw))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats = {
            'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
            'blcf': '17:32:18.690, +57.37.28.536, I, 1.42087e+09Hz',
            'max': numpy.array([ 21.92034912]),
            'maxpos': numpy.array([59, 21,  0,  0], dtype=numpy.int32),
            'maxposf': '17:10:00.642, +58.42.19.808, I, 1.42087e+09Hz',
            'mean': numpy.array([ 7.84297146]),
            'min': numpy.array([ 3.36271787]),
            'minpos': numpy.array([71, 50,  0,  0], dtype=numpy.int32),
            'minposf': '17:04:49.308, +60.07.45.791, I, 1.42087e+09Hz',
            'npts': numpy.array([ 4217.]),
            'rms': numpy.array([ 8.70721651]),
            'sigma': numpy.array([ 3.7824345]),
            'sum': numpy.array([ 33073.81065345]),
            'sumsq': numpy.array([ 319714.46711966]),
            'trc': numpy.array([74, 74,  0,  0], dtype=numpy.int32),
            'trcf': '17:03:03.151, +61.19.10.757, I, 1.42087e+09Hz'
        }
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=False,
                             projection=projection)
        
    def test_projection_TAN(self):
        """test_projection_TAN: create image with TAN (Gnomonic) projection"""
        projection = 'TAN'
        spw = '0'
        self.task_param.update(dict(projection=projection, spw=spw))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats = {
            'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
            'blcf': '17:32:17.872, +57.37.35.824, I, 1.42087e+09Hz',
            'max': numpy.array([21.91863632]),
            'maxpos': numpy.array([59, 21,  0,  0], dtype=numpy.int32),
            'maxposf': '17:10:00.782, +58.42.20.655, I, 1.42087e+09Hz',
            'mean': numpy.array([ 7.84080757]),
            'min': numpy.array([ 3.36540604]),
            'minpos': numpy.array([71, 50,  0,  0], dtype=numpy.int32),
            'minposf': '17:04:49.729, +60.07.44.771, I, 1.42087e+09Hz',
            'npts': numpy.array([ 4222.]),
            'rms': numpy.array([ 8.7050746]),
            'sigma': numpy.array([ 3.78198999]),
            'sum': numpy.array([ 33103.88957095]),
            'sumsq': numpy.array([ 319936.08330953]),
            'trc': numpy.array([74, 74,  0,  0], dtype=numpy.int32),
            'trcf': '17:03:04.170, +61.19.04.235, I, 1.42087e+09Hz'
        }
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=False,
                             projection=projection)

    def test_projection_CAR(self):
        """test_projection_CAR: create image with CAR (Plate Caree) projection"""
        projection = 'CAR'
        spw = '0'
        self.task_param.update(dict(projection=projection, spw=spw))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats = {
            'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
            'blcf': '17:32:18.122, +57.37.29.933, I, 1.42087e+09Hz',
            'max': numpy.array([21.91925812]),
            'maxpos': numpy.array([59, 21,  0,  0], dtype=numpy.int32),
            'maxposf': '17:10:00.722, +58.42.19.922, I, 1.42087e+09Hz',
            'mean': numpy.array([ 7.84154849]),
            'min': numpy.array([ 3.36489725]),
            'minpos': numpy.array([71, 50,  0,  0], dtype=numpy.int32),
            'minposf': '17:04:49.481, +60.07.45.807, I, 1.42087e+09Hz',
            'npts': numpy.array([ 4219.]),
            'rms': numpy.array([ 8.70603491]),
            'sigma': numpy.array([ 3.78266474]),
            'sum': numpy.array([ 33083.49308872]),
            'sumsq': numpy.array([ 319779.29008623]),
            'trc': numpy.array([74, 74,  0,  0], dtype=numpy.int32),
            'trcf': '17:03:03.803, +61.19.09.870, I, 1.42087e+09Hz'
        }
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=False,
                             projection=projection)
        
    def test_projection_SFL(self):
        """test_projection_SFL: create image with SFL (Sanson-Flamsteed) projection"""
        projection = 'SFL'
        spw = '0'
        self.task_param.update(dict(projection=projection, spw=spw))
        outshape = (self.imsize[0],self.imsize[1],1,self.nchan)
        refstats = {
            'blc': numpy.array([0, 0, 0, 0], dtype=numpy.int32),
            'blcf': '17:32:18.553, +57.37.29.753, I, 1.42087e+09Hz',
            'max': numpy.array([21.91932678]),
            'maxpos': numpy.array([59, 21,  0,  0], dtype=numpy.int32),
            'maxposf': '17:10:00.673, +58.42.19.909, I, 1.42087e+09Hz',
            'mean': numpy.array([ 7.84234172]),
            'min': numpy.array([ 3.36329484]),
            'minpos': numpy.array([71, 50,  0,  0], dtype=numpy.int32),
            'minposf': '17:04:49.429, +60.07.45.787, I, 1.42087e+09Hz',
            'npts': numpy.array([ 4218.]),
            'rms': numpy.array([ 8.70668658]),
            'sigma': numpy.array([ 3.78252027]),
            'sum': numpy.array([ 33078.99737787]),
            'sumsq': numpy.array([ 319751.35842591]),
            'trc': numpy.array([74, 74,  0,  0], dtype=numpy.int32),
            'trcf': '17:03:03.322, +61.19.09.669, I, 1.42087e+09Hz'
        }
        self.run_test_common(self.task_param, refstats, outshape,
                             compstats=self.keys, ignoremask=False,
                             projection=projection)
    
"""
# utility for sdimaging_test_mapextent
# commented out since sd tool is no longer available in CASA (CAS-10301)
def get_mapextent(infile, scan=None):
    s = sd.scantable(infile, average=False)
    outfile = infile.rstrip('/') + '.tmp'
    try:
        s.save(outfile)
        if scan is None:
            with tbmanager(outfile) as tb:
                dir = tb.getcol('DIRECTION')
        else:
            with table_selector(outfile, taql='SCANNO==16') as tb:
                dir = tb.getcol('DIRECTION')
        rad2deg = lambda x: x * 180.0 / numpy.pi
        xmin = rad2deg(dir[0].min())
        xmax = rad2deg(dir[0].max())
        ymin = rad2deg(dir[1].min())
        ymax = rad2deg(dir[1].max())
        return numpy.array([xmax, ymin]), numpy.array([xmin, ymax])
    finally:
        if os.path.exists(outfile):
            shutil.rmtree(outfile)

def get_mapextent_ephemeris(infiles):
    mapcenter = None
    xmin = None
    xmax = None
    ymin = None
    ymax = None
    for infile in infiles:
        blc, trc = get_mapextent(infile)
        if mapcenter is None:
            mapcenter = 0.5 * (blc + trc)
        if xmin is None:
            xmin = trc[0]
        else:
            xmin = min(xmin, trc[0])
        if xmax is None:
            xmax = blc[0]
        else:
            xmax = max(xmax, blc[0])
        if ymin is None:
            ymin = blc[1]
        else:
            ymin = min(ymin, blc[1])
        if ymax is None:
            ymax = trc[1]
        else:
            ymax = max(ymax, trc[1])
    return numpy.array([xmax, ymin]), numpy.array([xmin, ymax])
"""

def str_to_deg(s):
    return qa.quantity(s)['value']

def calc_statistics(imagename):
    with toolmanager(imagename, 'ia') as ia:
        s = ia.statistics()
    return s
    
def calc_mapproperty(statistics):
    ra_in_deg = lambda x: qa.quantity(x.split(',')[0])['value']
    dec_in_deg = lambda x: qa.quantity(x.split(',')[1])['value']
    blcf = statistics['blcf']
    trcf = statistics['trcf']
    blcra = ra_in_deg(blcf)
    if blcra > 180.0:
        blcra -= 360.0
    blcdec = dec_in_deg(blcf)
    trcra = ra_in_deg(trcf)
    if trcra > 180.0:
        trcra -= 360.0
    trcdec = dec_in_deg(trcf)
    npix = statistics['trc'][:2] + 1
    dra = abs((trcra - blcra) * numpy.cos(0.5 * (blcdec + trcdec)))
    ddec = abs(trcdec - blcdec)
    return {'extent': numpy.array([dra, ddec]), 'npix': npix,
            'blc': numpy.array([blcra, blcdec]), 'trc': numpy.array([trcra, trcdec])}


def suite():
    return [sdimaging_test0,sdimaging_test1,
            sdimaging_test2,sdimaging_test3,
            sdimaging_autocoord,sdimaging_test_selection,
            sdimaging_test_flag,sdimaging_test_polflag,
            sdimaging_test_mslist,
            sdimaging_test_restfreq, 
            sdimaging_test_mapextent,
            sdimaging_test_interp,
            sdimaging_test_clipping,
            sdimaging_test_projection
            ]
