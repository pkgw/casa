import os
import sys
import shutil
import copy
from collections import namedtuple
import numpy
from scipy.optimize import curve_fit
from tasks import *
from taskinit import *
from casa_stack_manip import stack_frame_find
import unittest
from matplotlib.transforms import offset_copy

datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/sdsidebandsplit/'

# Gaussian fit
def gauss_func(x, *p):
    amp, center, width, offset = p
    y = amp*numpy.exp(-(x-center)**2/(2. * width**2)) + offset
    return y

def gauss_fit(x, y):
    # initial guess
    o = y.mean()
    a = numpy.abs(y - o).max()
    c = x[numpy.where( numpy.abs(y - o)==a)[0]]
    w = numpy.abs(y - o).sum()/a
    #print("initial guess: (%f, %f, %f, %f)" % (a,c,w,o))
    return curve_fit(gauss_func, x, y, p0=(a,c,w,o))

# a named tuple to store spectral information.
# start, end: start and end of a channel range
# max, min: max and min value of the channel range
# peak, center, width: gaussian fit parameters of the channel range
SpectralInfo = namedtuple('SpectralInfo',
                          ['start', 'end', 'max', 'min', 'peak', 'center', 'width', 'offset']) 


class sdsidebandsplitTestBase(unittest.TestCase):
    standard_param = dict(
        imagename = ['onepix_noiseless_shift0.image', 'onepix_noiseless_shift-102.image',
                     'onepix_noiseless_shift8.image', 'onepix_noiseless_shift62.image',
                     'onepix_noiseless_shift88.image', 'onepix_noiseless_shift100.image'],
        outfile='separated.image',
        overwrite = False,
        signalshift = [0.0, -102, +8, +62, +88, +100],
        imageshift = [0.0, 102, -8, -62, -88, -100],
        getbothside = False,
        refchan = 0.0,
        refval = '805GHz',
        otherside = False,
        threshold = 0.2
        )
    
    def update_task_param(self, new_param={}):
        """
        Overwrite standard task parameter and return a new dictionary
        with updated parameters.
        Note this task does not check validity of parameter names in
        the input parameter.
        
        Parameter
            new_param : a dictionary of parameter names (key) and
                        values (value) to overwrite standard task
                        execution parameters in tests.
        """
        if type(new_param) is not dict:
            raise TypeError('The input should be a dictionary')
        updated_param = copy.deepcopy(self.standard_param)
        updated_param.update(new_param)
        return updated_param
    
    def setUp(self):
        # copy input images
        for name in self.standard_param['imagename']:
            if os.path.exists(name):
                shutil.rmtree(name)
            shutil.copytree(datapath+name, name)
        # remove output files of previous run
        prefix = self.standard_param['outfile']
        for suffix in ['.signalband', '.imageband']:
            if os.path.exists(prefix+suffix):
                shutil.rmtree(prefix+suffix)

    def tearDown(self):
        # remove input images
        for name in self.standard_param['imagename']:
            if os.path.exists(name):
                shutil.rmtree(name)
                
        # remove output files
        prefix = self.standard_param['outfile']
        for suffix in ['.signalband', '.imageband']:
            if os.path.exists(prefix+suffix):
                shutil.rmtree(prefix+suffix)    

    def run_test(self, reference, **new_param):
        """
        Run sdsidebandsplit with given parameters and test result
        
        Arguments
            reference : a reference to compare results.
                        A dictionary with keys, 'signal' and 'image',
                        for signal and image sideband, respectively.
                        The data structure of 'signal' and 'image' values
                        depend on the implementation of test,
                        e.g., compare_image_data method.
            other key word arguments : test specific prameters to run tests
        """
        # Run task
        task_param = self.update_task_param(new_param)
        ret = sdsidebandsplit(**task_param)
        self.assertEqual(ret, None, 'The return value of task should be None')
        # Test results
        template_image = task_param['imagename'][0]
        self.assertTrue(os.path.exists(template_image),
                        "Could not find template image '%s'" % template_image)
        refcsys, refshape = self.get_image(template_image)
        self.assertTrue('signal' in reference,
                        'Internal Error: No valid reference value for signal sideband')
        # test signal band image
        imagename = task_param['outfile']+'.signalband'
        self.check_result(imagename, refcsys, refshape, reference['signal'])
        # test image band image
        imagename = task_param['outfile']+'.imageband'
        if task_param['getbothside']:
            self.assertTrue('image' in reference,
                            'Internal Error: No valid reference value for image sideband')
            # modify refcsys for image sideband
            spid = refcsys.findaxisbyname('spectral')
            refcsys.setreferencepixel(task_param['refchan'], 'spectral')
            myqa = qatool()
            refcsys.setreferencevalue(myqa.convert(task_param['refval'],
                                                   refcsys.units()[spid])['value'],
                                      'spectral')
            inc = refcsys.increment(format='n', type='spectral')['numeric'][0]
            refcsys.setincrement(-inc, 'spectral')
            
            self.check_result(imagename, refcsys, refshape, reference['image'])

    def get_image(self, imagename, getdata=False, getmask=False):
        """
        Returns image coordinate system object, shape.
        Optionally returns image pixel and mask values.
        Return values are in the order of
            csys, shape, data (optional), mask (optional).
        
        Parameters:
            imagename : the name of image
            getdata   : if True, returns image pixel values
            getmask   : if True, return image pixel mask
        """
        self.assertTrue(os.path.exists(imagename),
                        "Could not find image '%s'" % imagename)
        myia = gentools(['ia'])[0]
        myia.open(imagename)
        try:
            imcsys = myia.coordsys()
            imshape = myia.shape()
            if getdata:
                imdata = myia.getchunk()
            if getmask:
                immask = myia.getchunk(getmask=True)
        finally:
            myia.close()
        retval = [imcsys, imshape]
        if getdata: retval.append(imdata)
        if getmask: retval.append(immask)
        return retval

    def check_result(self, imagename, ref_csys, ref_shape, ref_value):
        """
        Compare an image with reference coordinate system, shape, and values.
        Details of tests shold be defined by methods called from this method,
        i.e., compare_image_coordinate and compare_image_data.
        
        Arguments
            imagename : the name of image to be tested
            ref_csys  : the reference coordinate system
            ref_shape : the reference of image shape
            ref_value : the data structure which defines image data
        """
        self.assertTrue(os.path.exists(imagename),
                        "Output image '%s' does not exist." % imagename)
        mycsys, myshape, mydata = self.get_image(imagename, getdata=True)
        self.compare_image_coordinate(mycsys, myshape, ref_csys, ref_shape)
        self.compare_image_data(mydata, ref_value)

    def compare_image_coordinate(self, csys, shape, ref_csys, ref_shape):
        """
        This method compares a coordinate system and shape with reference ones.
        The order of axes should be the same.
        Tested items:
        - dimension of shape
        - shape of each dimension
        - coordinate types of axes
        - axes units
        - reference pixel ids
        - reference values
        - increments
        """       
        # dimension
        self.assertEqual(len(shape), len(ref_shape), 'Dimension of shape differs from reference.')
        # dimension of csys
        self.assertEqual(ref_csys.naxes(), csys.naxes(),
                         'Dimension of coordinate system differs from reference')
        # confirm dimension of csys and shape
        self.assertEqual(len(shape), csys.naxes(),
                         'Dimention mismatch between shape and coordinate system')
        for i in range(len(ref_shape)):
            # shape of each dimension
            self.assertEqual(shape[i], ref_shape[i],
                             'Shape in %d-th dimension differs' % i)
            # axis type
            self.assertEqual(csys.axiscoordinatetypes()[i],
                             ref_csys.axiscoordinatetypes()[i],
                             'Axis coordinate type does not match in dimension %d' % i)
            # axis unit
            self.assertEqual(csys.units()[i], ref_csys.units()[i],
                             'Axis unit does not match in dimension %d' % i)
            # axis reference pixel
            self.assertAlmostEqual(csys.referencepixel()['numeric'][i],
                                   ref_csys.referencepixel()['numeric'][i],
                                   'Reference pixel does not match in dimension %d' % i)
            # axis reference value
            self.assertAlmostEqual(csys.referencevalue()['numeric'][i],
                                   ref_csys.referencevalue()['numeric'][i],
                                   ' does not match in dimension %d' % i)
            # axis increment
            self.assertAlmostEqual(csys.increment()['numeric'][i],
                                   ref_csys.increment()['numeric'][i],
                                   ' does not match in dimension %d' % i)

class failureTestCase(sdsidebandsplitTestBase):
    """
    A class to test invalid task parameters to run sdsidebandsplit.
    Implemented based on test case table attached to CAS-8091
    """
    def setUp(self):
        self.g = stack_frame_find()
        if '__rethrow_casa_exceptions' in self.g:
            self.rethrow_backup = self.g['__rethrow_casa_exceptions']
        else:
            self.rethrow_backup = None
        self.g['__rethrow_casa_exceptions'] = True
        super(failureTestCase, self).setUp()
    
    def tearDown(self):
        if self.rethrow_backup is None:
            self.g.pop('__rethrow_casa_exceptions')
        else:
            self.g['__rethrow_casa_exceptions'] = self.rethrow_backup
        super(failureTestCase, self).tearDown()
        del self.g

    def run_exception(self, ref_message, **new_param):
        """
        Run task and compare 
        """
        task_param = self.update_task_param(new_param)
        self.assertRaisesRegex(Exception, ref_message, sdsidebandsplit, **task_param)            

    # T-001
    def test_imagename_1image(self):
        """test failure: len(imagename)<2"""
        imagename = [ self.standard_param['imagename'][0] ]
        ref_message = 'At least two valid input data are required for processing'
        self.run_exception(ref_message, imagename = imagename)

    # T-005
    def test_imagename_invalidname(self):
        """test failure: len(imagename)==2 but includes an invalid imagename"""
        invalid_name = 'invalid.image'
        imagename = self.standard_param['imagename'][:-2] + [invalid_name]
        ref_message = 'Could not find %s' % invalid_name
        self.run_exception(ref_message, imagename = imagename)
        
    # T-006
    def test_outfile_undefined(self):
        """test failure: outfile is empty"""
        ref_message = 'Output file name is undefined.'
        self.run_exception(ref_message, outfile = '')

    # T-008, T-009
    def test_outfile_exists(self):
        """test failure: overwrite=F and outfile already exists."""
        for sideband in ('signalband', 'imageband'):
            print(('Test %s' % sideband))
            name = self.standard_param['outfile'] + '.' + sideband
            os.mkdir(name)
            ref_message = 'Image %s already exists.' % name
            param = dict(getbothside=(sideband == 'imageband'))
            self.run_exception(ref_message, **param)
            shutil.rmtree(name)

    # T-012
    def test_shifts_undefined(self):
        """test failure: both signalshift and imageshift are undefined"""
        ref_message = 'Channel shift was not been set.'
        self.run_exception(ref_message, signalshift=[], imageshift=[])
        
    # T-014, T-015, T-017, T-018
    def test_shift_wrong_length(self):
        """test failure: lengh of signalshift or imageshift does not match len(imagename)"""
        ref_message = "The number of shift should match that of images"
        for sideband in ['signalshift', 'imageshift']:
            myshift = self.standard_param[sideband] + [50]
            for shift in (myshift[:5], myshift):
                print(('Test len(%s)=%d' % (sideband, len(shift))))
                param = {sideband: shift}
                self.run_exception(ref_message, **param)   

    # T-022, T-023, T-024
    def test_refval_invalid(self):
        """test failure: refval is invalid (empty, a negative freqency or not a frequency)"""
        ref_message = ('refval is not a dictionary',
                       'Frequency should be positive',
                       'From/to units not consistent.')
        for refval, message in zip(('', '-100GHz', '300K'), ref_message):
            print(("Test refval='%s'" % refval))
            self.run_exception(message, refval=refval, getbothside=True)

    # T-027, T-031
    def test_threshold_outofrange(self):
        """test failure: threshold = 0.0, 1.0"""
        ref_message = 'Rejection limit should be > 0.0 and < 1.0'
        for thres in (0.0, 1.0):
            print(('Test threshold=%f' % thres))
            self.run_exception(ref_message, threshold=thres)
        
class standardTestCase(sdsidebandsplitTestBase):
    """
    A class to test valid task parameters to run sdsidebandsplit.
    Implemented based on test case table attached to CAS-8091
    The input images are synthesized spectra of 
    1 x 1 pixel, stokes I, 4080 channels.
    """

    standard_reference = dict(signal=(SpectralInfo(0, 1500, 4.06522, 0.99518, 2.96347, 898.4841, 30.48852, 1.08165),
                                      SpectralInfo(1700, 2700, 6.05933, 1.02671, 4.92205, 2297.6872, 19.75842, 1.14450)),
                              image=(SpectralInfo(1000, 2000, 8.07553, 1.05448, 6.94301, 1600.1052, 19.95953, 1.13069),
                                     SpectralInfo(2500, 3500, 3.06859, 1.00263, 1.99147, 2999.9953, 9.92538, 1.07988)))
    
    def assertAlmostEqual2(self, first, second, eps=1.0e-7, msg=None):
        if second == 0:
            self.assertLessEqual(abs(first), eps, msg)
        else:
            reldiff = abs((first - second) / second)
            self.assertLessEqual(reldiff, eps, msg)
 
    def compare_image_data(self, data, reference):
        """
        Compare image data with reference.

        Arguments
            data      : pixel value of image
            reference : a list of namedtuple, SpectralInfo,
                        which defines spectral feature of segments of spectrum
                        See the begining of this code about definition of SpectralInfo.
        """
        self.assertEqual(data.shape, (1,1,1,4080), 'Data shape is not expected one')
        for seg in reference:
            sp = data[0,0,0,seg.start:seg.end]
            x = list(range(seg.start, seg.end))
            #print('Max: ref {0} val {1}'.format(seg.max, sp.max()))
            #print('Min: ref {0} val {1}'.format(seg.min, sp.min()))
            self.assertAlmostEqual2(sp.max(), seg.max, 1e-3, 'Max comparison failed')
            self.assertAlmostEqual2(sp.min(), seg.min, 0.01, 'Min comparison failed')
            # compare gaussian fit
            fitp, _dummy = gauss_fit(x, sp)
            #print('Peak: ref {0} val {1}'.format(seg.peak, fitp[0]))
            #print('Peak Pos: ref {0} val {1}'.format(seg.center, fitp[1]))
            #print('Width: ref {0} val {1}'.format(seg.width, fitp[2]))
            #print('Offset: ref {0} val {1}'.format(seg.offset, fitp[3]))
            self.assertAlmostEqual2(fitp[0], seg.peak, 1e-3, 'Peak comparison failed')
            self.assertAlmostEqual2(fitp[1], seg.center, 1e-3, 'Peak position comparison failed')
            self.assertAlmostEqual2(numpy.abs(fitp[2]), numpy.abs(seg.width), 1e-3, 'Width comparison failed')
            self.assertAlmostEqual2(fitp[3], seg.offset, 1e-3, 'Offset comparison failed')
    
    # T-002
    def test_imagename_2images(self):
        """len(imagename)==2"""
        reference = dict(signal=[SpectralInfo(0, 1500, 4.0, 1.0, 2.86439, 898.70730, 30.33598, 1.09944),
                                 SpectralInfo(1500, 3000, 6.0, 1.0, 4.91510, 2297.94793, 19.555875, 1.10176954)])
        imagename = self.standard_param['imagename'][:2]
        signalshift = self.standard_param['signalshift'][:2]
        imageshift  = self.standard_param['imageshift'][:2]
        self.run_test(reference, imagename=imagename,
                              signalshift=signalshift, imageshift=imageshift)

    # T-007
    def test_imagename_6images(self):
        """standard run: valid outfile, len(imagename)==6"""
        self.run_test(self.standard_reference)

    # T-010
    def test_imageband_exists_signalonly(self):
        """imageband image exists but only signal band is solved (must succeed)"""
        imageband = self.standard_param['outfile'] + '.imageband'
        os.mkdir(imageband)
        self.assertTrue(os.path.exists(imageband), "Failed to create '%s'" % imageband)
        self.run_test(self.standard_reference, getbothside=False, overwrite=False)

    # T-011
    def test_overwrite(self):
        """overwrite = True"""
        for sideband in ['.signalband', '.imageband']:
            name = self.standard_param['outfile'] + sideband
            os.mkdir(name)
            self.assertTrue(os.path.exists(name), "Failed to create '%s'" % name)
        self.run_test(self.standard_reference, overwrite=True, getbothside=True)

    # T-013
    def test_signalshift_from_imageshift(self):
        """obtain signalshift from imageshift"""
        self.run_test(self.standard_reference, signalshift=[])

    # T-016
    def test_imageshift_from_signalshift(self):
        """obtain imageshift from signalshift"""
        self.run_test(self.standard_reference, imageshift=[])

    # T-019
    def test_getbothside(self):
        """getbothside = True"""
        self.run_test(self.standard_reference, getbothside=True)

    # T-020
    def test_refchan_negative(self):
        """refchan = -1.0"""
        self.run_test(self.standard_reference, getbothside=True, refchan=-1.0)

    # T-021
    def test_refchan_large(self):
        """refchan > nchan"""
        self.run_test(self.standard_reference, getbothside=True, refchan=5000.0)

    # T-025
    def test_otherside(self):
        """otherside = True"""
        reference = dict(signal=(SpectralInfo(0, 1500, 2.77938, -0.198638, 2.90884, 898.4125, 29.57543, -0.12153),
                                 SpectralInfo(1500, 3000, 4.74292, -0.20648, 4.88534, 2298.0446, 19.75160, -0.13152)),
                         image=(SpectralInfo(1000, 2200, 6.67820, -0.24841, 6.84674, 1599.9131, 19.75747, -0.15339),
                                SpectralInfo(2500, 3500, 1.88367, -0.11315, 1.96069, 2999.8335, 9.94222, -0.07489)))
        for doboth in [False, True]:
            print(('getbothside = %s'% str(doboth)))
            self.run_test(reference, otherside=True, getbothside=doboth, overwrite=True)

    # T-028, T-029, T-030
    def test_threshold(self):
        """various threshold values"""
        ref_small = dict(signal=(SpectralInfo(0, 1500, 4.09385, 1.08961, 2.99224, 898.04215, 29.99680, 1.09845),
                                 SpectralInfo(1500, 3000, 6.08972, 1.08962, 4.99175, 2297.9978, 19.98419, 1.09853)))
        ref_mid = dict(signal=(SpectralInfo(0, 1500, 4.06263, 0.96078, 2.97665, 898.26048, 29.81733, 1.07312),
                               SpectralInfo(1500, 3000, 6.02824, 0.99354, 4.87862, 2297.8214, 19.10228, 1.19242)))
        ref_large = dict(signal=(SpectralInfo(0, 1500, 3.99807, 0.97891, 2.96490, 898.00416, -29.48243, 1.04387),
                               SpectralInfo(1500, 3000, 5.99822, 0.96876, 4.83709, 2298.0035, -18.79443, 1.21969)))
        for val, ref in zip((0.0001, 0.5, 0.9999), (ref_small, ref_mid, ref_large)):
            print(('Threshold=%f' % val))
            self.run_test(ref, threshold=val, overwrite=True)

class MultiPixTestCase(sdsidebandsplitTestBase):
    """
    A class to test sdsidebandsplit with multi-pixel images.
    Implemented based on test case table attached to CAS-8091 (T-032)
    """
    standard_param = dict(
        imagename = ['multipix_noiseless_shift0.image', 'multipix_noiseless_shift-102.image',
                     'multipix_noiseless_shift8.image', 'multipix_noiseless_shift62.image',
                     'multipix_noiseless_shift88.image', 'multipix_noiseless_shift100.image'],
        outfile='separated.image',
        overwrite = False,
        signalshift = [0.0, -102, +8, +62, +88, +100],
        imageshift = [0.0, 102, -8, -62, -88, -100],
        getbothside = False,
        refchan = 0.0,
        refval = '805GHz',
        otherside = False,
        threshold = 0.2
        )

    def compare_image_data(self, data, reference):
        """
        Compare image data with reference.

        Arguments
            data      : pixel value of image
            reference : a reference image name
        """
        self.assertTrue(os.path.exists(reference),
                        "Could not find reference image '%s'" % reference)
        myia = gentools(['ia'])[0]
        myia.open(reference)
        try:
            ref_data = myia.getchunk()
        finally:
            myia.close()
        self.assertEqual(data.shape, ref_data.shape, 'Image shape comparison failed')
        self.assertAlmostEqual(data.max(), ref_data.max(), 3, 'Max comparison failed')
#         self.assertEqual(numpy.where(data==data.max()),
#                          numpy.where(ref_data==ref_data.max()),
#                          'Max position comparison failed')
        self.assertAlmostEqual(data.min(), ref_data.min(), 3, 'Min comparison failed')
#         self.assertEqual(numpy.where(data==data.min()),
#                          numpy.where(ref_data==ref_data.min()),
#                          'Max position comparison failed')
        self.assertAlmostEqual(data.std(), ref_data.std(), 3, 'StdDev comparison failed')

    # T-032
    def test_multi_pixels(self):
        """images with 10x10 spatial pixel"""
        reference = dict(signal = datapath+'ref_multipix.signalband',
                         image = datapath+'ref_multipix.imageband')
        self.run_test(reference, getbothside=True)

def suite():
    return [failureTestCase, standardTestCase, MultiPixTestCase]