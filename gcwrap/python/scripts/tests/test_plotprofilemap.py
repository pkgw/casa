import os
import sys
import shutil
import numpy
import re
import imghdr
from __main__ import default
from taskinit import gentools
import unittest
import matplotlib
import pylab as pl

try:
    import testutils
except:
    import tests.testutils as testutils

from plotprofilemap import plotprofilemap
from exportfits import exportfits

(myia, myrg,) = gentools(['ia', 'rg'])

class plotprofilemap_test(unittest.TestCase):
    """
    This is a test suite for plotprofilemap task. 
    
    List of Tests:
        test_image_not_exist: input image does not exist (causes error)
        test_not_overwrite: output image already exists (causes error)
        test_pol_not_out_of_range: pol index is out of range (causes error)
        test_plotmasked_invalid: unsupported plotmasked value (causes error)
        test_numpanel_5x5: standard test (5x5 panels)
        test_numpanel_10x10: standard test (10x10 panels)
        test_plotmasked_empty: plotmasked is empty
        test_plotmasked_zero: plotmasked is zero
        test_plotmasked_text: plotmasked is text
        test_plotmasked_plot: plotmasked is plot
        test_export_image: test export the plot to PNG file
        test_fits_image: input image is FITS cube
    """
    # Data path of input/output
    datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/imregrid/'
    
    imagename_ref = 'expected.im'
    prefix = 'plotprofilemap_test'
    imagename = prefix + '.im'
    fitsimage = prefix + '.fits'
    figfile = prefix + '.png'
    
    standard_task_param = {'imagename': imagename,
                           'numpanels': '5,5',
                           'pol': 0,
                           'plotmasked': 'empty'}
    
    def setUp(self):
        # turn off interactive mode
        self.is_interactive = matplotlib.is_interactive()
        if self.is_interactive:
            pl.ioff()
        
        # copy input image
        testutils.copytree_ignore_subversion(self.datapath, self.imagename_ref, self.imagename)
        
        # make parameters default
        default(plotprofilemap)
        
        # copy standard task parameter set
        self.task_param = self.standard_task_param.copy()
        
        # make mask
        self.make_mask(self.imagename)
        
        # initial check
        self.assertTrue(os.path.exists(self.imagename))
        self.assertFalse(matplotlib.is_interactive())
        
    def tearDown(self):
        # restore original state
        if self.is_interactive:
            pl.ion()
        else:
            pl.ioff()
                
        # delete files 
        for f in [self.fitsimage, self.figfile]:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(self.imagename):
            shutil.rmtree(self.imagename)
            
    def make_mask(self, imagename):
        """
        make mask so that there is at least one blank panel when 
        profile map with 5x5 panels is created.
        """
        myia.open(imagename)
        try:
            imageshape = myia.shape()
            nx = imageshape[0] / 5 - 1
            ny = imageshape[1] / 5 - 1
            ns = imageshape[2] - 1
            print 'masked region: blc=[0,0,0], trc=[{0},{1},{2}]'.format(
                nx, ny, ns)
            region = myrg.box(blc=[0,0,0], trc=[nx,ny,ns])
            myia.set(pixelmask=False, region=region)
            msk = myia.getchunk(getmask=True)
            self.assertTrue(numpy.all(msk[:nx,:ny,:] == False))
        finally:
            myia.close()
        
    def run_task(self, **kwargs):
        self.task_param.update(kwargs)
        print self.task_param
        res = plotprofilemap(**self.task_param)
        return res
    
    def verify(self, numpanels, plotmasked=None):
        # get figure object
        figure = pl.gcf()
        
        # get list of axes object
        alist = figure.axes
        
        # parse numpanels
        s = numpanels.split(',')
        if len(s) == 1:
            nx = int(s[0])
            ny = int(s[0])
        elif len(s) == 2:
            nx = int(s[0])
            ny = int(s[1])
        else:
            self.assertFalse(True)
        
        # expected number of axes are 
        #   nx * ny (number of plot panels)
        #   + (nx + ny) (number of position labels)
        #   + 2 (number of position titles)
        expected_num_panels = nx * ny + nx + ny + 2
        self.assertEqual(len(alist), expected_num_panels)
        
        if plotmasked is not None:
            pass
        
    def verify_figfile(self, figfile):
        # figfile must exist
        self.assertTrue(os.path.exists(figfile))
        
        # figfile must be PNG format
        self.assertEqual(imghdr.what(figfile), 'png')
        
    def test_image_not_exist(self):
        """test_image_not_exist: input image does not exist (causes error)"""
        imagename = 'blabla.im'
        res = self.run_task(imagename=imagename)
        self.assertFalse(res)

    def test_not_overwrite(self):
        """test_not_overwrite: output image already exists (causes error)"""
        figfile = self.figfile
        
        # create figfile
        os.system('echo "" > {0}'.format(self.figfile))
        
        with self.assertRaises(RuntimeError) as cm:
            res = self.run_task(figfile=figfile, overwrite=False)
            self.assertFalse(res)
        the_exception = cm.exception
        #print 'Exception reported: "{0}"'.format(str(the_exception))
        self.assertTrue(re.match('^overwrite is False and output file exists:', str(the_exception)) is not None)

    def test_pol_not_out_of_range(self):
        """test_pol_not_out_of_range: pol index is out of range (causes error)"""
        # pol index is out of range
        pol = 1  
        
        with self.assertRaises(RuntimeError) as cm:
            res = self.run_task(pol=1)
            self.assertFalse(res)
        the_exception = cm.exception
        #print 'Exception reported: "{0}"'.format(str(the_exception))
        self.assertTrue(re.match('^pol {0} is out of range'.format(pol), str(the_exception)) is not None)

    def test_plotmasked_invalid(self):
        """test_plotmasked_invalid: unsupported plotmasked value (causes error)"""
        # invalid plotmasked value
        plotmasked = 'shadow'
        res = self.run_task(plotmasked=plotmasked)
        self.assertFalse(res)

    def test_numpanel_5x5(self):
        """test_numpanel_5x5: standard test (5x5 panels)"""
        numpanels = '5,5'
        
        res = self.run_task(numpanels=numpanels)
        
        self.verify(numpanels)

    def test_numpanel_10x10(self):
        """test_numpanel_10x10: standard test (10x10 panels)"""
        numpanels = '10,10'
        
        res = self.run_task(numpanels=numpanels)
        
        self.verify(numpanels)

    def test_plotmasked_empty(self):
        """test_plotmasked_empty: plotmasked is empty"""
        numpanels = '5,5'
        plotmasked = 'empty'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_zero(self):
        """test_plotmasked_zero: plotmasked is zero"""
        numpanels = '5,5'
        plotmasked = 'zero'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_text(self):
        """test_plotmasked_text: plotmasked is text"""
        numpanels = '5,5'
        plotmasked = 'text'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_plot(self):
        """test_plotmasked_plot: plotmasked is plot"""
        numpanels = '5,5'
        plotmasked = 'plot'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_export_image(self):
        """test_export_image: test export the plot to PNG file"""
        numpanels = '5,5'
        figfile = self.figfile
        
        res = self.run_task(numpanels=numpanels, figfile=figfile)
        
        self.verify(numpanels)
        self.verify_figfile(figfile)

    def test_fits_image(self):
        """test_fits_image: input image is FITS cube"""
        # convert input image to FITS
        self.assertFalse(os.path.exists(self.fitsimage))
        exportfits(imagename=self.imagename, fitsimage=self.fitsimage)
        self.assertTrue(os.path.exists(self.fitsimage))
        
        imagename=self.fitsimage
        numpanels = '5,5'
        
        res = self.run_task(imagename=imagename, numpanels=numpanels)
        
        self.verify(numpanels)


def suite():
    return [plotprofilemap_test]
