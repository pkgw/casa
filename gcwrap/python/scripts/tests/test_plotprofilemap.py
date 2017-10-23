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
    from . import testutils
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
        test_plotmasked_none: plotmasked is none
        test_export_image: test export the plot to PNG file
        test_fits_image: input image is FITS cube
        test_title: put title to the plot
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
            print('masked region: blc=[0,0,0], trc=[{0},{1},{2}]'.format(
                nx, ny, ns))
            region = myrg.box(blc=[0,0,0], trc=[nx,ny,ns])
            myia.set(pixelmask=False, region=region)
            msk = myia.getchunk(getmask=True)
            self.assertTrue(numpy.all(msk[:nx,:ny,:] == False))
        finally:
            myia.close()
        
    def run_task(self, **kwargs):
        self.task_param.update(kwargs)
        #print self.task_param
        res = plotprofilemap(**self.task_param)
        return res
    
    def _verify_axes_for_text(self, a, text=None):
        self.assertFalse(a.axison)
        lines = a.get_lines()
        self.assertEqual(len(lines), 0)
        texts = a.texts
        self.assertEqual(len(texts), 1)
        if text is not None:
            self.assertEqual(texts[0].get_text(), text)
            
    def _verify_plotmasked(self, a, plotmasked):
        axison = a.axison
        lines = a.get_lines()
        texts = a.texts
        
        if plotmasked == 'empty':
            # empty
            # show empty panel
            self.assertTrue(axison)
            self.assertEqual(len(lines), 0)
            self.assertEqual(len(texts), 0)
        elif plotmasked == 'zero':
            # zero
            # plot zero level
            self.assertTrue(axison)
            self.assertEqual(len(lines), 1)
            self.assertEqual(len(texts), 0)
            xdata, ydata = lines[0].get_data()
            self.assertTrue(numpy.all(ydata == 0.0))
        elif plotmasked == 'none':
            # none
            # show nothing
            self.assertFalse(axison)
            self.assertEqual(len(lines), 0)
            self.assertEqual(len(texts), 0)
        elif plotmasked == 'text':
            # text
            # show text 'NO DATA'
            self.assertTrue(axison)
            self.assertEqual(len(lines), 0)
            self.assertEqual(len(texts), 1)
            text = texts[0].get_text()
            self.assertEqual(text, 'NO DATA')
        elif plotmasked == 'plot':
            # plot
            # plot data with different color
            self.assertTrue(axison)
            self.assertEqual(len(lines), 1)
            self.assertEqual(len(texts), 0)
            xdata, ydata = lines[0].get_data()
            self.assertFalse(numpy.all(ydata == 0.0))
        else:
            self.fail('Invalid plotmasked value {0}'.format(plotmasked))
            
    
    def verify(self, numpanels, plotmasked=None, title=None):
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
        if title is not None:
            # if title is specified, the task adds additional axes for title
            expected_num_panels += 1
        self.assertEqual(len(alist), expected_num_panels)
        
        # verify axes
        # 0~nx-1: axes for horizontal axis label (right to left)
        # nx~nx+ny-1: axes for vertical axis label (bottom to top)
        # nx+ny: axes for horizontal axis title
        # nx+ny+1: axes for vertical axis title
        # nx+ny+2: axes for plot title (if title is specified)
        # subsequent: axes for plots (bottom right -> top right 
        #                              -> bottom of next column -> top of next column 
        #                              -> ... -> bottom left -> top left)
        # NOTE: (ny * (nx - 1) + 1)-th panel is empty 
        index = 0
        
        # axes for horizontal axis label
        #   - axison should be False
        #   - no line 
        #   - one text entry
        for i in range(nx):
            self._verify_axes_for_text(alist[index])
            index += 1
            
        # axes for vertical axis label
        #   - axison should be False
        #   - no line 
        #   - one text entry
        for i in range(ny):
            self._verify_axes_for_text(alist[index])
            index += 1
            
        # axes for horizontal axis title
        #   - axison should be False
        #   - no line 
        #   - one text entry
        #   - text should be 'Right Ascension (J2000)'
        for i in range(1):
            text = 'Right Ascension (J2000)'
            self._verify_axes_for_text(alist[index], text=text)
            index += 1

        # axes for vertical axis title
        #   - axison should be False
        #   - no line 
        #   - one text entry
        #   - text should be 'Declination (J2000)'
        for i in range(1):
            text = 'Declination (J2000)'
            self._verify_axes_for_text(alist[index], text=text)
            index += 1
        
        # axes for vertical axis title
        #   - axison should be False
        #   - no line 
        #   - one text entry
        #   - text should be title  
        if title is not None:
            for i in range(1):
                self._verify_axes_for_text(alist[index], text=title)
                index += 1
            
        # axes for profile map
        #   - axison should be True
        #   - one line (if not empty)
        #   - no text (if not empty)
        empty_panel = ny * (nx - 1) 
        for i in range(nx * ny):
            #print 'verify axes {0}'.format(index)
            a = alist[index]
            if plotmasked is None or i != empty_panel:
                self.assertTrue(a.axison)
                lines = a.get_lines()
                texts = a.texts
                self.assertEqual(len(lines), 1)
                self.assertEqual(len(texts), 0)
            else:
                print('verifying plotmasked parameter: axes {0} plotmasked {1}'.format(index, plotmasked))
                # verify plotmasked parameter 
                self._verify_plotmasked(a, plotmasked)
            index += 1

        
    def verify_figfile(self, figfile):
        # figfile must exist
        self.assertTrue(os.path.exists(figfile))
        
        # figfile must be PNG format
        self.assertEqual(imghdr.what(figfile), 'png')
        
    def skip_if_darwin(self):
        sysname = os.uname()[0]
        if sysname == 'Darwin':
            self.skipTest('Skip test_export_image on OS X since it may cause segfault')
            return
        
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
        self.skip_if_darwin()

        numpanels = '5,5'
        
        res = self.run_task(numpanels=numpanels)
        
        self.verify(numpanels)

    def test_numpanel_10x10(self):
        """test_numpanel_10x10: standard test (10x10 panels)"""
        self.skip_if_darwin()

        numpanels = '10,10'
        
        res = self.run_task(numpanels=numpanels)
        
        self.verify(numpanels)

    def test_plotmasked_empty(self):
        """test_plotmasked_empty: plotmasked is empty"""
        self.skip_if_darwin()

        # make mask
        self.make_mask(self.imagename)
        
        numpanels = '5,5'
        plotmasked = 'empty'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_zero(self):
        """test_plotmasked_zero: plotmasked is zero"""
        self.skip_if_darwin()

        # make mask
        self.make_mask(self.imagename)
        
        numpanels = '5,5'
        plotmasked = 'zero'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_text(self):
        """test_plotmasked_text: plotmasked is text"""
        self.skip_if_darwin()

        # make mask
        self.make_mask(self.imagename)
        
        numpanels = '5,5'
        plotmasked = 'text'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_plot(self):
        """test_plotmasked_plot: plotmasked is plot"""
        self.skip_if_darwin()

        # make mask
        self.make_mask(self.imagename)
        
        numpanels = '5,5'
        plotmasked = 'plot'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_plotmasked_none(self):
        """test_plotmasked_plot: plotmasked is none"""
        self.skip_if_darwin()

        # make mask
        self.make_mask(self.imagename)
        
        numpanels = '5,5'
        plotmasked = 'none'
        
        res = self.run_task(numpanels=numpanels, plotmasked=plotmasked)
        
        self.verify(numpanels, plotmasked)

    def test_export_image(self):
        """test_export_image: test export the plot to PNG file"""
        self.skip_if_darwin()
        
        numpanels = '5,5'
        figfile = self.figfile
        
        res = self.run_task(numpanels=numpanels, figfile=figfile)
        
        self.verify(numpanels)
        self.verify_figfile(figfile)

    def test_fits_image(self):
        """test_fits_image: input image is FITS cube"""
        self.skip_if_darwin()

        # convert input image to FITS
        self.assertFalse(os.path.exists(self.fitsimage))
        exportfits(imagename=self.imagename, fitsimage=self.fitsimage)
        self.assertTrue(os.path.exists(self.fitsimage))
        
        imagename=self.fitsimage
        numpanels = '5,5'
        
        res = self.run_task(imagename=imagename, numpanels=numpanels)
        
        self.verify(numpanels)
        
    def test_title(self):
        """test_title: put title to the plot"""
        self.skip_if_darwin()
        
        numpanels = '5,5'
        title = 'This is test image'
        
        res = self.run_task(numpanels=numpanels, title=title)
        
        self.verify(numpanels, title=title)


def suite():
    return [plotprofilemap_test]
