import os
import sys
import shutil
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import fnmatch
import sha
import time
import matplotlib.pyplot as plt
import numpy as np

# Paths for data
datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/plotms/"
altdatapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/setjy/"
calpath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/gaincal/"
overlaypath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/mstransform/"

# Pick up alternative data directory to run tests on MMSs
if os.environ.has_key('TEST_DATADIR'):   
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/plotms/'
    if os.path.isdir(DATADIR):
        datapath = DATADIR         

print 'plotms tests will use data from '+ datapath

class plotms_test_base(unittest.TestCase):

    testms  = "pm_ngc5921.ms"
    testms2 = "ngc5921.ms"
    testms3 = "sun.subset.pentagon.ms"
    testms4 = "split_ddid_mixedpol_CAS-12283.ms"
    testcaltable = 'ngc5921.ref1a.gcal'
    outputDir="/tmp/" + str(os.getpid()) + "/"
    plotfile_jpg = "/tmp/myplot.jpg"
    display = os.environ.get("DISPLAY")
    ms = os.path.join(outputDir, testms)
    ms2 = os.path.join(outputDir, testms2)
    ms3 = os.path.join(outputDir, testms3)
    ms4 = os.path.join(outputDir, testms4)
    caltable = os.path.join(outputDir, testcaltable)

    def cleanUp(self):
        if os.path.exists(self.outputDir):
            shutil.rmtree(self.outputDir)

    def setUpData(self):
        res = None
        default(plotms)
        if not os.path.exists(self.ms):
            shutil.copytree(os.path.join(datapath,self.testms), 
                    self.ms, symlinks=True)

    def tearDownData(self):
        self.cleanUp()
        pm.setPlotMSFilename("")

    def setUpAltData(self):
        if not os.path.exists(self.ms2):
            shutil.copytree(os.path.join(altdatapath,self.testms2),
                    self.ms2, symlinks=True)

    def setUpCalData(self):
        res = None
        default(plotms)
        if not os.path.exists(self.ms2):
            shutil.copytree(os.path.join(calpath,self.testms2), 
                    self.ms2, symlinks=True)
        testcaltable = os.path.join(self.outputDir, self.caltable)
        if not os.path.exists(self.caltable):
            shutil.copytree(os.path.join(calpath, self.testcaltable),
                    self.caltable, symlinks=True)

    def setUpPointingData(self):
        if not os.path.exists(self.ms3):
            shutil.copytree(os.path.join(datapath,self.testms3),
                    self.ms3, symlinks=True)

    def setUpOverlayData(self):
        if not os.path.exists(self.ms4):
            shutil.copytree(os.path.join(overlaypath,self.testms4),
                    self.ms4, symlinks=True)

    def checkPlotfile(self, plotfileName, minSize, maxSize=None):
        self.assertTrue(os.path.isfile(plotfileName), "Plot was not created")
        plotSize = os.path.getsize(plotfileName)
        print plotfileName, 'file size is', plotSize
        self.assertGreater(plotSize, minSize)
        if maxSize:
            self.assertLess(plotSize, maxSize)
    
    # plotms can return True if it catches an error, but makes no plot file
    def checkNoPlotfile(self, plotfileName):
        self.assertFalse(os.path.isfile(plotfileName), "Plot was created")

    def getFilecount(self, dirName, namePattern ):
        nameTarget = namePattern + '*'
        count = 0
        for  file in os.listdir( dirName ):
            if fnmatch.fnmatch( file, nameTarget):
                count = count + 1
        return count
    
    def removeFiles(self, dirName, namePattern ):
        nameTarget = namePattern + '*'
        for  file in os.listdir( dirName ):
            if fnmatch.fnmatch( file, nameTarget):
                os.remove(os.path.join(dirName,file))

    def removePlotfile(self, plotfile=None):
        try:
            if not plotfile:
                plotfile = self.plotfile_jpg
            os.remove(plotfile)
        except OSError:  # "No such file or directory"
            pass
            
    def checkDisplay(self):
        self.assertGreater(len(self.display), 0, 'DISPLAY not set, cannot run test')

# ------------------------------------------------------------------------------
 
class test_basic(plotms_test_base):
    ''' Test basic single plots and overplots '''

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        # alternate data needed for plotting two different MSes
        self.setUpAltData()
        
    def tearDown(self):
        self.tearDownData()
            
    def test_basic_plot(self):
        '''test_basic_plot: Basic plot'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic01.jpg")
        self.removePlotfile()
        time.sleep(5)
        # default axes
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        self.removePlotfile()
        # default xaxis only
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat="jpg", 
                     yaxis='freq', showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 240000)
        self.removePlotfile()
        # default yaxis only
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat="jpg", 
                     xaxis='scan', showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        print
        
    def test_basic_blankplot(self):               
        '''test_basic_blankplot: Blank plot running plotms with no arguments'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic02.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms( showgui=False, plotfile=self.plotfile_jpg, expformat='jpg',
                highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 22000)
        print    
           
    def test_basic_overwrite(self):
        '''test_basic_overwrite: Check overwrite functionality'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic03.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Create plotfile
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)

        # Next, overwrite is False so the save should fail.
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, highres=True)   
        self.assertFalse(res)

        # Next, overwrite is True so the save should succeed.
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        print

    def test_basic_overplot2MS(self):
        '''test_basic_overplot2MS: Overplot two data sets on one plot.'''
        if os.path.exists(os.path.join(self.outputDir, self.ms2)):
            self.plotfile_jpg = os.path.join(self.outputDir, "testBasic04.jpg")
            self.removePlotfile()
            time.sleep(5)
            # Plot first MS scan vs time
            res = plotms(vis=self.ms, 
                         showgui=False, yaxis='scan',
                         showlegend=True, legendposition='lowerRight',
                         customsymbol=[True], symbolshape=['diamond'], symbolsize=[3],
                         symbolcolor=['ff0000'], symbolfill=['mesh3'], highres=True)
            self.assertTrue(res)
            # Plot second MS field vs time and export it
            res = plotms(vis=self.ms2, plotfile=self.plotfile_jpg, expformat='jpg', 
                         showgui=False, yaxis='field',
                         rowindex=0, colindex=0, plotindex=1, clearplots=False,
                         showlegend=True, legendposition='lowerRight',
                         customsymbol=[True], symbolshape=['circle'], symbolsize=[3],
                         symbolcolor=['00FF00'], symbolfill=['mesh3'], highres=True)   
            self.assertTrue(res) 
            self.checkPlotfile(self.plotfile_jpg, 45000)
        else:
            print "Skipping test, no path to alternate MS"
        print

    def test_basic_overplot2MS_freq(self):
        '''test_basic_overplot2MS_freq: CAS-6975 overplotting problem'''
        if os.path.exists(os.path.join(self.outputDir, self.ms2)):
            self.plotfile_jpg = os.path.join(self.outputDir, "testBasic05.jpg")
            self.removePlotfile()
            time.sleep(5)
            # Create the first plot
            res = plotms(vis=self.ms, showgui=False, highres=True,
                         xaxis="freq", yaxis="phase", avgchannel="63")
            self.assertTrue(res)
            # Do an overplot with a different file
            res = plotms(vis=self.ms2, showgui=False, plotfile=self.plotfile_jpg,
                         xaxis="freq", yaxis="phase", avgchannel="63",
                         plotindex=1, clearplots=False, highres=True) 
            self.assertTrue(res)
            self.checkPlotfile(self.plotfile_jpg, 30000)
        else:
            print "Skipping test, no path to alternate MS"
        print

    def test_basic_overplot2colors(self):
        '''test_basic_overplot2colors: CAS-7043 Create overplot with different color for each'''
        if os.path.exists(os.path.join(self.outputDir, self.ms2)):
            self.plotfile_jpg = os.path.join(self.outputDir, "testBasic06.jpg")
            self.removePlotfile()
            time.sleep(5)
            # Create the first plot
            res = plotms(vis=self.ms, showgui=False, highres=True,
                         xaxis="time", yaxis="amp", antenna="1", avgchannel="1000",
                         customsymbol=True, symbolshape='diamond', 
                         symbolsize=5,symbolcolor='00ff00')
            self.assertTrue(res)
            # Do an overplot with a different file
            res = plotms(vis=self.ms2, showgui=False, plotfile=self.plotfile_jpg,
                         xaxis="time", yaxis="amp", antenna="!1",avgchannel="1000",
                         plotindex=1, clearplots=False,
                         customsymbol=True, symbolshape='diamond',
                         symbolsize=1, symbolcolor='0000ff', highres=True)
            self.assertTrue(res)
            self.checkPlotfile(self.plotfile_jpg, 40000)
        else:
            print "Skipping test, no path to alternate MS"
        print

    def xtest_basic_screenExport(self):
        '''test_basic_screenExport: Export plot in screen resolution'''
        # not working yet...
        self.plotfile_jpg = os.path.join(self.outputDir, "testBasic07.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=False)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 30000)
        print

    def xtest_basic_pngExport(self):
        '''test_basic_pngExport: Export plot in png format'''
        plotfile_png = os.path.join(self.outputDir, "testBasic08.png")
        self.removePlotfile(plotfile_png)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_png, expformat="png", 
                     showgui=False, highres=False)   
        self.assertTrue(res)
        self.checkPlotfile(plotfile_png, 30000)
        self.removePlotfile(plotfile_png)
        print

    def xtest_basic_pdfExport(self):
        '''test_basic_pdfExport: Export plot in pdf format'''
        plotfile_pdf = os.path.join(self.outputDir, "testBasic09.pdf")
        self.removePlotfile(plotfile_pdf)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_pdf, expformat="pdf", 
                     showgui=False, highres=False)   
        self.assertTrue(res)
        self.checkPlotfile(plotfile_pdf, 30000)
        self.removePlotfile(plotfile_pdf)
        print



# ------------------------------------------------------------------------------

class test_averaging(plotms_test_base):
    ''' test plotms averaging options '''

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()
    
    def test_averaging_time(self):
        '''test_averaging_time: Average time'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging01.jpg")
        self.removePlotfile()
        time.sleep(5)
        # interval = 30s
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg,
                     showgui=False, avgtime='60', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000, 310000)
        print
    
    def test_averaging_timescan(self):
        '''test_averaging_timescan: Average time over scans'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging02.jpg")
        self.removePlotfile()
        time.sleep(5)
        # interval = 30s
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgtime='120', avgscan=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 140000, 210000)
        print

    def test_averaging_timefield(self):
        '''test_averaging_timefield: Average time over fields'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging03.jpg")
        self.removePlotfile()
        time.sleep(5)
        # interval = 30s
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgtime='120', avgfield=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 140000, 210000)
        print

    def test_averaging_chan(self):
        '''test_averaging_chan: Average channel'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging04.jpg")
        self.removePlotfile()
        time.sleep(5)
        # nchan = 63
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgchannel='7')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 150000, 230000)
        print

    def test_averaging_baseline(self):
        '''test_averaging_baseline: Average over baseline'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging05.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgbaseline=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 140000, 210000)
        print

    def test_averaging_antenna(self):
        '''test_averaging_antenna: Average per antenna'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging06.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgantenna=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 160000, 310000)
        print

    def test_averaging_blnant(self):
        '''test_averaging_blnant: Average over baseline and per antenna (should fail)'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging07.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgbaseline=True, avgantenna=True)
        self.assertFalse(res)
        print

    def test_averaging_spw(self):
        '''test_averaging_spw: Average over spw'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAveraging08.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, avgspw=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000, 310000)
        print

# ------------------------------------------------------------------------------

class test_axis(plotms_test_base):
    ''' test axis and datacolumn options '''

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        self.setUpPointingData()
        
    def tearDown(self):
        self.tearDownData()

    def test_axis_list(self):
        '''test_axis_list: Overplot scan and field y-axis with time x-axis.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis01.jpg")
        self.removePlotfile()
        time.sleep(5)
        '''Create a (scan & field)/time plot'''
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='time', yaxis=['scan','field'],
                     showlegend=True, legendposition='lowerRight',
                     customsymbol=[True,True], symbolshape=['diamond','circle'], 
                     symbolsize=[5,5], symbolcolor=['ff0000','00ff00'], 
                     symbolfill=['mesh3','mesh3'], highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 40000)
        print
        
    def test_axis_twoAxes(self):
        '''test_axis_twoAxes: Single plot with two y axes/locations.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis02.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, gridrows=1, gridcols=1, 
                     showgui=False, yaxis=['amp','scan'], 
                     yaxislocation=['left','right'],xaxis='time',
                     plotfile=self.plotfile_jpg, expformat='jpg',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 200000)
        self.removePlotfile()
        
        print

    def test_axis_wtamp(self):
        '''test_axis_wtamp: Test that wt*amp works for x-and y-axis choices.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis03.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Use a 2x1 grid and plot wt*amp vs time on with one graph
        # and time vs wt*amp on other
        res = plotms(vis=self.ms, gridrows=2, gridcols=1,
                     showgui=False, yaxis='wt*amp', xaxis='time')
        self.assertTrue(res)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='wt*amp', yaxis='time',
                     plotindex=1, rowindex=1, colindex=0, clearplots=False,
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 210000)
        print

    def test_axis_nodatacol(self):
        '''test_axis_nodatacol: Test non-visibility data (no datacolumn) x vs. y'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis04.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, xaxis='elevation', yaxis='azimuth')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 30000)
        print

    def test_axis_datacolumns(self):
        '''test_axis_datacolumns: Test datacolumn options'''
        datacols = ['data', 'corrected', 'model', 'residual',
                    'corrected-model', 'corrected-model_vector', 'corrected-model_scalar',
                    'data-model', 'data-model_vector', 'data-model_scalar',
                    'corrected/model', 'corrected/model_vector', 'corrected/model_scalar',
                    'data/model', 'data/model_vector', 'data/model_scalar']
        for datacol in datacols:
            filename = "testAxis05_datacolumn.jpg"
            plotfile = os.path.join(self.outputDir, filename)
            self.removePlotfile(plotfile)
            res = plotms(vis=self.ms, plotfile=plotfile, highres=True,
                         showgui=False, ydatacolumn=datacol)
            self.assertTrue(res)
            minSize = 180000
            if datacol is 'model':
                minSize = 30000
            self.checkPlotfile(plotfile, minSize)
            self.removePlotfile(plotfile)
        print

    def test_axis_baddatacolumn(self):
        '''test_axis_baddatacolumn: Test invalid datacolumn'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis06.jpg")
        self.removePlotfile()
        time.sleep(5)
        try:
            res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, ydatacolumn='corr/model')
            self.assertFalse(res)
        except RuntimeError:  # "Parameter verification failed" rethrown on test server
            pass
        print

    def test_axis_datacolumnNoFloat(self):
        '''test_axis_datacolumnNoFloat: Choose 'float' datacolumn that does not exist'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis07.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, ydatacolumn='float')
        self.assertFalse(res)
        print

    def test_axis_synonyms(self):
        '''test_axis_synonym: Test shortened form of axis names'''
        # Just test a commonly-used subset
        syns = ['chan', 'freq', 'vel', 'ant', 'ant1', 'ant2',
                'imag', 'vel', 'wtsp']
        for syn in syns:
            filename = "testAxis08_" + syn + ".jpg"
            plotfile = os.path.join(self.outputDir, filename)
            self.removePlotfile(plotfile)
            res = plotms(vis=self.ms, plotfile=plotfile, highres=True,
                         showgui=False, yaxis=syn)
            self.assertTrue(res)
            self.checkPlotfile(plotfile, 30000)
            self.removePlotfile(plotfile)
        # test when yaxis is list
        res = plotms(vis=self.ms, plotfile=plotfile, highres=True,
            showgui=False, yaxis=['chan','freq'])
        self.checkPlotfile(plotfile, 30000)
        self.removePlotfile(plotfile)
        print

    def test_axis_syn_bad(self):
        '''test_axis_syn_bad: Test invalid axis synonym'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis09.jpg")
        self.removePlotfile()
        # now test bad synonym
        try:
            res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, highres=True,
                     showgui=False, yaxis='veloc')
            self.assertFalse(res)
        except RuntimeError:  # "Parameter verification failed" rethrown on test server
            pass
        print

    # Removed from regular testing but available for manual testing!
    def xtest_axis_ALL(self):
        '''test_axis_ALL: Test ALL axis names'''
        axes = ['scan', 'field', 'time', 'interval', 'spw', 'chan',
                'freq', 'vel', 'corr', 'ant1', 'ant2', 'baseline',
                'row', 'amp', 'phase', 'real', 'imag', 'wt', 'wtsp',
                'sigma', 'sigmasp', 'flag', 'flagrow', 'uvdist', 
                'uvwave', 'u', 'v', 'w', 'uwave', 'vwave', 'wwave', 
                'azimuth', 'elevation', 'hourang', 'parang', 
                'ant', 'ant-azimuth', 'ant-elevation', 'ant-parang', 'ant-ra', 'ant-dec',
                'observation', 'intent']
        for axis in axes:
            filename = "testAxis10_" + axis + ".jpg"
            plotfile = os.path.join(self.outputDir, filename)
            self.removePlotfile(plotfile)
            axis_vis = self.ms
            if axis in ['ant-ra','ant-dec']:
                axis_vis = self.ms3
            res = plotms(vis=axis_vis, plotfile=plotfile, highres=True,
                         showgui=False, yaxis=axis)
            self.assertTrue(res)
            self.checkPlotfile(plotfile, 40000)
            self.removePlotfile(plotfile)
        print

    def test_axis_twoDataAxes(self):
        '''test_axis_twoDataAxes: Single plot with two y data axes.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testAxis11.jpg")
        self.removePlotfile()
        time.sleep(5)
        # cannot plot with duplicate y-axes (same axis, same data)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg',
                     overwrite=True, showgui=False, yaxis=['amp','amp'], 
                     yaxislocation=['left','right'], xaxis='time', highres=True)
        self.assertFalse(res) 
        # Plot amp:data vs time and amp:model vs time. 
        res = plotms(vis=self.ms, gridrows=1, gridcols=1, 
                     showgui=False, yaxis=['amp','amp'],
                     ydatacolumn=['data','model'],
                     yaxislocation=['left','right'],xaxis='time',
                     showlegend=True, legendposition='lowerRight',
                     customsymbol=[True,True], symbolshape=['diamond','circle'], 
                     symbolsize=[5,5], symbolcolor=['ff0000','00ff00'], 
                     symbolfill=['mesh3','mesh3'], 
                     plotfile=self.plotfile_jpg, expformat='jpg',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 200000)
        self.removePlotfile()
        # test ms has no corrected data, should revert to data
        # Plot amp:corrected(->data) vs time and amp:model vs time. 
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg',
                     overwrite=True, showgui=False, yaxis=['amp','amp'], 
                     ydatacolumn=['corrected','model'],
                     yaxislocation=['left','right'], xaxis='freq', 
                     showlegend=True, legendposition='lowerRight',
                     customsymbol=[True,True], symbolshape=['diamond','circle'], 
                     symbolsize=[5,5], symbolcolor=['ff0000','00ff00'], 
                     symbolfill=['mesh3','mesh3'], highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 120000)
        self.removePlotfile()
        # test two y-axes with default x-axis
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg',
                     overwrite=True, showgui=False, yaxis=['amp','amp'], 
                     ydatacolumn=['corrected','model'],
                     yaxislocation=['left','right']) 
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        print

    def test_axis_radec_params(self,debug=False):
        '''test_axis_radec_params: Test ant-ra/ant-dec parameters'''
        yx_axes = [('ant-ra','time'),
                   ('ant-dec','time'),
                   ('ant-dec','ant-ra')
                  ]
        interp_methods = ['nearest','cubic spline']
        # Sub-plots grid
        grid_cols = len(yx_axes)
        grid_rows = len(interp_methods)
        n_plots = grid_rows*grid_cols
        # Create 1 plot file per supported reference frame
        ref_frames = {'icrs':30000,'j2000':30000,'b1950':30000,'galactic':30000,'azelgeo':30000}
        for ref_frame, plot_min_size in ref_frames.iteritems():
            # Plot file
            plot_filename = "testAxis13_radec_" + ref_frame + ".png"
            plot_dir = self.outputDir if not debug else '/tmp'
            plot_path = os.path.join(plot_dir, plot_filename)
            self.removePlotfile(plot_path)
            # Create sub-plots, export plot when plotting last sub-plot
            plot_index = 0
            for row, interp_method in enumerate(interp_methods):
                for col, (y_axis,x_axis) in enumerate(yx_axes):
                    is_first_plot = (plot_index == 0)
                    is_last_plot = ( (plot_index + 1) == n_plots)
                    title_fmt = 'ref={ref_frame}, interp={interp}'
                    title = title_fmt.format(ref_frame=ref_frame,
                                             interp=interp_method)
                    res = plotms(
                              vis = self.ms3,
                              #
                              title = title,
                              titlefont = 10,
                              #
                              gridrows = grid_rows,
                              gridcols = grid_cols,
                              #
                              rowindex = row,
                              colindex = col,
                              plotindex = plot_index,
                              #
                              xaxis = x_axis,
                              xframe = ref_frame,
                              xinterp = interp_method,
                              #
                              yaxis = y_axis,
                              yframe = ref_frame,
                              yinterp = interp_method,
                              #
                              coloraxis = 'spw',
                              #
                              plotfile = plot_path if is_last_plot else '',
                              width = 1024,
                              height = 768,
                              highres = True,
                              overwrite = True,
                              #
                              showgui = False,
                              #
                              clearplots = is_first_plot
                          )
                    self.assertTrue(res)
                    plot_index = plot_index + 1
                    # Note: last plotms call is blocking : plotms waits 
                    #        until all plots are drawn before exporting
            self.checkPlotfile(plot_path, plot_min_size)
            if not debug:
                self.removePlotfile(plot_path)
        print

# ------------------------------------------------------------------------------

class test_calibration(plotms_test_base):
    ''' test plotms callib parameter '''

    def setUp(self):
        self.checkDisplay()
        self.setUpCalData()
        
    def tearDown(self):
        self.tearDownData()

    def test_calibration_callib(self):
        '''test_calibration_callib: CAS-3034, CAS-7502 callib parameter for OTF calibration'''
        if os.path.exists(calpath):
            self.plotfile_jpg = os.path.join(self.outputDir, "testCalibration01.jpg")
            self.removePlotfile()
            time.sleep(5)

            # callib is a string not a filename
            callibStr = "caltable='" + self.caltable + "' calwt=True tinterp='nearest'"
            res = plotms(vis=self.ms2, plotfile = self.plotfile_jpg,
                         ydatacolumn="corrected", xaxis="frequency",
                         showgui=False, callib=callibStr, highres=True)
            self.assertTrue(res)
            self.checkPlotfile(self.plotfile_jpg, 240000)
        else:
            print "Skipping test, no path to calibration tables"
        print 

    def test_calibration_badcallib(self):
        '''test_calibration_badcallib: callib file does not exist'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalibration02.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms2, plotfile = self.plotfile_jpg, 
                     ydatacolumn="corrected", xaxis="frequency",
                     showgui=False, callib='/tmp/nocallib.txt',
                     highres=True)
        self.assertFalse(res)
        print 

# ------------------------------------------------------------------------------
 
class test_calplots(plotms_test_base):
    ''' Test basic single plots and overplots '''

    def setUp(self):
        self.checkDisplay()
        # cal table for plotting
        self.setUpCalData()
        
    def tearDown(self):
        self.tearDownData

    def test_calplot_basic(self):
        '''test_calplot_basic: Basic plot of caltable'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalPlot01.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.caltable, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 30000)
        self.removePlotfile()
        print
 
    def test_calplot_axes(self):
        '''test_calplot_axes: Basic plot of caltable with non-default axes'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalPlot02.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.caltable, xaxis='scan',
                     plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 20000)
        self.removePlotfile()
        res = plotms(vis=self.caltable, yaxis='phase',
                     xaxis='baseline', overwrite=True,
                     plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 60000)
        self.removePlotfile()
        print

    def test_calplot_iteration(self):
        '''test_calplot_iteration: caltable with corr iteraxis'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalPlot03.jpg")
        plotfile1 = os.path.join(self.outputDir, "testCalPlot03_Poln1_2.jpg")
        self.removeFiles(self.outputDir, "testCalPlot03_")
        time.sleep(5)
        res = plotms(vis=self.caltable, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True, iteraxis='corr', exprange='all',
                     overwrite=True)   
        self.assertTrue(res)
        fileCount = self.getFilecount( self.outputDir, "testCalPlot03_" )
        self.assertEqual(fileCount,2)
        self.checkPlotfile(plotfile1, 30000)
        self.removeFiles(self.outputDir, "testCalPlot03_")
        print

    def test_calplot_selection(self):
        '''test_calplot_selection: caltable with corr selection'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalPlot04.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.caltable, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True, correlation='R', overwrite=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 30000)
        self.removePlotfile()
        print

    def test_calplot_ratioplot(self):
        '''test_calplot_ratioplot: caltable with corr selection'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testCalPlot05.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.caltable, plotfile=self.plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=True, correlation='/', overwrite=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        self.removePlotfile()
        print

# ------------------------------------------------------------------------------

class PlotmsPageHeader:
    ''' Analyze PlotMS page header from png image of PlotMS Plot 
        Assumptions: graphical area has a black frame, header background color is light
    '''
    def __init__(self,png_path,debug=False):
        self.png_path = png_path
        self.debug = debug
        self.color_img = plt.imread(png_path) # shape = (height,width,4)
        self.img_height , self.img_width , depth = self.color_img.shape
        self.height = -1
        self.height_ratio = -1.0
        self.rows = -1
        self._analyze()

    def _analyze(self):
        gray_img = self.color_img.min(axis=2)
        # Binarize
        is_light = gray_img > 0.75
        bin_img = np.zeros(gray_img.shape,dtype=gray_img.dtype)
        bin_img[is_light] = 1.0
        # Project on X (vertical) axis
        gray_xproj_bin = bin_img.min(axis=1)
        # White to black transitions
        (steps_down,) = np.where(np.diff(gray_xproj_bin) == -1.0 )
        if steps_down.size > 0 :
            self.height = steps_down[-1]
            self.rows = steps_down.size - 1
        else:
            self.height = 0
            self.rows = 0
        if self.img_height > 0:
            self.height_ratio = float(self.height) / self.img_height
        if self.debug:
            gray_xmask = np.outer(gray_xproj_bin,np.ones(self.img_width,dtype=np.float))
            path_parts = os.path.splitext(self.png_path)
            gray_xmask_path = ''.join([path_parts[0],'.xmask',path_parts[1]])
            plt.imsave(gray_xmask_path,gray_xmask,cmap=plt.cm.gray)

    def empty(self):
        return self.height <= 0

    def hasCorrectHeightRatio(self):
        return 0.10 <= self.height_ratio <= 0.30

class test_pageheader(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()

    def tearDown(self):
        if not self.debug:
            self.tearDownData()

    def checkPageHeader(self,expected_rows):
        page_header = PlotmsPageHeader(self.plotfile_png)
        # Empty or not
        if expected_rows == 0:
            err_msg = 'Page header: is unexpected, has a height of: '
            err_msg += '{:d} pixels'.format(page_header.height)
            self.assertTrue(page_header.empty(),err_msg)
            return
        else:
            err_msg = 'Page header: is missing'
            self.assertTrue(not page_header.empty(),err_msg)
        # Height ratio
        err_msg = 'Page header: has incorrect height ratio: '
        err_msg += '{:2.0f}% not in [10%,30%]'.format(100*page_header.height_ratio)
        self.assertTrue(page_header.hasCorrectHeightRatio(), err_msg)
        # Number of rows
        err_msg = 'Page header: has wrong number of rows: '
        err_msg += '{:d} rows detected, {:d} rows expected'.format(page_header.rows,expected_rows)
        self.assertTrue(page_header.rows == expected_rows, err_msg)

    def test_pageheader_none(self):
        '''test_pageheader_none: no page header if no header items'''
        self.plotfile_png = os.path.join(self.outputDir, "testPageHeader01.png")
        self.removePlotfile(self.plotfile_png)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_png, expformat='png', 
                     antenna='0&2',
                     showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_png, 15000)
        self.checkPageHeader(expected_rows=0)

    def test_pageheader_items01(self):
        '''test_pageheader_items01: filename,telescope,projid,observer'''
        self.plotfile_png = os.path.join(self.outputDir, "testPageHeader02.png")
        self.removePlotfile(self.plotfile_png)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_png, expformat='png', 
                     antenna='0&2',headeritems='filename,telescope,projid,observer',
                     showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_png, 20000)
        self.checkPageHeader(expected_rows=2)

    def test_pageheader_items02(self):
        '''test_pageheader_items02: targdir,telescope,targname,observer,ycolumn'''
        self.plotfile_png = os.path.join(self.outputDir, "testPageHeader03.png")
        self.removePlotfile(self.plotfile_png)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_png, expformat='png', 
                     antenna='0&2',headeritems='targdir,telescope,targname,observer,ycolumn',
                     showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_png, 20000)
        self.checkPageHeader(expected_rows=3)

# ------------------------------------------------------------------------------

class test_display(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()
       
    def test_display_symbol(self):
        '''test_display_symbol: Set a custom plotting symbol'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay01.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Test diamond shape
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, highres=True,
                     customsymbol=True, symbolshape='diamond', symbolsize=5,
                     symbolcolor='00ff00', symbolfill='mesh3')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        # Test pixel shape
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, highres=True, 
                     customsymbol=True, symbolshape='pixel',
                     symbolcolor='00aa00')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 84000)
        print

    def test_display_flsymbol(self):
        '''test_display_flsymbol: Set a custom plotting flaggedsymbol'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay02.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, customflaggedsymbol=True, highres=True,
                     flaggedsymbolshape='diamond', flaggedsymbolsize=5, 
                     flaggedsymbolcolor='00ff00', flaggedsymbolfill='mesh3')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 40000)
        print
        
    def test_display_flsymbol_false( self ):
        '''test_display_flsymbol_false:  CAS-7046:  customflaggedsymbol=False '''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay03.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testDisplay03_2.jpg")
        self.removePlotfile()
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
        # Set customflaggedsymbol=True
        res = plotms(vis=self.ms, plotfile = self.plotfile_jpg, showgui=False,
                     customflaggedsymbol=True,flaggedsymbolshape='diamond',
                     flaggedsymbolsize=5,flaggedsymbolcolor='00ff00',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 40000)
        # Set customflaggedsymbol=False
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, showgui=False,
                     customflaggedsymbol=False,flaggedsymbolshape='diamond',
                     flaggedsymbolsize=5,flaggedsymbolcolor='00ff00',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 40000)
        self.removePlotfile(plotfile2_jpg)
        print

    def test_display_legend(self):
        '''test_display_legend: Place a legend on a plot.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay04.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Place a legend in the upper right corner of the plot
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, showlegend=True, legendposition='upperRight',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        print

    def test_display_legend_overplot(self):
        '''test_display_legend_overplot: Test that legend works with overplots'''
        # Note when testing this, don't just check the file size, but look at the
        # plot and make sure there is a legend there.
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay05.jpg")
        self.removePlotfile()
        time.sleep(5)
        # First plot: scan vs time
        res = plotms(vis=self.ms, showgui=False, yaxis='scan', highres=True,
                     plotindex=0, showlegend=True, legendposition='lowerRight',
                     customsymbol=[True], symbolshape=['diamond'], symbolsize=[3],
                     symbolcolor=['ff0000'], symbolfill=['mesh3'])
        self.assertTrue(res)
        # Overplot: field vs time.
        res = plotms(vis=self.ms, showgui=False, yaxis='field',
                     plotindex=1, clearplots=False,
                     showlegend=True, legendposition='lowerRight',
                     customsymbol=[True], symbolshape=['circle'], symbolsize=[3],
                     symbolcolor=['00FF00'], symbolfill=['mesh3'],
                     plotfile=self.plotfile_jpg, highres=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 45000)
        print
         
    def test_display_coloraxis(self):
        '''test_display_coloraxis: Colorize by time on an elevation x amp plot.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay06.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Colorize by time.
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='elevation', yaxis='amp',
                     coloraxis='time', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)

        # Colorize by synonym, see CAS-6921.
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, highres=True,
                     xaxis='elevation', yaxis='amp',
                     coloraxis='chan')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 180000)

        # Colorize by averaged time on an elevation x amp plot.
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, highres=True,
                     xaxis='elevation', yaxis='amp',
                     coloraxis='time', averagedata=True, avgtime='3600')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 40000)
        print    

    def test_display_yaxisloc(self):
        '''test_display_yaxisloc: Overplot scan on left y-axis and field on right y-axis.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay07.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='time', yaxis=['scan','field'],
                     yaxislocation=['left','right'], highres=True,
                     showlegend=True, legendposition='exteriorTop',
                     customsymbol=[True,True], symbolshape=['diamond','circle'],
                     symbolsize=[5,5], symbolcolor=['ff0000','00ff00'], 
                     symbolfill=['mesh3','mesh3'])
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        print    
        
    def test_display_yaxisloc_iter(self):
        '''test_display_yaxisloc_iter: Overplot (scan and field) vs time and iterate over antenna.'''
        plotfile_jpg = os.path.join(self.outputDir, "testDisplay08.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testDisplay08_Antenna1@VLA:N7,2@VLA:W1,3@VLA:W2,4@VLA:E1.jpg")
        self.removePlotfile(plotfile_jpg)
        self.removePlotfile(plotfile1_jpg)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='time', yaxis=['scan','field'],
                     yaxislocation=['left','right'], highres=True,
                     gridrows=2, gridcols=2, iteraxis="antenna",
                     showlegend=True, legendposition='lowerLeft',
                     customsymbol=[True,True], symbolshape=['diamond','circle'],
                     symbolsize=[5,5], symbolcolor=['ff0000','00ff00'], 
                     symbolfill=['mesh3','mesh3'])
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 150000)
        self.removePlotfile(plotfile1_jpg)
        print       
      
    # Removing from test suite until tarball seg fault fixed
    def xtest_display_sharedaxis(self):
        '''test_display_sharedaxis: CAS-7074 sharedaxis on 2x2 grid'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay09.jpg")
        self.removeFiles(self.outputDir, "testDisplay09_" )
        # Create the plot and check that there are 19 iterations
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, exprange='all', 
                     scan='3', antenna='1~3', xaxis="time", yaxis="amp",
                     showgui=False, iteraxis='baseline',
                     gridrows=2, gridcols=2, highres=True,
                     xselfscale=True, yselfscale=True, 
                     xsharedaxis=True, ysharedaxis=True)
        self.assertTrue(res)
        fileCount = self.getFilecount( self.outputDir, "testDisplay09_" )
        self.assertEqual(fileCount,19)
        print

    def test_display_sharedaxis_scale(self):
        '''test_display_sharedaxis_scale: CAS-7074 xsharedaxis needs to be a subparameter of global scale'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay10.jpg")
        self.removePlotfile()
        time.sleep(5)
        # xsharedaxis without xselfscale!
        res = plotms(vis=self.ms, plotfile = self.plotfile_jpg, exprange='all',
                     gridrows=2, gridcols=2, iteraxis='antenna',
                     showgui=False, xsharedaxis=True, highres=True)
        self.assertFalse(res)
        print

    def test_display_labels(self):
        '''test_display_labels: test custom title and axis labels'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay11.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, highres=True,
                     title='NGC5921', xlabel='x axis', ylabel='y axis')
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        print

    def test_display_gridlines(self):
        '''test_display_gridlines: show major and minor grids'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testDisplay12.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, showmajorgrid=True, showminorgrid=True,
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 280000)
        print

# ------------------------------------------------------------------------------

class test_grid(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()

    def test_grid_location(self):
        '''test_grid_location: Set rowindex & colindex to place a plot in a particular location'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testGrid01.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Use a grid with 2 rows and 3 columns. 
        # Put the plot in the second row, second col
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, gridrows=2, gridcols=3, xaxis="time",
                     rowindex=1, colindex=1, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)
        print

    def test_grid_fill(self):
        '''test_grid_fill: Set grid and fill each location with a plot'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testGrid02.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testGrid02_2.jpg")
        self.removePlotfile()
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
        # Use a grid with 2 rows and 2 columns. Fill all the plots in the grid
        res = plotms(vis=self.ms, plotindex=0, title='Plot A',
                     showgui=False, gridrows=2, gridcols=2,
                     rowindex=0, colindex=0)
        self.assertTrue(res)
        res = plotms(vis=self.ms, plotindex=1, title='Plot B',
                     showgui=False, clearplots=False, 
                     rowindex=0, colindex=1)
        self.assertTrue(res)
        res = plotms(vis=self.ms, plotindex=2, title='Plot C',
                     showgui=False, clearplots=False, 
                     rowindex=1, colindex=0)
        self.assertTrue(res)
        res = plotms(vis=self.ms, plotindex=3, title='Plot D',
                     showgui=False, clearplots=False,
                     rowindex=1, colindex=1,
                     plotfile=self.plotfile_jpg, expformat='jpg',
                     highres=True)  
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 50000)

        # Now use a grid with 1 row and 2 cols.  Fill in all the plots in the grid.
        res = plotms(vis=self.ms, plotindex=0, title='Plot A',
                     showgui=False,  gridrows=1, gridcols=2,
                     rowindex=0, colindex=0)
        res = plotms(vis=self.ms, plotindex=1, title='Plot B',
                     showgui=False, clearplots=False,
                     rowindex=0, colindex=1,
                     plotfile=plotfile2_jpg, expformat='jpg',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 50000)
        self.removePlotfile(plotfile2_jpg)
        print

    def test_grid_fill_overplot(self):
        '''test_grid_fill_overplot: 2x2 multiplot display with single plots and overplots'''
        plotFile = os.path.join(self.outputDir, "testGrid03.jpg")
        plotFiles = [os.path.join(self.outputDir, "testGrid03_Scan1,2.jpg"),
                     os.path.join(self.outputDir, "testGrid03_Scan3,4,5,6_2.jpg"),
                     os.path.join(self.outputDir, "testGrid03_Scan7_3.jpg")]
        
        for  i in range(0, len(plotFiles)):
            self.removePlotfile(plotFiles[i])
        time.sleep(5)
        
        # First put the multiplot in, scan,field vs time
        # Plot scan vs time.
        res = plotms(vis=self.ms, gridrows=2, gridcols=2,  
                     showgui=False, yaxis='scan',
                     rowindex=0, colindex=0, plotindex=0,
                     customsymbol=[True], symbolshape=['diamond'], symbolsize=[3],
                     symbolcolor=['ff0000'], symbolfill=['mesh3'])
        self.assertTrue(res)
        print 'Added overplot 1'
        # Overplot field vs time in same plot
        res = plotms(vis=self.ms, gridrows=2, gridcols=2,
                     showgui=False, yaxis='field',
                     rowindex=0, colindex=0, plotindex=1, clearplots=False,
                     showlegend=True, legendposition='lowerRight',
                     customsymbol=[True], symbolshape=['circle'], symbolsize=[3],
                     symbolcolor=['00FF00'], symbolfill=['mesh3'])   
        self.assertTrue(res)
        print 'Added overplot 2'  
        # Add a single basic plot amp vs time in the first row, second column
        res = plotms(vis=self.ms, showgui=False, gridrows=2, gridcols=2, 
                     rowindex=0, colindex=1, plotindex=2, clearplots=False)
        self.assertTrue(res)
        print 'Added single plot'
        # Add in an iteration plot and export it.
        res = plotms(vis=self.ms, plotfile=plotFile, expformat='jpg',
                     showgui=False, iteraxis='scan', exprange='all',
                     gridrows=2, gridcols=2, highres=True,
                     rowindex=1, colindex=0, plotindex=3, clearplots=False)
        self.assertTrue(res)
        print 'Added iteration plot, now exporting'
        self.checkPlotfile(plotFiles[0], 180000)
        self.removePlotfile(plotFiles[0])
        self.checkPlotfile(plotFiles[1], 150000)
        self.removePlotfile(plotFiles[1])
        self.checkPlotfile(plotFiles[2], 30000)
        self.removePlotfile(plotFiles[2])
        print

    def test_grid_badplotindex(self):
        '''test_grid_badplotindex: Set plotindex out of range'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testGrid04.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Fail with nonzero plotindex when clearplots=True (first plot)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, expformat='jpg', 
                     showgui=False, iteraxis='baseline',
                     plotindex=1, rowindex=0, colindex=1, 
                     gridrows=2, gridcols=2, highres=True) 
        self.assertFalse(res)
        print

    def test_grid_badindex(self):
        '''test_grid_badindex: use row/col index larger than gridrow/gridcol'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testGrid05.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, exprange='all',
                     showgui=False, iteraxis='scan', highres=True,
                     gridrows=2, gridcols=2, rowindex=2, colindex=2)
        self.assertTrue(res)
        self.checkNoPlotfile(self.plotfile_jpg)
        print

    def xtest_grid_screenExport(self):
        '''test_grid_screenExport: Export grid plot in screen resolution'''
        # not working yet so removed this test
        plotfile_jpg = os.path.join(self.outputDir, "testGrid06.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testGrid06_Scan1,2,3,4.jpg")
        self.removePlotfile(plotfile1_jpg)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=False, iteraxis='scan',
                     gridrows=2, gridcols=2)   
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 130000)
        self.removePlotfile(plotfile1_jpg)
        # with external axes
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=False, iteraxis='scan',
                     gridrows=2, gridcols=2, xsharedaxis=True,
                     ysharedaxis=True, xselfscale=True, yselfscale=True)   
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 70000)
        self.removePlotfile(plotfile1_jpg)
        # with right yaxis
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat="jpg", 
                     showgui=False, highres=False, iteraxis='scan',
                     gridrows=2, gridcols=2, xsharedaxis=True,
                     ysharedaxis=True, xselfscale=True, yselfscale=True,
                     yaxislocation='right')   
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 70000)
        self.removePlotfile(plotfile1_jpg)
        print

# ------------------------------------------------------------------------------

class test_iteration(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()
    
    def test_iteration_scan(self):
        '''test_iteration_scan: Iterate by scan and export all'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration01.jpg")
        plotFiles = [os.path.join(self.outputDir, "testIteration01_Scan1.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan2_2.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan3_3.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan4_4.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan5_5.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan6_6.jpg"),
                     os.path.join(self.outputDir, "testIteration01_Scan7_7.jpg")]
        for i in range(0, len(plotFiles)):
            self.removePlotfile(plotFiles[i])
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat='jpg', 
                     showgui=False, iteraxis='scan', exprange='all',
                     highres=True)
        self.assertTrue(res)
        # Check each page got saved
        for  i in range(0, len(plotFiles)):
            self.checkPlotfile(plotFiles[i], 40000)
            self.removePlotfile(plotFiles[i])
        print

    def test_iteration_antenna(self):
        '''test_iteration_antenna: Iterate by antenna and export all'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration02.jpg")
        # Create the plot and check that there are 27 iterations
        res = plotms(vis=self.ms, plotfile = plotfile_jpg, exprange='all',
                     showgui=False, iteraxis='antenna', overwrite=True,
                     highres=True)
        self.assertTrue(res)
        fileCount = self.getFilecount(self.outputDir, "testIteration02_")
        # no Antenna23
        self.assertEqual(fileCount, 27)
        self.removeFiles(self.outputDir, "testIteration02_")
        print

    def test_iteration_time(self):
        '''test_iteration_time: Iterate over time'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration03.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration03_Time09:18:59.9998.jpg")
        self.removePlotfile(plotfile1_jpg)
        time.sleep(5)
        # Export first plot only 
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, showgui=False,
                     xaxis='elevation', yaxis='amp', iteraxis='time',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 30000)
        self.removePlotfile(plotfile1_jpg)
        print

    def test_iteration_timeavg(self):
        '''test_iteration_time: Iterate over averaged time'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration04.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration04_Time09:19:15.0000.jpg")
        self.removePlotfile(plotfile1_jpg)
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat='jpg', 
                     showgui=False, xaxis='elevation', yaxis='amp', avgtime='60',
                     iteraxis='time', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 30000)
        self.removePlotfile(plotfile1_jpg)
        print

    def test_iteration_grid(self):
        '''test_iteration_grid: Iterate by scan on square grid.'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration05.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration05_Scan1,2,3,4.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testIteration05_Scan5,6,7_2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg) 
        time.sleep(5)
        # Make 2 pages of 2x2 iteration plots over scan
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat='jpg',
                     showgui=False, iteraxis='scan', exprange='all',
                     gridrows=2, gridcols=2, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 80000)
        self.removePlotfile(plotfile1_jpg)
        self.checkPlotfile(plotfile2_jpg, 55000)
        self.removePlotfile(plotfile2_jpg) 
        print    
        
    def test_iteration_grid2(self):
        '''test_iteration_grid2: Iterate by scan with right axis and non-square grid.'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration06.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration06_Scan1,2,3,4,5,6.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testIteration06_Scan7_2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg) 
        time.sleep(5)
        # Make 2 pages of 3x2 iteration plots over scan
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, expformat='jpg',
                     showgui=False, iteraxis='scan', exprange='all',
                     gridrows=3, gridcols=2, highres=True,
                     yaxis=['amp'], yaxislocation=['right'])
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 170000)
        self.removePlotfile(plotfile1_jpg)
        self.checkPlotfile(plotfile2_jpg, 25000)
        self.removePlotfile(plotfile2_jpg) 
        print


    def test_iteration_selection( self ):
        '''test_iteration_selection: CAS-7050:(Pipeline) Iteration with selection'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration07.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration07_Antenna1@VLA:N7.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testIteration07_Antenna2@VLA:W1_2.jpg")
        plotfile3_jpg = os.path.join(self.outputDir, "testIteration07_Antenna3@VLA:W2_3.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        self.removePlotfile(plotfile3_jpg)
        time.sleep(5)
        # Select 3: check that there are 3 iteration plots
        res = plotms(vis=self.ms, plotfile=plotfile_jpg, exprange='all',
                     showgui=False, xaxis="time", yaxis="amp",
                     antenna='1~3&&&', iteraxis='antenna', highres=True)
        self.assertTrue(res)
        fileCount = self.getFilecount( self.outputDir, "testIteration07_" )
        self.assertEqual(fileCount, 3)
        self.checkPlotfile(plotfile1_jpg, 190000)
        self.removePlotfile(plotfile1_jpg)
        self.checkPlotfile(plotfile2_jpg, 170000)
        self.removePlotfile(plotfile2_jpg)
        self.checkPlotfile(plotfile3_jpg, 160000)
        self.removePlotfile(plotfile3_jpg)
        print
        
    def test_iteration_select1( self ):
        '''test_iteration_select1: CAS-7050 (Pipeline) Iteration with selection of 1'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration08.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration08_Antenna28@VLA:W7.jpg")
        self.removePlotfile(plotfile1_jpg)
        time.sleep(5)
        # One valid selection: check that there is only 1 iteration plot
        res = plotms(vis=self.ms, plotfile = plotfile_jpg, exprange='all',
                     xaxis="time", yaxis="amp", antenna='28~31',
                     showgui=False, iteraxis='antenna', highres=True)
        self.assertTrue(res)
        fileCount = self.getFilecount(self.outputDir, "testIteration08_" )
        self.assertEqual(fileCount, 1) 
        self.checkPlotfile(plotfile1_jpg, 60000)
        self.removePlotfile(plotfile1_jpg)
        print
       
    def test_iteration_select0( self ):
        '''test_iteration_select0: CAS-7050 (Pipeline) Iteration with empty selection'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testIteration09.jpg")
        self.removePlotfile()
        time.sleep(5)
        # Create plot with empty antenna selection and check the result is false'''
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, exprange='all',
                     xaxis="time", yaxis="amp", antenna="100,101,102",
                     showgui=False, iteraxis='antenna', highres=True)
        self.assertFalse(res)
        print

    def test_iteration_badselection( self ):
        '''test_iteration_badselection:  CAS-7050:  (Pipeline) Iteration with bad selection, you should be able to tell what plots are skipped.'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration10.jpg")
        plotFiles = [os.path.join(self.outputDir, "testIteration10_Scan1.jpg"),
                     os.path.join(self.outputDir, "testIteration10_Scan2_2.jpg"),
                     os.path.join(self.outputDir, "testIteration10_Scan3_3.jpg"),
                     os.path.join(self.outputDir, "testIteration10_Scan4_4.jpg")]
        for i in range(0, len(plotFiles)):
            self.removePlotfile(plotFiles[i])
        time.sleep(5)
        '''The expectation is the plot with the bad scan, 100, will be skipped'''
        res = plotms(vis=self.ms, plotfile = plotfile_jpg, exprange='all',
                     showgui=False, scan="1,2,100,3,4", iteraxis='scan',
                     highres=True)
        self.assertTrue(res)
        fileCount = self.getFilecount(self.outputDir, "testIteration10_" )
        self.assertEqual(fileCount, 4) 
        for  i in range(0, len(plotFiles)):
            self.checkPlotfile(plotFiles[i], 40000)
            self.removePlotfile(plotFiles[i])
        print

    def test_iteration_corr(self):
        '''test_iteration_corr: Iterate by correlation and export all'''
        plotfile_jpg = os.path.join(self.outputDir, "testIteration11.jpg")
        plotfile1_jpg = os.path.join(self.outputDir, "testIteration11_CorrRR.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testIteration11_CorrLL_2.jpg")
        # Create the plot and check that there are 27 iterations
        res = plotms(vis=self.ms, plotfile = plotfile_jpg, exprange='all',
                     showgui=False, iteraxis='corr', overwrite=True,
                     highres=True)
        self.assertTrue(res)
        fileCount = self.getFilecount(self.outputDir, "testIteration11_")
        self.assertEqual(fileCount, 2)
        self.checkPlotfile(plotfile1_jpg, 230000)
        self.removePlotfile(plotfile1_jpg)
        self.checkPlotfile(plotfile2_jpg, 220000)
        self.removePlotfile(plotfile2_jpg)
        print

# ------------------------------------------------------------------------------

class test_multi(plotms_test_base):
    ''' tests of multiple plotms arguments '''

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()

    def test_multi_cookbook(self):
        '''test_cookbook: Juergen's cookbook Plotting Multiple Data Sets example'''
        plotfile1_jpg = os.path.join(self.outputDir, "testMulti01.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testMulti01_2.jpg")
        plotfile3_jpg = os.path.join(self.outputDir, "testMulti01_3.jpg")
        plotfile4_jpg = os.path.join(self.outputDir, "testMulti01_4.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        self.removePlotfile(plotfile3_jpg)
        self.removePlotfile(plotfile4_jpg)
        time.sleep(5)

        # Plot in the second column, first row, plotindex=0
        print 'Test plot 1'
        res = plotms(vis=self.ms, gridrows=2, gridcols=2,
                     rowindex=0, colindex=1, highres=True,
                     showgui=False, plotfile=plotfile1_jpg,
                     customsymbol=True, symbolshape='diamond', symbolsize=5,
                     symbolcolor='ff0000')
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 50000)
        self.removePlotfile(plotfile1_jpg)
        # Overplot in the same panel, plotindex=1
        print 'Test plot 2'
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg,
                     showgui=False, clearplots=False, 
                     plotindex=1, rowindex=0, colindex=1,
                     gridrows=2, gridcols=2, yaxislocation='right', 
                     customsymbol=True, symbolshape='circle', symbolsize=5,
                     symbolcolor='00ff00', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 50000)
        self.removePlotfile(plotfile2_jpg)
        # Define a second plot plotindex=2, in the lower right corner
        print 'Test plot 3'
        res = plotms(vis=self.ms, clearplots=False, 
                     plotindex=2, rowindex=1, colindex=1, 
                     gridrows=2, gridcols=2,
                     showgui=False, plotfile=plotfile3_jpg,
                     customsymbol=False, yaxislocation='', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile3_jpg, 50000)
        self.removePlotfile(plotfile3_jpg)
        
        print 'Test plot 4'
        # Move the plot with the overplot one panel to the left. 
        # Clear the plots and rerun the script.
        res = plotms(vis=self.ms, showgui=False,
                     gridrows=2, gridcols=2,
                     rowindex=0, colindex=0,
                     customsymbol=True, symbolshape='diamond', symbolsize=5,
                     symbolcolor='ff0000')
        res = plotms(vis=self.ms, showgui=False, clearplots=False,
                     gridrows=2, gridcols=2,
                     plotindex=1, rowindex=0, colindex=0,
                     yaxislocation='right', 
                     customsymbol=True, symbolshape='circle', symbolsize=5,
                     symbolcolor='00ff00')
        res = plotms(vis=self.ms, showgui=False, clearplots=False,
                     gridrows=2, gridcols=2,
                     plotindex=2, rowindex=1, colindex=1, 
                     plotfile=plotfile4_jpg, highres=True,
                     customsymbol=False, yaxislocation='')
        self.assertTrue(res)
        self.checkPlotfile(plotfile4_jpg, 50000)
        self.removePlotfile(plotfile4_jpg)
        print

    def test_multi_args(self):
        '''test_multi_args: CAS-6662 Pipeline no plot scenario with multiple arguments.'''
        # Note: cannot tell what the issue/fix was from this ticket
        # Seems to be part GUI and part cache issues.

        plotfile1_jpg = os.path.join(self.outputDir, "testMulti02.jpg")
        plotfile2_jpg = os.path.join(self.outputDir, "testMulti02_2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
    
        res = plotms(vis=self.ms, showgui=False, plotfile=plotfile1_jpg,
                     xaxis='uvdist', yaxis='amp',ydatacolumn='model',
                     spw='0', scan='2,4,6,8', coloraxis='spw',
                     highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile1_jpg, 20000)
        self.removePlotfile(plotfile1_jpg)
        
        res = plotms(vis=self.ms, ydatacolumn="corrected", field="1", 
                     scan="2,3", correlation="LL,RR", coloraxis="antenna2",
                     avgtime="1e8", avgscan=True, veldef="RADIO",
                     customsymbol=True, symbolshape="autoscaling", symbolsize=2,
                     symbolcolor="0000ff",symbolfill="fill",symboloutline=False,
                     customflaggedsymbol=False, flaggedsymboloutline=False,
                     flaggedsymbolshape="nosymbol", flaggedsymbolsize=2,
                     flaggedsymbolcolor="ff0000", flaggedsymbolfill="fill", 
                     plotfile=plotfile2_jpg, showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 30000)
        self.removePlotfile(plotfile2_jpg)
        print

# ------------------------------------------------------------------------------

class test_overlays(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        self.setUpOverlayData()
        
    def tearDown(self):
        self.tearDownData()

    def test_atm_overlays(self):
        '''test_atm_overlays: showatm and showtsky overplots'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testOverlay01.jpg")
        self.removePlotfile()
        time.sleep(5)
        # basic plot with showatm, xaxis chan
        res = plotms(vis=self.ms, xaxis='chan', plotfile=self.plotfile_jpg,
            expformat="jpg", showgui=False, highres=True, showatm=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 80000)
        self.removePlotfile()
        # basic plot with showtsky, xaxis freq
        self.plotfile_jpg = os.path.join(self.outputDir, "testOverlay02.jpg")
        self.removePlotfile()
        res = plotms(vis=self.ms, xaxis='freq', plotfile=self.plotfile_jpg,
            expformat="jpg", showgui=False, highres=True, showtsky=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 80000)
        self.removePlotfile()
        # plotfile without overlay: xaxis must be chan or freq
        # so ignores showatm/tsky
        self.plotfile_jpg = os.path.join(self.outputDir, "testOverlay03.jpg")
        self.removePlotfile()
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg,
            expformat="jpg", showgui=False, highres=True, showatm=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        self.removePlotfile()
        print

    def test_image_overlay(self):
        '''test_image_overlay: atm and image overplots'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testOverlay04.jpg")
        self.removePlotfile()
        time.sleep(5)

        # showimage fails but should still plot showatm
        res = plotms(vis=self.ms, xaxis='chan', plotfile=self.plotfile_jpg,
            expformat="jpg", showgui=False, highres=True, showatm=True, showimage=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 80000)
        self.removePlotfile()

        # showimage succeeds
        self.plotfile_jpg = os.path.join(self.outputDir, "testOverlay05.jpg")
        self.removePlotfile()
        res = plotms(vis=self.ms4, xaxis='freq', spw="58", plotfile=self.plotfile_jpg,
            expformat="jpg", showgui=False, highres=True, showatm=True, showimage=True)   
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 190000)
        self.removePlotfile()
        print

# ------------------------------------------------------------------------------

class test_selection(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()
 
    def test_selection_scan(self):
        '''test_selection_scan: Check scan invalid/valid selections'''
        # Will check max size to ensure selection was done
        plotfile1_jpg = os.path.join(self.outputDir, "testSelection01_scan1.jpg") #FAIL
        plotfile2_jpg = os.path.join(self.outputDir, "testSelection01_scan2.jpg")
        plotfile3_jpg = os.path.join(self.outputDir, "testSelection01_scan3.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        self.removePlotfile(plotfile3_jpg)
        time.sleep(5)
        # Fail with invalid scan
        res = plotms(vis=self.ms, plotfile=plotfile1_jpg, expformat='jpg', 
                     showgui=False, scan='8', highres=True)
        self.assertFalse(res)
        # Succeed with valid scans
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, scan='2,4', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 60000, 130000)
        self.removePlotfile(plotfile2_jpg)
        # Succeed with different scan selection (CAS-6813)
        res = plotms(vis=self.ms, plotfile=plotfile3_jpg, expformat='jpg', 
                     overwrite=True, showgui=False, scan='5,7', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile3_jpg, 60000, 130000)
        self.removePlotfile(plotfile3_jpg)
        print

    def test_selection_spw(self):
        '''test_selection_spw: Check spw invalid/valid selections'''
        # Will check max size to ensure selection was done
        plotfile1_jpg = os.path.join(self.outputDir, "testSelection02_spw1.jpg") #FAIL
        plotfile2_jpg = os.path.join(self.outputDir, "testSelection02_spw2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
        # Fail with invalid spw
        res = plotms(vis=self.ms, plotfile=plotfile1_jpg, expformat='jpg', 
                     showgui=False, spw='500', highres=True)
        self.assertFalse(res)
        # Succeed with valid spw
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, expformat='jpg', 
                     showgui=False, spw='0', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 190000, 310000)
        self.removePlotfile(plotfile2_jpg)
        print
        
    def test_selection_ant(self):
        # Will check max size to ensure selection was done
        '''test_selection_ant: Check antenna invalid/valid selections'''
        plotfile1_jpg = os.path.join(self.outputDir, "testSelection03_ant1.jpg") #FAIL
        plotfile2_jpg = os.path.join(self.outputDir, "testSelection03_ant2.jpg")
        plotfile3_jpg = os.path.join(self.outputDir, "testSelection03_ant3.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        self.removePlotfile(plotfile3_jpg)
        time.sleep(5)
        # Fail with invalid antenna
        res = plotms(vis=self.ms, plotfile=plotfile1_jpg, expformat='jpg', 
                     showgui=False, selectdata=True, antenna='ea22&&*',
                     highres=True)  
        self.assertFalse(res)
        # Succeed without antenna selection (make sure it cleared)
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, expformat='jpg', 
                     showgui=False, highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 190000, 310000)
        self.removePlotfile(plotfile2_jpg)
        # Succeed with valid antenna 
        res = plotms(vis=self.ms, plotfile=plotfile3_jpg, expformat='jpg', 
                     showgui=False, antenna='0~1', highres=True) 
        self.assertTrue(res)
        self.checkPlotfile(plotfile3_jpg, 50000, 110000)
        self.removePlotfile(plotfile3_jpg)
        print

    def test_selection_field(self):
        '''test_selection_field: Check field invalid/valid selections'''
        # Will check max size to ensure selection was done
        plotfile1_jpg = os.path.join(self.outputDir, "testSelection04_field1.jpg") #FAIL
        plotfile2_jpg = os.path.join(self.outputDir, "testSelection04_field2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
        # Fail with invalid field
        res = plotms(vis=self.ms, plotfile=plotfile1_jpg, expformat='jpg', 
                     showgui=False, field='3', highres=True)
        self.assertFalse(res)
        # Succeed with valid field
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, expformat='jpg', 
                     showgui=False, field='1', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 70000, 160000)
        self.removePlotfile(plotfile2_jpg)
        print
       
    def test_selection_corr(self):
        '''test_selection_corr: Check corr invalid/valid selections'''
        # Will check max size to ensure selection was done
        plotfile1_jpg = os.path.join(self.outputDir, "testSelection04_corr1.jpg") #FAIL
        plotfile2_jpg = os.path.join(self.outputDir, "testSelection04_corr2.jpg")
        self.removePlotfile(plotfile1_jpg)
        self.removePlotfile(plotfile2_jpg)
        time.sleep(5)
        # Fail with invalid corr
        res = plotms(vis=self.ms, plotfile=plotfile1_jpg, expformat='jpg', 
                     showgui=False, correlation='XX', highres=True)
        self.assertFalse(res)
        # Succeed with valid corr
        res = plotms(vis=self.ms, plotfile=plotfile2_jpg, expformat='jpg', 
                     showgui=False, correlation='RR', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(plotfile2_jpg, 190000, 310000)
        self.removePlotfile(plotfile2_jpg)
        print

    # Test MS has no STATE table for intent selection

# ------------------------------------------------------------------------------

class test_transform(plotms_test_base):

    def setUp(self):
        self.checkDisplay()
        self.setUpData()
        
    def tearDown(self):
        self.tearDownData()
        
    def test_transform_freqframe(self):
        '''test_transform_freqframe: Test frequency frames'''
        frames = ['LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO',
                  'GALACTO', 'LGROUP', 'CMB'] 
        for frame in frames:
            filename = "testTransform01_" + frame + ".jpg"
            plotfile = os.path.join(self.outputDir, filename)
            self.removePlotfile(plotfile)
            res = plotms(vis=self.ms, plotfile=plotfile, yaxis='freq', 
                         showgui=False, freqframe=frame, highres=True)
            self.assertTrue(res)
            self.checkPlotfile(plotfile, 180000)
            self.removePlotfile(plotfile)
        print

    def test_transform_badframe(self):
        '''test_transform_badframe: Test that invalid freqframe fails.'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testTransform02.jpg")
        self.removePlotfile()
        time.sleep(5)
        try:
            res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, yaxis='freq', 
                         showgui=False, freqframe='J2000', highres=True)
            self.assertFalse(res)
        except RuntimeError:  # "Parameter verification failed" rethrown on test server
            pass
        print
        
    def test_transform_veldef(self):
        '''test_transform_veldef: Test velocity definitions'''
        vels = ['RADIO', 'OPTICAL', 'TRUE']
        for vel in vels:
            filename = "testTransform03_" + vel + ".jpg"
            plotfile = os.path.join(self.outputDir, filename)
            self.removePlotfile(plotfile)
            res = plotms(vis=self.ms, plotfile=plotfile, yaxis='freq', 
                         showgui=False, veldef=vel, highres=True)
            self.assertTrue(res)
            self.checkPlotfile(plotfile, 260000)
            self.removePlotfile(plotfile)
        print

    def test_transform_restfreq(self):
        '''test_transform_restfreq: Test rest frequency'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testTransform04.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, yaxis='freq', 
                     showgui=False, restfreq='1420', highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 260000)
        print
        
    def test_transform_shift(self):
        '''test_transform_shift: Test phase shift'''
        self.plotfile_jpg = os.path.join(self.outputDir, "testTransform05.jpg")
        self.removePlotfile()
        time.sleep(5)
        res = plotms(vis=self.ms, plotfile=self.plotfile_jpg, yaxis='phase', 
                     showgui=False, shift=[-15, -15], highres=True)
        self.assertTrue(res)
        self.checkPlotfile(self.plotfile_jpg, 90000)
        print

# ------------------------------------------------------------------------------
 
def suite():
    print 'Tests may fail due to DBUS timeout if the version of Qt is not at least 4.8.5'
    print
    return [test_basic,
            test_averaging,
            test_axis,
            test_calibration,
            test_calplots,
            test_pageheader,
            test_display,
            test_grid,
            test_iteration,
            test_multi,
            test_overlays,
            test_selection,
            test_transform
           ]
 
