import os
import string
import sys
import shutil
import unittest
from __main__ import default
from tasks import *
from taskinit import *

'''
Unit tests for task plotweather. It tests the following parameters:
    vis:                wrong and correct values
    seasonal_weight:    default (0.5) and other values
    doPlot:             default (True) and False
    plotName:           if output is created; test formats

    return value:       [opacity] (type='list')
'''

datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/listobs/'

# Read the data sets from another directory
if 'TEST_DATADIR' in os.environ:  
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/listobs/'
    if os.path.isdir(DATADIR):
        datapath = DATADIR
    else:
        print('WARN: directory '+DATADIR+' does not exist')

print('plotweather tests will use data from '+datapath)

class plotweather_test(unittest.TestCase):

    # Input MS, must have WEATHER table
    msfile = 'nep2-shrunk.ms'
    msNoWeatherfile = 'ngc5921_ut.ms'
    # output plots
    fig = '/tmp/plotweathertest.png'
    defaultFig = msfile + ".plotweather.png"

    def setUp(self):
        default(plotweather)
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)
        shutil.copytree(datapath+self.msfile, self.msfile)
    
    def tearDown(self):
        if (os.path.exists(self.msfile)):
            shutil.rmtree(self.msfile)
        if (os.path.exists(self.fig)):
            os.remove(self.fig)
        if (os.path.exists(self.defaultFig)):
            os.remove(self.defaultFig)
        
    def test0(self):
        '''Test 0: Default parameters'''
        opac = plotweather()
        self.assertIsNone(opac)
       
    def test1(self):
        '''Test 1: Bad input file'''
        badmsfile = 'badfile.ms'
        opac = plotweather(vis=badmsfile)
        self.assertIsNone(opac)
        
    def test2(self):
        '''Test 2: ms with no weather, no plot '''
        if (os.path.exists(self.msNoWeatherfile)):
            shutil.rmtree(self.msNoWeatherfile)
        shutil.copytree(datapath+self.msNoWeatherfile, self.msNoWeatherfile)

        opac = plotweather(vis=self.msNoWeatherfile, plotName=self.fig)
        self.assertIsNotNone(opac)
        self.assertAlmostEqual(opac[0], 0.0054234724819465846)
        self.assertFalse(os.path.exists(self.fig))
        if (os.path.exists(self.msNoWeatherfile)):
            shutil.rmtree(self.msNoWeatherfile)
        
    def test3(self):
        '''Test 3: Good input file and output exists'''
        res = plotweather(vis=self.msfile, plotName=self.fig)
        self.assertIsNotNone(res)
        opac = res[0]/1e55
        self.assertAlmostEqual(opac, 1.3867727940788754)
        self.assertTrue(os.path.exists(self.fig))

    def test4(self):
        '''Test 4: Good input file and no output plot exists'''
        res = plotweather(vis=self.msfile, doPlot=False)
        self.assertIsNotNone(res)
        opac = res[0]/1e55
        self.assertAlmostEqual(opac, 1.3867727940788754)
        defaultFig = self.msfile + ".plotweather.png"
        self.assertFalse(os.path.exists(defaultFig))

    def test5(self):
        '''Test 5: seasonal_weight'''
        res = plotweather(vis=self.msfile, seasonal_weight=0.75, plotName=self.fig)
        self.assertIsNotNone(res)
        opac = res[0]/1e54
        self.assertAlmostEqual(opac, 6.9338639703943761)
        self.assertTrue(os.path.exists(self.fig))

    def test6(self):
        '''Test 6: pdf output format'''
        plot = '/tmp/plotweathertest.pdf'
        opac = plotweather(vis=self.msfile, plotName=plot)
        self.assertTrue(os.path.exists(plot))
        os.remove(plot)

    def test7(self):
        '''Test 7: ps output format'''
        plot = '/tmp/plotweathertest.ps'
        opac = plotweather(vis=self.msfile, plotName=plot)
        self.assertTrue(os.path.exists(plot))
        os.remove(plot)

    def test8(self):
        '''Test 8: eps output format'''
        plot = '/tmp/plotweathertest.eps'
        opac = plotweather(vis=self.msfile, plotName=plot)
        self.assertTrue(os.path.exists(plot))
        os.remove(plot)
        
    def test9(self):
        '''Test 9: svg output format'''
        plot = '/tmp/plotweathertest.svg'
        opac = plotweather(vis=self.msfile, plotName=plot)
        self.assertTrue(os.path.exists(plot))
        os.remove(plot)
        
def suite():
    return [plotweather_test]

        
        
        
        
        

