import os
import string
import sys
import shutil
import unittest
from __main__ import default
from tasks import *
#from taskinit import *
from __casac__ import tableplot

'''
Unit tests for task plotants. It tests the following parameters:
    vis:           wrong and correct values
    figfile:       if output is created
'''
tp = tableplot.tableplot()
class plotants_test(unittest.TestCase):
    # Input and output names
    msfile = 'ic2233_1.ms'
    res = None
    fig = 'plotantstest.png'
    #tp = tableplot.tableplot()

    def setUp(self):
        self.res = None
        default(plotants)
        # Switch off the displaying of the GUI
        tp.setgui(gui=False)

        # It is not necessary to copy it for all tests
        if (not os.path.exists(self.msfile)):
            datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/ic2233/'
            shutil.copytree(datapath+self.msfile, self.msfile)

    def tearDown(self):
        if (os.path.exists(self.msfile)):
            os.system('rm -rf ' + self.msfile)

        os.system('rm -rf ' + self.fig)
        # Switch GUI back on
        tp.setgui(gui=True)


    def test1(self):
       '''Test 1: Default parameters'''
       self.res = plotants()
       self.assertFalse(self.res)

    def test2(self):
        '''Test 2: Bad input file'''
        msfile = 'badfile'
        self.res = plotants(vis=msfile)
        self.assertFalse(self.res)

    def test3(self):
        '''Test 3: Good input file and output exists'''
        self.res = plotants(vis=self.msfile, figfile=self.fig)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test4(self):
        '''Test 4: Label antenna IDs'''
        self.res = plotants(vis=self.msfile, figfile=self.fig, antindex=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test5(self):
        '''Test 5: Logarithmic antenna positions'''
        self.res = plotants(vis=self.msfile, figfile=self.fig, logpos=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test6(self):
        '''Test 6: Exclude antenna positions'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            exclude='1,5,19,14,10,13')
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test7(self):
        '''Test 7: checkbaselines'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            checkbaselines=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test8(self):
        '''Test 8: exclude checkbaselines'''
        # antenna (name) 11 is already excluded by checkbaselines
        # (warning)
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            exclude='11', checkbaselines=True)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test9(self):
        '''Test 9: Title'''
        self.res = plotants(vis=self.msfile, figfile=self.fig,
            title='IC2233')
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test10(self):
        '''Test 10: All arguments'''
        self.res = plotants(self.msfile, self.fig, True, True, '1,3,5,7,9',
            True, "IC2233")
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

def suite():
    return [plotants_test]

