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


    def test0(self):
       '''Test 0: Default parameters'''
       self.res = plotants()
       self.assertFalse(self.res)

    def test1(self):
        '''Test 1: Bad input file'''
        msfile = 'badfile'
        self.res = plotants(vis=msfile)
        self.assertFalse(self.res)

    def test2(self):
        '''Test 2: Good input file and output exists'''
        self.res = plotants(vis=self.msfile, figfile=self.fig)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

    def test3(self):
        '''Test 3: Logarithmic antenna positions'''
        self.res = plotants(vis=self.msfile, logpos=True, figfile=self.fig)
        self.assertEqual(self.res,None)
        self.assertTrue(os.path.exists(self.fig))

def suite():
    return [plotants_test]

