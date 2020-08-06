import os
import shutil
import testhelper as th
import numpy as np
from __main__ import default
from tasks import predictcomp 
from taskinit import *
import unittest
import exceptions


''' Python unit tests for the predictcomp task

 - tests the following parameters:
 objname: error for unsupported object vs supported object
          non-visible case
 standard: wrong standard vs correct standard
 minfreq/maxfreq: wrong unit vs correct unit
 output: check for the cl file
 antennalist: use of the configuration file to plot 'observed' visibility
              amplitudes vs uvdist. GUI is turned off.
 

'''

datapath = os.environ.get('CASAPATH').split()[0] +\
                            '/data/alma/simmos/'


class predictcomp_test(unittest.TestCase):

    def setUp(self):
        self.res=None
        default(predictcomp) 

    def tearDown(self):
        #pass
        os.system('rm -rf *.cl')

        
    def test_default(self):
        '''predictcomp: test defaults'''
        self.res=predictcomp() 
        self.assertIsNone(self.res)
 
    def test_invalid_objname(self): 
        '''predictcomp: invalid objname'''
        self.res=predictcomp(objname='Moon', minfreq='100GHz',maxfreq='120GHz') 
        self.assertIsNone(self.res)
        
    def test_valid_objname(self):
        '''predictcomp: valid objname'''
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                             standard='Butler-JPL-Horizons 2012') 
        print("type(self.res) = ",type(self.res))
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
             
    def test_invalid_freqrange(self):
        '''predictcomp: invalud freqrange'''
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100',maxfreq='120')
        self.assertIsNone(self.res)

    def test_predicted_visplot(self):
        '''predictcomp: generate visibility plot for a given array configuration''' 
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',maxfreq='120GHz', 
                 standard='Butler-JPL-Horizons 2012',antennalist=datapath+'alma.cycle5.1.cfg', showplot=False,savefig='visplot.png') 
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
        self.assertTrue(os.path.exists('visplot.png'))
          
    def test_valid_but_not_visible_objname(self):
        '''predictcomp: valid but not visible objname'''
        self.res=predictcomp(objname='Mars', epoch='2018/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                             antennalist=datapath+'alma.cycle5.1.cfg',showplot=False) 
        self.assertIsNone(self.res)

def suite():
    return [predictcomp_test]

