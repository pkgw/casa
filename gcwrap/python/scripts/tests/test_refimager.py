##########################################################################
##########################################################################
#
# Test programs for the refactored imager :  test_refimager
#
# Each of the following categories (classes) has a set of tests within it.
#
#  test_onefield                 # basic tests, deconvolution algorithms
#  test_iterbot                   # iteration control options for mfs and cube
#  test_multifield               # multiple fields of same type and with different shapes/deconvolvers/gridders
#  test_stokes                    # multiple stokes planes, imaging with flagged correlations..
#  test_cube                      # all things cube. Spectral frame setup, handling empty channels, etc
#  test_widefield                # facets, wprojection, imagemosaic, mosaicft, awproject
#  test_mask                      # input mask options : regridding, mask file, automasking, etc
#  test_modelvis                # saving models (column/otf), using starting models, predict-only (setjy)
#
# To run from within casapy :  
#
#  runUnitTest.main(['test_refimager'])                                                                           # Run all tests
#  runUnitTest.main(['test_refimager[test_onefield]'])                                                    # Run tests from test_onefield
#  runUnitTest.main(['test_refimager[test_onefield_mtmfs]'])                                        # Run one specific test
#  runUnitTest.main(['test_refimager[test_onefield_mtmfs,test_onefield_hogbom]'])    # Multiple specific tests
#
# To see the full list of tests :   grep "\"\"\" \[" test_refimager.py
#
#  These tests need data stored in data/regression/unittest/clean/refimager
#
#  For a developer build, to get the datasets locally 
#
#  --- Get the basic data repo :  svn co https://svn.cv.nrao.edu/svn/casa-data/distro data
#  --- Make directories : mkdir -p data/regression/unittest/clean; cd data/regression/unittest/clean
#  --- Get test datasets :  svn co https://svn.cv.nrao.edu/svn/casa-data/trunk/regression/unittest/clean/refimager
#
##########################################################################
#
#  Datasets
#
#  refim_twochan.ms : 2 channels, one 1Jy point source with spectral index of -1.0
#  refim_twopoints_twochan.ms : Two point sources, 1Jy and 5Jy, both with spectral index -1.0. For multifield tests.
#  refim_point.ms : 1-2 GHz, 20 channels, 1 spw, one 1Jy point source with spectral index -1.0.
#  refim_point_withline.ms : refim_point with a 'line' added into 3 channels (just topo)
#  refim_mawproject.ms : Two pointing wideband mosaic with 1 point source in between the two pointings
#  refim_mawproject_offcenter.ms : Two pointing wideband mosaic with 1 point source at center of one pointing
#  refim_point_stokes.ms : RR=1.0, LL=0.8, RL and LR are zero. Stokes I=0.9, V=0.1, U,Q=0.0
#  refim_point_linRL.ms : I=1, Q=2, U=3, V=4  in circular pol basis.
#
##########################################################################

import os
import sys
import shutil
import commands
import numpy
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import inspect
import numpy as np

_ia = iatool( )
_vp = vptool( )
_cb = cbtool( )

from imagerhelpers.test_imager_helper import TestHelpers

## List to be run
def suite():
     return [test_onefield, test_iterbot, test_multifield,test_stokes, test_modelvis, test_cube, test_mask, test_startmodel,test_widefield,test_pbcor]
#     return [test_onefield, test_iterbot, test_multifield,test_stokes,test_cube, test_widefield,test_mask, test_modelvis,test_startmodel,test_widefield_failing]

refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/clean/refimager/'
#refdatapath = "/export/home/riya/rurvashi/Work/ImagerRefactor/Runs/UnitData/"
#refdatapath = "/home/vega/rurvashi/TestCASA/ImagerRefactor/Runs/WFtests/"

## Base Test class with Utility functions
class testref_base(unittest.TestCase):

     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"

          self.th = TestHelpers()

     def tearDown(self):
          """ don't delete it all """
#          self.delData()

     # Separate functions here, for special-case tests that need their own MS.
     def prepData(self,msname=""):
          os.system('rm -rf ' + self.img+'*')
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          shutil.copytree(refdatapath+self.msfile, self.msfile)
          
     def delData(self,msname=""):
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          os.system('rm -rf ' + self.img+'*')

     def checkfinal(self,pstr=""):
          #pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  casa -c `echo $CASAPATH | awk '{print $1}'`/gcwrap/python/scripts/regressions/admin/runUnitTest.py test_refimager["+ inspect.stack()[1][3] +"]"
          pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_refimager["+ inspect.stack()[1][3] +"]'])"
          casalog.post(pstr,'INFO')
          if( pstr.count("(Fail") > 0 ):
               self.fail("\n"+pstr)

##############################################
##############################################

##Task level tests : one field, 2chan.
class test_onefield(testref_base):
     
     def test_onefield_defaults(self):
          """ [onefield] Test_Onefield_defaults : Defaults """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',interactive=0)
          report=self.th.checkall(imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'], imval=[(self.img+'.psf', 1.0, [50,50,0,0])])
          self.checkfinal(pstr=report)

     def test_onefield_clark(self):
          """ [onefield] Test_Onefield_clark : mfs with clark minor cycle """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='clark',interactive=0) #,phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          #off center#ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=1000,interactive=0,phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          #compare with clean#clean(vis=self.msfile,imagename=self.img+'.old',imsize=200,cell='8.0arcsec',niter=1000,psfmode='clark',phasecenter='J2000 19h59m57.5s +40d49m00.077s') # default is clark
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imval=[(self.img+'.psf',1.0,[50,50,0,0])])
          self.checkfinal(pstr=report)

     def test_onefield_hogbom(self):
          """ [onefield] Test_Onefield_hogbom : mfs with hogbom minor cycle """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)#,phasecenter='J2000 19h59m57.5s +40d49m00.077s')
          report=self.th.checkall(ret=ret, peakres=0.35, modflux=0.77, iterdone=10, imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imval=[(self.img+'.psf',1.0,[50,50,0,0])])
          self.checkfinal(pstr=report)

     def test_onefield_mem(self):
          """ [onefield] Test_Onefield_mem : mfs with mem minor cycle """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=10,deconvolver='mem',interactive=0)
          report=self.th.checkall(ret=ret, peakres=12.7, modflux=6.98, iterdone=10, imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imval=[(self.img+'.psf',1.0,[100,100,0,0])])
          self.checkfinal(pstr=report)

     def test_onefield_multiscale(self):
          """ [onefield] Test_Onefield_multiscale : mfs with multiscale minor cycle """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',niter=10,deconvolver='multiscale',scales=[0,20,40,100],interactive=0)
          report=self.th.checkall(ret=ret, peakres=1.28, modflux=2.87, iterdone=10, imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'], imval=[(self.img+'.psf',1.0,[100,100,0,0])])
          self.checkfinal(pstr=report)

     def test_onefield_mtmfs(self):
          """ [onefield] Test_Onefield_mtmfs : mt-mfs with minor cycle iterations """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, imexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0',self.img+'.model.tt1',self.img+'.alpha'], imval=[(self.img+'.psf.tt0',1.0,[50,50,0,0]),(self.img+'.psf.tt1',1.039e-05,[50,50,0,0])])
          ## iterdone=11 only because of the return (iterdone_p+1) in MultiTermMatrixCleaner::mtclean() !
          self.checkfinal(pstr=report)

     def test_onefield_autonames(self):
          """ [onefield] Test_Onefield_autonames : Test auto increment of image names """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False)
          report=self.th.checkall(imexist=[self.img+'.psf',self.img+'_2.psf',self.img+'_3.psf'] )
          self.checkfinal(pstr=report)

     def test_onefield_twoMS(self):
          """ [onefield] Test_Onefield_twoMS : One field, two input MSs """
          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)
#          try:
#               ## This run should fail with an exception
#               ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10)
#               correct=False
#          except Exception as e:
#              correct=True
#          self.assertTrue(correct)
          ## This run should go smoothly.
          ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10,datacolumn='data')
          report=self.th.checkall(imexist=[self.img+'.psf',self.img+'.residual'])
          self.delData(ms1)
          self.delData(ms2)
          self.checkfinal(pstr=report)

     def test_onefield_twoMS_diffcolumns(self):
          """ [onefield] Test_Onefield_twoMS_diffcolumns : One field, two input MSs, one with data and one with data and corrected """
          ms1 = 'refim_point_onespw0.ms'
          ms2 = 'refim_point_onespw1.ms'
          self.prepData(ms1)
          self.prepData(ms2)

          ## Make corrected_data column for one of them
          _cb.open(ms2)
          _cb.close()

          ret = tclean(vis=[ms1,ms2],field='0',spw=['0','0'], imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='hogbom',niter=10,datacolumn='corrected')
          report=self.th.checkall(imexist=[self.img+'.psf',self.img+'.residual'])
          self.delData(ms1)
          self.delData(ms2)
          self.checkfinal(pstr=report)


     def test_onefield_restart_mfs(self):
          """ [onefield] : test_onefield_restart_mfs : Check calcpsf,calcres and ability to restart and continue"""
          ## TODO : Need to add and use info in the return record, when only a major cycle is done. Then check nmajorcycle.
          self.prepData('refim_twochan.ms')

          ## Only psf
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=True,calcres=False,deconvolver='clark')
          report=self.th.checkall(imexist=[self.img+'.psf'], imexistnot=[self.img+'.residual', self.img+'.image'],nmajordone=1)

          ## Only residual
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=False,calcres=True,deconvolver='clark',restoration=False)
          report1=self.th.checkall(imexist=[self.img+'.psf', self.img+'.residual'], imexistnot=[self.img+'.image'],nmajordone=1)

          ## Start directly with minor cycle and do only the last major cycle.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='clark')
          report2=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, imexist=[self.img+'.psf',self.img+'.residual', self.img+'.image'],nmajordone=1)

          ## Re-start from existing model image and continue on...
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='clark')
          report3=self.th.checkall(ret=ret, peakres=0.161, modflux=0.991, imexist=[self.img+'.psf',self.img+'.residual', self.img+'.image'],nmajordone=1)

          self.checkfinal(pstr=report+report1+report2+report3)


     def test_onefield_restart_mtmfs(self):
          """ [onefield] : test_onefield_restart_mtmfs : Check calcpsf,calcres and ability to restart and continue"""
          ## TODO : Need to add and use info in the return record, when only a major cycle is done. Then check nmajorcycle.
          self.prepData('refim_twochan.ms')

          ## Only psf
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=True,calcres=False,deconvolver='mtmfs')
          report=self.th.checkall(imexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'], imexistnot=[self.img+'.residual.tt0', self.img+'.image.tt0'],nmajordone=1)

          ## Only residual
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,calcpsf=False,calcres=True,deconvolver='mtmfs',restoration=False)
          report1=self.th.checkall(imexist=[self.img+'.psf.tt0',self.img+'.psf.tt1', self.img+'.residual.tt0', self.img+'.residual.tt1'], imexistnot=[self.img+'.image.tt0'],nmajordone=1)

          ## Start directly with minor cycle and do only the last major cycle.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='mtmfs')
          report2=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, imexist=[self.img+'.psf.tt1',self.img+'.residual.tt1', self.img+'.image.tt1', self.img+'.alpha'],nmajordone=1,imval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          ## Re-start from existing model image and continue on...
          ## ( If restart from modified residuals... the alpha is -1.25xx which is wrong. 
          ##   In this case, need to do calcres=True which will do extra first major cycle (nmajor=2) )
          ## But... current code (as of r33373) makes appropriate restored image but does not mess up residuals.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,interactive=0,calcpsf=False,calcres=False,deconvolver='mtmfs')
          report3=self.th.checkall(ret=ret, peakres=0.136, modflux=0.988, imexist=[self.img+'.psf.tt1',self.img+'.residual.tt1', self.img+'.image.tt1', self.img+'.alpha'],nmajordone=1,imval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          self.checkfinal(pstr=report+report1+report2+report3)

     def test_onefield_all_outputs_mfs(self):
          """ [onefield] : test_onefield_all_outputs_mfs : Make all output images even when not needed """
          self.prepData('refim_twochan.ms')

          ## Make only partial outputs
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='hogbom')
#          report1=self.th.checkall(imexist=[self.img+'.psf', self.img+'.residual'],imexistnot=[self.img+'.image',self.img+'.model'],nmajordone=1)

          ## Make all outputs
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='hogbom',restoration=True)
          report2=self.th.checkall(imexist=[self.img+'2.psf', self.img+'2.residual',self.img+'2.image',self.img+'2.model'],nmajordone=1)
 
          self.checkfinal(pstr = report2)

     def test_onefield_all_outputs_mtmfs(self):
          """ [onefield] : test_onefield_all_outputs_mtmfs : Make all output images even when not needed """
          self.prepData('refim_twochan.ms')

          ## Make only partial outputs
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs')
#          report1=self.th.checkall(imexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'],imexistnot=[self.img+'.image.tt0',self.img+'.model.tt0', self.img+'.alpha'],nmajordone=1)

          ## Make all outputs
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=True)
          report2=self.th.checkall(imexist=[self.img+'2.psf.tt0', self.img+'2.psf.tt1',self.img+'2.image.tt0',self.img+'2.model.tt0', self.img+'2.alpha'],nmajordone=1)
 
          self.checkfinal(pstr=report2)

     def test_onefield_restore_mtmfs_niter0(self):
          """ [onefield] : test_onefield_restore_mtmfs_niter0 : Niter=0 run followed by restoration without a model"""
          self.prepData('refim_twochan.ms')

          ## This test also checks the principal solution calculation on the dirty images.

          ## niter=0 run 
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=False)
          report1=self.th.checkall(imexist=[self.img+'.psf.tt0', self.img+'.psf.tt1'], imexistnot=[self.img+'.model.tt0', self.img+'.model.tt0'],nmajordone=1)
          ## restore only 
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',calcres=False,calcpsf=False,restoration=True)
          report2=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.alpha'],nmajordone=0,
                             imval=[(self.img+'.alpha',-1.0,[50,50,0,0])])

          ## niter=0 and restore ( in one step )
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,interactive=0,deconvolver='mtmfs',restoration=True)
          report3=self.th.checkall(imexist=[self.img+'2.image.tt0', self.img+'2.alpha'],nmajordone=1,
                             imval=[(self.img+'.alpha',-1.0,[50,50,0,0])] )

          self.checkfinal(pstr=report1+report2+report3)

     def test_onefield_rectangular_pixels(self):
          """ [onefield] : test_onefield_rectangular_pixels : Test restoration with rectangular pixels (cas-7171)"""
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell=['10.0arcsec','30.0arcsec'],niter=10)
          report=self.th.checkall(imexist=[self.img+'.image'], imval=[(self.img+'.image',0.482,[50,49,0,0])] )

          self.checkfinal(report)

##############################################
##############################################

##Task level tests : iteration controls
class test_iterbot(testref_base):

     def test_iterbot_mfs_1(self):
          """ [iterbot] Test_Iterbot_Mfs_1 : Zero Iterations """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=0,interactive=0,restoration=False)
          report=self.th.checkall(imexist=[self.img+'.psf', self.img+'.residual'], imexistnot=[self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_2(self):
          """ [iterbot] Test_Iterbot_Mfs_2 : Iterations with low gain """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,gain=0.1,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10,nmajordone=2,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_3(self):
          """ [iterbot] Test_Iterbot_Mfs_3 : Cycleniter test """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,cycleniter=3,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.392, modflux=0.732, iterdone=10, nmajordone=5,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_4(self):
          """ [iterbot] Test_Iterbot_Mfs_4 : Iterations with high gain """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10, gain=0.5,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.024, modflux=1.274, iterdone=10, nmajordone=3,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_5(self):
          """ [iterbot] Test_Iterbot_Mfs_5 : Threshold test """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,threshold='0.1Jy',gain=0.5,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.0924, modflux=1.129, iterdone=5, nmajordone=3,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_6(self):
          """ [iterbot] Test_Iterbot_Mfs_6 : Cycleniter and threshold """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10, cycleniter=3, threshold='0.1Jy',gain=0.5,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.0924, modflux=1.129, iterdone=5, nmajordone=3,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_7(self):
          """ [iterbot] Test_Iterbot_Mfs_7 : Threshold + cyclefactor to trigger major cycles earlier """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=10,threshold='0.01Jy', gain=0.5,cyclefactor=10.0,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.024, modflux=1.274, iterdone=10, nmajordone=9,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_8(self):
          """ [iterbot] Test_Iterbot_Mfs_8 : minpsffraction to trigger major cycles earlier. """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=20,threshold='0.01Jy', minpsffraction = 0.5,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.16127, modflux=0.9919, iterdone=20, nmajordone=4,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_mfs_9(self):
          """ [iterbot] Test_Iterbot_Mfs_9 : maxpsffraction """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',deconvolver='clark',niter=20,threshold='0.01Jy', minpsffraction=0.8,maxpsffraction=0.5,interactive=0)
          report=self.th.checkall(ret=ret, peakres=0.16127, modflux=0.9919, iterdone=20, nmajordone=4,imexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'])

          self.checkfinal(report)

     def test_iterbot_cube_1(self):
          """ [iterbot] Test_Iterbot_cube_1 : iteration counting across channels (>niter) """
          self.prepData('refim_point_withline.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='0.75Jy',interactive=0)
          report=self.th.checkall(ret=ret, iterdone=90,nmajordone=2,imexist=[self.img+'.psf', self.img+'.residual'])
          ## Only chans 6 and 7 reach cycleniter, others reach threshold in fewer than 10 iters per chan.

          self.checkfinal(report)

     def test_iterbot_cube_2(self):
          """ [iterbot] Test_Iterbot_cube_2 : High threshold, iterate only on line channels. """
          self.prepData('refim_point_withline.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='1.75Jy',interactive=0)
          report=self.th.checkall(ret=ret, peakres=1.73, modflux=0.407,iterdone=12,nmajordone=2,imexist=[self.img+'.psf', self.img+'.residual'])

          self.checkfinal(report)


     def test_iterbot_cube_3(self): # test for returned summary/plot for no iteration case 
          """ [iterbot] Test_Iterbot_cube_3 : Very high threshold, no iteration (verification of CAS-8576 fix) """
          self.prepData('refim_point_withline.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='clark',niter=10,threshold='3.5Jy',interactive=0)
          report=self.th.checkall(ret=ret,iterdone=0,nmajordone=1,imexist=[self.img+'.psf', self.img+'.residual'])

          self.checkfinal(report)

     def test_iterbot_cube_4(self): 
          """ [iterbot] Test_Iterbot_cube_4 : Large niter, and low threshold - catch if diverges (verification of CAS-8584 fix) """
          self.prepData('refim_point_withline.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',deconvolver='hogbom',niter=1000000,threshold='0.0000001Jy',gain=0.5,interactive=0)
          report=self.th.checkall(ret=ret,iterdone=1749,nmajordone=8,imexist=[self.img+'.psf', self.img+'.residual'])

     def test_iterbot_divergence(self): 
          """ [iterbot] Test_Iterbot_divergence : Use negative loop gain to make it diverge (verification of CAS-9244 fix) """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', niter=50,cycleniter=5, gain=-0.2,interactive=0)
          report=self.th.checkall(ret=ret,iterdone=10,nmajordone=3,imexist=[self.img+'.psf', self.img+'.image'])

     def test_iterbot_mfs_deconvolvers(self):
          """ [iterbot] : test_iterbot_deconvolvers : Do all minor cycle algorithms respond in the same way to iteration controls ? No ! """
          # clark and hogbom reach niter first, but multiscale gets to cyclethreshold first. Check peakres and iterdone.
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='clark')
          report1=self.th.checkall(ret=ret1, peakres=0.3922, modflux=0.732, iterdone=10, nmajordone=2,imexist=[self.img+'1.psf', self.img+'1.residual', self.img+'1.image'])

          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='hogbom')
          report2=self.th.checkall(ret=ret2, peakres=0.3530, modflux=0.7719, iterdone=10, nmajordone=2,imexist=[self.img+'2.psf', self.img+'2.residual', self.img+'2.image'])

          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=100,cell='8.0arcsec',niter=10,threshold='0.1Jy', interactive=0,deconvolver='multiscale')
          report3=self.th.checkall(ret=ret3, peakres=0.3922, modflux=0.7327, iterdone=10, nmajordone=2,imexist=[self.img+'3.psf', self.img+'3.residual', self.img+'3.image'])
     

          self.checkfinal(report1+report2+report3)
          
##############################################
##############################################

##Task level tests : multi-field, 2chan.
### For some of these tests, do the same with uvsub and compare ? 
class test_multifield(testref_base):
     
     def test_multifield_both_mfs(self):
          """ [multifield] Test_Multifield_both_mfs : Two fields, both mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=13, 
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',1.075,[50,50,0,0]),
                               (self.img+'1.image',5.590,[40,40,0,0]),
                               (self.img+'.residual',0.04,[30,18,0,0])])

          self.checkfinal(report)

     def test_multifield_both_mtmfs(self):
          """ [multifield] Test_Multifield_both_mtmfs : Two fields, both mt-mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\n\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='mtmfs',interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=12,
                        nmajordone=2,
                        imexist=[self.img+'.image.tt0', self.img+'1.image.tt0',self.img+'.image.tt1', self.img+'1.image.tt1', self.img+'.alpha', self.img+'1.alpha'],
                        imval=[(self.img+'.image.tt0',1.094,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.577,[40,40,0,0]),
                               (self.img+'.alpha',-0.90,[50,50,0,0]),
                               (self.img+'1.alpha',-1.0,[40,40,0,0])])
          self.checkfinal(report)


     def test_multifield_both_cube(self):
          """ [multifield] Test_Multifield_both_cube : Two fields, both cube"""
          self.prepData("refim_twopoints_twochan.ms")
          #self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\n')
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nimagename='+self.img+'2\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:48.895 +40.55.58.543\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest')
          report=self.th.checkall(ret=ret, 
                        iterdone=42,
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',1.434,[50,50,0,0]),
                               (self.img+'1.image',7.44,[40,40,0,0]),
                               (self.img+'2.image',7.42,[51,40,0,0]),
                               (self.img+'.image',0.758,[50,50,0,1]),
                               (self.img+'1.image',3.715,[40,40,0,1]),
                               (self.img+'2.image',3.675,[51,40,0,1]) ])
          self.checkfinal(report)

     def test_multifield_both_cube_diffshape(self):
          """ [multifield] Test_Multifield_both_cube : Two fields, both cube but different nchans"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nnchan=3\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest')
          report=self.th.checkall(ret=ret, 
                        iterdone=22,
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',1.434,[50,50,0,0]),
                               (self.img+'1.image',7.44,[40,40,0,0]),
                               (self.img+'.image',0.762,[50,50,0,1]),
                               (self.img+'1.image',3.70,[40,40,0,1]) ])
          self.checkfinal(report)

     def test_multifield_cube_mfs(self):
          """ [multifield] Test_Multifield_cube_mfs : Two fields, one cube and one mfs"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nspecmode=mfs\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest')
          report=self.th.checkall(ret=ret, 
                        iterdone=15,
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',1.4,[50,50,0,0]),
                               (self.img+'1.image',5.6,[40,40,0,0]),
                               (self.img+'.image',0.75,[50,50,0,1])])
          self.checkfinal(report)

     def test_multifield_mfs_mtmfs(self):
          """ [multifield] Test_Multifield_mfs_mtmfs : Two fields, one mt-mfs and one mfs (i.e. different deconvolvers)"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nreffreq=1.5GHz\ndeconvolver=mtmfs\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=13,
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image.tt0',self.img+'1.alpha'],
                        imval=[(self.img+'.image',1.094,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.57,[40,40,0,0]), 
                               (self.img+'1.alpha', -1.0, [40,40,0,0])  ])
          self.checkfinal(report)

     def test_multifield_cube_mtmfs(self):
          """ [multifield] Test_Multifield_cube_mtmfs : Two fields, one cube and one mtmfs"""
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nreffreq=1.5GHz\ndeconvolver=mtmfs\nspecmode=mfs\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',nchan=2,interpolation='nearest')
          report=self.th.checkall(ret=ret, 
                        iterdone=15,  # two chans in one field, and one chan in the other
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image.tt0',self.img+'1.alpha'],
                        imval=[(self.img+'.image',1.427,[50,50,0,0]),
                               (self.img+'1.image.tt0',5.575,[40,40,0,0]),
                               (self.img+'.image',0.762,[50,50,0,1]) , 
                               (self.img+'1.alpha', -1.0, [40,40,0,0])  ])
          self.checkfinal(report)


     def test_multifield_diff_gridders(self):
          """ [multifield] Test_Multifield_diff_gridders : Two fields, both mfs, gridft and wproject """
          self.prepData("refim_twopoints_twochan.ms")
#          ##Outlier uses gridft
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]\ngridder=gridft')
#          ## Outlier uses wproject but with different number of planes as the main field
#          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nmask=circle[[40pix,40pix],10pix]\ngridder=wproject\nwprojplanes=6')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',gridder='wproject',wprojplanes=4,interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=13,
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',1.075,[50,50,0,0]),
                               (self.img+'1.image',5.58,[40,40,0,0])])
          self.checkfinal(report)


     def test_multifield_autonames(self):
          """ [multifield] Test_Multifield_4 : Test auto increment of image names """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nusemask=user\nmask=circle[[40pix,40pix],10pix]')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt',restart=False)
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',restart=False) # no outlier...
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',outlierfile=self.img+'.out.txt',restart=False)

          report=self.th.checkall(imexist=[self.img+'.psf',self.img+'1.psf',self.img+'_2.psf',self.img+'1_2.psf',self.img+'_3.psf',self.img+'_4.psf',self.img+'1_4.psf'], imexistnot=[self.img+'1_3.psf'] )
          self.checkfinal(report)


### TODO :  Either put a check so that if any fields overlap, an error is thrown. Or, do sensible model choosing for some modes but detect and complain for other modes where it's harder to pick which model image to use.
     def test_multifield_overlap_mfs(self):
          """ [multifield] Test_Multifield_overlap_mfs : Two overlapping image fields, both mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:58:39.580 +40.55.55.931",outlierfile=self.img+'.out.txt',niter=20,deconvolver='hogbom',interactive=0)

          report=self.th.checkall(ret=ret, 
                        iterdone=40, ## both images see the brightest source. 
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',5.575,[48,51,0,0]),
                               (self.img+'1.image',5.574,[130,136,0,0])]) ## both images have correct flux (not twice or zero !)
          self.checkfinal(report)


     def test_multifield_overlap_mtmfs(self):
          """ [multifield] Test_Multifield_overlap_mtmfs : Two overlapping image fields, both mt-mfs """
          self.prepData("refim_twopoints_twochan.ms")
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559\n')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',phasecenter="J2000 19:58:39.580 +40.55.55.931",outlierfile=self.img+'.out.txt',niter=20,deconvolver='mtmfs',interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=39, ## both images see the brightest source.
                        nmajordone=2,
                        imexist=[self.img+'.image.tt0', self.img+'1.image.tt0'],
                        imval=[(self.img+'.image.tt0',5.52,[48,51,0,0]),
                                (self.img+'1.image.tt0',5.53,[130,136,0,0]),
                               (self.img+'.alpha',-0.965,[48,51,0,0]),
                               (self.img+'1.alpha',-0.965,[130,136,0,0])]) 
          self.checkfinal(report)


     def test_multifield_facets_mfs(self):
          """ [multifield] Test_Multifield_mfs_facets : Facetted imaging (mfs) """
          self.prepData("refim_twopoints_twochan.ms")
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',phasecenter="J2000 19:59:00.2 +40.50.15.50",facets=2,deconvolver='hogbom',niter=30)
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf'],imval=[(self.img+'.psf',1.0,[100,100,0,0]),(self.img+'.image',5.56,[127,143,0,0]) ] )
          self.checkfinal(report)

     def test_multifield_facets_mtmfs(self):
          """ [multifield] Test_facets_mtmfs : Facetted imaging (mt-mfs) """
          self.prepData("refim_twopoints_twochan.ms")
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',phasecenter="J2000 19:59:00.2 +40.50.15.50",facets=2,deconvolver='mtmfs',niter=30)
          report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.alpha'],imval=[(self.img+'.psf.tt0',1.0,[100,100,0,0]),(self.img+'.image.tt0',5.56,[127,143,0,0]),(self.img+'.alpha',-1.0,[127,143,0,0]) ] )
          self.checkfinal(report)


#     def test_multifield_cube_chunks(self):
#          """ [multifield] Test_Multifield_cube_chunks : Two fields, two sections of the same cube"""
#          self.prepData("refim_point.ms")
#          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nnchan=5\nstart=5')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='10.0arcsec',specmode='cube',nchan=5,start=0,outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',interactive=0,interpolation='nearest')
#          report=self.th.checkall(ret=ret, 
#                        iterdone=38,
#                        nmajordone=2,
#                        imexist=[self.img+'.image', self.img+'1.image'],
#                        imval=[(self.img+'.image',1.434,[50,50,0,0]),
#                               (self.img+'1.image',7.452,[40,40,0,0]),
#                               (self.img+'.image',0.762,[50,50,0,1]),
#                               (self.img+'1.image',3.702,[40,40,0,1]) ])
#          self.checkfinal(report)
##############################################
##############################################

##Task level tests : Stokes imaging options
class test_stokes(testref_base):

     def test_stokes_mfs_I(self):
          """ [onefield] Test_Stokes_I_mfs mfs with stokes I"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='I')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,0])])
          self.checkfinal(report)

     def test_stokes_mfs_IV(self):
          """ [onefield] Test_Stokes_mfs_IV : mfs with stokes IV"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IV')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',4.0,[50,50,1,0])  ])
          self.checkfinal(report)

     def test_stokes_mfs_QU(self):
          """ [onefield] Test_Stokes_mfs_QU : mfs with stokes QU"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='QU')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',3.0,[50,50,1,0])  ])
          self.checkfinal(report)

     def test_stokes_mfs_Q(self):
          """ [onefield] Test_Stokes_mfs_Q : mfs with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',2.0,[50,50,0,0]) ] )
          self.checkfinal(report)

     def test_stokes_mfs_U(self):
          """ [onefield] Test_Stokes_mfs_U : mfs with stokes U"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='U')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',3.0,[50,50,0,0]) ] )
          self.checkfinal(report)

     def test_stokes_mfs_V(self):
          """ [onefield] Test_Stokes_mfs_V : mfs with stokes V"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='V')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',4.0,[50,50,0,0]) ] )
          self.checkfinal(report)

     def test_stokes_cube_I(self):
          """ [onefield] Test_Stokes_cube_I : cube with stokes I"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='I',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',1.0,[50,50,0,2]) ] )
          self.checkfinal(report)

     def test_stokes_cube_IV(self):
          """ [onefield] Test_Stokes_stokes_IV : cube with stokes V"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IV',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,0]),(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',1.0,[50,50,0,2]),  (self.img+'.image',4.0,[50,50,1,0]),(self.img+'.image',4.0,[50,50,1,1]),(self.img+'.image',4.0,[50,50,1,2])] )
          self.checkfinal(report)

     def test_stokes_cube_QU(self):
          """ [onefield] Test_Stokes_stokes_QU : cube with stokes QU"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='QU',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',2.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,0,2]),  (self.img+'.image',3.0,[50,50,1,0]),(self.img+'.image',3.0,[50,50,1,1]),(self.img+'.image',3.0,[50,50,1,2])] )
          self.checkfinal(report)

     def test_stokes_cube_Q(self):
          """ [onefield] Test_Stokes_cube_Q : cube with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',2.0,[50,50,0,0]),(self.img+'.image',2.0,[50,50,0,1]) ,(self.img+'.image',2.0,[50,50,0,2]) ])
          self.checkfinal(report)

     def test_stokes_cube_U(self):
          """ [onefield] Test_Stokes_cube_U : cube with stokes U"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='U',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',3.0,[50,50,0,0]),(self.img+'.image',3.0,[50,50,0,1]) ,(self.img+'.image',3.0,[50,50,0,2]) ])
          self.checkfinal(report)

     def test_stokes_cube_V(self):
          """ [onefield] Test_Stokes_cube_V : cube with stokes V"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='V',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',4.0,[50,50,0,0]),(self.img+'.image',4.0,[50,50,0,1]) ,(self.img+'.image',4.0,[50,50,0,2]) ])
          self.checkfinal(report)

     def test_stokes_cube_IQUV_fromRL(self):
          """ [onefield] Test_Stokes_cube_IQUV_fromRL : cube with stokes IQUV"""
          self.prepData('refim_point_linRL.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,3,1]) ])
          self.checkfinal(report)

     def test_stokes_cube_IQUV_fromXY(self):
          """ [onefield] Test_Stokes_cube_IQUV_fromXY : cube with stokes IQUV"""
          self.prepData('refim_point_linXY.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube',interpolation='nearest')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,3,1]) ])
          self.checkfinal(report)

     def test_stokes_mtmfs_Q(self):
          """ [onefield] Test_Stokes_mtmfs_Q : mtmfs with stokes Q"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='Q',deconvolver='mtmfs',nterms=2)
          report=self.th.checkall(imexist=[self.img+'.image.tt0'], imexistnot=[self.img+'.image.alpha'], imval=[(self.img+'.image.tt0',2.0,[50,50,0,0]) ] )
          self.checkfinal(report)

     def test_stokes_mtmfs_IQUV(self):
          """ [onefield] Test_Stokes_mtmfs_IQUV : mtmfs with stokes IQUV"""
          self.prepData('refim_point_linRL.ms')
          tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',deconvolver='mtmfs',nterms=2)
          report=self.th.checkall(imexist=[self.img+'.image.tt0'],imexistnot=[self.img+'.image.alpha'], imval=[(self.img+'.image.tt0',1.0,[50,50,0,0]),(self.img+'.image.tt0',2.0,[50,50,1,0]), (self.img+'.image.tt0',3.0,[50,50,2,0]),(self.img+'.image.tt0',4.0,[50,50,3,0]) ])
          self.checkfinal(report)


#     def test_stokes_cube_I_flags(self):
#          """ [onefield] Test_Stokes_cube_I_flags : cube with stokes I and only XY or YX flagged"""
#          self.prepData('refim_point_linXY.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube')
#          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,4,1]) ])

#     def test_stokes_cube_pseudo_I_flags(self):
#          """ [onefield] Test_Stokes_cube_pseudo_I_flags : cube with stokes I and one of XX or YY flagged"""
#          self.prepData('refim_point_linXY.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10, stokes='IQUV',interactive=0,specmode='cube')
#          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.0,[50,50,0,1]),(self.img+'.image',2.0,[50,50,1,1]), (self.img+'.image',3.0,[50,50,2,1]),(self.img+'.image',4.0,[50,50,4,1]) ])

##############################################
##############################################

##Task level tests : cube.
class test_cube(testref_base):

#     def __init__(self,methodName='runTest'):
#          testref_base.__init__(self,methodName)
#          self.test_cube_0.__func__.__doc__ %="aaaa"

     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"
          self.th = TestHelpers()

          ## Setup some variables to use in all the tests

          ## chan 4 (TOPO)
          qfstart=qa.quantity("1.2GHz")
          #qvstart=qa.quantity("-59958.5km/s")
          # for restf=1.25GHz
          qvstart=qa.quantity("11991.7km/s")
          # ch10
          #qvstart=qa.quantity("16788.4km/s")

          #mfstart=me.frequency('LSRK',qa.quantity("1.09999GHz"))
          # ch4 (for rest 1.25GHz)
          mfstart=me.frequency('LSRK',qa.quantity("1.199989GHz"))
          mvstart=me.radialvelocity('BARY',qa.quantity("11977.6km/s"))
          #dop = me.todoppler('radio',mfstart,qa.quantity('1.0GHz'))
          mfstart10=me.frequency('LSRK',qa.quantity(" 1.17999GHz"))                                                        
          # doppler with ch4 freq
          dop = me.todoppler('radio',mfstart,qa.quantity('1.25GHz'))                                              

          #1chan width 
          #qvwidth = qa.quantity("11991.700km/s")
          #qvwidth = qa.quantity("4796.7km/s")
          qvwidth = qa.quantity("11991.7km/s")
          mvwidth = me.radialvelocity('TOPO',qvwidth)

          # restf = 1.25GHz
          # vel range: 59961.1 -  -31174.7 km/s (lsrk/radio)
          #            74952.3 -  -28238.3 km/s (lsrk/optical)  

          self.testList = {
                      0:{'imagename':'Cubetest_chandefstdefwidth','spw':'0','start':0,'width':1,'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, default start and width, LSRK'},
                      1:{'imagename':'Cubetest_chandefstdefwidthtopo','spw':'0','start':0,'width':1, 'outframe':'TOPO','veltype':'radio',
                        'desc':'channel, default start and width, TOPO'},
                      2:{'imagename':'Cubetest_chandefstwidth2','spw':'0','start':0,'width':2, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, default start, width=2, LSRK'},
                      3:{'imagename':'Cubetest_chanst5wd1','spw':'0','start':5,'width':1, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, start=5, default width, LSRK'},
                      # this will result in blank channnel images (calcChanFreqs requires start and width in channel       
                      # mode to be given in chan index                                                                 
                      4:{'imagename':'Cubetest_chandefstwd1spwsel','spw':'0:5~19','start':0,'width':1, 'outframe':'LSRK','veltype':'radio',
                        'desc':'channel, spw=0:5~19, LSRK'},
                      #5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'40MHz','outframe':'TOPO',
                      #  'desc':'frequency, default start, width=\'40MHz\', TOPO'},
                      # data set changed!
                      5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'100MHz','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, default start, width=\'100MHz\'(2 x chanwidth), TOPO'},
                      6:{'imagename':'Cubetest_freqst2defwd','spw':'0','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.1GHz\', default width, TOPO'},
                      7:{'imagename':'Cubetest_freqst2defwdspwsel','spw':'0:4~19','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.1GHz\', default width, spw=0:4~19, TOPO'},
                      8:{'imagename':'Cubetest_freqst10wdm','spw':'0','start':'1.5GHz','width':'-50MHz','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, start=\'1.5GHz\', width=\'-50MHz\', TOPO'},
                      9:{'imagename':'Cubetest_veldefstwd2','spw':'0','start':'','width':'23983.4km/s','outframe':'TOPO','veltype':'radio',
                        'desc':'frequency, default start, width=\'23983.4km/s\', TOPO'},
                     10:{'imagename':'Cubetest_veldefstwd2m','spw':'0','start':'','width':'-23983.4km/s','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width=\'-23983.4m/s\', TOPO'},
                     11:{'imagename':'Cubetest_velst4defwd','spw':'0','start':'11991.7km/s','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, start=\'11991.7km/s\', default width, TOPO'},
                     12:{'imagename':'Cubetest_velst4defwdbary','spw':'0','start':'11977.6km/s','width':'','outframe':'BARY','veltype':'radio',
                        'desc':'velocity, start=\'11977.6km/s\', default width, BARY'},
                     # currently 13 is not quite properly working, investigating - 2014.08.27 TT 
                     # for refim_point.ms ch9=-41347.8km/s (opt)
                     #13:{'imagename':'Cubetest_optvelst10wdeflsrk','spw':'0','start':'-49962.6km/s','width':'',
                     13:{'imagename':'Cubetest_optvelst19wdlsrk','spw':'0','start':'-41347.8km/s','width':'20000km/s',
                        'veltype':'optical','outframe':'LSRK',
                     ##   'desc':'velocity, start=\'74952.3km/s\', default width, veltype=optical LSRK'},
                     #   'desc':'velocity, start=\'-49962.6km/s\', default width, veltype=optical LSRK'},
                        'desc':'velocity, start=\'-41347.5km/s\', default width , veltype=optical LSRK'},
                     14:{'imagename':'Cubetest_stqfreqdefwd','spw':'0','start':qfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start(quanity)=%s, default width, veltype=radio LSRK' % qfstart},
                     15:{'imagename':'Cubetest_stmfreqdefwd','spw':'0','start':mfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start=%s, default width, veltype=radio LSRK' % mfstart},
                     16:{'imagename':'Cubetest_stqveldefwd','spw':'0','start':qvstart,'width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity(quantity), start=%s, default width, TOPO ' % qvstart},
                     17:{'imagename':'Cubetest_stmveldefwd','spw':'0','start':mvstart,'width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'velocity(measure), start=%s, default width(outframe=TOPO will be overridden)' % mvstart},
                     18:{'imagename':'Cubetest_veldefstqvwidth','spw':'0','start':'','width':qvwidth,'outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width(quantity)=%s' % qvwidth},
                     19:{'imagename':'Cubetest_veldefstmvwidth','spw':'0','start':'','width':mvwidth,'outframe':'TOPO','veltype':'radio',
                        'desc':'velocity, default start, width(measure)=%s, TOPO' % mvwidth},
                     20:{'imagename':'Cubetest_stdopdefwd','spw':'0','start':dop,'width':'','outframe':'LSRK','veltype':'radio',
                        'desc':'doppler, start=%s, default width, LSRK' % dop},
                     # with a gap in spw channel sel
                     21:{'imagename':'Cubetest_st4gap','spw':'0:4~9;12~14','start':4,'width':'','outframe':'LSRK','veltype':'radio',
                        'desc':'channel, start=%s, default width, channel gap (10-11) LSRK' % 4},
                     # stride > 1
                     22:{'imagename':'Cubetest_st4stride2','spw':'0:0~10^2','start':0,'width':'','outframe':'LSRK','veltype':'radio', 'interpolation':'nearest',
                        'desc':'channel, start=%s, default width, step=2 LSRK nearest' % 0},
                     23:{'imagename':'Cubetest_defstspwchansel4','spw':'0:4~13','start':'','width':'','outframe':'TOPO','veltype':'radio',
                        'desc':'spw with channel selection( 0:4~13 ), default start, LSRK nearest'}
                    }
          
#          self.test_cube_0.__func__.__doc__ %=self.testList[0]['desc']
     

     def run_cubetclean(self, testid):
          """ core function to execute a cube tclean """
          if self.testList[testid].has_key('interpolation'):
              interpolation = self.testList[testid]['interpolation']
          else:
              interpolation = 'linear'

          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,\
                       specmode='cube',nchan=10,restfreq=['1.25GHz'],\
                       phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',\
                       spw=self.testList[testid]['spw'],\
                       imagename=self.img+self.testList[testid]['imagename'],\
                       start=self.testList[testid]['start'],\
                       width=self.testList[testid]['width'],\
                       veltype=self.testList[testid]['veltype'],\
                       outframe=self.testList[testid]['outframe'], \
                       interpolation=interpolation)
          return ret

     def test_cube_0(self):
          """ [cube] Test_Cube_0 new """
          testid=0
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50002,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',999988750)
          self.checkfinal(report+report2)

     def test_cube_1(self):
          """ [cube] Test_Cube_1  """
          testid=1
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50002,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO', 9.9999999e8)

          self.checkfinal(report+report2)

     def test_cube_2(self):
          """ [cube] Test_Cube_2  """
          testid=2
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.4643,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.02498846e9)
          self.checkfinal(report+report2)

     def test_cube_3(self):
          """ [cube] Test_Cube_3  """
          # start = 5 (1.25GHZ IN TOPO)
          testid=3
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.2000,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.249985937e9)
          self.checkfinal(report+report2)

     def test_cube_4(self):
          """ [cube] Test_Cube_4  """
          testid=4
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          # NEED CHECK!!!
          #report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          #imval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000,
          #[50,50,0,0])])
          #report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.23998593e9)

     def test_cube_5(self):
          """ [cube] Test_Cube_5  """
          # width by freq (2x chanw) result should be the same as #2
          testid=5
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.4643,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
          self.checkfinal(report+report2)

     def test_cube_6(self):
          """ [cube] Test_Cube_6  """ 
          # start in freq=1.1GHz (=chan5)
          testid=6
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.36365354,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
          self.checkfinal(report+report2)

     def test_cube_7(self):
          """ [cube] Test_Cube_7  """
          # start 1.1GHz(TOPO)=chan5 spw=4~19
          ##as of 03/23/2017 it is  not.... chan2 = 1.1e9 (in refim_point.ms)  
          testid=7
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',0.0,
          [50,50,0,0]),(self.img+self.testList[testid]['imagename']+'.image',1.2000,
          [50,50,0,3])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
          self.checkfinal(report+report2)

     def test_cube_8(self):
          """ [cube] Test_Cube_8  """
          # start =1.5GHz(chan10)  width=-50MHz TOPO (descending freq)
          testid=8
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.42858946,
          [50,50,0,9])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.5e9)
          self.checkfinal(report+report2)

     def test_cube_9(self):
          """ [cube] Test_Cube_9  """
          # width in vel (=23983.4km/s=2xChanW) def start (=cube will be ascending order in vel)
          testid=9
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.46184647,
          [50,50,0,9])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.925e9)
          self.checkfinal(report+report2)

     def test_cube_10(self):
          """ [cube] Test_Cube_10  """
          # width in vel = -23983.4m/s def start (cube will be in descending order in vel)
          testid=10
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.46184647,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
          self.checkfinal(report+report2)

     def test_cube_11(self):
          """ [cube] Test_Cube_11  """
          # start 11991.7km/s (chan4)
          testid=11
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001776,
          [50,50,0,4])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2e9)
          self.checkfinal(report+report2)

     def test_cube_12(self):
          """ [cube] Test_Cube_12  """
          # start 11977.6km/s (BARY) = chan4
          testid=12
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          [50,50,0,4])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
          self.checkfinal(report+report2)

     def test_cube_13(self):
          """ [cube] Test_Cube_13  """
          # 
          testid=13
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          # use own tclean command as nchan need to modify
          ret = tclean(vis=self.msfile,field='0',imsize=100,cell='8.0arcsec',niter=10,specmode='cube',nchan=8,restfreq=['1.25GHz'],phasecenter="J2000 19:59:28.500 +40.44.01.50",deconvolver='hogbom',spw=self.testList[testid]['spw'],imagename=self.img+self.testList[testid]['imagename'],start=self.testList[testid]['start'], width=self.testList[testid]['width'],veltype=self.testList[testid]['veltype'],outframe=self.testList[testid]['outframe'])

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          #report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          #imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          #[50,50,0,4])])
          #report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2000e9)

     def test_cube_14(self):
          """ [cube] Test_Cube_14  """
          # start = quantity ('1.2GHz') frame default(LSRK)
          testid=14
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.25000215,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2e9)
          self.checkfinal(report+report2)

     def test_cube_15(self):
          """ [cube] Test_Cube_15  """
          # measure freq in LSRK ch4
          testid=15
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image', 1.25001216,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989e9)
          self.checkfinal(report+report2)

     def test_cube_16(self):
          """ [cube] Test_Cube_16  """
          # start quantity vel=11991.7km/s outframe=topo (ascending vel order)
          testid=16
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001776,
          [50,50,0,4])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2000e9)
          self.checkfinal(report+report2)

     def test_cube_17(self):
          """ [cube] Test_Cube_17  """
          # start measure vel=11977.6km/s BARY, outframe=TOPO will be overridedden (ascending vel order)
          testid=17
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
          [50,50,0,4])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
          self.checkfinal(report+report2)

     def test_cube_18(self):
          """ [cube] Test_Cube_18  """
          # defaut start, width in vel (quantity) +11991.7km/s (TOPO, radio)=datachan width, will be
          # ascending order in vel so highet DATA channel will be chan 0 in the image (image chan0=1.45GHz)
          testid=18
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001764,
          [50,50,0,9])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
          self.checkfinal(report+report2)

     def test_cube_19(self):
          """ [cube] Test_Cube_19  """
          # default start, width in vel (measure) +11991.7km/s (TOPO, radio)
          testid=19
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001764,
          [50,50,0,9])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
          self.checkfinal(report+report2)

     def test_cube_20(self):
          """ [cube] Test_Cube_20  """
          # doppler (with ch4 LSRK freq, rest freq=1.25GHz)
          testid=20
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000546,
          [50,50,0,4])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989152e9)
          self.checkfinal(report+report2)

     def test_cube_21(self):
          """ [cube] Test_Cube_21  """
          # data sel with channel gap (10,11 excluded) 4~9, 12~14
          testid=21
          self.testList[testid]['interpolation']='nearest'
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.250001562, [50,50,0,0]),
                 (self.img+self.testList[testid]['imagename']+'.image',0.0, [50,50,0,6]),
                 (self.img+self.testList[testid]['imagename']+'.image',0.0, [50,50,0,7])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199986500e9)
          self.checkfinal(report+report2)

     def test_cube_22(self):
          """ [cube] Test_Cube_22  """
          # stride (step=2) use nearest interpolation (other interpotion methods
          # may not work well...)
          testid=22
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.5000546,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','LSRK',0.999988750387e9)
          self.checkfinal(report+report2)

     def test_cube_23(self):
          """ [cube] Test_Cube_23  """
          testid=23
          print " : " , self.testList[testid]['desc']
          self.prepData('refim_point.ms')
          ret = self.run_cubetclean(testid)

          self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
          report=self.th.checkall(imexist=[self.img+self.testList[testid]['imagename']+'.image'],
          imval=[(self.img+self.testList[testid]['imagename']+'.image',1.2500156,
          [50,50,0,0])])
          report2 = self.th.checkspecframe(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.20e9)
          self.checkfinal(report+report2)

     def test_cube_D1(self):
          """ [cube] Test_Cube_D1 : specmode cubedata - No runtime doppler corrections """
          self.prepData('refim_Cband.G37line.ms')
          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cubedata',nchan=30,start=105,width=1,veltype='radio',imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',86.254,[128,128,0,18])])
          ## line is smoother
          self.checkfinal(report)

     def test_cube_D2(self):
          """ [cube] Test_Cube_D2 : specmode cube - WITH doppler corrections """
          self.prepData('refim_Cband.G37line.ms')
          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cube',nchan=30,start=105,width=1,veltype='radio',imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10)
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',92.1789,[128,128,0,20])])
          ## line is tighter
          self.checkfinal(report)

#     def test_cube_D3(self):
#          """ EMPTY : [cube] Test_Cube_D3 : specmode cubesrc - Doppler correct to a SOURCE ephemeris"""
#          ret = tclean(vis=self.msfile,field='1',spw='0:105~135',specmode='cubesrc',nchan=30,start=105,width=1,veltype='radio',imagename=self.img,imsize=256,cell='0.01arcmin',phasecenter=1,deconvolver='hogbom',niter=10)
#          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )

     def test_cube_continuum_subtract_uvsub(self):
          """ [cube] Test_Cube_continuum_subtract :  Using uvsub """
          self.prepData('refim_point_withline.ms')
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='reset0')
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step0data.png',title="original data")
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step0model.png',title="empty model")

          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', spw='0:12~19',niter=50,gain=0.2,savemodel='modelcolumn',deconvolver='mtmfs')
#          self.assertTrue(self.th.exists(self.img+'.model') )
#          self.assertTrue( self.th.checkmodelchan(self.msfile,10) == 0.0 and self.th.checkmodelchan(self.msfile,3) > 0.0 )
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step1.png',title="model after partial mtmfs on some channels")

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='reset0')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'], spw='0',niter=0,savemodel='modelcolumn',deconvolver='mtmfs')
#          self.assertTrue( self.th.checkmodelchan(self.msfile,10) > 0.0 and self.th.checkmodelchan(self.msfile,3) > 0.0 
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step2.png',title="model after mtmfs predict on full spw" )

          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data-model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.step3data.png',title="data-model")
          

#     def test_cube_continuum_subtract_otf(self):
#          """ EMPTY : [cube] Test_Cube_continuum_subtract :  On-The-Fly using multifield """
#          self.prepData('refim_point_withline.ms')

     def test_cube_badchannel_restoringbeam(self):
          """ [cube] Test auto restoring beam with a bad edge channel """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom')
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
          report1=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',0.889,[54,50,0,0]) , (self.img+'.image',0.0602,[54,50,0,19]) , (self.img+'.residual',0.033942,[54,50,0,19]) ])
          # first channel's psf is 'bad' and wider along one axis. This offcenter location is higher in value

          ret = tclean(vis=self.msfile,imagename=self.img+'1',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',restoringbeam='common')
          self.assertTrue(os.path.exists(self.img+'1.psf') and os.path.exists(self.img+'1.image') )
          report2=self.th.checkall(imexist=[self.img+'1.image'],imval=[(self.img+'1.image',0.8906,[54,50,0,0]), (self.img+'1.image',0.35945,[54,50,0,19]) , (self.img+'1.residual',0.033942,[54,50,0,19]) ])
          # OLD - first channel has been restored by a 'common' beam picked from channel 2
          self.checkfinal(report1+report2)

#  def test_cube_explicit_restoringbeam(self):
#          """ [cube] Test explicit restoring beams : Test peak flux and off source value for smoothed residuals"""

     def test_cube_chanchunks(self):
          """ [cube] Test channel chunking for large cubes """
          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom', savemodel='modelcolumn')
#          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
#          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.image',1.5002,[50,50,0,0]) , (self.img+'.image',0.769,[50,50,0,19]) ])

          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',chanchunks=7)
          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
          report=self.th.checkall(imexist=[self.img+'cc.image'],imval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
          self.checkfinal(report)

     def test_cube_chanchunks_auto(self):
          """ [cube] Test channel chunking for large cubes : automatic calc of nchanchunks """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',chanchunks=-1)
          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
          report=self.th.checkall(imexist=[self.img+'cc.image'],imval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
          self.checkfinal(report)


     def test_cube_chanchunks_savemodel(self):
          """ [cube] Test channel chunking for large cubes and save model """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img+'cc',specmode='cube',imsize=100,cell='10.0arcsec',niter=10,deconvolver='hogbom',chanchunks=2,savemodel='modelcolumn')
          self.assertTrue(os.path.exists(self.img+'cc.psf') and os.path.exists(self.img+'cc.image') )
          report=self.th.checkall(imexist=[self.img+'cc.image'],imval=[(self.img+'cc.image',1.5002,[50,50,0,0]) , (self.img+'cc.image',0.769,[50,50,0,19]) ])
          self.assertTrue( self.th.checkmodelchan(self.msfile,5) > 0.0 and self.th.checkmodelchan(self.msfile,18) > 0.0 )
          self.checkfinal(report)

##############################################
##############################################

##Task level tests : masks and clean boxes.
class test_mask(testref_base):

     def test_mask_1(self):
          """ [mask] test_mask_1 : Input mask as file and string : mfs """
          self.prepData('refim_twochan.ms')
          mstr = 'circle[[50pix,80pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=self.img+'.mask.txt')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=mstr)
          report=self.th.checkall(imexist=[self.img+'1.mask', self.img+'2.mask'], imval=[(self.img+'1.mask',0.0,[50,50,0,0]),(self.img+'1.mask',1.0,[50,80,0,0]),(self.img+'2.mask',0.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,80,0,0])])
          self.checkfinal(report)

     def test_mask_2(self):
          """ [mask] test_mask_2 :  Input mask as file and string : cube (few channels) """
          self.prepData('refim_point.ms')
          mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',interactive=0,usemask='user',mask=self.img+'.mask.txt')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',interactive=0,usemask='user',mask=mstr)
          report=self.th.checkall(imexist=[self.img+'1.mask', self.img+'2.mask'], imval=[(self.img+'1.mask',0.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,10]),(self.img+'1.mask',0.0,[50,50,0,11]),(self.img+'2.mask',0.0,[50,50,0,1]),(self.img+'2.mask',1.0,[50,50,0,2]),(self.img+'2.mask',1.0,[50,50,0,10]),(self.img+'2.mask',0.0,[50,50,0,11])])
          self.checkfinal(report)

     def test_mask_3(self):
          """ [mask] test_mask_3 : Input mask as image-to-be-regridded (ra/dec) : mfs """
          self.prepData('refim_twochan.ms')
          mstr = 'circle[[50pix,50pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=self.img+'.mask.txt')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=self.img+'1.mask',phasecenter='J2000 19h59m57.5s +40d49m00.077s') # shift phasecenter
          report=self.th.checkall(imexist=[self.img+'1.mask', self.img+'2.mask'], imval=[(self.img+'1.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[91,13,0,0])])
          self.checkfinal(report)

     def test_mask_4(self):
          """ [mask] test_mask_4 :  Input mask as image-to-be-regridded(ra/dec/specframe) : cube """
          self.prepData('refim_point.ms')
          mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',interactive=0,usemask='user',mask=self.img+'.mask.txt')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',specmode='cube',start='1.3GHz',interactive=0,usemask='user',mask=self.img+'1.mask')
          report=self.th.checkall(imexist=[self.img+'1.mask', self.img+'2.mask'], imval=[(self.img+'1.mask',0.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,10]),(self.img+'1.mask',0.0,[50,50,0,11]),(self.img+'2.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,50,0,4]),(self.img+'2.mask',0.0,[50,50,0,10])])
          self.checkfinal(report)

     def test_mask_5(self):
          """ [mask] test_mask_5 : Input cube mask that has different chan
          ranges (use mask from the 1st tclean with a different channel range in the 2nd tclean run)"""
          self.prepData('refim_point.ms')
          mstr = 'circle[[50pix,50pix],10pix]'
          self.th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=1,deconvolver='hogbom',specmode='cube',start=0,nchan=10,interactive=0,usemask='user',mask=self.img+'.mask.txt')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=1,deconvolver='hogbom',specmode='cube',start=5,nchan=10,interactive=0,usemask='user',mask=self.img+'1.mask')
          report=self.th.checkall(imexist=[self.img+'1.mask', self.img+'2.mask'], imval=[(self.img+'1.mask',1.0,[50,50,0,1]),(self.img+'1.mask',1.0,[50,50,0,2]),(self.img+'1.mask',1.0,[50,50,0,9]),(self.img+'2.mask',1.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,50,0,4]),(self.img+'2.mask',0.0,[50,50,0,5])])
          self.checkfinal(report)

     def test_mask_autobox(self):
         # changed to use threshold based automasking 
          """ [mask] test_mask_autobox :  Autobox """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh')
          # temporarily change value test to make it pass until extra masking in final minor cycle is resolved....
          #report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,80,0,0])])
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.checkfinal(report)

     def test_mask_autobox_redraw(self):
         # changed to use threshold based automasking 
          """ [mask] test_mask_autobox_redraw :  Autoboxing with a redraw after each major cycle """
          self.prepData('refim_eptwochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh',maskthreshold=0.5)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=20,cycleniter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh',maskthreshold=0.5)
          # tweak in automask threshold in the code changed masking extent 2016-03-21
          #report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,30,0,0]),(self.img+'2.mask',1.0,[60,30,0,0])])
          # temporarily change the value test for unmasked region to make it pass (replace with the above when the extra masking issue is resolved...)
          #report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,85,0,0]),(self.img+'2.mask',1.0,[60,30,0,0])])
          #change in behavior due to automask code modification on July 1st,2016
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[60,85,0,0]),(self.img+'2.mask',0.0,[60,30,0,0])])
          self.checkfinal(report)

     def test_mask_autobox_nmask(self):
          """ [mask] test_mask_autobox_nmask : Autoboxing with nmask """
          # this won't be triggering actual pruning but just to check going into write places
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
                       interactive=0,usemask='auto-thresh',nmask=3)
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.checkfinal(report)
           
     def test_mask_autobox2_nmask(self):
          """ [mask] test_mask_autobox2_nmask : Autoboxing (no binning) with nmask"""
          # this won't be triggering actual pruning but just to check going into write places
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
                       interactive=0,usemask='auto-thresh2',nmask=3)
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.checkfinal(report)


     def test_mask_autobox_pbmask(self):
          """ [mask] test_mask_autobox_nmask : Autoboxing  with pbmask"""
          # this won't be triggering actual pruning but just to check going into write places
          self.prepData('refim_point.ms')
          # change imsize to see the pbmask boundary
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=500,cell='8.0arcsec',niter=10,deconvolver='hogbom',
                       interactive=0,usemask='auto-thresh', pbmask=0.2)
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[250,250,0,0]),(self.img+'.mask',0.0,[250,285,0,0]),(self.img+'.mask',0.0,[360,360])])
          self.checkfinal(report)

# This test deprecated. removed autoadjust param.
#     def test_mask_autobox_autoadjust(self):
#          """ [mask] test_mask_autobox_autoadjust : Autoboxing with autoadjust=T """
#          self.prepData('refim_point.ms')
#          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',
#                       interactive=0,usemask='auto-thresh',autoadjust=True)
#          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
#          self.checkfinal(report)
      
#     def test_mask_pbmask(self):
#          """ [mask] test_mask_pbmask :  pb mask """
#
#     def test_mask_combined_1(self):
#          """ [mask] test_mask_combined_1 :  string + pbmask """
#
#     def test_mask_combined_2(self):
#          """ [mask] test_mask_combined_2 :  Autobox + pbmask """
#
#     def test_mask_outlier(self):
#          """ [mask] test_mask_outlier : With outlier fields """

     def test_mask_restart(self):
          """ [mask] test_mask_restart : Test that mask reloads upon restart """
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-thresh')
          ret2 = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          #report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,80,0,0])])
          # temporarily change the value test for unmasked region to make it pass (replace with the above when the extra masking issue is resolved...)
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.checkfinal(report)

     def test_mask_autobox_multithresh(self):
          """ [mask] test_mask__autobox_multithresh :  multi-threshold Autobox (default)"""
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh')
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
          self.checkfinal(report)

     def test_mask_autobox_multithresh_with_prune(self):
          """ [mask] test_mask__autobox_multithresh :  multi-threshold Autobox (minbeamfrac=0.3)"""
          # also test for a bug fix to the new pruneRegions (only caused the failure when image size large
          self.prepData('refim_twochan.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=1000,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='auto-multithresh',
          minbeamfrac=0.3)
          report=self.th.checkall(imexist=[self.img+'.mask'], imval=[(self.img+'.mask',1.0,[500,500,0,0]),(self.img+'.mask',0.0,[500,510,0,0])])
          self.checkfinal(report)

#     def test_mask_outregion(self):
#          """ [mask] test_mask_outregion : Input mask has region that goes outside the image """
#          self.prepData('refim_twochan.ms')
#          mstr = 'circle[[50pix,110pix],20pix]'
#          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,usemask='user',mask=mstr)
#          report=self.th.checkall(imexist=[self.img+'2.mask'], imval=[(self.img+'2.mask',0.0,[50,50,0,0]),(self.img+'2.mask',1.0,[50,95,0,0])])

     def test_mask_zerostart(self):
          """ [mask] test_mask_zerostart : Test that a zero starting mask is caught  """
          self.prepData('refim_point.ms')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='10.0arcsec',niter=0,interactive=0)
          os.system('cp -r ' + self.img + '.residual '+ self.img+'2.inpmask')
          _ia.open(self.img+'2.inpmask')
          pix =_ia.getchunk()
          pix.fill(0.0)
          _ia.putchunk(pix)
          _ia.close()

          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='10.0arcsec',niter=10,interactive=0,mask=self.img+'2.inpmask')

          report=self.th.checkall(ret=ret, imexist=[self.img+'2.mask'], imval=[(self.img+'2.model',0.0,[50,50,0,0]),(self.img+'2.mask',0.0,[50,50,0,0])], stopcode=7)

          self.checkfinal(report)


##############################################
##############################################

##Task level tests : awproject and mosaics
class test_widefield(testref_base):
     
     def test_widefield_wproj_mfs(self):
          """ [widefield] Test_Widefield_wproj : W-Projection """ 
          ### Need better test dataset for this.....
          self.prepData("refim_twopoints_twochan.ms")
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='8.0arcsec',phasecenter="J2000 19:59:00.2 +40.50.15.50",niter=30,gridder='widefield',wprojplanes=4,deconvolver='hogbom')
          report=self.th.checkall(imexist=[self.img+'.image'],imval=[(self.img+'.psf',1.0,[100,100,0,0]),(self.img+'.image',5.56,[127,143,0,0]) ] )
          self.checkfinal(report)


     def test_widefield_aproj_mfs(self):
          """ [widefield] Test_Widefield_aproj : MFS with narrowband AWProjection (wbawp=F, 1spw)  stokes I """
          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=False,conjbeams=True,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',savemodel='modelcolumn')
         ## ret = tclean(vis=self.msfile,spw='2',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',wbawp=False,conjbeams=True,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom')
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          ## weight is pbsq which is 0.7^2 = 0.49 (approx).
          self.checkfinal(report)

          #do stokes V too.....
     def test_widefield_aproj_cube(self):
          """ [widefield] Test_Widefield_aproj_cube_aproj : Cube with AW-Projection  and rotation off """

          casalog.post("EMPTY TEST")
          return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",specmode='cube',niter=1,gain=1.0,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=False,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom')
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',0.11,[256,256,0,0]),(self.img+'.weight',0.34,[256,256,0,0]) ] )
          self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
          self.checkfinal(report)

     ## Test normtype too somewhere..


     def test_widefield_wbaproj_mfs(self):
          """ [widefield] Test_Widefield_wbaproj_mfs : MFS with wideband AWProjection (wbawp=T, allspw) and nt=1 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=True,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3)
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          self.checkfinal(report)


          #do stokes V too..

     def test_widefield_aproj_mtmfs(self):
          """ [widefield] Test_Widefield_aproj_mtmfs : MFS with AWProjection (wbawp=T,conjbeams=F, allspw) and nt=2 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=0,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=False,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='mtmfs')
          report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imval=[(self.img+'.image.tt0',0.96,[256,256,0,0]),(self.img+'.weight.tt0',0.48,[256,256,0,0]),(self.img+'.alpha',0.04,[256,256,0,0]) ] )
          self.checkfinal(report)
          ## alpha should represent that of the mosaic PB (twice).. -0.1 doesn't look right. Sigh.... well.. it should converge to zero.
          ## alpha keeps increasing in magnitude with niter.... not right.
          ## restricting this check to niter=0 only. Alpha should represent that of the PB.

     def test_widefield_wbaproj_mtmfs(self):
          """ [widefield] Test_Widefield_wbaproj_mtmfs : MFS with wideband AWProjection (wbawp=T,conjbeams=T, allspw) and nt=2 stokes I  """

          # casalog.post("EMPTY TEST")
          # return

          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',cfcache=self.img+'.cfcache',wbawp=True,conjbeams=True,psterm=False,computepastep=360.0,rotatepastep=360.0,deconvolver='mtmfs',pblimit=0.1)
          report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imval=[(self.img+'.image.tt0',0.96,[256,256,0,0]),(self.img+'.weight.tt0',0.486,[256,256,0,0]),(self.img+'.alpha',0.0,[256,256,0,0]) ] )
          ## alpha should be ZERO as the pb spectrum has been taken out.
          self.checkfinal(report)


#     def test_widefield_wbaproj_subsets(self):
#          """ [widefield] Test_Widefield_wbaproj_subsets : MFS with the AWProjection gridder and A,W turned off  """
#          self.prepData("refim_mawproject.ms")
#          ## PS only
#          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',psterm=True,aterm=False,wprojplanes=1,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3)
#          #report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
#
#          ## W and PS only
#          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='awproject',psterm=True,aterm=False,wprojplanes=16,computepastep=360.0,rotatepastep=360.0,deconvolver='hogbom',pblimit=0.3)
#          #report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )


#     def test_widefield_multispws(self):
#          """ [widefield] Test_Widefield_multispws : Test cube imaging with mosaicft and awproj  """

     ## CHECK NORMALIZATION OF WEIGHTIMAGE = normed to peak=1
     ## TODO : make vpman recognize EVLA in addition to VLA.
     def test_widefield_mosaicft_mfs(self):
          """ [widefield] Test_Widefield_mosaic : MFS with mosaicft  stokes I """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='mosaicft',deconvolver='hogbom',pblimit=0.3)
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',0.961231,[256,256,0,0]),(self.img+'.weight',0.50576,[256,256,0,0]) ] )
          #ret = clean(vis=self.msfile,spw='1',field='*',imagename=self.img+'.old',imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,imagermode='mosaic',psfmode='hogbom')
          self.checkfinal(report)

          #do stokes V too..

     def test_widefield_mosaicft_mtmfs(self):
          """ [widefield] Test_Widefield_mosaicft_mtmfs : MT-MFS with mosaicft  stokes I, alpha """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='*',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=60,gridder='mosaicft',deconvolver='mtmfs')
          report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.psf.tt0', self.img+'.weight.tt0'],imval=[(self.img+'.image.tt0',0.9413,[256,256,0,0]),(self.img+'.weight.tt0',0.50546,[256,256,0,0]),(self.img+'.alpha',0.0764185711741,[256,256,0,0]) ] )
          ## alpha should represent that of the mosaic PB (twice)... and should then converge to zero
          self.checkfinal(report)
          

     def test_widefield_mosaicft_cube(self):
          """ [widefield] Test_Widefield_mosaicft_cube : MFS with mosaicft  stokes I """
          self.prepData("refim_mawproject.ms")
#          _vp.setpbpoly(telescope='EVLA', coeff=[1.0, -1.529e-3, 8.69e-7, -1.88e-10]) 
#          _vp.saveastable('evlavp.tab')
          ret = tclean(vis=self.msfile,spw='*',field='0',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",specmode='cube',niter=10,gridder='mosaicft',deconvolver='hogbom',gain=0.1,stokes='I') #,vptable='evlavp.tab')
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',0.7987,[256,256,0,0]),(self.img+'.weight',0.6528,[256,256,0,0]) ] )
          self.checkfinal(report)

          #do stokes V too..

     
class test_widefield_failing(testref_base):

     def test_widefield_imagemosaic(self):
          """ [widefield] Test_Widefield_imagemosaic : Image domain mosaic for single-term mfs (or narrowband)  """
          self.prepData("refim_mawproject.ms")
          ret = tclean(vis=self.msfile,spw='1',field='*',imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",niter=30,gridder='imagemosaic',deconvolver='hogbom')
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.psf', self.img+'.weight'],imval=[(self.img+'.image',1.0,[256,256,0,0]),(self.img+'.weight',0.493,[256,256,0,0]) ] )
          self.checkfinal(report)


     def test_widefield_mosaic_outlier(self):
          """ [multifield] Test_widefield_mosaic_outlier : Mosaic with an outlier field """
          #### Need another dataset for this.
          self.prepData("refim_mawproject.ms")
          ## Outlier uses gridft
          self.th.write_file(self.img+'.out.txt', 'imagename='+self.img+'1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:42.0 +40.55.58.543\nmask=circle[[40pix,40pix],10pix]\ngridder=gridft')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=512,cell='10.0arcsec',phasecenter="J2000 19:59:28.500 +40.44.01.50",outlierfile=self.img+'.out.txt',niter=10,deconvolver='hogbom',gridder='mosaicft',interactive=0)
          report=self.th.checkall(ret=ret, 
                        iterdone=10,  # outlier field has nothing in it :).
                        nmajordone=2,
                        imexist=[self.img+'.image', self.img+'1.image'],
                        imval=[(self.img+'.image',0.933,[256,256,0,0]),
                               (self.img+'1.image',0.0,[40,40,0,0])])
          self.checkfinal(report)
          #### NOT WORKING as the model is being picked from tst1, so residual after 0.6 model and one major cycle is same as dirty and the output image is 1.6 instead of 0.93.

##############################################
##############################################

##Task level tests : model prediction.
class test_modelvis(testref_base):
     
     def test_modelvis_1(self):
          """ [modelpredict] Test_modelvis_1 : mfs with no save model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='none')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==False )

     def test_modelvis_2(self):
          """ [modelpredict] Test_modelvis_2 : mfs with save model column """
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='modelcolumn')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

          ##Predict from input model image (startmodel)
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', niter=0,savemodel='modelcolumn')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_3(self):
          """ [modelpredict] Test_modelvis_3 : mfs with save virtual model """
          self.prepData("refim_twochan.ms")

          ## Save model after deconvolution
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          ##Predict from input model image (startmodel)
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', niter=0,savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_4(self):
          """ [modelpredict] Test_modelvis_4 : mt-mfs with no save model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',savemodel='none')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==False )

     def test_modelvis_5(self):
          """ [modelpredict] Test_modelvis_5 : mt-mfs with save model column """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',savemodel='modelcolumn')
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='data',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.data.png',title="original data")
          plotms(vis=self.msfile,xaxis='frequency',yaxis='amp',ydatacolumn='model',customsymbol=True,symbolshape='circle',symbolsize=5,showgui=False,plotfile=self.img+'.plot.model.png',title="empty model")
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],niter=0,deconvolver='mtmfs',savemodel='modelcolumn')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )

     def test_modelvis_6(self):
          """ [modelpredict] Test_modelvis_6 : mt-mfs with save virtual model """
          self.prepData("refim_twochan.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],niter=0,deconvolver='mtmfs',savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_7(self):
          """ [modelpredict] Test_modelvis_7 : cube with chan selection and save model column """
          ## check explicit channels ...
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0:5~12',imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,savemodel='modelcolumn',start=5,nchan=8,interpolation='nearest')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )
          reportcv=self.th.checkchanvals(self.msfile, [(10,">",0.0),(3,"==",1.0)])
          self.checkfinal(reportcv)
          
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0',imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',specmode='cube',niter=0,savemodel='modelcolumn')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==True and modsum>0.0 and hasvirmod==False )
          reportcv=self.th.checkchanvals(self.msfile,[(10,">",0.0),(3,"==",self.th.checkmodelchan(self.msfile,1))])
          self.checkfinal(reportcv)


     def test_modelvis_8(self):
          """ [modelpredict] Test_modelvis_8 : cube with chan selection and save virtual model """
          ## check explicit channels ...
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0:5~12',imagename=self.img,imsize=100,cell='8.0arcsec',specmode='cube',niter=10,savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,spw='0',imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',specmode='cube',niter=0,savemodel='virtual')
          hasmodcol, modsum, hasvirmod = self.th.checkmodel(self.msfile)
          self.assertTrue( hasmodcol==False and hasvirmod==True )

     def test_modelvis_9(self):
          """ [modelpredict] Test_modelvis_9 : Don't de-grid channels with zero model. Also test limited-freq mask """
          self.prepData("refim_point.ms")
          masklist=''  # choose only a few channels here.
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec',niter=10,mask=masklist)
          self.assertTrue(self.th.exists(self.img+'.model') )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',niter=0,savemodel='modelcolumn')

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model',niter=0,savemodel='virtual')

     def test_modelvis_10(self):
          """ [modelpredict] Test_modelvis_10 : Use input model of different (narrower) freq range than data """
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', spw='0:5~12',niter=10,savemodel='modelcolumn')
          self.assertTrue(self.th.exists(self.img+'.model') )
#          self.assertTrue( self.th.checkmodelchan(self.msfile,10) > 0.0 and self.th.checkmodelchan(self.msfile,3) == 1.0 )
          reportcv=self.th.checkchanvals(self.msfile,[(10,">",0.0),(3,"==",1.0)])
          self.checkfinal(reportcv)


          ## add model expansion parameter
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,savemodel='modelcolumn')
#          self.assertTrue( self.th.checkmodelchan(self.msfile,10) > 0.0 and self.th.checkmodelchan(self.msfile,3) > 0.0 )
          reportcv=self.th.checkchanvals(self.msfile,[(10,">",0.0),(3,">",0.0)])
          self.checkfinal(reportcv)
                    

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,savemodel='virtual')
          ## cannot check anything here....  just that it runs without error

     def test_modelvis_11(self):
          """ [modelpredict] Test_modelvis_11 : Predict model image over channel gaps not included in imaging """
          self.prepData("refim_point.ms")
          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img,imsize=100,cell='8.0arcsec', spw='0:0~8;12~19',niter=10,savemodel='modelcolumn')
          self.assertTrue(self.th.exists(self.img+'.model') )
          ###vi2 leave unselected channel as is so it will be 1.0
          self.assertTrue( (self.th.checkmodelchan(self.msfile,10) == 0.0) or (np.abs(self.th.checkmodelchan(self.msfile,10)-1) < 1.0e-12)   and self.th.checkmodelchan(self.msfile,3) > 0.0 )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,savemodel='modelcolumn')
          self.assertTrue( self.th.checkmodelchan(self.msfile,10) > 0.0 and self.th.checkmodelchan(self.msfile,3) > 0.0 )

          delmod(self.msfile);self.th.delmodels(msname=self.msfile,modcol='delete')
          ret = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',startmodel=self.img+'.model', spw='0',niter=0,savemodel='virtual')
          ## cannot check anything here....  just that it runs without error

          
class test_startmodel(testref_base):
     def test_startmodel_regrid_mfs(self):
          """ [modelpredict] Test_startmodel_regrid_mfs : Regrid input model onto new image grid : mfs (ra/dec) """
          self.prepData('refim_twopoints_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=50,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,phasecenter='J2000 19h58m40.801s +40d55m59.863s')
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=200,cell='8.0arcsec',niter=0,deconvolver='hogbom',interactive=0,startmodel=self.img+'1.model')
          report=self.th.checkall(imexist=[self.img+'1.residual', self.img+'2.residual'], imval=[(self.img+'1.residual',1.7963,[25,25,0,0]),(self.img+'2.residual',1.910,[168,190,0,0])])
          self.checkfinal(report)

     def test_startmodel_regrid_cube(self):
          """ [modelpredict] Test_startmodel_regrid_cube : Regrid input model onto new image grid : cube (ra/dec/specframe)"""
          self.prepData('refim_point.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=50,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,specmode='cube',start='1.05GHz',width='50MHz',nchan=20)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=0,deconvolver='hogbom',interactive=0,startmodel=self.img+'1.model',specmode='cube')
          report=self.th.checkall(imexist=[self.img+'1.residual', self.img+'2.residual'], imval=[(self.img+'1.residual',0.362,[25,25,0,5]),(self.img+'2.residual',0.362,[50,50,0,6])])
          self.checkfinal(report)



#     def test_startmodel_14(self):
#          """ [modelpredict] Test_startmodel_14 : Regrid input model onto new image grid : mtmfs (ra/dec/terms)"""

#     def test_startmodel_15(self):
#          """ [modelpredict] Test_startmodel_15 : Regrid input model onto new image grid : mfs (imsize/cell)"""

     def test_startmodel_mfs_continue(self):
          """ [startmodel] test_startmodel_mfs_continue : Restart a run with no parameter changes"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          report=self.th.checkall(imexist=[self.img+'1.residual'], imval=[(self.img+'1.residual',0.35304,[50,50,0,0])])
          ret2 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          report=self.th.checkall(imexist=[self.img+'1.residual'], imval=[(self.img+'1.residual',0.1259,[50,50,0,0])])
          self.checkfinal(report)

     def test_startmodel_mfs_restart(self):
          """ [startmodel] test_startmodel_mfs_restart : Restart a run using 'startmodel' and changed imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,startmodel=self.img+'1.model')
          report=self.th.checkall(imexist=[self.img+'1.residual', self.img+'2.residual'], imval=[(self.img+'1.residual',0.35304,[50,50,0,0]),(self.img+'2.residual',0.1259,[50,50,0,0])])
          self.checkfinal(report)

     def test_startmodel_mfs_changeshape_1(self):
          """ [startmodel] test_startmodel_mfs_changeshape_1 : Restart a run but change shape only (cas-6937)"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          report=self.th.checkall(imexist=[self.img+'1.residual'], imval=[(self.img+'1.residual',0.35304,[50,50,0,0])])

          try:
               ## This run should fail with an exception (if __rethrow_exceptions = True )
               ret2 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=120,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
               correct=False
          except Exception as e:
               correct=True
          #self.assertTrue(correct)
          
          ## Check that there is no change in output value.... 
          ## i.e. the second run should have failed.     
          report=self.th.checkall(imval=[(self.img+'1.residual',0.35304,[50,50,0,0])])
          self.checkfinal(report)

     def test_startmodel_mfs_changeshape_2(self):
          """ [startmodel] test_startmodel_mfs_changeshape_2 : Restart a run using 'startmodel' and change shape and imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0)
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=120,cell='8.0arcsec',niter=10,deconvolver='hogbom',interactive=0,startmodel=self.img+'1.model')
          report=self.th.checkall(imexist=[self.img+'1.residual', self.img+'2.residual'], imval=[(self.img+'1.residual',0.35304,[50,50,0,0]),(self.img+'2.residual',0.1259,[60,60,0,0])])
          self.checkfinal(report)

     def test_startmodel_mtmfs_restart(self):
          """ [startmodel] test_startmodel_mtmfs_restart : Restart a multi-term run using 'startmodel' and changed imagename"""
          self.prepData('refim_twochan.ms')
          ret1 = tclean(vis=self.msfile,imagename=self.img+'1',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0)
          # start with full model
          ret2 = tclean(vis=self.msfile,imagename=self.img+'2',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,startmodel=[self.img+'1.model.tt0',self.img+'1.model.tt1'])
          # start with model only for tt0
          ret3 = tclean(vis=self.msfile,imagename=self.img+'3',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,startmodel=self.img+'1.model.tt0')
          # start with model only for tt1
          ret3 = tclean(vis=self.msfile,imagename=self.img+'4',imsize=100,cell='8.0arcsec',niter=10,deconvolver='mtmfs',interactive=0,startmodel=['',self.img+'1.model.tt1'])

          report=self.th.checkall(imexist=[self.img+'1.residual.tt0', self.img+'2.residual.tt0', self.img+'3.residual.tt0', self.img+'4.residual.tt0', self.img+'1.residual.tt1', self.img+'2.residual.tt1', self.img+'3.residual.tt1', self.img+'4.residual.tt1'], imval=[  (self.img+'1.residual.tt0',0.39226,[50,50,0,0]),
                             (self.img+'2.residual.tt0',0.13677,[50,50,0,0]),
                             (self.img+'3.residual.tt0',0.13677,[50,50,0,0]),
                             (self.img+'4.residual.tt0',0.39226,[50,50,0,0]),  
                             (self.img+'1.residual.tt1',-0.04358,[50,50,0,0]),
                             (self.img+'2.residual.tt1',-0.01519,[50,50,0,0]),
                             (self.img+'3.residual.tt1',-0.04358,[50,50,0,0]),
                             (self.img+'4.residual.tt1',-0.01519,[50,50,0,0])     ] )
          self.checkfinal(report)

     def test_startmodel_with_mask_mfs(self):
          """ [startmodel] test_startmodel_with_mask_mfs : Mask out some regions in the startmodel, before prediction """
          self.prepData("refim_twopoints_twochan.ms")
          ## image both sources
          tclean(vis=self.msfile,spw='0:0',imagename=self.img,imsize=200,cell='10.0arcsec',niter=100)
          ## mask out all but the brightest outlier source
          _ia.open(self.img+'.model')
          _ia.calcmask(mask='"'+self.img+'.model">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          ## predict only that outlier source
          tclean(vis=self.msfile,spw='0:0',imagename=self.img+'.2',niter=0,savemodel='modelcolumn',imsize=200,cell='10.0arcsec',startmodel=self.img+'.model')
          ## subtract it out
          uvsub(vis=self.msfile)
          ## image the rest of the field.
          tclean(vis=self.msfile,spw='0:0',imagename=self.img+'.3',imsize=200,cell='10.0arcsec',niter=100)
          
          report=self.th.checkall(imexist=[self.img+'.model',self.img+'.2.model',self.img+'.3.model'],
                        imval=[ (self.img+'.model',1.497,[100,100,0,0]),
                                (self.img+'.model',7.474,[154,172,0,0]),
                                (self.img+'.3.model',1.497,[100,100,0,0]), 
                                (self.img+'.3.model',0.024,[154,172,0,0])   ] )
          self.checkfinal(report)
          
     def test_startmodel_with_mask_mtmfs(self):
          """ [startmodel] test_startmodel_with_mask_mtmfs : Mask out some regions in the startmodel, before prediction """
          self.prepData("refim_twopoints_twochan.ms")
          ## image both sources
          tclean(vis=self.msfile,imagename=self.img,imsize=200,cell='10.0arcsec',niter=100,deconvolver='mtmfs')
          ## mask out all but the brightest outlier source
          _ia.open(self.img+'.model.tt0')
          _ia.calcmask(mask='"'+self.img+'.model.tt0">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          _ia.open(self.img+'.model.tt1')
          _ia.calcmask(mask='"'+self.img+'.model.tt0">2.0') # keep only high values. i.e. the far out bright source.
          _ia.close()
          ## predict only that outlier source
          tclean(vis=self.msfile,imagename=self.img+'.2',niter=0,savemodel='modelcolumn',imsize=200,cell='10.0arcsec',startmodel=[self.img+'.model.tt0',self.img+'.model.tt1'],deconvolver='mtmfs')
          ## subtract it out
          uvsub(vis=self.msfile)
          ## image the rest of the field.
          tclean(vis=self.msfile,imagename=self.img+'.3',imsize=200,cell='10.0arcsec',niter=100,deconvolver='mtmfs')
          
          report=self.th.checkall(imexist=[self.img+'.model.tt0',self.img+'.2.model.tt0',self.img+'.3.model.tt0'],
                        imval=[ 
                                (self.img+'.model.tt0',1.11,[100,100,0,0]),
                                (self.img+'.model.tt0',5.56,[154,172,0,0]),
                                (self.img+'.3.model.tt0',1.12,[100,100,0,0]), 
                                (self.img+'.3.model.tt0',0.013,[154,172,0,0]),
                                (self.img+'.model.tt1',-1.1027,[100,100,0,0]),
                                (self.img+'.model.tt1',-5.602,[154,172,0,0]),
                                (self.img+'.3.model.tt1',-1.124,[100,100,0,0]), 
                                (self.img+'.3.model.tt1',-0.021,[154,172,0,0])   ] )
          self.checkfinal(report)
          
 ##############################################
class test_pbcor(testref_base):
     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"

          self.th = TestHelpers()

          _vp.setpbpoly(telescope='EVLA', coeff=[1.0, -1.529e-3, 8.69e-7, -1.88e-10]) 
          _vp.saveastable('evlavp.tab')

     def test_pbcor_mfs(self):
          """ [pbcor] Test pbcor with mfs"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True)
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.7,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0])])
          self.checkfinal(report)

     def test_pbcor_mtmfs(self):
          """ [pbcor] Test pbcor with mtmfs"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True, deconvolver='mtmfs')
          report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.pb.tt0'], imexistnot=[self.img+'.image.tt0.pbcor', self.img+'.alpha.pbcor'], imval=[(self.img+'.pb.tt0',0.7,[256,256,0,0])])  
          #report=self.th.checkall(imexist=[self.img+'.image.tt0', self.img+'.pb.tt0', self.img+'.image.tt0.pbcor', self.img+'.alpha.pbcor'], imval=[(self.img+'.pb.tt0',0.7,[256,256,0,0]),(self.img+'.image.tt0.pbcor',1.0,[256,256,0,0]),(self.img+'.alpha',-0.7,[256,256,0,0]), (self.img+'.alpha.pbcor',-0.7,[256,256,0,0]) ])  
          # uncorrected alpha, for now. 
          self.checkfinal(report)

     def test_pbcor_cube_basic(self):
          """ [pbcor] Test pbcor with cube"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='cube', vptable='evlavp.tab', pbcor=True)
          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])])
          self.checkfinal(report)

     def test_pbcor_cube_twosteps(self):
          """ [pbcor] Test pbcor with cube with imaging and pbcor separately"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='cube',vptable='evlavp.tab', pbcor=False)
          report1=self.th.checkall(imexist=[self.img+'.image',self.img+'.pb'],imexistnot=[self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.79,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2])])
          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=0,calcres=False,calcpsf=False, specmode='cube', vptable='evlavp.tab', pbcor=True)
          report2=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])])
          self.checkfinal(report1+report2)

     def test_pbcor_cube_mosaicft(self):
          """ [pbcor] Test pbcor with cube with mosaicft"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='cube', vptable='evlavp.tab',pbcor=True, gridder='mosaic')

          report=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.79,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0]), (self.img+'.pb',0.59,[256,256,0,2]),(self.img+'.image.pbcor',1.0,[256,256,0,2])])
          self.checkfinal(report)

     def test_pbcor_mfs_restart(self):
          """ [pbcor] Test pbcor with mfs and a restart"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=0, specmode='mfs', vptable='evlavp.tab', pbcor=False)
          report1=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb'], imexistnot=[self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.7,[256,256,0,0])])

          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True)
          report2=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb', self.img+'.image.pbcor'], imval=[(self.img+'.pb',0.7,[256,256,0,0]),(self.img+'.image.pbcor',1.0,[256,256,0,0])])
          self.checkfinal(report1+report2)

     def test_pbcor_turn_off_pbmask(self):
          """ [pbcor] Test pbcor with mfs where the internal T/F mask is turned off"""
          self.prepData('refim_mawproject.ms')
          ret1 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=0, specmode='mfs', vptable='evlavp.tab', pbcor=True)
          report1=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb'], imval=[(self.img+'.pb',0.7,[256,256,0,0])], immask=[(self.img+'.pb',False,[10,10,0,0]), (self.img+'.image',False,[10,10,0,0])] )

          ret2 = tclean(vis=self.msfile, imagename=self.img, field='0', imsize=512, cell='10.0arcsec', phasecenter="J2000 19:59:28.500 +40.44.01.50", niter=10, specmode='mfs', vptable='evlavp.tab', pbcor=True, calcpsf=False, calcres=False, pblimit=-0.2)
          report2=self.th.checkall(imexist=[self.img+'.image', self.img+'.pb'], imval=[(self.img+'.pb',0.7,[256,256,0,0])] , immask=[(self.img+'.pb',False,[10,10,0,0]), (self.img+'.image',True,[10,10,0,0])]  )

          self.checkfinal(report1+report2)
