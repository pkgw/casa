import os
import sys
import shutil
import subprocess
import numpy
from __main__ import default
from tasks import *
from taskinit import *
import unittest
import inspect


#################################
## Tool level test script for the refactored imager framework
#################################
from refimagerhelper import PySynthesisImager, PyParallelContSynthesisImager,PyParallelCubeSynthesisImager, PyParallelDeconvolver, PyParallelImagerHelper,ImagerParameters

import subprocess
#### Specify parameters.
def getparams(testnum=1,testid=0, parallelmajor=False,parallelminor=False,parallelcube=False):

     # Iteration parameters - common to all tests below
     niter=200
     cycleniter=100
     threshold=0.001
     loopgain=0.1

     restoringbeam=[]
     #restoringbeam='30.0arcsec'

     # Interaction ON or OFF
     interactive=False
     #mask='ttt.mask'
     mask=''
     #mask = 'circle[[50pix,50pix],10pix]'

     if(testnum==23):  ## Cube with AWP (and mosft)
          casalog.post("==================================");
          casalog.post("Test 23 Cube with AWP (and mosft)");
          casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/reg_mawproject.ms', #_offcenter.ms',
                                       field='*',scan='',
                                       spw='*',\
                                       usescratch=True,readonly=True,\
                                       mode='cube',\
                                       imagename='mytest0', nchan=-1,\
                                       imsize=[512,512],\
                                       cell=['10.0arcsec','10.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='awprojectft',
                                       startmodel='', weighting='natural',\
                                       aterm=True, psterm=False, mterm=True,\
                                       wbawp = False,
                                       cfcache = "perm.mytest0.cube.cfcache.mos",
                                       dopointing = False, dopbcorr = True, conjbeams = True,
                                       computepastep =360.0, rotatepastep =5.0,\
                                       deconvolver='hogbom',\
                                       #pblimit=0.1,normtype='flatsky',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==22):  ## 22 image-field, mfs --- readonly/savevirtualmodel/savemodelcolumn.
          casalog.post("==================================");
          casalog.post("Test 22 image-field, mfs --- readonly/savevirtualmodel/savemodelcolumn.");
          casalog.post("==================================");

          msname = 'DataTest/point_twospws.ms'
          resetmodelcol(msname)  ## Set model column to zero
          delmod(msname)  ## Get rid of OTF model
          delmodkeywords(msname) ## Get rid of extra OTF model keywords that sometimes persist...

          print("At start : ", checkmodel(msname))

          testList = {
               ## readonly
               0:{'readonly':True,  'usescratch':True, 'deconvolver':'hogbom', 'ntaylorterms':1, 'mtype':'default'},
               ## readonly
               1:{'readonly':True,  'usescratch':False, 'deconvolver':'hogbom', 'ntaylorterms':1, 'mtype':'default'},
               ## save model column in last major cycle
               2:{'readonly':False,  'usescratch':True, 'deconvolver':'hogbom', 'ntaylorterms':1, 'mtype':'default'},
               ## save virtual model in last major cycleexit
               3:{'readonly':False,  'usescratch':False, 'deconvolver':'hogbom', 'ntaylorterms':1, 'mtype':'default'},
               ## Multi-term test : save model column in last major cycle
               4:{'readonly':False,  'usescratch':True, 'deconvolver':'mtmfs', 'ntaylorterms':2, 'mtype':'multiterm'},
               ## Multi-term test : save virtual model in last major cycle
               5:{'readonly':False,  'usescratch':False, 'deconvolver':'mtmfs', 'ntaylorterms':2, 'mtype':'multiterm'}
               }

          ###  Note : 2 or 3 has a bug. When trying to read the virtual model (in plotms), the plotted values are
          ###           different when using a modelcolumn or a virtual model. Needs more checking.

          if testid > 5:
               print('No such test.')
               return

          paramList = ImagerParameters(msname=msname,field='0',spw='0',\
                                       usescratch=testList[testid]['usescratch'],readonly=testList[testid]['readonly'],\
                                       mode='mfs',\
                                       ntaylorterms=testList[testid]['ntaylorterms'],mtype=testList[testid]['mtype'],
                                       imagename='mytest0', nchan=1,\
                                       imsize=[110,110],\
                                       cell=['8.0arcsec','8.0arcsec'],stokes='I',\
                                       phasecenter=0,
                                       ftmachine='gridft', startmodel='', weighting='briggs',\
                                       deconvolver=testList[testid]['deconvolver'],\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)
          print('Run to check model  :  checkmodel("' + msname + '")')

     if(testnum==21):  ## 21 image-field, mfs --- Multiple Stokes planes -- Clark
          casalog.post("==================================");
          casalog.post("Test 21 image-field, mfs --- Multiple Stokes planes -- Clark.");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/point_stokes.ms',field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[110,110],\
                                       cell=['8.0arcsec','8.0arcsec'],stokes='I',\
                                       phasecenter=0,
                                       ftmachine='gridft', startmodel='', weighting='briggs',\
                                       deconvolver='clarkstokes',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==20):  ## 1 image-field, mfs, multiscale
          casalog.post("==================================");
          casalog.post("Test 20 MFS 1 term + multiscale");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/eptwochan.ms',field='0',spw='0:0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[200,200],\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       ### center
                                       #phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ### offset
                                       #phasecenter="J2000 19:59:23.591 +40.44.01.50",\
                                       phasecenter=0,
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='multiscale',scales=[0,20,40,100],\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==19): ## mode=cubesrc (For interface test purpose only, the mode is not fully implemeted yet)
          casalog.post("==================================");
          casalog.post("Test 19 image-field, cubesrc --- Imaging with various cube parameter specifications");
          #casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/point_twospws.ms', field='0',\
                                       spw='0',\
                                       usescratch=True,readonly=True,\
                                       imagename="CubesrcTest", mode='cubesrc',\
                                       nchan=10,\
                                       start="1.1GHz",\
                                       width="",\
                                       veltype='radio',\
                                       sysvel="-3000km/s",\
                                       sysvelframe="LSRK",\
                                       restfreq=["1.0GHz"],\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==18): ## mode=cube (with doppler correction)
          casalog.post("==================================");
          casalog.post("Test 18 image-field, cube --- With doppler corrections");
          #casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/Cband.G37line.DopplerTest.ms', field='1',\
                                       spw='0:105~135',\
                                       usescratch=True,readonly=True,\
                                       imagename="mytest0", mode='cube',\
                                       nchan=30,\
                                       start=105,\
                                       width=1,\
                                       veltype='radio',\
                                       imsize=[256,256],\
                                       cell='0.01arcmin',\
                                       phasecenter=1,
                                       ftmachine='gridft', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==17): ## mode=cubedata
          casalog.post("==================================");
          casalog.post("Test 17 image-field, cubedata --- No runtime dopper corrections");
          #casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/Cband.G37line.DopplerTest.ms', field='1',\
                                       spw='0:105~135',\
                                       usescratch=True,readonly=True,\
                                       imagename="mytest0", mode='cubedata',\
                                       nchan=30,\
                                       start=105,\
                                       width=1,\
                                       veltype='radio',\
                                       imsize=[256,256],\
                                       cell='0.01arcmin',\
                                       phasecenter=1,
                                       ftmachine='gridft', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==16):  ## 2 image-fields, mfs - one with nterms 1 and one with nterms 2
          casalog.post("==================================");
          casalog.post("Test 16 image-fields, mfs (one with nterms=1. one with nterms=2).");
          casalog.post("==================================");

          write_file('out16.txt', 'imagename=mytest1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nntaylorterms=2\nmtype=multiterm\nreffreq=1.5GHz\ndeconvolver=mtmfs')
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',\
                                       field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       outlierfile='out16.txt',\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                        restoringbeam=restoringbeam,
                                      interactive=interactive,mask=mask)

     if(testnum==15):  ## 2 image-fields, mfs, Overlapping models. Both multi-term
          casalog.post("==================================");
          casalog.post("Test 15 two overlapping image-fields, mfs, both multi term");
          casalog.post("==================================");

          write_file('out15.txt', 'imagename=mytest1\nnchan=1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559\n')
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',\
                                       field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       outlierfile='out15.txt',\
                                       mode='mfs',\
                                       ntaylorterms=2,mtype='multiterm',restfreq=['1.5GHz'],\
                                       deconvolver='mtmfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:58:39.580 +40.55.55.931",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==14):  ## 2 image-fields, mfs, Overlapping models. Single term
          casalog.post("==================================");
          casalog.post("Test 14 two overlapping image-fields, mfs, single term");
          casalog.post("==================================");

          write_file('out14.txt', 'imagename=mytest1\nnchan=1\nimsize=[200,200]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:59:02.426 +40.51.14.559')
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',\
                                       field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       outlierfile='out14.txt',\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:58:39.580 +40.55.55.931",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==13): ## 1 image-field, cube --- Imaging with various cube parameter specifications
          casalog.post("==================================");
          casalog.post("Test 13 image-field, cube --- Imaging with various cube parameter specifications");
          #casalog.post("==================================");
          ## chan 5 (TOPO)
          qfstart=qa.quantity("1.1GHz")
          #qvstart=qa.quantity("-59958.5km/s")
          # for restf=1.25GHz
          #qvstart=qa.quantity("35975.1km/s")
          # ch10
          qvstart=qa.quantity("16788.4km/s")

          ##mfstart=me.frequency('LSRK',qa.quantity("1.1GHz"))
          mfstart=me.frequency('LSRK',qa.quantity("1.09999GHz"))
          #mvstart=me.radialvelocity('BARY',qa.quantity("-59976.1km/s"))
          mvstart=me.radialvelocity('BARY',qa.quantity("11977.6km/s"))
          #dop = me.todoppler('radio',mfstart,qa.quantity('1.0GHz'))
          mfstart10=me.frequency('LSRK',qa.quantity(" 1.17999GHz"))
          dop = me.todoppler('radio',mfstart10,qa.quantity('1.25GHz'))
          #1chan width
          #qvwidth = qa.quantity("11991.700km/s")
          qvwidth = qa.quantity("4796.7km/s")
          mvwidth = me.radialvelocity('TOPO',qvwidth)

          # restf = 1.25GHz
          # vel range: 59961.1 -  -31174.7 km/s (lsrk/radio)
          #            74952.3 -  -28238.3 km/s (lsrk/optical)

          testList = {
                      0:{'imagename':'Cubetest_chandefstdefwidth','spw':'0','start':0,'width':1,'outframe':'LSRK',
                        'desc':'channel, default start and width, LSRK'},
                      1:{'imagename':'Cubetest_chandefstdefwidthtopo','spw':'0','start':0,'width':1, 'outframe':'TOPO',
                        'desc':'channel, default start and width, TOPO'},
                      2:{'imagename':'Cubetest_chandefstwidth2','spw':'0','start':0,'width':2, 'outframe':'LSRK',
                        'desc':'channel, default start, width=2, LSRK'},
                      3:{'imagename':'Cubetest_chanst5wd1','spw':'0','start':5,'width':1, 'outframe':'LSRK',
                        'desc':'channel, start=5, default width, LSRK'},
                      # this will result in blank channnel images (calcChanFreqs requires start and width in channel
                      # mode to be given in chan index
                      4:{'imagename':'Cubetest_chandefstwd1spwsel','spw':'0:5~19','start':0,'width':1, 'outframe':'LSRK',
                        'desc':'channel, spw=0:5~19, LSRK'},
                      5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'40MHz','outframe':'TOPO',
                        'desc':'frequency, default start, width=\'40MHz\', TOPO'},
                      6:{'imagename':'Cubetest_freqst5defwd','spw':'0','start':'1.1GHz','width':'','outframe':'TOPO',
                        'desc':'frequency, start=\'1.1GHz\', default width, TOPO'},
                      7:{'imagename':'Cubetest_freqst5defwdspwsel','spw':'0:4~19','start':'1.1GHz','width':'','outframe':'TOPO',
                        'desc':'frequency, start=\'1.1GHz\', default width, spw=0:6~19, TOPO'},
                      8:{'imagename':'Cubetest_freqst10wdm','spw':'0','start':'1.2GHz','width':'-20MHz','outframe':'TOPO',
                        'desc':'frequency, start=\'1.2GHz\', width=\'-20MHz\', TOPO'},
                      9:{'imagename':'Cubetest_veldefstwd2','spw':'0','start':'','width':'9593.4km/s','outframe':'TOPO',
                        'desc':'frequency, default start, width=\'9593.4km/s\', TOPO'},
                     10:{'imagename':'Cubetest_veldefstwd2m','spw':'0','start':'','width':'-9593.40km/s','outframe':'TOPO',
                        'desc':'velocity, default start, width=\'-9593.40m/s\', TOPO'},
                     11:{'imagename':'Cubetest_velst10defwd','spw':'0','start':'11991.7km/s','width':'','outframe':'TOPO',
                        'desc':'velocity, start=\'11991.7km/s\', default width, TOPO'},
                     12:{'imagename':'Cubetest_velst10defwdbary','spw':'0','start':'11977.6km/s','width':'','outframe':'BARY',
                        'desc':'velocity, start=\'11977.6km/s\', default width, BARY'},
                     # currently 13 is not quite properly working, investigating - 2014.08.27 TT
                     #13:{'imagename':'Cubetest_optvelst10defwdlsrk','spw':'0','start':'-49962.6km/s','width':'4234.40km/s',
                     ##13:{'imagename':'Cubetest_optvelst10defwdlsrk','spw':'0','start':'-49962.6km/s','width':'5500.0km/s',
                     ##13:{'imagename':'Cubetest_optvelst10defwdlsrk','spw':'0','start':'74952.3km/s','width':'5500.0km/s',
                     ##   'veltype':'optical','outframe':'LSRK',
                     ##   'desc':'velocity, start=\'74952.3km/s\', default width, veltype=optical LSRK'},
                     ### smallest w in vopt = 4792.2km/s (ch19-20), largest w in vopt: 8817km/s (ch0-1)
                     ##13:{'imagename':'Cubetest_optvelst10wdeflsrk','spw':'0','start':'74952.3km/s','width':'',
                     13:{'imagename':'Cubetest_optvelst10wdeflsrk','spw':'0','start':'12494.8km/s','width':'',
                     #13:{'imagename':'Cubetest_optvelst0defwdlsrk','spw':'0','start':'26072.5km/s','width':'8817km/s',
                     #13:{'imagename':'Cubetest_optvelst2defwdlsrk','spw':'0','start':'132605km/s','width':'-8817km/s',
                        'veltype':'optical','outframe':'LSRK',
                     ##   'desc':'velocity, start=\'74952.3km/s\', default width, veltype=optical LSRK'},
                        'desc':'velocity, start=\'12494.8km/s\', default width, veltype=optical LSRK'},
                     14:{'imagename':'Cubetest_stqfreqdefwd','spw':'0','start':qfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start=%s, default width, veltype=radio TOPO' % qfstart},
                     15:{'imagename':'Cubetest_stmfreqdefwd','spw':'0','start':mfstart,'width':'', 'veltype':'radio','outframe':'',
                        'desc':'frequency, start=%s, default width, veltype=radio LSRK' % mfstart},
                     16:{'imagename':'Cubetest_stqveldefwd','spw':'0','start':qvstart,'width':'','outframe':'TOPO',
                        'desc':'velocity, start=%s, default width, TOPO' % qvstart},
                     17:{'imagename':'Cubetest_stmveldefwd','spw':'0','start':mvstart,'width':'','outframe':'TOPO',
                        'desc':'velocity, start=%s, default width, BARY' % mvstart},
                     18:{'imagename':'Cubetest_veldefstqvwidth','spw':'0','start':'','width':qvwidth,'outframe':'TOPO',
                        'desc':'velocity, default start, width=%s, TOPO' % qvwidth},
                     19:{'imagename':'Cubetest_veldefstmvwidth','spw':'0','start':'','width':mvwidth,'outframe':'TOPO',
                        'desc':'velocity, default start, width=%s, TOPO' % mvwidth},
                     # this also fails, investigating - 2014.08.27 TT
                     20:{'imagename':'Cubetest_stdopdefwd','spw':'0','start':dop,'width':'','outframe':'TOPO',
                        'desc':'doppler, start=%s, default width, LSRK' % dop}
                    }

          casalog.post("Sub-test "+str(testid)+":"+testList[testid]['desc']);
          casalog.post("==================================");
          if 'veltype' in testList[testid]:
               inveltype=testList[testid]['veltype']
          else:
               inveltype='radio'
          paramList = ImagerParameters(msname='DataTest/point_twospws.ms', field='0',\
                                       spw=testList[testid]['spw'],\
                                       usescratch=True,readonly=True,\
                                       imagename=testList[testid]['imagename'], mode='cube',\
                                       #nchan=3,\
                                       nchan=10,\
                                       #nchan=20,\
                                       restfreq=['1.25GHz'],\
                                       #restfreq=['1.5GHz'],\
                                       start=testList[testid]['start'],\
                                       width=testList[testid]['width'],\
                                       veltype=inveltype,\
                                       outframe=testList[testid]['outframe'],\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)



     if(testnum==12):  ## 1 image-field, mfs --- WB AWP ( multi term mosaic )
          casalog.post("==================================");
          casalog.post("Test 12 image-field, mfs --- WB AWP( multi term )");
          casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/reg_mawproject.ms',
                                       field='*',scan='',
                                       spw='*',\
                                       usescratch=True,readonly=True,\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[512,512],\
                                       cell=['10.0arcsec','10.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ntaylorterms=2,mtype='multiterm',restfreq=['1.5GHz'],\
                                       ftmachine='awprojectft',
                                       startmodel='', weighting='natural',\
                                       aterm=True, psterm=False, mterm=True,\
                                       wbawp = True,
                                       cfcache = "perm.wb.nt2.mytest0.cfcache",\
                                       dopointing = False, dopbcorr = True, conjbeams = True,
                                       computepastep =360.0, rotatepastep =360.0,#5.0,\
                                       deconvolver='mtmfs',\
                                       pblimit=0.1,normtype='flatnoise',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==11):  ## 1 image-field, mfs --- WB AWP ( single term )
          casalog.post("==================================");
          casalog.post("Test 11 image-field, mfs --- WB AWP (single term)");
          casalog.post("==================================");
          ## Note : For this test, if doMajor peak comes to 0.006 instead of 0.69,
          ##           it means unnormalized weight is being used (refilled perhaps but not normed).
          paramList = ImagerParameters(msname='DataTest/reg_mawproject.ms',
                                       field='1',scan='',
                                       spw='*',\
                                       usescratch=True,readonly=True,\
                                       imagename='mytest0', nchan=1,start='1.5GHz', width='1.0GHz',\
                                       imsize=[512,512],\
                                       cell=['10.0arcsec','10.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='awprojectft',
                                       startmodel='', weighting='natural',\
                                       aterm=True, psterm=False, mterm=True,\
                                       wbawp = True,
                                       cfcache = "perm.wb.mytest0.cfcache",\
                                       dopointing = False, dopbcorr = True, conjbeams = True,
                                       computepastep =360.0, rotatepastep =5.0,\
                                       deconvolver='hogbom',\
                                       pblimit=0.1,normtype='flatnoise',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==10):  ## 1 image-field, mfs --- Narrow-band AWP
          casalog.post("==================================");
          casalog.post("Test 10 image-field, mfs --- WB AWP");
          casalog.post("==================================");
          paramList = ImagerParameters(msname='DataTest/reg_mawproject.ms', #_offcenter.ms',
                                       field='*',scan='',
                                       spw='1',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[512,512],\
                                       cell=['10.0arcsec','10.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       stokes='V',\
                                       ftmachine='awprojectft',
                                       startmodel='', weighting='natural',\
                                       aterm=True, psterm=False, mterm=True,\
                                       wbawp = False,
                                       cfcache = "perm.mytest0.cfcache.mos",#.offcenter",\
                                       dopointing = False, dopbcorr = True, conjbeams = True,
                                       computepastep =360.0, rotatepastep =5.0,\
                                       deconvolver='hogbom',\
                                       #pblimit=0.1,normtype='flatsky',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==9):  ## Image domain mosaic for single-term (narrow band )
          casalog.post("==================================");
          casalog.post("Test 9 : Image domain mosaic for single-term");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/reg_mawproject.ms',\
                                       field='',spw='1',scan='',\
#                                       usescratch=True,readonly=True,\
                                       savemodel='none',\
                                       specmode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       ntaylorterms=1,gridder='imagemosaic',restfreq=['1.5GHz'],\
                                       imsize=[512,512],\
                                       cell=['10.0arcsec','10.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       #phasecenter=1,\
                                       #ftmachine='gridft',
                                       startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       #pblimit=0.1,normtype='flatsky',
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)

     if(testnum==8):  # MTMFS + Facets
          casalog.post("==================================");
          casalog.post("Test 8 : 1 image-field, mfs, nt=1, 2x2 facets ");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       ntaylorterms=2,mtype='multiterm',restfreq=['1.5GHz'],\
                                       imsize=[200,200], facets=2,\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       phasecenter="J2000 19:59:00.2 +40.50.15.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='mtmfs',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==7):  # MTMFS
          casalog.post("==================================");
          casalog.post("Test 7 : 1 image-field, mfs, With ntaylorterms=2 ");
          casalog.post("==================================");

#          paramList = ImagerParameters(msname='DataTest/eptwochan.ms',field='0',spw='*',\
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       ntaylorterms=2,mtype='multiterm',
                                       reffreq='1.6GHz',\
                                       #restfreq=['1.5GHz'],\
                                       imsize=[200,200],\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       #phasecenter="J2000 19:59:00.2 +40.50.15.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='mtmfs',\
#                                       scales=[0,20,40,100],\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==6): ## Facetted imaging
          casalog.post("==================================");
          casalog.post("Test 6 : 1 image-field, mfs, nt=1, 2x2 facets ");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       ntaylorterms=1,mtype='default',restfreq=['1.5GHz'],\
                                       imsize=[200,200], facets=2,\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       phasecenter="J2000 19:59:00.2 +40.50.15.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==5):  ## 1 image-field, mfs, multiple input MSs --- Imaging.
          casalog.post("==================================");
          casalog.post("Test 5 image-field, mfs, multiple input MSs --- Imaging.");
          casalog.post("==================================");
          paramList = ImagerParameters(msname=['DataTest/point_onespw0.ms','DataTest/point_onespw1.ms'],\
                                       field='0',spw=['0','0'],\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==4):  ## 2 image-fields, one cube, one mfs --- Imaging.
          casalog.post("==================================");
          casalog.post("Test 4 image-fields, one cube, one mfs --- Imaging.");
          casalog.post("==================================");

          write_file('out4.txt', 'imagename=mytest1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nmode=mfs\nstart=1.0GHz\nwidth=2.0GHz')
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',\
                                       field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       outlierfile='out4.txt',\
                                       mode='cube',\
                                       imagename='mytest0', nchan=2,\
                                       start='1.0GHz', width='1.0GHz',
                                       #start='1.0GHz', width='2.0GHz',
                                       #start=0, width=1,
                                       outframe='TOPO',\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)



     if(testnum==3):  ## 2 image-fields, mfs --- Imaging.
          casalog.post("==================================");
          casalog.post("Test 3 image-fields, mfs --- Imaging.");
          casalog.post("==================================");

          if len(mask)==0:
               write_file('out3.txt', 'imagename=mytest1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543')
          else:
               write_file('out3.txt', 'imagename=mytest1\nnchan=1\nimsize=[80,80]\ncell=[8.0arcsec,8.0arcsec]\nphasecenter=J2000 19:58:40.895 +40.55.58.543\nmask=circle[[40pix,40pix],10pix]')
          paramList = ImagerParameters(msname='DataTest/twopoints_twochan.ms',\
                                       field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       outlierfile='out3.txt',\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==2):  ## 1 image-field, cube --- Imaging.
          casalog.post("==================================");
          casalog.post("Test 2 image-field, cube --- Imaging.");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/point_twospws.ms', field='0',spw='',\
                                       usescratch=True,readonly=True,\
                                       mode='cube',\
                                       imagename='mytest0', nchan=10,start='1.0GHz', width='40MHz',\
                                       imsize=[100,100],\
                                       cell=['8.0arcsec','8.0arcsec'],\
                                       phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ftmachine='gridft', startmodel='', weighting='natural',\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)


     if(testnum==1):  ## 1 image-field, mfs --- Imaging.
          casalog.post("==================================");
          casalog.post("Test 1 image-field, mfs --- Imaging.");
          casalog.post("==================================");

          paramList = ImagerParameters(msname='DataTest/point_twospws.ms',field='0',spw='0',\
                                       usescratch=True,readonly=True,\
                                       mode='mfs',\
                                       imagename='mytest0', nchan=1,\
                                       imsize=[110,110],\
                                       cell=['8.0arcsec','8.0arcsec'],stokes='I',\
                                       ### center
                                       #phasecenter="J2000 19:59:28.500 +40.44.01.50",\
                                       ### offset
                                       #phasecenter="J2000 19:59:23.591 +40.44.01.50",\
                                       phasecenter=0,
                                       ftmachine='gridft', startmodel='', weighting='briggs',robust=0.5,\
                                       deconvolver='hogbom',\
                                       niter=niter,cycleniter=cycleniter,\
                                       threshold=threshold,loopgain=loopgain,\
                                       restoringbeam=restoringbeam,
                                       interactive=interactive,mask=mask)



     ### Check input parameters, and parse outlier files.
     #if paramList.checkParameters() == False:
     #   return [None, "", False, False,False]

     #paramList.printParameters()

     return [ paramList , parallelmajor, parallelminor, parallelcube ]

#####################################################


#######################################
###    Function to run the Imaging.
#######################################
def doClean( params = [None,False,False,False] , doplot=True ):

    os.system('rm -rf mytest*')

    pmajor = params[1]
    pcube = params[3]

    if pmajor==False and pcube==False:
         imager = PySynthesisImager(params[0])
    elif pmajor==True:
         imager = PyParallelContSynthesisImager(params[0])
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params[0])
    else:
         print('Invalid parallel combination in doClean.')
         return

    ### Set up Imagers, Deconvolvers, IterControl, and ParallelSync.
    imager.initializeImagers()
    imager.initializeDeconvolvers()
    imager.initializeNormalizers()
    imager.initializeIterationControl()

    ### Run it.
    imager.makePSF()
    imager.runMajorMinorLoops()

    imager.restoreImages()

    if( doplot == True ):
         imager.getSummary();

    imager.deleteTools()
########################################


########################################
########################################
#  Run only Major Cycle
########################################
def doMajor( params = [None,False,False,False] , doplot=True , tomake='both'):

    os.system('rm -rf mytest*')

    pmajor=params[1]
    pcube = params[3]

    if pmajor==False and pcube==False:
         imager = PySynthesisImager(params[0])
    elif pmajor==True:
         imager = PyParallelContSynthesisImager(params[0])
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params[0])
    else:
         print('Invalid parallel combination in doClean.')
         return


    ### Set up Imagers and ParallelSync.
    imager.initializeImagers()
    imager.initializeNormalizers()

    ### Run it.
    if tomake=='both' or tomake=='psf':
         imager.makePSF()
    if tomake=='both':
         imager.runMajorCycle()

    imager.deleteTools()

########################################
########################################
###   Run only the minor cycle....
########################################
def doMinor( params = [None,False,False,False] , doplot=True ):

##    os.system('rm -rf mytest*')

    pminor=params[2]

    if pminor==False:
         imager = PySynthesisImager(params[0])
    else:
         imager = PyParallelDeconvolver(params[0])

    ### Set up Deconvolvers and IterControl
    imager.initializeDeconvolvers()
    imager.initializeIterationControl()

    # Run it
    while ( not imager.hasConverged() ):
        imager.runMinorCycle()

    imager.restoreImages()

    if( doplot == True ):
         imager.getSummary();

    imager.deleteTools()
########################################

########################################
#
#   Predict Model only (for use by setjy too)
#
########################################
########################################
#
#   Restore model image, Inputs : model, residual, psf.
#
########################################
def doRestore( params = [None,False,False,False] ):

    pminor=params[2]
    imager = PySynthesisImager(params[0])
    imager.initializeDeconvolvers()
    imager.restoreImages()
    imager.deleteTools()

########################################
def doSetjy( params = [None,False,False,False] ):

    pminor=params[2]
    imager = PySynthesisImager(params[0])
    imager.initializeDeconvolvers()

    imager.predictModel()

    imager.deleteTools()


########################################

def write_file(filename,str_text):
    '''Save the string in a text file'''
    inp = filename
    cmd = str_text
    # remove file first
    if os.path.exists(inp):
        os.system('rm -f '+ inp)
    # save to a file
    with open(inp, 'w') as f:
        f.write(cmd)
    f.close()
    return


########################################

## Test function for live ImageStore transport between tools.
def toolTestMajorCycle( testnum=1 ):
     params = getparams( testnum=testnum )[0]

     SItool = casac.synthesisimager()
     SItool.initializemajorcycle( params.getSelPars(), params.getImagePars(), params.getGridPars() )

     PStool = casac.synthesisnormalizer()
     normpars = {'imagename':params.getImagePars()['0']['imagename']}
     PStool.setupnormalizer( normpars )

     PStool.scattermodel()
     SItool.makepsf()
     SItool.executemajorcycle()

     PStool.gatherpsfweight( )
     PStool.gatherresidual( )

     PStool.dividepsfbyweight()
     PStool.divideresidualbyweight()

     SItool.done()
     PStool.done()


###################################################

## Test function for parameter-list interface for the major cycle.
def toolTestMajorCycle2( testnum=1 ):

     os.system('rm -rf mytest*')

     params = getparams( testnum=testnum )[0]
     allselpars = params.getSelPars()
     allimagepars = params.getImagePars()
     allgridpars = params.getGridPars()

     SItool = casac.synthesisimager()

     for mss in sorted(allselpars.keys()):
          SItool.selectdata( **(allselpars[ mss ] ) )

     for fld in sorted(allimagepars.keys()):
          print('Setup imaging and image for field ' + str(fld))
          SItool.defineimage( **(allimagepars[fld]) )

     SItool.setweighting( **allgridpars )

     PStools = []
     nfld = len( list(allimagepars.keys()) )
     for fld in range(0, nfld ):
          PStool.append(casac.synthesisnormalizer())
          normpars = {'imagename':allimagepars[ (list(allimagepars.keys()))[fld] ]['imagename']}
          PStool.setupnormalizer( normpars )

     for fld in range(0, nfld):
          PStool[fld].scattermodel()

     SItool.makepsf()
     SItool.executemajorcycle()

     for fld in range(0, nfld):
          PStool[fld].gatherpsfweight( )
          PStool[fld].gatherresidual( )

          PStool[fld].dividepsfbyweight()
          PStool[fld].divideresidualbyweight()

     SItool.done()

     for fld in range(0, nfld):
          PStool[fld].done()


###################################################


def checkDataPartitioningCode():

     ## Make parameter lists.
     #paramList = ImagerParameters(msname=['x1.ms','x2.ms'], field='0',spw=['0','2'], usescratch=True)
     ## Sync input lists to the same size.
     #paramList.checkParameters()

     params = getparams( testnum=2 ,parallelmajor=True )
     paramList = params[0]

     selpars = paramList.getSelPars()
     impars = paramList.getImagePars()

     synu = casac.synthesisutils()
     print(synu.contdatapartition( selpars , 2))

     print(synu.cubedatapartition( selpars, 2))

     print(synu.cubeimagepartition( impars, 2))

     synu.done()

#     ppar = PyParallelImagerHelper()
#
#     print 'Selpars : ', selpars
#     newselpars = ppar.partitionCubeDataSelection( selpars )
#
#     print 'Impars : ', impars
#     newimpars = ppar.partitionCubeDeconvolution( impars )
#
#     ppar.takedownCluster()
#

     # The output dictionary should be indexed as follows ( for 2 nodes, and 2 MSs )
     #{ '0' : { 'ms0' : { 'msname':xxx1, 'spw':yyy1 } ,
     #            'ms1' : { 'msname':xxx2, 'spw':yyy1 }
     #         }
     #  '1' : { 'ms0' : { 'msname':xxx1, 'spw':yyy2 } ,
     #            'ms1' : { 'msname':xxx2, 'spw':yyy2 }
     #        }
     # }
     #


def checkPars():

     params = getparams( testnum=7 )
     paramList = params[0]

     selpars = paramList.getSelPars()
     impars = paramList.getImagePars()

     synu = casac.synthesisutils()

#     print selpars['ms0']
#     fixrec = synu.checkselectionparams( selpars['ms0'] )
#     print fixrec
#     print synu.checkselectionparams( fixrec )


     print(impars['0'])
     fixrec = synu.checkimageparams( impars['0'] )
     print("---------------------------------------------")
     print(fixrec)
     print("---------------------------------------------")
     print(synu.checkimageparams( fixrec ))

     synu.done()

def testImageCoordinates( testnum=1, testid=0):

    multitest=False
    if testnum==13: multitest=True
    if multitest and testid==-1: testid=list(range(21))
    if type(testid)==int: testid=[testid]
    for tst in testid:
      params = getparams( testnum, testid=tst )
      paramList = params[0]
      impars=paramList.getImagePars()['0']
      if os.path.exists(impars['imagename']+'.residual'):
        os.system('rm -rf '+impars['imagename']+'*')

      doMajor(params)

      ok = ia.open(impars['imagename']+'.residual')
      if ok:
        selpars = paramList.getSelPars()
        outcsys=ia.coordsys().torecord()
        imsummary=ia.summary()
        ia.done()
        #print selpars
        nchan = impars['nchan']
        if (nchan!=imsummary['shape'][3]):
          print("No. of image channel =%s while nchan asked is %s" % (imsummary['shape'][3], nchan))
        checkimcoord(outcsys,selpars,impars)

def checkimcoord(csys,selpars,impars):
    """
    check the output image from testImageCoordinates if it is consistent
    with the input data and image parameters (selpars, impars)
    """

    errs=[]
    # freq frame dict
    specframes={0:'REST',1:'LSRK',2:'LSRD',3:'BARY',4:'GEO',5:'TOPO',6:'GALACTO',7:'LGROUP',8:'CMB',64:'Undefined'}

    msname=selpars['ms0']['msname']
    selspw=selpars['ms0']['spw']
    selfield=selpars['ms0']['field']
    if selspw=='': selspw='*'
    ms.open(msname)
    sel={}
    sel['spw']=selspw
    sel['field']=selfield
    ms.msselect(sel)
    selected = ms.msselectedindices()
    #print "selected=",selected
    tb.open(msname+'/SPECTRAL_WINDOW')
    chanfreqs=tb.getcell('CHAN_FREQ',selected['spw'])
    chanw=tb.getcell('CHAN_WIDTH',selected['spw'])
    msframe=tb.getcell('MEAS_FREQ_REF',selected['spw'][0])
    tb.close()
    dataframe=specframes[msframe]
    #print "dataframe=",dataframe
    chansel = selected['channel'][0][1]
    chanselend = selected['channel'][0][2]
    #print "start freq=",chanfreqs[chansel]
    #print "csys=",csys
    # kower edge freq
    descendingfreq = False
    if chanw[0] < 0:
      descendingfreq = True

    floedge = chanfreqs[chansel]-chanw[chansel]/2.0
    fhiedge = chanfreqs[chanselend]+chanw[chansel]/2.0
    #print "floedge=",floedge, "fhiedge=",fhiedge
    outframe=impars['outframe']
    veltype=impars['veltype']
    start=impars['start']
    width=impars['width']
    restf=impars['restfreq'][0]

    reversechanorder=False
    # default start handling with width specification
    # - also depends on freq ordering in the data
    if type(start)==str and start=="":
      widthval=0
      if type(width)==dict:
        if 'type' in width and width['type']=='radialvelocity':
          widthval = width['m0']['value']
        elif 'unit' in width and width['unit'].count('m/s'):
          widthval = width['value']
        if (descendingfreq and widthval < 0) or (not descendingfreq and widthval > 0):
          reversechanorder=True
      elif width.count('m/s')!=0:
        widthval = qa.quantity(width)['value']
        if (descendingfreq and widthval < 0) or (not descendingfreq and widthval > 0):
          reversechanorder=True

    #print "descendingfreq=",descendingfreq
    #print "reversechanorder=",reversechanorder
    needfreqconv=False
    if dataframe!=outframe:
       needfreqconv=True

    # from output image --------------------------------
    # freqframe
    if 'spectral2' in csys:
      specCoord = csys['spectral2']
      if 'tabular' in specCoord:
        imch0=specCoord['tabular']['crval'][0]
        iminc=specCoord['tabular']['cdelt'][0]
      else:
        imch0=specCoord['wcs']['crval']
        iminc=specCoord['wcs']['cdelt']

    imframe=specCoord['system']
    # conversion layer- this should be LSRK by default
    convframe=specCoord['conversion']['system']
    obsdate=csys['obsdate']
    # if image channel reversed, check last channel
    if reversechanorder:
      imch0=imch0+iminc*(impars['nchan']-1)
    # check frames
    if convframe!="LSRK":
      errs.append('Error: conversion layer is not set to LSRK!')
    if imframe!=outframe:
      errs.append('Error: inconsistent output spec frame:%s, should be %s' % (imframe,outframe))

    startfreq = 0
    imch0conv = imch0
    vmode = False
    if type(start)==int:
      if start>=chansel:
        if type(width)==int and width > 0:
          range = start + width - 1
          startfreq = (chanfreqs[range] + chanfreqs[start]) /2.0
        else:
          startfreq=chanfreqs[start]
      else:
        # start chan idx less than chansel by spw. Will results in some blank channel images
        startfreq = chanfreqs[start]
    else:
      if type(start)==str:
      # start parameter
        if start=="":
          freq0q = qa.quantity(str(chanfreqs[chansel])+'Hz')
        elif start.count('Hz') !=0:
          freq0q = qa.quantity(start)
        elif start.count('m/s') !=0:
          freq0q = vfconv(start,outframe,restf,veltype)['m0']
          vmode = True
      elif type(start)==dict:
        # measure freq/vel format
        if 'type' in start:
          mframe = start['refer']
          if start['type']=='frequency':
            freq0q = start['m0']
          else:
            freq0q = vfconv(str(start['m0']['value'])+start['m0']['unit'], mframe, restf, veltype)['m0']
            vmode = True
        else:
          if start['unit'].count('Hz') != 0:
            freq0q = start
          elif start['unit'].count('m/s') != 0:
            freq0q = vfconv(str(start['value'])+start['unit'],outframe,restf,veltype)['m0']
            vmode = True

      # step(width)
      qhalfw = 0
      if type(width)==str and width!="":
        if width.count('Hz') != 0:
          qhalfw = qa.div(qa.quantity(width),2)
        elif width.count('m/s') != 0:
          qhalfw = qa.div(qa.quantity(width),2)
          vmode = True

      if qhalfw:
        # add width part to the starting frequency
        if vmode:
          qvst = qa.quantity(start)
          if qvst['value']==0.0:
            qvst = vfconv(str(freq0q['value'])+freq0q['unit'],outframe, restf, veltype)['m0']
          qvstedge = qa.sub(qvst,qhalfw)
          qfstedge = vfconv(str(qvstedge['value'])+qvstedge['unit'], outframe, restf, veltype)['m0']
          fstedge = qa.convert(qfstedge, 'Hz')['value']

          if reversechanorder:
            if (fstedge > floedge and not descendingfreq) or (fstedge < floedge and descendingfreq):
              qvelcen = qa.sub(vfconv(str(floedge)+'Hz',outframe, restf,veltype)['m0'],qhalfw)
            else:
              qvelcen = qvst
          else:
            if (fstedge < floedge and not descendingfreq) or (fstedge > floedge and descendingfreq):
              # shift start to fit into data....
              #qvelcen = qa.add(qa.quantity(start),qhalfw)
              qvelcen = qa.add(vfconv(str(floedge)+'Hz',outframe, restf,veltype)['m0'],qhalfw)
            else:
              qvelcen = qvst
          qstfreq = vfconv(str(qvelcen['value'])+qvelcen['unit'], outframe, restf, veltype)['m0']
          startfreq = qa.convert(qstfreq,'Hz')['value']
        else:
          qfstedge = qa.sub(freq0q,qhalfw)
          fstedge = qa.convert(qfstedge, 'Hz')['value']
          if (fstedge < floedge and not descendingfreq) or (fstedge > floedge and descendingfreq):
            qfrqcen = qa.add(qa.quantity(floedge,'Hz'),qhalfw)
            qstfreq = qa.convert(qfrqcen,'Hz')
          else:
            qstfreq = freq0q
          startfreq = qa.convert(qstfreq,'Hz')['value']
      else:
        startfreq = qa.convert(freq0q,'Hz')['value']

    #print "startfreq=",startfreq
    if needfreqconv:
      #print "converting image ch0 freq to freq in dataframe=",dataframe
      me.doframe(me.epoch(obsdate['refer'],str(obsdate['m0']['value'])+obsdate['m0']['unit']))
      me.doframe(csys['telescopeposition'])
      me.doframe(me.direction(csys['direction0']['conversionSystem'],str(csys['direction0']['crval'][0])+'rad',str(csys['direction0']['crval'][1])+'rad'))
      print("******************\n")
      #print me.showframe()
      #print "imch0=", imch0
      imch0conv=me.measure(me.frequency(imframe,str(imch0)+'Hz'),dataframe)['m0']['value']
    #print "dataframe=",dataframe, " outframe=",outframe, " imframe=",imframe
    if (abs(imch0conv-startfreq)/startfreq < 0.0001):
      if (reversechanorder):
        print("OK, image last chan freq match with the data selections and start specification: start freq=", startfreq)
      else:
        print("OK, image ch0 freq match with the data selections and start specification: start freq=", startfreq)
    else:
      print("Error, image ch0 freq does not match with the data selections+start specification, imch0(in %s)=%s start freq=%s" % (dataframe,imch0conv,startfreq))


def vfconv(velorfreqstr, frame, restfstr, veltype):
    """
    vel to freq or freq to vel convesion
    velorfreq - input vel of freq in a string with unit
    """
    if (qa.quantity(velorfreqstr)['unit'].find('m/s') > -1 ):
      dop = me.doppler(veltype, qa.quantity(velorfreqstr))
      retvf = me.tofrequency(frame, dop, qa.quantity(restfstr))
    elif (qa.quantity(velorfreqstr)['unit'].find('Hz') > -1 ):
      mfreq = me.frequency(frame, qa.quantity(velorfreqstr))
      retvf = me.todoppler(veltype, mfreq, qa.quantity(restfstr))
    else:
      print("Error! cannot convert ", velorfreqstr)
    return retvf


#def testmakeimage():
#     params = getparams( testnum=1 )
#     paramList = params[0]
#     impars = (paramList.getImagePars())['0']
#
#     os.system('rm -rf ' + impars['imagename'])
#
#     synu = casac.synthesisutils()
#     #synu.makeimage( impars , 'DataTest/twopoints_twochan.ms')
#     synu.makeimage( impars , 'DataTest/point_twospws.ms')
#     synu.done()


def delmodkeywords(msname=""):
     delmod(msname)
     tb.open( msname+'/SOURCE', nomodify=False )
     keys = tb.getkeywords()
     for key in keys:
          tb.removekeyword( key )
     tb.close()

def resetmodelcol(msname=""):
     tb.open( msname, nomodify=False )
     dat = tb.getcol('MODEL_DATA')
     dat.fill( complex(1.0,0.0) )
     tb.putcol('MODEL_DATA', dat)
     tb.close();

def delmodels(msname=""):
     resetmodelcol(msname)  ## Set model column to zero
     delmod(msname)  ## Get rid of OTF model
     delmodkeywords(msname) ## Get rid of extra OTF model keywords that sometimes persist...

def checkmodel(msname=""):
     tb.open( msname )
     hasmodcol = (  (tb.colnames()).count('MODEL_DATA')>0 )
     modsum=0.0
     if hasmodcol:
          dat = tb.getcol('MODEL_DATA')
          modsum=dat.sum()
     tb.close()
     tb.open( msname+'/SOURCE' )
     keys = tb.getkeywords()
     if len(keys)>0:
          hasvirmod=True
     else:
          hasvirmod=False
     tb.close()
     print(msname , ": modelcol=", hasmodcol, " modsum=", modsum, " virmod=", hasvirmod)
