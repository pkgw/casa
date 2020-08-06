###############################################
#                                             #
# Regression/Benchmarking Script for NGC 7538 #
#                                             #
###############################################

import time
import os

os.system('rm -rf ngc7538d* 1328.* 2229.cont2* ap314.* ngc7538*.ms')

startTime = time.time()
startProc = time.clock()

datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/ATST1/NGC7538/'

print('--Import--')
default('importvla')
importvla(archivefiles=[datapath+'AP314_A950519.xp1',
                        datapath+'AP314_A950519.xp2',
                        datapath+'AP314_A950519.xp3'],
	  vis='ngc7538.ms', bandname='K', frequencytol=10000000.0)
importtime = time.time() 
print('--Observation summary--')
listobs(vis='ngc7538.ms')
#listtime = time.time()
print('--Flag auto-correlations--')
default('flagdata')
flagdata(vis='ngc7538.ms', mode='manual', autocorr=True)
flagtime = time.time()
print('--Setjy--')
default('setjy')
setjy(vis='ngc7538.ms',field='0',standard='Perley-Taylor 99',scalebychan=False) #set flux density for 1331+305 (3C286)
setjytime = time.time()

print('--Gencal(opac)--')
default('gencal')
gencal(vis='ngc7538.ms', caltable='ap314.opac',
       caltype='opac',parameter=[0.08])
gencaltime = time.time()

print('--Gaincal--')
default('gaincal')
gaincal(vis='ngc7538.ms', caltable='ap314.gcal',
	field='<2', spw='0~1:2~56', gaintype='G',
	solint='inf', combine='', refant='VA19',
	gaintable=['ap314.opac'])
gaintime = time.time()
print('--Bandpass--')
default('bandpass')
bandpass(vis='ngc7538.ms', caltable='1328.bcal',
	 field='0',
	 gaintable=['ap314.opac','ap314.gcal'], interp=['','nearest'],
	 refant='VA19')
bptime = time.time()
print('--Fluxscale--')
default('fluxscale')
fluxscale(vis='ngc7538.ms', caltable='ap314.gcal', fluxtable='ap314.fluxcal',
	  reference=['1328+307'], transfer=['2229+695'])
fstime = time.time()
print('--Apply Cal--')
default('applycal')
applycal(vis='ngc7538.ms',
	 field='1~5',
	 gaintable=['ap314.opac','ap314.fluxcal', '1328.bcal'],
	 gainfield=['','1'])
	 
correcttime = time.time()

print('--Split (fluxcal data)--')
default('split')
split(vis='ngc7538.ms', outputvis='ngc7538_cal.split.ms',
      field='0',spw='0:0~61', datacolumn='model')

print('--Split (continuum)--')
default('split')
# This _averages_ 
split(vis='ngc7538.ms', outputvis='ngc7538d.cont.ms',
      field='3',spw='0:2~56',width=[55], datacolumn='corrected')

print('--Split (mf cont,)--')
default('split')
split(vis='ngc7538.ms', outputvis='ngc7538.cont.ms',
	field='3,4,5',spw='0:2~56',width=[55], datacolumn='corrected')
print('--Split (bandcal data)--')
default('split')
split(vis='ngc7538.ms', outputvis='2229.cont2.ms',
      field='1', spw='0:2~56,1:2~56',width=[55,55], datacolumn='corrected')
splitcaltime = time.time()
default('split')
split(vis='ngc7538.ms',outputvis='ngc7538d.line.ms',
      field='3',spw='0:2~56,1:2~56',datacolumn='corrected')
splitsrctime = time.time()

print('--Clean cal--')
default('clean')
clean(vis='2229.cont2.ms',imagename='2229.cont2',mode='channel',
      psfmode='hogbom',niter=6000,gain=0.1,threshold=8.,mask='',
      nchan=1,start=0,width=1,field='0',spw='0',
      imsize=[256,256],cell=[0.5,0.5],
      weighting='briggs',robust=0.5,imagermode='')
cleantime1 = time.time()
print('--Clean src cont--')
default('clean')
clean(vis='ngc7538d.cont.ms',imagename='ngc7538d.cont',mode='channel',
      psfmode='hogbom',niter=5000,gain=0.1,threshold=3.,mask='',
      nchan=1,start=0,width=1,field='0',spw='0',
      imsize=[1024,1024],cell=[0.5,0.5],
      weighting='briggs',robust=0.5,imagermode='')
cleantime2 = time.time()
print('--Clean src line--')
default('clean')
clean(vis='ngc7538d.line.ms',imagename='ngc7538d.cube',mode='channel',
      psfmode='hogbom',niter=5000,gain=0.1,threshold=30.,mask='',
      nchan=48,start=2,width=1,field='0',spw='0',
      imsize=[128,128],cell=[4.,4.],
      weighting='briggs',robust=2.,imagermode='',
      uvtaper=True, outertaper=['12.0arcsec','12.0arcsec', '0deg'])
cleantime3 = time.time()
# -- Not done in old regression but should be
#print '--Contsub (image plane)--'
#ia.open('ngc7538d.cube.image')
#myim=ia.continuumsub('ngc7538d_subed.line.im','ngc7538d_res.cont.im',channels=range(0,48),fitorder=1)
#ia.close()
#myim.close()
#contsubtime=time.time()
#print '--View image--'
#viewer('ngc5921_task.image')

endProc = time.clock()
endTime = time.time()

# Regression

test_name_cal = 'NGC7538--Calibrater maximum amplitude test'
test_name_src = 'NGC7538--Source maximum amplitude test'
test_name_immax = 'NGC7538--Image maximum test'
test_name_imrms = 'NGC7538--Image rms test'

ms.open('ngc7538_cal.split.ms')
thistest_cal=max(ms.range(["amplitude"]).get('amplitude'))
ms.close()
ms.open('ngc7538d.line.ms')
thistest_src=max(ms.range(["amplitude"]).get('amplitude'))
ms.close()
ia.open('ngc7538d.cube.image')
# ia.statistics returns dictionary with 'return','statsout'
# get the second value in the dictionary (statsout)
statistics=ia.statistics(list=True, verbose=True)
ia.close()
# note thistest_immax will be a list with one value 
thistest_immax=statistics['max'][0]
# note thistest_imrms will be a list with one value 
thistest_imrms=statistics['rms'][0]

cal_max=2.413
src_max=18.3638
im_max=0.2606
im_rms=0.0124


diff_cal=abs((cal_max-thistest_cal)/cal_max)
diff_src=abs((src_max-thistest_src)/src_max)
diff_immax=abs((im_max-thistest_immax)/im_max)
diff_imrms=abs((im_rms-thistest_imrms)/im_rms)

import datetime
datestring=datetime.datetime.isoformat(datetime.datetime.today())
outfile='ngc7538.'+datestring+'.log'
logfile=open(outfile,'w')

print('', file=logfile)
print('********** Data Summary *********', file=logfile)
print('* Observer: unavailable     Project: AP314                                  *', file=logfile)
print('* Observation: VLA(27 antennas)                                             *', file=logfile)
print('*  Telescope Observation Date Observer       Project                        *', file=logfile)
print('*  VLA       [                4.30759e+09, 4.30759e+09]unavailable    AP314 *', file=logfile)
print('*  VLA       [                4.30759e+09, 4.30762e+09]unavailable    AP314 *', file=logfile)
print('*  VLA       [                4.30762e+09, 4.30763e+09]unavailable    AP314 *', file=logfile)
print('*Data records: 838404       Total integration time = 36000 seconds          *', file=logfile)
print('*   Observed from   09:23:45   to   19:23:45                                *', file=logfile)
print('*Fields: 6                                                                  *', file=logfile)
print('*  ID   Name          Right Ascension  Declination   Epoch                  *', file=logfile)
print('*  0    1328+307      13:31:08.29      +30.30.33.04  J2000                  *', file=logfile)
print('*  1    2229+695      22:30:36.48      +69.46.28.00  J2000                  *', file=logfile)
print('*  2    NGC7538C      23:14:02.48      +61.27.14.86  J2000                  *', file=logfile)
print('*  3    NGC7538D      23:13:43.82      +61.27.00.18  J2000                  *', file=logfile)
print('*  4    NGC7538E      23:13:34.64      +61.27.26.44  J2000                  *', file=logfile)
print('*  5    NGC7538F      23:13:35.76      +61.28.33.66  J2000                  *', file=logfile)
print('* Data descriptions: 2 (2 spectral windows and 2 polarization setups)       *', file=logfile)
print('*   ID  #Chans Frame Ch1(MHz)    Resoln(kHz) TotBW(kHz)  Ref(MHz)    Corrs  *', file=logfile)
print('*   0       63 TOPO  23691.4682  118.164062  6152.34375  23694.4955  RR     *', file=logfile)
print('*   1       63 TOPO  23719.6063  118.164062  6152.34375  23722.6336  LL     *', file=logfile)
print('*********************************', file=logfile)
print('', file=logfile)
print('********** Regression ***********', file=logfile)
print('*                               *', file=logfile)
passfail = {True: '* Passed',
            False: '* FAILED'}
print(passfail[diff_cal < 0.05], 'cal max amplitude test *', file=logfile)
print('* Cal max amp '+str(thistest_cal), file=logfile)
print(passfail[diff_src < 0.05], 'src max amplitude test *', file=logfile)
print('* Src max amp '+str(thistest_src), file=logfile)
print(passfail[diff_immax < 0.05], 'image max test         *', file=logfile)
print('* Image max '+str(thistest_immax), file=logfile)
print(passfail[diff_imrms < 0.05], 'image rms test         *', file=logfile)
print('* Image rms '+str(thistest_imrms), file=logfile)
if ((diff_src<0.05) & (diff_cal<0.05) & (diff_immax<0.05) & (diff_imrms<0.05)): 
	regstate=True
	print('---', file=logfile)
	print('Passed Regression test for NGC7538', file=logfile)
	print('---', file=logfile)
else: 
	regstate=False
	print('----FAILED Regression test for NGC7538', file=logfile)
print('*********************************', file=logfile)

print('', file=logfile)
print('', file=logfile)
print('********* Benchmarking *****************', file=logfile)
print('*                                      *', file=logfile)
print('Total wall clock time was: '+str(endTime - startTime), file=logfile)
print('Total CPU        time was: '+str(endProc - startProc), file=logfile)
print('Processing rate MB/s  was: '+str(240.3/(endTime - startTime)), file=logfile)
print('* Breakdown:                           *', file=logfile)
print('*   import       time was: '+str(importtime-startTime), file=logfile)
print('*   flagautocorr time was: '+str(flagtime-importtime), file=logfile)
print('*   setjy        time was: '+str(setjytime-flagtime), file=logfile)
print('*   gencal       time was: '+str(gencaltime-setjytime), file=logfile)
print('*   gaincal      time was: '+str(gaintime-gencaltime), file=logfile)
print('*   bandpass     time was: '+str(bptime-gaintime), file=logfile)
print('*   fluxscale    time was: '+str(fstime-bptime), file=logfile)
print('*   applycal     time was: '+str(correcttime-fstime), file=logfile)
print('*   split-cal    time was: '+str(splitcaltime-correcttime), file=logfile)
print('*   split-src    time was: '+str(splitsrctime-splitcaltime), file=logfile)
print('*   clean-cal    time was: '+str(cleantime1-splitsrctime), file=logfile)
print('*   clean-src-c  time was: '+str(cleantime2-cleantime1), file=logfile)
print('*   clean-src-l  time was: '+str(cleantime3-cleantime2), file=logfile)
#print '*   contsub      time was: ',contsubtime-cleantime3
print('*****************************************', file=logfile)
print('basho (test cpu) time was: 500 seconds', file=logfile)

logfile.close()

