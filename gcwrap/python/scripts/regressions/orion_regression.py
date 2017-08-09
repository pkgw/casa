###############################################
#                                             #
# Regression/Benchmarking Script for Orion    #
#        Single Dish + Synthesis              #
###############################################

import time
import os
import shutil

os.system('rm -rf orion_t* gbt_gau.im')
os.system('rm -rf orion.ms')
os.system('rm -rf orion.t*')

pathname=os.environ.get('CASAPATH')
datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/ATST3/Orion/'
#datapath='/home/casa/data/trunk/regression/ATST3/Orion/'

print('--Copy data to local directory--')
mspath='cp -r '+datapath+'orion.ms .'
os.system(mspath)
os.system('chmod -R a+wx orion.ms')

startTime = time.time()
startProc = time.clock()


print('--Feather--')
# Starting from:
#    VLA Orion mosaic image : orion_vlamem.im
#    GBT OTF image : orion.gbt.im
feather('orion_tfeather.im',datapath+'orion_vlamem.im',datapath+'orion.gbt.im')
feathertime = time.time()
#GBT:   Max:5.129806e+01        Flux:2.425065e+02 Jy    rms:1.277546e+01
#VLA:   Max:8.340111e-01        Flux:1.891523e+02 Jy    rms:1.099514e-01

print('--Feather - create synth image--')
# Starting from:
#    VLA Orion mosaic image : orion.ms
#    GBT OTF image : orion.gbt.im
#    Some details about the mosaic:
#    primary beam at X-band = 5.4',
#    mosaic field spacing = 2.5'
#    total mosaic size = approx. 9.5' = 570"
#    synthesized beam size = 8.4" in D config at 3.6 cm, 8.3 GHz
#    cell size = 2" and nx,ny = 300 (600" field size)
#    phase center = center field: J2000 05:35:17.470, -005.23.06.790
#    NOTE: field 10 is outside of the 9 point primary mosaic (sitting
#     on M43 -- but the flux is resolved out so there is no use to
#     add it to the mosaic.  The script below leaves it out.
shutil.rmtree('orion.task.model', True)
default('mosaic')
mosaic('orion.ms',
       'orion.task',
       'mfs',
       'entropy',
       niter=29,
       sigma='4mJy',
       targetflux='240Jy',
       mask=datapath+'orion.mask6',
       field=[2, 3, 4, 5, 6, 7, 8, 9,10],
       spw=[0, 1],
       stokes='I',
       cell=[2, 2],
       imsize=[300, 300],
       weighting='briggs',
       mosweight=True,
       rmode='norm',
       robust=-1,
       cyclefactor=5,
       cyclespeedup=500,
       prior='',
       phasecenter=6,
       ftmachine='ft',
       minpb=0.1,
       scaletype='PBCOR')

feather('orion_tfeather2.im','orion.task.image',datapath+'orion.gbt.im')
#GBT:   Max:5.129806e+01        Flux:2.425065e+02 Jy    rms:1.277546e+01
#VLA:   Max:8.340111e-01        Flux:1.891523e+02 Jy    rms:1.099514e-01
feathersynthtime = time.time()

print('--Single Dish as Model (multi-scale)--')
## Starting from:
##    VLA calibrated visibilities: orion.ms
##    GBT OTF cube: orion.gbt.im
shutil.rmtree('orion_tsdms.model', True)
shutil.rmtree('orion_tsdms.mask', True)
default('clean')
clean(vis='orion.ms',
      imagename='orion_tsdms',
      field='2~10',
      spw='0,1',
      mode='mfs',
      niter=5000,
      gain=0.4,
      threshold='10.0mJy',
      imagermode='mosaic',
      mosweight=True,
      ftmachine='ft',
      cyclefactor=4,
      cyclespeedup=300,
      multiscale=[0,2,8,16,32,64],
      negcomponent=-1,
      flatnoise=False,
      mask='',  #datapath+'orion.mask6',
      modelimage=datapath+'orion.gbt.im',
      imsize=[300,300],
      cell=['2.0arcsec','2.0arcsec'],
      phasecenter=6,
      stokes='I',
      weighting='briggs',
      robust=-1.0,
      pbcor=True,
      minpb=0.1)
sdmodelmstime = time.time()
###combo: Max:1.195286e+00        Flux:2.873779e+02 Jy    rms:9.069330e-02
###GBT:   Max:5.129806e+01        Flux:2.425065e+02 Jy    rms:1.277546e+01
###VLA:   Max:8.340111e-01        Flux:1.891523e+02 Jy    rms:1.099514e-01
##
print('--Single Dish as Model (MEM)--')
### Starting from:
###    VLA calibrated visibilities: orion.ms
###    GBT OTF cube: orion.gbt.im
shutil.rmtree('orion_tsdmem.model', True)
default('mosaic')
mosaic('orion.ms', 'orion_tsdmem', 'mfs', 'entropy', niter=3, sigma='4mJy',
       targetflux='240Jy', mask=datapath+'orion.mask6',
       field=[2, 3, 4, 5, 6, 7, 8, 9, 10], spw=[0, 1], stokes='I',
       cell=[2, 2], imsize=[300, 300],
       weighting='briggs', mosweight=True, rmode='norm', robust=-1,
       cyclefactor=4, cyclespeedup=500,
       phasecenter=6,
       modelimage='orion_tsdmem', sdimage=datapath+'orion.gbt.im',
       ftmachine='ft',
       prior='orion_tsdmem',
       minpb=0.1, scaletype='PBCOR')
sdmodelmemtime=time.time()
###GBT:   Max:5.129806e+01        Flux:2.425065e+02 Jy    rms:1.277546e+01
###VLA:   Max:8.340111e-01        Flux:1.891523e+02 Jy    rms:1.099514e-01
##
######
#### Create synthesis (BIMA) data cube, deconvolve single dish cube
#### DO joint deconvolution
######
def joint_deconvolve(datapath):
        print('--Joint deconvolution --')

        #Regrid GBT image onto synth imaging coordinates
        ia.open('orion_tsdms.image')
        csys = ia.coordsys()
        ia.close()
        ia.open(datapath+'orion.gbt.im')
        ib=ia.regrid(outfile='orion_tgbt_regrid.im',shape=[300,300,1,1], axes=[0,1],
                  csys=csys.torecord(),overwrite=True, asvelocity=False)
        ib.setcoordsys(csys.torecord())
        ia.close()
        ib.done()

        #Deconvolve GBT image
        # Sigh.  dc.open will warn about the lack of a PSF, but I can't seem to
        # define a PSF before calling dc.open.
        dc.open('orion_tgbt_regrid.im', psf='', warn=False)
        #make gaussian for PSF (best guess for GBT beam based on beamsize
        #report in GBT image)
        dc.makegaussian('gbt_gau.im',bmaj='55arcsec',bmin='55arcsec',
                        bpa='0deg',normalize=False)
        dc.close()
        os.system("rm -rf orion_tjoint3")
        dc.open('orion_tgbt_regrid.im',psf='gbt_gau.im')
        dc.setscales(scalemethod='uservector',uservector=[30.,50.,100.])
        dc.clean(algorithm='fullmsclean',model='orion_tjoint3',niter=30,
                 gain=0.3,mask=datapath+'orion.mask6',threshold='0.5Jy')
        dc.close()
        #default('clean')
        #ia.open('orion_tjoint3')
        #ia.calc(pixels='orion_tjoint3*"'+datapath+'orion.mask6"')
        #ia.close()

        im.open('orion.ms')
        im.selectvis(field=[2,3,4,5,6,7,8,9,10],spw=[0,1])
        im.defineimage(nx=300,cellx='2arcsec',phasecenter=6,spw=[0,1])
        im.setvp(dovp=True)
        #im.setscales(scalemethod='uservector',uservector=[0,3,10,30,100])
        im.setscales(scalemethod='uservector',uservector=[0,3,10,30,100])
        ###if clean component for large scale goes negative continue to use
        ##that scale
        im.setoptions(ftmachine='ft', padding=1.2)
        im.setmfcontrol(stoplargenegatives=-1, cyclefactor=5.0, cyclespeedup=500, flatnoise=False)
        im.weight(type='briggs',rmode='norm',robust=-1,mosaic=True)
        im.clean(algorithm='mfmultiscale', model='orion_tjoint3',
                 image='orion_tjoint3.image', gain=0.3, niter=4000,
                 mask='', threshold='4mJy')
        im.close()
        return time.time()

jointmemtime = joint_deconvolve(datapath)

endProc = time.clock()
endTime = time.time()

# Regression
import datetime
datestring = datetime.datetime.isoformat(datetime.datetime.today())
outfile = 'orion.' + datestring + '.log'
logfile = open(outfile, 'w')

test_results = {}
for k, fn in (('Feather 1',           'orion_tfeather.im'),
              ('Feather 2',           'orion_tfeather2.im'),
              ('SD Model (MS)',       'orion_tsdms.image'),
              ('SD Model (MEM)',      'orion_tsdmem.image'),
              ('Joint Deconvolution', 'orion_tjoint3.image')):
        if ia.open(fn):
                test_results[k] = ia.statistics(list=True, verbose=True)
                ia.close()
        else:
                print("Could not open", fn, "for reading!", file=logfile)

print('', file=logfile)
print('********** Data Summary *********', file=logfile)
print('           GBT image ', file=logfile)
print('Image name       : orion.gbt.im', file=logfile)
print('Object name      :', file=logfile)
print('Image type       : PagedImage', file=logfile)
print('Image quantity   : Intensity', file=logfile)
print('Pixel mask(s)    : mask0', file=logfile)
print('Region(s)        : None', file=logfile)
print('Image units      : Jy/beam', file=logfile)
print('Restoring Beam   : 98.0547 arcsec, 98.0547 arcsec, 79.25 deg', file=logfile)
print('Direction reference : J2000', file=logfile)
print('Spectral  reference : LSRK', file=logfile)
print('Velocity  type      : RADIO', file=logfile)
print('Rest frequency      : 1 Hz', file=logfile)
print('Telescope           : GBT', file=logfile)
print('Observer            : UNKNOWN', file=logfile)
print('Date observation    : UNKNOWN', file=logfile)
print('Axis Coord Type      Name             Proj Shape Tile   Coord value at pixel    Coord incr Units', file=logfile)
print('------------------------------------------------------------------------------------------------', file=logfile)
print('1    1     Direction Right Ascension   SIN   300  300  05:35:17.470   151.00 -2.000000e+00 arcsec', file=logfile)
print('2    1     Direction Declination       SIN   300  300 -05.23.06.790   151.00  2.000000e+00 arcsec', file=logfile)
print('3    2     Stokes    Stokes                    1    1             I', file=logfile)
print('4    3     Spectral  Frequency                 1    1    8.4351e+09     1.00  6.050000e+07 Hz', file=logfile)
print('                     Velocity                                  -inf     1.00           nan km/s', file=logfile)
print('   Observer: unavailable     Project: DSTST', file=logfile)
print('Observation: VLA(26 antennas)', file=logfile)
print('Data records: 1093716       Total integration time = 8545 seconds', file=logfile)
print('   Observed from   11:15:48   to   13:38:13', file=logfile)
print('Fields: 12', file=logfile)
print('  ID   Name          Right Ascension  Declination   Epoch', file=logfile)
print('  0    0518+165      05:21:09.89      +16.38.22.04  J2000', file=logfile)
print('  1    0539-057      05:41:38.09      -05.41.49.43  J2000', file=logfile)
print('  2    ORION1        05:35:07.42      -05.25.36.07  J2000', file=logfile)
print('  3    ORION2        05:35:17.42      -05.25.36.79  J2000', file=logfile)
print('  4    ORION3        05:35:27.42      -05.25.37.52  J2000', file=logfile)
print('  5    ORION4        05:35:27.47      -05.23.07.52  J2000', file=logfile)
print('  6    ORION5        05:35:17.47      -05.23.06.79  J2000', file=logfile)
print('  7    ORION6        05:35:07.47      -05.23.06.07  J2000', file=logfile)
print('  8    ORION7        05:35:07.52      -05.20.36.07  J2000', file=logfile)
print('  9    ORION8        05:35:17.52      -05.20.36.80  J2000', file=logfile)
print('  10   ORION9        05:35:27.52      -05.20.37.52  J2000', file=logfile)
print('  11   ORION10       05:35:32.61      -05.16.07.88  J2000', file=logfile)
print('Spectral Windows:  (2 unique spectral windows and 1 unique polarization setups)', file=logfile)
print('  SpwID  #Chans Frame Ch1(MHz)  Resoln(kHz) TotBW(kHz) Ref(MHz) Corrs', file=logfile)
print('  0           1 TOPO  8435.1    50000       50000      8435.1   RR  RL  LR  LL', file=logfile)
print('  1           1 TOPO  8485.1    50000       50000      8485.1   RR  RL  LR  LL', file=logfile)
print('*********************************', file=logfile)
print('', file=logfile)
print('********** Regression ***********', file=logfile)
print('*                               *', file=logfile)

#              Test name          Stat type Expected  Label irregularities
test_descs = (('Feather 1',           'max',  0.780,  ' '),
              ('Feather 2',           'max',  0.978,  ' '),
              ('SD Model (MS)',       'max',  1.14),
              ('SD Model (MEM)',      'max',  0.906),
              ('Joint Deconvolution', 'max',  1.10, '', 'Joint Decon1'), # 1.014
              ('Feather 1',           'flux', 242.506,  ' '),
              ('Feather 2',           'flux', 242.506,  ' '),
              ('SD Model (MS)',       'flux', 347, ' ', 'SD Model (MS)', 'SD Model MS'),
              ('SD Model (MEM)',      'flux', 286, '', 'SD Model (MEM)', 'Joint Deconvolution'),
              ('Joint Deconvolution', 'flux', 272, '', 'Joint Decon2')) # 360.468

def log_test_result(test_results, testdesc, logfile):
        """Append testdesc to logfile and return whether or not the test was
        successful."""
        result = test_results[testdesc[0]][testdesc[1]][0]
        reldiff = abs(1.0 - result / testdesc[2])
        if reldiff < 0.05:
                print('* Passed', end=' ', file=logfile)
                print('* Alright', end=' ')
                retval = True
        else:
                print('! FAILED', end=' ', file=logfile)
                print('#####FAILED')
                retval = False

        # RR 4/18/2009: I think this complication might stem from a bug in the
        # original version of the script, but since it is a regression script I
        # am hesitant to change the output.
        if len(testdesc) > 5:
                title1 = testdesc[5]
        else:
                title1 = testdesc[0]
        if len(testdesc) > 4:
                title2 = testdesc[4]
        else:
                title2 = testdesc[0]
        title2 += ':'
        if testdesc[1] == 'max':
                title1 += ' image max'
                title2 += ' Image max'
        else:
                title1 += ' ' + testdesc[1]
                title2 += ' ' + testdesc[1].title()
        if len(testdesc) > 3:
                title2 += testdesc[3]

        print(title1, 'test', file=logfile)
        print(title1, 'test')
        print('*--  ' + title2 + str(result) + ',' + str(testdesc[2]), file=logfile)
        print('*--  ' + title2 + str(result) + ',' + str(testdesc[2]))
        return retval

regstate = True
for td in test_descs:
        regstate &= log_test_result(test_results, td, logfile)

if regstate:
        print('---', file=logfile)
        print('Passed Regression test for Orion', file=logfile)
        print('---', file=logfile)
        print('')
        print('Regression PASSED')
        print('')
else:
        print('----FAILED Regression test for Orion', file=logfile)
        print('')
        print('Regression FAILED')
        print('')

print('*********************************', file=logfile)

print('', file=logfile)
print('', file=logfile)
print('********* Benchmarking *****************', file=logfile)
print('*                                      *', file=logfile)
print('Total wall clock time was: '+str(endTime - startTime), file=logfile)
print('Total CPU        time was: '+str(endProc - startProc), file=logfile)
print('Processing rate MB/s  was: '+str(760./(endTime - startTime)), file=logfile)
print('* Breakdown:                           *', file=logfile)
print('*   feather      time was: '+str(feathertime-startTime), file=logfile)
print('*   feathersynth time was: '+str(feathersynthtime-feathertime), file=logfile)
print('*   sdmodel (MS) time was: '+str(sdmodelmstime-feathersynthtime), file=logfile)
print('*   sdmodel(MEM) time was: '+str(sdmodelmemtime-sdmodelmstime), file=logfile)
print('*   joint decon  time was: '+str(jointmemtime-sdmodelmemtime), file=logfile)
print('*****************************************', file=logfile)
#
logfile.close()
