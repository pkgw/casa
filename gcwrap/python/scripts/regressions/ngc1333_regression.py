#############################################################################
#                                                                           #
# Test Name:                                                                #
#    Regression/Benchmarking Script for NGC 1333                            #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    This script illustrates data reduction of a mosaic.                    #
#      The dataset is a VLA SiO Mosaic of a YSO in NGC 1333 consisting of   #
#    10 source fields in an almost linear array oriented SE-NW.             #
#    Two spectral windows are centered on red & blue-shifted emission peaks.#
#    Central field also observed in 2 spectral widows centered on red and   #
#    blue-shifted emission peaks.                                           #
#      The NGC 1333 Mosaic was observed on 2 separate days: 2may03 & 8may03 #
#    The quality of the 2nd dataset is significantly worse and doesn't add  #
#    much to the image.                                                     #
#      Observations are of shocked SiO(1-0) emission in a bipolar outflow   #
#    (~ 43 GHz, VLA Q band receivers).  There is no detectable continuum    #
#    emission in this data.  There are 2 overlapping spectral windows with  #
#    starting frequencies of 43.4197 GHz & 43.4224 GHz.                     #
#    Each calibrator is observed in both spectral windows:                  #
#      0336+323 = Gain calibrator (Sv ~ 2.1 Jy on 2may; 2.8 Jy on 8may)     #
#      0542+498 = Flux density calibrator (3c147, Sv = 0.9144 Jy)           #
#         This flux calibrator is somewhat resolved. you can derive gain    #
#         solutions with uvrange=[0,50] (details are in the cookbook).      #
#      0319+415 = Bandpass calibrator (3c84, Sv ~ 4.9 Jy on 2may;           #
#         7.4 Jy on 8may)                                                   #
#    A good reference antenna for gain & bandpass calibration = ANT 26 for  #
#    both days.                                                             #
#    The test data are in uvfits files called N1333_1.UVFITS and            #
#    N1333_2.UVFITS.                                                        #
#                                                                           #
#    The NGC 1333 data on 2 May suffers from some correlator glitches,      #
#    a correlator re-boot, and some antennas were not fringing. The         #
#    correlator glitches are easy to find - amplitudes                      #
#    are 100s of Janskys during a glitch. Non-fringing antennas just        #
#    look like low-amplitude noise and are difficult to find during         #
#    initial inspection of the data - but gain                              #
#    solution plots show the amplitudes are near 0 for a given              #
#    timestamp. Even after you get out obvious bad data, you may still      #
#    want to clip the source data above a threshold.                        #
#                                                                           #
#    The 8 May data has some non-fringing antennas and is taken under       #
#    very poor weather conditions. There are major chunks of data on        #
#    some antennas that were flagged on-line so it is best to delete        #
#    what remains on these antennas.  Clipping of the final calibrated      #
#    source data is a must.                                                 #
#                                                                           #
# Features Tested:                                                          #
#    This script illustrates end-to-end processing with CASA as depicted    #
#    in the following flow-chart.                                           #
#                                                                           #
# Input Data             Process          Output Data                       #
#                                                                           #
# NGC1333_1.UVFITS---> importuvfits ----> n1333_1.ms                        #
#                           |                                               #
#                           v                                               #
#                        listobs    ----> casa.log                          #
#                           |                                               #
#                           v                                               #
#                       flagdata    ----> n1333_1.ms                        #
#                           |                                               #
#                           v                                               #
#                         setjy     ----> n1333_1.ms                        #
#                           |                                               #
#                           v                                               #
#                        gaincal    ----> n1333_1.[12].gcal                 #
#                           |                                               #
#                           v                                               #
#                       bandpass    ----> n1333_1.[12].bcal                 #
#                           |                                               #
#                           v                                               #
#                       fluxscale   ----> n1333_1.[12].fcal                 #
#                           |                                               #
#                           v                                               #
#                       applycal    ----> n1333_1.ms                        #
#                           |                                               #
#                           v                                               #
#                         split     ----> n1333_1.src.t2may.ms +            #
#                           |             n1333_1.gcal_t[12].2may.ms        #
#                           v                                               #
# NGC1333_2.UVFITS---> importuvfits ----> n1332_2.ms                        #
#                           |                                               #
#                           v                                               #
#                        listobs    ----> casa.log                          #
#                           |                                               #
#                           v                                               #
#                       flagdata    ----> n1333_2.ms                        #
#                           |                                               #
#                           v                                               #
#                         setjy     ----> n1333_2.ms                        #
#                           |                                               #
#                           v                                               #
#                        gaincal    ----> n1333_2.[12].gcal                 #
#                           |                                               #
#                           v                                               #
#                       bandpass    ----> n1333_2.[12].bcal                 #
#                           |                                               #
#                           v                                               #
#                       fluxscale   ----> n1333_2.[12].fcal                 #
#                           |                                               #
#                           v                                               #
#                       applycal    ----> n1333_2.ms                        #
#                           |                                               #
#                           v                                               #
#                         split     ----> n1333_2.src.t8may.ms +            #
#                           |             n1333_2.gcal_t[12].8may.ms        #
#                           v                                               #
#                        concat     ----> n1333_both.ms                     #
#                           |                                               #
#                           v                                               #
#                         mosaic    ----> n1333_both.model +                #
#                           |             n1333_both.residual +             #
#                           |             n1333_both.image +                #
#                           |             n1333_both.flux                   #
#                           v                                               #
#                       immoments   ----> n1333_both.src.tom0.red +         #
#                           |             n1333_both.src.tom0.blu +         #
#                           |             n1333_both.src.tom0.all +         #
#                           |             n1333_both.src.tom1.all           #
#                           v                                               #
#                 Get various statistics                                    #
#                                                                           #
# Success/failure criteria:                                                 #
#  Compare statistics with saved results from past trials.  Fail            #
#  when there is a significant variation from previous runs.                #
#                                                                           #
#############################################################################
# Updated by RRusk 2008-02-08 from tool-based ngc1333_regression.py         #
#############################################################################
#                                                                           #
# UPDATED BY gmoellen 2008-11-17:                                           #
#                                                                           #
#   There have long been several missteps in the flagging & calibration     #
#   which I have repaired:                                                  #
#                                                                           #
#   o clip of Amp>0.2 was being applied to all fields, and thus spuriously  #
#      flagging a very large fraction of 0319.                              #
#   o 0319 appears to be highly resolved.  Now use 0336 (integrated over    #
#      the whole observation) instead.  (0319 is now dropped entirely)      #
#   o added use of canned image model for 0542 (3C147_Q.im)                 #
#   o There *is* an opacity tip for 08May!  Using opacity=0.062 for 08May   #
#      and opacity=0.06 for 02May (from tip on 01May)                       #
#   o Use gaincurve for *all* calibration solve and apply, not just some    #
#   o The phase stability is not too bad, but scan-based solutions          #
#      (solint='inf') in gaincal are too long, especially for the long      #
#      scans on the flux density calibrator.  Using solint='int' yields a   #
#      much more consistent set of flux density results for 0336, among the #
#      spws, and (marginally) between the two dates.  The net effect is to  #
#      decrease the flux density of 0336 and the ngc1333.                   #
#   o I've added a few a priori flag commands to deal with bad calibrator   #
#      scans, based on poor solutions and/or calibrated data                #
#   o The regression tests now do channel-averaging (to avoid goofy edge    #
#      channels), and the test values have been revised accordingly. The    #
#      final peak SNR (~9) is better than before (~6), even though the      #
#      peak flux density is less (due to better f.d. calibration).          #
#   o Added a variable to optionally disable the imaging/moments portion    #
#                                                                           #
#############################################################################
import os
import sys
import time
import regression_utility as tstutl

# Enable benchmarking?
benchmarking = True
usedasync = False

# Optionally turn off imaging (which is very time-consuming)
doimage = True

#
# Set up some useful variables
#
# Get path to CASA home directory by stipping name from '$CASAPATH'
pathname=os.environ.get('CASAPATH').split()[0]

# This is where the NGC1333 UVFITS data will be
datapath=pathname+'/data/regression/ATST2/NGC1333/'
fitsdata1=datapath+'N1333_1.UVFITS'
fitsdata2=datapath+'N1333_2.UVFITS'

# 3C147 model image for setjy
modelim=pathname+'/data/nrao/VLA/CalModels/3C147_Q.im'

# The testdir where all output files will be kept
testdir='ngc1333_regression'

# The prefix to use for output files.
prefix=testdir+"/"+'n1333_both'
prefix1=testdir+"/"+'n1333_1'
prefix2=testdir+"/"+'n1333_2'

# Make new test directory
# (WARNING! Removes old test directory of the same name if one exists)
tstutl.maketestdir(testdir)

# Start benchmarking
if benchmarking:
    startTime = time.time()
    startProc = time.clock()

#
#=====================================================================
#
# Import the data from FITS to MS
#
print('*** 02 MAY ***')
print('--Import--')

# Safest to start from task defaults
default('importuvfits')

# Set up the MS filename and save as new global variable
msfile1 = prefix1 + '.ms'

# Use task importuvfits
fitsfile = fitsdata1
vis = msfile1
antnamescheme="new"
importuvfits()

# Record import time
if benchmarking:
    importtime = time.time()

#
#=====================================================================
#
# List a summary of the MS
#
print('--Listobs--')

# Don't default this one.  Make use of the previous setting of
# vis.  Remember, the variables are GLOBAL!

# You may wish to see more detailed information, like the scans.
# In this case use the verbose = True option
verbose = True

listobs()

# Record listing completion time
if benchmarking:
    listtime = time.time()

#
#=====================================================================
#  Flagging
#
print('--Flagdata--')

#
# The following information on bad data comes from a test report
# created by Debra Shepherd.  It is currently available at
# http://aips2.nrao.edu/projectoffice/almatst1.1/TST1.1.data.description.pdf
#
# NGC 1333: bad data on 2 May:
# o Correlator glitches
#   - Spwid 1, field 1: 02-May-2003/21:44:00 to 21:50:00
#   - Spwid 1, fields 2,3,4,5,6: 02-May-2003/21:40:00 to 21:51:00
#   - Spwid 1, field 7: 02-May-2003/21:56:00 to 21:58:00
#   - Spwid 2, fields 8,9,10,11,12: 02-May-2003/21:55:00 to 22:20:00
# o Non-fringing antennas:
#   - Ant 9 - everything: all spwids, all fields, all times
#   - Ant 14 - spwid 1, 02-May-2003/18:36:00 to 18:53:00
#   - Ant 14 - spwid 2, 02-May-2003/18:52:00 to 19:10:00
# o Bad bandpass solutions:
#   - Ant 22 - spwid 2 (pathological bandpass solution,
#     don't trust this antenna).
# o End channels:
#   - Flag channels 1,2,3 and 61, 62, 63 (very noisy)
#     Flag them to prevent higher noise in some image planes.
# o Even after all this, the calibrated source data still has amplitudes
# that vary from a maximum of 10Jy/channel to 30Jy/channel.
# The uv weighting will properly downweight the poor data.
#





# Use flagdata() in vector mode
default('flagdata')
vis = msfile1
mode = 'list'


##################################################
#
# Get rid of the autocorrelations from the MS
#
# Flag antenna with antennaid 8
#
# Flag all data whose amplitude  are not in range [0.0,2.0] on the
# parallel hands       ( not done )
#
# Flag data (which is bad) in a time range
#
# Flag all antenna 14, 15 data in the time ranges stated
#
# Sequential flagdata execution:
#
# flagdata(vis=msfile1, autocorr=True, mode='manualflag')
# flagdata(vis=msfile1,antenna='9', mode='manualflag')
# flagdata(vis=msfile1,mode='manualflag',clipexpr='ABS RR',
#          clipminmax=[0.0,2.0],clipoutside=True)
#
# flagdata(vis=msfile1,mode='manualflag',clipexpr='ABS LL',
#          clipminmax=[0.0,2.0],clipoutside=True)
# flagdata(vis=msfile1,mode='manualflag',timerange='2003/05/02/21:40:58~2003/05/02/22:01:30')
#
# flagdata(vis=msfile1,mode='manualflag', antenna='14',
#          timerange='2003/05/02/18:50:50~2003/05/02/19:13:30')
# flagdata(vis=msfile1,mode='manualflag', antenna='15', spw='0',
#          timerange='2003/05/02/22:38:49~2003/05/02/22:39:11')

# Parallel flagdata execution:
flagcmds = ["autocorr=True",
       "autocorr=False antenna='VA09'",
       "autocorr=False timerange='2003/05/02/21:40:58~2003/05/02/22:01:30'",
       "autocorr=False antenna='VA14' timerange='2003/05/02/18:50:50~2003/05/02/19:13:30'",
       "autocorr=False antenna='VA15' timerange='2003/05/02/22:38:49~2003/05/02/22:39:11' spw='0'"]

inpfile=flagcmds
#autocorr  = [true, false , false, false , false ]
#antenna   = [''  , 'VA09', ''   , 'VA14', 'VA15']
#timerange = [''  , ''    , '2003/05/02/21:40:58~2003/05/02/22:01:30', \
#                           '2003/05/02/18:50:50~2003/05/02/19:13:30', \
#                           '2003/05/02/22:38:49~2003/05/02/22:39:11'] 
#spw       = [''  , ''    , ''   , ''    , '0'   ]


###################################################
# Finally, apply all the flag specifications
#

flagdata()


# Record flagging completion time
if benchmarking:
    flagtime = time.time()


#
#=====================================================================
#
# Set the fluxes of the primary calibrator(s)
#
print('--Setjy--')
default('setjy')

setjy(vis=msfile1,field='0542+498_1',modimage=modelim,scalebychan=False,standard='Perley-Taylor 99') 
setjy(vis=msfile1,field='0542+498_2',modimage=modelim,scalebychan=False,standard='Perley-Taylor 99')

# Record setjy completion time
if benchmarking:
    setjytime = time.time()


#
#=====================================================================
#
# Prior calibrations (opacity and gaincurve)
#
print('--Gencal (opacity, gaincurve)--')
default('gencal')

opactable = prefix1 + '.opac'
gencal(vis=msfile1,caltable=opactable,
       caltype='opac',
       parameter=[0.06])

gcurvetable = prefix1 + '.gaincurve'
gencal(vis=msfile1,caltable=gcurvetable,
       caltype='gc')

# Record gencal completion time
if benchmarking:
    gencaltime = time.time()


#
#=====================================================================
#
# Gain calibration
#
print('--Gaincal--')
default('gaincal')

gtable1 = prefix1 + '.1.gcal'
gaincal(vis=msfile1,caltable=gtable1,
	field='0,12,14',spw='0:4~58', gaintype='G',
	solint='int',combine='',refant='VA27',minsnr=2.,
        gaintable=[opactable,gcurvetable],interp=['',''])

gtable2 = prefix1 + '.2.gcal'
gaincal(vis=msfile1,caltable=gtable2,
	field='6,13,15',spw='1:4~58', gaintype='G',
	solint='int',combine='',refant='VA27',
        gaintable=[opactable,gcurvetable],interp=['',''])

# gaincal calibration completion time
if benchmarking:
    gaintime = time.time()

#
#=====================================================================
#
# Bandpass calibration
#
print('--Bandpass--')
default('bandpass')

btable1 = prefix1 + '.1.bcal'
bandpass(vis=msfile1,caltable=btable1,
	 field='0',spw='0',
	 gaintable=[opactable,gcurvetable,gtable1],interp=['','','nearest'],
	 refant='VA27',solint='inf',combine='scan')
btable2 = prefix1 + '.2.bcal'
bandpass(vis=msfile1,caltable=btable2,
	 field='6',spw='1',
	 gaintable=[opactable,gcurvetable,gtable2],interp=['','','nearest'],
	 refant='VA27',solint='inf',combine='scan')

# bandpass calibration completion time
if benchmarking:
    bptime = time.time()

#
#=====================================================================
#
# Bootstrap flux scale
#   Transfer the flux density  from flux calibrater to gain calibraters
#
print('--Fluxscale--')
default('fluxscale')

ftable1 = prefix1 + '.1.fcal'
fluxscale(vis=msfile1,caltable=gtable1,fluxtable=ftable1,
	  reference='0542+498_1',transfer=['0336+323_1', '0319+415_1'])
ftable2 = prefix1 + '.2.fcal'
fluxscale(vis=msfile1,caltable=gtable2,fluxtable=ftable2,
	  reference='0542+498_2',transfer=['0336+323_2', '0319+415_2'])

# Record fluxscale completion time
if benchmarking:
    fstime = time.time()

#
#=====================================================================
#
# Apply our calibration solutions to the data
# (This will put calibrated data into the CORRECTED_DATA column)
#
print('--ApplyCal--')
default('applycal')

applycal(vis=msfile1,
	 field='0~5',spw='0',
	 gaintable=[opactable,gcurvetable,ftable1,btable1],
	 gainfield=['','','0'],calwt=False)
applycal(vis=msfile1,
	 field='6~11',spw='1',
	 gaintable=[opactable,gcurvetable,ftable2,btable2],
	 gainfield=['','','6'],calwt=False)

# Record applycal completion time
if benchmarking:
    correcttime = time.time()


#
#=====================================================================
#
# Split the target and gain calibrater data
#
print('--Split (target and cals) --')
default('split')

# Split out the corrected data of target
splitms1 = prefix1 + '.src.t2may.ms'
split(vis=msfile1,outputvis=splitms1,
      field='1~5,7~11',spw='0;1:0~62',
      datacolumn='corrected')

# Record split src data completion time
if benchmarking:
    splitsrctime = time.time()

# Split out calibraters
calsplitms1 = prefix1 + '.gcal_t1.2may.ms'
split(vis=msfile1,outputvis=calsplitms1,
      field='0',spw='0:0~62',
      datacolumn='corrected')
calsplitms2 = prefix1 + '.gcal_t2.2may.ms'
split(vis=msfile1,outputvis=calsplitms2,
      field='6',spw='1:0~62',
      datacolumn='corrected')

# Record split cal completion time
if benchmarking:
    splitcaltime = time.time()

#
#=====================================================================
#
# Import the data from FITS to MS
#
print('*** 08 MAY ***')
print('--Import--')

# Safest to start from task defaults
default('importuvfits')

# Set up the MS filename and save as new global variable
msfile2 = prefix2 + '.ms'

# Use task importuvfits
fitsfile = fitsdata2
vis = msfile2
antnamescheme='new'
importuvfits()

# Note that there will be a ngc5921.ms.flagversions
# containing the initial flags as backup for the main ms flags.

# Record import time
if benchmarking:
    importtime2 = time.time()

#
#=====================================================================
#
# List a summary of the MS
#
print('--Listobs--')

# Don't default this one.  Make use of the previous setting of
# vis.  Remember, the variables are GLOBAL!

# You may wish to see more detailed information, like the scans.
# In this case use the verbose = True option
verbose = False

listobs()

# Record listing completion time
if benchmarking:
    listtime2 = time.time()

#
##################################################
#
#  Flagging
#

print('--Flagdata--')

#
# Bad data as identified in report by Debra Shepherd at
# http://aips2.nrao.edu/projectoffice/almatst1.1/TST1.1.data.description.pdf
#
# NGC 1333: bad data on 8 May:
# o Non-fringing antennas:
#   - Ant 9 - everything: all spwids, all fields, all times
#   - Ant 15 - spwid 1, 08-May-2003/16:32:10 to 16:46:10
#   - Ant 14 - spwid 2, 08-May-2003/16:55:10 to 17:10:00
# o Antennas that were mostly flagged on-line - flag the rest of the data:
#   - Ants 10 & 19 - everything: all spwids, all fields, all times
# o End channels:
#   - Flag channels 1,2,3 and 61, 62, 63 (very noisy)
# o Even after all this, the calibrated source data still has amplitudes that
#   vary from a maximum of 40Jy/channel to 140Jy/channel (this data was taken
#   under significantly worse conditions than on 2may). I suggest you flag
#   the worst timestamps:
#     - 08-May-2003/17:03:00 to 17:10:00
#     - 08-May-2003/18:50:00 to 18:57:00
# o The uv weighting will properly down-weight the remaining poor data
#   (and down-weight the 8may data relative to the better 2may data).
#   Now you know why this 2nd day of data doesn't add much to the
#   final image RMS.
#

default('flagdata')
vis = msfile2
mode = 'list'

#
#
##################################################
#
# Get rid of the autocorrelations from the MS
#
#autocorr[0] = true

#
# Flag all data from antennas 8,9,18
#

#flagdata(vis=msfile2,antenna='9,10,19', mode='manualflag')

#antenna[1] = 'VA09,VA10,VA19'

#
# Flag antenna 14 for timerange.  ANTENNAID 14 is ANTENNANAME 15
#

#default('flagdata')
#flagdata(vis=msfile2,mode='manualflag',
#	 antenna='15',
#	 timerange='2003/05/08/00:00:00~2003/05/08/20:00:00')

#
#  antenna '22' has a bad scan 144
#flagdata(vis=msfile2,mode='manualflag',
#	 antenna='6,22',
#	 scan='144')

#
#
# Flag the last data channels
#

#default('flagdata')
#flagdata(vis=msfile2, spw='*:59;60;61;62', mode='manualflag')

#spw[4] = '*:59;60;61;62'

flagcmds = ["autocorr=True",
            "antenna='VA09,VA10,VA19'",
            "antenna='VA15' timerange='2003/05/08/00:00:00~2003/05/08/20:00:00'",
            "antenna='VA06,VA22' scan='144'",
            "spw='*:59;60;61;62'"]

inpfile=flagcmds
#
###################################################
# Finally, apply all the flag specifications
#

flagdata()

# Record flagging completion time
if benchmarking:
    flagtime2 = time.time()


#
#=====================================================================
#
# Set the fluxes of the primary calibrator(s)
#
print('--Setjy--')
default('setjy')

setjy(vis=msfile2,field='0542+498_1',modimage=modelim,scalebychan=False,standard='Perley-Taylor 99')
setjy(vis=msfile2,field='0542+498_2',modimage=modelim,scalebychan=False,standard='Perley-Taylor 99')

# Record setjy completion time
if benchmarking:
    setjytime2 = time.time()


#
#=====================================================================
#
# Prior calibrations (opacity and gaincurve)
#
print('--Gencal (opacity, gaincurve)--')
default('gencal')

opactable2 = prefix2 + '.opac'
gencal(vis=msfile1,caltable=opactable2,
       caltype='opac',
       parameter=[0.062])

gcurvetable2 = prefix2 + '.gaincurve'
gencal(vis=msfile1,caltable=gcurvetable2,
       caltype='gc')

# gencal2 completion time
if benchmarking:
    gencaltime2 = time.time()

#
#=====================================================================
#
# Gain calibration
#
print('--Gaincal--')
default('gaincal')

gtable2_1 = prefix2 + '.1.gcal'
gaincal(vis=msfile2,caltable=gtable2_1,
	field='0,12,14',spw='0:4~58', gaintype='G',
	solint='int',combine='',refant='VA27',
        gaintable=[opactable2,gcurvetable2],interp=['',''])
gtable2_2 = prefix2 + '.2.gcal'
gaincal(vis=msfile2,caltable=gtable2_2,
	field='6,13,15',spw='1:4~58', gaintype='G',
	solint='int',combine='',refant='VA27',
        gaintable=[opactable2,gcurvetable2],interp=['',''])

# gaincal calibration completion time
if benchmarking:
    gaintime2 = time.time()


#
#=====================================================================
#
# Bandpass calibration
#
print('--Bandpass--')
default('bandpass')

btable2_1 = prefix2 + '.1.bcal'
bandpass(vis=msfile2,caltable=btable2_1,
	 field='0',spw='0',
	 gaintable=[opactable2,gcurvetable2,gtable2_1],interp=['','','nearest'],
	 refant='VA27',
	 solint='inf',combine='scan')
btable2_2 = prefix2 + '.2.bcal'
bandpass(vis=msfile2,caltable=btable2_2,
	 field='6',spw='1',
	 gaintable=[opactable2,gcurvetable2,gtable2_2],interp=['','','nearest'],
	 refant='VA27',
	 solint='inf',combine='scan')

# bandpass calibration completion time
if benchmarking:
    bptime2 = time.time()

#
#=====================================================================
#
# Bootstrap flux scale
#   Transfer the flux density  from flux calibrater to gain calibraters
#
print('--Fluxscale--')
default('fluxscale')

ftable2_1 = prefix2 + '.1.fcal'
fluxscale(vis=msfile2, caltable=gtable2_1, fluxtable=ftable2_1,
	  reference='0542+498_1',transfer=['0336+323_1', '0319+415_1'])
ftable2_2 = prefix2 + '.2.fcal'
fluxscale(vis=msfile2, caltable=gtable2_2, fluxtable=ftable2_2,
	  reference='0542+498_2',transfer=['0336+323_2', '0319+415_2'])

# Record fluxscale completion time
if benchmarking:
    fstime2 = time.time()

#
#=====================================================================
#
# Apply our calibration solutions to the data
# Do the correction for the above solutions of bandpass and gain,
# including gain curve.
#
print('--ApplyCal--')
default('applycal')

applycal(vis=msfile2,
	 field='0~5',spw='0',
	 gaintable=[opactable2,gcurvetable2,ftable2_1,btable2_1],
	 gainfield=['','','0'],calwt=False)
applycal(vis=msfile2,
	 field='6~11',spw='1',
	 gaintable=[opactable2,gcurvetable2,ftable2_2,btable2_2],
	 gainfield=['','','6'],calwt=False)


# Record applycal completion time
if benchmarking:
    correcttime2 = time.time()

#
#=====================================================================
#
# Split the target and gain calibrater data
#
print('--Split (target and cals) --')
default('split')

# Split out the corrected data of target
splitms2 = prefix2 + '.src.t8may.ms'
split(vis=msfile2,outputvis=splitms2,
      field='1~5,7~11',spw='0;1:0~62',
      datacolumn='corrected')

# Record split src data completion time
if benchmarking:
    splitsrctime2 = time.time()

# Split out calibraters
calsplitms2_1 = prefix2 + '.gcal_t1.8may.ms'
split(vis=msfile2,outputvis=calsplitms2_1,
      field='0',spw='0:0~62',datacolumn='corrected')
calsplitms2_2 = prefix2 + '.gcal_t2.8may.ms'
split(vis=msfile2,outputvis=calsplitms2_2,
      field='6',spw='1:0~62',datacolumn='corrected')

# Record split cal completion time
if benchmarking:
    splitcaltime2 = time.time()


#
#=====================================================================
#
# Merge the 2 corrected data sets
#
print('--Concatenate data sets--')
default('concat')

msfileboth = prefix + '.ms'
shellcmd = 'cp -r '+splitms1+' '+msfileboth
os.system(shellcmd)
vis = [msfileboth, splitms2]
concatvis = msfileboth
freqtol='10MHz'
dirtol='1arcsec'
concat()

# Record concatenation completion time
if benchmarking:
    concattime = time.time()

#
#=====================================================================
# Imaging using mosaic
#
# 63 channels in each spectral window (spwid) with about 30 channel
# overlap so about 96 total independent channels.  Imager will figure
# out how the windows overlap and will place them properly on the grid.
#
# Set total mosaic size with nx, ny; starting with the third channel (2),
# make 18 channels total with 5-channel averaging.  The mosaic phase
# centre is fieldid=0.  Use both spectral windows.
#

print('--Image data--')

"""
im.open(thems=msfileboth)
im.selectvis(field=range(0,10),spw=[0,1],nchan=63,start=0,step=1)
im.defineimage(nx=800,ny=800,cellx='0.5arcsec',stokes='I',mode='channel',nchan=18,start=2,step=5,phasecenter=0,spw=[0,1])
### natural weighting
im.weight(type='natural')
### Use default primary beam
im.setvp(dovp=True,usedefaultvp=True,dosquint=False)
### Use the newer mosaic gridder where mosaicing is done in convolution
im.setoptions(padding=1.0,ftmachine='mosaic')
#### Define up to which point in the beam image will be shown
imfluxscale = prefix+'.src.task.fluxscale'
im.setmfcontrol(scaletype='PBCOR',minpb=0.07,constpb=0.4, fluxscale=imfluxscale)
immodel = prefix+'.src.tall.cln.model'
imimage = prefix+'.src.task.image'
imresidual = prefix+'.src.tall.cln.resid'
im.clean(algorithm='mfclark',niter=1,gain=0.1,threshold=5.,model=immodel,
         mask='',image=imimage,residual=imresidual,npercycle=30)
im.close()
"""

imagetime=0.0
momenttime=0.0
if doimage:
    default('clean')

    vis = msfileboth
    imagename = prefix

    mode = 'channel'
### image centered on field 0, making 18 channels
### averaging 5 data channels to 1 image channel
    nchan = 18
    start = 4
    width = 5
    psfmode = 'clark'
    imsize = [800,800]
    cell = ['0.5arcsec','0.5arcsec']
    stokes='I'
#### clean with niter of 1 ...basically doing a just a dirty mosaic
    niter = 1
    gain = 0.1
    threshold = '5.0mJy'
    mask = ''
#cleanbox = []
###selecting data from 10 fields 0-9 all channels of the 2 spectral windows
    field = '0~9'
    spw = '0~1'
    timerange = ''
    restfreq = ''
    phasecenter='0'
    modelimage = ''
    weighting = 'natural'
    imagermode='mosaic'
    mosweight = False
    ftmachine = 'mosaic'
    flatnoise=True
    cyclefactor = 1.5
    cyclespeedup = -1
    scaletype = 'SAULT'
    pbcor=False
    minpb = 0.1
    sigma = '0.001Jy'
    targetflux = '1.0Jy'
    constrainflux = False
    prior = ['']
    negcomponent = 2
    scales = [0, 3, 10]
    interpolation = 'nearest'
    async = False
    
    clean()

# Record imager completion time
    if benchmarking:
        imagetime = time.time()

#
#=====================================================================
# Do moments
#

    print('--Calculate moments--')

    default('immoments')

    imimage=prefix+'.image'
    mom0redoutfile = prefix+'.src.tmom0.red'
    mom0blueoutfile = prefix+'.src.tmom0.blu'
    mom0alloutfile = prefix+'.src.tmom0.all'
    mom1alloutfile = prefix+'.src.tmom1.all'
    
    imagename = imimage
    moments = [0]
    axis = 'spec'
    chans = '2~8'
    includepix=[0.003,100.0]
    excludepix = [-1]
    outfile = mom0redoutfile
    async = False
    immoments()
    outfile = mom0blueoutfile
    chans='9~15'
    immoments()

    outfile = mom0alloutfile
    chans='2~15'
    immoments()

    outfile = mom1alloutfile
    moments = 1
    includepix=[0.02,100.0]
    immoments()
# Record moment estimation completion time
    if benchmarking:
        momenttime = time.time()

"""
ia.open(infile=imimage)
ia.moments(outfile=mom0redoutfile,
           moments=0,axis=3,mask='indexin(3,[2:9])',includepix=[0.003,100.0])
ia.moments(outfile=mom0blueoutfile,
           moments=0,axis=3,mask='indexin(3,[9:15])',includepix=[0.003,100.0])
ia.moments(outfile=mom0alloutfile,
           moments=0,axis=3,mask='indexin(3,[2:15])',includepix=[0.003,100.0])
ia.moments(outfile=mom1alloutfile,
           moments=1,axis=3,mask='indexin(3,[2:15])',includepix=[0.02,100.0])
ia.close()
"""


endProc = time.clock()
endTime = time.time()

#
#=====================================================================
# Do regression test
#

ms.open(thems=splitms1)
ms.selectchannel(1,2,59,1);
src_2may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close

ms.open(thems=splitms2)
ms.selectchannel(1,2,59,1);
src_8may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close

ms.open(thems=calsplitms2_1)
ms.selectinit(0,T)
ms.selectchannel(1,2,59,1);
gcal1_8may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close()

ms.open(thems=calsplitms2_2)
ms.selectchannel(1,2,59,1);
gcal2_8may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close()

ms.open(thems=calsplitms1)
ms.selectchannel(1,2,59,1);
gcal1_2may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close()

ms.open(thems=calsplitms2)
ms.selectchannel(1,2,59,1);
gcal2_2may=max(ms.range(items=["amplitude"]).get('amplitude'))
ms.close()

thistest_immax=0.0
thistest_imrms=0.0
if doimage:
    ia.open(infile=mom0alloutfile)
    statistics=ia.statistics(list=True, verbose=True)
    thistest_immax=statistics['max'][0]
    thistest_imrms=statistics['rms'][0]
    ia.close()
# 7499.6601 (n1333_both.ms)

# test values 
cal1_2may=4.56 # (channel averaged)
cal2_2may=4.00 # (channel averaged)
cal1_8may=5.72 # (channel averaged)
cal2_8may=6.70 # (channel averaged)
src2may=3.13   # (channel averaged)
src8may=8.18   # (channel averaged)
immax=0.70
imrms=0.069

diff_cal1_2may=abs((cal1_2may-gcal1_2may)/cal1_2may)
diff_cal2_2may=abs((cal2_2may-gcal2_2may)/cal2_2may)
diff_cal1_8may=abs((cal1_8may-gcal1_8may)/cal1_8may)
diff_cal2_8may=abs((cal2_8may-gcal2_8may)/cal2_8may)
diff_src_2may=abs((src2may-src_2may)/src2may)
diff_src_8may=abs((src8may-src_8may)/src8may)
diff_immax=abs((immax-thistest_immax)/immax)
diff_imrms=abs((imrms-thistest_imrms)/imrms)


if not benchmarking:
    print('')
    print('--- Done ---')
else:
    import datetime
    datestring=datetime.datetime.isoformat(datetime.datetime.today())
    outfile='ngc1333.'+datestring+'.log'
    logfile=open(outfile,'w')

    for x in [sys.stdout, logfile]: print('', file=x)
    for x in [sys.stdout, logfile]: print('********** Data Summary *********', file=x)
    for x in [sys.stdout, logfile]: print('*********************************', file=x)
    for x in [sys.stdout, logfile]: print('', file=x)
    for x in [sys.stdout, logfile]: print('********** Regression ***********', file=x)
    for x in [sys.stdout, logfile]: print('*                               *', file=x)
    regstate = True
    if (diff_cal1_2may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed cal1 max amplitude test (2may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed cal1 max amplitude test (2may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Cal1 max amp (2may) '+str(gcal1_2may)+' ('+str(cal1_2may)+')', file=x)
    if (diff_cal2_2may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed cal2 max amplitude test (2may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed cal2 max amplitude test (2may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Cal2 max amp (2may) '+str(gcal2_2may)+' ('+str(cal2_2may)+')', file=x)
    if (diff_cal1_8may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed cal1 max amplitude test (8may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed cal1 max amplitude test (8may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Cal1 max amp (8may) '+str(gcal1_8may)+' ('+str(cal1_8may)+')', file=x)
    if (diff_cal2_8may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed cal2 max amplitude test (8may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed cal2 max amplitude test (8may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Cal2 max amp (8may) '+str(gcal2_8may)+' ('+str(cal2_8may)+')', file=x)
    if (diff_src_2may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed src max amplitude test (2may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed src max amplitude test (2may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Src max amp (2may) '+str(src_2may)+' ('+str(src2may)+')', file=x)
    if (diff_src_8may < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed src max amplitude test (8may) *', file=x)
    else:
        for x in [sys.stdout, logfile]: print('* Failed src max amplitude test (8may) *', file=x)
        regstate = False
    for x in [sys.stdout, logfile]: print('*   Src max amp (8may) '+str(src_8may)+' ('+str(src8may)+')', file=x)
    if (diff_immax < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed image max test                *', file=x)
    else:
        if doimage:
            for x in [sys.stdout, logfile]: print('* Failed image max test                *', file=x)
            regstate = False
        else:
            for x in [sys.stdout, logfile]: print('* Did not do image max test', file=x)
            
    for x in [sys.stdout, logfile]: print('*   Image max '+str(thistest_immax)+' ('+str(immax)+')', file=x)
    if (diff_imrms < 0.05):
        for x in [sys.stdout, logfile]: print('* Passed image rms test                *', file=x)
    else:
        if doimage:
            for x in [sys.stdout, logfile]: print('* Failed image rms test                *', file=x)
            regstate = False
        else:
            for x in [sys.stdout, logfile]: print('* Did not do image rms test', file=x)
    for x in [sys.stdout, logfile]: print('*   Image rms '+str(thistest_imrms)+' ('+str(imrms)+')', file=x)
    #if ((diff_cal1_2may<0.05) & (diff_cal2_2may<0.05) &
    #    (diff_cal1_8may<0.05) & (diff_cal2_8may<0.05) &
    #    (diff_src_2may<0.05) & (diff_src_8may<0.05) &
    #    (diff_immax<0.05) & (diff_imrms<0.05)): 
    if regstate:
        for x in [sys.stdout, logfile]: print('---', file=x)
	for x in [sys.stdout, logfile]: print('Passed Regression test for NGC1333', file=x)
	for x in [sys.stdout, logfile]: print('---', file=x)
        tstutl.note("Passed Regression test for NGC1333","NORMAL")
    else: 
	for x in [sys.stdout, logfile]: print('----FAILED Regression test for NGC1333', file=x)
        tstutl.note("FAILED Regression test for NGC1333","SEVERE")
    for x in [sys.stdout, logfile]: print('*********************************', file=x)
    for x in [sys.stdout, logfile]: print('', file=x)
    for x in [sys.stdout, logfile]: print('', file=x)
    for x in [sys.stdout, logfile]: print('********* Benchmarking *****************', file=x)
    for x in [sys.stdout, logfile]: print('*                                      *', file=x)
    for x in [sys.stdout, logfile]: print('Total wall clock time was: '+str(endTime - startTime), file=x)
    for x in [sys.stdout, logfile]: print('Total CPU        time was: '+str(endProc - startProc), file=x)
    for x in [sys.stdout, logfile]: print('Processing rate MB/s  was: '+str(240.3/(endTime - startTime)), file=x)
    for x in [sys.stdout, logfile]: print('* Breakdown:                           *', file=x)
    for x in [sys.stdout, logfile]: print('*   import       time was: '+str(importtime-startTime), file=x)
    for x in [sys.stdout, logfile]: print('*   listobs      time was: '+str(listtime-importtime), file=x)
    for x in [sys.stdout, logfile]: print('*   flagdata     time was: '+str(flagtime-listtime), file=x)
    for x in [sys.stdout, logfile]: print('*   setjy        time was: '+str(setjytime-flagtime), file=x)
    for x in [sys.stdout, logfile]: print('*   gencal       time was: '+str(gencaltime-setjytime), file=x)
    for x in [sys.stdout, logfile]: print('*   gaincal      time was: '+str(gaintime-gencaltime), file=x)
    for x in [sys.stdout, logfile]: print('*   bandpass     time was: '+str(bptime-gaintime), file=x)
    for x in [sys.stdout, logfile]: print('*   fluxscale    time was: '+str(fstime-bptime), file=x)
    for x in [sys.stdout, logfile]: print('*   correct      time was: '+str(correcttime-fstime), file=x)
    for x in [sys.stdout, logfile]: print('*   split        time was: '+str(splitcaltime-correcttime), file=x)
    for x in [sys.stdout, logfile]: print('*   import       time was: '+str(importtime2-splitcaltime), file=x)
    for x in [sys.stdout, logfile]: print('+   listobs      time was: '+str(listtime2-importtime2), file=x)
    for x in [sys.stdout, logfile]: print('*   flagdata     time was: '+str(flagtime2-listtime2), file=x)
    for x in [sys.stdout, logfile]: print('*   setjy        time was: '+str(setjytime2-flagtime2), file=x)
    for x in [sys.stdout, logfile]: print('*   gencal       time was: '+str(gencaltime2-setjytime2), file=x)
    for x in [sys.stdout, logfile]: print('*   gaincal      time was: '+str(gaintime2-gencaltime2), file=x)
    for x in [sys.stdout, logfile]: print('*   bandpass     time was: '+str(bptime2-gaintime2), file=x)
    for x in [sys.stdout, logfile]: print('*   fluxscale    time was: '+str(fstime2-bptime2), file=x)
    for x in [sys.stdout, logfile]: print('*   correct      time was: '+str(correcttime2-fstime2), file=x)
    for x in [sys.stdout, logfile]: print('*   split        time was: '+str(splitcaltime2-correcttime2), file=x)
    for x in [sys.stdout, logfile]: print('*   concatenate  time was: '+str(concattime-splitcaltime2), file=x)
    for x in [sys.stdout, logfile]: print('*   image        time was: '+str(imagetime-concattime), file=x)
    for x in [sys.stdout, logfile]: print('*   moments      time was: '+str(momenttime-imagetime), file=x)
    for x in [sys.stdout, logfile]: print('*****************************************', file=x)

    #
    logfile.close()
