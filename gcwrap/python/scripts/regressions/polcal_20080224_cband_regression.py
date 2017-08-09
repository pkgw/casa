##########################################################################
#                                                                        #
# Use Case Script for POLCAL 6cm Data                                    #
# Using POLCA data 20080224 BnC-config C-band                            #
#                                                                        #
# Last Updated STM 2008-05-23 (Beta Patch 2)                             #
# Updated      STM 2008-06-11 (Beta Patch 2.0)                           #
#    Uses new clean task                                                 #
# Updated      STM 2008-06-18 (Beta Patch 2.0) Format as regression      #
# Updated      STM 2008-09-17 (Beta Patch 3.0) Use imval task            #
# Updated      GAM 2013-05-15 Removed redundant gaincurve/opacity refs   #
# Updated      GAM 2013-11-25 Removed (deprecated) plotxy clause         #
# Update       GAM 2017-02-24 Use task command syntax throughout         #
#                             removed flagmanager calls(flagdata does it)#
#                             various fixes to bring up-to-date          #
#                                                                        #
##########################################################################

import time
import os
import pickle
import pylab as pl

pathname=os.environ.get('CASAPATH').split()[0]

#
#=====================================================================
# SET UP THE SCRIPT CONTROL PARAMETERS HERE
#=====================================================================
#
# Set up some useful variables to control subsequent actions:

# This script may have some interactive commands: scriptmode = True
# if you are running it and want it to stop during interactive parts.
scriptmode = True

# Enable benchmarking?
benchmarking = True

# This name will prefix all output files
prefix = 'polcal_20080224.cband.regression'

# This name is the prefix all input (script,regression) files
#regressdir=pathname + '/data/tutorial/VLA/Polcal/
regressdir = './'
scriptprefix = 'polcal_20080224_cband_regression'

#=====================================================================
#
# Clean up old files
os.system('rm -rf '+prefix+'*')

print('Regression Script for VLA POLCAL C-Band Data')
print('Will do: import, flagging, calibration, imaging, analysis')
print('')

#=====================================================================

# Import data from export or use already existing MS?  Or UVFITS?
importmode = 'vla'               # 'vla','fits','ms'
# This is the name of the datafile used in import
# or the name of a previously made ms that will be copied
# NOTE: if an ms name must be different than prefix + '.ms'
#datafile = 'polcal_20080224.cband.edited.ms'
#datafile = '20080224C.UVF'
#
# NOTE: This file may be obtained from the CASA repository:
# http://casa.nrao.edu/Data/VLA/Polcal/POLCA_20080224_1
datafile = ['POLCA_20080224_1']
datafile = pathname + '/data/regression/polcal/20080224_cband/POLCA_20080224_1'

#
# If from export set these:
exportproject = 'POLCA'
exportband = 'C'
#
# Spectral windows to use in ms (usually 0,1)
usespw = ''
usespwlist = ['0','1']

# The ms will have this name
msfile = prefix + '.ms'

# These are names of calibration tables
gtable = prefix + '.gcal'
ftable = prefix + '.fluxscale'
ptable = prefix + '.pcal'
xtable = prefix + '.polx'

# Flagging:
#myquackinterval = 14.0        # if >0 then quack scan beginnings
myquackinterval = 20.0        # if >0 then quack scan beginnings
# Flagging these antennas (if blank then no flagging)
# NOTE: This script uses NEW names, so VLA ants are VAxx
flagants = 'EA04'
#flagants = 'EA*'             # keep only VLA antennas
#flagants = 'VA*'             # keep only EVLA antennas

#
# List of sources in ms
#
#         0    A    1924-292      19:24:51.06      -29.14.30.12  J2000
#         1    A    1743-038      17:43:58.86      -03.50.04.62  J2000
#         2    A    2202+422      22:02:43.29      +42.16.39.98  J2000
#         3    A    2253+161      22:53:57.75      +16.08.53.56  J2000
#         4    B    2136+006      21:36:38.59      +00.41.54.21  J2000
#         5    B    0137+331      01:37:41.30      +33.09.35.13  J2000
#         6    A    2355+498      23:55:09.46      +49.50.08.34  J2000
#         7    B    0319+415      03:19:48.16      +41.30.42.10  J2000
#         8    B    0359+509      03:59:29.75      +50.57.50.16  J2000
#
# These sources are the gain calibrators
gaincalfield = ['0137+331','2202+422','1743-038','1924-292','2136+006','2253+161','2355+498','0319+415','0359+509']
#
# These sources will have calibration transferred from srclist
targets = []

# Assemble field strings from lists
fieldgain = ','.join(gaincalfield)
transgain = ','.join(gaincalfield[1:])    # gain transfer to all but 0137+331
fieldtargets = ','.join(targets)

#
# This list is used for final clean and stats
srclist = gaincalfield + targets

# Location of Cal Models
# e.g. for MacOSX
#fluxcaldir = '/opt/casa/data/nrao/VLA/CalModels/'
# or standard distro
fluxcaldir = pathname + '/data/nrao/VLA/CalModels/'
# or in place
#fluxcaldir = './'

# Calibration parameters:
fluxcalfield = '0137+331'    # primary calibrator for setjy
fluxcalmodel = '3C48_C.im'   # if non-blank use this model image
gaincalfield = ''            # names of gain calibrators (''=all fields)
calrefant = 'VA15'           # reference antenna name for calibration (VA15,EA19)
gainsolint = 20.0            # 20s for gaincal solutions
polcalfield = '2202+422'     # polarization (D-term) calibrator
polcalmode = 'D+QU'          # polarization (D-term) calibration mode
polduvrange = ''             # uvrange for polcal D
setpolmodel = True           # if true then use setjy to set pol model
polxfield = '0137+331'       # polarization angle (X) calibrator
polxuvrange = ''             # uvrange for polcal X
#
setjymode = 'set'            # mode for fluxcal setyjy: 'set', 'flux', 'ft'

# This is the name of the split file for corrected data
srcsplitms = prefix + '.split.ms'

#
# Set up general clean parameters

# This is BnC-config VLA 6cm (4.85GHz) obs
# Check the observational status summary
# Primary beam FWHM = 45'/f_GHz = 557"
# Synthesized beam for VLA/EVLA at C-Band:
#      A-config FWHM = 0.4"
#      B-config FWHM = 1.2"
#      C-config FWHM = 3.9"
#      D-config FWHM = 14.0"
# RMS in 10min (600s) = 0.06 mJy (thats now, but close enough)
#
# Set the output image size and cell size (arcsec)
# 0.4" will give 3x oversampling at least
# clean will say to use a composite integer (e.g.288) for efficiency
#clnalg = 'clark'
clnalg = 'hogbom'
usecsclean = False
clnimsize = 288
clncell = 0.4
# Fix maximum number of iterations
clniter = 200
# Also set flux residual threshold (0.04 mJy)
# Our scans are around 120s
# With rms of 0.06 mJy in 600s ==> rms = 0.13 mJy
# Set to 10x thermal rms
clthreshold = 1.3

# Set up a clean box in the center (1/8 of image)
clncenter = clnimsize/2
clnblc = clncenter - clnimsize/8
clntrc = clncenter + clnimsize/8
# For poor uv coverage, use tigher box (6 x SynthBeam = 18xcell)
clnblc = clncenter - 10
clntrc = clncenter + 10
centerbox = [clnblc,clnblc,clntrc,clntrc]

myclnbox = centerbox
# Can also force interactive cleaning
#myclnbox = 'interactive'

#
#=====================================================================
#
# Polarization of X angle calibrator 0137+331
# If setpolmodel = True
#
# Set up fluxcalmodel
#
fcalmodel = {}
#
# The flux model for 0137+331 (C-band)
fcalfield = {}
# NOTE: you must have entries for all spw in usespwlist
# I,Q,U,V
fcalfield['0'] = [5.405,0,0,0]
fcalfield['1'] = [5.458,0,0,0]
fcalmodel['0137+331'] = fcalfield
# Put in 2202+422
# These values from AIPS (http://www.vla.nrao.edu/astro/calib/polar/2004/)
fcalfield = {}
fcalfield['0'] = [2.465,0,0,0]
fcalfield['1'] = [2.461,0,0,0]
fcalmodel['2202+422'] = fcalfield
#
# Set up pcalmodel
#
pcalmodel = {}
#
# The polarization model for 0137+331
pcalfield = {}
# NOTE: you must have entries for all spw in usespwlist
# From calibrator manual: C-band RLPD=-148deg P/I=0.041
# IPOL,FPOL,RLPHASE
pcalfield['0'] = [5.405,0.041,-148.0]
pcalfield['1'] = [5.458,0.041,-148.0]
pcalmodel['0137+331'] = pcalfield
# Put in 2202+422 (with effective flux of 1.0 before fluxscale)
# These values from AIPS (http://www.vla.nrao.edu/astro/calib/polar/2004/)
pcalfield = {}
pcalfield['0'] = [1.0,0.072,-55.00]
pcalfield['1'] = [1.0,0.072,-55.00]
pcalmodel['2202+422'] = pcalfield
#
# Set the polmodel from pcalmodel
#
print('--Setting up Polarization models--')

polmodel = {}
for field in list(pcalmodel.keys()) :
    spwmodel = {}
    # the RLPD is atan2(U,Q) so Q=I*P/I*cos(RLPD)  U=I*P/I*sin(RLPD)
    for spw in usespwlist:
        ipol = pcalmodel[field][spw][0]
        fpol = pcalmodel[field][spw][1]
        rlpd_deg = pcalmodel[field][spw][2]
        rlpd = rlpd_deg*pl.pi/180.0
        ppol = ipol*fpol
        qpol = ppol*cos(rlpd)
        upol = ppol*sin(rlpd)
        fluxdensity=[ipol,qpol,upol,0.0]

        pmodel = {}
        pmodel['rlpd_deg'] = rlpd_deg
        pmodel['rlpd'] = rlpd
        pmodel['fpol'] = fpol

        fmodel = {}
        fmodel['flux'] = fluxdensity
        fmodel['poln'] = pmodel
        spwmodel[spw] = fmodel

    polmodel[field] = spwmodel

print("Created polmodel dictionary")
print(polmodel)
#
#=====================================================================
# Start processing
#=====================================================================
#
if benchmarking:
    startTime=time.time()
    startProc=time.clock()

#
#=====================================================================
# Data Import and List
#=====================================================================
#
if ( importmode == 'vla' ):
    #
    # Import the data from VLA Export to MS
    #
    print('--ImportVLA--')
    print("Use importvla to read VLA Export and make an MS")

    importvla(archivefiles=datafile,
              vis=msfile,
              bandname=exportband,
              autocorr=False,
              antnamescheme='new',
              project=exportproject)

elif ( importmode == 'fits' ):
    #
    # Import the data from VLA Export to MS
    #
    print('--ImportUVFITS--')
    print("Use importuvfits to read UVFITS and make an MS")

    importuvfits(fitsfile=datafile,
                 vis=msfile,
                 async=False)

else:
    #
    # Copy from msfile
    #
    print('--MS Copy--')
    print("Copying "+datafile+" to "+msfile)
    os.system('cp -r '+datafile+' '+msfile)


if benchmarking:
    import2time=time.time()

#
#=====================================================================
#
print('--Listobs--')

print("List summary of MS")

listobs(vis=msfile)

###############################################
###  Begin Task: listobs  ###
#
# MeasurementSet Name:
#     /home/sandrock/smyers/Testing/2008-03/polcal_20080224/polcal_20080224.cband.raw.ms
# MS Version 2
#
#          Observer: unavailable     Project: POLCA
#       Observation: VLA
#   Data records: 318708       Total integration time = 9836.67 seconds
#          Observed from   17:10:52   to   19:54:48
#
#          ObservationID = 0         ArrayID = 0
#         Date        Timerange                Scan  FldId FieldName      SpwIds
#         24-Feb-2008/17:10:51.7 - 17:12:08.3     1      0 1924-292       [0, 1]
#                     17:21:01.7 - 17:22:18.3     2      1 1743-038       [0, 1]
#                     17:34:31.7 - 17:35:48.3     3      2 2202+422       [0, 1]
#                     17:45:01.7 - 17:46:18.3     4      3 2253+161       [0, 1]
#                     17:55:11.7 - 17:56:28.3     5      4 2136+006       [0, 1]
#                     18:08:01.7 - 18:09:18.3     6      5 0137+331       [0, 1]
#                     18:22:11.7 - 18:23:58.3     7      6 2355+498       [0, 1]
#                     18:32:51.7 - 19:07:58.3     8      2 2202+422       [0, 1]
#                     19:20:51.7 - 19:22:18.3     9      5 0137+331       [0, 1]
#                     19:32:11.7 - 19:33:48.3    10      7 0319+415       [0, 1]
#                     19:42:01.7 - 19:43:18.3    11      8 0359+509       [0, 1]
#                     19:53:31.7 - 19:54:48.3    12      2 2202+422       [0, 1]
#   Fields: 9
#         ID   Code Name          Right Ascension  Declination   Epoch
#         0    A    1924-292      19:24:51.06      -29.14.30.12  J2000
#         1    A    1743-038      17:43:58.86      -03.50.04.62  J2000
#         2    A    2202+422      22:02:43.29      +42.16.39.98  J2000
#         3    A    2253+161      22:53:57.75      +16.08.53.56  J2000
#         4    B    2136+006      21:36:38.59      +00.41.54.21  J2000
#         5    B    0137+331      01:37:41.30      +33.09.35.13  J2000
#         6    A    2355+498      23:55:09.46      +49.50.08.34  J2000
#         7    B    0319+415      03:19:48.16      +41.30.42.10  J2000
#         8    B    0359+509      03:59:29.75      +50.57.50.16  J2000
#   Spectral Windows:  (2 unique spectral windows and 1 unique polarization setups)
#         SpwID  #Chans Frame Ch1(MHz)    ChanWid(kHz)TotBW(kHz)  Ref(MHz)    Corrs
#         0    1 TOPO  4885.1      50000       50000       4885.1      RR  RL  LR  LL
#         1    1 TOPO  4835.1      50000       50000       4835.1      RR  RL  LR  LL
#   Feeds: 27: printing first row only
#         Antenna   Spectral Window     # Receptors    Polarizations
#         1         -1                  2              [         R, L]
#   Antennas: 27:
#         ID   Name  Station   Diam.    Long.         Lat.
#         0    EA24  VLA:W12   25.0 m   -107.37.37.4  +33.53.44.2
#         1    EA16  VLA:W6    25.0 m   -107.37.15.6  +33.53.56.4
#         2    EA01  VLA:W10   25.0 m   -107.37.28.9  +33.53.48.9
#         3    EA19  VLA:W4    25.0 m   -107.37.10.8  +33.53.59.1
#         4    VA08  VLA:W16   25.0 m   -107.37.57.4  +33.53.33.0
#         5    EA17  VLA:W14   25.0 m   -107.37.46.9  +33.53.38.9
#         6    VA06  VLA:W8    25.0 m   -107.37.21.6  +33.53.53.0
#         7    VA22  VLA:W2    25.0 m   -107.37.07.4  +33.54.00.9
#         8    EA04  UNKNOWN   25.0 m   -107.37.41.3  +33.53.42.0
#         9    VA20  VLA:E12   25.0 m   -107.36.31.7  +33.53.48.5
#         10   VA15  VLA:E4    25.0 m   -107.37.00.8  +33.53.59.7
#         11   VA28  VLA:E6    25.0 m   -107.36.55.6  +33.53.57.7
#         12   VA10  VLA:E8    25.0 m   -107.36.48.9  +33.53.55.1
#         13   EA14  VLA:E16   25.0 m   -107.36.09.8  +33.53.40.0
#         14   EA11  VLA:E10   25.0 m   -107.36.40.9  +33.53.52.0
#         15   VA03  VLA:E14   25.0 m   -107.36.21.3  +33.53.44.5
#         16   EA23  VLA:E18   25.0 m   -107.35.57.2  +33.53.35.1
#         17   EA21  VLA:E2    25.0 m   -107.37.04.4  +33.54.01.1
#         18   VA12  VLA:N4    25.0 m   -107.37.06.5  +33.54.06.1
#         19   VA02  VLA:N20   25.0 m   -107.37.13.2  +33.55.09.5
#         20   EA13  VLA:N16   25.0 m   -107.37.10.9  +33.54.48.0
#         21   EA26  VLA:N32   25.0 m   -107.37.22.0  +33.56.33.6
#         22   EA25  VLA:N24   25.0 m   -107.37.16.1  +33.55.37.7
#         23   VA09  VLA:N8    25.0 m   -107.37.07.5  +33.54.15.8
#         24   EA18  VLA:N12   25.0 m   -107.37.09.0  +33.54.30.0
#         25   VA07  VLA:N36   25.0 m   -107.37.25.6  +33.57.07.6
#         26   VA27  VLA:N28   25.0 m   -107.37.18.7  +33.56.02.5
#
#
#       Tables:
#          MAIN                  318708 rows
#          ANTENNA                   27 rows
#          DATA_DESCRIPTION           2 rows
#          DOPPLER                    2 rows
#          FEED                      27 rows
#          FIELD                      9 rows
#          FLAG_CMD             <empty>
#          FREQ_OFFSET         <absent>
#          HISTORY                    6 rows
#          OBSERVATION                1 row
#          POINTING             <empty>
#          POLARIZATION               1 row
#          PROCESSOR            <empty>
#          SOURCE                     9 rows
#          SPECTRAL_WINDOW            2 rows
#          STATE                <empty>
#          SYSCAL              <absent>
#          WEATHER             <absent>
#
###  End Task: listobs  ###
###############################################

# Note that the antennas are out of order as loaded by importvla

if benchmarking:
    list2time=time.time()

#
#=====================================================================
# Data Flagging if needed
#=====================================================================
#
if ( myquackinterval > 0.0 ):
    #
    # First quack the data
    #
    print('--Flagdata (scan starts)--')
    print("Quacking scan beginnings using interval "+str(myquackinterval))
    flagdata(vis=msfile,
             correlation='',
             field='',
             antenna='',
             spw=usespw,
             mode='quack',
             quackinterval=myquackinterval)

#
if (flagants != '' and not flagants.isspace() ):
    print('--Flagdata (antennas)--')
    print("Flag all data to AN "+flagants)

    flagdata(vis=msfile,
             correlation='',
             field='',
             spw=usespw,
             mode='manual',
             antenna=flagants)

flagtimes = '19:06:50~19:06:57,19:21:17~19:21:20'
#
if (flagtimes != '' and not flagtimes.isspace() ):
    print('--Flagdata (timerange)--')
    print("Flag timeranges "+flagtimes)

    flagdata(vis=msfile,
             correlation='',
             field='',
             spw=usespw,
             mode='manual',
             antenna='',
             timerange=flagtimes)

if benchmarking:
    flag2time=time.time()

#
#=====================================================================
# Calibration
#=====================================================================
#
# Set the fluxes of the primary calibrator(s)
#
if ( setjymode == 'flux' ):
    print('--Setjy--')

    print("Use setjy to set flux of "+fluxcalfield+" to point model")

    # Loop over spw
    for spw in usespwlist:
        print("Setting SPW "+spw+" to "+str(fluxdensity))
        setjy(vis=msfile,
              field=fluxcalfield,
              spw=usespw,
              modimage=fluxcaldir+fluxcalmodel,  # If we need a model then put this here
              standard='Perley-Taylor 99',        # enforce older standard
              scalebychan=False,                  # old default
              fluxdensity=fcalmodel[fluxcalfield][spw])


elif ( setjymode == 'ft' ):
    print('--FT--')

    for spw in usespwlist:
        model = fluxcaldir + fluxcalmodel+'_'+spw+'_IQUV.model'
        print("Use FT to set model "+model)
        ft(vis=msfile,
           field=fluxcalfield,
           model=model)

else:
    print('--Setjy--')
    default('setjy')

    print("Use setjy to set flux of "+fluxcalfield)

    setjy(vis=msfile,
          field=fluxcalfield,
          spw=usespw,
          modimage=fluxcaldir+fluxcalmodel,
          standard='Perley-Taylor 99',
          scalebychan=False)

    #
    # You should see something like this in the logger and casa.log file:
    #
    # 0137+331  spwid=  0  [I=5.405, Q=0, U=0, V=0] Jy, (Perley-Taylor 99)
    # 0137+331  spwid=  1  [I=5.458, Q=0, U=0, V=0] Jy, (Perley-Taylor 99)

    # cf. AIPS
    #  SETJY     '0137+331        ' IF =  1 FLUX = 5.4054 (Jy calcd)
    #  SETJY     '0137+331        ' IF =  2 FLUX = 5.4585 (Jy calcd)

    print("Look in logger for the fluxes (should be 5.405 and 5.458 Jy)")

if benchmarking:
    setjy2time=time.time()

#=====================================================================
#
# Initial gain calibration
#
print('--Gaincal--')

print("Solve for antenna gains on sources "+gaincalfield)
print("We have 2 single-channel continuum spw")
print("Output gain table name is "+gtable)
print("Calibrating using fields "+field)

gaincal(vis=msfile,
        caltable=gtable,
        field=fieldgain,
        spw=usespw,
        parang=False,
        gaintype='G',
        solint=gainsolint,
        calmode='ap',
        refant=calrefant,
        minsnr=3)

# use plotcal to view or listcal to list

if benchmarking:
    gaincal2time=time.time()

#=====================================================================
#
# List gain calibration
#
print('--Listcal--')

listfile=caltable + '.list'
print("Listing calibration to file "+listfile)

listcal(vis=msfile,
        caltable=gtable,
        listfile=listfile)

if benchmarking:
    listgcal2time=time.time()

#
#=====================================================================
#
# Bootstrap flux scale
#
print('--Fluxscale--')
print("Use fluxscale to rescale gain table to make new one")

ftable = prefix + '.fluxscale'
print("Output scaled gain cal table is "+ftable)

fluxscale(vis=msfile,
          fluxtable=ftable,
          caltable=gtable,
          reference=fluxcalfield,
          transfer=transgain)

# You should see in the logger something like:
# Found reference field(s): 0137+331
# Found transfer field(s): 1924-292 1743-038 2202+422 2253+161 2136+006 2355+498 0319+415 0359+509
# Flux density for 1924-292 in SpW=0 is: 8.25145 +/- 0.00988121 (SNR = 835.065, nAnt= 13)
# Flux density for 1924-292 in SpW=1 is: 8.22457 +/- 0.0140951 (SNR = 583.505, nAnt= 13)
# Flux density for 1743-038 in SpW=0 is: 5.31336 +/- 0.00603626 (SNR = 880.239, nAnt= 13)
# Flux density for 1743-038 in SpW=1 is: 5.3184 +/- 0.00480634 (SNR = 1106.54, nAnt= 13)
# Flux density for 2202+422 in SpW=0 is: 2.46545 +/- 0.00335055 (SNR = 735.833, nAnt= 13)
# Flux density for 2202+422 in SpW=1 is: 2.46072 +/- 0.00353799 (SNR = 695.512, nAnt= 13)
# Flux density for 2253+161 in SpW=0 is: 8.74607 +/- 0.0142334 (SNR = 614.474, nAnt= 13)
# Flux density for 2253+161 in SpW=1 is: 8.77219 +/- 0.0102289 (SNR = 857.587, nAnt= 13)
# Flux density for 2136+006 in SpW=0 is: 9.97863 +/- 0.013815 (SNR = 722.303, nAnt= 13)
# Flux density for 2136+006 in SpW=1 is: 9.99001 +/- 0.0170089 (SNR = 587.339, nAnt= 13)
# Flux density for 2355+498 in SpW=0 is: 1.29395 +/- 0.00181169 (SNR = 714.221, nAnt= 13)
# Flux density for 2355+498 in SpW=1 is: 1.29893 +/- 0.00217214 (SNR = 597.995, nAnt= 13)
# Flux density for 0319+415 in SpW=0 is: 13.5742 +/- 0.0221722 (SNR = 612.218, nAnt= 13)
# Flux density for 0319+415 in SpW=1 is: 13.5481 +/- 0.0230828 (SNR = 586.932, nAnt= 13)
# Flux density for 0359+509 in SpW=0 is: 5.13982 +/- 0.00906505 (SNR = 566.993, nAnt= 13)
# Flux density for 0359+509 in SpW=1 is: 5.10322 +/- 0.00990264 (SNR = 515.339, nAnt= 13)
# Storing result in polcal_20080224.cband.vla_3c84.fluxscale
# Writing solutions to table: polcal_20080224.cband.vla_3c84.fluxscale

if benchmarking:
    fluxscale2time=time.time()

#=====================================================================
#
# List fluxscale table
#
print('--Listcal--')
print("Listing calibration to file "+listfile)

listcal(vis=msfile,
        caltable=ftable,
        listfile=caltable + '.list')

if benchmarking:
    listfcal2time=time.time()

#=====================================================================
#
# Plot final gain calibration
#
print('--Plotcal--')

figfile = caltable + '.plot.amp.png'
print("Plotting calibration to file "+figfile)

plotcal(caltable=ftable,
        xaxis='time',
        yaxis='amp',
        showgui=False,
        figfile=figfile)

figfile = caltable + '.plot.phase.png'
print("Plotting calibration to file "+figfile)

plotcal(caltable=ftable,
        xaxis='time',
        yaxis='phase',
        showgui=False,
        figfile=figfile)

figfile = caltable + '.plot.antamp.png'
print("Plotting calibration to file "+figfile)

plotcal(caltable=ftable,
        xaxis='antenna',
        yaxis='amp',
        showgui=False,
        figfile=figfile)

if benchmarking:
    plotcal2time=time.time()

dosetpoljy=False
if ( setpolmodel and polcalmode.count('X') > 0 ):
    dosetpoljy=True
    #
    # =====================================================================
    #
    # Now run setjy to (re)set model for polxfield
    #
    print('--Setjy--')
    default('setjy')

    print("Use setjy to set IQU fluxes of "+polxfield)

    for spw in usespwlist:
        setjy(vis=msfile,
              field=polxfield,
              scalebychan=False,
              fluxdensity=polmodel[field][spw]['flux'])

if benchmarking:
    setpoljy2time=time.time()

#=====================================================================
#
# Polarization (D-term) calibration
#
print('--PolCal--')
default('polcal')

print("Polarization D-term Calibration (linear approx) on "+polcalfield)

ptable = prefix + '.pcal'

polcal(vis=msfile,
       caltable=ptable,
       field=polcalfield,
       spw=usespw,
       selectdata=True,
       uvrange=polduvrange,
       poltype=polcalmode,
       solint='inf',
       refant=calrefant,
       minsnr=3,
       gaintable=gtable,
       gainfield='nearest')

# You should see something like:
# Fractional polarization solution for 2202+422 (spw = 0):
# : Q = 0.00356182, U = 0.0717148  (P = 0.0718032, X = 43.5783 deg)
# Fractional polarization solution for 2202+422 (spw = 1):
# : Q = -0.00561314, U = -0.0720833  (P = 0.0723015, X = -47.2263 deg)

if benchmarking:
    polcal2time=time.time()

#=====================================================================
#
# List polcal solutions
#
print('--Listcal--')

listfile = caltable+'.list'
print("Listing calibration to file "+listfile)

listcal(vis=msfile,
        caltable=ptable,
        listfile=listfile)

if benchmarking:
   listpcal2time=time.time()

#=====================================================================
#
# Plot polcal solutions
#
print('--Plotcal--')

figfile = caltable + '.plot.reim.png'
print("Plotting calibration to file "+figfile)
plotcal(caltable=ptable,
        xaxis='real',
        yaxis='imag',
        showgui=False,
        figfile=figfile)

figfile = caltable + '.plot.antamp.png'
print("Plotting calibration to file "+figfile)
plotcal(caltable=ptable,
        xaxis='antenna',
        yaxis='amp',
        showgui=False,
        figfile=figfile)

figfile = caltable + '.plot.antphase.png'
print("Plotting calibration to file "+figfile)
plotcal(caltable=ptable,
        xaxis='antenna',
        yaxis='phase',
        showgui=False,
        figfile=figfile)

figfile = caltable + '.plot.antsnr.png'
print("Plotting calibration to file "+figfile)
plotcal(caltable=ptable,
        xaxis='antenna',
        yaxis='snr',
        showgui=False,
        figfile=figfile)

if benchmarking:
    plotpcal2time=time.time()

#=====================================================================
# Do Chi (X) pol angle calibration if possible
#=====================================================================
#
dopolx = False
if ( polxfield in pcalmodel ):
    dopolx = True

    if ( setpolmodel and not polcalmode.count('X') > 0 ):
        #
        # =============================================================
        #
        # Now run setjy if we havent already
        #

        print('--Setjy--')
        print("Use setjy to set IQU fluxes of "+polxfield)

        for spw in usespwlist:
            setjy(vis=msfile,
                  field=polxfield,
                  scalebychan=False,
                  standard='manual',
                  fluxdensity=polmodel[polxfield][spw]['flux'])

    if benchmarking:
        setxjy2time=time.time()

    #
    # =====================================================================
    #
    # Polarization (X-term) calibration
    #
    print('--PolCal--')
    default('polcal')

    print("Polarization R-L Phase Calibration (linear approx)")
    xtable = prefix + '.polx'

    polcal(vis=msfile,
           caltable=xtable,
           field=polxfield,
           spw=usespw,
           selectdata=True,
           uvrange=polxuvrange,
           poltype='Xf',
           solint='inf',
           refant=calrefant,
           minsnr=3,
           gaintable=[gtable,ptable],
           gainfield=['nearest',''])

    # You should get something like:
    # Position angle offset solution for 0137+331 (spw = 0) = 72.437 deg.
    # Position angle offset solution for 0137+331 (spw = 1) = -21.0703 deg.

    if benchmarking:
        xpolcal2time=time.time()

    #
    # =====================================================================
    #
    # List polcal solutions
    #
    #print '--Listcal--'

    #listfile = caltable + '.list'

    #print "Listing calibration to file "+listfile

    #listcal()

    #
    # =====================================================================
    #
    # Plot polcal solutions
    #
    print('--Plotcal--')

    figfile = caltable + '.plot.png'
    print("Plotting calibration to file "+figfile)

    plotcal(caltable=xtable,
            xaxis='antenna',
            yaxis='phase',
            showgui=False,
            figfile=figfile)


    if benchmarking:
        plotxcal2time=time.time()

else:
    if (polxfield != '' and not polxfield.isspace() ):
        print("DO NOT HAVE PCALMODEL FOR "+polxfield)
        print("PCALMODEL = ",pcalmodel)

    if benchmarking:
        setxjy2time=time.time()
        xpolcal2time=time.time()
        plotxcal2time=time.time()


#=====================================================================
#
# Correct the data
# (This will put calibrated data into the CORRECTED_DATA column)
#
# First using gaincalfield
#
print('--ApplyCal--')
default('applycal')

print("This will apply the calibration to the DATA")
print("Fills CORRECTED_DATA")

# Start with the fluxscaled G table, the D table, and the X table
if (dopolx):
    gaintable = [ftable,ptable,xtable]
else:
    gaintable = [ftable,ptable]

print("Applying calibration to all fields.")

applycal(vis=msfile,
         parang=True,
         gaintable=gaintable,
         gainfield=['nearest','',''])

if benchmarking:
    correct2time=time.time()

#
#=====================================================================
#
# Now write out the corrected data
#
print('--Split--')
default('split')

srcsplitms=prefix + '.split.ms'
print("Split CORRECTED_DATA into DATA in new ms "+srcsplitms)

split(vis=msfile,
      outputvis=srcsplitms,
      datacolumn='corrected')

if benchmarking:
    split2time=time.time()


#
#=====================================================================
# CLEAN the sources
#=====================================================================

clnmodel = {}
#
#=====================================================================
# Loop over sources and spw
# Set up for new clean in patch 2
#

imagermode=''
if usecsclean:
    imagermode='csclean'

for src in srclist:

    srcmodel = {}

    for spwid in usespwlist:

        print('-- Clean '+src+' spw '+spwid+' --')
        default('clean')

        imname1 = prefix + '.' + src + '.' + spwid + '.clean'
        print("  Output images will be prefixed with "+imname1)

        clean(vis=srcsplitms,
              imagename=imname1,
              field=src,
              spw=spwid,
              mode='mfs',
              stokes='IQUV',
              psfmode=clnalg,
              imagermode=imagermode,
              imsize=[clnimsize,clnimsize],
              cell=[clncell,clncell],
              gain=0.1,
              niter=clniter,
              threshold=clthreshold,
              weighting='natural',
              mask=myclnbox)

        # Set up variables
        clnimage1 = imname1+'.image'
        clnmodel1 = imname1+'.model'
        clnresid1 = imname1+'.residual'
        clnmask1  = imname1+'.mask'
        clnpsf1   = imname1+'.psf'
        clnflux1  = imname1+'.flux'

        #
        # =====================================================================
        #
        # Get some statistics of the clean image
        #

        field = src
        spw = spwid

        # Use the clean box
        mybox = str(clnblc)+','+str(clnblc)+','+str(clntrc)+','+str(clntrc)

        spwmodel = {}

        spwstats = {}
        spwfluxes = {}
        spwsum = {}
        spwmod = {}

        for stokes in ['I','Q','U','V']:

            # Use the clean image
            xstat=imstat(imagename=clnimage1,
                         box=mybox,
                         stokes=stokes)

            spwstats[stokes] = xstat

            # Peak (max or min) in box
            xmax = xstat['max'][0]
            xmin = xstat['min'][0]
            if( abs(xmin) > abs(xmax) ):
                xpol = xmin
            else:
                xpol = xmax

            spwfluxes[stokes]= xpol

            # Integrated flux in box
            xsum = xstat['flux'][0]
            spwsum[stokes]= xsum

            # Use the clean model and no box
            xstat=imstat(imagename=clnmodel1,
                         box='',
                         stokes=stokes)

            # Integrated flux in image
            xmod = xstat['sum'][0]
            spwmod[stokes]= xmod

        # Done with stokes

        spwmodel['stat'] = spwstats
        spwmodel['flux'] = spwfluxes
        spwmodel['integ'] = spwsum
        spwmodel['model'] = spwmod

        # Use imval task or ia tool for pixel values in the restored image
        imagename = clnimage1
        # Get image values at the reference pixel
        spwref = {}
        #ia.open(imagename)
        #
        # Stokes I
        # Get reference pixel using imhead in list mode
        myhead = imhead(imagename,'list')
        xref = int(myhead['crpix1'])
        yref = int(myhead['crpix2'])
        #ipix = ia.pixelvalue()
        #iflx = ipix['value']['value']
        refbox = str(xref)+','+str(yref)+','+str(xref)+','+str(yref)
        ipix = imval(imagename,box=refbox,stokes='I')
        iflx = ipix['data'][0]
        spwref['I'] = iflx
        #
        # Stokes Q
        #qpix = ia.pixelvalue([xref,yref,1,0])
        #qflx = qpix['value']['value']
        qpix = imval(imagename,box=refbox,stokes='Q')
        qflx = qpix['data'][0]
        spwref['Q'] = qflx
        #
        # Stokes U
        #upix = ia.pixelvalue([xref,yref,2,0])
        #uflx = upix['value']['value']
        upix = imval(imagename,box=refbox,stokes='U')
        uflx = upix['data'][0]
        spwref['U'] = uflx
        #
        # Stokes V
        #vpix = ia.pixelvalue([xref,yref,3,0])
        #vflx = vpix['value']['value']
        vpix = imval(imagename,box=refbox,stokes='V')
        vflx = vpix['data'][0]
        spwref['V'] = vflx
        #
        # Polarization quantities
        pflx = sqrt( qflx**2 + uflx**2 )
        fflx = pflx/iflx
        xflx = atan2(uflx,qflx)*180.0/pi
        spwref['P'] = pflx
        spwref['F'] = fflx
        spwref['X'] = xflx
        spwref['xref'] = xref
        spwref['yref'] = yref
        #

        # Now the values at the maximum of I
        spwmax = {}
        #
        # Pull the maxpos of I
        xref = spwstats['I']['maxpos'][0]
        yref = spwstats['I']['maxpos'][1]
        refbox = str(xref)+','+str(yref)+','+str(xref)+','+str(yref)
        #
        # Stokes I
        iflx = spwstats['I']['max'][0]
        spwmax['I'] = iflx
        #
        # Stokes Q
        #qpix = ia.pixelvalue([xref,yref,1,0])
        #qflx = qpix['value']['value']
        qpix = imval(imagename,box=refbox,stokes='Q')
        qflx = qpix['data'][0]
        spwmax['Q'] = qflx
        #
        # Stokes U
        #upix = ia.pixelvalue([xref,yref,2,0])
        #uflx = upix['value']['value']
        upix = imval(imagename,box=refbox,stokes='U')
        uflx = upix['data'][0]
        spwmax['U'] = uflx
        #
        # Stokes V
        #vpix = ia.pixelvalue([xref,yref,3,0])
        #vflx = vpix['value']['value']
        vpix = imval(imagename,box=refbox,stokes='V')
        vflx = vpix['data'][0]
        spwmax['V'] = vflx

        spwmax['xref'] = xref
        spwmax['yref'] = yref
        # Done with ia tool
        #ia.close()

        spwmodel['refval'] = spwref
        spwmodel['maxval'] = spwmax

        srcmodel[spwid] = spwmodel

    # Done with spw

    clnmodel[src] = srcmodel

if benchmarking:
    clean2time=time.time()

# Done with srcs
#
if benchmarking:
    endProc=time.clock()
    endTime=time.time()

#
#=====================================================================
# Previous results to be used for regression

regression = {}
regressmodel = {}
regressfile = regressdir + scriptprefix + '.pickle'

try:
    fr = open(regressfile,'r')
except:
    print("No previous regression results file "+regressfile)
else:
    u = pickle.Unpickler(fr)
    regression = u.load()
    fr.close()
    print("Previous regression results filled from "+regressfile)

    if 'results' in regression:
        print("  on "+regression['host']+" for "+regression['version']+" at "+regression['date'])
        regressmodel = regression['results']
    else:
        # Older version
        regressmodel = regression

#=====================================================================
# Report Final Stats
#=====================================================================
#
print('Results for '+prefix+' :')
print("")

import datetime
datestring=datetime.datetime.isoformat(datetime.datetime.today())

outfile = 'out.'+prefix+'.'+datestring+'.log'
logfile=open(outfile,'w')

# Some date and version info
myvers = casalog.version()
myuser = os.getenv('USER')
myhost = str( os.getenv('HOST') )
mycwd = os.getcwd()
myos = os.uname()

# Print version to outfile
print('Running '+myvers+' on host '+myhost, file=logfile)
print('at '+datestring, file=logfile)
print('', file=logfile)

print('Results for '+prefix+' :', file=logfile)
print("", file=logfile)

if ( polxfield in polmodel ):
    # Check RL phase offset on X calibrator
    print("R-L phase residual from image of "+polxfield)
    print("")
    print("R-L phase residual from image of "+polxfield+" :", file=logfile)
    print("", file=logfile)

    src = polxfield
    rlcor = {}

    for spwid in usespwlist:
        ipol = clnmodel[src][spwid]['flux']['I']
        qpol = clnmodel[src][spwid]['flux']['Q']
        upol = clnmodel[src][spwid]['flux']['U']
        vpol = clnmodel[src][spwid]['flux']['V']
        rlpd = atan2(upol,qpol)
        rlpdcal = polmodel[src][spwid]['poln']['rlpd']
        rlpcor = rlpdcal - rlpd
        scor = sin(rlpcor); ccor = cos(rlpcor); rlpcor = atan2(scor,ccor)
        rlcor[spwid] = rlpcor
        rlpcor_deg = rlpcor*180.0/pl.pi

        print("R-L Phase Correction SPW "+spwid+" = %7.2f deg" % rlpcor_deg)
        print("R-L Phase Correction SPW "+spwid+" = %7.2f deg" % rlpcor_deg, file=logfile)


if 'results' in regression:
    print("", file=logfile)
    print("Regression versus "+regression['version']+" on host "+regression['host']+" at "+regression['date'], file=logfile)

#
#=====================================================================
#
# Loop over sources and spw
#
print("")
print("Final Stats:")
print("")

print("", file=logfile)
print("Final Stats:", file=logfile)
print("", file=logfile)

new_regression = {}

# Save info in regression dictionary
new_regression['date'] = datestring
new_regression['version'] = myvers
new_regression['user'] = myuser
new_regression['host'] = myhost
new_regression['cwd'] = mycwd
new_regression['os'] = myos

new_regression['dataset'] = 'VLA POLCAL 20080224'

outpolmodel = {}

passfail = True
for src in srclist:

    print("Source "+src+" :")
    print("Source "+src+" :", file=logfile)

    outpolsrc = {}
    for spwid in usespwlist:

        field = src
        spw = spwid

        # Get fluxes from images

        ipol = clnmodel[src][spwid]['flux']['I']
        qpol = clnmodel[src][spwid]['flux']['Q']
        upol = clnmodel[src][spwid]['flux']['U']
        vpol = clnmodel[src][spwid]['flux']['V']

        # Now get polarization results

        ppol = sqrt(qpol**2 + upol**2)
        fpol = ppol/ipol
        rlpd = atan2(upol,qpol)
        rlpd_deg = rlpd*180.0/pl.pi

        outpolspw = {}
        outpolspw['ipol']=ipol
        outpolspw['ppol']=ppol
        outpolspw['fpol']=fpol
        outpolspw['rlpd']=rlpd
        outpolspw['rlpd_deg']=rlpd_deg

        outpolsrc[spwid] = outpolspw

        #print '  spw %s CASA I = %7.3f Q = %7.3f U = %7.3f V = %7.3f ' %\
        #      (spwid,ipol,qpol,upol,vpol)
        print('  spw %s CASA I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg' %\
              (spwid,ipol,ppol,fpol,rlpd_deg))
        print('  spw %s CASA I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg' %\
              (spwid,ipol,ppol,fpol,rlpd_deg), file=logfile)

        if (src in regressmodel):
            iflx = regressmodel[src][spwid]['ipol']
            fflx = regressmodel[src][spwid]['fpol']
            rlreg = regressmodel[src][spwid]['rlpd']
            rlreg_deg = regressmodel[src][spwid]['rlpd_deg']

            pflx = iflx*fflx
            qflx = pflx*cos(rlreg)
            uflx = pflx*sin(rlreg)
            vflx = 0.0

            print('  spw %s PREV I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg' %\
                  (spwid,iflx,pflx,fflx,rlreg_deg))
            print('  spw %s PREV I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg' %\
                  (spwid,iflx,pflx,fflx,rlreg_deg), file=logfile)

            ipol_diff = ipol - iflx
            ppol_diff = ppol - pflx
            fpol_diff = fpol - fflx
            rldiff = rlpd - rlreg
            rlpd_diff = atan2( sin(rldiff), cos(rldiff) )
            rlpd_diff_deg = rlpd_diff*180.0/pl.pi

            test = (abs(ipol_diff/ipol) < 0.08) & (abs(fpol_diff) < 0.008) & (abs(rlpd_diff < 0.08))
            if test:
                teststr = 'PASS'
            else:
                teststr = 'FAIL'

            passfail = passfail & test

            print('  spw %s DIFF I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg %s ' %\
                  (spwid,ipol_diff,ppol_diff,fpol_diff,rlpd_diff_deg,teststr))
            print('  spw %s PREV I = %7.3f P = %7.3f F = %7.4f X = %7.2f deg %s ' %\
                  (spwid,ipol_diff,ppol_diff,fpol_diff,rlpd_diff_deg,teststr), file=logfile)
            print('')
            print("", file=logfile)

    # Done with spw
    outpolmodel[src] = outpolsrc

    if (src in regressmodel):
        pass
    else:
        print("")
        print("", file=logfile)

# Done with src
if (src in regressmodel):
    if passfail:
        passfailstr = 'PASSED'
    else:
        passfailstr = 'FAILED'

    print('POLCAL Regression '+passfailstr)
    print('POLCAL Regression '+passfailstr, file=logfile)

    print('')
    print('Regression '+passfailstr)
    print('')

new_regression['results'] = outpolmodel

# Should see something like:
#
# R-L phase residual from image of 0137+331
#
# R-L Phase Correction SPW 0 =    0.28 deg
# R-L Phase Correction SPW 1 =    0.33 deg
# No regression results file polcal_20080224_cband_regression.polmodel.regress.pickle
#
# Final Stats:
#
# Source 0137+331 :
#   spw 0 CASA I =   5.263 P =   0.243 F =  0.0462 X = -148.28 deg
#   spw 1 CASA I =   5.313 P =   0.222 F =  0.0418 X = -148.33 deg
#
# Source 2202+422 :
#   spw 0 CASA I =   2.521 P =   0.181 F =  0.0720 X =  -58.04 deg
#   spw 1 CASA I =   2.516 P =   0.184 F =  0.0732 X =  -53.20 deg
#
# Source 1743-038 :
#   spw 0 CASA I =   5.437 P =   0.085 F =  0.0156 X =    6.59 deg
#   spw 1 CASA I =   5.425 P =   0.068 F =  0.0126 X =   -2.26 deg
#
# Source 1924-292 :
#   spw 0 CASA I =   8.079 P =   0.086 F =  0.0106 X =    2.57 deg
#   spw 1 CASA I =   8.012 P =   0.053 F =  0.0066 X =   15.01 deg
#
# Source 2136+006 :
#   spw 0 CASA I =  10.284 P =   0.139 F =  0.0135 X = -160.31 deg
#   spw 1 CASA I =  10.300 P =   0.149 F =  0.0145 X = -167.58 deg
#
# Source 2253+161 :
#   spw 0 CASA I =   8.937 P =   0.492 F =  0.0550 X =   -1.40 deg
#   spw 1 CASA I =   8.905 P =   0.530 F =  0.0595 X =    8.29 deg
#
# Source 2355+498 :
#   spw 0 CASA I =   1.320 P =   0.006 F =  0.0046 X =  140.29 deg
#   spw 1 CASA I =   1.327 P =   0.002 F =  0.0015 X = -130.23 deg
#
# Source 0319+415 :
#   spw 0 CASA I =  13.911 P =   0.066 F =  0.0047 X = -136.14 deg
#   spw 1 CASA I =  13.932 P =   0.024 F =  0.0017 X =  -73.09 deg
#
# Source 0359+509 :
#   spw 0 CASA I =   5.253 P =   0.100 F =  0.0190 X = -126.80 deg
#   spw 1 CASA I =   5.222 P =   0.083 F =  0.0158 X = -128.69 deg
#
#
#=====================================================================
#
# Benchmarking results
#
if benchmarking:
    print('', file=logfile)
    print('********* Benchmarking *****************', file=logfile)
    print('*                                      *', file=logfile)
    print('Total wall clock time was: '+str(endTime - startTime), file=logfile)
    print('Total CPU        time was: '+str(endProc - startProc), file=logfile)
    print('Processing rate MB/s  was: '+str(300./(endTime - startTime)), file=logfile)
    print('* Breakdown:                           *', file=logfile)
    print('*   import       time was: '+str(import2time-startTime), file=logfile)
    print('*   listobs      time was: '+str(list2time-import2time), file=logfile)
    print('*   flagdata     time was: '+str(flag2time-list2time), file=logfile)
    print('*   setjy        time was: '+str(setjy2time-flag2time), file=logfile)
    print('*   gaincal      time was: '+str(gaincal2time-setjy2time), file=logfile)
    print('*   listcal(G)   time was: '+str(listgcal2time-gaincal2time), file=logfile)
    print('*   fluxscale    time was: '+str(fluxscale2time-listgcal2time), file=logfile)
    print('*   listcal(F)   time was: '+str(listfcal2time-fluxscale2time), file=logfile)
    print('*   plotcal(F)   time was: '+str(plotcal2time-listfcal2time), file=logfile)

    if dosetpoljy:
        print('*   setjy(D)     time was: '+str(setpoljy2time-plotcal2time), file=logfile)
        print('*   polcal(D)    time was: '+str(polcal2time-setpoljy2time), file=logfile)
    else:
        print('*   polcal(D)    time was: '+str(polcal2time-plotcal2time), file=logfile)

    print('*   listcal(D)   time was: '+str(listpcal2time-polcal2time), file=logfile)
    print('*   plotcal(D)   time was: '+str(plotpcal2time-listpcal2time), file=logfile)

    if dopolx:
        print('*   setjy(X)     time was: '+str(setxjy2time-plotpcal2time), file=logfile)
        print('*   polcal(X)    time was: '+str(xpolcal2time-setxjy2time), file=logfile)
        print('*   plotcal(X)   time was: '+str(plotxcal2time-xpolcal2time), file=logfile)
        print('*   applycal     time was: '+str(correct2time-plotxcal2time), file=logfile)
    else:
        print('*   applycal     time was: '+str(correct2time-plotpcal2time), file=logfile)

    print('*   split        time was: '+str(split2time-correct2time), file=logfile)
    print('*   clean/stat   time was: '+str(clean2time-split2time), file=logfile)
    print('*****************************************', file=logfile)
    print('sandrock (2008-06-17) wall time was: 255 seconds', file=logfile)
    print('sandrock (2008-06-17) CPU  time was: 233 seconds', file=logfile)

logfile.close()

print('')
if benchmarking:
    print('Total wall clock time was: '+str(endTime - startTime))
    print('Total CPU        time was: '+str(endProc - startProc))
    print('Processing rate MB/s  was: '+str(300./(endTime - startTime)))
    print('')

print("Done with POLCAL Regression")
#
# Done
#
logfile.close()
print("Results are in "+outfile)

#
#=====================================================================
#
# Now save stat dictionaries using Pickle
#
pickfile = 'out.'+prefix + '.polmodel.'+datestring+'.pickle'
f = open(pickfile,'w')
p = pickle.Pickler(f)
# The regression results for this run
p.dump(new_regression)
# Now the clean results and the input pol models
p.dump(clnmodel)
p.dump(polmodel)
f.close()
print("")
print("Dictionaries new_regression,clnmodel,polmodel saved in "+pickfile)
print("")
print("Use Pickle to retrieve these")
print("")

# e.g.
# f = open(pickfile)
# u = pickle.Unpickler(f)
# regression = u.load()
# clnmodel = u.load()
# polmodel = u.load()
# f.close()

print("")
print("Completed POLCAL Regression")
if passfail:
    print("Regression PASSED")
else:
    print("Regression FAILED")

