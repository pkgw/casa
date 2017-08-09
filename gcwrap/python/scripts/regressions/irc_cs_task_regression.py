###############################
#
# IRC+10216 13CS Reduction Script
# using new sd tasks
# Nod data
# 
# tasks used
# sdlistold
# sdcalold
# sdsmoothold
# sdbaselineold
# sdplotold
# sdstatold
# sdsaveold
###############################
import time
import os

os.system('rm -rf IRC+10216_rawACSmod IRC+10216_rawACSmod_cal IRC+10216_rawACSmod_cal_sm IRC+10216_rawACSmod_cal_sm_bs  irc_cs_reducedSCAN0_CYCLE0_BEAM0_IF0.txt irc_cs_reduced.eps irc_cs_fit.txt')

#enable/disable plotting
doplot = False

casapath = os.environ['CASAPATH'].split()[0]
datapath = casapath+'/data/regression/ATST5/IRC+10216/IRC+10216_rawACSmod'
copystring = 'cp -r '+datapath+' .'
os.system(copystring)

startTime = time.time()
startProc = time.clock()

#   Project: AGBT06A_018_01
#Observation: GBT(1 antennas)
#
#Tue Jan 30 20:24:29 2007    NORMAL ms::summary:
#Data records: 4992       Total integration time = 4244.15 seconds
#   Observed from   10:31:01   to   11:41:45
#
#Tue Jan 30 20:24:29 2007    NORMAL ms::summary:
#Fields: 3
#  ID   Name          Right Ascension  Declination   Epoch
#  0    IRC+10216     09:47:57.38      +13.16.43.70  J2000
#  1    IRC+10216     09:47:57.38      +13.16.43.70  J2000
#  2    IRC+10216     09:47:57.38      +13.16.43.70  J2000
#
#Tue Jan 30 20:24:29 2007    NORMAL ms::summary:
#Spectral Windows:  (13 unique spectral windows and 1 unique polarization setups)
#  SpwID  #Chans Frame Ch1(MHz)    Resoln(kHz) TotBW(kHz)  Ref(MHz)    Corrs
#  0        4096 LSRK  48263.5858  12.2067945  49999.0303  48288.5853  RR  LL H2CO
#  1        4096 LSRK  45899.1486  12.2067945  49999.0303  45924.1481  RR  LL H213CO
#  2        4096 LSRK  43857.1942  12.2067945  49999.0303  43882.1937  RR  LL 13CS
#  3        4096 LSRK  48970.0031  12.2067945  49999.0303  48995.0026  RR  LL CS
#  7        4096 LSRK  46226.6452  12.2067949  49999.0319  46251.6447  RR  LL 
#  17       4096 LSRK  45469.0638  12.2067955  49999.0343  45494.0633  RR  LL HC3N
#  18       4096 LSRK  45280.5214  12.2067955  49999.0343  45305.5209  RR  LL HCC13CN
#  19       4096 LSRK  44054.8331  12.2067955  49999.0343  44079.8326  RR  LL CH3OH
#  20       4096 LSRK  44146.1063  12.2067955  49999.0343  44171.1058  RR  LL HCCC15N
#  27       4096 LSRK  43941.6784  12.206796   49999.0365  43966.6779  RR  LL HN13CO
#  28       4096 LSRK  42625.1528  12.206796   49999.0365  42650.1523  RR  LL H15NCO
#  29       4096 LSRK  41574.856   12.206796   49999.0365  41599.8556  RR  LL HNC18O
#  30       4096 LSRK  43402.4488  12.206796   49999.0365  43427.4483  RR  LL SiO

# Scans: 230-231,233-236,253-256 Setup 1 H2CO et al
# Scans: 237-240,249-252 Setup 2 HC3N et al
# Scans: 241-248         Setup 3 SiO et al


asap_init()                             #load ASAP module

					#scan numbers (zero-based) as compared to GBTIDL
					#changes made to get to IRC+10216_rawACSmod
					#  -- merge spectral windows with tolerance

# summary
#default(sdlistold)
#infile = 'IRC+10216_rawACSmod'
#sdlistold()

if doplot:
   localplotlevel = 1
else:
   localplotlevel = 0

# calibartion and averaging
# calibrate nod scans for CS line (IF=3)
default(sdcalold)
infile = 'IRC+10216_rawACSmod'
fluxunit = 'K'
calmode = 'nod'
#scanlist = [230,231]
scan = '230,231'
#iflist = [3]
spw = '3'
tau = 0.09         # do opacity correction
overwrite = True
plotlevel = localplotlevel
sdcalold()
# output
localoutfile = infile+'_cal'

#  averaging and smoothing
default(sdaverageold)
infile = localoutfile
# do time and polarization average
timeaverage = True # average in time
tweight = 'tintsys' # weighted by integ time and Tsys for time averaging
scanaverage = False
polaverage = True  # average polarization
pweight = 'tsys'   # weighted by Tsys for pol. averaging
#smoothing
# do boxcar smoothing with channel width=5
kernel = 'boxcar'
kwidth = 5
overwrite = True
plotlevel = localplotlevel
sdaverageold()
localoutfile = infile+'_ave'

#fit and remove baselines
# do baseline fit with cubic spline with one knot (npiece=2)
# 3-sigma clipping plus 1 iteration applied.
# automatically detect lines to exclude from fitting
default(sdbaselineold)
infile = localoutfile
maskmode = 'auto'
#edge = [50]
thresh = 5
avg_limit = 4
#blfunc = 'poly'
#order = 1
blfunc = 'cspline'
npiece = 2
clipthresh = 3.0
clipniter = 1
overwrite = True
plotlevel = localplotlevel
sdbaselineold()
localoutfile = infile+'_bs'

#plotting the reslut
#plot the spectrum and save to a postscript file
if doplot:
   default(sdplotold)
   infile = localoutfile
   specunit = 'GHz'
   outfile = 'irc_cs_reduced.eps'
   #sd.plotter.set_histogram(hist=True)     # draw spectrum using histogram                 # histogram
   #sd.plotter.axhline(color='r',linewidth=2) # zline                                       # zline
   sdplotold()
else:
   print("Plotting the result is skipped.")

# statistics
default(sdstatold)
# select line free regions to get rms
infile = localoutfile
#masklist = [800,1500]
spw = '*:800~1500'
xstat = sdstatold()
curr_rms = xstat['rms']
#rms=
#
# select the line region
#masklist = [1850,2300]
spw = '*:1850~2300'
xstat = sdstatold()
xstat
curr_max = xstat['max']
curr_sum = xstat['sum']
curr_median = xstat['median']
curr_mean = xstat['mean']

# Save the spectrum
# in different formats
default(sdsaveold)
infile = localoutfile
outfile = 'irc_cs_reduced'
outform = 'ASCII'
overwrite = True
sdsaveold()
#outfile = 'irc_cs_reduced.ms'
#outform = 'MS2'
#sdsaveold()

#
endProc = time.clock()
endTime = time.time()

# --- end of irc cs script

#irc_max=3.3
#irc_rms=0.147
#irc_sum=627.5
# Regression values of CASA 2.3(#6654)+ASAP 2.2.0(#1448)
# on 64bit REL5.2 (2008/12/01)
#irc_max=3.325
#irc_rms=0.1473
#irc_sum=627.9
# Regression values of CASA 3.1+ASAP3 on 64bit RHEL5.5 (2010/10?)
# Regression values of CASA 3.2(#14238)+ASAP3(#2031) on 64bit RHEL5.5
#   using cspline with (npiece,clipthresh,clipniter)=(2,3.0,1) in sdbaseline (2011/03/25 WK)
irc_max = 3.3270
irc_rms = 0.1472
irc_sum = 630.75

diff_max = abs( (irc_max - curr_max) / irc_max )
diff_rms = abs( (irc_rms - curr_rms) / irc_rms )
diff_sum = abs( (irc_sum - curr_sum) / irc_sum )

import datetime
datestring = datetime.datetime.isoformat(datetime.datetime.today())
outfile = 'irc.cs.task'+datestring+'.log'
logfile = open(outfile,'w')

print('', file=logfile)
print('********** Regression ***********', file=logfile)
print('*                               *', file=logfile)
if (diff_max < 0.05): print('* Passed spectrum max test ', file=logfile)
print('*  Spectrum max '+str(curr_max), file=logfile)
if (diff_rms < 0.05): print('* Passed spectrum rms test ', file=logfile)
print('*  Spectrum rms '+str(curr_rms), file=logfile)
if (diff_sum < 0.05): print('* Passed spectrum (line) sum test', file=logfile)
print('*  Line integral '+str(curr_sum), file=logfile)
if ((diff_max<0.05) & (diff_rms<0.05) & (diff_sum<0.05)):
	regstate = True
        print('')
        print('Regression PASSED')
        print('')
        print('---', file=logfile)
        print('Passed Regression test for IRC-CS', file=logfile)
        print('---', file=logfile)
else:
	regstate = False
        print('')
        print('Regression FAILED')
        print('')
        print('----FAILED Regression test for IRC-CS', file=logfile)
print('*********************************', file=logfile)

print('', file=logfile)
print('', file=logfile)
print('********* Benchmarking *****************', file=logfile)
print('*                                      *', file=logfile)
print('Total wall clock time was: '+str(endTime - startTime), file=logfile)
print('Total CPU        time was: '+str(endProc - startProc), file=logfile)
#print 'Processing rate MB/s  was: ', 35.1/(endTime - startTime)

logfile.close()
