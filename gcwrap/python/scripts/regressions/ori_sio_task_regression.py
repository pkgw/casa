###############################
#
# ORION-S SiO Reduction Script
# using SD tasks
# Position-Switched data
#
# used tasks
# sdlistold
# sdcalold
# sdaverageold
# sdbaselineold
# sdplotold
# sdstatold
# sdfitold
# sdsaveold
###############################
#
import time
import os

os.system('rm -rf OrionS_rawACSmod OrionS_rawACSmod_cal OrionS_rawACSmod_cal_sm OrionS_rawACSmod_cal_sm_bs orions_sio_reducedSCAN0_CYCLE0_BEAM0_IF0.txt orions_sio_reduced.eps orions_sio_fit.txt')

#enable/disable plotting
doplot = False

casapath = os.environ['CASAPATH']
datapath = casapath.split()[0]+'/data/regression/ATST5/OrionS/OrionS_rawACSmod'
copystring = 'cp -r '+datapath+' .'
os.system(copystring)

startTime = time.time()
startProc = time.clock()

#           MeasurementSet Name:  /home/rohir3/jmcmulli/SD/OrionS_rawACSmod      MS Version 2
#
# Project: AGBT06A_018_01
# Observation: GBT(1 antennas)
#
#Data records: 256       Total integration time = 1523.13 seconds
#   Observed from   01:45:58   to   02:11:21
#
#Fields: 4
#  ID   Name          Right Ascension  Declination   Epoch
#  0    OrionS        05:15:13.45      -05.24.08.20  J2000
#  1    OrionS        05:35:13.45      -05.24.08.20  J2000
#  2    OrionS        05:15:13.45      -05.24.08.20  J2000
#  3    OrionS        05:35:13.45      -05.24.08.20  J2000
#
#Spectral Windows:  (8 unique spectral windows and 1 unique polarization setups)
#  SpwID  #Chans Frame Ch1(MHz)    Resoln(kHz) TotBW(kHz)  Ref(MHz)    Corrs
#  0        8192 LSRK  45464.3506  6.10423298  50005.8766  45489.3536  RR  LL HC3N
#  1        8192 LSRK  45275.7825  6.10423298  50005.8766  45300.7854  RR  LL HN15CO
#  2        8192 LSRK  44049.9264  6.10423298  50005.8766  44074.9293  RR  LL CH3OH
#  3        8192 LSRK  44141.2121  6.10423298  50005.8766  44166.2151  RR  LL HCCC15N
#  12       8192 LSRK  43937.1232  6.10423356  50005.8813  43962.1261  RR  LL HNCO
#  13       8192 LSRK  42620.4173  6.10423356  50005.8813  42645.4203  RR  LL H15NCO
#  14       8192 LSRK  41569.9768  6.10423356  50005.8813  41594.9797  RR  LL HNC18O
#  15       8192 LSRK  43397.8198  6.10423356  50005.8813  43422.8227  RR  LL SiO

# Scans: 21-24  Setup 1 HC3N et al
# Scans: 25-28  Setup 2 SiO et al

asap_init()                             #load ASAP module
                                        #Orion-S (SiO line reduction only)
                                        #Notes:
                                        #scan numbers (zero-based) as compared to GBTIDL

                                        #changes made to get to OrionS_rawACSmod
                                        #modifications to label sig/ref positions
#os.environ['CASAPATH'] = casapath


# summary
#default(sdlistold)
#infile = 'OrionS_rawACSmod'
#sdlistold()

if doplot:
   localplotlevel = 1
else:
   localplotlevel = 0

# calibartion
# calibrate position-switched CH3OH scans (IF=15)
default(sdcalold)
infile = 'OrionS_rawACSmod'
fluxunit = 'K'
calmode = 'ps'
#scanlist = [25,26,27,28]
scanlist = '25,26,27,28'
#iflist = [15]
spw = '15'
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
tweight = 'tintsys' #weighted by integ time and Tsys for time averaging
scanaverage = False
polaverage = True  # average polarization
pweight = 'tsys'   # weighted by Tsys for pol averaging
#smoothing
# do boxcar smoothing with channel width=5
kernel = 'boxcar'
kwidth = 10
overwrite = True
plotlevel = localplotlevel
sdaverageold()
localoutfile = infile+'_ave'

#fit and remove baselines
# do baseline fit with polynomial order of 2
# automatically detect lines to exclude from fitting
default(sdbaselineold)
infile = localoutfile
maskmode = 'list'
#masklist = [[500,3500],[5000,7500]]
spw = '*:500~3500;5000~7500'
blfunc = 'poly'
order = 5
overwrite = True
plotlevel = localplotlevel
sdbaselineold()
localoutfile = infile+'_bs'
#sd.plotter.plot(spave)                 # plot                                          # baseline

#plotting the reslut
#plot the spectrum and save to a postscript file
if doplot:
   default(sdplotold)
   infile = localoutfile
   specunit = 'GHz'
   outfile = 'orions_sio_reduced.eps'
   #sd.plotter.set_histogram(hist=True)     # draw spectrum using histogram                 # histogram
   #sd.plotter.axhline(color='r',linewidth=2) # zline                                       # zline
   sdplotold()
else:
   print("Plotting the result is skipped.")

# statistics
default(sdstatold)
# select line free regions to get rms
infile = localoutfile
#masklist = [1000,3000]
spw = '*:1000~3000'
xstat = sdstatold()
curr_rms = xstat['rms']
#rms= 0.037199538201093674 [CASA 2.3(#6654)+ASAP 2.2.0(#1448)]
#
# select the line region
#masklist = [3900,4300]
spw = '*:3900~4300'
xstat = sdstatold()
xstat
curr_max = xstat['max']
curr_sum = xstat['sum']
curr_median = xstat['median']
curr_mean = xstat['mean']

# fitting
default(sdfitold)
infile = localoutfile
#sd.plotter.plot(spave)                 # plot spectrum
fitmode = 'list'
#maskline = [3900,4300] # create region around line
spw = '*:3900~4300'     # create region around line
nfit = 1
plotlevel = localplotlevel
outfile = 'orions_sio_fit.txt'
xstat = sdfitold()
xstat  # print fit statistics

# Save the spectrum
# in different formats
default(sdsaveold)
infile = localoutfile
outfile = 'orions_sio_reduced'
outform = 'ASCII'
overwrite = True
sdsaveold()
outfile = 'orions_sio_reduced.ms'
outform = 'MS2'
sdsaveold()

endProc = time.clock()
endTime = time.time()
#
# --- end of orion-s ch3oh script
# Regression

#ori_max=0.366
#ori_rms=0.037
#ori_sum=49.9
# Regression values of CASA 2.3(#6654)+ASAP 2.2.0(#1448)
# on 64bit REL5.2 (2008/12/01)
ori_max = 0.3662
ori_rms = 0.03720
ori_sum = 49.95
diff_max = abs( (ori_max - curr_max) / ori_max )
diff_rms = abs( (ori_rms - curr_rms) / ori_rms )
diff_sum = abs( (ori_sum - curr_sum) / ori_sum )

import datetime
datestring = datetime.datetime.isoformat(datetime.datetime.today())
outfile = 'ori.sio.task.'+datestring+'.log'
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
        print('---', file=logfile)
        print('Passed Regression test for OrionS-SiO', file=logfile)
        print('---', file=logfile)
        print('')
        print('Regression PASSED')
        print('')
else:
        regstate = False
        print('----FAILED Regression test for OrionS-SiO', file=logfile)
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
#print 'Processing rate MB/s  was: ', 35.1/(endTime - startTime)

logfile.close()
