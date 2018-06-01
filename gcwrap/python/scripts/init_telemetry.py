###
### Telemetry shutdown hook
###
from casa_shutdown import add_shutdown_hook
from telemetry import telemetry
import __casac__
import fnmatch
import os
from casac import casac
import TelemetryLogMonitor
from casac import casac

telemetryhelper = casac.telemetryhelper()
hostid = telemetryhelper.getUniqueId()

logdir = casa['dirs']['rc']
logpattern = 'casastats-*.log'

def submitStatistics():
    if (casa['state']['telemetry-enabled'] == True):
        casalog.poststat("Stop CASA")
    casalog.post("Checking telemetry submission interval")
    mytelemetry = telemetry(logdir + '/', logpattern, casalog)
    mytelemetry.createStampFile()
    if (mytelemetry.isSubmitInterval()):
        mytelemetry.send('https://casa.nrao.edu/cgi-bin/crash-report.pl')
        mytelemetry.refreshStampFile()

casa['state']['telemetry-enabled'] = False

casa_util = __casac__.utils.utils()
rcTelemetryFlag = str.upper(casa_util.getrc("EnableTelemetry"))

# Size limit for the telemetry logs
tLogSizeLimit = 10000
# File size check interval
tLogSizeInterval = 60
try:
    tLogSizeLimit = int(casa_util.getrc("TelemetryLogLimit"))
    tLogSizeInterval = int(casa_util.getrc("TelemetryLogSizeInterval"))
except:
    pass

if ( casa['flags'].telemetry or
    (os.environ.has_key('CASA_ENABLE_TELEMETRY') and
     os.environ['CASA_ENABLE_TELEMETRY'].upper( ) == 'TRUE') or
     rcTelemetryFlag == 'TRUE'):

     # Enable telemetry
     casa['state']['telemetry-enabled'] = True
     casa['files']['telemetry-logfile'] = "casastats.log"

     # Create a stats file if one doesn't exist. Append to an existing file otherwise.
     logfiles = []

     for file in os.listdir(logdir):
         if fnmatch.fnmatch(file, logpattern):
             #print "Matched: " + file
             logfiles.append(file)

     logfiles.sort(reverse=True)

     # Size of the existing (non-active) logfiles
     inactiveTLogSize = 0

     if (logfiles and logfiles[0] != None):
         print "Found an existing telemetry logfile. Appending telemetry data to " + casa['dirs']['rc'] + "/" + logfiles[0]
         casa['files']['telemetry-logfile'] = casa['dirs']['rc'] + "/" + logfiles[0]
         for i in range(1, len(logfiles)):
             inactiveTLogSize = inactiveTLogSize + os.path.getsize(casa['dirs']['rc'] + "/" + logfiles[i])/1024
         #print "Inactive log size: " + str(inactiveTLogSize)
     else :
         print "Creating a new telemetry file"
         myUtils = casac.utils()
         ver = myUtils.version()
         casaver = "CASA-" + str(ver[0])+ str(ver[1]) + str(ver[2])+ "-" + str(ver[3])
         casa['files']['telemetry-logfile'] = casa['dirs']['rc'] + '/casastats-' + casaver +'-'  + hostid + "-" + time.strftime("%Y%m%d-%H%M%S", time.gmtime()) + '.log'

     # Submit statistics at shutdown
     add_shutdown_hook(submitStatistics)

     # Subtract the inactive log sizes from the total log file size limit
     tLogSizeLimit = tLogSizeLimit - inactiveTLogSize
     if (tLogSizeLimit <= 0):
         print "Telemetry log size limit exceeded. Disabling telemetry."
         casa['state']['telemetry-enabled'] = False
     else :
         tLogMonitor = TelemetryLogMonitor.TelemetryLogMonitor()
         tLogMonitor.start(casa['files']['telemetry-logfile'],tLogSizeLimit, tLogSizeInterval, casa)
         print "Telemetry initialized."
