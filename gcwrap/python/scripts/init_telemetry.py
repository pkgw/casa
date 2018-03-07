###
### Telemetry shutdown hook
###
from casa_shutdown import add_shutdown_hook
from telemetry import telemetry
import __casac__
import fnmatch
import os

logdir = casa['dirs']['rc']
logpattern = 'casastats-*.log'
def submitStatistics():
    casalog.poststat("Stop CASA")
    print "Submitting telemetry"
    mytelemetry = telemetry(logdir + '/', logpattern, casalog)
    mytelemetry.createStampFile()
    if (mytelemetry.isSubmitInterval()):
        mytelemetry.send('https://casa.nrao.edu/cgi-bin/crash-report.pl')
        mytelemetry.refreshStampFile()

casa['state']['telemetry-enabled'] = False

casa_util = __casac__.utils.utils()
rcTelemetryFlag = str.upper(casa_util.getrc("EnableTelemetry"))

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
             print "Matched: " + file
             logfiles.append(file)

     logfiles.sort(reverse=True)

     if (logfiles and logfiles[0] != None):
         print "Found an existing telemetry logfile. Appending telemetry data to " + casa['dirs']['rc'] + "/" + logfiles[0]
         casa['files']['telemetry-logfile'] = casa['dirs']['rc'] + "/" + logfiles[0]
     else :
         print "Creating a new telemetry file"
         casa['files']['telemetry-logfile'] = casa['dirs']['rc'] + '/casastats-' +time.strftime("%Y%m%d-%H%M%S", time.gmtime()) + '.log'

     # Submit statistics at shutdown
     add_shutdown_hook(submitStatistics)

     print "Telemetry initialized."
