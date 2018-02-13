###
### Telemetry shutdown hook
###
from casa_shutdown import add_shutdown_hook
from telemetry import telemetry

def submitStatistics():
    casalog.poststat("Stop CASA")
    print "Submitting telemetry"
    mytelemetry = telemetry(casa['dirs']['rc'] + '/', 'casastats-*.log', casalog)
    mytelemetry.createStampFile()
    if (mytelemetry.isSubmitInterval()):
        mytelemetry.send('https://casa.nrao.edu/cgi-bin/crash-report.pl')
        mytelemetry.refreshStampFile()

casa['state']['telemetry-enabled'] = False

if ( casa['flags'].telemetry or
    (os.environ.has_key('CASA_ENABLE_TELEMETRY') and
     os.environ['CASA_ENABLE_TELEMETRY'].upper( ) == 'TRUE')):
     casa['state']['telemetry-enabled'] = True
     add_shutdown_hook(submitStatistics)
