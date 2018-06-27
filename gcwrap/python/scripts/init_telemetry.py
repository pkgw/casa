###
### Telemetry shutdown hook
###
from casa_shutdown import add_shutdown_hook
from telemetry import telemetry
import fnmatch
import os
from casac import casac
import __casac__

def logShutdown():
    if (casa['state']['telemetry-enabled'] == True):
        casalog.poststat("Stop CASA")

casa['state']['telemetry-enabled'] = False

casa_util = __casac__.utils.utils()
rcTelemetryFlag = str.upper(casa_util.getrc("EnableTelemetry"))

if ( casa['flags'].telemetry or
    (os.environ.has_key('CASA_ENABLE_TELEMETRY') and
     os.environ['CASA_ENABLE_TELEMETRY'].upper( ) == 'TRUE') or
     rcTelemetryFlag == 'TRUE'):

     # Enable telemetry
     casa['state']['telemetry-enabled'] = True

     add_shutdown_hook(logShutdown)

     casatelemetry = telemetry(casa)
