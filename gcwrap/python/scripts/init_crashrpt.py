###
### Initialize the crash dump feature
###

import os
import signal
import tempfile

casa['state']['crash-reporter'] = False
if ( casa['flags'].crash_report and
     'CASA_USE_CRASH_REPORTER' in os.environ and
     os.environ['CASA_USE_CRASH_REPORTER'].upper( ) == 'TRUE'):
    try:
        temporaryDirectory = tempfile.gettempdir()
        posterApp = casa['helpers']['crashPoster']
        if posterApp is None: posterApp = "" # handle case where it wasn't found
        postingUrl = "https://casa.nrao.edu/cgi-bin/crash-report.pl"
        theLogFile = casa['files']['logfile']
        message = casac.utils()._crash_reporter_initialize(temporaryDirectory, posterApp, postingUrl, theLogFile)
        if len (message) > 0:
            if message != "no-op":
                print(("***\n*** Crash reporter failed to initialize: " + message))
        else:
            casa['state']['crash-reporter'] = True
    except Exception as e:
        print("***\n*** Crash reporter initialization failed.\n***")
        print("*** exception={0}\n***".format (e))
