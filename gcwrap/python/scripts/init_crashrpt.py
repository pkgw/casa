###
### Initialize the crash dump feature
###

import os
import signal
import crashrpt_conf

casa['state']['crash-reporter'] = False
if ( casa['flags'].crash_report or (
     'CASA_USE_CRASH_REPORTER' in os.environ and
     os.environ['CASA_USE_CRASH_REPORTER'].upper( ) == 'TRUE')):
    try:
        systemTempDir = crashrpt_conf.systemTempDir
        temporaryDirectoryCommon = crashrpt_conf.temporaryDirectoryCommon
        temporaryDirectory = crashrpt_conf.temporaryDirectory
        # Create common temporary directory for multi-user environments
        # Configure permissions based on system temporary directory
        if not os.path.exists(temporaryDirectoryCommon):
                os.makedirs(temporaryDirectoryCommon)
                st = os.stat(systemTempDir)
                perm = st.st_mode
                os.chmod(temporaryDirectoryCommon, perm)
        # Create temporary directory for user
        if not os.path.exists(temporaryDirectory):
                os.makedirs(temporaryDirectory)
        posterApp = casa['helpers']['crashPoster']
        if posterApp is None: posterApp = "" # handle case where it wasn't found
        postingUrl = "https://casa.nrao.edu/cgi-bin/crash-report.pl"
        if 'CASA_CRASHREPORT_URL' in os.environ :
                postingUrl = os.environ['CASA_CRASHREPORT_URL']
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
