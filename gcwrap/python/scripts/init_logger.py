from mpi4casa.MPIEnvironment import MPIEnvironment


def casalogger(logfile=''):
    """
    Spawn a new casalogger using logfile as the filename.
    You should only call this if the casalogger dies or you close it
    and restart it again.

    Note: if you changed the name of the log file using casalog.setlogfile
    you will need to respawn casalogger with the new log filename. Eventually,
    we will figure out how to signal the casalogger with the new name but not
    for a while.
    """
    from casa_system import procmgr

    if logfile == '':
        if 'files' in casa and 'logfile' in casa['files'] :
            logfile = casa['files']['logfile']
        else:
            casa['files']['logfile'] = os.getcwd( ) + '/casa.log'
            logfile = 'casa.log'

    if (os.uname()[0]=='Darwin'):

        if casa['helpers']['logger'] == 'console':
            procmgr.create("logger",['/usr/bin/open','-a','console', logfile])
        else:
            procmgr.create("logger",[casa['helpers']['logger'],logfile])

    elif (os.uname()[0]=='Linux'):
        procmgr.create("logger",[casa['helpers']['logger'],logfile])

    else:
        print('Unrecognized OS: No logger available')


deploylogger = True

if not os.access('.', os.W_OK):
    if casa['flags'].nologfile == False:
        print()
        print("********************************************************************************")
        print("Warning: no write permission in current directory, no log files will be written.")
        print("********************************************************************************")
    deploylogger = False

if casa['flags'].nologger :
    deploylogger = False

if casa['flags'].nogui :
    deploylogger = False

## do not start logger gui on MPI clients...
if MPIEnvironment.is_mpi_enabled and not MPIEnvironment.is_mpi_client:
    deploylogger = False

if deploylogger and casa['files']['logfile'] != '/dev/null':
    casalogger(casa['files']['logfile'])

if (casa['state']['telemetry-enabled']):
    casalog = casac.logsink(casa['files']['logfile'], True, casa['files']['telemetry-logfile'])
    casatelemetry.setCasaLog(casalog)
    if not MPIEnvironment.is_mpi_enabled or (MPIEnvironment.is_mpi_enabled and MPIEnvironment.is_mpi_client):
        casatelemetry.submitStatistics()
else :
    casalog = casac.logsink(casa['files']['logfile'])


processor_origin = MPIEnvironment.processor_origin
casalog.processorOrigin(processor_origin)

casalog.showconsole((MPIEnvironment.is_mpi_enabled and MPIEnvironment.log_to_console) or casa['flags'].log2term)

###
### For reasons unknown, setglobal causes the global logger to steal the
### log sink associated with this casalog object we've created, making
### our casalog object useless for the purposes of configuration. As a
### result, it's best to call this last.
###
casalog.setglobal(True)

try:
    casalog.post("CASA Version " + casa['build']['version'])
except:
    print("Error: the logfile is not writable")
    sys.exit(1)
