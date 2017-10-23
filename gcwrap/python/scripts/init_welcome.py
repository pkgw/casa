import atexit
import os
import sys
import traceback
import platform
import datetime

from init_welcome_helpers import redirect_argv, immediate_exit_with_handlers

if (casa['state']['telemetry-enabled'] == True):
    casa['state']['telemetry-starttime'] = str(datetime.datetime.now())
    casalog.origin("CASAStart")
    casalog.poststat("Starting CASA at: " + casa['state']['telemetry-starttime'] + " Version " + casa['build']['version'] + " Platform: " + platform.platform() +  " Variant: " + casa['variant'])
    # Set back to "casa" so that the current logging is not altered
    casalog.origin("")

if casa['flags'].execute:
    import os.path

    if '/' in casa['flags'].execute[0]:
        ## qualified path
        __paths_to_check = [ '' ]
    else:
        ## non-qualified path
        __paths_to_check = [ "./", casa['dirs']['python'] + '/' ]

    __candidates = filter( os.path.isfile, map(lambda dir: dir + casa['flags'].execute[0], __paths_to_check) )

    if len(__candidates) > 0:
        # Run file with filename given in the command line
        _err = 0
        try:
            with redirect_argv(casa['flags'].execute):
                execfile(__candidates[0])

        except NameError, err:
            _err = 1
            if str(err) == "name 'T' is not defined" or \
               str(err) == "name 'F' is not defined" or \
               str(err) == "name 'true' is not defined" or \
               str(err) == "name 'false' is not defined" :
                print "------------------------------------------------------------------------------"
                print "Warning: CASA no longer defines T/true and F/false as synonyms for True/False"
                print "------------------------------------------------------------------------------"
                traceback.print_exc()
            else:
                traceback.print_exc()

        except Exception, err:
            _err = 1
            traceback.print_exc()

        sys.stdout.flush()
        sys.stderr.flush()
        immediate_exit_with_handlers(_err)

    else:
        # python command provided on the command line...
        _err = 0
        try:
            exec(casa['flags'].execute[0])
        except Exception, err:
            _err = 1
            traceback.print_exc()

        immediate_exit_with_handlers(_err)

else:

    ###
    ### Revisit if we ever get rid of plotcal al a matplotlib...
    ###
    ### ----------------------------------------------------------------------
    ### without this we get errors when our customized TkInter backend
    ### is being used in an open matplotlib window like:
    ### ----------------------------------------------------------------------
    ### Error in sys.exitfunc:
    ### Traceback (most recent call last):
    ###  File "/opt/casa/02/lib/python2.7/atexit.py", line 24, in _run_exitfuncs
    ###    func(*targs, **kargs)
    ###  File "/opt/casa/02/lib/python2.7/site-packages/matplotlib/_pylab_helpers.py", line 82, in destroy_all
    ###    manager.destroy()
    ###  File "/opt/casa/02/lib/python2.7/site-packages/matplotlib/backends/backend_tkagg.py", line 453, in destroy
    ###    self.window.destroy()
    ###  File "/opt/casa/02/lib/python2.7/lib-tk/Tkinter.py", line 1860, in destroy
    ###    self.tk.call('destroy', self._w)
    ###_tkinter.TclError: can't invoke "destroy" command:  application has been destroyed
    class ___protect_exit(object):
        "direct shutdown errors to /dev/null"
        def __init__( self, handler ):
            self.__handler = handler

        def __call__( self ):
            try:
                self.handler( )
            except:
                pass

    sys.exitfunc = ___protect_exit(sys.exitfunc)

    from casa_builtin import enable_builtin_protection, register_builtin

    ###
    ### backward compatibility at the command line...
    ### removed because we now create casa specific error messages about T/F removal...
    ###
    #T = True
    #F = False
    #true = True
    #false = False
    #register_builtin("T")
    #register_builtin("F")
    #register_builtin("true")
    #register_builtin("false")

    register_builtin("casa")
    register_builtin("cu")
    register_builtin(["viewer", "imview", "msview"])

    enable_builtin_protection()
    _blue = '\033[94m'
    _end = '\033[0m'
    print "Enter " + _blue + "doc('start')" + _end + " for help getting started with CASA..."
    #print "CASA Version " + casa['build']['version'] + "\n  Compiled on: " + casa['build']['time']
