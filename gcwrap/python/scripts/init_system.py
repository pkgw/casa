import os
import sys
import time
import argparse
import multiprocessing
from IPython.terminal.prompts import Prompts, Token

try:
    from casac import casac
except ImportError, e:
    print "failed to load casa:\n", e
    os._exit(1)

try:
    import matplotlib
except ImportError, e:
    print "failed to load matplotlib:\n", e
    print "sys.path =", "\n\t".join(sys.path)

from asap_init import *
from casa_system import casa

if not os.environ.has_key('OMP_NUM_THREADS'):
    # if OMP_NUM_THREADS is not set, set it to max(1,N_CPU-2)
    os.environ['OMP_NUM_THREADS'] = str(max(1,multiprocessing.cpu_count()-2))

class _Prompt(Prompts):
     def in_prompt_tokens(self, cli=None):
         return [(Token.Prompt, 'CASA <'),
                 (Token.PromptNum, str(self.shell.execution_count)),
                 (Token.Prompt, '>: ')]

_ip = get_ipython()
_ip.prompts = _Prompt(_ip)

###
### provide extra context for T/F errors...
###
def true_false_handler(self, etype, value, tb, tb_offset=None):
    if type(etype) is type(NameError):
        if str(value) == "name 'T' is not defined" or \
           str(value) == "name 'F' is not defined" or \
           str(value) == "name 'true' is not defined" or \
           str(value) == "name 'false' is not defined" :
            print "------------------------------------------------------------------------------"
            print "Warning: CASA no longer defines T/true and F/false as synonyms for True/False"
            print "------------------------------------------------------------------------------"
    return self.showtraceback()

_ip.set_custom_exc((BaseException,), true_false_handler)


##
## toplevel frame marker
##
_casa_top_frame_ = True

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
## set up casa root
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
if os.environ.has_key('CASAPATH') :
    __casapath__ = os.environ['CASAPATH'].split(' ')[0]
    __casaarch__ = os.environ['CASAPATH'].split(' ')[1]
    if not os.path.exists(__casapath__ + "/data") :
        print "DEBUG: CASAPATH = %s" % (__casapath__)
        print "Unable to find the data repository directory in your CASAPATH. Please fix."
        os._exit(1)
    else :
        casa['dirs']['root'] = __casapath__
        casa['dirs']['data'] = __casapath__ + "/data"

        if os.path.exists(__casapath__ + "/lib/python2.7/start_casa.py"):
            casa['dirs']['python'] = __casapath__ + "/lib/python2.7"
        elif os.path.exists(__casapath__ + "/Resources/python/start_casa.py"):
            casa['dirs']['python'] = __casapath__ + "/Resources/python"
        elif os.path.exists(__casapath__ + "/" + __casaarch__ + "/lib/python2.7/start_casa.py"):
            casa['dirs']['python'] = __casapath__ + "/" + __casaarch__ + "/lib/python2.7"

        if casa['dirs']['python'] is not None:
            casa['dirs']['recipes'] = casa['dirs']['python'] + "/recipes"

        if os.path.exists(__casapath__ + "/" + __casaarch__ + "/xml"):
            casa['dirs']['xml'] = __casapath__ + "/" + __casaarch__ + "/xml"
        elif os.path.exists(__casapath__ + "/xml"):
            casa['dirs']['xml'] = __casapath__ + "/xml"
        else:
            raise RuntimeError, "Unable to find the XML constraints directory in your CASAPATH"

        casa['dirs']['doc'] = None
        if os.path.exists(__casapath__ + "/share/doc"):
            casa['dirs']['doc'] = __casapath__ + "/share/doc"
        elif os.path.exists(__casapath__ + "/doc"):
            casa['dirs']['doc'] = __casapath__ + "/doc"
        elif os.path.exists(__casapath__ + "/Resources/doc"):
            casa['dirs']['doc'] = __casapath__ + "/Resources/doc"

else :
    __casapath__ = casac.__file__
    while __casapath__ and __casapath__ != "/" :
        if os.path.exists( __casapath__ + "/data") :
            break
        __casapath__ = os.path.dirname(__casapath__)
    if not os.path.exists(__casapath__ + "/data") :
        raise RuntimeError, "casa path could not be determined"
    else :
        casa['dirs']['root'] = __casapath__
        casa['dirs']['data'] = __casapath__ + "/data"
        if os.path.exists(__casapath__ + "/" + __casaarch__ + "python/2.7/assignmentFilter.py"):
            casa['dirs']['python'] = __casapath__ + "/" + __casaarch__ + "/python/2.7"
        elif os.path.exists(__casapath__ + "/lib/python2.7/assignmentFilter.py"):
            casa['dirs']['python'] = __casapath__ + "/lib/python2.7"
        elif os.path.exists(__casapath__ + "/Resources/python/assignmentFilter.py"):
            casa['dirs']['python'] = __casapath__ + "/Resources/python"

        if casa['dirs']['python'] is not None:
            casa['dirs']['recipes'] = casa['dirs']['python'] + "/recipes"

        if os.path.exists(__casapath__ + "/" + __casaarch__ + "/xml"):
            casa['dirs']['xml'] = __casapath__ + "/" + __casaarch__ + "/xml"
        elif os.path.exists(__casapath__ + "/xml"):
            casa['dirs']['xml'] = __casapath__ + "/xml"
        else:
            raise RuntimeError, "Unable to find the XML constraints directory in your CASAPATH"

        casa['dirs']['doc'] = None
        if os.path.exists(__casapath__ + "/share/doc"):
            casa['dirs']['doc'] = __casapath__ + "/share/doc"
        elif os.path.exists(__casapath__ + "/doc"):
            casa['dirs']['doc'] = __casapath__ + "/doc"
        elif os.path.exists(__casapath__ + "/Contents/Resources/doc"):
            casa['dirs']['doc'] = __casapath__ + "/Contents/Resources/doc"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
## try to set casapyinfo path...
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
if os.path.exists( __casapath__ + "/bin/casa-config") :
    casa['helpers']['info'] = __casapath__ + "/bin/casa-config"

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
##     first try to find executables using casapyinfo...
##            (since system area versions may be incompatible)...
##     next try likely system areas...
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
##
##   note:  hosts which have dbus-daemon-1 but not dbus-daemon seem to have a broken dbus-daemon-1...
##
for info in [ (['dbus-daemon'],'dbus'),
              (['CrashReportPoster'],'crashPoster'),
              (['ipcontroller','ipcontroller-2.6'], 'ipcontroller'),
              (['ipengine','ipengine-2.6'], 'ipengine') ]:
    exelist = info[0]
    entry = info[1]
    for exe in exelist:
        if casa['helpers']['info']:
            casa['helpers'][entry] = (lambda fd: fd.readline().strip('\n'))(os.popen(casa['helpers']['info'] + " --exec 'which " + exe + "'"))
        if casa['helpers'][entry] and os.path.exists(casa['helpers'][entry]):
            break
        else:
            casa['helpers'][entry] = None

        ### first look in known locations relative to top (of binary distros) or known casa developer areas
        for srchdir in [ __casapath__ + '/MacOS', __casapath__ + '/lib/casa/bin', '/usr/lib64/casa/01/bin', '/opt/casa/01/bin' ] :
            dd = srchdir + os.sep + exe
            if os.path.exists(dd) and os.access(dd,os.X_OK) :
                casa['helpers'][entry] = dd
                break
        if casa['helpers'][entry] is not None:
            break

    ## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    ##     next search through $PATH for executables
    ## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    if casa['helpers'][entry] is None:
        for exe in exelist:
            for srchdir in os.getenv('PATH').split(':') :
                dd = srchdir + os.sep + exe
                if os.path.exists(dd) and os.access(dd,os.X_OK) :
                    casa['helpers'][entry] = dd
                    break
            if casa['helpers'][entry] is not None:
                break

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
## try to set pipeline path...
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
if os.path.exists(casa['dirs']['root']+"/pipeline"):
    casa['dirs']['pipeline'] = casa['dirs']['root']+"/pipeline"

# initialize/finalize Sakura library
if hasattr(casac,'sakura'):
    #casalog.post('Managing Sakura lifecycle', priority='DEBUG')
    casac.sakura().initialize_sakura()
    import atexit
    atexit.register(lambda: __import__('casac').casac.sakura().cleanup_sakura())

class iArgumentParser(argparse.ArgumentParser):
    '''iPython thinks it knows that a user would never want to
    exit in any way other than typing "exit" at the command line'''
    def exit(self, status=0, message=None):
        if message:
            self._print_message(message, sys.stderr)
        os._exit(status)

argparser = iArgumentParser(prog="casa",description='Start CASA (Common Astronomy Software Applications)')

argparser.add_argument( '--rcdir',dest='rcdir',default=casa['dirs']['rc'],
                        help='location for startup files' )
argparser.add_argument( '--logfile',dest='logfile',default=casa['files']['logfile'],
                        help='path to log file' )
argparser.add_argument( "--maclogger",dest='maclogger',action='store_const',const='console',
                        default=__casapath__+'/'+__casaarch__+'/apps/casalogger.app/Contents/MacOS/casalogger',
                        help='logger to use on Apple systems' )
argparser.add_argument( "--log2term",dest='log2term',action='store_const',const=True,default=False,
                        help='direct output to terminal' )
argparser.add_argument( "--nologger",dest='nologger',action='store_const',const=True,default=False,
                        help='do not start CASA logger' )
argparser.add_argument( "--nologfile",dest='nologfile',action='store_const',const=True,default=False,
                        help='do not create a log file' )
argparser.add_argument( "--nogui",dest='nogui',action='store_const',const=True,default=False,
                        help='avoid starting GUI tools' )
argparser.add_argument( '--colors', dest='prompt', default='NoColor',
                        help='prompt color', choices=['NoColor', 'Linux', 'LightBG'] )
argparser.add_argument( "--trace",dest='trace',action='store_const',const=True,default=False,
                        help='list imported modules' )
argparser.add_argument( "--pipeline",dest='pipeline',action='store_const',const=True,default=False,
                        help='start CASA pipeline run' )
argparser.add_argument( "--agg",dest='agg',action='store_const',const=True,default=False,
                        help='startup without tkagg' )
argparser.add_argument( '--iplog',dest='ipython_log',default=False,
                          const=True,action='store_const',
                          help='create ipython log' )
argparser.add_argument( "-c",dest='execute',default=[],nargs=argparse.REMAINDER,
                        help='python eval string or python script to execute' )

casa['flags'], casa['args'] = argparser.parse_known_args( )
#### must keep args in sync with 'casa' state...
casa['files']['logfile'] = casa['flags'].logfile
casa['dirs']['rc'] = casa['flags'].rcdir

#### pipeline requires the Agg backend; any use of
#### matplotlib before 'init_pipeline.py' is loaded
#### would affect the ability to set the backend...
if casa['flags'].pipeline or casa['flags'].agg:
    matplotlib.use('Agg')

### provide details about what is being imported:
### before the module is imported:
###
###      importer => importee
###
### is printed. After the module has been imported:
###
###      ---> importee: <path to importee>
###
### is printed...
if casa['flags'].trace:
    import inspect
    import __builtin__
    _savimp = __builtin__.__import__

    def _newimp(name, *x):
        caller = inspect.currentframe( ).f_back
        print "%s => %s" % (caller.f_globals.get('__name__'), name)
        result = _savimp(name, *x)
        print "---> %s: %s" % (name, result.__file__ if hasattr(result,'__file__') else '?')
        return result

    __builtin__.__import__ = _newimp

print "CASA %s -- Common Astronomy Software Applications\n" % casa['build']['version']
