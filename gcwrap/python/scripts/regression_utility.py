import string
import inspect
import sys
import os
import shutil
from casa_stack_manip import stack_frame_find

myf=stack_frame_find( )
casalog=myf['casalog']

#
# Utils
#
def note(message, priority="INFO", origin="regression_utility", ntime=None, postcli='F'):
    #if not ntime:  #if (ntime==None):
    #    ntime=time.asctime()
    #print ntime, priority, origin, message
    if postcli: print(message)
    casalog.postLocally(message, priority, origin)
###
def info(message):
    #note(message,origin='regression_utility')
    #print message
    casalog.postLocally(message, priority="INFO", origin="regression_utility")

def fail(message=""):
    casalog.postLocally(message, priority="SEVERE", origin='regression_utility')
    #print message
    raise RuntimeError(message)

###
def stop(message=""):
    note(message ,priority='SEVERE', origin='regression_utility')
    raise RuntimeError(message)

###
def cleanup(dir):
    if (os.path.isdir(dir)):
        info("Cleaning up directory "+dir)
        def errFunc(raiser, problemPath, excInfo):
            #print raiser.__name__,'failed on',problemPath
            note(raiser.__name__+'failed on'+problemPath,"SEVERE")
            raise RuntimeError("Cleanup of " + dir + " fails!")
        shutil.rmtree(dir,0,errFunc)
    return True

###
def maketestdir(testdir):
    if not cleanup(testdir):
        note("Cleanup of "+testdir+" failed","SEVERE")
        return False
    try:
        os.mkdir(testdir)
    except IOError as e:
        note(e, "SEVERE")
        raise RuntimeError("mkdir " + testdir + " fails!")
