import os
from taskinit import casalog, cbtool, mstool
from mstools import write_history
from parallel.parallel_task_helper import ParallelTaskHelper
from parallel.parallel_data_helper import ParallelDataHelper

_ms = mstool( )
_cb = cbtool( )

def delmod(vis=None,otf=None,field=None,scr=None):

    casalog.origin('delmod')

    # Do the trivial parallelization
    if ParallelDataHelper.isMMSAndNotServer(vis):
        helper = ParallelTaskHelper('delmod', locals())
        helper.go()
        return


    #Python script
    try:

        # only if vis exists...
        if ((type(vis)==str) & (os.path.exists(vis))):
            # ... and we are asked to do something...
            # open without adding anything!
            _cb.open(vis,addcorr=False,addmodel=False)
            _cb.delmod(otf=otf,field=field,scr=scr)
            _cb.close()
        else:
            raise Exception('Visibility data set not found - please verify the name')

        # Write history to MS
        try:
            param_names = delmod.__code__.co_varnames[:delmod.__code__.co_argcount]
            param_vals = [eval(p) for p in param_names]
            write_history(mstool(), vis, 'delmod', param_names,
                          param_vals, casalog)
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                         'WARN')

    except Exception as instance:
        print('*** Error ***',instance)
