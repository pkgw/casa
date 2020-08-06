import numpy
import os
from taskinit import casalog, gentools
from mstools import write_history
import sdutil
ms,sdms,tb = gentools(['ms','sdms','tb'])

def sdsmooth(infile=None, datacolumn=None, antenna=None, 
              field=None, spw=None, timerange=None, scan=None, 
              pol=None, intent=None, reindex=None,
              kernel=None, kwidth=None,
              outfile=None, overwrite=None):

    casalog.origin('sdsmooth')

    try:
        if len(outfile) == 0:
            errmsg = 'outfile is empty.'
            raise_exception(errmsg)
        
        if (os.path.exists(outfile)) and (not overwrite):
            errmsg = outfile+' exists.'
            raise_exception(errmsg)

        sdms.open(infile)
        sdms.set_selection(spw=spw, field=field, 
                           antenna=antenna,
                           timerange=timerange, scan=scan,
                           polarization=pol, intent=intent,
                           reindex=reindex)
        sdms.smooth(type=kernel, width=kwidth, datacolumn=datacolumn, outfile=outfile)
        
        # Write to HISTORY of outfile MS
        param_names = sdsmooth.__code__.co_varnames[:sdsmooth.__code__.co_argcount]
        param_vals = [eval(p) for p in param_names]
        write_history(ms, outfile, 'sdsmooth', param_names,
                      param_vals, casalog)

    except Exception as instance:
        raise Exception(instance)
    finally:
        sdms.close()

def raise_exception(errmsg):
    casalog.post(errmsg, priority='SEVERE')
    raise Exception(errmsg)
