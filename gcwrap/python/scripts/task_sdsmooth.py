import numpy
import os
from taskinit import *
import sdutil
ms,sdms,tb = gentools(['ms','sdms','tb'])

def sdsmooth(infile=None, datacolumn=None, antenna=None,
              field=None, spw=None, timerange=None, scan=None,
              pol=None, intent=None, kernel=None, kwidth=None,
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
                           polarization=pol, intent=intent)
        sdms.smooth(type=kernel, width=kwidth, datacolumn=datacolumn, outfile=outfile)
    except Exception as instance:
        raise Exception(instance)
    finally:
        sdms.close()

def raise_exception(errmsg):
    casalog.post(errmsg, priority='SEVERE')
    raise Exception(errmsg)
