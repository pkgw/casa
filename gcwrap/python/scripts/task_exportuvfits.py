import os
from taskinit import *

def exportuvfits(
    vis, fitsfile, datacolumn, field, spw, antenna, time,
    writesyscal, multisource, combinespw, 
    writestation, padwithflags, overwrite
):
    casalog.origin('exportuvfits')
    try:
        myms = mstool()
        if ((type(vis)==str) & (os.path.exists(vis))):
            myms.open( vis, lock=True )
        else:
            raise Exception('Visibility data set not found - please verify the name')
        writesyscal=False #until ms syscal table defined
        res = myms.tofits(
            fitsfile=fitsfile,
            column=datacolumn,
            field=field, spw=spw,
            baseline=antenna, time=time,
            writesyscal=writesyscal,
            multisource=multisource,
            combinespw=combinespw,
            writestation=writestation,
            padwithflags=padwithflags,
            overwrite=overwrite
        )
        if res:
            return True
        else:
            raise Exception("exportuvfits failed")
    except Exception as instance:
        casalog.post( '*** Error ***'+str(instance), 'SEVERE' )
        raise
    finally:
        if myms:
            myms.done()

