from taskinit import *
from ialib import write_image_history

def specsmooth(
    imagename, outfile, box, chans, stokes, region,  mask,
    overwrite, stretch, axis, function, width, dmethod
):
    casalog.origin('specsmooth')
    myia = iatool()
    myia.dohistory(False)
    outia = None
    try:
        if (not myia.open(imagename)):
            raise Exception("Cannot create image analysis tool using " + imagename)
        if (len(outfile) == 0):
            raise Exception("outfile must be specified.")
        function = function.lower()
        drop = len(dmethod) > 0
        if (function.startswith("b")):
            outia = myia.boxcar(
                outfile=outfile, region=region, mask=mask, overwrite=overwrite,
                stretch=stretch, axis=axis, width=width, drop=drop, dmethod=dmethod
            )
        elif (function.startswith("h")):
            outia = myia.hanning(
                outfile=outfile, region=region, mask=mask, overwrite=overwrite,
                stretch=stretch, axis=axis, drop=drop, dmethod=dmethod
            )
        else:
            raise Exception("Unsupported convolution function " + function)
        try:
            param_names = specsmooth.__code__.co_varnames[:specsmooth.__code__.co_argcount]
            param_vals = [eval(p) for p in param_names]   
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        return True
    except Exception as instance:
        casalog.post( str( '*** Error ***') + str(instance), 'SEVERE')
        raise
    finally:
        if myia:
            myia.done()
        if outia:
            outia.done()
        
