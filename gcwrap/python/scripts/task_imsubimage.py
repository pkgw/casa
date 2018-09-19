from taskinit import *
from ialib import write_image_history

_rg = rgtool()

def imsubimage(
    imagename, outfile, box, region, chans, stokes, mask, dropdeg,
    overwrite, verbose, stretch, keepaxes
):
    casalog.origin('imsubimage')
    myia = iatool()
    myia.dohistory(False)
    outia = None
    try:
        if (not myia.open(imagename)):
            raise Exception("Cannot create image analysis tool using " + imagename)
        if (len(outfile) == 0):
            raise Exception("outfile must be specified.")
        xregion = region
        if (type(region) != type({})):
            xregion = _rg.frombcs(
                csys=myia.coordsys().torecord(), shape=myia.shape(), box=box,
                chans=chans, stokes=stokes, stokescontrol="a", region=region
            )
        outia = myia.subimage(
            outfile=outfile, region=xregion, mask=mask, dropdeg=dropdeg,
            overwrite=overwrite, list=verbose, stretch=stretch, keepaxes=keepaxes
        )
        try:
            param_names = imsubimage.__code__.co_varnames[:imsubimage.__code__.co_argcount]
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
        
