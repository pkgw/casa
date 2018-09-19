from taskinit import *

from ialib import write_image_history

def imdev(
    imagename, outfile, region, box, chans,
    stokes, mask, overwrite, stretch,
    grid, anchor, xlength, ylength, interp, stattype, statalg,
    zscore, maxiter
):
    _myia = iatool()
    _myrg = rgtool()
    _mycs = cstool()
    try:
        casalog.origin('imdev')
        _myia.open(imagename)
        _mycs = _myia.coordsys()
        csrec = _mycs.torecord()
        shape =  _myia.shape()
        reg = _myrg.frombcs(
            csrec, shape,
            box, chans, stokes, "a", region
        )
        outia = _myia.deviation(
            outfile=outfile, region=reg, mask=mask,
            overwrite=overwrite, stretch=stretch, grid=grid,
            anchor=anchor, xlength=xlength, ylength=ylength,
            interp=interp, stattype=stattype, statalg=statalg,
            zscore=zscore, maxiter=maxiter
        )
        try:
            param_names = imdev.__code__.co_varnames[:imdev.__code__.co_argcount]
            param_vals = [eval(p) for p in param_names]   
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        outia.done() 
        return True
    except Exception as instance:
        casalog.post( '*** Error ***'+str(instance), 'SEVERE' )
        raise
    finally:
        _myia.done()
        _myrg.done()
        _mycs.done()
        if outia:
            outia.done()
