from taskinit import *

_rg = rgtool( )

def imsubimage(
    imagename, outfile, box, region, chans, stokes, mask, dropdeg,
    overwrite, verbose, stretch, keepaxes
):
    casalog.origin('imsubimage')
    myia = iatool()
    outia = None
    try:
        if (not myia.open(imagename)):
            raise Exception("Cannot create image analysis tool using " + imagename)
        if (len(outfile) == 0):
            raise Exception("outfile must be specified.")
        if (type(region) != type({})):
            region = _rg.frombcs(
                csys=myia.coordsys().torecord(), shape=myia.shape(), box=box,
                chans=chans, stokes=stokes, stokescontrol="a", region=region
            )
        outia = myia.subimage(
            outfile=outfile, region=region, mask=mask, dropdeg=dropdeg,
            overwrite=overwrite, list=verbose, stretch=stretch, keepaxes=keepaxes
        )
        return True
    except Exception as instance:
        casalog.post( str( '*** Error ***') + str(instance), 'SEVERE')
        raise
    finally:
        if myia:
            myia.done()
        if outia:
            outia.done()
        
