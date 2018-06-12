# sd task for imaging
from casac import casac
from taskinit import *

# def sdsidebandsplit(infiles, outfile, overwrite, field, spw, antenna, scan, intent,
#                     imageshift, getbothside, lo1, loframe, reftime, refdir, threshold,
#                     mode, nchan, start, width, veltype, outframe,
#                     gridfunction, convsupport,truncate, gwidth, jwidth,
#                     imsize, cell, phasecenter, ephemsrcname, pointingcolumn,
#                     restfreq, stokes, minweight, clipminmax):
def sdsidebandsplit(imagename, outfile, overwrite, signalshift, imageshift,
                    getbothside, refchan, refval, useother, threshold):
    casalog.origin('sdsidebandsplit')
    
    separator = casac.sidebandseparator()
    try:
        separator.open(imagename)
        separator.setshift(signalshift, True)
        if len(imageshift) > 0:
            separator.setshift(imageshift, False)
        separator.setlimit(threshold)
        separator.setboth(getbothside)
        if getbothside:
            if refval == '':
                qrefval = -1.0
            else:
                myqa = qatool()
                qrefval = myqa.quantity(refval)
            separator.set_imageband_frequency(refchan, qrefval)
        separator.setsolveother(useother)
        separator.separate(outfile, overwrite)
    finally:
        separator.close()
    
    
    """
    min_casarev = 23268
    try:
        # checking for CASA revision
        a=inspect.stack()
        for k in range(len(a)):
            if (string.find(a[k][1], 'ipython console') > 0):
                stacklevel=k
                break
        myf=sys._getframe(stacklevel).f_globals
        pkgrev = int(myf['casa']['build']['number'])
        if pkgrev < min_casarev:
            casalog.post("CASA version >= %d required for this interface"\
                         % min_casarev, "ERROR")


        sbsep = sd.sbseparator(infiles)
        
        # Set frequency information to select data
        if not (qa.compare(freqtol, "Hz") or qa.quantity(freqtol)['unit'] == ''):
            del sbsep
            casalog.post("Invalid unit for frequency tolerance", "ERROR")
        sbsep.set_frequency(ifno, freqtol, frame=frame)

        if len(imageshift) > 0:
            # The format is checked cpp class
            sbsep.set_shift(imageshift=imageshift)

        # Set direction tolerance to grid data
        sbsep.set_dirtol(dirtol)

        # Set rejection limit of solution
        sbsep.set_limit(threshold)
        # Set which solution to get
        sbsep.set_both(getbothside)
        # Set information to get frequency information of image sideband
        if getbothside and len(lo1) > 0:
            if isinstance(reftime, str):
                myme = metool()
                mereftime = myme.epoch('TAI', reftime)
                reftime = mereftime['m0']['value']
            sbsep.set_lo1(lo1, loframe, reftime, refdir)
        
        sbsep.set_solve_other(False)
        
        # Invoke separation supression
        sbsep.separate(outfile, overwrite = overwrite)

        del sbsep

    except Exception, instance:
        #print '***Error***',instance
        casalog.post( str(instance), priority = 'ERROR' )
        raise Exception, instance
        return
    """
