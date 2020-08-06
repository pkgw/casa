from taskinit import mstool, tbtool, casalog
from mstools import write_history
import flaghelper

def statwt(
    vis, selectdata, field, spw, intent, array, observation, scan,
    combine, timebin, slidetimebin, chanbin, minsamp, statalg,
    fence, center, lside, zscore, maxiter, fitspw, excludechans,
    wtrange, flagbackup, preview, datacolumn
):
    casalog.origin('statwt')
    if not selectdata:
        # CAS-10761, requirement provided by Urvashi
        if field or spw or intent or array or observation:
            casalog.post(
                "selectdata=False, any explicitly set data "
                + "selection parameters will be ignored",
                "WARN"
            )
        field = ""
        spw = ""
        intent = ""
        array = ""
        observation = ""
        scan = ""
    try:
        if (flagbackup):
            if (preview):
                casalog.post(
                    "Running in preview mode. No flags will be "
                    + "modified, so existing flags will not be backed up."
                )
            else:
                casalog.post('Backup original flags before applying new flags')
                flaghelper.backupFlags(aflocal=None, msfile=vis, prename='statwt')
        myms = mstool()
        myms.open(vis, nomodify=preview)
        sel = {}
        sel['spw'] = spw
        #sel['time'] = timerange
        sel['field'] = field
        #sel['baseline'] = antenna
        sel['scan'] = scan
        sel['scanintent'] = intent
        #sel['polarization'] = correlation
        #sel['uvdist'] = uvrange
        sel['observation'] = str(observation)
        sel['array'] = array
        #sel['feed'] = feed
        # Select the data. Only-parse is set to false.
        myms.msselect(sel, False)
        
        rval = None
        rval = myms.statwt(
            combine=combine, timebin=timebin,
            slidetimebin=slidetimebin, chanbin=chanbin,
            minsamp=minsamp, statalg=statalg, fence=fence,
            center=center, lside=lside, zscore=zscore,
            maxiter=maxiter, fitspw=fitspw, excludechans=excludechans,
            wtrange=wtrange, preview=preview, datacolumn=datacolumn
        )
        
        # Write to HISTORY of MS
        if rval != None and preview == False:
            # Write history to MS
            try:
                param_names = statwt.__code__.co_varnames[:statwt.__code__.co_argcount]
                param_vals = [eval(p) for p in param_names]
                write_history(mstool(), vis, 'statwt', param_names,
                              param_vals, casalog)
            except Exception as instance:
                casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                             'WARN')            
            
        return rval
    except Exception as instance:
        casalog.post( '*** Error ***'+str(instance), 'SEVERE' )
        raise
    finally:
        myms.done()
