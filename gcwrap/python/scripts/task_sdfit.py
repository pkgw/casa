import numpy
import os
import time, datetime
from taskinit import gentools, casalog
import sdutil
ms, sdms, tb = gentools(['ms','sdms','tb'])


from numpy import ma, array, logical_not, logical_and


@sdutil.sdtask_decorator
def sdfit(infile=None, datacolumn=None, antenna=None, field=None, spw=None,
           timerange=None, scan=None, pol=None, intent=None,
           timebin=None, timespan=None,
           polaverage=None, 
           fitfunc=None, fitmode=None, nfit=None, thresh=None, avg_limit=None,
           minwidth=None, edge=None, outfile=None, overwrite=None):
    casalog.origin('sdfit')

    try:
        if os.path.exists(outfile):
            if overwrite:
                os.system('rm -rf %s' % outfile)
            else:
                raise ValueError(outfile + ' exists.')
        if (fitmode not in  ['list', 'auto']):
            raise ValueError("fitmode='%s' is not supported yet" % fitmode)
        if (spw == ''): spw = '*'

        selection = ms.msseltoindex(vis=infile, spw=spw, field=field, 
                                    baseline=antenna, time=timerange, 
                                    scan=scan)

        sdms.open(infile)
        sdms.set_selection(spw=sdutil.get_spwids(selection), field=field, 
                           antenna=antenna, timerange=timerange, 
                           scan=scan, polarization=pol, intent=intent)

        tempfile = 'temp_sdfit_'+str(datetime.datetime.fromtimestamp(time.time())).replace('-','').replace(' ','').replace(':','')
        if os.path.exists(tempfile):
            tempfile += str(datetime.datetime.fromtimestamp(time.time())).replace('-','').replace(' ','').replace(':','')
            if os.path.exists(tempfile):
                raise Exception('temporary file ' + tempfile + ' exists...')
        tempoutfile = tempfile + '_temp_output_ms'
        if os.path.exists(tempoutfile):
            tempoutfile += str(datetime.datetime.fromtimestamp(time.time())).replace('-','').replace(' ','').replace(':','')
            if os.path.exists(tempoutfile):
                raise Exception('temporary ms file ' + tempoutfile + ' exists...')
        
        if fitmode=='auto': nfit = [-1]
        num_fit_str = str(',').join(map(str, nfit))

        sdms.fit_line(datacolumn=datacolumn, spw=spw, pol=pol, 
                      timebin=timebin, timespan=timespan,
                      polaverage=polaverage, 
                      fitfunc=fitfunc, nfit=num_fit_str,
                      linefinding=(fitmode=='auto'), threshold=thresh,
                      avg_limit=avg_limit, minwidth=minwidth, edge=edge,
                      tempfile=tempfile, tempoutfile=tempoutfile)

        if os.path.exists(tempfile):
            return get_results(tempfile, fitfunc, nfit, outfile)
        else:
            raise Exception('temporary file was unexpectedly not created.')

    except Exception as instance:
        raise Exception(instance)
    finally:
        if 'tempfile' in locals() and os.path.exists(tempfile):
            os.system('rm -f %s' % tempfile)
        if 'tempoutfile' in locals() and os.path.exists(tempoutfile):
            os.system('rm -rf %s' % tempoutfile)

def get_results(tempfile, fitfunc, nfit, outfile):
    res = {'cent':[], 'peak':[], 'fwhm':[], 'nfit':[]}
    if (fitfunc == 'gaussian'):
        func = 'gauss'
    elif (fitfunc == 'lorentzian'):
        func = 'lorentz'
    ncomp = sum(nfit)
    iline = 0
    with open(tempfile, 'r') as f:
        outfile_exists = (outfile != '')
        if outfile_exists:
            fout = open(outfile, 'a')
            s = '#SCAN\tTIME\t\tANT\tBEAM\tSPW\tPOL\tFunction\tP0\t\tP1\t\tP2\n'
            fout.write(s)
        
        for line in f:
            component = line.strip().split(':')   # split into each component
            if (ncomp > 0): 
                assert(len(component) == ncomp)
            res['cent'].append([])
            res['peak'].append([])
            res['fwhm'].append([])
            res['nfit'].append(ncomp if ncomp>=0 else len(component))
            for icomp in range(res['nfit'][-1]):
                fit_result = component[icomp].split(',')   # split into each parameter
                num_ids = 6 # scan, time, ant, beam, spw, pol
                assert(len(fit_result) == 2*(len(list(res.keys()))-1)+num_ids)
                res['cent'][iline].append([float(fit_result[num_ids+0]), float(fit_result[num_ids+1])])
                res['peak'][iline].append([float(fit_result[num_ids+2]), float(fit_result[num_ids+3])])
                res['fwhm'][iline].append([float(fit_result[num_ids+4]), float(fit_result[num_ids+5])])

                if outfile_exists:
                    s =  fit_result[0] + '\t'       #scanID
                    s += fit_result[1] + '\t'       #time
                    s += fit_result[2] + '\t'       #antennaID
                    s += fit_result[3] + '\t'       #beamID
                    s += fit_result[4] + '\t'       #spwID
                    s += fit_result[5] + '\t'       #polID
                    s += func + str(icomp) + '\t\t' #function
                    s += fit_result[8] + '\t'       #P0 (peak)
                    s += fit_result[6] + '\t'       #P1 (cent)
                    s += fit_result[10] + '\t'      #P2 (fwhm)
                    s += '\n'
                    fout.write(s)
            iline += 1
        
        if outfile_exists: fout.close()
    
    return res
