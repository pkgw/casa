import os, re
import shutil
import string
import copy
import math
import time
from taskinit import mttool, mstool, tbtool, casalog, qa
from mstools import write_history
from parallel.parallel_data_helper import ParallelDataHelper
import flaghelper as fh
from update_spw import update_spwchan
from callibrary import callibrary

"""
The following code is based on the mstransform code, with 
task name and some task parameters modified. 
To minimise code modification, the parameters used by 
mstransform but not by sdpolaverage are kept and the 
default values for mstransform are given to them.
(CAS-12083, 2019/1/22 WK)
"""

def sdpolaverage(
             infile,
             datacolumn,
             antenna, 
             field,
             spw, 
             timerange, 
             scan,
             intent,
             polaverage,
             outfile):

    # followings are parameters of mstransform but not used by sdpolaverage.
    # just putting default values
    vis = infile             # needed for ParallelDataHelper
    outputvis = outfile      # needed for ParallelDataHelper
    createmms = False
    separationaxis = "auto"
    numsubms = "auto"
    tileshape = [0]
    correlation = ""
    array = ""
    uvrange = ""
    observation = ""
    feed = ""
    realmodelcol = False
    keepflags = True
    usewtspectrum = False
    combinespws = False
    chanaverage = False
    chanbin = 1
    hanning = False
    regridms = False
    mode = "channel"
    nchan = -1
    start = 0
    width = 1
    nspw = 1
    interpolation = "linear"
    phasecenter = ""
    restfreq = ""
    outframe = ""
    veltype = "radio"
    preaverage = False
    timeaverage = False
    timebin = "0s"
    timespan = ""
    maxuvwdistance = 0.0
    docallib = False
    callib = ""
    douvcontsub = False
    fitspw = ""
    fitorder = 0
    want_cont = False
    denoising_lib = True
    nthreads = 1
    niter = 1
    disableparallel = False
    ddistart = -1
    taql = ""
    monolithic_processing = False
    reindex = True
    
    casalog.origin('sdpolaverage')
    
    # Initialize the helper class
    pdh = ParallelDataHelper('sdpolaverage', locals())
    
    # When dealing with MMS, process in parallel or sequential
    # disableparallel is a hidden parameter. Only for debugging purposes!
    if disableparallel:
        pdh.bypassParallelProcessing(1)
    else:
        pdh.bypassParallelProcessing(0)
    
    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Process the input Multi-MS
    if ParallelDataHelper.isMMSAndNotServer(infile) == True and monolithic_processing == False:
        '''
        retval{'status': True,  'axis':''}         --> can run in parallel        
        retval{'status': False, 'axis':'value'}    --> treat MMS as monolithic MS, set new axis for output MMS
        retval{'status': False, 'axis':''}         --> treat MMS as monolithic MS, create an output MS
        '''
        
        retval = pdh.validateInputParams()
        # Cannot create an output MMS.
        if retval['status'] == False and retval['axis'] == '':
            casalog.post('Cannot process MMS with the requested transformations','WARN')
            casalog.post('Use task listpartition to see the contents of the MMS')
            casalog.post('Will create an output MS','WARN')
            createmms = False
            
        # MMS is processed as monolithic MS. 
        elif retval['status'] == False and retval['axis'] != '':
            createmms = True
            pdh.override__args('createmms', True)
            pdh.override__args('monolithic_processing', True)
            separationaxis = retval['axis']
            pdh.override__args('separationaxis', retval['axis'])
            casalog.post("Will process the input MMS as a monolithic MS",'WARN')
            casalog.post("Will create an output MMS with separation axis \'%s\'"%retval['axis'],'WARN')
            
        # MMS is processed in parallel
        else:
            createmms = False
            try:
                pdh.override__args('createmms', False)
                pdh.setupCluster('sdpolaverage')
                pdh.go()
            except Exception as instance:
                casalog.post('%s'%instance,'ERROR')
                return False
            
            return True

    # Create an output Multi-MS
    if createmms == True:
        
        # Check the heuristics of separationaxis and the requested transformations
        pval = pdh.validateOutputParams()
        if pval == 0:
            raise Exception('Cannot create MMS using separationaxis=%s with some of the requested transformations.'\
                            %separationaxis)
                             
        try:
            pdh.setupCluster('sdpolaverage')
            pdh.go()
            monolithic_processing = False
        except Exception as instance:
            casalog.post('%s'%instance,'ERROR')
            return False
        
        return True
                    
        
    # Create a local copy of the MSTransform tool
    mtlocal = mttool()
    mslocal = mstool()
        
    try:
        # Gather all the parameters in a dictionary.
        config = {}
        
        if keepflags:
            taqlstr = ''
        else:
            taqlstr = "NOT (FLAG_ROW OR ALL(FLAG))"
        
        # MMS taql selection
        if taql != '' and taql != None:
            if not keepflags:
                taqlstr = taqlstr + " AND "+taql
            else:
                taqlstr = taql
        
        config = pdh.setupParameters(inputms=infile, outputms=outfile, field=field, 
                    spw=spw, array=array, scan=scan, antenna=antenna, correlation=correlation,
                    uvrange=uvrange,timerange=timerange, intent=intent, observation=str(observation),
                    feed=feed, taql=taqlstr)
        
        # ddistart will be used in the tool when re-indexing the spw table
        config['ddistart'] = ddistart
        
        # re-index parameter is used by the pipeline to not re-index any sub-table and the associated IDs
        config['reindex'] = reindex        
        
        config['datacolumn'] = datacolumn
        dc = datacolumn.upper()            
        # Make real a virtual MODEL column in the output MS
        if "MODEL" in dc or dc == 'ALL':
            config['realmodelcol'] = realmodelcol

        config['usewtspectrum'] = usewtspectrum
        
        # Add the tile shape parameter
        if tileshape.__len__() == 1:
            # The only allowed values are 0 or 1
            if tileshape[0] != 0 and tileshape[0] != 1:
                raise ValueError('When tileshape has one element, it should be either 0 or 1.')
                
        elif tileshape.__len__() != 3:
            # The 3 elements are: correlations, channels, rows
            raise ValueError('Parameter tileshape must have 1 or 3 elements.')
            
        config['tileshape'] = tileshape                

        if combinespws:
            casalog.post('Combine spws %s into new output spw'%spw)
            config['combinespws'] = True
            
        # Only parse chanaverage if chanbin is valid
        if chanaverage and isinstance(chanbin, int) and chanbin <= 1:
            raise Exception('Parameter chanbin must be > 1 to do channel averaging')
            
        # Validate the case of int or list chanbin
        if chanaverage and pdh.validateChanBin():
            casalog.post('Parse channel averaging parameters')
            config['chanaverage'] = True
            
            # convert numpy types, until CAS-6493 is not fixed
            chanbin = fh.evaluateNumpyType(chanbin)
            config['chanbin'] = chanbin
            
        if hanning:
            casalog.post('Apply Hanning smoothing')
            config['hanning'] = True
            
        if regridms:
            casalog.post('Parse regridding parameters')            
            config['regridms'] = True
            # Reset the defaults depending on the mode
            # Only add non-empty string parameters to config dictionary
            start, width = pdh.defaultRegridParams()
            config['mode'] = mode
            config['nchan'] = nchan
            if start != '':
                config['start'] = start
            if width != '':
                config['width'] = width
            if nspw > 1:
                casalog.post('Separate MS into %s spws'%nspw)
            config['nspw'] = nspw
            config['interpolation'] = interpolation
            if restfreq != '':
                config['restfreq'] = restfreq
            if outframe != '':
                config['outframe'] = outframe
            if phasecenter != '':
                config['phasecenter'] = phasecenter
            config['veltype'] = veltype
            config['preaverage'] = preaverage
            
        # Only parse timeaverage parameters when timebin > 0s
        if timeaverage:
            tb = qa.convert(qa.quantity(timebin), 's')['value']
            if not tb > 0:
                raise Exception("Parameter timebin must be > '0s' to do time averaging")
                       
        if timeaverage:
            casalog.post('Parse time averaging parameters')
            config['timeaverage'] = True
            config['timebin'] = timebin
            config['timespan'] = timespan
            config['maxuvwdistance'] = maxuvwdistance
            
        polaverage_ = polaverage.strip()
        if polaverage_ != '':
            config['polaverage'] = True
            config['polaveragemode'] = polaverage_

        if docallib:
            casalog.post('Parse docallib parameters')
            mycallib = callibrary()
            mycallib.read(callib)
            config['calibration'] = True
            config['callib'] = mycallib.cld
            
        if douvcontsub:
            casalog.post('Parse uvcontsub parameters')
            config['uvcontsub'] = True
            uvcontsub_config = {}
            uvcontsub_config['fitspw'] = fitspw
            uvcontsub_config['fitorder'] = fitorder
            uvcontsub_config['want_cont'] = want_cont
            uvcontsub_config['denoising_lib'] = denoising_lib   
            uvcontsub_config['nthreads'] = nthreads            
            uvcontsub_config['niter'] = niter                 
            config['uvcontsublib'] = dict(uvcontsub_config)
        
        # Configure the tool and all the parameters
        casalog.post('%s'%config, 'DEBUG')
        mtlocal.config(config)
        
        # Open the MS, select the data and configure the output
        mtlocal.open()
        
        # Run the tool
        casalog.post('Apply the transformations')
        mtlocal.run()        

        mtlocal.done()
                    
    except Exception as instance:
        mtlocal.done()
        casalog.post('%s'%instance,'ERROR')
        return False

    # Update the FLAG_CMD sub-table to reflect any spw/channels selection
    # If the spw selection is by name or FLAG_CMD contains spw with names, skip the updating    
    
    if ((spw != '') and (spw != '*')) or chanaverage == True:
        isopen = False

        try:
            mytb = tbtool()
            mytb.open(outfile + '/FLAG_CMD', nomodify=False)
            isopen = True
            nflgcmds = mytb.nrows()
            
            if nflgcmds > 0:
                updateFlagCmd = False

                # If spw selection is by name in FLAG_CMD, do not update, CAS-7751
                mycmd = mytb.getcell('COMMAND', 0)
                cmdlist = mycmd.split()
                for cmd in cmdlist:
                    # Match only spw indices, not names
                    if cmd.__contains__('spw'):
                        cmd = cmd.strip("spw=")
                        spwstr = re.search('^[^a-zA-Z]+$', cmd)
                        if spwstr != None and spwstr.string.__len__() > 0:
                            updateFlagCmd = True
                            break                
                

                if updateFlagCmd:
                    mademod = False
                    cmds = mytb.getcol('COMMAND')
                    widths = {}
                    #print "width =", width
                    if hasattr(chanbin, 'has_key'):
                        widths = chanbin
                    else:
                        if hasattr(chanbin, '__iter__') and len(chanbin) > 1:
                            for i in range(len(chanbin)):
                                widths[i] = chanbin[i]
                        elif chanbin != 1:
    #                        print 'using ms.msseltoindex + a scalar width'
                            numspw = len(mslocal.msseltoindex(vis=infile,
                                                         spw='*')['spw'])
                            if hasattr(chanbin, '__iter__'):
                                w = chanbin[0]
                            else:
                                w = chanbin
                            for i in range(numspw):
                                widths[i] = w
    #                print 'widths =', widths 
                    for rownum in range(nflgcmds):
                        # Matches a bare number or a string quoted any way.
                        spwmatch = re.search(r'spw\s*=\s*(\S+)', cmds[rownum])
                        if spwmatch:
                            sch1 = spwmatch.groups()[0]
                            sch1 = re.sub(r"[\'\"]", '', sch1)  # Dequote
                            # Provide a default in case the split selection excludes
                            # cmds[rownum].  update_spwchan() will throw an exception
                            # in that case.
                            cmd = ''
                            try:
                                #print 'sch1 =', sch1
                                sch2 = update_spwchan(infile, spw, sch1, truncate=True,
                                                      widths=widths)
                                #print 'sch2 =', sch2
                                ##print 'spwmatch.group() =', spwmatch.group()
                                if sch2:
                                    repl = ''
                                    if sch2 != '*':
                                        repl = "spw='" + sch2 + "'"
                                    cmd = cmds[rownum].replace(spwmatch.group(), repl)
                            #except: # cmd[rownum] no longer applies.
                            except Exception as e:
                                casalog.post(
                                    "Error %s updating row %d of FLAG_CMD" % (e,
                                                                              rownum),
                                             'WARN')
                                casalog.post('sch1 = ' + sch1, 'DEBUG1')
                                casalog.post('cmd = ' + cmd, 'DEBUG1')
                            if cmd != cmds[rownum]:
                                mademod = True
                                cmds[rownum] = cmd
                    if mademod:
                        casalog.post('Updating FLAG_CMD', 'INFO')
                        mytb.putcol('COMMAND', cmds)

                else:
                    casalog.post('FLAG_CMD table contains spw selection by name. Will not update it!','DEBUG')
                
            mytb.close()
            
        except Exception as instance:
            if isopen:
                mytb.close()
            mslocal = None
            mytb = None
            casalog.post("*** Error \'%s\' updating FLAG_CMD" % (instance),
                         'SEVERE')
            return False

    mytb = None

    # Write history to output MS, not the input ms.
    try:
        param_names = sdpolaverage.__code__.co_varnames[:sdpolaverage.__code__.co_argcount]
        param_vals = [eval(p) for p in param_names]
        write_history(mslocal, outfile, 'sdpolaverage', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                     'WARN')
        return False

    mslocal = None
    
    return True
    
 
    
