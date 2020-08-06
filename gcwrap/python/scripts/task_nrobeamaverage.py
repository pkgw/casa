import os, re
import shutil
import string
import copy
import math
import time
import datetime
import sdutil
import contextlib
from taskinit import mttool, mstool, tbtool, casalog, qa, gentools
from mstools import write_history
from parallel.parallel_data_helper import ParallelDataHelper
import flaghelper as fh
from update_spw import update_spwchan
from callibrary import callibrary

"""
The following code is based on the mstransform code, with 
task name and some task parameters modified. 
To minimise code modification, the parameters used by 
mstransform but not by nrobeamaverage are kept as much as 
possible, and the default values for mstransform are given 
to them.
(CAS-12475, 2019/6/7 WK)
"""

@contextlib.contextmanager
def open_table(table, nomodify=True):
    (tb,) = gentools(['tb'])
    tb.open(table, nomodify=nomodify)
    try:
        yield tb
    finally:
        tb.close()

def nrobeamaverage(
             infile,
             datacolumn,
             field,
             spw, 
             timerange, 
             scan,
             beam,
             timebin,
             outfile):

    casalog.origin('nrobeamaverage')
    try:
        #set temporary data name
        tmpfile = 'tmp-nrobeamaverage-' + os.path.basename(infile.rstrip('/')) + '-' + "{0:%Y%m%d%H%M%S.%f}".format(datetime.datetime.now()) + '.ms'

        #data selection
        do_mst(infile=infile, datacolumn=datacolumn, 
               field=field, spw=spw, timerange=timerange, scan=scan, 
               timebin='0s', outfile=tmpfile)

        #open tmpfile and rewrite antenna column of the ON spectra using beam
        idx_on = None
        with open_table(os.path.join(tmpfile, 'STATE')) as tb:
            ocol = tb.getcol('OBS_MODE')
            for i in range(len(ocol)):
                if ocol[i] == 'OBSERVE_TARGET#ON_SOURCE':
                    idx_on = i
                    break
        if idx_on is None: raise Exception('ON_SOURCE data not found.')

        with open_table(os.path.join(tmpfile, 'ANTENNA')) as tb:
            num_beams = len(tb.getcol('NAME'))
            _beam, min_beamid = get_beamid(beam, num_beams)

        with open_table(tmpfile, nomodify=False) as tb:
            acol = tb.getcol('ANTENNA1')
            scol = tb.getcol('STATE_ID')
            for i in range(len(acol)):
                if (acol[i] in _beam) and (scol[i] == idx_on):
                    acol[i] = min_beamid
            tb.putcol('ANTENNA1', acol)
            tb.putcol('ANTENNA2', acol)

        #time averaging
        do_mst(infile=tmpfile, datacolumn=datacolumn, 
               field='', spw='', timerange='', scan='', 
               timebin=timebin, outfile=outfile)

        add_history(casalog=casalog, infile=infile, datacolumn=datacolumn,
                    field=field, spw=spw, timerange=timerange, scan=scan,
                    timebin=timebin, beam=beam, outfile=outfile)

    except Exception as e:
        casalog.post('Exception from task_nrobeamaverage : ' + str(e), "SEVERE", origin=origin)
    finally:
        #delete tmpfile
        if os.path.isdir(tmpfile): shutil.rmtree(tmpfile)


def get_beamid(beam, num_beams):
    _beam = beam
    try:
        if isinstance(_beam, str):
            _beam = _beam.strip().split(",")
        elif not isinstance(_beam, list):
            raise ValueError('the parameter beam must be list or string.')

        #the default case (beam='')
        if (len(_beam) == 1) and (_beam[0] == ''):
            _beam = []
            for i in range(num_beams): _beam.append(i)
        else:
            for i in range(len(_beam)): _beam[i] = int(_beam[i])
    except Exception as e:
        casalog.post("Error \'%s\' input beam ID is invalid" % (e))
    
    min_beamid = _beam[0]
    for i in range(len(_beam)):
        if _beam[i] < min_beamid: min_beamid = _beam[i]

    return _beam, min_beamid
    


def do_mst(infile, datacolumn, field, spw, timerange, scan, timebin, outfile):
    # followings are parameters of mstransform but not used by nrobeamaverage.
    # just putting default values
    vis = infile             # needed for ParallelDataHelper
    outputvis = outfile      # needed for ParallelDataHelper
    tileshape = [0]
    antenna = ""
    intent = ""
    correlation = ""
    array = ""
    uvrange = ""
    observation = ""
    feed = ""
    realmodelcol = False
    usewtspectrum = False
    chanbin = 1
    mode = "channel"
    start = 0
    width = 1
    timeaverage = False
    timespan = "scan"
    maxuvwdistance = 0.0
    ddistart = -1
    reindex = True
    
    # Initialize the helper class
    pdh = ParallelDataHelper('nrobeamaverage', locals())
    pdh.bypassParallelProcessing(0)
    
    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Create a local copy of the MSTransform tool
    mtlocal = mttool()
    mslocal = mstool()
    
    try:
        # Gather all the parameters in a dictionary.
        config = {}
        
        config = pdh.setupParameters(inputms=infile, outputms=outfile, field=field, 
                    spw=spw, array=array, scan=scan, antenna=antenna, correlation=correlation,
                    uvrange=uvrange, timerange=timerange, intent=intent, observation=str(observation),
                                     feed=feed, taql='')
        
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
        config['tileshape'] = tileshape

        # Only parse timeaverage parameters when timebin > 0s
        """
        if timeaverage:
            tb = qa.convert(qa.quantity(timebin), 's')['value']
            if not tb > 0:
                raise Exception, "Parameter timebin must be > '0s' to do time averaging"
        """
        tbin = qa.convert(qa.quantity(timebin), 's')['value']
        if tbin < 0:
            raise Exception("Parameter timebin must be > '0s' to do time averaging")
        timeaverage = (tbin > 0)
                       
        if timeaverage:
            casalog.post('Parse time averaging parameters')
            config['timeaverage'] = True
            config['timebin'] = timebin
            config['timespan'] = timespan
            config['maxuvwdistance'] = maxuvwdistance
            
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
    
    if (spw != '') and (spw != '*'):
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
    mslocal = None

    return True


def add_history(casalog, infile, datacolumn, field, spw, timerange, scan, timebin, beam, outfile):
    mslocal = mstool()
    # Write history to output MS, not the input ms.
    try:
        param_names = nrobeamaverage.__code__.co_varnames[:nrobeamaverage.__code__.co_argcount]
        param_vals = [eval(p) for p in param_names]
        write_history(mslocal, outfile, 'nrobeamaverage', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                     'WARN')
        return False

    mslocal = None
    
    return True
    
 
    
