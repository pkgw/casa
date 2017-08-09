import os
import re
import string
import time
import shutil
import numpy as np
from taskinit import *
from update_spw import update_spwchan
import flaghelper as fh
from parallel.parallel_data_helper import ParallelDataHelper

def split(vis,
          outputvis,
          keepmms,
          field,
          spw,
          scan,
          antenna,
          correlation,
          timerange,
          intent,
          array,
          uvrange,
          observation,
          feed,
          datacolumn,
          keepflags,
          width,
          timebin,
          combine
          ):

    """Create a visibility subset from an existing visibility set"""

    casalog.origin('split')

    # Initialize the helper class
    pdh = ParallelDataHelper("split", locals())

    # Validate input and output parameters
    try:
        pdh.setupIO()
    except Exception as instance:
        casalog.post('%s'%instance,'ERROR')
        return False

    # Input vis is an MMS
    if pdh.isParallelMS(vis) and keepmms:

        retval = pdh.validateInputParams()
        if not retval['status']:
            raise Exception('Unable to continue with MMS processing')

        pdh.setupCluster('split')

        # Execute the jobs
        try:
            pdh.go()
        except Exception as instance:
            casalog.post('%s'%instance,'ERROR')
            return False

        return True


    # Create local copy of the MSTransform tool
    mtlocal = mttool()

    try:

        # Gather all the parameters in a dictionary.
        config = {}

        if keepflags:
            taqlstr = ''
        else:
            taqlstr = "NOT (FLAG_ROW OR ALL(FLAG))"

        if type(correlation) == list:
            correlation = ', '.join(correlation)
            correlation = correlation.upper()

        config = pdh.setupParameters(inputms=vis, outputms=outputvis, field=str(field),
                    spw=str(spw), array=str(array), scan=str(scan), antenna=str(antenna), correlation=correlation,
                    uvrange=uvrange,timerange=timerange, intent=intent, observation=str(observation),
                    feed=str(feed), taql=taqlstr)

        config['datacolumn'] = datacolumn

        # Channel averaging
        chanaverage = False
        chanbin = width

        # String type
        if isinstance(width, str):
            if width.isdigit():
                chanbin = string.atoi(width)
            else:
                casalog.post('Parameter width is invalid. Using 1 as default', 'WARN')
                chanbin = width = 1

            if chanbin > 1:
                chanaverage = True

        # List type
        elif isinstance(width, list):
            if isinstance(width[0], str):
                if width[0].isdigit():
                    chanbin = list(map(int,width))
                else:
                    casalog.post('Parameter width is invalid. Using 1 as default', 'WARN')
                    chanbin = width = 1

            # If any chanbin in list is > 1, chanaverage=True
            testbin = [i for i in chanbin if i > 1]
            if len(testbin) > 0:
                chanaverage = True

        # Any other type
        if not isinstance(chanbin,str) and not isinstance(chanbin,list):
            casalog.post('Original type(width) is %s'%type(chanbin),'DEBUG')
            if chanbin > 1:
                chanaverage = True

        if chanaverage:
            casalog.post('Parse channel averaging parameters')
            config['chanaverage'] = True
            # verify that the number of spws is the same of the number of chanbin
            pdh.validateChanBin()
            # convert numpy types, until CAS-6493 is not fixed
            chanbin = fh.evaluateNumpyType(chanbin)
            casalog.post('Converted type(width) is %s'%type(chanbin),'DEBUG')
            config['chanbin'] = chanbin

        # Time averaging
        timeaverage = False
        tb = qa.convert(qa.quantity(timebin), 's')['value']
        if tb > 0:
            timeaverage = True

        if timeaverage:
            casalog.post('Parse time averaging parameters')
            config['timeaverage'] = True
            config['timebin'] = timebin
            config['timespan'] = combine
            config['maxuvwdistance'] = 0.0

        # Configure the tool
        casalog.post('%s'%config, 'DEBUG1')
        mtlocal.config(config)

        # Open the MS, select the data and configure the output
        mtlocal.open()

        # Run the tool
        mtlocal.run()

        mtlocal.done()

    except Exception as instance:
        mtlocal.done()
        casalog.post('%s'%instance,'ERROR')
        return False

    # Local copy of ms tool
    mslocal = mstool()

    # Update the FLAG_CMD sub-table to reflect any spw/channels selection
    # If the spw selection is by name or FLAG_CMD contains spw with names, skip the updating

    if ((spw != '') and (spw != '*')) or chanaverage == True:
        isopen = False

        try:
            mytb = tbtool()
            mytb.open(outputvis + '/FLAG_CMD', nomodify=False)
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
                            numspw = len(mslocal.msseltoindex(vis=vis,
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
                                sch2 = update_spwchan(vis, spw, sch1, truncate=True,
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
        param_names = split.__code__.co_varnames[:split.__code__.co_argcount]
        param_vals = [eval(p) for p in param_names]
        write_history(mslocal, outputvis, 'split', param_names,
                      param_vals, casalog)
    except Exception as instance:
        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                     'WARN')
        return False

    mslocal = None

    return True

