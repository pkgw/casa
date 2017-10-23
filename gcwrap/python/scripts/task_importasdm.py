import os
import shutil
from taskinit import *
import flaghelper as fh
import testhelper as th
from casac import casac
from parallel.parallel_data_helper import ParallelDataHelper
import recipes.ephemerides.convertephem as ce

def importasdm(
    asdm=None,
    vis=None,
    createmms=None,
    separationaxis=None,
    numsubms=None,    
    corr_mode=None,
    srt=None,
    time_sampling=None,
    ocorr_mode=None,
    compression=None,
    lazy=None,
    asis=None,
    wvr_corrected_data=None,
    scans=None,
    ignore_time=None,
    process_syspower=None,
    process_caldevice=None,
    process_pointing=None,
    process_flags=None,
    tbuff=None,
    applyflags=None,
    savecmds=None,
    outfile=None,
    flagbackup=None,
    verbose=None,
    overwrite=None,
    showversion=None,
    useversion=None,
    bdfflags=None,
    with_pointing_correction=None,
    remove_ref_undef=None,
    convert_ephem2geo=None,
    polyephem_tabtimestep=None
    ):
    """Convert an ALMA Science Data Model observation into a CASA visibility file (MS) or single-dish data format (Scantable).
           The conversion of the ALMA SDM archive format into a measurement set.  This version
           is under development and is geared to handling many spectral windows of different
           shapes.

           Keyword arguments:
           asdm -- Name of input ASDM file (directory)
               default: none; example: asdm='ExecBlock3'

       vis       -- Root ms or scantable name, note a prefix (.ms or .asap) is NOT appended to this name
           default: none
           
       createmms  -- Create a Multi-MS
           default: False
           
       corr_mode -- correlation mode to be considered on input. Could
            be one or more of the following, ao, co, ac, or all
           default: all

       srt       -- spectral resolution type. Could be one or more of
                    the following, fr, ca, bw, or all
           default: all

       time_sampling -- specifies the time sampling, INTEGRATION and/or
                            SUBINTEGRATION. could be one or more of the following
                            i, si, or all.
           default: all

       ocorr_mode    -- output data for correlation mode AUTO_ONLY 
                            (ao) or CROSS_ONLY (co) or CROSS_AND_AUTO (ca)
           default: ca

      compression  -- produces comrpressed columns in the resulting measurement set.
                 default: False

       lazy         -- Make the MS DATA column read the ASDM Binary data directly
                       (faster import, smaller MS)
                 default: False

       asis         --  creates verbatim copies of the ASDM tables in 
                        the output measurement set. The value given to
                    this option must be a list of table names separated
                    by space characters; the wildcard character '*' is 
                            allowed in table names.

       wvr_corrected_data -- specifies wich values are considered in the 
                      ASDM binary data to fill the DATA column in 
                      the MAIN table of the MS. Expected values for 
                      this option are 'no' for the uncorrected data 
                      (this is the default), 'yes' for the corrected
                      data and 'both' for corrected and uncorrected 
                      data. In the latter case, two measurement sets
                      are created, one containing the uncorrected 
                      data and the other one, whose name is suffixed
                      by '-wvr-corrected', containing the corrected 
                      data.

       scans --  processes only the scans specified in the option's value. This value is a semicolon 
                 separated list of scan specifications. A scan specification consists in an exec bock index 
                 followed by the character ':' followed by a comma separated list of scan indexes or scan 
                 index ranges. A scan index is relative to the exec block it belongs to. Scan indexes are 
                 1-based while exec blocks's are 0-based. "0:1" or "2:2~6" or "0:1,1:2~6,8;2:,3:24~30" "1,2" 
                 are valid values for the option. "3:" alone will be interpreted as 'all the scans of the 
                 exec block#3'. An scan index or a scan index range not preceded by an exec block index will
                 be interpreted as 'all the scans with such indexes in all the exec blocks'.  By default 
                 all the scans are considered.

       ignore_time -- All the rows of the tables Feed, History, Pointing, Source, SysCal, CalDevice, SysPower,
                      and Weather are processed independently of the time range of the selected exec block / scan.

       process_syspower -- The SysPower table is processed if and only if this parameter is set to True.
              default: True

       process_caldevice -- The CalDevice table is processed if and only if this parameter is set to True.
              default: True

       process_pointing -- The Pointing table is processed if and only if this parameter is set to True.
                       If the parameter is set to False the resulting MS will have an empty POINTING table.
              default: True

       process_flags -- Process the online flags and save them to the FLAG_CMD sub-table.
              default: True

            &gt;&gt;&gt; process_flags expandable parameter
                 tbuff -- Time padding buffer (in seconds).
                    default: 0.0

                 applyflags -- Apply the online flags to the MS.
                    default: False

                 savecmds -- Save the online flags to an ASCII file.
                    default: False
                    
                 outfile -- Filename to save the online flags.
                    default: ''

       flagbackup -- Backup the FLAG column in the .flagversions.
              default: True

       verbose     -- produce log output as asdm2MS is being run.

       overwrite -- Over write an existing MS.

       showversion -- report the version of the asdm2MS being used.

       useversion -- Selects the version of asdm2MS to be used (presently only \'v3\' is available).
                     default: v3
                     
       bdfflags -- Set the MS FLAG column according to the ASDM _binary_ flags
                   default: false

       with_pointing_correction -- add (ASDM::Pointing::encoder - ASDM::Pointing::pointingDirection)
                 to the value to be written in MS::Pointing::direction 
                   default: false

       remove_ref_undef -- if set to True then apply fixspwbackport on the resulting MSes.

       convert_ephem2geo -- if True, convert any attached ephemerides to the GEO reference frame

       polyephem_tabtimestep -- Timestep (days) for the tabulation of polynomial ephemerides. A value <= 0 disables tabulation.
                   Presently, VLA data can contain polynomial ephemerides. ALMA data uses tabulated values.
                   default: 0.          

        """

    # Python script
    
    # make agentflagger tool local
    aflocal = casac.agentflagger()

    # make table tool local
    tblocal = casac.table()

    try:
        casalog.origin('importasdm')
        viso = ''
        visoc = ''  # for the wvr corrected version, if needed
        if len(vis) > 0:
            viso = vis
            tmps = vis.rstrip('.ms')
            if tmps == vis:
                visoc = vis + '-wvr-corrected'
            else:
                visoc = tmps + '-wvr-corrected.ms'
        else:
            viso = asdm.rstrip("/") + '.ms'
            visoc = asdm.rstrip("/") + '-wvr-corrected.ms'
            vis = asdm.rstrip("/")



        useversion = 'v3'
        theexecutable = 'asdm2MS'

        execute_string = theexecutable + ' --icm "' + corr_mode \
            + '" --isrt "' + srt + '" --its "' + time_sampling \
            + '" --ocm "' + ocorr_mode + '" --wvr-corrected-data "' \
            + wvr_corrected_data + '" --asis "' + asis \
            + '" --logfile "' + casalog.logfile() + '"'

        if len(scans) > 0:
            execute_string = execute_string + ' --scans ' + scans
        if ignore_time:
            execute_string = execute_string + ' --ignore-time'
        if useversion == 'v3':
            if not process_syspower:
                execute_string = execute_string + ' --no-syspower'
            if not process_caldevice:
                execute_string = execute_string + ' --no-caldevice'
            if not process_pointing:
                execute_string = execute_string + ' --no-pointing'

        if compression:
            execute_string = execute_string + ' --compression'
        elif lazy:
            execute_string = execute_string + ' --lazy'
            
        if verbose:
            execute_string = execute_string + ' --verbose'
#         if not overwrite and os.path.exists(viso):
#             raise Exception, \
#                 'You have specified an existing MS and have indicated you do not wish to overwrite it'

        # Compression
        if compression:
                   # viso = viso + '.compressed'
            viso = viso.rstrip('.ms') + '.compressed.ms'
            visoc = visoc.rstrip('.ms') + '.compressed.ms'

        vistoproc = [] # the output MSs to post-process
        if wvr_corrected_data == 'no' or wvr_corrected_data == 'both':
            vistoproc.append(viso)
        if (wvr_corrected_data == 'yes' or wvr_corrected_data == 'both') : 
            vistoproc.append(visoc)

        for ff in vistoproc:
            if not overwrite and os.path.exists(ff):
                raise Exception('You have specified an existing MS and have indicated you do not wish to overwrite it: %s'%ff)

        # If viso+".flagversions" then process differently depending on the value of overwrite..
        #
        if flagbackup:
            for myviso in vistoproc:
                dotFlagversion = myviso + '.flagversions'
                if os.path.exists(dotFlagversion):
                    if overwrite:
                        casalog.post("Found '" + dotFlagversion
                                     + "' . It'll be deleted before running the filler."
                                     )
                        os.system('rm -rf %s' % dotFlagversion)
                    else:
                        casalog.post("Found '%s' but can't overwrite it."
                                     % dotFlagversion)
                        raise Exception("Found '%s' but can't overwrite it." \
                            % dotFlagversion)
               
        # Make outfile always a list             
        if isinstance(outfile, str):
            if outfile == '': 
                outfile = []
            else:
                noutfile = [outfile]
                outfile = noutfile
            
        if savecmds:
            if len(outfile) == 0:
                # Create default names for the online flags
                for myviso in vistoproc:
                    outfile.append(myviso.replace('.ms','_cmd.txt'))
            elif len(outfile) != len(vistoproc):
                casalog.post('List of outfile names does not match list of MSs','WARN')
                casalog.post('Will save online flags to temporary filenames', 'WARN')
                outfile = []
                for myviso in vistoproc:
                    online_file = myviso.replace('.ms','_TEMP_cmd.txt')
                    outfile.append(online_file)
                                     
            if not overwrite:
                for of in outfile:
                    if os.path.exists(of):
                        raise Exception("Cannot overwrite online flags file '%s'; overwrite is set to False."% of)
                
            
        execute_string = execute_string + ' ' + asdm + ' ' + viso

        if showversion:
            casalog.post("You set option \'showversion\' to True. Will just output the version information and then terminate."
                         , 'WARN')
            execute_string = theexecutable + ' --revision'

        if with_pointing_correction:
            execute_string = execute_string + ' --with-pointing-correction'

        if (polyephem_tabtimestep!=None) and (type(polyephem_tabtimestep)==int or type(polyephem_tabtimestep)==float):
            if polyephem_tabtimestep>0:
                casalog.post('Will tabulate all attached polynomial ephemerides with a time step of '
                             +str(polyephem_tabtimestep)+' days.')
                if polyephem_tabtimestep>1.:
                    casalog.post('A tabulation timestep of <= 1 days is recommended.', 'WARN')
                execute_string = execute_string + ' --polyephem-tabtimestep '+str(polyephem_tabtimestep)

        casalog.post('Running ' + theexecutable
                     + ' standalone invoked as:')
        # print execute_string
        casalog.post(execute_string)
        exitcode = os.system(execute_string)
        if exitcode != 0:
            if not showversion:
                casalog.post(theexecutable
                             + ' terminated with exit code '
                             + str(exitcode), 'SEVERE')
                raise Exception('ASDM conversion error. Please check if it is a valid ASDM and that data/alma/asdm is up to date.')

        if showversion:
            return
        
        #
        # Possibly remove the name of the measurement set expected to contain the corrected data from the list of of produced measurement
        # sets if it appears the filler did not find any corrected data.
        #
        if not os.path.exists(visoc):
            vistoproc = [myviso for myviso in vistoproc if myviso != visoc]

        # CAS-7369. HISTORY should be written after createmms is tested
        #
        # Populate the HISTORY table of the MS with information about the context in which it's been created
        #
        try: 
            mslocal = mstool() 
            param_names = importasdm.__code__.co_varnames[:importasdm.__code__.co_argcount] 
            param_vals = [eval(p) for p in param_names]

            for myviso in vistoproc:
                write_history(mslocal, myviso, 'importasdm', param_names, 
                              param_vals, casalog) 

        except Exception as instance: 
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 
                         'WARN')
            return False 

        if mslocal:
            mslocal = None 
            
        # 
        # Do we apply fixspwbackport
        if remove_ref_undef :
            casalog.post('remove_ref_undef=True: fixspwbackport will be applied ...')
            
            for myviso in vistoproc:
                cmd = 'fixspwbackport ' + myviso
                casalog.post('Running fixspwbackport standalone invoked as:')
                casalog.post(cmd)
                cmdexitcode = os.system(cmd)

                if cmdexitcode != 0:
                    casalog.post(cmd
                                 + ' terminated with exit code '
                                 + str(cmdexitcode), 'SEVERE')
                    raise Exception('fixspwbackport error.')

        # Binary Flag processing
        if bdfflags:
            
            casalog.post('Parameter bdfflags==True: flags from the ASDM binary data will be used to set the MS flags ...')
            
            bdffexecutable = 'bdflags2MS '
            bdffexecstring_base = bdffexecutable + ' -f ALL' + ' --ocm "' + ocorr_mode \
            + '" --logfile "' + casalog.logfile() + '"'
 
            if len(scans) > 0:
                bdffexecstring_base = bdffexecstring_base + ' --scans ' + scans

            if lazy and not compression:
                bdffexecstring_base = bdffexecstring_base + ' --lazy=true'

            for myviso in vistoproc:
                if myviso.find("wvr-corrected") != -1:
                    options = " --wvr-corrected=True " 
                else:
                    options = " "

                bdffexecstring = bdffexecstring_base + options + asdm + ' ' + myviso

                casalog.post('Running '+bdffexecutable+' standalone invoked as:')
                casalog.post(bdffexecstring)

                bdffexitcode = os.system(bdffexecstring)
                if bdffexitcode != 0:
                    casalog.post(bdffexecutable
                                 + ' terminated with exit code '
                                 + str(bdffexitcode), 'SEVERE')
                    raise Exception('ASDM binary flags conversion error. Please check if it is a valid ASDM and that data/alma/asdm is up to date.')


        theephemfields = ce.findattachedephemfields(myviso,field='*')
        if len(theephemfields)>0: 
            # until asdm2MS does this internally: recalc the UVW coordinates for ephem fields
            imt = imtool()
            imt.open(myviso, usescratch=False)
            imt.calcuvw(theephemfields, refcode='J2000', reuse=False)
            imt.close()

        if convert_ephem2geo:
            for myviso in vistoproc:
                ce.convert2geo(myviso, '*') # convert any attached ephemerides to GEO
        
        if len(theephemfields)>0: 
            # also set the direction column in the SOURCE table
            tblocal.open(myviso+'/FIELD', nomodify=False)
            sourceids = tblocal.getcol('SOURCE_ID')
            ftimes = tblocal.getcol('TIME')
            ftimekw = tblocal.getcolkeywords('TIME')
            tmpa = tblocal.getcol('PHASE_DIR')
            origphasedir = tmpa

            affectedsids = []
            thesamplefields = []
            for fld in theephemfields: # determine all source ids used by the ephem fields
                if not (sourceids[fld] in affectedsids): # this source id wasn't handled yet
                    affectedsids.append(sourceids[fld])
                    thesamplefields.append(fld)
                    # need to temporarily change the offset (not all mosaics have an element at (0,0))
                    tmpa[0][0][fld]=0.
                    tmpa[1][0][fld]=0.
                #endif
            #endfor
            tblocal.putcol('PHASE_DIR', tmpa)
            tblocal.close()

            tblocal.open(myviso+'/SOURCE')
            sourceposref = tblocal.getcolkeywords('DIRECTION')['MEASINFO']['Ref']
            tblocal.close()

            directions = []
            msmdlocal = casac.msmetadata()
            msmdlocal.open(myviso)
            
            for fld in thesamplefields:
                thedirmeas = msmdlocal.phasecenter(fld)
                if thedirmeas['refer']!=sourceposref:
                    casalog.post('Ephemeris is in '+thedirmeas['refer']+' instead of '+sourceposref
                                 +' frame.\nEntry in SOURCE table will be converted to '+sourceposref, 'WARN')
                    melocal = metool()
                    melocal.doframe(thedirmeas)
                    thedirmeas = melocal.measure(thedirmeas, sourceposref)

                directions.append([thedirmeas['m0']['value'], thedirmeas['m1']['value']])
                thetime = me.epoch(v0=str(ftimes[fld])+'s', rf=ftimekw['MEASINFO']['Ref'])
                casalog.post("Will set SOURCE direction for SOURCE_ID "+str(sourceids[fld])
                             +" to ephemeris phase center for time "+str(thetime['m0']['value'])+" "+thetime['m0']['unit']+" "+thetime['refer']) 
            #endfor
            msmdlocal.close()
             
            # restore original PHASE_DIR
            tblocal.open(myviso+'/FIELD', nomodify=False)
            tblocal.putcol('PHASE_DIR', origphasedir)
            tblocal.close()

            # write source directions
            tblocal.open(myviso+'/SOURCE', nomodify=False)
            ssourceids = tblocal.getcol('SOURCE_ID')
            sdirs = tblocal.getcol('DIRECTION')
            for row in range(0,len(ssourceids)):
                for i in range(0,len(affectedsids)):
                    if ssourceids[row]==affectedsids[i]:
                        sdirs[0][row] = directions[i][0]
                        sdirs[1][row] = directions[i][1]
                        break
                #endfor
            #endfor
            tblocal.putcol('DIRECTION', sdirs) # write back corrected directions
            tblocal.close()
                
        #end if        

        ##############################################################################################3
        # CAS-7369 - Create an output Multi-MS (MMS)
        if createmms:
            # Get the default parameters of partition
            from tasks import partition
            fpars = partition.parameters
            for mypar in list(fpars.keys()):
                fpars[mypar] = partition.itsdefault(mypar)
                
            # Call the cluster for each MS
            for myviso in vistoproc:
                casalog.origin('importasdm')
                
                # Move original MS to tempdir
                tempname = myviso+'.temp.ms'
                outputmms = myviso
                shutil.move(myviso, tempname)
                
                # Get the proper column
                datacolumn = 'DATA'
                dcols = ['DATA', 'FLOAT_DATA']
                for dc in dcols:
                    if len(th.getColDesc(tempname, dc)) > 0:
                        datacolumn = dc
                        break
                    
                fpars['datacolumn'] = datacolumn
                    
                casalog.post('Will create a Multi-MS for: '+myviso)
                
                fpars['vis'] =  tempname
                fpars['flagbackup'] =  False 
                fpars['outputvis'] = outputmms
                fpars['separationaxis'] = separationaxis
                fpars['numsubms'] = numsubms
                pdh = ParallelDataHelper('partition', fpars) 
            
                # Get a cluster
                pdh.setupCluster(thistask='partition')
                try:
                    pdh.go()
                    
                    # Remove original MS
                    shutil.rmtree(tempname)

                except Exception as instance:
                    # Restore MS in case of error in MMS creation
                    shutil.move(tempname, myviso)
                    casalog.post('%s'%instance,'ERROR')
                    return False
                
            casalog.origin('importasdm')

        # Create a .flagversions for the MS or MMS
        if flagbackup:
            for myviso in vistoproc:
                if os.path.exists(myviso):
                    aflocal.open(myviso)
                    aflocal.saveflagversion('Original',
                            comment='Original flags at import into CASA',
                            merge='save')
                    aflocal.done()
                
        # Importasdm Flag Parsing
        if os.access(asdm + '/Flag.xml', os.F_OK):
            # Find Flag.xml
            casalog.post('Found Flag.xml in SDM')
            
            # Find Antenna.xml
            if os.access(asdm + '/Antenna.xml', os.F_OK):
                casalog.post('Found Antenna.xml in SDM')

            else:
                raise Exception('Failed to find Antenna.xml in SDM')
            
            # Find SpectralWindow.xml
            if os.access(asdm + '/SpectralWindow.xml', os.F_OK):
                casalog.post('Found SpectralWindow.xml in SDM')

            else:
                raise Exception('Failed to find SpectralWindow.xml in SDM')
                    
            #
            # Parse Flag.xml into flag dictionary
            #
            if process_flags:
                flagcmds = fh.parseXML(asdm, float(tbuff))
                onlinekeys = list(flagcmds.keys())
                nflags = onlinekeys.__len__()
                                
                # Apply flags to the MS
                if nflags > 0:
                    idx = 0
                    for myviso in vistoproc:
                        if applyflags:
                            # Open the MS and attach it to the tool
                            aflocal.open(myviso)
                            # Select the data
                            aflocal.selectdata()
                            # Setup the agent's parameters
                            fh.parseAgents(aflocal, flagcmds, [], True, True, '')
                            # Initialize the agents
                            aflocal.init()
                            # Run the tool
                            aflocal.run(True, True)
                            casalog.post('Applied %s flag commands to %s'%(nflags,myviso))
                            # Destroy the tool and de-attach the MS
                            aflocal.done()
                            # Save to FLAG_CMD table. APPLIED is set to True.
                            fh.writeFlagCommands(myviso, flagcmds, True, '', '', True)       
                        else:
                            casalog.post('Will not apply flags to %s (apply_flags=False), use flagcmd to apply'%myviso)

                            # Write to FLAG_CMD, APPLIED is set to False
                            fh.writeFlagCommands(myviso, flagcmds, False, '', '', True)
                    
                        # Save the flag cmds to an ASCII file
                        if savecmds:
                            # Save to standard filename
                            fh.writeFlagCommands(myviso, flagcmds, False, '', outfile[idx], False)
                            casalog.post('Saved %s flag commands to %s'%(nflags,outfile[idx]))
                            idx += 1
                    
                else:
                    casalog.post('There are no flag commands to process')
                
        else:
            casalog.post('There is no Flag.xml in ASDM', 'WARN')

        
        return
    
    except Exception as instance:

        print('*** Error ***', instance)


