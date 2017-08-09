from __future__ import absolute_import, division, print_function
from taskinit import *
import glob

# Concatenation of ephemeris tables
#
# Example:
#      import recipes.ephemerides.concatephem as cce
#      cce.concatephem(ephems=['Neptune1.tab', 'Neptune2.tab', 'Neptune3.tab'],
#                              outputephem='NeptuneAll.tab')


def concatephem(ephems=[], outputephem=''):
    """
    concatephem

    Concatenate the given ephemeris tables into a single output table
    filling gaps with dummy rows and removing overlap rows.
    Before concatenation, test that basic conditions are fulfilled:
    same time grid, list of input ephemerides is ordered in time.

    ephems - List of the ephemeris tables to be concatenated
            default:
    outputephem - Name of the output ephemeris to be created from the concatenation
            If empty, concatephem will only perform a dryrun.
            default: ''

    """
    mytb = tbtool()
    starts = []
    ends = []
    stepsizes = []
    hasoverlap_with_previous = []
    gap_to_previous = [] # negative entry means overlap
    shouldconcat_with_previous = []
    canconcat_with_previous = []

    stepsize_rel_tolerance = 1E-6 # relative difference in time step size must be smaller than this
    stepsize_abs_tolerance_d = 0.1/86400. # absolute difference in time step size (days) must be smaller than this

    for myephem in ephems:
        mytb.open(myephem)
        mynrows = mytb.nrows()
        if mynrows==0:
            mytb.close()
            casalog.post('Ephemeris '+myephem+' has no rows.', 'SEVERE')
            return
        mystart = mytb.getcell('MJD',0)
        myend = mytb.getcell('MJD', mynrows-1)
        if mynrows<2:
            mystepsize = 0
            casalog.post('Ephemeris '+myephem+' has only one row.', 'WARN')
        else:
            mystepsize =  mytb.getcell('MJD',1) - mystart

        mytb.close()

        starts.append(mystart)
        ends.append(myend)
        stepsizes.append(mystepsize)

    if os.path.exists(outputephem):
        casalog.post('Output ephemeris table '+outputephem+' exists already.', 'SEVERE')
        return

    casalog.post('Ephem no., startMJD, endMJD, step size (d)', 'INFO')
    for i in range(0,len(starts)):
        casalog.post(str(i)+', '+str(starts[i])+', '+str(ends[i])+', '+str(stepsizes[i]), 'INFO')

    backstep = 0
    for i in range(0,len(starts)):
        shouldconcat_with_previous.append(False)
        canconcat_with_previous.append(False)
        hasoverlap_with_previous.append(False)
        gap_to_previous.append(0)
        if i>0:
            if (abs(stepsizes[i] - stepsizes[i-1-backstep])/stepsizes[i] < stepsize_rel_tolerance) \
                    and abs(stepsizes[i] - stepsizes[i-1-backstep]) < stepsize_abs_tolerance_d:
                casalog.post( 'Ephemerides '+str(i-1-backstep)+' and '+str(i)+' have same step size.', 'INFO')
                if starts[i-1-backstep] <= starts[i]:
                    casalog.post( 'Ephemeris '+str(i-1-backstep)+' begins before '+str(i), 'INFO')
                    if ends[i-1-backstep] < ends[i]:
                        casalog.post( 'Ephemeris '+str(i-1-backstep)+' ends before '+str(i), 'INFO')
                        shouldconcat_with_previous[i] = True
                        numsteps_to_add = (starts[i]-ends[i-1-backstep])/stepsizes[i] - 1
                        gap_to_previous[i] = int(round(numsteps_to_add))
                        if abs(round(numsteps_to_add) - numsteps_to_add)<1E-4:
                            casalog.post( 'Gap between ephemerides '+str(i-1-backstep)+' and '+str(i)+' is '+str(gap_to_previous[i])+' steps', 'INFO')
                            canconcat_with_previous[i] = True
                            backstep = 0
                            if numsteps_to_add<0.:
                                casalog.post(  'Ephemerides '+str(i-1-backstep)+' and '+str(i)+' have overlap of '+str(-gap_to_previous[i])+' steps', 'INFO')
                                hasoverlap_with_previous[i]=True
                        else:
                            casalog.post( 'Gap between ephemerides '+str(i-1-backstep)+' and '+str(i)+' is not an integer number of steps:', 'WARN')
                            casalog.post( str(round(numsteps_to_add) - numsteps_to_add), 'INFO')
                            canconcat_with_previous[i] = False
                            backstep += 1
                    else:
                        casalog.post( 'Ephemeris '+str(i-1-backstep)+' does not end before '+str(i), 'INFO')
                        shouldconcat_with_previous[i] = False
                        canconcat_with_previous[i] = False
                        backstep += 1
                else:
                    casalog.post( 'Ephemeris '+str(i-1-backstep)+' does not begin before or at the same time as '+str(i), 'INFO')
                    casalog.post('Ephemerides are in wrong order.', 'SEVERE')
                    return
            else:
                casalog.post( 'Ephemerides '+str(i-1-backstep)+' and '+str(i)+' do not have same step size.', 'INFO')
                casalog.post('Ephemerides have inhomogeneous steps sizes.', 'SEVERE')
                return
        #end if


    if outputephem=='' or len(starts)<2:
        return
    else:
        casalog.post( 'Creating output ephemeris ...', 'INFO')

        os.system('cp -R '+ephems[0]+' '+outputephem)

        mytb2 = tbtool()

        try:
            mytb.open(outputephem, nomodify=False)
        except:
            casalog.post('Error opening table '+outputepehem+' for writing.', 'SEVERE')
            return


        for i in range(1,len(starts)):

            if shouldconcat_with_previous[i] and canconcat_with_previous[i]:

                mynrows = mytb.nrows()

                mytb2.open(ephems[i])
                mynrows2 = mytb2.nrows()
                startrow = 0
                if(hasoverlap_with_previous[i]):
                    startrow = -gap_to_previous[i]
                elif(gap_to_previous[i]>0): # fill gap with dummy rows
                    mytb.addrows(gap_to_previous[i])
                    startmjd = mytb.getcell('MJD', mynrows-1)
                    for j in range(mynrows, mynrows+gap_to_previous[i]):
                        newmjd = startmjd+stepsizes[i]*(j-mynrows+1)
                        mytb.putcell('MJD', j, newmjd)
                        mytb.putcell('RA', j, 0.)
                        mytb.putcell('DEC', j, 0.)
                        mytb.putcell('Rho', j, 0.)
                        mytb.putcell('RadVel', j, 0.)
                        mytb.putcell('diskLong', j, 0.)
                        mytb.putcell('diskLat', j, 0.)

                    casalog.post( str(gap_to_previous[i])+' dummy rows appended to fill gap', 'INFO')

                mynrows = mytb.nrows()

                mytb.addrows(mynrows2-startrow)
                for j in range(mynrows, mynrows+mynrows2-startrow):
                    fromrow = j - mynrows + startrow
                    mytb.putcell('MJD', j, mytb2.getcell('MJD', fromrow))
                    mytb.putcell('RA', j, mytb2.getcell('RA', fromrow))
                    mytb.putcell('DEC', j, mytb2.getcell('DEC', fromrow))
                    mytb.putcell('Rho', j, mytb2.getcell('Rho', fromrow))
                    mytb.putcell('RadVel', j, mytb2.getcell('RadVel', fromrow))
                    mytb.putcell('diskLong', j, mytb2.getcell('diskLong', fromrow))
                    mytb.putcell('diskLat', j, mytb2.getcell('diskLat', fromrow))

                casalog.post( str(mynrows2-startrow)+' rows copied over', 'INFO')

            else:

                casalog.post( 'Ephemeris '+str(i)+' skipped', 'INFO')

            #endif
        #endfor
    #endif

    return


def findephems(vis=[], field=''):

    """
       findephems

       Search the MSs given in the list "vis" for ephemerides for
       a given field and return the list of paths to the ephemeris tables
       in the same order as vis.

       vis - list of the MSs to search for ephemerides
             default: []

       field - field for which to seach ephemerides
             default:
    """

    if vis==[] or field=='':
        return []

    if type(vis) == str:
        vis = [vis]

    mytb=tbtool()

    rval = [] # list of ephemeris tables to be returned

    for visname in vis:
        if not os.path.exists(visname+'/FIELD'):
            casalog.post('MS '+visname+' does not exist.', 'SEVERE')
            return []

        mytb.open(visname+'/FIELD')
        fnames = mytb.getcol('NAME')

        thefield = field
        if type(thefield) == int:
            thefield = str(fnames[thefield])
        if type(thefield) != str:
            casalog.post('Parameter field must be str or int.', 'SEVERE')
            mytb.close()
            return []

        if thefield in fnames:
            colnames = mytb.colnames()
            if not 'EPHEMERIS_ID' in colnames:
                casalog.post('MS '+visname+' has no ephemerides attached.', 'WARN')
                rval.append('')
                mytb.close()
                continue

            ephids = mytb.getcol('EPHEMERIS_ID')
            mytb.close()
            theephid = ephids[list(fnames).index(thefield)]

            thetabs = glob.glob(visname+'/FIELD/EPHEM'+str(theephid)+'_*')

            if len(thetabs)==0:
                casalog.post('MS '+visname+' has no ephemerides for field '+thefield+' attached.', 'WARN')
                rval.append('')
                continue

            if len(thetabs)>1:
                casalog.post('MS '+visname+' has more than one ephemeris for field '+thefield+' attached.', 'SEVERE')
                return[]

            rval.append(thetabs[0])

        else:
             casalog.post('MS '+visname+' has no field '+thefield, 'WARN')
             rval.append('')
             continue
    #endfor

    return rval


def concatreplaceephem(vis=[], field='', commontime=False):

    """
       Search the MSs given in the list "vis" for ephemerides for
       a given field, concatenate the ephemerides, and replace
       all of the original ephemerides with the concatenated one.

       vis - list of the MSs to search for ephemerides
             default: []

       field - field for which to seach ephemerides
             default:

       commontime - set the FIELD table reference time in all input MSs to a common time.
             If True, the common time is taken from the time of the first MS.
             If a float value is provided, it is interpreted as a time in seconds in the native
             CASA time representation.
             default: False - do not modify the FIELD table TIME column

    """

    ephemfield = field
    thetabs = findephems(vis, ephemfield)

    if not (type(commontime)==bool) and not (type(commontime)==float):
        raise Exception('ERROR: parameter commontime must be bool or float')

    if thetabs != [] and not ('' in thetabs):
        tmptab = os.path.basename(thetabs[0])+'.concattmp'
        concatephem(thetabs, tmptab)
        if os.path.exists(tmptab):
            for targettab in thetabs:
                if not os.path.exists(targettab):
                    raise Exception('Internal ERROR: ephemeris '+targettab+' does not exist')
                os.system('rm -rf '+targettab)
                os.system('cp -R '+tmptab+' '+targettab)
        else:
            casalog.post('ERROR while concatenating ephemerides for field '+str(ephemfield), 'SEVERE')
            raise Exception('Concatenation of ephemerides for field '+str(ephemfield)+' failed.')

        os.system('rm -rf '+tmptab)

        if type(commontime)==float or commontime==True:
            casalog.post('Will set FIELD table reference time for field '+str(ephemfield)+' to common value:', 'INFO')
            thecommontime = commontime
            if type(commontime)==float:
                casalog.post('   '+str(thecommontime), 'INFO')
            mytb = tbtool()
            for thevis in vis:
                mytb.open(thevis+'/FIELD', nomodify=False)
                thenames = mytb.getcol('NAME')
                thetimes = mytb.getcol('TIME')

                for i in range(0,len(thenames)):
                    if (thenames[i]==ephemfield):
                        if thecommontime==True and thevis==vis[0]:
                            thecommontime = thetimes[i] # memorize the common time
                            casalog.post('   '+str(thecommontime), 'INFO')
                        else:
                            thetimes[i] = thecommontime
                #end for
                mytb.putcol('TIME', thetimes)
            #end for
            mytb.close()

    else:
        casalog.post('ERROR while searching for ephemerides for field '+str(ephemfield), 'SEVERE')
        raise Exception('Cannot find ephemerides for field '+str(ephemfield)+' in all input MSs.')

    casalog.post('All ephemerides for field '+str(ephemfield)+' replaced by concatenated ephemeris.', 'INFO')

    return True
