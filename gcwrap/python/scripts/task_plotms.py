import os
import time
from taskinit import *
from casa_system import procmgr

def plotms(vis=None, 
           gridrows=None, gridcols=None,
           rowindex=None,colindex=None,
           plotindex=None,
           xaxis=None, xdatacolumn=None, xframe=None, xinterp=None,
           yaxis=None, ydatacolumn=None, yframe=None, yinterp=None, yaxislocation=None,
           selectdata=None, field=None, spw=None,
           timerange=None, uvrange=None, antenna=None, scan=None,
           correlation=None, array=None, observation=None, 
           intent=None, feed=None, msselect=None,
           averagedata=None,
           avgchannel=None, avgtime=None, avgscan=None, avgfield=None,
           avgbaseline=None, avgantenna=None, avgspw=None, scalar=None,
           transform=None,
           freqframe=None,restfreq=None,veldef=None,shift=None,
           extendflag=None,
           extcorr=None, extchannel=None,
           iteraxis=None,xselfscale=None,yselfscale=None,
           xsharedaxis=None, ysharedaxis=None,
           customsymbol=None, symbolshape=None, symbolsize=None,
           symbolcolor=None, symbolfill=None, symboloutline=None,
           coloraxis=None,
           customflaggedsymbol=None, flaggedsymbolshape=None,
           flaggedsymbolsize=None, flaggedsymbolcolor=None,
           flaggedsymbolfill=None, flaggedsymboloutline=None,
           xconnector=None, timeconnector=False,
           plotrange=None,
           title=None, titlefont=None, 
           xlabel=None, xaxisfont=None, ylabel=None, yaxisfont=None,
           showmajorgrid=None, majorwidth=None, majorstyle=None,  majorcolor=None,    
           showminorgrid=None, minorwidth=None, minorstyle=None,  minorcolor=None, 
           showlegend=None, legendposition=None,   
           plotfile=None, expformat=None, verbose=True, exprange=None,
           highres=None, dpi=None, width=None, height=None, overwrite=None,
           showgui=None, clearplots=None,
           callib=None, headeritems=None, showatm=None, showtsky=None, showimage=None):
# we'll add these later
#           extspw=None, extantenna=None,
#           exttime=None, extscans=None, extfield=None,

    """
    
            Task for plotting and interacting with visibility data.  A variety
        of axes choices (including data column) along with MS selection and
        averaging options are provided for data selection.  Flag extension
        parameters are also available for flagging operations in the plotter.
        
            All of the provided parameters can also be set using the GUI once
        the application has been launched.  Additional and more specific
        operations are available through the GUI and/or through the plotms
        tool (pm).
        

    Keyword arguments:
    vis -- input visibility dataset
           default: ''
    
    gridrows -- Number of subplot rows.
                    default: 1
    gridcols -- Number of subplot columns.
                    default: 1 
    rowindex -- Row location of the subplot (0-based).
                    default: 0
    colindex -- Column location of the subplot (0-based).
                    default: 0          
    plotindex -- Index to address a subplot (0-based).
                    default: 0            
    xaxis, yaxis -- what to plot on the two axes
                    default: '' (uses PlotMS defaults/current set).
        &gt;&gt;&gt; xaxis, yaxis expandable parameters
        xdatacolumn, 
        ydatacolumn -- which data column to use for data axes
                       default: '' (uses PlotMS default/current set).
        xframe,
        yframe      -- which coordinates frame to use for ant-ra,ant-dec axes
                       default: '' (uses PlotMS default/current set).
        xinterp,
        yinterp     -- which interpolation method to use for ant-ra,ant-dec axes
                       default: '' (uses PlotMS default/current set).
        yaxislocation -- whether the data should be plotted using the left or right y-axis
                       default: '' (uses PlotMS default).
    iteraxis -- what axis to iterate on when doing iteration plots
                default: ''
              &gt;&gt;&gt; xsharedaxis, ysharedaxis, xselfscale, yselfscale expandable parameters 
        xselfscale -- If true, iterated plots should share a common x-axis label per column.
                       default: False.
        yselfscale -- If true, iterated plots should share a common y-axis label per row.
                       default: False.
        xsharedaxis -- use a common x-axis for vertically aligned plots (must also set xselfscale=True)
                        default: False.
        ysharedaxis -- use a common y-axis for horizontally aligned plots (must also set yselfscale=True)
                        default: False.
    selectdata -- data selection parameters flag
                  (see help par.selectdata for more detailed information)
                  default: False
      &gt;&gt;&gt; selectdata expandable parameters
        field -- select using field ID(s) or field name(s)
                 default: '' (all).
        spw -- select using spectral window/channels
               default: '' (all)
        timerange -- select using time range
                     default: '' (all).
        uvrange -- select using uvrange
                   default: '' (all).
        antenna -- select using antenna/baseline
                   default: '' (all).
        scan -- select using scan number
                default: '' (all).
        correlation -- select using correlations
                       default: '' (all).
        array -- select using (sub)-array range
                 default: '' (all).
        observation -- select by observation ID(s).
                 default: '' (all).
        intent -- select observing intent  
                  default: ''  (no selection by intent)  
                  intent='*BANDPASS*'  (selects data labelled with  
                                        BANDPASS intent)
        feed   -- select feed ID
                  default: '' (all)
                  feed='1~2'
        msselect -- TaQL selection expression
                    default: '' (all).
    
    averagedata -- data averaging parameters flag
                   default: False.
      &gt;&gt;&gt; averagedata expandable parameters
        avgchannel -- average over channel?  either blank for none, or a value
                      in channels.
                      default: '' (none).
        avgtime -- average over time?  either blank for none, or a value in
                   seconds.
                   default: '' (none).
        avgscan -- average over scans?  only valid if time averaging is turned
                   on.
                   default: False.
        avgfield -- average over fields?  only valid if time averaging is
                    turned on.
                    default: False.
        avgbaseline -- average over all baselines?  mutually exclusive with
                       avgantenna.
                       default: False.
        avgantenna -- average by per-antenna?  mutually exclusive with
                      avgbaseline.
                      default: False.
        avgspw -- average over all spectral windows?
                  default: False.
    
    extendflag -- have flagging extend to other data points?
                  default: False.
      &gt;&gt;&gt; extendflag expandable parameters
        extcorr -- extend flags based on correlation?  blank = none.
                          default: ''.
        extchannel -- extend flags based on channel?
                      default: False.
        extspw -- extend flags based on spw?
                  default: False.
        extantenna -- extend flags based on antenna?  should be either blank,
                      'all' for all baselines, or an antenna-based value.
                      default: ''.
        exttime -- extend flags based on time (within scans)?
                   default: False.
        extscans -- extend flags based on scans?  only valid if time extension
                    is turned on.
                    default: False.
        extfield -- extend flags based on field?  only valid if time extension
                    is turned on.
                    default: False.

    coloraxis -- which axis to use for colorizing
                     default: ''  (ignored - same as colorizing off)              
    
    title  -- title along top of plot (called "canvas" in some places)
    titlefont -- plot title font size
                 default: 0 (autosize depending on grid)
    exprange -- Export all iteration plots ('all') or only the current one.
                    default: '' (only export the current iteration plot)
    xlabel, ylabel -- text to label horiz. and vert. axes, with formatting (%% and so on)
    xaxisfont, yaxisfont -- int for axis font size
    
    showlegend -- show a legend on the plot
                    default: False
    legendposition -- position for the legend.  Legends can be interior or exterior to the plot
                    Interior legends can be located in the upper right, lower right, upper left, or lower left.
                    Exterior legends can be located on the right, left, top, or bottom.
                    default: 'upperright'
    showgui -- Whether or not to display the plotting GUI
                    default: True; example showgui=False
    clearplots -- clear existing plots so that the new ones coming in can replace them.                 
    callib -- calibration library string, list of strings, or filename for on-the-fly calibration
    headeritems -- string of comma-separated page header items keywords
    showatm -- show atmospheric transmission curve
    showtsky -- show sky temperature curve
    showimage -- show image sideband curve

    """
    # Check if DISPLAY environment variable is set.
    if os.getenv('DISPLAY') == None:
        casalog.post('ERROR: DISPLAY environment variable is not set! Cannot run plotms.', 'SEVERE')
        return False

    # using procmgr?
    usingprocmgr = False
    if casa['state']['init_version'] > 0:
        usingprocmgr = True
        if not procIsRunning("dbus"):
            casalog.post('ERROR: dbus-daemon has stopped, cannot run plotms. Please restart casa.', 'SEVERE')
            return False
 
    # check arguments
    # check plotfile for export
    if plotfile:
        if not os.path.dirname(plotfile):
            # CAS-7148: Use dir that user cd'ed to in casapy session
            # instead of dir that 
            plotfile = os.path.join(os.getcwd(), plotfile)
        if (os.path.exists(plotfile) and not overwrite):
            casalog.post("Plot file " + plotfile + " exists and overwrite is false, cannot write the file", "SEVERE")
            return False

    # Define axis synonyms
    # format is:  synonym['new_term'] = 'existing_term'
    # existing_term in PlotMSConstants.h
    # CAS-8532: match capitalization in axis names in GUI
    if True:
        synonyms = {}
        synonyms['Scan'] = 'scan'
        synonyms['Field'] = 'field'
        synonyms['Time'] = 'time'
        synonyms['timeinterval'] = synonyms['timeint'] = synonyms['time_interval'] = synonyms['Interval'] = 'interval'
        synonyms['Spw'] = 'spw'
        synonyms['chan'] = synonyms['Channel'] = 'channel'
        synonyms['freq'] = synonyms['Frequency'] = 'frequency'
        synonyms['vel'] = synonyms['Velocity'] = 'velocity'
        synonyms['correlation'] = synonyms['Corr'] = 'corr'
        synonyms['ant1'] = synonyms['Antenna1'] = 'antenna1'
        synonyms['ant2'] = synonyms['Antenna2'] = 'antenna2'
        synonyms['Baseline'] = 'baseline'
        synonyms['Row'] = 'row'
        synonyms['Observation'] = 'observation'
        synonyms['Intent'] = 'intent'
        synonyms['Feed1'] = 'feed1'
        synonyms['Feed2'] = 'feed2'
        synonyms['amplitude'] = synonyms['Amp'] = 'amp'
        synonyms['Phase'] = 'phase'
        synonyms['Real'] = 'real'
        synonyms['imaginary'] = synonyms['Imag'] = 'imag'
        synonyms['weight'] = synonyms['Wt'] = synonyms['Weight'] = 'wt'
        synonyms['wtamp'] = synonyms['Wt*Amp'] = 'wtamp'
        synonyms['weightspectrum'] = synonyms['WtSp'] = synonyms['WeightSpectrum'] = 'wtsp'
        synonyms['Sigma'] = 'sigma'
        synonyms['sigmaspectrum'] = synonyms['SigmaSpectrum'] = synonyms['SigmaSp'] = 'sigmasp'
        synonyms['Flag'] = 'flag'
        synonyms['FlagRow'] = 'flagrow'
        synonyms['UVdist'] = 'uvdist'
        synonyms['uvdistl'] = synonyms['uvdist_l']=synonyms['UVwave'] = 'uvwave'
        synonyms['U'] = 'u'
        synonyms['V'] = 'v'
        synonyms['W'] = 'w'
        synonyms['Uwave'] = 'uwave'
        synonyms['Vwave'] = 'vwave'
        synonyms['Wwave'] = 'wwave'
        synonyms['Azimuth'] = 'azimuth'
        synonyms['Elevation'] = 'elevation'
        synonyms['hourang'] = synonyms['HourAngle'] = 'hourangle'
        synonyms['parang'] = synonyms['parallacticangle'] = synonyms['ParAngle'] = 'parangle'
        synonyms['ant'] = synonyms['Antenna'] = 'antenna'
        synonyms['Ant-Azimuth'] = 'ant-azimuth'
        synonyms['Ant-Elevation'] = 'ant-elevation'
        synonyms['Ant-Ra'] = synonyms['Ant-RA'] = 'ant-ra'
        synonyms['Ant-Dec'] = synonyms['Ant-DEC'] = 'ant-dec'
        synonyms['ant-parallacticangle']=synonyms['ant-parang'] = synonyms['Ant-ParAngle'] = 'ant-parangle'
        synonyms['gamp']=synonyms['gainamp']=synonyms['GainAmp']='Gain Amp'
        synonyms['gphase']=synonyms['gainphase']=synonyms['GainPhase']='Gain Phase'
        synonyms['greal']=synonyms['gainreal']=synonyms['GainReal']='Gain Real'
        synonyms['gimag']=synonyms['gainimag']=synonyms['GainImag']='Gain Imag'
        synonyms['del']=synonyms['delay']=synonyms['Delay']='delay'
        synonyms['swp']=synonyms['swpower']=synonyms['switchedpower']=synonyms['SwPower']=synonyms['spgain']='swpower'
        synonyms['tsys']=synonyms['Tsys']=synonyms['TSYS']='tsys'
        synonyms['opac']=synonyms['opacity']=synonyms['Opac']='opac'
        synonyms['snr']=synonyms['SNR']='SNR'
        synonyms['antpos']='Antenna Positions'
        synonyms['radialvelocity']= synonyms['Radial Velocity'] = 'Radial Velocity [km/s]'
        synonyms['rho']=synonyms['Distance']='Distance (rho) [km]'
        # data columns: unspecified residuals default to vector
        synonyms['residual']=synonyms['corrected-model']='corrected-model_vector'
        synonyms['data-model']='data-model_vector'
        synonyms['corrected/model']='corrected/model_vector'
        synonyms['data/model']='data/model_vector'

    if True:  # ant-ra/ant-dec axes parameters: Python/C++ parameters maps
        # Reference Frames
        cpp_radec_frame = {}
        for py_radec_ref_frame in ['icrs','j2000','b1950','azelgeo','galactic']:
            cpp_radec_frame[py_radec_ref_frame] = py_radec_ref_frame.upper()
        # Interpolation Methods
        cpp_radec_interp = {}
        for py_radec_interp in ['nearest','cubic spline']:
            cpp_radec_interp[py_radec_interp] = ' '.join(
                [word.capitalize() for word in py_radec_interp.split()]  )
        cpp_radec_interp['spline']=cpp_radec_interp['cubic spline']
    try:
        # Do preliminary checks on argument values
        # Set synonyms to existing_terms
        if(xaxis in synonyms):
            xaxis = synonyms[xaxis]
        if isinstance(yaxis, str):
            if yaxis in synonyms:
                yaxis = synonyms[yaxis]
        elif isinstance(yaxis, list):
            for index,axis in enumerate(yaxis):
                if axis in synonyms:
                    yaxis[index] = synonyms[axis]

        if isinstance(coloraxis, str):
            if coloraxis in synonyms:
                coloraxis = synonyms[coloraxis]

        if(xdatacolumn in synonyms):
            xdatacolumn = synonyms[xdatacolumn]
        if isinstance(ydatacolumn, str):
            if ydatacolumn in synonyms:
                yaxis = synonyms[ydatacolumn]
        elif isinstance(ydatacolumn, list):
            for index,col in enumerate(ydatacolumn):
                if col in synonyms:
                    ydatacolumn[index] = synonyms[col]

        for param_name in ['xframe','yframe']:
            param_py_value = eval(param_name)
            if isinstance(param_py_value, str):
                if param_py_value in cpp_radec_frame:
                    exec('{p_name} = cpp_radec_frame[param_py_value]'.format(
                            p_name=param_name)
                    )
            elif isinstance(param_py_value, list):
                for index,frame in enumerate(param_py_value):
                    if frame in cpp_radec_frame:
                        exec('{p_name}[{index}] = cpp_radec_frame[frame]'.format(
                            p_name=param_name,
                            index=index)
                        )

        for param_name in ['xinterp','yinterp']:
            param_py_value = eval(param_name)
            if isinstance(param_py_value, str):
                if param_py_value in cpp_radec_interp:
                    exec('{p_name} = cpp_radec_interp[param_py_value]'.format(
                            p_name=param_name)
                    )
            elif isinstance(param_py_value, list):
                for index,interp in enumerate(param_py_value):
                    if interp in cpp_radec_interp:
                        exec('{p_name}[{index}] = cpp_radec_interp[interp]'.format(
                            p_name=param_name,
                            index=index)
                        )

        # check vis exists
        vis = vis.strip()
        if len(vis) > 0:
            vis = os.path.abspath(vis)
            if not os.path.exists(vis):
                casalog.post('\n'.join(['Input file not found:',vis]),"SEVERE")
                return False

        # check plotindex
        if not plotindex:
            plotindex = 0  
        if plotindex < 0:
            casalog.post("A negative plotindex is not valid.", "SEVERE")
            return False  
        if clearplots and plotindex > 0:   
            casalog.post("A nonzero plotindex is not valid when clearing plots.", "SEVERE")
            return False
       
        # start plotms with the procmgr, use casa logfile
        if usingprocmgr:
            if not procIsRunning("plotms"):
                # set plotmsApp
                plotmsApp = 'casaplotms'
                for dir in os.getenv('PATH').split(':'):
                    dd = dir + os.sep + plotmsApp
                    if os.path.exists(dd) and os.access(dd, os.X_OK):
                        plotmsApp = dd
                        break
                # set logfile
                try:
                    logfile = casa['files']['logfile']
                except KeyError:
                    logfile = ""
                # start process with procmgr
                procmgr.create("plotms", [plotmsApp, "--nogui", "--nopopups",
                    "--casapy", "--logfilename="+logfile])
            if procIsRunning("plotms"):
                # plotms is running, but is it connected? CAS-11306
                procmgrPid = procmgr.fetch('plotms').pid
                plotmsPid = pm.getPlotMSPid()
                if plotmsPid != procmgrPid:
                    # connect future pm calls to procmgr plotms
                    pm.setPlotMSPid(procmgrPid)
            else:
                casalog.post("plotms failed to start", "SEVERE")
                return False

        # Determine whether this is going to be a scripting client or 
        # a full GUI supporting user interaction.  This must be done 
        # before any other properties are set because it affects the
        # constructor of plotms.
        if casa['flags'].nogui:
            showgui = False
        pm.setShowGui( showgui )

        if clearplots:
            # Clear any existing plots unless still drawing last one
            if pm.isDrawing():
                casalog.post("Plotms is running in GUI mode and cannot be run again until the current drawing completes.", "SEVERE")
                return False
            pm.clearPlots()

        # set grid
        gridChange = False    
        if gridrows > 0 or gridcols > 0:
            gridChange = True
            if not gridrows:
                gridrows = 1
            if not gridcols:
                gridcols = 1
        if gridChange:
            pm.setGridSize( gridrows, gridcols )

        # set vis filename
        pm.setPlotMSFilename(vis, False, plotindex )

        # set yaxis defaults as needed
        if isinstance(yaxis, tuple):
            # yaxis default from xml is str or list: yaxis=('', [])
            # set it to the empty string
            yaxis = yaxis[0]
        if not yaxis or isinstance(yaxis, str):
            if not yaxis:
                yaxis = ''
            if yaxis == 'ant-ra' or yaxis == 'ant-dec':
                # Handle empty lists as empty strings
                if isinstance(yinterp, list) and not yinterp:
                    yinterp = ''
                if isinstance(yframe, list) and not yframe:
                    yframe = ''
                if isinstance(yinterp, str) and isinstance(yframe, str):
                    # For now, ignore cases where xinterp or xframe is a list
                    if isinstance(xframe, list):
                        msg_fmt = "Assuming xframe={assumed} instead of xframe={org}"
                        assumed_xframe = '' if not xframe else xframe[0]
                        msg = msg_fmt.format(assumed=assumed_xframe, org=xframe)
                        casalog.post(msg,'WARN','set_axes')
                        xframe = assumed_xframe
                    if isinstance(xinterp, list):
                        msg_fmt = "Assuming xinterp={assumed} instead of xinterp={org}"
                        assumed_xinterp = '' if not xinterp else xinterp[0]
                        msg = msg_fmt.format(assumed=assumed_xinterp, org=xinterp)
                        casalog.post(msg,'WARN','set_axes')
                        xinterp = assumed_xinterp
                    if isinstance(yaxislocation, list):
                        yaxislocation= 'left' if not yaxislocation else yaxislocation[0]
                    if not isinstance(yaxislocation, str):
                        yaxislocation= 'left'
                    xdatacolumn = ydatacolumn = ''
                    pm.setPlotAxes(xaxis, yaxis, xdatacolumn, ydatacolumn,
                    xframe, yframe, xinterp, yinterp,
                    yaxislocation,
                    False, plotindex, 0)
                else:
                    # Handle yinterp, yframe and yaxislocation as parallel lists, which
                    # 1. must have the same length
                    # 2. must NOT contain empty strings, otherwise C++ vectors won't have the same length
                    if isinstance(yinterp, list):
                        if isinstance(yframe, str):
                            # Allow usage: plotms(yinterp=['nearest','spline'],yframe='')
                            if not yframe:
                                yframe = 'icrs'
                            yframe = [yframe for i in yinterp]
                        else:
                            if len(yframe) != len(yinterp):
                                msg_fmt = "Length mismatch: yframe={0} and yinterp={1}"
                                msg = msg_fmt.format(yframe,yinterp)
                                casalog.post(msg,'ERROR',set_axes)
                                return False
                        if isinstance(yaxislocation, str):
                            # Allow usage: plotms(yinterp=['nearest','spline'],yaxislocation='')
                            if not yaxislocation:
                                yaxislocation = 'left'
                            yaxislocation = [yaxislocation for i in yinterp]
                        else:
                            if len(yaxislocation) != len(yinterp):
                                msg_fmt = "Length mismatch: yaxislocation={0} and yinterp={1}"
                                msg = msg_fmt.format(yaxislocation,yinterp)
                                casalog.post(msg,'ERROR',set_axes)
                                return False
                        # For now: enforce xframe=yframe, xinterp=yinterp in this case
                        casalog.post('Enforcing xframe=yframe, xinterp=yinterp','WARN','set_axes')
                        xdatacolumn = ydatacolumn = 'data'
                        if not xaxis:
                            xaxis = 'time'
                        for dataindex, (frame,interp,yaxisloc) in enumerate(zip(yframe,yinterp,yaxislocation)):
                            pm.setPlotAxes(xaxis, yaxis, xdatacolumn, ydatacolumn,
                                           frame, frame, interp, interp,
                                           yaxisloc,
                                           False, plotindex, dataindex)
                    else:
                        casalog.post('Not yet implemented: yframe=list','SEVERE','set_axes')
                        return False
            else:
                if not yaxislocation or not isinstance(yaxislocation, str):
                    yaxislocation='left'
                if not ydatacolumn or not isinstance(ydatacolumn, str):
                    ydatacolumn=''
                pm.setPlotAxes(xaxis, yaxis, xdatacolumn, ydatacolumn,
                    xframe, yframe, xinterp, yinterp,
                    yaxislocation,
                    False, plotindex, 0)
        else:
            # make ydatacolumn and yaxislocation same length as yaxis
            # and check that no duplicate y axes
            yAxisCount = len(yaxis)
            yDataCount = 0
            if ydatacolumn!=['']:
                yDataCount = len(ydatacolumn)
            yLocationCount = 0
            if yaxislocation!=['']:
                yLocationCount = len(yaxislocation)
            '''Make sure all the y-axis values are unique.'''
            uniqueY = True
            for i in range( yAxisCount ):
                yDataColumnI = ''
                if  i < yDataCount :
                    yDataColumnI = ydatacolumn[i]
                for j in range(i):
                    if yaxis[j] == yaxis[i] :
                        yDataColumnJ = ''
                        if j < yDataCount:
                            yDataColumnJ = ydatacolumn[j]
                        if yDataColumnI == yDataColumnJ :
                            # same axis, same datacolumn!
                            uniqueY = False
                            break
                if not uniqueY :
                    break
            if ( uniqueY ):
                for i in range(yAxisCount):
                    yDataColumn=''
                    if i < yDataCount:
                        yDataColumn = ydatacolumn[i]
                    yAxisLocation = 'left'
                    if i < yLocationCount:
                        yAxisLocation = yaxislocation[i]
                    if xaxis in ['ant-ra','ant-dec'] or yaxis[i]  in ['ant-ra','ant-dec']:
                        raise Exception('Currently not supported: multiple y-axes involving ant-ra or ant-dec')
                    # Always make C++ ra/dec parameters vectors the same length as yaxis
                    xframe = yframe = 'icrs'
                    xinterp = yinterp = 'nearest'
                    pm.setPlotAxes(xaxis, yaxis[i], xdatacolumn, yDataColumn, 
                        xframe, yframe, xinterp, yinterp,
                        yAxisLocation,
                        False, plotindex, i)
            else :
                raise Exception('Please remove duplicate y-axes.')

        if not showatm:
            showatm = False
        if not showtsky:
            showtsky = False
        if not showimage:
            showimage = False
        if showatm and showtsky:
            casalog.post('You have selected both showatm and showtsky.  Defaulting to showatm=True only.', "WARN")
            showtsky = False
        if showatm or showtsky:  # check that xaxis is None, chan, or freq
            validxaxis = not xaxis or xaxis in ["channel", "frequency"]
            if not validxaxis:
                casalog.post('showatm and showtsky are only valid when xaxis is channel or frequency', 'SEVERE')
                return False
        if showimage and (not showatm and not showtsky):
            casalog.post('Defaulting to showimage=False because showatm and showtsky are False.', "WARN")
            showimage = False
        pm.setShowCurve(showatm, showtsky, showimage, False, plotindex)

        # Set selection
        if selectdata:
            pm.setPlotMSSelection(field, spw, timerange, uvrange, antenna, scan,
                                  correlation, array, str(observation), intent,
                                  feed, msselect, False, plotindex)
        else:
            pm.setPlotMSSelection('', '', '', '', '', '', '', '', '', '', '',
                                  '', False, plotindex)
       
        # Set averaging
        if not averagedata:
            avgchannel = avgtime = ''
            avgscan = avgfield = avgbaseline = avgantenna = avgspw = False
            scalar = False
        if avgbaseline and avgantenna:
            casalog.post('Averaging over baselines is mutually exclusive with per-antenna averaging.', "SEVERE")
            return False
        if avgchannel and (float(avgchannel) < 0.0):
            casalog.post('Cannot average negative number of channels', "SEVERE")
            return False
        try:
            if avgtime and (float(avgtime) < 0.0):
                casalog.post('Cannot average negative time value', "SEVERE")
                return False
        except ValueError:
            casalog.post('avgtime value must be numerical string in seconds (no units)', "SEVERE")
            return False
        pm.setPlotMSAveraging(avgchannel, avgtime, avgscan, avgfield, avgbaseline, 
                              avgantenna, avgspw, scalar, False, plotindex)

        # Set transformations
        if not transform:
            freqframe=''
            restfreq=''
            veldef='RADIO'
            shift=[0.0,0.0]
        pm.setPlotMSTransformations(freqframe,veldef,restfreq,shift[0],shift[1],
                                    False, plotindex)

        # Set calibration: None or string (filename)
        useCallib = False
        callibString = ''
        if isinstance(callib, str):
            # Determine if filename or string of params
            if '=' in callib:
                useCallib = True
                callibString = callib
            else:
                callibFile = callib.strip()
                if len(callibFile) > 0:
                    callibFile = os.path.abspath(callib)
                    if os.path.exists(callibFile):
                        useCallib = True
                        callibString = callibFile
                    else:
                        casalog.post("Callib file does not exist", "SEVERE")
                        return False
        elif isinstance(callib, list): # default is callib=['']
            if len(callib[0]) > 0:  # argument set to list of strings
                useCallib = True
                callibString = ",".join(callib)
        pm.setPlotMSCalibration(useCallib, callibString, False, plotindex) 

        # Set flag extensions; for now, some options here are not available
        # pm.setFlagExtension(extendflag, extcorrelation, extchannel, extspw,
        #    extantenna, exttime, extscans, extfield)
        extcorrstr=''
        if extcorr:
            extcorrstr='all'
        pm.setFlagExtension(extendflag, extcorrstr, extchannel)
        
        # Export range
        if not exprange or exprange == "":
            exprange='current'
        pm.setExportRange(exprange)
        # for pm.save:
        if not dpi:
            dpi = -1
        if not width:
            width = -1
        if not height:
            height = -1

        # Set additional axes (iteration, colorization, etc.)
        # (Iteration)
        if not iteraxis:
            iteraxis = ""
        if iteraxis=="":
            xselfscale = yselfscale = False
            xsharedaxis = ysharedaxis = False
        if not rowindex:
            rowindex = 0
        if not colindex:
            colindex = 0
        if not xselfscale:
            xselfscale = False
        if not yselfscale:
            yselfscale = False
        if not xsharedaxis:
            xsharedaxis = False
        if not ysharedaxis:
            ysharedaxis = False
        if not xselfscale and xsharedaxis:
            casalog.post("Plots cannot share an x-axis unless they use the same x-axis scale.", "ERROR")
            return False
        if not yselfscale and ysharedaxis:
            casalog.post( "Plots cannot share a y-axis unless they use the same y-axis scale.", "ERROR")
            return False
        if xsharedaxis and gridrows < 2:
            casalog.post( "Plots cannot share an x-axis when gridrows=1.", "WARN")
            xsharedaxis=False
        if ysharedaxis and gridcols < 2:
            casalog.post( "Plots cannot share a y-axis when gridcols=1.", "WARN")
            ysharedaxis=False
        pm.setPlotMSIterate(iteraxis,rowindex,colindex,
                            xselfscale,yselfscale,
                            xsharedaxis,ysharedaxis,False,plotindex);
        
        # (Colorization)
        if coloraxis:
            pm.setColorAxis(coloraxis,False,plotindex)

        # Set custom symbol
        # Make the custom symbol into a list if it is not already.
        if type(customsymbol) is bool and customsymbol:
            customSymbolValue = customsymbol
            customsymbol=[customSymbolValue]
            
        if type(symbolshape) is str:
            symbolValue = symbolshape
            symbolshape=[symbolValue]
            
        if type(symbolsize) is int:
            symbolValue = symbolsize
            symbolsize=[symbolValue]    
        
        if type(symbolcolor) is str:
            symbolValue = symbolcolor
            symbolcolor=[symbolValue]  
            
        if type(symbolfill) is str:
            symbolValue = symbolfill
            symbolfill=[symbolValue]
            
        if type(symboloutline) is bool:
            symbolValue = symboloutline
            symboloutline=[symbolValue]                   
        
        if type(customsymbol) is list:
            customSymbolCount = len(customsymbol)
            for i in range(0,customSymbolCount):
                
                if  i >= len(symbolshape) or not symbolshape[i]:
                    symbolShapeI = 'autoscaling'
                else:
                    symbolShapeI = symbolshape[i]
                symbolShape = symbolShapeI 
                
                if customsymbol[i]:
                    if i >=len(symbolsize) or not symbolsize[i]:
                        symbolSizeI = 2
                    else:
                        symbolSizeI = symbolsize[i]
                    symbolSize = symbolSizeI
                    
                    if i>=len(symbolcolor) or not symbolcolor[i]:
                        symbolColorI = '0000ff'
                    else:
                        symbolColorI = symbolcolor[i]
                    symbolColor = symbolColorI
                    
                    if i>=len(symbolfill) or not symbolfill[i]:
                        symbolFillI = 'fill'
                    else:
                        symbolFillI = symbolfill[i]
                    symbolFill = symbolFillI
                    
                    if type( symboloutline) is bool:
                        symbolOutlineI = symboloutline
                    elif type(symboloutline) is list:
                        if i>=len(symboloutline) or not symboloutline[i]:
                            symbolOutlineI=False
                        else:
                            symbolOutlineI = symboloutline[i]
                    else:
                        symbolOutlineI = False        
                    symbolOutline = symbolOutlineI
                    
                else:
                    symbolSize = 2
                    symbolColor = '0000ff'
                    symbolFill = 'fill'
                    symbolOutline = False
                pm.setSymbol(symbolShape, symbolSize, symbolColor,
                     symbolFill, symbolOutline, False,plotindex,i)
            
        # Set custom flagged symbol
        if type(customflaggedsymbol) is bool and customflaggedsymbol:
            customSymbolValue = customflaggedsymbol
            customflaggedsymbol=[customSymbolValue]
            
        if type(flaggedsymbolshape) is str:
            symbolValue = flaggedsymbolshape
            flaggedsymbolshape=[symbolValue]
            
        if type(flaggedsymbolsize) is int:
            symbolValue = flaggedsymbolsize
            flaggedsymbolsize=[symbolValue]    
        
        if type(flaggedsymbolcolor) is str:
            symbolValue = flaggedsymbolcolor
            flaggedsymbolcolor=[symbolValue]  
            
        if type(flaggedsymbolfill) is str:
            symbolValue = flaggedsymbolfill
            flaggedsymbolfill=[symbolValue]
            
        if type(flaggedsymboloutline) is bool:
            symbolValue = flaggedsymboloutline
            flaggedsymboloutline=[symbolValue]  
        
        if type(customflaggedsymbol) is list:
            customSymbolCount = len(customflaggedsymbol)
            for i in range(0,customSymbolCount):
                if i>=len(flaggedsymbolshape) or not flaggedsymbolshape[i]:
                    flaggedSymbolShapeI = 'nosymbol'
                else:
                    flaggedSymbolShapeI = flaggedsymbolshape[i]
                flaggedSymbolShape = flaggedSymbolShapeI
                
                if customflaggedsymbol[i]:
                    
                    if i >=len(flaggedsymbolsize) or not flaggedsymbolsize[i]:
                        flaggedSymbolSizeI = 2
                    else:
                        flaggedSymbolSizeI = flaggedsymbolsize[i]
                    flaggedSymbolSize = flaggedSymbolSizeI
                    
                    if i >=len(flaggedsymbolcolor) or not flaggedsymbolcolor[i]:
                        flaggedSymbolColorI = 'ff0000'
                    else:
                        flaggedSymbolColorI = flaggedsymbolcolor[i]
                    flaggedSymbolColor = flaggedSymbolColorI
                    
                    if i>=len(flaggedsymbolfill) or not flaggedsymbolfill[i]:
                        flaggedSymbolFillI = 'fill'
                    else:
                        flaggedSymbolFillI = flaggedsymbolfill[i]
                    flaggedSymbolFill = flaggedSymbolFillI
                    
                    if type(flaggedsymboloutline) is bool:
                        flaggedSymbolOutlineI = flaggedsymboloutline
                    elif type(flaggedsymboloutline) is list:
                        if i>=len(flaggedsymboloutline) or not flaggedsymboloutline[i]:
                            flaggedSymbolOutlineI = False
                        else:
                            flaggedSymbolOutlineI = flaggedsymboloutline[i]
                    else:
                        flaggedSymbolOutlineI = False
                    flaggedSymbolOutline = flaggedSymbolOutlineI
                else:
                    flaggedSymbolSize = 2
                    flaggedSymbolColor = 'ff0000'
                    flaggedSymbolFill = 'fill'
                    flaggedSymbolOutline = False    
                pm.setFlaggedSymbol(flaggedSymbolShape, flaggedSymbolSize,
                            flaggedSymbolColor, flaggedSymbolFill,
                            flaggedSymbolOutline, False, plotindex, i)
 
        # connect the dots
        if not xconnector:
            xconnector = 'none'
        if not timeconnector:
            timeconnector = False
        pm.setConnect(xconnector, timeconnector, False, plotindex)

        # Legend
        if not showlegend:
            showlegend = False
        if not legendposition:
            legendposition = 'upperRight' 
        pm.setLegend( showlegend, legendposition, False, plotindex )          
        
        # Set various user-directed appearance parameters
        pm.setTitle(title,False,plotindex)
        pm.setTitleFont(titlefont,False,plotindex)
        pm.setXAxisLabel(xlabel,False,plotindex)
        pm.setXAxisFont(xaxisfont,False,plotindex)
        pm.setYAxisLabel(ylabel,False,plotindex)
        pm.setYAxisFont(yaxisfont,False,plotindex)
        pm.setGridParams(showmajorgrid, majorwidth, majorstyle, majorcolor,
                         showminorgrid, minorwidth, minorstyle, minorcolor, False, plotindex)

        # Plot ranges
        if len(plotrange) == 0:
            plotrange=[0.0, 0.0, 0.0, 0.0]
        elif len(plotrange) != 4:
            casalog.post('plotrange parameter has incorrect number of elements.', 'SEVERE')
            return False
        else:
            try:
                for i,val in enumerate(plotrange):
                    plotrange[i] = float(val)
            except (TypeError, ValueError) as e:
                casalog.post("plotrange elements must be numeric", 'SEVERE')
                return False
        xranges = plotrange[1] - plotrange[0]
        yranges = plotrange[3] - plotrange[2]
        pm.setXRange((xranges<=0.0), plotrange[0], plotrange[1], False, plotindex)
        pm.setYRange((yranges<=0.0), plotrange[2], plotrange[3], False, plotindex)
        
        # Page Header Items
        # Python keywords for specifying header items are defined in CAS-8082, 
        # Erik's comment dated 9-jun-2016
        # Python / C++ header items keywords map
        # format is header_cpp_kw['python_keyword'] = 'c++_keyword', where 
        # the c++ keyword is what's coded in PlotMSPageHeaderParam.h
        header_cpp_kw = {}
        header_cpp_kw['filename'] = 'filename'
        header_cpp_kw['ycolumn']  = 'y_columns'
        header_cpp_kw['obsdate']  = 'obs_start_date'
        header_cpp_kw['obstime']  = 'obs_start_time'
        header_cpp_kw['observer'] = 'obs_observer'
        header_cpp_kw['projid']   = 'obs_project'
        header_cpp_kw['telescope'] = 'obs_telescope_name'
        header_cpp_kw['targname'] = 'target_name'
        header_cpp_kw['targdir']  = 'target_direction'
        
        if type(headeritems) is str:
            cpp_headeritems = []
            for headeritem_word in headeritems.split(','):
                py_headeritem = headeritem_word.strip()
                if py_headeritem == "":
                    continue
                if py_headeritem in header_cpp_kw:
                    ccp_headeritem = header_cpp_kw[py_headeritem]
                    cpp_headeritems.append(ccp_headeritem)
                else:
                    casalog.post("Ignoring invalid page header item: " + py_headeritem ,"WARN")
    
            pm.setPlotMSPageHeaderItems(','.join(cpp_headeritems), False, plotindex)

        # Update - ready to plot!
        plotUpdated = pm.update()
        if not plotUpdated:
            casalog.post( "There was a problem updating the plot.", "ERROR")
        else:
            # write file if requested
            if(plotfile != ""):
                casalog.post("Plot file " + plotfile, 'NORMAL')
                # kluge: isDrawing checks if *any* thread is running, could be cache
                # thread or drawing thread! Give it time for cache to finish...
                time.sleep(0.5)
                if (pm.isDrawing()):
                    casalog.post("Will wait until drawing of the plot has completed before exporting it",'NORMAL')
                    while (pm.isDrawing()):
                        time.sleep(1.0)
                casalog.post("Exporting the plot.",'NORMAL')
                plotUpdated = pm.save( plotfile, expformat, verbose, highres, dpi, width, height)

    except Exception as instance:
        plotUpdated = False
        print("Exception during plotms task: ", instance)
        
    if not plotUpdated:
        checkProcesses() # see if something crashed, log failure
    return plotUpdated

def procIsRunning(procname):
    procrun = procmgr.running(procname)
    if not procrun:
        time.sleep(2) # for slow startups
        procrun = procmgr.running(procname)
    if procrun:
        try:
            process = procmgr.fetch(procname)
            if not process.is_alive():  # crash!
                process.stop()  # let procmgr know it crashed
                procrun = False
        except AttributeError:  # fetch returned None: None.is_alive()
            pass
    return procrun

def checkProcesses():
    if not procIsRunning('plotms'):
        casalog.post( "plotms has stopped running. Check logs for error and run again.", "SEVERE")

    if not procIsRunning('dbus'):
        casalog.post( "dbus-daemon has stopped running.  Please restart casa.", "SEVERE")
