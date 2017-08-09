import os
import time
from taskinit import *
from casa_stack_manip import stack_frame_find

def plotms(vis=None,
           gridrows=None, gridcols=None,
           rowindex=None,colindex=None,
           plotindex=None,
           xaxis=None, xdatacolumn=None,
           yaxis=None, ydatacolumn=None, yaxislocation=None,
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
           plotrange=None,
           title=None, titlefont=None,
           xlabel=None, xaxisfont=None, ylabel=None, yaxisfont=None,
           showmajorgrid=None, majorwidth=None, majorstyle=None,  majorcolor=None,
           showminorgrid=None, minorwidth=None, minorstyle=None,  minorcolor=None,
           showlegend=None, legendposition=None,
           plotfile=None, expformat=None, exprange=None,
           highres=None, dpi=None, width=None, height=None, overwrite=None,
           showgui=None, clearplots=None,
           callib=None
):

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
        showgui -- Whether or not to display the plotting GUI
                  default: True; example showgui=False

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
    clearplots -- clear existing plots so that the new ones coming in can replace them.
    callib -- calibration library string, list of strings, or filename for on-the-fly calibration

    """
    # Check if DISPLAY environment variable is set.
    if os.getenv('DISPLAY') == None:
        casalog.post('ERROR: DISPLAY environment variable is not set!  Cannot open plotms.', 'SEVERE')
        return False

    if plotfile:
        if not os.path.dirname(plotfile):
            # CAS-7148: Use dir that user cd'ed to in casapy session
            # instead of dir that
            plotfile = os.path.join(os.getcwd(), plotfile)
        if (os.path.exists(plotfile) and not overwrite):
            casalog.post("Plot file " + plotfile + " exists and overwrite is false, cannot write the file", "SEVERE")
            return False


    try:
        # Check synonyms
        # format is:  synonym['new_term'] = 'existing_term', with
        # the existing term being what's coded in PlotMSConstants.h  (case insensitive)
        # CAS-8532 match capitalization in axis names in GUI
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
        synonyms['ant-parallacticangle']=synonyms['ant-parang'] = synonyms['Ant-ParAngle'] = 'ant-parangle'
        synonyms['gamp']=synonyms['gainamp']=synonyms['GainAmp']='Gain Amplitude'
        synonyms['gphase']=synonyms['gainphase']=synonyms['GainPhase']='Gain Phase'
        synonyms['greal']=synonyms['gainreal']=synonyms['GainReal']='Gain Real'
        synonyms['gimag']=synonyms['gainimag']=synonyms['GainImag']='Gain Imag'
        synonyms['del']=synonyms['delay']=synonyms['Delay']='delay'
        synonyms['swp']=synonyms['swpower']=synonyms['switchedpower']=synonyms['SwPower']=synonyms['spgain']='swpower'
        synonyms['tsys']=synonyms['Tsys']=synonyms['TSYS']='tsys'
        synonyms['opac']=synonyms['opacity']=synonyms['Opac']='opac'
        synonyms['snr']=synonyms['SNR']='SNR'
        synonyms['radialvelocity']= synonyms['Radial Velocity'] = 'Radial Velocity [km/s]'
        synonyms['rho']=synonyms['Distance']='Distance (rho) [km]'

        if(xaxis in synonyms): xaxis = synonyms[xaxis]
        if type(yaxis) is str:
            if(yaxis in synonyms): yaxis = synonyms[yaxis]
        elif type(yaxis) is list:
            for index,axis in enumerate(yaxis):
                if (axis in synonyms): yaxis[index] = synonyms[axis]
        if type(coloraxis) is str:
            if(coloraxis in synonyms): coloraxis = synonyms[coloraxis]


        # synonyms for data columns (only one, so just hardcode it)
        if (xdatacolumn=='cor' or xdatacolumn=='corr'):  xdatacolumn='corrected'
        if (ydatacolumn=='cor' or ydatacolumn=='corr'):  ydatacolumn='corrected'

        vis = vis.strip()
        if len(vis) > 0:
            vis = os.path.abspath(vis)

        if not plotindex:
            plotindex = 0
        if plotindex < 0:
            casalog.post("A negative plotindex is not valid.", "SEVERE")
            return False
        if clearplots and plotindex > 0:
            casalog.post("A nonzero plotindex is not valid when clearing plots.", "SEVERE")
            return False

        # start plotms with the procmgr, use logfile
        myframe = stack_frame_find()
        if myframe['casa']['state']['init_version'] > 0:
            from casa_system import procmgr
            pmsrun = procmgr.running("plotms")
            pms = procmgr.fetch("plotms")
            try:
                if pmsrun and not pms.is_alive():  # crash!
                    pms.stop()
                    pmsrun = False
            except AttributeError:  # fetch failed: pms=None
                pass
            if not pmsrun:
                # start (or restart)
                plotmsApp = 'casaplotms'
                for dir in os.getenv('PATH').split(':'):
                    dd = dir + os.sep + plotmsApp
                    if os.path.exists(dd) and os.access(dd, os.X_OK):
                        plotmsApp = dd
                        break
                try:
                    logfile = myframe['casa']['files']['logfile']
                except KeyError:
                    logfile = ""
                procmgr.create("plotms", [plotmsApp, "--nogui", "--nopopups",
                    "--casapy", "--logfilename="+logfile])
                if procmgr.running("plotms"):
                    # connect future calls to this plotms
                    procpid = procmgr.fetch('plotms').pid
                    pm.setPlotmsPid(procpid)

        # Determine whether this is going to be a scripting client or
        # a full GUI supporting user interaction.  This must be done
        # before other properties are set because it affects the
        # constructor of plotms.
        pm.setShowGui( showgui )

        if pm.isDrawing() and clearplots:
            casalog.post("Plotms is running in GUI mode and cannot be run again until the current drawing completes.", "SEVERE")
            return False

        #Clear any existing plots.
        if clearplots:
            pm.clearPlots()

        gridChange = False
        if gridrows > 0 or gridcols > 0:
            gridChange = True
            if not gridrows:
                gridrows = 1
            if not gridcols:
                gridcols = 1

        if gridChange:
            pm.setGridSize( gridrows, gridcols )
        pm.setPlotMSFilename(vis, False, plotindex )

        if type(yaxis) is tuple:
            yaxis = yaxis[0]

        if not yaxis or type(yaxis) is str:
            if not yaxislocation or not type(yaxislocation) is str:
                yaxislocation='left'
            if not ydatacolumn or not type(ydatacolumn) is str:
                ydatacolumn=''
            if not yaxis:
                yaxis = ''

            pm.setPlotAxes( xaxis, yaxis, xdatacolumn, ydatacolumn, yaxislocation, False, plotindex, 0)

        else:
            yAxisCount = len(yaxis)

            yDataCount = 0
            if ydatacolumn!=['']:
                yDataCount = len(ydatacolumn)
            yLocationCount = 0
            if yaxislocation!=['']:
                yLocationCount = len(yaxislocation)

            '''Make sure all the y-axis values are unique.'''
            uniqueY = True
            for i in range(0, yAxisCount ):
                yDataColumnI = ''
                if  i < yDataCount :
                    yDataColumnI = ydatacolumn[i]
                for j in range(0, i):
                    if yaxis[j] == yaxis[i] :
                        yDataColumnJ = ''
                        if j < yDataCount:
                            yDataColumnJ = ydatacolumn[j]
                        if yDataColumnI == yDataColumnJ :
                            uniqueY = False
                            break
                if not uniqueY :
                    break
            if ( uniqueY ):
                for i in range(0,yAxisCount):
                    yDataColumn=''
                    if i < yDataCount:
                        yDataColumn = ydatacolumn[i]
                    yAxisLocation = 'left'
                    if i < yLocationCount:
                        yAxisLocation = yaxislocation[i]
                    pm.setPlotAxes(xaxis, yaxis[i], xdatacolumn, yDataColumn, yAxisLocation, False, plotindex, i)
            else :
                raise Exception('Please remove duplicate y-axes.')

        # Set selection
        if (selectdata and os.path.exists(vis)):
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

        # Set calibration: None, string (filename), dictionary
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
                        casalog.post("Callib file does not exist")
                        raise RuntimeError("Callib file does not exist")
        elif isinstance(callib, list):
            if len(callib[0]) > 0:  # no param is a list (default in plotms.xml?)
                useCallib = True
                callibString = ",".join(callib)
        pm.setPlotMSCalibration(useCallib, callibString, False, plotindex)

        # Set flag extension
        # for now, some options here are not available:
        # pm.setFlagExtension(extendflag, extcorrelation, extchannel, extspw, extantenna, exttime, extscans, extfield)
        extcorrstr=''
        if extcorr:
            extcorrstr='all'
        pm.setFlagExtension(extendflag, extcorrstr, extchannel)

        # Export range
        if not exprange or exprange == "":
            exprange='current'
        pm.setExportRange(exprange)
        if not dpi:
            dpi = -1
        if not width:
            width = -1
        if not height:
            height = -1

        # Set stuff that informs the plot on additional axes
        #  (iteration, colorization, etc.)
        # (Iteration)
        if not iteraxis:
            iteraxis = ""
        if iteraxis=="":
            xselfscale=yselfscale=False
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



          #Determine if there should be a legend.
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



        #Plot range
        if (len(plotrange)!=4):
            if (len(plotrange)==0):
                plotrange=[0.0,0.0,0.0,0.0]
            else:
                raise Exception('plotrange parameter has incorrect number of elements.')

        xrange=plotrange[1]-plotrange[0]
        yrange=plotrange[3]-plotrange[2]

        pm.setXRange((xrange<=0.), plotrange[0],plotrange[1], False, plotindex)
        pm.setYRange((yrange<=0.), plotrange[2],plotrange[3], False, plotindex)

        # Update
        plotUpdated = pm.update()
        if not plotUpdated:
            casalog.post( "There was a problem updating the plot.")
        else:
            # write file if requested
            if(plotfile != ""):
                casalog.post("Plot file " + plotfile, 'NORMAL')
                time.sleep(0.5)
                if (pm.isDrawing()):
                    casalog.post("Will wait until drawing of the plot has completed before exporting it",'NORMAL')
                    while (pm.isDrawing()):
                        time.sleep(1.0)
                casalog.post("Exporting the plot.",'NORMAL')
                casalog.post("Calling pm.save", 'NORMAL')
                plotUpdated = pm.save( plotfile, expformat, highres, dpi, width, height)

    except Exception as instance:
        plotUpdated = False
        print("Exception during plotms task: ", instance)

    return plotUpdated
