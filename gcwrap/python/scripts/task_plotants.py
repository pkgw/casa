import os
import numpy as np
import pylab as pl
from textwrap import wrap
from taskinit import msmdtool, gentools, qatool, casalog, mstool

def plotants(vis=None, figfile=None, 
        antindex=None, logpos=None, 
        exclude=None, checkbaselines=None,
        title=None):
    """Plot the antenna distribution in the local reference frame:

    The location of the antennas in the MS will be plotted with
    X-toward local east; Y-toward local north.  The name of each 
    antenna is shown next to its respective location.

    Keyword arguments:
    vis -- Name of input visibility file.
            default: none. example: vis='ngc5921.ms'

    figfile -- Save the plotted figure in this file.
            default: ''. example: figfile='myFigure.png'

    antindex -- Label antennas with name and antenna ID
            default: False. example: antindex=True

    logpos -- Produce a logarithmic position plot.
            default: False. example: logpos=True

    exclude -- antenna IDs or names to exclude from plotting
            default: []. example: exclude=[2,3,4], exclude='DV15'

    checkbaselines -- Only plot antennas in the MAIN table.
          This can be useful after a split.  WARNING:  Setting 
          checkbaselines to True will add to runtime in
          proportion to the number of rows in the dataset.
            default: False. example: checkbaselines=True

    title -- Title written along top of plot
            default: ''

    You can zoom in by pressing the magnifier button (bottom,
    third from right) and making a rectangular region with
    the mouse.  Press the home button (left most button) to
    remove zoom.

    A hard-copy of this plot can be obtained by pressing the
    button on the right at the bottom of the display. A file 
    dialog will allow you to choose the directory, filename,
    and format of the export.
    """

    try:
        # remove trailing / for title basename
        if vis[-1]=='/':
            vis = vis[:-1]
        if isinstance(exclude, tuple):
            exclude = list(exclude)
            if not all(isinstance(ant, (int, str)) for ant in exclude):
                casalog.post("'exclude' list must contain str or int", "ERROR")
                return
        # make exclude a list
        if isinstance(exclude, int):
            exclude = [exclude]
        elif isinstance(exclude, str):
            myms = mstool()
            try:
                exclude = myms.msseltoindex(vis, baseline=exclude)['antenna1'].tolist()
            except RuntimeError as rterr:  # MSSelection failed
                mesg = "Antenna plot failed: exclude " + str(rterr)
                casalog.post(mesg, "ERROR")
                return

        pl.clf()
        telescope, names, ids, xpos, ypos = getAntennaInfo(vis, exclude, 
            checkbaselines)
        if not names:
            casalog.post("No antennas selected.  Exiting plotants.", "ERROR")
            return
        if logpos:
            plotAntennasLog(telescope, names, ids, xpos, ypos, antindex)
        else:
            plotAntennas(telescope, names, ids, xpos, ypos, antindex)
        if not title:
            msname = os.path.basename(vis)
            title = "Antenna Positions for "
            if len(msname) > 55:
                title += '\n'
            title += msname
        pl.title(title, {'fontsize':12})
        if figfile:
            pl.savefig(figfile)
    except Exception as instance:
        casalog.post("Antenna plot failed: " + str(instance), "ERROR")

def getAntennaInfo(msname, exclude, checkbaselines):

    tb, me = gentools(['tb', 'me'])
    qa = qatool()

    telescope, arrayPos = getObservatoryInfo(msname)
    arrayWgs84 = me.measure(arrayPos, 'WGS84')
    arrayLon, arrayLat, arrayAlt = [arrayWgs84[i]['value'] 
        for i in ['m0','m1','m2']]

    # Open the ANTENNA subtable to get the names of the antennas in this MS and
    # their positions.  Note that the entries in the ANTENNA subtable are pretty
    # much in random order, so antNames translates between their index and name
    # (e.g., index 11 = STD155).  We'll need these indices for later, since the
    # main data table refers to the antennas by their indices, not names.

    anttabname = msname + '/ANTENNA'
    tb.open(anttabname)
    # Get antenna names from antenna table
    colname = 'NAME'
    if telescope == 'VLBA':
        colname = 'STATION'
    antNames = np.array(tb.getcol(colname)).tolist()
    allAntIds = range(len(antNames))
    # Get antenna positions from antenna table
    antPositions = np.array([me.position('ITRF', qa.quantity(x, 'm'),
        qa.quantity(y, 'm'), qa.quantity(z, 'm'))
        for (x, y, z) in tb.getcol('POSITION').transpose()])
    tb.close()

    if checkbaselines:
        # Get antenna ids from main table; this will add to runtime
        tb.open(msname)
        ants1 = tb.getcol('ANTENNA1')
        ants2 = tb.getcol('ANTENNA2')
        tb.close()
        antIdsUsed = set(np.append(ants1, ants2))
    else:
        # use them all!
        antIdsUsed = allAntIds

    # handle exclude -- remove from antIdsUsed
    for ant in exclude:
        if isinstance(ant, int):
            if ant not in allAntIds:
                casalog.post("Cannot exclude antenna id " + str(ant) + 
                        ": does not exist", "WARN")
            else:
                try:
                    antIdsUsed.remove(ant)
                except KeyError: # id exists but not used anyway
                    casalog.post("Cannot exclude antenna " + str(ant) + 
                        ": not in main table", "WARN")
        if isinstance(ant, str):
            if ant:  # default is empty string
                start = 0
                count = antNames.count(ant)
                if count is 0:
                    casalog.post("Cannot exclude antenna " + ant + 
                            ": does not exist", "WARN")
                for i in range(count):
                    try:
                        idx = antNames[start:].index(ant)
                        antIdsUsed.remove(start + idx)
                        start += idx + 1
                    except ValueError: # id not in list of used ids
                        casalog.post("Cannot exclude antenna " + ant + 
                            ": not in main table", "WARN")

    # apply antIdsUsed mask
    ids = list(antIdsUsed)
    antNames = [antNames[i] for i in ids]
    antPositions = [antPositions[i] for i in ids]

    nAnts = len(ids)
    print "Number of points being plotted:", nAnts
    casalog.post("Number of points being plotted: " + str(nAnts))
    if nAnts == 0: # excluded all antennas
        return telescope, antNames, [], [], []
    
    # Get the names, indices, and lat/lon/alt coords of "good" antennas.
    antWgs84s = np.array([me.measure(pos, 'WGS84') for pos in antPositions])
    antLons, antLats, antAlts = [np.array( [pos[i]['value'] 
        for pos in antWgs84s]) for i in ['m0','m1','m2']]
    
    # Convert from lat, lon, alt to X, Y, Z (unless VLBA)
    # where X is east, Y is north, Z is up,
    # and 0, 0, 0 is the center of the LWA1.  
    # Note: this conversion is NOT exact, since it doesn't take into account
    # Earth's ellipticity!  But it's close enough.
    radE = 6370000.
    antXs = (antLons - arrayLon) * radE * np.cos(arrayLat)
    antYs = (antLats - arrayLat) * radE
    antZs = antAlts - arrayAlt
    
    return telescope, antNames, ids, antXs, antYs

def getObservatoryInfo(msname):
    metadata = msmdtool()
    metadata.open(msname)
    telescope = metadata.observatorynames()[0]
    arrayPos = metadata.observatoryposition()
    metadata.close()
    return telescope, arrayPos

def plotAntennasLog(telescope, names, ids, xpos, ypos, antindex):
    fig = pl.figure(1)
    # Set up subplot.
    ax = fig.add_subplot(1, 1, 1, polar=True, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    # Do not show azimuth labels.
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Do not show grid.
    ax.grid(False)

    # code from pipeline summary.py 
    # PlotAntsChart draw_polarlog_ant_map_in_subplot
    if telescope in ('VLA', 'EVLA'):
        # For (E)VLA, set a fixed local center position that has been
        # tuned to work well for its array configurations (CAS-7479).
        xcenter, ycenter = -32, 0
        rmin_min, rmin_max = 12.5, 350
    else:
        # For non-(E)VLA, take the median of antenna offsets as the
        # center for the plot.
        xcenter = np.median(xpos)
        ycenter = np.median(ypos)
        rmin_min, rmin_max = 3, 350

    # Derive radial offset w.r.t. center position.
    r = ((xpos-xcenter)**2 + (ypos-ycenter)**2)**0.5
    # Set rmin, clamp between a min and max value, ignore station
    # at r=0 if one is there.
    rmin = min(rmin_max, max(rmin_min, 0.8*np.min(r[r > 0])))
    # Update r to move any points below rmin to r=rmin.
    r[r <= rmin] = rmin
    rmin = np.log(rmin)
    # Set rmax.
    rmax = np.log(1.5*np.max(r))
    # Derive angle of offset w.r.t. center position.
    theta = np.arctan2(xpos-xcenter, ypos-ycenter)

    # Draw circles at specific distances from the center.
    angles = np.arange(0, 2.01*np.pi, 0.01*np.pi)
    show_circle = True
    circles = [30, 100, 300, 1000, 3000, 10000]
    if telescope == "VLBA":
        circles = [1e5, 3e5, 1e6, 3e6, 1e7]
    for cr in circles:

        # Only draw circles outside rmin.
        if cr > np.min(r) and show_circle:

            # Draw the circle.
            radius = np.ones(len(angles))*np.log(cr)
            ax.plot(angles, radius, 'k:')

            # Draw tick marks on the circle at 1 km intervals.
            inc = 0.1*10000/cr
            if telescope == "VLBA":
                inc = 0.1*1e8/cr
            if cr > 100:
                for angle in np.arange(inc/2., 2*np.pi+0.05, inc):
                    ax.plot([angle, angle],
                               [np.log(0.95*cr), np.log(1.05*cr)], 'k-')

            # Add text label to circle to denote distance from center.
            va = 'top'
            circle_label_angle = -20.0 * np.pi / 180.
            if cr >= 1000:
                if np.log(cr) < rmax:
                    ax.text(circle_label_angle, np.log(cr),
                               '%d km' % (cr/1000), size=8, va=va)
                    ax.text(circle_label_angle + np.pi, np.log(cr),
                               '%d km' % (cr / 1000), size=8, va=va)
            else:
                ax.text(circle_label_angle, np.log(cr), '%dm' % (cr), 
                    size=8, va=va)
                ax.text(circle_label_angle + np.pi, np.log(cr), 
                        '%dm' % (cr), size=8, va=va)

        # Find out if most recently drawn circle was outside all antennas, 
        # if so, no more circles will be drawn.
        if np.log(cr) > rmax:
            show_circle = False

    # plot points and antenna names/ids
    for i, antenna in enumerate(names):
        ax.plot(theta[i], np.log(r[i]), 'ko', ms=5, mfc='r')
        if antindex:
            antenna += ' (' + str(ids[i]) + ')'
        ax.text(theta[i], np.log(r[i]), '   '+antenna, size=8, va='top')

    # Set minimum and maximum radius.
    ax.set_rmax(rmax)
    ax.set_rmin(rmin)

    # Make room for 2-line title
    pl.subplots_adjust(top=0.88)

def plotAntennas(telescope, names, ids, xpos, ypos, antindex):
    fig = pl.figure(1)
    ax = fig.add_subplot(111)

    # use m or km units
    units = ' (m)'
    if np.median(xpos) > 1e6 or np.median(ypos) > 1e6:
        xpos /= 1e3
        ypos /= 1e3
        units = ' (km)'

    # plot points and antenna names/ids
    for i, (x, y, name) in enumerate(zip(xpos, ypos, names)):
        ax.plot(x, y, 'ro')
        if antindex:
            name += ' (' + str(ids[i]) + ')'
        ax.text(x, y+2, '  ' + name, size=8, va='top')
        fig.show()

    pl.xlabel('X' + units)
    pl.ylabel('Y' + units)
    pl.margins(0.1, 0.1)
