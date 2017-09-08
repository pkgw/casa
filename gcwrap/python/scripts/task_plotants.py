import os
import numpy as np
import pylab as pl
from taskinit import gentools, qatool

def plotants(vis=None, figfile=''):
    """Plot the antenna distribution in the local reference frame:

    The location of the antennas in the MS will be plotted with
    X-toward local east; Y-toward local north.

    Keyword arguments:
    vis -- Name of input visibility file.
            default: none. example: vis='ngc5921.ms'

    figfile -- Save the plotted figure in this file.
            default: ''. example: figfile='myFigure.png'

            The name of each antenna (egs. vla=antenna number) is
               shown next to its respective location.

            DO NOT use the buttons on the Mark Region line.  These are
               not implemented yet and might abort CASA.

            You can zoom in by pressing the magnifier button (bottom,
               third from left) and making a rectangular region with
               the mouse.  Press the home button (left most button) to
               remove zoom.

            A hard-copy of this plot can be obtained by pressing the
               button on the right at the bottom of the display.  This
               produces a png format file.
    """

    try:
        antNames, antXs, antYs = getAntennaInfo(vis)
        pl.clf()
        plotAntennas(antXs, antYs, antNames)
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')
        pl.margins(0.1, 0.1)
        if figfile:
            pl.savefig(figfile)
    except Exception, instance:
          print '*** Error ***',instance

def getAntennaInfo(msname):

    msmd, me, tb = gentools(['msmd', 'me', 'tb'])
    qa = qatool()

    msmd.open(msname)
    arrPos = msmd.observatoryposition()
    msmd.close()
    arrWgs84 = me.measure(arrPos, 'WGS84')
    arrLon, arrLat, arrAlt = [arrWgs84[i]['value'] for i in ['m0','m1','m2']]

    # Open the ANTENNA subtable to get the names of the antennas in this MS and
    # their positions.  Note that the entries in the ANTENNA subtable are pretty
    # much in random order, so allAnts translates between their index and name
    # (e.g., index 11 = STD155).  We'll need these indices for later, since the
    # main data table refers to the antennas by their indices, not names.
    anttabname = msname + '/ANTENNA'
    tb.open(anttabname)
    antNames = np.array(tb.getcol('NAME'))
    antPoss = np.array([me.position('ITRF', qa.quantity(x, 'm'),
                                    qa.quantity(y, 'm'), qa.quantity(z, 'm'))
                        for (x, y, z) in tb.getcol('POSITION').transpose()])
    tb.close()
    
    # Get the names, indices, and lat/lon/alt coords of our "good" antennas.
    antWgs84s = np.array([me.measure(pos, 'WGS84') for pos in antPoss])
    antLons, antLats, antAlts = [np.array( [pos[i]['value'] 
        for pos in antWgs84s]) for i in ['m0','m1','m2']]
    
    # Convert from lat, lon, alt to X, Y, Z, 
    # where X is east, Y is north, Z is up,
    # and 0, 0, 0 is the center of the LWA1.  
    # Note: this conversion is NOT exact, since it doesn't take into account
    # Earth's ellipticity!  But it's close enough.
    radE = 6370000.
    antXs = (antLons - arrLon) * radE * np.cos(arrLat)
    antYs = (antLats - arrLat) * radE
    antZs = antAlts - arrAlt
    
    return antNames, antXs, antYs

def plotAntennas(antXs, antYs, antNames):
    fig = pl.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(antXs, antYs, 'ro')
    for (x, y, name) in zip(antXs, antYs, antNames):
        ax.text(x + 0.5, y + 0.5, name)
        fig.show()
