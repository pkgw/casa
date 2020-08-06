import os
import numpy as np
import pylab as pl
from textwrap import wrap
from taskinit import msmdtool, gentools, qatool, casalog, mstool
from casa_system import casa

def plotants(vis=None, figfile=None, 
        antindex=None, logpos=None, 
        exclude=None, checkbaselines=None,
        title=None, showgui=None):
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
	showplot = showgui and not casa['flags'].nogui and not casa['flags'].agg and \
		not casa['flags'].pipeline  # pipeline sets backend to 'agg'
	if not showplot:
		pl.close()
		pl.ioff()
	else:
		pl.show()
		pl.ion()
	pl.clf()

	try:
		# remove trailing / for title basename
		if vis[-1]=='/':
			vis = vis[:-1]
		myms = mstool()
		try:
			exclude = myms.msseltoindex(vis, baseline=exclude)['antenna1'].tolist()
		except RuntimeError as rterr:  # MSSelection failed
			errmsg = str(rterr)
			errmsg = errmsg.replace('specificion', 'specification')
			errmsg = errmsg.replace('Antenna Expression: ', '')
			casalog.post("Exclude selection error: " + errmsg, "ERROR")
			return

		telescope, names, ids, xpos, ypos, stations = getPlotantsAntennaInfo(vis,
			logpos, exclude, checkbaselines)
		if not names:
			casalog.post("No antennas selected.  Exiting plotants.", "ERROR")
			return

		if not title:
			msname = os.path.basename(vis)
			title = "Antenna Positions for "
			if len(msname) > 55:
				title += '\n'
			title += msname
		pl.title(title, {'fontsize':12})

		if logpos:
			plotAntennasLog(telescope, names, ids, xpos, ypos, antindex, stations)
		else:
			plotAntennas(telescope, names, ids, xpos, ypos, antindex, stations, showplot)
		if figfile:
			pl.savefig(figfile)

	except Exception as instance:
		casalog.post("Error: " + str(instance), "ERROR")

def getPlotantsAntennaInfo(msname, log, exclude, checkbaselines):

	tb, me = gentools(['tb', 'me'])
	qa = qatool()

	telescope, arrayPos = getPlotantsObservatoryInfo(msname)
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
	antNames = np.array(tb.getcol("NAME")).tolist()
	stationNames = np.array(tb.getcol("STATION")).tolist()
	if telescope == 'VLBA':  # names = ant@station
		antNames = ['@'.join(antsta) for antsta in zip(antNames,stationNames)]
	# Get antenna positions from antenna table
	antPositions = np.array([me.position('ITRF', qa.quantity(x, 'm'),
		qa.quantity(y, 'm'), qa.quantity(z, 'm'))
		for (x, y, z) in tb.getcol('POSITION').transpose()])
	tb.close()

	allAntIds = list(range(len(antNames)))
	if checkbaselines:
		# Get antenna ids from main table; this will add to runtime
		tb.open(msname)
		ants1 = tb.getcol('ANTENNA1')
		ants2 = tb.getcol('ANTENNA2')
		tb.close()
		antIdsUsed = list(set(np.append(ants1, ants2)))
	else:
		# use them all!
		antIdsUsed = allAntIds

	# handle exclude -- remove from antIdsUsed
	for antId in exclude:
		try:
			antNameId = antNames[antId] + " (id " + str(antId) + ")"
			antIdsUsed.remove(antId)
			casalog.post("Exclude antenna " + antNameId)
		except ValueError:
			casalog.post("Cannot exclude antenna " + antNameId + ": not in main table", "WARN")

	# apply antIdsUsed mask
	antNames = [antNames[i] for i in antIdsUsed]
	antPositions = [antPositions[i] for i in antIdsUsed]
	stationNames = [stationNames[i] for i in antIdsUsed]

	nAnts = len(antIdsUsed)
	print("Number of points being plotted:", nAnts)
	casalog.post("Number of points being plotted: " + str(nAnts))
	if nAnts == 0: # excluded all antennas
		return telescope, antNames, [], [], []

	# Get the names, indices, and lat/lon/alt coords of "good" antennas.
	antWgs84s = np.array([me.measure(pos, 'WGS84') for pos in antPositions])

	# Convert from lat, lon, alt to X, Y, Z (unless VLBA)
	# where X is east, Y is north, Z is up,
	# and 0, 0, 0 is the center
	# Note: this conversion is NOT exact, since it doesn't take into account
	# Earth's ellipticity!  But it's close enough.
	if telescope == 'VLBA' and not log:
		antLons, antLats = [[pos[i] for pos in antWgs84s] for i in ['m0','m1']]
		antXs = [qa.convert(lon, 'deg')['value'] for lon in antLons]
		antYs = [qa.convert(lat, 'deg')['value'] for lat in antLats]
	else:
		antLons, antLats = [np.array( [pos[i]['value']
			for pos in antWgs84s]) for i in ['m0','m1']]
		radE = 6370000.
		antXs = (antLons - arrayLon) * radE * np.cos(arrayLat)
		antYs = (antLats - arrayLat) * radE
	return telescope, antNames, antIdsUsed, antXs, antYs, stationNames

def getPlotantsObservatoryInfo(msname):
	metadata = msmdtool()
	metadata.open(msname)
	telescope = metadata.observatorynames()[0]
	arrayPos = metadata.observatoryposition()
	metadata.close()
	return telescope, arrayPos

def getAntennaLabelProps(telescope, station, log=False):
	# CAS-7120 make plots more readable (rotate labels)
	vAlign = 'center'
	hAlign = 'left'
	rotAngle = 0
	if station and "VLA" in telescope:
		# these have non-standard format:
		# strip off VLA: or VLA:_ prefix if any
		if 'VLA:' in station:
			station = station[4:]
		if station[0] == '_':
			station = station[1:]
		if station in ['W01', 'E01', 'W1', 'E1']:
			vAlign = 'top'
			hAlign = 'center'
		else:
			vAlign = 'bottom'
			if 'W' in station or 'MAS' in station:
				hAlign = 'right'
				rotAngle = -35
			elif 'E' in station or ('N' in station and not log):
				rotAngle = 35
	return vAlign, hAlign, rotAngle

def plotAntennasLog(telescope, names, ids, xpos, ypos, antindex, stations):
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
	if 'VLA' in telescope:
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
	for i, (name, station) in enumerate(zip(names, stations)):
		if station and 'OUT' not in station:
			ax.plot(theta[i], np.log(r[i]), 'ko', ms=5, mfc='r')
			if antindex:
				name += ' (' + str(ids[i]) + ')'
			# set alignment and rotation angle (for VLA)
			valign, halign, angle = getAntennaLabelProps(telescope, station, log=True)
			# adjust so text is not on the circle:
			yoffset = 0
			if halign is 'center':
				yoffset = 0.1
			ax.text(theta[i], np.log(r[i])+yoffset, ' '+name, size=8, va=valign, ha=halign,
				rotation=angle, weight='semibold')

	# Set minimum and maximum radius.
	ax.set_rmax(rmax)
	ax.set_rmin(rmin)

	# Make room for 2-line title
	pl.subplots_adjust(top=0.88)

def plotAntennas(telescope, names, ids, xpos, ypos, antindex, stations, showplot):
	fig = pl.figure(1)
	ax = fig.add_subplot(111)

	if telescope == 'VLBA':
		labelx = 'Longitude (deg)'
		labely = 'Latitude (deg)'
	else:
		# use m or km units
		units = ' (m)'
		if np.median(xpos) > 1e6 or np.median(ypos) > 1e6:
			xpos /= 1e3
			ypos /= 1e3
			units = ' (km)'
		labelx = 'X' + units
		labely = 'Y' + units
	if "VLA" in telescope:
		spacer = ' '
	else:
		spacer = '  '

	# plot points and antenna names/ids
	for i, (x, y, name, station) in enumerate(zip(xpos, ypos, names, stations)):
		if station and 'OUT' not in station:
			ax.plot(x, y, 'ro')
			if antindex:
				name += ' (' + str(ids[i]) + ')'
			# set alignment and rotation angle (for VLA)
			valign, halign, angle = getAntennaLabelProps(telescope, station)
			# adjust so text is not on the circle:
			if halign is 'center':
				y -= 10
			ax.text(x, y, ' '+name, size=8, va=valign, ha=halign, rotation=angle,
				weight='semibold')
			if showplot:
				fig.show()

	pl.xlabel(labelx)
	pl.ylabel(labely)
	pl.margins(0.1, 0.1)

