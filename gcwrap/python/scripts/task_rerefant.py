import os
from taskinit import *

_cb = cbtool( )

def rerefant(vis,tablein,caltable,refantmode,refant):
	""" Smooth calibration solution(s) derived from one or more sources:

	Keyword arguments:
	vis -- Name of input visibility file (MS)
		default: none; example: vis='ngc5921.ms'
	tablein -- Input calibration table (any type)
		default: none; example: tablein='ngc5921.gcal'
	caltable -- Output calibration table (re-refant-ed)
		default: none
	refantmode -- Reference antenna application mode (TBD)
	refant -- Reference antenna name(s)
	"""

	#Python script
	try:
		casalog.origin('smoothcal')
                if ((type(vis)==str) & (os.path.exists(vis))):
                        _cb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
			       
                else:
                        raise Exception('Visibility data set not found - please verify the name')

		_cb.rerefant(tablein=tablein,tableout=caltable,refantmode=refantmode,refant=refant);
		_cb.close()

	except Exception as instance:
		print('*** Error ***',instance)
		_cb.close()
		raise Exception(instance)


