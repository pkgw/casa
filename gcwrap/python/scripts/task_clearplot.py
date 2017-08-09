import os
from taskinit import *
import pylab as pl

def clearplot():
	"""Clear the matplotlib plotter and all layers:

	"""
	try:
		pl.ion()
		print("Calling tp.clearplot()")
		ok=tp.clearplot(0,0,0)
	except Exception as instance:
		print('*** Error ***',instance)
