import os
import numpy as np
from taskinit import *

from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

def apparentsens(vis=None,
		 field=None,spw=None,intent=None,
		 selectdata=None,timerange=None,uvrange=None,antenna=None,scan=None,
		 observation=None, #  msselect=None,
		 imsize=None,cell=None,
		 stokes=None,specmode=None,
		 weighting=None,
		 robust=None,npixels=None,uvtaper=None):

	#Python script
        casalog.origin('apparentsens')

	try: 

		# insist on continuum calculation only (for now)
		if specmode!='mfs':
			casalog.post( "The apparentsens task only supports specmoe='mfs' at this time", "WARN", "task_apparentsens" )

		imname=vis+'.apparentsens'

		# remove the image files that imager makes below
		os.system('rm -Rf '+imname+'.*')

		# fill the relevant parameters
		paramList = ImagerParameters(
			msname=vis,
			field=field,
			spw=spw,
			timestr=timerange,
			uvdist=uvrange,
			antenna=antenna,
			scan=scan,
			obs=observation,
			state=intent,
			imagename=imname,
			imsize=imsize,
			cell=cell,
			stokes=stokes,
			specmode=specmode,
			nchan=-1,
			weighting=weighting,
			robust=robust,
			npixels=npixels,
			uvtaper=uvtaper
			)

		imager = PySynthesisImager(params=paramList)

		imager.initializeImagers()
		imager.initializeNormalizers()
		imager.setWeighting()

		out = imager.calcVisAppSens()

#		imager.makePSF()
		
		imager.deleteTools()

		return out

	except Exception as instance:
		print('*** Error ***', instance)
		casalog.post("Error in apparentsens: %s" % str(instance), "SEVERE")
		raise Exception("Error in apparentsens: "+str(instance))

