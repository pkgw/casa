import os
import numpy as np
from taskinit import *

from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

def apparentsens(vis=None,
		 field=None,spw=None,intent=None,
		 selectdata=None,timerange=None,uvrange=None,antenna=None,scan=None,
		 observation=None, msselect=None,
		 imsize=None,cell=None,
		 stokes=None,specmode=None,
		 weighting=None,
		 robust=None,npixels=None,uvtaper=None):

	#Python script
        casalog.origin('apparentsens')

	try: 

		paramList = ImagerParameters(
			msname=vis,
			field=field,
			scan=scan,
			spw=spw,
			imagename=vis+'.apparentsens',
			imsize=imsize,
			cell=cell,
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
		imager.makePSF()
		
		imager.deleteTools()

	except Exception, instance:
		print '*** Error ***', instance
		casalog.post("Error in apparentsens: %s" % str(instance), "SEVERE")
		raise Exception, "Error in apparentsens: "+str(instance)

