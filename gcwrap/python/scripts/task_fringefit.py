import os
import numpy as np
from callibrary import *
from taskinit import *
# For stack frames in debugging
import sys
def fringefit(vis=None,caltable=None,
	      field=None,spw=None,intent=None,
	      selectdata=None,timerange=None,antenna=None,scan=None,
	      observation=None, msselect=None,
	      solint=None,combine=None,refant=None,
	      minsnr=None,zerorates=None,globalsolve=None,niter=None,
              delaywindow=None,ratewindow=None,append=None,
	      docallib=None, callib=None, gaintable=None,gainfield=None,interp=None,spwmap=None,
	      parang=None):
	#Python script
	casalog.origin('fringefit')
	try: 
		mycb = cbtool()
		oldness = False
		casalog.post("Setting vi to {}".format(oldness))
		mycb.setvi(oldness) # Use new calibration framework
		if ((type(vis)==str) & (os.path.exists(vis))):
			mycb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
		else:
			raise Exception('Visibility data set not found - please verify the name')
		# Do data selection according to selectdata
		if (selectdata):
			casalog.post("Selecting data")
			# pass all data selection parameters in as specified
			mycb.selectvis(time=timerange,spw=spw, scan=scan, field=field,
				     intent=intent, observation=str(observation),
				     baseline=antenna, chanmode='none',
				     msselect=msselect)
		else:
			# selectdata=F, so time,scan,baseline,msselect=''
			# using spw and field specifications only
			mycb.selectvis(time='',spw=spw,scan='',field=field,
				       observation='', baseline='', 
				       chanmode='none', msselect='')
			
		if docallib:
			# by cal library from file
			mycallib=callibrary()
			mycallib.read(callib)
			mycb.setcallib(mycallib.cld)
		else:

			# by traditional parameters
			ngaintab = 0;
			if (gaintable!=['']):
				ngaintab=len(gaintable)

			ngainfld = len(gainfield)
			nspwmap = len(spwmap)
			ninterp = len(interp)

			# handle list of list issues with spwmap
			if (nspwmap>0):
				if (type(spwmap[0])!=list):
					# first element not a list, only one spwmap specified
					# make it a list of list
					spwmap=[spwmap];
					nspwmap=1;

			for igt in range(ngaintab):
				if (gaintable[igt]!=''):

					# field selection is null unless specified
					thisgainfield=''
					if (igt<ngainfld):
						thisgainfield=gainfield[igt]

					# spwmap is null unless specifed
					thisspwmap=[-1]
					if (igt<nspwmap):
						thisspwmap=spwmap[igt];

					# interp is 'linear' unless specified
					thisinterp='linear'
					if (igt<ninterp):
						if (interp[igt]==''):
							interp[igt]=thisinterp;
						thisinterp=interp[igt];

					mycb.setapply(t=0.0,table=gaintable[igt],field=thisgainfield,
						      calwt=True,spwmap=thisspwmap,interp=thisinterp)
                        if len(delaywindow) != 2:
                                delaywindow = [-1e6, 1e6]
                        if len(ratewindow) != 2:
                                ratewindow = [-1e6, 1e6]
		# ...and now the specialized terms
		# (BTW, interp irrelevant for these, since they are evaluated)
		
		# Apply parallactic angle, if requested
		if parang: mycb.setapply(type='P')

		# Set up for solving; only support one gaintype
		mycb.setsolve(type="FRINGE",t=solint,refant=refant,preavg=0.01,
			      minsnr=minsnr,combine=combine,
			      zerorates=zerorates,
                              globalsolve=globalsolve,
                              niter=niter,
                              delaywindow=delaywindow,
                              ratewindow=ratewindow,
                	      table=caltable,append=append)
		mycb.solve()
		reportsolvestats(mycb.activityrec());
		mycb.close()

	except Exception as instance:
		print('*** Error ***', instance)
		mycb.close()
		exc_type, exc_obj, exc_tb = sys.exc_info()
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		casalog.post("Error in fringefit: %s" % str(instance), "SEVERE")
		casalog.post("{}:{}:{}".format(exc_type, fname, exc_tb.tb_lineno))
		raise Exception("Error in fringefit: "+str(instance))

def reportsolvestats(rec):
	if (list(rec.keys()).count('origin')==1 and
	    rec['origin']=='Calibrater::genericGatherAndSolve'):
		casalog.post("Calibration solve statistics per spw:  (expected/attempted/succeeded):")
		nexp=rec['nExpected']
		natt=rec['nAttempt']
		nsuc=rec['nSucceed']
		for ispw in range(len(nexp)):
			solstatstr="  Spw "+str(ispw)+": "
			solstatstr+=str(nexp[ispw])+"/"
			solstatstr+=str(natt[ispw])+"/"
			solstatstr+=str(nsuc[ispw])
			casalog.post(solstatstr)
