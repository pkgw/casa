import os
import numpy as np
from callibrary import *
from taskinit import *

def accor(vis=None,caltable=None,
	  field=None,spw=None,intent=None,
	  selectdata=None,timerange=None,antenna=None,scan=None,
	  observation=None, msselect=None,
	  solint=None,combine=None,append=None,
	  docallib=None,callib=None,
	  gaintable=None,gainfield=None,interp=None,spwmap=None):

	#Python script
        casalog.origin('accor')

	try: 
                mycb = cbtool()

                if ((type(vis)==str) & (os.path.exists(vis))):
                        mycb.open(filename=vis,compress=False,addcorr=False,addmodel=False)
                else:
                        raise Exception('Visibility data set not found - please verify the name')

		# Do data selection according to selectdata
		casalog.post("NB: accor automatically excludes crosso-correlations.")
		if (selectdata):
			# insist only CCs
			if len(msselect)>0:
				msselect='('+msselect+') && ANTENNA==ANTENNA2'
			else:
				msselect='ANTENNA1==ANTENNA2'

			# pass all data selection parameters in as specified
			mycb.selectvis(time=timerange,spw=spw, scan=scan, field=field,
				     intent=intent, observation=str(observation),
				     baseline=antenna,uvrange='',chanmode='none',
				     msselect=msselect);
		else:
			# selectdata=F, so time,scan,baseline,uvrange,msselect=''
			# using spw and field specifications only
			# also insist only CCs
			mycb.selectvis(time='',spw=spw,scan='',field=field,intent=intent,
                                     observation='', baseline='', uvrange='',
                                     chanmode='none', msselect='ANTENNA1==ANTENNA2')

		# Arrange applies....

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

		# Set up for solving:  
		mycb.setsolve(type='ACCOR',t=solint,combine=combine,
			      table=caltable, append=append)
		mycb.solve()
		mycb.close()

	except Exception as instance:
		print('*** Error ***', instance)
		mycb.close()
		casalog.post("Error in accor: %s" % str(instance), "SEVERE")
		raise Exception("Error in accor: "+str(instance))
