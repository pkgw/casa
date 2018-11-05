import os
import time
import datetime
import regression_utility as regutl
import sys, traceback
import numpy as np
import pipeline

pathname = os.environ.get('CASAPATH').split()[0]
rootdatapath = pathname+'/data/regression/pipeline/'


'''Initial VLA pipeline regression
   B. Kent, May 2015
   Update September 1, 2015
   Update April 20, 2018
   Update June  01, 2018
'''

THISHOME  = "working/"
startTime=0.0
endTime=0.0
startProc=0.0
endProc=0.0
regstate = True
standard_file = rootdatapath + 'VLApipeline44-standard'
# MIN_CASA_REVISION = 36095



def load_context(filename):
    with open(filename, 'rb') as picklefile:
        return pipeline.infrastructure.utils.pickle_load(picklefile)

#
EPS       = 1e-5  # Logical "zero"
#
#--------------------------------------------------------------
#
def pipeline_regression():
    global regstate
    # global MIN_CASA_REVISION
        
    
    #revision = int(casadef.subversion_revision)
    #if MIN_CASA_REVISION > revision:
    #    msg = ('Minimum CASA revision for the pipeline is r%s, '
    #           'got CASA %s (r%s).' % (MIN_CASA_REVISION,
    #            cu.version_info( ),
    #            casadef.subversion_revision))
    #    print msg
    #    regstate = False
    #    raise EnvironmentError(msg)
    
    
    
    # ASDM      = "/export/home/icarus_2/awells/CASA_stable/data/regression/foo/vla_pipeline_data/rawdata/13A-537.sb24066356.eb24324502.56514.05971091435"
    ASDM = rootdatapath  + "13A-537.sb24066356.eb24324502.56514.05971091435"
    try:
        import pipeline.recipes.hifv as hifv
    except ImportError, e:
        print(e)
        
    
    # Check to see if the ASDM exists
    if not os.path.exists(ASDM):
        print("Unable to open ASDM ", ASDM)
	regstate=False
        raise IOError
    else:
        print("Using ", ASDM)
    
    # Run the CASA VLA Pipeline standard recipe
    try:
        hifv.hifv([ASDM], importonly=False)
    except Exception as e:
        print(e)
        regstage=False



def run():
    '''Timing and running of main regression function
    '''
    global startTime, endTime, startProc, endProc, regstate;
    startTime = time.time()
    startProc = time.clock()
    regstate = False
    pipeline_regression()
    endTime = time.time()
    endProc = time.clock()
    print("Run Time = ",endTime-startTime,endProc-startProc)

def stats():
    global startTime, endTime, startProc, endProc, regstate, standard_context_file
    
    datestring = datetime.datetime.isoformat(datetime.datetime.today())
    outfile = 'vlapipeline.'+datestring+'.log'
    logfile = open(outfile,'w')
    
    try:
        import pipeline
    except ImportError, e:
        print e
        print >>logfile, "Unable to import the CASA pipeline"
        regstate = False
    
    try:
        # Open context
        context = pipeline.Pipeline(context='last').context
        regstate = True
        print >>logfile,"VLA pipeline context created."
        print("VLA pipeline context created.")
        print >>logfile,"Context verification - VLA pipeline regression PASSED."
        print("Context verification - VLA pipeline regression PASSED.")
        
    except Exception as e:
        regstate = False
        print >>logfile,"VLA pipeline context NOT created."
        print("VLA pipeline context NOT created.")
        print >>logfile,"Context verification - VLA pipeline regression FAILED."
        print("Context verification - VLA pipeline regression FAILED.")
        print >>logfile, e
        print e

    # Test fluxscale values
    # Test that hifv_fluxboot stage was the 13th stage run (index 12)
    try:
        # Open context of regression pipeline run
        context = pipeline.Pipeline(context='last').context
        fluxlist = context.results[12].read()[0].flux_densities
        rtol = 1.0e-5  #Relative Tolerance
        atol = 1.0e-5  #Absolute Tolerance

        # value_compare43 = 0.6934577681171487
        # value_compare = 0.716367318068  # CASA 4.5
        # value_compare = 0.717857716108  # CASA 4.6
        # value_compare = 0.716364780148  # CASA 4.6.144 r36095, pipeline r36209 (trunk)
        # value_compare = 0.718457023749  # CASA-test 5.0.11-DEV (r38177) 
        # value_compare = 0.71857779577   # CASA-prerelease 5.0.0-187, pipeline r40156
        # value_compare = 0.718551806637  # CASA-prerelease 5.1.0-34, pipeline r40615
        # value_compare = 0.718519030402  # CASA-prerelease 5.3.0-26, pipeline r40909
        # value_compare =  0.718342661383 # CASA-prerelease 5.3.0-114, pipeline r41386
        value_compare =  0.717434107414   # CASA-prerelease 5.4.0-3, pipeline r41527
        
        result_bool = np.isclose(fluxlist[0][0], value_compare, rtol=rtol, atol=atol, equal_nan=False)
        
        print >>logfile, "Accepted test value is: ", value_compare, " from CASA-prerelease 5.4.0-3, pipeline r41527 (trunk)"
        print "Accepted test value is: ", value_compare, " from CASA-prerelease 5.4.0-3, pipeline r41527 (trunk)"
        print >>logfile, "Regression generated value is: ", fluxlist[0][0]
        print "Regression generated value is: ", fluxlist[0][0]
        
        if (result_bool):
            regstate=True
            print >>logfile,"hifv_fluxboot values match within relative tolerance of 1.0e-05 and absolute tolerance of 1.0e-05"
            print "hifv_fluxboot values match within relative tolerance of 1.0e-05 and absolute tolerance of 1.0e-05"
            print >>logfile,"hifv_fluxboot test PASSED."
            print "hifv_fluxboot test PASSED."
        else:
            regstate=False
            print >>logfile,"hifv_fluxboot values are not within tolerances."
            print "hifv_fluxboot values are not within tolerances."
            print >>logfile,"hifv_fluxboot test FAILED."
            print "hifv_fluxboot test FAILED."
        
    except Exception as e:
        regstate=False
        print >>logfile,"hifv_fluxboot values are not within tolerances."
        print "hifv_fluxboot values are not within tolerances."
        print >>logfile,"hifv_fluxboot test FAILED."
        print "hifv_fluxboot test FAILED."
        print >>logfile, e
        print e

    logfile.close()

def main():
    try:
        run()
        stats()
    except KeyboardInterrupt:
        print("Interrupt requested...exiting")
    except Exception:
        traceback.print_exc(file=sys.stdout)

if __name__ == "__main__":
    main()
    print "Regstate:" , regstate
    if regstate:
        print("Regression PASSED")
    else:
        print("Regression FAILED")
