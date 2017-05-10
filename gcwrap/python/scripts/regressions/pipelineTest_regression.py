import os
import time
import datetime
import regression_utility as regutl
import sys, traceback
import numpy as np
import pipeline



pathname=os.environ.get('CASAPATH').split()[0]
rootdatapath = pathname+'/data/regression/pipeline/'


'''Initial VLA pipeline regression
   B. Kent, May 2015
   Last update September 1, 2015
'''

THISHOME  = "working/"
startTime=0.0
endTime=0.0
startProc=0.0
endProc=0.0
regstate = True
standard_file = rootdatapath + 'VLApipeline44-standard'
#MIN_CASA_REVISION = 36095



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
    #global MIN_CASA_REVISION
        
    
    #revision = int(casadef.subversion_revision)
    #if MIN_CASA_REVISION > revision:
    #    msg = ('Minimum CASA revision for the pipeline is r%s, '
    #           'got CASA %s (r%s).' % (MIN_CASA_REVISION,
    #            cu.version_info( ),
    #            casadef.subversion_revision))
    #    print msg
    #    regstate = False
    #    raise EnvironmentError(msg)
    
    
    
    #ASDM      = "/export/home/icarus_2/awells/CASA_stable/data/regression/foo/vla_pipeline_data/rawdata/13A-537.sb24066356.eb24324502.56514.05971091435"
    ASDM = rootdatapath  + "13A-537.sb24066356.eb24324502.56514.05971091435"
    try:
        import pipeline.recipes.hifv as hifv
    except ImportError, e:
        print e
        
    
    #Check to see if the ASDM exists
    if not os.path.exists(ASDM):
        print "Unable to open ASDM " + ASDM
	regstate=False
        raise IOError
    else:
        print "Using " + ASDM
    
    #Run the CASA VLA Pipeline standard recipe
    try:
        hifv.hifv([ASDM], importonly=False)
    except Exception, e:
        print e
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
    print "Run Time = ",endTime-startTime,endProc-startProc

def stats():
    global startTime, endTime, startProc, endProc, regstate, standard_context_file
    
    datestring=datetime.datetime.isoformat(datetime.datetime.today())
    outfile='vlapipeline.'+datestring+'.log'
    logfile=open(outfile,'w')
    
    try:
        import pipeline
    except ImportError, e:
        print e
        print >>logfile, "Unable to import the CASA pipeline"
        regstate=False
    
    try:
        #Open context
        context = pipeline.Pipeline(context='last').context
        regstate=True
        print >>logfile,"VLA pipeline context created."
        print "VLA pipeline context created."
        print >>logfile,"Context verification - VLA pipeline regression PASSED."
        print "Context verification - VLA pipeline regression PASSED."
        
    except Exception, e:
        regstate=False
        print >>logfile,"VLA pipeline context NOT created."
        print "VLA pipeline context NOT created."
        print >>logfile,"Context verification - VLA pipeline regression FAILED."
        print "Context verification - VLA pipeline regression FAILED."
        print >>logfile, e
        print e

    #Test fluxscale values
    #Test that hifv_fluxboot stage was the 13th stage run (index 12)
    try:
        #Open context of regression pipeline run
        context = pipeline.Pipeline(context='last').context
        fluxlist = context.results[12].read()[0].flux_densities
        rtol=1.0e-5  #Relative Tolerance
        atol=1.0e-8  #Absolute Tolerance
        
        #Open context "standard" for comparison
        #Test commit
        #standard_context = pipeline.Pipeline(context=standard_file+'.context', path_overrides={'name':standard_file, 'output_dir':os.getcwd()}).context
        #standard_fluxlist = standard_context.results[12].read()[0].flux_densities
        #value_compare43 = 0.6934577681171487
        #value_compare = 0.716367318068  #CASA 4.5
        #value_compare = 0.717857716108  #CASA 4.6
        #value_compare = 0.716364780148  #CASA 4.6.144 r36095, pipeline r36209 (trunk)
        #value_compare = 0.718457023749  #CASA-test 5.0.11-DEV (r38177) 
        value_compare = 0.71857779577  #CASA-prerelease 5.0.0-187, pipeline r40156
        
        #result_bool = np.isclose(fluxlist[0][0], standard_fluxlist[0][0], rtol=rtol, atol=atol, equal_nan=False)
        result_bool = np.isclose(fluxlist[0][0], value_compare, rtol=rtol, atol=atol, equal_nan=False)
        
        print >>logfile, "Accepted test value is: ", value_compare, " from CASA 4.6.144 r36095 and pipeline r36150 (trunk)"
        print "Accepted test value is: ", value_compare, " from CASA 4.6.144 r36095 and pipeline r36209 (trunk)"
        print >>logfile, "Regression generated value is: ", fluxlist[0][0]
        print "Regression generated value is: ", fluxlist[0][0]
        
        if (result_bool):
            regstate=True
            print >>logfile,"hifv_fluxboot values match within relative tolerance of 1.0e-05 and absolute tolerance of 1.0e-08"
            print "hifv_fluxboot values match within relative tolerance of 1.0e-05 and absolute tolerance of 1.0e-08"
            print >>logfile,"hifv_fluxboot test PASSED."
            print "hifv_fluxboot test PASSED."
        else:
            regstate=False
            print >>logfile,"hifv_fluxboot values are not within tolerances."
            print "hifv_fluxboot values are not within tolerances."
            print >>logfile,"hifv_fluxboot test FAILED."
            print "hifv_fluxboot test FAILED."
            print >>logfile, e
            print e
        
    except Exception, e:
        regstate=False
        print >>logfile,"hifv_fluxboot values are not within tolerances."
        print "hifv_fluxboot values are not within tolerances."
        print >>logfile,"hifv_fluxboot test FAILED."
        print "hifv_fluxboot test FAILED."
        print >>logfile, e
        print e

    '''
    #Test flagging values
    #Test flagging statistics from the final applycal stage 
    try:
        #Open context of regression pipeline run
        context = pipeline.Pipeline(context='last').context
        flagsummary= context.results[14].read()[0].flagsummary
        
        #Open context "standard" for comparison
        
        standard_context = pipeline.Pipeline(context=standard_file+'.context', path_overrides={'name':standard_file, 'output_dir':os.getcwd()}).context
        standard_flagsummary = standard_context.results[14].read()[0].flagsummary
        try:
            assert cmp(flagsummary, standard_flagsummary) == 0
            regstate=True
            print >>logfile,"hifv_applycals flagsummary dictionaries match."
            print "hifv_applycals flagsummary dictionaries match."
            print >>logfile,"hifv_applycals test PASSED."
            print "hifv_applycals test PASSED."
        except Exception, e:
            regstate=False
            print >>logfile,"hifv_applycals flagsummary dictionaries do NOT match."
            print "hifv_applycals flagsummary dictionaries do NOT match."
            print >>logfile,"hifv_applycals test FAILED."
            print "hifv_applycals test FAILED."
            print >>logfile, e
            print e
            

    except Exception, e:
        regstate=False
        print >>logfile,"hifv_fluxboot values are not within tolerances."
        print "hifv_fluxboot values are not within tolerances."
        print >>logfile,"hifv_fluxboot test FAILED."
        print "hifv_fluxboot test FAILED."
        print >>logfile, e
        print e
    '''

    logfile.close()

def main():
    try:
        run()
        stats()
    except KeyboardInterrupt:
        print "Interrupt requested...exiting"
    except Exception:
        traceback.print_exc(file=sys.stdout)
    #sys.exit(0)

if __name__ == "__main__":
    main()
    print "Regstate:" , regstate
    if regstate:
    	print "Regression PASSED"
    else:
    	print "Regression FAILED"
