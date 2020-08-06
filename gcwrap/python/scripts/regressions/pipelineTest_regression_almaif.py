import os
import time
import datetime
import regression_utility as regutl
import sys, traceback
import numpy as np
import pipeline
import shutil

pathname = os.environ.get('CASAPATH').split()[0]
rootdatapath = pathname+'/data/regression/pipeline/almaif/'

'''Initial ALMA pipeline regression
    March 2019
'''

THISHOME  = "working/"
startTime = 0.0
endTime = 0.0
startProc = 0.0
endProc = 0.0
EPS = 1e-5  # Logical "zero"
regstate = True

def printmsg(logfile, msg):
    print(msg, file=logfile)
    print(msg)
    return

def load_context(filename):
    with open(filename, 'rb') as picklefile:
        return pipeline.infrastructure.utils.pickle_load(picklefile)

def pipeline_regression():
    global regstate
    # global MIN_CASA_REVISION

    ASDM = rootdatapath  + "uid___A002_Xc46ab2_X15ae"
    fluxcsv_file = 'flux.csv'
    fluxcsv = rootdatapath + fluxcsv_file
    
    try:
        import pipeline.recipes.hifacal as hifacal
    except ImportError as e:
        print(e)
        
    
    # Check to see if the ASDM exists
    if not os.path.exists(ASDM):
        print(("Unable to open ASDM: ", ASDM))
	regstate=False
        raise IOError
    else:
        print(("Using {}".format(ASDM)))
        
        
    # copy flux.csv file 
    try:
        shutil.copy2( fluxcsv, fluxcsv_file )
    except:
        print(( "Could not copy {0} file from {1}".format(fluxcsv_file, fluxcsv) ))
        regstate=False
        raise IOError

    # verify if flux.csv exists
    if not os.path.exists( fluxcsv ):
        print(( "flux.csv file {} not found.".format(fluxcsv) ))
        regstate=False
        raise IOError
    
    # Run the CASA ALMA Pipeline standard calibration recipe
    try:
        hifacal.hifacal([ASDM], importonly=False, dbservice=False)
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
    print(("Run Time = ",endTime-startTime,endProc-startProc))

def stats():
    global startTime, endTime, startProc, endProc, regstate, standard_context_file
    
    datestring = datetime.datetime.isoformat(datetime.datetime.today())
    outfile = 'almapipeline.' + datestring + '.log'
    logfile = open(outfile,'w')
    
    try:
        import pipeline
    except ImportError as e:
        printmsg(logfile, "Unable to import the CASA pipeline.")
        printmsg(logfile, e)
        regstate = False
    
    try:
        # Open context
        context = pipeline.Pipeline(context='last').context
        regstate = True
        printmsg(logfile, "ALMA pipeline context created.")
        printmsg(logfile,"Context verification - ALMA pipeline regression PASSED.")        
    except Exception as e:
        regstate = False
        printmsg(logfile,"ALMA pipeline context NOT created.")
        printmsg("Context verification - ALMA pipeline regression FAILED.")
        printmsg(logfile, e)
        
    # Test imaging statistics
    # Test results of ALMA calibrator imaging from 19th stage run (index 18)
    try:
        # Open context of regression pipeline run
        context = pipeline.Pipeline(context='last').context
        imlist = context.results[18].read().results
        rtol = 1.0e-5  #Relative Tolerance
        atol = 1.0e-5  #Absolute Tolerance
        casaversion = 'CASA-prerelease 5.5.0-111'
        pipelinerevision = 'pipeline r42335 (trunk)'
        
        # Test image RMS
        value_compare = [0.00099723566425690248, 0.00075701982013768481, 0.0009601672821003053, 0.00094808680062035723, 0.0017550242121000309, 0.0020991357404819321, 0.0021395170847661651, 0.0024631465885457502]
        result_bool = [(np.isclose(imlist[i].image_rms, value_compare[i], rtol=rtol, atol=atol, equal_nan=False), os.path.basename(imlist[i].image)) for i in range(0,len(value_compare))]

        printmsg(logfile, "Accepted test RMS values are: {!s} from {!s}, {!s}".format(str(value_compare), casaversion, pipelinerevision))
        printmsg(logfile, "Regression generated values are: {!s}".format(str([imlist[i].image_rms for i in range(0,len(value_compare))])))
        
        for (rb, imagename) in result_bool:
            printmsg(logfile, ' ')
            printmsg(logfile, imagename)
            if (rb):
                regstate=True
                printmsg(logfile, "hif_makeimages RMS values match within relative tolerance of {!s} and absolute tolerance of {!s}".format(rtol,atol))
                printmsg(logfile, "hif_makeimages test PASSED.")
            else:
                regstate=False
                printmsg(logfile, "hif_makeimages RMS values are not within tolerances.")
                printmsg(logfile, "hif_makeimages test FAILED.")
    except Exception as e:
        regstate=False
        printmsg(logfile, "hif_makeimages pipeline error occurred.")
        printmsg(logfile, "hif_makeimages test FAILED.")
        printmsg(logfile, e)

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
    print(("Regstate:" , regstate))
    if regstate:
        print("Regression PASSED")
    else:
        print("Regression FAILED")
