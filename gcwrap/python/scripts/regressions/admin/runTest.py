####################            Imports
import os
import sys
import getopt
import traceback
import unittest
import string
import re
import shutil
import pprint
import nose
from taskinit import casalog
from casa_stack_manip import *
import testwrapper
from testwrapper import *
import StringIO
import contextlib
import logging
import argparse
import time

####################            Temp Constants
regressionsDisabled = True

####################            Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)
# create a file handler
handler = logging.FileHandler('runTest-debug.log',mode='w')
handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(handler)


####################            Functions

def use_pybot():
    raise Exception, "Pybot option not avaliable"

def remove_duplicates(values):
    output = []
    seen = set()
    for value in values:
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output

def haslist(name):
    '''Check if specific list of tests have been requested'''
    n0 = name.rfind('[')
    n1 = name.rfind(']')
    if n0 == -1:
        return False
    return True

def getname(testfile):
    '''Get the test name from the command-line Ex: from test_clean[test1], returns test_clean'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        return testfile[:n0]

def gettests(testfile):
    '''Get the list of specific tests from the command-line Ex: from test_clean[test1,test3] returns [test1,test3]'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        temp = testfile[n0+1:n1]
        tests = temp.split(',')
        return tests

def settestdir(datadir):
    logger.debug("Start def settestdir()")
    '''Set an environmental variable for the data directory'''
    absdatadir = os.path.abspath(datadir)
    os.environ.__setitem__('TEST_DATADIR',absdatadir)
    #print "Setting TEST_DATADIR Environmental variable"
    logger.info("Setting TEST_DATADIR Environmental variable: %s",absdatadir)
    return

def getclasses(testnames):
    logger.debug("Start def getclasses()")
    '''Get the classes of a test script It will copy the test script to /tmp and remove it afterwards'''
    here = os.getcwd()
    tmpdir = '/tmp'
    for filename in testnames:
        if not filename.startswith("test_"):
            print "Cannot Get Classes for Regression Test"
            logger.error("Cannot Get Classes for Regression Test: %s",filename)
            return
    try:
        os.chdir(tmpdir)
        
        for filename in testnames:
            tt = UnitTest(filename)
            tt.copyTest(copyto=tmpdir)
            classes = tt.getTestClasses(filename)
            for c in classes:
                pprint.pprint('Class '+c.__name__)
                for attr, value in c.__dict__.iteritems():
                    if len(attr) >= len("test") and attr[:len("test")] == "test":
                        print '\t%s'%c(attr)
            os.remove(filename+'.py')
            os.remove(filename+'.pyc')
        os.chdir(here)
    except:
        print '--> ERROR: Cannot copy script to %s'%tmpdir
        logger.error('Failed to open file', exc_info=True)
        return

def readJSON(FILE):
    logger.debug("Start Function: readJSON()")
    logger.debug("Filename: %s", FILE)
    import json
    List = []
    if(not os.path.exists(FILE)):
        print 'ERROR: List of tests does not exist'
        logger.error('ERROR: List of tests does not exist')
        return []
    tests = json.load(open(FILE))
    return tests['testlist']

def list_tests():
    logger.info("Start Function: list_tests()")
    print 'Full list of tests'
    print '-----------------------'
    for test in readJSON(LISTofTESTS):
        print test['testType'], ":" ,test["testScript"]
    print '-----------------------'
    print '-----------------------'
    print 'Full list of tags'
    tmpArray = []    
    for test in readJSON(LISTofTESTS):
        tmpArray.append(test['tag'])
    listtags = remove_duplicates(tmpArray)
    listtags.append('functionalTest')
    if not regressionsDisabled:
        listtags.append('regression')
    listtags.append('TS1')
    listtags.append('TS2')
    #listtags.append('TS3') #'Disabled Due to CAS-10844'
    for listtag in listtags:
        print listtag

@contextlib.contextmanager
def stdoutIO(stdout=None):
    old = sys.stdout
    if stdout is None:
        stdout = StringIO.StringIO()
    sys.stdout = stdout
    yield stdout
    sys.stdout = old

# Creates A testcase for regression tests to be executed in nose
def test_dummy():
    '''Executing Regression Test'''
    print "Regression: %s"%(test_name)
    #test_result = regression_test.run(True)
    #TODO: Import regression setup and regression execution (if both are avaliable)
    with stdoutIO() as regressionResult:
        test_result = regression_test.run(False) # Some regressions are True/ False. This should be set within the test script
    print regressionResult.getvalue()
    if "PASS" in str(regressionResult.getvalue()).upper(): # Looks for PASS, PASSED
        assert True
    elif "FAIL" in str(regressionResult.getvalue()).upper():# Looks for FAIL, FAILED
        assert False,"Regression: %s Failed"%(test_name)

def _find_unit_path():
    TESTS_DIR = CASA_DIR + "/" + os.environ["CASAPATH"].split()[1] + '/lib/python' + PYVER + '/tests/'
    if not os.access(TESTS_DIR, os.F_OK):
        if os.access(CASA_DIR+'/lib64', os.F_OK):
            TESTS_DIR = CASA_DIR+'/lib64/python' + PYVER + '/tests/'
        elif os.access(CASA_DIR+'/lib', os.F_OK):
            TESTS_DIR = CASA_DIR+'/lib/python'+ PYVER +'/tests/'
        else:            #Mac release
            TESTS_DIR = CASA_DIR+'/Resources/python/tests/'
    return TESTS_DIR

def checkForMPI():
    logger.debug("Start Function: checkForMPI()")
    if (not os.environ.has_key("OMP_NUM_THREADS")) or (not os.environ.has_key("OMPI_COMM_WORLD_SIZE")): 
        print "MPICASA is not enabled"
        return False
    print "MPICASA is enabled"
    return True

def checkForPIPELINE():
    logger.debug("Start Function: def checkForPIPELINE()")
    try:
        import pipeline
    except ImportError, e:
        print e
        print "Unable to import the CASA pipeline"
        return False
    return True

def checkForCASAGUIDE_DATA():
    logger.debug("Start Function: checkForCASAGUIDE_DATA()")
    if not os.path.isdir(os.environ["CASAPATH"].split()[0] + '/data/casaguidedata/'):
        print "Could Not Find CasaGuide Data"
        return False
    return True

def run(testnames=[]):
    logger.info("Start Function: run(%s)",testnames)
    # Global variable used by regression framework to determine pass/failure status
    global regstate  
    global regressionResult
    global regression_test
    global test_name

    regressionResult = ''
    regstate = False
    listtests = testnames
    whichtests = 0

    if listtests == 'all':
        logger.info("Executing All Tests")
        if (HAVE_MEMTEST and MEM):
            logging.critical('Cannot Execute All Tests in MEM mode')
            raise Exception, 'Cannot Execute All Tests in MEM mode'
        whichtests = 0
        # Get the full list of tests from file
        skip_mpi = False
        skip_casaguide = False
        listtests = []
        tests = readJSON(LISTofTESTS)
        mpiEnabled = checkForMPI()
        pipelineEnabled = checkForPIPELINE()
        casaguideEnabled = checkForCASAGUIDE_DATA()
        logger.debug("mpiEnabled: %s",mpiEnabled)
        logger.debug("pipelineEnabled: %s",pipelineEnabled)
        logger.debug("casaguideEnabled: %s",casaguideEnabled)

        for test in tests:
            logger.debug("Test: %s", test['testScript'])
            if not mpiEnabled and "mpi" in test["tag"]:
                #print "Skipping ", test['testScript'], "MPI Not Enabled"
                logger.debug("Skipping %s %s", test['testScript'], "MPI Not Enabled")
                continue
            if not pipelineEnabled and "pipeline" in test["tag"]:
                logger.debug("Skipping %s %s", test['testScript'], "Pipeline Not Enabled")
                continue
            if not casaguideEnabled and "casaguide" in test["tag"]:
                logger.debug("Skipping %s %s", test['testScript'], "Casaguide")
                continue
            listtests.append(str(test['testScript']))
            logger.debug("To Execute Test: %s", test['testScript'])

        if listtests == []:
            logging.critical('List of tests \"%s\" is empty or does not exist', LISTofTESTS)
            raise Exception, 'List of tests \"%s\" is empty or does not exist'%LISTofTESTS

    elif listtests == 'subset':
        whichtests = 1
        listtests = []
        tests = readJSON(LISTofTESTS)
        mpiEnabled = checkForMPI()
        pipelineEnabled = checkForPIPELINE()
        casaguideEnabled = checkForCASAGUIDE_DATA()
        logger.debug("mpiEnabled: %s",mpiEnabled)
        logger.debug("pipelineEnabled: %s",pipelineEnabled)
        logger.debug("casaguideEnabled: %s",casaguideEnabled)
        logger.debug("TestTag: %s",testTag)
        # Begin Run Test Suite

        for test in tests:
            # Run MPI tests. Check for MPI enabled 
            if ('mpi'in testTag):
                 if not mpiEnabled:
                    logger.critical("MPICASA Not Enabled. Cannot Generate Suite for 'mpi' Tag")
                    raise Exception, 'MPICASA Not Enabled'

            # Run Pipeline tests. Check for Pipeline enabled 
            if ('pipeline'in testTag):
                if not pipelineEnabled:
                    logger.critical("Pipeline Not Enabled. Cannot Generate Suite for 'pipeline' Tag")
                    raise Exception, 'Pipeline Not Enabled'

            # Run casaguide tests. Check for casaguide data 
            if ('casaguide'in testTag):
                if not casaguideEnabled:
                    logger.critical('Cannot Find Casa Guide Data in: ' + os.environ["CASAPATH"].split()[0] + '/data/casaguidedata')
                    raise Exception, 'Cannot Find Casa Guide Data in: ' + os.environ["CASAPATH"].split()[0] + '/data/casaguidedata'

            # Test Tag to run all Functional Tests or all regression tags
            if ("functionalTest" in testTag) or ("regression" in testTag):
                # BEGIN TODO TEMP
                if regressionsDisabled:
                    if 'regression' in testTag: 
                        raise Exception, 'Disabled Regression Tests'
                # END TODO TEMP
                if not mpiEnabled and ("mpi" in test["tag"]): continue
                if not pipelineEnabled and ("pipeline" in test["tag"]): continue
                if not casaguideEnabled and ("casaguide" in test["tag"]): continue
                elif (testTag in test['testType']): 
                    listtests.append(str(test['testScript']))

            # Run Priority Tests ( Originally mapped to TS1, TS2 , etc)
            elif testTag.startswith("TS"):

                # BEGIN TODO TEMP
                if 'TS3' in testTag: # TS3 has almost as many tests as --all option. Disabled till CAS-10844 is resolved
                    raise Exception, 'Disabled Due to CAS-10844.'
                # END TODO TEMP
                if not mpiEnabled and ("mpi" in test["tag"]): 
                    continue
                if not pipelineEnabled and ("pipeline" in test["tag"]):
                    continue
                if not casaguideEnabled and ("casaguide" in test["tag"]): 
                    continue
                for testPriority in range(1,int(testTag[2:])+1): # run test suites inclusively (i.e. 'TS1' runs TS1, 'TS2' runs TS1 tests & TS2 tests, etc. 
                    if (str(testPriority) in test['priority']): 


                        # ****************** BEGIN TODO TEMP NOT TO RUN REGRESSIONS *******************
                        # When regressions are fix, delete this section and uncomment last line

                        if regressionsDisabled:
                            if "regression" in test["testType"]:
                                continue
                            else:
                                listtests.append(str(test['testScript']))

                        # ****************** END TODO TEMP NOT TO RUN REGRESSIONS *******************
                        #listtests.append(test['testScript'])

            elif testTag in test["tag"]:
                # ****************** BEGIN TODO TEMP NOT TO RUN REGRESSIONS *******************
                # When regressions are fix, delete this section and uncomment last line
                if regressionsDisabled:
                    if "regression" in test["testType"]:
                        continue
                    else:
                        listtests.append(str(test['testScript']))

                # ****************** END TODO TEMP NOT TO RUN REGRESSIONS *******************
                #listtests.append(test['testScript'])

        if listtests == []:
            raise Exception, 'List of tests is empty'

    elif (type(testnames) != type([])):                
        if (os.path.isfile(testnames)):
            # How to prevent it from opening a real test???
            whichtests = 1
            listtests = readJSON(testnames)
            if listtests == []:
                raise Exception, 'List of tests is empty'
        else:
            raise Exception, 'List of tests does not exist'
            
    else:
        # run specific tests
        whichtests = 1

    # Directories
    PWD = os.getcwd()
    WDIR = PWD+'/nosedir/'

    # Create a working directory
    workdir = WDIR
    #print 'Creating work directory '+ workdir
    logger.info("Creating work directory: %s", workdir)
    if os.access(workdir, os.F_OK) is False:
        os.makedirs(workdir)
    else:
        shutil.rmtree(workdir)
        os.makedirs(workdir)
    logger.debug("Working Directory: %s", workdir)

    # Move to working dir
    os.chdir(workdir)
    
    # Create a directory for nose's xml files
    xmldir = WDIR+'xml/'
    if os.access(xmldir, os.F_OK) is False:
        os.makedirs(xmldir)
    else:
        shutil.rmtree(xmldir)
        os.makedirs(xmldir)
    
    #print "Starting tests for %s: " %(listtests)
    logger.info("Starting Tests for: %s", listtests)
    # ASSEMBLE and RUN the TESTS
    if not whichtests:
        '''Run all tests'''
        logger.info("Executing all tests")
        list = []
        suiteList = []
        for f in listtests:
            suite = unittest.TestSuite()
            try:
                if f.startswith("test_"):
                    tests = UnitTest(f).getUnitTest()
                    if DRY_RUN:
                        tests = UnitTest(f).getUnitTest()
                        suiteList = suiteList + tests
                    else:
                        for test in tests:
                            suite.addTest(test)
                        logger.debug("\nTestName: %s: Suite: %s", f, suite)
                        suiteList.append(suite)
                else:
                    logger.debug("Disabled Regression Test: %s not compatiable with runTest.py",f)
                    pass
            except:
                traceback.print_exc()
        list = suiteList

    # memTest.MemTest() plugin does not work with suites, requires 1 list of tests
    elif (whichtests == 1): 
        list = []
        suiteList = []
        suite = unittest.TestSuite()
        for f in listtests:
            try:
                if f.startswith("test_"):
                    tests = UnitTest(f).getUnitTest() 
                    list = list+tests
                else:
                    logger.debug("Disabled Regression Test: %s not compatiable with runTest.py",f)
                    pass
            except:
                traceback.print_exc()

    if (len(list) == 0):
        os.chdir(PWD)
        logger.critical("ERROR: There are no valid tests to run")
        raise Exception, 'ERROR: There are no valid tests to run'

    # Run all tests and create a XML report
    xmlfile = xmldir+'nose.xml'
    if not os.environ.get('NOSE_XUNIT_FILE'):
        os.environ['NOSE_XUNIT_FILE'] = xmlfile
        logger.debug("Setting os.environ['NOSE_XUNIT_FILE']=%s",xmlfile)

    try:
        if (HAVE_MEMTEST and MEM):
            logger.debug('Executing nose.run(argv=[sys.argv[0],"-d","-s","--with-memtest","--verbosity=2","--memtest-file="%s"], suite=%s, addplugins=[memTest.MemTest()])',xmlfile,list)
            regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-memtest","--verbosity=2","--memtest-file="+xmlfile], suite=list, addplugins=[memTest.MemTest()])
            logger.debug("Regstate: %s", regstate)
            logger.info("Mem Test XML File: %s", xmlfile)

        elif (HAVE_COVTEST and COV):
            logger.warn("Due to the directory structure and various dependencies, Coverage Reporting May not Be Accurate")
            time.sleep(5) # Give 5 Seconds to View this warning
            logger.debug('Executing nose.run(argv=[sys.argv[0],"-d","-s","--with-coverage","--verbosity=2","--cover-xml-file="%s","--cover-html-dir="%s","--cover-html","--cover-erase","--cover-inclusive","--cover-branches"], suite=%s)',xmlfile,xmldir,list)
            regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-coverage","--verbosity=2","--cover-xml-file="+xmlfile,"--cover-html-dir="+xmldir,"--cover-html","--cover-erase","--cover-inclusive","--cover-branches"], suite=list)
            logger.debug("Regstate: %s", regstate)
            logger.info("Cov Test XML Directory: %s", xmldir)

        else:
            cmd = [sys.argv[0],"-d","-s","--with-xunit","--verbosity=2","--xunit-file="+xmlfile]
            if DRY_RUN:
                cmd = cmd + ["--collect-only"]
            logger.debug('Executing nose.run(argv=%s, suite=%s)',cmd,list)
            regstate = nose.run(argv=cmd, suite=list)
            logger.debug("Regstate: %s", regstate)
            logger.info("XUNIT File: %s", xmlfile)

    except:
        print "Failed to run one or more tests"
        logger.critical("Failed to run one or more tests", exc_info=True)
        traceback.print_exc()

    if logger.isEnabledFor(logging.DEBUG):
        logger.info("Debug Log: runTest-debug.log")

    os.chdir(PWD)

####################            Constants

# Python Version
PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

# CASA Directory
CASA_DIR = os.environ["CASAPATH"].split()[0]

# mem mode variables
HAVE_MEMTEST=True
MEM = 0
try:
    import memTest
except:
    HAVE_MEMTEST = False

# cov mode variables
HAVE_COVTEST=True
COV = 0
try:
    import coverage
except:
    HAVE_COVTEST = False

# Dry run of Tests
DRY_RUN = False

# Define which tests to run
whichtests = 0

# List of Tests
LISTofTESTS = _find_unit_path() + "uTest_list.json"

global testTag
testTag = ''

####################            Main

if __name__ == "__main__":
    logger.info('Start of Main')
    logger.info("Disabled Regression Testing. FunctionalTest will be executed.")
    original_datapath = casa['dirs']['data']
    # List of tests to run
    testnames = []

    # Get command line arguments
    if "-c" in sys.argv:
        # If called with ... -c runUnitTest.py from the command line,then parse the command line parameters
        i = sys.argv.index("-c")
        if len(sys.argv) >= i + 2 and re.compile("runTest\.py$").search(sys.argv[i + 1]):

            parser = argparse.ArgumentParser()
            parser.add_argument("-a", "--all", action='store_true',help='run all tests defined in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
            parser.add_argument("-d", "--datadir",help='set an env. variable to a directory, TEST_DATADIR=<dir> that can be used inside the tests.')
            parser.add_argument("-f", "--file",nargs='?', type=argparse.FileType('r'),help='run the tests defined in an ASCII file <list>; one test per line')
            parser.add_argument("-g", "--debug",action='store_true',help='Set casalog.filter to DEBUG. Set runTest-debug.log to DEBUG')
            parser.add_argument("-l", "--list",action='store_true',help='print the list of tests & tags defined in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
            parser.add_argument("-m", "--mem",action='store_true',help='show the memory used by the tests and the number of files left open.')
            parser.add_argument("-p", "--pybot",action='store_true',help=argparse.SUPPRESS) # TODO

            tmpArray = []
            for test in readJSON(LISTofTESTS):
                tmpArray.append(str(test['tag']))
            listtags = remove_duplicates(tmpArray)
            listtags.append('functionalTest')
            if not regressionsDisabled:
                listtags.append('regression')
            listtags.append('TS1')
            listtags.append('TS2')
            listtags.append('TS3') #'Disabled Due to CAS-10844'

            parser.add_argument("-t", "--subset", choices=listtags,metavar='tag',help='run a suite of tests defined from tags in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
            parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script (those returned by suite())')
            parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
            parser.add_argument("-z", "--coverage",action='store_true',help='show the coverage of the tests')

            args, unknownArgs = parser.parse_known_args(sys.argv[i+2:])

            # No option print help
            if sys.argv[i+2:] == []:
                parser.print_help()
                os._exit(0)

            if args.dry_run:
                # Some Parameters do not work with some plugins
                if args.coverage:
                    raise Exception, "Cannot have a Dry Run in Coverage Mode"
                if args.mem:
                    raise Exception, "Cannot have a Dry Run in Mem Mode"
                DRY_RUN = True

            # Run All Test
            if args.all:
                logger.info('Running All of Tests')
                whichtests = 0
                testnames = 'all'

            # Run Subset of Test
            if args.subset is not None:
                    logger.info('Running Subset of Tests')
                    testTag = args.subset
                    if args.all or testTag == 'TS3': # TS3 tag is the same as running --all flag
                        logger.info('Running All Tests')
                        whichtests = 0
                        testnames = 'all'
                    else:
                        testnames = 'subset'

            # List All Avaliable Tests
            if args.list:
                list_tests()
                os._exit(0)

            # Print the classes from a test script
            if args.classes is not None:
                getclasses(args.classes)
                os._exit(0)


            # Options to change test Execution Modes
            if args.mem: # run specific tests in mem mode
                logger.info('Setting Mem Mode')
                MEM = 1 

            if args.coverage: # run specific tests in cov mode
                if (HAVE_COVTEST):
                    logger.info('Setting Cov Mode')
                    COV = 1
                else:
                    print "You don't have module coverage installed: See https://pypi.python.org/pypi/coverage"
                    os._exit(0)

            # set casalog.filter to DEBUG. set runTest.log to DEBUG
            if args.debug:
                logger.info('Setting Debug Mode')
                casalog.filter('DEBUG')
                logger.setLevel(logging.DEBUG)
                logging.basicConfig(level=logging.DEBUG)
                handler.setLevel(logging.DEBUG)
                logger.debug("Starting Debug Mode")

            if args.file is not None:
                logger.info('Reading Test List from %s: ', args.file)
                for line in args.file:
                    try:
                        logger.debug("Adding Test %s from file %s",re.sub(r'[\n\r]+', '',line),args.file)
                        testnames.append(re.sub(r'[\n\r]+', '',line))
                    except:
                        raise Exception, " The list should contain one test per line."

            if args.pybot:
                # TODO Define use_pybot function to run pybot and produce a .html output file
                use_pybot()
                os._exit(0)

            if args.datadir is not None:
                '''This will create an environmental variable called TEST_DATADIR that can be read by the tests to use
                   an alternative location for the data. This is used to test tasks with MMS data directory with test data'''

                if not os.path.isdir(args.datadir):
                    raise Exception, 'Value of --datadir is not a directory -> '+ args.datadir
                logger.debug("Data Path: %s")
                # Set an environmental variable for the data directory
                # Also, overwrite casa paramaters dictionary to point to new data

                logger.info("Setting Test Data Dir: %s",args.datadir)
                casa['dirs']['data'] = args.datadir
                settestdir(args.datadir)
                if not os.environ.has_key('TEST_DATADIR'):
                    logging.critical("Could not create environmental variable TEST_DATADIR")
                    raise Exception, 'Could not create environmental variable TEST_DATADIR'

            # Deal with other arguments
            for arg in unknownArgs:
                if arg.startswith(("-", "--")):
                    raise ValueError('unrecognized argument: %s'%(arg))
                #
                else:
                    testnames.append(arg)

            #If no tests are given, no subet tag or --all option
            if testnames == []:
                raise Exception, 'List of tests is empty'

    if logger.isEnabledFor(logging.DEBUG):
        # Debug Information
        logger.debug("PYVER: %s",PYVER)
        logger.debug("CASA_DIR: %s",CASA_DIR)
        logger.debug("HAVE_MEMTEST: %s",HAVE_MEMTEST)
        logger.debug("HAVE_COVTEST: %s",HAVE_COVTEST)
        logger.debug("DRY_RUN: %s",DRY_RUN)
        logger.debug("whichtests: %s",whichtests)
        logger.debug("Tests: %s",testnames)
        logger.debug("TestTag: %s",testTag)
        logger.debug("Arguments: %s",args)

    try:
        run(testnames)
        casa['dirs']['data'] = original_datapath #Set Datapath to original datapath if '-d' option was given
    except:
        traceback.print_exc()

    logging.shutdown()
