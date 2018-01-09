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


regressionsDisabled = True
####################            Functions
def usage():
    print '========================================================================='
    print '\nRunTest will execute Python test(s) of CASA tasks.'
    print 'Usage:\n'
    print 'casa [casapy-options] -c runTest.py [options] test_name\n'
    print 'Options:'
    print '  no option              print this message and exit.'
    #print '  -a or --all            run all tests defined in '
    #print '                         trunk/gcwrap/python/scripts/tests/uTest_list.json'
    print '  <test_name>            run only <test_name> (more tests are separated by spaces).'
    print '  -f or --file <list>    run the tests defined in an ASCII file <list>; one test per line.'
    print '  -d or --datadir <dir>  set an env. variable to a directory, TEST_DATADIR=<dir> '
    print '                         that can be used inside the tests.'
    print '  -m or --mem            show the memory used by the tests and the number of files left open.'
    print '  -g or --debug          set casalog.filter to DEBUG.'
    print '  -l or --list           print the list of tests from & print the list of tags from '
    print '                         trunk/gcwrap/python/scripts/tests/uTest_list.json.'
    print '  -s or --classes        print the classes from a test script (those returned by suite()).'
    print '  -H or --Help           print this message and exit.\n'
    print 'NOTE: it will look for tests in the install directory, which usually is \r'
    print '      <casa_install_dir>/python/2.7/tests'
    print 'See documentation in: http://www.eso.org/~scastro/ALMA/CASAUnitTests.htm\n'
    print '=========================================================================='

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
    '''Set an environmental variable for the data directory'''
    absdatadir = os.path.abspath(datadir)
    os.environ.__setitem__('TEST_DATADIR',absdatadir)
    print "Setting TEST_DATADIR Environmental variable"
    return

def getclasses(testnames):
    '''Get the classes of a test script It will copy the test script to /tmp and remove it afterwards'''
    here = os.getcwd()
    tmpdir = '/tmp'
    for filename in testnames:
        if not filename.startswith("test_"):
            print "Cannot Get Classes for Regression Test"
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
        return

def readJSON(FILE):
    import json
    List = []
    #print FILE
    if(not os.path.exists(FILE)):
        print 'ERROR: List of tests does not exist'
        return []
    tests = json.load(open(FILE))
    return tests['testlist']

def list_tests():
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

def _find_regression_path():
    TESTS_DIR = CASA_DIR + "/" + os.environ["CASAPATH"].split()[1] + '/lib/python' + PYVER + '/regressions/'
    _potential_data_directories = ( "/opt/casa/data",
                                    "/home/casa/data",
                                    "/home/casa/data/trunk",
                                    "/home/casa/data/master",
                                    "/opt/casa/data/master",
                                    "/export/data/casa" )
    REGRESSION_DATA = filter(lambda x: os.access(x,os.F_OK),map(lambda y: y+"/regression",_potential_data_directories))
    if not os.access(TESTS_DIR, os.F_OK):
        if os.access(CASA_DIR+'/lib64', os.F_OK):
            TESTS_DIR = CASA_DIR+'/lib64/python' + PYVER + '/regressions/'
        elif os.access(CASA_DIR+'/lib', os.F_OK):
            TESTS_DIR = CASA_DIR+'/lib/python' + PYVER + '/regressions/'
        else:            #Mac release
            TESTS_DIR = CASA_DIR+'/Resources/python/regressions/'
    stack_frame_find()['TESTS_DIR']=TESTS_DIR
    return TESTS_DIR

def checkForMPI():
    if (not os.environ.has_key("OMP_NUM_THREADS")) or (not os.environ.has_key("OMPI_COMM_WORLD_SIZE")): 
        print "MPICASA is not enabled"
        return False
    print "MPICASA is enabled"
    return True

def checkForPIPELINE():
    try:
        import pipeline
    except ImportError, e:
        #print e
        print "Unable to import the CASA pipeline"
        return False
    return True

def checkForCASAGUIDE_DATA():
    if not os.path.isdir(os.environ["CASAPATH"].split()[0] + '/data/casaguidedata/'):
        print "Could Not Find CasaGuide Data"
        return False
    return True

def main(testnames=[]):

    # Global variable used by regression framework to determine pass/failure status
    global regstate  
    global regressionResult
    regressionResult = ''
    regstate = False
        
    listtests = testnames
    if listtests == '--Help' or listtests == []:
        usage()
        sys.exit()
        
       
    if listtests == 'all':
        raise Exception, 'Disabled Due to CAS-10844.'
        whichtests = 1
        # Get the full list of tests from file
        skip_mpi = False
        skip_casaguide = False
        listtests = []
        tests = readJSON(LISTofTESTS)
        mpiEnabled = checkForMPI()
        pipelineEnabled = checkForPIPELINE()
        casaguideEnabled = checkForCASAGUIDE_DATA()

        # Check For MPI
        for test in tests:
            if not mpiEnabled and "mpi" in test["tag"]:
                print "Skipping ", test['testScript'], "MPI Not Enabled"
                continue
            if not pipelineEnabled and "pipeline" in test["tag"]:
                print "Skipping ", test['testScript'], "MPI Not Enabled"
                continue
            if not casaguideEnabled and "casaguide" in test["tag"]:
                print "Skipping ", test['testScript'], "Casaguide"
                continue
            listtests.append(test['testScript'])


        if listtests == []:
            raise Exception, 'List of tests \"%s\" is empty or does not exist'%LISTofTESTS

    elif listtests == '-t' or listtests == '--subset':
        whichtests = 1
        listtests = []
        tests = readJSON(LISTofTESTS)
        mpiEnabled = checkForMPI()
        pipelineEnabled = checkForPIPELINE()
        casaguideEnabled = checkForCASAGUIDE_DATA()

        # Begin Run Test Suite

        for test in tests:
            # Run MPI tests. Check for MPI enabled 
            if ('mpi'in testTag):
                 if not mpiEnabled: raise Exception, 'MPICASA Not Enabled'

            # Run Pipeline tests. Check for Pipeline enabled 
            if ('pipeline'in testTag):
                if not pipelineEnabled:raise Exception, 'Pipeline Not Enabled'

            # Run Pipeline tests. Check for Pipeline enabled 
            if ('casaguide'in testTag):
                if not casaguideEnabled:raise Exception, 'Cannot Find Casa Guide Data in: ' + os.environ["CASAPATH"].split()[0] + '/data/casaguidedata'

            # Test Tag to run all Functional Tests or all regression tags
            if ("functionalTest" in testTag) or ("regression" in testTag):

                # BEGIN TODO TEMP
                if 'functionalTest' in testTag: 
                    raise Exception, 'Disabled Due to CAS-10844.'
                # END TODO TEMP

                # BEGIN TODO TEMP
                if regressionsDisabled:
                    if 'regression' in testTag: 
                        raise Exception, 'Disabled Regression Tests'
                # END TODO TEMP
                if not mpiEnabled and ("mpi" in test["tag"]): continue
                if not pipelineEnabled and ("pipeline" in test["tag"]): continue
                if not casaguideEnabled and ("casaguide" in test["tag"]): continue
                elif (testTag in test['testType']): 
                    listtests.append(test['testScript'])

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
                                listtests.append(test['testScript'])

                        # ****************** END TODO TEMP NOT TO RUN REGRESSIONS *******************
                        #listtests.append(test['testScript'])

            elif testTag in test["tag"]:
                # ****************** BEGIN TODO TEMP NOT TO RUN REGRESSIONS *******************
                # When regressions are fix, delete this section and uncomment last line
                if regressionsDisabled:
                    if "regression" in test["testType"]:
                        continue
                    else:
                        listtests.append(test['testScript'])

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
    
    # BEGIN TODO TEMP
    if regressionsDisabled:
        for listtest in listtests:
            if not listtest.startswith("test_"):
                raise Exception, 'Disabled Regression Tests'
    # END TODO TEMP
    print "Tests to run: ", listtests
    # Directories
    PWD = os.getcwd()
    WDIR = PWD+'/nosedir/'
    
    # Create a working directory
    workdir = WDIR
    print 'Creating work directory '+ workdir
    if os.access(workdir, os.F_OK) is False:
        os.makedirs(workdir)
    else:
        shutil.rmtree(workdir)
        os.makedirs(workdir)
    
    # Move to working dir
    os.chdir(workdir)
    
    # Create a directory for nose's xml files
    xmldir = WDIR+'xml/'
    if os.access(xmldir, os.F_OK) is False:
        os.makedirs(xmldir)
    else:
        shutil.rmtree(xmldir)
        os.makedirs(xmldir)
    
    print "Starting tests for %s: " %(listtests)
    global regression_test
    global test_name
    # ASSEMBLE and RUN the TESTS
    if not whichtests:
        '''Run all tests'''
        list = []
        print "RUNNING ALL TESTS"
        try:
            for f in listtests:
                TESTDIR = WDIR + '/' + str(f) + '/'
                if os.access(TESTDIR, os.F_OK) is False:
                    os.makedirs(TESTDIR)
                else:
                    shutil.rmtree(TESTDIR)
                    os.makedirs(TESTDIR)
                os.chdir(TESTDIR)
                if f.startswith("test_"):
                    print "-------------------------------------Test: %s---------------------------------------------"%f
                    _find_unit_path()
                    if not haslist(f):                
                        testcases = UnitTest(f).getUnitTest()
                        list = list+testcases
                else:
                    print "-------------------------------------Regression: %s---------------------------------------------"%f
                    name = f
                    _find_regression_path( )
                    #print "Copying %s to %s"%(TESTS_DIR + "tests/" + name + ".py",os.getcwd()+'/'+name + ".py")
                    print "------------------------------------------------------------------------------------------------------------------------"
                    print "starting test %s (%s)" % (name,TESTS_DIR + "tests/" + name + ".py")
                    print "------------------------------------------------------------------------------------------------------------------------"
                    shutil.copy2(TESTS_DIR + "tests/" + name + ".py",os.getcwd()+'/'+name + ".py")
                    global regression_test
                    global test_name
                    test_name = name
                    regression_test = __import__(name)
                    list = list+[unittest.FunctionTestCase(test_dummy)]

        except:
            traceback.print_exc()
                    
    elif (whichtests == 1):
        '''Run specific tests'''
        list = []
        # Get one list of tests
        tmp = []
        for i in range(0,len(listtests)):
            if haslist(listtests[i]):
                tmp = tmp + gettests(listtests[i])
        if len(tmp) > 0:
            listtests = tmp

        # Remost Duplicates
        #print "Removing Duplicates"
        listtests = remove_duplicates(listtests)

        for f in listtests:
            TESTDIR = WDIR + '/' + str(f) + '/'
            if os.access(TESTDIR, os.F_OK) is False:
                os.makedirs(TESTDIR)
            else:
                shutil.rmtree(TESTDIR)
                os.makedirs(TESTDIR)
            os.chdir(TESTDIR)

            if f.startswith("test_"):
                print "-------------------------------------Test: %s---------------------------------------------"%f
                _find_unit_path()
                if not haslist(f):                
                    testcases = UnitTest(f).getUnitTest()
                    list = list+testcases

                else:
                    ff = getname(f)
                    tests = gettests(f)
                    print "-------------------------------------Test: %s---------------------------------------------"%f
                    # allow splitting of large test groups into smaller chunks
                    # large long running test groups make parallel test scheduling
                    # a bit more complicated so splitting them to smaller groups
                    # helps
                    # syntax: [testsplit:chunk_index-number_of_chunks]
                    if len(tests) == 1 and tests[0].startswith('testsplit:'):
                        import math
                        testcases = UnitTest(ff).getUnitTest()
                        chk, nchk = map(int, tests[0].split(':')[1].split('-'))
                        if chk > nchk or chk < 1:
                            raise ValueError('testsplit chunk must be 1 <= nchunks')
                        nchk = min(len(testcases), nchk)
                        chksz = int(math.ceil(len(testcases) / float(nchk)))
                        offset = (chk - 1) * chksz
                        print 'running tests %d to %d' % \
                            (offset, min(offset + chksz, len(testcases)))
                        testcases = testcases[offset:offset + chksz]
                    else:
                        testcases = UnitTest(ff).getUnitTest(tests)
                    list = list+testcases                
            else:
                print "-------------------------------------Regression: %s---------------------------------------------"%f
                name = f
                _find_regression_path( )
                #print "Copying %s to %s"%(TESTS_DIR + "tests/" + name + ".py",os.getcwd()+'/'+name + ".py")
                print "------------------------------------------------------------------------------------------------------------------------"
                print "starting test %s (%s)" % (name,TESTS_DIR + "tests/" + name + ".py")
                print "------------------------------------------------------------------------------------------------------------------------"
                shutil.copy2(TESTS_DIR + "tests/" + name + ".py",os.getcwd()+'/'+name + ".py")
                test_name = name # Passing test_name as a global to test_dummy. Passing regression_test as a global to test_dummy. Import it as a test instance
                regression_test = __import__(name)
                list = list+[unittest.FunctionTestCase(test_dummy)]

                
        if (len(list) == 0):
            os.chdir(PWD)
            raise Exception, 'ERROR: There are no valid tests to run'
                                                                     
                
    # Run all tests and create a XML report
    xmlfile = xmldir+'nose.xml'
    #print list
    #sys.exit()

    # BEGIN TODO TEMP
    if len(list) > 1000:
        raise Exception, 'Disabled Due to CAS-10844.'
    # END TODO TEMP
    #sys.exit()
    try:
        if (HAVE_MEMTEST and MEM):
            if not os.environ.get('NOSE_XUNIT_FILE'):
                os.environ['NOSE_XUNIT_FILE'] = xmlfile
            regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-memtest","--verbosity=2","--memtest-file="+xmlfile], suite=list, addplugins=[memTest.MemTest()])
            print "Mem Test XML File: ", xmlfile
        else:
            regstate = nose.run(argv=[sys.argv[0],"-d","-s","--with-xunit","--verbosity=2","--xunit-file="+xmlfile], suite=list)
            print "Xunit File: ",xmlfile

        os.chdir(PWD)
    except:
        print "Failed to run one or more tests"
        traceback.print_exc()
    else:
        os.chdir(PWD)

# ------------------ NOTE ---------------------------------------------
# Once CASA moves to Python 2.7, the getpopt module should be replaced
# by argparse. The next section will need to be updated accordingly
# ---------------------------------------------------------------------

####################            Constants

# Python Version
PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

# CASA Directory
CASA_DIR = os.environ["CASAPATH"].split()[0]

# memory mode variable
MEM = 0

# Define which tests to run    
whichtests = 0

# List of Tests

LISTofTESTS = _find_unit_path() + "uTest_list.json"
HAVE_MEMTEST=True
try:
    import memTest
except:
    HAVE_MEMTEST = False

global testTag
testTag = ''
####################            Main
if __name__ == "__main__":

    original_datapath = casa['dirs']['data']
    # Get command line arguments
    if "-c" in sys.argv:
        # If called with ... -c runUnitTest.py from the command line,
        # then parse the command line parameters
        i = sys.argv.index("-c")
        if len(sys.argv) >= i + 2 and \
               re.compile("runTest\.py$").search(sys.argv[i + 1]):
            
        
            try:
                # Get only this script options
                opts,args=getopt.getopt(sys.argv[i+2:], "Hpalt:mgs:f:d:", ["Help","pybot","all","list","subset=","mem",
                                                                     "debug","classes=","file=",
                                                                     "datadir="])
                
            except getopt.GetoptError, err:
                # Print help information and exit:
                print str(err) # will print something like "option -a not recognized"
                usage()
                os._exit(2)
                
            # List of tests to run
            testnames = []
            
            # Boolean for file with tests.
            # One could allow the use of --file with specific tests given in
            # the command line by removing this option and appending to the
            # testnames list in the args handling
            hasfile = False
            alltests = False
            
            #If no option is given, show the Help page
            if opts == [] and args == []:
                usage()
                os._exit(0)
                
            #print args
            # All other options       
            for o, a in opts:
                #print opts
                if o in ("-H", "--Help"):
                    usage()
                    os._exit(0) 
                if o in ("-p", "--pybot"):
                    raise Exception, "Pybot option not avaliable"
                    use_pybot()
                    os._exit(0) 
                if o in ("-l", "--list"):
                    list_tests()
                    os._exit(0)
                if o in ("-s", "--classes"): 
                    testnames.append(a)
                    getclasses(testnames)
                    os._exit(0)
                if o in ("-m", "--mem"):
                    # run specific tests in mem mode            
                    MEM = 1
                elif o in ("-t", "--subset"):
                    testTag = a
                    testnames = o
                elif o in ("-g", "--debug"):
                    #Set the casalog to DEBUG
                    casalog.filter('DEBUG')
                elif o in ("-d", "--datadir"):
                    # This will create an environmental variable called
                    # TEST_DATADIR that can be read by the tests to use
                    # an alternative location for the data. This is used 
                    # to test tasks with MMS data
                    # directory with test data
                    datadir = a
                    if not os.path.isdir(datadir):                            
                        raise Exception, 'Value of --datadir is not a directory -> '+datadir  
                    
                    # Set an environmental variable for the data directory
                    # Also, overwrite casa paramaters dictionary to point to new data
                    print "Setting Test Data Dir: ", datadir
                    casa['dirs']['data'] = datadir
                    settestdir(datadir)
                    if not os.environ.has_key('TEST_DATADIR'):    
                        raise Exception, 'Could not create environmental variable TEST_DATADIR'                        
                        
                elif o in ("-a", "--all"):
                    alltests = True
                    whichtests = 0
                    testnames = 'all'
                    break
                elif o in ("-f", "--file"):
                    hasfile = True
                    if  os.path.isfile(a):
                        f = open(a,"r") 
                        for line in f:
                            try:
                                testnames.append(re.sub(r'[\n\r]+', '',line))
                            except:
                                raise Exception, " The list should contain one test per line."
                    else:
                        raise Exception, str(a) + " is Not a Valid input File"
                    
                else:
                    assert False, "unhandled option"

            # Deal with other arguments
            if args != [] and not hasfile and not alltests:
                testnames = args
                                        
        else:
            testnames = []
        
    else:
        # Not called with -c (but possibly execfile() from iPython)
        testnames = []

    try:
        main(testnames)
        casa['dirs']['data'] = original_datapath #Set Datapath to original datapath if '-d' option was given
    except:
        traceback.print_exc()


