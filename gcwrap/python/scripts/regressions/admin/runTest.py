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
import io
import contextlib
import logging
import argparse
import time

##
## testwrapper.py depends upon the current directory being in the path because
## it changes to the directory where the test is located and then imports it.
## CASA no longer leaves empty strings in sys.path to avoid confusion when
## stray files are in the current directory.
##
sys.path.insert(0,'')

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

def use_pybot(tests):
    logger.debug("Start def use_pybot(%s)",tests)
    tests = tests
    pybotWorkspace = os.environ["CASAPATH"].split()[0] + "/pybotWorkspace"
    if not os.path.isdir(pybotWorkspace):
        os.makedirs(pybotWorkspace)
    else:
        shutil.rmtree(pybotWorkspace)
        os.makedirs(pybotWorkspace)
    logger.debug("Pybot Workspace Directory: %s",pybotWorkspace)
    for test in tests:
        generateHTML(pybotWorkspace,test)
    # Fetch Testing Dir From Casa-pkg
    if logger.isEnabledFor(logging.DEBUG): verbose = "--verbose"
    else: verbose = "--quiet"

    if not os.path.isdir("casa-pkg"):
        subprocess.call(["/opt/local/bin/git","clone","https://open-bitbucket.nrao.edu/scm/casa/casa-pkg.git",verbose])
    else:
        shutil.rmtree("casa-pkg")
        subprocess.call(["/opt/local/bin/git","clone","https://open-bitbucket.nrao.edu/scm/casa/casa-pkg.git",verbose])

    shutil.copyfile(os.getcwd() + "/casa-pkg/testing/pybot-regression/__init__.html",pybotWorkspace+"/__init__.html")
    logger.debug("Copying %s to %s", os.getcwd() + "/casa-pkg/testing/pybot-regression/__init__.html",pybotWorkspace+"/__init__.html")

    shutil.copytree(os.getcwd() + "/casa-pkg/testing/pybot-regression/lib",pybotWorkspace + "/lib")
    logger.debug("Copying %s to %s", os.getcwd() + "/casa-pkg/testing/pybot-regression/lib",pybotWorkspace + "/lib")

    shutil.copytree(os.getcwd() + "/casa-pkg/testing/pybot-regression/bin",pybotWorkspace + "/bin")
    logger.debug("Copying %s to %s", os.getcwd() + "/casa-pkg/testing/pybot-regression/bin",pybotWorkspace + "/bin")

    os.environ['CASAROOT'] = os.environ["CASAPATH"].split()[0]
    return pybotWorkspace

def generateHTML(outdir,test):
    logger.debug("Start def generateHTML()")
    # Open Template File
    fp = open(HTML_Template, "rb" )
    msg = fp.read( )
    fp.close( )

    # Open Testlists
    tests = readJSON(LISTofTESTS)
    logger.debug("Generating HTML For: %s", test)
    for localtest in tests:
        if localtest['testScript'] == test:
            testNumber      = localtest['testNumber']
            testName        = localtest['testName']
            forceTag        = "TS" + str(localtest['priority'])
            testPyName      = localtest['testScript']
            maintainer      = localtest['Maintainer']
            maintainerEmail = localtest['MaintainerEmail']
            tag             = localtest['tag']
    templateInfo = msg % (testNumber, testName, testNumber, testName, forceTag, testPyName, maintainer, maintainerEmail,tag)
    filename = outdir + "/%s_%s.html"%(testNumber, testName.replace(" ", "_"))
    file = open(filename, "w")
    for line in  templateInfo:
        file.write(line)
    file.close()
    logger.debug("Generated %s", filename)

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
            #print "Cannot Get Classes for Regression Test"
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
                for attr, value in c.__dict__.items():
                    if len(attr) >= len("test") and attr[:len("test")] == "test":
                        print('\t%s'%c(attr))
            os.remove(filename+'.py')
            os.remove(filename+'.pyc')
        os.chdir(here)
    except:
        print('--> ERROR: Cannot copy script to %s'%tmpdir)
        logger.error('Failed to open file', exc_info=True)
        return

def getsubtests(filename,list=[]):
    logger.debug("Start def getsubtests()")
    logger.debug("Test: %s", filename)
    f = filename
    testlist_to_execute = []
    here = os.getcwd()
    tmpdir = '/tmp'
    os.chdir(tmpdir)
    tt = UnitTest(f)
    tt.copyTest(copyto=tmpdir)
    classes = tt.getTestClasses(f)
    logger.debug("Checking For Attribute: %s",ATTR_VAL)
    for c in classes:
        # Check if class has @attr(tag=''). Note: class @attr takes priority over func attr
        if 'tag' in c.__dict__:
            if c.tag == ATTR_VAL:
                for attr, value in c.__dict__.items():
                    if len(attr) >= len("test") and attr[:len("test")] == "test":
                        testlist_to_execute.append([attr,value.__module__])
        else:
            # Check if functions within each class has @attr(tag = '') or func.tag = ''
            for attr, value in c.__dict__.items():
                if len(attr) >= len("test") and attr[:len("test")] == "test":
                    if hasattr(value,'tag'):
                        if value.tag == ATTR_VAL:
                          testlist_to_execute.append([attr,value.__module__])
    os.remove(f+'.py')
    os.remove(f+'.pyc')
    os.chdir(here)

    return testlist_to_execute

def readJSON(FILE):
    logger.debug("Start Function: readJSON()")
    logger.debug("Reading JSON File: %s", FILE)
    import json
    List = []
    if(not os.path.exists(FILE)):
        #print 'ERROR: List of tests does not exist'
        logger.error('ERROR: List of tests does not exist')
        return []
    tests = json.load(open(FILE))
    return tests['testlist']

def list_tests():
    logger.info("Start Function: list_tests()")
    print('Full list of tests')
    print('-----------------------')
    for test in readJSON(LISTofTESTS):
        print(test['testType'], ":" ,test["testScript"])
    print('-----------------------')
    print('-----------------------')
    print('Full list of tags')
    tmpArray = []    
    for test in readJSON(LISTofTESTS):
        tmpArray.append(test['tag'])
    listtags = remove_duplicates(tmpArray)
    listtags.append('functionalTest')
    listtags.append('regression')
    listtags.append('TS1')
    listtags.append('TS2')
    #listtags.append('TS3') #'Disabled Due to CAS-10844'
    for listtag in listtags:
        print(listtag)

@contextlib.contextmanager
def stdoutIO(stdout=None):
    old = sys.stdout
    if stdout is None:
        stdout = io.StringIO()
    sys.stdout = stdout
    yield stdout
    sys.stdout = old

class RegressionTestCase(unittest.TestCase):
# Creates A testcase for regression tests to be executed in nose


    def test_dummy(self):
        '''Executing Regression Test'''
        from casac import casac
        ## flush output
        #sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        #sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

        PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

        CASA_DIR = os.environ["CASAPATH"].split()[0]
        TESTS_DIR = CASA_DIR + "/" + os.environ["CASAPATH"].split()[1] + '/lib/python' + PYVER + '/regressions/'

        _potential_data_directories = ( "/opt/casa/data","/home/casa/data","/home/casa/data/trunk","/home/casa/data/master","/opt/casa/data/master","/export/data/casa" )

        REGRESSION_DATA = [x for x in [y+"/regression" for y in _potential_data_directories] if os.access(x,os.F_OK)]

        if not os.access(TESTS_DIR, os.F_OK):
            if os.access(CASA_DIR+'/lib64', os.F_OK):
                TESTS_DIR = CASA_DIR+'/lib64/python' + PYVER + '/regressions/'
            elif os.access(CASA_DIR+'/lib', os.F_OK):
                TESTS_DIR = CASA_DIR+'/lib/python' + PYVER + '/regressions/'
            else:            #Mac release
                TESTS_DIR = CASA_DIR+'/Resources/python/regressions/'
        stack_frame_find()['TESTS_DIR']=TESTS_DIR

        if os.access(TESTS_DIR+name+".py",os.F_OK):
            path = TESTS_DIR+name+".py"
        elif os.access(TESTS_DIR+name+"_regression.py",os.F_OK):
            path = TESTS_DIR+name+"_regression.py"
        elif os.access(TESTS_DIR+name+"-regression.py",os.F_OK):
            path = TESTS_DIR+name+"-regression.py"
        elif os.access(TESTS_DIR+name+"_regression1.py",os.F_OK):
            path = TESTS_DIR+name+"_regression1.py"
        elif os.access(TESTS_DIR+name+"-regression1.py",os.F_OK):
            path = TESTS_DIR+name+"-regression1.py"
        elif os.access(TESTS_DIR+name+"_regression2.py",os.F_OK):
            path = TESTS_DIR+name+"_regression2.py"
        elif os.access(TESTS_DIR+name+"-regression2.py",os.F_OK):
            path = TESTS_DIR+name+"-regression2.py"
        elif os.access(TESTS_DIR+"test_"+name+".py",os.F_OK):
            path = TESTS_DIR+"test_"+name+".py"
        else:
            raise RuntimeError("task %s not found" % name)
        print()
        print("------------------------------------------------------------------------------------------------------------------------")
        print("starting test %s (%s)" % (name,path))
        print("------------------------------------------------------------------------------------------------------------------------")
        with stdoutIO() as regressionResult:
            #execfile(path,globals()) 
            exec(open(path).read(),globals())
        print(regressionResult.getvalue())
        if "PASS" in str(regressionResult.getvalue()).upper(): # Looks for PASS, PASSED
            assert True
        elif "FAIL" in str(regressionResult.getvalue()).upper():# Looks for FAIL, FAILED
            assert False #,"Regression: %s Failed"%(test_name)

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
    if ("OMP_NUM_THREADS" not in os.environ) or ("OMPI_COMM_WORLD_SIZE" not in os.environ): 
        print("MPICASA is not enabled")
        return False
    print("MPICASA is enabled")
    return True

def checkForPIPELINE():
    logger.debug("Start Function: def checkForPIPELINE()")
    try:
        import pipeline
    except ImportError as e:
        print(e)
        print("Unable to import the CASA pipeline")
        return False
    return True

def checkForCASAGUIDE_DATA():
    logger.debug("Start Function: checkForCASAGUIDE_DATA()")
    if not os.path.isdir(os.environ["CASAPATH"].split()[0] + '/data/casaguidedata/'):
        print("Could Not Find CasaGuide Data")
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
            raise Exception('Cannot Execute All Tests in MEM mode')
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
            raise Exception('List of tests \"%s\" is empty or does not exist'%LISTofTESTS)

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
                    raise Exception('MPICASA Not Enabled')

            # Run Pipeline tests. Check for Pipeline enabled 
            if ('pipeline'in testTag):
                if not pipelineEnabled:
                    logger.critical("Pipeline Not Enabled. Cannot Generate Suite for 'pipeline' Tag")
                    raise Exception('Pipeline Not Enabled')

            # Run casaguide tests. Check for casaguide data 
            if ('casaguide'in testTag):
                if not casaguideEnabled:
                    logger.critical('Cannot Find Casa Guide Data in: ' + os.environ["CASAPATH"].split()[0] + '/data/casaguidedata')
                    raise Exception('Cannot Find Casa Guide Data in: ' + os.environ["CASAPATH"].split()[0] + '/data/casaguidedata')

            # Test Tag to run all Functional Tests or all regression tags
            if ("functionalTest" in testTag) or ("regression" in testTag):
                if not mpiEnabled and ("mpi" in test["tag"]): continue
                if not pipelineEnabled and ("pipeline" in test["tag"]): continue
                if not casaguideEnabled and ("casaguide" in test["tag"]): continue
                elif (testTag in test['testType']): 
                    listtests.append(str(test['testScript']))

            # Run Priority Tests ( Originally mapped to TS1, TS2 , etc)
            elif testTag.startswith("TS"):
                # BEGIN TODO TEMP
                if 'TS3' in testTag: # TS3 has almost as many tests as --all option. Disabled till CAS-10844 is resolved
                    raise Exception('Disabled Due to CAS-10844.')
                # END TODO TEMP
                if not mpiEnabled and ("mpi" in test["tag"]): 
                    continue
                if not pipelineEnabled and ("pipeline" in test["tag"]):
                    continue
                if not casaguideEnabled and ("casaguide" in test["tag"]): 
                    continue
                for testPriority in range(1,int(testTag[2:])+1): # run test suites inclusively (i.e. 'TS1' runs TS1, 'TS2' runs TS1 tests & TS2 tests, etc. 
                    if (str(testPriority) in test['priority']): 
                        listtests.append(test['testScript'])
            elif testTag in test["tag"]:
                listtests.append(test['testScript'])

        if listtests == []:
            raise Exception('List of tests is empty')

    elif (type(testnames) != type([])):                
        if (os.path.isfile(testnames)):
            # How to prevent it from opening a real test???
            whichtests = 1
            listtests = readJSON(testnames)
            if listtests == []:
                raise Exception('List of tests is empty')
        else:
            raise Exception('List of tests does not exist')
            
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
    #sys.exit()
    global name
    if not (HAVE_ROBOT and USE_PYBOT):
        if not whichtests:
            '''Run all tests'''
            logger.info("Executing all tests")
            suiteList = []
            testlist_to_execute= []
            for f in listtests:
                suite = unittest.TestSuite()
                try:
                    if f.startswith("test_"):
                        tests = UnitTest(f).getUnitTest()
                        if RUN_SUBTEST:
                            testlist_to_execute = testlist_to_execute + getsubtests(f,tests)
                        if DRY_RUN:
                            tests = UnitTest(f).getUnitTest()
                            suiteList = suiteList + tests
                        else:
                            for test in tests:
                                suite.addTest(test)
                            logger.debug("\nTestName: %s: Suite: %s", f, suite)
                            suiteList.append(suite)
                    else:
                        name = f
                        suite = unittest.TestSuite()
                        suite.addTest(RegressionTestCase('test_dummy'))
                        list = list + [suite]
                except:
                    traceback.print_exc()
            list = suiteList

        # memTest.MemTest() plugin does not work with suites, requires 1 list of tests
        elif (whichtests == 1): 
            list = []
            suiteList = []
            testlist_to_execute= []
            suite = unittest.TestSuite()
            for f in listtests:
                try:
                    # Functional Tests
                    if f.startswith("test_"):
                        tests = UnitTest(f).getUnitTest() 
                        list = list+tests
                        if RUN_SUBTEST:
                            testlist_to_execute = testlist_to_execute + getsubtests(f,list)
                    # Regression Tests
                    else:
                        name = f
                        suite = unittest.TestSuite()
                        suite.addTest(RegressionTestCase('test_dummy'))
                        list = list + [suite]
                except: 
                        traceback.print_exc()
    elif (HAVE_ROBOT and USE_PYBOT):
        list = listtests

    if RUN_SUBTEST:
        if len(testlist_to_execute) == 0:
            raise ValueError("Cannot Find Tests with Attribute:'%s'"%(ATTR_VAL))
        if not whichtests:
            for i in range(0,len(list)):
                tmp = []
                for item in list[i]:
                    if [item._testMethodName,item.__module__] in testlist_to_execute:
                        tmp.append(item)
                list[i] =  unittest.TestSuite(tmp)
        else:
            tmp = []
            for item in list:
                if [item._testMethodName,item.__module__] in testlist_to_execute:
                    tmp.append(item)
            list = tmp
    if (len(list) == 0):
        os.chdir(PWD)
        logger.critical("ERROR: There are no valid tests to run")
        raise Exception('ERROR: There are no valid tests to run')

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

        elif (HAVE_ROBOT and USE_PYBOT):
            """Programmatic entry point for running tests.

            :param tests: Paths to test case files/directories to be executed similarly
                as when running the ``robot`` command on the command line.
            :param options: Options to configure and control execution. Accepted
                options are mostly same as normal command line options to the ``robot``
                command. Option names match command line option long names without
                hyphens so that, for example, ``--name`` becomes ``name``.

            Most options that can be given from the command line work. An exception
            is that options ``--pythonpath``, ``--argumentfile``, ``--escape`` ,
            ``--help`` and ``--version`` are not supported.

            Options that can be given on the command line multiple times can be
            passed as lists. For example, ``include=['tag1', 'tag2']`` is equivalent
            to ``--include tag1 --include tag2``. If such options are used only once,
            they can be given also as a single string like ``include='tag'``.

            Options that accept no value can be given as Booleans. For example,
            ``dryrun=True`` is same as using the ``--dryrun`` option.

            Options that accept string ``NONE`` as a special value can also be used
            with Python ``None``. For example, using ``log=None`` is equivalent to
            ``--log NONE``.

            ``listener``, ``prerunmodifier`` and ``prerebotmodifier`` options allow
            passing values as Python objects in addition to module names these command
            line options support. For example, ``run('tests', listener=MyListener())``.

            To capture the standard output and error streams, pass an open file or
            file-like object as special keyword arguments ``stdout`` and ``stderr``,
            respectively.

            A return code is returned similarly as when running on the command line.
            Zero means that tests were executed and no critical test failed, values up
            to 250 denote the number of failed critical tests, and values between
            251-255 are for other statuses documented in the Robot Framework User Guide.

            Example::

                from robot import run

                run('path/to/tests.robot')
                run('tests.robot', include=['tag1', 'tag2'], splitlog=True)
                with open('stdout.txt', 'w') as stdout:
                    run('t1.robot', 't2.robot', name='Example', log=None, stdout=stdout)

            Equivalent command line usage::

                robot path/to/tests.robot
                robot --include tag1 --include tag2 --splitlog tests.robot
                robot --name Example --log NONE t1.robot t2.robot > stdout.txt
            """
            if args.debug:
                robot.run(use_pybot(list),dryrun=DRY_RUN,xuint=xmlfile,outputdir=xmldir,console='verbose')
            else:
                robot.run(use_pybot(list),dryrun=DRY_RUN,xuint=xmlfile,outputdir=xmldir)

        else:
            cmd = [sys.argv[0],"-d","-s","--with-xunit","--verbosity=2","--xunit-file="+xmlfile]
            if DRY_RUN:
                cmd = cmd + ["--collect-only"]
            logger.debug('Executing nose.run(argv=%s, suite=%s)',cmd,list)
            regstate = nose.run(argv=cmd, suite=list)
            logger.debug("Regstate: %s", regstate)
            logger.info("XUNIT File: %s", xmlfile)

    except:
        print("Failed to run one or more tests")
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

# pybot mode variables
HAVE_ROBOT = True
USE_PYBOT = 0
try:
    import robot
except:
    HAVE_ROBOT = False


# Use Nose attribute Functionality 
RUN_SUBTEST = False

# Dry run of Tests
DRY_RUN = False

# Define which tests to run
whichtests = 0

# List of Tests
LISTofTESTS = _find_unit_path() + "uTest_list.json"

# HTML Template 
HTML_Template =  _find_unit_path() + "uTest_Template.dat"

# Define Generate HTML

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
            parser.add_argument("-p", "--pybot",action='store_true',help="BETA: Select and run tests using pybot (Robot Framework)") # TODO
            parser.add_argument("-r", "--attr",action='append',help='BETA: Select and run tests with attribute from <testname>.py')

            tmpArray = []
            for test in readJSON(LISTofTESTS):
                tmpArray.append(str(test['tag']))
            listtags = remove_duplicates(tmpArray)
            listtags.append('functionalTest')
            listtags.append('regression')
            listtags.append('TS1')
            listtags.append('TS2')
            listtags.append('TS3') #'Disabled Due to CAS-10844'
            listtags.append('TS4') #'Disabled Due to CAS-10844'

            parser.add_argument("-t", "--subset", choices=listtags,metavar='tag',help='run a suite of tests defined from tags in trunk/gcwrap/python/scripts/tests/uTest_list.json.')
            parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script (those returned by suite())')
            parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
            parser.add_argument("-z", "--coverage",action='store_true',help='BETA: show the coverage of the tests')

            parser.add_argument("--generate-html",action='store_true',help='Generate HTML files for use with pybot')

            args, unknownArgs = parser.parse_known_args(sys.argv[i+2:])

            # No option print help
            if sys.argv[i+2:] == []:
                parser.print_help()
                os._exit(0)

            # set casalog.filter to DEBUG. set runTest.log to DEBUG
            if args.debug:
                logger.info('Setting Debug Mode')
                casalog.filter('DEBUG')
                logger.setLevel(logging.DEBUG)
                logging.basicConfig(level=logging.DEBUG)
                handler.setLevel(logging.DEBUG)
                logger.debug("Starting Debug Mode")

            if args.generate_html:
                logger.info('Generating HTML Suite')
                fp = open(HTML_Template, "rb" )
                logger.debug("Reading HTML Template File: %s",HTML_Template)
                msg = fp.read( )
                fp.close( )
                if not os.path.isdir(os.getcwd()+ "/HTML"):
                    os.makedirs(os.getcwd()+ "/HTML")
                outdir = os.getcwd()+ "/HTML"
                logger.debug("Save Directory: %s",outdir)
                for test in readJSON(LISTofTESTS):
                    testNumber      = test['testNumber']
                    testName        = test['testName']
                    forceTag        = "TS" + str(test['priority'])
                    testPyName      = test['testScript']
                    maintainer      = test['Maintainer']
                    maintainerEmail = test['MaintainerEmail']
                    tag             = test['tag']
                    templateInfo = msg % (testNumber, testName, testNumber, testName, forceTag, testPyName, maintainer, maintainerEmail,tag)
                    filename = outdir + "/%s_%s.html"%(testNumber, testName.replace(" ", "_"))
                    logger.debug("Creating File: %s_%s.html",testNumber, testName.replace(" ", "_"))
                    file = open(filename, "w")
                    for line in  templateInfo:
                        file.write(line)
                    file.close()
                os._exit(0)

            # List All Avaliable Tests
            if args.list:
                list_tests()
                os._exit(0)

            # Print the classes from a test script
            if args.classes is not None:
                getclasses(args.classes)
                os._exit(0)

            if args.dry_run:
                # Some Parameters do not work with some plugins
                if args.coverage:
                    raise Exception("Cannot have a Dry Run in Coverage Mode")
                if args.mem:
                    raise Exception("Cannot have a Dry Run in Mem Mode")
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

            if args.attr: 
                logger.info('Setting Attr')
                RUN_SUBTEST = True
                if len(args.attr) != 1: raise Exception("Using multiple attributes not yet avaliable.")
                ATTR_VAL = args.attr[0]

            # Options to change test Execution Modes
            if args.mem: # run specific tests in mem mode
                logger.info('Setting Mem Mode')
                MEM = 1 

            if args.coverage: # run specific tests in cov mode
                if (HAVE_COVTEST):
                    logger.info('Setting Cov Mode')
                    COV = 1
                else:
                    logger.critical("You don't have module coverage installed: See https://pypi.python.org/pypi/coverage")
                    os._exit(0)

            if args.file is not None:
                logger.info('Reading Test List from %s: ', args.file)
                for line in args.file:
                    try:
                        logger.debug("Adding Test %s from file %s",re.sub(r'[\n\r]+', '',line),args.file)
                        testnames.append(re.sub(r'[\n\r]+', '',line))
                    except:
                        raise Exception(" The list should contain one test per line.")

            if args.pybot: #run specific tests in pybot mode
                if (HAVE_ROBOT):
                    logger.info('Setting Pybot Mode')
                    USE_PYBOT  = 1
                else:
                    logger.critical("You don't have module robot installed")
                    os._exit(0)

            if args.datadir is not None:
                '''This will create an environmental variable called TEST_DATADIR that can be read by the tests to use
                   an alternative location for the data. This is used to test tasks with MMS data directory with test data'''

                if not os.path.isdir(args.datadir):
                    raise Exception('Value of --datadir is not a directory -> '+ args.datadir)
                logger.debug("Data Path: %s",args.datadir)
                # Set an environmental variable for the data directory
                # Also, overwrite casa paramaters dictionary to point to new data

                logger.info("Setting Test Data Dir: %s",args.datadir)
                casa['dirs']['data'] = args.datadir
                settestdir(args.datadir)
                if 'TEST_DATADIR' not in os.environ:
                    logging.critical("Could not create environmental variable TEST_DATADIR")
                    raise Exception('Could not create environmental variable TEST_DATADIR')

            # Deal with other arguments
            for arg in unknownArgs:
                if arg.startswith(("-", "--")):
                    raise ValueError('unrecognized argument: %s'%(arg))
                #
                else:
                    testnames.append(arg)

            #If no tests are given, no subet tag or --all option
            if testnames == []:
                raise Exception('List of tests is empty')

    if logger.isEnabledFor(logging.DEBUG):
        # Debug Information
        logger.debug("PYVER: %s",PYVER)
        logger.debug("CASA_DIR: %s",CASA_DIR)
        logger.debug("HAVE_MEMTEST: %s",HAVE_MEMTEST)
        logger.debug("HAVE_COVTEST: %s",HAVE_COVTEST)
        logger.debug("USE_PYBOT: %s",USE_PYBOT)
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
