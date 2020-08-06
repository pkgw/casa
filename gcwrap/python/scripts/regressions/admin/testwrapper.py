""" Class to wrap a test to use with unittest framework."""

import os
import subprocess
import sys
import shutil
import inspect
import re
import string
import traceback
import casac
import unittest

import tests.selection_syntax as selection_syntax
import importlib
SyntaxTestClasses = [selection_syntax.SelectionSyntaxTest]

PYVER = str(sys.version_info[0]) + "." + str(sys.version_info[1])

AIPS_DIR = os.environ["CASAPATH"].split()[0]
DATA_DIR = AIPS_DIR+'/data'

SCRIPT_REPOS=AIPS_DIR + "/" + os.environ["CASAPATH"].split()[1] + '/lib/python' + PYVER + '/tests/'

if not os.access(SCRIPT_REPOS, os.F_OK):
    if os.access(AIPS_DIR+'/lib64', os.F_OK):
        SCRIPT_REPOS = AIPS_DIR+'/lib64/python' + PYVER + '/tests/'
    elif os.access(AIPS_DIR+'/lib', os.F_OK):
        SCRIPT_REPOS = AIPS_DIR+'/lib/python' + PYVER + '/tests/'
    else:            #Mac release
        SCRIPT_REPOS = AIPS_DIR+'/Resources/python/tests/'

class Helper():
    # This class is called when a test is not found. It will
    # raise an exception and make nose fail. This way, the not
    # found test will be counted as an error and won't be ignored.
    def __init__(self, name):
        self.tname = name
    
    def test_dummy(self):
        '''Helper function'''
        raise Exception("Cannot find test %s"%self.tname)

class UnitTest:
    def __init__(self,testname=''):
        """Take the name of a test file (without .py), wrap it and run"""
        self.testname = testname
        self.workdir = testname+'_work'
        self.scriptdir = SCRIPT_REPOS       
        self.datadir = [DATA_DIR]
        self.dataFiles = []
    
    def funcdesc(self):
        '''Name of test for FunctionTestCase'''
        return 'Test '+self.testname
    
    def funcSetup(self):        
        """Copy data files to local working directory"""
        
        dataFiles = self.dataFiles
        print('Searching for input data in %s'%(self.datadir))
        for datafile in dataFiles: 
            file = self.locatedata(datafile, self.datadir)
            #if data already exist, remove them
            if(file != ''):
                os.system('rm -rf '+ self.workdir+'/'+datafile)
            if(os.path.isdir(file)):
                shutil.copytree(file, self.workdir+'/'+datafile)
            if(os.path.isfile(file)):
                shutil.copy(file, self.workdir+'/'+datafile)
        
    def funcTeardown(self):
        """Remove data files from working directory"""
        
        dataFiles = self.dataFiles
        for datafile in dataFiles:
            file = self.workdir+'/'+datafile
            os.system('rm -rf ' + file)
        
        self.dataFiles = []

    def getFuncTest(self):
        print('-------------- Unit Test for %s ---------------'%self.testname)
        """Wrap a script using unittest"""
        testscript = self.searchscript(self.testname, self.scriptdir)
        
        # avoid creating a _work directory
        if (testscript == ""):
            return

        # copy test to local directory
        self.workdir = os.getcwd()
        self.getTest(testscript, self.testname,self.scriptdir, self.workdir)

        # import the test
        mytest = __import__(self.testname)
        importlib.reload(mytest)

        #get the data
        try:
            self.dataFiles = mytest.data()
        except:
            print('No data needed or found')
              
        # Wrap the test, funcSetup and funcTeardown in a FunctionTestCase and return it
        testcase = (unittest.FunctionTestCase(mytest.run,setUp=self.funcSetup,
                                              tearDown=self.funcTeardown,
                                              description=self.funcdesc()))
        
        return testcase


    def getUnitTest(self,list=[]):
        """Set up a unit test script to run wit nose"""    
        print('-------------- Unit Test for %s ---------------'%self.testname)
        if list:
            print('List of specific tests %s'%(list))
            
        # search for script in repository
        testscript = self.searchscript(self.testname, self.scriptdir)

        if (testscript == ""):
            testlist = []
            # Create a dummy list and return it so that nose
            # includes this test in the list of erroneous tests
            # instead of ignoring it.
            t = unittest.FunctionTestCase(Helper(self.testname).test_dummy)
            return [t]

        # copy test to local directory
        self.workdir = os.getcwd()
        self.getTest(testscript, self.testname,self.scriptdir, self.workdir)

        # import the test
        mytest = __import__(self.testname)
        importlib.reload(mytest)
        
        # get the classes
        classes = mytest.suite()
        testlist = []
        
        # Check if specific tests/classes were requested
        if not list:
#            print "no list"
            for c in classes:
                for attr, value in c.__dict__.items():                
#                    print attr, " = ", value
                    if len(attr) >= len("test") and \
                        attr[:len("test")] == "test" : \
                        testlist.append(c(attr))
                for base in SyntaxTestClasses:
                    if issubclass(c, base):
                        test_name_list = [t._testMethodName for t in testlist]
                        for attr, value in base.__dict__.items():
                            if len(attr) >= len("test") and \
                                attr[:len("test")] == "test" :
                                if attr not in test_name_list:
                                    testlist.append(c(attr))                        
        else:
            # verify if list contains classes and/or methods
            for input in list:
                for c in classes:
#                    print c
                    if self.isaclass(c,input):
#                        print "it is a class"
                        # It is a class. Get all its methods
                        for attr, value in c.__dict__.items():
                            if len(attr) >= len("test") and \
                                attr[:len("test")] == "test" : \
                                # append each test method to the list    
                                testlist.append(c(attr))
                        for base in SyntaxTestClasses:
                            if issubclass(c, base):
                                test_name_list = [t._testMethodName for t in testlist]
                                for attr, value in base.__dict__.items():
                                    if len(attr) >= len("test") and \
                                        attr[:len("test")] == "test" :
                                        if attr not in test_name_list:
                                            testlist.append(c(attr))
                    else:
                        # maybe it is a method. Get only this one
#                        print "maybe it is a method"
                        for attr, value in c.__dict__.items():
                            if input == attr:
                                testlist.append(c(attr))
                        for base in SyntaxTestClasses:
                            if issubclass(c, base):
                                test_name_list = [t._testMethodName for t in testlist]
                                for attr, value in base.__dict__.items():
                                    if input == attr:
                                        if attr not in test_name_list:
                                            testlist.append(c(attr))
                                                                                        
#            for attr, value in c.__dict__.iteritems():
#                print attr, value
#                if list:
#                    print 'There is a list'
#                    for test in list:
#                        print test
#                        if test == attr:
#                            testlist.append(c(attr))
#                else:
#                    print attr, " = ", value
#                   if len(attr) >= len("test") and \
#                        attr[:len("test")] == "test" : \
##                        attr.rfind('test') != -1 :
#                        testlist.append(c(attr))
#        print testlist
        return testlist
            

    def cleanup(self,workdir):
        # for safety, avoid removing the local directory
        if (workdir == '.'):
            workdir = '/tmp/utests'
        
        if os.path.isdir(workdir):
            print('Cleaning up '+ workdir)
            shutil.rmtree(workdir)


    def createDir(self, workdir):
        """Create a working directory"""
        if os.access(workdir, os.F_OK) is False:
            print(workdir+' does not exist, creating it')
            os.makedirs(workdir)


    def locatedata(self, datafile, datadir):
        
        for repository in datadir :

            #Skip hidden directories
            filter_hidden = ' | grep -vE "^\\."  | grep -vE "/\\."'
            
            #See if find understands -L or -follow (depends on find version)
            (err, a) = subprocess.getstatusoutput('find -L ' + repository+'/ 1>/dev/null 2>&1')
            if not err:
                findstr='find -L '+repository+'/ -name '+datafile+' -print 2>/dev/null' + filter_hidden
            else:
                findstr='find '+repository+'/ -follow -name '+datafile+' -print 2>/dev/null' + filter_hidden
            # A '/' is appended to the directory name; otherwise sometimes find doesn't find.
            #Also, ignore error messages such as missing directory permissions
            
            (find_errorcode, a)=subprocess.getstatusoutput(findstr)   # stdout and stderr
            #if find_errorcode != 0:
            #    print >> sys.stderr, "%s failed: %s" % (findstr, a)
            retval=''
            b=['']
            if(a!=''):
                b=string.split(a, '\n')
                retval=b[len(b)-1]
                if(len(b) > 1):
                    print('More than 1 file found with name '+datafile)
                print('Will use', retval)
                return retval
        raise Exception('Could not find datafile %s in the repository directories %s' \
              % (datafile, datadir))
 

    def searchscript(self, testname, scriptdir):
        """Search for the script"""
        print("Searching for script %s in %s" %(testname,scriptdir))                  
#        TestName=string.lower(testname)
        TestName = testname

        theScript=''
        numOfScript=0

        # search for DIR/<name>.py
        if os.path.isdir(scriptdir):
            allScripts=os.listdir(scriptdir)
        else:
            allScripts=[]
#        print "allScripts = ", allScripts
        for scr in allScripts:
#            print scriptdir, scr, testname
            if (scr == TestName + '.py'):
                theScript = scr
                numOfScript += 1  
                      
        if numOfScript == 0:
#            raise Exception, 'Could not find test %s' %TestName 
            print('ERROR: Could not find test %s' %TestName)
            return ""  
            
        if( numOfScript > 1) :
            print('More than 1 scripts found for name '+TestName)
            print('Using the following one '+ theScript)
            
        print("Found", theScript)
        return theScript

    def getTest(self, testnamek, testName, scriptdir, workdir):
        print('Copy the script to the working dir')
        if testnamek[0:6] == 'tests/':
            shutil.copy(scriptdir+'/'+testnamek,
                        workdir+'/'+testName+'.py')
        else:    
            shutil.copy(scriptdir+'/'+testnamek, \
                    workdir+'/')

    def isaclass(self,myclass,input):
        '''Check if input equals the class name'''
#        print "input=%s myclass=%s"%(input,myclass.__name__)
        if input == myclass.__name__:
            return True
            
        return False

    def copyTest(self, copyto='/tmp'):
        """Copy the test script to the given directory
           It will search for the test script defined in
           self.testname from the location defined in
           self.scriptdir. It raises an exception if it 
           cannot copy to destination."""
        
        TestName = self.testname

        theScript=''
        numOfScript=0

        if os.path.isdir(self.scriptdir):
            allScripts=os.listdir(self.scriptdir)
        else:
            allScripts=[]
            
        for scr in allScripts:
            if (scr == TestName + '.py'):                
                fscript = os.path.join(self.scriptdir, scr)
                numOfScript += 1  
                      
        if numOfScript == 0:
            print('ERROR: Could not find test %s' %TestName)
            return False  
            
        if( numOfScript > 1) :
            print('More than 1 scripts found for name '+TestName)
            print('Copying the following one '+ fscript)
                    
        try:
            shutil.copy(fscript, copyto)
        except:
            raise
        
        return True

    def getTestClasses(self, atest=''):
        '''Return the classes contained in a test script.
           It will return only the classes returned by the
           suite() function. If atest is not given it will
           get the classes from the self.testname script.
           atest --> a test name (without the .py extension)'''
        
        if atest == '':
            atest = self.testname
        # import the test
        mytest = __import__(atest)
        importlib.reload(mytest)
        
        # get the classes from the suite function
        classes = mytest.suite()
        return classes

class ExecTest(unittest.TestCase,UnitTest):
    """Wraps scripts to run with execfile"""
    
    def setup(self):
        self.workdir = self.testname+'_work'
        self.scriptdir = SCRIPT_REPOS       
        
        self.testscript = self.searchscript(self.testname,self.scriptdir)
        # avoid creating a _work directory
        if (self.testscript == ""):
            return
        
        # create a working directory
        self.cleanup(self.workdir)
        self.createDir(self.workdir)
        
        # copy test to workdir
        self.getTest(self.testscript, self.testname, self.scriptdir, self.workdir)
        thisDir = os.getcwd()
        os.chdir(self.workdir)    
           
    def testrun(self):
        #self.testname is defined in the calling function
        # run self.setup before calling this function
        """Run test with execfile"""
        
        # run the test
        a=inspect.stack()
        stacklevel=0
        for k in range(len(a)):
            if (string.find(a[k][1], 'ipython console') > 0):
                stacklevel=k
                break
        gl=sys._getframe(stacklevel).f_globals
        
        exec(compile(open(self.testscript, "rb").read(), self.testscript, 'exec'), gl)   
        

