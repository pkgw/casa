#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for ASDM version 1.2, 1.3 import to MS          #
#    and the "inverse filler" task exportasdm 
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The conversion of ASDM to MS and back needs to be verified.            #
#                                                                           # 
# Features tested:                                                          #
#    1) Is the import performed without raising exceptions                  #
#    2) Do all expected tables exist                                        #
#    3) Can the MS be opened                                                #
#    4) Do the tables contain expected values                               #
#    5) Is exportasdm performed without raising exceptions                  #
#    6) Is the created ASDM well-formed (XML) and complete                  #
#    7) Can the resulting ASDM be reimported without raising exceptions     #
#    8) Does it have the same number of integrations as the original        #
#    9) Can a mixed mode ASDM be imported in lazy mode                      #
#   10) Does the lazy mode when used with auto-only and the FLOAT_DATA      #
#       column produce an equivalent MS as the non-lazy mode does with the  #
#       same data selection                                                 #
#                                                                           #
# Input data:                                                               #
#     one dataset for the filler of ASDM 1.0                                #
#     one simulated MS dataset                                              #
#                                                                           #
#############################################################################
import os
import sys
import shutil
import numpy
from __main__ import default
from tasks import importasdm, flagdata, exportasdm, flagcmd
from taskinit import mstool, tbtool, cbtool, qatool, aftool, casalog
import testhelper as th
import flaghelper as fh
import unittest
import partitionhelper as ph
from parallel.parallel_data_helper import ParallelDataHelper
import recipes.ephemerides.convertephem as ce

myname = 'test_importasdm'

# default ASDM dataset name
myasdm_dataset_name = 'uid___X5f_X18951_X1'
myms_dataset_name = 'M51.ms'

# name of the resulting MS
msname = myasdm_dataset_name+'.ms'

# name of the exported ASDM
asdmname = myms_dataset_name+'.asdm'

# name of the reimported MS
reimp_msname = 'reimported-'+myms_dataset_name

# make local copies of the tools
tblocal = tbtool()
mslocal = mstool()

def checktable(themsname, thename, theexpectation):
    global myname
    tblocal.open(themsname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print(myname, ": comparing ", mycell)
        value = tblocal.getcell(mycell[0], mycell[1])
        # see if value is array
        try:
            isarray = value.__len__
        except:
            # it's not an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement = (value == mycell[2])
            else:
                in_agreement = ( abs(value - mycell[2]) < mycell[3]) 
        else:
            # it's an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement =  (value == mycell[2]).all() 
            else:
                in_agreement = (abs(value - mycell[2]) < mycell[3]).all() 
        if not in_agreement:
            print(myname, ":  Error in MS subtable", thename, ":")
            print("     column ", mycell[0], " row ", mycell[1], " contains ", value)
            print("     expected value is ", mycell[2])
            tblocal.close()
            return False
    tblocal.close()
    print(myname, ": table ", thename, " as expected.")
    return True

#########################

def verify_asdm(asdmname, withPointing):
    print("Verifying asdm ", asdmname)
    if(not os.path.exists(asdmname)):
        print("asdm ", asdmname, " doesn't exist.")
        raise Exception
    # test for the existence of all obligatory tables
    allTables = [ "Antenna.xml",
                  "ASDM.xml",
                 # "CalData.xml",
                 # "CalDelay.xml",
                 # "CalReduction.xml",
                  "ConfigDescription.xml",
                  "CorrelatorMode.xml",
                  "DataDescription.xml",
                  "ExecBlock.xml",
                  "Feed.xml",
                  "Field.xml",
                 #"FocusModel.xml",
                 #"Focus.xml",
                  "Main.xml",
                  "PointingModel.xml",
                  "Polarization.xml",
                  "Processor.xml",
                  "Receiver.xml",
                  "SBSummary.xml",
                  "Scan.xml",
                  "Source.xml",
                  "SpectralWindow.xml",
                  "State.xml",
                  "Station.xml",
                  "Subscan.xml",
                  "SwitchCycle.xml"
                  ]
    isOK = True
    for fileName in allTables:
        filePath = asdmname+'/'+fileName
        if(not os.path.exists(filePath)):
            print("ASDM table file ", filePath, " doesn't exist.")
            isOK = False
        else:
            # test if well formed
            rval = os.system('xmllint --noout '+filePath)
            if(rval !=0):
                print("Table ", filePath, " is not a well formed XML document.")
                isOK = False

    print("Note: xml validation not possible since ASDM DTDs (schemas) not yet online.")
        
    if(not os.path.exists(asdmname+"/ASDMBinary")):
        print("ASDM binary directory "+asdmname+"/ASDMBinary doesn't exist.")
        isOK = False

    if(withPointing and not os.path.exists(asdmname+"/Pointing.bin")):
        print("ASDM binary file "+asdmname+"/Pointing.bin doesn't exist.")
        isOK = False

    if (not isOK):
        raise Exception

# Setup for different data importing
class test_base(unittest.TestCase):
    def setUp_m51(self):
        res = None
        if(os.path.exists(myasdm_dataset_name)):
            shutil.rmtree(myasdm_dataset_name)

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/asdm-import/input/'
        shutil.copytree(datapath + myasdm_dataset_name, myasdm_dataset_name)
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/exportasdm/input/'
        shutil.copytree(datapath + myms_dataset_name, myms_dataset_name)
        default(importasdm)

    def setUp_xosro(self):
        self.asdm = 'X_osro_013.55979.93803716435'  #EVLA SDM
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/flagdata/'
        if(not os.path.lexists(self.asdm)):
            os.system('ln -s '+datapath+self.asdm +' '+self.asdm)
            
        default(importasdm)
        default(flagdata)

    def setUp_polyuranus(self):
        self.asdm = 'polyuranus'  # EVLA SDM with ephemeris
        datapath = os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/importevla/'
        if (not os.path.lexists(self.asdm)):
            os.system('ln -s '+datapath+self.asdm+' '+self.asdm)
        default(importasdm)
        default(flagdata)

    def setUp_autocorr(self):
        self.asdm = 'AutocorrASDM'  # ALMA 
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/importasdm/'
        if(not os.path.lexists(self.asdm)):
            os.system('ln -s '+datapath+self.asdm +' '+self.asdm)
            
        default(importasdm)

    def setUp_acaex(self):
        res = None
        myasdmname = 'uid___A002_X72bc38_X000' # ACA example ASDM with mixed pol/channelisation

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/asdm-import/input/'
        os.system('ln -sf '+datapath+myasdmname)
        default(importasdm)

    def setUp_12mex(self):
        res = None
        myasdmname = 'uid___A002_X71e4ae_X317_short' # 12m example ASDM with mixed pol/channelisation

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/asdm-import/input/'
        os.system('ln -sf '+datapath+myasdmname)
        default(importasdm)

    def setUp_eph(self):
        res = None
        myasdmname = 'uid___A002_X997a62_X8c-short' # 12m example ASDM with ephemerides

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/asdm-import/input/'
        os.system('ln -sf '+datapath+myasdmname)
        default(importasdm)

    def setUp_flags(self):
        res = None
        myasdmname = 'test_uid___A002_X997a62_X8c-short' # Flag.xml is modified

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/importasdm/'
        os.system('ln -sf '+datapath+myasdmname)
        default(importasdm)

    def setUp_SD(self):
        res = None
        myasdmname = 'uid___A002_X6218fb_X264' # Single-dish ASDM

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/alma-sd/M100/'
        os.system('ln -sf '+datapath+myasdmname)
        default(importasdm)

    def setUp_numbin(self):
        res = None
        # need full copies as this test involves renaming some xml files
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/asdm-import/input/'
        for this_asdm_name in ['alma_numbin_mixed','evla_numbin_2','evla_numbin_4']:
            if (os.path.exists(this_asdm_name)):
                shutil.rmtree(this_asdm_name)
            shutil.copytree(datapath + this_asdm_name, this_asdm_name)
        default(importasdm)


###########################
# beginning of actual test 
class asdm_import1(test_base):
    
    def setUp(self):
        self.setUp_m51()
        
    def tearDown(self):
        os.system('rm -rf '+myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        os.system('rm -rf reimported-M51.ms*')
        os.system('rm -rf myinput.ms')
        os.system('rm -rf '+asdmname)
                
    def test1(self):
        '''Asdm-import: Test good v1.2 input with filler v3 and inverse filler v3 '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importasdm(myasdm_dataset_name, useversion='v3')
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 1.008, 0],
                         ['DATA',      42, [ [10.5526886+0.0j] ], 1E-7]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
    
            expected = [
    # old values using TAI     ['UVW',       638, [-65.07623467,   1.05534109, -33.65801386], 1E-8],
                         ['UVW',       638, [-65.14758508, 1.13423277, -33.51712451], 1E-7],
                         ['EXPOSURE',  638, 1.008, 0],
                         ['DATA',      638, [ [0.00362284+0.00340279j] ], 1E-8]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [2202242.5520, -5445215.1570, -2485305.0920], 0.0001],
                         ['DISH_DIAMETER',1, 12.0, 0]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table ANTENNA failed'
            
            name = "POINTING"
            expected = [ ['DIRECTION',       10, [[ 1.94681283],[ 1.19702955]], 1E-8],
                         ['INTERVAL',        10, 0.048, 0],
                         ['TARGET',          10, [[ 1.94681283], [ 1.19702955]], 1E-8],
                         ['TIME',            10, 4758823736.016000, 1E-6],
                         ['TIME_ORIGIN',     10, 0., 0],
                         ['POINTING_OFFSET', 10, [[ 0.],[ 0.]], 0],
                         ['ENCODER',         10, [ 1.94851533,  1.19867576], 1E-8 ]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table POINTING failed'
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])
                
        myvis = myms_dataset_name
        os.system('rm -rf exportasdm-output.asdm myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        default('exportasdm')
        try:
            print("\n>>>> Test of exportasdm v3: input MS  is ", myvis)
            print("(a simulated input MS with pointing table)")
            rval = exportasdm(
                vis = 'myinput.ms',
                asdm = 'exportasdm-output.asdm',
                archiveid="S002",
                apcorrected=False,
                useversion='v3'
                )
            if not rval:
                raise Exception
            os.system('rm -rf '+asdmname+'; mv exportasdm-output.asdm '+asdmname)
            verify_asdm(asdmname, True)
        except:
            print(myname, ': *** Unexpected error exporting MS to ASDM, regression failed ***')   
            raise
            
        try:
            print("Reimporting the created ASDM (v3)....")
            importasdm(asdm=asdmname, vis=reimp_msname, wvr_corrected_data='no', useversion='v3')
            print("Testing existence of reimported MS ....")
            if(not os.path.exists(reimp_msname)):
                print("MS ", reimp_msname, " doesn't exist.")
                raise Exception
            print("Testing equivalence of the original and the reimported MS.")
            tblocal.open(myms_dataset_name)
            nrowsorig = tblocal.nrows()
            print("Original MS contains ", nrowsorig, "integrations.")
            tblocal.close()
            tblocal.open(reimp_msname)
            nrowsreimp = tblocal.nrows()
            tblocal.close()
            print("Reimported MS contains ", nrowsreimp, "integrations.")
            if(not nrowsreimp==nrowsorig):
                print("Numbers of integrations disagree.")
                raise Exception
        except:
            print(myname, ': *** Unexpected error reimporting the exported ASDM, regression failed ***')   
            raise

class asdm_import2(test_base):
    
    def setUp(self):
        self.setUp_m51()
        
    def tearDown(self):
        shutil.rmtree(myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        shutil.rmtree('myinput.ms', ignore_errors=True)
        os.system('rm -rf reimported-M51.ms*')
        shutil.rmtree('M51.ms.asdm', ignore_errors=True)
                
    def test_import2(self):
        '''Asdm-import: Test good v1.2 input with filler v3 and inverse filler v3 '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importasdm(myasdm_dataset_name, useversion='v3')
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
#                            "table.f0",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [ 0., 0., 0. ], 1E-8],
                         ['EXPOSURE',  42, 1.008, 0],
                         ['DATA',      42, [ [10.5526886+0.0j] ], 1E-7]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
            else:
                retValue['success']=True
    
            expected = [
    # old values using TAI     ['UVW',       638, [-65.07623467,   1.05534109, -33.65801386], 1E-8],
                         ['UVW',       638, [-65.14758508, 1.13423277, -33.51712451], 1E-7],
                         ['EXPOSURE',  638, 1.008, 0],
                         ['DATA',      638, [ [0.00362284+0.00340279j] ], 1E-8]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
            else:
                retValue['success']=True
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [2202242.5520, -5445215.1570, -2485305.0920], 0.0001],
                         ['DISH_DIAMETER',1, 12.0, 0]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table ANTENNA failed'
            else:
                retValue['success']=True
            
            name = "POINTING"
            expected = [ ['DIRECTION',       10, [[ 1.94681283],[ 1.19702955]], 1E-8],
                         ['INTERVAL',        10, 0.048, 0],
                         ['TARGET',          10, [[ 1.94681283], [ 1.19702955]], 1E-8],
                         ['TIME',            10, 4758823736.016000, 1E-6],
                         ['TIME_ORIGIN',     10, 0., 0],
                         ['POINTING_OFFSET', 10, [[ 0.],[ 0.]], 0],
                         ['ENCODER',         10, [ 1.94851533,  1.19867576], 1E-8 ]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table POINTING failed'
            else:
                retValue['success']=True
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])
                
        myvis = myms_dataset_name
        os.system('rm -rf exportasdm-output.asdm myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        default('exportasdm')
        try:
            print("\n>>>> Test of exportasdm v3: input MS  is ", myvis)
            print("(a simulated input MS with pointing table)")
            rval = exportasdm(
                vis = 'myinput.ms',
                asdm = 'exportasdm-output.asdm',
                archiveid="S002",
                apcorrected=False,
                useversion='v3'
                )
            if not rval:
                raise Exception
            os.system('rm -rf '+asdmname+'; mv exportasdm-output.asdm '+asdmname)
            verify_asdm(asdmname, True)
        except:
            print(myname, ': *** Unexpected error exporting MS to ASDM, regression failed ***')   
            raise
            
        try:
            print("Reimporting the created ASDM (v3)....")
            importasdm(asdm=asdmname, vis=reimp_msname, wvr_corrected_data='no', useversion='v3')
            print("Testing existence of reimported MS ....")
            if(not os.path.exists(reimp_msname)):
                print("MS ", reimp_msname, " doesn't exist.")
                raise Exception
            print("Testing equivalence of the original and the reimported MS.")
            tblocal.open(myms_dataset_name)
            nrowsorig = tblocal.nrows()
            print("Original MS contains ", nrowsorig, "integrations.")
            tblocal.close()
            tblocal.open(reimp_msname)
            nrowsreimp = tblocal.nrows()
            tblocal.close()
            print("Reimported MS contains ", nrowsreimp, "integrations.")
            if(not nrowsreimp==nrowsorig):
                print("Numbers of integrations disagree.")
                raise Exception
        except:
            print(myname, ': *** Unexpected error reimporting the exported ASDM, regression failed ***')   
            raise
        
class asdm_import3(test_base):
    # most of these (except test_CAS4532) were formerly found in test_importevla
    
    def setUp(self):
        # multiple SDMs used here, do not trust self.asdm to be correct, this just make a link to the data
        self.setUp_xosro()
        self.setUp_polyuranus()
        
    def tearDown(self):
        for myasdmname in ['X_osro_013.55979.93803716435', 'polyuranus']:
            os.system('rm -rf '+myasdmname)  # a link
            shutil.rmtree(myasdmname+".ms",ignore_errors=True)
            shutil.rmtree(myasdmname+".ms.flagversions",ignore_errors=True)
        os.system('rm -rf X_osro_013.55979.93803716435_cmd.txt')

    # functions to duplicate what importevla did after filling the MS - shadow and zero-level flagging
    def getmsmjds(self,msname):
        # Get start and end times from MS, return in mjds
        # this might take too long for large MS
        # NOTE: could also use values from OBSERVATION table col TIME_RANGE
        mslocal2 = mstool()
        success = True
        ms_startmjds = None
        ms_endmjds = None
        try:
            mslocal2.open(msname)
            timd = mslocal2.range(['time'])
            mslocal2.close()
        except:
            success = False
            print('Error opening MS ' + vis)
        if success: 
            ms_startmjds = timd['time'][0]
            ms_endmjds = timd['time'][1]
        else:
            print('WARNING: Could not open vis as MS to find times')
        return (ms_startmjds, ms_endmjds)

    def flagzero(self,startMJDsec, endMJDsec, flagpol):
        # returns a dictionary of flag commands related to zero data values
        # over the time interval given in the arguments
        # startMJDsec : start of interval, MJD, seconds
        # endMJDsec : end of interval, MJD, seconds
        # if flagpol = True, then flag all polarizations, else just RR and LL
        # there will be either a single command (flagpol=True) or
        # two commands (flagpol=False)

        result = {}

        flagz = {}
        flagz['time'] = 0.5*(startMJDsec+endMJDsec)
        flagz['interval'] = endMJDsec - startMJDsec 
        flagz['level'] = 0
        flagz['severity'] = 0
        flagz['type'] = 'FLAG'
        flagz['applied'] = False
        flagz['antenna'] = ''
        flagz['mode'] = 'clip'
        
        if flagpol:
            flagz['reason'] = 'CLIP_ZERO_ALL'
            command = {}
            command['mode'] = 'clip'
            command['clipzeros'] = True
            command['correlation'] = 'ABS_ALL'
            flagz['command'] = command
            flagz['id'] = 'ZERO_ALL'
            result[len(result)] = flagz.copy()
        else:
            flagz['reason'] = 'CLIP_ZERO_RR'
            command = {}
            command['mode'] = 'clip'
            command['clipzeros'] = True
            command['correlation'] = 'ABS_RR'
            flagz['command'] = command
            flagz['id'] = 'ZERO_RR'
            result[len(result)] = flagz.copy()
            
            flagz['reason'] = 'CLIP_ZERO_LL'
            command = {}
            command['mode'] = 'clip'
            command['clipzeros'] = True
            command['correlation'] = 'ABS_LL'
            flagz['command'] = command
            flagz['id'] = 'ZERO_LL'
            result[len(result)] = flagz.copy()

        return result

    def flagshadow(self,startMJDsec, endMJDsec, tolerance, addantenna):
        # create a set of flag commands related to shadowed data
        # the tolerance is given by tolerance
        # addantenna : may be a filename or a string
    
        result = {}
        command = {}
        # flag shadowed data
        result['time'] = 0.5 * (startMJDsec + endMJDsec)
        result['interval'] = endMJDsec - startMJDsec
        result['level'] = 0
        result['severity'] = 0
        result['type'] = 'FLAG'
        result['applied'] = False
        result['antenna'] = ''
        result['mode'] = 'shadow'
        result['reason'] = 'SHADOW'
    
        command['mode'] = 'shadow'
        command['tolerance'] = tolerance

        if type(addantenna) == str:
            if addantenna != '':
                # it's a filename, create a dictionary
                antdict = fh.readAntennaList(addantenna)
                command['addantenna'] = antdict
                    
        elif type(addantenna) == dict:
            if addantenna != {}:
                command['addantenna'] = addantenna
            
        result['command'] = command
        result['id'] = 'SHADOW'

        return result

    def applyflags(self,allflags, mspath):
        # aftool is agentflagger
        aflocal = aftool()
        if (len(allflags)) > 0:
            aflocal.open(mspath)
            aflocal.selectdata()
            fh.parseAgents(aflocal,allflags,[],True,True,'')
            aflocal.init()
            stats = aflocal.run(True,True)
            aflocal.done()
            # writes these to the FLAG_CMD table, append is True
            fh.writeFlagCommands(mspath,allflags,True,'','',True)
       
    def test_CAS4532(self):
        '''importasdm CAS-4532: white spaces on Antenna.xml'''
        # The X_osro_scan1/Antenna.xml and SpectralWindow.xml 
        # contain white spaces between some of the contents and 
        # the tags. This should not cause any error in the XML 
        # parser from fh.readXML
        import flaghelper as fh

        self.asdm = 'X_osro_013.55979.93803716435'
        
        flagcmddict = fh.readXML(self.asdm, 0.0)
        self.assertTrue(flagcmddict, 'Some XML file may contain white spaces not handled by readXML')
        
        self.assertEqual(list(flagcmddict.keys()).__len__(),214)

    def test_evlatest1(self):
        '''test of importing evla data, test1: Good input asdm'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }  

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm+".ms"

        # note that ocorr_mode defaulted to "co" for importevla and with_pointing_correction defaulted to True, so both must be set here
        # polyephem_tabtimestep defaulted to 0.001 in importevla, but that's only relevant for ephemeral objects, not relevant here
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2',ocorr_mode='co',with_pointing_correction=True)
        print(myname, ": Success! Now checking output ...")
        # this is probably not the best way to test sucess, too dependent on table system details, which might change without invalidating the result
        # but this is how it was written. These are the components it currently produces.
        # not a complete list of all of the table.* files in each directory
        mscomponents = set(["table.dat",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            # importevla copied these over, importasdm does not
                            #"Antenna.xml",
                            #"Flag.xml",
                            #"SpectralWindow.xml",
                            "ANTENNA/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.dat",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.dat",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.dat",
                            "FEED/table.f0",
                            "FIELD/table.dat",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.dat",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.dat",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.dat",
                            "OBSERVATION/table.f0",
                            "POINTING/table.dat",
                            "POINTING/table.f0",
                            "POLARIZATION/table.dat",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.dat",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.dat",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.dat",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.dat",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "SYSPOWER/table.f0",
                            "WEATHER/table.dat",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       42, [  1607.50778695, -1241.40287976 , 584.50368163], 1E-8],
                         ['EXPOSURE',  42, 1.0, 0]
                       ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
            else:
                retValue['success']=True
    
            expected = [
                         ['UVW',       638, [14.20193237, 722.59606805 , 57.57988905], 1E-8],
                         ['EXPOSURE',  638, 1.0, 0]
                       ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table MAIN failed'
            else:
                retValue['success']=True
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ -4.80000000e-12,  0.,  0.], 0],
                         ['POSITION',     1, [-1599644.8611, -5042953.6623, 3554197.0332], 0.0001],
                         ['DISH_DIAMETER',1, 25.0, 0]
                         ]
            results = checktable(msname, name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table ANTENNA failed'
            else:
                retValue['success']=True
            
        self.assertTrue(retValue['success'])

    def test_evlatest2(self):
        '''test of importing evla data, test2: Good input asdm with polynomial ephemeris'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.asdm = 'polyuranus'
        msname = self.asdm+".ms"

        # note that ocorr_mode defaulted to "co" for importevla and with_pointing_correction defaulted to True, so both must be set here
        # polyephem_tabtimestep defaulted to 0.001 must also be set here since this includes ephemeris data (unsure what the default is otherwise)
        # importevla did not have the convert_ephem2geo step so it must be disabled here to duplicate that behavior
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='0:5',ocorr_mode='co',with_pointing_correction=True,polyephem_tabtimestep=0.001,convert_ephem2geo=False)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
                            "table.f1",
                            "table.f2",
                            "table.f3",
                            "table.f4",
                            "table.f5",
                            "table.f6",
                            "table.f7",
                            "table.f8",
                            # importevla copied these over, importasdm does not
                            #"Antenna.xml",
                            #"Flag.xml",
                            #"SpectralWindow.xml",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "WEATHER/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FIELD/EPHEM0_URANUS_57545.8.tab/table.f0",  # the ephemeris table!
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            mslocal.open(msname)
            mssum = mslocal.summary()
            mslocal.close()

            if(mssum['scan_5']['0']['FieldName']=='URANUS' and mssum['field_2']['direction']['m0']['value']==0.37832756704958537):
                print(myname, ": MS summary as expected.")
                retValue['success']=True
            else:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Unexepcted mssummary.'
                
        self.assertTrue(retValue['success'])

    def test_evla_apply1(self):
        '''test of importing evla data, apply all flags, and save to file'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)
            
        # this applies and saves the online commands, which is all importasdm can do
        # note that there is no ephemeris data here so those arguments are not used in this call
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=True, applyflags=True, savecmds=True, flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. doing flagzero and shadow flagging")

                (startMJDs,endMJDs) = self.getmsmjds(msname)
                flagz = self.flagzero(startMJDs,endMJDs,True)
                flagh = self.flagshadow(startMJDs,endMJDs,0.0,'')
                allflags = flagz.copy()
                allflags[len(allflags)] = flagh
                self.applyflags(allflags,msname)
                fh.writeFlagCommands(msname,allflags,False,'',cmdfile,True)
                print(myname," : Checking flags")
            
                # Check flags
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],2446080)
        
                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
        
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist),216,'unexpected number of flag commands in saved ascii file : %s'%str(len(cmdlist)))
            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])
        
    def test_evla_apply2(self):
        '''test of importing evla data and applying the online flags'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
            
        # importasdm can everything importevla did in this case
        self.res  = importasdm(asdm=self.asdm, vis=msname, scans='2',ocorr_mode='co', with_pointing_correction=True,
                               process_flags=True, applyflags=True, savecmds=False, flagbackup=False)
        
        # Check flags
        res = flagdata(vis=msname, mode='summary')
        self.assertEqual(res['flagged'],2446080)
        self.assertEqual(res['scan']['2']['flagged'],2446080)
        
    def test_evla_apply3(self):
        '''test of importing evla data, do not apply online flags; post fill - apply clip zeros on RR and LL and save to file'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)
            
        # do NOT use the online flags here
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2,13', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=False,flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK, doing flagzero with flagpol=False")
                (startMJDs,endMJDs) = self.getmsmjds(msname)
                flagz = self.flagzero(startMJDs,endMJDs,False)
                self.applyflags(flagz,msname)
                fh.writeFlagCommands(msname,flagz,False,'',cmdfile,True)

                print(myname," : Checking flags")

                # Check flags - not the most useful test case
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'There are no zeros in this data set')
                self.assertEqual(res['scan']['2']['flagged'],0,'No flags should have been applied')
                self.assertEqual(res['scan']['13']['flagged'],0,'No flags should have been applied')
        
                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
        
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist),2,'Only clip zeros should be saved to file (2) : %s'%str(len(cmdlist)))

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])

    def test_evla_apply4(self):
        '''test of importing evla data: Save online flags to FLAG_CMD and file; do not apply'''

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)
            
        # importasdm can do all of this
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=True, applyflags=False,savecmds=True, flagbackup=False)
        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        ok = True
        try:
            mslocal.open(msname)
            mslocal.close()
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            ok = False

        self.assertTrue(ok,'Error opening MS')

        # No flags were applied
        res = flagdata(vis=msname, mode='summary')
        self.assertEqual(res['flagged'],0)
        
        # Apply only row 213 using flagcmd
        # The command in row 213 is the following:
        # antenna='ea06' timerange='2012/02/22/22:30:55.200~2012/02/22/22:35:08.199' 
        # spw='EVLA_X#A0C0#0' correlation='LL,LR,RL
        flagcmd(vis=msname, action='apply', tablerows=213)
        
        # Check flags. RR should no be flagged
        res = flagdata(vis=msname, mode='summary')
        self.assertEqual(res['correlation']['RR']['flagged'],0,'RR should not be flagged')
        self.assertEqual(res['correlation']['LL']['flagged'],29440)
        self.assertEqual(res['correlation']['LR']['flagged'],29440)
        self.assertEqual(res['correlation']['RL']['flagged'],29440)
        self.assertEqual(res['antenna']['ea06']['flagged'],88320)
        self.assertEqual(res['antenna']['ea07']['flagged'],3840,'Only a few baselines should be flagged')
        self.assertEqual(res['antenna']['ea08']['flagged'],3840,'Only a few baselines should be flagged')
        
        # Check output file existence       
        self.assertTrue(os.path.exists(cmdfile))
        
        # Check file content
        cmdlist = open(cmdfile,'r').readlines()
        self.assertEqual(len(cmdlist), 214, 'Only Online cmds should have been saved to file')

    def test_evla_apply5(self):
        '''test of importing evla data: Apply only shadow flags'''

        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        if os.path.exists(msname):
            os.system('rm -rf '+msname)
 
        self.res = importasdm(asdm=self.asdm, vis=msname, ocorr_mode='co', with_pointing_correction=True,
                              process_flags=False, flagbackup=False)
        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. shadow flagging")

                (startMJDs,endMJDs) = self.getmsmjds(msname)
                flagh = self.flagshadow(startMJDs,endMJDs,0.0,'')
                allflags = {}
                allflags[0] = flagh
                self.applyflags(allflags,msname)

                # This data set doesn't have shadow. Not a very useful sdm.
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'There are shadowed antenna in this data set')

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'

        self.assertTrue(retValue['success'])

    def test_evla_savepars(self):
        '''test importing evla data: save the flag commands and do not apply'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)

        self.res = importasdm(asdm=self.asdm, vis=msname,scans='11~13', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=True, applyflags=False, savecmds=True, flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. doing flagzero and shadow flagging")

                (startMJDs,endMJDs) = self.getmsmjds(msname)
                flagz = self.flagzero(startMJDs,endMJDs,True)
                flagh = self.flagshadow(startMJDs,endMJDs,0.0,'')
                allflags = flagz.copy()
                allflags[len(allflags)] = flagh

                # do NOT apply, but only append them to the command file
                fh.writeFlagCommands(msname,allflags,False,'',cmdfile,True)

                # Check flags - non should be applied
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'No flags should have been applied')

                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
        
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist), 216, 'Online, shadow and clip zeros should be saved to file')

                # NOW apply the flags
                # Apply flags using flagdata
                flagdata(vis=msname, mode='list', inpfile=cmdfile)
        
                # and check that they've been applied as expected
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],6090624)

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])
                
    def test_evla_apply1_flagdata(self):
        '''test of importing evla data, apply onlineflags, use flagdata to do shadow and zero level clipping'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)
            
        # this applies and saves the online commands, which is all importasdm can do
        # note that there is no ephemeris data here so those arguments are not used in this call
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=True, applyflags=True, savecmds=True, outfile=cmdfile, flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. doing flagzero and shadow flagging")

                # workaround here. flagdata cannot append to outfile, just overwrite it
                # create 2 independent cmdfiles and then append them to the one produced by importasdm
                clipCmdfile = msname.replace('.ms','_clip_cmd.txt')
                shadowCmdfile = msname.replace('.ms','_shadow_cmd.txt')

                flagdata(vis=msname, mode='clip',clipzeros=True,savepars=True,outfile=clipCmdfile,cmdreason='CLIP_ZERO_ALL')
                flagdata(vis=msname, mode='shadow',savepars=True,outfile=shadowCmdfile,cmdreason='SHADOW')

                # concatenate them
                os.system('cat %s >> %s' % (clipCmdfile,cmdfile))
                os.system('cat %s >> %s' % (shadowCmdfile,cmdfile))
                # remove the single line files
                os.system('rm %s' % clipCmdfile)
                os.system('rm %s' % shadowCmdfile)
                
                print(myname," : Checking flags")
            
                # Check flags
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],2446080)
        
                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
                
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist),216,'unexpected number of flag commands in saved ascii file : %s'%str(len(cmdlist)))
            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])
        
    def test_evla_apply3_flagdata(self):
        '''test of importing evla data, do not apply online flags; post fill - save to file, use flagdata'''
        # flagdata can not just clip the auto-correlation polationarizations. Can not duplicate intent of test_evla_apply3
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)
            
        # do NOT use the online flags here
        self.res = importasdm(asdm=self.asdm, vis=msname, scans='2,13', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=False,flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK, doing flagzero")  # equivalent to flagpol=False
                flagdata(vis=msname, mode='clip', clipzeros=True, correlation="ABS_RR,ABS_LL", savepars=True, cmdreason='CLIP_ZERO_AUTO_ONLY', outfile=cmdfile)
                
                print(myname," : Checking flags")

                # Check flags - not the most useful test case
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'There are no zeros in this data set')
                self.assertEqual(res['scan']['2']['flagged'],0,'No flags should have been applied')
                self.assertEqual(res['scan']['13']['flagged'],0,'No flags should have been applied')
        
                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
        
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist),1,'Only clip zeros should be saved to file (1) : %s'%str(len(cmdlist)))

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])

    def test_evla_apply5_flagdata(self):
        '''test of importing evla data: Apply only shadow flags; using flagdata'''

        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        if os.path.exists(msname):
            os.system('rm -rf '+msname)
 
        self.res = importasdm(asdm=self.asdm, vis=msname, ocorr_mode='co', with_pointing_correction=True,
                              process_flags=False, flagbackup=False)
        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. shadow flagging")
                flagdata(vis=msname, mode='shadow', savepars=True)

                # This data set doesn't have shadow. Not a very useful sdm.
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'There are shadowed antenna in this data set')

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'

        self.assertTrue(retValue['success'])

    def test_evla_savepars_flagdata(self):
        '''test importing evla data: save the flag commands and do not apply; using flagdata'''
        retValue = {'success':True, 'msg':"", 'error_msgs':''}

        self.asdm = 'X_osro_013.55979.93803716435'
        msname = self.asdm + ".ms"
        cmdfile = msname.replace('.ms','_cmd.txt')

        if os.path.exists(msname):
            os.system('rm -rf '+msname)
        if os.path.exists(cmdfile):
            os.system('rm -rf '+cmdfile)

        self.res = importasdm(asdm=self.asdm, vis=msname,scans='11~13', ocorr_mode='co', with_pointing_correction=True,
                              process_flags=True, applyflags=False, savecmds=True, outfile=cmdfile, flagbackup=False)

        print(myname, ": importasdm success. Checking that filled MS can be opened as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success'] = False
            retValue['error_msgs'] = retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()

        if (retValue['success']):
            try:
                print(myname," : OK. doing flagzero and shadow flagging")
                # workaround here. flagdata cannot append to outfile, just overwrite it
                # create 2 independent cmdfiles and then append them to the one produced by importasdm
                clipCmdfile = msname.replace('.ms','_clip_cmd.txt')
                shadowCmdfile = msname.replace('.ms','_shadow_cmd.txt')
                
                # no NOT apply
                flagdata(vis=msname, mode='clip',clipzeros=True,action='calculate',savepars=True,outfile=clipCmdfile,cmdreason='CLIP_ZERO_ALL')
                flagdata(vis=msname, mode='shadow',action='calculate',savepars=True,outfile=shadowCmdfile,cmdreason='SHADOW')

                # concatenate them
                os.system('cat %s >> %s' % (clipCmdfile,cmdfile))
                os.system('cat %s >> %s' % (shadowCmdfile,cmdfile))
                # remove the single line files
                os.system('rm %s' % clipCmdfile)
                os.system('rm %s' % shadowCmdfile)

                # Check flags - non should be applied
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],0,'No flags should have been applied')

                # Check output file existence
                self.assertTrue(os.path.exists(cmdfile))
        
                # Check file content
                cmdlist = open(cmdfile,'r').readlines()
                self.assertEqual(len(cmdlist), 216, 'Online, shadow and clip zeros should be saved to file')

                # NOW apply the flags
                # Apply flags using flagdata
                flagdata(vis=msname, mode='list', inpfile=cmdfile)
        
                # and check that they've been applied as expected
                res = flagdata(vis=msname, mode='summary')
                self.assertEqual(res['flagged'],6090624)

            except AssertionError as error:
                print(myname,' : assertion error while testing flags after filling: ' + str(error))
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msg']+str(error)
            except:
                print(myname," : post fill flagging and checking failed.")
                retValue['success'] = False
                retValue['error_msg']=retValue['error_msgs']+'Unexpected post-fill flagging result'
        self.assertTrue(retValue['success'])
                
class asdm_import4(test_base):
    
    def setUp(self):        self.setUp_autocorr()
        
    def tearDown(self):
        os.system('rm -rf '+self.asdm)
        os.system('rm -rf x54.ms*')
        os.system('rm -rf scan3.ms* autocorr.ms*')
        os.system('rm -rf scan3flags.txt')
        os.system('rm -rf scan3flags1.txt')
        os.system('rm -rf fbackup1.ms*')
        
    def test_autocorr(self):
        '''importasdm: auto-correlations should be written to online flags'''
        outfile='scan3flags.txt'
        importasdm(asdm=self.asdm, vis='x54.ms', scans='3', savecmds=True, outfile=outfile)
        self.assertTrue(os.path.exists(outfile))
        ff = open(outfile,'r')
        cmds = ff.readlines()
        self.assertEqual(cmds.__len__(), 2832)
        
        # auto-correlation should have been written to online flags               
        self.assertTrue(cmds[0].__contains__('&&*'))

    def test_autocorr_no_append(self):
        '''importasdm: online flags should not be appended to existing file'''
        outfile='scan3flags1.txt'
        importasdm(asdm=self.asdm, vis='x54.ms', scans='3', savecmds=True, outfile=outfile)
        self.assertTrue(os.path.exists(outfile))
        ff = open(outfile,'r')
        cmds = ff.readlines()
        self.assertEqual(cmds.__len__(), 2832)
        
        # Run again and verify that the online flags are not overwritten
        os.system('rm -rf x54.ms*')
        importasdm(asdm=self.asdm, vis='x54.ms', scans='3', savecmds=True, outfile=outfile, overwrite=False)
        print('Expected Error!')
        ff = open(outfile,'r')
        cmds = ff.readlines()
        self.assertEqual(cmds.__len__(), 2832)
        
       
    def test_flagautocorr1(self):
        '''importasdm: test that auto-correlations from online flags are correctly flagged'''        
        importasdm(asdm=self.asdm, vis='scan3.ms', scans='3', applyflags=True)
        res = flagdata('scan3.ms',mode='summary', basecnt=True)
        self.assertEqual(res['flagged'], 298)
        self.assertEqual(res['baseline']['DA44&&DA44']['flagged'], 76)
        self.assertEqual(res['baseline']['PM03&&PM03']['flagged'], 16)

    # CAS-4698, CSV-2576
    def test_flagautocorr2(self):
        '''importasdm: apply auto-correlations from online flags in flagcmd''' 
        outputms = 'autocorr.ms'       
        importasdm(asdm=self.asdm, vis=outputms, scans='3', applyflags=False)
        
        # Flag with flagcmd
        flagcmd(vis=outputms, action='apply', flagbackup=False)
        res = flagdata(vis=outputms, mode='summary', basecnt=True)
        self.assertEqual(res['flagged'], 298)
        self.assertEqual(res['baseline']['DA44&&DA44']['flagged'], 76)
        self.assertEqual(res['baseline']['PM03&&PM03']['flagged'], 16)

    # CAS-5286
    def test_flagautocorr3(self):
        '''importasdm: apply auto-correlations using autocorr=true in flagdata''' 
        outputms = 'autocorr.ms'       
        importasdm(asdm=self.asdm, vis=outputms, scans='3', process_flags=False,
                   flagbackup=False)
        
        # Applys with flagdata. Because all the data has PROCESSOR TYPE RADIOMETER
        # none of the auto-correlations should be flagged
        flagdata(vis=outputms, mode='manual', autocorr=True, flagbackup=False)
        res = flagdata(vis=outputms, mode='summary', basecnt=True)
        self.assertEqual(res['flagged'], 0)
        self.assertEqual(res['baseline']['DA44&&DA44']['flagged'], 0)
        self.assertEqual(res['baseline']['PM03&&PM03']['flagged'], 0)

    def test_flagbackup(self):
        '''importasdm: Create a flagbackup by default''' 
        outputms = 'fbackup1.ms'
        fbackup = outputms+'.flagversions' 
        
        # Do not create a flagbackup
        importasdm(asdm=self.asdm, vis=outputms, scans='3', flagbackup=False)
        self.assertFalse(os.path.exists(fbackup))

        # Create a backup by default
        importasdm(asdm=self.asdm, vis=outputms, scans='3', overwrite=True)
        self.assertTrue(os.path.exists(fbackup))
        
        # Delete it and create again using the parameter
        os.system('rm -rf '+outputms+'*')
        importasdm(asdm=self.asdm, vis=outputms, scans='3', flagbackup=True)
        self.assertTrue(os.path.exists(fbackup))
        

class asdm_import5(test_base):
    
    def setUp(self):
        self.setUp_m51()
        
    def tearDown(self):
        shutil.rmtree(myasdm_dataset_name)
        shutil.rmtree(myms_dataset_name)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        os.system('rm -rf reference.ms* reimported-M51.ms*')
        
                
    def test1_lazy1(self):
        '''Asdm-import: Test good v1.2 input with default filler in lazy mode, with_pointing_correction=True'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importasdm(myasdm_dataset_name, lazy=True, with_pointing_correction=True)
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.f17asdmindex",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", msname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+msname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm=myasdm_dataset_name, vis='reference.ms', lazy=False, overwrite=True)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.compTables('reference.ms', msname, ['FLAG', 'FLAG_CATEGORY', 'DATA','WEIGHT_SPECTRUM'], 
                                                    0.001) 
                retValue['success'] = th.compVarColTables('reference.ms', msname, 'DATA',1E-5) 
                retValue['success'] = th.compVarColTables('reference.ms', msname, 'FLAG') 

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 #"POINTING", # expect difference since with_pointing_correction=True
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                      msname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                        
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            msname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)

                print("\n*** Subtable POINTING")
                try:
                    retValue['success'] = retValue['success'] and not (th.compTables('reference.ms/POINTING', # expect difference since with_pointing_correction=True
                                                                                     msname+'/POINTING', [], 0.0))
                except:
                    retValue['success'] = False
                    print("ERROR: POINTING tables should differ in this test.")
            


                
        self.assertTrue(retValue['success'],retValue['error_msgs'])


class asdm_import6(test_base):

    def setUp(self):
        self.setUp_acaex()
        
    def tearDown(self):
        myasdmname = 'uid___A002_X72bc38_X000'
        os.system('rm '+myasdmname) # a link
        shutil.rmtree(myasdmname+".ms",ignore_errors=True)
        shutil.rmtree(myasdmname+".ms.flagversions",ignore_errors=True)
        os.system('rm -rf reference.ms*')
               
    def test6_lazy1(self):
        '''Asdm-import: Test good ACA ASDM with mixed pol/channelisation input with default filler in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X72bc38_X000'
        themsname = myasdmname + ".ms"

        self.res = importasdm(myasdmname, vis=themsname, lazy=True, scans='0:1~3') # only the first 3 scans to save time
        self.assertEqual(self.res, None)

        #test that scratch columns can be created from lazy import
        cblocal = cbtool()
        cblocal.open(themsname)
        cblocal.close()
        tblocal = tbtool()
        tblocal.open(themsname)
        self.assertIn('CORRECTED_DATA', tblocal.getdesc())
        self.assertIn('MODEL_DATA', tblocal.getdesc())
        tblocal.close()

        print(myname, ": Success! Now checking output ...")
        # crude test of success ...
        # the following tables should exist as directies containing table.dat and table.f0
        expectedTabs = ["ANTENNA","CALDEVICE","DATA_DESCRIPTION","FEED","FIELD",
                        "FLAG_CMD","HISTORY","OBSERVATION","POINTING","POLARIZATION",
                        "PROCESSOR","SOURCE","SPECTRAL_WINDOW","STATE","SYSCAL",
                        "SYSPOWER","WEATHER"]
        for name in expectedTabs:
            for fname in ["table.dat","table.f0"]:
                tabname = name+"/"+fname
                thisName = themsname+"/"+tabname
                if not os.access(thisName, os.F_OK):
                    print(myname, ": Error  ", thisName+"/"+name, "doesn't exist ...")
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+tabname+' does not exist'
                else:
                    print(myname, ": ", tabname+"/"+fname, "present.")
        if (retValue['success']):
            print(myname, ": MS exists. All tables present. Try opening as MS ...")
            try:
                mslocal.open(themsname)
            except:
                print(myname, ": Error  Cannot open MS table", themsname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
            else:
                mslocal.close()
                print(myname, ": OK. Checking tables in detail ...")
    
                importasdm(asdm=myasdmname, vis='reference.ms', lazy=False, overwrite=True, scans='0:1~3')

                if(os.path.exists('reference.ms')):
                    retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                           +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                    if not retValue['success']:
                        print("ERROR: DATA does not agree with reference.")
                    else:
                        print("DATA columns agree.")
                    retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                   +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                    if not retValueTmp:
                        print("ERROR: FLAG does not agree with reference.")
                    else:
                        print("FLAG columns agree.")

                    retValue['success'] = retValue['success'] and retValueTmp

                    for subtname in expectedTabs:
                        print("\n*** Subtable ",subtname)
                        excllist = []
                        if subtname=='SOURCE':
                            excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                        if subtname=='SYSCAL':
                            excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                        if subtname=='CALDEVICE':
                            excllist=["NOISE_CAL","CAL_EFF"]
                        if subtname=='HISTORY':
                            excllist=['MESSAGE']
                        if subtname=="SYSCAL":
                            exclist=["TANT_SPECTRUM","TANT_TSYS_SPECTRUM"]
                        if subtname=='SPECTRAL_WINDOW':
                            excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION']
                            for colname in excllist: 
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                        
                        try:    
                            retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                                themsname+'/'+subtname, excllist, 
                                                                0.01) and retValue['success']
                        except:
                            retValue['success'] = False
                            print("ERROR for table ", subtname)

                    if retValue['success']:
                        # For this MS, the WEATHER/DEW_POINT_FLAG is all True and
                        # the other *FLAG columns are all false, just check PRESSURE_FLAG
                        tblocal.open(themsname+'/WEATHER')
                        dewptFlag = tblocal.getcol('DEW_POINT_FLAG')
                        pressFlag = tblocal.getcol('PRESSURE_FLAG')
                        tblocal.close()
                        retValue['success'] = (dewptFlag.sum()==len(dewptFlag)) and (pressFlag.sum()==0)
                        if not (retValue['success']):
                            print("ERROR for WEATHER table. Expected values of DEW_POINT_FLAG and PRESSURE_FLAG are not seen.")
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
class asdm_import7(test_base):

    def setUp(self):
        self.setUp_12mex()
        self.setUp_eph()
        self.setUp_SD()
        
    def tearDown(self):
        for myasdmname in ['uid___A002_X71e4ae_X317_short', 'uid___A002_X997a62_X8c-short', 'uid___A002_X6218fb_X264']:
            os.system('rm -f '+myasdmname) # a link
            shutil.rmtree(myasdmname+".ms",ignore_errors=True)
            shutil.rmtree(myasdmname+'.ms.flagversions',ignore_errors=True)
        shutil.rmtree("reference.ms",ignore_errors=True)
        shutil.rmtree("reference.ms.flagversions",ignore_errors=True)
        shutil.rmtree("uid___A002_X997a62_X8c-short.interp.ms",ignore_errors=True)
               
    def test7_lazy1(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, lazy=True, scans='0:1~4') # only the first 4 scans to save time
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "WEATHER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
            mslocal.close()
            print(myname, ": MS can be opened. Now testing the changing of the asdmref ...")
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)
            
            mslocal.open(themsname)
            
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm="moved_"+myasdmname, vis='reference.ms', lazy=False, overwrite=True, scans='0:1~4')

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                       +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                               +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")
                    
                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "CALDEVICE",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL",
                                 "SYSPOWER",
                                 "WEATHER"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='CALDEVICE':
                        excllist=['NOISE_CAL','CAL_EFF']
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)

                try:
                    # test that the PRESSURE column has the expected units
                    wxcalOK = tblocal.open(themsname+'/WEATHER')
                    if wxcalOK:
                        wxcalOK = tblocal.getcolkeyword("PRESSURE","QuantumUnits")[0] == 'hPa'
                        tblocal.close()
                    retValue['success'] = wxcalOK and retValue['success']
                    if not wxcalOK:
                        print("PRESSURE column in WEATHER table is missing or has incorrect units")
                except:
                    retValue['success'] = False
                    print("ERROR getting units of PRESSURE column in WEATHER table.")

                try:
                    # test that the SDM_WINDOW_FUNCTION column exists and has the exepcted values
                    winFuncOK = tblocal.open(themsname+'/SPECTRAL_WINDOW')
                    if winFuncOK:
                        winFuncCol = tblocal.getcol('SDM_WINDOW_FUNCTION')
                        tblocal.close()
                        # test values here
                        # expect 55 rows, rows 1:24 are HANNING, the rest are UNIFORM
                        indx = numpy.arange(len(winFuncCol))
                        winFuncOK = winFuncOK and (numpy.array_equal(indx[winFuncCol=="HANNING"],(numpy.arange(24)+1)))
                        winFuncOK = winFuncOK and (len(indx[winFuncCol=="UNIFORM"])==31)

                    retValue['success'] = winFuncOK and retValue['success']
                    if not winFuncOK:
                        print("SDM_WINDOW_FUNCTION column in the SPECTRAL_WINDOW table is missing or has incorrect values")
                except:
                    retValue['success'] = False
                    print("ERROR checking the value of the SDM_WINDOW_FUNCTION column in the SPECTRAL_WINDOW table.")

                try:
                    # test that the SDM_NUM_BIN column exists and has the exepcted values
                    numBinOK = tblocal.open(themsname+'/SPECTRAL_WINDOW');
                    if numBinOK:
                        numBinCol = tblocal.getcol('SDM_NUM_BIN');
                        tblocal.close()
                        # all values are 1
                        numBinOK = numpy.all(numBinCol==1)

                    retValue['success'] = numBinOK and retValue['success']
                    if not numBinOK:
                        print("SDM_NUM_BIN column in the SPECTRAL_WINDOW table is missing or has incorrect values")
                except:
                    retValue['success'] = False
                    print("ERROR checking the value of the SDM_NUM_BIN column in the SPECTRAL_WINDOW table.")

        os.system("mv moved_"+myasdmname+" "+myasdmname)
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test7_lazy2(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode with reading the BDF flags'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, lazy=True, bdfflags=True) 
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm=myasdmname, vis='reference.ms', lazy=True, overwrite=True, bdfflags=False)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")
                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") != 0
                if not retValueTmp:
                    print("ERROR: FLAG columns do agree with reference but they shouldn't.")
                else:
                    print("FLAG columns do not agree as expected.")

                retValue['success'] = retValue['success'] and retValueTmp

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)
            
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test7_lazy3(self):
        '''Asdm-import: Test good 12 m ASDM with Ephemeris table in lazy mode'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X997a62_X8c-short'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, lazy=True, convert_ephem2geo=True, 
                              process_pointing=False, flagbackup=False) 
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["FIELD/table.dat",
                            "FIELD/EPHEM0_Mars_57034.9.tab",
                            "FIELD/EPHEM1_Titania_57034.9.tab"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All relevant tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname


        print(myname, " :  testing FIELD values in ms.summary()")
        try:
            mssum = mslocal.summary()
            # only Mars appears here because this short SDM only contains a single scan and that uses Mars
            self.assertEqual(mssum['scan_1']['0']['FieldName'],'Mars')
            self.assertAlmostEqual(mssum['field_0']['direction']['m0']['value'],-0.4770797859505159,15)
            self.assertAlmostEqual(mssum['field_0']['direction']['m1']['value'],-0.2154815444753364,15)
        except:
            print(myname, ": Error ms summary has an unexpect source or direction value")
            retValue['success']=False
            retValue['error_msg']=retValue['err_msg']+'Unexpected source or direction value in ms summary '+thismsname + '\n'


        mslocal.close()

        ephems = []
        # values from indivual rows were picked for no particular reason except verify they've not changed
        ephems.append({'name':"FIELD/EPHEM0_Mars_57034.9.tab",
                       'nrows':27,
                       'rows':[{'row':10,
                                'values':{'MJD':57035.041666666664,
                                          'RA':332.7140437500001,
                                          'DEC':-12.327346944444447,
                                          'Rho':2.024609480125507,
                                          'RadVel':723729.77502873}},
                               {'row':22,
                                'values':{'MJD':57035.208333333336,
                                          'RA':332.8387870833333,
                                          'DEC':-12.2793975,
                                          'Rho':2.0254053468626436,
                                          'RadVel':705588.202526264}}
                               ]
                       }
                      )

        ephems.append({'name':"FIELD/EPHEM1_Titania_57034.9.tab",
                       'nrows':45,
                       'rows':[{'row':17,
                                'values':{'MJD':57035.055555555555,
                                          'RA':11.813166666666666,
                                          'DEC':4.365749999999999,
                                          'Rho':20.150883673698488,
                                          'RadVel':2730048.0839084117}},
                               {'row':40,
                                'values':{'MJD':57035.21527777778,
                                          'RA':11.816041666666667,
                                          'DEC':4.3661111111111115,
                                          'Rho':20.153736461701364,
                                          'RadVel':2711142.1699538543}}
                               ]
                       }
                      )

        for ephem in ephems:
            print(myname,": Testing various things in ephemeris ", ephem['name'], " ...")

            tblocal.open(themsname+"/"+ephem['name'])
            kw = tblocal.getkeywords()
            nrows = tblocal.nrows()
            if not nrows==ephem['nrows']:
                print(myname,": Error. unexpected number of rows in ephemeris :",ephem['name'])
                retValue['success']=False
                retValue['error_msg']=retValue['error_msgs']+' Unexpected number of rows in ephemeris table :'+ ephem['name'] + '\n'

            for row in ephem['rows']:
                thisRow = row['row']
                for colname in row['values']:
                    thisVal = tblocal.getcell(colname,thisRow)
                    self.assertAlmostEqual(thisVal,row['values'][colname],10)

            # unfilled columns
            self.assertEqual((tblocal.getcol('diskLong') != 0.0).sum(),0)
            self.assertEqual((tblocal.getcol('diskLat') != 0.0).sum(),0)

            tblocal.close()
            geodist = kw['GeoDist'] # (km above reference ellipsoid)
            geolat = kw['GeoLat'] # (deg)
            geolong = kw['GeoLong'] # (deg)
            if not (geodist==geolat==geolong==0.):
                print(myname, ": ERROR.")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' Ephemeris was not converted to GEO for '+themsname+'\n'
            prsys = kw['posrefsys']
            if not (prsys=="ICRF/ICRS"):
                print(myname, ": ERROR.")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' posrefsys keyword is not ICRF/ICRS '+themsname+'\n'

        # fill and request an interpolated table.  Tests asdm2MS directly as this option isn't 
        # available in importasdm

        print(myname," filling an interpolated version of the same ephemeris")
        themsname_interp = myasdmname+".interp.ms"
        execute_string = "asdm2MS --no-pointing --interpolate-ephemeris 'yes' " + myasdmname + ' ' + themsname_interp
        print(myname, ' executing : ', execute_string)
        exitcode = os.system(execute_string)
        self.assertEqual(exitcode,0)
        ce.convert2geo(themsname_interp, '*') # convert the ephemeris to GEO
        # note that the recalculation of UVW and the adjustment of the SOURCE table are not
        # done here the way they would be done if filled via importasdm
        print(myname, ": Success! Now checking output ...")
        for name in mscomponents:
            if not os.access(themsname_interp+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname_interp+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname_interp+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All relevant tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname_interp)
        except:
            print(myname, ": Error  Cannot open MS table", themsname_interp)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname_interp
        print(myname, " :  testing FIELD values in ms.summary()")
        try:
            mssum = mslocal.summary()
            # only Mars appears here because this short SDM only contains a single scan and that uses Mars
            self.assertEqual(mssum['scan_1']['0']['FieldName'],'Mars')
            # difference here is < 0".0004 of the above, non-interpolated value
            self.assertAlmostEqual(mssum['field_0']['direction']['m0']['value'],-0.4770797877079177,15)
            # difference here is < 0".00005 of the above, non-interpolated value
            self.assertAlmostEqual(mssum['field_0']['direction']['m1']['value'],-0.2154815442529733,15)
        except:
            print(myname, ": Error ms summary has an unexpect source or direction value")
            retValue['success']=False
            retValue['error_msg']=retValue['err_msg']+'Unexpected source or direction value in ms summary '+thismsname + '\n'

        mslocal.close()
        ephems = []
        # values from indivual rows were picked for no particular reason except verify they've not changed
        # these rows 
        ephems.append({'name':"FIELD/EPHEM0_Mars_57034.9.tab",
                       'nrows':361,
                       'rows':[{'row':105,
                                'values':{'MJD':57035.008000000001630,
                                          'RA':332.688983703339886,
                                          'DEC':-12.337033046664128,
                                          'Rho':2.024447666954067,
                                          'RadVel':722270.482337458524853}},
                               {'row':320,
                                'values':{'MJD':57035.222999999998137,
                                          'RA':332.849807533339913,
                                          'DEC':-12.275182948886352,
                                          'Rho':2.025474163348287,
                                          'RadVel':702390.653405150165781}}
                               ]
                       }
                      )

        ephems.append({'name':"FIELD/EPHEM1_Titania_57034.9.tab",
                       'nrows':306,
                       'rows':[{'row':95,
                                'values':{'MJD':57035.033000000003085,
                                          'RA':11.812802333333494,
                                          'DEC':4.365715333333369,
                                          'Rho':20.150479182212013,
                                          'RadVel':2725243.279430569149554}},
                               {'row':250,
                                'values':{'MJD':57035.188000000001921,
                                          'RA':11.815509000000159,
                                          'DEC':4.366057555555590,
                                          'Rho':20.153251583732832,
                                          'RadVel':2721431.250284913461655}}
                               ]
                       }
                      )

        for ephem in ephems:
            print(myname,": Testing various things in ephemeris ", ephem['name'], " ...")

            tblocal.open(themsname_interp+"/"+ephem['name'])
            kw = tblocal.getkeywords()
            nrows = tblocal.nrows()
            if not nrows==ephem['nrows']:
                print(myname,": Error. unexpected number of rows in ephemeris :",ephem['name'])
                retValue['success']=False
                retValue['error_msg']=retValue['error_msgs']+' Unexpected number of rows in ephemeris table :'+ ephem['name'] + '\n'

            for row in ephem['rows']:
                thisRow = row['row']
                for colname in row['values']:
                    thisVal = tblocal.getcell(colname,thisRow)
                    self.assertAlmostEqual(thisVal,row['values'][colname],10)

            # unfilled columns
            self.assertEqual((tblocal.getcol('diskLong') != 0.0).sum(),0)
            self.assertEqual((tblocal.getcol('diskLat') != 0.0).sum(),0)

            tblocal.close()
            geodist = kw['GeoDist'] # (km above reference ellipsoid)
            geolat = kw['GeoLat'] # (deg)
            geolong = kw['GeoLong'] # (deg)
            if not (geodist==geolat==geolong==0.):
                print(myname, ": ERROR.")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' Ephemeris was not converted to GEO for '+themsname_interp+'\n'
            prsys = kw['posrefsys']
            if not (prsys=="ICRF/ICRS"):
                print(myname, ": ERROR.")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+' posrefsys keyword is not ICRF/ICRS '+themsname_interp+'\n'

        self.assertTrue(retValue['success'],retValue['error_msgs'])
        print(myname, ": OK.")


    def test7_lazy4(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler in lazy mode selecting only AUTO data, writing to FLOAT_DATA'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, ocorr_mode="ao", lazy=True, scans='0:1~4') # only the first 4 scans to save time
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
            mslocal.close()
            print(myname, ": MS can be opened. Now testing the changing of the asdmref ...")
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)
            
            mslocal.open(themsname)
            
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm="moved_"+myasdmname, vis='reference.ms', ocorr_mode="ao", lazy=False, overwrite=True, scans='0:1~3')

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.FLOAT_DATA,t2.FLOAT_DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")
                    
                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)

        os.system("mv moved_"+myasdmname+" "+myasdmname)
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test7_lazy5(self):
        '''Asdm-import: Test TP asdm with default filler in lazy mode selecting only AUTO data, writing to FLOAT_DATA'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X6218fb_X264'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, ocorr_mode="ao", bdfflags=True, applyflags=True, lazy=True)
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["ANTENNA/table.dat",
                            "CALDEVICE/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "SYSPOWER/table.dat",
                            "WEATHER/table.dat",
                            "ANTENNA/table.f0",
                            "CALDEVICE/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0",
                            "SYSPOWER/table.f0",
                            "WEATHER/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
            mslocal.close()
            print(myname, ": MS can be opened. Now testing the changing of the asdmref ...")
            mslocal.open(themsname)
            mslocal.asdmref("./moved_"+myasdmname)
            mslocal.close()
            os.system("mv "+myasdmname+" moved_"+myasdmname)
            
            mslocal.open(themsname)
            
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm="moved_"+myasdmname, vis='reference.ms', ocorr_mode="ao", lazy=False, bdfflags=True, applyflags=True, overwrite=True)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.FLOAT_DATA,t2.FLOAT_DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")

                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.WEIGHT,t2.WEIGHT, 1.e-06)))") == 0
                if not retValueTmp:
                    print("ERROR: WEIGHT does not agree with reference.")
                else:
                    print("WEIGHT columns agree.")
                    
                retValueTmp2 = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") == 0
                if not retValueTmp2:
                    print("ERROR: FLAG does not agree with reference.")
                else:
                    print("FLAG columns agree.")

                retValue['success'] = retValue['success'] and retValueTmp and retValueTmp2

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL",
                                 "WEATHER"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)

        os.system("mv moved_"+myasdmname+" "+myasdmname)
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test7_skiprows1(self):
        '''Asdm-import: Test TP asdm, comparing output when duplicate DATA rows are skipped versus not-skipped, lazy and regular, with bdflagging on'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X6218fb_X264'
        themsname = myasdmname+".ms"

        # the same tests are done for lazy being False and True
        for lazy in [False, True]:
            # always start with a clean slate
            shutil.rmtree(themsname,True)
            shutil.rmtree('referemce.ms',True)
            # use importasdm, which always looks for and skips duplicate DATA rows.
            self.res = importasdm(myasdmname, vis=themsname, ocorr_mode="ao", bdfflags=True, lazy=lazy, overwrite=True)
            self.assertEqual(self.res, None)
            print(myname,": Success! Now checking output ...")
            mscomponents = set(["ANTENNA/table.dat",
                                "CALDEVICE/table.dat",
                                "DATA_DESCRIPTION/table.dat",
                                "FEED/table.dat",
                                "FIELD/table.dat",
                                "FLAG_CMD/table.dat",
                                "HISTORY/table.dat",
                                "OBSERVATION/table.dat",
                                "POINTING/table.dat",
                                "POLARIZATION/table.dat",
                                "PROCESSOR/table.dat",
                                "SOURCE/table.dat",
                                "SPECTRAL_WINDOW/table.dat",
                                "STATE/table.dat",
                                "SYSCAL/table.dat",
                                "SYSPOWER/table.dat",
                                "WEATHER/table.dat",
                                "ANTENNA/table.f0",
                                "CALDEVICE/table.f0",
                                "DATA_DESCRIPTION/table.f0",
                                "FEED/table.f0",
                                "FIELD/table.f0",
                                "FLAG_CMD/table.f0",
                                "HISTORY/table.f0",
                                "OBSERVATION/table.f0",
                                "POINTING/table.f0",
                                "POLARIZATION/table.f0",
                                "PROCESSOR/table.f0",
                                "SOURCE/table.f0",
                                "SPECTRAL_WINDOW/table.f0",
                                "STATE/table.f0",
                                "SYSCAL/table.f0",
                                "SYSPOWER/table.f0",
                                "WEATHER/table.f0"
                                ])
            for name in mscomponents:
                if not os.access(themsname+"/"+name, os.F_OK):
                    print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
                else:
                    print(myname, ": ", name, "present.")
            print(myname, ": MS exists. All tables present. Try opening as MS ...")
            try:
                mslocal.open(themsname)
                print(myname, ": MS can be opened")
                mslocal.close()
            
            except:
                print(myname, ": Error  Cannot open MS table", themsname)
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
            else:
                print(myname, ": OK. Generating a reference MS with first integration checking turned off")

                # this must be done using asdm2MS and bdflags2MS directly
                asdm2MScmd = 'asdm2MS --ocm "ao" --checkdupints false'
                if lazy:
                    asdm2MScmd = asdm2MScmd + " --lazy"
                asdm2MScmd = asdm2MScmd + " " + myasdmname + " reference.ms"
                print(myname,'Running asdm2MS standalone invoked as:')
                print(asdm2MScmd)
                exitcode = os.system(asdm2MScmd)
                if exitcode != 0:
                    print(myname,"asdm2MS terminated with exit code ",exitcode)
                    retValue['success'] = False
                    retValue['error_msgs']=retValue['error_msgs']+' standalone execution of asdm2MS failed'
                    # this should break out of the main loop over lazy values
                    break
                
                bdflags2MScmd = 'bdflags2MS -f ALL --ocm "ao" --checkdupints false'
                if lazy:
                    bdflags2MScmd = bdflags2MScmd + " --lazy=true"
                bdflags2MScmd = bdflags2MScmd + " " + myasdmname + " reference.ms"
                print(myname,'Running bdflags2MS standalone invoked as:')
                print(bdflags2MScmd)
                exitcode = os.system(bdflags2MScmd)
                if exitcode != 0:
                    print(myname,"bdflags2MS terminated with exit code ",exitcode)
                    retValue['success'] = False
                    retValue['error_msgs']=retValue['error_msgs']+' standalone execution of bdflags2MS failed'
                    # this should break out of the main loop over lazy values
                    break                    

                # at this point, reference.ms should exist, with all of the auto rows (ao), including duplicates, with BDF flags applied
                if(os.path.exists('reference.ms')):
                    # and the test here is to make sure that the two MSs are identical in the Main table. The same values in 
                    # the same order except for gaps where rows from the larger table are skipped. 

                    # expected size (rows)
                    msSize = 31972
                    refSize = 32040
                    
                    # expected gaps start at these rows and are always 4 rows long, row numbers in reference.ms
                    gaps = [2280,3176,5328,7120,9276,10172,11432,12328,14480,15376,18752,19648,21800,26636,27896,28792,30944]

                    mstb = tbtool()
                    mstb.open(themsname)
                    if mstb.nrows() != msSize:
                        print(myname,'MS size is not of the expected number of rows : ',mstb.nrows(),' != ',msSize)
                        retValue['success'] = False
                        retValue['error_msgs'] = retValue['error_msgs'] + 'bad size for MS'
                    
                    reftb = tbtool()
                    reftb.open('reference.ms')
                    if reftb.nrows() != refSize:
                        print(myname,'Reference MS size is not of the expected number of rows : ',reftb.nrows(),' != ',refSize)
                        retValue['success'] = False
                        retValue['error_msgs'] = retValue['error_msgs'] + 'bad size for reference MS'

                    if retValue['success']:
                        refrow = 0
                        msrow = 0
                        gapStart = -1
                        gapIndex = -1
                        # assumes set of column names is the same
                        cols = mstb.colnames()

                        while((msrow < mstb.nrows()) and (refrow < reftb.nrows())):
                            matchFound = True
                            for col in cols:
                                if mstb.iscelldefined(col,msrow) != reftb.iscelldefined(col,refrow):
                                    # defined in one cell, not defined in the other
                                    matchFound = False
                                else:
                                    if (mstb.iscelldefined(col,msrow)):
                                        msVal = mstb.getcell(col,msrow)
                                        refVal = reftb.getcell(col,refrow)
                                        # assumes the type is the same for the same column in both tables
                                        if isinstance(msVal,numpy.ndarray):
                                            matchFound = numpy.array_equal(msVal,refVal)
                                        else:
                                            matchFound = msVal == refVal
                                if not matchFound:
                                    # no point in checking the other columns
                                    break
                            if matchFound:
                                if gapStart >= 0:
                                    # a gap has ended, verify that it was exactly 4 rows long
                                    if (refrow-gapStart) != 4:
                                        print(myname,'Unexpected gap length not equal to 4 rows. Gap length = ',(refrow-gapStart),' starting at row ',refRow)
                                        retValue['success'] = False
                                        retValue['error_msg'] = 'Unexpected gap length not equal to 4 rows'
                                    gapStart = -1
                                msrow = msrow + 1
                            else:
                                # do not increment msrow in this case, keep looking for a match
                                if gapStart < 0:
                                    # new gap, increment index and verify it's at the expected place
                                    gapIndex = gapIndex+1
                                    gapStart = refrow
                                    if gapIndex > len(gaps):
                                        print(myname,'Unexpected gap seen past end of known gaps. Starting at row ',refrow)
                                        retValue['success'] = False
                                        retValue['error_msg'] = 'Unexpected gap after end of known gaps'
                                    else:
                                        if gapStart != gaps[gapIndex]:
                                            print(myname,'Unexpected gap start at row ',gapStart,' expected at row ',gaps[gapIndex])
                                            retValue['success'] = False
                                            retValue['error_msg'] = 'Unexpected gap start row'
                            # refrow is always incremented
                            refrow = refrow + 1

                            # bail out on failure
                            if not retValue['success']:
                                break
                    
                    reftb.close()
                    mstb.close()
            # break out of the lazy loop on failure
            if not retValue['success']:
                break
                    
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        
    def test7_bdflags1(self):
        '''Asdm-import: Test good 12 m ASDM with mixed pol/channelisation input with default filler selecting "co" on output and using the BDF flags'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        myasdmname = 'uid___A002_X71e4ae_X317_short'
        themsname = myasdmname+".ms"

        self.res = importasdm(myasdmname, vis=themsname, ocorr_mode="co", bdfflags=True) 
        self.assertEqual(self.res, None)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "FEED/table.dat",
                            "FIELD/table.dat",
                            "FLAG_CMD/table.dat",
                            "HISTORY/table.dat",
                            "OBSERVATION/table.dat",
                            "POINTING/table.dat",
                            "POLARIZATION/table.dat",
                            "PROCESSOR/table.dat",
                            "SOURCE/table.dat",
                            "SPECTRAL_WINDOW/table.dat",
                            "STATE/table.dat",
                            "SYSCAL/table.dat",
                            "ANTENNA/table.f0",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.f0",
                            "FIELD/table.f0",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.f0",
                            "POINTING/table.f0",
                            "POLARIZATION/table.f0",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.f0",
                            "SPECTRAL_WINDOW/table.f0",
                            "STATE/table.f0",
                            "SYSCAL/table.f0"
                            ])
        for name in mscomponents:
            if not os.access(themsname+"/"+name, os.F_OK):
                print(myname, ": Error  ", themsname+"/"+name, "doesn't exist ...")
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+themsname+'/'+name+' does not exist'
            else:
                print(myname, ": ", name, "present.")
        print(myname, ": MS exists. All tables present. Try opening as MS ...")
        try:
            mslocal.open(themsname)
        except:
            print(myname, ": Error  Cannot open MS table", themsname)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+themsname
        else:
            mslocal.close()
            print(myname, ": OK. Checking tables in detail ...")
    
            importasdm(asdm=myasdmname, vis='reference.ms', overwrite=True, ocorr_mode="co", bdfflags=False)

            if(os.path.exists('reference.ms')):
                retValue['success'] = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                                    +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1.DATA,t2.DATA, 1.e-06)))") == 0
                if not retValue['success']:
                    print("ERROR: DATA does not agree with reference.")
                else:
                    print("DATA columns agree.")
                retValueTmp = th.checkwithtaql("select from [select from reference.ms orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "
                                            +themsname+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1.FLAG==t2.FLAG)) ") != 0
                if not retValueTmp:
                    print("ERROR: FLAG columns do agree with reference but they shouldn't.")
                else:
                    print("FLAG columns do not agree as expected.")

                retValue['success'] = retValue['success'] and retValueTmp

                for subtname in ["ANTENNA",
                                 "DATA_DESCRIPTION",
                                 "FEED",
                                 "FIELD",
                                 "FLAG_CMD",
                                 "OBSERVATION",
                                 "POLARIZATION",
                                 "PROCESSOR",
                                 "SOURCE",
                                 "SPECTRAL_WINDOW",
                                 "STATE",
                                 "SYSCAL"]:
                    
                    print("\n*** Subtable ",subtname)
                    excllist = []
                    if subtname=='SOURCE':
                        excllist=['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL']
                    if subtname=='SYSCAL':
                        excllist=['TANT_SPECTRUM', 'TANT_TSYS_SPECTRUM']
                    if subtname=='SPECTRAL_WINDOW':
                        excllist=['CHAN_FREQ', 'CHAN_WIDTH', 'EFFECTIVE_BW', 'RESOLUTION', 'ASSOC_SPW_ID', 'ASSOC_NATURE']
                        for colname in excllist:
                            if colname!='ASSOC_NATURE':
                                retValue['success'] = th.compVarColTables('reference.ms/SPECTRAL_WINDOW',
                                                                          themsname+'/SPECTRAL_WINDOW', colname, 0.01) and retValue['success']
                    if subtname=='POLARIZATION':
                        excllist=['CORR_TYPE', 'CORR_PRODUCT']
                        for colname in excllist: 
                            retValue['success'] = th.compVarColTables('reference.ms/POLARIZATION',
                                                                      themsname+'/POLARIZATION', colname, 0.01) and retValue['success']
                    try:    
                        retValue['success'] = th.compTables('reference.ms/'+subtname,
                                                            themsname+'/'+subtname, excllist, 
                                                            0.01) and retValue['success']
                    except:
                        retValue['success'] = False
                        print("ERROR for table ", subtname)
            
                
        self.assertTrue(retValue['success'],retValue['error_msgs'])

class asdm_import8(test_base):
    # these are more like unit tests, difficult to test without invoking all of importasdm
    # currently this is just tests on SDM_NUM_BIN
    
    def setUp(self):
        self.setUp_numbin()

    def tearDown(self):
        pass
        #for this_asdm_name in ['alma_numbin_mixed','evla_numbin_2','evla_numbin_4']:
        #    os.system('rm -rf '+this_asdm_name+"*")

    def doNumTest(self, testName, asdm_name, ms_name, spWin_name, execBlock_name, expWinFunCol, expNumBinCol, expResCol):
        retValue = {'success': True, 'error_msgs': '' } 
        print(testName,": testing SDM columns in",asdm_name)

        originalSpWin = None
        originalExecBlock = None

        if spWin_name is not None:
            print(testName,": using",spWin_name,"for SpectralWindow.xml")
            spWin_path = asdm_name+"/"+spWin_name
            originalSpWin = asdm_name+"/SpectralWindow.xml.original"
            if not os.path.exists(spWin_path):
                msg = spWin_path+" does not exist"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            if os.path.exists(originalSpWin):
                msg = originalSpWin+" already exists, will not overwrite"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            shutil.move(asdm_name+'/SpectralWindow.xml',originalSpWin)
            shutil.copyfile(spWin_path,asdm_name+'/SpectralWindow.xml')

        if execBlock_name is not None:
            print(testName,": using",execBlock_name,"for ExecBlock.xml")
            execBlock_path = asdm_name+"/"+execBlock_name
            originalExecBlock = asdm_name+"/ExecBlock.xml.original"
            if not os.path.exists(execBlock_path):
                msg = execBlock_path+" does not exist"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            if os.path.exists(originalExecBlock):
                msg = originalExecBlock+" already exists, will not overwrite"
                retValue['success'] = False
                retValue['error_msgs'] = msg
                return retValue
            shutil.move(asdm_name+'/ExecBlock.xml',originalExecBlock)
            shutil.copyfile(execBlock_path,asdm_name+'/ExecBlock.xml')
 
        importasdm(asdm=asdm_name,vis=ms_name,lazy=True,process_syspower=False,process_caldevice=False,process_pointing=False,process_flags=False)

        # the only table this test cares about is SPECTRAL_WINDOW
        spwName = ms_name + "/SPECTRAL_WINDOW"
        if not os.access(spwName,os.F_OK):
            print(testName,": Error ", spwName, "doesn't exist ...")
            retValue['success'] = False
            retValue['error_msgs']=spwName+' does not exist'
        else:
            ok = tblocal.open(spwName)
            if (ok):
                try:
                    winFunCol = tblocal.getcol('SDM_WINDOW_FUNCTION')
                    if not numpy.all(winFunCol==expWinFunCol):
                        retValue['success'] = False
                        msg = "ERROR Unexpected SDM_WINDOW_FUNCTION values when filling "+asdm_name
                        retValue['error_msgs']=msg
                        print(testName,":", msg)  
                except:
                     retValue['success'] = False
                     msg = "ERROR getting/testing SDM_WINDOW_FUNCTION column in "+spwName
                     retValue['error_msgs']=msg
                     print(testName,":",msg)

                try:
                    numBinCol = tblocal.getcol('SDM_NUM_BIN')
                    if not numpy.all(numBinCol==expNumBinCol):
                        retValue['success'] = False
                        msg = "ERROR Unexpected SDM_NUM_BIN values when filling "+asdm_name
                        # there may already be messages in error_msgs
                        if len(retValue['error_msgs']>0):
                            retValue['error_msgs']=retValue['error_msgs']+'\n'
                        retValue['error_msgs']=retValue['error_msgs']+msg
                        print(testName,":",msg)             
                except:
                    retValue['success'] = False
                    msg = "ERROR getting/testing SDM_NUM_BIN column in "+spwName
                    # there may already be messages in error_msgs
                    if len(retValue['error_msgs'])>0:
                        retValue['error_msgs']=retValue['error_msgs']+'\n'
                    retValue['error_msgs']=retValue['error_msgs']+msg
                    print(testName,":",msg)

                # only test RESOLUTION values if expResCol is not None
                if expResCol is not None:
                    try:
                        resCol = tblocal.getcol('RESOLUTION')
                        # only test first value in each row, assumes all values in a row are equal
                        resCol = resCol[0,:]
                        if not numpy.all(resCol==expResCol):
                            retValue['success'] = False
                            msg = "ERROR Unexpected RESOLUTION values when filling "+asdm_name
                            # there may already be messages in error_msgs
                            if len(retValue['error_msgs']>0):
                                retValue['error_msgs']=retValue['error_msgs']+'\n'
                            retValue['error_msgs']=retValue['error_msgs']+msg
                            print(testName,":",msg)             
                    except:
                        retValue['success'] = False
                        msg = "ERROR getting/testing RESOLUTION column in "+spwName
                        # there may already be messages in error_msgs
                        if len(retValue['error_msgs'])>0:
                            retValue['error_msgs']=retValue['error_msgs']+'\n'
                        retValue['error_msgs']=retValue['error_msgs']+msg
                        print(testName,":",msg)
                tblocal.close()

            else:
                msg = "ERROR opening",spwName
                retValue['success'] = False
                retValue['error_msgs'] = msg
                print(testName,":",msg)

        if originalSpWin is not None:
            os.remove(asdm_name+'/SpectralWindow.xml')
            shutil.move(originalSpWin,asdm_name+'/SpectralWindow.xml')
            print(testName,": restored original SpectralWindow.xml")

        if originalExecBlock is not None:
            os.remove(asdm_name+'/ExecBlock.xml')
            shutil.move(originalExecBlock,asdm_name+'/ExecBlock.xml')
            print(testName,": restored original ExecBlock.xml")

        return retValue 

    def test_alma_numbin(self):
        retValue = {'success': True, 'error_msgs': '' } 

        # original SpectralWindow.xml and inferred numBin  values
        asdm_name = 'alma_numbin_mixed'
        ms_name = asdm_name+".ms"
        expWinFunCol = numpy.array(['UNIFORM']*5 + ['HANNING']*8 + ['UNIFORM']*4 + ['HANNING']*18 + ['UNIFORM']*42)
        # expected values, 8 @ 24, 2 @ 27,29,31,33, rest are 1
        expNumBinCol = numpy.ones(77,dtype=numpy.int32)
        expNumBinCol[25] = 8
        for indx in [27,29,31,33]:
            expNumBinCol[indx] = 2
        res = self.doNumTest(myname,asdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with appropriate numBin values, should yield same column values
        ms_name = asdm_name+".numbin.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        # SpectralWindow.xml with faked resolution and expectedBw values but no numBin, tests other inferred numBin values
        expNumBinCol[5] = 4
        expNumBinCol[7] = 16
        ms_name = asdm_name+".faked.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.faked',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        # SpectralWindow.xml with faked resolution and expectedBw values and added numBin values, same expected values as previous test
        ms_name = asdm_name+".faked.numBin.ms"
        res = self.doNumTest(myname,asdm_name,ms_name,'SpectralWindow.xml.faked.numBin',None,expWinFunCol,expNumBinCol,None)
        retValue['success'] = retValue['success'] and res['success']
        retValue['error_msgs'] = retValue['error_msgs'] + res['error_msgs']

        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test_evla_numbin(self):
        retValue = {'success':True, 'error_msgs':''}

        # numbin=2 related tests
        sdm_name = 'evla_numbin_2'

        # original SpectralWindow.xml and inferred numBin values, all equal to 2
        expWinFunCol = numpy.array(['UNIFORM']*16)
        expNumBinCol = numpy.empty(16,dtype=numpy.int32)
        expNumBinCol.fill(2)
        # also should alter resolution to these expected values
        expResCol = numpy.empty(16)
        expResCol.fill(4000000.)
        ms_name = sdm_name+".ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with numBin field and appropriately modified resolution, same expected values
        ms_name = sdm_name+".numBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # SpectralWindow.xml with mostly numBin and alterned resolution, but one row has the original values, same expected values
        ms_name = sdm_name+".mixed.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.mixed',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml but with one non-physical (bad) value of resolution leading to the algorithm giving up and numBin=1 and resolution the original bad value
        ms_name = sdm_name+".bad.ms"
        expNumBinCol[0] = 1
        expResCol[0] = 9000000.
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.bad',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # numbin=4 related tests
        sdm_name = 'evla_numbin_4'
        
        # original SpectralWindow.xml and inferred numBin values, all equal to 4
        expWinFunCol = numpy.array(['UNIFORM']*16)
        expNumBinCol = numpy.empty(16,dtype=numpy.int32)
        expNumBinCol.fill(4)
        # also should alter resolution to these expected values
        expResCol = numpy.empty(16)
        expResCol.fill(8000000.)
        ms_name = sdm_name+".ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml with numBin field and altered resolution, same expected values
        ms_name = sdm_name+".numBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.numBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml with numBin field and original resolution
        # expected numBin is the same, expected resolution is now the original values
        expResCol /= 4.0
        ms_name = sdm_name+".onlyNumBin.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,'SpectralWindow.xml.onlyNumBin',None,expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        # original SpectralWindow.xml and altered ExecBlock so that the telescope is UNKNOWN
        # numBin is all 1 and expected resolution is the original resolution
        expNumBinCol.fill(1)
        ms_name = sdm_name+".unknownTel.ms"
        res = self.doNumTest(myname,sdm_name,ms_name,None,'ExecBlock.xml.unknownTel',expWinFunCol,expNumBinCol,expResCol)
        retValue['success'] = res['success']
        retValue['error_msgs'] = res['error_msgs']

        self.assertTrue(retValue['success'],retValue['error_msgs'])


def suite():
    return [asdm_import1, 
            asdm_import2, 
            asdm_import3, 
            asdm_import4,
            asdm_import5,
            asdm_import6,
            asdm_import7,
            asdm_import8]
        
    
