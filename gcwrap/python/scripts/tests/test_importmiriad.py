#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for the Miriad uv data import to MS             #
#    
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The conversion of Miriad uv data to MS needs to be verified.           #
#                                                                           # 
# Features tested:                                                          #
#    1) Is the import performed without raising exceptions                  #
#    2) Do all expected tables exist                                        #
#    3) Can the MS be opened                                                #
#    4) Do the tables contain expected values                               #
#                                                                           #
# Input data:                                                               #
#                                                                           #
#############################################################################
import os
import sys
import shutil
from __main__ import default
from tasks import *
from taskinit import *
import unittest

myname = 'importmiriad-unit-test'

# default dataset name
my_dataset_names = ['1934.uv']

# name of the resulting MS
msname = my_dataset_names[0].split('.')[0]+'.ms'

def checktable(thename, theexpectation, dataslice=[]):
    global msname, myname
    tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print(myname, ": comparing ", mycell)
        if mycell[0]=="DATA" or mycell[0]=="CHAN_WIDTH" or mycell[0]=="CHAN_FREQ":
	    value = tb.getcellslice(mycell[0], mycell[1],dataslice[0],dataslice[1],dataslice[2])
	else:
	    value = tb.getcell(mycell[0], mycell[1])
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
            if isinstance(value, str):
                in_agreement = value == mycell[2]
            else:
                # it's an array
                # zero tolerance?
                if mycell[3] == 0:
                    in_agreement =  (value == mycell[2]).all() 
                else:
                    try:
                        in_agreement = (abs(value - mycell[2]) < mycell[3]).all()
                    except:
                        in_agreement = False
        if not in_agreement:
            print(myname, ":  Error in MS subtable", thename, ":")
            print("     column ", mycell[0], " row ", mycell[1], " contains ", value)
            print("     expected value is ", mycell[2])
            tb.close()
            return False
    tb.close()
    print(myname, ": table ", thename, " as expected.")
    return True


###########################
# beginning of actual test 

class test_importmiriad(unittest.TestCase):
    
    def setUp(self):
        res = None

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/importmiriad/'
        for fname in my_dataset_names:
            if(os.path.exists(fname)):
                shutil.rmtree(fname)
            shutil.copytree(datapath + fname, fname)
        default(importmiriad)
        
    def tearDown(self):
        for fname in my_dataset_names:
            shutil.rmtree(fname)
        shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        
    def test1(self):
        '''miriad-import: Test good input'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importmiriad(my_dataset_names[0], msname,tsys=False)
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
                            "table.f0",
			    "table.f0i",
                            "table.f1",
                            "table.f1_TSM1",
                            "table.f2",
                            "table.f2_TSM1",
                            "table.f3",
                            "table.f3_TSM0",
                            "ANTENNA/table.dat",
                            "DATA_DESCRIPTION/table.dat",
                            "DOPPLER/table.dat",
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
                            "DOPPLER/table.f0",
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
                            "FEED/table.f0i",
                            "FIELD/table.f0i",
                            "POINTING/table.f0i",
                            "POLARIZATION/table.f0i",
                            "SOURCE/table.f0i",
                            "SPECTRAL_WINDOW/table.f0i",
                            "SYSCAL/table.f0i"
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
            ms.open(msname)
        except:
            print(myname, ": Error  Cannot open MS table", tablename)
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            ms.close()
            print(myname, ": OK. Checking tables in detail ...")
            retValue['success']=True
    
            # check main table first
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                         ['UVW',       0, [ 167.85437179, 11.69765597, 6.01248136], 1E-8],
                         ['EXPOSURE',  0, 9.856, 1E-4],
                         ['DATA',      0,[[ 40.19867706 -8.43087101j],
      					  [ -0.65212655 +0.32635808j],
      					  [  0.77692640 -0.51819152j],
      					  [ 47.81456757+12.74779224j]], 1E-8]
                         ]
	    dataslice=[[0,1025],[3,1025],[1,1]]
            results = checktable(name, expected, dataslice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            expected = [
                         ['UVW',       209, [3725.85784932, 228.183798, 117.68756041], 1E-8],
                         ['EXPOSURE',  209, 9.856, 1E-4],
                         ['DATA',      209,[[ 54.56027985-2.18005967j],
       				            [ -0.06506564-0.72650862j],
      				            [  0.65535468+0.22498345j],
       					    [ 45.48416138+1.93808842j]], 1E-8]
                         ]
            results = checktable(name, expected, dataslice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [-4751589.52238021,  2791757.53976021, -3200482.25099623], 0.0001],
                         ['DISH_DIAMETER', 1, 22.0, 0.0] 
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "SPECTRAL_WINDOW"
            expected = [ ['NUM_CHAN',        0, 2049, 0],
                         ['TOTAL_BANDWIDTH', 0, 2049e6, 0.01e6],
                         ['CHAN_WIDTH',      0, [ -1e6, -1e6, -1e6], 0.1],
                         ['CHAN_FREQ',       0, [  3124e6,   3123e6,  3122e6], 0.01E6]
                         ]
	    freqSlice=[[0],[2],[1]]
            results = checktable(name, expected,freqSlice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])
                
    
def suite():
    return [test_importmiriad]
