#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for the ATCA RPFITS data import to MS           #
#    
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The conversion of RPFITS data to MS needs to be verified.              #
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
import numpy

myname = 'importatca-unit-test'

# default dataset name
my_dataset_names = ['2012-10-25_0707-002.C2728']

# name of the resulting MS
msname = my_dataset_names[0].split('.')[1].lower()+'.ms'

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

class test_importatca(unittest.TestCase):
    
    def setUp(self):
        res = None

        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/importatca/'
        for fname in my_dataset_names:
            if(os.path.exists(fname)):
                os.remove(fname)
            shutil.copy(datapath + fname, fname)
        default(importatca)
        
    def tearDown(self):
        for fname in my_dataset_names:
            os.remove(fname)
        #shutil.rmtree(msname,ignore_errors=True)
        shutil.rmtree(msname+'.flagversions',ignore_errors=True)
        
    def test1(self):
        '''atca-rpfits-import: Test good input'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }    

        self.res = importatca(my_dataset_names[0], msname, options="noac")
        print(myname, ": Success! Now checking output ...")
        mscomponents = set(["table.dat",
                            "table.f0",
			    "table.f0i",
                            "table.f1",
                            "table.f2",
                            "table.f2_TSM1",
                            "table.f2_TSM2",
                            "table.f2_TSM3",
                            "table.f3",
                            "table.f3_TSM1",
                            "table.f3_TSM2",
                            "table.f3_TSM3",
                            "table.f4",
                            "table.f4_TSM1",
                            "table.f4_TSM2",
                            "table.f4_TSM3",
                            "table.f5",
                            "table.f5_TSM1",
                            "table.f5_TSM2",
                            "table.f5_TSM3",
                            "ANTENNA/table.dat",
                            "ANTENNA/table.f0",
                            "ATCA_SCAN_INFO/table.dat",
                            "ATCA_SCAN_INFO/table.f0",
                            "ATCA_SCAN_INFO/table.f1",
                            "DATA_DESCRIPTION/table.dat",
                            "DATA_DESCRIPTION/table.f0",
                            "FEED/table.dat",
                            "FEED/table.f0",
                            "FEED/table.f0i",
                            "FIELD/table.dat",
                            "FIELD/table.f0",
                            "FIELD/table.f0i",
                            "FLAG_CMD/table.dat",
                            "FLAG_CMD/table.f0",
                            "HISTORY/table.dat",
                            "HISTORY/table.f0",
                            "OBSERVATION/table.dat",
                            "OBSERVATION/table.f0",
                            "POINTING/table.dat",
                            "POINTING/table.f0",
                            "POINTING/table.f0i",
                            "POINTING/table.f1",
                            "POLARIZATION/table.dat",
                            "POLARIZATION/table.f0",
                            "POLARIZATION/table.f0i",
                            "PROCESSOR/table.dat",
                            "PROCESSOR/table.f0",
                            "SOURCE/table.dat",
                            "SOURCE/table.f0",
                            "SOURCE/table.f0i",
                            "SPECTRAL_WINDOW/table.dat",
                            "SPECTRAL_WINDOW/table.f0",
                            "SPECTRAL_WINDOW/table.f0i",
                            "STATE/table.dat",
                            "STATE/table.f0",
                            "SYSCAL/table.dat",
                            "SYSCAL/table.f0",
                            "SYSCAL/table.f1",
                            "SYSCAL/table.f1_TSM0"
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
                         ['UVW',       0, [ 37.08396912, 110.68421173, 121.33368683], 1E-8],
                         ['EXPOSURE',  0, 9.856, 1E-4],
                         ['DATA',      0,[[-0.19179995-0.27503902j],
					  [-0.07504284-0.07763118j],
					  [-0.13403413+0.04170537j],
					  [ 0.24940215-0.23159257j]], 1E-8]
                         ]
	    dataslice=[[0,1025],[3,1025],[1,1]]
            results = checktable(name, expected, dataslice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
    
            expected = [
                         ['UVW',       854, [870.02380371, 2447.86230469, 2683.08862305], 1E-8],
                         ['EXPOSURE',  854, 9.856, 1E-4],
                         ['DATA',      854,[[ 0.69642496-0.57600111j],
					    [ 0.01779159+0.27384654j],
					    [-0.07370736-0.30479908j],
					    [ 0.45262983+0.29306653j]], 1E-8]
                         ]
            results = checktable(name, expected, dataslice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "ANTENNA"
            expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                         ['POSITION',     1, [-4751591.759, 2791758.855, -3200483.75], 0.001],
                         ['DISH_DIAMETER', 1, 22.0, 0.0] 
                         ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            
            name = "SPECTRAL_WINDOW"
            expected = [ ['NUM_CHAN',        0, 2049, 0],
                         ['TOTAL_BANDWIDTH', 0, 2048e6, 0.01e6],
                         ['CHAN_WIDTH',      0, [-1e6, -1e6, -1e6], 0.1],
                         ['CHAN_FREQ',       0, [ 3124e6,   3123e6,  3122e6], 0.01E6],
                         ['CHAN_FREQ',       2, [ 1411e6, 1410.96875e6, 1410.9375e6], 10.]
                         ]
	    freqSlice=[[0],[2],[1]]
            results = checktable(name, expected,freqSlice)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
                
        self.assertTrue(retValue['success'])
        
        # test atm pressure
        self._check_atm_pressure(msname)
                
    def _check_atm_pressure(self, vis):
        weather_table = os.path.join(vis, 'WEATHER')
        tb.open(weather_table)
        try:
            # PRESSURE column should exist
            self.assertTrue('PRESSURE' in tb.colnames())
            
            # unit should be hPa
            colkeys = tb.getcolkeywords('PRESSURE')
            self.assertTrue('QuantumUnits' in colkeys)
            pressure_unit = colkeys['QuantumUnits'][0]
            print(('Pressure unit is {0}'.format(pressure_unit)))
            self.assertEqual(pressure_unit, 'hPa')
            
            # value should be in reasonable range
            pressure_min = 400.0
            pressure_max = 1100.0
            pressure_value = tb.getcol('PRESSURE')
            self.assertTrue(numpy.all(pressure_min < pressure_value))
            self.assertTrue(numpy.all(pressure_value < pressure_max))
        finally:
            tb.close()
    
def suite():
    return [test_importatca]
