//# PointingDirectionCalculator_GTest.cc: this defines unit tests of
//# PointingDirectionCalculator  using google test framework
//#
//# Copyright (C) 2018
//# National Astronomical Observatory of Japan
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

// G_Test Include //

#include "gtest/gtest.h"

// include file for TABLE Manupilation //

#include <casacore/casa/OS/Directory.h>
#include <casacore/casa/OS/Path.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <casacore/tables/Tables/ScaColData.h>

// org include //

#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <synthesis/Utilities/PointingDirectionCalculator.h>

#include <casa/aipstype.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/BasicSL/String.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Containers/Block.h>
#include <casa/Utilities/BinarySearch.h>
#include <casa/Logging/LogIO.h>
#include <tables/TaQL/ExprNode.h>
#include <ms/MSSel/MSSelection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>

using namespace casacore;
using namespace std;

//+
// Additional CASACORE include files
//-

#include <cstdio>
#include <casa/OS/EnvVar.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

namespace casa {

//******************************************************************************
//  Varisous Types of MaeasurementSet defnition
//
//  - Some MSs are not compatible with the required, which may cause Exception.
//  - Expected Exception are describged in the definition
//******************************************************************************


//+
// MeasurementSet NameList class
//-

class MSNameList {

public:
    typedef struct _MSDef
    {
        bool   ExThrow;  // True = cause Exeption
        String name;     // MS name, with relative path from CASAPATH
    } MSDef;

    // Get File Name by Number //
    const String  getName(uInt No ) {
        const String msg = "Internal Bugcheck:";
        if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());
        return TestMSList[No].name;
    }

    // True is this access makes Exception 
    bool isExceptionActivated(uInt No) { 
        const String msg = "Internal Bugcheck:";
        if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());
        return TestMSList[No].ExThrow;
   }

    uInt count() { return TestMSList.size();  }

private:

    std::vector<MSDef> TestMSList 
    {
        // Exeption(bool)  , Filename //
 
        {true, "./sdimaging-t.ms"   		},
        {false, "sdimaging/sdimaging.ms"        },
        {false, "listobs/uid___X02_X3d737_X1_01_small.ms" },

        // Following 2 MS are affected assert(), cannot run on UT
        // Release EXE must throws Excepton. 
#if 0
        {true,  "concat/input/A2256LC2_4.5s-1.ms"               },
        {true,  "concat/input/A2256LC2_4.5s-2.ms"               },
#endif 
        {false, "sdimaging/Uranus1.cal.Ant0.spw34.ms"   },
        {false, "sdimaging/Uranus2.cal.Ant0.spw34.ms"   },
        {false, "sdimaging/azelpointing.ms"      	},
        {false, "sdimaging/clipping_1row.ms"            },
        {false, "sdimaging/clipping_2rows.ms"           },
        {false, "sdimaging/clipping_3rows.ms"           },
        {false, "sdimaging/clipping_3rows_2chans.ms"    },
        {false, "sdimaging/clipping_3rows_suprious.ms"  },
        {false, "sdimaging/pointing6.ms"                },
        {false, "sdimaging/sdimaging_flagtest.ms"       },
        {false, "sdimaging/selection_intent.ms"         },
        {false, "sdimaging/selection_misc.ms"           },
        {false, "sdimaging/selection_spw.ms"            },
        {false, "sdimaging/selection_spw_unifreq.ms"    },
        {true, ".sis14_twhya_calibrated_flagged.ms"     },
        {true, ".sis14_twhya_calibrated_flagged-t.ms"   },
        {true,  "sdimaging/hogehoge.ms"                 },
        {true,  "sdimaging/hogehoge.ms"                 },  // Extra Hoge for self-degbug .
 
        // Any additional definition can be written here as you want. // 
    };

};

//+
// NEW: Under construction
// Direction Column (List and other service)
//-
class DirectionColumnList 
{
public:
        string dirName(uint n) { return DirList[n]; }
private:
        std::vector<string> DirList 
         = {"DIRECTION", "TARGET", "POINTING_OFFSET", "SOURCE_OFFSET", "ENCODER" };

};

//+
// DEBUG Tentative 
//-

bool use_spline = false;


//****************************************
// Execution / Running Enviromnent
// CASAPATH is the base directory
//****************************************

class RunEnv
{
public:

    RunEnv()
    {
        //+
        // Path and directory Position by Env. Variable.
        //-
        
        CasaPath        = GetCasaPath( "CASAPATH" );
        CasaMasterPath  = CasaPath + "/data/regression/unittest/"; 
        UnitTestTMasterFileName = CasaMasterPath + "sdimaging/sdimaging.ms";

        printf("RunEnv:: Environment Variable Information -----\n" );
        printf("CASAPATH      :%s \n", CasaPath.c_str());
        printf("CasaMasterPath:%s \n", CasaMasterPath.c_str());
    }


    const String getCasaPath()
    {
        return CasaPath;
    }

    const String getCasaMasterPath()
    {
        return CasaMasterPath;
    }
     

private:
    
    String GetCasaPath(const String pathname )
    {
        if (casacore::EnvironmentVariable::isDefined(pathname)) 
        {
            string casapath = casacore::EnvironmentVariable::get(pathname);
            size_t endindex = casapath.find(" ");
            if (endindex != string::npos)
            {
                string casaroot = casapath.substr(0, endindex);
                cout << pathname << "=" << casaroot << endl;
                return (casaroot);
            } 
            else 
            {
                cout << "hit npos" << endl;
                return "/hoge/";
            }
        } 
        else 
        {
            cout << "ERROR: Specified path " << pathname  << " is not defined" << endl;
            return "";
        }
    }

    String CasaPath;            // translated from CASAPATH 
    String CasaMasterPath;
    String UnitTestTMasterFileName;

};




//************************************** 
// Base TestClass
//**************************************

class BaseClass : public ::testing::Test
{

public:

     RunEnv       env;

     //+
     // MS for Test
     //-
     void CopyDefaultMStoWork();
     void DeleteWorkingMS();

     //+
     // Console Message of test progress
     //-

     void TestDescription( const String &Title );
     void FunctionalDescription(const String &Title, const String &Param);
     void Description(const String &Title, const String &Param);

    uInt    ExpectedNrow = 0;   // C++11 feature //

    BaseClass()  { }

    ~BaseClass() { }

    void SetUp() { }
    void TearDown() { }


private:

     // Test MS copy and delete flag //
     bool copyOperation    = true;
     bool deleteOperation  = false;
 
};

//*
// Global Constant
//*
const String  DefaultLocalMsName = "./sdimaging-t.ms";


//******************************************************
//  Log Title Output Functions for readable text 
//  of this UT.
//******************************************************

void BaseClass::TestDescription( const String &Title )
{
    printf( "///////////////////////////////////////////////////////////////////////////// \n");
    printf( " %s  \n",Title.c_str() );
    printf( "///////////////////////////////////////////////////////////////////////////// \n");
}

void BaseClass::FunctionalDescription(const String &Title, const String &Param)
{
    printf("=============================================\n");
    printf("# %s \n",Title.c_str());
    printf("# [%s] \n",Param.c_str());
    printf("=============================================\n");
}

void BaseClass::Description(const String &Title, const String &Param)
{
    printf("+---------------------------------------------------\n");
    printf("| %s   [%s] \n", Title.c_str(), Param.c_str());
    printf("+---------------------------------------------------\n");
}

//**************************************************************************
//  Copying template MS from Master Repository.
//
//   This is for Some testing items which must contain test data in the MS.
//   After copying, Modifying fuction for MS is executed depending on
//   Test cases/items.
//***************************************************************************
void BaseClass::CopyDefaultMStoWork()
{
    //  Environment //

        RunEnv env;

    // Src/Dst Path (string) 

        const String src = env.getCasaMasterPath() + "sdimaging/sdimaging.ms";
        const String dst = DefaultLocalMsName;

    // Src/Dst Path (Path) 

        casacore::Path        sourcePath(src);
        casacore::Path        targetPath(dst);       
        casacore::Directory   dir_ctrl(sourcePath);

        Description( "Copying Default MeasurementSet for modifed use to; " ,dst.c_str() );

        printf( " - src filespec  : %s \n", src.c_str() );
        printf( " - dest filespec : %s \n", dst.c_str() );

    //+
    // Copy File, use copy() method.
    //   WARNING; Destination directory must be full described.
    //-

        if(copyOperation)
        {
            dir_ctrl.copy( targetPath,
                           True,    // Overwrite 
                           True  ); // Users permisssion 
        }
}

//+
// The Working File is to be deleted 
//  whenever One Test Fixture ends. 
//-

void BaseClass::DeleteWorkingMS()
{
    String dst         = DefaultLocalMsName;

    casacore::Path        path(dst);
    casacore::Directory   dir_ctrl(path);

    Description( "Deleting Working MeasurementSet for modifed use." ,dst.c_str() );

    // Delete File (Recursively done) 
    // NOTE: for debug use, please change true-> falase )

    if (deleteOperation)
    {
         dir_ctrl. removeRecursive(false /*keepDir=False */ );
    }
}


//********************************************************
//  INTERPOLATION  Generation 
//   for Interporation Verification TEST 
//  - Generate testing trajectry in Direction
//  -  Interval(POINTING, MAIN)  (tunable in initialize() ).
//********************************************************

class CurveFunction;
class EvaluateInterporation /* : public BaseClass */ 
{
public:

    EvaluateInterporation() { printf("EvaluateInterporation constructor \n" ); /* Minimum SetUp only */  }

    // Interval Time //

        Double  getPointingTableInterval() { return pointingIntervalSec;}
        Double  getMainTableInterval()     { return mainIntervalSec;}

    // Curve Function select //

        void setCurveFunctionNo(uInt num) { printf("EvaluateInterporation:: selected Curve Function = %u \n", num ); 
                                            curveFunctionNo = defaultCurveFunctionNo = num; }
    // Calculation Error Limit //
    
         void setInterpolationErrorLimit(Double val)
                                          { printf("EvaluateInterporation:: set interpolation err. limit = %e \n", val );
                                            interpolationErrorLimit = defaultInterpolationErrorLimit = val; }

    // Init and Define Parameters 

        void Initialize();
        void Initialize( Double, Double);
        void setMainRowCount( uint n ) { requiredMainTestingRow =  defultMainTestingRow       = n; }  

    // Offset between POINTNG and MAIN //

        /* function reserved */
 
    // Pseudo Trace(Direction) for the Test

        casacore::Vector<Double>  pseudoDirInfoPointing(Double tn);
        casacore::Vector<Double>  pseudoDirInfoMain    (Double tn);

    // available POINTING TABLE count //

       uInt getAvailablePointingTestingRow() { return availablePointingTestingRow; }

    // required MAIN TABLE count 

       uInt getRequiredMainTestingRow()      { return requiredMainTestingRow; }

    //+
    // Numerical Error Statictic 
    //-
 
        Double getInterpolationErrorLimit() { return interpolationErrorLimit;  } ;

    // Row count (to add) 

        uInt getAddInerpolationTestPointingTableRow() {return addInerpolationTestPointingTableRow; };
        uInt getAddInerpolationTestMainTableRow()     {return addInerpolationTestMainTableRow; };
	
    // Antenna Count //

        uInt getNumberOfAntenna() { return numberOfAntenna; };

        bool checkExtendAvailable() 
        {
            if(availablePointingTestingRow < defInerpolationTestPointingTableRow) return true;
            else return false; 
        }

private:

       // initialize  Commmon //

       void init();

        casacore::Vector<Double>  pseudoDirInfo(Double tn);

        //+
        //  Relative Time (r_time) and Total Time
        //   for Curve Function
        //-

        Double  Interval  = 0.0; // Initilal ..
        uInt    nRow      = 0;
        Double  r_time     = 0.0;

       // Pre-located row , use tables with extended. See MS (sdimaging.ms) by tool //

        const uInt defInerpolationTestPointingTableRow   = 3843;
        const uInt defInerpolationTestMainTableRow       = 3843;

        // Row Count to execute //
        uInt requiredMainTestingRow     = 0;    //   MUST BE SET 
        uInt defultMainTestingRow       = 5000; //   copied to required when Initialize() 

        Double  requiredPointingTestingRow = 0;    //   internally calculated
        Double  availablePointingTestingRow = 0;    //   internally calculated   

        Double addInerpolationTestPointingTableRow ;   // calculated when start
        Double addInerpolationTestMainTableRow ;       // calculated when start

       // Error Limit (threshold) in GoogleTest Macro //

        Double interpolationErrorLimit ;
        Double defaultInterpolationErrorLimit = 2.0e-03 ;

        //+
        // Testing Function Select [TENTATIVE]
        //   21-JAN-2018: Migrating to isolate CurveFunction.
        //   class CurveFunction is under constrution.
        //-

        uInt curveFunctionNo;                  // Testing Function Mo/
        uInt defaultCurveFunctionNo = 0;       // Testing Function default //

    // Interval Second.

         Double pointingIntervalSec =0.0;            // Interval Time to set in POINTING 
         Double mainIntervalSec     =0.0;            // Interval Time to set in MAIN 

    // Number of Antenna 

        Double dayOffset=0       ;      // Day offset   (Day)  ** NOT USED *** 

    // Antena to USE (CAS-8418)

        const uInt numberOfAntenna =3;
};

//+
// Initialize
//-

void EvaluateInterporation::init()
{
        // Bugcheck //
        assert(requiredMainTestingRow!=0);

        printf( "EvaluateInterporation::init()::  Pointing Interval = %f\n", pointingIntervalSec);
        printf( "EvaluateInterporation::init()::  Main     Interval = %f\n", mainIntervalSec);
        printf( "EvaluateInterporation::init()::  Main Testing Row  = %d\n", requiredMainTestingRow );

        //+
        //  Define Row Count in POINTING and Main
        //   - the reuired minimum numbers are up to Interval ration.
        //-

            Double TotalTime =  requiredMainTestingRow * mainIntervalSec ;

            availablePointingTestingRow = requiredPointingTestingRow  = TotalTime / pointingIntervalSec;
 
           // if POINTING table already have sufficient length 
            if (requiredPointingTestingRow < defInerpolationTestPointingTableRow){
                    availablePointingTestingRow = requiredPointingTestingRow;
                    requiredPointingTestingRow  = defInerpolationTestPointingTableRow;
            }

         //+
         // Interpolation Error Limit 
         //-

         interpolationErrorLimit =  defaultInterpolationErrorLimit ;

        //+
        // Optimize Row Count
        //   - when required row is insufficient, set adding count and expand MS later. 
        //   - when multiple-antenna is used number of tables are increased as below.(CAS-8418)
        //-

        if ( requiredPointingTestingRow > defInerpolationTestPointingTableRow )
        {
            addInerpolationTestPointingTableRow = requiredPointingTestingRow * numberOfAntenna 
                                                   - defInerpolationTestPointingTableRow;
        }
        else
        {
            addInerpolationTestPointingTableRow = 0;
        }

        if ( requiredMainTestingRow > defInerpolationTestMainTableRow )
        {
            addInerpolationTestMainTableRow = requiredMainTestingRow * numberOfAntenna
                                                - defInerpolationTestMainTableRow;
        }
        else
        {
            addInerpolationTestMainTableRow  =0 ;
        }

        printf( "EvaluateInterporation::init()::File Size: Pointing, required =%f, adding size = %f \n", 
                requiredPointingTestingRow ,addInerpolationTestPointingTableRow);
        printf( "EvaluateInterporation::init()::File Size: MAIN    , required =%u, adding size = %f \n", 
                requiredMainTestingRow ,addInerpolationTestMainTableRow);

}

void EvaluateInterporation::Initialize(Double p_interval, Double m_interval )
{
        printf("EveInterp::initialize()::Setting Up INTERVALs " );

        pointingIntervalSec         =  p_interval;
        mainIntervalSec             =  m_interval;

        /*Tunable*/ 
        requiredMainTestingRow =  defultMainTestingRow;

        /*Tunable*/ 
        curveFunctionNo     =  defaultCurveFunctionNo;

        /*Tunable*/
        interpolationErrorLimit  = defaultInterpolationErrorLimit;

        //+
        //  Offset Time between in MAIN and in POITING 
        //    Positive value forward the time (MAIN).
        //-  

        /*Tunable*/ // second and day offset on POINTING,  UNDER CONSTRUCTION // 
 
        // common init //
        init();
}

void EvaluateInterporation::Initialize( )
{

        Double p_int = 0.048;
        Double m_int = 1.008;

        setMainRowCount(5000);
        Initialize(p_int, m_int);  // Give Pointing and Main Intervals. //
}

//*****************************************************
// Interpolation Testing Curve lcass
//  
//****************************************************

typedef void (*FUNCTYPE)(Double, Double&, Double&);

//+
// Local function definition
//     Normalized time ::     0 <=  r_time <= 1.0 
//     func(t) = {x(t), y(t) }:   (0 <= t <= 1)
//-

void Function_SimpleLinear( Double r_time, Double &X, Double &Y )
{
    X = -1.0 + 2.0 * r_time;
    Y = -1.0 + 2.0 * r_time;

    return;
}

void Function_NormalizedLinear( Double r_time, Double &X, Double &Y )
{
    // Normalized time :: | Rel_Time | < 1.0 , this case [0,1.0] is used //

    X =  (r_time * 2.0 - 1.0 ) * M_PI ;
    Y =  (r_time * 2.0 - 1.0 ) * (M_PI / 2.0);

    return;
}

void Function_sinusoid_slow( Double r_time, Double& X, Double& Y)
{
 
    X = 1.0 * cos( 2.0*M_PI  * r_time );
    Y = 1.0 * sin( 2.0*M_PI  * r_time );

    return;
}

void Function_sinusoid_quick( Double r_time, Double& X, Double& Y)
{   
    
    X = 2.0 * cos( 10* 2.0*M_PI  * r_time );
    Y = 1.0 * sin( 10*  2.0*M_PI  * r_time );
   
    return;
}   

void Function_sinusoid_hasty( Double r_time, Double& X, Double& Y)
{
    double FREQ= 100.0; 
    X = 2.0 * cos( FREQ*  2.0*M_PI  * r_time );
    Y = 1.0 * sin( FREQ*  2.0*M_PI  * r_time );

    return;
}

void Function_harmonics_sinusoid( Double r_time, Double& X, Double& Y)
{        
    const Double Amp1 = 0.5;
    const Double Amp2 = 0.6;
    const Double Omega = 2.0 * M_PI ; 
         
    Double x1  = Amp1 * cos( Omega * r_time );
    Double y1  = Amp2 * sin( Omega * r_time );

    Double x4  = Amp1/1.5 * cos( 4.0 * Omega * r_time );
    Double y4  = Amp2/1.5 * sin( 4.0 * Omega * r_time );        

    X = x1 + x4;
    Y = y1 + y4;

    return;
}

void Function_gauss( Double r_time, Double& X, Double& Y)
{

    Double t   = r_time - 0.5;
    Double A = 50;

    Double gauss  = exp (-A*t*t);

    X = gauss;
    Y = gauss;

    return;
}

void Function_Err(Double r_time, Double& X, Double& Y)
{
    X = r_time;
    Y = r_time;

    throw;
}

//****
// UNDER CoNSTRUCTION
//****
class CurveFunction
{
public:
    // Curve Type //
    enum CurveFuncType {
       SimpleLinear,
       NormalizedLinear,
       SinusoidSlow,
       SinusoidQuick,
       SinusoidHasty,
       HarmonicsSinusoid,
       Gauss,
       Err 
    };

    // constructor //
    CurveFunction(uInt no=0){ curveFunctionNo = no; }

    void calc(Double r_time, Double& X, Double& Y) {
           (*fpCurvefunc[curveFunctionNo])( r_time, X, Y ); return; }

    void setFunctionType(uInt no ) {
           if( no < fpCurvefunc.size() )  curveFunctionNo = no;}

private:

    uInt curveFunctionNo =0;

    std::vector<FUNCTYPE>  fpCurvefunc
    {    
        Function_SimpleLinear,        // 0 
        Function_NormalizedLinear,    // 1 
        Function_sinusoid_slow,       // 2 
        Function_sinusoid_quick,      // 3 
        Function_sinusoid_hasty,      // 4 
        Function_harmonics_sinusoid,  // 5 
        Function_gauss,               // 6   (new 12/11)
        Function_Err
    };   

};

//+
//  Generate Pseuo Direction / Time Infomation
//  both for Pointing and Main.
//-


FUNCTYPE fpCurvefunc[]  = 
{
    Function_SimpleLinear,        // 0
    Function_NormalizedLinear,    // 1
    Function_sinusoid_slow,       // 2
    Function_sinusoid_quick,      // 3
    Function_sinusoid_hasty,      // 4
    Function_harmonics_sinusoid,  // 5
    Function_gauss,               // 6   (new 12/11)
    Function_Err

};

casacore::Vector<Double> EvaluateInterporation::pseudoDirInfo(Double delta)
{
        casacore::Vector<Double> point;
        point.resize(5);

        //+
        //  relative time limit (r_time)
        //   (on private)
        //-

        if (r_time > 1.0) {
            printf( "r_time::Exceeded 1.0: %f \n",r_time );
            assert(r_time > 1.0);
        }
 

        //+
        //  Determin TIME 
        //    dd : in day.
        //-

        Double time  = delta * Interval;
        Double dd    =  (22 *3600.0 
                         +  5*60 +  41.5 
                         + time  
                       ) / (3600*24) ;
     
        casacore::MVTime  basetime (2003,11,12 ,dd);     
 
        //+
        // Designed Function of POINITNG location
        //-

        Double X = 0.0;
        Double Y = 0.0;

        //+
        //  Choose Function
        //  (In near future, change to Templated Function in C++ )
        //-

// rsv://        curv_func.calc( r_time, X, Y);

        (*fpCurvefunc[curveFunctionNo])( r_time, X, Y );
        
        if(false)    printf("PSD_DEBUG time=%f, r_time=%f, XY( %f,%f) \n", time, r_time, X, Y);

        // Probe the range //

	assert( abs(X) <=  M_PI );
        assert( abs(Y) <=  M_PI/2.0 );

        // Direction Values //

        point[0] = X;
        point[1] = Y;

        // Time and Interval //

        point[2] = basetime.second();                 // Time  (sec) 
        point[3] = Interval;               // Interval (sec) from Initial def.
        point[4] = r_time;                 // Relative Time r_time: [0,1]
      
        //+
        // Time on MAIN
        //   basically with offset from main.
        //   in orregular case , Time range in MAIN is not included in that of POINTNG,
        //   or the starting time in MAIN is not included Pointing Time range.
        // 
        //   X following statemment can add time shift on MAIN time.
        //-
       
        return point;

}

casacore::Vector<Double> EvaluateInterporation::pseudoDirInfoPointing(Double delta)
{
        Interval =   pointingIntervalSec;
        nRow     =   availablePointingTestingRow;
        r_time   =   delta/availablePointingTestingRow;
             
        return(pseudoDirInfo(delta));
}

casacore::Vector<Double> EvaluateInterporation::pseudoDirInfoMain(Double delta)
{
        Interval = mainIntervalSec;
              
        if(pointingIntervalSec <= mainIntervalSec ) // Ordinary case
        {
            nRow=  requiredMainTestingRow;
            r_time =   delta / nRow ;
        }   
        else   // This case may happen that not all the rows in Pointng are used.
        {
            nRow = max ( requiredMainTestingRow ,  defInerpolationTestMainTableRow );
            r_time =   delta / nRow   ;
        }
   
        return(pseudoDirInfo(delta));
}

//+ 
//  Interpolation Error statistics 
//-
class ErrorMax 
{
public:
    ErrorMax() {  MaxErr.resize(2); MaxErr[0]=0.0; MaxErr[1]=0.0;}
    void put( Vector<Double> newval )
    {
        Double v0 = abs(newval[0]);
        Double v1 = abs(newval[1]);
        MaxErr[0] = max(v0, MaxErr[0]); 
        MaxErr[1] = max(v1, MaxErr[1]);
    }
    std::vector<Double> get() { return MaxErr; }
private:

    std::vector<Double> MaxErr;

};

//************************************************
// [NEW C++]
// Access Service class 
//   for  POINTING table 
//************************************************

class PointingTableAccess
{
public:
      // Constructor //
      PointingTableAccess(String const &MsName, bool WriteAccess =false ) 
      {
//          printf("MS File [%s] \n", MsName.c_str());
          auto option = casacore::Table::TableOption::Old;
          if(WriteAccess) option = casacore::Table::TableOption::Update;
 
          MeasurementSet ms_t(MsName, option );
          ms = std::move(ms_t);
 
          init();           
          prepareColumns();

      }
      // Destructor //
      ~PointingTableAccess() { }

      // Flush and close //
      void flush() 
      {
          ms.flush();
          ms.resync();
      }
      // Row //
      uInt getNrow()   { return (nrow = hPointing.nrow()) ;  };
      void appendRow(uInt AddCnt)  { hPointing.addRow(AddCnt);}
      void removeRow(uInt row)     { hPointing.removeRow(row);}

      // Duplicate Column //
      void duplicateColumns()
      { 
          TableDesc  tblDsc = hPointing.tableDesc();

          ColumnDesc  OrgColumnDesc = tblDsc.columnDesc ( "DIRECTION" ) ;    
          ColumnDesc  RevColumnDesc1 = ColumnDesc(OrgColumnDesc);
          ColumnDesc  RevColumnDesc2 = ColumnDesc(OrgColumnDesc);
          ColumnDesc  RevColumnDesc3 = ColumnDesc(OrgColumnDesc);

          String colname1 = "POINTING_OFFSET";
          String colname2 = "SOURCE_OFFSET";
          String colname3 = "ENCODER";

          RevColumnDesc1.setName( colname1 );
          RevColumnDesc2.setName( colname2 );
          RevColumnDesc3.setName( colname3 );

          printf( "Adding 3 Columns(if not exist)  on Pointing Table \n" );

          if( ! tblDsc.isColumn( colname1 ) ){
              hPointing.addColumn(RevColumnDesc1);
              printf("Column[%s]was created.\n",colname1.c_str());
          }
          if( ! tblDsc.isColumn( colname2 ) ){
              hPointing.addColumn(RevColumnDesc2);
              printf("Column[%s]was created.\n",colname2.c_str());
          }
          if( ! tblDsc.isColumn( colname3 ) ){
              hPointing.addColumn(RevColumnDesc3); 
              printf("Column[%s]was created.\n",colname3.c_str());
          }

           pointingPointingOffset = columnPointing ->pointingOffset();
           pointingSourceOffset   = columnPointing ->sourceOffset();
           pointingEncoder        = columnPointing ->encoder();

          flush();
      }
      void fillNewColumns()
      {   printf ("Filling Data on New columns.(nrow=%d) \n",nrow);
          //+
          // Fill empty data
          //-
             IPosition Ipo = pointingDirection.shape(0);        
             /* dummy data */ 
             Array<Double> init_data1( Ipo, -0.1);
             Array<Double> init_data2( Ipo, -0.2);
             Array<Double> init_data3( Ipo, -0.3);
             for (uInt row=0; row<nrow; row++)
             {
                 // set Shape of New added Colum //
                 pointingPointingOffset. setShape(row, Ipo);
                 pointingSourceOffset.   setShape(row, Ipo);
                 pointingEncoder.        setShape(row, Ipo);
 
                 // put initial data to Column //
                 pointingPointingOffset. put( row, init_data1 );
                 pointingSourceOffset.   put( row, init_data2 );
                 pointingEncoder.        put( row, init_data3 );
              }

      }

      // Column check //
      bool checkColumn(String const &columnName ) 
      {
          String columnNameUpcase = columnName;
          columnNameUpcase.upcase();
          if (true == (ms.pointing().tableDesc().isColumn(columnNameUpcase))) return true;
          else return false;
      } 

      //+
      // Read (get) method  
      //-

      /* Vector */ 
        Array<Double> getDirection     (uInt row) { return pointingDirection      .get(row); }
        Array<Double> getTarget        (uInt row) { return pointingTarget         .get(row); }
        Array<Double> getPointingOffset(uInt row) { return pointingPointingOffset .get(row); }
        Array<Double> getSourceOffset  (uInt row) { return pointingSourceOffset   .get(row); }
        Array<Double> getEncoder       (uInt row) { return pointingEncoder        .get(row); }
 
      /* Scalor */
        Double getTime     (uInt row)       { return pointingTime.      get(row); }
        Double getInterval (uInt row)       { return pointingInterval.  get(row); }

      /// Write(put) ///
        /* Vector type */
          void putDirection      (uInt row, Array<Double> dir ) { pointingDirection.      put(row, dir ); }
          void putTarget         (uInt row, Array<Double> dir ) { pointingTarget.         put(row, dir ); }
          void putPointingOffset (uInt row, Array<Double> dir ) { pointingPointingOffset. put(row, dir ); }
          void putSourceOffset   (uInt row, Array<Double> dir ) { pointingSourceOffset.   put(row, dir ); }
          void putEncoder        (uInt row, Array<Double> dir ) { pointingEncoder.        put(row, dir ); }

        /* Scolar type */
          void putTime(uInt row, Double dd)      {  pointingTime.      put(row, dd ); }
          void putInterval(uInt row, Double dd)  {  pointingInterval.  put(row, dd ); }

private:

      // MS assign /
      void init()
      {      
           hPointing = ms.pointing();
           nrow   = hPointing.nrow();

           unique_ptr<casacore::ROMSPointingColumns>  colPt( new casacore::ROMSPointingColumns(hPointing) );
           columnPointing = std::move(colPt);
      }

      // Column Handle 
      void prepareColumns() 
      {
          pointingAntennaId      = columnPointing ->antennaId();
          pointingTime           = columnPointing ->time();
          pointingInterval       = columnPointing ->interval();
          pointingName           = columnPointing ->name();
          pointingNumPoly        = columnPointing ->numPoly();
          pointingTimeOrigin     = columnPointing ->timeOrigin();
   
          pointingDirection      = columnPointing ->direction();
          pointingTarget         = columnPointing ->target();
   
          pointingPointingOffset = columnPointing ->pointingOffset();
          pointingSourceOffset   = columnPointing ->sourceOffset();
          pointingEncoder        = columnPointing ->encoder();
      }

   // local var. //
    casacore::MeasurementSet     ms;
    casacore::MSPointing         hPointing;
    casacore::uInt               nrow;

    std::unique_ptr<casacore::ROMSPointingColumns>   columnPointing;

    ROScalarColumn<Int>    pointingAntennaId      ;
    ROScalarColumn<Double> pointingTime           ;
    ROScalarColumn<Double> pointingInterval       ;
    ROScalarColumn<String> pointingName           ;
    ROScalarColumn<Int>    pointingNumPoly        ;
    ROScalarColumn<Double> pointingTimeOrigin     ;

    ROArrayColumn<Double>  pointingDirection      ;
    ROArrayColumn<Double>  pointingTarget         ;

    ROArrayColumn<Double>  pointingPointingOffset ;
    ROArrayColumn<Double>  pointingSourceOffset   ;
    ROArrayColumn<Double>  pointingEncoder        ;

};


//************************************************
// [NEW C++]
// Access Service class 
//   for ANTENNA table 
//************************************************
// Antenna Data (for Wtiting individual record) 

typedef  struct AntTblBuf_ {
    
    String name;
    String station;
    String type;
    String mount;
    
    Vector<Double> position;
    Vector<Double> offset;
    
    Double  dish_diameter;
    Int     orbit_id;
    Double  mean_orbit[6];

} ANTENNADataBuff ;

class AntennaTableAccess {

public:
      // Constructor //
      AntennaTableAccess(String const &MsName, bool WriteAccess =false )
      {
//          printf("MS File [%s] \n", MsName.c_str());
          auto option = casacore::Table::TableOption::Old;
          if(WriteAccess) option = casacore::Table::TableOption::Update;

          MeasurementSet ms_t(MsName, option );

          ms = std::move(ms_t);

          init();
          prepareColumns();
      }

      // Destructor //
      ~AntennaTableAccess() { }

      void flush() 
      {
          ms.flush();
          ms.resync();
      }
      // Row //
      uInt getNrow()               { return (nrow = hAntenna.nrow()) ;  };
      void appendRow(uInt AddCnt)  { hAntenna.addRow(AddCnt); }
      void removeRow(uInt row)      { hAntenna.removeRow(row); }

      // block Read  //
      ANTENNADataBuff  getRowData(uInt row)
      {
        ANTENNADataBuff out;     
   
        out.name     = antennaName.    get(row);
        out.station  = antennaStation. get(row);
        out.type     = antennaType.    get(row);
        out.mount    = antennaMount.   get(row);

        out.position = antennaPosition.get(row);
        out.offset   = antennaOffset.   get(row);

        out.dish_diameter = antennaDishDiameter. get(row);

        return out;
      }

      // block Write //
      void putRowData(uInt row, ANTENNADataBuff &data)
      {

        antennaName.         put(row, data.name);

        antennaStation.      put(row, data.station);
        antennaType.         put(row, data.type);
        antennaMount.        put(row, data.mount);

        // Vector data 
        Array<Double>        data_1(ipoPosition, 0.0 ); data_1 = data.position;
        antennaPosition.     put(row, data_1);

        Array<Double>        data_2(ipoOffset, 0.0 );data_2 = data.offset;
        antennaOffset.       put(row, data_2);

        antennaDishDiameter. put(row, data.dish_diameter);

      }

private:

     void init()
     {
         hAntenna = ms.antenna();
         nrow     = hAntenna.nrow();

         unique_ptr<casacore::MSAntennaColumns>  colAnt( new MSAntennaColumns( hAntenna ) );
         columnAntenna = std::move(colAnt);
 
   //      columnAntenna = new casacore::MSAntennaColumns( hAntenna );
     }
 
     void prepareColumns()
     {
         antennaName      = columnAntenna->name();
         antennaStation   = columnAntenna->station();
         antennaType      = columnAntenna->type();
         antennaMount     = columnAntenna->mount();
         antennaPosition  = columnAntenna->position();
         antennaOffset    = columnAntenna->offset();   
         antennaDishDiameter    = columnAntenna->dishDiameter();
 
         ipoPosition = antennaPosition.shape(0);
         ipoOffset   = antennaOffset.shape(0);
     }

     casacore::MeasurementSet     ms;
     casacore::MSAntenna          hAntenna;
     casacore::uInt               nrow;

     unique_ptr<casacore::MSAntennaColumns>  columnAntenna;

     ROScalarColumn<String>  antennaName          ;
     ROScalarColumn<String>  antennaStation       ;
     ROScalarColumn<String>  antennaType          ;
     ROScalarColumn<String>  antennaMount         ;
     ROArrayColumn<Double>   antennaPosition      ;
     ROArrayColumn<Double>   antennaOffset        ;
     ROScalarColumn<Double>  antennaDishDiameter  ;

     IPosition ipoPosition;
     IPosition ipoOffset;
};


//*******************************************************
// MeasurementSet Edit Class (MsEdit)
//  - Modifying test-MS
//  - Addinng artificial(pseudo) data onto MS
//*******************************************************
class PointingTableAccess;
class MsEdit         
{
public:
 
    EvaluateInterporation  evgen;

//rsv//    CurveFunction          cvfunc(0);

    MsEdit() {     }   // Current 

    MsEdit(const String &ms, bool mode=false )
    {     
        // CAS-8418
        unique_ptr<PointingTableAccess>  ptTemp( new PointingTableAccess(ms,mode));
        unique_ptr<AntennaTableAccess>   anTemp( new AntennaTableAccess(ms,mode));
        ptTblAcc_ = std::move(ptTemp);
        anTblAcc_ = std::move(anTemp);

    }

    // Add Row //

        uInt appendRowOnAntennaTable(uInt n);
        uInt appendRowOnPointingTable(uInt n);
        uInt appendRowOnMainTable(uInt n);
    
    // Write Data on Antenna Table // 
    //            on Pointing Table //

        void writeDataToAntennaTable( uInt Row =0 );
        void writeDataToPointingTable(String MsName =DefaultLocalMsName );

    // Write (generated) Test Data on Pointing Table //

        void writeInterpolationTestDataOnPointingTable(Double dt, String MsName =DefaultLocalMsName );

    // Write (generated) Test Data in MAIN Table //

        void writeInterpolationTestDataOnMainTable(Double dt, String MsName =DefaultLocalMsName);

    // Add or Remove Column (Pointing) ////

        void duplicateNewColumnsFromDirection();

    //+
    // Default File Name
    //-

        const String MsName = DefaultLocalMsName;  // default (C++11) //

    //+
    // Buff between table and local buff to write.
    //-

        ANTENNADataBuff  AntennaData;   // for Read 
        ANTENNADataBuff  AntennaData1;  // for Write

private:

   //+
   // CAS-8418 (Under construction) 
   //  these object will provide various service , same as the current ones 
   //-
    std::unique_ptr<casa::PointingTableAccess>     ptTblAcc_;
    std::unique_ptr<casa::AntennaTableAccess>      anTblAcc_;

};


//+
// Add one row on Antanna Table
//  returns latest nrow.
//-

uInt  MsEdit::appendRowOnAntennaTable(uInt addCnt)
{
    // CAS-8418 new code //

    AntennaTableAccess ata(MsName,true);
    ata.appendRow(addCnt);
    ata.flush();
    uInt nrow = ata.getNrow();

    return nrow;
}

//+
//  Write Data to Antenna Table 
//-

void MsEdit::writeDataToAntennaTable( uInt Row )
{
    // Open MS by Update mode //

        String MsName =DefaultLocalMsName;
        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String AntennaTableName = ms0.antennaTableName(); 
        printf("MsEdit:Antanna  Table name \n" );
        printf(" [%s] \n",AntennaTableName.c_str());

    // Prepeare Handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Get current row count //

        uInt nrow_a = hAntennaTable.nrow();

        printf( "MsEdit:Antenna Table nrow  =%d \n",nrow_a);

    //
    // Get Column handle from Table  (Antenna)
    //

        std::unique_ptr<casacore::MSAntennaColumns> 
                columnAntenna( new casacore::MSAntennaColumns( hAntennaTable ));

    //+
    // LIST
    //-

    // Special Column //
        
        ROScalarColumn<String> antennaName      = columnAntenna->name(); 
        ROScalarColumn<String> antennaStation   = columnAntenna->station();
        ROScalarColumn<String> antennaType      = columnAntenna->type();
        ROScalarColumn<String> antennaMount     = columnAntenna->mount();
        ROArrayColumn<Double>  antennaPosition  = columnAntenna->position();
        ROArrayColumn<Double>  antennaOffset    = columnAntenna->offset();

        ROScalarColumn<Double>  antennaDishDiameter    = columnAntenna->dishDiameter();

        // Set Up Locate data  to be copied to Columns. //

        printf( "setting on Local Structure \n" );
        
        AntennaData1.name            =  "Z80A";
        AntennaData1.station         =  "SN123";
        AntennaData1.type            =  "KUMIKO-based";
        AntennaData1.mount           =  "IBM-PC";
        
        AntennaData1.position.resize(3);
        AntennaData1.position[0]    =  100051.01 ;
        AntennaData1.position[1]    =  -100052.02 ;
        AntennaData1.position[2]    =  100053.03 ;

        AntennaData1.offset.resize(3);
        AntennaData1.offset[0]    =   33.3 ;
        AntennaData1.offset[1]    =   44.4 ;
        AntennaData1.offset[2]    =   55.5 ;

        AntennaData1.dish_diameter  = 2.23620679;


        // Put Data to Column. //

        printf( "putting Local Structure to Column.\n" );

        antennaName.         put(Row, AntennaData1.name);
        antennaStation.      put(Row, AntennaData1.station);
        antennaType.         put(Row, AntennaData1.type);
        antennaMount.        put(Row, AntennaData1.mount);

        antennaPosition.     put(Row, AntennaData1.position);
        antennaOffset.       put(Row, AntennaData1.offset);

        antennaDishDiameter. put(Row, AntennaData1.dish_diameter);

        // Flush //
        
        ms0.flush();

}


//+
// Add rows by specified count on Pointing Table Table
//-

uInt  MsEdit::appendRowOnPointingTable(uInt AddCount )
{
    // CAS-8418 new code //

    PointingTableAccess pta(MsName,true);
    pta.appendRow(AddCount);
    pta.flush();

    uInt nrow = pta.getNrow();
   
    return nrow;
}

//+
// Write Columns on Pointing Table  
//-
void MsEdit::writeDataToPointingTable(String MsName )
{

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

        printf("MsEdit:Adding Data on new columns in POINTING table[%s]. \n",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        uInt  nrow_p = hPointingTable.nrow();
        printf( "MsEdit:Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        // create the Smart Pointer in use. //

        std::unique_ptr<casacore::ROMSPointingColumns> 
                columnPointing( new casacore::ROMSPointingColumns( hPointingTable ));

    //+
    // Listing  (Pointing) 
    //-

        ROScalarColumn<Int>    pointingAntennaId      = columnPointing ->antennaId();
        ROScalarColumn<Double> pointingTime           = columnPointing ->time();
        ROScalarColumn<Double> pointingInterval       = columnPointing ->interval();
        ROScalarColumn<String> pointingName           = columnPointing ->name();
        ROScalarColumn<Int>    pointingNumPoly        = columnPointing ->numPoly();
        ROScalarColumn<Double> pointingTimeOrigin     = columnPointing ->timeOrigin();

        ROArrayColumn<Double>  pointingDirection      = columnPointing ->direction();
        ROArrayColumn<Double>  pointingTarget         = columnPointing ->target();

        ROArrayColumn<Double>  pointingPointingOffset = columnPointing ->pointingOffset();
        ROArrayColumn<Double>  pointingSourceOffset   = columnPointing ->sourceOffset();
        ROArrayColumn<Double>  pointingEncoder        = columnPointing ->encoder();

        // Matrix Shape //

        IPosition Ipo;
        Ipo  = pointingDirection.shape(0);
        printf(" - Shape of Direction.[%ld, %ld] \n", Ipo[0], Ipo[1] );

        printf( "attempting to add Data on Pointing Table. Nrow=%d \n", nrow_p  );

        // prepare initial value //

        Array<Double> init_data1( Ipo, -0.1);
        Array<Double> init_data2( Ipo, -0.2);
        Array<Double> init_data3( Ipo, -0.3);

        for (uInt row=0; row<nrow_p; row++)
        {
            // set Shape of New added Colum //

            pointingPointingOffset. setShape(row, Ipo); 
            pointingSourceOffset.   setShape(row, Ipo);
            pointingEncoder.        setShape(row, Ipo);

            // put initial data to Column //
            
            pointingPointingOffset. put( row, init_data1 );
            pointingSourceOffset.   put( row, init_data2 );
            pointingEncoder.        put( row, init_data3 );

        }

        // Flush // 
        ms0.flush();
 
}

//+
//   Create new columns
//   copied by Direction Table 
//  
//    - This is for testing MovingSourceCorrection, 
//      only OFFSET tables are applied. 
//-

void MsEdit::duplicateNewColumnsFromDirection()
{
        String MsName = DefaultLocalMsName;
        printf( "MsEdit:Cpoied from Direction and Adding 3 Columns on Pointing Table[%s]\n", MsName.c_str() );

    // CAS-8418 new code //

        PointingTableAccess ata1(MsName,true);
         ata1.duplicateColumns();

        PointingTableAccess ata2(MsName,true);
         ata2.fillNewColumns();

}

//*****************************************************************
// Wtite Test Data on Direction Column in Pointing Table
//  - Values are got by sub fuction above.
//  - SetUp() in class TestDirection calls this.
//  - see also TEST Fixture  
//****************************************************************


void  MsEdit::writeInterpolationTestDataOnPointingTable(Double dt, String MsName)
{
    printf( "MsEdit:Writing Test Data on Direcotion Column in Pointing Table [%s]\n",MsName.c_str()  );

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        uInt nrow_p = hPointingTable.nrow();

        printf( "MsEdit:Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        // create the Smart Pointer in use. //

        std::unique_ptr<casacore::ROMSPointingColumns> 
                columnPointing( new casacore::ROMSPointingColumns( hPointingTable ));

    //+
    //   Antenna ID
    //-
        ROScalarColumn<Int>    pointingAntennaId      = columnPointing ->antennaId();
    //+
    // Time Info
    //-
        ROScalarColumn<Double> pointingTime           = columnPointing ->time();
        ROScalarColumn<Double> pointingInterval       = columnPointing ->interval();
    //+
    // Listing Columns(Direction related) on Pointing 
    //-

        ROArrayColumn<Double>  pointingDirection              = columnPointing ->direction();
        ROArrayColumn<Double>  pointingTarget                 = columnPointing ->target();

        ROArrayColumn<Double>  pointingPointingOffset         = columnPointing ->pointingOffset();
        ROArrayColumn<Double>  pointingSourceOffset           = columnPointing ->sourceOffset();
        ROArrayColumn<Double>  pointingEncoder                = columnPointing ->encoder();

        IPosition Ipo = pointingDirection.shape(0);
        printf(" - Shape of pointingDirection.[%ld, %ld] \n", Ipo[0], Ipo[1] );


    //+
    //  Loop for each Row,
    //  (NOTE) In particular case , LoopCnt requires ONE more.
    //          In ordinary case +1 causes Exceeding row count.
    //-
        uInt LoopCnt = evgen.getAvailablePointingTestingRow();

        if(evgen.checkExtendAvailable() )
        { 
            printf( "- Extend access to Pointing table = Available.\n" );
            LoopCnt += 5;
        }
   
        for (uInt ant_id=0; ant_id < evgen.getNumberOfAntenna() ; ant_id++ )
        {
            printf( "- [CAS-8418]creating  Antenna %d data. Loop=%d \n", ant_id, LoopCnt );
            for (uInt row=0; row < LoopCnt; row++)
            {
                uInt  rowA = row + (ant_id * LoopCnt);

                //+
                // Experiment   the following cannot be compiled.
                //  Set AZEL on this Colun 
                //-
 
                //+
                // DIRECTION  
                //   CAS-8418::   1-Feb-2019
                //   updated to make indivisual value on Direction Columns
                //-
                    Double delta = (Double)row + dt ;    // delta must be [0,1]

                    Array<Double> direction1(Ipo, 0.0);   // IP shape and initial val // 
                    Array<Double> direction2(Ipo, 0.0);   // IP shape and initial val // 
                    Array<Double> direction3(Ipo, 0.0);   // IP shape and initial val // 
                    Array<Double> direction4(Ipo, 0.0);   // IP shape and initial val // 
                    Array<Double> direction5(Ipo, 0.0);   // IP shape and initial val // 

                    Vector<Double>  psd_data  
                        = evgen.pseudoDirInfoPointing(delta); // generated pseudo data. (Pointing) //

                           
                    direction1[0][0] = psd_data[0];
                    direction1[0][1] = psd_data[1];

                // write access //
                // (CAS-8418: Multiple Antenna supported) use rowA instead of row//

                    pointingDirection.           put(rowA, direction1 );
                    pointingTarget.              put(rowA, direction1 );

                    pointingPointingOffset.      put(rowA, direction1 );
                    pointingSourceOffset.        put(rowA, direction1 );
                    pointingEncoder.             put(rowA, direction1 );

               //+ 
               // New Time   (intentionally activates interporation)
               //  basically FIXED values.
               //-
            
                    pointingTime.           put(rowA, psd_data[2] ); // Time
                    pointingInterval.       put(rowA, psd_data[3] ); // Interval

                    pointingAntennaId.      put(rowA, ant_id ) ;      // AntennaID 

                //+
                // Show (option)
                //-

                if(false)
                {
                    printf( "%d, Time-Interval-dirX- dirY,", row);
                    printf( "%f, %f , ", psd_data[2], psd_data[3] );
                    printf( "%f, %f \n", psd_data[0], psd_data[1] );
                }

            }//end row
        }// end antid

        // Flush //
        
        ms0.flush();
}

//+
// Add rows by specified count 
// on Main Table Table
//-

uInt  MsEdit::appendRowOnMainTable(uInt AddCount )
{
    // Measurement Set (use default name) //

        MeasurementSet ms0( MsName.c_str(), casacore::Table::TableOption:: Update );

    // Add Row //

        printf( "MsEdit:INTERPOLATION::Attempt to append [%d] new rows on Pointing Table. \n",AddCount );

        ms0.addRow(AddCount);
        uInt nrow = ms0.nrow();

        printf( "MsEdit:INTERPOLATION::   New nrow count is %d \n", nrow ); 

       ms0.flush();
       ms0.resync();

       return nrow;
}

//+
// Wtite Test Data on Direction Column in MAIN TABLE
//  -  sample shift  [ 0 < deleta < 1 ]
//  -  time diffrence = delta * interval 
//-

void  MsEdit::writeInterpolationTestDataOnMainTable(Double delta_shift, String MsName)
{
    printf( "MsEdit:INTERPOLATION::writeInterpolationTestDataOnMainTable ,1) Writing Time in MAIN Table [%s]\n",
             MsName.c_str()  );
    printf( "ddddddddddddddddddddddddddddddddddddd\n");
    printf( "   delta_shift =, %f \n", delta_shift );
    printf( "ddddddddddddddddddddddddddddddddddddd\n");

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Get current row count //

        uInt nrow_ms = ms0.nrow();
        printf( "MsEdit:INTERPOLATION::Main Table nrow =%d \n",nrow_ms);

    //+
    // Column::
    //   Time Info
    //-

        ROScalarColumn<Int>    antenna_col       ;

        ROScalarColumn<Double> mainTime_col           ;
        ROScalarColumn<Double> mainInterval_col       ;
 
    // Attach ..//

        antenna_col       .attach( ms0  , "ANTENNA1");

        mainTime_col      .attach( ms0 , "TIME");
        mainInterval_col  .attach( ms0 , "INTERVAL");

    // WRite .. //

        uInt LoopCnt = evgen.getRequiredMainTestingRow();

        for (uInt ant_id=0; ant_id < evgen.getNumberOfAntenna() ; ant_id++ )
        {
            for (uInt row=0; row < LoopCnt; row++)
            {
                uInt  rowA = row + (ant_id * LoopCnt);

                // Pseudo Data (TEST DATA );

                   Vector<Double>  psd_data 
                       = evgen.pseudoDirInfoMain( (Double)row); // generated pseudo data. (Main table) //

                // Time Set  //
                    Double interval = psd_data[3];
                    Double time     = psd_data[2];
                
                    Double SetTime = time + (delta_shift * interval) ; 

                    antenna_col.            put(rowA, ant_id   );      // AntennaID (CAS-8418)
                    mainTime_col.           put(rowA, SetTime  );      // Time     (( REvised 10.26))
                    mainInterval_col.       put(rowA, interval );      // Interval (( Revised 10.26))

                // Antenna //
                // Show Time //

                    if (false){
                            printf( "Main Tbl, %d SetTime=%f,  Interval=%f, \n",
                            row, SetTime, interval );
                }
            }
        }

        // Flush //
        
        ms0.flush();

}







//+
// PointingDirectionCalculator Class 
// Constructor
//-
class TestMeasurementSet : public BaseClass
{

public:

protected:

        TestMeasurementSet() { }
        ~TestMeasurementSet() { }

        virtual void SetUp()
        {
            BaseClass::SetUp();
        }

        virtual void TearDown()
        {
            BaseClass::TearDown();
        }

        // Test Fixture Sub //
        void test_constructor(String const name );
};

/*---------------------------------------
  attempt to open vaious Measurement Set
 ----------------------------------------*/
void TestMeasurementSet::test_constructor(String const name)
{
    Description("Testing Construcror." , name );

    // CONSTRUCTOR  //

        MeasurementSet ms0( name  );      
        PointingDirectionCalculator calc(ms0);

    // Initial brief Inspection //
   
        printf("# Constuctor Initial Check.  [%s] \n", name.c_str()  );
        printf("  detected nrow : %d \n", calc.getNrowForSelectedMS() );
      
        EXPECT_NE((uInt)0, calc.getNrowForSelectedMS() );
}

/*--------------------------------------------
 *  Constructor test by Vaisous Measurment Set 
 *
 *  - prepare some MSs.
 *  - Most of them works normal, some throws
 *    Exception.
 *  - MS descrition is defines in this file.
 * --------------------------------------------*/

TEST_F(TestMeasurementSet, variousConstructor )
{
    // MS name database //
    MSNameList  MsList;
    TestDescription( "CALC Constructor by various MS " );

    for(uInt n=0; n< MsList.count(); n++)
    {
        FunctionalDescription( "CALC Constructor by various MS", 
                               std::to_string(n)+". " + MsList.getName(n).c_str()  );
        String name = env.getCasaMasterPath()+ MsList.getName(n);
        if ( MsList.isExceptionActivated(n))
        {
            EXPECT_ANY_THROW( test_constructor(name) );
        }
        else
        {
          EXPECT_NO_THROW( test_constructor(name) );
        }
    }
}

/*-----------------------------------------------------
  Inspect RowId consistency

 Step 1: make SelectedMS by selectData()
         a pair of Antena is the key (6 combiniations)
 Step 2:  execute getRowId() 
 Step 3:  record the obtained RowId
 Step 4: Insect that 1) All the ID appear, 
          2) No duplication
 -----------------------------------------------------*/

TEST_F(TestMeasurementSet, RowId_inMS )
{

    TestDescription( "RowId functions  by various MS"  );

    String MsName = "listobs/uid___X02_X3d737_X1_01_small.ms";
     FunctionalDescription("Testing RowId functions." , MsName );

    String name = env.getCasaMasterPath()+MsName;

   // Measurment Set and Constructor //

     MeasurementSet ms0( name  );       
     PointingDirectionCalculator calc(ms0);
     uInt nrow0 = calc.getNrowForSelectedMS(); 

   // Prepare  Arraay ..//

    Vector<int> rowIdCheckList(nrow0,0);

    std::vector<String> selectList = {
      "DV01&&DV01",
      "DV02&&DV02",
      "PM03&&PM03",
      "DV01&&DV02",
      "DV01&&PM03",
      "DV02&&PM03",
    };

    for(uInt sw=0; sw <selectList.size() ; sw++)
    {
         //+
         // setlectData() key 
         //  by Antenna Pair
         //-

         String ant_sel = selectList[sw];

         Description("Testing RowId functions." , "case="+std::to_string(sw) );

        // getRowIdForOriginalMS //

          Vector<uInt> vRowIdOrgMS = calc.getRowIdForOriginalMS();

         printf( " calling selectData() \n");     
         calc.selectData( ant_sel,  "","","","","","","","","" );

        // Nrow from Selected //

          uInt nrow = calc.getNrowForSelectedMS();
          printf("Selected nrow =%d\n",nrow);

        // Vecrtor<uInt> getRowID() 
        
          Description("(1) Vector<uInt> getRowId() ", name );
          Vector<uInt> vRowId = calc.getRowId();

        // getRowId( int ) //
 
          Description("(2) uInt getRowId() ", name );

        // Show and Verify //
          
          printf( "Num Row (Org) = %d \n" , nrow0 );
          printf( "Num Row (Sel) = %d \n" , nrow );
          printf( "    key=[%s]\n", ant_sel.c_str() ); 

          for (uInt k=0; k < nrow; k++)
          {
              uInt RowId = calc.getRowId(k);
              printf( "RowID ,%d, checked \n", RowId ); 
              
              // Check List //
             rowIdCheckList[RowId] ++;
          }

    } // end for 
    
    //*
    // RowId List 
    //  - checck all the ID have been appeared once,
    //  - Duplicattion IS NOT Allowed
    //* 

    for (uInt i=0; i < nrow0; i++)
    {
        EXPECT_EQ( rowIdCheckList[i],1 );    
    }
}

/*---------------------------------------
  getDirection  base class
 ----------------------------------------*/
class TestDirection : public BaseClass
{

public:
        uInt const InterpolationDivCount = 10;
protected:

        //+
        // NEW: CAS-8418 related
        //
        
        // select Direction column (Internal Test) //
        
         uInt  GTestSelectDirectionColumn = 0;

        // MeasurementSet Editting
     
            casa::MsEdit  msedit;

        // Add 3 OFFSET Colums (basically reserved) //
          void addColumnsOnPointing();

        // Add avalable data on OFFSET columns  //
          void addColumnDataOnPointing();

        // Sub-function of TEST_F(TestDirection....)

        std::vector<Double>  testDirectionForDeltaTime(Double dt);
        vector<Double>       testDirectionByInterval(Double p, Double m);


        TestDirection()
        {
        }

        ~TestDirection()
        {
        }

        virtual void SetUp()
        {
            BaseClass::SetUp();

            //+
            // Copy and add culumns, 
            //  init columns, 
            //  generate test dta
            //-
            
            CopyDefaultMStoWork();

            addColumnsOnPointing();
            addColumnDataOnPointing();   // FILL DATA

        }

        virtual void TearDown()
        {
            BaseClass::TearDown();

            // Delete Working MS 

             DeleteWorkingMS();
        }

};
/*--------------------------------------
   Add 3 new Columns and locate data 
   on PPOINTING table
  --------------------------------------*/

void TestDirection::addColumnsOnPointing()
{
    msedit.duplicateNewColumnsFromDirection();
}


void TestDirection::addColumnDataOnPointing()
{
#if 0
    msedit.writeDataToPointingTable( DefaultLocalMsName );
#endif 
}
 
/*---------------------------------------
  setDirectionColumn
  - check thr accsessor 
  - Must not make exception.
 ----------------------------------------*/

//+
// Sub Function (reserved) 
//  - Method: getAccessor() must be coded in the header 
//    to return the accessor address
//  - At the moment, getAccessor() is not implemented
//-  

static void assert_accessor( PointingDirectionCalculator  &calc )
{

    uInt Id = calc.getCurretAccessorId();
    printf( "# Currnet Accessor ID = %d\n",Id);
  
}


TEST_F(TestDirection, setDirectionColumn  )
{

    TestDescription( "setDirectionColumn (String Fram)" );
    const String MsName = DefaultLocalMsName;

    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    // setDirectionColumm()  calls 
    //
    //  (note) Checcking 'accessor' address reqires to change header.
    //         This is not urgent unless, any address error happens.
    //-
        String ColName;

    for( uInt n=0; n<2;n++ ) 	// 2 Times. run .../
    {     
        ColName = "DIRECTION";
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);

        ColName = "TARGET";
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);       
 
        ColName = "POINTING_OFFSET";    // NEED to ADD Table in advance  //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);

        ColName = "SOURCE_OFFSET"; // NEED to Add Table in advance //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);
 
        ColName = "ENCODER";      // NEED to add Table in advance  //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);

        ColName = "hogehoge";
        Description("Column Name" , ColName );
        EXPECT_ANY_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor(calc);

    }    
}

/*----------------------------------------------------------

  getDirection() and MovingSourceCorrection

 - Test performing MovingSource Correction.
 - when "POINTING_OFFSET" or "SOURCE_OFFSET" is specified,
  ignore this functonality.

  xx detail description of this fixture , in progress.
------------------------------------------------------------*/
TEST_F(TestDirection, MovingSourceCorrection  )
{

    TestDescription( "performMovingSourceCorrection and setDirectionColumns" );
    const String MsName = DefaultLocalMsName;     

    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    // setDirectionListMatrixShape           //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    // setFrame()  //

        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //+
    // setDirectionColumm()  calls 
    //-
        const vector<String> ColName
        {
            "DIRECTION",       
             "TARGET",
             "POINTING_OFFSET",  // *** NEED to ADD in advance  //
             "SOURCE_OFFSET",    // *** NEED to Add in advance //
             "ENCODER",          // *** NEED to add in advance  //
             "DIRECTION"         // extra //
        };

        //+
        // Normal Seq. with setMovingSoure Convert.
        //  using String , deffault values are set.
        //-

	FunctionalDescription("Normal Seq.", "Selectve Convert");

        for(size_t k=0; k < ColName.size(); k++)
        {
            Description("Column Name" , ColName[k] );
            EXPECT_NO_THROW( calc.setDirectionColumn( ColName[k] ) );

            if(true)  // Selective:: TESTING dependency how  setMovingSource() relates. //
            {
                String src = "SUN";
                Description("setMovingSource()", src);
                EXPECT_NO_THROW( calc.setMovingSource( src ) );
            }    

            Description("calling  getDirection() ", ColName[k] );
            Matrix<Double>  DirList;
            EXPECT_NO_THROW( DirList= calc.getDirection() );

            uInt  n_row    = DirList.nrow();

            printf( "Number of Row = %d \n", n_row );
            EXPECT_EQ( n_row, ExpectedNrow);

        }

        //+
        // No SetMovingSouce executution Sequence. 
        //-
        
        FunctionalDescription("Normal Seq.", "Always call setDirectionColumn");

        for(uInt k=0; k < ColName.size(); k++)
        {
            Description("Column Name" , ColName[k] );
            EXPECT_NO_THROW( calc.setDirectionColumn( ColName[k] ) );

            String src = "SUN";
            Description("setMovingSource()", src);
            EXPECT_NO_THROW( calc.setMovingSource( src ) );
            
            Description("calling  getDirection() ", ColName[k] );
            Matrix<Double>  DirList;   
            EXPECT_NO_THROW( DirList = calc.getDirection() );

            uInt  n_row = DirList.nrow();

            printf( "Number of Row = %d \n", n_row);
            EXPECT_EQ( n_row, ExpectedNrow);

        }

}

/*-------------------------------------------
  Verification Test of CAS-11818
  - If you use old source, this test causes
    core dump.

   XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   XXXX   DO NOT UPDATE !!  XXXX
   XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
---------------------------------------------*/
TEST_F(TestDirection, VerifyCAS11818 )
{

    TestDescription( "configureMovingSourceCorrection(CAS11818) Test" );
    const String MsName = DefaultLocalMsName;    //  

    // MS name for this Test //

        String name =   MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
   
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
        printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
        ExpectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, ExpectedNrow );


    // setFrame, setDirectionListMatrixshape    //
    // setMovingSource Convert                  //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //
    // setDirectionColumm()  calls 
    //
        String Name[6];
        
        Name[0] = "DIRECTION";        
        Name[1] = "TARGET";
        Name[2] = "POINTING_OFFSET";    // *** NEED to ADD in advance  //
        Name[3] = "SOURCE_OFFSET";      // *** NEED to Add in advance //
        Name[4] = "ENCODER";            // *** NEED to add in advance  //
        Name[5] = "hogehoge";         /// NOT USED ///

         
        // No SetMovingSouce executution 

        FunctionalDescription("BUG-TEST by assert()" , "setDirectionColumn+ getDirection, without setMovingSource");

        for(uInt k=0; k<5; k++)
        {
            Description("Column Name" , Name[k] );
            EXPECT_NO_THROW( calc.setDirectionColumn( Name[k] ) );

            Description("calling  getDirection() ", Name[k] );
            Matrix<Double>  DirList;   
#if 1
            EXPECT_NO_THROW( DirList= calc.getDirection() );
#else
            DirList= calc.getDirection();
#endif 
            uInt  N_Row;
            N_Row  = DirList.nrow();

            printf( "Number of Row = %d \n", N_Row );
            EXPECT_EQ( N_Row, ExpectedNrow);

        }

}

/*-----------------------------------------------
  UnsetMovingSource and 
   with the Conveniation of setMovingSource 
 - 5 different (C++11)senario, with set/unset / No op.
------------------------------------------------*/
TEST_F(TestDirection, setMovingSource  )
{

    TestDescription( "performMovingSourceCorrection and setDirectionColumns" );
    const String MsName = DefaultLocalMsName;    //  

     // List all info on Pointing Table. //
     //   _List series. was Removed. 
     //   Replaace to Dump series.

    // Create Object //
    
        MeasurementSet ms( MsName ); 
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
        printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
        ExpectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, ExpectedNrow );

    // setDirectionListMatrixShape()    //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    // setFrame() //
    
        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //
    // setDirectionColumm()  calls 
    //
        const std::vector<String> ColName
        {
            "DIRECTION",        
            "TARGET",
            "POINTING_OFFSET",  // *** NEED to ADD in advance  //
            "SOURCE_OFFSET",    // *** NEED to Add in advance //
            "ENCODER",          // *** NEED to add in advance  //
//          "hogehoge"          // => This makes Exception 
        };

    //+
    //  Four senarios are executed.
    //-  	

    for(uInt senario = 0; senario <= 3; senario++ )	// senario [0,1,2,3] //
    {   
          FunctionalDescription("Senario" , std::to_string(senario).c_str() ); 

          for(size_t k=0; k < ColName.size(); k++)
          {
              Description("- Column Name" , ColName[k] );
              EXPECT_NO_THROW( calc.setDirectionColumn( ColName[k] ) );

              if( senario==0 )     // Always No call //
              {
                /*  nothing */
              }
              else
              if( senario==1 )   // Allways execute setMovingSource() //
              {
                 String src = "SUN";
                 Description("- setMovingSource()", src);
                 EXPECT_NO_THROW( calc.setMovingSource( src ) );
              }
              else
              if( senario==2 )   // Allways execute unsetMovngSource() //
              {
                 String src = "SUN";
                 Description("- unsetMovingSource()", src);
                 EXPECT_NO_THROW( calc.unsetMovingSource() );
              }
              else 
              if( senario==3 )   // Once set and unset MovingSource() //
              {
                 String src = "SUN";
                 Description("- setMovingSource() consequently call unsetMovingSource()", src);

                 EXPECT_NO_THROW( calc.setMovingSource( src ) );
                 EXPECT_NO_THROW( calc.unsetMovingSource() );
              }

              Matrix<Double>  DirList;

              //+
              // getDirection()
              //-

              Description("- getDirection() ", ColName[k] );
              EXPECT_NO_THROW( DirList  = calc.getDirection() );

              uInt n_row;
              n_row  = DirList.nrow();

              printf( "- Number of Row = %d \n", n_row );
              EXPECT_EQ( n_row, ExpectedNrow);

          }
    }

}

/*------------------------------------------
   Matric Shape ( COLUMN or ROW )
   - Inspect  Coumn/Row Number shoud be 
     in the specification. 
  ------------------------------------------*/
TEST_F(TestDirection, Matrixshape )
{

    TestDescription( "setDirectionListMatrixShape()" );
    const String MsName = DefaultLocalMsName;    //  
    
    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
    
       PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    // A set of API Call   
    //    - setDirectionListMatrixShape() and getDirection() 
    //    - getDirection return Matrix<Double> and the shape is determined in 
    //      getDirection by IPPosition(p,q,r)
    //-

        Matrix<Double> DirList1;
        Matrix<Double> DirList2;
        uInt N_Col;  
        uInt N_Row ;  

        // COLUMN //
    
        Description("setDirectionListMatrixShape", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );
        
        DirList1  = calc.getDirection();
        N_Col    = DirList1.ncolumn();
        N_Row    = DirList1.nrow();
        printf ("# NCol = %d , NRow = %d \n", N_Col, N_Row );

        EXPECT_EQ( N_Col, (uInt)2);

        // ROW  //

        Description("setDirectionListMatrixShape", "ROW_MAJOR");
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::ROW_MAJOR) );

        DirList2  = calc.getDirection();
        N_Col    = DirList2.ncolumn();
        N_Row    = DirList2.nrow();
        printf ("# NCol = %d , NRow = %d \n", N_Col, N_Row );

        EXPECT_EQ( N_Row, (uInt)2);

}

/*------------------------------------
  getDirection () and MovingSource()
 -------------------------------------*/
TEST_F(TestDirection, getDirection1 )
{

    TestDescription( "getDirection (J2000) No data selection" );
    const String MsName = DefaultLocalMsName;    

    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    // setDirectionColumn() 
    //-
  
        String ColName = "DIRECTION"; 
        Description( "getDirectionColumn()", ColName  );

        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

    //+
    //  MatrixShape (COLUMN_MAJOR) 
    //-
        Description("calling setDirectionListMatrixShape()" ,"Column Major" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    //+
    //  setMovingSourceDirection() 
    //-
        String src = "MOON";
        Description("calling setMovingSource()", src);

        EXPECT_NO_THROW( calc.setMovingSource( src ) ); 

    //+
    //  getDirection()
    //-

        Description("calling  getDirection() ","" );

        Matrix<Double>  DirList1  = calc.getDirection();
        uInt  n_row    = DirList1.nrow();

        printf( "Number of Row = %d \n", n_row );
        EXPECT_EQ( n_row, ExpectedNrow);

    //+
    // Dump Matrix
    //-

    if (true) 
    {
        Description("Dump obtined Direction info. ","" );

        for (uInt row=0; row< n_row; row++)
        {
            // Direction //

            Double Val_1 = DirList1(row,0);
            Double Val_2 = DirList1(row,1);

            casacore::MDirection  MovDir  = calc.getMovingSourceDirection();
            String strMovDir = MovDir.toString();

            printf(    "Dir at, %d, %f,%f, [Mov:%s]  \n",  
                    row, Val_1, Val_2, strMovDir.c_str() );
         }

    }

}

/*--------------------------------------
   Interpolation Test in getDirection
 ---------------------------------------*/

void listPointingTable(String MsName)
{
    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str());

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();
        printf("Pointing Table name [%s] \n",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        uInt nrow_p = hPointingTable.nrow();

        printf( "Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        std::unique_ptr<casacore::ROMSPointingColumns> 
                columnPointing( new casacore::ROMSPointingColumns( hPointingTable ));

    //+
    // Listing  (Pointing) 
    //-

        ROScalarColumn<Double> pointingTime           = columnPointing ->time();
        ROScalarColumn<Double> pointingInterval       = columnPointing ->interval();
        ROArrayColumn<Double>  pointingDirection      = columnPointing ->direction();
        ROArrayColumn<Double>  pointingTarget         = columnPointing ->target();

        printf( "================================\n");
        printf( " Pointing (Time and Dir) \n");
        printf( "-----------------+--------------\n");

        for (uInt row=0; row<nrow_p; row++)
        {

            Vector<Double> valDirection  =   pointingDirection. get(row);
            Vector<Double> valTarget     =   pointingTarget.    get(row);
	    Double         valTime       =    pointingTime.     get(row);

            printf( "Pointing,(pos time dir target), %d , %f, , %f, %f, ,%f, %f)  \n", 
                    row,  valTime, 
                    valDirection[0],      valDirection[1],
                    valTarget[0],         valTarget[1] );

        }

}
    

/*---------------------------------------------
    Interporation Test (sub) 
     dt : displacement betweeen measure point
     [0 <= dt <= 1]   dt=0; X[n],  dt=1, X[n+1]
  ---------------------------------------------*/

std::vector<Double>  TestDirection::testDirectionForDeltaTime(Double dt )
{
    TestDescription( "testDirectionForDeltaTime(dt)" ); 
    printf( " dt = %f (sec)\n",dt ); 

    const String MsName = DefaultLocalMsName;

    //  MS (Dump) ** develper option ** //

        if(false)  listPointingTable(MsName);

    // Create Object //

        MeasurementSet ms( MsName.c_str() );
        PointingDirectionCalculator calc(ms);   

    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    //  MatrixShape (COLUMN_MAJOR) 
    //-
        Description("calling setDirectionListMatrixShape()" ,"Column Major" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) ); 

    //+
    // setFrame()
    //-

        String FrameName= "J2000";
        EXPECT_NO_THROW( calc.setFrame( FrameName ));

    //+
    // selectData
    //    0=DA61, 1=PM03, 2=PM04
    //-

        if(true)
        {
              const String AntSel = "GBT&&&";
            calc.selectData( AntSel,  "","","","","","","","","" );
 
        }

    //+
    // setDirectionColumn() 
    //-
  
        String DColName[] = {"DIRECTION", "TARGET", "POINTING_OFFSET", "SOURCE_OFFSET", "ENCODER" }; 
        
        Description( "getDirectionColumn()", DColName[GTestSelectDirectionColumn]  );

        EXPECT_NO_THROW( calc.setDirectionColumn( DColName[GTestSelectDirectionColumn] ) );

    //  Dump (before)

        if(false)
        {
            Description("Dump Pointing (before getDirection) ","" );
            listPointingTable( MsName );
        }

    //+
    //  InterPolation mode (Foece Unuse)
    //-
 
        calc.setSplineInterpolation(use_spline);

    //+
    //  getDirection()
    //-

        Description("Calling  getDirection() ","" );

        Matrix<Double>  DirList1  = calc.getDirection();
        size_t   n_col    = DirList1.ncolumn();
        size_t   n_row    = DirList1.nrow();

        printf( "Number of Row = %zu \n", n_row );
        printf( "Number of Col = %zu \n", n_col );

    //+
    // Dump Matrix
    //-

    Description("Inspecting  Direction Inerpolation starts.","" ); 

    Double maxErr_1 = 0.0;
    Double maxErr_2 = 0.0;             
    Double absErr_1 = 0.0;
    Double absErr_2 = 0.0;

    if(true)
    {
        printf("Row, cal_x, cal_y, gen_x. gen_y, err_x, err_y \n");

        uInt LoopCnt = n_row-2;
        for (uInt row=0; row < LoopCnt; row++)   
        {
            // Direction(1) by getDirection //

              Double calculated_1 = DirList1(row,0);
              Double calculated_2 = DirList1(row,1);

            // Direction(2) by generated/estimated //

              casacore::Vector<Double>  gen_out 
                      = msedit.evgen.pseudoDirInfoMain(  (Double)row   
                                                          + dt         );    // dt:Interpolation offset  (sec)
            //+
            // Error calculation 
            //-

                Double generated_1 = gen_out[0]; 
                Double generated_2 = gen_out[1]; 
                
                Double Err_1 = calculated_1 - generated_1 ;   
                Double Err_2 = calculated_2 - generated_2 ;

                absErr_1 = abs(Err_1);
                absErr_2 = abs(Err_2);
                
                if( maxErr_1 < absErr_1 ) maxErr_1 = absErr_1;
                if( maxErr_2 < absErr_2 ) maxErr_2 = absErr_2;

            //+ 
            // Google Test
            //-

		EXPECT_LE( absErr_1, msedit.evgen.getInterpolationErrorLimit()  ); 
                EXPECT_LE( absErr_2, msedit.evgen.getInterpolationErrorLimit()  ); 

            // Output List (in One line)//

                if(false) 
                {
                    printf( "Evaluation,");
                    printf( "%6d, %13.10f,%13.10f,", row,  calculated_1, calculated_2 );
                    printf( "%13.10f,%13.10f,",      generated_1,  generated_2 );
                    printf( "%12.4e,%12.4e \n",      Err_1,     Err_2);
                }
        }    
    }
   
   std::vector<Double> vv(2);
   vv[0] = maxErr_1;   
   vv[1] = maxErr_2;

   return vv; 
}

//+
// TestSub
//  - Pointing table interval and Main table interval are specified
//  - execute numerical error test by activating Interpolation
//    by using slided time in Main Table.
//-
std::vector<Double> TestDirection::testDirectionByInterval(Double p_int, Double m_int)
{
    printf("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
    printf("      INTERPOLATION:: testing [%f,%f ] \n", p_int,m_int );
    printf("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");

    // Max Error ///

      std::vector<Double> reterr;
      ErrorMax            maxerr;
    
    // Add INTERPOLATION TEST DATA 

      printf ("calling TestwriteInterpolationTestDataOnPointingTable\n"); 
      msedit.writeInterpolationTestDataOnPointingTable( 0.0 );  // Pointing //

    // Test Loop  
    //   - nDiv is given, internd to spacify number of separating count 
    //     between one exact point and the next,  

    uInt nDiv = InterpolationDivCount; 

    for (uInt loop=0; loop < nDiv; loop ++ )
    {
         //+
         // SetUp Testing  MeasurmentSet
         //-

           Description("INTERPOLATION:: Making MeasurementSet",
                        "INTERPOLATION::k="+to_string( (Double)loop/(Double)nDiv  ));

           msedit.writeInterpolationTestDataOnMainTable( (Double)loop/(Double)nDiv  );

        // excuttion ..// 
        
          Description("Execution starts. ","" );
          reterr = testDirectionForDeltaTime( (Double)loop/(Double)nDiv );

          printf( "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");
          printf( " Max Error =, %e, %e \n", reterr[0], reterr[1] );
          printf( "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");

          maxerr.put(reterr);
    }
    std::vector<Double> e_max = maxerr.get();

    printf("----------------------------------------\n");
    printf("INTERPOLATION:: testing [%f,%f] END     \n", p_int, m_int );
    printf("----------------------------------------\n");

    return e_max; 
}
/*-----------------------------------------------------------------------
   Interporatio Test in getDirection()   Multiple-Combiniation

 - Use Combiniation of Pointing table Interval and Main table Interval
 - Capable of selecting curve fucntion for simulated pointing trajectry
  ----------------------------------------------------------------------*/
 
TEST_F(TestDirection, InterpolationFull )
{
    // Combiniation List of Pointing Interval and Main Interval //

    vector<bool>   InterpolationMode     = { false, true };
    vector<Double> Main_IntervalList     = { 1.0 };
    vector<Double> Pointing_IntervalList = { 1.0, 0.5, 0.2, 0.1, 0.05  };

    ErrorMax  maxerr;
    std::vector<Double> r_err = {0.0}; 

  for(uInt s=0; s<InterpolationMode.size(); s++)
  {
    use_spline = InterpolationMode[s];

   // Error Limit 
    msedit.evgen.    setInterpolationErrorLimit( 1e-04 );

   // Combiniation Loop 
   for( uint m=0; m < Main_IntervalList.size(); m++)
    {
        for( uint p=0; p < Pointing_IntervalList.size(); p++)
        {
             Double p_i = Pointing_IntervalList[p];
             Double m_i = Main_IntervalList[m];
  
             SetUp();
             msedit.evgen  .    setCurveFunctionNo(0);   // set Curve Fuction
             msedit.evgen.      setMainRowCount(5000);
             msedit.evgen.      Initialize(p_i, m_i) ;

           // Prepate Antenna (For Multiple-set)  //

             msedit.appendRowOnAntennaTable(2);   // 1 + 2 more
             msedit.writeDataToAntennaTable(1);
             msedit.writeDataToAntennaTable(2);

           // Increase Row on MS for large-file.

              msedit.appendRowOnPointingTable (msedit.evgen.  getAddInerpolationTestPointingTableRow() );
              msedit.appendRowOnMainTable     (msedit.evgen. getAddInerpolationTestMainTableRow() );

              addColumnDataOnPointing();   // FILL DATA 
 
           //
           // Execute Main-Body , get error info //
           //
              r_err = TestDirection::testDirectionByInterval( p_i, m_i );
              maxerr.put(r_err);

        }
    }

    std::vector<Double>  e = maxerr.get();
    printf ( "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n");
    printf ( "   THE WORST Error = %e, %e \n", e[0], e[1] );
    printf ( "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n");

  }
}

/*-----------------------------------------------------------------------
  Interporatio Test in getDirection()   SINGLE mode.

- Set up one set of Pointing table Interval and Main table Interval
- Table sizes to be created is automatically tuned.
- Capable of selecting curve fucntion for simulated pointing trajectry
 ----------------------------------------------------------------------*/

TEST_F(TestDirection, InterpolationSingle )
{

    TestDescription( "Interpolation test in getDirection() SINGLE-parameter mode" );

    //+
    //  Select Function
    //    - See EvaluateInterporation class
    //
    //  Set Testing Count
    //    - define test count. some rows are automatically added
    //-

      use_spline = true;

      msedit.evgen.    setCurveFunctionNo(0);   // set Curve Fuction
      msedit.evgen.    setMainRowCount   (5000);  // aprox. 1-2H 
      msedit.evgen.      Initialize( 2.99827,     // Pointing Interval
                                     2.99827 ) ;  // Main Interval
 
      msedit.evgen.    setInterpolationErrorLimit( 1.0E-04 );

    // Prepate Antenna (for Multple-set) //

      msedit.appendRowOnAntennaTable(2);   // 1 + 2 more
      msedit.writeDataToAntennaTable(1);
      msedit.writeDataToAntennaTable(2);

    // Increase(Append)  Row on MS for large-file.:

      msedit.appendRowOnPointingTable ( msedit.evgen.  getAddInerpolationTestPointingTableRow() );
      msedit.appendRowOnMainTable     ( msedit.evgen.  getAddInerpolationTestMainTableRow() );
      addColumnDataOnPointing();    

    //+
    //   Pointing TBL and MAIN TBL
    //-
 
      std::vector<Double> r_err = {0.0};

    // get currect interval time.. //
    
      Double p_interval = msedit.evgen.getPointingTableInterval();
      Double m_interval = msedit.evgen.getMainTableInterval();

     // Execute and get numerical error info 

       r_err = TestDirection::testDirectionByInterval(p_interval, m_interval);

       printf( "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");
       printf( " Total Max Error = %e, %e \n", r_err[0], r_err[1] );
       printf( "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");


}

TEST_F(TestDirection, CompareInterpolation )
{
  TestDescription( "Interpolation Result Comparing.(Linear vs Spline)" );

    std::vector<Double>   r_err1 ;
    std::vector<Double>   r_err2 ;
   
    TestDescription( "getDirection (J2000) with selected data. uvw available" );
 
    // Use DefaultMS as a simple sequence.
    // Use the below to show uv valuses from MS.

    std::vector<String> MsList = {
      "sdimaging/Uranus1.cal.Ant0.spw34.ms",
      "sdimaging/Uranus2.cal.Ant0.spw34.ms"
    };

    for(uInt fno=0; fno<MsList.size(); fno++)
    {
        // Direction 
        Matrix<Double>  DirList1; // for Linear
        Matrix<Double>  DirList2; // for Spline

        // selected MS name //
        String name = env.getCasaMasterPath()+MsList[fno]; 
        printf( "MS[%s] is used. \n",name.c_str() );

        // Create Object //
        MeasurementSet ms0( name );
        PointingDirectionCalculator calc(ms0);

        //+
        // setDirectionColumn() 
        //-

        String ColName = "DIRECTION";
        Description( "getDirectionColumn()", ColName  );
        calc.setDirectionColumn( ColName ) ;

        //+
        //  MatrixShape (COLUMN_MAJOR) 
        //-

        Description("calling setDirectionListMatrixShape()" ,"Column Major" );
        calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) ;

        //+
        // setFrame()
        //-

        if(true)
        {
            String FrameName= "AZELGEO";
            calc.setFrame( FrameName );
        }

        if(true)
        {
     
            // Set Interporation Mode //
         
            calc.setSplineInterpolation(false);
        
            //+
            //  getDirection()
            //-
          
            Description("calling getDirection()" ,"#1" );
         
            DirList1  = calc.getDirection();
            size_t   n_col    = DirList1.ncolumn();
            size_t   n_row    = DirList1.nrow();
         
            printf( "getDirection()::Number of Row = %zu \n", n_row );
            printf( "getDirection()::Number of Col = %zu \n", n_col );
    
        }

        if(true)
        {

            // Set Interporation Mode //

            calc.setSplineInterpolation(true);

            //+
            //  getDirection()
            //-

            Description("calling getDirection()" ,"#2" );
 
            DirList2  = calc.getDirection();
            size_t   n_col    = DirList2.ncolumn();
            size_t   n_row    = DirList2.nrow();
 
            printf( "getDirection()::Number of Row = %zu \n", n_row );
            printf( "getDirection()::Number of Col = %zu \n", n_col );
        }

        uInt size = DirList1.nrow();
        for( uint row=0; row<size; row++)
        {
            Double er1 = DirList2(row,0) - DirList1(row,0);
            Double er2 = DirList2(row,1) - DirList1(row,1);

            printf("row[%4d] ", row );
            printf("Dir1=, %e,%e  ,", DirList1(row,0),DirList1(row,1) ); 
            printf("Dir2=, %e,%e  ,", DirList2(row,0),DirList2(row,1) );
            printf("Err =, %e,%e \n", er1, er2 );
        }
    }
}

/*---------------------------------------------------
    getDirection  with uvw data dump,
     - Ordinary (standard sequence) 
  --------------------------------------------------*/

TEST_F(TestDirection, getDirectionExtended )
{

    TestDescription( "getDirection (J2000) with selected data. uvw available" );

    // Use DefaultMS as a simple sequence.
    // Use the below to show uv valuses from MS.

    const String MsName = env.getCasaMasterPath() + "listobs/uid___X02_X3d737_X1_01_small.ms";

    // Create Object //
    
        MeasurementSet ms0( MsName.c_str() );
    
        PointingDirectionCalculator calc(ms0);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );


    //+
    // * Option *
    // selectData()   
    //  Please change if() to true,
    //-
        if(true)
        {
            const String AntSel = "DV01&&DV02";
            calc.selectData( AntSel,  "","","","","","","","","" );

            ExpectedNrow = calc.getNrowForSelectedMS();
        }

    //+
    // setDirectionColumn() 
    //-
  
        String ColName = "DIRECTION"; 
        Description( "getDirectionColumn()", ColName  );

        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

    //+
    //      MatrixShape (COLUMN_MAJOR) 
    //-
        Description("calling setDirectionListMatrixShape()" ,"Column Major" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    //+
    //  setMovingSourceDirection() 
    //-
        String src = "MOON";
        Description("calling setMovingSource()", src);
         EXPECT_NO_THROW( calc.setMovingSource( src ) ); 

    //+
    //  getDirection()
    //-

        Description("calling  getDirection() ","" );

         Matrix<Double>  DirList1  = calc.getDirection();
         uInt  n_row    = DirList1.nrow();

         printf( "Number of Row = %d \n", n_row );
         EXPECT_EQ( n_row, ExpectedNrow);


    //+
    // uv value 
    //  (first attach to Column e )
    //-

        casacore::ROArrayColumn<casacore::Double> uvwColumn;	
        uvwColumn .attach( ms0 , "UVW");

    //+
    // Dump Matrix
    //-

#if 1

    Description("Dump obtined Direction info. ","" );

    for (uInt row=0; row< n_row; row++)
    {
        // (u,v,w) //

        uInt RowId = calc.getRowId(row); // point the row in Original. //
        casacore::Vector<Double>  val_uvw = uvwColumn.get(RowId);

        Double u = val_uvw[0];
        Double v = val_uvw[1];
        Double w = val_uvw[2];

        // Direction //

        Double Val_1 = DirList1(row,0);
        Double Val_2 = DirList1(row,1);

        casacore::MDirection  MovDir  = calc.getMovingSourceDirection();
        String strMovDir = MovDir.toString();

        printf(    "Dir at, %d, %f,%f, [Mov:%s], uv=, %f,%f,%f  \n",  
                row, Val_1, Val_2, strMovDir.c_str(), 
                u, v, w );

    }

#endif 

}
/*------------------------------------------
   selectData ( <various types of keys>  )

   - Inspect that expected number of selected
     data comes out.
  ------------------------------------------*/

class TestSelectData : public BaseClass
{

public:

protected:
 
        String  AntSel              = "";      // [C++11] Should be Null, to go throgh GetTEN 
        String  SpwSel              = "";     
        String  FieldSel            = "";
        String  TimeSel             = "";
        String  ScanSel             = "";
        String  FeedSel             = "";
        String  IntentSel           = "";
        String  ObservationSel      = "";
        String  UVRangeSel          = "";
        String  MSSelect            = "";

        String DefMsName = "listobs/uid___X02_X3d737_X1_01_small.ms";

        TestSelectData()
        {
        }

        ~TestSelectData()
        {
        }


        virtual void SetUp()
        {
            BaseClass::SetUp();
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
       }

        void test_selectdata(PointingDirectionCalculator & calc);

};

/*------------------------------------------------------
   Execution of selectData(...)
    - args are given from external variables
    - Note that 'calc' MUST be given by reference.
  ------------------------------------------------------*/

void TestSelectData::test_selectdata(PointingDirectionCalculator& calc)
{

      calc.selectData( AntSel,
                         SpwSel,
                         FieldSel,
                         TimeSel,
                         ScanSel,
                         FeedSel,
                         IntentSel,
                         ObservationSel,
                         UVRangeSel,
                         MSSelect
                );

}

/*------------------------------------------
   selectData( Antenna )
   - Inspect the result with Antenna key. 
  ------------------------------------------*/

TEST_F(TestSelectData, Antenna )
{
    TestDescription( "selectData (key=Antenna)" );

    // MS name for this Test //
   
       String name = env.getCasaMasterPath() + DefMsName;
       printf( " Used MS is [%s] \n", name.c_str() );
  
    // Create Object //
    
       MeasurementSet ms( name.c_str() );
  
       PointingDirectionCalculator calc(ms);

    // Initial brief Inspection //
   
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n"); 
       ExpectedNrow = calc.getNrowForSelectedMS();
       printf("Expected nrow =%d\n",ExpectedNrow);

       EXPECT_NE( (uInt)0, ExpectedNrow );
       
    //+
    // TEST SET 
    //  1) Specify the condition
    //  2) probe and confirm Exception
    //  3) Result check , expected Nrow was selected 
    //-

    uInt  nrow ;

      AntSel = "";
        FunctionalDescription("Antenna: by NULL  (Matches) ",AntSel);

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow); 
        EXPECT_EQ (ExpectedNrow, nrow);     // see MS in detail //

      AntSel = "hoge&&&";
        FunctionalDescription("Testing Abnormal Name = hoge (No Matches))",AntSel );

        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ (ExpectedNrow, nrow);

      AntSel = "DV01&&&";
        FunctionalDescription("Antenna: Normal specific Name. (Matches) ",AntSel);

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ ((uInt)180, nrow);     // see MS in detail //

      AntSel = "DV02&&&";
        FunctionalDescription("Antenna: Normal specific Name (No Matches)",AntSel );

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ ((uInt)180, nrow);

      AntSel = "PM03&&&";
        FunctionalDescription("Antenna: Normal specific Name (No Matches)",AntSel );

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ ((uInt)180, nrow);
 
      AntSel = "DV*&&&";

        FunctionalDescription("Antenna: Normal with Wild card char. (Matches)",AntSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        printf("=> Calling getNrowForSelectedMS() after oelectData() called.\n");
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ ( (uInt)360, nrow);

        //+
        // Other detail cases are to be programmed here if needed. 
        //-

        return;
        
}

/*------------------------------------------
   selectData( Spw )
   - Inspect the result with Spw. 
  ------------------------------------------*/

TEST_F(TestSelectData, Spw )
{
    TestDescription( "selectData (key=Spw)" );

    // MS name for this Test //

        String name = env.getCasaMasterPath() + DefMsName;
        printf( " Used MS is [%s] \n", name.c_str() );

    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      SpwSel = "";
        FunctionalDescription( "Spw:Nothig specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

    
      SpwSel = "*";
        FunctionalDescription( "Spw: Wildcard *  specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      SpwSel = "hoge";
        FunctionalDescription( "Spw: abnormal letters  specified.",SpwSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      SpwSel = "0:13~20";
        FunctionalDescription( "Spw: spw=0, ch=13~20 specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ( (uInt)270, nrow);
#if 0
//
// Once execution OK, the next Error query makes unexpected Return ?
//
      SpwSel = "0:13=20";
        Description( "Spw: spw=0, ch=13~20 specified.(Syntax ERROR)",SpwSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

#endif 

      SpwSel = "1:13~15";
        FunctionalDescription( "Spw: spw=1, ch=13~20 specified.(None)",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ((uInt)810, nrow);


}

/*------------------------------------------
   selectData( Field )
   - Inspect the result with Field. 
  ------------------------------------------*/

TEST_F(TestSelectData, Field )
{ 
    TestDescription( "selectData (key=Field)" );

    // Using MS //
    
        const String MsName = "sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
      printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE( (uInt)0, ExpectedNrow );

      uInt nrow;
      FieldSel = "";
        FunctionalDescription( "Files:Nothig specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      FieldSel = "*";
        FunctionalDescription( "Field: Wildcard *  specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      FieldSel = "hoge";
        FunctionalDescription( "Fieled: abnormal letters  specified.",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      FieldSel = "0";  // Field ID 
        FunctionalDescription( "Field: ID=0 specified.(No exits)",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ((uInt)0, nrow);       // On-Table , None on MAIN .. //

      FieldSel = "1";  // Field ID 
        FunctionalDescription( "Field: ID=1 specified.(exits)",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      FieldSel = "9";  // Out of Range
        FunctionalDescription( "Field: ID=9 Not on the FIELD TABLE",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

}

/*------------------------------------------
   selectData( Time )
   - Inspect the result with Time range. 
  ------------------------------------------*/

TEST_F(TestSelectData, Time )
{ 
    TestDescription( "selectData (key=Time)" );

    // Using MS //
    
        const String MsName = "sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
      printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      TimeSel = "";
        FunctionalDescription( "Time: Nothig specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      TimeSel = "hoge";
        FunctionalDescription( "Time: abnormal letters  specified.",TimeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      TimeSel = ">1900/01/01/00:00:00";
        FunctionalDescription( "Time: since 1900/01/01/00:00:00  specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      TimeSel = ">2018/01/01/00:00:00";
        FunctionalDescription( "Time: since 2018/01/01/00:00:00  specified. None matches",TimeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ((uInt)0, nrow);

      TimeSel = "<2014/12/4/00:39:25";  //  Time Condition  
        FunctionalDescription( "Field: a specifc Time, Limited Number  match are expected",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        printf( "# Actually detected Nrow = %d \n" , nrow); 
        EXPECT_GT ((uInt)20, nrow);      // On-Table  //
        EXPECT_EQ ((uInt)5,   nrow);    

      TimeSel = "2014/12/4/00:40:00~2014/12/4/00:40:10";  //  Combined  10 sec 
        FunctionalDescription( "Field: a specifc Time, Limited Number  match are expected",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        printf( "# Actually detected Nrow = %d \n" , nrow);
        EXPECT_GT ((uInt)20, nrow);      // On-Table  //
        EXPECT_EQ ((uInt)17, nrow);    // SEE ACTUAL COUNT on Browser //

}

/*------------------------------------------
   selectData( Feed )
   - Inspect the result with Feed key. 
  ------------------------------------------*/

TEST_F(TestSelectData, Feed )
{ 
    TestDescription( "selectData (key=Feed)" );

    // Using MS //
    
        const String MsName = "/sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //

       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      FeedSel = "";
        FunctionalDescription( "Feed: Nothig specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      FeedSel = "0";    // NO DATA /
        FunctionalDescription( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

      FeedSel = "1";    // NO DATA /
        FunctionalDescription( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

}

/*------------------------------------------
   selectData( Intent )
   - Inspect the result with Intent key. 
  ------------------------------------------*/

TEST_F(TestSelectData, Intent )
{ 
    TestDescription( "selectData (key=Intent)" );

    // Using MS //
    
        const String MsName = "sdimaging/selection_intent.ms";    // 
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //

       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      IntentSel = "";
        FunctionalDescription( "Intent: Nothig specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);
 
       // Usage : In what way this term is used. 
        //   and where to exist in MS.
          
      IntentSel = "*HOGE*";     // 
        FunctionalDescription( "Scan: ID specified.",IntentSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

      IntentSel = "*CAL*";      // 
        FunctionalDescription( "Scan: ID specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)512, nrow);

      IntentSel = "*BAND*";      // 
        FunctionalDescription( "Scan: ID specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)512, nrow);

}

/*------------------------------------------
   selectData( Observation )
   - Inspect the result with Observation info. 
  ------------------------------------------*/

TEST_F(TestSelectData, Observation )
{ 
    TestDescription( "selectData (key=Observation)" );

    // Using MS //
    
        const String MsName = "/sdimaging/selection_spw.ms";    //    Three Observation entries.
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      ObservationSel = "";
        FunctionalDescription( "Observation: Nothig specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

     ObservationSel = "hoge";     
        FunctionalDescription( "Observation: abnormal expr..",ObservationSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

      ObservationSel = "0";     // one DATA /
        FunctionalDescription( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)256, nrow);

      ObservationSel = "1";    // second DATA _
        FunctionalDescription( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uint)256, nrow);

      ObservationSel = "2";    // 3rd.Data 
        FunctionalDescription( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

         EXPECT_EQ ((uint)225, nrow);

      ObservationSel = "9";    // No Data (err) 
        FunctionalDescription( "Observation: ID specified.",ObservationSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)0, nrow);

}

/*------------------------------------------
   selectData( UVRange )
   - Inspect the result with UVRange. 
  ------------------------------------------*/

TEST_F(TestSelectData, UVRange )
{ 
    TestDescription( "selectData (key=UV Range)" );

    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + DefMsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      UVRangeSel = "";
        FunctionalDescription( "UVrange: Nothig specified.",UVRangeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      UVRangeSel = "hoge";     
        FunctionalDescription( "UVrange: abnormal expr..",UVRangeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

      UVRangeSel = ">1.0lambda";        //  Exprecasacore::Table::TableOption:: Updatession Unknown...../
        FunctionalDescription( "UVrange: ID specified.",UVRangeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)540, nrow);



}

/*------------------------------------------
   selectData( MS select ) - TaQL query -
   - Inspect the result selected by TaQL form. 
  ------------------------------------------*/

TEST_F(TestSelectData, MSselect )
{ 
    TestDescription( "selectData (key=MS Select)" );

    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + DefMsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

      uInt nrow;
      MSSelect = "";
        FunctionalDescription( "MSselect: Nothig specified.",MSSelect);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);

      MSSelect = "hoge";     
        FunctionalDescription( "MSselect: abnormal expr..",MSSelect);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, nrow);

      MSSelect = "FIELD_ID==0";        //  Query description is here. //

        FunctionalDescription( "MSselect: ID specified.",MSSelect);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)600, nrow);

      MSSelect = "FIELD_ID==1";        //  Query description is here. //

        FunctionalDescription( "MSselect: ID specified.",MSSelect);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)360, nrow);

      MSSelect = "FIELD_ID==2";        //  Query description is here. //

        FunctionalDescription( "MSselect: ID specified.",MSSelect);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)120, nrow);

      MSSelect = "FIELD_ID==3";        //  ERROR :Query description is here. //

        FunctionalDescription( "MSselect: ID specified.",MSSelect);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ ((uInt)0, nrow);
}


class TestSetFrame : public BaseClass
{

public:
typedef struct _FrameName {
    bool   available;
    String name;
} FrameTypeList;

protected:

const std::vector<FrameTypeList> DefinedFrametypes
    = {
        {true,  "J2000"},  
        {true,  "JMEAN"}, 
        {true,  "JTRUE"}, 
        {true,  "APP"}, 
        {true,  "B1950"}, 
        {true,  "B1950_VLA"}, 
        {true,  "BMEAN"}, 
        {true,  "BTRUE"}, 
        {true,  "GALACTIC"}, 
        {true,  "HADEC"}, 
        {true,  "AZEL"}, 
        {true,  "AZELSW"}, 
        {true,  "AZELGEO"}, 
        {true,  "AZELSWGEO"}, 
        {true,  "JNAT"}, 
        {true,  "ECLIPTIC"}, 
        {true,  "MECLIPTIC"}, 
        {true,  "TECLIPTIC"}, 
        {true,  "SUPERGAL"}, 
        {true,  "ITRF"}, 
        {true,  "TOPO"}, 
        {true,  "ICRS"}, 
        {false, "N_Types"}, 
        {true,  "MERCURY"}, 
        {true,  "VENUS"}, 
        {true,  "MARS"}, 
        {true,  "JUPITER"}, 
        {true,  "SATURN"}, 
        {true,  "URANUS"}, 
        {true,  "NEPTUNE"}, 
        {true,  "PLUTO"}, 
        {true,  "SUN"}, 
        {true,  "MOON"}, 
        {true,  "COMET"}, 
        {false,  "N_Planets"}, 
        {false,  "EXTRA"}, 
        {false,  "DEFAULT"}, 
        {false,  "AZELNE"}, 
        {false,  "AZELNEGEO"},

        {false,  "Hoge"},               // Bad name //

    };

        TestSetFrame () { }

        ~TestSetFrame () { }

        virtual void SetUp()
        {
            BaseClass::SetUp();
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
       }
      
        // Test Fixture Sub //

        void check_direction_info(PointingDirectionCalculator& calc, uInt n_frame );
};

/*-----------------------------------------------
  check_direction_into()

   1) setFrame( given name by ARG )
   2) getDirectionType()
   - 1)and 2) MUST BE SAME.
   3) output (converted to string) must be same 

  Internally use: MDirection getDirectionType()
  -----------------------------------------------*/

void TestSetFrame::check_direction_info(PointingDirectionCalculator& calc, uInt n_frame )

{
    //+
    // Test setFrame (Str) 
    //   No Exception is expected all the time 
    //-

      EXPECT_NO_THROW( calc.setFrame( DefinedFrametypes[n_frame].name ));

    //+
    // Durection Type 
    //  converting to string form to compare.
    //  output text must contain the specified KEY.
    //-  

         // Direction Type //
  
       casacore::MDirection DirType  = calc.getDirectionType();
 
          printf( "#   MDirection: [%s] \n",  DirType.toString().c_str()  ) ;
          printf( "#   Given String [%s] \n", DefinedFrametypes[n_frame].name.c_str() );

        
        String converted = DirType.toString();
        String sub_str   = DefinedFrametypes[n_frame].name;

        // GTEST :: Check SubString // 
        
    Description( "Checking frame sub-string. ",DirType.toString().c_str() );

    // Some of them are not supported and throw Exception //

    if(  DefinedFrametypes[n_frame].available == true)       
        EXPECT_TRUE( converted.find(sub_str) !=std::string::npos); 
    else
        EXPECT_FALSE( converted.find(sub_str) !=std::string::npos);

}

/*-----------------------------------------------
  setFrama()
  - Strings are defined above.
  - Some keywords intentionally  makes Exception. 
  -----------------------------------------------*/

TEST_F(TestSetFrame, setFrame )
{ 
    TestDescription( "setFrame (String FrameName)" );

    // Using MS //
    
        const String MsName = "listobs/uid___X02_X3d737_X1_01_small.ms";    
    
    // MS name for this Test //

        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );

    // Create Object //
    
        MeasurementSet ms( name.c_str() );
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
        printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
        ExpectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, ExpectedNrow );

    // Various Frame Type (String) //

        for(uInt k = 0; k < DefinedFrametypes.size() ; k++  )
        {
            FunctionalDescription( "setFrame(Rsved Name) ", DefinedFrametypes[k].name.c_str() );

            // Execute and check Exception and other requirements//
        
            check_direction_info( calc, k ) ;
        }
}

class MyDebug : public BaseClass
{

public:

protected:

        MyDebug () { }

        ~MyDebug () { }

        virtual void SetUp()
        {
            BaseClass::SetUp();
            CopyDefaultMStoWork();
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
            DeleteWorkingMS();
       }
      
        // Test Fixture Sub //

};

//+
// DEBUG THIS UNIT TEST 
//-

// Construcor and Column check//
TEST_F(MyDebug, Debug1_Col )
{

    MSNameList   ms;

    for (uInt m=0; m< ms.count(); m++ ) { 
        if(ms.isExceptionActivated(m) == true ) continue ; // Ignore un-available ones. //

        String name = env.getCasaMasterPath() + ms.getName(m);

        PointingTableAccess   pta1(name);
        printf("======%d  %s ======= \n", m, ms.getName(m).c_str() );

        printf("check Column[Direction]      =%d \n", pta1.checkColumn("Direction"));
        printf("check Column[Target]         =%d \n", pta1.checkColumn("TarGet"));
        printf("check Column[PointingOffset] =%d \n", pta1.checkColumn("Pointing_Offset"));
        printf("check Column[SourceOffset]   =%d \n", pta1.checkColumn("Source_Offset"));
        printf("check Column[Encoder]        =%d \n", pta1.checkColumn("ENcoder"));

    }

}

// Construcor and Column check//
void debug1_sub(const String& name )
{
    PointingTableAccess   poAcc(name);
    AntennaTableAccess    atAcc(name);
}
// Object alloc/free Test , 1000 times. // 
TEST_F(MyDebug, Debug1_mem1 )
{
    const String MsName = "./sdimaging-t.ms";
    String name = /* env.getCasaMasterPath()+ */  MsName;

    for (uInt m=0; m<1000;m++ ) {

        debug1_sub(name);
    }
}

// Construcor and Column check//
void debug2_sub(MeasurementSet &ms )
{
     PointingDirectionCalculator calc(ms);
}
// Object alloc/free Test , 1000 times. // 
TEST_F(MyDebug, Debug1_mem2 )
{
    const String MsName = "./sdimaging-t.ms";
    String name = /* env.getCasaMasterPath()+ */  MsName;

    MeasurementSet ms0(name);
    for (uInt m=0; m<100;m++ ) {

        debug2_sub(ms0);
    }
}

// Copied MS for modify access //
TEST_F(MyDebug, Debug2_checkColumn )
{
    const String MsName = "./sdimaging-t.ms";
    String name =  MsName;

    PointingTableAccess   pta2(name,true);

    printf("check Column[Direction]      =%d \n", pta2.checkColumn("Direction"));
    printf("check Column[Target]         =%d \n", pta2.checkColumn("TarGet"));
    printf("check Column[PointingOffset] =%d \n", pta2.checkColumn("Pointing_Offset"));
    printf("check Column[SourceOffset]   =%d \n", pta2.checkColumn("Source_Offset"));
    printf("check Column[Encoder]        =%d \n", pta2.checkColumn("ENcoder"));

    pta2.duplicateColumns();

    printf("check Column[Direction]      =%d \n", pta2.checkColumn("Direction"));
    printf("check Column[Target]         =%d \n", pta2.checkColumn("TarGet"));
    printf("check Column[PointingOffset] =%d \n", pta2.checkColumn("Pointing_Offset"));
    printf("check Column[SourceOffset]   =%d \n", pta2.checkColumn("Source_Offset"));
    printf("check Column[Encoder]        =%d \n", pta2.checkColumn("ENcoder"));

// Row operation //
}

// Pointing Row
TEST_F(MyDebug, Debug3_Row )
{
    const String MsName = "./sdimaging-t.ms";
    PointingTableAccess   pta(MsName,true);

    printf( "1) nrow = %u\n", pta.getNrow() );
    printf( "2) Adding 3 row \n" );

    pta.appendRow(3);

    printf( "3) nrow = %u\n", pta.getNrow() );

    printf( "4) Removing 3 rows \n" );
    
    pta.removeRow(3843);
    pta.removeRow(3843);
    pta.removeRow(3843);
   
    printf( "5) nrow = %u\n", pta.getNrow() );

}

TEST_F(MyDebug, Debug4_AddColumn)
{
    const String MsName = "./sdimaging-t.ms";
    PointingTableAccess   pta2(MsName,true);

    pta2.duplicateColumns();

    printf("check Column[Direction]      =%d \n", pta2.checkColumn("Direction"));
    printf("check Column[Target]         =%d \n", pta2.checkColumn("TarGet"));
    printf("check Column[PointingOffset] =%d \n", pta2.checkColumn("Pointing_Offset"));
    printf("check Column[SourceOffset]   =%d \n", pta2.checkColumn("Source_Offset"));
    printf("check Column[Encoder]        =%d \n", pta2.checkColumn("ENcoder"));
    printf("check Column[HOGEHOGE]       =%d \n", pta2.checkColumn("HogeHOGE"));

}

// Pointing Dump //
TEST_F(MyDebug, Debug5_Dump )
{
    const String MsName = "./sdimaging-t.ms";
    PointingTableAccess   pta(MsName,false);

    uint nrow = pta.getNrow();
    for (uInt row=0; row<nrow; row++)
    {
        Double time              = pta.getTime(row);
        Double interval          = pta.getInterval(row);
        Vector<Double> direction = pta.getDirection(row);
        Vector<Double> target    = pta.getTarget(row);

        printf( "%u, Time=,%f,Interval=,%f, Dir=,%f,%f, Target=,%f,%f     \n",
                 row, time, interval,direction[0], direction[1], target[0], target[1] );

    }
}

// Antenna: Row //
TEST_F(MyDebug, Debug10_Row )
{
    const String MsName = "./sdimaging-t.ms";
    AntennaTableAccess   ata(MsName,true);

    printf( "nrow = %u\n", ata.getNrow() );

    printf( "Adding 3 row \n" );

    ata.appendRow(3);

    printf( "nrow = %u\n", ata.getNrow() );

    printf( "Removing 3 rows \n" );
    
    ata.removeRow(3);
    ata.removeRow(2);
    ata.removeRow(1);
   
    printf( "nrow = %u\n", ata.getNrow() );
}

// Antenna WRite and Dump //
TEST_F(MyDebug, Debug11_Dump )
{
    const String MsName = "./sdimaging-t.ms";

    AntennaTableAccess   ata(MsName,true);
    uint nrow = ata.getNrow();

    for (uInt row=0; row<nrow; row++)
    {
        ANTENNADataBuff data = ata.getRowData(row);
        String  name             = data.name;
        String station           = data.station;
        Vector<Double>  position = data.position;
        Vector<Double>  offset   = data.offset;

        printf( "AntID=%u, Name=[%s], Pos=(%f,%f,%f)  \n",
                 row, name.c_str(), position[0],position[1],position[2]  );

    }

    printf( " Writing 2 rows \n" );

    ANTENNADataBuff write_data;
      write_data.name            =  "...";
      write_data.station         =  "MyStation";
      write_data.type            =  "KUMIKO-based";
      write_data.mount           =  "MyMount";

    ata.appendRow(2);

    write_data.name = "DA01";
    write_data.position.resize(3);
    write_data.offset.resize(3);
    write_data.position[0] = 11.0;
    write_data.position[1] = 11.1;
    write_data.position[2] = 11.2;

      ata.putRowData(1, write_data);

    write_data.name = "DA02";

    write_data.position[0] = 222.0;
    write_data.position[1] = 222.1;
    write_data.position[2] = 222.2;

      ata.putRowData(2, write_data);


     nrow = ata.getNrow();
     for (uInt row=0; row<nrow; row++)
     {
     
         ANTENNADataBuff data = ata.getRowData(row);
         String  name             = data.name;
         String station           = data.station;
         Vector<Double>  position = data.position;
         Vector<Double>  offset   = data.offset;

         printf( "AntID=%u, Name=[%s], Pos=(%f,%f,%f)  \n",
                  row, name.c_str(), position[0],position[1],position[2]  );
  
     } 

}

//+
// Antenna List. 
//-
TEST_F(MyDebug, Debug12_Dump )
{
    const String MsName = "./sdimaging-t.ms";

    std::vector<String>   MsList = {
      "sdimaging/Uranus1.cal.Ant0.spw34.ms",
      "sdimaging/Uranus2.cal.Ant0.spw34.ms"   };

    for (uInt m=0; m<MsList.size(); m++)
    {
        String name = env.getCasaMasterPath() + MsList[m];
        AntennaTableAccess   ata(name);
        uint nrow = ata.getNrow();

        printf("Listing File[%s] \n", name.c_str());
        for (uInt row=0; row<nrow; row++)
        {
            ANTENNADataBuff data = ata.getRowData(row);
            String  name             = data.name;
            String station           = data.station;
            Vector<Double>  position = data.position;
            Vector<Double>  offset   = data.offset;

            printf( " - AntID=%u, Name=[%s], Pos=(%f,%f,%f)  \n",
                     row, name.c_str(), position[0],position[1],position[2]  );
        }
    }

}

TEST_F(MyDebug, SplineConstructor )
{
    const String MsName = "./sdimaging-t.ms";
    String name = /* env.getCasaMasterPath()+ */  MsName;

    MeasurementSet ms0(name);
    PointingDirectionCalculator calc(ms0); 

    casa::SplineInterpolation  *spH =  calc.getCurrentSplineObj();

    calc.setDirectionColumn("DIRECTION" );
      spH =  calc.getCurrentSplineObj();
    calc.setDirectionColumn("DIRECTION" );  
      spH =  calc.getCurrentSplineObj();
    calc.setDirectionColumn("TARGET" );
      spH =  calc.getCurrentSplineObj();
    calc.setDirectionColumn("TARGET" );
      spH =  calc.getCurrentSplineObj();

}

}  // END namespace

/************************************************
   Unit Test Main (google test)
    - Based on instructed Template for GTest.
    - Such minimum statements are recommended.
(Note) Interpolation Test completed 
(History) 
-  12/7 Merged master (to get new CMakeList)
-  12/7 Added initialize (CAS-12114,old 11427-2)
 *************************************************/

int main (int nArgs, char * args [])
 {

   // Initialize //

    ::testing::InitGoogleTest(& nArgs, args);

   // Run Test //

    return (RUN_ALL_TESTS()) ;

}

