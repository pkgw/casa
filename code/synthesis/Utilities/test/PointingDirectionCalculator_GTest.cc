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
    typedef struct _MSDef   {
        bool   ExThrow;  // True = cause Exeption
        String name;     // MS name, with relative path from CASAPATH
    } MSDef;

    // Get File Name by Number //
    const String  name(uInt No ) {
        const String msg = "Internal Bugcheck:";
        if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());
        return TestMSList[No].name;
    }

    // True is this access makes Exception 
    bool isException(uInt No) { 
        const String msg = "Internal Bugcheck:";
        if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());
        return TestMSList[No].ExThrow;
   }
    uInt count() { return TestMSList.size();  }

private:

    uInt  pos =0;  
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
 
    };
};

//+
// CAS-8418 added
// Direction Column List
//-
class PointingColumnList 
{
public:
        size_t size() { return DirList.size();}
        string name(uint n) { return DirList[n]; }
private:
        std::vector<string> DirList 
         = {"DIRECTION", "TARGET", "POINTING_OFFSET", "SOURCE_OFFSET", "ENCODER" };
};


//***************************
// DEBUG flags
//***************************

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
    
    String GetCasaPath(const String &pathname )
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
// Programme Wide constant 
//**************************************

const String  DefaultLocalMsName =  "./sdimaging-t.ms";

//************************************** 
// Base TestClass 
//  for TEST FIXTURE
//**************************************
class BaseClass : public ::testing::Test
{
public:
     // Running Environment //
     RunEnv       env;

     //+
     // MS copy/delete to make Test-MS
     //-
     void CopyDefaultMStoWork();
     void DeleteWorkingMS();

     //+
     // Console Message of test progress
     //-

     void TestDescription( const String &Title );
     void FunctionalDescription(const String &Title, const String &Param);
     void Description(const String &Title, const String &Param);

    uInt    expectedNrow = 0;   // C++11 feature //

    BaseClass()  { }

    ~BaseClass() { }

    void SetUp() { }
    void TearDown() { }

private:

     // Test MS copy and delete flag //
     bool fgCopyMS    = true;
     bool fgDeleteMS  = false;
 
};

//+
//  Log Title Output Functions for readable text 
//  of this UT.
//-
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
//   This is for Some testing items which must contain planned Data in the MS.
//   After copying, This UT programme moify the MS for each purpose.
//***************************************************************************
void BaseClass::CopyDefaultMStoWork()
{
    // Src/Dst Path (string) 
        const String src = env.getCasaMasterPath() + "sdimaging/sdimaging.ms";
        const String dst = DefaultLocalMsName;

    // Src/Dst Path (Path) 
        casacore::Path        sourcePath(src);
        casacore::Path        targetPath(dst);       
        casacore::Directory   dir_ctrl(sourcePath);

        printf( "Copying Default MeasurementSet for modifed use to [%s] \n" ,dst.c_str() );
        printf( " - src filespec  : %s \n", src.c_str() );
        printf( " - dest filespec : %s \n", dst.c_str() );

    // Copy File   //
    if(fgCopyMS)
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

    printf( "Deleting Working MeasurementSet[%s]. \n" ,dst.c_str() );

    // Delete File (Recursively done) 
    if (fgDeleteMS)
    {
         dir_ctrl. removeRecursive(false /*keepDir=False */ );
    }
}

//*****************************************************
// Interpolation Testing Trajectory lcass
// (implemented by Singleton) 
//****************************************************

class TrajectoryFunction
{
    typedef void (*FUNCTYPE)(Double, Double&, Double& );

public:

    typedef enum _Tr_ {
        Simple_linear,        // 0
        Normalized_linear,    // 1
        Sinusoid_slow,       // 2
        Sinusoid_quick,      // 3
        Sinusoid_hasty,      // 4
        Harmonics_sinusoid,  // 5
        Gauss,               // 6   (new 12/11)
        Zero,                // 7
        Const                // 8
    } Type;

    size_t size() { return fpTrajectoryfunc.size(); }
    
    // calculation (default) //
    void calc(Double r_time, Double& X, Double& Y) {
           (*fpTrajectoryfunc[trajFunctionNo])( r_time, X, Y ); return; }

    // calculation (specified) //
    void calc(Double r_time, Double& X, Double& Y, uInt Fno) {
           (*fpTrajectoryfunc[Fno])( r_time, X, Y ); return; }

    // set type //
    void setType(uInt no ) {
           if( no < fpTrajectoryfunc.size() )  trajFunctionNo = no;}

    // get instance (first call) //
    static TrajectoryFunction &getInstance() { static TrajectoryFunction inst; return inst; }


private:

    TrajectoryFunction() {}
    ~TrajectoryFunction(){}
    TrajectoryFunction(const TrajectoryFunction&);
    TrajectoryFunction& operator=(const TrajectoryFunction&);

    // Selected Function //
     uInt trajFunctionNo =0;

static void Function_SimpleLinear( Double r_time, Double &X, Double &Y )
{
    X = -1.0 + 2.0 * r_time;
    Y = -1.0 + 2.0 * r_time;
    return;
} 

static void Function_NormalizedLinear( Double r_time, Double &X, Double &Y )
{
    // Normalized time :: | Rel_Time | < 1.0 , this case [0,1.0] is used //

    X =  (r_time * 2.0 - 1.0 ) * M_PI ;
    Y =  (r_time * 2.0 - 1.0 ) * (M_PI / 2.0);
    return;
}

static void Function_sinusoid_slow( Double r_time, Double& X, Double& Y)
{
 
    X = 1.0 * cos( 2.0*M_PI  * r_time );
    Y = 1.0 * sin( 2.0*M_PI  * r_time );
    return;
}

static void Function_sinusoid_quick( Double r_time, Double& X, Double& Y)
{   
    
    X = 2.0 * cos( 10* 2.0*M_PI  * r_time );
    Y = 1.0 * sin( 10*  2.0*M_PI  * r_time );
    return;
}   

static void Function_sinusoid_hasty( Double r_time, Double& X, Double& Y)
{
    double FREQ= 100.0; 
    X = 2.0 * cos( FREQ*  2.0*M_PI  * r_time );
    Y = 1.0 * sin( FREQ*  2.0*M_PI  * r_time );
    return;
}

static void Function_harmonics_sinusoid( Double r_time, Double& X, Double& Y)
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

static void Function_gauss( Double r_time, Double& X, Double& Y)
{

    Double t   = r_time - 0.5;
    Double A = 50;

    Double gauss  = exp (-A*t*t);

    X = gauss;
    Y = gauss;
    return;
}

static void Function_zero(Double r_time, Double& X, Double& Y)
{
    X = 0.0 + 0.0*r_time;
    Y = 0.0 + 0.0*r_time;
    return;
}
static void Function_const(Double r_time, Double& X, Double& Y)
{
    X = 1.0 + 0.0*r_time;
    Y = 1.0 + 0.0*r_time;
    return;
}

    // Function Table //
    std::vector<FUNCTYPE>  fpTrajectoryfunc
    {    
        Function_SimpleLinear,        // 0
        Function_NormalizedLinear,    // 1
        Function_sinusoid_slow,       // 2
        Function_sinusoid_quick,      // 3
        Function_sinusoid_hasty,      // 4
        Function_harmonics_sinusoid,  // 5
        Function_gauss,               // 6   (new 12/11)
        Function_zero,                 // 7
        Function_const                 // 8
    };   

 

}; // end class
//********************************************************
//  INTERPOLATION  Generation 
//   for Interporation Verification TEST 
//  - Generate testing trajectry in Direction
//  -  Interval(POINTING, MAIN)  (tunable in initialize() ).
//********************************************************

class TuneMSParam  /*: public BaseClass */ 
{
public:

    typedef struct _Pointing_ 
    {
        Double time;
        Double interval;
        Double relativeTime;
    
        std::pair<Double,Double> position[5];

    } PseudoPointingData;
 
    TuneMSParam() { printf("TuneMSParam constructor \n" ); /* Minimum SetUp only */  }

    // Interval Time //

        Double  getPointingTableInterval() { return pointingIntervalSec_;}
        Double  getMainTableInterval()     { return mainIntervalSec_;}

    // Calculation Error Limit //
    
        void setInterpolationErrorLimit(Double val)
             { printf("TuneMSParam:: set interpolation err. limit = %e \n", val );
               errorLimit_ = defaultInterpolationErrorLimit_ = val; }

    // Init and Define Parameters 

        void Initialize();
        void Initialize( Double, Double);
        void setMainRowCount( uint n ) { requiredMainTestingRow_ =  defultMainTestingRow       = n; }  

    // Offset between POINTNG and MAIN //

        /* function reserved */
 
    // Pseudo Trace(Direction) for the Test

        PseudoPointingData     pseudoPointingInfoPointing2(Double tn);
        PseudoPointingData     pseudoPointingInfoMain2    (Double tn);

    // available POINTING TABLE count //

       uInt getAvailablePointingTestingRow() { return availableNrowInPointing_; }

    // required MAIN TABLE count 

       uInt getRequiredMainTestingRow()      { return requiredMainTestingRow_; }

    //+
    // Numerical Error Statictic 
    //-
 
        Double getInterpolationErrorLimit() { return errorLimit_;  } ;

    // Row count (to add) 

        uInt getAddInerpolationTestPointingTableRow() {return extraNrowInPointing_; };
        uInt getAddInerpolationTestMainTableRow()     {return extraNrowInMain; };
	
    // Antenna Count //

        uInt getMaxOfAntenna() { return prepareMaxAntenna_; }
        void setMaxAntenna(uInt n) { prepareMaxAntenna_ = n;  }
   

    // Special TEST (new : 8-MAr-29190)

        bool ifCoeffLocTest() { return fgCoeffLocationTest; }
        void setCoeffLocTest(bool val) { fgCoeffLocationTest = val; }

private:
        // initialize  Commmon //
          void init();

        //+
        // create pseudo Pointing Info by--rajectory-function 
        // return value contains various info 
        //  => see PseudoPointing
        //-

          PseudoPointingData        pseudoPointingBaseInfo(Double deltaTime);

        //+
        //  Relative Time (r_time) and Total Time
        //   for Trajectory Function
        //-

          Double  Interval  = 0.0; // Initilal ..
          uInt    nRow      = 0;
          Double  r_time     = 0.0;

       // Pre-located row , use tables with extended. See MS (sdimaging.ms) by tool //

        uInt defultMainTestingRow       = 5040; 

        const uInt defInerpolationTestPointingTableRow_   = 3843;
        const uInt defInerpolationTestMainTableRow_       = 3843;

        // Row Count to execute //
        uInt requiredMainTestingRow_     = 0;    //   MUST BE SET 

        Double  requiredNrowInPointing_ = 0;    //   internally calculated
        Double  availableNrowInPointing_ = 0;    //   internally calculated   

        Double  extraNrowInPointing_ ;   // calculated when start
        Double extraNrowInMain ;       // calculated when start

       // Error Limit (threshold) in GoogleTest Macro //

        Double errorLimit_ ;
        Double defaultInterpolationErrorLimit_ = 1.0e-06 ;


    // Interval Second.

         Double pointingIntervalSec_ =0.0;            // Interval Time to set in POINTING 
         Double mainIntervalSec_     =0.0;            // Interval Time to set in MAIN 

    // Number of Antenna (CAS-8418)  

        uInt prepareMaxAntenna_ = 3;
        bool fgCoeffLocationTest  = false; 

};

//+
// Initialize
//-

void TuneMSParam::init()
{
        // Bugcheck //
          assert(requiredMainTestingRow_!=0);

        // Interpolation Error Limit 
          errorLimit_ =  defaultInterpolationErrorLimit_ ;

        //+
        // Measurment Set Size Set-Up.
        //- 
          printf( "TuneMSParam::init()::  Pointing Interval = %f\n", pointingIntervalSec_);
          printf( "TuneMSParam::init()::  Main     Interval = %f\n", mainIntervalSec_);
          printf( "TuneMSParam::init()::  Main Testing Row  = %d\n", requiredMainTestingRow_ );

        //+
        //  Define Row Count in POINTING and Main
        //   - the reuired minimum numbers are up to Interval ration.
        //-
       
            Double TotalTime =  requiredMainTestingRow_ * mainIntervalSec_ ;
            availableNrowInPointing_ = requiredNrowInPointing_  
                                     = TotalTime / pointingIntervalSec_;

           // if POINTING table already have sufficient length 
            if (requiredNrowInPointing_ < defInerpolationTestPointingTableRow_)
            {
                printf(  "&&&&& Nrow in Pointing, forced to be set. \n" );
                requiredNrowInPointing_  = requiredMainTestingRow_ *  mainIntervalSec_ / pointingIntervalSec_ ;
            }

        //+
        // Optimize Row Count
        //   - when required row is insufficient, set adding count and expand MS later. 
        //   - when multiple-antenna is used number of tables are increased as below.(CAS-8418)
        // (bug fix)
        //    7-MAR-2019 prepareMaxAntenna_ was applied for exact row control.
        //-

        if ( (prepareMaxAntenna_ * requiredNrowInPointing_) > defInerpolationTestPointingTableRow_ )
        {
            extraNrowInPointing_ = requiredNrowInPointing_ * prepareMaxAntenna_ 
                                 - defInerpolationTestPointingTableRow_;
        }
        else
        {
            extraNrowInPointing_ = 0;
        }

        //+
        // Expanding MAIN
        //-
        if ( (prepareMaxAntenna_ * requiredMainTestingRow_) > defInerpolationTestMainTableRow_ )
        {
            extraNrowInMain = requiredMainTestingRow_ * prepareMaxAntenna_
                            - defInerpolationTestMainTableRow_;
        }
        else
        {
            extraNrowInMain  = 0 ;
        }

        printf( "DBG::TuneMSParam::init()::File Size: Pointing, required =%f, adding size = %f \n", 
                requiredNrowInPointing_ ,extraNrowInPointing_);
        printf( "DBG::TuneMSParam::init()::File Size: MAIN    , required =%u, adding size = %f \n", 
                requiredMainTestingRow_ ,extraNrowInMain);

}

void TuneMSParam::Initialize(Double p_interval, Double m_interval )
{
        pointingIntervalSec_         =  p_interval;
        mainIntervalSec_             =  m_interval;

        // Number of Row count used in TEST // 
        requiredMainTestingRow_ =  defultMainTestingRow;

        // Inerpolation Error limit (currently active)
        errorLimit_  = defaultInterpolationErrorLimit_;

        // common init //
        init();
}

void TuneMSParam::Initialize( )
{
        Double p_int = 0.048;
        Double m_int = 1.008;

        setMainRowCount(defultMainTestingRow);
        Initialize(p_int, m_int);  // Give Pointing and Main Intervals. //
}


//+
//  Generate Pseuo Direction / Time Infomation
//  both for Pointing and Main.
//-

TuneMSParam::PseudoPointingData  TuneMSParam::pseudoPointingBaseInfo(Double deltaTime)
{
        casacore::Vector<Double> point;
        point.resize(5);

        PseudoPointingData  point2; 

        //+
        //  relative time limit (r_time)
        //   (on private)
        //-

        if (r_time > 1.0) {
            printf( "r_time::Exceeded 1.0: %f \n",r_time );
            assert( r_time < 1.0); 
        }
 
        //+
        //  Determin TIME upon Base Date. 
        //    dd : in day.
        //-

        Double time  = deltaTime * Interval;
        Double dd    =  (22 *3600.0 
                         +  5*60 +  41.5 
                         + time  
                       ) / (3600*24) ;
     
        casacore::MVTime  basetime (2003,11,12 ,dd);     
 
        //+
        // Designed Function of POINITNG location
        //-

        Double X2[5] = {}; // all clear
        Double Y2[5] = {}; // all clear 
       
        //  Trajectory Function execution
        //    (memo) debug function should be build in this class.
        //-

        // prepare five sets // 

        TrajectoryFunction::getInstance().calc( r_time, X2[0], Y2[0] );
        TrajectoryFunction::getInstance().calc( r_time, X2[1], Y2[1] );
        TrajectoryFunction::getInstance().calc( r_time, X2[2], Y2[2] );
        TrajectoryFunction::getInstance().calc( r_time, X2[3], Y2[3] );
        TrajectoryFunction::getInstance().calc( r_time, X2[4], Y2[4] );

        // Intentional Offset (for identify 5 sets)

        if(false)
        {
            for(uInt n=0;n<5;n++)
            {
                X2[n] += 0.01 * (Double)n;
                Y2[n] += 0.01 * (Double)n;
            }
        }

        // Probe the range //

        for(uInt n=0;n<5;n++) {
	    assert( abs(X2[n]) <=  M_PI );
            assert( abs(Y2[n]) <=  M_PI/2.0 );
        }

        // Direction Values //
        // CAS-8418 New  Interface //

        point2.time         = basetime.second();
        point2.interval     = Interval;
        point2.relativeTime = r_time;  

        // Five different output values. //
        for(uInt n=0;n<5;n++) {
            point2.position[n] = make_pair(X2[n],Y2[n]);
        }        

        return point2; 
}


TuneMSParam::PseudoPointingData TuneMSParam::pseudoPointingInfoPointing2(Double deltaTime)
{
    // privide local conditon on private variables //
      Interval =   pointingIntervalSec_;
      nRow     =   availableNrowInPointing_; 
      r_time   =   deltaTime/availableNrowInPointing_;
             
      return(pseudoPointingBaseInfo(deltaTime));
}
             
TuneMSParam::PseudoPointingData  TuneMSParam::pseudoPointingInfoMain2(Double deltaTime)
{
    // privide local conditon on private variables //
 
    Interval = mainIntervalSec_;

    //+
    // Determine number of row of Pointing and Main table.
    // this depends on which total time is longer.
    //-              
    if(pointingIntervalSec_ <= mainIntervalSec_ ) // Ordinary case
    {
        nRow=  requiredMainTestingRow_;
        r_time =   deltaTime / nRow ;
    }   
    else   // This case may happen that not all the rows in Pointng are used.
    {
        nRow = max ( requiredMainTestingRow_ ,  defInerpolationTestMainTableRow_ );
        r_time =   deltaTime / nRow   ;
    } 

    return(pseudoPointingBaseInfo(deltaTime));

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
 
      /* IPO */

        IPosition getIpo() {return  pointingDirection.shape(0); }

      /* Scalor */
        uInt   getAntennaId(uInt row)       { return pointingAntennaId. get(row); }
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
          void putAntennaId(uInt row, Double dd) {  pointingAntennaId.  put(row, dd ); }

     /// Debug Utility ///

          void dump(String fname)
          {
              FILE *fp=fopen( fname.c_str(), "w");
              Double prevTime =0.0;
              Double intervalD =0.0;
              for(uInt row=0;row<nrow;row++)
              {
                  uInt   antId    = getAntennaId(row);
                  Double time     = getTime(row);
                  Double interval = getInterval(row); 
               
                  Vector<Double> Dir = getDirection(row);
                  Vector<Double> Tar = getTarget(row);
                  intervalD = time - prevTime;
                  prevTime  = time;
                  if(intervalD <  interval * 1.1)  intervalD = 0.0;
 
                  fprintf(fp, "%d,|,%d,%f,%f,|,%f,%f,|,%f,%f,,%f \n",
                              row, antId, time, interval,
                              Dir[0], Dir[1], Tar[0],Tar[1],intervalD  );
              }
              fclose(fp);
          }
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


//************************************************
// [NEW C++]
// Access Service class
//   for MAIN table
//************************************************

class MainTableAccess {

public:
    // constructor //
    MainTableAccess(String const &MsName, bool WriteAccess =false )
    {
        auto option = casacore::Table::TableOption::Old;
        if(WriteAccess) option = casacore::Table::TableOption::Update;

        MeasurementSet ms_t(MsName, option );
        ms = std::move(ms_t);

        init();
        prepareColumns();
    }

    // Destructor //
    ~MainTableAccess() { }
    void flush()
    {
        ms.flush();
        ms.resync();
    }
    // Row //
    uInt getNrow()               { return (nrow = ms.nrow()) ;  };
    void appendRow(uInt AddCnt)  { ms.addRow(AddCnt); }
    void removeRow(uInt row)      { ms.removeRow(row); }

    // Write (put) //

    void putAntenna  (uInt row, Double dd) { antenna1_col.     put(row, dd );}
    void putAntenna2 (uInt row, Double dd) { antenna2_col.     put(row, dd );}

    void putTime    (uInt row, Double dd) { mainTime_col.    put(row, dd );}
    void putInterval(uInt row, Double dd) { mainInterval_col.put(row, dd );}

private:

    void init()
    {
        nrow     = ms.nrow();
    }

    void prepareColumns()
    {
        antenna1_col       .attach( ms  , "ANTENNA1");
        antenna2_col       .attach( ms  , "ANTENNA2");
        mainTime_col      .attach( ms , "TIME");
        mainInterval_col  .attach( ms , "INTERVAL");
    }


    // handle (in MS, directly connects to Columns)
     casacore::MeasurementSet     ms;
     casacore::uInt               nrow;
    // Columns //
     ROScalarColumn<Int>    antenna1_col  ;
     ROScalarColumn<Int>    antenna2_col  ;
     ROScalarColumn<Double> mainTime_col ;
     ROScalarColumn<Double> mainInterval_col ;

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
 
    TuneMSParam  tuneMS;

    MsEdit() {     }   // Current 

    MsEdit(const String &ms, bool mode=false )
    {     
        // CAS-8418
        unique_ptr<PointingTableAccess>  ptTemp( new PointingTableAccess(ms,mode));
        unique_ptr<AntennaTableAccess>   anTemp( new AntennaTableAccess(ms,mode));
        ptTblAcc_ = std::move(ptTemp);
        anTblAcc_ = std::move(anTemp);

    }

    // Add or Remove Column (Pointing) ////
   
        void duplicateNewColumnsFromDirection();

    // Add Row //

        uInt appendRowOnAntennaTable(uInt n);
        uInt appendRowOnPointingTable(uInt n);
        uInt appendRowOnMainTable(uInt n);
    
    // Write Data on Antenna Table // 

        void writeDataToAntennaTable( uInt Row =0 );
        void prepareDataToAntennaTable();

    // Write (generated) Test Data on Pointing Table //
    // Write (generated) Test Data on MAIN Table //

        void writePseudoOnPointing ( );
        void writePseudoOnMainTable     (Double dt, String MsName =DefaultLocalMsName);

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
   // CAS-8418 
   //  Prividing easier access on this UT to construct functions
   //  in this class.  
   //-
    std::unique_ptr<casa::PointingTableAccess>     ptTblAcc_;
    std::unique_ptr<casa::AntennaTableAccess>      anTblAcc_;

};


//+
// CAS-8418 Add one row on Antanna Table
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

void setData(ANTENNADataBuff &data, uInt id)
{
    data.name            =  "ZZ0"+std::to_string(id);

    printf( "writeDataToAntennaTable():setData  name =[%s]\n",data.name.c_str() );

    data.station         =  "Station";
    data.type            =  "lovely-KUMIKO";
    data.mount           =  "Mount Fuji";
}


void MsEdit::prepareDataToAntennaTable( )
{
    String MsName =DefaultLocalMsName;
    AntennaTableAccess ata(MsName,true);

    uInt nrow_a = ata.getNrow();
    printf( "MsEdit: prepareDataToAntennaTable () current nrow=%d \n", nrow_a);

    // Buffer //
    ANTENNADataBuff dt0;

    // Position, Offset  3dim. //
    dt0.position.resize(3);
    dt0.offset.resize(3);
    for(uInt i=0; i<3; i++)
    {
        dt0.position[i] = (Double)i * 0.1 + 100050.0 ;
        dt0.offset[i]   = (Double)i * 0.1 + 50.0 ;
    }

    // Write records as defined. //
    uInt N = tuneMS.getMaxOfAntenna();
    printf( "DBG:: number of Antena  = %d ready \n",N );
    for(uInt ant = 0; ant < N; ant++)
    {
        // set  //
        setData(dt0, ant);
        // Put Data //
        ata.putRowData(ant, dt0);
    }

    // flush and close // 
    ata.flush();

}

//+
// Add rows by specified count on Pointing Table Table
//-

uInt  MsEdit::appendRowOnPointingTable(uInt AddCount )
{
    // CAS-8418 new code //
    {
      PointingTableAccess pta(MsName,true);
      pta.appendRow(AddCount);
      pta.flush();
    }
    
    PointingTableAccess pta(MsName,true);
    uInt nrow = pta.getNrow();
   
    return nrow;
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
    // CAS-8418 new code //
    { 
        String MsName = DefaultLocalMsName;
        PointingTableAccess ata(MsName,true);
        ata.duplicateColumns();
    } /* Once close */ 

    {
        PointingTableAccess ata(MsName,true);
        ata.fillNewColumns();
    }
} 

//*****************************************************************
// Wtite Test Data on Direction Column in Pointing Table
//  - Values are got by sub fuction above.
//  - SetUp() in class TestDirection calls this.
//  - see also TEST Fixture  
//****************************************************************

void  MsEdit::writePseudoOnPointing()
{
    // CAS-8418 CODE 
    String MsName =DefaultLocalMsName;
    PointingTableAccess pT( MsName, true);

    //+
    //  Loop for each Row,
    //  (NOTE) In particular case , LoopCnt requires ONE more.
    //          In ordinary case +1 causes Exceeding row count.
    //-

    uInt LoopCnt =   tuneMS.getAvailablePointingTestingRow();

    //+
    // For all Antenna and Row 
    //+
    uInt N = pT.getNrow();
    printf( "DBG:writePseudoOnPointing:: Creating Data on Pointing Table. (Currently nrow=%d)\n",N);
    for (uInt ant_id=0; ant_id < tuneMS.getMaxOfAntenna() ; ant_id++ )
    {
            printf( "DBG:writePseudoOnPointing::creating row data ( Ant=%d) data. Loop=%d \n", 
                     ant_id, LoopCnt );

            for (uInt row=0; row < LoopCnt; row++)
            {
                uInt  rowA = row + (ant_id * LoopCnt);

                //+
                // DIRECTION  
                //   CAS-8418::   1-Feb-2019
                //   updated to make indivisual value on Direction Columns
                //-
                Double deltaTime = (Double)row  ;    // deltaTime must be [0,1]
                IPosition Ipo = pT.getIpo();
                Array<Double> direction(Ipo, 0.0);   // IP shape and initial val // 

                Vector< Array<Double>  > Dir5;
                Dir5.resize(5);

                // Calculate Pseudo-Direction as test data //

                TuneMSParam::PseudoPointingData  psd_data2  
                        = tuneMS.pseudoPointingInfoPointing2(deltaTime); // generated pseudo data. (Pointing) //
 
                if(tuneMS.ifCoeffLocTest() )
                {
                    //+
                    // Test Pattern Data 
                    //-
                    for(uInt cno=0; cno<5; cno++)
                    {   
                        // const value ..// 
                          direction[0][0] = 0.1 + 0.1 * (Double)cno  ;
                          direction[0][1] = 0.1 + 0.1 * (Double)ant_id;
                          Dir5[cno] = direction;
                    }
                }
                else  // Ordinary cases... // 
                {
                    for(uInt cno=0; cno<5; cno++)
                    {
                       // Add offset 
                          direction[0][0] = psd_data2.position[cno].first  ;
                          direction[0][1] = psd_data2.position[cno].second ;
                          Dir5[cno] = direction;
                    } 

                }

                // write access //
                // (CAS-8418: Multiple Antenna supported) use rowA instead of row//
                pT.putDirection      (rowA, Dir5[0] );
                pT.putTarget         (rowA, Dir5[1] );
                pT.putPointingOffset (rowA, Dir5[2] );
                pT.putSourceOffset   (rowA, Dir5[3] );
                pT.putEncoder        (rowA, Dir5[4] );

                //+ 
                // New Time (shifted)  (this  activates interporation)
                //  basically FIXED values.
                //-

                pT.putTime      (rowA, psd_data2.time );     // Time
                pT.putInterval  (rowA, psd_data2.interval ); // Interval

                pT.putAntennaId (rowA, ant_id ) ;     // AntennaID 


            }//end row
    }// end antid


    pT.flush();
    printf( "DBG:writePseudoOnPointing:: end\n");
}

//+
// CAS-8418 Add rows by specified count 
// on Main Table Table
//-

uInt  MsEdit::appendRowOnMainTable(uInt AddCount )
{

// CAS-8418 NEW CODE //

     MainTableAccess   mta(MsName,true);
     mta.appendRow(AddCount);
     uInt nRow = mta.getNrow();
     mta.flush();

     return nRow;
}

//+
// Wtite Test Data on Direction Column in MAIN TABLE
//  -  sample shift  [ 0 < deleta < 1 ]
//  -  time diffrence = deltaTime * interval 
//-

void  MsEdit::writePseudoOnMainTable(Double deltaTime, String MsName)
{
    printf( "   deltaTime =, %f \n", deltaTime );

//******************
// CAS-8418 CODE 
//******************
    MainTableAccess   mta(MsName,true);

    uInt nrow_ms = mta.getNrow();
    uInt LoopCnt = tuneMS.getRequiredMainTestingRow();

    printf("writing to MAIN, nrow=%d, number of data on each antenna=%d \n", nrow_ms,LoopCnt );
    for (uInt ant_id=0; ant_id < tuneMS.getMaxOfAntenna() ; ant_id++ )
    {
        printf("DBG:: Ant=%d \n",ant_id);
        for (uInt row=0; row < LoopCnt; row++)
        {
            uInt  rowA = row + (ant_id * LoopCnt);

            // Pseudo Data (TEST DATA );

            TuneMSParam::PseudoPointingData  psd_data2
                   = tuneMS.pseudoPointingInfoMain2( (Double)row); // generated pseudo data. (Main table) //


            // Time Set  //

            Double interval = psd_data2.interval;
            Double time     = psd_data2.time;

            Double SetTime = time + (deltaTime * interval) ; 

            mta.  putAntenna  (rowA, ant_id   );      // AntennaID1 (CAS-8418)
            mta.  putAntenna2 (rowA, 0   );           // AntennaID2  always fixed = 0 

            mta.  putTime    (rowA, SetTime  );      // Time     (( REvised 10.26))
            mta.  putInterval(rowA, interval );      // Interval (( Revised 10.26))

        }
    }

     mta.flush(); 
     printf( "DBG:return writePseudoOnMainTable  \n");

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
    // CONSTRUCTOR  //

        MeasurementSet ms0( name  );      
        PointingDirectionCalculator calc(ms0);

    // Initial brief Inspection //
   
        printf("# Constuctor Initial Check.  [%s] row =%d\n", name.c_str(),calc.getNrowForSelectedMS()  );
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
                               std::to_string(n)+". " + MsList.name(n).c_str()  );
        String name = env.getCasaMasterPath()+ MsList.name(n);
        if ( MsList.isException(n))
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

        //+ 
        // TestCondition const
        // (Programmer Tunable)
        //-

            // Resources (Tunable) //
            uInt const InterpolationDivCount_ = 10;
            uInt const numAntenna_            = 3;
            uInt const numPointingColumn_     = 5;

            // Listing option
            bool       fgResultListing = false;
            // Convertion option (by setFrame)
            bool       fgConversion    = false;
            // Google Test On/OFF
            bool       fgGoogleTest    = true;

        // Prepared Antenna in MS //

        void setMaxAntenna(uInt n) { msedit.tuneMS.setMaxAntenna(n);    }

        // Pointing Colum List (common definition) //

        PointingColumnList pColLis_;

        // NEW: CAS-8418 relatedly commonized

        unique_ptr<casa::PointingDirectionCalculator> calc0;
        casa::PointingDirectionCalculator             *pdc;

        void start(const String msName)
        {
            MeasurementSet ms( msName.c_str() );

            // create PDCalc obj. on this class to share //
            unique_ptr<casa::PointingDirectionCalculator> 
                   calcTmp( new casa::PointingDirectionCalculator(ms));
            calc0 = std::move(calcTmp);
            pdc   = calc0.get();
        }

        void setCondition(uInt numRow, Double pointingInterval, Double mainInterval, Double errLimit)
        {
            msedit.tuneMS.    setMainRowCount   (numRow);       // aprox. 1-2H 
            msedit.tuneMS.    Initialize( pointingInterval,     // Pointing Interval
                                         mainInterval ) ;      // Main Interval

            msedit.tuneMS.    setInterpolationErrorLimit( errLimit );
        }
        void selectTrajectory( uInt no )
        {
            TrajectoryFunction::getInstance(). setType(no);   // Trajectory Function Type // 
        }

        void prepareAntenna()
        {
            uInt N =  msedit.tuneMS.getMaxOfAntenna();

            msedit.appendRowOnAntennaTable(N-1);   // add  more (#0 is ready)
            msedit.prepareDataToAntennaTable();
        }

        void prepareRows()
        {
            msedit.appendRowOnPointingTable ( msedit.tuneMS. getAddInerpolationTestPointingTableRow() );
            msedit.appendRowOnMainTable     ( msedit.tuneMS. getAddInerpolationTestMainTableRow() );
        }

        Double getPointingInterval() { return(msedit.tuneMS.getPointingTableInterval()); }
        Double getMainInterval()     { return(msedit.tuneMS.getMainTableInterval()); }
 

        Double getErrorLimit() { return (msedit.tuneMS.getInterpolationErrorLimit());}

        // Current Pointing Column //
        void setCurrentPointingColumn(uInt no ) {currentPointingColumn = no;       }
        uInt getCurrentPointingColumn()         { return currentPointingColumn ;     }

        // Current Antenna Select //
        void setSelectedAntenna(uInt no) {curretSelectedAntenna  = no;       }
        uInt getSelectedAntenna()        { return curretSelectedAntenna;     }

        // Easy Access // 
        void writeOnPointing()
        {
            msedit.writePseudoOnPointing () ;
        }
 
        void writeOnMain(Double dt, String MsName =DefaultLocalMsName)
        {
            msedit.writePseudoOnMainTable (dt, MsName);
        }

protected:

        // MeasurementSet Editting  //
     
          casa::MsEdit  msedit;

        // Handle Spline Object , and Coeff Table//
        SplineInterpolation  *sp ;           // =  calc.getCurrentSplineObj();
        SplineInterpolation::COEFF coeff ;   // =  sp->getCoeff();

        // Add 3 OFFSET Colums ,copied from DIRECTION column. //

        void addColumnsOnPointing() { msedit.duplicateNewColumnsFromDirection(); }

        //*
        // Sub-function of TEST_F(TestDirection....)  
        //*

        std::vector<Double>  testDirectionForDeltaTime2(Double dt, uInt p, uInt a);// Extended(CAS-8418)
        vector<Double>       testDirectionByInterval(Double p, Double m, uInt pc, uInt a);

        TestDirection(){ }
        ~TestDirection() { }

        virtual void SetUp()
        {
            BaseClass::SetUp();

            //+
            // Copy and add columns,  init columns, 
            //-

            CopyDefaultMStoWork();
            addColumnsOnPointing();

        }

        virtual void TearDown()
        {
            BaseClass::TearDown();

            // Delete Working MS 

             DeleteWorkingMS();
        }
private:
 
           uInt currentPointingColumn = 0; 
           uInt curretSelectedAntenna = 0;
};

//+
//  check Accessor ID (So far, this is no Assetion)
//-  

static void inspectAccessor( PointingDirectionCalculator  &calc )
{
    uInt Id = calc.getCurretAccessorId();
    SplineInterpolation *sp =  calc.getCurrentSplineObj();

    // Inspect Coefficient Table //
    ASSERT_EQ( sp->isCofficientReady(), true );

    printf( "Spline Cofficient in Pointing Column[%d] GTest OK.\n",Id);
}


//------------------------------
// Revised Edition (CAS-8418)
//  (replaced from old to new)
//------------------------------
TEST_F(TestDirection, setDirectionColumn  )
{
    start(DefaultLocalMsName);
    expectedNrow = pdc->getNrowForSelectedMS();
    EXPECT_NE((uInt)0, expectedNrow );

    for( uInt n=0; n<2;n++ ) 	// 2 Times. run .../
    {
        for(size_t k=0; k < pColLis_.size(); k++)
        {     
            String ColName = pColLis_.name(k);
            Description("Column Name" , ColName );

            // UNIT TEST //
            EXPECT_NO_THROW( pdc->setDirectionColumn( ColName ) );
            inspectAccessor(*pdc);
        }
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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

    // setDirectionListMatrixShape           //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    // setFrame()  //

        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //+
    // setDirectionColumm()  calls 
    //-

        //+
        // Normal Seq. with setMovingSoure Convert.
        //  using String , deffault values are set.
        //-

	FunctionalDescription("Normal Seq.", "Selectve Convert");


        for(size_t k=0; k < pColLis_.size(); k++)
        {
            Description("Column Name" , pColLis_.name(k) );
            EXPECT_NO_THROW( calc.setDirectionColumn( pColLis_.name(k) ) );

            if(true)  // Selective:: TESTING dependency how  setMovingSource() relates. //
            {
                String src = "SUN";
                Description("setMovingSource()", src);
                EXPECT_NO_THROW( calc.setMovingSource( src ) );
            }    

            Description("calling  getDirection() ", pColLis_.name(k) );
            Matrix<Double>  DirList;
            EXPECT_NO_THROW( DirList= calc.getDirection() );

            uInt  n_row    = DirList.nrow();

            printf( "Number of Row = %d \n", n_row );
            EXPECT_EQ( n_row, expectedNrow);

        }

        //+
        // No SetMovingSouce executution Sequence. 
        //-
        
        FunctionalDescription("Normal Seq.", "Always call setDirectionColumn");

        for(uInt k=0; k < pColLis_.size(); k++)
        {
            Description("Column Name" , pColLis_.name(k) );
            EXPECT_NO_THROW( calc.setDirectionColumn( pColLis_.name(k) ) );

            String src = "SUN";
            Description("setMovingSource()", src);
            EXPECT_NO_THROW( calc.setMovingSource( src ) );
            
            Description("calling  getDirection() ", pColLis_.name(k) );
            Matrix<Double>  DirList;   
            EXPECT_NO_THROW( DirList = calc.getDirection() );

            uInt  n_row = DirList.nrow();

            printf( "Number of Row = %d \n", n_row);
            EXPECT_EQ( n_row, expectedNrow);

        }

}

/*-------------------------------------------
  Verification Test of CAS-11818
  - If you use old source, this test causes
    core dump.
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
        expectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, expectedNrow );


    // setFrame, setDirectionListMatrixshape    //
    // setMovingSource Convert                  //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //
    // setDirectionColumm()  calls 
    //
         
        // No SetMovingSouce executution 

        FunctionalDescription("BUG-TEST by assert()" , "setDirectionColumn+ getDirection, without setMovingSource");


        for(uInt k=0; k<5; k++)
        {
            Description("Column Name" , pColLis_.name(k) );
            EXPECT_NO_THROW( calc.setDirectionColumn( pColLis_.name(k) ) );

            Description("calling  getDirection() ", pColLis_.name(k) );
            Matrix<Double>  DirList;   
#if 1
            EXPECT_NO_THROW( DirList= calc.getDirection() );
#else
            DirList= calc.getDirection();
#endif 
            uInt  N_Row;
            N_Row  = DirList.nrow();

            printf( "Number of Row = %d \n", N_Row );
            EXPECT_EQ( N_Row, expectedNrow);

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
        expectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, expectedNrow );

    // setDirectionListMatrixShape()    //

        Description("setDirectionListMatrixShape()", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    // setFrame() //
    
        Description("setFrame()", "J2000" );
        EXPECT_NO_THROW( calc.setFrame( "J2000" ) );

    //+
    //  Four senarios are executed.
    //-  	

    for(uInt senario = 0; senario <= 3; senario++ )	// senario [0,1,2,3] //
    {   
          FunctionalDescription("Senario" , std::to_string(senario).c_str() ); 

          for(size_t k=0; k < pColLis_.size(); k++)
          {
              Description("- Column Name" , pColLis_.name(k) );
              EXPECT_NO_THROW( calc.setDirectionColumn( pColLis_.name(k) ) );

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

              // getDirection //
              Matrix<Double>  DirList;
              Description("- getDirection() ", pColLis_.name(k) );
              EXPECT_NO_THROW( DirList  = calc.getDirection() );

              uInt n_row  = DirList.nrow();

              printf( "- Number of Row = %d \n", n_row );
              EXPECT_EQ( n_row, expectedNrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

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
        EXPECT_EQ( n_row, expectedNrow);

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

//-------------------------------------------------
//    Interporation Test (sub) 
//     dt : displacement betweeen measure point
//     [0 <= dt <= 1]   dt=0; X[n],  dt=1, X[n+1]
//------------------------------------------------
//+
// BEING REVISED (CAS-8418)   22-FEB-2019 
//    Extend multiple Direction Colums, multiple Antenna.
//-
std::vector<Double>  TestDirection::testDirectionForDeltaTime2(Double dt, uInt colNo, uInt ant )
{
    const String MsName = DefaultLocalMsName;

    // Create Object //
        MeasurementSet ms( MsName.c_str() );
        PointingDirectionCalculator calc(ms);   

    //  MatrixShape (COLUMN_MAJOR)
        calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) ; 

    // setFrame()
        String FrameName= "J2000";
        calc.setFrame( FrameName );

    //+
    // selectData on generated MS
    //  Antenna name is ANT0, ANT1, AT2
    //-
   
        const String AntSel = "ZZ0" + std::to_string(ant)+ "&&ZZ00"; 

        calc.selectData( AntSel,  "","","","","","","","","" );
        /* uInt match_row = calc.getNrowForSelectedMS();    */
     
    //  InterPolation mode (Foece Unuse)
        calc.setSplineInterpolation(use_spline);

    //+
    // setDirectionColumn()
    //
    //   NOTES: If multiple loop is intened.
    //       setDirectionColumn() must be working with psd data.  
    //-
#if 0
        uInt colNo = getCurrentPointingColumn() ; 
#endif 
        printf( "setDirectionColumn(%s)\n", pColLis_.name( colNo ).c_str() );

        calc.setDirectionColumn( pColLis_.name(colNo) ) ;

    //+
    // setFrame call for "conversion"
    // in CAS-8418, try to examine converted result.
    //-
        if(fgConversion)
        {
            String frameName = "AZELGEO" ;
            calc.setFrame( frameName );
        }

    //+
    //  getDirection()
    //-
        Matrix<Double>  DirList1  = calc.getDirection();
        size_t   n_row    = DirList1.nrow();

    //+
    // Direction Interpolation
    //-
    Double maxErr_1 = 0.0;
    Double maxErr_2 = 0.0;             
    Double absErr_1 = 0.0;
    Double absErr_2 = 0.0;


    uInt LoopCnt = n_row-2;
    for (uInt row=0; row < LoopCnt; row++)   
    {
        // Direction(1) by getDirection //

          Double calculated_1 = DirList1(row,0);
          Double calculated_2 = DirList1(row,1);

        // Direction by generated/estimated //

          TuneMSParam::PseudoPointingData  gen_out2
                  = msedit.tuneMS.pseudoPointingInfoMain2 ( (Double)row  + dt );  // dt:Interpolation offset  (sec)

          Double generated_1 = gen_out2.position[colNo].first;
          Double generated_2 = gen_out2.position[colNo].second;

        //+
        // Error calculation 
        //-

          Double Err_1 = calculated_1 - generated_1 ;   
          Double Err_2 = calculated_2 - generated_2 ;

          absErr_1 = abs(Err_1);
          absErr_2 = abs(Err_2);
                
          if( maxErr_1 < absErr_1 ) maxErr_1 = absErr_1;
          if( maxErr_2 < absErr_2 ) maxErr_2 = absErr_2;

          // Google Test On/Off
          if(fgGoogleTest)
          {
	      EXPECT_LE( absErr_1, msedit.tuneMS.getInterpolationErrorLimit()  ); 
              EXPECT_LE( absErr_2, msedit.tuneMS.getInterpolationErrorLimit()  ); 
          }

          // Output List (in One line) On/OFF   //
          if(fgResultListing) 
          {
              printf( "Evaluation,");
              printf( "%6d, %13.10f,%13.10f,", row,  calculated_1, calculated_2 );
              printf( "%13.10f,%13.10f,",      generated_1,  generated_2 );
              printf( "%12.4e,%12.4e \n",      Err_1,     Err_2);
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
std::vector<Double> TestDirection::testDirectionByInterval(Double p_int, Double m_int,
                                                           uInt p_col, uInt antenna )
{
    // Max Error ///

      std::vector<Double> reterr;
      ErrorMax            maxerr;
    
    // Add INTERPOLATION TEST DATA 

      writeOnPointing();

    // Test Loop  
    //   - nDiv is given, internd to spacify number of separating count 
    //     between one exact point and the next,  

    uInt nDiv = InterpolationDivCount_; 
    for (uInt loop=0; loop < nDiv; loop ++ )
    {
        //+
        // SetUp Testing  MeasurmentSet
        //-  

        writeOnMain((Double)loop/(Double)nDiv);

        // Execution //
        reterr = testDirectionForDeltaTime2( (Double)loop/(Double)nDiv, p_col, antenna );

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
   Interporatio Test in getDirection()  
   FULL combiniation of Interval Times and spline mode.

 - Use Combiniation of Pointing table Interval and Main table Interval
 - Capable of selecting curve fucntion for simulated pointing trajectry
  ----------------------------------------------------------------------*/
 
TEST_F(TestDirection, InterpolationFull )
{
  TestDescription( "Interpolation Full-combiniation 1) mode,2) Interval" );
    // Combiniation List of Pointing Interval and Main Interval //

    vector<bool>   InterpolationMode     = { false, true};
    vector<Double> Pointing_IntervalList = { 0.1, 0.05, 0.01  };
    vector<Double> Main_IntervalList     = { 0.2, 0.1, 0.01,  };

    ErrorMax  maxerr;
    std::vector<Double> r_err = {0.0}; 

    //+
    // Use following paramters
    //-
      uInt usingColumn  = 0;
      uInt usingAntenna = 0;

    // Multiple Loop //
    for(uInt s=0; s<InterpolationMode.size(); s++)
    {  
        // Spline OFF, ON //
        use_spline = InterpolationMode[s];

       // Combiniation Loop
       //  5-MAR-2919 ; changed out/in loop var. 
        for( uint p=0; p < Pointing_IntervalList.size(); p++)
        {
            for( uint m=0; m < Main_IntervalList.size(); m++)
            {
                // Interval // 
                 Double p_i = Pointing_IntervalList[p];
                 Double m_i = Main_IntervalList[m];

                // Error Limit
                 Double ErrLimit =0.0;
               
                 if( use_spline==true) { ErrLimit = 1E-06; }
                 else                  {  ErrLimit = 1E-04;}
       
                 if(p_i > m_i )        {  ErrLimit = 1E-02;}

                 printf( "======================================================\n");
                 printf( " Mode[%d] Interval (P=%f,M=%f)  Limit =%f             \n", 
                          use_spline , p_i, m_i , ErrLimit );
                 printf( "======================================================\n");

                 // Copy Template MS //
                   SetUp();

                 //+
                 // set Examination Condition (revised by CAS-8418) //
                 //-
    
                 selectTrajectory( TrajectoryFunction::Type::Simple_linear ); 
                 setCondition( 5040,   /*numinTestingRow */      //number of row
                               p_i,    // Pointing Interval
                               m_i,    // Main Interval
                               ErrLimit );  // Error limit 
 
                 // Prepate Antenna (for Multple-set) //
                   prepareAntenna();
 
                 // Increase(Append)  Row on MS for large-file.:
                   prepareRows();
          
                 //+
                 // Execute Main-Body , get error info //
                 //-
                   r_err = TestDirection::testDirectionByInterval( p_i, m_i, usingColumn, usingAntenna );
                   maxerr.put(r_err);
             }
        } // End Interval combiniation

        std::vector<Double>  e = maxerr.get();
        printf ( "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n");
        printf ( "   THE WORST Error = %e, %e \n", e[0], e[1] );
        printf ( "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n");

    } // End Interpolation Mode 
}

/*-----------------------------------------------------------------------
  Interporatio Test in getDirection()   
  Specific COMBIIATION Parameter mode
 - Set of testing parameters are given
------------------------------------------------------------------------*/ 


typedef struct Parm {
    bool   use_spline;
    Double p_interval;
    Double m_interval;
    Double errLimit;
} ParamList;

std::vector<ParamList>
 pList =
{
    {true, 0.01,  0.05,  5.0E-05 },
    {true, 0.05,  0.05,  5.0E-05 },
    {true, 0.10,  0.05,  5.0E-05 },
    {true, 0.15,  0.05,  5.0E-04 },
    {true, 0.20,  0.05,  5.0E-04 },
    {true, 0.25,  0.05,  5.0E-04 },
    {true, 0.30,  0.05,  5.0E-04 },
    {true, 0.35,  0.05,  5.0E-04 },

    {true, 0.1,    0.1,   5.0E-08},
    {true, 0.048,  0.001, 4.0E-05},
    {true, 0.050,  2.5,   4.0E-05}
};

TEST_F(TestDirection, InterpolationCombiniation )
{
  TestDescription( "Interpolation by Listed condition." );
    // Combiniation List of Pointing Interval and Main Interval //

    ErrorMax  maxerr;
    std::vector<Double> r_err = {0.0}; 

    //+
    // Commonly use following paramters
    //-
      uInt usingColumn  = 0;
      uInt usingAntenna = 0;

    // Loop //
    for(uInt n=0; n<pList.size();n++)
    {
          printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n"   ); 
          printf("DBG:: parameter Set[%d]  starts. \n",n );
          printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n"   );

          use_spline       = pList[n].use_spline;
          Double p_i       = pList[n].p_interval;
          Double m_i       = pList[n].m_interval;
          Double err_limit = pList[n].errLimit;

          printf( "DBG:: Interval (Poinitng, Main) = (%f,%f) \n", p_i, m_i );

          // Copy Template MS //
          SetUp();

          //+
          // set Examination Condition (revised by CAS-8418) //
          //-
            selectTrajectory( TrajectoryFunction::Type::Simple_linear ); 
            setCondition( 5040,   /*numinTestingRow */      //number of row
                          p_i,    // Pointing Interval
                          m_i,    // Main Interval
                          err_limit );  // Error limit 
 
          // Prepate Antenna (for Multple-set) //
            prepareAntenna();
 
          // Increase(Append)  Row on MS for large-file.:
            prepareRows();
          
          //+
          // Execute Main-Body , get error info //
          //-
            r_err = TestDirection::testDirectionByInterval( p_i, m_i, usingColumn, usingAntenna );
            maxerr.put(r_err);
    }
}

/*-----------------------------------------------------------------------
  Interporation Test in getDirection()   
  SINGLE Parameter mode

- Set up one set of Pointing table Interval and Main table Interval
- Table sizes to be created is automatically tuned.
- Capable of selecting curve fucntion for simulated pointing trajectry
 ----------------------------------------------------------------------*/

TEST_F(TestDirection, InterpolationSingle )
{
    TestDescription( "Interpolation test in getDirection() SINGLE-parameter mode" );

      use_spline = true;

    // define Number of Antenna prepeared in MS //

      setMaxAntenna(5);

    // set Examination Condition (revised by CAS-8418) //

      selectTrajectory( TrajectoryFunction::Type::Simple_linear );// Trajectory(Curve) Function
      setCondition( 5040,       // number of row
                    0.1,          // Pointing Interval
                    0.05,         // Main Interval
                    1E-07  );  // Error limit 


    // Prepate Antenna (for Multple-set) //
      prepareAntenna();

    // Increase(Append)  Row on MS for large-file.:
      prepareRows();

    //+
    // For all Pointing Columns  and,
    //   For all regisetered Antenna 
    //-

    for(uInt pcol=0; pcol <numPointingColumn_; pcol++)  // THIS IS NOT A SECURE CODE // 
    {
        for(uInt ant=0;ant< numAntenna_ ; ant++)
        {   
            String info = "Col="+std::to_string(pcol)+", Ant="+to_string(ant);
            Description( "Single mode", info );

            // Indicate antenna select. //
              setSelectedAntenna( ant );

            // get currect interval time.. //
              Double p_interval = getPointingInterval();
              Double m_interval = getMainInterval();

            // Execute and get numerical error info 
             std::vector<Double> r_err = testDirectionByInterval(p_interval, m_interval,
                                                                 pcol,  ant );
             printf( " Total Max Error = %e, %e \n", r_err[0], r_err[1] );
        }
    }
}

//*****************************************
// Interporation Test in getDirection()
// Examine Coeff Table
//    with Pointing Column and AntennaId 
//*****************************************

TEST_F(TestDirection, CoefficientOnColumnAndAntenna )
{
      use_spline = true;

    // set Examination Condition  //

      selectTrajectory(TrajectoryFunction::Type::Zero); // Trajectory(Curve) Function
      setCondition( 5000,       // number of row
                    1.0,        // Pointing Interval
                    1.0,        // Main Interval
                    5.0E-06 );  // Error limit 

    // Prepate Antenna (for Multple-set) //
      prepareAntenna();

    // Increase(Append)  Row on MS for large-file.:
      prepareRows();

    // NEW: 8-MAR-2019 Set Pseudo control in Special mode 
      msedit.tuneMS.setCoeffLocTest( true ) ; 

    //+
    // Create MS for all Pointing Columns  and,
    //   for all regisetered Antenna 
    //-

    for(uInt pcol=0; pcol <numPointingColumn_; pcol++)  // THIS IS NOT A SECURE CODE // 
    {
        printf("DBG:: Direction Column = %d \n", pcol );

        // Internally indicate Pointing Column to be used
        setCurrentPointingColumn( pcol );

        // Write Data on Pointing TAble  
        writeOnPointing();
    }

    //+
    // setDirectionCplumn() call to create Spline Object.
    //   creating 5 spline objects with specified Pointing Column.
    //-

    const String MsName = DefaultLocalMsName;

    // Create Object //
    MeasurementSet ms( MsName.c_str() );
    PointingDirectionCalculator calc(ms);
   
    //+
    // Poinying-Column and Antenna Loop
    //- 
    PointingColumnList pList;
    for(uInt pcol=0; pcol <numPointingColumn_; pcol++) 
    {
        String name =pColLis_.name(pcol);
        printf( "DBG::calling setDirectionColumn(Pointing Column = %s) \n\n",name.c_str() );

        calc.setDirectionColumn(name);

        // Check out Spline Object ..//
          sp =  calc.getCurrentSplineObj();
          coeff = sp->getCoeff();

        //*
        // TENTATIVE (working) //
        // Check if the coeff is expected value.
        // at the moment only showing values.
        //*
        for(uInt ant=0;ant<numAntenna_;ant++)
        {
            for(uInt i =0; i<10; i++)
            {
                uInt Index = i;
                Double a0 = coeff[ant][Index][0][0];
                Double a1 = coeff[ant][Index][0][1];
                Double a2 = coeff[ant][Index][0][2];
                Double a3 = coeff[ant][Index][0][3];

                Double b0 = coeff[ant][Index][1][0];
                Double b1 = coeff[ant][Index][1][1];
                Double b2 = coeff[ant][Index][1][2];
                Double b3 = coeff[ant][Index][1][3];

                printf( "COEFF[%d][%d][XY][0-3] = %f,%f,%f,%f,|,",ant, Index, a0,a1,a2,a3  );
                printf( " %f,%f,%f,%f \n",b0,b1,b2,b3  );
            }
            printf( "\n");
        }
    }


}


TEST_F(TestDirection, CompareInterpolation )
{
  TestDescription( "Interpolation Result Comparing.(Linear vs Spline)" );

    std::vector<Double>   r_err1 ;
    std::vector<Double>   r_err2 ;

    Double            errorLimit = 5e-04;
   
    TestDescription( "getDirection (J2000) with selected data. uvw available" );
 
    // Use the below MS to evaluate the difference of
    //  interpolation results.

    std::vector<String> MsList = {
      "sdimaging/Uranus1.cal.Ant0.spw34.ms",
      "sdimaging/Uranus2.cal.Ant0.spw34.ms"
    };

    // TEST LOOP //

    for(uInt fno=0; fno<MsList.size(); fno++)
    {
        // Direction 
        Matrix<Double>  DirList1; // for Linear
        Matrix<Double>  DirList2; // for Spline

        // selected MS name //
        String name = env.getCasaMasterPath()+MsList[fno]; 
        printf( "MS[%s] is used. \n",name.c_str() );

        if(true){
           PointingTableAccess pta(name);
           pta.dump("Pointing_"+std::to_string(fno)+ ".csv" );
        }
        
        // Create Object //
        MeasurementSet ms0( name );
        PointingDirectionCalculator calc(ms0);

        //+
        // setDirectionColumn() 
        //-

        String ColName = "DIRECTION";
        calc.setDirectionColumn( ColName ) ;

        //+
        //  MatrixShape (COLUMN_MAJOR) 
        //-

        calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) ;

        //+
        // setFrame()
        //-

        if(true)
        {
            String FrameName= "AZELGEO";
            calc.setFrame( FrameName );
        }

        //******************
        // Examine Linear
        //******************
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
        //******************
        // Examine Spline
        //******************
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

        //****************************
        // List, Compare and Examine
        //****************************
        uInt size = DirList1.nrow();
        for( uint row=0; row<size; row++)
        {
            Double er1 = DirList2(row,0) - DirList1(row,0);
            Double er2 = DirList2(row,1) - DirList1(row,1);

            EXPECT_LT( er1, errorLimit );
            EXPECT_LT( er2, errorLimit );
       
            printf("row,%d, ", row );
            printf("Dir1, %e,%e,", DirList1(row,0),DirList1(row,1) ); 
            printf("Dir2, %e,%e,", DirList2(row,0),DirList2(row,1) );
            printf("Err, %e,%e \n", er1, er2 );
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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );


    //+
    // * Option *
    // selectData()   
    //  Please change if() to true,
    //-
        if(true)
        {
            const String AntSel = "DV01&&DV02";
            calc.selectData( AntSel,  "","","","","","","","","" );

            expectedNrow = calc.getNrowForSelectedMS();
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
         EXPECT_EQ( n_row, expectedNrow);

    //+
    // uv value 
    //  (first attach to Column e )
    //-
        casacore::ROArrayColumn<casacore::Double> uvwColumn;	
        uvwColumn .attach( ms0 , "UVW");

    //+
    // Dump Matrix
    //-
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

        // Select Keywords //
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
       expectedNrow = calc.getNrowForSelectedMS();
       printf("Expected nrow =%d\n",expectedNrow);

       EXPECT_NE( (uInt)0, expectedNrow );
       
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
        EXPECT_EQ (expectedNrow, nrow);     // see MS in detail //

      AntSel = "hoge&&&";
        FunctionalDescription("Testing Abnormal Name = hoge (No Matches))",AntSel );

        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);
        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      SpwSel = "";
        FunctionalDescription( "Spw:Nothig specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

    
      SpwSel = "*";
        FunctionalDescription( "Spw: Wildcard *  specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      SpwSel = "hoge";
        FunctionalDescription( "Spw: abnormal letters  specified.",SpwSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

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
        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE( (uInt)0, expectedNrow );

      uInt nrow;
      FieldSel = "";
        FunctionalDescription( "Files:Nothig specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      FieldSel = "*";
        FunctionalDescription( "Field: Wildcard *  specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      FieldSel = "hoge";
        FunctionalDescription( "Fieled: abnormal letters  specified.",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      FieldSel = "0";  // Field ID 
        FunctionalDescription( "Field: ID=0 specified.(No exits)",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ((uInt)0, nrow);       // On-Table , None on MAIN .. //

      FieldSel = "1";  // Field ID 
        FunctionalDescription( "Field: ID=1 specified.(exits)",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      FieldSel = "9";  // Out of Range
        FunctionalDescription( "Field: ID=9 Not on the FIELD TABLE",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      TimeSel = "";
        FunctionalDescription( "Time: Nothig specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      TimeSel = "hoge";
        FunctionalDescription( "Time: abnormal letters  specified.",TimeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      TimeSel = ">1900/01/01/00:00:00";
        FunctionalDescription( "Time: since 1900/01/01/00:00:00  specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      FeedSel = "";
        FunctionalDescription( "Feed: Nothig specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      FeedSel = "0";    // NO DATA /
        FunctionalDescription( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

      FeedSel = "1";    // NO DATA /
        FunctionalDescription( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      IntentSel = "";
        FunctionalDescription( "Intent: Nothig specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);
 
       // Usage : In what way this term is used. 
        //   and where to exist in MS.
          
      IntentSel = "*HOGE*";     // 
        FunctionalDescription( "Scan: ID specified.",IntentSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      ObservationSel = "";
        FunctionalDescription( "Observation: Nothig specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

     ObservationSel = "hoge";     
        FunctionalDescription( "Observation: abnormal expr..",ObservationSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      UVRangeSel = "";
        FunctionalDescription( "UVrange: Nothig specified.",UVRangeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      UVRangeSel = "hoge";     
        FunctionalDescription( "UVrange: abnormal expr..",UVRangeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

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
       expectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, expectedNrow );

      uInt nrow;
      MSSelect = "";
        FunctionalDescription( "MSselect: Nothig specified.",MSSelect);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (expectedNrow, nrow);

      MSSelect = "hoge";     
        FunctionalDescription( "MSselect: abnormal expr..",MSSelect);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (expectedNrow, nrow);

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
        expectedNrow = calc.getNrowForSelectedMS();
        EXPECT_NE((uInt)0, expectedNrow );

    // Various Frame Type (String) //

        for(uInt k = 0; k < DefinedFrametypes.size() ; k++  )
        {
            FunctionalDescription( "setFrame(Rsved Name) ", DefinedFrametypes[k].name.c_str() );

            // Execute and check Exception and other requirements//
        
            check_direction_info( calc, k ) ;
        }
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

    casa::MSNameList   nameList;


    ::testing::InitGoogleTest(& nArgs, args);

   // Run Test //

    return (RUN_ALL_TESTS()) ;

}

