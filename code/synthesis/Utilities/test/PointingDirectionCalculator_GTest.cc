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
 
// MS name and expected Exception //

struct _MSDef
{
    bool   ExThrow;  // True = cause Exeption
    String name;     // MS name, with relative path from CASAPATH
};

// MS name-list //

std::vector<struct _MSDef> TestMSList 
{
    // Exeption(bool)  , Filename //
 
        {true, "./sdimaging-t.ms"   		},
        {false, "sdimaging/sdimaging.ms"        },
        {false, "listobs/uid___X02_X3d737_X1_01_small.ms" },

    // Following 2 MS are affected assert(), cannot run on UT
    // Release EXE must throws Excepton. 
    //      {true,  "concat/input/A2256LC2_4.5s-1.ms"               },
    //      {true,  "concat/input/A2256LC2_4.5s-2.ms"               },

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

// Get File Name by Number //
const String  getMsNameFromList(uInt No )
{
    const String msg = "Internal Bugcheck:";
    if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());

    return TestMSList[No].name;
}

// Get Exception information //
bool  getMSThrowFromList(uInt  No )
{
    const String msg = "Internal Bugcheck:";
    if( !(TestMSList.size() >  No) ) throw AipsError(msg.c_str());
    
    return TestMSList[No].ExThrow;
}

// Get Count of MS List //
size_t getMSCountFromList()
{
    return TestMSList.size();
}

//******************************************************
//  Log Title Output Functions for readable text 
//  of this UT.
//******************************************************

void TestDescription( const String &Title )
{
    printf( "///////////////////////////////////////////////////////////////////////////// \n");
    printf( " %s  \n",Title.c_str() );
    printf( "///////////////////////////////////////////////////////////////////////////// \n");
    sleep(2);
}

void FunctionalDescription(const String &Title, const String &Param)
{
    printf("=============================================\n");
    printf("# %s \n",Title.c_str());
    printf("# [%s] \n",Param.c_str());
    printf("=============================================\n");
}

void Description(const String &Title, const String &Param)
{
    printf("+-----------------------------\n");
    printf("| %s \n",Title.c_str());
    printf("| [%s] \n",Param.c_str());
    printf("+-----------------------------\n");
}

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


//**************************************************************************
//  Copying template MS from Master Repository.
//
//   This is for Some testing items which must contain test data in the MS.
//   After copying, Modifying fuction for MS is executed depending on
//   Test cases/items.
//***************************************************************************

static const String DefaultLocalMsName = "./sdimaging-t.ms";

void CopyDefaultMStoWork()
{
    //  Environment //

        RunEnv env;

    // Src/Dsr Path (string) 

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

        if(true)
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

void DeleteWorkingMS()
{
    String dst         = DefaultLocalMsName;

    casacore::Path        path(dst);
    casacore::Directory   dir_ctrl(path);

    Description( "Deleting Working MeasurementSet for modifed use." ,dst.c_str() );

    // Delete File (Recursively done) 
    // NOTE: for debug use, please change true-> falase )

    if (true)
    {
         dir_ctrl. removeRecursive(false /*keepDir=False */ );
    }
}

//+
// Editting local MS for detail Test
// (Tables) ANTENNA and POINTING . 
//-

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

//+
// MeasurementSet Edit Class (MsEdit)
//  - Functions are defines.
//-

class MsEdit 
{
public:

    MsEdit()        { };

    // Add or Remove Row  (Antenna) //

        uInt AntennaTable_AppendRow();
        void AntennaTable_RemoveRow(uInt NRow );
        
    // List Table Contents. //

        void AntennaTable_List(String MsName =DefaultLocalMsName);
        void PointingTable_List(String MsName =DefaultLocalMsName, bool showAll=false );

    // Write Data on Antenna Table // 
     
        void AntennaTable_WriteData(String MsName =DefaultLocalMsName, uInt Row =0 );

    // Write new Columns and init data // 

        void PointingTable_WriteData(String MsName =DefaultLocalMsName );

    // Write (generated) Test Data on Pointing Table //

        void WriteTestDataOnPointingTable(String MsName =DefaultLocalMsName );

    // Write (generated) Test Data in MAIN Table //

        void WriteTestDataOnMainTable(String MsName =DefaultLocalMsName);

    // Add or Remove Column (Pointing) ////

        void CreateNewColumnsFromDirection();

    // Generating artificial POINTING values

        casacore::Vector<double> GeneratePointngForInterporation(double tm);

    // SetUp Evaluation parameters for POINTING data //
    
        void SetUpEvaluationParametersInPointing(Int param, Int N );

    //+
    // Default File Name
    //-

        const String MsName = DefaultLocalMsName;  // default (C++11) //

    //+
    // Buff between table and local var.
    //-

        ANTENNADataBuff  AntennaData;   // for Read 
        ANTENNADataBuff  AntennaData1;  // for Write

    //+
    // Interpolation Test
    //-

        bool useRealData  = false;   // True = use real data in the MS. False = use generated data

        double ErrorWorst_1 = 0.0;
        double ErrorWorst_2 = 0.0;

    //+
    // Generating Parameter 
    //   (assume AZEL)
    //    -PI <  AZ < PI,  -PI/2 < EL < PI/2
    //-

        Double Dir_X_start      = -3.14;
        Double Dir_X_increment  = 0.001;

        Double Dir_Y_start      = -0.157;
        Double Dir_Y_increment  = 0.0002;

    // Time Offset and Expected Result //

        Double  ExpectedOffset_X = 0.0;    // 0.0005;     
        Double  ExpectedOffset_Y = 0.0;    // 0.0005;

    // Parameters  (Time)
    
        Double slideOffset   = 0;      // Sliding Time (sec)
        Double dayOffset     = 0;      // Day offset   (Day)

    ////////////////////
    
        // Interval and shift 

        Double pointingInterval = 1.0 ;	        // Interval Time to set in POINTING
        Double TestOffset       = 0.49;          //   0 < TestOffset < Interval 

        // Interporation Error Limit //
    
        Double interpolationErrorLimit = 5.0E-9;
    
    ////////////////////

};

//+
// Add onmsedit. erow on Antanna Table
//  returns latest nrow.
//-

uInt  MsEdit::AntennaTable_AppendRow()
{
    // Measurement Set (use default name) //

        MeasurementSet ms0( MsName.c_str(), casacore::Table::TableOption:: Update );

    // Table handle //

        MSAntenna   hAntennaTable = ms0.antenna();

    // Add Row //

        hAntennaTable.addRow();
        uInt nrow = hAntennaTable.nrow();

        return nrow;
}

//+
// Remove specified row from Antanna Table
//-

void MsEdit::AntennaTable_RemoveRow(uInt nrow )
{
    // Measurment Set (use default name ) 

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Table handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Remove Row //
      
         hAntennaTable.removeRow( nrow );

}

//+
// List all the data from Antanna Table
//  of specified MS.
//- 

void MsEdit::AntennaTable_List(String MsName )
{

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String AntennaTableName = ms0.antennaTableName(); 
        printf("Antanna  Table name \n" );
        printf(" [%s] \n",AntennaTableName.c_str());

    // Prepeare Handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Get current row count //

        uInt nrow_a = hAntennaTable.nrow();

        printf( "Antenna Table nrow  =%d \n",nrow_a);

    //+
    // Get Column handle from Table  (Antenna)
    //  NEW FEASURE by C++11: Use samrt pointer.
    //-

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

        printf( "=====================================================\n");
        for (uInt row=0; row<columnAntenna->nrow(); row++)
        {

            AntennaData.name            =  antennaName.         get(row);
            AntennaData.station         =  antennaStation.      get(row);
            AntennaData.type            =  antennaType.         get(row);
            AntennaData.mount           =  antennaMount.        get(row);
            AntennaData.position        =  antennaPosition.     get(row);
            AntennaData.offset          =  antennaOffset.       get(row);

            AntennaData.dish_diameter   =  antennaDishDiameter. get(row);

            printf( "Antenna[%2d]: name     [%s]\n",row,  AntennaData.name.    c_str() );
            printf( "Antenna[%2d]: station  [%s]\n",row,  AntennaData.station. c_str() );
            printf( "Antenna[%2d]: type     [%s]\n",row,  AntennaData.type.    c_str() );
            printf( "Antenna[%2d]: mount    [%s]\n",row,  AntennaData.mount.   c_str() );

            printf( "Antenna[%2d]: position [%f,%f,%f] \n",row,  
                                 AntennaData.position[0],
                                 AntennaData.position[1], 
                                 AntennaData.position[2]  );

            printf( "Antenna[%2d]: offset   [%f,%f,%f] \n",row,
                                 AntennaData.offset[0],
                                 AntennaData.offset[1],
                                 AntennaData.offset[2]  );


 
            printf( "Antenna[%2d]: dish diameter  [%f]\n",row,  AntennaData.dish_diameter );

            printf( "------------------\n");

        } // end for

}

//+
//  Write Data to Antenna Table 
//-

void MsEdit::AntennaTable_WriteData(String MsName, uInt Row )
{
    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String AntennaTableName = ms0.antennaTableName(); 
        printf("Antanna  Table name \n" );
        printf(" [%s] \n",AntennaTableName.c_str());

    // Prepeare Handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Get current row count //

        uInt nrow_a = hAntennaTable.nrow();

        printf( "Antenna Table nrow  =%d \n",nrow_a);

    //
    // Get Column handle from Table  (Antenna)
    //

        std::unique_ptr<casacore::MSAntennaColumns> 
                columnAntenna( new casacore::MSAntennaColumns( hAntennaTable ));

    //+
    // LIST
    //-

    // Special Cplun //
        
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
 // List ALL THE DATA  from Pointing Table
 //  of specified MS.
 //- 

void MsEdit::PointingTable_List(String MsName, bool showAll)
{
    // Open MS by Update mode //

//        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );
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

        printf("Loop Start.\n" );
        for (uInt row=0; row<nrow_p; row++)
        {
            if(showAll)
            {
                printf( "Pointing: Antenna ID  [%d] = %d \n", row, pointingAntennaId.    get(row)  );
                printf( "Pointing: Time        [%d] = %f \n", row, pointingTime.         get(row)  );
                printf( "Pointing: Interval    [%d] = %f \n", row, pointingInterval.     get(row)  );
                printf( "Pointing: Name        [%d] = \"%s\" \n", row, pointingName.       get(row).c_str()  );
                printf( "Pointing: Num Poly    [%d] = %d \n", row, pointingNumPoly.      get(row)  );
                printf( "Pointing: Time Origin [%d] = %f \n", row, pointingTimeOrigin.   get(row)  );

            }

            Vector<Double> valDirection  =   pointingDirection. get(row);
            Vector<Double> valTarget     =   pointingTarget. get(row);
#if 0
            Vector<Double> valPointingOffset     =   pointingPointingOffset. get(row);
            Vector<Double> valSourceOffset       =   pointingSourceOffset.   get(row);
            Vector<Double> valEncoder            =   pointingEncoder.        get(row);
#endif 
            printf( "Pointing: Direction        [%d] = (%f,%f)  \n", row, valDirection[0],      valDirection[1] );
            printf( "Pointing: Target           [%d] = (%f,%f)  \n", row, valTarget[0],         valTarget[1] );
#if 0
            printf( "Pointing: Pointing Offset  [%d] = (%f,%f)  \n", row, valPointingOffset[0], valPointingOffset[1] );
            printf( "Pointing: Source   Offset  [%d] = (%f,%f)  \n", row, valSourceOffset[0],   valSourceOffset[1] );
            printf( "Pointing: encoder          [%d] = (%f,%f)  \n", row, valEncoder[0],        valEncoder[1] );
#endif 

#if 0
            Vector<Double> valPointingOffset     =   pointingPointingOffset. get(row);
            Vector<Double> valSourceOffset       =   pointingSourceOffset  . get(row); 
         
            printf( "Pointing: Pointing Offset [%d] = (%f,%f)  \n", row, valPointingOffset[0],    valPointingOffset[1] );
            printf( "Pointing: Source   Offset [%d] = (%f,%f)  \n", row, valSourceOffset[0],      valSourceOffset[1] );
#endif 
            printf( "------------------\n"); 

        }

}



//+
// Write Columns on Pointing Table  
//-

void MsEdit::PointingTable_WriteData(String MsName )
{

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

        Description("Adding Data on new columns in POINTING table. ",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        uInt  nrow_p = hPointingTable.nrow();
        printf( "Pointing Table nrow =%d \n",nrow_p);

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

        Description( "attempting to add Data on Pointing Table.", "Nrow="+std::to_string(nrow_p)  );

        for (uInt row=0; row<nrow_p; row++)
        {
            // set Shape of New added Colum //

            pointingPointingOffset. setShape(row, Ipo); 
            pointingSourceOffset.   setShape(row, Ipo);
            pointingEncoder.        setShape(row, Ipo);

            // prepare initial value //

            Array<Double> init_data1( Ipo, -1.0);
            Array<Double> init_data2( Ipo, -2.0);
            Array<Double> init_data3( Ipo, -3.0);

            // put initial data to Column //
            
            pointingPointingOffset. put( row, init_data1 );
            pointingSourceOffset.   put( row, init_data2 );
            pointingEncoder.        put( row, init_data3 );

        }

        // Flush // 
        ms0.flush();
 
}

//   copied by Direction Table 
//  
//    - This is for testing MovingSourceCorrection, 
//      only OFFSET tables are applied. 
//-

void MsEdit::CreateNewColumnsFromDirection()
{

 
    // Open MS by Update mode //

        String MsName = DefaultLocalMsName;
        
        Description( "Cpoied from Direction and Adding 3 Columns on Pointing Table ", MsName.c_str() );

        String name =  MsName;

        MeasurementSet ms0( name.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();
        printf("Pointing Table name is [%s] \n",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Prepare Column //

        // create the Smart Pointer in use. //

        std::unique_ptr<casacore::MSPointingColumns> 
                columnPointing( new casacore::MSPointingColumns( hPointingTable ));

    // each Column.. used in setDirectionColumn() //

       ArrayColumn< Double > colDirection       =  columnPointing->direction ();
       ArrayColumn< Double > colTarget          =  columnPointing->target ();
       ArrayColumn< Double > colPointingOffset  =  columnPointing->pointingOffset ();
       ArrayColumn< Double > colSourceOffset    =  columnPointing->sourceOffset ();
       ArrayColumn< Double > colEncoder         =  columnPointing->encoder ();


    //+
    // Attempt to create new Column.
    //-

    // Table Desc linked from MS //

       TableDesc  tblDsc = hPointingTable.tableDesc();

    //+
    // Show current Columns. 
    //-
        Description("=== Current Columns ======","");

        Vector<String> col_name = tblDsc.columnNames();
        
        for (uInt i=0 ; i<col_name.size() ; i++)
        {
                printf( "%s \n", col_name[i].c_str() );
        }

    // Copy Column Descriptor //

        ColumnDesc  OrgColumnDesc = tblDsc.columnDesc ( "DIRECTION" ) ;    

        ColumnDesc  RevColumnDesc1 = ColumnDesc(OrgColumnDesc);
        ColumnDesc  RevColumnDesc2 = ColumnDesc(OrgColumnDesc);
        ColumnDesc  RevColumnDesc3 = ColumnDesc(OrgColumnDesc);

    // Coumns Name def //
    
        String colname1 = "POINTING_OFFSET";
        String colname2 = "SOURCE_OFFSET";
        String colname3 = "ENCODER";

        RevColumnDesc1.setName( colname1 );
        RevColumnDesc2.setName( colname2 );
        RevColumnDesc3.setName( colname3 );

        Description( "Adding 3 Columns on Pointing Table ", MsName.c_str() );

        if( ! tblDsc.isColumn( colname1 ) )
        {
            printf(" Intended Table Not exist, Attempt to add Column.\n" );
      
            // Add Revied Column to TABLE //

            hPointingTable.addColumn(RevColumnDesc1);
            hPointingTable.addColumn(RevColumnDesc2);
            hPointingTable.addColumn(RevColumnDesc3); 
        }

#if 0
      // Remove Colume (by Name )

          hPointingTable.removeColumn(colname1);
          hPointingTable.removeColumn(colname2);
          hPointingTable.removeColumn(colname3);
#endif 
  
     // Flush ..//

         ms0.flush(); 
 
     //+
     //  Option:: Access data (shown by Array) 
     //   - If need to inspect , write a code here.
     //-

#if 0  
         Array<Double> dir  = colDirection .get(0);
         Array<Double> tar  = colTarget    .get(0);
         Array<Double> ptOff  = colPointingOffset .get(0);
         Array<Double> scOff  = colSourceOffset    .get(0);
         Array<Double> enc    = colEncoder        .get(0);
#endif 

      printf(" Intended Change completed.\n" );

}

//+
// Generate Pseudo Data on Direction
//   returns pointng info. by Vecror.
//-

casacore::Vector<double> MsEdit::GeneratePointngForInterporation(Double tn )

{
        casacore::Vector<Double> point;
        point.resize(6);


        Double d =  (22 *3600.0 
                     +  5*60 +  41.5 
                     + (tn * pointingInterval) + dayOffset )/(3600*24) 
                     + slideOffset;
     
        casacore::MVTime  tm00(2003,11,12 ,d);     
 
        //+
        // Designed Function of POINITNG location
        //-

        switch(1)
        {
            case 1:    /// 1.st order ///
        
                point[0] = Dir_X_start + Dir_X_increment * tn;     // Direction
                point[1] = Dir_Y_start + Dir_Y_increment * tn;    // Direction
                break;

            case 2:       /// 2.nd oder /// 
   
                point[0] = Dir_X_start + pow((Dir_X_increment * tn), 2) /(2.0*M_PI);     // Direction
                point[1] = Dir_Y_start + pow((Dir_Y_increment * tn), 2) /(1.0*M_PI);    // Direction
                break;

            case 3:       /// 3.rd. order ///
                
                point[0] = pow((Dir_X_increment * (tn-1400)) , 3)/(2.0*M_PI*M_PI);  
                point[1] = pow((Dir_Y_increment * (tn-1400)) , 3)/(2.0*M_PI*M_PI);

                break;

            case 0:
            default:
                point[0] = 0.0;
                point[1] = 0.0;
                break;
        }

        // Time and Interval //

        point[2] = tm00.second();                  // Time  (sec) 
        point[3] = pointingInterval;               // Interval (sec)
        
        //+
        // (Reserved)
        //-
        
        point[4] = 0.0 ;
        point[5] = 0.0 ;

        return point;
}

void MsEdit::SetUpEvaluationParametersInPointing(Int Param, Int Num )
{
        Description("Setting Up parameters ","" );
  
        //+
        // Generating Parameter 
        //   (assume AZEL)
        //    -PI <  AZ < PI,  -PI/2 < EL < PI/2
        //-

            Dir_X_start      = -3.14;
            Dir_X_increment  = 0.0016;

            Dir_Y_start      = -1.57;
            Dir_Y_increment  = 0.0008;

        //+
        // Define Interval and Interpolation Condition // 
        //-
            pointingInterval             = 1.0 ;                          // Interval Time to set in POINTING
            TestOffset                   = (double)Param / (double) Num;  // Offset bwtween 2 Samples. ( 0 < TestOffset < Interval) 

            interpolationErrorLimit      = 1.0e-02 ;  // 1.2e-8;        // Error Limit

        //+
        // SetUp Threashold  Condition and MeasurmentSet
        //-

          ExpectedOffset_X =  Dir_X_increment * ( TestOffset / pointingInterval );
          ExpectedOffset_Y =  Dir_Y_increment * ( TestOffset / pointingInterval );

}

//+
// Wtite Test Data on Direction Column in Pointing Table
//  - Values are got by sub fuction above.
//  - SetUp() in class TestDirection calls this.
//  - see also TEST Fixture  
//-

void  MsEdit::WriteTestDataOnPointingTable(String MsName)
{
    Description( "Writing Test Data on Direcotion Column in Pointing Table", 
                  MsName.c_str()  );

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        uInt nrow_p = hPointingTable.nrow();

        printf( "Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        // create the Smart Pointer in use. //

        std::unique_ptr<casacore::ROMSPointingColumns> 
                columnPointing( new casacore::ROMSPointingColumns( hPointingTable ));

    //+
    // Time Info
    //-

        ROScalarColumn<Double> pointingTime           = columnPointing ->time();
        ROScalarColumn<Double> pointingInterval       = columnPointing ->interval();
 
    //+
    // Listing Columns(Direction related) on Pointing 
    //-

        ROArrayColumn<Double>  pointingDirection      = columnPointing ->direction();
        ROArrayColumn<Double>  pointingTarget         = columnPointing ->target();

        IPosition Ipo = pointingDirection.shape(0);
        printf(" - Shape of pointingDirection.[%ld, %ld] \n", Ipo[0], Ipo[1] );

        for (uInt row=0; row<nrow_p; row++)
        {

            // DIRECTION  //

                Array<Double> direction(Ipo, 0.0);   // IP shape and initial val // 

                Vector<Double>  psd_data  = GeneratePointngForInterporation( (Double)row ); // generated pseudo data. //

                direction[0][0] = psd_data[0];
                direction[0][1] = psd_data[1];

            // write access //

                pointingDirection.   put(row, direction );
                pointingTarget.      put(row, direction );

            // Time Info. (current) //
 
                Double curTime      = pointingTime.     get(row);
                Double curInterval  = pointingInterval. get(row);

            // Show Time (OPTION)//

            if(false) { 
                printf( "[%d] Curr / Generated Time, %f ,%f \n", 
                        row, curTime, psd_data[2] ); 
                printf( "     Generated Direction (%f ,%f) \n",
                        psd_data[0], psd_data[1] );
            }

            // Time Set  //

            if(useRealData)
            {
                pointingTime.           put(row, curTime     ); // copy curr. Time
                pointingInterval.       put(row, curInterval ); // copy curr. Interval 
            }
            else
            {
                pointingTime.           put(row, psd_data[2] + TestOffset); // Time
                pointingInterval.       put(row, psd_data[3] ); // Interval
            }
        }

        // Flush //
        
        ms0.flush();
        sleep(2);
}

//+
// Wtite Test Data on Direction Column in MAIN TABLE
//  - Values are got by common function.
//-

void  MsEdit::WriteTestDataOnMainTable(String MsName)
{
    Description( "Writing Test Data (Time) in MAIN Table", 
                  MsName.c_str()  );

    // Open MS by Update mode //

        MeasurementSet ms0( MsName.c_str(),casacore::Table::TableOption:: Update );

    // Get current row count //

        uInt nrow_ms = ms0.nrow();

        printf( "Main Table nrow =%d \n",nrow_ms);

    //+
    // Column::
    //   Time Info
    //-

        ROScalarColumn<Double> mainTime           ;
        ROScalarColumn<Double> mainInterval       ;
 
    // Attach ..//

        mainTime      .attach( ms0 , "TIME");
        mainInterval  .attach( ms0 , "INTERVAL");

    // WRite .. //

        for (uInt row=0; row<nrow_ms; row++)
        {
            // Pseudo Data (TEST DATA )

               Vector<Double>  psd_data  = GeneratePointngForInterporation( (Double)row ); // generated pseudo data. //

            // Time Info. (current) //
 
                Double curTime      = mainTime.     get(row);
                Double curInterval  = mainInterval. get(row);

            // Show Time //
  
                if (false){  
                        printf( "[%d] Curr / Interval, %f ,%f \n", 
                            row, curTime, curInterval ); 
                }

            // Time Set  //

            if(useRealData)
            {
                mainTime.           put(row, curTime ); // copy curr. Time
                mainInterval.       put(row, curInterval ); // copy curr. Interval 
            }
            else
            {
                mainTime.           put(row, psd_data[2] ); // Time
                mainInterval.       put(row, psd_data[3] ); // Interval
            }   

        }

        // Flush //
        
        ms0.flush();

}

// 
// Base TestClass
//-

class BaseClass : public ::testing::Test
{

public:

     RunEnv       env;

protected:

    uInt    ExpectedNrow = 0;   // C++11 feature //

    BaseClass()  { }

    ~BaseClass() { }

    virtual void SetUp()
    {
    }

    virtual void TearDown()
    {
    }


private:
 

};

//+
// PointingDirectionCalculator Class 
// Constructor
//-
class TestMeasurementSet : public BaseClass
{

public:

protected:

        TestMeasurementSet()
        {
        }

        ~TestMeasurementSet()
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

        void test_constructor(uInt num );

};

/*---------------------------------------
  attempt to open vaious Measurement Set
 ----------------------------------------*/

void TestMeasurementSet::test_constructor(uInt num )
{

    String name = env.getCasaMasterPath()+getMsNameFromList( num );
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
    TestDescription( "CALC Constructor by various MS" );

    size_t max_count = getMSCountFromList();
    printf("- %zu MeasruementSet for this test are ready.\n" , max_count); 

    for(uInt k=0; k< max_count ; k++)
    {
       FunctionalDescription( "CALC Constructor by various MS",getMsNameFromList(k).c_str()  );
 
        if (getMSThrowFromList(k))
        {
          EXPECT_ANY_THROW( test_constructor( k ) );
        }
        else
        {
          EXPECT_NO_THROW( test_constructor( k ) );
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

    Vector<int> RowIdList(nrow0,0);

    for(int sw=1; sw <=6; sw++)
    {
         //+
         // setlectData() key 
         //  by Antenna Pair
         //-

         String AntSel ="";
         if (sw==1)
         { 
             AntSel = "DV01&&DV01";
         }
         if (sw==2)
         { 
             AntSel = "DV02&&DV02";
         }
         if (sw==3)
         { 
             AntSel = "PM03&&PM03";
         }
         if (sw==4)
         { 
             AntSel = "DV01&&DV02";
         }
         if (sw==5)
         { 
             AntSel = "DV01&&PM03";
         }
         if (sw==6)
         { 
             AntSel = "DV02&&PM03";
         }

         Description("Testing RowId functions." , "case="+std::to_string(sw) );

        // getRowIdForOriginalMS //

          Vector<uInt> vRowIdOrgMS = calc.getRowIdForOriginalMS();

         printf( " calling selectData() \n");     
         if( sw !=0 ) 
	 {
             calc.selectData( AntSel,
                              "","","","","","","","","" );
         }

        // Nrow from Selected //

        uInt nrow = calc.getNrowForSelectedMS();
        printf("Selected nrow =%d\n",nrow);

        // Vecrtor<uInt> getRowID() 
        
          Description("(1) Vector<uInt> getRowId() ", name );
          Vector<uInt> vRowId = calc.getRowId();


        // getRowId( int ) //
 
          Description("(2) uInt getRowId() ", name );

        // Show and Verify //
        //  
          printf( "Num Row (Org) = %d \n" , nrow0 );
          printf( "Num Row (Sel) = %d \n" , nrow );
          printf( "    key=[%s]\n", AntSel.c_str() ); 

          for (uInt k=0; k < nrow; k++)
          {
              uInt RowId = calc.getRowId(k);
              printf( "RowID = [%d] \n", RowId ); 
              
              // Check List //
             RowIdList[RowId] ++;
          }
    }
    
    //*
    // RowId List 
    //  - checck all the ID have been appeared once,
    //  - Duplicattion IS NOT Allowed
    //* 

    for (uInt i=0; i < nrow0; i++)
    {
        // printf( " RowIdList[%d]  = %d \n" ,i, RowIdList[i] );
        EXPECT_EQ( RowIdList[i],1 );    
    }
}

/*---------------------------------------
  getDirection  base class
 ----------------------------------------*/
class TestDirection : public BaseClass
{

public:

protected:

        casa::MsEdit  msedit;

        // Add 3 OFFSET Colums (basically reserved) //
        void addColumnsOnPointing();

        // Add avalable data on OFFSET columns  //
        void addColumnDataOnPointing();

        // Add Testing Data(generated) to direction on POINTING //
        void addTestDataForGetDirection();

        // subfunction of TEST_F(TestDirection....)
        void subTestDirection(Double dt);


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
            addColumnDataOnPointing();
            addTestDataForGetDirection();

        }

        virtual void TearDown()
        {

            BaseClass::TearDown();

            // Delete Working MS 

            DeleteWorkingMS();

        }

};

void TestDirection::addTestDataForGetDirection()
{
    // Pointing //

    msedit.WriteTestDataOnPointingTable();

    // MAIN //

    msedit.WriteTestDataOnMainTable();
 
}

void TestDirection::addColumnsOnPointing()
{
    msedit.CreateNewColumnsFromDirection();
}

void TestDirection::addColumnDataOnPointing()
{
    msedit.PointingTable_WriteData( DefaultLocalMsName );
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

static void assert_accessor()
{
#if 0
    void *  pAccessor = calc.getAccessor();
    printf( "#   Accessor [%p] \n", pAccessor );
     
#endif

}


TEST_F(TestDirection, setDirectionColumn  )
{

    TestDescription( "setDirectionColumn (String Fram)" );
    const String MsName = DefaultLocalMsName;

    if(false){
        printf( "Listing all POINTING TABLE  \n");
        msedit.PointingTable_List( MsName, false );
    }
 
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

        ColName = "DIRECTION";
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();

        ColName = "TARGET";
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();       
 
        ColName = "POINTING_OFFSET";    // NEED to ADD Table in advance  //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();

        ColName = "SOURCE_OFFSET"; // NEED to Add Table in advance //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();
 
        ColName = "ENCODER";      // NEED to add Table in advance  //
        Description("Column Name" , ColName );
        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();

        ColName = "hogehoge";
        Description("Column Name" , ColName );
        EXPECT_ANY_THROW( calc.setDirectionColumn( ColName ) );

        assert_accessor();
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

    // List all info on Pointing Table. //

    if(true) {
          printf( "Listing all POINTING TABLE  \n");
          msedit.PointingTable_List( MsName, false );
    }

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
                Description("setMovingSourceConvert()", src);
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

        for(unsigned int k=0; k < ColName.size(); k++)
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
   
    // List all the point ..//
#if 0
        printf( "Listing all POINTING TABLE  \n");
        msedit.PointingTable_List( MsName, false );
#endif 
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

    // List all the point ..//
     
        if(false)
        {
            printf( "Listing all POINTING TABLE  \n");
            msedit.PointingTable_List( MsName, false );
        }

    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
    
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

#if 1

    Description("Dump obtined Direction info. ","" );

    for (uInt row=0; row< n_row; row++)
    {
        // Direction //

        double Val_1 = DirList1(row,0);
        double Val_2 = DirList1(row,1);

        casacore::MDirection  MovDir  = calc.getMovingSourceDirection();
        String strMovDir = MovDir.toString();

        printf(    "Dir at, %d, %f,%f, [Mov:%s]  \n",  
                row, Val_1, Val_2, strMovDir.c_str() );

    }

#endif 

}

/*--------------------------------------
   Interpolation Test in getDirection
 ---------------------------------------*/

//+
// Sub Fucntion
//-


void DumpPointingTable(String MsName)
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
    
void TestDirection::subTestDirection( Double dt )
{

    TestDescription( "Interpolation Test in getDirection ()" );

    const String MsName = DefaultLocalMsName;

    //  MS (Dump) ** develper option ** //

        if(false)  DumpPointingTable(MsName);

    // Create Object //
    
        MeasurementSet ms( MsName.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

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
  
        String ColName = "DIRECTION"; 
        Description( "getDirectionColumn()", ColName  );

        EXPECT_NO_THROW( calc.setDirectionColumn( ColName ) );

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

    //  Dump (before)

        if(false)
        {
            Description("Dump Pointing (before getDirection) ","" );
            DumpPointingTable( MsName );
        }

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

    Description("Inspecting  Direction Inerpolation ","" ); 

    if(true)
    {
        for (uInt row=1; row< n_row; row++)  // ACTUNG !!! start from 1 /// 
        {
            // Direction //

              double calculated_1 = DirList1(row,0);
              double calculated_2 = DirList1(row,1);

    
            // Generated Data for Test //
            
                casacore::Vector<double>  gen_out 
                    = msedit.GeneratePointngForInterporation((Double)row - dt);
            
            //+
            // Google Test 
            //-


                double generated_1 = gen_out[0]; 
                double generated_2 = gen_out[1]; 
                
                double Err_1 = calculated_1 - generated_1 ;   
                double Err_2 = calculated_2 - generated_2 ;

                double absErr_1 = abs(Err_1);
                double absErr_2 = abs(Err_2);
 
                if(absErr_1 > msedit.ErrorWorst_1 ) msedit.ErrorWorst_1 = absErr_1;
                if(absErr_2 > msedit.ErrorWorst_2 ) msedit.ErrorWorst_2 = absErr_2;

#if 1
		EXPECT_LE( absErr_1, msedit.interpolationErrorLimit  ); 
                EXPECT_LE( absErr_2, msedit.interpolationErrorLimit  ); 
#else
                ASSERT_LE( absErr_1, msedit.interpolationErrorLimit  ); 
                ASSERT_LE( absErr_2, msedit.interpolationErrorLimit  ); 
#endif
            // Output //

                printf( "----\n");
                printf( "Main Table Dir [%6d], %12.9f,%12.9f \n",   row,  calculated_1, calculated_2 );
                printf( "Generated data [%6d], %12.9f,%12.9f \n",   row,  generated_1,  generated_2 );
                printf( "Numerical error[%6d],%5.2e,%5.2e \n",      row,  absErr_1,     absErr_2);
        }    
    }
   
    printf( "Error Worst 1 = %e \n", msedit.ErrorWorst_1 );
    printf( "Error Worst 2 = %e \n", msedit.ErrorWorst_2 ); 

}

TEST_F(TestDirection, Interpolation0 )
{

    TestDescription( "Interpolation test in getDirection() " );

    // Reset Statictic //

    msedit.ErrorWorst_1 = 0.0;
    msedit.ErrorWorst_2 = 0.0;

    // Test Loop 
    for (uInt loop=1; loop <= 9; loop ++ )
    {

        // Set Up Parapeters 

           msedit.SetUpEvaluationParametersInPointing( loop, 10 );

        // SetUp Testing  MeasurmentSet

           Description("Making MeasurementSet","k="+to_string( (double)loop/10.0) );
           addTestDataForGetDirection();

        // Executtion ..// 
        
          Description("Execution starts. ","" );
          subTestDirection( (double)loop/10.0 );

    }

}

TEST_F(TestDirection, getDirectionExtended )
{

    TestDescription( "getDirection (J2000) with selected data. uvw available" );

    // Use DefaultMS as a simple sequence.
    // Use the below to show uv valuses from MS.

#if 0
    const String MsName = DefaultLocalMsName;    
#else
    const String MsName = env.getCasaMasterPath() + "listobs/uid___X02_X3d737_X1_01_small.ms";
#endif 

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
        casacore::Vector<double>  val_uvw = uvwColumn.get(RowId);

        double u = val_uvw[0];
        double v = val_uvw[1];
        double w = val_uvw[2];

        // Direction //

        double Val_1 = DirList1(row,0);
        double Val_2 = DirList1(row,1);

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

//+
// Execution of selectData(...)
//  - args are given from external variables
//  - Note that 'calc' MUST be given by reference.
//-

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

//
// Frame definitioan for  setFrame()
// - Strings and a bool var. wheather it is available in casacore.
// - Some of them causes Exeption Message from inside.
//

struct _FrameName {
    bool   available;
    String name;
};

class TestSetFrame : public BaseClass
{

public:

protected:

const std::vector<struct _FrameName> DefinedFrametypes
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

        TestSetFrame ()
        {
        }

        ~TestSetFrame ()
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



}  // END namespace


/************************************************
   Unit Test Main (google test)
    - Based on instructed Template for GTest.
    - Such minimum statements are recommended.
 *************************************************/

int main (int nArgs, char * args [])
 {

   // Initialize //

    ::testing::InitGoogleTest(& nArgs, args);

   // Run Test //

    return (RUN_ALL_TESTS()) ;

}

