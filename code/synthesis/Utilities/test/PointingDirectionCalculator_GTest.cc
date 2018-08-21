//# SDPosInterpolator_GTest.cc: this defines unit tests of
//# SDPosInterpolator using google test framework
//#
//# Copyright (C) 2016
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

// Nishie include for TABLE Manupilation //

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
using namespace casacore;
using namespace std;


//+
// Additional CASACORE
// Added by Nishie //
//-


#include <cstdio>
#include <casa/OS/EnvVar.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

namespace casa {

//+
// Environment 
//-   
 
struct _MSDef
{
    bool   ExThrow;
    String name;
};

const struct _MSDef
    MSNames[30] = {
        
        // Exeption(bool)  , Filename //
        {true, "./sdimaging-t.ms"                    },
        {false, "sdimaging/sdimaging.ms"                      },
//        {true, ".sis14_twhya_calibrated_flagged.ms"        },
//        {true, ".sis14_twhya_calibrated_flagged-t.ms"      },
        {false, "listobs/uid___X02_X3d737_X1_01_small.ms" },
        {true,  "concat/input/A2256LC2_4.5s-1.ms"               },
        {true,  "concat/input/A2256LC2_4.5s-2.ms"               },
        {false, "sdimaging/Uranus1.cal.Ant0.spw34.ms" },
        {false, "sdimaging/Uranus2.cal.Ant0.spw34.ms" },
        {false, "sdimaging/azelpointing.ms"                   },
        {false, "sdimaging/clipping_1row.ms"          },
        {false, "sdimaging/clipping_2rows.ms"         },
        {false, "sdimaging/clipping_3rows.ms"         },
        {false, "sdimaging/clipping_3rows_2chans.ms"  },
        {false, "sdimaging/clipping_3rows_suprious.ms"        },
/*15*/  {false, "sdimaging/pointing6.ms"                      },
        {false, "sdimaging/sdimaging_flagtest.ms"             },
        {false, "sdimaging/selection_intent.ms"               },
        {false, "sdimaging/selection_misc.ms"         },
/*21*/  {false, "sdimaging/selection_spw.ms"          },
        {false, "sdimaging/selection_spw_unifreq.ms"  },
        {true,  "sdimaging/hogehoge.ms"  },
        {false, "",  }
    };

//+

//-

const String  getMSName(int No )
{
    return MSNames[No].name;
}
bool  getMSThrow(int No )
{
    return MSNames[No].ExThrow;
}

//+
//  Log Title
//-

void TestTitle( const String &Title )
{
    printf( "====================================== \n");
    printf( " %s  \n",Title.c_str() );
    printf( "====================================== \n");
}

void Title(const String &Title, const String &Param)
{
    printf("##########################\n");
    printf("# %s \n",Title.c_str());
    printf("# [%s] \n",Param.c_str());
    printf("##########################\n");
}


//+ 
// Enviromnent
//-

class MyEnv
{

public:

    MyEnv()
    {
        //+
        // Path and directory Position by Env. Variable.
        //-
        
        CasaPath        = GetCasaPath( "CASAPATH" );
        CasaMasterPath  = CasaPath + "regression/unittest/"; 
        CasaMaserMSFullName = CasaMasterPath + "sdimaging/sdimaging.ms";

        printf("MyEnv:: Environment Variable Information -----\n" );
        printf("CASAPATH      :%s \n", CasaPath.c_str());
        printf("CasaMasterPath:%s \n", CasaMasterPath.c_str());
    }
    
    String GetCasaPath(const String pathname )
    {
        if (casacore::EnvironmentVariable::isDefined(pathname)) {
            string casapath = casacore::EnvironmentVariable::get(pathname);
            size_t endindex = casapath.find(" ");
            if (endindex != string::npos) {
                string casaroot = casapath.substr(0, endindex);
                cout << pathname << "=" << casaroot << endl;
                return (casaroot);
             } else {
                cout << "hit npos" << endl;
                return "/hoge/";
             }
        } else {
            cout << "ERROR: Specified path " << pathname  << " is not defined" << endl;
            return "";
        }
    }

    String getCasaPath()
    {
        return CasaPath;
    }
    String getCasaMasterPath()
    {
        return CasaMasterPath;
    }


private:

    String CasaPath;            // translated from CASAPATH 
    String CasaMasterPath;

    String CasaMaserMSFullName;

};


//+
//  Copying File from Master Repository.
//-

void CopyDefaultMStoWork()
{
//  Environment //

    MyEnv env;

// TEST of Directory.copy() //

String src = env.getCasaMasterPath() + "sdimaging/sdimaging.ms";
String dst = "./sdimaging-t.ms";

     casacore::Path        sourcePath(src);
     casacore::Path        targetPath(dst);
       
     casacore::Directory   dir_ctrl(sourcePath);

     Title( "Copying Default MeasurementSet for modifed use to; " ,dst.c_str() );

     printf( " - src filespec  : %s \n", src.c_str() );
     printf( " - dest filespec : %s \n", dst.c_str() );


#if 1
     dir_ctrl.copy  (  targetPath,
                        True,    // Overwrite 
                        True  ); // Users permisssion
#endif 


}

void DeleteWorkingMS()
{
    String dst         = "./sdimaging-t.ms";
    casacore::Path        path(dst);
    casacore::Directory   dir_ctrl(path);

     Title( "Deleting Working MeasurementSet for modifed use." ,dst.c_str() );

    // Delete File (Revursively done) //

    if (false) dir_ctrl. removeRecursive(false);

}



typedef  struct ChgAntennaTable_ {

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

class MSEdit  : public MyEnv
{

public:
        MSEdit()        { };

        // Env: File Path //
        
        MyEnv env;
        
        // Add or Remove Row  (Antenna) //

        int  AntennaTable_Add_Row();
        void AntennaTable_Remove_Row(int NRow );
        
        // List Table //

        void ListAntennaTable(String MSname ="./sdimaging-t.ms");
        void ListPointingTable(String MSname ="./sdimaging-t.ms", bool showAll=false );

        // Write Data on Table // 
     
        void WriteDataAntennaTable(String MSname ="./sdimaging-t.ms", int Row =0 );

        // Write new Columns and init data // 

        void WriteColumnsPointingTable(String MSname ="./sdimaging-t.ms" );

        // Write (generated) Test Data on Pointing Table //

        void WriteTestDataOnDirection(String MSname ="./sdimaging-t.ms" );

        // Add or Remove Column (Pointing) ////

        void AddOffsetColumnsToPointingTable();

private:

        String MSname = "./sdimaging-t.ms";  // default //

        //+
        // Buff between table and local var.
        //-

         ANTENNADataBuff  AntennaData;   // for Read 
         ANTENNADataBuff  AntennaData1;  // for Write

};

int  MSEdit::AntennaTable_Add_Row()
{
    // Measurement Set (use default name) //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Table handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Add Row //

         hAntennaTable.addRow();

        int nrow = hAntennaTable.nrow();

        return nrow;
}


void MSEdit::AntennaTable_Remove_Row(int NRow )
{
    // Measurment Set (use default name ) //
        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Table handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Add Row //
      
         hAntennaTable.removeRow( NRow );

}





void MSEdit::ListAntennaTable(String MSname )
{

    // Open MS by Update mode //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String AntennaTableName = ms0.antennaTableName(); 
        printf("Antanna  Table name \n" );
        printf(" [%s] \n",AntennaTableName.c_str());

    // Prepeare Handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Get current row count //

        int nrow_a = hAntennaTable.nrow();

        printf( "Antenna Table nrow  =%d \n",nrow_a);

    //
    // Get Column handle from Table  (Antenna)
    //

        casacore::MSAntennaColumns     * columnAntenna;
        columnAntenna = new casacore::MSAntennaColumns( hAntennaTable );

    // LIST
    //  as follows.  Number of Rows = nrow()
    //  get() privides Native data in C   
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
        for (int i=0; i<(int)columnAntenna->nrow(); i++)
        {

            AntennaData.name            =  antennaName.         get(i);
            AntennaData.station         =  antennaStation.      get(i);
            AntennaData.type            =  antennaType.         get(i);
            AntennaData.mount           =  antennaMount.        get(i);
            AntennaData.position        =  antennaPosition.     get(i);
            AntennaData.offset          =  antennaOffset.       get(i);

            AntennaData.dish_diameter   =  antennaDishDiameter. get(i);

            printf( "Antenna[%2d]: name     [%s]\n",i,  AntennaData.name.    c_str() );
            printf( "Antenna[%2d]: station  [%s]\n",i,  AntennaData.station. c_str() );
            printf( "Antenna[%2d]: type     [%s]\n",i,  AntennaData.type.    c_str() );
            printf( "Antenna[%2d]: mount    [%s]\n",i,  AntennaData.mount.   c_str() );

            printf( "Antenna[%2d]: position [%f,%f,%f] \n",i,  
                                 AntennaData.position[0],
                                 AntennaData.position[1], 
                                 AntennaData.position[2]  );

            printf( "Antenna[%2d]: offset   [%f,%f,%f] \n",i,
                                 AntennaData.offset[0],
                                 AntennaData.offset[1],
                                 AntennaData.offset[2]  );


 
            printf( "Antenna[%2d]: dish diameter  [%f]\n",i,  AntennaData.dish_diameter );


          printf( "------------------\n");

        }



}


//+
//  Write Data to Antenna Table 
//-

void MSEdit::WriteDataAntennaTable(String MSname, int Row )
{
    // Open MS by Update mode //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String AntennaTableName = ms0.antennaTableName(); 
        printf("Antanna  Table name \n" );
        printf(" [%s] \n",AntennaTableName.c_str());

    // Prepeare Handle //

        MSAntenna   hAntennaTable  = ms0.antenna();

    // Get current row count //

        int nrow_a = hAntennaTable.nrow();

        printf( "Antenna Table nrow  =%d \n",nrow_a);

    //
    // Get Column handle from Table  (Antenna)
    //

        casacore::MSAntennaColumns     * columnAntenna;
        columnAntenna = new casacore::MSAntennaColumns( hAntennaTable );

    // LIST
    //  as follows.  Number of Rows = nrow()
    //  get() privides Native data in C   
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
// Write Antenna Table 
//-

void MSEdit::ListPointingTable(String MSname, bool showAll)
{

    // Open MS by Update mode //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();
        printf("Pointing Table name [%s] \n",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        int nrow_p = hPointingTable.nrow();

        printf( "Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        casacore::ROMSPointingColumns  *columnPointing;
        columnPointing =  new casacore::ROMSPointingColumns( hPointingTable );

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

        for (int i=0; i<nrow_p; i++)
        {
          if(showAll)
          {
                printf( "Pointing: Antenna ID  [%d] = %d \n", i, pointingAntennaId.    get(i)  );
                printf( "Pointing: Time        [%d] = %f \n", i, pointingTime.         get(i)  );
                printf( "Pointing: Interval    [%d] = %f \n", i, pointingInterval.     get(i)  );
                printf( "Pointing: Name        [%d] = \"%s\" \n", i, pointingName.       get(i).c_str()  );
                printf( "Pointing: Num Poly    [%d] = %d \n", i, pointingNumPoly.      get(i)  );
                printf( "Pointing: Time Origin [%d] = %f \n", i, pointingTimeOrigin.   get(i)  );

           }

           Vector<Double> valDirection  =   pointingDirection. get(i);
           Vector<Double> valTarget     =   pointingTarget. get(i);

           Vector<Double> valPointingOffset     =   pointingPointingOffset. get(i);
           Vector<Double> valSourceOffset       =   pointingSourceOffset.   get(i);
           Vector<Double> valEncoder            =   pointingEncoder.        get(i);


           printf( "Pointing: Direction        [%d] = (%f,%f)  \n", i, valDirection[0], valDirection[1] );
           printf( "Pointing: Target           [%d] = (%f,%f)  \n", i, valTarget[0],    valTarget[1] );
           printf( "Pointing: Pointing Offset  [%d] = (%f,%f)  \n", i, valPointingOffset[0], valPointingOffset[1] );
           printf( "Pointing: Source   Offset  [%d] = (%f,%f)  \n", i, valSourceOffset[0],   valSourceOffset[1] );
           printf( "Pointing: encoder          [%d] = (%f,%f)  \n", i, valEncoder[0],        valEncoder[1] );


#if 0
           Vector<Double> valPointingOffset     =   pointingPointingOffset. get(i);

           Vector<Double> valSourceOffset       =   pointingSourceOffset  . get(i); 
         
           printf( "Pointing: Pointing Offset [%d] = (%f,%f)  \n", i, valPointingOffset[0],    valPointingOffset[1] );
           printf( "Pointing: Source   Offset [%d] = (%f,%f)  \n", i, valSourceOffset[0],      valSourceOffset[1] );
#endif 
          printf( "------------------\n"); 
        }

}

//+
// 
// Write Columns on Pointing Table  
//-

void MSEdit::WriteColumnsPointingTable(String MSname )
{

    // Open MS by Update mode //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

        Title("Adding Data on new columns in POINTING table. ",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        int nrow_p = hPointingTable.nrow();

        printf( "Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        casacore::ROMSPointingColumns  *columnPointing;
        columnPointing =  new casacore::ROMSPointingColumns( hPointingTable );


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


        IPosition Ipo;
          Ipo  = pointingDirection.shape(0);
          printf(" - Shape of Direction.[%d, %d] \n", Ipo[0], Ipo[1] );


        for (int row=0; row<nrow_p; row++)
        {

//          printf( " - attempting to add a data [%d] ( on Testing) \n", row );

            // set Shape of New added Colum //

            pointingPointingOffset. setShape(row, Ipo); 
            pointingSourceOffset.   setShape(row, Ipo);
            pointingEncoder.        setShape(row, Ipo);

            // initial calue //
#if 0
             Array<Double> init_data0( Ipo,  0.1);
#endif 
            Array<Double> init_data1( Ipo, -1.0);
            Array<Double> init_data2( Ipo, -2.0);
            Array<Double> init_data3( Ipo, -3.0);


            // put initial data //
#if 0
            pointingDirection.      put( row, init_data0 );
            pointingTarget.         put( row, init_data0 );
#endif 
            pointingPointingOffset. put( row, init_data1 );
            pointingSourceOffset.   put( row, init_data2 );
            pointingEncoder.        put( row, init_data3 );

        }


        // Flush //
            Array<Double> init_data( Ipo, -1.0);
        
        ms0.flush();

}
 

//+
//   Add and Remove Column 
//    to/from POINTNG Table 
//-

void MSEdit::AddOffsetColumnsToPointingTable()
{

 
    // Open MS by Update mode //

        String MsName ="./sdimaging-t.ms";
        
        Title( "Adding 3 Columns on Pointing Table ", MsName.c_str() );

        String name =  MsName;

        MeasurementSet ms0( name.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();
        printf("Pointing Table name [%s] \n",PointingTableName.c_str());

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Prepare Column //

       casacore::MSPointingColumns     * columnPointing;
       columnPointing = new casacore::MSPointingColumns( hPointingTable );

    // each Column.. used in setDirectionColumn() //

       ArrayColumn< Double > colDirection       =  columnPointing->direction ();
       ArrayColumn< Double > colTarget          =  columnPointing->target ();
       ArrayColumn< Double > colPointingOffset  =  columnPointing->pointingOffset ();
       ArrayColumn< Double > colSourceOffset    =  columnPointing->sourceOffset ();
       ArrayColumn< Double > colEncoder         =  columnPointing->encoder ();


     //
     // Attempt to create new Column.
     //

    // Table Desc linked from MS //

    TableDesc  tblDsc = hPointingTable.tableDesc();


    //+
    // Show current Columns. 
    //-
        printf("=== Current Columns ====== \n");

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

       Title( "Adding 3 Columns on Pointing Table ", MsName.c_str() );


    if( ! tblDsc.isColumn( colname1 ) )
    {
      printf(" ! Intended Table Not exist, Attempt to add Column.\n" );
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


     //
     //  Access data (shown by Array) 
     //

       Array<Double> dir  = colDirection .get(0);
       Array<Double> tar  = colTarget    .get(0);
#if 0
       Array<Double> ptOff  = colPointingOffset .get(0);
       Array<Double> scOff  = colSourceOffset    .get(0);
       Array<Double> enc    = colEncoder        .get(0);
#endif 

      printf(" ! Intended Change completed.\n" );

}

casacore::Vector<double>  generatePseudoDirection(int i)
{
    casacore::Vector<Double> point;
    point.resize(4);

    point[0] = 0.0 + 0.01 * (double)i;      // Direction
    point[1] = 0.0 - 0.01 * (double)i;

    //+
    // THis value can change Interporation behabior
    //   - use 0, +10000.0, -10000.0
    //-  

    Double intentionalShift = 10000.0;
    Double interval = 2.99827;

    Double        d = (22 *3600.0 +  5*60 +  41.5 + (double)i * interval + intentionalShift )/(3600*24); 
                      

    casacore::MVTime  tm(2003,11,12 ,d);     
 
 
    point[2] = tm.second();    // Time 
    point[3] = interval;           // Interval 

        
    return point;
}

void  MSEdit::WriteTestDataOnDirection(String MSname)
{
        Title( "Writing Test Data on Direcotion Column in Pointing Table", 
                   MSname.c_str()  );

    // Open MS by Update mode //

        MeasurementSet ms0( MSname.c_str(),casacore::Table::TableOption:: Update );

    // Tables Name //

        String PointingTableName = ms0.pointingTableName();

    // Prepeare Handle //

        MSPointing  hPointingTable = ms0.pointing();

    // Get current row count //

        int nrow_p = hPointingTable.nrow();

        printf( "Pointing Table nrow =%d \n",nrow_p);

    //
    // Get Column handle from Table  (Pointing)
    //  

        casacore::ROMSPointingColumns  *columnPointing;
        columnPointing =  new casacore::ROMSPointingColumns( hPointingTable );

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

        ROArrayColumn<Double>  pointingPointingOffset = columnPointing ->pointingOffset();
        ROArrayColumn<Double>  pointingSourceOffset   = columnPointing ->sourceOffset();
        ROArrayColumn<Double>  pointingEncoder        = columnPointing ->encoder();


        IPosition Ipo = pointingDirection.shape(0);
        printf(" - Shape of pointingDirection.[%d, %d] \n", Ipo[0], Ipo[1] );

        for (int row=0; row<nrow_p; row++)
        {

             // DIRECTION  //

               Array<Double> direction(Ipo, -1);   // IP shape and initial val // 

                Vector<Double>  psd_data  = generatePseudoDirection( row );

                direction[0][0] = psd_data[0];
                direction[0][1] = psd_data[1];
 
             // write access //

                if(false)       pointingDirection.   put(row, direction );
                if(false)       pointingTarget.      put(row, direction );

             // Time Info.  //
             
                if(true)       pointingTime.           put(row, psd_data[2] ); // Time
                if(false)       pointingInterval.       put(row, psd_data[3] ); // Interval 
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

        uInt    ExpectedNrow = 0;

protected:

        BaseClass()
        {
        }

        ~BaseClass()
        {
        }


        virtual void SetUp()
        {
        }

        virtual void TearDown()
       {
       }

    MyEnv       env;

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
            printf ("TestMeasurementSet::SetUp:: called \n");
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
            printf ("TestMeasurementSet::TearDown:: called \n");
       }

       void test_constructor(int num );

};




/*---------------------------------------
  attempt to open vaious Measurement Set
 ----------------------------------------*/

void TestMeasurementSet::test_constructor(int num )
{

    String name = env.getCasaMasterPath()+getMSName( num );
    Title("Testing Construcror." , name );

    // CONSTRUCTOR  //

      MeasurementSet ms0( name  );      

      PointingDirectionCalculator calc(ms0);

    // Initial brief Inspection //
   

      printf("# Constuctor Initial Check.  [%s] \n", name.c_str()  );
      printf("   detected nrow : %d \n", (int)calc.getNrowForSelectedMS() );
      EXPECT_NE((int)0, (int)calc.getNrowForSelectedMS() );
    
}


TEST_F(TestMeasurementSet, variousConstructor )
{
    TestTitle( "CALC Constructor by various MS" );
 
    for(int k=0;k<30;k++)
    {
        if(getMSName(k)=="")    break;

        if (getMSThrow(k))
        {
          EXPECT_ANY_THROW( test_constructor( k ) );
        }
        else
        {
          EXPECT_NO_THROW( test_constructor( k ) );
        }

    }


}

class TestDirection : public BaseClass
{

public:

        casa::MSEdit  msedit;

        // Add 3 OFFSET Colums (basically reserved) //
        void addColumnsOnPointing();

        // Add avalable data on OFFSET columns  //
        void addColumnDataOnPointing();

        // Add Testing Data(generated) to direction on POINTING //
        void addTestDataOnDirection();

protected:

        TestDirection()
        {
        }

        ~TestDirection()
        {
        }


        virtual void SetUp()
        {
            BaseClass::SetUp();
            printf ("TestDirection::SetUp:: called \n");

            //+
            // Copy and add culumns, 
            //  init columns, 
            //  generate test dta
            //-
            
            CopyDefaultMStoWork();

            addColumnsOnPointing();
 
            addColumnDataOnPointing();
           
            addTestDataOnDirection();

        }

        virtual void TearDown()
       {

            BaseClass::TearDown();

            DeleteWorkingMS();

            printf ("TestDirection::TearDown:: called \n");
       }
    

};

void TestDirection::addTestDataOnDirection()
{
    msedit.WriteTestDataOnDirection();
}

void TestDirection::addColumnsOnPointing()
{
    msedit.AddOffsetColumnsToPointingTable();
}

void TestDirection::addColumnDataOnPointing()
{
    String MsName = "./sdimaging-t.ms";    //  

    msedit.WriteColumnsPointingTable( MsName );
}

TEST_F(TestDirection, setDirectionColumn  )
{

    TestTitle( "setDirectionColumn (String FrameName)" );
    String MsName = "./sdimaging-t.ms";    //  

            printf( "Listing all POINTING TABLE  \n");
            msedit.ListPointingTable( MsName, false );

 
    // MS name for this Test //
   
        String name =   MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );


    //
    // setDirectionColumm()  calls 
    //
        String Name;

        Name = "DIRECTION";
        Title("Column Name" , Name );
        EXPECT_NO_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }   
   
        
        Name = "TARGET";
        Title("Column Name" , Name );
        EXPECT_NO_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }


        Name = "POINTING_OFFSET";    // *** NEED to ADD in advance  //
        Title("Column Name" , Name );
        EXPECT_NO_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }

        Name = "SOURCE_OFFSET"; // *** NEED to Add in advance //
        Title("Column Name" , Name );
        EXPECT_NO_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }

        Name = "ENCODER";      // *** NEED to add in advance  //
        Title("Column Name" , Name );
        EXPECT_NO_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }

        Name = "hogehoge";
        Title("Column Name" , Name );
        EXPECT_ANY_THROW( calc.setDirectionColumn( Name ) );

        { void *  pAccessor = calc.getAccessor();
          printf( "#   Accessor [%p] \n", pAccessor );
        }
}


TEST_F(TestDirection, Matrixshape )
{

    TestTitle( "setDirectionListMatrixShape()" );
    String MsName = "./sdimaging-t.ms";    //  
    
    // MS name for this Test //
        String name =  MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    

       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    // A set of API Call   
    //    - setDirectionListMatrixShape() and getDirection() 
    //    - getDirection return Matrix<Double> and the shape is determined in 
    //      getDirection by IPPosition(p,q,r)
    //-

        Matrix<Double> DirList1;
        Matrix<Double> DirList2;
        int N_Col;  
        int N_Row ;  

    // COLUMN //
    
        Title("setDirectionListMatrixShape", "COLUMN_MAJOR" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );
        
        DirList1  = calc.getDirection();
        N_Col    = DirList1.ncolumn();
        N_Row    = DirList1.nrow();
        printf ("# NCol = %d , NRow = %d \n", N_Col, N_Row );

        EXPECT_EQ( N_Col, 2);

    // ROW  //

        Title("setDirectionListMatrixShape", "ROW_MAJOR");
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::ROW_MAJOR) );

        DirList2  = calc.getDirection();
        N_Col    = DirList2.ncolumn();
        N_Row    = DirList2.nrow();
        printf ("# NCol = %d , NRow = %d \n", N_Col, N_Row );

        EXPECT_EQ( N_Row, 2);

}

TEST_F(TestDirection, getDirection )
{

    TestTitle( "getDirection (J2000)" );
    String MsName = "./sdimaging-t.ms";    // 


    // MS name for this Test //
        String name =  MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((uInt)0, ExpectedNrow );

    //+
    // setDirectionColumn() 
    //-

        calc.setDirectionColumn("DIRECTION");

    //+
    //      getDirection
    //-
        Title("calling setDirectionListMatrixShape()" ,"Column Major" );
        EXPECT_NO_THROW( calc.setDirectionListMatrixShape(PointingDirectionCalculator::COLUMN_MAJOR) );

    //+
    //  getDirection() call
    //-

        Title("calling  getDirection() ","" );

        Matrix<Double>  DirList  = calc.getDirection();
        int  N_Row    = DirList.nrow();

        printf( "Number of Row = %d \n", N_Row );
        EXPECT_EQ( N_Row, ExpectedNrow);

    //+
    // Dump Matrix
    //-

#if 1
    FILE *fp = fopen( "Dir.csv", "w" );

    for (int row=0; row<N_Row; row++)
    {
        double Val_1 = DirList(row,0);
        double Val_2 = DirList(row,1);

         printf(    "Dir at, [%d], %f,%f \n",   row, Val_1, Val_2 );
        fprintf( fp,"Dir at, [%d], %f,%f \n",   row, Val_1, Val_2 );
 
    }

    fclose(fp);
#endif 



}


class TestSelectData : public BaseClass
{

public:
 
    String  AntSel              = "";      // Should be Null, to go throgh GetTEN 
    String  SpwSel              = "";     
    String  FieldSel            = "";
    String  TimeSel             = "";
    String  ScanSel             = "";
    String  FeedSel             = "";
    String  IntentSel           = "";
    String  ObservationSel      = "";
    String  UVRangeSel          = "";
    String  MSSelect            = "";

    String MsName = "listobs/uid___X02_X3d737_X1_01_small.ms";


protected:

        TestSelectData()
        {
        }

        ~TestSelectData()
        {
        }


        virtual void SetUp()
        {
            BaseClass::SetUp();
            printf ("TestMeasurementSet::SetUp:: called \n");
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
            printf ("TestMeasurementSet::TearDown:: called \n");
       }

        void test_selectdata(PointingDirectionCalculator & calc);

};


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




TEST_F(TestSelectData, Antenna )
{
    TestTitle( "selectData (key=Antenna)" );


    // MS name for this Test //
   
       String name = env.getCasaMasterPath() + MsName;
       printf( " Used MS is [%s] \n", name.c_str() );
  
    // Create Object //
    
       MeasurementSet ms( name.c_str() );
  
       PointingDirectionCalculator calc(ms);

    // Initial brief Inspection //
   
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n"); 
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE( (int)0, ExpectedNrow );
       
    //+
    // TEST SET 
    //  1) Specify the condition
    //  2) probe and confirm Exception
    //  3) Result check , expected Nrow was selected 
    //-

    int nrow ;

      AntSel = "";
        Title("Antenna: by NULL  (Matches) ",AntSel);

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, nrow);     // see MS in detail //

      AntSel = "hoge&&&";
        Title("Testing Abnormal Name = hoge (No Matches))",AntSel );

        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      AntSel = "DV01&&&";
        Title("Antenna: Normal specific Name. (Matches) ",AntSel);

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (180, (int)nrow);     // see MS in detail //

      AntSel = "DV02&&&";
        Title("Antenna: Normal specific Name (No Matches)",AntSel );

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (180, (int)nrow);

      AntSel = "PM03&&&";
        Title("Antenna: Normal specific Name (No Matches)",AntSel );

        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (180, (int)nrow);
#if 0
//+  
// THis hoge search may makes enexpected resut.
// Compare the first execution of HOGE search.
//-
      AntSel = "hoge&&&";
        Title("Testing Abnormal Name = hoge (No Matches))",AntSel );

        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow); 
#endif 
      AntSel = "DV*&&&";
        Title("Antenna: Normal with Wild card char. (Matches)",AntSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        printf("=> Calling getNrowForSelectedMS() after selectData() called.\n");
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ( 360, (int)nrow);

        //+
        // Other detail cases are programmed if needed. 
        //-
        
        
        return;
        
}


TEST_F(TestSelectData, Spw )
{
    TestTitle( "selectData (key=Spw)" );

    // MS name for this Test //
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );

    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      SpwSel = "";
        Title( "Spw:Nothig specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

    
      SpwSel = "*";
        Title( "Spw: Wildcard *  specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      SpwSel = "hoge";
        Title( "Spw: abnormal letters  specified.",SpwSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      SpwSel = "0:13~20";
        Title( "Spw: spw=0, ch=13~20 specified.",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ ( 270, (int)nrow);
#if 0
//
// Once execution OK, the next Error query makes unexpected Return ?
//
      SpwSel = "0:13=20";
        Title( "Spw: spw=0, ch=13~20 specified.(Syntax ERROR)",SpwSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

#endif 

      SpwSel = "1:13~15";
        Title( "Spw: spw=1, ch=13~20 specified.(None)",SpwSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (810, (int)nrow);


}

TEST_F(TestSelectData, Field )
{ 
    TestTitle( "selectData (key=Field)" );

    // Using MS //
    
        String MsName = "sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
      printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      FieldSel = "";
        Title( "Files:Nothig specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FieldSel = "*";
        Title( "Field: Wildcard *  specified.",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FieldSel = "hoge";
        Title( "Fieled: abnormal letters  specified.",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FieldSel = "0";  // Field ID 
        Title( "Field: ID=0 specified.(No exits)",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (0, (int)nrow);       // On-Table , None on MAIN .. //

      FieldSel = "1";  // Field ID 
        Title( "Field: ID=1 specified.(exits)",FieldSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FieldSel = "9";  // Out of Range
        Title( "Field: ID=9 Not on the FIELD TABLE",FieldSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

}



TEST_F(TestSelectData, Time )
{ 
    TestTitle( "selectData (key=Time)" );

    // Using MS //
    
        String MsName = "sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
      printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      TimeSel = "";
        Title( "Time: Nothig specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      TimeSel = "hoge";
        Title( "Time: abnormal letters  specified.",TimeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      TimeSel = ">1900/01/01/00:00:00";
        Title( "Time: since 1900/01/01/00:00:00  specified.",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      TimeSel = ">2018/01/01/00:00:00";
        Title( "Time: since 2018/01/01/00:00:00  specified. None matches",TimeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (0, (int)nrow);


      TimeSel = "<2014/12/4/00:39:25";  //  Time Condition  
        Title( "Field: a specifc Time, Limited Number  match are expected",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        printf( "# Actually detected Nrow = %d \n" , nrow); 
        EXPECT_GT (20, (int)nrow);      // On-Table  //
        EXPECT_EQ (5,   (int)nrow);    

      TimeSel = "2014/12/4/00:40:00~2014/12/4/00:40:10";  //  Combined  10 sec 
        Title( "Field: a specifc Time, Limited Number  match are expected",TimeSel);
        EXPECT_NO_THROW( test_selectdata(calc) ); // getTEN makes.
        nrow = calc.getNrowForSelectedMS();
        printf( "# Actually detected Nrow = %d \n" , nrow);
        EXPECT_GT (20, (int)nrow);      // On-Table  //
        EXPECT_EQ (17,   (int)nrow);    // SEE ACTUAL COUNT on Browser //




}


TEST_F(TestSelectData, Feed )
{ 
    TestTitle( "selectData (key=Feed" );

    // Using MS //
    
        String MsName = "/sdimaging/Uranus1.cal.Ant0.spw34.ms";    // One definition MEAS_FREQ_REF =5
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //

       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      FeedSel = "";
        Title( "Feed: Nothig specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FeedSel = "0";    // NO DATA /
        Title( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, (int)nrow);

      FeedSel = "1";    // NO DATA /
        Title( "Feed: ID specified.",FeedSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, (int)nrow);

}


TEST_F(TestSelectData, Intent )
{ 
    TestTitle( "selectData (key=Intent)" );

    // Using MS //
    
        String MsName = "sdimaging/selection_intent.ms";    // 
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //

       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE(0, ExpectedNrow );

      int nrow;
      IntentSel = "";
        Title( "Intent: Nothig specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);
 
       // Usage : In what way this term is used. 
        //   and where to exist in MS.
          
      IntentSel = "*HOGE*";     // 
        Title( "Scan: ID specified.",IntentSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, (int)nrow);

      IntentSel = "*CAL*";      // 
        Title( "Scan: ID specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (512, (int)nrow);

      IntentSel = "*BAND*";      // 
        Title( "Scan: ID specified.",IntentSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (512, (int)nrow);

}


TEST_F(TestSelectData, Observation )
{ 
    TestTitle( "selectData (key=Observation)" );

    // Using MS //
    
        String MsName = "/sdimaging/selection_spw.ms";    //    Three Observation entries.
    
    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      ObservationSel = "";
        Title( "Observation: Nothig specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

     ObservationSel = "hoge";     
        Title( "Observation: abnormal expr..",ObservationSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, (int)nrow);

      ObservationSel = "0";     // one DATA /
        Title( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (256, (int)nrow);

      ObservationSel = "1";    // second DATA _
        Title( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (256, (int)nrow);

      ObservationSel = "2";    // 3rd.Data 
        Title( "Observation: ID specified.",ObservationSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();


      ObservationSel = "9";    // No Data (err) 
        Title( "Observation: ID specified.",ObservationSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (0, (int)nrow);



}


TEST_F(TestSelectData, UVRange )
{ 
    TestTitle( "selectData (key=UV Range)" );

    // MS name for this Test //
    
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );
    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

      int nrow;
      UVRangeSel = "";
        Title( "UVrange: Nothig specified.",UVRangeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();
        EXPECT_EQ (ExpectedNrow, (int)nrow);

      UVRangeSel = "hoge";     
        Title( "UVrange: abnormal expr..",UVRangeSel);
        EXPECT_ANY_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (ExpectedNrow, (int)nrow);

      UVRangeSel = ">1.0lambda";        //  Exprecasacore::Table::TableOption:: Updatession Unknown...../
        Title( "UVrange: ID specified.",UVRangeSel);
        EXPECT_NO_THROW( test_selectdata(calc) );
        nrow = calc.getNrowForSelectedMS();

        EXPECT_EQ (540, (int)nrow);



}


struct FrameName {
  bool   available;
  String name;
};

class TestSetFrame : public BaseClass
{

public:

struct FrameName  MyFrameTypes[50]
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
        {false,  "HogeHogeHoge"},       // Bad name //
        {false,  "" }                   // terminator //

};


protected:

        TestSetFrame ()
        {
        }

        ~TestSetFrame ()
        {
        }


        virtual void SetUp()
        {
            BaseClass::SetUp();
            printf ("TestSetFrame::SetUp:: called \n");
        }

        virtual void TearDown()
       {
            BaseClass::TearDown();
            printf ("TestSetFrame::TearDown:: called \n");
       }
        
        void check_direction_info(PointingDirectionCalculator& calc, int n_frame );


};

void TestSetFrame::check_direction_info(PointingDirectionCalculator& calc, int n_frame )

{
    //+
    // Test setFrame (Str) 
    //   No Exception is expected all the time 
    //-

      EXPECT_NO_THROW( calc.setFrame( MyFrameTypes[n_frame].name ));


#if 0 
    //
    // reference frame 
    //
      casacore::MeasFrame refframe = calc.getReferenceFrame(); 

    // Get the epoch pointer (0 if not present)
         const Measure* epoch   = refframe .epoch();
   
    // Get the position pointer (0 if not present)
         const Measure* position        = refframe .position();
   
    // Get the direction pointer (0 if not present)
         const Measure* direction   = refframe .direction(); 
                          
    // Direction Column  Name //
         casacore::String  DirColName  = calc. getDirectionColumnName(); 
#endif 

   
    // Durection Type //
         casacore::MDirection DirType  = calc.getDirectionType();
 
          printf( "#   MDirection: [%s] \n",  DirType.toString().c_str()  ) ;
          printf( "#   Given String [%s] \n", MyFrameTypes[n_frame].name.c_str() );

        
        String converted = DirType.toString();
        String sub_str   = MyFrameTypes[n_frame].name;

        // GTEST :: Check SubString //

    if(  MyFrameTypes[n_frame].available == true)       
        EXPECT_TRUE( converted.find(sub_str) !=std::string::npos); 
    else
        EXPECT_FALSE( converted.find(sub_str) !=std::string::npos);


}

TEST_F(TestSetFrame, setFrame )
{ 
    TestTitle( "setFrame (String FrameName)" );

    // Using MS //
    
        String MsName = "listobs/uid___X02_X3d737_X1_01_small.ms";    
    
    // MS name for this Test //
        String name = env.getCasaMasterPath() + MsName;
        printf( " Used MS is [%s] \n", name.c_str() );

    
    // Create Object //
    
        MeasurementSet ms( name.c_str() );
    
        PointingDirectionCalculator calc(ms);
    
    // Initial brief Inspection //
    
       printf("=> Calling getNrowForSelectedMS() in Initial Inspection\n");
       ExpectedNrow = (int)calc.getNrowForSelectedMS();
       EXPECT_NE((int)0, ExpectedNrow );

    // Various Frame Type (String) //

    for( int i = 0; MyFrameTypes[i].name !="" ; i++  )
    {
        Title( "setFrame: Reserved Name.", MyFrameTypes[i].name.c_str() );

        // Execute and check Exception and other requirements//
        
        check_direction_info( calc, i ) ;
  


    }


}



//+
//   Environment Set  (UNDER CONSTRUCTION)
//
//   - copy file(s) from Master Directort
//-

void InitEnvironment()
{
    MyEnv env;
    String cmd = " " ;

    //+
    // Copy Master MS from standard MS data stsrted from CASAPATH
    //-

    printf( "cmd will be [%s]\n", cmd.c_str() );

//    system ( cmd.c_str() );

}


}  // END namespace



//--------------------------------------
//   Unit Test Main (google test)
//--------------------------------------



int main (int nArgs, char * args [])
 {

    casa::MSEdit msedit ;  ////  TENTATIME /////


    ::testing::InitGoogleTest(& nArgs, args);

   // Run Test //

    printf ("Going to RUN_TEST \n" );
    return (RUN_ALL_TESTS()) ;

}

