//# PointingDirectionCalculator.h: Does for MSes various fixes which do not involve calibrating.
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
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
//#

#ifndef _SYNTHESIS_POINTING_DIRECTION_CALCULATOR_H_
#define _SYNTHESIS_POINTING_DIRECTION_CALCULATOR_H_

#include <casa/aipstype.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/MVDirection.h>
#include <casa/Utilities/CountedPtr.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSPointing.h>
#include <ms/MeasurementSets/MSPointingColumns.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasFrame.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/TableMeasures/ArrayMeasColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <synthesis/Utilities/SDPosInterpolator.h>

namespace casa {

///
// <summary>Calculates Pointing Direction Efficiently</summary>
//
//<use visibility=export>
//
// <reviewed reviewer="UNKNOWN" date="before2018/05/22" tests="" demos="">
// </reviewed>
//
// The PointingDirectionCalculator calculates pointing directions for selected MS rows.
// The calculation is based on the information recorded in MS POINTING table.
// In general, sampling intervals for spectral (autocorrelation) data (MS MAIN table)
// and pointing direction information (MS POINTING table) are independent. Thus,
// any interpolation in time domain should be necessary to calculate the pointing
// direction from the MS POINTING table. In addition, a direction reference frame
// in MS POINTING table may be different from the desired one. In that case,
// the frame conversion is also required. The class is intended to automate those
// complicated calculation. User can obtain a list of pointing directions by
// calling a few configuration method to select data and/or specify output direction
// reference frame. The class supports pointing direction calculation for moving
// objects. If the flag for moving source is turn on, returned direction will be
// shifted to the first timestamp.
//
// <prerequisite>
// </prerequisite>
//
// <etymology>
// </etymology>
//
// <synopsis>
// The PointingDirectionCalculator calculates pointing directions for selected
// MS rows. The calculation is based on the information recorded in MS POINTING
// table.
// </synopsis>
//
// <example>
// Here is a typical example for ALMA data. In ALMA, POINTING information is
// recorded in AZELGEO reference. However, there is a case that we need
// pointing direction in celestial coordinate such as J2000. In such case,
// the code will look like the following.
// <srcblock>
// MeasurementSet ms; // any MS
// PointingDirectionCalculator calculator(ms);
// calculator.setDirectionListMatrixShape(PointingDirectionCalculator::ROW_MAJOR);
// calculator.selectData(spw="19", antenna="PM01&&&"); // select spw 19 and antnena PM01
// calculator.setDirectionColumn("DIRECTION");
// calculator.setFrame("J2000");
// Matrix directions = calculator.getDirection();
// </srcblock>
// </example>
//
// <motivation>
// Pointing direction is one of the crucial information on single-dish data
// reduction. However, such information is not directly associated with
// spectral (autocorrelation) data. Instead, pointing directions must be
// calculated from the records in MS POINTING table. That calculation is
// sometimes complicated. It may contain an interpolation in time as well
// as a conversion to desired direction reference frame. This class is
// designed to automate such calculation. The implementation is somewhat
// specialized for the usecase of single-dish data reduction.
// </motivation>
//
//#! <templating arg=T>
//#! </templating>
//
// <thrown>
// </thrown>
//
// <todo asof="before2018/05/22">
//   <li> Define unit test
//   <li> Implement more sophisticated interpolation method such as cubic spline
// </todo>
///


//+
//  CAS-8418:
//    typedef of accessor_ and 
//    Direction column types.
//-
typedef 
casacore::MDirection (*ACCESSOR)(  casacore::ROMSPointingColumns &pointingColumns,
                                       casacore::uInt rownr);
typedef
enum DC_ { Undefined, DIRECTION, TARGET, POINTING_OFFSET, SOURCE_OFFSET, ENCODER } DirectionColumnID;


class SplineInterpolation;  // CAS-8418::Forward Refference //
class PointingDirectionCalculator {
public:
    // Enumerations for memory layout of the output pointing direction array.
    // User should select the layout according to an usercase of this class.
    enum MatrixShape {
        // Memory layout is "column major"
        COLUMN_MAJOR,
        // Memory layout is "row major"
        ROW_MAJOR
    };

    //+
    // Constructor
    //-
   
    //+
    //  CAS-8418 Extension.
    //   Create Spline five Objects on PointingDirectionCalculator.
    //
    //   - The setDirectionColumn(name) performs all the initialization and make Coefficient table.
    //   - In this constructor, init() and setDirectionColumn("DIRECTION") are called,
    //     soon or later, setDirectionColumn will be called, whenever new Direction column is used
    //     and interporation is made.
    //   - The spline object is created by a column, which contains all antenna_Id.
    //   - Cofficient table is obtained from SDPosInterpolation. 
    //   - SplineInterpolation class in this file provides calculation but it is sort of customize 
    //     for PoinitngDirectionCalculator.
    //   - To see the initialization, please start from setDirectionColumn().  
    //- 

    PointingDirectionCalculator(casacore::MeasurementSet const &ms);

    // Destructor
    ~PointingDirectionCalculator();

    // Select data in the given MS.
    // Each selection parameters accept MS data selection syntax.
    void selectData(casacore::String const &antenna = "", 
		casacore::String const &spw ="",
		casacore::String const &field = "", 
		casacore::String const &time = "",
            	casacore::String const &scan = "", 
		casacore::String const &feed = "",
            	casacore::String const &intent = "", 
		casacore::String const &observation = "", 
		casacore::String const &uvrange = "",
            	casacore::String const &msselect = "");

    //+
    // Select which POINTING column to use for pointing direction calculation.
    // Possible values are "DIRECTION" (default), "TARGET", "POINTING_OFFSET",
    // "SOURCE_OFFSET", and "ENCODER". These values are all case-sensitive.
    //-
    
    void setDirectionColumn(casacore::String const &columnName = "DIRECTION");

    //+
    // Set output direction reference frame. This accepts reference strings which
    // <linkto class=MDirection>MDirection</linkto> can recognize.
    // If given string is invalid, the frame will be set to "J2000".
    //-
 
    void setFrame(casacore::String const frameType);

    //+
    // Set output direction matrix shape. If "ROW_MAJOR" is given, the shape will be
    // (2, nrow) where nrow is a number of selected rows in MS. If "COLUMN_MAJOR" is
    // set, the shape will be (nrow, 2). User can choose appropriate shape according
    // to the access pattern for the output direction matrix.
    //-

    void setDirectionListMatrixShape(
            PointingDirectionCalculator::MatrixShape const shape);
    
    //+
    // <group>
    // Set source name for the moving source.
    // The method accepts source names which <linkto class=MDirection>MDirection</linkto>
    // can recognize (e.g. "Moon", "Jupiter").
    // If given string is invalid, exception will be thrown. User can specify the moving source
    // using a string or the MDirection instance.
    
    void setMovingSource(casacore::String const sourceName);
    void setMovingSource(casacore::MDirection const &sourceDirection);
    
    // </group>

    //+
    // Clear the moving source setting.
    //-
    
    void unsetMovingSource();

    //+
    // Return number of rows for selected MS.
    //-

    casacore::uInt getNrowForSelectedMS() {return selectedMS_->nrow();}

    //+
    // Return direction type as a <linkto class=MDirection>MDirection::Types</linkto>
    // enum.
    //-

    casacore::MDirection::Types const &getDirectionType() {return directionType_;}

    //+
    // Return an information on the moving source as a <linkto class=MDirection>MDirection</linkto>
    // instance.
    //-
    
    casacore::MDirection const &getMovingSourceDirection() {return *movingSource_;}

    //+
    // Return pointing direction matrix. Its shape depends on the user set the shape to
    // "ROW_MAJOR" or "COLUMN_MAJOR". Returned directions are interpolated to timestamps
    // recorded in the selected MS, and are converted to desired direction reference
    // frame if necessary.
    //-     

    casacore::Matrix<casacore::Double> getDirection();

    //+
    // Return pointing direction for specified row.
    // If irow is larger than or equal to the number of rows for selected MS,
    // exception will be thrown.
    //-

    casacore::Vector<casacore::Double> getDirection(casacore::uInt irow);

    // <group>
    // Return a list of row ids for selected rows. The getRowId will return the ids
    // in selected MS. On the other hand, the getRowIdForOriginalMS will return
    // the ids in original MS.

    casacore::Vector<casacore::uInt> getRowId();
    casacore::Vector<casacore::uInt> getRowIdForOriginalMS();

    // </group>

    //+
    // Return a row id for specified row.
    // If irow is larger than or equal to the number of rows for selected MS,
    // exception will be thrown.
    //-

    casacore::uInt getRowId(casacore::uInt irow);

    //*************************************
    // CAS-8418  Spline Interpolation API 
    //*************************************

    void setSplineInterpolation(bool mode) {useSplineInterpolation_ = mode;};

    // Spline Object handle
    casa::SplineInterpolation      *getCurrentSplineObj() { return currSpline_; }

    // Curret Direction column (=accessor in this source )
    casacore::uInt getCurretAccessorId()  { return  accessorId_ ; };

    // Spline device status (on Debug) //

    bool isInitializeDone (casacore::uInt n)   {return initializeReady_ [n]; }
    bool isCofficientReady(casacore::uInt n)   {return coefficientReady_[n]; }

private:

    void init();
    void initPointingTable(casacore::Int const antennaId);
    void resetAntennaPosition(casacore::Int const antennaId);
    void resetTime(casacore::Double const timestamp);
    void inspectAntenna();
    void configureMovingSourceCorrection();
    casacore::Vector<casacore::Double> doGetDirection(casacore::uInt irow);

    // table access stuff

    casacore::CountedPtr<casacore::MeasurementSet> 		originalMS_;
    casacore::CountedPtr<casacore::MeasurementSet> 		selectedMS_;
    casacore::CountedPtr<casacore::MSPointing> 			pointingTable_;
    casacore::CountedPtr<casacore::ROMSPointingColumns> 	pointingColumns_;
    casacore::ROScalarMeasColumn<casacore::MEpoch> 		timeColumn_;
    casacore::ROScalarColumn<casacore::Double> 			intervalColumn_;
    casacore::ROScalarColumn<casacore::Int> 			antennaColumn_;
    casacore::String 						directionColumnName_;

    casacore::MDirection (*accessor_)(	casacore::ROMSPointingColumns &pointingColumns,
            				casacore::uInt rownr);
    // conversion stuff

    casacore::MPosition 			antennaPosition_;
    casacore::MEpoch 				referenceEpoch_;
    casacore::MeasFrame 			referenceFrame_;
    casacore::CountedPtr<casacore::MDirection::Convert> directionConvert_;
    casacore::MDirection::Types 			directionType_;
    casacore::CountedPtr<casacore::MDirection> 		movingSource_;
    casacore::CountedPtr<casacore::MDirection::Convert> movingSourceConvert_;

    void (*movingSourceCorrection_)(
            casacore::CountedPtr<casacore::MDirection::Convert> &convertToAzel,
            casacore::CountedPtr<casacore::MDirection::Convert> &convertToCelestial,
            casacore::Vector<casacore::Double> &direction);

    // other
   
    casacore::Vector<casacore::uInt> 		antennaBoundary_;
    casacore::uInt 				numAntennaBoundary_;
    casacore::Vector<casacore::Double> 		pointingTimeUTC_;
    casacore::Double 				lastTimeStamp_;
    casacore::Int 				lastAntennaIndex_;
    casacore::uInt 				pointingTableIndexCache_;
    PointingDirectionCalculator::MatrixShape 	shape_;

    // privatize  default constructor

    PointingDirectionCalculator();

    //***************************
    // CAS-8418 Spline Extended
    //***************************

     // Spline object and flag (for each Direction Column)
      
        bool useSplineInterpolation_    = true;      // default: Use Spline if TRUE
        bool useOldInterpolationModule_ = false;     // default: Use origial doGetDirection 

        casa::SplineInterpolation    *currSpline_;      // Currently Active Object (traditiona Pointer)

     // Spline Object for each Direction-Column //

        // Hold the spline object. switched and  assigned to 'currSpline_.  pointer.//

        std::unique_ptr<casa::SplineInterpolation>     splineObj_[5];


        // Internal conditions to check limitted service.  
        casacore::Vector<bool>                         initializeReady_  ;
        casacore::Vector<bool>                         coefficientReady_ ;

     // Accessor ID (See typedef above. ) 

//      DirectionColumnID     accessorId_ ; 
        casacore::uInt        accessorId_ ;
  
     // creating temporary Spline object (in construction)

        bool checkColumn(casacore::MeasurementSet const &ms,
                         casacore::String const &columnName );

        bool activateSplinefromDirectionColumn (casacore::MeasurementSet const &ms, 
                                                casacore::String const &colname, 
                                                ACCESSOR acc);
        bool activateSplinefromDirectionColumn(casacore::MeasurementSet const &ms,
                                               casacore::uInt DirColNo, bool makeActive);

     //+
     // [new] doGetDirection(uint row)   CAS-8418
     // Changes:   
     //   - Performs Spline-Interpolation
     //   - In case number of dat is insfficient, alternatively uses Linea-Interpolation
     //   - Logic is mostly same as the old one, new calculation was inserted.
     // Limitatioon:
     //   - provided calcuation() requires deviding point (as written dt in src.) 
     //     between the given fixed points.
     // Migration:
     //   - Old module remains until the compatibility is validated.
     //   - see doGetDirection() in detail.  
     //-
     
     // The Old and the New // 
        casacore::Vector<casacore::Double> doGetDirectionOrg(casacore::uInt irow); // Org::Lenear only 
        casacore::Vector<casacore::Double> doGetDirectionNew(casacore::uInt irow); // New::Spline/inear
};

//+
// Interpolation Base class (reserved for future)
//-
class Interpolation {
public:

    // Interpolation TYpe
      enum InterpolationType {
        LINEAR,
        SPLINE
   };

private:
      casacore::uInt reservedInformaiton; 
};

//+
// CAS-8418
// Spline Interpolation Class
//-

class SplineInterpolation :public Interpolation  {

public:
        // Constructor 
          SplineInterpolation(casacore::MeasurementSet const &ms, ACCESSOR acc );

         ~SplineInterpolation() { };

        // Calculate function //
        casacore::Vector<casacore::Double>   calculate(casacore::uInt row,
                                                       casacore::Double dt,
                                                       casacore::uInt AntennaID =0);
        // Internal stauts (for inspection)  //

        bool isCofficientReady()   {return stsCofficientReady; }

private:
        //  default constructor 

         SplineInterpolation();

        //  constructor sub. 

          void init( casacore::MeasurementSet const &ms, ACCESSOR const my_accessor);

          casacore::Vector<casacore::Vector<casacore::Vector<casacore::Vector<casacore::Double> > > > 
          coeff_;
       
        // Interal Staus (one Spline status)//
 
          bool stsCofficientReady = false;
  
        // debug routune //
          void dumpCsvCoeff();

        // Programmers Debug Flag //
       
          bool verifySDP          = false;
          bool dumpCoeffientTable = false;

};

//+
// CAS-8418
// Antenna Boundary Class
//-
class AntennaBoundary {
public:
         AntennaBoundary(casacore::MeasurementSet const &ms) ;
        ~AntennaBoundary() { };

        std::pair<casacore::uInt, casacore::uInt>  getAntennaBoundary( casacore::uInt n );
    
        casacore::uInt  getNumOfAntenna() {return numAntennaBoundary_ - 1;}

        casacore::MSPointing  getPointingHandle() { return hPointing_; };

private:
       //  AntennaBoundary on Pointing Tablle 
         casacore::Vector<casacore::uInt>            antennaBoundary_;
         casacore::uInt                              numAntennaBoundary_;
     
       // Pointing Table handle
         casacore::MSPointing hPointing_; 
 
       // Programmer's Debug //
         void showTable() ;
         bool showAntennaTable =false;
};


} //# NAMESPACE CASA - END

#endif /* _SYNTHESIS_POINTINGi_DIRECTION_CALCULATOR_H_ */
