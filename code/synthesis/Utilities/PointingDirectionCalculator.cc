//# PointingDirectionCalculator.cc: Implementation of PointingDirectionCalculator.h
//# All helper functions of imager moved here for readability
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
// NEW //
#include <synthesis/Utilities/SDPosInterpolator.h>

using namespace casacore;
using namespace casacore;
using namespace std;

// Debug Message Handling
// if DIRECTIONCALC_DEBUG is defined, the macro debuglog and
// debugpost point standard output stream (std::cout and
// std::endl so that debug messages are sent to standard
// output. Otherwise, these macros basically does nothing.
// "Do nothing" behavior is implemented in NullLogger
// and its associating << operator below.
//
// Usage:
// Similar to standard output stream.
//
//   debuglog << "Any message" << any_value << debugpost;
//
  
//   #define DIRECTIONCALC_DEBUG

namespace {
struct NullLogger {
};

template<class T>
inline NullLogger &operator<<(NullLogger &logger, T /*value*/) {
    return logger;
}

#ifndef DIRECTIONCALC_DEBUG
NullLogger nulllogger;
#endif
}

#ifdef DIRECTIONCALC_DEBUG
#define debuglog cout << "PointingDirectionCalculator::DEBUG "
#define debugpost endl
#else
#define debuglog nulllogger
#define debugpost 0
#endif

namespace {
#define ARRAY_DIRECTION(ColumnName) \
inline MDirection ColumnName ## Accessor(ROMSPointingColumns &pointingColumns, uInt rownr) { \
    return pointingColumns.ColumnName ## Meas(rownr); \
}

#define SCALAR_DIRECTION(ColumnName) \
inline MDirection ColumnName ## Accessor(ROMSPointingColumns &pointingColumns, uInt rownr) { \
    return pointingColumns.ColumnName ## Meas()(rownr); \
}

ARRAY_DIRECTION(direction)
ARRAY_DIRECTION(target)
ARRAY_DIRECTION(pointingOffset)
ARRAY_DIRECTION(sourceOffset)
SCALAR_DIRECTION(encoder)

// working function for moving source correction
// convertToAzel must be configured with moving source direction and
// proper reference frame. Also, convertToCelestial must refer proper
// reference frame.
inline void performMovingSourceCorrection(
        CountedPtr<MDirection::Convert> &convertToAzel,
        CountedPtr<MDirection::Convert> &convertToCelestial,
        Vector<Double> &direction) {
    // moving source handling
    // If moving source is specified, output direction list is always
    // offset from reference position of moving source

        debuglog << "MovingSourceCorrection <Working>." << debugpost;

    // DEBUG (CAS-11818) //
    assert( convertToCelestial != nullptr );
    assert( convertToAzel != nullptr );
         
    MDirection srcAzel = (*convertToAzel)();
    MDirection srcDirection = (*convertToCelestial)(srcAzel);
    Vector<Double> srcDirectionVal = srcDirection.getAngle("rad").getValue();
    direction -= srcDirectionVal;
}

inline void skipMovingSourceCorrection(
        CountedPtr<MDirection::Convert> &/*convertToAzel*/,
        CountedPtr<MDirection::Convert> &/*convertToCelestial*/,
        Vector<Double> &/*direction*/) {

        debuglog << "MovingSourceCorrection <NO ACTION>" << debugpost;

    // do nothing
}
} // anonymous namespace

using namespace casacore;
namespace casa {
PointingDirectionCalculator::PointingDirectionCalculator(
        MeasurementSet const &ms) :
        originalMS_(new MeasurementSet(ms)), selectedMS_(), pointingTable_(), pointingColumns_(), timeColumn_(), intervalColumn_(), antennaColumn_(), directionColumnName_(), accessor_(
        NULL), antennaPosition_(), referenceEpoch_(), referenceFrame_(
                referenceEpoch_, antennaPosition_), directionConvert_(
        NULL), directionType_(MDirection::J2000), movingSource_(NULL), movingSourceConvert_(
        NULL), movingSourceCorrection_(NULL), antennaBoundary_(), numAntennaBoundary_(
                0), pointingTimeUTC_(), lastTimeStamp_(-1.0), lastAntennaIndex_(
                -1), pointingTableIndexCache_(0), shape_(
                PointingDirectionCalculator::COLUMN_MAJOR), 
/*CAS-8418*/    allAntennaBoundary_(0),allNumAntennaBoundary_(0)
{
    accessor_ = directionAccessor;

    Block<String> sortColumns(2);
    sortColumns[0] = "ANTENNA1";
    sortColumns[1] = "TIME";
    selectedMS_ = new MeasurementSet(originalMS_->sort(sortColumns));
//+
// CAS-8418: In ideal, Just here, Creating Interpolation should be performed
//   Making Table MUST BE done only ONE-TIME. 
//-
    initializeSplineInterpolation();

    init();

    // set default output direction reference frame
    setFrame("J2000");

    // set default direction column name
    setDirectionColumn("DIRECTION");
}
void PointingDirectionCalculator::init() {
    // attach column
    timeColumn_.attach(*selectedMS_, "TIME");
    intervalColumn_.attach(*selectedMS_, "INTERVAL");
    antennaColumn_.attach(*selectedMS_, "ANTENNA1");

    // initial setup
    debuglog << "inspectAntenna" << debugpost;
    inspectAntenna();
    debuglog << "done" << debugpost;

    resetAntennaPosition(antennaColumn_(0));
}

void PointingDirectionCalculator::selectData(String const &antenna,
        String const &spw, String const &field, String const &time,
        String const &scan, String const &feed, String const &intent,
        String const &observation, String const &uvrange,
        String const &msselect) {
    // table selection
    MSSelection thisSelection;
    thisSelection.setAntennaExpr(antenna);
    thisSelection.setSpwExpr(spw);
    thisSelection.setFieldExpr(field);
    thisSelection.setTimeExpr(time);
    thisSelection.setScanExpr(scan);
    thisSelection.setStateExpr(intent);
    thisSelection.setObservationExpr(observation);
    thisSelection.setUvDistExpr(uvrange);
    thisSelection.setTaQLExpr(msselect);

    TableExprNode exprNode = thisSelection.getTEN(&(*originalMS_));

    // sort by ANTENNA1 and TIME for performance reason
    Block<String> sortColumns(2);
    sortColumns[0] = "ANTENNA1";
    sortColumns[1] = "TIME";
    if (exprNode.isNull()) {
        debuglog << "NULL selection" << debugpost;
        selectedMS_ = new MeasurementSet(originalMS_->sort(sortColumns));
    } else {
        debuglog << "Sort of selection" << debugpost;
        MeasurementSet tmp = (*originalMS_)(exprNode);
        selectedMS_ = new MeasurementSet(tmp.sort(sortColumns));
    }
    debuglog << "selectedMS_->nrow() = " << selectedMS_->nrow() << debugpost;
    if (selectedMS_->nrow() == 0) {
        stringstream ss;
        ss << "Selected MS is empty for given selection: " << endl;
        if (!antenna.empty()) {
            ss << "\tantenna \"" << antenna << "\"" << endl;
        }
        if (!spw.empty()) {
            ss << "\tspw \"" << spw << "\"" << endl;
        }
        if (!field.empty()) {
            ss << "\tfield \"" << field << "\"" << endl;
        }
        if (!time.empty()) {
            ss << "\ttime \"" << time << "\"" << endl;
        }
        if (!scan.empty()) {
            ss << "\tscan \"" << scan << "\"" << endl;
        }
        if (!feed.empty()) {
            ss << "\tfeed \"" << feed << "\"" << endl;
        }
        if (!intent.empty()) {
            ss << "\tintent \"" << intent << "\"" << endl;
        }
        if (!observation.empty()) {
            ss << "\tobservation \"" << observation << "\"" << endl;
        }
        if (!uvrange.empty()) {
            ss << "\tuvrange \"" << uvrange << "\"" << endl;
        }
        if (!msselect.empty()) {
            ss << "\tmsselect \"" << msselect << "\"" << endl;
        }

        throw AipsError(ss.str());
    }

    init();

    debuglog << "done selectdata" << debugpost;
}

void PointingDirectionCalculator::configureMovingSourceCorrection() {
#if 0   // Causes Crash //
    if (!movingSource_.null() || directionColumnName_.contains("OFFSET")) {
#else   // FIXED //
    if ( !movingSource_.null() && !directionColumnName_.contains("OFFSET") )
    {
#endif
        debuglog << "configureMovingSourceCorrection::Perfrom." << debugpost;
        movingSourceCorrection_ = performMovingSourceCorrection;
    } else {
        debuglog << "configureMovingSourceCorrection::Skip." << debugpost;
        movingSourceCorrection_ = skipMovingSourceCorrection;
    }
}

void PointingDirectionCalculator::setDirectionColumn(String const &columnName) {
    String columnNameUpcase = columnName;
    columnNameUpcase.upcase();
    if (!(originalMS_->pointing().tableDesc().isColumn(columnNameUpcase))) {
        stringstream ss;
        ss << "Column \"" << columnNameUpcase
                << "\" doesn't exist in POINTING table.";
        throw AipsError(ss.str());
    }

    directionColumnName_ = columnNameUpcase;

    if (directionColumnName_ == "DIRECTION") {
        accessor_ = directionAccessor;
    } else if (directionColumnName_ == "TARGET") {
        accessor_ = targetAccessor;
    } else if (directionColumnName_ == "POINTING_OFFSET") {
        accessor_ = pointingOffsetAccessor;
    } else if (directionColumnName_ == "SOURCE_OFFSET") {
        accessor_ = sourceOffsetAccessor;
    } else if (directionColumnName_ == "ENCODER") {
        accessor_ = encoderAccessor;
    } else {
        stringstream ss;
        ss << "Column \"" << columnNameUpcase << "\" is not supported.";
        throw AipsError(ss.str());
    }

    configureMovingSourceCorrection();
}

void PointingDirectionCalculator::setFrame(String const frameType) {
    Bool status = MDirection::getType(directionType_, frameType);
    if (!status) {
        LogIO os(LogOrigin("PointingDirectionCalculator", "setFrame", WHERE));
        os << LogIO::WARN << "Conversion of frame string \"" << frameType
                << "\" into direction type enum failed. Use J2000."
                << LogIO::POST;
        directionType_ = MDirection::J2000;
    }

    // create conversion engine

    // Accessor 
    MDirection nominalInputMeasure = accessor_(*pointingColumns_, 0);

    // RefFrame
    MDirection::Ref outReference(directionType_, referenceFrame_);

    // Conversion 
    directionConvert_ = new MDirection::Convert(nominalInputMeasure,
            outReference);
    // Epoch 
    const MEpoch *e = dynamic_cast<const MEpoch *>(referenceFrame_.epoch());
    const MPosition *p =
            dynamic_cast<const MPosition *>(referenceFrame_.position());
    debuglog << "Conversion Setup: Epoch "
            << e->get("s").getValue() << " " << e->getRefString() << " Position "
            << p->get("m").getValue() << " " << p->getRefString()
            << debugpost;
}

void PointingDirectionCalculator::setDirectionListMatrixShape(
        PointingDirectionCalculator::MatrixShape const shape) {
    shape_ = shape;
}

void PointingDirectionCalculator::setMovingSource(String const sourceName) {
    MDirection sourceDirection(Quantity(0.0, "deg"), Quantity(90.0, "deg"));
    sourceDirection.setRefString(sourceName);
    setMovingSource(sourceDirection);
}

void PointingDirectionCalculator::setMovingSource(
        MDirection const &sourceDirection) {
    movingSource_ = dynamic_cast<MDirection *>(sourceDirection.clone());

    // create conversion engine for moving source
    MDirection::Ref refAzel(MDirection::AZEL, referenceFrame_);
    movingSourceConvert_ = new MDirection::Convert(*movingSource_, refAzel);

    configureMovingSourceCorrection();
}

void PointingDirectionCalculator::unsetMovingSource() {
  if (!movingSource_.null()) {
    movingSource_ = nullptr;
  }
}

//+
// NEW for SPLINE Interpolation
//-

void PointingDirectionCalculator::splineInit()
{

    //+
    // get Antenna count from  allAntennaBoundary 
    //-

        uInt numAnt = allNumAntennaBoundary_ -1;

    // prepere MS handle
    //  to access time and direction

        MSPointing hPointing  = selectedMS_->pointing();
        std::unique_ptr<casacore::ROMSPointingColumns>
                columnPointing( new casacore::ROMSPointingColumns( hPointing ));
        
    // Prepare Time and direction//

      Vector<Vector<Double> >          tmp_time;
      Vector<Vector<Vector<Double> > > tmp_dir;

    // Resize (top level) //
    
      tmp_time.        resize(numAnt);
      tmp_dir.         resize(numAnt);

    for(uInt ant=0; ant <numAnt; ant++)
    {
/*SN*/  printf("splineInit()::ant=%d\n",ant);

        uInt startPos = allAntennaBoundary_[ant];
        uInt endPos   = allAntennaBoundary_[ant+1];
        int size = endPos - startPos;

        // define size of each antenna  
          tmp_dir [ant]. resize(size);
          tmp_time[ant]. resize(size);

        // for each row // 
        for (uInt row = startPos; row < endPos; row++) 
        {
            uInt index = row - startPos;

            // resizei (for Dir) //
            tmp_dir[ant][index].resize(2);

            // values //
#if 0
            Double time           = pointingTimeUTC_[index];
            MDirection     dir    = accessor_(*pointingColumns_, index);
#else 
            Double        time    = (columnPointing->time()).get(index);
            MDirection     dir    = accessor_(*columnPointing, index);
#endif 
            Vector<Double> dirVal = dir.getAngle("rad").getValue();

            // set on Vector //
            tmp_time[ant][index] = time;
            tmp_dir [ant][index] = dirVal;

            if(false)
            {
                printf("SDP arg index=%d ant=%d, time=%f, dir=[%f,%f]\n", index, ant, time, dirVal[0],dirVal[1]);
            }
        }
    }

    //+
    // Minimum Condition Inspection
    // [TENTATIVE]
    //-

    for(uInt ant=0; ant <numAnt; ant++)
    {
        if(tmp_time[ant].size() <= 4) return;
        if(tmp_dir [ant].size() <= 4) return;
    }
    
    //+
    // SDPosInterpolator Objct 
    //   - calulate Coefficient Table - 
    //-

      SDPosInterpolator  sdp_ (tmp_time, tmp_dir);
   
    // Obtain Coeff , save//

      splineCoeff_ = sdp_.getSplineCoeff();

    // Dump //

    if(false)
    {
      FILE* fp = fopen( "coeff.csv","w" );
      for(uInt ant=0; ant < splineCoeff_.size(); ant++ )
      {
        uInt size2 = splineCoeff_[ant].size();

        for(uInt i=0; i< size2; i++)
        {
            Double x_c0 = splineCoeff_[ant][i][0][0];
            Double x_c1 = splineCoeff_[ant][i][0][1];
            Double x_c2 = splineCoeff_[ant][i][0][2];
            Double x_c3 = splineCoeff_[ant][i][0][3];

            Double y_c0 = splineCoeff_[ant][i][1][0];
            Double y_c1 = splineCoeff_[ant][i][1][1];
            Double y_c2 = splineCoeff_[ant][i][1][2];
            Double y_c3 = splineCoeff_[ant][i][1][3];

            if(true) // File (csv) 
            {   
                fprintf(fp,"Spline::COEFF[%d],%4d,", ant, i );
                fprintf(fp, "X, %-12.5e, %-12.5e, %-12.5e, %-12.5e,|,",
                        x_c0, x_c1, x_c2, x_c3 );
                
                fprintf(fp, "Y, %-12.5e, %-12.5e, %-12.5e, %-12.5e \n",
                      y_c0, y_c1, y_c2, y_c3 );
            }
        }
      }
      fclose(fp);
    }
 
}

//+
//  Coefficient Build API ( related with CAS-8418)
//  RESERVED
//- 

void  PointingDirectionCalculator::initializeSplineInterpolation()
{
     printf("iniializeSplineInterpolation() called. \n");

    if(true)  // ALWAYS TRUE, if false -> spline coeff is NOT READY./
    {
        //+
        // CAS-8418 (19-DEC-2018) 
        //     AntennaBounday for  Interporation was transferd from inspectAntenna() 
        //     No longer the original AntennaBounday_ ... are used.
        // 
        // Here, new AntenaBoundary is created, however directly accessing 
        // Pointing Table. Existing AntenaBoundary info is still untracable.
        // detail behavior is partl vague.   (20-Dec-2018)    
        //-

         printf("initializeSplineInterpolation()::START making  antennaList on Pointing. \n" ); 

        // (1) Do the same, but see PointingTable.antenna 

        allAntennaBoundary_.resize(selectedMS_->antenna().nrow() + 1);
        allAntennaBoundary_ = -1;

        Int countP = 0;
        allAntennaBoundary_[countP] = 0;
        ++countP;

        // Look at Pointing Table, here using original CASACORE call. //

        MSPointing hPointing  = selectedMS_->pointing();
        std::unique_ptr<casacore::ROMSPointingColumns>
                columnPointing( new casacore::ROMSPointingColumns( hPointing ));

        ROScalarColumn<casacore::Int>  antennaColumnP = columnPointing->antennaId();  
        Vector<Int> antennaListP =  antennaColumnP.getColumn();

        uInt nrowP = antennaListP.nelements();

        Int lastAntP = antennaListP[0];

        //+
        //(PROGRAM ROBUSTNESS)
        //
        // Antenna List Scan as done in original code //
        //  - This section must be Robust.
        //  - In case irregular table structure is detected ,exit and...
        //     disable SPLINE and optionally inform this event. 
        //- 
     
        for (uInt i = 0; i < nrowP; ++i) 
        {
            if(false) printf( "AntennaList[%d/%d] = %d \n",i, nrowP, antennaListP[i]  );

            if (antennaListP[i] > lastAntP) 
            {
                allAntennaBoundary_[countP] = i;
                ++countP;
                lastAntP = antennaListP[i];
            }
            else if (antennaListP[i] < lastAntP )
            {
                printf( "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n" );
                printf( "xx This antennaID alignment is not supported. xxxxxx \n" );
                printf( "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n" );
                fgSpline = false; // FORCE DISABLE//

                return; // ABORT and Ignore //
            }
        }

        allAntennaBoundary_[countP] = nrowP;
        ++countP;
        allNumAntennaBoundary_ = countP;


        // Show Boundary List .. //
        
        printf( "Created Antenna Boundary Info, allNumAntennaBoundary_  = %u \n", allNumAntennaBoundary_  );
        for (uInt b=0; b< allAntennaBoundary_.size(); b++)
        {
            printf( "Created antennaBoundary_[%d] = %d \n", b, allAntennaBoundary_[b] );
        }   
 
        // update flag //
        doneAntenaBoundaryCreate = true;
        printf( "inspectantenna():: END making AtennaBoundary and splineInit() ... \n");

        //+
        // Spline Make Table 
        //  - Only once !
        //-

        splineInit();   // In this module, all access to MS is performed.  //
    }

}

Vector<Double> PointingDirectionCalculator::splineCalulate(uInt index, Double dt,uInt antID )
{
    uInt arraySize = splineCoeff_[antID].size();

    Vector<Double> outval(2);  // Local work for return //

    // Coeffcient //

    Double a0=0.0;
    Double a1=0.0;
    Double a2=0.0;
    Double a3=0.0;

    Double b0=0.0;
    Double b1=0.0;
    Double b2=0.0;
    Double b3=0.0;

    if(  index < arraySize)
    {
        a0 = splineCoeff_[antID][index][0][0];
        a1 = splineCoeff_[antID][index][0][1];
        a2 = splineCoeff_[antID][index][0][2];
        a3 = splineCoeff_[antID][index][0][3];

        b0 = splineCoeff_[antID][index][1][0];
        b1 = splineCoeff_[antID][index][1][1];
        b2 = splineCoeff_[antID][index][1][2];
        b3 = splineCoeff_[antID][index][1][3];
    }

//+
// Spline
//-
    double Xs =  (((0* dt + a3)*dt + a2)*dt + a1)*dt + a0;
    double Ys =  (((0* dt + b3)*dt + b2)*dt + b1)*dt + b0;

// Return //

    // Spline interpolated//
      outval[0] = Xs;
      outval[1] = Ys;
 
    return outval;

}

Matrix<Double> PointingDirectionCalculator::getDirection() {
    assert(!selectedMS_.null());

    uInt const nrow = selectedMS_->nrow();
    debuglog << "selectedMS_->nrow() = " << nrow << debugpost;
    Vector<Double> outDirectionFlattened(2 * nrow);
    // column major data offset and increment for outDirectionFlattened,
    // and output matrix shape
    uInt offset = nrow;
    uInt increment = 1;
    // matrix shape: number of rows is nrow and number of columns is 2
    IPosition outShape(2, nrow, 2);
    if (shape_ == PointingDirectionCalculator::ROW_MAJOR) {
        // column major specific offset, increment and output shape
        offset = 1;
        increment = 2;
        // matrix shape: number of rows is 2 and number of columns is nrow
        outShape = IPosition(2, 2, nrow);
    }

    for (uInt i = 0; i < numAntennaBoundary_ - 1; ++i) {
        uInt start = antennaBoundary_[i];
        uInt end = antennaBoundary_[i + 1];
        uInt currentAntenna = antennaColumn_(start);

        printf( "getDirection::resetAntennaPosition(id=%d) calls. in %d Antennas.  \n",
                 currentAntenna,numAntennaBoundary_ );

        resetAntennaPosition(currentAntenna);

        debuglog << "antenna " << currentAntenna << " start " << start
                << " end " << end << debugpost;
        uInt const nrowPointing = pointingTimeUTC_.nelements();
        debuglog << "nrowPointing = " << nrowPointing << debugpost;
        debuglog << "pointingTimeUTC = " << min(pointingTimeUTC_) << "~"
        << max(pointingTimeUTC_) << debugpost;

        for (uInt j = start; j < end; ++j) {
            debuglog << "start index " << j << debugpost;

            // doGetDirection call //
            Vector<Double> direction = doGetDirection(j);

            debuglog << "index for lat: " << (j * increment)
                    << " (cf. outDirectionFlattened.nelements()="
                    << outDirectionFlattened.nelements() << ")" << debugpost;
            debuglog << "index for lon: " << (offset + j * increment)
                    << debugpost;
            outDirectionFlattened[j * increment] = direction[0];
            outDirectionFlattened[offset + j * increment] = direction[1];
        }
        debuglog << "done antenna " << currentAntenna << debugpost;
    }
    debuglog << "done getDirection" << debugpost;
    return Matrix < Double > (outShape, outDirectionFlattened.data());
}
Vector<Double> PointingDirectionCalculator::doGetDirection(uInt irow) {
    debuglog << "doGetDirection(" << irow << ")" << debugpost;
    Double currentTime =
            timeColumn_.convert(irow, MEpoch::UTC).get("s").getValue();
    resetTime(currentTime);

    // search and interpolate if necessary
    Bool exactMatch;
    uInt const nrowPointing = pointingTimeUTC_.nelements();
    // pointingTableIndexCache_ is not so effective in terms of performance
    // simple binary search may be enough,
    Int index = binarySearch(exactMatch, pointingTimeUTC_, currentTime,
            nrowPointing, 0);
    debuglog << "binarySearch result " << index << debugpost;
//    uInt n = nrowPointing - pointingTableIndexCache_;
//    Int lower = pointingTableIndexCache_;
//    debuglog << "do binarySearch n=" << n << " lower=" << lower
//            << " nrowPointing=" << nrowPointing << " cache="
//            << pointingTableIndexCache_ << debugpost;
//    Int index = binarySearch(exactMatch, pointingTimeUTC_, currentTime, n,
//            lower);
//    debuglog << "binarySearch result " << index << debugpost;
//    if (!exactMatch && lower > 0 && index == lower) {
//        // maybe out of range
//        n = nrowPointing - n;
//        lower = 0;
//        debuglog << "do second binarySearch n=" << n << " lower=" << lower
//                << " nrowPointing=" << nrowPointing << " cache="
//                << pointingTableIndexCache_ << debugpost;
//        index = binarySearch(exactMatch, pointingTimeUTC_, currentTime, n,
//                lower);
//        debuglog << "second binarySearch result " << index << debugpost;
//    }
//    pointingTableIndexCache_ = (uInt) max((Int) 0, (Int) (index - 1));
    debuglog << "Time " << setprecision(16) << currentTime << " idx=" << index
            << debugpost;
    MDirection direction;
    assert(accessor_ != NULL);
    if (exactMatch) {
        debuglog << "exact match" << debugpost;
        direction = accessor_(*pointingColumns_, index);
    } else if (index <= 0) {
        debuglog << "take 0th row" << debugpost;
        direction = accessor_(*pointingColumns_, 0);
    } else if (index > (Int) (nrowPointing - 1)) {
        debuglog << "take final row" << debugpost;
        direction = accessor_(*pointingColumns_, nrowPointing - 1);
//            } else if (currentInterval > pointingIntervalColumn(index)) {
//                // Sampling rate of pointing < data dump rate
//                // nearest interpolation
//                debuglog << "nearest interpolation" << debugpost;
//                Double dt1 = abs(pointingTimeUTC[index] - currentTime);
//                Double dt2 = abs(currentTime - pointingTimeUTC[index - 1]);
//                if (dt1 >= dt2) {
//                    // midpoint takes the value at index - 1
//                    direction = accessor_(*pointingColumns_, index - 1);
//                } else {
//                    direction = accessor_(*pointingColumns_, index);
//                }
    } else {
        debuglog << "linear interpolation " << debugpost;
        // Sampling rate of pointing > data dump rate (fast scan)
        // linear interpolation
        Double t0 = pointingTimeUTC_[index - 1];
        Double t1 = pointingTimeUTC_[index];
        Double dt = t1 - t0;
        debuglog << "Interpolate between " << setprecision(16) << index - 1
                << " (" << t0 << ") and " << index << " (" << t1 << ")"
                << debugpost;
        MDirection dir1 = accessor_(*pointingColumns_, index - 1);
        MDirection dir2 = accessor_(*pointingColumns_, index);
        String dirRef1 = dir1.getRefString();
        String dirRef2 = dir2.getRefString();
        MDirection::Types refType1, refType2;
        MDirection::getType(refType1, dirRef1);
        MDirection::getType(refType2, dirRef2);
        debuglog << "dirRef1 = " << dirRef1 << " ("
                << MDirection::showType(refType1) << ")" << debugpost;
        if (dirRef1 != dirRef2) {
            MeasFrame referenceFrameLocal((pointingColumns_->timeMeas())(index),
                    *(referenceFrame_.position()));
            dir2 = MDirection::Convert(dir2,
                    MDirection::Ref(refType1, referenceFrameLocal))();
        }

        // Get Original Direction //
        Vector<Double> dirVal1 = dir1.getAngle("rad").getValue();
        Vector<Double> dirVal2 = dir2.getAngle("rad").getValue();

        // SPLINE:: Preserve original calculation //
          Vector<Double> scanRate;
          Vector<Double> interpolated(2);

        if(fgSpline) // New CAS-8418 //
        { 
            //+
            // NEW Spline Interpolation
            //   using original var. see above for t0,t1,dt and nrowPointing.
            //-

            uInt antID = 0; // TENTATIVE //
            Double dtime =  (currentTime - t0) ;

            // determin section 
            //  please refer  exact retuen specification of binarySearch() 
            
            uInt uIndex;
            if( index >=1 )  uIndex = index-1;
            else if (index > (Int)(nrowPointing-1) )   uIndex = nrowPointing-1;
            else { printf( "BUGCHECK\n");  throw; } 
 
            Vector<Double> ttDir = splineCalulate(uIndex, dtime, antID );

            interpolated[0] = ttDir[0]; // 3rd. order Spline
            interpolated[1] = ttDir[1];

            if(false) {
                printf( "Nishie:: index=%d,  dtime=%f,", index, dtime );
                printf ("Dir=, %f, %f \n", ttDir[0],ttDir[1]);
            }
        }
        else // Original //
        {
            //+
            // Original Linear Interpolation
            //-

              scanRate = dirVal2 - dirVal1;
              interpolated = dirVal1 + scanRate * (currentTime - t0) / dt;
        }

        // Convert the interpolated diretion from MDirection to Vector //
          direction = MDirection(Quantum<Vector<Double> >(interpolated, "rad"),refType1);
        
    }
    debuglog << "direction = "
            << direction.getAngle("rad").getValue() << " (unit rad reference frame "
            << direction.getRefString()
            << ")" << debugpost;
    Vector<Double> outVal(2);
    if (direction.getRefString() == MDirection::showType(directionType_)) {
        outVal = direction.getAngle("rad").getValue();
    } else {
        MDirection converted = (*directionConvert_)(direction);
        outVal = converted.getAngle("rad").getValue();
        debuglog << "converted = " << outVal << "(unit rad reference frame "
                << converted.getRefString() << ")" << debugpost;
    }

    // moving source correction
    assert(movingSourceCorrection_ != NULL);
    movingSourceCorrection_(movingSourceConvert_, directionConvert_, outVal);

    return outVal;
}

Vector<Double> PointingDirectionCalculator::getDirection(uInt i) {
    if (i >= selectedMS_->nrow()) {
        stringstream ss;
        ss << "Out of range row index: " << i << " (nrow for selected MS "
                << getNrowForSelectedMS() << ")" << endl;
        throw AipsError(ss.str());
    }
    debuglog << "start row " << i << debugpost;
    Int currentAntennaIndex = antennaColumn_(i);
    debuglog << "currentAntennaIndex = " << currentAntennaIndex
            << " lastAntennaIndex_ = " << lastAntennaIndex_ << debugpost;
    Double currentTime =
            timeColumn_.convert(i, MEpoch::UTC).get("s").getValue();
    resetAntennaPosition(currentAntennaIndex);
    debuglog << "currentTime = " << currentTime << " lastTimeStamp_ = "
            << lastTimeStamp_ << debugpost;
    if (currentTime != lastTimeStamp_) {
        resetTime(i);
    }
    debuglog << "doGetDirection" << debugpost;
    Vector<Double> direction = doGetDirection(i);
    return direction;
}

Vector<uInt> PointingDirectionCalculator::getRowId() {
    return selectedMS_->rowNumbers();
}

Vector<uInt> PointingDirectionCalculator::getRowIdForOriginalMS() {
    return selectedMS_->rowNumbers(*originalMS_, True);
}

uInt PointingDirectionCalculator::getRowId(uInt i) {
    return selectedMS_->rowNumbers()[i];
}


void PointingDirectionCalculator::inspectAntenna() {
    // selectedMS_ must be sorted by ["ANTENNA1", "TIME"]
    antennaBoundary_.resize(selectedMS_->antenna().nrow() + 1);
    antennaBoundary_ = -1;
    Int count = 0;
    antennaBoundary_[count] = 0;
    ++count;

    Vector<Int> antennaList = antennaColumn_.getColumn();
    uInt nrow = antennaList.nelements();
    Int lastAnt = antennaList[0];

    for (uInt i = 0; i < nrow; ++i) {
        if (antennaList[i] != lastAnt) {
            antennaBoundary_[count] = i;
            ++count;
            lastAnt = antennaList[i];
        }
    }
    antennaBoundary_[count] = nrow;
    ++count;
    numAntennaBoundary_ = count;
    debuglog << "antennaBoundary_=" << antennaBoundary_ << debugpost;
    debuglog << "numAntennaBoundary_=" << numAntennaBoundary_ << debugpost;
}

void PointingDirectionCalculator::initPointingTable(Int const antennaId) {
    if (!pointingTable_.null() && !pointingColumns_.null()
            && pointingTable_->nrow() > 0
            && pointingColumns_->antennaId()(0) == antennaId) {
        // no need to update
        return;
    }
    debuglog << "update pointing table for antenna " << antennaId << debugpost;
    MSPointing original = selectedMS_->pointing();
    MSPointing selected = original(original.col("ANTENNA_ID") == antennaId);
    if (selected.nrow() == 0) {
        debuglog << "no rows for antenna " << antennaId << " try -1"
                << debugpost;
        // try ANTENNA_ID == -1
        selected = original(original.col("ANTENNA_ID") == -1);
 
#if 0   // follwing assert() invalidate throw AipsError() //
        assert(selected.nrow() > 0);
#endif 
        if (selected.nrow() == 0) {
            stringstream ss;
            ss << "Internal Error: POINTING table has no entry for antenna "
                    << antennaId << "." << endl;
            throw AipsError(ss.str());
        }
    }
    debuglog << "selected pointing rows " << selected.nrow() << debugpost;
    pointingTable_ = new MSPointing(selected.sort("TIME"));

    // attach columns
    pointingColumns_ = new ROMSPointingColumns(*pointingTable_);

    // initialize pointingTimeUTC_
    uInt const nrowPointing = pointingTable_->nrow();
    pointingTimeUTC_.resize(nrowPointing);
    ROScalarMeasColumn<MEpoch> pointingTimeColumn =
            pointingColumns_->timeMeas();
    for (uInt i = 0; i < nrowPointing; ++i) {
        MEpoch e = pointingTimeColumn(i);
        if (e.getRefString() == MEpoch::showType(MEpoch::UTC)) {
            pointingTimeUTC_[i] = e.get("s").getValue();
        } else {
            pointingTimeUTC_[i] =
                    MEpoch::Convert(e, MEpoch::UTC)().get("s").getValue();
        }
    }

    // reset index cache for pointing table
    pointingTableIndexCache_ = 0;

    debuglog << "done initPointingTable" << debugpost;
}

void PointingDirectionCalculator::resetAntennaPosition(Int const antennaId) {
    MSAntenna antennaTable = selectedMS_->antenna();
    uInt nrow = antennaTable.nrow();
    if (antennaId < 0 || (Int) nrow <= antennaId) {
        stringstream ss;
        ss << "Internal Error: Invalid ANTENNA_ID is specified (" << antennaId
                << ")." << endl;
        throw AipsError(ss.str());
    } else if (antennaId != lastAntennaIndex_ || lastAntennaIndex_ == -1) {
        ScalarMeasColumn < MPosition
                > antennaPositionColumn(antennaTable, "POSITION");
        antennaPosition_ = antennaPositionColumn(antennaId);
        debuglog << "antenna position: "
                << antennaPosition_.getRefString() << " "
                << setprecision(16) << antennaPosition_.get("m").getValue() << debugpost;
        referenceFrame_.resetPosition(antennaPosition_);

        printf("initPointingTalbe(id=%d) calls. \n",antennaId);
        initPointingTable(antennaId);

        lastAntennaIndex_ = antennaId;
    }
}

void PointingDirectionCalculator::resetTime(Double const timestamp) {
    debuglog << "resetTime(Double " << timestamp << ")" << debugpost;
    debuglog << "lastTimeStamp_ = " << lastTimeStamp_ << " timestamp = "
            << timestamp << debugpost;
    if (timestamp != lastTimeStamp_ || lastTimeStamp_ < 0.0) {
        referenceEpoch_ = MEpoch(Quantity(timestamp, "s"), MEpoch::UTC);
        referenceFrame_.resetEpoch(referenceEpoch_);

        lastTimeStamp_ = timestamp;
    }
}

}  //# NAMESPACE CASA - END
