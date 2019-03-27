///i/# PointingDirectionCalculator.cc: Implementation of PointingDirectionCalculator.h
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

//+
// Original Constructor
//- 
PointingDirectionCalculator::PointingDirectionCalculator(
        MeasurementSet const &ms) :
        originalMS_(new MeasurementSet(ms)), selectedMS_(), pointingTable_(), 
        pointingColumns_(), timeColumn_(), intervalColumn_(), antennaColumn_(), 
        directionColumnName_(), accessor_(NULL), antennaPosition_(), referenceEpoch_(),
        referenceFrame_(referenceEpoch_, antennaPosition_), 
        directionConvert_(NULL), directionType_(MDirection::J2000), movingSource_(NULL), 
        movingSourceConvert_(NULL), movingSourceCorrection_(NULL), 
        antennaBoundary_(), numAntennaBoundary_(0), pointingTimeUTC_(), lastTimeStamp_(-1.0),
        lastAntennaIndex_(-1), pointingTableIndexCache_(0), 
        shape_(PointingDirectionCalculator::COLUMN_MAJOR),
      /*CAS-8418*/ useSplineInterpolation_(true),	// Set when Spline is used. 
      /*CAS-8418*/ initializeReady_(5,false),           // Spline initialization Ready
      /*CAS-8418*/ coefficientReady_(5,false),          // Spline Coefficient Ready
      /*CAS-8418*/ accessorId_(DIRECTION)               // specify default accessor ID
{ 

// -- original code -- //
    accessor_ = directionAccessor;

    Block<String> sortColumns(2);
    sortColumns[0] = "ANTENNA1";
    sortColumns[1] = "TIME";
    selectedMS_ = new MeasurementSet(originalMS_->sort(sortColumns));

    debuglog << "Calling init()" << debugpost;
    init();

    // set default output direction reference frame
    debuglog << "Calling setFrame(J2000)" << debugpost;
    setFrame("J2000");

    // set default direction column name
    setDirectionColumn("DIRECTION");
}


PointingDirectionCalculator::~PointingDirectionCalculator()
{
  // Nothing at the moment. //
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

#if 0   // Old code //

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

#else

//+  
// New code reuired by CAS-8418  
//   - When setDirectionColumn is called , 
//   - new spline coeffcient table corresoinding to specified Direction column is generated. 
//      Once generated and from next time, lately created Obj. is pointed and used.
//   - from init(), this function is called.
//   - See indetal in;
//         initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );  
//-

    if (directionColumnName_ == "DIRECTION") {
        accessor_ = directionAccessor;
        accessorId_ =  DIRECTION;
        initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );

    } else if (directionColumnName_ == "TARGET") {
        accessor_ = targetAccessor;
        accessorId_ = TARGET;
        initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );
 
    } else if (directionColumnName_ == "POINTING_OFFSET") {
        accessor_ = pointingOffsetAccessor;
        accessorId_ = POINTING_OFFSET;
        initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );

    } else if (directionColumnName_ == "SOURCE_OFFSET") {
        accessor_ = sourceOffsetAccessor;
        accessorId_ = SOURCE_OFFSET;
        initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );

    } else if (directionColumnName_ == "ENCODER") {
        accessor_ = encoderAccessor;
        accessorId_ = ENCODER;
        initializeSplinefromPointingColumn(*originalMS_, accessorId_, true );

    } else {
        stringstream ss;
        ss << "Column \"" << columnNameUpcase << "\" is not supported.";
        throw AipsError(ss.str());
    }

    //+
    // Limmited service, force to use traditional Linear Interpolation. 
    //-

     if (getCurrentSplineObj()->isCoefficientReady() == false )
     {
         LogIO os(LogOrigin("PointingDirectionCalculator", "doGetDirection(i)", WHERE));
         os << LogIO::WARN << "INSUFFICIENT NUMBER OF POINTING DATA,  \n"
                           << "forced to Linear Interpolation " << LogIO::POST;
         useSplineInterpolation_ = false;
     }

#endif 
// ---org code ---
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
////////////////////////////////////////
// Reforming goGetdiretion [TENTATIVE]
// 
// - Switch to old doGetdirection() for tracable test
//   and to new doGetdirection() that includes Spline.
// - Accuracy(=Uncertainity) Test must be possible.
// - At the moment, this is trial.
//
///////////////////////////////////////

Vector<Double> PointingDirectionCalculator::doGetDirection(uInt irow)
{
    // In case Old source is needed, please locate Org.source and
    // make a link here.

    return (doGetDirectionNew(irow));
}

//----------------------------------
// CAS-8418 NEW doGetDirection(i)
//----------------------------------
Vector<Double> PointingDirectionCalculator::doGetDirectionNew(uInt irow) {
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
    debuglog << "Time " << setprecision(16) << currentTime << " idx=" << index
            << debugpost;
    //+
    // Check section on series of pointing data 
    // by 'Binary match' wth Pointing Table //
    // -
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
    } else {
        debuglog << "linear interpolation " << debugpost;

        //+
        // Copied from original (same logic)
        // - checcking Reference type.
        // - If needed , make conversion. 
        //-

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

        //+
        // CAS-8418::  Spline Interpolation
        //   using original var. see above for t0,t1,dt and nrowPointing.
        //-

        uInt antID = lastAntennaIndex_;   // WARNING: depends on resetAntennaPosition() // 
        Double dtime =  (currentTime - t0) ;

        // determin section on 'uIndex'
        //  please refer exact return code of binarySearch() 
            
        uInt uIndex;
        if( index >=1 )
        {  
            uIndex = index-1;
        }
        else if (index > (Int)(nrowPointing-1) )
        {
            uIndex = nrowPointing-1;
        }
        else // Out of Range //
        { 
            stringstream ss; ss  << "BUGCHECK, never come here. "  << endl;
            throw AipsError( ss.str() ); 
        }

        //+
        // Limittedly use Linear Interpolation.
        //  In any case when Spline is incapable of.
        //-
        Vector<Double> interpolated(2);
        if(useSplineInterpolation_)
        {
            // getCurrentSplineObj() referes, current Spline Obj. //
            interpolated = getCurrentSplineObj()-> calculate(uIndex, dtime, antID );
        }
        else // Linear (copied from original) //
        {
            Vector<Double> dirVal1 = dir1.getAngle("rad").getValue();
            Vector<Double> dirVal2 = dir2.getAngle("rad").getValue();
            Vector<Double> scanRate = dirVal2 - dirVal1;

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

//*********************************************************
//  CAS-8418::
//  Extended dunctions in PointingDirectionCalculator
//   for initializing Spline Interpolation
//*********************************************************

//+
// Direction Columns and 
// AccessorId and accessor_ (function pointer) 
//-

std::vector<string>   dirColList
= { "DIRECTION","TARGET","POINTING_OFFSET","SOURCE_OFFSET","ENCODER" };

std::vector<ACCESSOR> accList
= {  directionAccessor, targetAccessor,pointingOffsetAccessor,
   sourceOffsetAccessor, encoderAccessor };

// Column checck in Pointing Table //
bool PointingDirectionCalculator::checkColumn(MeasurementSet const &ms,String const &columnName )
{
    String columnNameUpcase = columnName;
    columnNameUpcase.upcase();
    if (true == (ms.pointing().tableDesc().isColumn(columnNameUpcase))) return true;
    else return false;
}
//+
// Activate Spline Interpolation
// - check specified Column if exist
// - create Spline object from specified Direction column
// - prepare Coeffient table for calulation. 
//-
bool PointingDirectionCalculator::initializeSplinefromPointingColumn(MeasurementSet const &ms,
                                                                    DirectionColumnID  DirColNo, 
                                                                    bool makeActive=true)
{
    debuglog << "initializeSplinefromPointingColumn, columNo=" << DirColNo << debugpost;

    String colName = dirColList[DirColNo] ;
    ACCESSOR acc   = accList[DirColNo] ;

    // Column Range check //
    if( DirColNo >= initializeReady_.size())
    {
        stringstream ss;
        ss << "Bugcheck. No column on Pointing Table." << endl;
        throw AipsError(ss.str());
        return false;  // Bad Param //
    }
    if( initializeReady_[DirColNo] == true )
    {
        debuglog << "initializeSplinefromPointingColumn, Normal,already active."  << debugpost;

        // Change Master pointer  // 
        if(makeActive) currSpline_ = splineObj_[DirColNo].get();

        return true;   // Servece already OK //
    }
    if(checkColumn(ms, colName  ))
    {
        debuglog << "Spline Obj:: attempt to construct by " << colName.c_str() << debugpost;

        // Temporary Obj. //
          unique_ptr<SplineInterpolation> spTemp( new SplineInterpolation(ms,acc));

        // Spline Available (N>4)
          coefficientReady_ [DirColNo] = spTemp-> isCoefficientReady();

        // move to Spline obj. //
          splineObj_[DirColNo] = std::move(spTemp);

        // copy to Master pointer, if this is desired.  // 
          if(makeActive) currSpline_ = splineObj_[DirColNo].get();

        // Obj. available //
           initializeReady_[DirColNo] = true;

        debuglog << "initializeSplinefromPointingColumn, Normal End." << debugpost;
        return true;
    } 
    else
    {
         stringstream ss;
         ss << "Bugcheck. No column on Pointing Table." << endl;
         throw AipsError(ss.str());
    }
}

//***************************************************
//  CAS-8418: 
//  Antenna Boundary (for Pointing Table ) methods
//  - create antenna information on Pointing Table.
//***************************************************

// Constructor //
AntennaBoundary::AntennaBoundary(MeasurementSet const &ms) 
{
        // Antenna Boundary body //
        antennaBoundary_.resize(ms.antenna().nrow() + 1);
        antennaBoundary_ = -1;

        Int count = 0;
        antennaBoundary_[count] = 0;
        ++count;

        // Pointing Table Handle and the Columns Handle//
        //   with sorting by AntennaID and Time 

        MSPointing hPointing_org  = ms.pointing();

        // Sort keys //
        Block<String> sortColumns(2);
        sortColumns[0]="ANTENNA_ID";
        sortColumns[1]="TIME";

        // Pointing Table handle //
        MSPointing hPointingTmp(hPointing_org.sort(sortColumns));
        hPointing_ = hPointingTmp;

        // Column Handle //
        std::unique_ptr<casacore::ROMSPointingColumns>
                columnPointing( new casacore::ROMSPointingColumns( hPointing_ ));

        // Antenna List //

        ROScalarColumn<casacore::Int>  antennaColumn = columnPointing->antennaId();
        Vector<Int> antennaList =  antennaColumn.getColumn();

        uInt nrow = antennaList.nelements();
        Int lastAnt = antennaList[0];

        for (uInt i = 0; i < nrow; ++i) 
        {    
            if (antennaList[i] > lastAnt) 
            {    
                antennaBoundary_[count] = i; 
                count++;
                lastAnt = antennaList[i];
            }    
            else if (antennaList[i] < lastAnt )
            { 
                 stringstream ss;
                 ss << "Bugcheck. Bad sort in creating antenna list." << endl;
                 throw AipsError(ss.str());
            }    
        }    

        antennaBoundary_[count] = nrow;
        ++count;
        numAntennaBoundary_ = count;

        if(showAntennaTable) showTable();  // DEBUG //

}
// getAntenaBoundary(start, end) //
std::pair<casacore::uInt, casacore::uInt> AntennaBoundary::getAntennaBoundary( casacore::uInt n )
{
    std::pair<casacore::uInt, casacore::uInt> pos(antennaBoundary_[n],antennaBoundary_[n+1]);
    return pos;
}

// progmmer's debug //
void AntennaBoundary::showTable() 
{
    debuglog << "AntennaBoundary::  numAntennaBoundary_ = " << numAntennaBoundary_  << debugpost;
    for (uInt b=0; b< antennaBoundary_.size(); b++)
    {
        printf( "new Created antennaBoundary_[%d] = %d \n", b, antennaBoundary_[b] );
    }
}

//***************************************************
//  Spline Inerpolation  methods
//  CAS-8418
//***************************************************

// constructor (for each accessor) //
SplineInterpolation::SplineInterpolation(MeasurementSet const &ms, ACCESSOR accessor ) 
{
    stsCofficientReady  = false;
  
    init(ms, accessor);
}

// initialize //
void SplineInterpolation::init(MeasurementSet const &ms, ACCESSOR const my_accessor)
{
    // Antenna Bounday //

        AntennaBoundary  antb(ms);
        uInt numAnt = antb.getNumOfAntenna();

    // prepere MS handle from selectedMS_
        MSPointing hPoint = antb.getPointingHandle();
        std::unique_ptr<casacore::ROMSPointingColumns>
                columnPointing( new casacore::ROMSPointingColumns( hPoint ));

    // Prepare Time and direction//

      Vector<Vector<Double> >          tmp_time;
      Vector<Vector<Vector<Double> > > tmp_dir;

    // Resize (top level) //
    
      tmp_time.        resize(numAnt);
      tmp_dir.         resize(numAnt);

      deltaTime_.          resize(numAnt);  
      splineInvalid_.      resize(numAnt);

    // Column handle (only time,direction are needed, others are reserved) //
    
      ROScalarColumn<Double> pointingTime           = columnPointing ->time();
      ROScalarColumn<Double> pointingInterval       = columnPointing ->interval();

    // Following columns are accessed by 'accessor_'   //
   
      ROArrayColumn<Double>  pointingDirection      = columnPointing ->direction();
      ROArrayColumn<Double>  pointingTarget         = columnPointing ->target();    
      ROArrayColumn<Double>  pointingPointingOffset = columnPointing ->pointingOffset();
      ROArrayColumn<Double>  pointingSourceOffset   = columnPointing ->sourceOffset();
      ROArrayColumn<Double>  pointingencoder        = columnPointing ->encoder();
 
    for(uInt ant=0; ant <numAnt; ant++)
    {
       // Antenna Bounday Pos(start,end) //
         std::pair<uInt,uInt> pos = antb.getAntennaBoundary(ant);
         uInt startPos = pos.first;
         uInt endPos   = pos.second;

        // define size of each antenna
          uInt size = endPos - startPos;  
          tmp_dir [ant]. resize(size);
          tmp_time[ant]. resize(size);

          // Discont time //
          deltaTime_     [ant]. resize(size);
          splineInvalid_ [ant]. resize(size);
          fill( splineInvalid_.begin(),
                splineInvalid_.end(),  0 );

          Double  prevTime = pointingTime.get(0);

        // for each row // 
        for (uInt row = startPos; row < endPos; row++) 
        {
            uInt index = row - startPos;

            // resizei (for Dir) //
            tmp_dir[ant][index].resize(2);

            Double        interval= pointingInterval.get(row);
            Double        time    = pointingTime.get(row);
            MDirection     dir    = my_accessor(*columnPointing, row);
            Vector<Double> dirVal = dir.getAngle("rad").getValue();

            // Diff Time //
            Double         dTime  = time - prevTime;
                           prevTime = time;

            // set on Vector //
            
            tmp_time[ant][index] = time;
            tmp_dir [ant][index] = dirVal;

            // Discont Time //
            deltaTime_[ant][index] = dTime;

            //+
            // NEW: CAS-8418: Jump discontiuous section.
            //  When dT is far greater than the interval, this is
            //  a discontinuous section. The currve on close to edges 
            //  may behave unexpected.
            //-

            if(false)
            {
                //+
                // Experiment:   set Invalid Flag
                //-

                if(dTime > (interval+0.001)  )
                {
                    printf( "Set MARK (%d, %d) around this point. | dt=%f, interval=%f \n",
                             ant,index, dTime, interval);
                    // Mark //
                    splineInvalid_[ant][index-4]= -4;
                    splineInvalid_[ant][index-3]= -3;
                    splineInvalid_[ant][index-2]= -2;
                    splineInvalid_[ant][index-1]= -1;
                    splineInvalid_[ant][index+0]= 1;
                    splineInvalid_[ant][index+1]= 2;
                    splineInvalid_[ant][index+2]= 3;
                    splineInvalid_[ant][index+3]= 4;
                }
            }

            // DBG //
            if(showSDPParam)
            {
                string  info = "-";
                if ( dTime > (interval + 0.001) ) info = "Jumped";

                printf("SDP ix=,%d a=,%d,tm=%f,dt=,%f, dir,%f,%f,%s\n", 
                       index, ant, time, dTime, dirVal[0],dirVal[1], info.c_str() );
            }
        }
    }

    //+
    // Minimum Condition Inspection
    // (N >=4) 
    //-

    for(uInt ant=0; ant <numAnt; ant++)
    {
        uInt st =tmp_time[ant].size();
        uInt s1 =tmp_dir [ant].size();
        if ((st < 4)||(s1 < 4))
        {
          // Warning .. //

            LogIO os(LogOrigin("SplineInterpolation", "init()", WHERE));
            os << LogIO::WARN << "INSUFFICIENT NUMBER OF POINTING DATA, must be ge. 4 \n" 
               << "Alternatively, Linear Interpolation will be used. " << LogIO::POST;

           stsCofficientReady = false; // initially in-usable ..
           return;
        }
    }
 
    //+
    // SDPosInterpolator Objct 
    //   - calulate Coefficient Table - 
    //-

      SDPosInterpolator  sdp (tmp_time, tmp_dir);
   
    // Obtain Coeff (copy object) //

      coeff_ = sdp.getSplineCoeff();

    // Table Active ..
      stsCofficientReady = true;

    // Programmers inspection 
      if(dumpCoeffientTable) dumpCsvCoeff();

}

// Dump Coefficient Table ..//
void SplineInterpolation::dumpCsvCoeff()
{
    FILE* fp = fopen( "coeff.csv","w" );

    for(uInt ant=0; ant < coeff_.size(); ant++ )
    {
        uInt size2 = coeff_[ant].size();
        printf( "AntID=%u, size=%u \n",ant, size2);
        for(uInt i=0; i< size2; i++)
        {
            Double x_c0 = coeff_[ant][i][0][0];
            Double x_c1 = coeff_[ant][i][0][1];
            Double x_c2 = coeff_[ant][i][0][2];
            Double x_c3 = coeff_[ant][i][0][3];

            Double y_c0 = coeff_[ant][i][1][0];
            Double y_c1 = coeff_[ant][i][1][1];
            Double y_c2 = coeff_[ant][i][1][2];
            Double y_c3 = coeff_[ant][i][1][3];

            fprintf(fp,"Spline::COEFF[%d],%4d,", ant, i );
            fprintf(fp, "X, %-12.5e, %-12.5e, %-12.5e, %-12.5e,|,",
                    x_c0, x_c1, x_c2, x_c3 );
                
            fprintf(fp, "Y, %-12.5e, %-12.5e, %-12.5e, %-12.5e \n",
                    y_c0, y_c1, y_c2, y_c3 );
        }
    }
    fclose(fp);

}

//+
// Interpolation Calculation
//-
casacore::Vector<casacore::Double> SplineInterpolation::calculate(uInt index,
                                                                  Double dt,
                                                                  uInt antID )
{
    debuglog << "SplineInterpolation::calculate()" << debugpost;

    // Error check //
    uInt arraySize = coeff_[antID].size();
    if(  index >= arraySize)
    {
        // Exception handling is to be here... //
        stringstream ss;
        ss << "Bugcheck. Requested Index is too large." << endl; 
        throw AipsError(ss.str());
    }

  
    // Coefficient //

    Double a0 = coeff_[antID][index][0][0];
    Double a1 = coeff_[antID][index][0][1];
    Double a2 = coeff_[antID][index][0][2];
    Double a3 = coeff_[antID][index][0][3];

    Double b0 = coeff_[antID][index][1][0];
    Double b1 = coeff_[antID][index][1][1];
    Double b2 = coeff_[antID][index][1][2];
    Double b3 = coeff_[antID][index][1][3];

   // Spline Calc //

    Double Xs =  (((0* dt + a3)*dt + a2)*dt + a1)*dt + a0;
    Double Ys =  (((0* dt + b3)*dt + b2)*dt + b1)*dt + b0;

    //+
    // CAS-8418 Discontinuous section 
    //-

    if((false) && ( splineInvalid_[antID][index] !=0) ) 
    {
#if 1
        printf("Warning:: Discontinuous Section was detected in calculation. (Ant=%d, index=%d),flag = %2d \n",
                antID, index, splineInvalid_[antID][index] ); 
        printf( "dt = %f \n",dt);
#endif 
        //+
        // Alter to LINEAR 
        //-

        Double a0p = coeff_[antID][index -1][0][0];
        Double b0p = coeff_[antID][index -1][1][0];
  
        if(true)
        {
            Xs = a0p;
            Ys = b0p;
        }
    }

    // Return //

    Vector<Double> outval(2); 
    outval[0] = Xs;
    outval[1] = Ys;

    debuglog << "SplineInterpolation::calculate() Normal return." << debugpost;

    return outval;

}
}  //# NAMESPACE CASA - END
