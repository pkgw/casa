//# CTSelection.h: Class to represent a selection on a CASA CalTable
//# Copyright (C) 1996,1997,1998,1999,2001
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
//#
//# $Id$

#ifndef SYNTHESIS_CTSELECTION_H
#define SYNTHESIS_CTSELECTION_H

#include <casa/aips.h>
#include <ms/MSSel/MSSelection.h>
#include <ms/MSSel/MSSelectableTable.h>
#include <synthesis/CalTables/NewCalTable.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary> 
// CTSelection: Class to represent a selection on a CASA CalTable
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="" tests="" demos="">

// <prerequisite>
//   <li> <linkto class="NewCalTable">NewCalTable</linkto> module
//   <li> <linkto class="casacore::MSSelection">MSSelection</linkto> module
// </prerequisite>
//
// <etymology>
// From "CalTable" and "selection".
// </etymology>
//
// <synopsis>
// The CTSelection class represents a selection on a CASA CalTable (CT).
//
// The purpose of this class is to provides a simple expression based
// selection mechanism to both the end-user and developer wishing to
// perform query operations over a measurement set.  This class is a
// specialization of the CASACORE casacore::MSSelection class in order to
// override the antenna selection to select ANTENNA1 only instead of
// baselines, then using MSSelection as usual with the revised antennaExpr
// and taqlExpr.
//
// For a complete list of the STaQL interface refer to the
// casacore::MeasurementSet Selection Syntax document at: <a
// href="http://casa.nrao.edu/other_doc.shtml">Data Selection</a>
//
// The sub-expressions are interpreted in the order which they were
// set.  The order however in not important - any dependency on the
// order in which the expressions are evaluated is handled internally.
// The result of parsing the expressions is casacore::TableExprNode (TEN).
// All TENs from sub-expressions are finally ANDed and the resultant TEN
// is used to select the rows of the NewCalTable.
//
// </synopsis>
//
// <example>
// <srcblock>
// </srcblock>
// </example>
//
// <motivation>
// </motivation>
//
// <todo asof="Aug/14/2009">
// </todo>

  class CTSelection
  {
  public:

    // Default null constructor, and destructor
    CTSelection();
    virtual ~CTSelection();
    
    // Construct using a NewCalTable and the various selection expressions to
    // be applied to the given CT.  By default, the expressions will
    // be parsed immediately.  With mode=PARSE_LATE, the parsing will
    // be done with a call to toTableExprNode().
    CTSelection(const NewCalTable& ct,
		const casacore::MSSelection::MSSMode& mode= 
            casacore::MSSelection::PARSE_NOW,
		const casacore::String& timeExpr="",
		const casacore::String& antennaExpr="",
		const casacore::String& fieldExpr="",
		const casacore::String& spwExpr="",
		// const String& uvDistExpr="",         // not supported
		const casacore::String& taqlExpr="",
		//const String& polnExpr="",,           // not supported
		const casacore::String& scanExpr="",
		//const String& arrayExpr="",i          // not supported
		const casacore::String& stateExpr="",
		const casacore::String& observationExpr="");
    
    // Construct from a record representing a selection item at the
    // CLI or user interface level.  This is functionally same as the
    // constructor above with mode=PARSE_LATE.
    CTSelection(const casacore::Record& selectionItem);
    
    // Copy constructor
    CTSelection(const CTSelection& other);
    
    // Assignment operator
    CTSelection& operator=(const CTSelection& other);

    casacore::TableExprNode toTableExprNode(
            casacore::MSSelectableTable* msLike);

    inline casacore::TableExprNode getTEN() { return msSelection_p->getTEN(); }

    // clear selections
    inline void clear(const casacore::MSSelection::MSExprType type=
            casacore::MSSelection::NO_EXPR) 
        { msSelection_p->clear(type); }

    // Expression setters.  The following set*Expr() methods only set
    // the expressions.  Parsing is done with a call to
    // toTableExprNode().
    inline casacore::Bool setFieldExpr(const casacore::String& fieldExpr) {
        return msSelection_p->setFieldExpr(fieldExpr); }
    inline casacore::Bool setSpwExpr(const casacore::String& spwExpr) {
        return msSelection_p->setSpwExpr(spwExpr); }
    inline casacore::Bool setScanExpr(const casacore::String& scanExpr) {
        return msSelection_p->setScanExpr(scanExpr); }
    inline casacore::Bool setTimeExpr(const casacore::String& timeExpr) {
        return msSelection_p->setTimeExpr(timeExpr); }
    inline casacore::Bool setStateExpr(const casacore::String& stateExpr) {
        return msSelection_p->setStateExpr(stateExpr); }
    inline casacore::Bool setObservationExpr(
            const casacore::String& observationExpr) {
        return msSelection_p->setObservationExpr(observationExpr); }
    inline casacore::Bool setAntennaExpr(const casacore::String& antennaExpr) {
        return msSelection_p->setAntennaExpr(antennaExpr); }
    inline casacore::Bool setTaQLExpr(const casacore::String& taqlExpr) {
        return msSelection_p->setTaQLExpr(taqlExpr); }

    // Accessors for items selected:

    // Accessor for the list of antenna-1 selected.
    // Antennas affected by the baseline negation operator have the
    // antenna IDs multiplied by -1.
    // TBD: does this work with taql?
    inline casacore::Vector<casacore::Int> getAntenna1List(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getAntenna1List(ms); }
    
    // Accessor for the list of antenna-2 of the selected baselines.
    // Antennas affected by the baseline negation operator have the
    // antenna IDs multiplied by -1.
    inline casacore::Vector<casacore::Int> getAntenna2List(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getAntenna2List(ms); }
    
    // Accessor for the list of selected field IDs.
    inline casacore::Vector<casacore::Int> getFieldList(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getFieldList(ms); }

    // Accessor for the list of the specified time range(s) as the
    // start and end MJD values.  The time ranges are stored as columns.
    // Change 5/21/17: returns [startTime, stopTime, dT] CAS-10142
    inline casacore::Matrix<casacore::Double> getTimeList(
            const casacore::MeasurementSet* ms=NULL)
    { return msSelection_p->getTimeList(ms); }
    
    // Accessor for the list of the selected Spectral Window IDs.
    inline casacore::Vector<casacore::Int> getSpwList(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getSpwList(ms); }

    // Accessor for the list of the selected Observation IDs.
    inline casacore::Vector<casacore::Int> getObservationList(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getObservationList(ms); }

    // Accessor for the list of the selected Scan IDs.
    inline casacore::Vector<casacore::Int> getScanList(
            const casacore::MeasurementSet* ms=NULL) 
    { return msSelection_p->getScanList(ms); }

    // // casacore::Matrix<casacore::Int> getChanList(
    // //               const casacore::MeasurementSet* ms=NULL, 
    // // 			    const casacore::Int defaultStep=1,
    // // 			    const casacore::Bool sorted=false);
    // // { return msSelection_p->getChanList(ms, defaultStep, sorted); }

    // // //
    // // // Same as getChanList, except that the channels and steps are in Hz.
    // // //    
    // // casacore::Matrix<casacore::Double> getChanFreqList(
    // //                  const casacore::MeasurementSet* ms=NULL, 
    // // 				   const casacore::Bool sorted=false);
    // // { return msSelection_p->getChanFreqList(ms, sorted); }

    // This version of reset() works with generic MSSelectableTable
    // object.  Accessing the services of the CTSelection module via
    // this interface is recommended over the version of reset() that
    // uses MeasurementSet.
    void reset(casacore::MSSelectableTable& msLike,
	       const casacore::MSSelection::MSSMode& mode = 
               casacore::MSSelection::PARSE_NOW,
	       const casacore::String& timeExpr        = "",
	       const casacore::String& antennaExpr     = "",
	       const casacore::String& fieldExpr       = "",
	       const casacore::String& spwExpr         = "",
	       // const String& uvDistExpr             = "", // not supported
	       const casacore::String& taqlExpr        = "",
	       // const String& polnExpr               = "", // not supported
	       const casacore::String& scanExpr        = "",
	       // const String& arrayExpr              = "", // not supported
	       const casacore::String& stateExpr       = "",
	       const casacore::String& observationExpr = "");

private:

    // Resets msSelection_p expressions for antenna selection:
    // use taqlExpr to select on ANTENNA1 only for antenna selection,
    // use antennaExpr for baseline selection
    void doCalAntennaSel(const casacore::String& antennaExpr,
        casacore::MSSelectableTable* msLike);
    // append baseline selection and set taql selection
    void setAntennaSelections(casacore::String antsel,
          casacore::MSSelectableTable* msLike);

    // Reference antenna:
    // get reference antenna ids from cal table ANTENNA2 column
    casacore::Vector<casacore::Int> getRefAntIds(
        casacore::MSSelectableTable* msLike);
    // check if antennaId is a reference antenna
    bool isRefAntenna(casacore::Int antennaId, 
        casacore::Vector<casacore::Int> refantIds);
    // make baseline strings for ref ant with all ref ants
    casacore::String getRefAntBaselines(casacore::Int antId,
        casacore::Vector<casacore::Int> refantIds, casacore::String neg);

    // Antenna ID 0:
    // check if zero is selected (else it is negated but there is no -0)
    bool zeroIsSelected(casacore::String antennaExpr,
        casacore::MSSelectableTable* msLike);

    casacore::MSSelection* msSelection_p;
  };

} //# NAMESPACE CASA - END

#endif


