//# CTSelection.cc: Implementation of CTSelection class
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2003
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
//----------------------------------------------------------------------------

#include <synthesis/CalTables/CTSelection.h>
#include <ms/MSSel/MSSelectionTools.h>
#include <casa/Utilities/Sort.h>

namespace casa { //# NAMESPACE CASA - BEGIN

  CTSelection::CTSelection() {
    msSelection_p = new casacore::MSSelection();
  };

  CTSelection::CTSelection(
			const NewCalTable& ct,
			const casacore::MSSelection::MSSMode& mode,
			const casacore::String& timeExpr,
			const casacore::String& antennaExpr,
			const casacore::String& fieldExpr,
			const casacore::String& spwExpr,
			const casacore::String& taqlExpr,
			const casacore::String& scanExpr,
			const casacore::String& stateExpr,
			const casacore::String& observationExpr) : nct_p(ct) {
    // TBD: do antenna->taql!
    msSelection_p = new casacore::MSSelection(ct, mode, timeExpr, antennaExpr,
        fieldExpr, spwExpr, "", taqlExpr, "", scanExpr, "", stateExpr,
        observationExpr);
  }
  
  CTSelection::CTSelection (const CTSelection& other) {
    if (this != &other) {
        nct_p = other.nct_p;
        msSelection_p = other.msSelection_p;
        fullTEN_p = other.fullTEN_p;
    }
  }

  CTSelection& CTSelection::operator= (const CTSelection& other) {
    if (this != &other) {
        nct_p = other.nct_p;
        msSelection_p = other.msSelection_p;
        fullTEN_p = other.fullTEN_p;
    }
    
    return *this;
  }

  CTSelection::~CTSelection() {
    if (msSelection_p) {
        delete msSelection_p;
        msSelection_p = NULL;
    }
  }

  //---------------------------------------------------------------------------

  void CTSelection::reset(NewCalTable ct,
			  const casacore::MSSelection::MSSMode& mode,
			  const casacore::String& timeExpr,
			  const casacore::String& antennaExpr,
			  const casacore::String& fieldExpr,
			  const casacore::String& spwExpr,
			  const casacore::String& taqlExpr,
			  const casacore::String& scanExpr,
			  const casacore::String& stateExpr,
			  const casacore::String& observationExpr) {
      nct_p = ct;
      // TBD: do antenna->taql!
      msSelection_p->reset(ct, mode, timeExpr, antennaExpr, fieldExpr,
        spwExpr, "", taqlExpr, "", scanExpr, "", stateExpr, observationExpr);
  };

  //---------------------------------------------------------------------------

  casacore::TableExprNode CTSelection::toTableExprNode(
        casacore::MSSelectableTable* msLike) {
    // Convert the CT selection to a TableExprNode object, 
    // representing a TaQL selection in C++.
    // Input:
    //    msLike           const MSSelectableTable&  CalTable to bind TaQL
    // Output:
    //    toTableExprNode  TableExprNode             Table expression node
    //
    // Interpret all expressions and produce a consolidated TEN.  
    //
    casacore::String antExpr = msSelection_p->getExpr(
        casacore::MSSelection::ANTENNA_EXPR);
    if (!antExpr.empty()) doCalAntennaSel(antExpr, msLike); 
    return msSelection_p->toTableExprNode(msLike);
  }

  //---------------------------------------------------------------------------

  void CTSelection::doCalAntennaSel(const casacore::String& antennaExpr,
          casacore::MSSelectableTable* msLike) {
    // Override antenna selection for antenna-based cal tables
    casacore::MSSelectableTable::MSSDataType type = msLike->dataType();
    if ((type==casacore::MSSelectableTable::PURE_ANTENNA_BASED) ||
        (type==casacore::MSSelectableTable::REF_ANTENNA_BASED)) {
        // select ANTENNA1 only
        setTaqlAntennaSelection(antennaExpr, msLike);
        // also select baselines (given or constructed)
        setBaselineSelection(antennaExpr, msLike);
    }
  }

  void CTSelection::setTaqlAntennaSelection(const casacore::String& antennaExpr,
          casacore::MSSelectableTable* msLike) {
      // Convert antenna ids to taql ANTENNA1 selection
      // First use MSSelection to get selected antenna ids:
      // handles ranges (~), name->id conversion, negation, etc.
      casacore::MSSelection mssel;
      mssel.setAntennaExpr(antennaExpr); 
      casacore::TableExprNode ten = mssel.toTableExprNode(msLike);
      casacore::Vector<casacore::Int> selAntIds = mssel.getAntenna1List();

      // separate antenna selection and negation
      casacore::String antstr, notantstr;
      for (casacore::uInt i=0; i<selAntIds.size(); ++i) {
          casacore::Int antId = selAntIds(i);
          casacore::String antIdStr = casacore::String::toString(abs(antId));
          if (antId >= 0) {
              if (!antstr.empty()) antstr += ",";
              antstr += antIdStr;
          } else {
              if (!notantstr.empty()) notantstr += ",";
              notantstr += antIdStr;
          }
      }
      // make taqlExpr for ANTENNA1
      casacore::String taqlExpr;
      if (!antstr.empty()) 
          taqlExpr += "ANTENNA1 IN [" + antstr + "]";
      if (!notantstr.empty()) {
          if (!taqlExpr.empty()) taqlExpr += " && ";
          taqlExpr += "ANTENNA1 NOT IN [" + notantstr + "]";
      }
      // append taqlExpr to user's taql selection
      casacore::String msSelTaQL = msSelection_p->getExpr(
          casacore::MSSelection::TAQL_EXPR);
      if (!msSelTaQL.empty()) msSelTaQL += " && ";
      msSelTaQL += taqlExpr;
      //cout << "**** setTaQLExpr=" << msSelTaQL << endl;
      msSelection_p->setTaQLExpr(msSelTaQL);
  } 

  void CTSelection::setBaselineSelection(const casacore::String& antennaExpr,
          casacore::MSSelectableTable* msLike) {
        // Set baseline selections in MSSelection antennaExpr
        // First, separate baseline and antenna selections
        casacore::String baselineSel(""), antennaSel("");
        casacore::Vector<casacore::String> selections;
        selections = split(antennaExpr, ';', selections);
        for (casacore::uInt i=0; i<selections.size(); ++i) {
            if (selections(i).contains('&')) {
                if (!baselineSel.empty()) baselineSel += ";";
                baselineSel += selections(i);
            } else {
                if (!antennaSel.empty()) antennaSel += ";";
                antennaSel += selections(i);
            }
        }
        if (!antennaSel.empty()) {
            // Construct baselines (cross/auto-corr) for each antenna;
            // get ids for antenna selection only
            casacore::MSSelection mssel;
            mssel.setAntennaExpr(antennaSel);
            casacore::TableExprNode ten = mssel.toTableExprNode(msLike);
            casacore::Vector<casacore::Int> selAntIds = mssel.getAntenna1List();
            // get reference antennas
            casacore::Vector<casacore::Int> refAntIds = getRefAntIds(msLike);
            // construct baselines for antennas
            for (casacore::uInt i=0; i<selAntIds.size(); ++i) {
                casacore::Int antId = selAntIds(i);
                casacore::String sep = (baselineSel.empty() ? "" : ";");
                casacore::String neg = (antId < 0 ? "!" : "");
                casacore::String antstr = casacore::String::toString(abs(antId));
                casacore::String suffix = (isRefAntenna(antId, refAntIds) ?
                    "&&&" : "");
                baselineSel += sep + neg + antstr + suffix;
            }
        }
        // set new antennaExpr
        //cout << "**** setAntennaExpr=" << baselineSel << endl;
        msSelection_p->setAntennaExpr(baselineSel);
  }

  casacore::Vector<casacore::Int> CTSelection::getRefAntIds(
          casacore::MSSelectableTable* msLike) {
      // get unique antenna ids from antenna2 column
      NewCalTable nct = NewCalTable(*msLike->table());
      ROCTMainColumns ctmain(nct);
      casacore::Vector<casacore::Int> ant2 = ctmain.antenna2().getColumn();
      casacore::uInt nval = genSort(ant2, casacore::Sort::Descending, 
          (casacore::Sort::QuickSort | casacore::Sort::NoDuplicates));
      ant2.resize(nval, casacore::True);
      return ant2;
  }

  bool CTSelection::isRefAntenna(casacore::Int antennaId,
        casacore::Vector<casacore::Int> refantIds) {
      bool isrefant(false);
      for (casacore::uInt i=0; i<refantIds.size(); ++i) {
          if (abs(antennaId) == refantIds(i)) {
              isrefant = true;
              break;
          }
      }
      return isrefant;
  }

} //# NAMESPACE CASA - END

