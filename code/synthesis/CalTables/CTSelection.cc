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
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/CalTables/CTInterface.h>
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
			const casacore::String& observationExpr) {
    msSelection_p = new casacore::MSSelection();
    setTimeExpr(timeExpr);
    setAntennaExpr(antennaExpr);
    setFieldExpr(fieldExpr);
    setSpwExpr(spwExpr);
    setTaQLExpr(taqlExpr);
    setScanExpr(scanExpr);
    setStateExpr(stateExpr);
    setObservationExpr(observationExpr);

    if (mode==casacore::MSSelection::PARSE_NOW) {
        CTInterface cti(ct);
        toTableExprNode(&cti);
    }
  }
  
  CTSelection::CTSelection (const CTSelection& other) {
    if (this != &other) {
        msSelection_p = other.msSelection_p;
    }
  }

  CTSelection& CTSelection::operator= (const CTSelection& other) {
    if (this != &other) {
        msSelection_p = other.msSelection_p;
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

  void CTSelection::reset(casacore::MSSelectableTable& msLike,
			  const casacore::MSSelection::MSSMode& mode,
			  const casacore::String& timeExpr,
			  const casacore::String& antennaExpr,
			  const casacore::String& fieldExpr,
			  const casacore::String& spwExpr,
			  const casacore::String& taqlExpr,
			  const casacore::String& scanExpr,
			  const casacore::String& stateExpr,
			  const casacore::String& observationExpr) {
    clear();
    setTimeExpr(timeExpr);
    setAntennaExpr(antennaExpr);
    setFieldExpr(fieldExpr);
    setSpwExpr(spwExpr);
    setTaQLExpr(taqlExpr);
    setScanExpr(scanExpr);
    setStateExpr(stateExpr);
    setObservationExpr(observationExpr);
    if (mode==casacore::MSSelection::PARSE_NOW)
        toTableExprNode(&msLike);
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
        // separate baseline and antenna selections
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
        msSelection_p->setAntennaExpr(baselineSel);
        // set antenna expr and taql expr for antennas
        setAntennaSelections(antennaSel, msLike);
    }
  }

  void CTSelection::setAntennaSelections(casacore::String antsel,
          casacore::MSSelectableTable* msLike) {
    casacore::String antIdstr, antstr(""), notantstr(""), baselines("");
    casacore::MSSelection mssel;
    casacore::TableExprNode ten;
    casacore::Vector<casacore::Int> ant1list;
    casacore::Int antId;

    // so the table is not left open CAS-11483
    const casacore::Table* tab = msLike->table();
    CTInterface cti(*tab);

    // get taql antenna1 selections (not negations) from baselines
    casacore::String antExpr = msSelection_p->getExpr(
        casacore::MSSelection::ANTENNA_EXPR);
    if (!antExpr.empty()) {
        mssel.setAntennaExpr(antExpr);
        ten = mssel.toTableExprNode(&cti);
        ant1list = mssel.getAntenna1List();
        for (casacore::uInt i=0; i<ant1list.size(); ++i) {
            antId = ant1list(i);
            if ((antId>0) || (antId==0 && zeroIsSelected(antExpr, msLike)))
                antstr += (antstr.empty() ? "" : ",") +
                    casacore::String::toString(antId);
        }
    }

    // For each antenna:
    // 1. taql selection: IN/NOT IN ANTENNA1 only
    // 2. baseline selection: auto-corr for ref ant, else cross
    if (!antsel.empty()) {
        mssel.clear();
        ant1list.resize();
        mssel.setAntennaExpr(antsel);
        ten = mssel.toTableExprNode(&cti);
        ant1list = mssel.getAntenna1List();
        // Get antenna2 ids (reference antennas) for suffix
        casacore::Vector<casacore::Int> refAntIds = getRefAntIds(msLike);
        casacore::String sep, neg;
        for (casacore::uInt i=0; i<ant1list.size(); ++i) {
          antId = ant1list(i);
          antIdstr = casacore::String::toString(abs(antId));
          // separate antenna selection and negation for taql
          if ((antId>0) || (antId==0 && zeroIsSelected(antsel,msLike))) {
            antstr += (antstr.empty() ? "" : ",") + antIdstr;
            neg = "";
          } else {
            notantstr += (notantstr.empty() ? "" : ",") + antIdstr;
            neg = "!";
          }
          // make baseline string
          sep = (baselines.empty() ? "" : ";");
          if (isRefAntenna(antId, refAntIds)) {
            baselines += sep + getRefAntBaselines(antId, refAntIds, neg);
          } else {
            baselines += sep + neg + antIdstr;
          }
        }
        // antenna selection first, then baselines
        sep = (antExpr.empty() ? "" : ";");
        antExpr = baselines + sep + antExpr;
        msSelection_p->setAntennaExpr(antExpr);
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
    msSelection_p->setTaQLExpr(msSelTaQL);
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

  casacore::String CTSelection::getRefAntBaselines(casacore::Int antId,
        casacore::Vector<casacore::Int> refantIds, casacore::String neg) {
      // make baseline strings for antId with all refantIds
      casacore::String baselineStr(""), sep, antIdStr;
      for (casacore::uInt i=0; i<refantIds.size(); ++ i) {
            sep = (baselineStr.empty() ? "" : ";");
            antIdStr = casacore::String::toString(abs(antId));
            if (abs(antId)==refantIds(i))  // auto-correlation
                baselineStr += sep + neg + antIdStr + "&&&";
            else if (refantIds(i) != -1) {  // cross-correlation
                baselineStr += sep + neg + antIdStr + "&" + 
                    casacore::String::toString(refantIds(i));
            }
      }
      return baselineStr;
  }

  bool CTSelection::zeroIsSelected(casacore::String antennaExpr,
        casacore::MSSelectableTable* msLike) {
    // check if id "0" or antenna0 name is in expression
    NewCalTable nct = NewCalTable(*msLike->table());
    ROCTColumns ctcol(nct);
    casacore::String ant0name = ctcol.antenna().name().getColumn()(0);
    casacore::Vector<casacore::String> selections, antennas;
    selections = split(antennaExpr, ';', selections);
    for (casacore::uInt i=0; i<selections.size(); ++i) {
        if (selections(i).contains('&')) {
          // for baseline we only care about antenna1
          if (selections(i).startsWith("0") || 
              selections(i).startsWith(ant0name))
            return true;
        } else {
            if (!selections(i).startsWith("!")) {
                // split possible antenna list
                antennas = split(selections(i), ',', antennas);
                for (casacore::uInt j=0; j<antennas.size(); ++j) {
                    if (antennas(j).startsWith("0") || 
                        antennas(j).startsWith(ant0name))
                      return true;
                }
            }
        }
    }
    return false;
  }

} //# NAMESPACE CASA - END

