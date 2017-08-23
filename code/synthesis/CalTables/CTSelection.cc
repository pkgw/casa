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
    // TBD: do antenna->taql!
    msSelection_p = new casacore::MSSelection(ct, mode, timeExpr, antennaExpr,
        fieldExpr, spwExpr, "", taqlExpr, "", scanExpr, "", stateExpr,
        observationExpr);
  }
  
  CTSelection::CTSelection (const CTSelection& other) {
    if (this != &other) {
        msSelection_p = other.msSelection_p;
        fullTEN_p = other.fullTEN_p;
    }
  }

  CTSelection& CTSelection::operator= (const CTSelection& other) {
    if (this != &other) {
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
      // TBD: do antenna->taql!
      msSelection_p->reset(msLike, mode, timeExpr, antennaExpr, fieldExpr,
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
    return msSelection_p->toTableExprNode(msLike);
  }

} //# NAMESPACE CASA - END

