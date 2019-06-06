//# MWCAnimator.cc: Animator for MultiWCHolder class
//# Copyright (C) 2000,2001,2003
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

#include <casa/aips.h>
#include <casa/Containers/Record.h>
#include <casa/Exceptions/Error.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Utilities/DataType.h>
#include <casa/BasicSL/String.h>
#include <display/Display/MultiWCHolder.h>
#include <display/Display/Attribute.h>
#include <display/DisplayEvents/MWCAnimator.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// Constructor which makes an isolated MWCAnimator.
	MWCAnimator::MWCAnimator() {
	}

// Constructor for a single MultiWCHolder.
	MWCAnimator::MWCAnimator(MultiWCHolder &mholder) {
		addMWCHolder(mholder);
	}

// Destructor.
	MWCAnimator::~MWCAnimator() {
	}

// Add/remove MultiWCHolder/s.
	void MWCAnimator::addMWCHolder(MultiWCHolder &mholder) {
		if (isAlreadyRegistered(mholder)) {
			return;
		}
		itsMWCHList.push_back(&mholder);
	}

	void MWCAnimator::removeMWCHolder(MultiWCHolder &mholder) {
		if (!isAlreadyRegistered(mholder)) {
			return;
		}
		std::list<MultiWCHolder*> orig = itsMWCHList;
		itsMWCHList.clear( );
		std::copy_if( orig.begin( ), orig.end( ), std::back_inserter(itsMWCHList),
					  [&](MultiWCHolder *h){return h != &mholder;} );
	}

	void MWCAnimator::removeMWCHolders() { itsMWCHList.clear( ); }

// Set linear restrictions
	void MWCAnimator::setLinearRestrictions(AttributeBuffer &restrictions,
	                                        const AttributeBuffer &increments) {
        for ( auto h : itsMWCHList ) h->setLinearRestrictions(restrictions, increments);
	}

// Set linear restriction
	void MWCAnimator::setLinearRestriction(const RecordInterface &rec) {
		static String nameString("name"), valueString("value"),
		       incrementString("increment");
		if (rec.isDefined(nameString) && rec.isDefined(valueString) &&
		        rec.isDefined(incrementString) &&
		        (rec.dataType(valueString) == rec.dataType(incrementString))) {
			String name;
			rec.get(nameString, name);
			switch(rec.dataType(valueString)) {
			case TpInt: {
				Int val, inc;
				Int tol = 0;
				rec.get(valueString, val);
				rec.get(incrementString, inc);
				setLinearRestriction(name, val, inc, tol);
				break;
			}
			case TpFloat: {
				Float val, inc;
				Float tol = 1e-5;
				rec.get(valueString, val);
				rec.get(incrementString, inc);
				setLinearRestriction(name, val, inc, tol);
				break;
			}
			case TpDouble: {
				Double val, inc;
				Double tol = 1e-7;
				rec.get(valueString, val);
				rec.get(incrementString, inc);
				setLinearRestriction(name, val, inc, tol);
				break;
			}
			case TpRecord: {
				/*
				Record valrec, increc;
				valrec = rec.asRecord(valueString);
				increc = rec.asRecord(incrementString);
				QuantumHolder valh, inch;
				String error;
				if (valh.fromRecord(error, valrec)) {
				if (inch.fromRecord(error, increc)) {
				  setLinearRestriction(name, valh.asQuantum(), inch.asQuantum());
				}
				     }
				     */
			}
			default:
				// unsupported type at this point in time...should we say something?
				break;
			}
		} else {
			// failure - should we say something?
			return;
		}
	}


// Remove a restriction (including a 'linear' (ramped) one).
	void MWCAnimator::removeRestriction(const String& name) {
        for ( auto h : itsMWCHList ) h->removeRestriction(name);
	}


// Do we already have this WorldCanvasHolder registered?
	Bool MWCAnimator::isAlreadyRegistered(const MultiWCHolder &mholder) {
        return std::any_of( itsMWCHList.begin( ), itsMWCHList.end( ),
                            [&](MultiWCHolder *h){return h == &mholder;} );
	}

} //# NAMESPACE CASA - END

