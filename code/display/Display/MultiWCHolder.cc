//# MultiWCHolder.cc: Holder of multiple WorldCanvasHolders for panelling
//# Copyright (C) 2000,2001,2002,2003
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
#include <display/Display/WorldCanvasHolder.h>
#include <display/Display/WorldCanvas.h>
#include <display/DisplayDatas/DisplayData.h>
#include <display/DisplayDatas/PrincipalAxesDD.h>
#include <display/Display/AttributeBuffer.h>
#include <display/Display/AttValTol.h>
#include <display/Display/MultiWCHolder.h>
#include <casa/BasicMath/Math.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// Default constructor.
	MultiWCHolder::MultiWCHolder() :
		itsBLength(0), itsBIndex(0),
		itsHoldCount(0),
		itsRefreshHeld(false) {
		setBIndexName();
	}

// Constructor for a single WorldCanvasHolder.
	MultiWCHolder::MultiWCHolder(WorldCanvasHolder &holder) :
		itsHoldCount(0),
		itsRefreshHeld(false) {
		setBIndexName();
		addWCHolder(holder);
	}

// Destructor.
	MultiWCHolder::~MultiWCHolder() {
	}

// Add/remove WorldCanvasHolder/s.
	void MultiWCHolder::addWCHolder(WorldCanvasHolder &holder) {
		if (isAlreadyRegistered(holder)) {
			return;
		}

        itsWCHList.push_back(&holder);

		for (Int i = 0; i < itsHoldCount; i++) {
			holder.worldCanvas()->hold();
		}

		installRestrictions(holder);
		addAllDisplayDatas(holder);
	}
	void MultiWCHolder::removeWCHolder(WorldCanvasHolder &holder) {
		if (!isAlreadyRegistered(holder)) {
			return;
		}
		removeAllDisplayDatas(holder,true);
		std::list<WorldCanvasHolder*> orig = itsWCHList;
		itsWCHList.clear( );
		std::copy_if( orig.begin( ), orig.end( ), std::back_inserter( itsWCHList ),
					  [&](WorldCanvasHolder *h) { return h != &holder; } );
	}
	void MultiWCHolder::removeWCHolders() {
		itsWCHList.clear( );
	}

// Add/remove DisplayData/s.
	void MultiWCHolder::addDisplayData(DisplayData &displaydata, int position) {
		if (isAlreadyRegistered(displaydata)) {
			return;
		}
		hold();

		if ( position < 0 || position >= static_cast<int>(itsDDList.size( )) ) {
			itsDDList.push_back(&displaydata);
		} else {
			auto dd = itsDDList.begin( );
			std::advance( dd, position );
			itsDDList.insert( dd, &displaydata );
		}

		addToAllWorldCanvasHolders(displaydata, position);

		// Add a 'bIndex' restriction to newly-added DD.  It can be used to
		// alternate display of the various DDs by placing a similar restriction
		// on the WCHs.	 The index should reflect its order in the list.  However,
		// contours, vectors, etc are not displayed separately so the index needs
		//to take that into account.
		if ( isBlinkDD(&displaydata) ){

			itsBLength++;
			itsBlinkDDs.resize( itsBLength, true );
		}

		int index = 0;
		for ( auto dd : itsDDList ) {

			if ( dd ){
				if ( isBlinkDD( dd) ) {
					itsBlinkDDs[index] = dd;

				}
				Attribute bIndexAtt(itsBIndexName, index );
				dd->setRestriction( bIndexAtt );
				if ( isBlinkDD(dd)){
					index++;
				}
			}

		}

		refresh();
		release();
	}


	void MultiWCHolder::removeDisplayData(DisplayData &displaydata) {

		if (!isAlreadyRegistered(displaydata)) return;

		hold();
		removeFromAllWorldCanvasHolders(displaydata);
		std::list<DisplayData*> orig = itsDDList;
		std::list<DisplayData*> removed;
		itsDDList.clear( );

		std::partition_copy( orig.begin( ), orig.end( ),
							 std::back_inserter( removed ),
							 std::back_inserter( itsDDList ),
							 [&](DisplayData *dd) { return dd == &displaydata; } );

		// No point in leaving blink restriction hanging on the dd.
		if ( removed.size( ) > 0 && isBlinkDD(&displaydata) ) {

			bool found=false;

			for( Int ddBIndex=0; ddBIndex<itsBLength; ddBIndex++ ) {
				DisplayData* searchDD = static_cast<DisplayData*>(itsBlinkDDs[ddBIndex]);

				if( searchDD == &displaydata ) {
					// dd found in blinkDD list--it will be removed.
					found=true;

					if ( itsBIndex > ddBIndex ){

						itsBIndex--;
					}
				} else if( found ) {

					// DDs past the one being removed move back in the blinkDD list.
					// Their bIndex restriction must also be decremented.


					Int newddBIndex=ddBIndex-1;

					itsBlinkDDs[newddBIndex]=searchDD;
					Attribute bIndexAtt(itsBIndexName, newddBIndex);
					searchDD->setRestriction(bIndexAtt);

				}

			}
			// itsBIndex is communicated to the animator, and becomes the
			// WCH blink restriction setting.  It should be decremented
			// if it was selecting a DD past the one deleted, in order
			// to continue selecting the same DD.

			if ( found ) {		// (should be true).
				itsBLength--;
				if ( itsBLength >= 0 ){
					itsBlinkDDs[itsBLength] = NULL;
				}
				itsBIndex = max(0, min(itsBLength-1, itsBIndex));
				// Assure itsBIndex is in proper range

			}

		}

		refresh();
		release();
	}


	void MultiWCHolder::removeDisplayDatas() {
		hold();

		for ( auto dd : itsDDList ) {
			removeFromAllWorldCanvasHolders(*dd);
			if(isBlinkDD(dd)) dd->removeRestriction(itsBIndexName);
		}
		itsDDList.clear( );

		itsBLength = itsBIndex = 0;
		refresh();
		release();
	}


// Install/remove restriction/s.
	void MultiWCHolder::setRestriction(const Attribute &restriction) {
		itsAttributes.set(restriction);
		distributeRestrictions();
	}
	void MultiWCHolder::setRestrictions(const AttributeBuffer &restrictions) {
		itsAttributes.set(restrictions);
		distributeRestrictions();
	}
	void MultiWCHolder::removeRestriction(const String &name) {
		String nm = (name=="bIndex")?  itsBIndexName : name;
		itsAttributes.remove(nm);
        for ( auto holder : itsWCHList ) holder->removeRestriction(nm);
	}

	void MultiWCHolder::removeRestrictions() {
		itsAttributes.clear();
		distributeRestrictions();
		// dk note: line above accomplishes nothing; restrictions remain on
		// WCHs at present, i.e. this routine doesn't work.
		// (Implementation was never finished; to be fixed).
	}

// Distribute restrictions linearly.
	void MultiWCHolder::setLinearRestrictions(AttributeBuffer &restrictions,
			const AttributeBuffer &increments) {

		AttributeBuffer rstrs=restrictions;
		adjustBIndexName(rstrs);
		AttributeBuffer incrs=increments;
		adjustBIndexName(incrs);
		// Same buffers, except with modified name of 'bIndex' attribute.

		Int bInd = 0;
		Bool BIExists = ( itsBLength>0 &&
						  rstrs.getValue(itsBIndexName, bInd) &&
						  bInd>=0 );
		// There are blink DDs to control, and a bIndex
		// restriction (with a reasonable value) exists.

		if(BIExists) itsBIndex=bInd;
		// Maintain internal record of animator bIndex setting.
		// When DDs are removed, its appropriate value may change,
		// and is communicated back to the animator.

		for ( auto holder : itsWCHList ) {
			holder->setRestrictions(rstrs);

			restrictions += increments;
			// to retain (dubious) semantics of 'restrictions' return value...
			rstrs += incrs;

			if(BIExists) {

				// Do a modulo-length adjustment to blink index, so that there
				// are no empty panels.	 (In my opinion, this should be done for
				// zIndex as well.	(dk)).

				rstrs.getValue(itsBIndexName, bInd);
				if(bInd<0 || bInd>=itsBLength) {
					bInd = max(0,bInd) % itsBLength;
					restrictions.set("bIndex", bInd);
					rstrs.set(itsBIndexName, bInd);
				}
			}
		}
		refresh();
	}

	void MultiWCHolder::hold() {
		itsHoldCount++;
		for ( auto holder : itsWCHList )
			holder->worldCanvas()->hold();
	}

	void MultiWCHolder::release() {
		itsHoldCount--;
		if (itsHoldCount <= 0) {
			itsHoldCount = 0;
			if (itsRefreshHeld) {
				refresh(itsHeldReason);
			}
			itsRefreshHeld = false;
		}
		for ( auto holder : itsWCHList )
			holder->worldCanvas()->release();
	}

	void MultiWCHolder::refresh(const Display::RefreshReason &reason) {
		if (itsHoldCount) {
			if (!itsRefreshHeld) { // store only first reason
				itsRefreshHeld = true;
				itsHeldReason = reason;
			}
		} else {
			clear();
			for ( auto holder : itsWCHList )
				holder->refresh(reason);
		}
	}

// Do we already have this WorldCanvasHolder/DisplayData registered?
	Bool MultiWCHolder::isAlreadyRegistered(const WorldCanvasHolder &holder) {
		return std::any_of( itsWCHList.begin( ), itsWCHList.end( ),
							[&](WorldCanvasHolder *h) { return h == &holder; } );
	}
	Bool MultiWCHolder::isAlreadyRegistered(const DisplayData &displaydata) {
		return std::any_of( itsDDList.begin( ), itsDDList.end( ),
							[&](DisplayData *dd) { return dd == &displaydata; } );
	}

// Add/remove all the DisplayDatas to/from a WorldCanvasHolder.
	void MultiWCHolder::addAllDisplayDatas(WorldCanvasHolder &holder) {
		for ( auto dd : itsDDList ) {
			holder.addDisplayData( dd, -1 );
		}
	}
	void MultiWCHolder::removeAllDisplayDatas(WorldCanvasHolder &holder, const Bool& /*permanent*/) {
		for ( auto dd : itsDDList ) {
			holder.removeDisplayData( *dd, true );
		}
	}

// Add/remove a DisplayData to/from all WorldCanvasHolders.
	void MultiWCHolder::addToAllWorldCanvasHolders(DisplayData &displaydata, int position) {
		for ( auto holder : itsWCHList ) {
			holder->addDisplayData(&displaydata, position);
		}
	}
	void MultiWCHolder::removeFromAllWorldCanvasHolders(DisplayData &displaydata) {
		for ( auto holder : itsWCHList ) {
			holder->removeDisplayData(displaydata);
		}
	}
// Distribute blinkMode to all WorldCanvasHolders.
	void MultiWCHolder::setBlinkMode( bool mode ) {
		for ( auto holder : itsWCHList ) {
			holder->setBlinkMode(mode);
		}
	}
// Distribute restrictions to all WorldCanvasHolders.
	void MultiWCHolder::distributeRestrictions() {
		for ( auto holder : itsWCHList ) {
			holder->setRestrictions(itsAttributes);
		}
	}
// Install restrictions on a specific WorldCanvasHolder.
	void MultiWCHolder::installRestrictions(WorldCanvasHolder &holder) {
		if (isAlreadyRegistered(holder)) {
			//holder.removeRestrictions();
			holder.setRestrictions(itsAttributes);
		}
	}
// This will return the maximum 'nelements' (Z axis length) of all
// dds compatible with current canvas coordinates.  (Continuum image
// can be viewed along with selected channel of spectral image, e.g.).
	uInt MultiWCHolder::zLength() {
		uInt length = 0;
		if (itsWCHList.size( ) > 0) {
			length = itsWCHList.front( )->nelements();
		}
		// Returns the value of the first wch (should be the same for
		// all of them).
		return length;
	}

// Determines which DDs will be restricted, which are always active.
// May need refinement later; for now, blink Raster PADDs only; do not
// restrict other DDs.  (Contour DDs will always show, e.g.).
// (Note that GTkPanelDisplay assumes that isBlinkDD() is false for
// GTkDrawingDDs, at present).
// (12/04: This should probably be a DD method instead, so MWCH doesn't
// need to know about various DD classes...).
	Bool MultiWCHolder::isBlinkDD(DisplayData *dd) {
		return  dd->classType() == Display::Raster   &&
		        dynamic_cast<PrincipalAxesDD*>(dd) != 0;
	}

// (permanently) sets itsBIndexName (below).  Called only in constructor.
	void MultiWCHolder::setBIndexName() {
		ostringstream os;
		os<<"bIndex"<<this;
		itsBIndexName=String(os);
	}

// Adjust "bIndex" Attribute's name to include ID of this MWCH.
	void MultiWCHolder::adjustBIndexName(AttributeBuffer& rstrs) {
		if(!rstrs.exists("bIndex")) return;
		Attribute bIndexAtt(itsBIndexName, *(rstrs.getAttributeValue("bIndex")));
		rstrs.remove("bIndex");
		rstrs.set(bIndexAtt);
	}

// Return number of blink DDs, current appropriate blink index.  Sent to
// animator (by GtkPanelDisplay, actually) when DDs are added, removed.
// The animator in turn actually orders the bIndex 'LinearRestrictions'
// to be set or removed on the WCHs.
	Int MultiWCHolder::bLength() {
		return itsBLength;
	}
	Int MultiWCHolder::bIndex() {
		return itsBIndex;
	}



	Bool MultiWCHolder::conforms(DisplayData* dd,
	                             Bool testRstrs, Bool testCS, Bool testZ,
	                             Int wchIndex) {
		// Test conformance of a DD to a WCH of this MWCH (by default, test the
		// first one (WCH 0) which always exists).  The three aspects of
		// conformance can be selectively tested.

        if ( (size_t) wchIndex >= itsWCHList.size( ) || dd==0 ) return false;
        auto wchs = itsWCHList.begin( );
        std::advance( wchs, wchIndex );

		auto wch = *wchs;
		if ( wch==0 ) return false;

		return (!testZ     || dd->conformsToZIndex(*wch->worldCanvas()))  &&
		       (!testCS    || dd->conformsToCS(*wch->worldCanvas()))      &&
		       (!testRstrs || dd->conformsToRstrs(*wch->worldCanvas()));
	}


} //# NAMESPACE CASA - END

