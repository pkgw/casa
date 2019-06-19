//# MultiWCTool.cc: base class for MultiWorldCanvas event-based tools
//# Copyright (C) 2000,2001,2002
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

#include <display/Display/WorldCanvas.h>
#include <display/Display/PanelDisplay.h>
#include <display/DisplayEvents/MultiWCTool.h>
#include <display/Display/PixelCanvas.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// (Required) default constructor.
	MultiWCTool::MultiWCTool() :
		DisplayTool(),
		itsCurrentWC(0),
		itsEventHandlersRegistered(false) {
	}

	MultiWCTool::MultiWCTool(const Display::KeySym &keysym, bool enable_events ) :
		DisplayTool(keysym),
		itsCurrentWC(0),
		itsEventHandlersRegistered(false) {
		if ( enable_events ) enable();
	}

	MultiWCTool::~MultiWCTool() {
		MultiWCTool::disable();
	}

// (Required) copy constructor.
	MultiWCTool::MultiWCTool(const MultiWCTool &other) :
		DisplayTool(other), WCPositionEH( ), WCMotionEH( ), WCRefreshEH( ) {  }

// (Required) copy assignment.
	MultiWCTool &MultiWCTool::operator=(const MultiWCTool &other) {
		if (this != &other) DisplayTool::operator=(other);
		return *this;
	}

	void MultiWCTool::addWorldCanvas(WorldCanvas &worldcanvas) {
		itsWCList.push_back( &worldcanvas );
		if (itsEventHandlersRegistered) {
			worldcanvas.addPositionEventHandler(*this);
			worldcanvas.addMotionEventHandler(*this);
			worldcanvas.addRefreshEventHandler(*this);
		}
	}

	void MultiWCTool::removeWorldCanvas(WorldCanvas &worldcanvas) {
		std::list<WorldCanvas*> orig = itsWCList;
		std::list<WorldCanvas*> removed;
		itsWCList.clear( );
		std::partition_copy( orig.begin( ), orig.end( ),
							 std::back_inserter(removed),
							 std::back_inserter(itsWCList),
							 [&](WorldCanvas *wc){return wc == &worldcanvas;} );
		if ( removed.size( ) > 0 && itsEventHandlersRegistered ) {
			worldcanvas.removePositionEventHandler(*this);
			worldcanvas.removeMotionEventHandler(*this);
			worldcanvas.removeRefreshEventHandler(*this);
		}
	}

	void MultiWCTool::addWorldCanvases(PanelDisplay* pdisp) {
		pdisp->wcsApply( [&]( WorldCanvas *wcanvas ) {
							itsWCList.push_back(wcanvas);
							if (itsEventHandlersRegistered) {
								wcanvas->addPositionEventHandler(*this);
								wcanvas->addMotionEventHandler(*this);
								wcanvas->addRefreshEventHandler(*this);
							}
			} );
	}

	void MultiWCTool::removeWorldCanvases(PanelDisplay* pdisp) {
		disable();

		std::list<WorldCanvas*> orig = itsWCList;
		itsWCList.clear( );
		for ( auto wc : orig ) {
			Bool found = false;
			pdisp->wcsApply( [&]( WorldCanvas *wcanvas ) { found = found || wcanvas == wc; } );
			if( ! found ) itsWCList.push_back(wc);
		}

		enable();
	}

	void MultiWCTool::enable() {
		if ( ! itsEventHandlersRegistered ) {
			itsEventHandlersRegistered = true;
			for ( auto wc : itsWCList ) {
				wc->addPositionEventHandler(*this);
				wc->addMotionEventHandler(*this);
				wc->addRefreshEventHandler(*this);
			}
		}
	}

	void MultiWCTool::disable() {
		if ( itsEventHandlersRegistered ) {
			itsEventHandlersRegistered = false;
			for ( auto wc : itsWCList ) {
				wc->removePositionEventHandler(*this);
				wc->removeMotionEventHandler(*this);
				wc->removeRefreshEventHandler(*this);
			}
		}
	}

	void MultiWCTool::operator()(const WCPositionEvent &ev) {
		if (ev.key() != getKey() || getKey( ) == Display::K_None ) {
			if (ev.keystate()) {
				otherKeyPressed(ev);
			} else {
				otherKeyReleased(ev);
			}
		} else {
			if (ev.keystate()) {
				keyPressed(ev);
			} else {
				keyReleased(ev);
			}
		}
	}

	void MultiWCTool::refresh() {
		if(itsCurrentWC==0) return;
		PixelCanvas *pc = itsCurrentWC->pixelCanvas();
		if(pc==0) return;
		pc->copyBackBufferToFrontBuffer();
		pc->setDrawBuffer(Display::FrontBuffer);
		pc->callRefreshEventHandlers(Display::BackCopiedToFront);
	}


	void MultiWCTool::operator()(const WCRefreshEvent &ev) {
		static viewer::region::region_list_type empty;
		if (	itsCurrentWC != 0 &&
		        ev.worldCanvas() == itsCurrentWC &&
		        ev.reason() == Display::BackCopiedToFront &&
		        itsCurrentWC->pixelCanvas()->drawBuffer()==Display::FrontBuffer  )
			draw(ev,empty);
	}

	void MultiWCTool::operator()(const WCMotionEvent &ev) {
		static viewer::region::region_list_type empty;
		moved(ev,empty);
	}

	void MultiWCTool::setClipToDrawArea() {
		WorldCanvas *wc = itsCurrentWC;
		if(wc==0) return;
		PixelCanvas *pc = wc->pixelCanvas();
		if(pc==0) return;
		Int x0 = wc->canvasXOffset() + wc->canvasDrawXOffset();
		Int x1 = x0 + wc->canvasDrawXSize() - 1;
		Int y0 = wc->canvasYOffset() + wc->canvasDrawYOffset();
		Int y1 = y0 + wc->canvasDrawYSize() - 1;
		pc->setClipWindow(x0,y0, x1,y1);
		pc->enable(Display::ClipWindow);
	}

	void MultiWCTool::setClipToWC() {
		WorldCanvas *wc = itsCurrentWC;
		if(wc==0) return;
		PixelCanvas *pc = wc->pixelCanvas();
		if(pc==0) return;
		Int x0 = wc->canvasXOffset();
		Int x1 = x0 + wc->canvasXSize() - 1;
		Int y0 = wc->canvasYOffset();
		Int y1 = y0 + wc->canvasYSize() - 1;
		pc->setClipWindow(x0,y0, x1,y1);
		pc->enable(Display::ClipWindow);
	}

	void MultiWCTool::resetClip() {
		WorldCanvas *wc = itsCurrentWC;
		if(wc==0) return;
		PixelCanvas *pc = wc->pixelCanvas();
		if(pc==0) return;
		pc->disable(Display::ClipWindow);
	}

// Callbacks: responses to events.  To be implemented by derived classes
// as needed.

	void MultiWCTool::keyPressed(const WCPositionEvent &/*ev*/) {  }
	void MultiWCTool::keyReleased(const WCPositionEvent &) {  }
	void MultiWCTool::otherKeyPressed(const WCPositionEvent &) {  }
	void MultiWCTool::otherKeyReleased(const WCPositionEvent &) {  }
	void MultiWCTool::moved(const WCMotionEvent & /*ev*/, const viewer::region::region_list_type & /*selected_regions*/) { }
	void MultiWCTool::draw(const WCRefreshEvent&/*ev*/, const viewer::region::region_list_type & /*selected_regions*/) {  }



} //# NAMESPACE CASA - END

