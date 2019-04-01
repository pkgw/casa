//# PanelDisplay.cc: A 3d Slice Application
//# Copyright (C) 2001,2002,2003
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

#include <casa/aips.h>
#include <casa/Containers/Record.h>
#include <display/Display/PixelCanvas.h>
#include <display/Display/WorldCanvas.h>
#include <display/Display/MultiWCHolder.h>
#include <display/Display/PanelDisplay.h>
#include <display/Display/SlicePanelDisplay.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	SliceEH::SliceEH(PanelDisplay* mwch, uInt axis, uInt daxis)
		: itsPD(mwch),
		  itsPositionKey(Display::K_Pointer_Button1),
		  itsMotionKey(Display::KM_Pointer_Button_1) {
		itsAxis = axis;
		itsDrawAxis = daxis;
		itsAni.addMWCHolder(*itsPD);
	}

	SliceEH::~SliceEH() {
		itsAni.removeMWCHolder(*itsPD);
	}

	void SliceEH::matchKey(Display::KeySym key) {
		if (key == Display::K_Pointer_Button1) {
			itsMotionKey =  Display::KM_Pointer_Button_1;
			itsPositionKey = Display::K_Pointer_Button1;
		} else if (key == Display::K_Pointer_Button2) {
			itsMotionKey =  Display::KM_Pointer_Button_2;
			itsPositionKey = Display::K_Pointer_Button3;
		} else if (key == Display::K_Pointer_Button3) {
			itsMotionKey =  Display::KM_Pointer_Button_3;
			itsPositionKey = Display::K_Pointer_Button3;
		} else {
			itsPositionKey = Display::K_None;
			itsMotionKey = (Display::KeyModifier)0x0000;
		}
	}

	void SliceEH::draw(const WCMotionEvent &ev) {
		static Vector<Double> linA(2);
		static Vector<Double> linB(2);
		ev.worldCanvas()->pixelCanvas()->setColor(drawColor());
		ev.worldCanvas()->pixelCanvas()->setLineWidth(lineWidth());

		ev.worldCanvas()->copyBackBufferToFrontBuffer();

		WorldCanvas* wc = 0;
		itsPD->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		if ( wc != 0 ) {
			switch(itsDrawAxis) {
			case 0: {
				linA(0) = ev.linX();
				linA(1) = wc->linYMin();
				linB(0) = ev.linX();
				linB(1) = wc->linYMax();
				break;
			}
			case 1: {
				linA(0) = wc->linXMin();
				linA(1) = ev.linY();
				linB(0) = wc->linXMax();
				linB(1) = ev.linY();
				break;
			}
			case 2: {
				linA(0) = ev.linY();
				linA(1) = wc->linYMin();
				linB(0) = ev.linY();
				linB(1) = wc->linYMax();
				break;
			}
			case 3: {
				linA(0) = wc->linXMin();
				linA(1) = ev.linX();
				linB(0) = wc->linXMax();
				linB(1) = ev.linX();
				break;
			}
			default:
				break;
			}
		}
		wc->copyBackBufferToFrontBuffer();
		wc->drawLine(linA,linB,true);
		//clear
		linA(0) = ev.linX();
		linA(1) = ev.worldCanvas()->linYMin();
		linB(0) = ev.linX();
		linB(1) = ev.worldCanvas()->linYMax();
		ev.worldCanvas()->drawLine(linA,linB,true);
		linA(0) = ev.worldCanvas()->linXMin();
		linA(1) = ev.linY();
		linB(0) = ev.worldCanvas()->linXMax();
		linB(1) = ev.linY();
		ev.worldCanvas()->drawLine(linA,linB,true);
	}

	void SliceEH::operator()(const WCPositionEvent &ev) {
		if (ev.key() == itsPositionKey) {
			Int x,y;
			x = Int(ev.pixX()+0.5);
			y = Int(ev.pixY()+0.5);
			if (!ev.worldCanvas()->inDrawArea(x, y)) {
				return;
			}
			Int value;
			if (itsAxis == 0) {
				value = static_cast<Int>(ev.linX()+0.5);
			} else if (itsAxis == 1) {
				value = static_cast<Int>(ev.linY()+0.5);
			} else {
				return;
			}
			itsAni.setLinearRestriction("zIndex",value,0,0);
			draw(ev);
		} else if (ev.key() == Display::K_Escape) {
			ev.worldCanvas()->pixelCanvas()->copyBackBufferToFrontBuffer();
		}
	}
	void SliceEH::operator()(const WCMotionEvent &ev) {
		if (!(ev.modifiers() & itsMotionKey)) {
			return;
		}
		Int x,y;
		x = Int(ev.pixX()+0.5);
		y = Int(ev.pixY()+0.5);
		if (!ev.worldCanvas()->inDrawArea(x, y)) {
			return;
		}
		Int value;
		if (itsAxis == 0) {
			value = static_cast<Int>(ev.linX()+0.5);
		} else if (itsAxis == 1) {
			value = static_cast<Int>(ev.linY()+0.5);
		} else {
			return;
		}
		itsAni.setLinearRestriction("zIndex",value,0,0);
		draw(ev);
	}
// end SliceEH

	SlicePanelDisplay::SlicePanelDisplay(PixelCanvas* pcanvas) :
		itsActivePD(0) {
		itsPanelDisplays["xy"] = new PanelDisplay(pcanvas,1,1,0.0,0.5,0.5,0.5);
		itsPanelDisplays["zy"] = new PanelDisplay(pcanvas,1,1,0.5,0.5,0.5,0.5);
		itsPanelDisplays["xz"] = new PanelDisplay(pcanvas,1,1,0.0,0.0,0.5,0.5);
		itsSliceEHs["ZYcXY"] =	new SliceEH(itsPanelDisplays["xy"],0,1);
		itsSliceEHs["XZcXY"] =	new SliceEH(itsPanelDisplays["xy"],1,0);
		itsSliceEHs["XYcZY"] =	new SliceEH(itsPanelDisplays["zy"],0,1);
		itsSliceEHs["XZcZY"] =	new SliceEH(itsPanelDisplays["zy"],0,2);
		itsSliceEHs["XYcXZ"] =	new SliceEH(itsPanelDisplays["xz"],1,0);
		itsSliceEHs["ZYcXZ"] =	new SliceEH(itsPanelDisplays["xz"],1,3);
		installEHs();
		itsActivePD = itsPanelDisplays["xy"];
	}

	SlicePanelDisplay::~SlicePanelDisplay() {
		removeEHs();
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			PanelDisplay* pdisp = iter->second;
			if (pdisp) delete pdisp;
		}
		for (auto iter = itsSliceEHs.begin( ); iter != itsSliceEHs.end( ); ++iter) {
			SliceEH* seh = iter->second;
			if (seh) delete seh;
		}
	}

	void SlicePanelDisplay::installEHs() {
		WCMotionEH* meh = 0;
		WCMotionEH* meh1 = 0;
		WCPositionEH* peh = 0;
		WCPositionEH* peh1 = 0;
		PanelDisplay* pdisp = 0;

		WorldCanvas* wc = 0;

		pdisp = itsPanelDisplays["xy"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["XYcZY"];
		peh = itsSliceEHs["XYcZY"];
		meh1 = itsSliceEHs["XYcXZ"];
		peh1 = itsSliceEHs["XYcXZ"];

		if ( wc != 0 ) {
			wc->addMotionEventHandler(*meh);
			wc->addPositionEventHandler(*peh);
			wc->addMotionEventHandler(*meh1);
			wc->addPositionEventHandler(*peh1);
			wc->addMotionEventHandler(*this);
		}

		pdisp = itsPanelDisplays["zy"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["ZYcXY"];
		peh = itsSliceEHs["ZYcXY"];
		meh1 = itsSliceEHs["ZYcXZ"];
		peh1 = itsSliceEHs["ZYcXZ"];

		if ( wc != 0 ) {
			wc->addMotionEventHandler(*meh);
			wc->addPositionEventHandler(*peh);
			wc->addMotionEventHandler(*meh1);
			wc->addPositionEventHandler(*peh1);
			wc->addMotionEventHandler(*this);
		}

		pdisp = itsPanelDisplays["xz"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["XZcXY"];
		peh = itsSliceEHs["XZcXY"];
		meh1 = itsSliceEHs["XZcZY"];
		peh1 = itsSliceEHs["XZcZY"];

		if ( wc != 0 ) {
			wc->addMotionEventHandler(*meh);
			wc->addPositionEventHandler(*peh);
			wc->addMotionEventHandler(*meh1);
			wc->addPositionEventHandler(*peh1);
			wc->addMotionEventHandler(*this);
		}
	}

	void SlicePanelDisplay::removeEHs() {
		WCMotionEH* meh = 0;
		WCMotionEH* meh1 = 0;
		WCPositionEH* peh = 0;
		WCPositionEH* peh1 = 0;
		PanelDisplay* pdisp = 0;
		WorldCanvas* wc = 0;

		pdisp = itsPanelDisplays["xy"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["XYcZY"];
		peh = itsSliceEHs["XYcZY"];
		meh1 = itsSliceEHs["XYcXZ"];
		peh1 = itsSliceEHs["XYcXZ"];

		if ( wc != 0 ) {
			wc->removeMotionEventHandler(*meh);
			wc->removePositionEventHandler(*peh);
			wc->removeMotionEventHandler(*meh1);
			wc->removePositionEventHandler(*peh1);
			wc->removeMotionEventHandler(*this);
		}

		pdisp = itsPanelDisplays["zy"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["ZYcXY"];
		peh = itsSliceEHs["ZYcXY"];
		meh1 = itsSliceEHs["ZYcXZ"];
		peh1 = itsSliceEHs["ZYcXZ"];

		if ( wc != 0 ) {
			wc->removeMotionEventHandler(*meh);
			wc->removePositionEventHandler(*peh);
			wc->removeMotionEventHandler(*meh1);
			wc->removePositionEventHandler(*peh1);
			wc->removeMotionEventHandler(*this);
		}

		pdisp = itsPanelDisplays["xz"];
		pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

		meh = itsSliceEHs["XZcXY"];
		peh = itsSliceEHs["XZcXY"];
		meh1 = itsSliceEHs["XZcZY"];
		peh1 = itsSliceEHs["XZcZY"];

		if ( wc != 0 ) {
			wc->removeMotionEventHandler(*meh);
			wc->removePositionEventHandler(*peh);
			wc->removeMotionEventHandler(*meh1);
			wc->removePositionEventHandler(*peh1);
			wc->removeMotionEventHandler(*this);
		}
	}

	void  SlicePanelDisplay::refresh(const Display::RefreshReason& reason) {
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			iter->second->refresh(reason);
		}
	}

	void SlicePanelDisplay::hold() {
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			iter->second->hold();
		}
	}

	void SlicePanelDisplay::release() {
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			iter->second->release();
		}
	}


	void SlicePanelDisplay::enableTools() {
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			PanelDisplay* pdisp = iter->second;
			pdisp->enableTools();
		}
	}

	void SlicePanelDisplay::disableTools() {
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			PanelDisplay* pdisp = iter->second;
			pdisp->disableTools();
		}
	}

	void SlicePanelDisplay::setToolKey(const String& toolname,
	                                   const Display::KeySym& keysym) {
		if (toolname == "slice") {
			for (auto iter = itsSliceEHs.begin( ); iter != itsSliceEHs.end( ); ++iter) {
				iter->second->matchKey(keysym);
			}
			return;
		}
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			PanelDisplay* pdisp = iter->second;
			pdisp->setToolKey(toolname, keysym);
		}
	}


	Bool SlicePanelDisplay::setOptions(const Record& rec, Record& out) {
		// distribute
		Bool needRefresh = false;
		PanelDisplay* pdisp = 0;
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			pdisp = iter->second;
			needRefresh = pdisp->setOptions(rec,out);
			//pdisp->refresh();
		}
		return needRefresh;
	}

	Record SlicePanelDisplay::getOptions() const {
		// driven by the XY panel
        auto ptr = itsPanelDisplays.find("xy");
		return ptr->second->getOptions();
	}

	PanelDisplay* SlicePanelDisplay::getPanelDisplay(const String& pdname) {
        auto pdptr = itsPanelDisplays.find(pdname);
        return pdptr != itsPanelDisplays.end( ) ? pdptr->second : 0;
	}

	void SlicePanelDisplay::operator()(const WCMotionEvent& ev) {
		PanelDisplay* pdisp = 0;
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {

			pdisp = iter->second;

			WorldCanvas* wc = 0;
			pdisp->wcsApply( [&](WorldCanvas *w) { if ( wc == 0 ) wc = w; } );

			if ( wc != 0 ) {
				if ( wc == ev.worldCanvas() ) {
					itsActivePD = pdisp;
					break;
				}
			}
		}
	}

	void SlicePanelDisplay::precompute() {
		PanelDisplay* pdisp = 0;
		MWCAnimator tmpAni;
		uInt length = 0;
		for (auto iter = itsPanelDisplays.begin( ); iter != itsPanelDisplays.end( ); ++iter) {
			pdisp = iter->second;
			tmpAni.addMWCHolder(*pdisp);
			length = pdisp->zLength();
			for (uInt k = 0; k < length; k++) {
				tmpAni.setLinearRestriction("zIndex",static_cast<Int>(k),0,0);
			}
			tmpAni.removeMWCHolder(*pdisp);
		}
	}


} //# NAMESPACE CASA - END

