//# PanelDisplay.cc: Provision of panelled displays for data
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
#include <casa/Exceptions.h>
#include <casa/Containers/Record.h>
#include <display/Display/Attribute.h>
#include <display/Display/AttributeBuffer.h>
#include <display/Display/PixelCanvas.h>
#include <display/Display/WorldCanvas.h>
#include <display/Display/WorldCanvasHolder.h>
#include <display/Display/PanelDisplay.h>
#include <display/DisplayEvents/MultiWCTool.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	const String PanelDisplay::X_ORIGIN="xorigin";
	const String PanelDisplay::Y_ORIGIN = "yorigin";
	const String PanelDisplay::X_SIZE = "xsize";
	const String PanelDisplay::Y_SIZE = "ysize";

	// get controlled access to world canvases shared among a number of objects
	// return value indicates if the operation was possible (in the future it
	// may be necessary to serialize access so this function may return false
	// if mutual exclusion prevents access to the world canvas list)
	bool PanelDisplay::wcsApply( std::function<void(WorldCanvas *)> apply ) {
		for ( auto wc : itsWCList ) apply(wc);
		return true;
	}

// Constructor.
	PanelDisplay::PanelDisplay(PixelCanvas* pixelcanvas,
	                           const Int nx, const Int ny,
	                           const Float xOrigin, const Float yOrigin,
	                           const Float xSize, const Float ySize,
	                           const Float dx, const Float dy,
	                           const PanelDisplay::FillOrder order) :
		MultiWCHolder(),
		itsPixelCanvas(pixelcanvas),
		itsGeometrySet(false) {

		itslpgm =10;
		itsrpgm = 1; //4
		itstpgm = 1; //4
		itsbpgm = 7;
		setGeometry(nx, ny, xOrigin, yOrigin, xSize, ySize, dx, dy, order);
	}

// Destructor.
	PanelDisplay::~PanelDisplay() {
		unSetupGeometry();
		// cleanup Tools
		while (itsMWCTools.size() > 0u) {
			String key = itsMWCTools.begin( )->first;
			removeTool(key);
		}
		itsMWCTools.clear();
		/*if (itsWCLI) {
		  delete itsWCLI;
		}
		if (itsWCHLI) {
		  delete itsWCHLI;
		}*/
	}

	void PanelDisplay::setAttributes(AttributeBuffer& at) {
		for ( auto wc : itsWCList ) {
			wc->setAttributes(at);
		}
	}

	void PanelDisplay::getAttributeValue(const String& name, Int& newValue) const {
		if ( itsWCList.size( ) > 0 )
			itsWCList.front( )->getAttributeValue(name, newValue);
	}

// Option handling functions.
	void PanelDisplay::setDefaultOptions() {
		String attString;
		AttributeBuffer attBuffer;
		attString = "leftMarginSpacePG";
		attBuffer.add(attString, itslpgm);
		attString = "rightMarginSpacePG";
		attBuffer.add(attString, itsrpgm);
		attString = "bottomMarginSpacePG";
		attBuffer.add(attString, itsbpgm);
		attString = "topMarginSpacePG";
		attBuffer.add(attString, itstpgm);
		setAttributes(attBuffer);

	}

	Record PanelDisplay::getOptions() const {
		Record rec;
		String attString;
		Int temp;

		Record leftmarginspacepg;
		leftmarginspacepg.define("dlformat", WorldCanvas::LEFT_MARGIN_SPACE_PG);
		leftmarginspacepg.define("listname", "Left margin space (PG chars)");
		leftmarginspacepg.define("ptype", "intrange");
		leftmarginspacepg.define("pmin", 0);
		leftmarginspacepg.define("pmax", 20);
		leftmarginspacepg.define("default", itslpgm);
		attString = "leftMarginSpacePG";
		getAttributeValue(attString, temp);
		leftmarginspacepg.define("value", temp);
		leftmarginspacepg.define("allowunset", false);
		leftmarginspacepg.define("context", "Margins");
		rec.defineRecord(WorldCanvas::LEFT_MARGIN_SPACE_PG, leftmarginspacepg);

		Record bottommarginspacepg;
		bottommarginspacepg.define("dlformat", WorldCanvas::BOTTOM_MARGIN_SPACE_PG);
		bottommarginspacepg.define("listname", "Bottom margin space (PG chars)");
		bottommarginspacepg.define("ptype", "intrange");
		bottommarginspacepg.define("pmin", 0);
		bottommarginspacepg.define("pmax", 20);
		bottommarginspacepg.define("default", itsbpgm);
		attString = "bottomMarginSpacePG";
		getAttributeValue(attString, temp);
		bottommarginspacepg.define("value", temp);
		bottommarginspacepg.define("allowunset", false);
		bottommarginspacepg.define("context", "Margins");
		rec.defineRecord(WorldCanvas::BOTTOM_MARGIN_SPACE_PG, bottommarginspacepg);

		Record rightmarginspacepg;
		rightmarginspacepg.define("dlformat", WorldCanvas::RIGHT_MARGIN_SPACE_PG);
		rightmarginspacepg.define("listname", "Right margin space (PG chars)");
		rightmarginspacepg.define("ptype", "intrange");
		rightmarginspacepg.define("pmin", 0);
		rightmarginspacepg.define("pmax", 20);
		rightmarginspacepg.define("default", itsrpgm);
		attString = "rightMarginSpacePG";
		getAttributeValue(attString, temp);
		rightmarginspacepg.define("value", temp);
		rightmarginspacepg.define("allowunset", false);
		rightmarginspacepg.define("context", "Margins");
		rec.defineRecord(WorldCanvas::RIGHT_MARGIN_SPACE_PG, rightmarginspacepg);

		Record topmarginspacepg;
		topmarginspacepg.define("dlformat", WorldCanvas::TOP_MARGIN_SPACE_PG);
		topmarginspacepg.define("listname", "Top margin space (PG chars)");
		topmarginspacepg.define("ptype", "intrange");
		topmarginspacepg.define("pmin", 0);
		topmarginspacepg.define("pmax", 20);
		topmarginspacepg.define("default", itstpgm);
		attString = "topMarginSpacePG";
		getAttributeValue(attString, temp);
		topmarginspacepg.define("value", temp);
		topmarginspacepg.define("allowunset", false);
		topmarginspacepg.define("context", "Margins");
		rec.defineRecord(WorldCanvas::TOP_MARGIN_SPACE_PG, topmarginspacepg);

		Record nxpanels;
		nxpanels.define("dlformat", "nxpanels");
		nxpanels.define("listname", "Number of panels in x");
		nxpanels.define("ptype", "intrange");
		nxpanels.define("pmin", Int(1));
		nxpanels.define("pmax", Int(5));
		nxpanels.define("default", Int(1));
		nxpanels.define("value", itsNX);
		nxpanels.define("allowunset", false);
		nxpanels.define("context", "Number_of_panels");
		rec.defineRecord("nxpanels", nxpanels);

		Record nypanels;
		nypanels.define("dlformat", "nypanels");
		nypanels.define("listname", "Number of panels in y");
		nypanels.define("ptype", "intrange");
		nypanels.define("pmin", Int(1));
		nypanels.define("pmax", Int(5));
		nypanels.define("default", Int(1));
		nypanels.define("value", itsNY);
		nypanels.define("allowunset", false);
		nypanels.define("context", "Number_of_panels");
		rec.defineRecord("nypanels", nypanels);

		Record xspacing;
		xspacing.define("dlformat", "xspacing");
		xspacing.define("listname", "X-Spacing of Panels");
		xspacing.define("ptype", "floatrange");
		xspacing.define("pmin", 0.0);
		xspacing.define("pmax", 0.49);
		xspacing.define("default", Float(0.0));
		xspacing.define("value", itsDX);
		xspacing.define("allowunset", false);
		xspacing.define("context", "Number_of_panels");
		rec.defineRecord("xspacing", xspacing);

		Record yspacing;
		yspacing.define("dlformat", "yspacing");
		yspacing.define("listname", "Y-Spacing of Panels");
		yspacing.define("ptype", "floatrange");
		yspacing.define("pmin", 0.0);
		yspacing.define("pmax", 0.49);
		yspacing.define("default", Float(0.0));
		yspacing.define("value", itsDY);
		yspacing.define("allowunset", false);
		yspacing.define("context", "Number_of_panels");
		rec.defineRecord("yspacing", yspacing);

		return rec;

	}

	Bool PanelDisplay::setOptions(const Record& rec, Record& ) {
		Bool ret = false, localchange = false;
		Bool error;

		String attString;
		AttributeBuffer attBuffer;

		Bool geometrychange = false;
		geometrychange = readOptionRecord(itsNX, error, rec,  "nxpanels")
		                 || geometrychange;
		geometrychange = readOptionRecord(itsNY, error, rec,  "nypanels")
		                 || geometrychange;
		geometrychange = readOptionRecord(itsDX, error, rec,  "xspacing")
		                 || geometrychange;
		geometrychange = readOptionRecord(itsDY, error, rec,  "yspacing")
		                 || geometrychange;

		if (geometrychange) {
			setGeometry(itsNX,itsNY,itsXOrigin,itsYOrigin,
			            itsXSize,itsYSize,itsDX,itsDY,itsOrder);
			localchange = true;
		}

		// set distributed options

		attString = "leftMarginSpacePG";
		localchange = readOptionRecord(itslpgm, error, rec, WorldCanvas::LEFT_MARGIN_SPACE_PG)
		              || localchange;
		attBuffer.add(attString, itslpgm);

		attString = "rightMarginSpacePG";
		localchange = readOptionRecord(itsrpgm, error, rec, WorldCanvas::RIGHT_MARGIN_SPACE_PG)
		              || localchange;
		attBuffer.add(attString, itsrpgm);

		attString = "bottomMarginSpacePG";
		localchange = readOptionRecord(itsbpgm, error, rec, WorldCanvas::BOTTOM_MARGIN_SPACE_PG)
		              || localchange;
		attBuffer.add(attString, itsbpgm);

		attString = "topMarginSpacePG";
		localchange = readOptionRecord(itstpgm, error, rec, WorldCanvas::TOP_MARGIN_SPACE_PG)
		              || localchange;
		attBuffer.add(attString, itstpgm);
		setAttributes(attBuffer);

		ret = ret || localchange;
		return ret;
	}

	void PanelDisplay::getGeometry(Int& nx, Int& ny, Float& xOrigin,
	                               Float& yOrigin, Float& xSize, Float& ySize,
	                               Float& dx, Float& dy,
	                               PanelDisplay::FillOrder& order) const {
		nx = itsNX;
		ny = itsNY;
		xOrigin = itsXOrigin;
		yOrigin = itsYOrigin;
		xSize = itsXSize;
		ySize = itsYSize;
		dx = itsDX;
		dy = itsDY;
		order = itsOrder;
	}

	void PanelDisplay::getGeometry(RecordInterface& rec) const {
		rec.define("nxpanels", itsNX);
		rec.define("nypanels", itsNY);
		rec.define(X_ORIGIN, itsXOrigin);
		rec.define(Y_ORIGIN, itsYOrigin);
		rec.define(X_SIZE, itsXSize);
		rec.define(Y_SIZE, itsYSize);
		rec.define("xspacing", itsDX);
		rec.define("yspacing", itsDY);
	}


	void PanelDisplay::setGeometry(const RecordInterface& rec) {
		if (rec.isDefined("nxpanels") && rec.isDefined("nypanels") &&
		        rec.isDefined(X_ORIGIN) && rec.isDefined(Y_ORIGIN) &&
		        rec.isDefined(X_SIZE) && rec.isDefined(Y_SIZE) &&
		        rec.isDefined("xspacing") && rec.isDefined("yspacing")) {

			Int nx,ny;
			Float xOrigin,yOrigin, xSize,ySize, dx,dy;

			rec.get("nxpanels", nx);
			rec.get("nypanels", ny);
			rec.get(X_ORIGIN, xOrigin);
			rec.get(Y_ORIGIN, yOrigin);
			rec.get(X_SIZE, xSize);
			rec.get(Y_SIZE, ySize);
			rec.get("xspacing", dx);
			rec.get("yspacing", dy);

			setGeometry(nx,ny, xOrigin,yOrigin, xSize,ySize, dx,dy, itsOrder);

		} else {
			throw(AipsError("In PanelDisplay, the geometry description record "
			                "was missing one or more fields"));
		}
	}




	void PanelDisplay::clear() {


		// clear within screen pixel extents of the PD on the PC)

		Int pw = pixelCanvas()->width(),
		    ph = pixelCanvas()->height();

		Int xmin = Int(itsXOrigin*pw +.5),
		    ymin = Int(itsYOrigin*ph +.5),
		    xmax = Int((itsXOrigin+itsXSize)*pw +.5) - 1,
		    ymax = Int((itsYOrigin+itsYSize)*ph +.5) - 1;

		pixelCanvas()->clear(xmin, ymin,  xmax, ymax);

		if(pixelCanvas()->drawMode()==Display::Compile) return;

		// In draw mode, clear both front and back buffers.

		Display::DrawBuffer origbuf = itsPixelCanvas->drawBuffer();
		Display::DrawBuffer altbuf = (origbuf==Display::BackBuffer)?
		                             Display::FrontBuffer : Display::BackBuffer;

		pixelCanvas()->setDrawBuffer(altbuf);
		pixelCanvas()->clear(xmin, ymin,  xmax, ymax);

		pixelCanvas()->setDrawBuffer(origbuf);
	}



	void PanelDisplay::setGeometry(const Int nx, const Int ny,
	                               const Float xOrigin, const Float yOrigin,
	                               const Float xSize, const Float ySize,
	                               const Float dx, const Float dy,
	                               const PanelDisplay::FillOrder order) {
		if ((nx < 1) || (ny < 1)) {
			throw(AipsError("In PanelDisplay, the number of panels in each direction "
			                "must be at least one"));
		}
		itsNX = nx;
		itsNY = ny;
		itsXOrigin = xOrigin;
		itsYOrigin = yOrigin;
		itsXSize = xSize;
		itsYSize = ySize;
		itsDX = dx;
		itsDY = dy;
		itsOrder = order;


		updateTools(true,false);	// (remove them temporarily).


		Float xPanelSize = (itsXSize - static_cast<Float>(itsNX - 1) * itsDX) /
		                   static_cast<Float>(itsNX);
		Float yPanelSize = (itsYSize - static_cast<Float>(itsNY - 1) * itsDY) /
		                   static_cast<Float>(itsNY);

		// Prepare to synchornize zoom windows and CS master of any new WCs with
		// existing ones.

		AttributeBuffer zoomwindow;

		WorldCanvas* wc0 = 0;
		WorldCanvasHolder* wch0 = 0;

		Bool oldWCexists= (itsWCList.size( ) != 0);

		if(oldWCexists) {

			wc0 = itsWCList.front( );
			wch0 = itsWCHList.front( );

			Vector<Double> zoomBlc(2), zoomTrc(2);
			zoomBlc[0]=wc0->linXMin();
			zoomBlc[1]=wc0->linYMin();
			zoomTrc[0]=wc0->linXMax();
			zoomTrc[1]=wc0->linYMax();

			zoomwindow.add("manualZoomBlc", zoomBlc);
			zoomwindow.add("manualZoomTrc", zoomTrc);
		}


		Float y = itsYOrigin + itsYSize - yPanelSize;
		std::list<WorldCanvas*>::iterator wc_iter = itsWCList.begin( );
        std::list<WorldCanvasHolder*>::iterator wch_iter = itsWCHList.begin( );
		for (Int i = 0; i < itsNY; i++) {
			Float x = itsXOrigin;
			for (Int j = 0; j < itsNX; j++) {

				if ( wc_iter == itsWCList.end( ) ) {

					// out of WC[H]s--create new ones

					WorldCanvas* wc = new WorldCanvas(itsPixelCanvas,
					                                  x,y, xPanelSize,yPanelSize);
					WorldCanvasHolder* wch = new WorldCanvasHolder(wc);

					itsWCHList.insert(wch_iter,wch);
					itsWCList.insert(wc_iter,wc);

					// (To be fixed): _two identical WCH lists_ are maintained
					// (an oversight, no doubt).
					// The statement above adds to the one on this level; the
					// statement below adds to the base-level (MWCH) list. (dk 12/04)

					addWCHolder(*wch);	// (this one loads the DDs too).

					if(oldWCexists) {
						wch->syncCSmaster(wch0);
						// Use same CS master dd on new WC[H]s as on the other ones.
						// This also makes sure the new CS master (if any) sets
						// initial WC coordinate state.

						wc->setAttributes(zoomwindow);
						// The new WC will duplicate the zoom window of
						// the others on its first refresh.
					}
				} else {

					// just recycle / reposition old ones.

					(*wc_iter)->setWorldCanvasPosition(x,y, xPanelSize,yPanelSize);

					++wc_iter;
                    ++wch_iter;
				}

				x += xPanelSize + itsDX;
			}
			y -= yPanelSize + itsDY;
		}

		// remove any leftover WC{H}s
		std::list<WorldCanvasHolder*> wch_del;
		std::list<WorldCanvas*> wc_del;
		if ( wch_iter != itsWCHList.end( ) ) {
			wch_del.splice( wch_del.begin( ),
							itsWCHList, wch_iter, itsWCHList.end( ) );
		}
		if ( wc_iter != itsWCList.end( ) ) {
			wc_del.splice( wc_del.begin( ),
						   itsWCList, wc_iter, itsWCList.end( ) );
		}

        for ( auto d : wch_del ) {
            // lo, our parent (MultiWCHolder) maintains a DUPLICATE list
            // of WorldCanvasHolders, and it's even called the same thing
            // so convienient...
            removeWCHolder(*d);
            delete d;
        }
        for ( auto w : wc_del ) delete w;

		updateTools(false,true);	// (restore mouse tools).
		setDefaultOptions();
		itsGeometrySet = true;

		// can't refresh before realize is called, this would be called in
		// the 'initialize' process of X -> commented out.
		//itsPixelCanvas->refresh();

	}


	void PanelDisplay::unSetupGeometry() {
		if (!itsGeometrySet) {
			return;
		}
		updateTools(true,false);
		// 1. remove the WorldCanvasHolders
		std::list<WorldCanvasHolder*> wch_orig = itsWCHList;
		itsWCHList.clear( );
		for ( auto wch : wch_orig ) delete wch;
		// 2. delete WorldCanvases.
		std::list<WorldCanvas*> wc_orig = itsWCList;
		itsWCList.clear( );
		for ( auto wc : wc_orig ) delete wc;
		itsGeometrySet = false;
		// we have remove WorldCanvases from the PixelCanvas, so we should
		// refresh the entire PixelCanvas.

		//itsPixelCanvas->refresh();
	}

	WorldCanvasHolder* PanelDisplay::wcHolder(WorldCanvas* wcanvas) const {
		for ( auto wch : itsWCHList ) {
			if (wch->worldCanvas( ) == wcanvas) {
				return wch;
			}
		}
		return 0;
	}

	void PanelDisplay::setCSmaster( DisplayData* dd ) {
		for ( auto wch : itsWCHList ) {
			wch->setCSMaster(dd);
		}
		itsPixelCanvas->refresh(Display::WorldCoordinateChange, true);
	}

	Bool PanelDisplay::isCSmaster(const DisplayData *dd) const {
		// Is the specified DisplayData the one in charge of coordinate
		// state of the Panel's WCs?
		return itsWCHList.size( ) > 0 && itsWCHList.front( )->isCSmaster(dd);
	}


	Bool PanelDisplay::hasTools() {
		if (itsMWCTools.size() > 0) {
			return true;
		} else {
			return false;
		}
	}


	void PanelDisplay::updateTools(Bool remove, Bool add) {
		if (itsMWCTools.size() == 0) {
			return;
		}
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			if (remove) {
				iter->second->removeWorldCanvases(this);
			}
			if (add) {
				iter->second->addWorldCanvases(this);
			}
		}
	}

	void PanelDisplay::disableTools() {
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			iter->second->disable( );
		}
	}

	void PanelDisplay::enableTools() {
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			iter->second->enable();
		}
	}

	void PanelDisplay::enableTool(const String& toolname) {
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			if (iter->first == toolname) {
				iter->second->enable();
				break;
			}
		}
	}

	void PanelDisplay::disableTool(const String& toolname) {
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			if (iter->first == toolname) {
				iter->second->disable();
				break;
			}
		}
	}

	void PanelDisplay::setToolKey(const String& toolname,
	                              const Display::KeySym& keysym) {
		for (auto iter = itsMWCTools.begin( ); iter != itsMWCTools.end( ); ++iter) {
			if (iter->first == toolname) {
				iter->second->setKey(keysym);
				break;
			}
		}
	}

	void PanelDisplay::addTool(const String& key, const std::shared_ptr<MultiWCTool> &value) {
        auto wcptr = itsMWCTools.find(key);
		if ( wcptr == itsMWCTools.end( ) ) {
			itsMWCTools[key] = value;
			value->addWorldCanvases(this);
		}
	}

	void PanelDisplay::removeTool(const String& key) {
        auto wcptr = itsMWCTools.find(key);
		if(wcptr == itsMWCTools.end( )) return;
		std::shared_ptr<MultiWCTool> tool = wcptr->second;
		itsMWCTools.erase(wcptr);
		if ( tool.get( ) == 0 ) return;
		tool->removeWorldCanvases(this);
	}

	const std::shared_ptr<MultiWCTool> PanelDisplay::getTool(const String& key) {
        auto wcptr = itsMWCTools.find(key);
		if(wcptr == itsMWCTools.end( )) return std::shared_ptr<MultiWCTool>( );
		return wcptr->second;
	}

	float PanelDisplay::getDrawUnit(  ) const {
		return itsWCHList.size( ) > 0 ? itsWCHList.front( )->getDrawUnit() : 0.0;
	}

	int PanelDisplay::getColumnCount( ) const {
		return itsNX;
	}

	int PanelDisplay::getRowCount() const {
		return itsNY;
	}

} //# NAMESPACE CASA - END

