//# Copyright (C) 2009
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

#include "FlagActionUtil.h"
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Client/Client.h>

using namespace casacore;
namespace casa {

FlagActionUtil::FlagActionUtil()
  : flaggedPlots(){
}

FlagActionUtil::~FlagActionUtil() {
}

void FlagActionUtil::addRedrawPlot( PlotMSPlot* plot ){
	flaggedPlots.push_back( plot );
}

void FlagActionUtil::redrawPlots(Client *client, PlotMSPlot* plot, vector<PlotCanvasPtr>& visibleCanv  ){
	// For a flag/unflag, need to tell the plots to redraw themselves,
	// and clear selected regions.
	bool hold = client->allDrawingHeld();
	if(!hold) client->holdDrawing();

	for(unsigned int i = 0; i < flaggedPlots.size(); i++) {
		flaggedPlots[i]->plotDataChanged();

		vector<PlotCanvasPtr> canv = plot->canvases();
		for(unsigned int j = 0; j < canv.size(); j++) {
			// Only apply to visible canvases.
			bool visible = false;
			for(unsigned int k = 0;
					!visible && k < visibleCanv.size(); k++)
				if(canv[j] == visibleCanv[k]) visible = true;
			if(!visible) continue;

			canv[j]->clearSelectedRects();
		}
	}

	if(!hold) client->releaseDrawing();
	flaggedPlots.clear();
}

} /* namespace casa */
