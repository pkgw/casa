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

#include <iostream>
#include <iterator>
#include <algorithm>
#include <utility>
#include <limits>

#include "ActionFlagAll.h"
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Client/Client.h>
#include <plotms/Actions/FlagActionUtil.h>

using namespace casacore;
namespace casa {

ActionFlagAll::ActionFlagAll( Client* client )
	: ActionTool( client ),
	  FlagActionUtil() {
	itsType_ = TOOL_FLAG_ALL;

	auto plots = client->getCurrentPlots();
  auto visibleCanv = client->currentCanvases();
	for (size_t i = 0; i < plots.size(); ++i) {
	  auto plot = plots[i];
	  if (plot == NULL) continue;

	  // save current iteration
	  auto const currentIter = plot->iter();

	  // examine number of unflagged points in each plots
	  std::vector<MaskedScatterPlotPtr> scatterPlots = plot->plots();
	  //std::cout << "Number of scatter plots: " << scatterPlots.size() << std::endl;
	  std::vector<std::pair<unsigned int, unsigned int> > numData;
	  std::transform(scatterPlots.begin(), scatterPlots.end(), std::back_inserter(numData),
	      [](MaskedScatterPlotPtr x) { auto m = x->maskedData();
	      return std::make_pair(m->size(), m->sizeUnmasked()); });

	  auto canvases = plot->canvases();
	  auto numCanvases = canvases.size();
	  //std::cout << "Number of canvases: " << numCanvases << " current iteration is " << currentIter << std::endl;
	  size_t indexOffset = currentIter;

	  for (size_t i = 0; i < numCanvases; ++i) {
	    auto canv = canvases[i];
	    auto scatterPlotId = i + indexOffset;

	    auto n = numData[scatterPlotId];
	    if (!canv->title().empty() && n.first != 0 && n.second == 0) {
	      //std::cout << "No unmasked data in canvas " << i << ". Canvas background must be changed." << std::endl;
	      canv->setAllFlagged();
	      canv->refresh();
	    }

      // Only apply to visible canvases.
      // TODO: need to examine
      bool visible = false;
      for(unsigned int k= 0; !visible && k < visibleCanv.size(); k++)
        if(canv == visibleCanv[k]) visible = true;
      if(!visible) {
//        std::cout << "canvas " << canv->title() << " is not visible. continue." << std::endl;
        continue;
      }

	  }
	}
}

bool ActionFlagAll::doTool(PlotMSApp* plotms) {

  if (!toolEnabled) {
    // Locate/Flag/Unflag on all visible canvases.
    const vector<PlotMSPlot*>& plots = plotms->getPlotManager().plots();
    vector<PlotCanvasPtr> visibleCanv = client->currentCanvases();

    PlotMSPlot* plot = nullptr;

    // flagging/unflagging operation
    for(unsigned int i = 0; i < plots.size(); i++) {
      plot = plots[i];
      if(plot == NULL) continue;

      // Get parameters.
      PlotMSPlotParameters& params = plot->parameters();

      // Detect if we are showing flagged/unflagged points (for locate)
      PMS_PP_Display* d = params.typedGroup<PMS_PP_Display>();
      Bool showUnflagged=(d->unflaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);
      Bool showFlagged=(d->flaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);

      vector<PlotCanvasPtr> canv = plot->canvases();
      for(unsigned int j = 0; j < canv.size(); j++) {

        // Only apply to visible canvases.
        bool visible = false;
        for(unsigned int k= 0; !visible && k < visibleCanv.size(); k++)
          if(canv[j] == visibleCanv[k]) visible = true;
        if(!visible) {
//          std::cout << "canvas " << j << " is not visible. continue." << std::endl;
          continue;
        }

        bool isCanvasMarkedForFlag = canv[j]->isMarkedForFlag();
        bool isCanvasMarkedForUnFlag = canv[j]->isMarkedForUnflag();


        vector<PlotRegion> regions(1);
        constexpr double kDoubleMin = std::numeric_limits<double>::lowest();
        constexpr double kDoubleMax = std::numeric_limits<double>::max();
        PlotCoordinate const upperLeft(kDoubleMin, kDoubleMax);
        PlotCoordinate const lowerRight(kDoubleMax, kDoubleMin);
        regions[0] = PlotRegion(upperLeft, lowerRight);
        if (isCanvasMarkedForFlag) {
          // flag plotted data
          //std::cout << "canvas " << j << " is marked for flag." << std::endl;
          //std::cout << "regions[0]: " << regions[0].bottom() << ", " << regions[0].left()
          //    << ", " << regions[0].top() << ", " << regions[0].right() << std::endl;
          FlagActionUtil::flagRange(client, plot, j, regions, showUnflagged, showFlagged);

          // If this plot was flagged/unflagged, add it to the redraw list.
          addRedrawPlot( plot );
        } else if (isCanvasMarkedForUnFlag) {
          // unflag plotted data
          //std::cout << "canvas " << j << " is marked for unflag." << std::endl;
          //std::cout << "regions[0]: " << regions[0].bottom() << ", " << regions[0].left()
          //    << ", " << regions[0].top() << ", " << regions[0].right() << std::endl;
          FlagActionUtil::unflagRange(client, plot, j, regions, showUnflagged, showFlagged);

          // If this plot was flagged/unflagged, add it to the redraw list.
          addRedrawPlot( plot );
        }
      }
    }

    // reset the background color for each visible canvases
    for(unsigned int i = 0; i < plots.size(); i++) {
      plot = plots[i];
      if(plot == NULL) continue;

      vector<PlotCanvasPtr> canv = plot->canvases();
      for(unsigned int j = 0; j < canv.size(); j++) {
        // reset background color
        auto const currentBackground = canv[j]->background();
        if (canv[j]->isBackgroundColorChanged()) {
          canv[j]->setBackground(canv[j]->defaultBackground());
        }

        // clear all marks
        canv[j]->clearMark();

        // Only apply to visible canvases.
        bool visible = false;
        for(unsigned int k= 0; !visible && k < visibleCanv.size(); k++)
          if(canv[j] == visibleCanv[k]) visible = true;
        if(!visible) {
//          std::cout << "canvas " << j << " is not visible. continue." << std::endl;
          continue;
        }
        if (visible) {
          canv[j]->refresh();
        }
      }
    }
    redrawPlots(client, plot, visibleCanv);
  }

  return true;
}

ToolCode ActionFlagAll::getToolCode() const {
	return FLAGALL_TOOL;
}



ActionFlagAll::~ActionFlagAll() {
}

using namespace casacore;
} /* namespace casa */
