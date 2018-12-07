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

#include "ActionFlagAll.h"
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Client/Client.h>

using namespace casacore;
namespace casa {

ActionFlagAll::ActionFlagAll( Client* client )
	: ActionTool( client ){
  std::cout << "ActionFlagAll instance is created" << std::endl;
	itsType_ = TOOL_FLAG_ALL;
}

bool ActionFlagAll::doTool(PlotMSApp* plotms) {
  std::cout << "ActionFlagAll::doTool" << std::endl;
  std::cout << "toolEnabled = " << toolEnabled << std::endl;

  if (!toolEnabled) {
    // Locate/Flag/Unflag on all visible canvases.
    const vector<PlotMSPlot*>& plots = plotms->getPlotManager().plots();
    vector<PlotCanvasPtr> visibleCanv = client->currentCanvases();

    // Get flagging parameters.
    //PlotMSFlagging flagging = client->getFlagging();

    // default background color
    // TODO: find a way to obtain default (or current) background color
    PlotAreaFillPtr defaultBackground = NULL;

    PlotMSPlot* plot;
    for(unsigned int i = 0; i < plots.size(); i++) {
      plot = plots[i];
      if(plot == NULL) continue;

      // Get parameters.
      PlotMSPlotParameters& params = plot->parameters();
      PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();

      // Detect if we are showing flagged/unflagged points (for locate)
      PMS_PP_Display* d = params.typedGroup<PMS_PP_Display>();
      Bool showUnflagged=(d->unflaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);
      Bool showFlagged=(d->flaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);

      vector<PlotCanvasPtr> canv = plot->canvases();
      for(unsigned int j = 0; j < canv.size(); j++) {

        // Only apply to visible canvases.
        // TODO: need to examine
        bool visible = false;
        for(unsigned int k= 0; !visible && k < visibleCanv.size(); k++)
          if(canv[j] == visibleCanv[k]) visible = true;
        if(!visible) {
          std::cout << "canvas " << j << " is not visible. continue." << std::endl;
          continue;
        }

        bool isCanvasMarkedForFlag = canv[j]->isMarkedForFlag();
        bool isCanvasMarkedForUnFlag = canv[j]->isMarkedForUnflag();


        if (isCanvasMarkedForFlag) {
          std::cout << "canvas " << j << " is marked for flag." << std::endl;
        } else if (isCanvasMarkedForUnFlag) {
          std::cout << "canvas " << j << " is marked for unflag." << std::endl;
        } else {
          // get default background
          // TODO: exclude all-flagged panels from the beginning
          if (defaultBackground.null()) {
            defaultBackground = canv[j]->background();
          }
        }
      }
    }

    for(unsigned int i = 0; i < plots.size(); i++) {
      plot = plots[i];
      if(plot == NULL) continue;

      // Get parameters.
      PlotMSPlotParameters& params = plot->parameters();
      PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();

      // Detect if we are showing flagged/unflagged points (for locate)
      PMS_PP_Display* d = params.typedGroup<PMS_PP_Display>();
      Bool showUnflagged=(d->unflaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);
      Bool showFlagged=(d->flaggedSymbol()->symbol()!=PlotSymbol::NOSYMBOL);

      vector<PlotCanvasPtr> canv = plot->canvases();
      for(unsigned int j = 0; j < canv.size(); j++) {
        // reset background color
        auto const currentBackground = canv[j]->background();
        std::cout << "BG color " << currentBackground->color() << std::endl;
        if (canv[j]->isBackgroundColorChanged()) {
          std::cout << "BG color changed" << std::endl;
          if (!defaultBackground.null()) {
            canv[j]->setBackground(defaultBackground);
          }
        }

        // clear all marks
        canv[j]->clearMark();

        // Only apply to visible canvases.
        bool visible = false;
        for(unsigned int k= 0; !visible && k < visibleCanv.size(); k++)
          if(canv[j] == visibleCanv[k]) visible = true;
        if(!visible) {
          std::cout << "canvas " << j << " is not visible. continue." << std::endl;
          continue;
        }
        if (visible) {
          canv[j]->refresh();
        }
      }
    }
  }

  return true;
}

ToolCode ActionFlagAll::getToolCode() const {
  std::cout << "return FLAGALL_TOOL" << std::endl;
	return FLAGALL_TOOL;
}



ActionFlagAll::~ActionFlagAll() {
}

using namespace casacore;
} /* namespace casa */
