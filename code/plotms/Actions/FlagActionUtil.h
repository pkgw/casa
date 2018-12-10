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


#ifndef FLAGACTIONUTIL_H_
#define FLAGACTIONUTIL_H_

#include <vector>
#include <iostream>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Client/Client.h>

namespace casa {

// FlagActionUtil is a set of flagging methods that are
// accessed from multiple classes. All the methods are
// static and (hopefully) inlined.
class FlagActionUtil {
public:
	static PlotLogMessage* flagRange(Client *client, PlotMSPlot* plot,
	    int canvasIndex, std::vector<PlotRegion>& regions,
	    bool /*showUnflagged*/, bool /*showFlagged*/) {
	  PlotMSFlagging flagging = client->getFlagging();
	  std::cout << "Flagging the range " << std::endl;
	  PlotLogMessage* m = plot->flagRange(canvasIndex, flagging, Vector<PlotRegion>(regions), true);
	  return m;
	}

	static PlotLogMessage* unflagRange(Client *client, PlotMSPlot* plot,
      int canvasIndex, std::vector<PlotRegion>& regions,
      bool /*showUnflagged*/, bool /*showFlagged*/) {
	  PlotMSFlagging flagging = client->getFlagging();
	  PlotLogMessage* m = plot->flagRange(canvasIndex, flagging, Vector<PlotRegion>(regions), false);
	  return m;
	}

	// constructor
	FlagActionUtil();

	// destructor
	virtual ~FlagActionUtil();

protected:
  virtual void redrawPlots(Client *client, PlotMSPlot* plot, vector<PlotCanvasPtr>& visibleCanv  );
  virtual void addRedrawPlot( PlotMSPlot* plot );

private:
  // Keep list of plots that have to be redrawn.
  vector<PlotMSPlot*> flaggedPlots;
};

} /* namespace casa */
#endif /* FLAGACTIONUTIL_H_ */
