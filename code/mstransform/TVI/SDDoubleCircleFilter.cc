//# SDDoubleCircleFilter.h: Filter for SDDoubleCircleGainCal
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the Implied warranty of MERCHANTABILITY or
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

#include <mstransform/TVI/SDDoubleCircleFilter.h>

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE vi - BEGIN

// constructor
SDDoubleCircleFilter::SDDoubleCircleFilter(casacore::MeasurementSet const &ms,
    casacore::Record const &configuration)
: ms_p(ms), configuration_p(configuration) {
  initFilter();
}

// destructor
//SDDoubleCircleFilter::~SDDoubleCircleFilter() {
//
//}

// isFiltrate returns true if given vb does pass through the filter
bool SDDoubleCircleFilter::isFiltrate(VisBuffer2 const *vb) {
  return true;
}

void SDDoubleCircleFilter::initFilter() {

}

} //# NAMESPACE vi - END

} //# NAMESPACE CASA - END

#endif  // _MSVIS_SD_DOUBLE_CIRCLE_FILTER_H_
