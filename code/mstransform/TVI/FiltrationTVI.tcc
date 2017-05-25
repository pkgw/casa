//# FilteringTVI.tcc: Template class for data filtering TVI
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

#include <mstransform/TVI/FilteringTVI.h>

#include <casacore/casa/Exceptions/Error.h>

#include <msvis/MSVis/VisibilityIterator2.h>

using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE vi - BEGIN

// forward declaration

// FilteringTVI implementation
template<class Filter>
FilteringTVI<Filter>::FilteringTVI(ViImplementation2 * inputVi, Filter *filter) :
    TransformingVi2(inputVi), filter_p(filter) {
}

template<class Filter>
FilteringTVI<Filter>::~FilteringTVI() {
}

template<class Filter>
void FilteringTVI<Filter>::origin() {
  auto const vii = getVii();
  vii->origin();

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // filtering
  filter();
}

template<class Filter>
void FilteringTVI<Filter>::next() {
  // filtering
  filter();
}

template<class Filter>
void FilteringTVI<Filter>::filter() {
  auto const vii = getVii();

  for (; vii->more() && filter_p->isResidue(vii->getVisBuffer()); vii->next()) {
    // Synchronize own VisBuffer
    configureNewSubchunk();
  }
}

} //# NAMESPACE vi - END

} //# NAMESPACE CASA - END
