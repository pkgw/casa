//# FilteringTVI.h: Template class for data filtering TVI
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
#include <msvis/MSVis/ViiLayerFactory.h>

using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE vi - BEGIN

// forward declaration
class SDDoubleCircleFilter;

// Factory
FilteringTVIFactory::FilteringTVIFactory(casacore::Record const &configuration,
    ViImplementation2 *inputVII) :
    inputVII_p(inputVII), configuration_p(configuration) {
}

FilteringTVIFactory::~FilteringTVIFactory() {
}

ViImplementation2 * FilteringTVIFactory::createVi() const {
  ViImplementation2 *vii = nullptr;

  if (!configuration_p.isDefined("type")) {
    throw AipsError(String("You have to specify Filtering type"));
  }

  Int const type_indicator = configuration_p.asInt("type");

  switch(type_indicator) {
  case FilteringType::SDDoubleCircleFilter:
    // new filter impl
    // filter impl requires MS
    MeasurementSet const &ms = inputVII_p->ms();
    SDDoubleCircleFilter filter = new SDDoubleCircleFilter(ms, configuration_p);

    // new filter
    vii = new FilteringTVI<SDDoubleCircleFilter>(inputVII_p, filter);
    break;
  default:
    // unsupported type, throw exception
    throw AipsError(String("Invalid FilteringType for FilteringTVI: ")+String::toString(type_indicator));
    break;
  }

  return vii;
}

// LayerFactory
FilteringTVILayerFactory::FilteringTVILayerFactory(Record const &configuration) :
    ViiLayerFactory(), configuration_p(configuration) {
}

ViImplementation2 * FilteringTVILayerFactory::createInstance(
    ViImplementation2* vii0) const {
  // Make the FilteringTVI, using supplied ViImplementation2, and return it
  FilteringTVIFactory factory(configuration_p, vii0);
  return factory.createVi();
}

} //# NAMESPACE vi - END

} //# NAMESPACE CASA - END
