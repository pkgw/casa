//# Copyright (C) 1998,1999,2000,2001
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

#include <imageanalysis/Images/test/ComponentListImageTest.h>

#include <components/ComponentModels/ConstantSpectrum.h>
#include <components/ComponentModels/Flux.h>
#include <components/ComponentModels/GaussianShape.h>

using namespace casacore;

using namespace std;

using namespace casa;

namespace test {

ComponentListImageTest::ComponentListImageTest() {
}

ComponentListImageTest::~ComponentListImageTest() {}

void ComponentListImageTest::SetUp() {}

void ComponentListImageTest::TearDown() {}

ComponentList ComponentListImageTest::oneGaussianCL() const {
    MDirection dir(Quantity(0, "deg"), Quantity(0, "deg"), MDirection::J2000);
    Quantity majorAxis(10, "arcmin");
    Quantity minorAxis(8, "arcmin");
    Quantity pa(45, "deg");
    GaussianShape g(dir, majorAxis, minorAxis, pa);
    Flux<Double> flux;
    ConstantSpectrum cs;
    SkyComponent sc(flux, g, cs);
    ComponentList cl;
    cl.add(sc);
    return cl;
}


}
