//# tImager.cc:  this tests Imager
//# Copyright (C) 1996,1997,1999,2001
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

#include <casa/iostream.h>
#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casa/Containers/Block.h>
#include <measures/Measures/MRadialVelocity.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <casa/Logging/LogIO.h>
#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/ComponentShape.h>
#include <components/ComponentModels/Flux.h>
#include <synthesis/MeasurementEquations/Imager.h>
#include <synthesis/MeasurementEquations/ImagerMultiMS.h>
#include <synthesis/TransformMachines2/test/MakeMS.h>
#include <casa/namespace.h>

int main()
{
  using namespace std;
using namespace casacore;
  using namespace casa;
  try{
    
    
    MDirection thedir(Quantity(20.0, "deg"), Quantity(20.0, "deg"));
    String msname("Test.ms");
    casa::test::MakeMS::makems(msname, thedir);
    //Bool compress = false;
    Bool useModel = True;
    MeasurementSet myms(msname, Table::Update);
    ImagerMultiMS* imgr = new ImagerMultiMS( );
    cout <<"--Imager created for MeasurementSet object. " << endl;

    LogSink logSink_p=LogSink(LogMessage::NORMAL, false);	  
    logSink_p.clearLocally();
    LogIO os(LogOrigin("tImager", "main()", WHERE), logSink_p);
    MRadialVelocity dummy;
    MFrequency dummy1;
    imgr->setDataPerMS(msname, "none", Vector<Int>(1,100), Vector<Int>(1,0), Vector<Int>(1,1), Vector<Int>(1,0), Vector<Int>(1,0),"", "","", Vector<Int>(0), "", "", "", "", "", "", useModel);
    imgr->defineImage(600, 600, Quantity(1.1, "arcsec"), Quantity(1.1, "arcsec"), "I",
		      thedir, 0, "channel", 1, 40, 1, dummy1, dummy, Quantity(1.0, "Hz"), Vector<Int>(1,0));
    imgr->setvp(True,True,"",False, Quantity(360.0, "deg") ,Quantity(180.0,"deg"));
    //imgr->summary();
    imgr->makeimage("pb", "pb.image");
    {
      PagedImage<Float> pi("pb.image");
      pi.set(0.0);
    }
    ////Let's make a component now...bang on a pixel centre
     MDirection thedir1(Quantity(19.949941666666668, "deg"), Quantity(19.95019027777778, "deg"));
    {
      ComponentList cl;
      SkyComponent otherPoint(ComponentType::POINT);
      otherPoint.flux() = Flux<Double>(1.0);
      otherPoint.shape().setRefDirection(thedir1);
      cl.add(otherPoint);
      cl.rename(String("p.cl"));
    }
    imgr->pb("pb.image", "applypb.image", "p.cl", "", "apply", thedir, Quantity(360.0, "deg"),"pb");
    imgr->ft(Vector<String>(0), "p.cl");
    imgr->makeimage("model", "model.image");
    PagedImage<Float> mod("model.image");
    DirectionCoordinate dc=mod.coordinates().directionCoordinate();
    Vector<Double> dirpix=dc.toPixel(thedir1);
    IPosition pixat(4, Int(std::round(dirpix[0])), Int(std::round(dirpix[1])), 0,0);
    //cerr <<"dirpix " << dirpix << "   " << pixat << endl;
    Float modelval=mod.getAt(pixat);
    PagedImage<Float> compim("applypb.image");
    Float compval=compim.getAt(pixat);
    //cerr << "Values " << modelval << "  " << compval << endl;
    AlwaysAssertExit(near(modelval, compval, 1.0e-2));
    cerr <<"OK" << endl;
    exit(0);
  }catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
    return(1);
  }
};
