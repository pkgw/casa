/*
 * tSkipChannel.cc
//# Copyright (C) 2013-2014
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
 *
 *  test skioChannel
 */

#include <casacore/images/Images/ImageUtilities.h>
#include <casa/iostream.h>
#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Block.h>
#include <measures/Measures/MRadialVelocity.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <casa/Logging/LogIO.h>
#include <synthesis/ImagerObjects/SynthesisImager.h>
#include <synthesis/ImagerObjects/SIImageStore.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>
#include <lattices/Lattices/LatticeConcat.h>
#include <images/Images/PagedImage.h>
#include <images/Images/ImageConcat.h>
#include <images/Images/SubImage.h>
#include <casa/namespace.h>
#include <images/Images/TempImage.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <ms/MSSel/MSSourceIndex.h>
#include <synthesis/ImagerObjects/SDMaskHandler.h>

int main(int argc, char **argv)
{
  using namespace std;
  using namespace casacore;
  using namespace casa;
  try{

    Timer timer;

    String curimname, previmname;
    String outMaskName("testSkipChannels.im");
    CoordinateSystem csys;
    IPosition shape;

    if (argc > 1) {
      curimname=argv[1];
      previmname=argv[2];
    }
    else {

      curimname = "testSkipChanCur.im";
      previmname = "testSkipChanPrev.im";
      shape=IPosition(4, 100, 100, 1, 5);
      csys=CoordinateUtil::defaultCoords4D();
       
      PagedImage<Float> curImage(TiledShape(shape),csys,curimname);
      PagedImage<Float> prevImage(TiledShape(shape),csys,previmname);
      // no mask case (all true)
      prevImage.set(1);
      curImage.setUnits(Unit("Jy/pixel"));
      curImage.set(0.0);
      //masked regions
      curImage.putAt(1.0, IPosition(4,40,50,0,0));
      curImage.putAt(1.0, IPosition(4,45,55,0,0));
      curImage.putAt(1.0, IPosition(4,45,56,0,0));
      curImage.putAt(1.0, IPosition(4,0,0,0,0));
      curImage.putAt(1.0, IPosition(4,0,99,0,0));
      curImage.putAt(1.0, IPosition(4,99,99,0,0));
      curImage.putAt(1.0, IPosition(4,46,56,0,3));
      curImage.putAt(1.0, IPosition(4,45,56,0,4));
    }

    Float fracChange = 0.0;
    CountedPtr<ImageInterface<Float> > im;
    im = ImageUtilities::openImage<Float>(curimname);
    IPosition shp=im->shape();
    Int nchan = shp(3);
    csys = im->coordinates();
    CountedPtr<ImageInterface<Float> > previm;
    previm = ImageUtilities::openImage<Float>(previmname);
    IPosition constshp=im->shape();
    SDMaskHandler maskhandler;
    Vector<Bool> chanflag(nchan);
    Vector<Bool> zeromask(nchan);
    Matrix<String> type(1,nchan);
    Bool isthresholdreached(false);
    type.set("noise");
    chanflag.set(false);
    maskhandler.skipChannels(fracChange, *previm, *im, type, isthresholdreached, chanflag, zeromask);
  }
  catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
  }
  return 0;
}
