/*
 * tBinaryDilation.cc
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
 *  test binaryDilation
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

    String imagename, constraintimagename;
    String outMaskName("testBinaryDilationOut.im");
    CoordinateSystem csys;
    IPosition shape;

    if (argc > 1) {
      imagename=argv[1];
      constraintimagename=argv[2];
    }
    else {

      imagename = "testBDilationIn.im";
      constraintimagename = "testBDilationInConstaint.im";
      shape=IPosition(4, 100, 100, 1, 5);
      csys=CoordinateUtil::defaultCoords4D();
       
      PagedImage<Float> InImage(TiledShape(shape),csys,imagename);
      PagedImage<Float> InConstraintImage(TiledShape(shape),csys,constraintimagename);
      // no mask case (all true)
      InConstraintImage.set(1);
      InImage.setUnits(Unit("Jy/pixel"));
      InImage.set(0.0);
      //masked regions
      InImage.putAt(1.0, IPosition(4,40,50,0,0));
      InImage.putAt(1.0, IPosition(4,45,55,0,0));
      InImage.putAt(1.0, IPosition(4,45,56,0,0));
      InImage.putAt(1.0, IPosition(4,0,0,0,0));
      InImage.putAt(1.0, IPosition(4,0,99,0,0));
      InImage.putAt(1.0, IPosition(4,99,99,0,0));
      InImage.putAt(1.0, IPosition(4,46,56,0,3));
      InImage.putAt(1.0, IPosition(4,45,56,0,4));
    }

    CountedPtr<ImageInterface<Float> > im;
    im = ImageUtilities::openImage<Float>(imagename);
    IPosition shp=im->shape();
    csys = im->coordinates();
    CountedPtr<ImageInterface<Float> > constraintim;
    constraintim = ImageUtilities::openImage<Float>(constraintimagename);
    IPosition constshp=im->shape();
    PagedImage<Float> outmaskimage(TiledShape(shp), csys, outMaskName);
    SDMaskHandler maskhandler;
    //Structure Element
    IPosition axislen(2, 3, 3);
    Array<Float> se(axislen);
    se.set(0);
    se(IPosition(2,1,0))=1.0;
    se(IPosition(2,0,1))=1.0;
    se(IPosition(2,1,1))=1.0;
    se(IPosition(2,2,1))=1.0;
    se(IPosition(2,1,2))=1.0;

    TempImage<Float> maskImage(TiledShape(constshp), csys);
    maskImage.copyData(*constraintim); 
    maskImage.attachMask( LatticeExpr<Bool>( iif(maskImage > 0, true, false)) );
    ArrayLattice<Bool> mask( maskImage.getMask());
    IPosition specshape(2,1,shp(3));
    Array<Bool> chanmask(specshape);
    chanmask.set(true);

    maskhandler.binaryDilation(*im, se, 75, mask, chanmask, outmaskimage);

  }
  catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
  }
  return 0;
}
