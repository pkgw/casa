/*
 * tSynthesisUtils.cc
 *demo of SynthesisUtil functionality
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
 *  test labelRegions and findBlobSize
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


    String imagename;
    if (argc > 1) {
      imagename=argv[1];
    }
    else {
      imagename="tLabelandFindRegions.temp.mask";
      Int imsize = 100;
      Int nchan = 2; 
      IPosition blc(4,45,50,0,0); 
      IPosition trc(4,55,55,0,0); 
      IPosition shape(4, imsize, imsize, 1, nchan);
      CoordinateSystem csys=CoordinateUtil::defaultCoords4D();
      PagedImage<Float> maskImage(TiledShape(shape), csys, imagename);
      maskImage.setUnits(Unit("Jy/pixel"));
      maskImage.set(0.0);
      // sanity check
      if (blc(0) <= trc(0) && blc(1) <= trc(1) && blc(3) < nchan) {
        Int dx = trc(0) - blc(0);
        Int dy = trc(1) - blc(1);
        for (uInt i=0+blc(3); i < trc(3)+1; i++) {
          for (uInt j=0; j < (uInt)dx+1; j++) {
            for (uInt k=0; k < (uInt)dy+1; k++) {
              IPosition loc(4,blc(0)+j, blc(1)+k, 0, Int(i));
              Float val = 1.0;
              maskImage.putAt(val, loc);
            }
          }
        }
      }
    }

    CountedPtr<ImageInterface<Float> > im;
    im = ImageUtilities::openImage<Float>(imagename);
    IPosition shp=im->shape();
    cerr<<"shp="<<shp<<endl;
    IPosition start(4, 0, 0, 0, 0);
    IPosition length(4, shp(0),shp(1),shp(2),1);
    Slicer sl(start, length);
    AxesSpecifier aspec(False);
    SubImage<Float>* subIm = new SubImage<Float>(*im, sl, aspec, True);
    cerr<<"subIm shape="<<subIm->shape()<<" slice="<<sl<<endl;
    SDMaskHandler smh;

    TempImage<Float>* tempIm = new TempImage<Float> (TiledShape(subIm->shape(), subIm->niceCursorShape()), subIm->coordinates() );
    tempIm->copyData(LatticeExpr<Float> (*subIm));
    cerr<<"tempIm shape="<<tempIm->shape()<<endl;

    cerr<<"setup a blobmap ..."<<endl;
    TempImage<Float>* blobMap = new TempImage<Float> (TiledShape(subIm->shape(), subIm->niceCursorShape()), subIm->coordinates() );
    blobMap->set(0);
    cerr<<"blobMap shape="<<blobMap->shape()<<endl;

    cerr<<"Calling labelRegions()..."<<endl;
    //Array<Float> subimarr;
    //Array<Float> blobarr;
    //smh.getarr(*subIm, False);
    //smh.getarr(*blobMap, True);
    smh.labelRegions(*subIm, *blobMap);
    cerr<<"Calling findBlobSize()..."<<endl;
    Vector<Float> bloblist=smh.findBlobSize(*blobMap);
    for (uInt iblob=0; iblob<bloblist.nelements(); ++iblob) {
        cerr<<iblob<<" : "<<bloblist[iblob]<<endl;
    } 

  }catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
  }
  return 0;
};
