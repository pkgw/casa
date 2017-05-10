//# SDMaskHandler_GTest.cc: implementation of SDMaskHandler google test
//#
//# Copyright (C) 2015
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
//# 51 Franklin Street, Fifth FloorBoston, MA 02110-1335, USA 
//#

#include <synthesis/ImagerObjects/test/SDMaskHandler_GTest.h>
#include <synthesis/ImagerObjects/SDMaskHandler.h>
#include <images/Regions/WCBox.h>
#include <images/Regions/ImageRegion.h>

using namespace casacore;
using namespace casa;
using namespace std;

namespace test {

SDMaskHandlerTest::SDMaskHandlerTest() {};
SDMaskHandlerTest::~SDMaskHandlerTest() { }; 

void SDMaskHandlerTest::SetUp() {
    outMaskName="";
}

void SDMaskHandlerTest::TearDown() {
     
}

void SDMaskHandlerTest::generateBoxMaskImage(String imagename, Int imsize, Int nchan, IPosition blc, IPosition trc)
{
    //blc and trc --- 4 elements vect.
    IPosition shape(4, imsize, imsize, 4, nchan);
    csys=CoordinateUtil::defaultCoords4D();
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
    else {
     throw(AipsError("incompatible blc,trc, and/or nchan combination"));
    } 
}

ImageInterfaceTest::ImageInterfaceTest() {}
ImageInterfaceTest::~ImageInterfaceTest() {}


void ImageInterfaceTest::testMakeMaskByThreshold() 
{
    cout <<" Test makeMask()"<<endl;
    outMaskName="testMakeMask.im"; 
    Double thresval = 2.0;
    Quantity threshold(thresval, "Jy");
    IPosition shape(4, 100, 100, 1, 5);
    csys=CoordinateUtil::defaultCoords4D();
    //TempImage<Float> templateImage(TiledShape(shape),csys);
    PagedImage<Float> templateImage(TiledShape(shape),csys, String("mytempim"));
    templateImage.setUnits(Unit("Jy/pixel"));
    templateImage.set(0.0);
    for (uInt i=0; i < 5; i++) {
      for (uInt j=0; j < 3; j++ ) {
        //Int selpx = 45 + j;
        IPosition loc(4, 50, 45+j, 0, Int(i));
        IPosition loc2(4, 51, 45+j, 0, Int(i));
        IPosition loc3(4, 52, 45+j, 0, Int(i));
        Float val = 2.0 + Float(i);
        templateImage.putAt(val, loc);
        templateImage.putAt(val, loc2);
        templateImage.putAt(val, loc3);
      }
    }
    SHARED_PTR<ImageInterface<Float> > outmaskimage;
    SDMaskHandler maskhandler;
    outmaskimage = maskhandler.makeMask(outMaskName, threshold, templateImage);
    cerr<<"maskimage.shape()="<< outmaskimage->shape()<<endl;
    for (uInt i=0; i < 5; i++) {
      for (uInt j=0; j < 3; j++) {
        IPosition loc(4, 50, 45+j, 0, Int(i));
        IPosition loc2(4, 51, 45+j, 0, Int(i));
        IPosition loc3(4, 52, 45+j, 0, Int(i));
        // mask image pix values are 1 or 0
        //chan = 0 plane should be zero by threshold
        if (i==0) {
          ASSERT_TRUE(outmaskimage->getAt(loc)==Float(0.0));
          ASSERT_TRUE(outmaskimage->getAt(loc2)==Float(0.0));
          ASSERT_TRUE(outmaskimage->getAt(loc3)==Float(0.0));
        }
        else {
          ASSERT_TRUE(outmaskimage->getAt(loc)==Float(1.0));
          ASSERT_TRUE(outmaskimage->getAt(loc2)==Float(1.0));
          ASSERT_TRUE(outmaskimage->getAt(loc3)==Float(1.0));
       }
     }
   }
}//testMakeMaskByThreshold

void ImageInterfaceTest::testRegionToMaskImage()
{
  cout <<" Test regionToImageMask()"<<endl;
  //static Bool regionToImageMask(const String& maskimage, Record* regionRec, Matrix<Quantity> & blctrcs,
  //                  Matrix<Float>& circles, const Float& value=1.0);
  outMaskName="testRegionToImageMask.im";
  IPosition shape(4, 100, 100, 1, 10);
  csys=CoordinateUtil::defaultCoords4D();
  PagedImage<Float> maskImage(TiledShape(shape),csys, outMaskName);
  //SDMaskHandler::cloneImShape(templateImage, maskname);
  // region record
  Vector<Quantum<Double> > qblc(4);
  Vector<Quantum<Double> > qtrc(4);
  // blctrcs box
  Vector<Quantum<Double> > qblcbox(2);
  Vector<Quantum<Double> > qtrcbox(2);
  Vector<Double> pBlc(4);
  Vector<Double> pTrc(4);
  Vector<Double> wBlc(4);
  Vector<Double> wTrc(4);
  Vector<Double> pBlcBox(4);
  Vector<Double> pTrcBox(4);
  Vector<Double> wBlcBox(4);
  Vector<Double> wTrcBox(4);
  // for regionRecord blc=(25,65,0,2) trc=(42,75,0,9)
  pBlc(0)=25.0; pBlc(1)=65.0; pBlc(2)=0; pBlc(3)=2.0;
  pTrc(0)=42.0; pTrc(1)=75.0; pTrc(2)=0; pTrc(3)=9.0;
  // for blctrc: box with (50,10,0,0) (65,40,0,0) 
  pBlcBox(0)=50.0; pBlcBox(1)=10.0; pBlcBox(2)=0; pBlcBox(3)=0;
  pTrcBox(0)=65.0; pTrcBox(1)=40.0; pBlcBox(2)=0; pBlcBox(3)=0;
  csys.toWorld(wBlc,pBlc);
  csys.toWorld(wTrc,pTrc);
  csys.toWorld(wBlcBox,pBlcBox);
  csys.toWorld(wTrcBox,pTrcBox);
  for (Int i = 0; i< 4; i++) {
    Int iaxis = csys.pixelAxisToWorldAxis(i);
    qblc(i) = Quantum<Double>(wBlc(i), csys.worldAxisUnits()(iaxis));
    qtrc(i) = Quantum<Double>(wTrc(i), csys.worldAxisUnits()(iaxis));
    if (i<2) {
      qblcbox(i) = Quantum<Double>(wBlcBox(i), csys.worldAxisUnits()(iaxis));
      qtrcbox(i) = Quantum<Double>(wTrcBox(i), csys.worldAxisUnits()(iaxis));
    }
  }
  Vector<Int> absRel; // null = abosolute
  WCBox wbox(qblc,qtrc,csys,absRel);
  LCBox box(pBlc,pTrc,shape);
  Record rec=wbox.toRecord("");
  Record *regionRec=0;
  regionRec = new Record();
  regionRec->assign(rec);
  // a box by blc and trc
  Matrix<Quantity> blctrcs(1,4);
  blctrcs(0,0)=qblcbox(0);
  blctrcs(0,1)=qblcbox(1);
  blctrcs(0,2)=qtrcbox(0);
  blctrcs(0,3)=qtrcbox(1);
  // define two circles with r=5 and r=12
  Matrix<Float> circles(2,3);
  circles(0,0) = 5;
  circles(0,1) = 10;
  circles(0,2) = 8;
  circles(1,0) = 12;
  circles(1,1) = 70;
  circles(1,2) = 65;
  //SDMaskHandler::regionToImageMask(maskname, regionRec, blctrcs, circles);
  SDMaskHandler::regionToImageMask(maskImage, regionRec, blctrcs, circles);
  //spot checks...
  //    --- checks at boundaries
  //a box by regionRec
  //inside the mask 
  ASSERT_TRUE(maskImage.getAt(IPosition(4,25,65,0,2))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,42,75,0,9))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,35,70,0,5))==Float(1.0));
  //outside the mask 
  ASSERT_TRUE(maskImage.getAt(IPosition(4,24,65,0,5))==Float(0.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,42,76,0,5))==Float(0.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,25,65,0,0))==Float(0.0));
  //a box by blctrcs
  //inside the mask
  ASSERT_TRUE(maskImage.getAt(IPosition(4,50,10,0,0))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,65,40,0,9))==Float(1.0));
  //outside the mask
  ASSERT_TRUE(maskImage.getAt(IPosition(4,40,10,0,0))==Float(0.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,65,50,0,9))==Float(0.0));
  //circles
  ASSERT_TRUE(maskImage.getAt(IPosition(4,10,8,0,0))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,15,8,0,0))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,70,65,0,0))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,70,77,0,0))==Float(1.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,10,14,0,0))==Float(0.0));
  ASSERT_TRUE(maskImage.getAt(IPosition(4,83,65,0,0))==Float(0.0));

  delete regionRec;
}//testRegionToMaskImage

void ImageInterfaceTest::testRegionText()
{
  cout <<" Test regionTextToImageRegion()"<<endl;
  outMaskName="testRegionText.im";
  IPosition shape(4, 100, 100, 1, 10);
  csys=CoordinateUtil::defaultCoords4D();
  //TempImage<Float> templateImage(TiledShape(shape),csys);
  PagedImage<Float> regionImage(TiledShape(shape),csys, outMaskName);
  //String crtfFile="/home/casa-dev-08/ttsutsum/testcrtf.txt";
  String crtfFile="box [[45pix,50pix],[85pix,65pix]]";
  ImageRegion* imageRegion=0;
  SDMaskHandler::regionTextToImageRegion(crtfFile, regionImage, imageRegion);
  if (imageRegion!=0)
    SDMaskHandler::regionToMask(regionImage,*imageRegion, Float(1.0));
  delete imageRegion;
}

void ImageInterfaceTest::testCopyMask()
{
  cout << " Test copyMask() "<<endl;
  String inMaskImageName("inmasktemp.im");
  IPosition blc(4,35,40,0,10);
  IPosition trc(4,55,50,0,15);
  // 100x100x1x50 cube with box mask between ch10-15
  generateBoxMaskImage(inMaskImageName,100,50,blc,trc);
  PagedImage<Float> inmask(inMaskImageName);
  TempImage<Float> outmask(inmask.shape(),inmask.coordinates());
  //PagedImage<Float> outmask(inmask.shape(),inmask.coordinates(),"outmasknew.im");
  SDMaskHandler maskhandler = SDMaskHandler();
  maskhandler.copyMask(inmask,outmask);
  ASSERT_TRUE(outmask.getAt(IPosition(4,35,40,0,10))==Float(1.0));
  ASSERT_TRUE(outmask.getAt(IPosition(4,35,40,0,15))==Float(1.0));
  ASSERT_TRUE(outmask.getAt(IPosition(4,55,50,0,15))==Float(1.0));
  ASSERT_TRUE(outmask.getAt(IPosition(4,55,50,0,15))==Float(1.0));
  ASSERT_TRUE(outmask.getAt(IPosition(4,45,45,0,9))==Float(0.0));
  ASSERT_TRUE(outmask.getAt(IPosition(4,45,45,0,16))==Float(0.0));
}

void ImageInterfaceTest::testMaskByPerPlaneThreshold()
{
   cout <<" Test makeMaskByPerChanThreshold()"<<endl;
    outMaskName="testPerChanMask.im";
    Vector<Float> thresval(5);
    thresval(0) = 2.0;
    thresval(1) = 1.0;
    thresval(2) = 0.0;
    thresval(3) = 2.0;
    thresval(4) = 1.0;

    IPosition shape(4, 100, 100, 1, 5);
    csys=CoordinateUtil::defaultCoords4D();
    PagedImage<Float> templateImage(TiledShape(shape),csys, String("testMaskPerChanInput.im"));
    templateImage.setUnits(Unit("Jy/pixel"));
    templateImage.set(0.0);

    templateImage.putAt(2.5,  IPosition(4, 50, 45, 0, 0));
    templateImage.putAt(1.0,  IPosition(4, 50, 44, 0, 0));
    templateImage.putAt(1.0,  IPosition(4, 50, 46, 0, 0));
    templateImage.putAt(1.0,  IPosition(4, 49, 45, 0, 0));
    templateImage.putAt(1.0,  IPosition(4, 49, 44, 0, 0));
    templateImage.putAt(1.0,  IPosition(4, 25, 25, 0, 1));
    templateImage.putAt(1.0,  IPosition(4, 25, 25, 0, 2));
    templateImage.putAt(1.0,  IPosition(4, 25, 25, 0, 3));
    templateImage.putAt(3.0,  IPosition(4, 45, 25, 0, 3));
    templateImage.putAt(2.0,  IPosition(4, 45, 25, 0, 4));
    
    PagedImage<Float> outmaskimage(TiledShape(shape), csys, String("testMaskPerChanOutput.im"));
    SDMaskHandler maskhandler;
    maskhandler.makeMaskByPerChanThreshold(templateImage, outmaskimage, thresval); 
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,50,45,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,50,46,0,0))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,25,25,0,1))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,25,25,0,2))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,25,25,0,3))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,45,25,0,3))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,45,25,0,4))==Float(1.0));
}

void ImageInterfaceTest::testBinaryDilation()
{
    cout <<" Test binaryDilation()"<<endl;
    outMaskName="testBDilationOut.im";

    IPosition shape(4, 100, 100, 1, 5);
    csys=CoordinateUtil::defaultCoords4D();
    PagedImage<Float> InImage(TiledShape(shape),csys, String("testBDilationIn.im"));
    //PagedImage<Float> dummyMaskImage(TiledShape(shape),csys, String("testBDilationDummyMask.im"));
    TempImage<Float> dummyMaskImage(TiledShape(shape),csys);
    dummyMaskImage.set(0);
    // 1 to execlude from the binary dilation
    dummyMaskImage.putAt(1.0, IPosition(4,0,0,0,0));
    dummyMaskImage.putAt(1.0, IPosition(4,99,99,0,0));
    // 1 -> false
    dummyMaskImage.attachMask( LatticeExpr<Bool> (iif(dummyMaskImage >  0, False, True)));
    ArrayLattice<Bool> mask(dummyMaskImage.getMask());
    Array<Bool> chanmask(IPosition(2,1,5));
    chanmask.set(true);
    InImage.setUnits(Unit("Jy/pixel"));
    InImage.set(0.0);
    InImage.putAt(1.0, IPosition(4,40,50,0,0));
    InImage.putAt(1.0, IPosition(4,45,55,0,0));
    InImage.putAt(1.0, IPosition(4,45,56,0,0));
    InImage.putAt(1.0, IPosition(4,0,0,0,0));
    InImage.putAt(1.0, IPosition(4,0,99,0,0));
    InImage.putAt(1.0, IPosition(4,99,99,0,0));
    InImage.putAt(1.0, IPosition(4,46,56,0,3));
    InImage.putAt(1.0, IPosition(4,45,56,0,4));
    PagedImage<Float> outmaskimage(TiledShape(shape), csys, outMaskName);
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
    
    maskhandler.binaryDilation(InImage, se, 1, mask, chanmask, outmaskimage); 

    //value test
    //chan0 
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,0,0,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,51,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,41,50,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,39,50,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,49,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,50,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,52,0,0))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,0,98,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,1,98,0,0))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,99,99,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,99,98,0,0))==Float(0.0));
    //chan1
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,50,0,1))==Float(0.0));
    //chan3
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,46,57,0,3))==Float(1.0));
    //chan4
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,45,55,0,4))==Float(1.0));

}// testBinaryDilation

// ToDO: combine with testBinaryDilation()?
void ImageInterfaceTest::testBinaryDilationIter()
{
    cout <<" Test binaryDilationIter"<<endl;
    outMaskName="testBDilationIter2Out.im";

    IPosition shape(4, 100, 100, 1, 5);
    csys=CoordinateUtil::defaultCoords4D();
    PagedImage<Float> InImage(TiledShape(shape),csys, String("testBDilationIn.im"));
    //PagedImage<Float> dummyMaskImage(TiledShape(shape),csys, String("testBDilationDummyMask.im"));
    TempImage<Float> dummyMaskImage(TiledShape(shape),csys);
    // no mask case (all true)
    dummyMaskImage.set(1);
    dummyMaskImage.attachMask( LatticeExpr<Bool>( iif(dummyMaskImage > 0, true, false)) );
    ArrayLattice<Bool> mask( dummyMaskImage.getMask());
    Array<Bool> chanmask(IPosition(2,1,5));
    chanmask.set(true);
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
    dummyMaskImage.set(1.0);
    PagedImage<Float> outmaskimage(TiledShape(shape), csys, outMaskName);
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
    
    maskhandler.binaryDilation(InImage, se, 2, mask, chanmask, outmaskimage); 
    //chan0 
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,48,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,53,0,0))==Float(0.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,1,1,0,0))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,1,2,0,0))==Float(0.0));
    //chan1
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,40,50,0,1))==Float(0.0));
    //chan3
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,46,54,0,3))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,46,53,0,3))==Float(0.0));
    //chan4
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,43,56,0,4))==Float(1.0));
    ASSERT_TRUE(outmaskimage.getAt(IPosition(4,46,58,0,4))==Float(0.0));
} //testBinaryDilationIter

void ImageInterfaceTest::testYAPruneRegions()
{
    cout <<" Test YAPruneRegions()"<<endl;

    // 4 dim image with 5 spectral channels
    IPosition shape(4, 100, 100, 1, 5);
    // 4 dim image with 1 spectral channels
    //IPosition shape(4, 100, 100, 1, 1);
    csys=CoordinateUtil::defaultCoords4D();
    PagedImage<Float> InImage(TiledShape(shape),csys, String("testYAPruneRegions.im"));
    //PagedImage<Float> dummyMaskImage(TiledShape(shape),csys, String("testBDilationDummyMask.im"));
    TempImage<Float> dummyMaskImage(TiledShape(shape),csys);
    dummyMaskImage.set(0);
    //PagedImage<Float> labelMap(TiledShape(shape),csys, String("testYAPruneRegions-labelmap.im"));
    //labelMap.set(0);
    // 1 to execlude from the binary dilation
    dummyMaskImage.putAt(1.0, IPosition(4,0,0,0,0));
    dummyMaskImage.putAt(1.0, IPosition(4,99,99,0,0));
    // 1 -> false
    dummyMaskImage.attachMask( LatticeExpr<Bool> (iif(dummyMaskImage >  0, False, True)));
    ArrayLattice<Bool> mask(dummyMaskImage.getMask());
    Array<Bool> chanmask(IPosition(2,1,5));
    chanmask.set(true);
    InImage.setUnits(Unit("Jy/pixel"));
    InImage.set(0.0);
    InImage.putAt(1.0, IPosition(4,40,50,0,0));
    InImage.putAt(1.0, IPosition(4,44,54,0,0));
    InImage.putAt(1.0, IPosition(4,45,55,0,0));
    InImage.putAt(1.0, IPosition(4,45,56,0,0));
    InImage.putAt(1.0, IPosition(4,46,55,0,0));
    InImage.putAt(1.0, IPosition(4,47,55,0,0));
    InImage.putAt(1.0, IPosition(4,47,54,0,0));
    InImage.putAt(1.0, IPosition(4,0,0,0,0));
    InImage.putAt(1.0, IPosition(4,0,99,0,0));
    InImage.putAt(1.0, IPosition(4,99,99,0,0));
    InImage.putAt(1.0, IPosition(4,46,56,0,3));
    InImage.putAt(1.0, IPosition(4,45,56,0,4));
    SDMaskHandler maskhandler;

    Double prunesize=2.0;
    Vector<Bool> pruned; 
    SHARED_PTR<ImageInterface<Float> > tempIm_ptr = maskhandler.YAPruneRegions(InImage,pruned,prunesize);
    PagedImage<Float> outMask(InImage.shape(), InImage.coordinates(), "testYAPruneRegions-out.mask");
    outMask.copyData(*(tempIm_ptr.get()) );
}
void ImageInterfaceTest::testYAPruneRegionsBigImage()
{
    cout <<" Test YAPruneRegionsBigImage()"<<endl;

    // 4 dim image with 1 spectral channels
    IPosition shape(4, 1500,1500, 1, 1);
    // 4 dim image with 1 spectral channels
    //IPosition shape(4, 100, 100, 1, 1);
    csys=CoordinateUtil::defaultCoords4D();
    PagedImage<Float> InImage(TiledShape(shape),csys, String("testYAPruneRegions.im"));
    //PagedImage<Float> dummyMaskImage(TiledShape(shape),csys, String("testBDilationDummyMask.im"));
    TempImage<Float> dummyMaskImage(TiledShape(shape),csys);
    dummyMaskImage.set(0);
    //PagedImage<Float> labelMap(TiledShape(shape),csys, String("testYAPruneRegions-labelmap.im"));
    //labelMap.set(0);
    // 1 to execlude from the binary dilation
    //dummyMaskImage.putAt(1.0, IPosition(4,0,0,0,0));
    //dummyMaskImage.putAt(1.0, IPosition(4,99,99,0,0));
    // 1 -> false
    //dummyMaskImage.attachMask( LatticeExpr<Bool> (iif(dummyMaskImage >  0, False, True)));
    //ArrayLattice<Bool> mask(dummyMaskImage.getMask());
    //Array<Bool> chanmask(IPosition(2,1,5));
    //chanmask.set(true);
    InImage.setUnits(Unit("Jy/pixel"));
    InImage.set(0.0);
    InImage.putAt(1.0, IPosition(4,45,55,0,0));
    InImage.putAt(1.0, IPosition(4,45,56,0,0));
    InImage.putAt(1.0, IPosition(4,46,55,0,0));
    InImage.putAt(1.0, IPosition(4,47,55,0,0));
    InImage.putAt(1.0, IPosition(4,47,54,0,0));
    SDMaskHandler maskhandler;

    Double prunesize=2.0;
    Vector<Bool> pruned; 
    SHARED_PTR<ImageInterface<Float> > tempIm_ptr = maskhandler.YAPruneRegions(InImage,pruned,prunesize);
    PagedImage<Float> outMask(InImage.shape(), InImage.coordinates(), "testYAPruneRegions2-out.mask");
    outMask.copyData(*(tempIm_ptr.get()) );
}


//methods in SDMaskHandler.h but no tests exist here


//methods in SDMaskHandler.h but no tests exist here
//resetMask(SHARED_PTR<SIImageStore> imstore)
//fillMask((SHARED_PTR<SIImageStore> imstore, Vector<String> maskStrings)
//fillMask((SHARED_PTR<SIImageStore> imstore, String maskStrings)
//
//called from regionToImageMask: boxRegionToImageRegion(), circleRegionToImageRegion(), recordRegionToImageRegion()
//
//copyAllMasks() - calls copyMask()
//expandMask()
//InMaskToImageRegion()
//maskInteractiveMask()
//test automask...
//makeAutoMask() - old simplistic a box around the peak
//autoMask()
//autoMaskByThreshold()
//makePBMask()
//auotMaskWithinPB()
//checkMaskImage()
//cloneImShape()


// Tests
TEST_F(ImageInterfaceTest, testMakeMaskByThreshold) {
   testMakeMaskByThreshold();
}

TEST_F(ImageInterfaceTest, testRegionToMaskImage) {
   testRegionToMaskImage();
}

TEST_F(ImageInterfaceTest, testRegionText) {
  testRegionText();
}

TEST_F(ImageInterfaceTest, testCopyMask) {
  testCopyMask();
}

TEST_F(ImageInterfaceTest, testMaskByPerPlaneThreshold) {
  testMaskByPerPlaneThreshold();
}

TEST_F(ImageInterfaceTest, testBinaryDilation) {
  testBinaryDilation();
}

TEST_F(ImageInterfaceTest, testBinaryDilationIter) {
  testBinaryDilationIter();
}

TEST_F(ImageInterfaceTest, testYAPruneRegions) {
  testYAPruneRegions();
}

TEST_F(ImageInterfaceTest, testYAPruneRegionsBigImage) {
  testYAPruneRegionsBigImage();
}
}//test

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
