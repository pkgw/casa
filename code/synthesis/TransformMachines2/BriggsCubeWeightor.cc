//# BriggsCubeWeight.cc: Implementation for Briggs weighting for cubes
//# Copyright (C) 2018-2019
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 3 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
//# License for more details.
//#
//# https://www.gnu.org/licenses/
//#
//# You should have received a copy of the GNU  General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Queries concerning CASA should be submitted at
//#        https://help.nrao.edu
//#
//#        Postal address: CASA Project Manager 
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>
#include<msvis/MSVis/VisImagingWeight.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include<synthesis/TransformMachines2/FTMachine.h>
#include<synthesis/TransformMachines2/BriggsCubeWeightor.h>


namespace casa{//# CASA namespace
namespace refim {//# namespace refactor imaging
  
using namespace casacore;
using namespace casa;
using namespace casacore;
using namespace casa::refim;
using namespace casacore;
using namespace casa::vi;

  BriggsCubeWeightor::BriggsCubeWeightor(): grids_p(0), ft_p(0), f2_p(0), d2_p(0), uscale_p(0), vscale_p(0), uorigin_p(0),vorigin_p(0), nx_p(0), ny_p(0), rmode_p(""), noise_p(0.0), robust_p(2), superUniformBox_p(0), multiField_p(False),initialized_p(False), refFreq_p(-1.0), freqInterpMethod_p(InterpolateArray1D<Double, Complex>::nearestNeighbour) {
    multiFieldMap_p.clear();
    
    
  }
  
  BriggsCubeWeightor::BriggsCubeWeightor( const String& rmode, const Quantity& noise, const Double robust, const Int superUniformBox, const Bool multiField)  : grids_p(0), ft_p(0), f2_p(0), d2_p(0), uscale_p(0), vscale_p(0), uorigin_p(0),vorigin_p(0), nx_p(0), ny_p(0), initialized_p(False), refFreq_p(-1.0),freqInterpMethod_p(InterpolateArray1D<Double, Complex>::nearestNeighbour) {

    rmode_p=rmode;
    noise_p=noise;
    robust_p=robust;
    superUniformBox_p=superUniformBox;
    multiField_p=multiField;
    multiFieldMap_p.clear();
  }


  BriggsCubeWeightor::BriggsCubeWeightor(vi::VisibilityIterator2& vi,
				       const String& rmode, const Quantity& noise,
				       const Double robust,
                                         const ImageInterface<Complex>& templateimage, const RecordInterface& inrec,
					 const Int superUniformBox, const Bool multiField){
    rmode_p=rmode;
    noise_p=noise;
    robust_p=robust;
    superUniformBox_p=superUniformBox;
    multiField_p=multiField;
    initialized_p=False;
    refFreq_p=-1.0;
    
    init(vi, templateimage,inrec);




  }
				       
 void BriggsCubeWeightor::init(vi::VisibilityIterator2& vi,
			       const ImageInterface<Complex>& templateimage, const RecordInterface& inRec)
 {
  LogIO os(LogOrigin("BriggsCubeWeightor", "constructor", WHERE));

  //freqInterpMethod_p=interpMethod;
  //freqFrameValid_p=freqFrameValid;
  //chanchunk may call the same object
  if(initialized_p && nx_p==templateimage.shape()(0) && ny_p==templateimage.shape()(1)){
    CoordinateSystem cs=templateimage.coordinates();
    Double freq=cs.toWorld(IPosition(4,0,0,0,0))[3];
    if(freq==refFreq_p)
      return;
  }
  //cerr << "in bgwt init " << endl;
  //Need to save previous wieght scheme of vi
  visWgt_p=vi.getImagingWeightGenerator();
    VisImagingWeight vWghtNat("natural");
  vi.useImagingWeight(vWghtNat);
  vi::VisBuffer2 *vb=vi.getVisBuffer();
  Int nIndices=0;
  for (vi.originChunks();vi.moreChunks();vi.nextChunk()) {
    for (vi.origin(); vi.more(); vi.next()) {
      String key=String::toString(vb->msId())+"_"+String::toString(vb->fieldId()(0));
      Int index=0;
      if(multiField_p){
	//find how many indices will be needed
	index=multiFieldMap_p.size();
	if(multiFieldMap_p.count(key) < 1)
	  multiFieldMap_p[key]=index;
	nIndices=multiFieldMap_p.size();
      }
      else{
	multiFieldMap_p[key]=0;
	nIndices=1;
      }

    }
  }
  //cerr << "nindices " << nIndices << endl;
  vi.originChunks();
  vi.origin();
  String key=String::toString(vb->msId())+"_"+String::toString(vb->fieldId()(0));
  IPosition shp=templateimage.shape();
  nx_p=shp[0];
  ny_p=shp[1];
  CoordinateSystem cs=templateimage.coordinates();
  refFreq_p=cs.toWorld(IPosition(4,0,0,0,0))[3];
  Vector<String> units = cs.worldAxisUnits();
  units[0]="rad"; units[1]="rad";
  cs.setWorldAxisUnits(units);
  Vector<Double> incr=cs.increment();
  uscale_p=(nx_p*incr[0]);
  vscale_p=(ny_p*incr[1]);
  uorigin_p=nx_p/2;
  vorigin_p=ny_p/2;
  ////TESTOO
  //IPosition shp=templateimage.shape();
  shp[3]=shp[3]+4; //add two channel at begining and end;
  Vector<Double> refpix=cs.referencePixel();
  refpix[3]+=2;
  cs.setReferencePixel(refpix);
  TempImage<Complex> newTemplate(shp, cs);
  ///////////////////////
  //ImageInterface<Complex>& newTemplate=const_cast<ImageInterface<Complex>& >(templateimage);
  Vector<Matrix<Double> > sumWgts(nIndices);
  
  for(int index=0; index < nIndices; ++index){
    initializeFTMachine(index, newTemplate, inRec);
    Matrix<Float> dummy;
    
    ft_p[index]->initializeToSky(newTemplate, dummy, *vb);
    Vector<Double> convFunc(2+superUniformBox_p, 1.0);
    //cerr << "superuniform box " << superUniformBox_p << endl;
    ft_p[index]->modifyConvFunc(convFunc, superUniformBox_p, 1);
    for (vi.originChunks();vi.moreChunks();vi.nextChunk()) {
      for (vi.origin(); vi.more(); vi.next()) {
	
	key=String::toString(vb->msId())+"_"+String::toString(vb->fieldId()(0));
      
	//cerr << "key and index "<< key << "   " << index << "   " << multiFieldMap_p[key] << endl; 
	if(multiFieldMap_p[key]==index){
	  ft_p[index]->put(*vb, -1, true, FTMachine::PSF);
	}
      
      }
    }
    Array<Float> griddedWeight;
    ft_p[index]->getGrid(griddedWeight);
    //cerr << index << " griddedWeight Shape " << griddedWeight.shape() << endl;
    grids_p[index]->put(griddedWeight.reform(newTemplate.shape()));
    sumWgts[index]=ft_p[index]->getSumWeights();
    //cerr << "sumweight " << sumWgts[index] << endl;
    //clear the ftmachine
    ft_p[index]->finalizeToSky();
  }
  ////Lets reset vi before returning 
  vi.originChunks();
  vi.origin();
  
  
  Int nchan=newTemplate.shape()(3);
  for (uInt index=0; index< f2_p.nelements();++index){
    //cerr << "rmode " << rmode_p << endl;
    
    for (uInt chan=0; chan < uInt(nchan); ++ chan){
      IPosition start(4,0,0,0,chan);
      IPosition shape(4,nx_p,ny_p,1,1);
      Array<Float> arr;
      grids_p[index]->getSlice(arr, start, shape, True);
      Matrix<Float> gwt(arr);
      if (rmode_p=="norm" && (sumWgts[index](0,chan)> 0.0)) {
	//os << "Normal robustness, robust = " << robust << LogIO::POST;
	Double sumlocwt = 0.;
	for(Int vgrid=0;vgrid<ny_p;vgrid++) {
	  for(Int ugrid=0;ugrid<nx_p;ugrid++) {
	    if(gwt(ugrid, vgrid)>0.0) sumlocwt+=square(gwt(ugrid,vgrid));
	  }
	}
	f2_p[index][chan] = square(5.0*pow(10.0,Double(-robust_p))) / (sumlocwt / sumWgts[index](0,chan));
	d2_p[index][chan] = 1.0;

      }
      else if (rmode_p=="abs") {
	//os << "Absolute robustness, robust = " << robust << ", noise = "
	//   << noise.get("Jy").getValue() << "Jy" << LogIO::POST;
	f2_p[index][chan] = square(robust_p);
	d2_p[index][chan] = 2.0 * square(noise_p.get("Jy").getValue());
	
      }
      else {
	f2_p[index][chan] = 1.0;
	d2_p[index][chan] = 0.0;
      }
      
      
      
    }//chan
    
    
    
  }

  initialized_p=True;
 }

  void BriggsCubeWeightor::weightUniform(Matrix<Float>& imweight, const vi::VisBuffer2& vb){

    if(multiFieldMap_p.size()==0)
      throw(AipsError("BroggsCubeWeightor has not been initialized"));
    String key=String::toString(vb.msId())+"_"+String::toString(vb.fieldId()(0));
    Int index=multiFieldMap_p[key];
    Vector<Int> chanMap=ft_p[0]->channelMap(vb);
    //cerr << "weightuniform chanmap " << chanMap << endl;
    ///No matching channels
    if(max(chanMap)==-1)
      return;
    Int nvischan=vb.nChannels();
    Int nRow=vb.nRows();
    Matrix<Double> uvw=vb.uvw();
    imweight.resize(nvischan, nRow);
    imweight.set(0.0);
    
    Matrix<Float> weight;
    VisImagingWeight::unPolChanWeight(weight, vb.weightSpectrum());
    Matrix<Bool> flag;
    cube2Matrix(vb.flagCube(), flag);

    Int nChanWt=weight.shape()(0);
    Double sumwt=0.0;
    Float u, v;
    IPosition pos(4,0);
    for (Int row=0; row<nRow; row++) {
	for (Int chn=0; chn<nvischan; chn++) {
	  if ((!flag(chn,row)) && (chanMap(chn) > -1)) {
	    pos(3)=chanMap(chn);
	    Float f=vb.getFrequency(0,chn)/C::c;
	    u=-uvw(0, row)*f;
	    v=-uvw(1, row)*f;
	    Int ucell=Int(std::round(uscale_p*u+uorigin_p));
	    Int vcell=Int(std::round(vscale_p*v+vorigin_p));
	    pos(0)=ucell; pos(1)=vcell;
	    ////TESTOO
	    //if(row==0){
	     
	    //  ofstream myfile;
	    //  myfile.open ("briggsLoc.txt", ios::out | ios::app | ios::ate );
	    //  myfile << vb.rowIds()(0) << " uv " << uvw.column(0) << " loc " << pos[0] << ", " << pos[1] << "\n"<< endl;
	    //  myfile.close();
  

	    //}
	    //////
	    imweight(chn,row)=0.0;
	    if((ucell>0)&&(ucell<nx_p)&&(vcell>0)&&(vcell<ny_p)) {
	      Float gwt=grids_p[index]->getAt(pos);
	      if(gwt >0){
		imweight(chn,row)=weight(chn%nChanWt,row);
		imweight(chn,row)/=gwt*f2_p[index][pos[3]]+d2_p[index][pos[3]];
		sumwt+=imweight(chn,row);
	      }
	    }
	    //else {
	    // imweight(chn,row)=0.0;
	      //ndrop++;
	    //}
	  }
	  else{
	    imweight(chn,row)=0.0;
	  }
    
	}
    }

   
    if(visWgt_p.doFilter()){
      visWgt_p.filter (imweight, flag, uvw, vb.getFrequencies(0), imweight);

    }
    
  }

void BriggsCubeWeightor::initializeFTMachine(const uInt index, const ImageInterface<Complex>& templateimage, const RecordInterface& inRec){
  Int nchan=templateimage.shape()(3);
  if(ft_p.nelements() <= index){
    ft_p.resize(index+1);
    grids_p.resize(index+1);
    f2_p.resize(index+1);
    d2_p.resize(index+1);
    f2_p[index]=Vector<Float>(nchan, 0.0);
    d2_p[index]=Vector<Float>(nchan, 0.0);
  }
  ft_p[index]=new refim::GridFT(Long(1000000), Int(200), "BOX",1.0, true, false);
  Int tmpInt;
  inRec.get("freqinterpmethod", tmpInt);
  freqInterpMethod_p=static_cast<InterpolateArray1D<Double, Complex >::InterpolationMethod>(tmpInt);
  ft_p[index]->setFreqInterpolation(freqInterpMethod_p);
  inRec.get("freqframevalid", freqFrameValid_p);
  ft_p[index]->setFrameValidity(freqFrameValid_p);
  String error;
  if(!(ft_p[index]->recoverMovingSourceState(error, inRec)))
    throw(AipsError("BriggsCubeWeightor could not get the state of the ftmachine:" + error));
  //remember to make the stokes I
  grids_p[index]=new TempImage<Float>(templateimage.shape(), templateimage.coordinates(), 0.0);
  
}
void BriggsCubeWeightor::cube2Matrix(const Cube<Bool>& fcube, Matrix<Bool>& fMat)
  {
	  fMat.resize(fcube.shape()[1], fcube.shape()[2]);
	  Bool deleteIt1;
	  Bool deleteIt2;
	  const Bool * pcube = fcube.getStorage (deleteIt1);
	  Bool * pflags = fMat.getStorage (deleteIt2);
	  for (uInt row = 0; row < fcube.shape()[2]; row++) {
		  for (Int chn = 0; chn < fcube.shape()[1]; chn++) {
			  *pflags = *pcube++;
			  for (Int pol = 1; pol < fcube.shape()[0]; pol++, pcube++) {
				  *pflags = *pcube ? *pcube : *pflags;
			  }
			  pflags++;
		  }
	  }
	  fcube.freeStorage (pcube, deleteIt1);
	  fMat.putStorage (pflags, deleteIt2);
  }


  }//# namespace refim ends
}//namespace CASA ends

