//# FTMachine.cc: Implementation of FTMachine class
//# Copyright (C) 1997-2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU General Public License
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
#include <cmath>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/UnitMap.h>
#include <casa/Quanta/UnitVal.h>
#include <measures/Measures/Stokes.h>
#include <casa/Quanta/Euler.h>
#include <casa/Quanta/RotMatrix.h>
#include <measures/Measures/MFrequency.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/Projection.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <casa/BasicSL/Constants.h>
#include <synthesis/TransformMachines2/FTMachine.h>
#include <synthesis/TransformMachines2/SkyJones.h>
#include <synthesis/TransformMachines2/VisModelData.h>
#include <scimath/Mathematics/RigidVector.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/StokesVector.h>
#include <msvis/MSVis/MSUtil.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/ImageUtilities.h>
#include <casa/Containers/Block.h>
#include <casa/Containers/Record.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/MatrixIter.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/BinarySearch.h>
#include <casa/Exceptions/Error.h>
#include <scimath/Mathematics/NNGridder.h>
#include <scimath/Mathematics/ConvolveGridder.h>
#include <measures/Measures/UVWMachine.h>

#include <casa/System/ProgressMeter.h>

#include <casa/OS/Timer.h>
#include <casa/sstream.h>
#include <casa/iostream.h>

using namespace casacore;
namespace casa{//# CASA namespace
namespace refim {//# namespace refactor imaging
  
using namespace casacore;
using namespace casa;
using namespace casacore;
using namespace casa::refim;
using namespace casacore;
using namespace casa::vi;
  FTMachine::FTMachine() : isDryRun(false), image(0), uvwMachine_p(0), 
			   tangentSpecified_p(false), fixMovingSource_p(false),
			   distance_p(0.0), lastFieldId_p(-1),lastMSId_p(-1), 
			   useDoubleGrid_p(false), 
			   freqFrameValid_p(false), 
			   freqInterpMethod_p(InterpolateArray1D<Double,Complex>::nearestNeighbour), 
			   pointingDirCol_p("DIRECTION"),
			   cfStokes_p(), cfCache_p(), cfs_p(), cfwts_p(), cfs2_p(), cfwts2_p(), 
			   canComputeResiduals_p(false), toVis_p(true), 
                           numthreads_p(-1), pbLimit_p(0.05),sj_p(0), cmplxImage_p( )
  {
    spectralCoord_p=SpectralCoordinate();
    isPseudoI_p=false;
    spwChanSelFlag_p=0;
    polInUse_p=0;
    pop_p = new PolOuterProduct;
  }
  
  FTMachine::FTMachine(CountedPtr<CFCache>& cfcache,CountedPtr<ConvolutionFunction>& cf):
    isDryRun(false), image(0), uvwMachine_p(0), 
    tangentSpecified_p(false), fixMovingSource_p(false),
    distance_p(0.0), lastFieldId_p(-1),lastMSId_p(-1), 
    useDoubleGrid_p(false), 
    freqFrameValid_p(false), 
    freqInterpMethod_p(InterpolateArray1D<Double,Complex>::nearestNeighbour), 
    pointingDirCol_p("DIRECTION"),
    cfStokes_p(), cfCache_p(cfcache), cfs_p(), cfwts_p(), cfs2_p(), cfwts2_p(),
    convFuncCtor_p(cf),canComputeResiduals_p(false), toVis_p(true), numthreads_p(-1), 
    pbLimit_p(0.05),sj_p(0), cmplxImage_p( )
  {
    spectralCoord_p=SpectralCoordinate();
    isPseudoI_p=false;
    spwChanSelFlag_p=0;
    polInUse_p=0;
    pop_p = new PolOuterProduct;
  }
  
  LogIO& FTMachine::logIO() {return logIO_p;};
  
  //---------------------------------------------------------------------- 
  FTMachine& FTMachine::operator=(const FTMachine& other)
  {
    if(this!=&other) {
      image=other.image;
      //generic selection stuff and state
      nAntenna_p=other.nAntenna_p;
      distance_p=other.distance_p;
      lastFieldId_p=other.lastFieldId_p;
      lastMSId_p=other.lastMSId_p;
      
      tangentSpecified_p=other.tangentSpecified_p;
      mTangent_p=other.mTangent_p;
      mImage_p=other.mImage_p;
      mFrame_p=other.mFrame_p;
      
      nx=other.nx;
      ny=other.ny;
      npol=other.npol;
      nchan=other.nchan;
      nvischan=other.nvischan;
      nvispol=other.nvispol;
      mLocation_p=other.mLocation_p;
      if(uvwMachine_p)
    	  delete uvwMachine_p;
      if(other.uvwMachine_p)
	uvwMachine_p=new casacore::UVWMachine(*other.uvwMachine_p);
      else
	uvwMachine_p=0;
      doUVWRotation_p=other.doUVWRotation_p;
      //Spectral and pol stuff 
      freqInterpMethod_p=other.freqInterpMethod_p;
      spwChanSelFlag_p.resize();
      spwChanSelFlag_p=other.spwChanSelFlag_p;
      freqFrameValid_p=other.freqFrameValid_p;
      //selectedSpw_p.resize();
      //selectedSpw_p=other.selectedSpw_p;
      imageFreq_p.resize();
      imageFreq_p=other.imageFreq_p;
      lsrFreq_p.resize();
      lsrFreq_p=other.lsrFreq_p;
      interpVisFreq_p.resize();
      interpVisFreq_p=other.interpVisFreq_p;
      //multiChanMap_p=other.multiChanMap_p;
      chanMap.resize();
      chanMap=other.chanMap;
      polMap.resize();
      polMap=other.polMap;
      nVisChan_p.resize();
      nVisChan_p=other.nVisChan_p;
      spectralCoord_p=other.spectralCoord_p;
      //doConversion_p.resize();
      //doConversion_p=other.doConversion_p;
      pointingDirCol_p=other.pointingDirCol_p;
      //moving source stuff
      movingDir_p=other.movingDir_p;
      fixMovingSource_p=other.fixMovingSource_p;
      firstMovingDir_p=other.firstMovingDir_p;
      
      //Double precision gridding for those FTMachines that can do
      useDoubleGrid_p=other.useDoubleGrid_p;
      cfStokes_p = other.cfStokes_p;
      polInUse_p = other.polInUse_p;
      cfs_p = other.cfs_p;
      cfwts_p = other.cfwts_p;
      cfs2_p = other.cfs2_p;
      cfwts2_p = other.cfwts2_p;
      canComputeResiduals_p = other.canComputeResiduals_p;

      pop_p = other.pop_p;
      toVis_p = other.toVis_p;
      spwFreqSel_p.resize();
      spwFreqSel_p = other.spwFreqSel_p;
      expandedSpwFreqSel_p = other.expandedSpwFreqSel_p;
      expandedSpwConjFreqSel_p = other.expandedSpwConjFreqSel_p;
      cmplxImage_p=other.cmplxImage_p;
      numthreads_p=other.numthreads_p;
      pbLimit_p=other.pbLimit_p;
      convFuncCtor_p = other.convFuncCtor_p;      
      sj_p.resize();
      sj_p=other.sj_p;
      isDryRun=other.isDryRun;
    };
    return *this;
  };
  
  FTMachine* FTMachine::cloneFTM(){
    Record rec;
    String err;
    if(!(this->toRecord(err, rec)))
      throw(AipsError("Error in cloning FTMachine"));
    return VisModelData::NEW_FT(rec);
  }

  //----------------------------------------------------------------------
  Bool FTMachine::changed(const vi::VisBuffer2&) {
      return false;
    }
  //----------------------------------------------------------------------
  FTMachine::FTMachine(const FTMachine& other)
  {
    operator=(other);
  }
  
  Bool FTMachine::doublePrecGrid(){
    return useDoubleGrid_p;
  }
  
  
  //----------------------------------------------------------------------
   void FTMachine::initPolInfo(const vi::VisBuffer2& vb)
   {
     //
     // Need to figure out where to compute the following arrays/ints
     // in the re-factored code.
     // ----------------------------------------------------------------
     {
       polInUse_p = 0;
       uInt N=0;
       for(uInt i=0;i<polMap.nelements();i++) if (polMap(i) > -1) polInUse_p++;
       cfStokes_p.resize(polInUse_p);
       for(uInt i=0;i<polMap.nelements();i++)
 	if (polMap(i) > -1) {cfStokes_p(N) = vb.correlationTypes()(i);N++;}
     }
   }
  //----------------------------------------------------------------------
    void FTMachine::initMaps(const vi::VisBuffer2& vb) {

      logIO() << LogOrigin("FTMachine", "initMaps") << LogIO::NORMAL;

      AlwaysAssert(image, AipsError);

      // Set the frame for the UVWMachine
      if(vb.isAttached())
	mFrame_p=MeasFrame(MEpoch(Quantity(vb.time()(0), "s"), ROMSColumns(vb.ms()).timeMeas()(0).getRef()), mLocation_p);
      else
	throw(AipsError("Cannot define some frame as no Visiter/MS is attached"));
      // First get the CoordinateSystem for the image and then find
      // the DirectionCoordinate
      casacore::CoordinateSystem coords=image->coordinates();
      Int directionIndex=coords.findCoordinate(Coordinate::DIRECTION);
      AlwaysAssert(directionIndex>=0, AipsError);
      DirectionCoordinate
        directionCoord=coords.directionCoordinate(directionIndex);

      // get the first position of moving source
      if(fixMovingSource_p){

        //First convert to HA-DEC or AZEL for parallax correction
        MDirection::Ref outref1(MDirection::AZEL, mFrame_p);
        MDirection tmphadec=MDirection::Convert(movingDir_p, outref1)();
        MDirection::Ref outref(directionCoord.directionType(), mFrame_p);
        firstMovingDir_p=MDirection::Convert(tmphadec, outref)();

      }


      // Now we need MDirection of the image phase center. This is
      // what we define it to be. So we define it to be the
      // center pixel. So we have to do the conversion here.
      // This is independent of padding since we just want to know
      // what the world coordinates are for the phase center
      // pixel
      {
        Vector<Double> pixelPhaseCenter(2);
        pixelPhaseCenter(0) = Double( image->shape()(0) / 2 );
        pixelPhaseCenter(1) = Double( image->shape()(1) / 2 );
        directionCoord.toWorld(mImage_p, pixelPhaseCenter);
      }

      // Decide if uvwrotation is not necessary, if phasecenter and
      // image center are with in one pixel distance; Save some
      //  computation time especially for spectral cubes.
      {
        Vector<Double> equal= (mImage_p.getAngle()-
  			     vb.phaseCenter().getAngle()).getValue();
        if((abs(equal(0)) < abs(directionCoord.increment()(0)))
  	 && (abs(equal(1)) < abs(directionCoord.increment()(1)))){
  	doUVWRotation_p=false;
        }
        else{
  	doUVWRotation_p=true;
        }
      }
      // Get the object distance in meters
      Record info(image->miscInfo());
      if(info.isDefined("distance")) {
        info.get("distance", distance_p);
        if(abs(distance_p)>0.0) {
  	logIO() << "Distance to object is set to " << distance_p/1000.0
  		<< "km: applying focus correction" << LogIO::POST;
        }
      }

      // Set up the UVWMachine.
      if(uvwMachine_p) delete uvwMachine_p; uvwMachine_p=0;
      String observatory;
      if(vb.isAttached())
	observatory=ROMSColumns(vb.ms()).observation().telescopeName()(0);
      else
	throw(AipsError("Cannot define frame because of no access to OBSERVATION table")); 
      if(observatory.contains("ATCA") || observatory.contains("DRAO")
         || observatory.contains("WSRT")){
        uvwMachine_p=new casacore::UVWMachine(mImage_p, vb.phaseCenter(), mFrame_p,
  				  true, false);
      }
      else{
        uvwMachine_p=new casacore::UVWMachine(mImage_p, vb.phaseCenter(), mFrame_p,
  				  false, tangentSpecified_p);
      }
      AlwaysAssert(uvwMachine_p, AipsError);
      
      lastFieldId_p=-1;
      
      lastMSId_p=vb.msId();
      phaseShifter_p=new UVWMachine(*uvwMachine_p);
      // Set up maps
      Int spectralIndex=coords.findCoordinate(Coordinate::SPECTRAL);
      AlwaysAssert(spectralIndex>-1, AipsError);
      
      spectralCoord_p=coords.spectralCoordinate(spectralIndex);

     
      //Store the image/grid channels freq values
      {
        Int chanNumbre=image->shape()(3);
        Vector<Double> pixindex(chanNumbre);
        imageFreq_p.resize(chanNumbre);
        Vector<Double> tempStorFreq(chanNumbre);
        indgen(pixindex);
        //    pixindex=pixindex+1.0;
        for (Int ll=0; ll< chanNumbre; ++ll){
  	if( !spectralCoord_p.toWorld(tempStorFreq(ll), pixindex(ll))){
  	  logIO() << "Cannot get imageFreq " << LogIO::EXCEPTION;

  	}
        }
        convertArray(imageFreq_p,tempStorFreq);
      }
      //Destroy any conversion layer Freq coord if freqframe is not valid
      if(!freqFrameValid_p){
        MFrequency::Types imageFreqType=spectralCoord_p.frequencySystem();
        spectralCoord_p.setFrequencySystem(imageFreqType);
        spectralCoord_p.setReferenceConversion(imageFreqType,
  					     MEpoch(Quantity(vb.time()(0), "s")),
  					     mLocation_p,
  					     mImage_p);
      }

      // Channel map: do this properly by looking up the frequencies
      // If a visibility channel does not map onto an image
      // pixel then we set the corresponding chanMap to -1.
      // This means that put and get must always check for this
      // value (see e.g. GridFT)

      nvischan  = vb.getFrequencies(0).nelements();
      interpVisFreq_p.resize();
      interpVisFreq_p=vb.getFrequencies(0);
      /*if(selectedSpw_p.nelements() < 1){
        Vector<Int> myspw(1);
        myspw[0]=vb.spectralWindows()(0);
        setSpw(myspw, freqFrameValid_p);
      }
      */

      //matchAllSpwChans(vb);
      
      chanMap.resize();
      matchChannel(vb);
      //chanMap=multiChanMap_p[vb.spectralWindows()(0)];
      if(chanMap.nelements() == 0)
        chanMap=Vector<Int>(vb.getFrequencies(0).nelements(), -1);

      {
        //logIO() << LogIO::DEBUGGING << "Channel Map: " << chanMap << LogIO::POST;
      }
      // Should never get here
      if(max(chanMap)>=nchan||min(chanMap)<-2) {
        logIO() << "Illegal Channel Map: " << chanMap << LogIO::EXCEPTION;
      }

      // Polarization map
      Int stokesIndex=coords.findCoordinate(Coordinate::STOKES);
      AlwaysAssert(stokesIndex>-1, AipsError);
      StokesCoordinate stokesCoord=coords.stokesCoordinate(stokesIndex);
      Vector<Stokes::StokesTypes> visPolMap(vb.getCorrelationTypesSelected());
      nvispol=visPolMap.nelements();
      AlwaysAssert(nvispol>0, AipsError);
      polMap.resize(nvispol);
      polMap=-1;
      Int pol=0;
      Bool found=false;
      // First we try matching Stokes in the visibilities to
      // Stokes in the image that we are gridding into.
      for (pol=0;pol<nvispol;pol++) {
        Int p=0;
        if(stokesCoord.toPixel(p, Stokes::type(visPolMap(pol)))) {
        	AlwaysAssert(p<npol, AipsError);
        	polMap(pol)=p;
        	found=true;
        }
      }
      // If this fails then perhaps we were looking to grid I
      // directly. If so then we need to check that the parallel
      // hands are present in the visibilities.
      if(!found) {
    	  Int p=0;
    	  if(stokesCoord.toPixel(p, Stokes::I)) {
    		  polMap=-1;
    		  if(vb.polarizationFrame()==MSIter::Linear) {
    			  p=0;
    			  for (pol=0;pol<nvispol;pol++) {
    				  if(Stokes::type(visPolMap(pol))==Stokes::XX)
    				  {polMap(pol)=0;p++;found=true;};
    				  if(Stokes::type(visPolMap(pol))==Stokes::YY)
    				  {polMap(pol)=0;p++;found=true;};
    			  }
        	}
        	else {
        		p=0;
        		for (pol=0;pol<nvispol;pol++) {
        			if(Stokes::type(visPolMap(pol))==Stokes::LL)
        			{polMap(pol)=0;p++;found=true;};
        			if(Stokes::type(visPolMap(pol))==Stokes::RR)
        			{polMap(pol)=0;p++;found=true;};
        		}
        	}
    		if(!found) {
    			logIO() <<  "Cannot find polarization map: visibility polarizations = "
    					<< visPolMap << LogIO::EXCEPTION;
    		}
    	else {
    		
    		//logIO() << LogIO::DEBUGGING << "Transforming I only" << LogIO::POST;
    	}
    	  };
      }
      //logIO() << LogIO::DEBUGGING << "Polarization map = "<< polMap
      //	    << LogIO::POST;

      initPolInfo(vb);
      Vector<Int> intpolmap(visPolMap.nelements());
      for (uInt kk=0; kk < intpolmap.nelements(); ++kk){
	intpolmap[kk]=Int(visPolMap[kk]);
      }
      pop_p->initCFMaps(intpolmap, polMap);

      //cerr << "initmaps polmap "<< polMap << endl;

    }

  FTMachine::~FTMachine() 
  {
    if(uvwMachine_p) delete uvwMachine_p; uvwMachine_p=0;
  }
  
  
  Bool FTMachine::interpolateFrequencyTogrid(const vi::VisBuffer2& vb,
  					     const Matrix<Float>& wt,
  					     Cube<Complex>& data,
  					     Cube<Int>& flags,
  					     Matrix<Float>& weight,
  					     FTMachine::Type type){
      Cube<Complex> origdata;
      Cube<Bool> modflagCube;
      // Read flags from the vb.
      setSpectralFlag(vb,modflagCube);
      Vector<Double> visFreq(vb.getFrequencies(0).nelements());
      //if(doConversion_p[vb.spectralWindows()[0]]){
      if(freqFrameValid_p){
        visFreq.resize(lsrFreq_p.shape());
        convertArray(visFreq, lsrFreq_p);
      }
      else{
        convertArray(visFreq, vb.getFrequencies(0));
        lsrFreq_p.resize();
        lsrFreq_p=vb.getFrequencies(0);
      }
      if(type==FTMachine::MODEL){
        origdata.reference(vb.visCubeModel());
      }
      else if(type==FTMachine::CORRECTED){
        origdata.reference(vb.visCubeCorrected());
      }
      else if(type==FTMachine::OBSERVED){
        origdata.reference(vb.visCube());
      }
      else if(type==FTMachine::PSF){
        // make sure its a size 0 data ...psf
        //so avoid reading any data from disk
        origdata.resize();

      }
      else{
        throw(AipsError("Don't know which column is being regridded"));
      }
      if((imageFreq_p.nelements()==1) || (freqInterpMethod_p== InterpolateArray1D<Double, Complex>::nearestNeighbour) || (vb.nChannels()==1)){
        data.reference(origdata);
        // do something here for apply flag based on spw chan sels
        // e.g.
        
        
        flags.resize(modflagCube.shape());
        flags=0;
        //flags(vb.flagCube())=true;
	
	flags(modflagCube)=true;
	
        weight.reference(wt);
        interpVisFreq_p.resize();
        interpVisFreq_p=lsrFreq_p;

        return false;
      }

      Cube<Bool>flag;

      //okay at this stage we have at least 2 channels
      Double width=fabs(imageFreq_p[1]-imageFreq_p[0])/fabs(visFreq[1]-visFreq[0]);
      //if width is smaller than number of points needed for interpolation ...do it directly
      //
      // If image chan width is more than twice the data chan width, make a new list of
      // data frequencies on which to interpolate. This new list is sync'd with the starting image chan
      // and have the same width as the data chans.
      if(((width >2.0) && (freqInterpMethod_p==InterpolateArray1D<Double, Complex>::linear)) ||
         ((width >4.0) && (freqInterpMethod_p !=InterpolateArray1D<Double, Complex>::linear))){
        Double minVF=min(visFreq);
        Double maxVF=max(visFreq);
        Double minIF=min(imageFreq_p);
        Double maxIF=max(imageFreq_p);
        if( ((minIF-fabs(imageFreq_p[1]-imageFreq_p[0])/2.0) > maxVF) ||
  	  ((maxIF+fabs(imageFreq_p[1]-imageFreq_p[0])/2.0) < minVF)){
  	//This function should not have been called with image
  	//being out of bound of data...but still
  	interpVisFreq_p.resize(imageFreq_p.nelements());
  	interpVisFreq_p=imageFreq_p;
  	chanMap.resize(interpVisFreq_p.nelements());
  	chanMap.set(-1);
        }
        else{ // Make a new list of frequencies.
  	Bool found;
  	uInt where=0;
  	Double interpwidth=visFreq[1]-visFreq[0];
  	if(minIF < minVF){ // Need to find the first image-channel with data in it
  	  where=binarySearchBrackets(found, imageFreq_p, minVF, imageFreq_p.nelements());
  	  if(where != imageFreq_p.nelements()){
  	    minIF=imageFreq_p[where];
  	  }
  	}

  	if(maxIF > maxVF){
  	   where=binarySearchBrackets(found, imageFreq_p, maxVF, imageFreq_p.nelements());
  	   if(where!= imageFreq_p.nelements()){
  	    maxIF=imageFreq_p[where];
  	   }

  	}

          // This new list of frequencies starts at the first image channel minus half image channel.
  	// It ends at the last image channel plus half image channel.
  	Int ninterpchan=(Int)ceil((maxIF-minIF+fabs(imageFreq_p[1]-imageFreq_p[0]))/fabs(interpwidth));
  	chanMap.resize(ninterpchan);
  	chanMap.set(-1);
  	interpVisFreq_p.resize(ninterpchan);
  	interpVisFreq_p[0]=(interpwidth > 0) ? minIF : maxIF;
  	interpVisFreq_p[0] -= fabs(imageFreq_p[1]-imageFreq_p[0])/2.0;
  	for (Int k=1; k < ninterpchan; ++k){
  	  interpVisFreq_p[k] = interpVisFreq_p[k-1]+ interpwidth;
  	}

  	for (Int k=0; k < ninterpchan; ++k){
  	  ///chanmap with width
  	  Double nearestchanval = interpVisFreq_p[k]- (imageFreq_p[1]-imageFreq_p[0])/2.0;
  	  where=binarySearchBrackets(found, imageFreq_p, nearestchanval, imageFreq_p.nelements());
  	  if(where != imageFreq_p.nelements())
  	    chanMap[k]=where;
  	}

        }// By now, we have a new list of frequencies, synchronized with image channels, but with data chan widths.
      }// end of ' if (we have to make new frequencies) '
      else{
        // Interpolate directly onto output image frequencies.
        interpVisFreq_p.resize(imageFreq_p.nelements());
        convertArray(interpVisFreq_p, imageFreq_p);
        chanMap.resize(interpVisFreq_p.nelements());
        indgen(chanMap);
      }

      

      if(type != FTMachine::PSF){ // Interpolating the data
   	//Need to get  new interpolate functions that interpolate explicitly on the 2nd axis
  	//2 swap of axes needed
  	Cube<Complex> flipdata;
  	Cube<Bool> flipflag;

          // Interpolate the data.
          //      Input flags are from the previous step ( setSpectralFlag ).
          //      Output flags contain info about channels that could not be interpolated
          //                                   (for example, linear interp with only one data point)
  	swapyz(flipflag,modflagCube);
  	swapyz(flipdata,origdata);
  	InterpolateArray1D<Double,Complex>::
  	  interpolate(data,flag,interpVisFreq_p,visFreq,flipdata,flipflag,freqInterpMethod_p);
  	flipdata.resize();
  	swapyz(flipdata,data);
  	data.resize();
  	data.reference(flipdata);
  	flipflag.resize();
  	swapyz(flipflag,flag);
  	flag.resize();
  	flag.reference(flipflag);
          // Note : 'flag' will get augmented with the flags coming out of weight interpolation
       }
      else
        { // get the flag array to the correct shape.
  	// This will get filled at the end of weight-interpolation.
           flag.resize(vb.nCorrelations(), interpVisFreq_p.nelements(), vb.nRows());
           flag.set(false);
      }
        // Now, interpolate the weights also.
        //   (1) Read in the flags from the vb ( setSpectralFlags -> modflagCube )
        //   (2) Collapse the flags along the polarization dimension to match shape of weight.
         Matrix<Bool> chanflag(wt.shape());
         AlwaysAssert( chanflag.shape()[0]==modflagCube.shape()[1], AipsError);
         AlwaysAssert( chanflag.shape()[1]==modflagCube.shape()[2], AipsError);
         chanflag=false;
         for(uInt pol=0;pol<modflagCube.shape()[0];pol++)
  	 chanflag = chanflag | modflagCube.yzPlane(pol);

         // (3) Interpolate the weights.
         //      Input flags are the collapsed vb flags : 'chanflag'
         //      Output flags are in tempoutputflag
         //            - contains info about channels that couldn't be interpolated.
         Matrix<Float> flipweight;
         flipweight=transpose(wt);
         Matrix<Bool> flipchanflag;
         flipchanflag=transpose(chanflag);
         Matrix<Bool> tempoutputflag;
         InterpolateArray1D<Double,Float>::
  	 interpolate(weight,tempoutputflag, interpVisFreq_p, visFreq,flipweight,flipchanflag,freqInterpMethod_p);
         flipweight.resize();
         flipweight=transpose(weight);
         weight.resize();
         weight.reference(flipweight);
         flipchanflag.resize();
         flipchanflag=transpose(tempoutputflag);
         tempoutputflag.resize();
         tempoutputflag.reference(flipchanflag);

         // (4) Now, fill these flags back into the flag cube
         //                 so that they get USED while gridding the PSF (and data)
         //      Taking the OR of the flags that came out of data-interpolation
         //                         and weight-interpolation, in case they're different.
         //      Expanding flags across polarization.  This will destroy any
         //                          pol-dependent flags for imaging, but msvis::VisImagingWeight
         //                          uses the OR of flags across polarization anyway
         //                          so we don't lose anything.

         AlwaysAssert( tempoutputflag.shape()[0]==flag.shape()[1], AipsError);
         AlwaysAssert( tempoutputflag.shape()[1]==flag.shape()[2], AipsError);
         for(uInt pol=0;pol<flag.shape()[0];pol++)
  	 flag.yzPlane(pol) = tempoutputflag | flag.yzPlane(pol);

         // Fill the output array of image-channel flags.
         flags.resize(flag.shape());
         flags=0;
         flags(flag)=true;

      return true;
    }

  void FTMachine::getInterpolateArrays(const vi::VisBuffer2& vb,
  				       Cube<Complex>& data, Cube<Int>& flags){


      if((imageFreq_p.nelements()==1) || (freqInterpMethod_p== InterpolateArray1D<Double, Complex>::nearestNeighbour)||  (vb.nChannels()==1)){
        Cube<Bool> modflagCube;
        setSpectralFlag(vb,modflagCube);
	
	 data.reference(vb.visCubeModel());
        //flags.resize(vb.flagCube().shape());
        flags.resize(modflagCube.shape());
        flags=0;
        //flags(vb.flagCube())=true;
        flags(modflagCube)=true;
        interpVisFreq_p.resize();
        interpVisFreq_p=vb.getFrequencies(0);
        return;
      }

      data.resize(vb.nCorrelations(), imageFreq_p.nelements(), vb.nRows());
      flags.resize(vb.nCorrelations(), imageFreq_p.nelements(), vb.nRows());
      data.set(Complex(0.0,0.0));
      flags.set(0);
      //no need to degrid channels that does map over this vb
      Int maxchan=max(chanMap);
      for (uInt k =0 ; k < chanMap.nelements() ; ++k){
        if(chanMap(k)==-1)
  	chanMap(k)=maxchan;
      }
      Int minchan=min(chanMap);
      if(minchan==maxchan)
        minchan=-1;


      for(Int k = 0; k < minchan; ++k)
        flags.xzPlane(k).set(1);

      for(uInt k = maxchan + 1; k < imageFreq_p.nelements(); ++k)
        flags.xzPlane(k).set(1);

      interpVisFreq_p.resize(imageFreq_p.nelements());
      convertArray(interpVisFreq_p, imageFreq_p);
      chanMap.resize(imageFreq_p.nelements());
      indgen(chanMap);
    }

  Bool FTMachine::interpolateFrequencyFromgrid(vi::VisBuffer2& vb,
  					       Cube<Complex>& data,
  					       FTMachine::Type type){

      Cube<Complex> *origdata;
      Vector<Double> visFreq(vb.getFrequencies(0).nelements());

      //if(doConversion_p[vb.spectralWindows()[0]]){
      if(freqFrameValid_p){
        convertArray(visFreq, lsrFreq_p);
      }
      else{
        convertArray(visFreq, vb.getFrequencies(0));
      }

      if(type==FTMachine::MODEL){
        origdata=const_cast <Cube<Complex>* > (&(vb.visCubeModel()));
      }
      else if(type==FTMachine::CORRECTED){
        origdata=const_cast<Cube<Complex>* >(&(vb.visCubeCorrected()));
      }
      else{
        origdata=const_cast<Cube<Complex>* >(&(vb.visCube()));
      }
    //
    // If visibility data (vb) has only one channel, or the image cube
    // has only one channel, resort to nearestNeighbour interpolation.
    // Honour user selection of nearestNeighbour.
    //
    if((imageFreq_p.nelements()==1) || 
       (vb.nChannels()==1) || 
       (freqInterpMethod_p== InterpolateArray1D<Double, Complex>::nearestNeighbour)){
        origdata->reference(data);
        return false;
      }

      //Need to get  new interpolate functions that interpolate explicitly on the 2nd axis
      //2 swap of axes needed
      Cube<Complex> flipgrid;
      flipgrid.resize();
      swapyz(flipgrid,data);

      Cube<Complex> flipdata((origdata->shape())(0),(origdata->shape())(2),
  			   (origdata->shape())(1)) ;
      flipdata.set(Complex(0.0));
      InterpolateArray1D<Double,Complex>::
        interpolate(flipdata,visFreq, imageFreq_p, flipgrid,freqInterpMethod_p);
      

      
      Cube<Bool>  copyOfFlag;
      //Vector<Int> mychanmap=multiChanMap_p[vb.spectralWindows()[0]];
      copyOfFlag.assign(vb.flagCube());
      for (uInt k=0; k< chanMap.nelements(); ++ k)
	if(chanMap(k) == -1)
	  copyOfFlag.xzPlane(k).set(true);
      flipgrid.resize();
      swapyz(flipgrid, copyOfFlag, flipdata);
      //swapyz(flipgrid,flipdata);
      vb.setVisCubeModel(flipgrid);

      return true;
    }


  void FTMachine::girarUVW(Matrix<Double>& uvw, Vector<Double>& dphase,
			 const vi::VisBuffer2& vb)
{
    
    
    
    //the uvw rotation is done for common tangent reprojection or if the 
    //image center is different from the phasecenter
    // UVrotation is false only if field never changes
  
   ROMSColumns mscol(vb.ms());
   if((vb.fieldId()(0)!=lastFieldId_p) || (vb.msId()!=lastMSId_p))
      doUVWRotation_p=true;
    if(doUVWRotation_p ||  fixMovingSource_p){
      
      mFrame_p.epoch() != 0 ? 
	mFrame_p.resetEpoch(MEpoch(Quantity(vb.time()(0), "s"))):
	mFrame_p.set(mLocation_p, MEpoch(Quantity(vb.time()(0), "s"), mscol.timeMeas()(0).getRef()));
      MDirection::Types outType;
      MDirection::getType(outType, mImage_p.getRefString());
      MDirection phasecenter=MDirection::Convert(vb.phaseCenter(), MDirection::Ref(outType, mFrame_p))();
      

      if(fixMovingSource_p){
       
      //First convert to HA-DEC or AZEL for parallax correction
	MDirection::Ref outref1(MDirection::AZEL, mFrame_p);
	MDirection tmphadec=MDirection::Convert(movingDir_p, outref1)();
	MDirection::Ref outref(mImage_p.getRef().getType(), mFrame_p);
	MDirection sourcenow=MDirection::Convert(tmphadec, outref)();
	//cerr << "Rotating to fixed moving source " << MVDirection(phasecenter.getAngle()-firstMovingDir_p.getAngle()+sourcenow.getAngle()) << endl;
	phasecenter.set(MVDirection(phasecenter.getAngle()+firstMovingDir_p.getAngle()-sourcenow.getAngle()));
	
    }


      // Set up the UVWMachine only if the field id has changed. If
      // the tangent plane is specified then we need a UVWMachine that
      // will reproject to that plane iso the image plane
      if((vb.fieldId()(0)!=lastFieldId_p) || (vb.msId()!=lastMSId_p) || fixMovingSource_p) {
	
	String observatory=mscol.observation().telescopeName()(0);
	if(uvwMachine_p) delete uvwMachine_p; uvwMachine_p=0;
	if(observatory.contains("ATCA") || observatory.contains("WSRT")){
		//Tangent specified is being wrongly used...it should be for a
	    	//Use the safest way  for now.
	    uvwMachine_p=new UVWMachine(phasecenter, vb.phaseCenter(), mFrame_p,
					true, false);
	    phaseShifter_p=new UVWMachine(mImage_p, phasecenter, mFrame_p,
					true, false);
	}
	else{
	  uvwMachine_p=new UVWMachine(phasecenter, vb.phaseCenter(),  mFrame_p,
				      false, false);
	  phaseShifter_p=new UVWMachine(mImage_p, phasecenter,  mFrame_p,
				      false, false);
	}
      }

      lastFieldId_p=vb.fieldId()(0);
	lastMSId_p=vb.msId();

      
      AlwaysAssert(uvwMachine_p, AipsError);
      
      // Always force a recalculation 
      uvwMachine_p->reCalculate();
      phaseShifter_p->reCalculate();
      
      // Now do the conversions
      uInt nrows=dphase.nelements();
      Vector<Double> thisRow(3);
      thisRow=0.0;
      //CoordinateSystem csys=image->coordinates();
      //DirectionCoordinate dc=csys.directionCoordinate(0);
      //Vector<Double> thePix(2);
      //dc.toPixel(thePix, phasecenter);
      //cerr << "field id " << vb.fieldId() << "  the Pix " << thePix << endl;
      //Vector<Float> scale(2);
      //scale(0)=dc.increment()(0);
      //scale(1)=dc.increment()(1);
      for (uInt irow=0; irow<nrows;++irow) {
	thisRow.assign(uvw.column(irow));
	//cerr << " uvw " << thisRow ;
	// This is for frame change
	uvwMachine_p->convertUVW(dphase(irow), thisRow);
	// This is for correlator phase center change
	MVPosition rotphase=phaseShifter_p->rotationPhase() ;
	//cerr << " rotPhase " <<  rotphase << " oldphase "<<  rotphase*(uvw.column(irow))  << " newphase " << (rotphase)*thisRow ;
	//	cerr << " phase " << dphase(irow) << " new uvw " << uvw.column(irow);
	//dphase(irow)+= (thePix(0)-nx/2.0)*thisRow(0)*scale(0)+(thePix(1)-ny/2.0)*thisRow(1)*scale(1);
	//Double pixphase=(thePix(0)-nx/2.0)*uvw.column(irow)(0)*scale(0)+(thePix(1)-ny/2.0)*uvw.column(irow)(1)*scale(1);
	//Double pixphase2=(thePix(0)-nx/2.0)*thisRow(0)*scale(0)+(thePix(1)-ny/2.0)*thisRow(1)*scale(1);
	//cerr << " pixphase " <<  pixphase <<  " pixphase2 " << pixphase2<< endl;
	//dphase(irow)=pixphase;
	RotMatrix rotMat=phaseShifter_p->rotationUVW();
	uvw.column(irow)(0)=thisRow(0)*rotMat(0,0)+thisRow(1)*rotMat(1,0);
	uvw.column(irow)(1)=thisRow(1)*rotMat(1,1)+thisRow(0)*rotMat(0,1);
	uvw.column(irow)(2)=thisRow(0)*rotMat(0,2)+thisRow(1)*rotMat(1,2)+thisRow(2)*rotMat(2,2);
	dphase(irow)+= rotphase(0)*uvw.column(irow)(0)+rotphase(1)*uvw.column(irow)(1);
      }
	
      
    }
}


  void FTMachine::rotateUVW(Matrix<Double>& uvw, Vector<Double>& dphase,
  			    const vi::VisBuffer2& vb)
    {



      //the uvw rotation is done for common tangent reprojection or if the
      //image center is different from the phasecenter
      // UVrotation is false only if field never changes
      if((vb.fieldId()(0)!=lastFieldId_p) || (vb.msId()!=lastMSId_p))
        doUVWRotation_p=true;
      if(doUVWRotation_p || tangentSpecified_p || fixMovingSource_p){
        ok();
	
        mFrame_p.epoch() != 0 ?
	  mFrame_p.resetEpoch(MEpoch(Quantity(vb.time()(0), "s"))):
	 
	  mFrame_p.set(mLocation_p, MEpoch(Quantity(vb.time()(0), "s"), ROMSColumns(vb.ms()).timeMeas()(0).getRef()));

        MDirection phasecenter=mImage_p;
        if(fixMovingSource_p){

        //First convert to HA-DEC or AZEL for parallax correction
  	MDirection::Ref outref1(MDirection::AZEL, mFrame_p);
  	MDirection tmphadec=MDirection::Convert(movingDir_p, outref1)();
  	MDirection::Ref outref(mImage_p.getRef().getType(), mFrame_p);
  	MDirection sourcenow=MDirection::Convert(tmphadec, outref)();
  	//cerr << "Rotating to fixed moving source " << MVDirection(phasecenter.getAngle()-firstMovingDir_p.getAngle()+sourcenow.getAngle()) << endl;
  	phasecenter.set(MVDirection(phasecenter.getAngle()+firstMovingDir_p.getAngle()-sourcenow.getAngle()));

      }


        // Set up the UVWMachine only if the field id has changed. If
        // the tangent plane is specified then we need a UVWMachine that
        // will reproject to that plane iso the image plane
        if((vb.fieldId()(0)!=lastFieldId_p) || (vb.msId()!=lastMSId_p) || fixMovingSource_p) {

  	String observatory=ROMSColumns(vb.ms()).observation().telescopeName()(0);
  	if(uvwMachine_p) delete uvwMachine_p; uvwMachine_p=0;
  	if(observatory.contains("ATCA") || observatory.contains("WSRT")){
  		//Tangent specified is being wrongly used...it should be for a
  	    	//Use the safest way  for now.
  	    uvwMachine_p=new UVWMachine(phasecenter, vb.phaseCenter(), mFrame_p,
  					true, false);
  	}
  	else{
  		uvwMachine_p=new UVWMachine(phasecenter, vb.phaseCenter(), mFrame_p,
  					false,tangentSpecified_p);
  	    }
       }

  	lastFieldId_p=vb.fieldId()(0);
  	lastMSId_p=vb.msId();


        AlwaysAssert(uvwMachine_p, AipsError);

        // Always force a recalculation
        uvwMachine_p->reCalculate();

        // Now do the conversions
        uInt nrows=dphase.nelements();
        Vector<Double> thisRow(3);
        thisRow=0.0;
        uInt irow;
        //#pragma omp parallel default(shared) private(irow,thisRow)
        {
  	//#pragma omp for
  	  for (irow=0; irow<nrows;++irow) {
  	    thisRow.reference(uvw.column(irow));
  	    convUVW(dphase(irow), thisRow);
  	  }

        }//end pragma
      }
    }
  void FTMachine::convUVW(Double& dphase, Vector<Double>& thisrow){
    //for (uInt i=0;i<3;i++) thisRow(i)=uvw(i,row);
    uvwMachine_p->convertUVW(dphase, thisrow);
    //for (uint i=0;i<3;i++) uvw(i,row)=thisRow(i);

  }


  void FTMachine::locateuvw(const Double*& uvw, const Double*& dphase,
			    const Double*& freq, const Int& nvchan,
			    const Double*& scale, const Double*& offset,  const Int& sampling, Int*& loc, Int*& off, Complex*& phasor, const Int& row, const bool& doW){
    
    Int rowoff=row*nvchan;
    Double phase;
    Double pos;
    Int nel= doW ? 3 : 2;
    for (Int f=0; f<nvchan; ++f){
      for (Int k=0; k <2; ++k){
	pos=(scale[k])*uvw[3*row+k]*(freq[f])/C::c+((offset[k])+1.0);
	loc[(rowoff+f)*nel+k]=std::lround(pos);
	off[(rowoff+f)*nel+k]=std::lround((Double(loc[(rowoff+f)*nel+k])-pos)*Double(sampling));
	//off[(rowoff+f)*2+k]=(loc[(rowoff+f)*2+k]-pos(k))*sampling;	
      }
      phase=-Double(2.0)*C::pi*dphase[row]*(freq[f])/C::c;
      phasor[rowoff+f]=Complex(cos(phase), sin(phase));
     
      ///This is for W-Projection
      if(doW){
	pos=sqrt(fabs(scale[2]*uvw[3*row+2]*(freq[f])/C::c))+offset[2]+1.0;
	loc[(rowoff+f)*nel+2]=std::lround(pos);
	off[(rowoff+f)*nel+2]=0;
      }
    }

    


  }

  void FTMachine::setnumthreads(Int num){
    numthreads_p=num;
  }
  Int FTMachine::getnumthreads(){
    return numthreads_p;
  }

  //
  // Refocus the array on a point at finite distance
  //
    // Refocus the array on a point at finite distance
    //
    void FTMachine::refocus(Matrix<Double>& uvw, const Vector<Int>& ant1,
  			  const Vector<Int>& ant2,
  			  Vector<Double>& dphase, const vi::VisBuffer2& vb)
    {

      ok();

      if(abs(distance_p)>0.0) {

        nAntenna_p=max(vb.antenna2())+1;

        // Positions of antennas
        Matrix<Double> antPos(3,nAntenna_p);
        antPos=0.0;
        Vector<Int> nAntPos(nAntenna_p);
        nAntPos=0;

        uInt aref = min(ant1);

        // Now find the antenna locations: for this we just reference to a common
        // point. We ignore the time variation within this buffer.
        uInt nrows=dphase.nelements();
        for (uInt row=0;row<nrows;row++) {
  	uInt a1=ant1(row);
  	uInt a2=ant2(row);
  	for(uInt dim=0;dim<3;dim++) {
  	  antPos(dim, a1)+=uvw(dim, row);
  	  antPos(dim, a2)-=uvw(dim, row);
  	}
  	nAntPos(a1)+=1;
  	nAntPos(a2)+=1;
        }

        // Now remove the reference location
        Vector<Double> center(3);
        for(uInt dim=0;dim<3;dim++) {
  	center(dim) = antPos(dim,aref)/nAntPos(aref);
        }

        // Now normalize
        for (uInt ant=0; ant<nAntenna_p; ant++) {
  	if(nAntPos(ant)>0) {
  	  for(uInt dim=0;dim<3;dim++) {
  	    antPos(dim,ant)/=nAntPos(ant);
  	    antPos(dim,ant)-=center(dim);
  	  }
  	}
        }

        // Now calculate the offset needed to focus the array,
        // including the w term. This must have the correct asymptotic
        // form so that at infinity no net change occurs
        for (uInt row=0;row<nrows;row++) {
  	uInt a1=ant1(row);
  	uInt a2=ant2(row);

  	Double d1=distance_p*distance_p-2*distance_p*antPos(2,a1);
  	Double d2=distance_p*distance_p-2*distance_p*antPos(2,a2);
  	for(uInt dim=0;dim<3;dim++) {
  	  d1+=antPos(dim,a1)*antPos(dim,a1);
  	  d2+=antPos(dim,a2)*antPos(dim,a2);
  	}
  	d1=sqrt(d1);
  	d2=sqrt(d2);
  	for(uInt dim=0;dim<2;dim++) {
  	  dphase(row)-=(antPos(dim,a1)*antPos(dim,a1)-antPos(dim,a2)*antPos(dim,a2))/(2*distance_p);
  	}
  	uvw(0,row)=distance_p*(antPos(0,a1)/d1-antPos(0,a2)/d2);
  	uvw(1,row)=distance_p*(antPos(1,a1)/d1-antPos(1,a2)/d2);
  	uvw(2,row)=distance_p*(antPos(2,a1)/d1-antPos(2,a2)/d2)+dphase(row);
        }
      }
    }

  void FTMachine::ok() {
    AlwaysAssert(image, AipsError);
    AlwaysAssert(uvwMachine_p, AipsError);
  }
  
  Bool FTMachine::toRecord(String& error, RecordInterface& outRecord, 
			   Bool withImage, const String diskimage) {
    // Save the FTMachine to a Record
    //
    outRecord.define("name", this->name());
    if(withImage){
      CoordinateSystem cs=image->coordinates();
      DirectionCoordinate dircoord=cs.directionCoordinate(cs.findCoordinate(Coordinate::DIRECTION));
      dircoord.setReferenceValue(mImage_p.getAngle().getValue());
      if(diskimage != ""){
	try{
	  PagedImage<Complex> imCopy(TiledShape(toVis_p ? griddedData.shape(): image->shape()), image->coordinates(), diskimage);
	  toVis_p ? imCopy.put(griddedData) : imCopy.copyData(*image);
	  ImageUtilities::copyMiscellaneous(imCopy, *image);
	  Vector<Double> pixcen(2);
	  pixcen(0)=Double(imCopy.shape()(0)/2); pixcen(1)=Double(imCopy.shape()(1)/2);
	  dircoord.setReferencePixel(pixcen);
	  cs.replaceCoordinate(dircoord, cs.findCoordinate(Coordinate::DIRECTION));
	  imCopy.setCoordinateInfo(cs);
	}
	catch(...){
	  throw(AipsError(String("Failed to save model image "+diskimage+String(" to disk")))); 
	}
	outRecord.define("diskimage", diskimage);
	
      }
      else{
	Record imrec;
	Vector<Double> pixcen(2);
	pixcen(0)=Double(griddedData.shape()(0)/2); pixcen(1)=Double(griddedData.shape()(1)/2);
	dircoord.setReferencePixel(pixcen);
	cs.replaceCoordinate(dircoord, cs.findCoordinate(Coordinate::DIRECTION));
	TempImage<Complex> imCopy(griddedData.shape(), cs);
	imCopy.put(griddedData) ;
	ImageUtilities::copyMiscellaneous(imCopy, *image);
	if(imCopy.toRecord(error, imrec))
	  outRecord.defineRecord("image", imrec);
      }
    }
    outRecord.define("nantenna", nAntenna_p);
    outRecord.define("distance", distance_p);
    outRecord.define("lastfieldid", lastFieldId_p);
    outRecord.define("lastmsid", lastMSId_p);
    outRecord.define("tangentspecified", tangentSpecified_p);
    saveMeasure(outRecord, String("mtangent_rec"), error, mTangent_p);
    saveMeasure(outRecord, "mimage_rec", error, mImage_p);
    //mFrame_p not necessary to save as it is generated from mLocation_p
    outRecord.define("nx", nx);
    outRecord.define("ny", ny);
    outRecord.define("npol", npol);
    outRecord.define("nchan", nchan);
    outRecord.define("nvischan", nvischan);
    outRecord.define("nvispol", nvispol);
    saveMeasure(outRecord, "mlocation_rec", error, mLocation_p);
    //no need to save uvwMachine_p
    outRecord.define("douvwrotation", doUVWRotation_p);
    outRecord.define("freqinterpmethod", static_cast<Int>(freqInterpMethod_p));
    outRecord.define("spwchanselflag", spwChanSelFlag_p);
    outRecord.define("freqframevalid", freqFrameValid_p);
    //outRecord.define("selectedspw", selectedSpw_p);
    outRecord.define("imagefreq", imageFreq_p);
    outRecord.define("lsrfreq", lsrFreq_p);
    outRecord.define("interpvisfreq", interpVisFreq_p);
    Record multichmaprec;
    //for (uInt k=0; k < multiChanMap_p.nelements(); ++k)
    //  multichmaprec.define(k, multiChanMap_p[k]);
    //outRecord.defineRecord("multichanmaprec", multichmaprec);
    outRecord.define("chanmap", chanMap);
    outRecord.define("polmap", polMap);
    outRecord.define("nvischanmulti", nVisChan_p);
    spectralCoord_p.save(outRecord, "spectralcoord");
    //outRecord.define("doconversion", doConversion_p);
    outRecord.define("pointingdircol", pointingDirCol_p);
    saveMeasure(outRecord, "movingdir_rec", error, movingDir_p);
    outRecord.define("fixmovingsource", fixMovingSource_p);
    saveMeasure(outRecord, "firstmovingdir_rec", error, firstMovingDir_p);
    outRecord.define("usedoublegrid", useDoubleGrid_p);
    outRecord.define("cfstokes", cfStokes_p);
    outRecord.define("polinuse", polInUse_p);
    outRecord.define("tovis", toVis_p);
    outRecord.define("sumweight", sumWeight);
    outRecord.define("numthreads", numthreads_p);
    //Need to serialized sj_p...the user has to set the sj_p after recovering from record
    return true;
  };
  
  Bool FTMachine::saveMeasure(RecordInterface& rec, const String& name, String& err, const Measure& meas){
    Record tmprec;
    MeasureHolder mh(meas);
    if(mh.toRecord(err, tmprec)){
      rec.defineRecord(name, tmprec);
      return true;
    }
    return false;
  }

  Bool FTMachine::fromRecord(String& error,
			     const RecordInterface& inRecord
			     ) {
    // Restore an FTMachine from a Record
    //
    uvwMachine_p=0; //when null it is reconstructed from mImage_p and mFrame_p
    //mFrame_p is not necessary as it is generated in initMaps from mLocation_p
    inRecord.get("nx", nx);
    inRecord.get("ny", ny);
    inRecord.get("npol", npol);
    inRecord.get("nchan", nchan);
    inRecord.get("nvischan", nvischan);
    inRecord.get("nvispol", nvispol);
    cmplxImage_p=NULL;
    inRecord.get("tovis", toVis_p);
    if(inRecord.isDefined("image")){
      cmplxImage_p=new TempImage<Complex>();
      image=&(*cmplxImage_p);
      
      const Record rec=inRecord.asRecord("image");
      if(!cmplxImage_p->fromRecord(error, rec))
	return false;   
      
    }
    else if(inRecord.isDefined("diskimage")){
      String theDiskImage;
      inRecord.get("diskimage", theDiskImage);
      try{
        File eldir(theDiskImage);
	if(! eldir.exists()){
	  String subPathname[30];
	  String sep = "/";
	  uInt nw = split(theDiskImage, subPathname, 20, sep);
	  String theposs=(subPathname[nw-1]);
	  Bool isExistant=File(theposs).exists();
	  if(isExistant) 
	    theDiskImage=theposs;
	  for (uInt i=nw-2 ; i>0; --i){
	    theposs=subPathname[i]+"/"+theposs;
	    File newEldir(theposs);
	    if(newEldir.exists()){
	      isExistant=true;
	      theDiskImage=theposs;
	    }
	  }
	  if(!isExistant)
	    throw(AipsError("Could not locate mage"));
	}
	cmplxImage_p=new PagedImage<Complex> (theDiskImage);
	image=&(*cmplxImage_p);
      }
      catch(...){
	throw(AipsError(String("Failure to load ")+theDiskImage+String(" image from disk")));
      }
    }
    if(toVis_p && !cmplxImage_p.null()) {
	griddedData.resize(image->shape());
	griddedData=image->get();
    }
    else if(!toVis_p){
      IPosition gridShape(4, nx, ny, npol, nchan);
      griddedData.resize(gridShape);
      griddedData=Complex(0.0);
    }

    nAntenna_p=inRecord.asuInt("nantenna");
    distance_p=inRecord.asDouble("distance");
    lastFieldId_p=inRecord.asInt("lastfieldid");
    lastMSId_p=inRecord.asInt("lastmsid");
    inRecord.get("tangentspecified", tangentSpecified_p);
    { const Record rec=inRecord.asRecord("mtangent_rec");
      MeasureHolder mh;
      if(!mh.fromRecord(error, rec))
	return false;
      mTangent_p=mh.asMDirection();
    }
    { const Record rec=inRecord.asRecord("mimage_rec");
      MeasureHolder mh;
      if(!mh.fromRecord(error, rec))
	return false;
      mImage_p=mh.asMDirection();
    }
    
   
    { const Record rec=inRecord.asRecord("mlocation_rec");
      MeasureHolder mh;
      if(!mh.fromRecord(error, rec))
	return false;
      mLocation_p=mh.asMPosition();
    }
    inRecord.get("douvwrotation", doUVWRotation_p);
   
    //inRecord.get("spwchanselflag", spwChanSelFlag_p);
    //We won't respect the chanselflag as the vister may have different selections
    spwChanSelFlag_p.resize();
    inRecord.get("freqframevalid", freqFrameValid_p);
    //inRecord.get("selectedspw", selectedSpw_p);
    inRecord.get("imagefreq", imageFreq_p);
    inRecord.get("lsrfreq", lsrFreq_p);
    inRecord.get("interpvisfreq", interpVisFreq_p);
    //const Record multichmaprec=inRecord.asRecord("multichanmaprec");
    //multiChanMap_p.resize(multichmaprec.nfields(), true, false);
    //for (uInt k=0; k < multichmaprec.nfields(); ++k)
    //  multichmaprec.get(k, multiChanMap_p[k]);
    inRecord.get("chanmap", chanMap);
    inRecord.get("polmap", polMap);
    inRecord.get("nvischanmulti", nVisChan_p);
    SpectralCoordinate *tmpSpec=SpectralCoordinate::restore(inRecord, "spectralcoord");
    if(tmpSpec){
      spectralCoord_p=*tmpSpec;
      delete tmpSpec;
    }
    //inRecord.get("doconversion", doConversion_p);
    inRecord.get("pointingdircol", pointingDirCol_p);
    { const Record rec=inRecord.asRecord("movingdir_rec");
      MeasureHolder mh;
      if(!mh.fromRecord(error, rec))
	return false;
      movingDir_p=mh.asMDirection();
    }
    inRecord.get("fixmovingsource", fixMovingSource_p);
    { const Record rec=inRecord.asRecord("firstmovingdir_rec");
      MeasureHolder mh;
      if(!mh.fromRecord(error, rec))
	return false;
      firstMovingDir_p=mh.asMDirection();
    }
    inRecord.get("usedoublegrid", useDoubleGrid_p);
    inRecord.get("cfstokes", cfStokes_p);
    inRecord.get("polinuse", polInUse_p);
    
    
    inRecord.get("sumweight", sumWeight);
    if(toVis_p){
      freqInterpMethod_p=InterpolateArray1D<Double, Complex>::nearestNeighbour;
    }
    else{
     Int tmpInt;
      inRecord.get("freqinterpmethod", tmpInt);
      freqInterpMethod_p=static_cast<InterpolateArray1D<Double, Complex >::InterpolationMethod>(tmpInt);
    }
    inRecord.get("numthreads", numthreads_p);
    return true;
  };
  
  // Make a plain straightforward honest-to-FSM image. This returns
  // a complex image, without conversion to Stokes. The representation
  // is that required for the visibilities.
  //----------------------------------------------------------------------
  void FTMachine::makeImage(FTMachine::Type type, 
			    vi::VisibilityIterator2& vi,
			    ImageInterface<Complex>& theImage,
			    Matrix<Float>& weight) {
    
    
    logIO() << LogOrigin("FTMachine", "makeImage0") << LogIO::NORMAL;
    
    // Loop over all visibilities and pixels
    vi::VisBuffer2* vb=vi.getVisBuffer();
    
    // Initialize put (i.e. transform to Sky) for this model
    vi.origin();
    
    if(vb->polarizationFrame()==MSIter::Linear) {
      StokesImageUtil::changeCStokesRep(theImage, StokesImageUtil::LINEAR);
    }
    else {
      StokesImageUtil::changeCStokesRep(theImage, StokesImageUtil::CIRCULAR);
    }
    
    initializeToSky(theImage,weight,*vb);
    Bool useCorrected= !(ROMSColumns(vi.ms()).correctedData().isNull());
    if((type==FTMachine::CORRECTED) && (!useCorrected))
      type=FTMachine::OBSERVED;
    Bool normalize=true;
    if(type==FTMachine::COVERAGE)
      normalize=false;

    // Loop over the visibilities, putting VisBuffers
    for (vi.originChunks();vi.moreChunks();vi.nextChunk()) {
      for (vi.origin(); vi.more(); vi.next()) {
	
	switch(type) {
	case FTMachine::RESIDUAL:
	  vb->setVisCube(vb->visCubeCorrected());
	  vb->setVisCube(vb->visCube()-vb->visCubeModel());
	  put(*vb, -1, false);
	  break;
	case FTMachine::MODEL:
	  put(*vb, -1, false, FTMachine::MODEL);
	  break;
	case FTMachine::CORRECTED:
	  put(*vb, -1, false, FTMachine::CORRECTED);
	  break;
	case FTMachine::PSF:
	  vb->setVisCube(Complex(1.0,0.0));
	  put(*vb, -1, true, FTMachine::PSF);
	  break;
	case FTMachine::COVERAGE:
	  vb->setVisCube(Complex(1.0));
	  put(*vb, -1, true, FTMachine::COVERAGE);
	  break;
	case FTMachine::OBSERVED:
	default:
	  put(*vb, -1, false, FTMachine::OBSERVED);
	  break;
	}
      }
    }
    finalizeToSky();
    // Normalize by dividing out weights, etc.
    getImage(weight, normalize);
  }
  
  
  
  
  Bool FTMachine::setFrameValidity(Bool validFrame){
    
    freqFrameValid_p=validFrame;
    return true;
  }
  
  /*
  Bool FTMachine::matchAllSpwChans(const vi::VisBuffer2& vb){

	  //////I have no clue how to get all the channel and data selection from all
	  ///spectral windows from Visbuffer2...
	  /// so this function is quite useless

	  Vector<Int>  elspw;
	  Vector<Int>  elstart;
	  Vector<Int>  elnchan;
	  Double elfstart, elfend, elfinc;
	  spectralCoord_p.toWorld(elfstart, 0.0);
	  spectralCoord_p.toWorld(elfend, Double(nchan));
	  if(elfend < elfstart){
		  Double tmpfreq=elfstart;
		  elfstart=elfend;
		  elfend=tmpfreq;
	  }
	  elfinc=(spectralCoord_p.increment()(0));

	  cerr << "elfstart " << elfstart << " elfend " << elfend << " elfinc "<< elfinc << endl;

	  MSUtil::getSpwInFreqRangeAllFields(elspw, elstart,
			  elnchan,vb.getVi()->ms(), elfstart,elfend,elfinc, MFrequency::LSRK);
	  selectedSpw_p.resize();
	  selectedSpw_p=elspw;
	  nVisChan_p.resize();
	  nVisChan_p=elnchan;
	  cerr << "elspw " << elspw << " elstart " << elstart  << " elnchan " << endl;

	  //doConversion_p.resize(max(selectedSpw_p)+1);
	  //doConversion_p.set(true);

      multiChanMap_p.resize(max(selectedSpw_p)+1, true);
      matchChannel(vb);
      /*Bool anymatchChan=false;
      Bool anyTopo=false;
      for (uInt k=0; k < selectedSpw_p.nelements(); ++k){
        Bool matchthis=matchChannel(selectedSpw_p[k], vb);
        anymatchChan= (anymatchChan || matchthis);
        anyTopo=anyTopo || ((MFrequency::castType(ROMSColumns(vb.getVi()->ms()).spectralWindow().measFreqRef()(selectedSpw_p[k]))==MFrequency::TOPO) && freqFrameValid_p);
      }

      // if TOPO and valid frame things may match later but not now  thus we'll go
      // through the data
      // hoping the user made the right choice
      if (!anymatchChan && !anyTopo){
        logIO() << "No overlap in frequency between image channels and selected data found for this FTMachine \n"
  	      << " Check your data selection and image parameters if you end up with a blank image"
  	      << LogIO::WARN << LogIO::POST;

      }
     //////////////////////
      return true;

    }
  */

  Bool FTMachine::matchChannel(const vi::VisBuffer2& vb){

    //Int spw=vb.spectralWindows()[0];
    nvischan  = vb.nChannels();
    chanMap.resize(nvischan);
    chanMap.set(-1);
    Vector<Double> lsrFreq(0);

      //cerr << "doConve " << spw << "   " << doConversion_p[spw] << " freqframeval " << freqFrameValid_p << endl;
//cerr <<"valid frame " << freqFrameValid_p << " polmap "<< polMap << endl;
     if(freqFrameValid_p)
    	 lsrFreq=vb.getFrequencies(0,MFrequency::LSRK);
     else
    	 lsrFreq=vb.getFrequencies(0);

     //cerr << "lsrFreq " << lsrFreq.shape() << " nvischan " << nvischan << endl;
     //     if(doConversion_p.nelements() < uInt(spw+1))
     //	 doConversion_p.resize(spw+1, true);
     // doConversion_p[spw]=freqFrameValid_p;

      if(lsrFreq.nelements() ==0){
        return false;
      }
      lsrFreq_p.resize(lsrFreq.nelements());
      lsrFreq_p=lsrFreq;
      Vector<Double> c(1);
      c=0.0;
      Vector<Double> f(1);
      Int nFound=0;


      //cout.precision(10);
      for (Int chan=0;chan<nvischan;chan++) {
        f(0)=lsrFreq[chan];
        if(spectralCoord_p.toPixel(c, f)) {
  	Int pixel=Int(floor(c(0)+0.5));  // round to chan freq at chan center
  	//cerr << " chan " << chan << " f " << f(0) << " pixel "<< c(0) << "  " << pixel << endl;
  	/////////////
  	//c(0)=pixel;
  	//spectralCoord_p.toWorld(f, c);
  	//cerr << "f1 " << f(0) << " pixel "<< c(0) << "  " << pixel << endl;
  	////////////////
  	if(pixel>-1&&pixel<nchan) {
  	  chanMap(chan)=pixel;
  	  nFound++;
  	  if(nvischan>1&&(chan==0||chan==nvischan-1)) {
  	    /*logIO() << LogIO::DEBUGGING
  		    << "Selected visibility channel : " << chan+1
  		    << " has frequency "
  		    <<  MFrequency(Quantity(f(0), "Hz")).get("GHz").getValue()
  		    << " GHz and maps to image pixel " << pixel+1 << LogIO::POST;
  	    */
  	  }
  	}
	else{
	  
	  if(nvischan > 1){
	    Double fwidth=lsrFreq[1]-lsrFreq[0];
	    Double limit=0;
	    Double where=c(0)*fabs(spectralCoord_p.increment()(0));
	    if( freqInterpMethod_p==InterpolateArray1D<Double,Complex>::linear)
	      limit=1;
	    else if( freqInterpMethod_p==InterpolateArray1D<Double,Complex>::cubic ||  freqInterpMethod_p==InterpolateArray1D<Double,Complex>::spline)
	      limit=2;
	    if(((pixel<0) && (where >= (0-limit*fabs(fwidth)))) )
	      chanMap(chan)=-2;
	    if((pixel>=nchan) ) {
	      where=f(0);
	      Double fend;
	      spectralCoord_p.toWorld(fend, Double(nchan));
	      if( ( (fwidth >0) &&where < (fend+limit*fwidth))  || ( (fwidth <0) &&where > (fend+limit*fwidth)) )
		chanMap(chan)=-2;
	    }
	  }


	}
        }
      }
      //cerr << "chanmap " << chanMap << endl;
      /* if(multiChanMap_p.nelements() < uInt(spw+1))
    	  multiChanMap_p.resize(spw+1);
      multiChanMap_p[spw].resize();
      multiChanMap_p[spw]=chanMap;
      */

      if(nFound==0) {
        /*
  	logIO()  << "Visibility channels in spw " << spw+1
  	<<      " of ms " << vb.msId() << " is not being used "
  	<< LogIO::WARN << LogIO::POST;
        */
        return false;
      }




      return true;

    }




  void FTMachine::gridOk(Int convSupport){
    
    if (nx <= 2*convSupport) {
      logIO_p 
	<< "number of pixels "
	<< nx << " on x axis is smaller that the gridding support "
	<< 2*convSupport   << " Please use a larger value " 
	<< LogIO::EXCEPTION;
    }
    
    if (ny <= 2*convSupport) {
      logIO_p 
	<< "number of pixels "
	<< ny << " on y axis is smaller that the gridding support "
	<< 2*convSupport   << " Please use a larger value " 
	<< LogIO::EXCEPTION;
    }
    
  }
  
  void FTMachine::setLocation(const MPosition& loc){
    
    mLocation_p=loc;
    
  }
  
  MPosition& FTMachine::getLocation(){
    
    return mLocation_p;
  }
  
  
  void FTMachine::setMovingSource(const String& sourcename){
    
    fixMovingSource_p=true;
    movingDir_p=MDirection(Quantity(0.0,"deg"), Quantity(90.0, "deg"));
    movingDir_p.setRefString(sourcename);
    
  }

  void FTMachine::setMovingSource(const MDirection& mdir){
    
    fixMovingSource_p=true;
    movingDir_p=mdir;
    
  }
  
  void FTMachine::setFreqInterpolation(const String& method){
    
    String meth=method;
    meth.downcase();
    if(meth.contains("linear")){
      freqInterpMethod_p=InterpolateArray1D<Double,Complex>::linear;
    }
    else if(meth.contains("splin")){
      freqInterpMethod_p=InterpolateArray1D<Double,Complex>::spline;  
    }	    
    else if(meth.contains("cub")){
      freqInterpMethod_p=InterpolateArray1D<Double,Complex>::cubic;
    }
    else{
      freqInterpMethod_p=InterpolateArray1D<Double,Complex>::nearestNeighbour;
    }
  
  }
  
  
  // helper function to swap the y and z axes of a Cube
  void FTMachine::swapyz(Cube<Complex>& out, const Cube<Complex>& in)
  {
    IPosition inShape=in.shape();
    uInt nxx=inShape(0),nyy=inShape(2),nzz=inShape(1);
    //resize breaks  references...so out better have the right shape 
    //if references is not to be broken
    if(out.nelements()==0)
      out.resize(nxx,nyy,nzz);
    Bool deleteIn,deleteOut;
    const Complex* pin = in.getStorage(deleteIn);
    Complex* pout = out.getStorage(deleteOut);
    uInt i=0, zOffset=0;
    for (uInt iz=0; iz<nzz; ++iz, zOffset+=nxx) {
      Int yOffset=zOffset;
      for (uInt iy=0; iy<nyy; ++iy, yOffset+=nxx*nzz) {
	for (uInt ix=0; ix<nxx; ++ix){ 
	  pout[i++] = pin[ix+yOffset];
	}
      }
    }
    out.putStorage(pout,deleteOut);
    in.freeStorage(pin,deleteIn);
  }

  void FTMachine::swapyz(Cube<Complex>& out, const Cube<Bool>& outFlag, const Cube<Complex>& in)
  {
    IPosition inShape=in.shape();
    uInt nxx=inShape(0),nyy=inShape(2),nzz=inShape(1);
    //resize breaks  references...so out better have the right shape 
    //if references is not to be broken
    if(out.nelements()==0)
      out.resize(nxx,nyy,nzz);
    Bool deleteIn,deleteOut, delFlag;
    const Complex* pin = in.getStorage(deleteIn);
    const Bool* poutflag= outFlag.getStorage(delFlag);
    Complex* pout = out.getStorage(deleteOut);
    uInt i=0, zOffset=0;
    for (uInt iz=0; iz<nzz; ++iz, zOffset+=nxx) {
      Int yOffset=zOffset;
      for (uInt iy=0; iy<nyy; ++iy, yOffset+=nxx*nzz) {
	for (uInt ix=0; ix<nxx; ++ix){ 
	  if(!poutflag[i])
	    pout[i] = pin[ix+yOffset];
	  ++i;
	}
      }
    }
    out.putStorage(pout,deleteOut);
    in.freeStorage(pin,deleteIn);
    outFlag.freeStorage(poutflag, delFlag);
  }

  // helper function to swap the y and z axes of a Cube
  void FTMachine::swapyz(Cube<Bool>& out, const Cube<Bool>& in)
  {
    IPosition inShape=in.shape();
    uInt nxx=inShape(0),nyy=inShape(2),nzz=inShape(1);
    if(out.nelements()==0)
      out.resize(nxx,nyy,nzz);
    Bool deleteIn,deleteOut;
    const Bool* pin = in.getStorage(deleteIn);
    Bool* pout = out.getStorage(deleteOut);
    uInt i=0, zOffset=0;
    for (uInt iz=0; iz<nzz; iz++, zOffset+=nxx) {
      Int yOffset=zOffset;
      for (uInt iy=0; iy<nyy; iy++, yOffset+=nxx*nzz) {
	for (uInt ix=0; ix<nxx; ix++) pout[i++] = pin[ix+yOffset];
      }
    }
    out.putStorage(pout,deleteOut);
    in.freeStorage(pin,deleteIn);
  }
  
  void FTMachine::setPointingDirColumn(const String& column){
    pointingDirCol_p=column;
    pointingDirCol_p.upcase();
    if( (pointingDirCol_p != "DIRECTION") &&(pointingDirCol_p != "TARGET") && (pointingDirCol_p != "ENCODER") && (pointingDirCol_p != "POINTING_OFFSET") && (pointingDirCol_p != "SOURCE_OFFSET")){
      
      //basically at this stage you don't know what you're doing...so you get the default
      
      pointingDirCol_p="DIRECTION";
      
    }    
  }
  
  String FTMachine::getPointingDirColumnInUse(){
    
    return pointingDirCol_p;
    
  }
  
  void FTMachine::setSpwChanSelection(const Cube<Int>& spwchansels) {
    spwChanSelFlag_p.resize();
    spwChanSelFlag_p=spwchansels;
  }
  
  void FTMachine::setSpwFreqSelection(const Matrix<Double>& spwFreqs) 
  {
    spwFreqSel_p.assign(spwFreqs);
    SynthesisUtils::expandFreqSelection(spwFreqs,expandedSpwFreqSel_p, expandedSpwConjFreqSel_p);
  }

  void FTMachine::setSpectralFlag(const vi::VisBuffer2& vb, Cube<Bool>& modflagcube){

      modflagcube.resize(vb.flagCube().shape());
      modflagcube=vb.flagCube();
      if(!isPseudoI_p){
	ArrayIterator<Bool> from(modflagcube, IPosition(1, 0));
	while(!from.pastEnd()){
	  if(anyTrue(from.array()))
	    from.array().set(True);
	  from.next();
	}
      }
      uInt nchan = vb.nChannels();
      uInt msid = vb.msId();
      uInt selspw = vb.spectralWindows()[0];
      Bool spwFlagIsSet=( (spwChanSelFlag_p.nelements() > 0) && (spwChanSelFlag_p.shape()(1) > selspw) &&
  			(spwChanSelFlag_p.shape()(0) > msid) &&
  			(spwChanSelFlag_p.shape()(2) >=nchan));
      //cerr << "spwFlagIsSet " << spwFlagIsSet << endl;
      for (uInt i=0;i<nchan;i++) {
        //Flag those channels that  did not get selected...
        //respect the flags from vb  if selected  or
        //if spwChanSelFlag is wrong shape
        if ((spwFlagIsSet) && (spwChanSelFlag_p(msid,selspw,i)!=1)) {
  	modflagcube.xzPlane(i).set(true);
        }
      }
    }

  Matrix<Double> FTMachine::negateUV(const vi::VisBuffer2& vb){
    Matrix<Double> uvw(vb.uvw().shape());
    for (Int i=0;i< vb.nRows() ; ++i) {
      for (Int idim=0;idim<2; ++idim) uvw(idim,i)=-vb.uvw()(idim, i);
      uvw(2,i)=vb.uvw()(2,i);
    }
    return uvw;
  }
  //-----------------------------------------------------------------------------------------------------------------
  //------------  Vectorized versions of initializeToVis, initializeToSky, finalizeToSky  
  //------------  that are called from CubeSkyEquation.
  //------------  They call getImage,getWeightImage, which are implemented in all FTMs
  //------------  Also, Correlation / Stokes conversions and gS/ggS normalizations.
 

  void FTMachine::setSkyJones(Vector<CountedPtr<casa::refim::SkyJones> >& sj){
    sj_p.resize();
    sj_p=sj;
    cout << "FTM : Set Sky Jones of length " << sj_p.nelements() << " and types ";
    for( uInt i=0;i<sj_p.nelements();i++) cout << sj_p[i]->type() << " ";
    cout << endl;
  }
  // Convert complex correlation planes to float Stokes planes
  void FTMachine::correlationToStokes(ImageInterface<Complex>& compImage, 
				      ImageInterface<Float>& resImage, 
				      const Bool dopsf)
  {
    // Convert correlation image to IQUV format
    AlwaysAssert(compImage.shape()[0]==resImage.shape()[0], AipsError);
    AlwaysAssert(compImage.shape()[1]==resImage.shape()[1], AipsError);
    AlwaysAssert(compImage.shape()[3]==resImage.shape()[3], AipsError);
    
    if(dopsf) 
      { 
	// For the PSF, choose only those stokes planes that have a valid PSF
	StokesImageUtil::ToStokesPSF(resImage,compImage);
      }
    else 
      {
	StokesImageUtil::To(resImage,compImage);
      }
  };
  
  // Convert float Stokes planes to complex correlation planes
  void FTMachine::stokesToCorrelation(ImageInterface<Float>& modelImage,
				      ImageInterface<Complex>& compImage)
  {
    /*
    StokesCoordinate stcomp=compImage.coordinates().stokesCoordinate(compImage.coordinates().findCoordinate(Coordinate::STOKES));
    StokesCoordinate stfloat = modelImage.coordinates().stokesCoordinate(modelImage.coordinates().findCoordinate(Coordinate::STOKES));

    cout << "Stokes types : complex : " << stcomp.stokes() << "    float : " << stfloat.stokes() << endl;
    cout << "Shapes : complex : " << compImage.shape() << "   float : " << modelImage.shape() << endl;
    */

    //Pol axis need not be same
    AlwaysAssert(modelImage.shape()[0]==compImage.shape()[0], AipsError);
    AlwaysAssert(modelImage.shape()[1]==compImage.shape()[1], AipsError);
    AlwaysAssert(modelImage.shape()[3]==compImage.shape()[3], AipsError);
    // Convert from Stokes to Complex
    StokesImageUtil::From(compImage, modelImage);
  };
  
  //------------------------------------------------------------------------------------------------------------------
  
  void FTMachine::normalizeImage(ImageInterface<Float>& skyImage,
				 Matrix<Float>& sumOfWts,
				 ImageInterface<Float>& sensitivityImage,
				 Bool dopsf, Float pblimit, Int normtype)
  {
    
    //Normalize the sky Image
    Int nXX=(skyImage).shape()(0);
    Int nYY=(skyImage).shape()(1);
    Int npola= (skyImage).shape()(2);
    Int nchana= (skyImage).shape()(3);
    
      IPosition pcentre(4,nXX/2,nYY/2,0,0);
    // IPosition psource(4,nXX/2+22,nYY/2,0,0);
    
    //    storeImg(String("norm_resimage.im") , skyImage);
    //    storeImg(String("norm_sensitivity.im"), sensitivityImage);
   
      /////    cout << "FTM::norm : pblimit : " << pblimit << endl; 
    
    // Note : This is needed because initial prediction has no info about sumwt.
    // Not a clean solution.  // ForSB -- if you see a better way, go for it.
    if(sumOfWts.shape() != IPosition(2,npola,nchana))
      {
	cout << "Empty SumWt.... resizing and filling with 1.0" << endl;
	sumOfWts.resize(IPosition(2,npola,nchana));
	sumOfWts=1.0;
      }
    
    //    if(dopsf)  cout << "*** FTM::normalizeImage : Image Center : " << skyImage.getAt(pcentre) << "  and weightImage : " << sensitivityImage.getAt(pcentre) << "  SumWt : " << sumOfWts[0,0];  
    //    else  cout << "*** FTM::normalizeImage : Source Loc : " << skyImage.getAt(psource) << "  and weightImage : " << sensitivityImage.getAt(psource) << "  SumWt : " << sumOfWts[0,0]; 
    
    
    
    IPosition blc(4,nXX, nYY, npola, nchana);
    IPosition trc(4, nXX, nYY, npola, nchana);
    blc(0)=0; blc(1)=0; trc(0)=nXX-1; trc(1)=nYY-1;
    //max weights per plane
    for (Int pol=0; pol < npola; ++pol){
      for (Int chan=0; chan < nchana ; ++chan){
	
	blc(2)=pol; trc(2)=pol;
	blc(3)=chan; trc(3)=chan;
	Slicer sl(blc, trc, Slicer::endIsLast);
	SubImage<Float> subSkyImage(skyImage, sl, false);
	SubImage<Float> subSensitivityImage(sensitivityImage, sl, false);
	SubImage<Float> subOutput(skyImage, sl, true);
	Float sumWt = sumOfWts(pol,chan);
	if(sumWt > 0.0){
	  switch(normtype)
	    {
	    case 0: // only sum Of Weights - FTM only (ForSB)
	      subOutput.copyData( (LatticeExpr<Float>) 
				  ((dopsf?1.0:-1.0)*subSkyImage/(sumWt)));
	      break;
	      
	    case 1: // only sensitivityImage   Ic/avgPB  (ForSB)
	      subOutput.copyData( (LatticeExpr<Float>) 
				  (iif(subSensitivityImage > (pblimit), 
				       (subSkyImage/(subSensitivityImage)),
				       (subSkyImage))));
				       // 0.0)));
	      break;
	      
	    case 2: // sum of Weights and sensitivityImage  IGridded/(SoW*avgPB) and PSF --> Id (ForSB)
	      subOutput.copyData( (LatticeExpr<Float>) 
				  (iif(subSensitivityImage > (pblimit), 
				       ((dopsf?1.0:-1.0)*subSkyImage/(subSensitivityImage*sumWt)),
				       //((dopsf?1.0:-1.0)*subSkyImage))));
				       0.0)));		
	      break;
	      
  	    case 3: // MULTIPLY by the sensitivityImage  avgPB
	      subOutput.copyData( (LatticeExpr<Float>) (subSkyImage * subSensitivityImage) );
	      break;
	      
  	    case 4: // DIVIDE by sqrt of sensitivityImage 
	      subOutput.copyData( (LatticeExpr<Float>) 
				  (iif((subSensitivityImage) > (pblimit), 
				       (subSkyImage/(sqrt(subSensitivityImage))),
				       (subSkyImage))));
				       //0.0)));
	      break;
	      
  	    case 5: // MULTIPLY by sqrt of sensitivityImage 
	      subOutput.copyData( (LatticeExpr<Float>) 
				  (iif((subSensitivityImage) > (pblimit), 
				       (subSkyImage * (sqrt(subSensitivityImage))),
				       (subSkyImage))));

	      break;

	    case 6: // divide by non normalized sensitivity image
	      {
		Float elpblimit=max( subSensitivityImage).getFloat() * pblimit;
		subOutput.copyData( (LatticeExpr<Float>) 
				    (iif(subSensitivityImage > (elpblimit), 
					 ((dopsf?1.0:-1.0)*subSkyImage/(subSensitivityImage)),
					 0.0)));
	      }
	      break;
	    default:
	      throw(AipsError("Unrecognized norm-type in FTM::normalizeImage : "+String::toString(normtype)));
	    }
	}
	else{
	  subOutput.set(0.0);
	}
      }
    }
    
    //if(dopsf)  cout << " Normalized (" << normtype << ") Image Center : " << skyImage.getAt(pcentre) << endl; 
     //     else  cout << " Normalized (" << normtype << ") Source Loc : " << skyImage.getAt(psource) << endl; 
    
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  /////  For use with the new framework 
  ///// (Sorry about these copies, but need to keep old system working)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  // Vectorized InitializeToVis
  void FTMachine::initializeToVisNew(const VisBuffer2& vb,
				     CountedPtr<SIImageStore> imstore)
  {
    AlwaysAssert(imstore->getNTaylorTerms(false)==1, AipsError);

    Matrix<Float> tempWts;

    // Convert from Stokes planes to Correlation planes
    stokesToCorrelation(*(imstore->model()), *(imstore->forwardGrid()));

    if(vb.polarizationFrame()==MSIter::Linear) {
      StokesImageUtil::changeCStokesRep(*(imstore->forwardGrid()),
					StokesImageUtil::LINEAR);
    }
    else {
      StokesImageUtil::changeCStokesRep(*(imstore->forwardGrid()),
					StokesImageUtil::CIRCULAR);
    }
   
    //------------------------------------------------------------------------------------
    // Image Mosaic only :  Multiply the input model with the Primary Beam
    if(sj_p.nelements() >0 ){
      for (uInt k=0; k < sj_p.nelements(); ++k){
	(sj_p(k))->apply(*(imstore->forwardGrid()), *(imstore->forwardGrid()), vb, 0, true);
      }
    }
    //------------------------------------------------------------------------------------

    // Call initializeToVis
    initializeToVis(*(imstore->forwardGrid()), vb); // Pure virtual
    
  };
  
  // Vectorized finalizeToVis is not implemented because it does nothing and is never called.
  
  // Vectorized InitializeToSky
  void FTMachine::initializeToSkyNew(const Bool dopsf, 
				     const VisBuffer2& vb,
				     CountedPtr<SIImageStore> imstore)
    
  {
    AlwaysAssert(imstore->getNTaylorTerms(false)==1, AipsError);
    
    // Make the relevant float grid. 
    // This is needed mainly for facetting (to set facet shapes), but is harmless for non-facetting.
    if( dopsf ) { imstore->psf(); } else { imstore->residual(); } 

    // Initialize the complex grid (i.e. tell FTMachine what array to use internally)
    Matrix<Float> sumWeight;
    initializeToSky(*(imstore->backwardGrid()) , sumWeight , vb);

  };
  
  // Vectorized finalizeToSky
  void FTMachine::finalizeToSkyNew(Bool dopsf, 
				   const VisBuffer2& vb,
				   CountedPtr<SIImageStore> imstore  )				   
  {
    // Check vector lengths. 
    AlwaysAssert( imstore->getNTaylorTerms(false)==1, AipsError);

    Matrix<Float> sumWeights;
    finalizeToSky(); 

    //------------------------------------------------------------------------------------
    // Straightforward case. No extra primary beams. No image mosaic
    if(sj_p.nelements() == 0 ) 
      {
	correlationToStokes( getImage(sumWeights, false) , ( dopsf ? *(imstore->psf()) : *(imstore->residual()) ), dopsf);
	
	if( useWeightImage() && dopsf ) { 
	  getWeightImage( *(imstore->weight())  , sumWeights); 
	  // Fill weight image only once, during PSF generation. Remember.... it is normalized only once
	  // during PSF generation.
	}
	
	// Take sumWeights from corrToStokes here....
	Matrix<Float> sumWeightStokes( (imstore->sumwt())->shape()[2], (imstore->sumwt())->shape()[3]   );
	StokesImageUtil::ToStokesSumWt( sumWeightStokes, sumWeights );

	AlwaysAssert( ( (imstore->sumwt())->shape()[2] == sumWeightStokes.shape()[0] ) && 
		      ((imstore->sumwt())->shape()[3] == sumWeightStokes.shape()[1] ) , AipsError );

	(imstore->sumwt())->put( sumWeightStokes.reform((imstore->sumwt())->shape()) );
	
      }
    //------------------------------------------------------------------------------------
    // Image Mosaic only :  Multiply the residual, and weight image by the PB.
    else 
      {
      
      // Take the FT of the gridded values. Writes into backwardGrid(). 
      getImage(sumWeights, false);

      // Multiply complex image grid by PB.
      if( !dopsf )
	{
	  for (uInt k=0; k < sj_p.nelements(); ++k){
	    (sj_p(k))->apply(*(imstore->backwardGrid()), *(imstore->backwardGrid()), vb, 0, true);
	  }
	}

      // Convert from correlation to Stokes onto a new temporary grid.
      SubImage<Float>  targetImage( ( dopsf ? *(imstore->psf()) : *(imstore->residual()) ) , true);
      TempImage<Float> temp( targetImage.shape(), targetImage.coordinates() );
      correlationToStokes( *(imstore->backwardGrid()), temp, false);

      // Add the temporary Stokes image to the residual or PSF, whichever is being made.
      LatticeExpr<Float> addToRes( targetImage + temp );
      targetImage.copyData(addToRes);
      
      // Now, do the same with the weight image and sumwt ( only on the first pass )
      if( dopsf )
	{
	  SubImage<Float>  weightImage(  *(imstore->weight()) , true);
	  TempImage<Float> temp(weightImage.shape(), weightImage.coordinates());
	  getWeightImage(temp, sumWeights);

	  for (uInt k=0; k < sj_p.nelements(); ++k){
	    (sj_p(k))->applySquare(temp,temp, vb, -1);
	  }

	  LatticeExpr<Float> addToWgt( weightImage + temp );
	  weightImage.copyData(addToWgt);
	  
	  AlwaysAssert( ( (imstore->sumwt())->shape()[2] == sumWeights.shape()[0] ) && 
			((imstore->sumwt())->shape()[3] == sumWeights.shape()[1] ) , AipsError );

	  SubImage<Float>  sumwtImage(  *(imstore->sumwt()) , true);
	  TempImage<Float> temp2(sumwtImage.shape(), sumwtImage.coordinates());
	  temp2.put( sumWeights.reform(sumwtImage.shape()) );
	  LatticeExpr<Float> addToWgt2( sumwtImage + temp2 );
	  sumwtImage.copyData(addToWgt2);
	  
	  //cout << "In finalizeGridCoreMos : sumwt : " << sumwtImage.get() << endl;
	  
	}

    }
    //------------------------------------------------------------------------------------


    
    return;
  };

  Bool FTMachine::changedSkyJonesLogic(const vi::VisBuffer2& vb, Bool& firstRow, Bool& internalRow)
  {
    firstRow=false;
    internalRow=false;

    if( sj_p.nelements()==0 ) 
      {throw(AipsError("Internal Error : Checking changedSkyJones, but it is not yet set."));}

    CountedPtr<SkyJones> ej = sj_p[0];
    if(ej.null())
      return false;
    if(ej->changed(vb,0))
      firstRow=true;
    Int row2temp=0;
    if(ej->changedBuffer(vb,0,row2temp)) {
      internalRow=true;
    }
    return (firstRow || internalRow) ;
  }
  void FTMachine::setCFCache(CountedPtr<CFCache>& /*cfc*/, const Bool /*loadCFC*/) 
  {
    throw(AipsError("FTMachine::setCFCache() directly called!"));
  }
  


  /*
  /// Move to individual FTMs............ make it pure virtual.
  Bool FTMachine::useWeightImage()
  {
    if( name() == "GridFT" || name() == "WProjectFT" )  
      { return false; }
    else
      { return true; }
  }
  */

  }//# namespace refim ends
}//namespace CASA ends

