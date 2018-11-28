//# HetArrayConvFunc.cc: Implementation for HetArrayConvFunc
//# Copyright (C) 2008-2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  General Public
//# License for more details.
//#
//# You should have received a copy of the GNU  General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Slice.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Containers/SimOrdMap.h>
#include <scimath/Mathematics/FFTServer.h>
#include <measures/Measures/MeasTable.h>
#include <scimath/Mathematics/MathFunc.h>
#include <scimath/Mathematics/ConvolveGridder.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/CompositeNumber.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/SubImage.h>
#include <images/Images/TempImage.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>

#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/Lattices/SubLattice.h>
#include <lattices/LRegions/LCBox.h>
#include <lattices/Lattices/LatticeConcat.h>
#include <lattices/LEL/LatticeExpr.h>
#include <lattices/Lattices/LatticeCache.h>
#include <lattices/LatticeMath/LatticeFFT.h>


#include <ms/MeasurementSets/MSColumns.h>

#include <msvis/MSVis/VisBuffer2.h>

#include <synthesis/TransformMachines2/Utils.h>
#include <synthesis/TransformMachines/PBMath1DAiry.h>
#include <synthesis/TransformMachines/PBMath1DNumeric.h>
#include <synthesis/TransformMachines/PBMath2DImage.h>
#include <synthesis/TransformMachines/PBMath.h>
#include <synthesis/TransformMachines2/HetArrayConvFunc.h>
#include <synthesis/MeasurementEquations/VPManager.h>

#include <casa/OS/Timer.h>
#include <iomanip>


using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

namespace refim {

using namespace casacore;
using namespace casa;
using namespace casacore;
using namespace casa::refim;

typedef unsigned long long ooLong;


  HetArrayConvFunc::HetArrayConvFunc() : convFunctionMap_p(0), nDefined_p(0), antDiam2IndexMap_p(-1),msId_p(-1), actualConvIndex_p(-1), vpTable_p("")
{
    calcFluxScale_p=true;
    init(PBMathInterface::AIRY);
}

HetArrayConvFunc::HetArrayConvFunc(const PBMathInterface::PBClass typeToUse, const String vpTable):
  convFunctionMap_p(0), nDefined_p(0), antDiam2IndexMap_p(-1),msId_p(-1), actualConvIndex_p(-1), vpTable_p(vpTable)
{
    calcFluxScale_p=true;
    init(typeToUse);

}

  HetArrayConvFunc::HetArrayConvFunc(const RecordInterface& rec, Bool calcFluxneeded):convFunctionMap_p(0), nDefined_p(0), antDiam2IndexMap_p(-1),msId_p(-1), actualConvIndex_p(-1) {
    String err;
    fromRecord(err, rec, calcFluxneeded);
}

HetArrayConvFunc::~HetArrayConvFunc() {
    //
}

void HetArrayConvFunc::init(const PBMathInterface::PBClass typeTouse) {
    doneMainConv_p=false;
    filledFluxScale_p=false;
    usePointingTable_p=False;
    pbClass_p=typeTouse;
    ///Generate the sincCache now as inside the multithread  part it may have race issues
    initSincCache();
    timer1_p=0.0;
    timer2_p=0.0;
    timer3_p=0.0;
    
      
}



void HetArrayConvFunc::findAntennaSizes(const vi::VisBuffer2& vb) {

      if(msId_p != vb.msId()) {
        msId_p=vb.msId();
        //ROMSColumns mscol(vb.ms());
        const ROMSAntennaColumns& ac=vb.subtableColumns().antenna();
        antIndexToDiamIndex_p.resize(ac.nrow());
        antIndexToDiamIndex_p.set(-1);
        Int diamIndex=antDiam2IndexMap_p.ndefined();
        Vector<Double> dishDiam=ac.dishDiameter().getColumn();
        Vector<String>dishName=ac.name().getColumn();
        String telescop=vb.subtableColumns().observation().telescopeName()(0);
        PBMath::CommonPB whichPB;
        if(pbClass_p==PBMathInterface::COMMONPB) {
            String band;
            String commonPBName;
            // This frequency is ONLY required to determine which PB model to use:
            // The VLA, the ATNF, and WSRT have frequency - dependent PB models
            Quantity freq( vb.subtableColumns().spectralWindow().refFrequency()(0), "Hz");


            PBMath::whichCommonPBtoUse( telescop, freq, band, whichPB, commonPBName );
            //Revert to using AIRY for unknown common telescope
            if(whichPB==PBMath::UNKNOWN)
                pbClass_p=PBMathInterface::AIRY;

        }
        if(pbClass_p== PBMathInterface::AIRY) {
	  LogIO os;
	os << LogOrigin("HetArrConvFunc", "findAntennaSizes")  << LogIO::NORMAL;
            ////////We'll be using dish diameter as key
            for (uInt k=0; k < dishDiam.nelements(); ++k) {
                if((diamIndex !=0) && antDiam2IndexMap_p.isDefined(String::toString(dishDiam(k)))) {
                    antIndexToDiamIndex_p(k)=antDiam2IndexMap_p(String::toString(dishDiam(k)));
                }
                else {
                    if(dishDiam[k] > 0.0) { //there may be stations with no dish on
                        antDiam2IndexMap_p.define(String::toString(dishDiam(k)), diamIndex);
                        antIndexToDiamIndex_p(k)=diamIndex;
                        antMath_p.resize(diamIndex+1);
                        if(pbClass_p== PBMathInterface::AIRY) {
			  //ALMA ratio of blockage to dish
                            Quantity qdiam= Quantity (dishDiam(k),"m");	
                            Quantity blockDiam= Quantity(dishDiam(k)/12.0*.75, "m");
			    Quantity support=Quantity(150, "arcsec");
                            ///For ALMA 12m dish it is effectively 10.7 m according to Todd Hunter
                            ///@ 2011-12-06
                            if((vb.subtableColumns().observation().telescopeName()(0) =="ALMA") || (vb.subtableColumns().observation().telescopeName()(0) =="ACA")){
			      Quantity fov(max(nx_p*dc_p.increment()(0), ny_p*dc_p.increment()(1)), dc_p.worldAxisUnits()(0));
			      if((abs(dishDiam[k] - 12.0) < 0.5)) {
				qdiam= Quantity(10.7, "m");
				blockDiam= Quantity(0.75, "m");
				support=Quantity(max(150.0, fov.getValue("arcsec")/5.0), "arcsec");
                                
			      }
			      else{
                                //2017 the ACA dishes are best represented by 6.25m:
                               
				qdiam= Quantity(6.25,"m");
				blockDiam = Quantity(0.75,"m");
				support=Quantity(max(300.0,fov.getValue("arcsec")/3.0) , "arcsec");
			      }
			    }
			     os << "Overriding PB with Airy of diam,blockage="<<qdiam<<","<<blockDiam<<" starting with antenna "<<k<<LogIO::POST; 
			    
			
		    

			antMath_p[diamIndex]=new PBMath1DAiry(qdiam, blockDiam,
							  support,
							  Quantity(100.0,"GHz"));
			}

		

                        //////Will no longer support this
                        /*else if(pbClass_p== PBMathInterface::IMAGE){
                          //Get the image name by calling code for the antenna name and array name
                          //For now hard wired to ALMA as this part of the code will not be accessed for non-ALMA
                          //see Imager::setMosaicFTMachine
                          // When ready to generalize then code that calls with telescope name, antenna name
                          //(via vb.msColumns) and/or diameter and frequency via vb.frequency (indexing will need to
                          //be upgraded to account for frequency too) should be done to return the
                          //right voltage pattern image.
                          String vpImageName="";
                          if (abs(dishDiam[k]-7.0) < 1.0)
                        Aipsrc::find(vpImageName, "alma.vp.7m", "");
                          else
                        Aipsrc::find(vpImageName, "alma.vp.12m", "") ;
                          //cerr << "first vpImagename " << vpImageName  << endl;
                          if(vpImageName==""){
                        String beamPath;
                        if(!MeasTable::AntennaResponsesPath(beamPath, "ALMA")){
                          throw(AipsError("Alma beam images requested cannot be found "));
                        }
                        else{
                          beamPath=beamPath.before(String("AntennaResponses"));
                          vpImageName= (abs(dishDiam[k]-7.0) < 1.0) ? beamPath
                            +String("/ALMA_AIRY_7M.VP") :
                            beamPath+String("/ALMA_AIRY_12M.VP");
                        }


                          }
                          //cerr << "Using the image VPs " << vpImageName << endl;
                          if(Table::isReadable(vpImageName))
                        antMath_p[diamIndex]=new PBMath2DImage(PagedImage<Complex>(vpImageName));
                          else
                        throw(AipsError(String("Cannot find voltage pattern image ") + vpImageName));
                        }
                        else{

                          throw(AipsError("Do not  deal with non airy dishes or images of VP yet "));
                        }
                        */
                        ++diamIndex;
	    
		    }

		}
	    }

	}
        else if(pbClass_p== PBMathInterface::IMAGE) {

            VPManager *vpman=VPManager::Instance();
            if(vpTable_p != String(""))
                vpman->loadfromtable(vpTable_p);
            ///else it is already loaded in the static object
            Vector<Record> recs;
            Vector<Vector<String> > antnames;

            if(vpman->imagepbinfo(antnames, recs)) {
                Vector<Bool> dishDefined(dishName.nelements(), false);
                Int nbeams=antnames.nelements();
                ///will be keying on file image file name here
                for (uInt k=0; k < dishDiam.nelements(); ++k) {
                    String key;
                    Bool beamDone=false;
                    Int recordToUse=0;
                    for (Int j =0; j < nbeams; ++j) {
                        key=recs[j].isDefined("realimage") ? recs[j].asString("realimage") : recs[j].asString("compleximage");
                        if(antnames[j][0]=="*" || anyEQ(dishName[k], antnames[j])) {
                            dishDefined[k]=true;
                            recordToUse=j;

                            if((diamIndex !=0) && antDiam2IndexMap_p.isDefined(key)) {
                                antIndexToDiamIndex_p(k)=antDiam2IndexMap_p(key);
                                beamDone=true;
                            }
                        }
                    }
                    if(!beamDone && dishDefined[k]) {
                        key=recs[recordToUse].isDefined("realimage") ? recs[recordToUse].asString("realimage") : recs[recordToUse].asString("compleximage");
                        antDiam2IndexMap_p.define(key, diamIndex);
                        antIndexToDiamIndex_p(k)=diamIndex;
                        antMath_p.resize(diamIndex+1);
                        if(recs[recordToUse].isDefined("realimage") && recs[recordToUse].isDefined("imagimage")) {
                            PagedImage<Float> realim(recs[recordToUse].asString("realimage"));
                            PagedImage<Float> imagim(recs[recordToUse].asString("imagim"));
                            antMath_p[diamIndex]=new PBMath2DImage(realim, imagim);
                        }
                        else {
                            antMath_p[diamIndex]=new PBMath2DImage(PagedImage<Complex>(recs[recordToUse].asString("compleximage")));
                        }
                        ++diamIndex;
                    }
                }
                if(!allTrue(dishDefined)) {
                    //cerr << "dishDefined " << dishDefined << endl;
                    throw(AipsError("Some Antennas in the MS did not have a VP defined"));
                }
            }
            else {
                throw(AipsError("Mosaic does not support non-image voltage patterns yet"));
            }

            //Get rid of the static class
            vpman->reset();
        }
	else if(vpTable_p != String("")){
	  ////When we get vpmanager to give beams on antenna name we
	  //should change this key to antenna name and loop over all antenna names
	  if((diamIndex !=0) && antDiam2IndexMap_p.isDefined(telescop+String("_")+String::toString(dishDiam(0)))) {
	    antIndexToDiamIndex_p.set(antDiam2IndexMap_p(telescop+String("_")+String::toString(dishDiam(0))));
	  }
	  else{
	    antDiam2IndexMap_p.define(telescop+"_"+String::toString(dishDiam(0)), diamIndex);
	    antIndexToDiamIndex_p.set(diamIndex);
	     VPManager *vpman=VPManager::Instance();
	     vpman->loadfromtable(vpTable_p);
	     Record rec;
	     vpman->getvp(rec, telescop);
	     antMath_p.resize(diamIndex+1);
	     antMath_p[diamIndex]=PBMath::pbMathInterfaceFromRecord(rec);
	     vpman->reset();
	  }
	  
	}
        else if(pbClass_p==PBMathInterface::COMMONPB) {
            //cerr << "Doing the commonPB thing" << endl;
            ///Have to use telescop part as string as in multims case different
            //telescopes may have same dish size but different commonpb
            // VLA and EVLA for e.g.
            if((diamIndex !=0) && antDiam2IndexMap_p.isDefined(telescop+String("_")+String::toString(dishDiam(0)))) {
                antIndexToDiamIndex_p.set(antDiam2IndexMap_p(telescop+String("_")+String::toString(dishDiam(0))));
            }
            else {
                antDiam2IndexMap_p.define(telescop+"_"+String::toString(dishDiam(0)), diamIndex);
                antIndexToDiamIndex_p.set(diamIndex);
                antMath_p.resize(diamIndex+1);
                antMath_p[diamIndex]=PBMath::pbMathInterfaceForCommonPB(whichPB, True);
            }
        }
        else {

            throw(AipsError("Mosaic  supports image based or Airy voltage patterns or known common pb   for now"));

        }

        //cerr << "antIndexTodiamIndex " << antIndexToDiamIndex_p << endl;
    }





}

void  HetArrayConvFunc::reset() {
    doneMainConv_p=false;
    convFunctions_p.resize(0, true);
    convWeights_p.resize(0, true);
    convSizes_p.resize(0, true);
    convSupportBlock_p.resize(0, true);
    convFunctionMap_p.resize(0);
    vbConvIndex_p.clear();
}


  void HetArrayConvFunc::findConvFunction2(const ImageInterface<Complex>& iimage,
                                        const vi::VisBuffer2& vb,
                                        const Int& convSamp, const Vector<Double>& visFreq,
                                        Array<Complex>& convFunc,
                                        Array<Complex>& weightConvFunc,
                                        Vector<Int>& convsize,
                                        Vector<Int>& convSupport,
                                        Vector<Int>& convFuncPolMap,
                                        Vector<Int>& convFuncChanMap,
                                        Vector<Int>& convFuncRowMap, Bool getConjConvFunc,
					const MVDirection& extraShift, const Bool useExtraShift)
  {
    //Double wtime0=omp_get_wtime();
    storeImageParams(iimage,vb);
    convFuncChanMap.resize(vb.nChannels());
    Vector<Double> beamFreqs;
    findUsefulChannels(convFuncChanMap, beamFreqs, vb, visFreq);
    Int nBeamChans=beamFreqs.nelements();
    /////For now not doing beam rotation or squints but to be enabled easily
    convFuncPolMap.resize(vb.nCorrelations());
    Int nBeamPols=1;
    convFuncPolMap.set(0);
    findAntennaSizes(vb);
    uInt ndish=antMath_p.nelements();
    Int isCached=checkPBOfField(vb, convFuncRowMap, extraShift, useExtraShift);
     if(isCached ==2) {
       //nominally pointings out of the image
        convFunc.resize();
        weightConvFunc.resize();
        convsize.resize();
        convSupport.resize();
        convFuncRowMap.resize();
        return;
    }
     actualConvIndex_p=convIndex(vb);
     //Double wtime0=omp_get_wtime();
     if(doneMainConv_p.shape()[0] < (actualConvIndex_p+1)) {
        //cerr << "resizing DONEMAIN " <<   doneMainConv_p.shape()[0] << endl;
        doneMainConv_p.resize(actualConvIndex_p+1, true);
        doneMainConv_p[actualConvIndex_p]=false;
        vpFunctions_p.resize(actualConvIndex_p+1);
        vpFunctions_p[actualConvIndex_p]=nullptr;
	pbFunctions_p.resize(actualConvIndex_p+1);
        pbFunctions_p[actualConvIndex_p]=nullptr;
	convFunctions_p.resize(actualConvIndex_p+1);
	convWeights_p.resize(actualConvIndex_p+1);
	convSupportBlock_p.resize(actualConvIndex_p+1);
	convFunctions_p[actualConvIndex_p]=nullptr;
	convWeights_p[actualConvIndex_p]=nullptr;
	convSupportBlock_p[actualConvIndex_p]=nullptr;
	convSizes_p.resize(actualConvIndex_p+1);
	convSizes_p[actualConvIndex_p]=nullptr;
	origConvSize_p.resize(actualConvIndex_p+1, True);
     }
     
       // Get the coordinate system
    CoordinateSystem coords(iimage.coordinates());
    Int directionIndex=coords.findCoordinate(Coordinate::DIRECTION);
    AlwaysAssert(directionIndex>=0, AipsError);
    // Set up the convolution function.
    Int nx=nx_p;
    Int ny=ny_p;
    Int support=max(nx_p, ny_p)/10;
    Int convSampling=1;
    //timer1_p+=omp_get_wtime()-wtime0;
    //wtime0=omp_get_wtime();
    if(!doneMainConv_p[actualConvIndex_p]) {
        for (uInt ii=0; ii < ndish; ++ii) {
            support=max((antMath_p[ii])->support(coords), support);
        }
	
        support=Int(min(Float(support), max(Float(nx_p), Float(ny_p)))*2.0)/2;
	///Multiply by 2 just to make sure in-baseline pointing can drift ...
	////If it is more than 2 HPFW that data is useless anyways.
        convSize_p=usePointingTable_p ? 2*support*convSampling : support*convSampling ;
        // Make this a nice composite number, to speed up FFTs
        CompositeNumber cn(uInt(convSize_p*2.0));
        convSize_p  = cn.nextLargerEven(Int(convSize_p));
        convSize_p=(convSize_p/16)*16;  // need it to be divisible by 8 in places
	origConvSize_p[actualConvIndex_p]=convSize_p;
	
	Vector<Double> sampling;
	DirectionCoordinate dc=dc_p;
       sampling = dc.increment();
       sampling*=Double(convSampling);
       sampling(0)*=Double(nx)/Double(convSize_p);
       sampling(1)*=Double(ny)/Double(convSize_p);
       dc.setIncrement(sampling);

        Vector<Double> unitVec(2);
        unitVec=convSize_p/2;
        dc.setReferencePixel(unitVec);
	MDirection fieldDir=vbutil_p->getPhaseCenter(vb);
        //make sure we are using the same units
        fieldDir.set(dc.worldAxisUnits()(0));
	
        dc.setReferenceValue(fieldDir.getAngle().getValue());
        coords.replaceCoordinate(dc, directionIndex);
	Int spind=coords.findCoordinate(Coordinate::SPECTRAL);
        SpectralCoordinate spCoord=coords.spectralCoordinate(spind);
        spCoord.setReferencePixel(Vector<Double>(1,0.0));
        spCoord.setReferenceValue(Vector<Double>(1, beamFreqs(0)));
        if(beamFreqs.nelements() >1)
	  spCoord.setIncrement(Vector<Double>(1, beamFreqs(1)-beamFreqs(0)));
	//cerr << "spcoord " ;
	//spCoord.print(std::cerr);
        coords.replaceCoordinate(spCoord, spind);
	CoordinateSystem conjCoord=coords;
	Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
	SpectralCoordinate conjSpCoord=spCoord;
		//cerr << "centreFreq " << centerFreq << " beamFreqs " << beamFreqs(0) << "  " << beamFreqs(1) << endl;
	conjSpCoord.setReferenceValue(Vector<Double>(1,SynthesisUtils::conjFreq(beamFreqs[0], centerFreq)));
	///Increment should go in the reverse direction
	////Do a tabular spectral coordinate for more than 1 channel 
	if(beamFreqs.nelements() >1){
	  Vector<Double> conjFreqs(beamFreqs.nelements());
	  for (uInt kk=0; kk< beamFreqs.nelements(); ++kk){
	    //conjFreqs[kk]=sqrt(2*centerFreq*centerFreq-beamFreqs(kk)*beamFreqs(kk));
	    conjFreqs[kk]=SynthesisUtils::conjFreq(beamFreqs[kk], centerFreq);
	  }
	  conjSpCoord=SpectralCoordinate(spCoord.frequencySystem(), conjFreqs, spCoord.restFrequency());
	  //conjSpCoord.setIncrement(Vector<Double>(1, beamFreqs(0)-beamFreqs(1)));
	}
	conjCoord.replaceCoordinate(conjSpCoord, spind);

	
	TempLattice<Complex> vpFuncTemp(TiledShape(IPosition(5, convSize_p, convSize_p, nBeamPols, nBeamChans, ndish), IPosition(5, convSize_p, convSize_p, 1, 1, 1)), 0);
        TempLattice<Complex> pbFuncTemp(TiledShape(IPosition(5, convSize_p, convSize_p, nBeamPols, nBeamChans, ndish), IPosition(5, convSize_p, convSize_p, 1, 1, 1)), 0);
	IPosition begin(5, 0, 0, 0, 0, 0);
        IPosition end(5, vpFuncTemp.shape()[0]-1,  vpFuncTemp.shape()[1]-1, nBeamPols-1, nBeamChans-1, 0);
	
	IPosition pbShape(4, convSize_p, convSize_p, 1, nBeamChans);
	Long memtot=HostInfo::memoryFree();
        Double memtobeused= Double(memtot)*1024.0;
        if(memtot <= 2000000)
	  memtobeused=0.0;
        TempImage<Complex> vpScreen(TiledShape(pbShape), coords, memtobeused/2.2);
	
        TempImage<Complex> pbScreen(TiledShape(pbShape), ((nchan_p==1) && getConjConvFunc) ?conjCoord : coords , memtobeused/2.2);

	for (uInt k=0; k < ndish; ++k) {
	  antMath_p[k]->setBandOrFeedName(bandName_p);
	  IPosition blcin(4, 0, 0, 0, 0);
	  IPosition trcin(4, convSize_p-1, convSize_p-1, 0, 0);
	  for (Int kk=0; kk < nBeamChans; ++kk) {
	    blcin[3]=kk;
	    trcin[3]=kk;
      		    //wtime0=omp_get_wtime();
	    Slicer slin(blcin, trcin, Slicer::endIsLast);
	    SubImage<Complex> subim(vpScreen, slin, true);
	    subim.set(Complex(1.0, 0.0));
	    (antMath_p[k])->applyVP(subim, subim, fieldDir);
	    SubImage<Complex> subim2(pbScreen, slin, true);
	    if((nchan_p==1) && getConjConvFunc){
	      subim2.copyData(subim);
	      //multiplySelfConjugate(subim2, beamFreqs[kk]);
	      (antMath_p[k])->applyVP(subim2, subim2, fieldDir);
	    }
	    else{
	      subim2.set(Complex(1.0,0.0));
	      (antMath_p[k])->applyPB(subim2, subim2, fieldDir);
	    }
	    ft_p.c2cFFTInDouble(subim);
	    ft_p.c2cFFTInDouble(subim2);
	  }
	  begin[4]=k;
	  end[4]=k;
	  Slicer slplane(begin, end, Slicer::endIsLast);
 
	  IPosition blcQ(4, 0, 0, 0, 0);
	  IPosition trcQ(4,  blcQ[0]+ pbShape(0)-1, blcQ[1]+pbShape(1)-1 , nBeamPols-1, nBeamChans-1);

                //cerr << "blcQ " << blcQ << " trcQ " << trcQ << " pbShape " << pbShape << endl;
                Slicer slQ(blcQ, trcQ, Slicer::endIsLast);
                {
                    SubImage<Complex>  vpSSub(vpScreen, slQ, false);
                    SubLattice<Complex> vpTempSub(vpFuncTemp,  slplane, true);
                    LatticeConcat<Complex> lc(4);
                    lc.setLattice(vpSSub);
                    vpTempSub.copyData(lc);
                }
                {
                    SubImage<Complex>  pbSSub(pbScreen, slQ, false);
                    SubLattice<Complex> pbTempSub(pbFuncTemp,  slplane, true);
                    LatticeConcat<Complex> lc(4);
                    lc.setLattice(pbSSub);
                    pbTempSub.copyData(lc);  
                }
	  
	}
	///Do the support size
	//Int maxSupport=findMaxsupport(pbFuncTemp);
	Int maxSupport=100;
	///
	Int lattSize=vpFuncTemp.shape()(0);
        //(*convSupportBlock_p[actualConvIndex_p])=maxSupport;
	//cerr << "convsupport " << convSupport_p << endl;

        if(maxSupport >lattSize) {
	  throw(AipsError("Programmer Error support is larger than FTbeam"));
	}
	IPosition blc(5, (lattSize/2)-(maxSupport/2),
		      (lattSize/2)-(maxSupport/2),0,0,0);
	IPosition trc(5, (lattSize/2)+(maxSupport/2-1),
		      (lattSize/2)+(maxSupport/2-1), nBeamPols-1, nBeamChans-1,ndish-1);
	IPosition shp(5, maxSupport, maxSupport, nBeamPols, nBeamChans, ndish);

	  vpFunctions_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, maxSupport, maxSupport, nBeamPols, nBeamChans, ndish ));
	  pbFunctions_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, maxSupport, maxSupport, nBeamPols, nBeamChans, ndish));
	  ///temporary
	  (*vpFunctions_p[actualConvIndex_p])=vpFuncTemp.getSlice(blc,shp);
	  (*pbFunctions_p[actualConvIndex_p])=pbFuncTemp.getSlice(blc, shp);
	  ////resample later
	  //	  (*vpFunctions_p[actualConvIndex_p])=resample(vpFuncTemp.getSlice(blc,shp),Double(convSamp)/Double(convSampling));
	  
	  //  (*pbFunctions_p[actualConvIndex_p])=resample(pbFuncTemp.getSlice(blc, shp),Double(convSamp)/Double(convSampling));
	  doneMainConv_p[actualConvIndex_p]=True;  
    }
    //timer2_p+=omp_get_wtime()-wtime0;
    //wtime0=omp_get_wtime();
    //Double wtime1=omp_get_wtime();
    // cerr << std::setprecision(10) << "init time "<< wtime1-wtime0 << endl;  
    diffPointingCorrection(vb, origConvSize_p[actualConvIndex_p], convFuncRowMap, convSamp, beamFreqs, (nchan_p == 1) && getConjConvFunc);
     if((nchan_p == 1) && getConjConvFunc) {
       convSupport_p.resize();
       convSupport_p=*convSupportBlock_p[actualConvIndex_p];
       Int conjsupp=conjSupport(beamFreqs) ;
       if(conjsupp > max(convSupport_p))
	 conjsupp=Int(1.5*Float(conjsupp));
       //cerr << "conjsupp " << conjsupp << " bef " << max(*convSupportBlock_p[actualConvIndex_p]) << " convsize " <<max(*convSizes_p[actualConvIndex_p]) << endl;
       ///Need to keep those that are zeros
       (*convSupportBlock_p[actualConvIndex_p])((*convSupportBlock_p[actualConvIndex_p]) != 0)=max(conjsupp, max( *convSupportBlock_p[actualConvIndex_p]));   
      
     }
    

     
     // cerr << "diffcorr time " << omp_get_wtime()-wtime1 << endl;
     convFunc.reference(*convFunctions_p[actualConvIndex_p]);
     weightConvFunc.reference(*convWeights_p[actualConvIndex_p]);
     convsize=(*convSizes_p[actualConvIndex_p]);
     convSupport= *convSupportBlock_p[actualConvIndex_p];
     //timer3_p+=omp_get_wtime()-wtime0;
     //cerr << "timers " << timer1_p << "    " << timer2_p <<"  " << timer3_p << " supp and size "<< max(convSupport) << "  " << max(convsize) << " field ID "<< vb.fieldId()(0) << " spw " << vb.spectralWindows()(0)<< endl;
  
 ////////////////TESTOOO
     //CoordinateSystem TMP = coords;
     //CoordinateUtil::addLinearAxes(TMP, Vector<String>(1,"gulu"), IPosition(1,convFunc.shape()(4))); 
     //PagedImage<Complex> SCREEN(TiledShape(convFunctions_p[actualConvIndex_p]->shape()), TMP, "convFunc"+String::toString(actualConvIndex_p));
     //SCREEN.put(*convFunctions_p[actualConvIndex_p]  );
    //	   PagedImage<Complex> SCREEN3(TiledShape(convWeights_p[actualConvIndex_p]->shape()), TMP, "FTWEIGHTVI2"+String::toString(actualConvIndex_p));
    //	  SCREEN3.put(*convWeights_p[actualConvIndex_p]  );
     
  
  }

  void HetArrayConvFunc::diffPointingCorrection(const vi::VisBuffer2& vb, const Int origSupportSize, Vector<Int>& convFuncRowMap, const Int convSamp, const Vector<Double>& freqs, const Bool doConj){
    if(vbutil_p.null())
      vbutil_p=new VisBufferUtil(vb);
   
    DirectionCoordinate dc=dc_p;
    IPosition shp=vpFunctions_p[actualConvIndex_p]->shape();
    //For now everyrow will get its own convfunction and weightfunction
    shp(4)=vb.nRows();
   
    //convFunctions_p[actualConvIndex_p]=new Array<Complex>(shp);
    //convWeights_p[actualConvIndex_p]=new Array<Complex>(shp);
    
    dc.setWorldAxisUnits(Vector<String>(2, "rad"));
    IPosition begin(5, 0, 0, 0, 0, 0);
    IPosition end(5, shp[0]-1, shp[1]-1, shp[2]-1, shp[3]-1, 0);
    ///not doing anything to economize on nrow right now
    convSupportBlock_p[actualConvIndex_p]=new Vector<Int>(shp[4]);
    //MDirection prevDir;
    Vector<Double>prevDirPix(2,-1);
    Int prevIndex=-1;
    IPosition prevBeg, prevEnd;
    MDirection dir1=vbutil_p->getPhaseCenter(vb);
    MDirection  dir2=dir1;
    Double wtime0;
    for (Int k=0; k< vb.nRows() ; ++k){
     
      begin[4]=k;
      end[4]=k;
      Int index1=antIndexToDiamIndex_p(vb.antenna1()(k));
      Int index2=antIndexToDiamIndex_p(vb.antenna2()(k));
      Double pixshift_x=0.0;
      Double pixshift_y=0.0;
      if(usePointingTable_p){
	wtime0=omp_get_wtime();
	dir1=vbutil_p->getPointingDir(vb, vb.antenna1()(k), k, dc_p.directionType());
	dir2=vbutil_p->getPointingDir(vb, vb.antenna2()(k), k, dc_p.directionType());
	timer1_p+=omp_get_wtime()-wtime0;
	wtime0=omp_get_wtime();
	roundShiftInPointing(pixshift_x, pixshift_y, origSupportSize, dir1, dir2, dc);
	//cerr << MDirection::showType(dc_p.directionType()) << " dirs " << dir1.toString() <<  " -- " <<  dir1.toString() << endl;
	timer2_p+=omp_get_wtime()-wtime0;
	
      }
      wtime0=omp_get_wtime();
      //cerr << "first bit " << omp_get_wtime()-wtime0 << endl;
      
      
      if(Int(abs(pixshift_x)) > origSupportSize/2 || Int(abs(pixshift_y)) > origSupportSize/2){

	//cerr << "antIDs (" << vb.antenna1()(k) << ", " << vb.antenna2()(k) << ") pixshifts "<< pixshift_x << "  ,  " << pixshift_y << " DIRS " << dir1.toString() << "  --->   "  << dir2.toString() << " fieldID " << vb.fieldId()(k)<< endl;
	///TESTOO
	/*{
	  Vector<Double> pxx(2);
	  dc_p.toPixel(pxx, dir1);
	  cerr << "dir1 pixel " << pxx;
	  dc_p.toPixel(pxx, dir2);
	  cerr << " dir2  " << pxx << endl;
	  }*/
	//////
	(*convSupportBlock_p[actualConvIndex_p])[k]=0;
	////Skip the remainder
	continue;
      }
      std::tuple<Int, Int,Int,Int,Int> bl(index1, index2, actualConvIndex_p, Int(pixshift_x), Int(pixshift_y));
      Int index=-1;
      //cerr << "bl index " << std::get<0>(bl) <<  ","<< std::get<1>(bl) <<  ","<< std::get<2>(bl) << "," <<  std::get<3>(bl) << "," << std::get<4>(bl)   << endl;
      if(baselinePointingDiffIndex_p.find(bl) == baselinePointingDiffIndex_p.end()){
	///fill diff pointing pb and pb2 for this pair of beams
	
	index=baselinePointingDiffIndex_p.size();
	baselinePointingDiffIndex_p[bl]=index;
	baslPBFunctions_p.resize(index+1, False, True);
	baslWeightFunctions_p.resize(index+1, False, True);
	baslSupport_p.resize(index+1, False, True);
	baslPBFunctions_p[index]=new Array<Complex>();
	baslWeightFunctions_p[index]=new Array<Complex>();
	baslSupport_p[index]=new Vector<Int>();
	rephaseBeams(*baslPBFunctions_p[index], *baslWeightFunctions_p[index],vb, k, origSupportSize, pixshift_x, pixshift_y);
	timer3_p+=omp_get_wtime() - wtime0;
	//support
	///////////////////////////
	convSupport_p.resize(1);
	convSupport_p=10000;
	Double fac=doConj ?  1.6*Double(conjSupport(freqs))/10000.0 : 1.0;
	///////////////
	if(fac <1.0)
	  fac=1.0;		   
	if(!supportAndNormalize(*baslPBFunctions_p[index], *baslWeightFunctions_p[index], *baslSupport_p[index], fac)){
	  (*convSupportBlock_p[actualConvIndex_p])[k]=0;
	  //cerr << "Failed suppport bl index " << std::get<0>(bl) <<  ","<< std::get<1>(bl) <<  ","<< std::get<2>(bl) << "," <<  std::get<3>(bl) << "," << std::get<4>(bl)   << endl;
	  //cerr << index << " maxes " << max(*baslPBFunctions_p[index]) << "      "  << max(*baslWeightFunctions_p[index]) << endl;
	  continue;
	}
	//	cerr << index << " afternorm " << sum(*baslPBFunctions_p[index]) << "      "  << max(*baslWeightFunctions_p[index]) << " support " << (*baslSupport_p[index]) << endl;
	baslPBFunctions_p[index]->assign(resample(*(baslPBFunctions_p[index]), Double(convSamp)));
	baslWeightFunctions_p[index]->assign(resample(*(baslWeightFunctions_p[index]), Double(convSamp)));
	

      }
      else{
	index=baselinePointingDiffIndex_p[bl];
	if(max(*baslSupport_p[index])==0){
	   (*convSupportBlock_p[actualConvIndex_p])[k]=0;
	   continue;
	}
      }
     
      
      Vector<Double>dirpix(2);
      dc_p.toPixel(dirpix, dir1);
      //dc_p.toPixel(prevDirPix, prevDir);
      dirpix[0]=dirpix[0]-nx_p/2;
      dirpix[1]=dirpix[1]-ny_p/2;
      //cerr << "index " << index << "actualConv " << actualConvIndex_p << " pix shift " << dirpix << endl;
      IPosition blc;
      IPosition trc;
      if(k < 1 || !(allEQ(prevDirPix, dirpix) && prevIndex==index)){
	
	Array<Complex> tmpArr;
	Array<Complex> tmpArr2;
	if(k >=1 && ((baslPBFunctions_p[index]->shape()[0]) !=  shp[0])){
	  if((baslPBFunctions_p[index]->shape()[0])> shp[0]){
	    /////This means that the first baseline VPs were wider apart than the latter ones.
	    ///Let's redo the size of the beams and beam^2s 
	    Int diff=(baslPBFunctions_p[index]->shape()[0])-shp[0];
	    blc=IPosition(5,  diff/2, diff/2, 0,0,0);
	    trc=IPosition(5, (baslPBFunctions_p[index]->shape()[0])-diff/2-1, (baslPBFunctions_p[index]->shape()[1])-diff/2-1, shp[2]-1, shp[3]-1, shp[4]-1);
	    ////Let us resize and copy back the smaller arrays
	    tmpArr.assign(*convFunctions_p[actualConvIndex_p]);
	    shp[0]=baslPBFunctions_p[index]->shape()[0];
	    shp[1]=baslPBFunctions_p[index]->shape()[1];
	    //cerr << "resizing the whole enchilada " << endl;
	    convFunctions_p[actualConvIndex_p]->resize(shp, False, ArrayInitPolicies::INIT);
	    //cerr << "SHAPES0 " <<  (*convFunctions_p[actualConvIndex_p])(blc, trc).shape() << "   " << tmpArr.shape() << endl;
	    (*convFunctions_p[actualConvIndex_p])(blc, trc)=tmpArr;
	    tmpArr.assign(*convWeights_p[actualConvIndex_p]);
	    convWeights_p[actualConvIndex_p]->resize(shp, False, ArrayInitPolicies::INIT);
	    (*convWeights_p[actualConvIndex_p])(blc, trc)=tmpArr;
	    ///readjust end
	    end[0]=shp[0]-1;
	    end[1]=shp[1]-1;
	    ////
	    tmpArr.assign((*baslPBFunctions_p[index]));
	    ///did the abs
	    tmpArr2.assign(*baslWeightFunctions_p[index]);
	    ////If the difference is so large as compared from 1st baseline.
	    ///should flag the first baseline i suspect
	    (*convSupportBlock_p[actualConvIndex_p])[k]=max(*baslSupport_p[index]);
	    
	    //cerr << "reshaped " << tmpArr.shape() << endl;
	  }
	  else{
	    trc=IPosition(5, shp[0], shp[1], shp[2], shp[3],1);
	    tmpArr.resize(trc);
	    tmpArr2.resize(trc);
	    tmpArr.set(0.0);
	    tmpArr2.set(0.0);
	    Int diff=shp[0]-(baslPBFunctions_p[index]->shape()[0]);
	    blc=IPosition(5,  diff/2, diff/2, 0,0,0);
	    ///Do the inner part
	    trc=IPosition(5, shp[0]-diff/2-1, shp[1]-diff/2-1, shp[2]-1, shp[3]-1, 0);
	    //cerr << tmpArr.shape() << " SHAPES " << tmpArr(blc,trc).shape() << "  ==   " << baslPBFunctions_p[index]->shape() <<  " blc trc "<< blc << "    " << trc  << endl;
	    tmpArr(blc,trc)=*baslPBFunctions_p[index];
	    tmpArr2(blc,trc)=*baslWeightFunctions_p[index];
	    //cerr << "reshaped2 " << tmpArr.shape() << endl;
	    (*convSupportBlock_p[actualConvIndex_p])[k]=max(*baslSupport_p[index]);


	    
	  }
	}
	else{
	  tmpArr.assign(*baslPBFunctions_p[index]);
	  tmpArr2.assign(*baslWeightFunctions_p[index]);
	  (*convSupportBlock_p[actualConvIndex_p])[k]=max(*baslSupport_p[index]);
	    //blc=IPosition(5,0,0,0,0,0);
	    //trc=tmpArr.shape()-1;
	}
	  //////TESTOO
	  /*if((*baslSupport_p[index]) > shp[0]/2){
	    cerr <<"TROUBLE " << "bl index " << std::get<0>(bl) <<  ","<< std::get<1>(bl) <<  ","<< std::get<2>(bl) << "," <<  std::get<3>(bl) << "," << std::get<4>(bl)   << endl;
	    cerr << "SHAPES " << shp <<  "    " << tmpArr.shape() << "   " << baslPBFunctions_p[index]->shape() << " blc trc " << blc << "     " << trc << endl;

	    }*/
	  /////////////
	//cerr << "k=" << k << " index " << index << " shapes " << tmpArr.shape() << "     " << shp << endl;
	//cerr << "bl index " << std::get<0>(bl) <<  ","<< std::get<1>(bl) <<  ","<< std::get<2>(bl) << "," <<  std::get<3>(bl) << "," << std::get<4>(bl)   << endl;
	if(k==0 ||  convFunctions_p[actualConvIndex_p].null()){
	  
	  shp[0]=tmpArr.shape()[0];
	  shp[1]=tmpArr.shape()[1];
	  if(convFunctions_p[actualConvIndex_p].null()){
	    ///Interestingly this malloc here takes a lot of time
	    // try to avoid when it can be.
	    convFunctions_p[actualConvIndex_p]=new Array<Complex>(shp);
	    convWeights_p[actualConvIndex_p]=new Array<Complex>(shp);
	  }
	  else if(convFunctions_p[actualConvIndex_p]->shape() != shp){
	    convFunctions_p[actualConvIndex_p]->resize(shp);
	    convWeights_p[actualConvIndex_p]->resize(shp);
	  }
	  end[0]=shp[0]-1;
	  end[1]=shp[1]-1;
	}
	//cerr <<k << " shape " << tmpArr.shape() << endl;
	//if((convFunctions_p[actualConvIndex_p])->shape()[0] != tmpArr.shape()[0]){
	 
	if(tmpArr.shape()[0] != shp[0]){
	  //cerr << "k=" << k << " shapes " << tmpArr.shape() << "     " << shp << endl;
	  
	  throw(AipsError("Tell the programmer something unexpected happened"));

	}
	  //shp[0]=tmpArr.shape()[0];
	  //shp[1]=tmpArr.shape()[1];
	 
	  //(convFunctions_p[actualConvIndex_p])->resize(shp);
	  //(convWeights_p[actualConvIndex_p])->resize(shp);
	  //convFunctions_p[actualConvIndex_p]=new Array<Complex>(shp);
	  //convWeights_p[actualConvIndex_p]=new Array<Complex>(shp);
	
	  //}
	  //end[0]=shp[0]-1;
	  //end[1]=shp[1]-1;

	////Have to do this properly per plane ...just testing now
	if(doConj){
	  //Array<Complex> backup;
	  //backup=*(convFunctions_p[actualConvIndex_p]);
	  IPosition blcfreq(tmpArr.shape().nelements(),0);
	  IPosition trcfreq=tmpArr.shape() - 1;
	  //cerr << "conj shapes " << tmpArr.shape() << endl;
	  for (uInt f=0; f < freqs.nelements(); ++f){
	    blcfreq[3]=f;
	    trcfreq[3]=f;
	    Array<Complex> arrSlice=tmpArr(blcfreq, trcfreq);
	    Array<Complex> outSlice(arrSlice.shape());
	    fillConjConvFunc(outSlice, arrSlice, freqs[f]);
	    arrSlice=outSlice;
	  }
	}
	
	//cerr << "third bit " << omp_get_wtime()-wtime0 << endl;
	applyPhaseGradient(tmpArr, dirpix[0], dirpix[1], nx_p, ny_p, convSamp);
	
	//(*(convFunctions_p[actualConvIndex_p]))(begin,end)=resample(tmpArr, Double(convSamp));
	(*(convFunctions_p[actualConvIndex_p]))(begin,end)=tmpArr;
	//tmpArr.assign((*baslWeightFunctions_p[index])(blc,trc));
	applyPhaseGradient( tmpArr2, dirpix[0], dirpix[1], nx_p, ny_p, convSamp);	
	//(*(convWeights_p[actualConvIndex_p]))(begin,end)=resample(tmpArr, Double(convSamp));
	(*(convWeights_p[actualConvIndex_p]))(begin,end)=tmpArr2;
      //Use only the max for this baseline for now.
	//DONE above for different conditions
	//(*convSupportBlock_p[actualConvIndex_p])[k]=max(*baslSupport_p[index]);


      }
      else{
	////previous values are the same
	(*(convFunctions_p[actualConvIndex_p]))(begin,end)=(*(convFunctions_p[actualConvIndex_p]))(prevBeg,prevEnd);
	(*(convWeights_p[actualConvIndex_p]))(begin,end)=	(*(convWeights_p[actualConvIndex_p]))(prevBeg,prevEnd);
	(*convSupportBlock_p[actualConvIndex_p])[k]=	(*convSupportBlock_p[actualConvIndex_p])[k-1]	;
      }
      prevDirPix=dirpix;
      //cerr << "prevDirpix " << prevDirPix << "  " << dirpix << endl;
      prevBeg=begin;
      prevEnd=end;
      prevIndex=index;
    }
    ///make row map
    convFuncRowMap.resize(vb.nRows());
    indgen(convFuncRowMap);
    ///For now letting convSize the same as original it came in
    convSizes_p[actualConvIndex_p]=new Vector<Int>(vb.nRows());
    convSizes_p[actualConvIndex_p]->set(shp[0]);
  }

  void HetArrayConvFunc::roundShiftInPointing(Double& xshift, Double& yshift, const Int origSupport, const MDirection& dir1, const MDirection& dir2, const DirectionCoordinate dc){

    Vector<Double>diff(2);
    Vector<Double>dirpix(2);
    dc.toPixel(dirpix, dir1);
    dc.toPixel(diff, dir2);
    diff-=dirpix;
    xshift=diff(0);
    yshift=diff(1);
    ////TESTOOO
    //if(abs(xshift)> 500){
    //  cerr << "DIRS " << dir2.toString() << "  -- " << dir1.toString() << endl;  
    //  cerr << "diff " <<  dir2.getAngle("rad").getValue("rad")-dir1.getAngle("rad").getValue("rad") << " diff0 " << diff << endl;
      
    //}
    /////////////
   
    Double ratio=std::round(Double(origSupport)/Double(vpFunctions_p[actualConvIndex_p]->shape()(0)));
    xshift=std::round(xshift/ratio)*ratio;
    yshift=std::round(yshift/ratio)*ratio;

  }
  
  void HetArrayConvFunc::rephaseBeams(Array<Complex>& basPB,
				      Array<Complex>& basPB2, const vi::VisBuffer2& vb,
				      const Int row, const Int origSize, const Double pixshift_x,
				      const Double pixshift_y){
       Int ind1=antIndexToDiamIndex_p(vb.antenna1()(row));
       Int ind2=antIndexToDiamIndex_p(vb.antenna2()(row));
       Int size=vpFunctions_p[actualConvIndex_p]->shape()(0);
       Int npol=vpFunctions_p[actualConvIndex_p]->shape()(2);
       Int nchan=vpFunctions_p[actualConvIndex_p]->shape()(3);
       Array<Complex> vp2;
       Array<Complex> pb2;
       IPosition begin(5, 0, 0, 0, 0, ind2);
       IPosition end(5, size-1, size-1, npol-1, nchan-1, ind2);
       vp2.assign((*(vpFunctions_p[actualConvIndex_p]))(begin, end));
       pb2.assign((*(pbFunctions_p[actualConvIndex_p]))(begin, end));
       //cerr << "maxes vp pb2 " << max(vp2) << "   " << max(pb2) << endl;
       applyPhaseGradient(vp2, pixshift_x, pixshift_y, origSize, origSize);
       applyPhaseGradient(pb2, pixshift_x, pixshift_y,origSize, origSize);
       //cerr << "after phgrad " << max(vp2) << "   " << max(pb2) << endl;
       Array<Complex> vp1;
       Array<Complex> pb1;
       begin[4]=ind1;
       end[4]=ind1;
       vp1.assign((*(vpFunctions_p[actualConvIndex_p]))(begin, end));
       pb1.assign((*(pbFunctions_p[actualConvIndex_p]))(begin, end));
       //cerr << "maxes vp pb2 " << max(vp1) << "   " << max(pb1) << endl;
       Bool st1, st2, stb1, stb2;
       Complex *pta1=vp1.getStorage(st1);
       Complex *pta2=vp2.getStorage(st2);
       Complex *ptb1=pb1.getStorage(stb1);
       Complex *ptb2=pb2.getStorage(stb2);
       FFT2D ft;
       ///FFT for all planes
       for (Int k=0; k < npol*nchan; ++k){
	 Complex* pt;
	 pt=pta1+(k*size*size);
	 ft.c2cFFT(pt, size, size, False);
	 pt=pta2+(k*size*size);
	 ft.c2cFFT(pt, size, size, False);
	 pt=ptb1+(k*size*size);
	 ft.c2cFFT(pt, size, size, False);
	 pt=ptb2+(k*size*size);
	 ft.c2cFFT(pt, size, size, False);
       }
       vp1.putStorage(pta1, st1);
       vp2.putStorage(pta2, st2);
       //cerr << "after fft " << max(vp1) << "   " << max(vp2) << endl;
       vp1=vp1*vp2;
       pb1.putStorage(ptb1, stb1);
       pb2.putStorage(ptb2, stb2);
       pb1=pb1*pb2;
       pta1=vp1.getStorage(st1);
       ptb1=pb1.getStorage(stb1);
       for (Int k=0; k < npol*nchan; ++k){
	 Complex *pt;
	 pt=pta1+(k*size*size);
	 ft.c2cFFT(pt, size, size, True);
	 pt=ptb1+(k*size*size);
	 ft.c2cFFT(pt, size, size, True);
       }
       vp1.putStorage(pta1, st1);
       pb1.putStorage(ptb1, stb1);
       //cerr << "after conv " << max(vp1) << "   " << max(pb1) << endl;
       //rescale to number of pixels
	Double sizecorr=Double(size)*Double(size)/Double(origSize)/Double(origSize);
       vp1= vp1/(sizecorr);
       pb1= pb1/(sizecorr);
       //cerr << "after size corr " << max(vp1) << "   " << max(pb1) << endl;
       basPB.reference(vp1);
       basPB2.reference(pb1);
       //cerr << "after reff " << max(basPB) << "   " << max(basPB2) << endl;
    }
  void HetArrayConvFunc::applyPhaseGradient(Array<Complex>& arr, const Double& xshift, const Double& yshift, const Int xsize, const Int ysize, const Int convSamp){
  
    
    IPosition shp=arr.shape();
    Int nx=shp(0);
    Int ny=shp(1);
    Double pshift_x=xshift*(-2.0*C::pi/Double(xsize*convSamp)) ;
    Double pshift_y=yshift*(-2.0*C::pi/Double(ysize*convSamp));
    ooLong nplanes=ooLong(shp(4))*ooLong(shp(3))*ooLong(shp(2));
    //cerr << nx << "   " << ny << " nplanes " << nplanes << " phaseshift "<< pshift_x << "    " << pshift_y << endl;
    Bool st;
    Complex* pta=arr.getStorage(st);
#ifdef _OPENMP
	omp_set_nested(0);
#endif
#pragma omp parallel for default(none) firstprivate(  pta, pshift_x,  pshift_y, nx, ny, nplanes)
    for(Int j=0; j < ny; ++j){
      Double cy, sy;
      SINCOS(Double(j-ny/2)*pshift_y, sy, cy);
      Complex phy(cy, sy);
      for(Int i=0; i< nx; ++i){
	Double cx, sx;
	SINCOS(Double(i-nx/2)*pshift_x, sx, cx);
	Complex phx(cx, sx);
	phx=phx*phy;
	for (ooLong k=0; k < nplanes; ++k){
	  ooLong indx=k*ooLong(nx)*ooLong(ny)+ooLong(j)*ooLong(nx)+ooLong(i);
	  pta[indx]  = pta[indx] *phx;
	}
      }      
    }
    arr.putStorage(pta, st);
    /*
    IPosition start(5,0,0,0,0,0);
    IPosition end(5,nx-1, ny-1, 0,0,0);
    for (Int k=0; k < shp[4]; ++k)
      for (Int l=0; l < shp[3];++l)
	for (Int m=0; m<shp[2]; ++m){
	  start[4]=k;
	  start[3]=l;
	  start[2]=m;
	  end[4]=k;
	  end[3]=l;
	  end[2]=m;
	  Matrix<Complex> pl(arr(start, end).reform(IPosition(2,nx,ny)));
	  for (Int j=0; j < ny; ++j){
	    Double cy, sy;
	    SINCOS(Double(j-ny/2)*pshift_y, sy, cy);
	    Complex phy(cy, sy);
	    for(Int i=0; i< nx; ++i){
	      Double cx, sx;
	      SINCOS(Double(i-nx/2)*pshift_x, sx, cx);
	      Complex phx(cx, sx);
	      pl(i,j)=pl(i,j)*phx*phy;
	    }
	  }

	}
    */
    


    
  }
  
void HetArrayConvFunc::findConvFunction(const ImageInterface<Complex>& iimage,
                                        const vi::VisBuffer2& vb,
                                        const Int& convSamp, const Vector<Double>& visFreq,
                                        Array<Complex>& convFunc,
                                        Array<Complex>& weightConvFunc,
                                        Vector<Int>& convsize,
                                        Vector<Int>& convSupport,
                                        Vector<Int>& convFuncPolMap,
                                        Vector<Int>& convFuncChanMap,
                                        Vector<Int>& convFuncRowMap, Bool getConjConvFunc,
					const MVDirection& extraShift, const Bool useExtraShift)
{

  if(usePointingTable_p){
    findConvFunction2(iimage, vb, convSamp, visFreq, convFunc, weightConvFunc, convsize, convSupport, convFuncPolMap, convFuncChanMap, convFuncRowMap, getConjConvFunc, extraShift, useExtraShift);
  return;
  }
  //////////////////////////////////////////////////////////////////
    storeImageParams(iimage,vb);
    convFuncChanMap.resize(vb.nChannels());
    Vector<Double> beamFreqs;
    findUsefulChannels(convFuncChanMap, beamFreqs, vb, visFreq);
    //cerr << "SPW " << vb.spectralWindow() << "   beamFreqs "<< beamFreqs <<  " chamMap " << convFuncChanMap << endl;
    Int nBeamChans=beamFreqs.nelements();
    /////For now not doing beam rotation or squints but to be enabled easily
    convFuncPolMap.resize(vb.nCorrelations());
    Int nBeamPols=1;
    convFuncPolMap.set(0);
    findAntennaSizes(vb);
    uInt ndish=antMath_p.nelements();
    if(ndish==0)
        throw(AipsError("Don't have dishsize"));
    Int ndishpair;
    if(ndish==1)
        ndishpair=1;
    else
        ndishpair=factorial(ndish)/factorial(ndish-2)/2 + ndish;

    convFunc.resize();
    weightConvFunc.resize();
    convFuncRowMap.resize();
    convsize.resize();
    convSupport.resize();

    Int isCached=checkPBOfField(vb, convFuncRowMap, extraShift, useExtraShift);
    //cout << "isCached " << isCached <<  endl;
    if(isCached==1 && (convFuncRowMap.shape()[0]==vb.nRows())) {
        /*convFunc.reference(convFunc_p);
        weightConvFunc.reference(weightConvFunc_p);
        convsize=*convSizes_p[actualConvIndex_p];
        convSupport=convSupport_p;
         return;
        */
    }
    else if(isCached ==2) {
        convFunc.resize();
        weightConvFunc.resize();
        convsize.resize();
        convSupport.resize();
        convFuncRowMap.resize();
        return;

    }
    actualConvIndex_p=convIndex(vb);
    //cerr << "actual conv index " << actualConvIndex_p << " doneMainconv " << doneMainConv_p << endl;
    if(doneMainConv_p.shape()[0] < (actualConvIndex_p+1)) {
        //cerr << "resizing DONEMAIN " <<   doneMainConv_p.shape()[0] << endl;
        doneMainConv_p.resize(actualConvIndex_p+1, true);
        doneMainConv_p[actualConvIndex_p]=false;
        convFunctions_p.resize(actualConvIndex_p+1);
        convFunctions_p[actualConvIndex_p]=nullptr;
    }
    ///// In multi ms mode ndishpair may change when meeting a new ms
    //// redo the calculation then
    if(msId_p != vb.msId())//doneMainConv_p[actualConvIndex_p] && ((convSupport_p.nelements() != uInt(ndishpair)) || convFunctions_p[actualConvIndex_p]->shape()[3] != nBeamChans))
    {
        doneMainConv_p[actualConvIndex_p]=false;
        //cerr << "invalidating doneMainConv " <<  convFunctions_p[actualConvIndex_p]->shape()[3] << " =? " << nBeamChans << " convsupp " << convSupport_p.nelements() << endl;
    }

    // Get the coordinate system
    CoordinateSystem coords(iimage.coordinates());
    Int directionIndex=coords.findCoordinate(Coordinate::DIRECTION);
    AlwaysAssert(directionIndex>=0, AipsError);
    // Set up the convolution function.
    Int nx=nx_p;
    Int ny=ny_p;
    Int support=max(nx_p, ny_p)/10;
    Int convSampling=1;
    if(!doneMainConv_p[actualConvIndex_p]) {
        for (uInt ii=0; ii < ndish; ++ii) {
            support=max((antMath_p[ii])->support(coords), support);
        }
	
        support=Int(min(Float(support), max(Float(nx_p), Float(ny_p)))*2.0)/2;
        convSize_p=support*convSampling;
        // Make this a nice composite number, to speed up FFTs
        CompositeNumber cn(uInt(convSize_p*2.0));
        convSize_p  = cn.nextLargerEven(Int(convSize_p));
        convSize_p=(convSize_p/16)*16;  // need it to be divisible by 8 in places
	
    }


    DirectionCoordinate dc=dc_p;
    //where in the image in pixels is this pointing
    Vector<Double> pixFieldDir(2);
    if(doneMainConv_p.shape()[0] < (actualConvIndex_p+1)) {
        //cerr << "resizing DONEMAIN " <<   doneMainConv_p.shape()[0] << endl;
        doneMainConv_p.resize(actualConvIndex_p+1, true);
        doneMainConv_p[actualConvIndex_p]=false;
    }
    //no need to call toPix here as its been done already above in checkOFPB
    //thus the values are still current.
    pixFieldDir=thePix_p;
    //toPix(pixFieldDir, vb);
    MDirection fieldDir=direction1_p;
    //shift from center
    pixFieldDir(0)=pixFieldDir(0)- Double(nx / 2);
    pixFieldDir(1)=pixFieldDir(1)- Double(ny / 2);
    //phase gradient per pixel to apply
    pixFieldDir(0)=-pixFieldDir(0)*2.0*C::pi/Double(nx)/Double(convSamp);
    pixFieldDir(1)=-pixFieldDir(1)*2.0*C::pi/Double(ny)/Double(convSamp);


    if(!doneMainConv_p[actualConvIndex_p]) {
      //cerr << "doneMainConv_p " << actualConvIndex_p << endl;

        Vector<Double> sampling;
	
        sampling = dc.increment();
	sampling*=Double(convSampling);
	sampling(0)*=Double(nx)/Double(convSize_p);
	sampling(1)*=Double(ny)/Double(convSize_p);
        dc.setIncrement(sampling);

        Vector<Double> unitVec(2);
        unitVec=convSize_p/2;
        dc.setReferencePixel(unitVec);
        //make sure we are using the same units
        fieldDir.set(dc.worldAxisUnits()(0));
        dc.setReferenceValue(fieldDir.getAngle().getValue());
        coords.replaceCoordinate(dc, directionIndex);
	Int spind=coords.findCoordinate(Coordinate::SPECTRAL);
        SpectralCoordinate spCoord=coords.spectralCoordinate(spind);
        spCoord.setReferencePixel(Vector<Double>(1,0.0));
        spCoord.setReferenceValue(Vector<Double>(1, beamFreqs(0)));
        if(beamFreqs.nelements() >1)
	  spCoord.setIncrement(Vector<Double>(1, beamFreqs(1)-beamFreqs(0)));
	//cerr << "spcoord " ;
	//spCoord.print(std::cerr);
        coords.replaceCoordinate(spCoord, spind);
	CoordinateSystem conjCoord=coords;
	Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
	SpectralCoordinate conjSpCoord=spCoord;
		//cerr << "centreFreq " << centerFreq << " beamFreqs " << beamFreqs(0) << "  " << beamFreqs(1) << endl;
	conjSpCoord.setReferenceValue(Vector<Double>(1,SynthesisUtils::conjFreq(beamFreqs[0], centerFreq)));
	///Increment should go in the reverse direction
	////Do a tabular spectral coordinate for more than 1 channel 
	if(beamFreqs.nelements() >1){
	  Vector<Double> conjFreqs(beamFreqs.nelements());
	  for (uInt kk=0; kk< beamFreqs.nelements(); ++kk){
	    //conjFreqs[kk]=sqrt(2*centerFreq*centerFreq-beamFreqs(kk)*beamFreqs(kk));
	    conjFreqs[kk]=SynthesisUtils::conjFreq(beamFreqs[kk], centerFreq);
	  }
	  conjSpCoord=SpectralCoordinate(spCoord.frequencySystem(), conjFreqs, spCoord.restFrequency());
	  //conjSpCoord.setIncrement(Vector<Double>(1, beamFreqs(0)-beamFreqs(1)));
	}
	conjCoord.replaceCoordinate(conjSpCoord, spind);
        IPosition pbShape(4, convSize_p, convSize_p, 1, nBeamChans);
        //TempImage<Complex> twoDPB(pbShape, coords);
	
	
        TempLattice<Complex> convFuncTemp(TiledShape(IPosition(5, convSize_p/4, convSize_p/4, nBeamPols, nBeamChans, ndishpair), IPosition(5, convSize_p/4, convSize_p/4, 1, 1, 1)), 0);
        TempLattice<Complex> weightConvFuncTemp(TiledShape(IPosition(5, convSize_p/4, convSize_p/4, nBeamPols, nBeamChans, ndishpair), IPosition(5, convSize_p/4, convSize_p/4, 1, 1, 1)), 0);
        //convFunc_p.resize(IPosition(5, convSize_p, convSize_p, nBeamPols, nBeamChans, ndishpair));

        // convFunc_p=0.0;
        //weightConvFunc_p.resize(IPosition(5, convSize_p, convSize_p, nBeamPols, nBeamChans, ndishpair));
        //weightConvFunc_p=0.0;
        IPosition begin(5, 0, 0, 0, 0, 0);
        IPosition end(5, convFuncTemp.shape()[0]-1,  convFuncTemp.shape()[1]-1, nBeamPols-1, nBeamChans-1, 0);
        //FFTServer<Float, Complex> fft(IPosition(2, convSize_p, convSize_p));
        //TempImage<Complex> pBScreen(TiledShape(pbShape, IPosition(4, convSize_p, convSize_p, 1, 1)), coords, 0);
        //TempImage<Complex> pB2Screen(TiledShape(pbShape, IPosition(4, convSize_p, convSize_p, 1, 1)), coords, 0);
        Long memtot=HostInfo::memoryFree();
        Double memtobeused= Double(memtot)*1024.0;
        if(memtot <= 2000000)
            memtobeused=0.0;
        TempImage<Complex> pBScreen(TiledShape(pbShape), coords, memtobeused/2.2);
		
        TempImage<Complex> pB2Screen(TiledShape(pbShape), ((nchan_p==1) && getConjConvFunc) ?conjCoord : coords  , memtobeused/2.2);
        IPosition start(4, 0, 0, 0, 0);
        convSupport_p.resize(ndishpair);
	//////////////////
	/*Double wtime0=0.0;
	Double wtime1=0.0;
	Double wtime2=0.0;
	wtime0=omp_get_wtime()
	*/;
	//////////////
        for (uInt k=0; k < ndish; ++k) {

            for (uInt j =k ; j < ndish; ++j) {
                //Timer tim;
                //Matrix<Complex> screen(convSize_p, convSize_p);
                //screen=1.0;
                //pBScreen.putSlice(screen, start);
                //cerr << "k " << k << " shape " << pBScreen.shape() <<  " direction1 " << direction1_p << " direction2 " << direction2_p << endl;

                //pBScreen.set(Complex(1.0, 0.0));
                //one antenna
		antMath_p[k]->setBandOrFeedName(bandName_p);
		antMath_p[j]->setBandOrFeedName(bandName_p);
                IPosition blcin(4, 0, 0, 0, 0);
                IPosition trcin(4, convSize_p-1, convSize_p-1, 0, 0);
                for (Int kk=0; kk < nBeamChans; ++kk) {
                    blcin[3]=kk;
                    trcin[3]=kk;
      		    //wtime0=omp_get_wtime();
                    Slicer slin(blcin, trcin, Slicer::endIsLast);
                    SubImage<Complex> subim(pBScreen, slin, true);
                    subim.set(Complex(1.0, 0.0));
                    (antMath_p[k])->applyVP(subim, subim, direction1_p);

                    //Then the other
                    (antMath_p[j])->applyVP(subim, subim, direction2_p);
                    //tim.show("After Apply ");
                    //tim.mark();
                    //pB2Screen.set(Complex(1.0,0.0));
                    SubImage<Complex> subim2(pB2Screen, slin, true);
					subim2.set(Complex(1.0,0.0));
                    
					if(nchan_p >1 || !getConjConvFunc){
						//one antenna
						(antMath_p[k])->applyPB(subim2, subim2, direction1_p);
						//Then the other
						(antMath_p[j])->applyPB(subim2, subim2, direction2_p);
					}
					else{
						//direct frequency PB
						//cerr << "orig coords " << subim.coordinates().toWorld(IPosition(4,0,0,0,0)) << " conj coords " <<  subim2.coordinates().toWorld(IPosition(4,0,0,0,0)) << endl;
						//cerr << "incr " << subim.coordinates().increment() << "   " << subim2.coordinates().increment() << endl;
						subim2.copyData(subim);
						//Now do the conjugate freq multiplication
						(antMath_p[k])->applyVP(subim2, subim2, direction1_p);

						//Then the other
						(antMath_p[j])->applyVP(subim2, subim2, direction2_p);
						
						/*
						//one antenna
						(antMath_p[k])->applyPB(subim2, subim2, direction1_p);
						//Then the other
						(antMath_p[j])->applyPB(subim2, subim2, direction2_p);
						*/
					}
                    //tim.show("After Apply2 ");
                    //tim.mark();
					//wtime1+=omp_get_wtime()-wtime0;
                    //subim.copyData((LatticeExpr<Complex>) (iif(abs(subim)> 5e-2, subim, 0)));
                    //subim2.copyData((LatticeExpr<Complex>) (iif(abs(subim2)> 25e-4, subim2, 0)));

					//wtime0=omp_get_wtime();
					ft_p.c2cFFTInDouble(subim);
					ft_p.c2cFFTInDouble(subim2);
					//ft_p.c2cFFT(subim);
					//ft_p.c2cFFT(subim2);
					//wtime2+=omp_get_wtime()-wtime0;
                    //  tim.show("after ffts ");


                }
		//cerr << "Apply " << wtime1 << "  fft " << wtime2 << endl;
                /*
                if(nBeamChans >1){
                  Slicer slin(blcin, trcin, Slicer::endIsLast);
                  SubImage<Complex> origPB(pBScreen, slin, false);
                  IPosition elshape= origPB.shape();
                  Matrix<Complex> i1=origPB.get(true);
                  SubImage<Complex> origPB2(pB2Screen, slin, false);
                  Matrix<Complex> i2=origPB2.get(true);
                  Int cenX=i1.shape()(0)/2;
                  Int cenY=i1.shape()(1)/2;

                  for (Int kk=0; kk < (nBeamChans-1); ++kk){
                    Double fratio=beamFreqs(kk)/beamFreqs(nBeamChans-1);
                    cerr << "fratio " << fratio << endl;
                    blcin[3]=kk;
                    trcin[3]=kk;
                    //Slicer slout(blcin, trcin, Slicer::endIsLast);
                    Matrix<Complex> o1(i1.shape(), Complex(0.0));
                    Matrix<Complex> o2(i2.shape(), Complex(0.0));
                    for (Int yy=0;  yy < i1.shape()(1); ++yy){
                      //Int nyy= (Double(yy-cenY)*fratio) + cenY;
                      Double nyy= (Double(yy-cenY)/fratio) + cenY;
                      Double cyy=ceil(nyy);
                      Double fyy= floor(nyy);
                      Int iy=nyy > fyy+0.5 ? cyy : fyy;
                      if(cyy <2*cenY && fyy >=0.0)
                 for(Int xx=0; xx < i1.shape()(0); ++ xx){
                   //Int nxx= Int(Double(xx-cenX)*fratio) + cenX;
                   Double nxx= Int(Double(xx-cenX)/fratio) + cenX;
                    Double cxx=ceil(nxx);
                    Double fxx= floor(nxx);
                    Int ix=nxx > fxx+0.5 ? cxx : fxx ;
                   if(cxx < 2*cenX && fxx >=0.0 ){
                     //Double dist=sqrt((nxx-cxx)*(nxx-cxx)+(nyy-cyy)*(nyy-cyy))/sqrt(2.0);
                     //o1(xx, yy)=float(1-dist)*i1(fxx, fyy)+ dist*i1(cxx,cyy);
                     o1(xx, yy)=i1( ix, iy);
                     //o2(xx, yy)=i2(nxx, nyy);
                     //o2(xx, yy)=float(1-dist)*i2(fxx, fyy)+ dist*i2(cxx,cyy);
                     o2(xx, yy)=i2(ix, iy);
                   }
                 }
                    }
                    pBScreen.putSlice(o1.reform(elshape), blcin);
                    pB2Screen.putSlice(o2.reform(elshape), blcin);
                  }

                }
                */

                //tim.show("after apply+apply2+masking+fft ");
                //tim.mark();
                //LatticeFFT::cfft2d(pBScreen);
                //LatticeFFT::cfft2d(pB2Screen);

                //Matrix<Complex> lala=pBScreen.get(true);
                //fft.fft0(lala, true);
                //fft.flip(lala, false, false);
                // pBScreen.put(lala.reform(IPosition(4, convSize_p, convSize_p, 1, 1)));
                //lala=pB2Screen.get(true);
                //fft.fft0(lala, true);
                //fft.flip(lala, false, false);
                //pB2Screen.put(lala.reform(IPosition(4, convSize_p, convSize_p, 1, 1)));

                //////////*****************
                /*if(0){
                  ostringstream os1;
                  os1 << "PB_field_" << Int(thePix_p[0]) << "_" << Int(thePix_p[1]) << "_antpair_" << k <<"_"<<j ;
                  PagedImage<Float> thisScreen(pbShape, coords, String(os1));
                  LatticeExpr<Float> le(abs(pBScreen));
                  thisScreen.copyData(le);
                  }*/
                ////////*****************/

                //tim.show("after FFT ");
                //tim.mark();
                Int plane=0;
                for (uInt jj=0; jj < k; ++jj)
                    plane=plane+ndish-jj-1;
                plane=plane+j;
                begin[4]=plane;
                end[4]=plane;
                Slicer slplane(begin, end, Slicer::endIsLast);
                //cerr <<  "SHAPES " << convFunc_p(begin, end).shape() << "  " << pBScreen.get(false).shape() << " begin and end " << begin << "    " << end << endl;
                //convFunc_p(begin, end).copyMatchingPart(pBScreen.get(false));
                //weightConvFunc_p(begin, end).copyMatchingPart(pB2Screen.get(false));
                IPosition blcQ(4, pbShape(0)/8*3, pbShape(1)/8*3, 0, 0);
                IPosition trcQ(4,  blcQ[0]+ pbShape(0)/4-1, blcQ[1]+pbShape(1)/4-1 , nBeamPols-1, nBeamChans-1);

                //cerr << "blcQ " << blcQ << " trcQ " << trcQ << " pbShape " << pbShape << endl;
                Slicer slQ(blcQ, trcQ, Slicer::endIsLast);
                {
                    SubImage<Complex>  pBSSub(pBScreen, slQ, false);
                    SubLattice<Complex> cFTempSub(convFuncTemp,  slplane, true);
                    LatticeConcat<Complex> lc(4);
                    lc.setLattice(pBSSub);
                    //cerr << "SHAPES " << cFTempSub.shape() << "   " <<  lc.shape() << endl;
                    cFTempSub.copyData(lc);
                    //cFTempSub.copyData(pBScreen);
                }
                {
                    SubImage<Complex>  pB2SSub(pB2Screen, slQ, false);
                    SubLattice<Complex> cFTempSub2(weightConvFuncTemp,  slplane, true);
                    LatticeConcat<Complex> lc(4);
                    lc.setLattice(pB2SSub);
                    cFTempSub2.copyData(lc);
                    // cFTempSub2.copyData(pB2Screen);
                    //weightConvFuncTemp.putSlice(pB2Screen.get(false), begin);

                }
                //	  supportAndNormalize(plane, convSampling);
                supportAndNormalizeLatt( plane, convSampling, convFuncTemp,  weightConvFuncTemp);



                // tim.show("After search of support ");
            }

        }


        doneMainConv_p[actualConvIndex_p]=true;
        convFunctions_p.resize(actualConvIndex_p+1);
        convWeights_p.resize(actualConvIndex_p+1);
        convSupportBlock_p.resize(actualConvIndex_p+1);
	//Using conjugate change support to be larger of either
	if((nchan_p == 1) && getConjConvFunc) {
	  Int conjsupp=conjSupport(beamFreqs) ;
	  if(conjsupp > max(convSupport_p)){
	    convSupport_p=Int(1.5*Float(conjsupp));
	  }

	}
        convFunctions_p[actualConvIndex_p]= new Array<Complex>();
        convWeights_p[actualConvIndex_p]= new Array<Complex>();
        convSupportBlock_p[actualConvIndex_p]=new Vector<Int>();
        Int newConvSize=2*(max(convSupport_p)+2)*convSampling;
        Int newRealConvSize=newConvSize* Int(Double(convSamp)/Double(convSampling));
        Int lattSize=convFuncTemp.shape()(0);
        (*convSupportBlock_p[actualConvIndex_p])=convSupport_p;
		//cerr << "convsupport " << convSupport_p << endl;

        if(newConvSize < lattSize) {
            IPosition blc(5, (lattSize/2)-(newConvSize/2),
                          (lattSize/2)-(newConvSize/2),0,0,0);
            IPosition trc(5, (lattSize/2)+(newConvSize/2-1),
                          (lattSize/2)+(newConvSize/2-1), nBeamPols-1, nBeamChans-1,ndishpair-1);
            IPosition shp(5, newConvSize, newConvSize, nBeamPols, nBeamChans, ndishpair);

            convFunctions_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, newRealConvSize, newRealConvSize, nBeamPols, nBeamChans, ndishpair ));
            convWeights_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, newRealConvSize, newRealConvSize, nBeamPols, nBeamChans, ndishpair ));
            (*convFunctions_p[actualConvIndex_p])=resample(convFuncTemp.getSlice(blc,shp),Double(convSamp)/Double(convSampling));
            convSize_p=newRealConvSize;
            (*convWeights_p[actualConvIndex_p])=resample(weightConvFuncTemp.getSlice(blc, shp),Double(convSamp)/Double(convSampling));
	    //cerr << "nchan " << nchan_p << " getconj " << getConjConvFunc << endl;
       
        }
        else {
            newRealConvSize=lattSize* Int(Double(convSamp)/Double(convSampling));
            convFunctions_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, newRealConvSize, newRealConvSize, nBeamPols, nBeamChans, ndishpair ));
            convWeights_p[actualConvIndex_p]= new Array<Complex>(IPosition(5, newRealConvSize, newRealConvSize, nBeamPols, nBeamChans, ndishpair ));

            (*convFunctions_p[actualConvIndex_p])=resample(convFuncTemp.get(),  Double(convSamp)/Double(convSampling));
            (*convWeights_p[actualConvIndex_p])=resample(weightConvFuncTemp.get(),  Double(convSamp)/Double(convSampling));
            convSize_p=newRealConvSize;
        }


	if((nchan_p == 1) && getConjConvFunc) {
	  fillConjConvFunc(beamFreqs);
	  /////////////////////////TESTOOO
	  /*PagedImage<Complex> SCREEN2(TiledShape(convFunctions_p[actualConvIndex_p]->shape()), TMP, "CONJU"+String::toString(actualConvIndex_p));
	  SCREEN2.put(*convFunctionsConjFreq_p[actualConvIndex_p]  );
	  */
	  ////////////////////////
	}

        convFunc_p.resize();
        weightConvFunc_p.resize();

    }
    else {
        convSize_p=max(*convSizes_p[actualConvIndex_p]);
        convSupport_p.resize();
        convSupport_p=*convSupportBlock_p[actualConvIndex_p];
    }
    /*
    rowMap.resize(vb.nRow());
    for (Int k=0; k < vb.nRow(); ++k){
      //plane of convfunc that match this pair of antennas
      rowMap(k)=antIndexToDiamIndex_p(vb.antenna1()(k))*ndish+
    antIndexToDiamIndex_p(vb.antenna2()(k));

    }
    */
    ////////////////TESTOOO
    //		 CoordinateSystem TMP = coords;
    //	  CoordinateUtil::addLinearAxes(TMP, Vector<String>(1,"gulu"), IPosition(1,nBeamChans)); 
    //	  PagedImage<Complex> SCREEN(TiledShape(convFunctions_p[actualConvIndex_p]->shape()), TMP, "NONCONJUVI2"+String::toString(actualConvIndex_p));
    //	  SCREEN.put(*convFunctions_p[actualConvIndex_p]  );
    //	   PagedImage<Complex> SCREEN3(TiledShape(convWeights_p[actualConvIndex_p]->shape()), TMP, "FTWEIGHTVI2"+String::toString(actualConvIndex_p));
    //	  SCREEN3.put(*convWeights_p[actualConvIndex_p]  );
	
    /////////////////

    makerowmap(vb, convFuncRowMap);
    ///need to deal with only the maximum of different baselines available in this
    ///vb
    ndishpair=max(convFuncRowMap)+1;

    convSupportBlock_p.resize(actualConvIndex_p+1);
    convSizes_p.resize(actualConvIndex_p+1);
    //convSupportBlock_p[actualConvIndex_p]=new Vector<Int>(ndishpair);
    //(*convSupportBlock_p[actualConvIndex_p])=convSupport_p;
    convSizes_p[actualConvIndex_p]=new Vector<Int> (ndishpair);

    /*    convFunctions_p[actualConvIndex_p]->resize(convSize_p, convSize_p, ndishpair);
    *(convFunctions_p[actualConvIndex_p])=convSave_p;
    convWeights_p[actualConvIndex_p]->resize(convSize_p, convSize_p, ndishpair);
    *(convWeights_p[actualConvIndex_p])=weightSave_p;
    */

    convFunc_p.resize();
	if((nchan_p == 1) && getConjConvFunc) {
	  // cerr << this << " recovering " << actualConvIndex_p <<  "   " <<convFunctionsConjFreq_p.size() << endl;
	  if(Int(convFunctionsConjFreq_p.size()) <= actualConvIndex_p){
	    fillConjConvFunc(beamFreqs);
	    
	  }
		convFunc_p=(*convFunctionsConjFreq_p[actualConvIndex_p]);
	}
	else{
		
		convFunc_p=(*convFunctions_p[actualConvIndex_p]);
	}
	
    weightConvFunc_p.resize();
    weightConvFunc_p=(*convWeights_p[actualConvIndex_p]);


    // cerr << "convfunc shapes " <<  convFunc_p.shape() <<  "   " << weightConvFunc_p.shape() << "  " << convSize_p << " pol " << nBeamPols << "  chan " << nBeamChans << " ndishpair " << ndishpair << endl;
    //convSupport_p.resize();
    //convSupport_p=(*convSupportBlock_p[actualConvIndex_p]);
    Bool delc;
    Bool delw;
    Double dirX=pixFieldDir(0);
    Double dirY=pixFieldDir(1);
    Complex *convstor=convFunc_p.getStorage(delc);
    Complex *weightstor=weightConvFunc_p.getStorage(delw);
    Int elconvsize=convSize_p;
    //cerr << "dirX Y "<< dirX << "   " << dirY << "  convsize " << convSize_p << endl; 
       #pragma omp parallel default(none) firstprivate(convstor, weightstor, dirX, dirY, elconvsize, ndishpair, nBeamChans, nBeamPols)
    {
       #pragma omp for
        for(Int iy=0; iy<elconvsize; ++iy) {
            applyGradientToYLine(iy,  convstor, weightstor, dirX, dirY, elconvsize, ndishpair, nBeamChans, nBeamPols);

        }
    }///End of pragma

    convFunc_p.putStorage(convstor, delc);
    weightConvFunc_p.putStorage(weightstor, delw);



    //For now all have the same size convsize;
    convSizes_p[actualConvIndex_p]->set(convSize_p);

    //We have to get the references right now
    //    convFunc_p.resize();
    //convFunc_p.reference(*convFunctions_p[actualConvIndex_p]);
    //weightConvFunc_p.resize();
    //weightConvFunc_p.reference(*convWeights_p[actualConvIndex_p]);

    convFunc.reference(convFunc_p);
    weightConvFunc.reference(weightConvFunc_p);
    convsize=*convSizes_p[actualConvIndex_p];
    convSupport=convSupport_p;
    //cerr << "convsize " << convsize << " convsupp " << convSupport << endl; 
     ////////////////TESTOOO
     //CoordinateSystem TMP = coords;
     //CoordinateUtil::addLinearAxes(TMP, Vector<String>(1,"gulu"), IPosition(1,convFunc.shape()(4))); 
     //PagedImage<Complex> SCREEN(TiledShape(convFunctions_p[actualConvIndex_p]->shape()), TMP, "convFunc_orig"+String::toString(actualConvIndex_p));
     //SCREEN.put(convFunc_p  );
    //	   PagedImage<Complex> SCREEN3(TiledShape(convWeights_p[actualConvIndex_p]->shape()), TMP, "FTWEIGHTVI2"+String::toString(actualConvIndex_p));
    //	  SCREEN3.put(*convWeights_p[actualConvIndex_p]  );
    //cerr  << " supp and size "<< max(convSupport) << "  " << max(convsize) << endl;

}


void HetArrayConvFunc::applyGradientToYLine(const Int iy, Complex*& convFunctions, Complex*& convWeights, const Double pixXdir, const Double pixYdir, Int convSize, const Int ndishpair, const Int nChan, const Int nPol) {
    Double cy, sy;

    SINCOS(Double(iy-convSize/2)*pixYdir, sy, cy);
    Complex phy(cy,sy) ;
    for (Int ix=0; ix<convSize; ix++) {
        Double cx, sx;
        SINCOS(Double(ix-convSize/2)*pixXdir, sx, cx);
        Complex phx(cx,sx) ;
        for (Int ipol=0; ipol< nPol; ++ipol) {
            //Int poloffset=ipol*nChan*ndishpair*convSize*convSize;
            for (Int ichan=0; ichan < nChan; ++ichan) {
                //Int chanoffset=ichan*ndishpair*convSize*convSize;
                for(Int iz=0; iz <ndishpair; ++iz) {
                    ooLong index=((ooLong(iz*nChan+ichan)*nPol+ipol)*ooLong(convSize)+ooLong(iy))*ooLong(convSize)+ooLong(ix);
                    convFunctions[index]= convFunctions[index]*phx*phy;
                    convWeights[index]= convWeights[index]*phx*phy;
                }
            }
        }

    }
}
Int  HetArrayConvFunc::conjSupport(const casacore::Vector<casacore::Double>& freqs){
  Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
  Double maxRatio=-1.0;
  //cerr << "freqs.nelem " << freqs.nelements() << "   " <<convSupport_p.nelements() << endl;
  for (Int k=0; k < freqs.shape()[0]; ++k) {
    //Double conjFreq=(centerFreq-freqs[k])+centerFreq;
    Double conjFreq=SynthesisUtils::conjFreq(freqs[k], centerFreq);
    //cerr << k << "  " << conjFreq/freqs[k]*convSupport_p[k] << "   " ;
    if(maxRatio < conjFreq/freqs[k] )      
      maxRatio=conjFreq/freqs[k];
  }
  //return  Int(max(convSupport_p)*sqrt(maxRatio)/2.0)*2;
  return  Int(max(convSupport_p)*(maxRatio)/2.0)*2;
}
void HetArrayConvFunc::fillConjConvFunc(const Vector<Double>& freqs) {
    //cerr << "Actualconv index " << actualConvIndex_p << endl;
    convFunctionsConjFreq_p.resize(actualConvIndex_p+1);
    //Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
    IPosition shp=convFunctions_p[actualConvIndex_p]->shape();
    convFunctionsConjFreq_p[actualConvIndex_p]=new Array<Complex>(shp, Complex(0.0));
    IPosition blc(5,0,0,0,0,0);
    IPosition trc=shp-1;
    for (Int k=0; k < freqs.shape()[0]; ++k) {
        blc[3]=k;
        trc[3]=k;
        Array<Complex> convSlice((*convFunctions_p[actualConvIndex_p])(blc, trc));
 
        Array<Complex> conjSlice=(*convFunctionsConjFreq_p[actualConvIndex_p])(blc, trc);
	fillConjConvFunc(conjSlice, convSlice, freqs[k]);
    }   
}

  void HetArrayConvFunc::fillConjConvFunc(Array<Complex>&out, const Array<Complex>& in, const Double& freq){

   Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
   Double conjFreq=SynthesisUtils::conjFreq(freq, centerFreq);
   Double ratio1= Double(Int(Double(in.shape()(0))*conjFreq/freq/2.0)*2)/Double(in.shape()(0));
   Double ratio2= Double(Int(Double(in.shape()(1))*conjFreq/freq/2.0)*2)/Double(in.shape()(1));

   IPosition shp=in.shape();
   Array<Complex> conjFreqSlice(resample(in, conjFreq/freq));
        if(conjFreq > freq) {
            IPosition end=shp-1;
            IPosition beg(5,0,0,0,0,0);
            beg(0)=(conjFreqSlice.shape()(0)-shp(0))/2;
            beg(1)=(conjFreqSlice.shape()(1)-shp(1))/2;
	    end(0)=beg(0)+shp(0)-1;
	    end(1)=beg(1)+shp(1)-1;
	    if(shp.nelements() >3)
	      end[3]=0;
            out=conjFreqSlice(beg, end);
        }
        else {
	  IPosition end=conjFreqSlice.shape()-1;
	  if(shp.nelements() >3)
	    end[3]=0;
	  IPosition beg(shp.nelements(),0);
	  beg(0)=(shp(0)-conjFreqSlice.shape()(0))/2;
	  beg(1)=(shp(1)-conjFreqSlice.shape()(1))/2;
	  end(0)+=beg(0);
	  end(1)+=beg(1);
	  out.set(0);
	  out(beg, end)=conjFreqSlice;
        }
	//cerr << "SUMS " << sum(real(in)) << "   new " << sum(real(out))/ratio1/ratio2 << endl; //" weight " << sum(real(weightSlice))/ratio1/ratio2<< endl;
	Complex renorm( 1.0/(ratio1*ratio2),0.0);
	out=out*renorm;

  }

  
Bool HetArrayConvFunc::toRecord(RecordInterface& rec) {

    try {
        rec.define("name", "HetArrayConvFunc");
        Int numConv=convFunctions_p.nelements();
        rec.define("ndefined", numConv);
        //rec.define("convfunctionmap", convFunctionMap_p);
        std::map<String, Int>::iterator it=vbConvIndex_p.begin();
        for (Int64 k=0; k < numConv; ++k) {
            rec.define("convfunctions"+String::toString(k), *(convFunctions_p[k]));
            rec.define("convweights"+String::toString(k), *(convWeights_p[k]));
            rec.define("convsizes"+String::toString(k), *(convSizes_p[k]));
            rec.define("convsupportblock"+String::toString(k), *(convSupportBlock_p[k]));
            rec.define(String("key")+String::toString(k), it->first);
            rec.define(String("val")+String::toString(k), it->second);
            it++;
        }
        rec.define("actualconvindex",  actualConvIndex_p);
        rec.define("donemainconv", doneMainConv_p);
        rec.define("vptable", vpTable_p);
        rec.define("pbclass", Int(pbClass_p));
	rec.define("usepointingtable", usePointingTable_p);

    }
    catch(AipsError& x) {
        return false;
    }
    return true;

}

Bool HetArrayConvFunc::fromRecord(String& err, const RecordInterface& rec, Bool calcfluxscale) {
    //Force pbmath stuff and saved image stuff
    nchan_p=0;
    msId_p=-1;
    try {
        if(!rec.isDefined("name") || rec.asString("name") != "HetArrayConvFunc") {
            throw(AipsError("Wrong record to recover HetArray from"));
        }
        nDefined_p=rec.asInt("ndefined");
        //rec.get("convfunctionmap", convFunctionMap_p);
        convFunctions_p.resize(nDefined_p, true, false);
        convSupportBlock_p.resize(nDefined_p, true, false);
        convWeights_p.resize(nDefined_p, true, false);
        convSizes_p.resize(nDefined_p, true, false);
        vbConvIndex_p.erase(vbConvIndex_p.begin(), vbConvIndex_p.end());
        for (Int64 k=0; k < nDefined_p; ++k) {
            convFunctions_p[k]=new Array<Complex>();
            convWeights_p[k]=new Array<Complex>();
            convSizes_p[k]=new Vector<Int>();
            convSupportBlock_p[k]=new Vector<Int>();
            rec.get("convfunctions"+String::toString(k), *(convFunctions_p[k]));
            rec.get("convweights"+String::toString(k), *(convWeights_p[k]));
            rec.get("convsizes"+String::toString(k), *(convSizes_p[k]));
            rec.get("convsupportblock"+String::toString(k), *(convSupportBlock_p[k]));
            String key;
            Int val;
            rec.get(String("key")+String::toString(k), key);
            rec.get(String("val")+String::toString(k), val);
            vbConvIndex_p[key]=val;
        }
        //Now that we are calculating all phase gradients on the fly ..
        ///we should clean up some and get rid of the cached variables

        convSize_p= nDefined_p > 0 ? (*(convSizes_p[0]))[0] : 0;
        //convSave_p.resize();
        //rec.get("convsave", convSave_p);
        //weightSave_p.resize();
        //rec.get("weightsave", weightSave_p);
        rec.get("vptable", vpTable_p);
        rec.get("donemainconv", doneMainConv_p);
	rec.get("usepointingtable", usePointingTable_p);
        //convSupport_p.resize();
        //rec.get("convsupport", convSupport_p);
        pbClass_p=static_cast<PBMathInterface::PBClass>(rec.asInt("pbclass"));
        calcFluxScale_p=calcfluxscale;
	///Generate the sincCache now as inside the multithread  part it may have race issues
	initSincCache();
    }
    catch(AipsError& x) {
        err=x.getMesg();
        return false;
    }

    return true;
}


void HetArrayConvFunc::supportAndNormalize(Int plane, Int convSampling) {

    LogIO os;
    os << LogOrigin("HetArrConvFunc", "suppAndNorm")  << LogIO::NORMAL;
    // Locate support
    Int convSupport=-1;
    IPosition begin(5, 0, 0, 0, 0, plane);
    IPosition end(5, convFunc_p.shape()[0]-1,  convFunc_p.shape()[1]-1, 0, 0, plane);
    Matrix<Complex> convPlane(convFunc_p(begin, end).reform(IPosition(2,convFunc_p.shape()[0], convFunc_p.shape()[1]))) ;
    Float maxAbsConvFunc=max(amplitude(convPlane));
    Float minAbsConvFunc=min(amplitude(convPlane));
    Bool found=false;
    Int trial=0;
    for (trial=convSize_p/2-2; trial>0; trial--) {
        //Searching down a diagonal
        if(abs(convPlane(convSize_p/2-trial,convSize_p/2-trial)) >  (1.0e-2*maxAbsConvFunc) ) {
            found=true;
            trial=Int(sqrt(2.0*Float(trial*trial)));
	   
            break;
        }
    }
    if(!found) {
        if((maxAbsConvFunc-minAbsConvFunc) > (1.0e-2*maxAbsConvFunc))
            found=true;
        // if it drops by more than 2 magnitudes per pixel
        trial=( (10*convSampling) < convSize_p) ? 5*convSampling : (convSize_p/2 - 4*convSampling);
    }


    if(found) {
        if(trial < 5*convSampling)
            trial= ( (10*convSampling) < convSize_p) ? 5*convSampling : (convSize_p/2 - 4*convSampling);
        convSupport=Int(0.5+Float(trial)/Float(convSampling))+1;
        //support is really over the edge
        if( (convSupport*convSampling) >= convSize_p/2) {
            convSupport=convSize_p/2/convSampling-1;
        }
    }
    else {
        /*
        os << "Convolution function is misbehaved - support seems to be zero\n"
           << "Reasons can be: \nThe image definition not covering one or more of the pointings selected \n"
           << "Or no unflagged data in a given pointing"

           << LogIO::EXCEPTION;
        */
        //OTF may have flagged stuff ...
        convSupport=0;
    }
    //cerr << "trial " << trial << " convSupport " << convSupport << " convSize " << convSize_p << endl;
    convSupport_p(plane)=convSupport;
    Double pbSum=0.0;
    /*
    Double pbSum1=0.0;

    for (Int iy=-convSupport;iy<=convSupport;iy++) {
      for (Int ix=-convSupport;ix<=convSupport;ix++) {
        Complex val=convFunc_p.xyPlane(plane)(ix*convSampling+convSize_p/2,
    					  iy*convSampling+convSize_p/2);

        pbSum1+=sqrt(real(val)*real(val)+ imag(val)*imag(val));
      }
    }

    */
    if(convSupport >0) {
        IPosition blc(2, -convSupport*convSampling+convSize_p/2, -convSupport*convSampling+convSize_p/2);
        IPosition trc(2, convSupport*convSampling+convSize_p/2, convSupport*convSampling+convSize_p/2);
        for (Int chan=0; chan < convFunc_p.shape()[3]; ++chan) {
            begin[3]=chan;
            end[3]=chan;
            convPlane.resize();
            convPlane.reference(convFunc_p(begin, end).reform(IPosition(2,convFunc_p.shape()[0], convFunc_p.shape()[1])));
            pbSum=real(sum(convPlane(blc,trc)))/Double(convSampling)/Double(convSampling);
            if(pbSum>0.0) {
                (convPlane)=convPlane*Complex(1.0/pbSum,0.0);
                convPlane.resize();
                convPlane.reference(weightConvFunc_p(begin, end).reform(IPosition(2,convFunc_p.shape()[0], convFunc_p.shape()[1])));
		 
                (convPlane) =(convPlane)*Complex(1.0/pbSum,0.0);
            }
            else {
                os << "Convolution function integral is not positive"
                   << LogIO::EXCEPTION;
            }
        }
    }
    else {
        //no valid convolution for this pointing
        for (Int chan=0; chan < convFunc_p.shape()[3]; ++chan) {
            begin[3]=chan;
            end[3]=chan;
            convFunc_p(begin, end).set(Complex(0.0));
            weightConvFunc_p(begin, end).set(Complex(0.0));
            //convFunc_p.xyPlane(plane).set(0.0);
            //weightConvFunc_p.xyPlane(plane).set(0.0);
        }
    }

}


  Bool HetArrayConvFunc::supportAndNormalize(Array<Complex>& ftPB, Array<Complex>& ftWt, Vector<Int>& support, const Double supportFactor){
    Int convSize=ftPB.shape()(0);
    Int npol=ftPB.shape()(2);
    Int nchan=ftPB.shape()(3);
    Int nbas=ftPB.shape()(4);
    support.resize(nbas);
    support=0;
    Bool retval=True;
    Float cutlevel=2.5e-2;
    //numeric needs a larger ft
    for (uInt k=0; k < antMath_p.nelements() ; ++k){
      if((antMath_p[k]->whichPBClass()) == PBMathInterface::NUMERIC)
	cutlevel=5e-3;
    }
    ArrayIterator<Complex> pbIt(ftPB, IPosition(2, 0,1));
    ArrayIterator<Complex> wtIt(ftWt,  IPosition(2, 0,1));
    while(!pbIt.pastEnd()){
      Int plane=wtIt.pos()(4);
      Bool found=false;
      Int trial=0;
      Matrix<Complex> wtplane(wtIt.array());
      Matrix<Complex> pbplane(pbIt.array());
      Float maxFunc, minFunc;
      IPosition minpos, maxpos;
      minMax(minFunc, maxFunc, minpos, maxpos, amplitude(wtplane));
      for (trial=0; trial< (convSize-max(maxpos.asVector())-2); ++trial) {
      ///largest along either axis
      //cerr << "rat1 " << abs(convPlane(maxpos[0]-trial,maxpos[1]))/maxAbsConvFunc << " rat2 " << abs(convPlane(maxpos[0],maxpos[1]-trial))/maxAbsConvFunc << endl;
	if((abs(wtplane(maxpos[0]-trial, maxpos[1])) <  (cutlevel*maxFunc)) &&(abs(wtplane(maxpos[0],maxpos[1]-trial)) <  (cutlevel*maxFunc)) )
	  {
            found=true;
            //trial=Int(sqrt(2.0*Float(trial*trial)));
	    
            break;
	  }
      }
      if(found){
	///might have to do it for every frequency channel
	Int suppThisPlane=(Int(0.5+Float(trial)))+1;
	//support is really over the edge
        if( (suppThisPlane >= convSize/2)) {
            suppThisPlane=convSize/2-1;
        }
	support(plane) =max(support(plane), suppThisPlane) ;
	//cerr << "convsupp " << convSupport << endl;
	IPosition blc(2,-suppThisPlane+convSize/2, -suppThisPlane+convSize/2);
	IPosition trc(2, suppThisPlane+convSize/2, suppThisPlane+convSize/2);
	Double pbSum=0.0;
	pbSum=real(sum(pbplane(blc,trc)));
	if(pbSum > 0.0){
	  pbplane=pbplane*Complex(1.0/pbSum, 0.0);
	  Double pbSum1=real(sum(wtplane(blc,trc)));
	  wtplane=wtplane*Complex(1.0/pbSum1, 0.0);
	}
	else {
	  //cerr << "PBSUM " << pbSum << " blc, trc " << blc <<  "   " << trc << endl;
	  //throw(AipsError("Convolution function integral is not positive"));
	  support(plane)=0;
	  retval=False;
	}
	//cerr << "pbplane sum " << sum(pbplane(blc,trc)) << "   "  << real(sum(pbplane(blc,trc))) << "   " <<(sum(pbIt.array()))<< endl;
      }
      else{
	support(plane)=0;
	retval=False;
      }
      
      pbIt.next();
      wtIt.next();
    }

    if(max(support)==0)
      return False;
    Int  newConvSize=std::round(min(Double(2*max(support+2))*supportFactor, Double(convSize))/2.0)*2;
    IPosition blc(5, convSize/2-newConvSize/2, convSize/2-newConvSize/2, 0, 0, 0);
    IPosition trc(5, convSize/2+newConvSize/2-1, convSize/2+newConvSize/2-1, npol-1, nchan-1, nbas-1);
    //cerr << "SUPP blc trc " << blc << "  " << trc << " szie " << ftPB.shape() << endl;
    Array<Complex> tempArr;
    tempArr.assign(ftPB(blc, trc));
    ftPB.assign(tempArr);
    tempArr.assign(ftWt(blc, trc));
    ftWt.assign(tempArr);
    return retval;
  }
  
void HetArrayConvFunc::supportAndNormalizeLatt(Int plane, Int convSampling, TempLattice<Complex>& convFuncLat,
        TempLattice<Complex>& weightConvFuncLat) {

    LogIO os;
    os << LogOrigin("HetArrConvFunc", "suppAndNorm")  << LogIO::NORMAL;
    // Locate support
    Int convSupport=-1;
    ///Use largest channel as highest freq thus largest conv func
    IPosition begin(5, 0, 0, 0, convFuncLat.shape()(3)-1, plane);
    IPosition shape(5, convFuncLat.shape()[0],  convFuncLat.shape()[1], 1, 1, 1);
    //Int convSize=convSize_p;
    Int convSize=shape(0);
    ///use FT weightconvlat as it is wider
    Matrix<Complex> convPlane=weightConvFuncLat.getSlice(begin, shape, true);
    Float maxAbsConvFunc, minAbsConvFunc;
    IPosition minpos, maxpos;
    minMax(minAbsConvFunc, maxAbsConvFunc, minpos, maxpos, amplitude(convPlane));
     Bool found=false;
    Int trial=0;
    Float cutlevel=2.5e-2;
    //numeric needs a larger ft
    for (uInt k=0; k < antMath_p.nelements() ; ++k){
      if((antMath_p[k]->whichPBClass()) == PBMathInterface::NUMERIC)
	cutlevel=5e-3;
    }

    for (trial=0; trial< (convSize-max(maxpos.asVector())-2); ++trial) {
      ///largest along either axis
      //cerr << "rat1 " << abs(convPlane(maxpos[0]-trial,maxpos[1]))/maxAbsConvFunc << " rat2 " << abs(convPlane(maxpos[0],maxpos[1]-trial))/maxAbsConvFunc << endl;
      if((abs(convPlane(maxpos[0]-trial, maxpos[1])) <  (cutlevel*maxAbsConvFunc)) &&(abs(convPlane(maxpos[0],maxpos[1]-trial)) <  (cutlevel*maxAbsConvFunc)) )
	{

            found=true;
            //trial=Int(sqrt(2.0*Float(trial*trial)));
	    
            break;
        }
    }
    if(!found) {
      if((maxAbsConvFunc-minAbsConvFunc) > (cutlevel*maxAbsConvFunc))
            found=true;
        // if it drops by more than 2 magnitudes per pixel
        //trial=( (10*convSampling) < convSize) ? 5*convSampling : (convSize/2 - 4*convSampling);
      trial=convSize/2 - 4*convSampling;
    }

    if(found) {
        if(trial < 5*convSampling)
            trial= ( (10*convSampling) < convSize) ? 5*convSampling : (convSize/2 - 4*convSampling);
        convSupport=(Int(0.5+Float(trial)/Float(convSampling)))+1 ;
	//cerr << "convsupp " << convSupport << endl;
        //support is really over the edge
        if( (convSupport*convSampling) >= convSize/2) {
            convSupport=convSize/2/convSampling-1;
        }
    }
    else {
        /*
        os << "Convolution function is misbehaved - support seems to be zero\n"
           << "Reasons can be: \nThe image definition not covering one or more of the pointings selected \n"
           << "Or no unflagged data in a given pointing"

           << LogIO::EXCEPTION;
        */
        //OTF may have flagged stuff ...
        convSupport=0;
    }
    convSupport_p(plane)=convSupport;
    Double pbSum=0.0;
    /*
    Double pbSum1=0.0;

    for (Int iy=-convSupport;iy<=convSupport;iy++) {
      for (Int ix=-convSupport;ix<=convSupport;ix++) {
        Complex val=convFunc_p.xyPlane(plane)(ix*convSampling+convSize_p/2,
    					  iy*convSampling+convSize_p/2);

        pbSum1+=sqrt(real(val)*real(val)+ imag(val)*imag(val));
      }
    }

    */
    //cerr << "convSize_p " << convSize_p <<  " convSize " << convSize << endl;
    if(convSupport >0) {
        IPosition blc(2, -convSupport*convSampling+convSize/2, -convSupport*convSampling+convSize/2);
        IPosition trc(2, convSupport*convSampling+convSize/2, convSupport*convSampling+convSize/2);
        for (Int chan=0; chan < convFuncLat.shape()[3]; ++chan) {
            begin[3]=chan;
            //end[3]=chan;
            convPlane.resize();
            convPlane=convFuncLat.getSlice(begin, shape, true);
            pbSum=real(sum(convPlane(blc,trc)))/Double(convSampling)/Double(convSampling);
            if(pbSum>0.0) {
                (convPlane)=convPlane*Complex(1.0/pbSum,0.0);
                convFuncLat.putSlice(convPlane, begin);
                convPlane.resize();
                convPlane=weightConvFuncLat.getSlice(begin, shape, true);
		Double pbSum1=0.0;
		pbSum1=real(sum(convPlane(blc,trc)))/Double(convSampling)/Double(convSampling);
                (convPlane) =(convPlane)*Complex(1.0/pbSum1,0.0);
                weightConvFuncLat.putSlice(convPlane, begin);
            }
            else {
                os << "Convolution function integral is not positive"
                   << LogIO::EXCEPTION;
            }
        }
    }
    else {
        //no valid convolution for this pointing
        for (Int chan=0; chan < convFuncLat.shape()[3]; ++chan) {
            begin[3]=chan;
            //end[3]=chan;
            convPlane.resize(shape[0], shape[1]);
            convPlane.set(Complex(0.0));
            convFuncLat.putSlice(convPlane, begin);
            weightConvFuncLat.putSlice(convPlane, begin);
            //convFunc_p.xyPlane(plane).set(0.0);
            //weightConvFunc_p.xyPlane(plane).set(0.0);
        }
    }

}

Int HetArrayConvFunc::factorial(Int n) {
    Int fact=1;
    for (Int k=1; k<=n; ++k)
        fact *=k;
    return fact;
}


Int HetArrayConvFunc::checkPBOfField(const vi::VisBuffer2& vb,
                                     Vector<Int>& /*rowMap*/, const MVDirection& extraShift, const Bool useExtraShift) {

  toPix(vb, extraShift, useExtraShift);
    Vector<Int> pixdepoint(2);
    convertArray(pixdepoint, thePix_p);
    if((pixdepoint(0) < 0) ||  pixdepoint(0) >= nx_p || pixdepoint(1) < 0 ||
            pixdepoint(1) >=ny_p) {
        //cout << "in pix de point off " << pixdepoint << endl;
        return 2;
    }
    String pointingid=String::toString(pixdepoint(0))+"_"+String::toString(pixdepoint(1));
    //Int fieldid=vb.fieldId();
    String msid=vb.msName(true);
    //If channel or pol length has changed underneath...then its time to
    //restart the map
    /*
    if(convFunctionMap_p.ndefined() > 0){
      if ((fluxScale_p.shape()[3] != nchan_p) || (fluxScale_p.shape()[2] != npol_p)){
    convFunctionMap_p.clear();
      }
    }

    */
    if(convFunctionMap_p.nelements() > 0) {
        if (calcFluxScale_p && ((fluxScale_p.shape()[3] != nchan_p) || (fluxScale_p.shape()[2] != npol_p))) {
            convFunctionMap_p.resize();
            nDefined_p=0;
        }
    }
    //String mapid=msid+String("_")+pointingid;
    /*
    if(convFunctionMap_p.ndefined() == 0){
      convFunctionMap_p.define(mapid, 0);
      actualConvIndex_p=0;
      fluxScale_p=TempImage<Float>(IPosition(4,nx_p,ny_p,npol_p,nchan_p), csys_p);
      filledFluxScale_p=false;
      fluxScale_p.set(0.0);
      return -1;
    }
    */
    if(convFunctionMap_p.nelements() == 0) {
        convFunctionMap_p.resize(nx_p*ny_p);
        convFunctionMap_p.set(-1);
        convFunctionMap_p[pixdepoint[1]*nx_p+pixdepoint[0]]=0;
        nDefined_p=1;
        actualConvIndex_p=0;
        if(calcFluxScale_p) {
            fluxScale_p=TempImage<Float>(IPosition(4,nx_p,ny_p,npol_p,nchan_p), csys_p);
            filledFluxScale_p=false;
            fluxScale_p.set(0.0);
        }
        return -1;
    }

    // if(!convFunctionMap_p.isDefined(mapid)){
    //  actualConvIndex_p=convFunctionMap_p.ndefined();
    //  convFunctionMap_p.define(mapid, actualConvIndex_p);
    if(convFunctionMap_p[pixdepoint[1]*nx_p+pixdepoint[0]] <0) {
        actualConvIndex_p=nDefined_p;
        convFunctionMap_p[pixdepoint[1]*nx_p+pixdepoint[0]]=nDefined_p;
        // ++nDefined_p;
        nDefined_p=1;
        return -1;
    }
    else {
        /*
        actualConvIndex_p=convFunctionMap_p[pixdepoint[1]*nx_p+pixdepoint[0]];
        convFunc_p.resize(); // break any reference
        weightConvFunc_p.resize();
        convSupport_p.resize();
        //Here we will need to use the right xyPlane for different PA range
        //and frequency may be
        convFunc_p.reference(*convFunctions_p[actualConvIndex_p]);
        weightConvFunc_p.reference(*convWeights_p[actualConvIndex_p]);
        //Again this for one time of antenna only later should be fixed for all
        // antennas independently
        //these are not really needed right now
        convSupport_p=(*convSupportBlock_p[actualConvIndex_p]);
        convSize_p=(*convSizes_p[actualConvIndex_p])[0];
        makerowmap(vb, rowMap);
        */
        actualConvIndex_p=0;
        return -1;
    }

    return 1;


}

void HetArrayConvFunc::makerowmap(const vi::VisBuffer2& vb,
                                  Vector<Int>& rowMap) {

    uInt ndish=antMath_p.nelements();
    rowMap.resize(vb.nRows());
    for (Int k=0; k < vb.nRows(); ++k) {
        Int index1=antIndexToDiamIndex_p(vb.antenna1()(k));
        Int index2=antIndexToDiamIndex_p(vb.antenna2()(k));
        if(index2 < index1) {
            index1=index2;
            index2=antIndexToDiamIndex_p(vb.antenna1()(k));
        }
        Int plane=0;
        for (Int jj=0; jj < index1; ++jj)
            plane=plane+ndish-jj-1;
        plane=plane+index2;
        //plane of convfunc that match this pair of antennas
        rowMap(k)=plane;

    }

}


  void HetArrayConvFunc::multiplySelfConjugate(ImageInterface<Complex>& im, const Double freq){
    Double centerFreq=SpectralImageUtil::worldFreq(csys_p, 0.0);
    Double conjFreq=SynthesisUtils::conjFreq(freq, centerFreq);
    Double factor=freq/conjFreq;
    Array<Complex> arr;
    arr=im.get();
    IPosition shp=arr.shape();
    Array<Complex>arr2=resample(arr, factor);
    if(factor > 1) {
            IPosition end=shp-1;
            IPosition beg(5,0,0,0,0,0);
            beg(0)=(arr2.shape()(0)-shp(0))/2;
            beg(1)=(arr2.shape()(1)-shp(1))/2;
	    end(0)=beg(0)+shp(0)-1;
	    end(1)=beg(1)+shp(1)-1;
	    if(shp.nelements() >3)
	      end[3]=0;
            arr *= arr2(beg, end);
        }
        else {
	  IPosition end=arr2.shape()-1;
	  if(shp.nelements() >3)
	    end[3]=0;
	  IPosition beg(shp.nelements(),0);
	  beg(0)=(shp(0)-arr2.shape()(0))/2;
	  beg(1)=(shp(1)-arr2.shape()(1))/2;
	  end(0)+=beg(0);
	  end(1)+=beg(1);
	  Array<Complex> out(shp);
	  out.set(0);
	  out(beg, end)=arr2;
	  arr*=out;
        }
    im.put(arr);

  }

Array<Complex> HetArrayConvFunc::resample(const Array<Complex>& inarray, const Double factor) {

  if(factor==1.0)
    return inarray;
    Double nx=Double(inarray.shape()(0));
    Double ny=Double(inarray.shape()(1));
    IPosition shp=inarray.shape();
    shp(0)=Int(nx*factor/2.0)*2;
    shp(1)=Int(ny*factor/2.0)*2;
    Int newNx=shp(0);
    Int newNy=shp(1);
    Double xinvfactor=Double(nx)/Double(newNx);
    Double yinvfactor=Double(ny)/Double(newNy);
    Array<Complex> out(shp, 0.0);
   // cerr << "SHP " << shp << endl;
    
   IPosition incursor=IPosition(inarray.shape().nelements(),1);
    incursor[0]=nx;
    incursor[1]=ny;
    IPosition outcursor=IPosition(inarray.shape().nelements(),1);
    outcursor[0]=shp[0];
    outcursor[1]=shp[1];
    ArrayIterator<Complex> inIt(inarray, IPosition(2,0,1), True);
    ArrayIterator<Complex> outIt(out, IPosition(2,0,1),True);
    inIt.origin();
    outIt.origin();
    //for (zzz=0; zzz< shp.(4); ++zzz){
    //  for(yyy=0; yyy< shp.(3); ++yyy){
    // for(xxx=0; xxx< shp.(2); ++xxx){
    while(!inIt.pastEnd()) {
       // cerr << "Iter shape " << inIt.array().shape() << endl;
        Matrix<Complex> inmat;
        inmat=inIt.array();    
        //Matrix<Float> leReal=real(Matrix<Complex>(inIt.array()));
        //Matrix<Float> leImag=imag(Matrix<Complex>(inIt.array()));
        Matrix<Float> leReal=real(inmat);
        Matrix<Float> leImag=imag(inmat);
        Bool leRealCopy, leImagCopy;
        Float *realptr=leReal.getStorage(leRealCopy);
        Float *imagptr=leImag.getStorage(leImagCopy);
        Bool isCopy;
        Matrix<Complex> outMat(outIt.array());
        Complex *intPtr=outMat.getStorage(isCopy);
        Float realval, imagval;
#ifdef _OPENMP
	omp_set_nested(0);
#endif
#pragma omp parallel for default(none) private(realval, imagval) firstprivate(intPtr, realptr, imagptr, nx, ny, newNx, newNy, xinvfactor, yinvfactor) 

        for (Int k =0; k < newNy; ++k) {
            Double y =Double(k)*yinvfactor;

            for (Int j=0; j < newNx; ++j) {
                //      Interpolate2D interp(Interpolate2D::LANCZOS);
                Double x=Double(j)*xinvfactor;
                //interp.interp(realval, where, leReal);
                realval=interpLanczos(x , y, nx, ny,
                                      realptr);
                imagval=interpLanczos(x , y, nx, ny,
                                      imagptr);
                //interp.interp(imagval, where, leImag);
                intPtr[k*Int(newNx)+j]=Complex(realval, imagval);
            }

        }
        outMat.putStorage(intPtr, isCopy);
        leReal.putStorage(realptr, leRealCopy);
        leImag.putStorage(imagptr, leImagCopy);
        inIt.next();
        outIt.next();
    }
    return out;
}
Matrix <Complex> HetArrayConvFunc::resample2(const Matrix<Complex>& inarray, const Double factor) {

    Double nx=Double(inarray.shape()(0));
    Double ny=Double(inarray.shape()(1));
    IPosition shp=inarray.shape();
    shp(0)=Int(nx*factor/2.0)*2;
    shp(1)=Int(ny*factor/2.0)*2;

    
    Matrix<Complex> outMat(shp, Complex(0.0));
    
    
   
     {
        //cerr << "Iter shape " << inarray.shape() << endl;
        
        Matrix<Float> leReal=real(inarray);
        Matrix<Float> leImag=imag(inarray);
        Bool leRealCopy, leImagCopy;
        Float *realptr=leReal.getStorage(leRealCopy);
        Float *imagptr=leImag.getStorage(leImagCopy);
        Bool isCopy;
        Complex *intPtr=outMat.getStorage(isCopy);
        Float realval, imagval;
#ifdef _OPENMP
//        omp_set_nested(0);
#endif
 //       #pragma omp parallel for default(none) private(realval, imagval) firstprivate(intPtr, realptr, imagptr, nx, ny) shared(leReal, leImag)

        for (Int k =0; k < shp(1); ++k) {
            Double y =Double(k)/Double(shp(1))*Double(ny);

            for (Int j=0; j < Int(nx*factor); ++j) {
                //      Interpolate2D interp(Interpolate2D::LANCZOS);
                Double x=Double(j)/Double(factor);
                //interp.interp(realval, where, leReal);
                realval=interpLanczos(x , y, nx, ny,
                                      realptr);
                imagval=interpLanczos(x , y, nx, ny,
                                      imagptr);
                //interp.interp(imagval, where, leImag);
                intPtr[k*Int(nx*factor)+j]=Complex(realval, imagval);
            }

        }
        outMat.putStorage(intPtr, isCopy);
        leReal.putStorage(realptr, leRealCopy);
        leImag.putStorage(imagptr, leImagCopy);
        
     
    }
    return outMat;
}
Float HetArrayConvFunc::sinc(const Float x)  {
  Int index=x*1000+4000;
 
    
  return sincCachePtr_p[index];

  

  
}
  void HetArrayConvFunc::initSincCache(){
    for (Float u=-4000; u<4000; ++u){ 
      Float ux=u/1000.0;
      if (ux == 0) {
        sincCache_p[u+4000]=1.0;
      }
      else{
	sincCache_p[u+4000]= sin(C::pi * ux) / (C::pi * ux);
      }
    }
    sincCachePtr_p=sincCache_p.data();
  }
  
Float HetArrayConvFunc::interpLanczos( const Double& x , const Double& y, const Double& nx, const Double& ny,   const Float* data, const Float a) {
    Double floorx = floor(x);
    Double floory = floor(y);
    Float result=0.0;
    if (floorx < a || floorx >= nx - a || floory < a || floory >= ny - a) {
        result = 0;
        return result;
    }
    for (Float i = floorx - a + 1; i <= floorx + a; ++i) {
        for (Float j = floory - a + 1; j <= floory + a; ++j) {
            result += Float(Double(data[Int(j*nx+i)]) * sinc(x - i)*sinc((x-i)/ a) * sinc(y - j)*sinc((y-j)/ a));
        }
    }
    return result;
}

ImageInterface<Float>&  HetArrayConvFunc::getFluxScaleImage() {
    if(!calcFluxScale_p)
        throw(AipsError("Programmer Error: flux image cannot be retrieved"));
    if(!filledFluxScale_p) {
        //The best flux image for a heterogenous array is the weighted coverage
        fluxScale_p.copyData(*(convWeightImage_p));
        IPosition blc(4,nx_p, ny_p, npol_p, nchan_p);
        IPosition trc(4, ny_p, ny_p, npol_p, nchan_p);
        blc(0)=0;
        blc(1)=0;
        trc(0)=nx_p-1;
        trc(1)=ny_p-1;

        for (Int j=0; j < npol_p; ++j) {
            for (Int k=0; k < nchan_p ; ++k) {

                blc(2)=j;
                trc(2)=j;
                blc(3)=k;
                trc(3)=k;
                Slicer sl(blc, trc, Slicer::endIsLast);
                SubImage<Float> fscalesub(fluxScale_p, sl, true);
                Float planeMax;
                LatticeExprNode LEN = max( fscalesub );
                planeMax =  LEN.getFloat();
                if(planeMax !=0) {
                    fscalesub.copyData( (LatticeExpr<Float>) (fscalesub/planeMax));

                }
            }
        }
        filledFluxScale_p=true;
    }


    return fluxScale_p;

}

void HetArrayConvFunc::sliceFluxScale(Int npol) {
    IPosition fshp=fluxScale_p.shape();
    if (fshp(2)>npol) {
        npol_p=npol;
        // use first npol planes...
        IPosition blc(4,0,0,0,0);
        IPosition trc(4,fluxScale_p.shape()(0)-1, fluxScale_p.shape()(1)-1,npol-1,fluxScale_p.shape()(3)-1);
        Slicer sl=Slicer(blc, trc, Slicer::endIsLast);
        //writeable if possible
        SubImage<Float> fluxScaleSub = SubImage<Float> (fluxScale_p, sl, true);
        SubImage<Float> convWeightImageSub = SubImage<Float> (*convWeightImage_p, sl, true);
        fluxScale_p = TempImage<Float>(fluxScaleSub.shape(),fluxScaleSub.coordinates());
        convWeightImage_p = new TempImage<Float> (convWeightImageSub.shape(),convWeightImageSub.coordinates());
        LatticeExpr<Float> le(fluxScaleSub);
        fluxScale_p.copyData(le);
        LatticeExpr<Float> le2(convWeightImageSub);
        convWeightImage_p->copyData(le2);
    }
}
} // namespace refim end
} //# NAMESPACE CASA - END




