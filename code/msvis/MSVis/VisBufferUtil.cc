//# VisBufferUtil.cc: VisBuffer Utilities
//# Copyright (C) 1996,1997,2001
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#include <casa/aips.h>

#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/OS/Timer.h>
#include <measures/Measures/UVWMachine.h>
#include <measures/Measures/MeasTable.h>
#include <ms/MeasurementSets/MSPointingColumns.h>
#include <msvis/MSVis/VisBufferUtil.h>
#include <msvis/MSVis/StokesVector.h>
#include <msvis/MSVis/VisibilityIterator.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/VisBuffer.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <casa/iostream.h>
#include <iomanip>
using namespace std;

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// <summary> 
// </summary>

// <reviewed reviewer="" date="" tests="tMEGI" demos="">

// <prerequisite>
// </prerequisite>
//
// <etymology>
// </etymology>
//
// <synopsis> 
// </synopsis> 
//
// <example>
// <srcblock>
// </srcblock>
// </example>
//
// <motivation>
// </motivation>
//
// <todo asof="">
// </todo>


VisBufferUtil::VisBufferUtil(): oldMSId_p(-1), timeAntIndex_p(0), cachedPointingDir_p(0){};


// Construct from a VisBuffer (sets frame info)
VisBufferUtil::VisBufferUtil(const VisBuffer& vb): oldMSId_p(-1), timeAntIndex_p(0), cachedPointingDir_p(0) {

  // The nominal epoch 
  MEpoch ep=vb.msColumns().timeMeas()(0);

  // The nominal position
  String observatory;
  MPosition pos;
  if (vb.msColumns().observation().nrow() > 0) {
    observatory = vb.msColumns().observation().telescopeName()
      (vb.msColumns().observationId()(0));
  }
  if (observatory.length() == 0 || 
      !MeasTable::Observatory(pos,observatory)) {
    // unknown observatory, use first antenna
    pos=vb.msColumns().antenna().positionMeas()(0);
  }
 
  // The nominal direction
  MDirection dir=vb.phaseCenter();

  // The nominal MeasFrame
  mframe_=MeasFrame(ep, pos, dir);

}

// Construct from a VisBuffer (sets frame info)
VisBufferUtil::VisBufferUtil(const vi::VisBuffer2& vb): oldMSId_p(-1) {
	if(vb.getVi() == NULL)
		ThrowCc("Programmer Error: used a detached Visbuffer when it should not have been so");
	ROMSColumns msc(vb.getVi()->ms());
  // The nominal epoch
  MEpoch ep=msc.timeMeas()(0);

  // The nominal position
  String observatory;
  MPosition pos;
  if (msc.observation().nrow() > 0) {
    observatory = msc.observation().telescopeName()
      (msc.observationId()(0));
  }
  if (observatory.length() == 0 ||
      !MeasTable::Observatory(pos,observatory)) {
    // unknown observatory, use first antenna
    pos=msc.antenna().positionMeas()(0);
  }

  // The nominal direction
  MDirection dir=vb.phaseCenter();

  // The nominal MeasFrame
  mframe_=MeasFrame(ep, pos, dir);

}
VisBufferUtil::VisBufferUtil(const vi::VisibilityIterator2& iter): oldMSId_p(-1) {

	ROMSColumns msc(iter.ms());
  // The nominal epoch
  MEpoch ep=msc.timeMeas()(0);

  // The nominal position
  String observatory;
  MPosition pos;
  if (msc.observation().nrow() > 0) {
    observatory = msc.observation().telescopeName()
      (msc.observationId()(0));
  }
  if (observatory.length() == 0 ||
      !MeasTable::Observatory(pos,observatory)) {
    // unknown observatory, use first antenna
    pos=msc.antenna().positionMeas()(0);
  }

  // The nominal direction
  //MDirection dir=iter.phaseCenter();
  MDirection dir=msc.field().phaseDirMeasCol()(0)(IPosition(1,0));
  // The nominal MeasFrame
  mframe_=MeasFrame(ep, pos, dir);

}
VisBufferUtil::VisBufferUtil(const MeasFrame& mframe): oldMSId_p(-1) {
	mframe_=mframe;

}
Bool VisBufferUtil::rotateUVW(const vi::VisBuffer2&vb, const MDirection& desiredDir,
				Matrix<Double>& uvw, Vector<Double>& dphase){

    Bool retval=true;
    mframe_.resetEpoch(vb.time()(0));
    UVWMachine uvwMachine(desiredDir, vb.phaseCenter(), mframe_,
			false, false);
    retval = !uvwMachine.isNOP();
    dphase.resize(vb.nRows());
    dphase.set(0.0);
    if(uvw.nelements() ==0)
      uvw=vb.uvw();
    for (Int row=0; row< vb.nRows(); ++row){
      Vector<Double> eluvw(uvw.column(row));
      uvwMachine.convertUVW(dphase(row), eluvw);
    }
    
    return retval;
  }

// Set the visibility buffer for a PSF
void VisBufferUtil::makePSFVisBuffer(VisBuffer& vb) {
  CStokesVector coh(Complex(1.0), Complex(0.0), Complex(0.0), Complex(1.0));
  vb.correctedVisibility()=coh;
}


Bool VisBufferUtil::interpolateFrequency(Cube<Complex>& data, 
					 Cube<Bool>& flag, 
					 const VisBuffer& vb,
					 const Vector<Float>& outFreqGrid, 	
					 const MS::PredefinedColumns whichCol, 
					 const MFrequency::Types freqFrame,
					 const InterpolateArray1D<Float,Complex>::InterpolationMethod interpMethod){

  Cube<Complex> origdata;
  // Convert the visibility frequency to the frame requested
  Vector<Double> visFreqD;
  convertFrequency(visFreqD, vb, freqFrame);
  //convert it to Float
  Vector<Float> visFreq(visFreqD.nelements());
  convertArray(visFreq, visFreqD);
  
  //Assign which column is to be regridded to origdata
  if(whichCol==MS::MODEL_DATA){
    origdata.reference(vb.modelVisCube());
  }
  else if(whichCol==MS::CORRECTED_DATA){
      origdata.reference(vb.correctedVisCube());
    }
  else if(whichCol==MS::DATA){
      origdata.reference(vb.visCube());
  }
  else{
    throw(AipsError("Don't know which column is being regridded"));
  }
  Cube<Complex> flipdata;
  Cube<Bool> flipflag;
  //The interpolator interpolates on the 3rd axis only...so need to flip the axes (y,z)
  swapyz(flipflag,vb.flagCube());
  swapyz(flipdata,origdata);

  //interpolate the data and the flag to the output frequency grid
  InterpolateArray1D<Float,Complex>::
    interpolate(data,flag, outFreqGrid,visFreq,flipdata,flipflag,interpMethod);
  flipdata.resize();
  //reflip the data and flag to be in the same order as in Visbuffer output
  swapyz(flipdata,data);
  data.resize();
  data.reference(flipdata);
  flipflag.resize();
  swapyz(flipflag,flag);
  flag.resize();     
  flag.reference(flipflag);

  return true;

}
  void VisBufferUtil::getFreqRange(Double& freqMin, Double& freqMax, vi::VisibilityIterator2& vi, MFrequency::Types freqFrame){
    vi.originChunks();
    vi.origin();

    Double freqEnd=0.0;
    Double freqStart=C::dbl_max;
    vi::VisBuffer2* vb=vi.getVisBuffer();
    for (vi.originChunks(); vi.moreChunks();vi.nextChunk())
    	{
	  for (vi.origin(); vi.more();vi.next())
    		{
		  Double localmax, localmin;
                  IPosition localmaxpos(1,0); 
                  IPosition localminpos(1,0);
		  Vector<Double> freqs=vb->getFrequencies(0, freqFrame);
		  if(freqs.nelements() ==0){
		    throw(AipsError("Frequency selection error" ));
		  }
                  minMax(localmin,localmax,localminpos, localmaxpos, freqs);
		  //localmax=max(freqs);
		  //localmin=min(freqs);
		  //freqEnd=max(freqEnd, localmax);
		  //freqStart=min(freqStart, localmin);

                  Int nfreq = freqs.nelements(); 
                  Vector<Int> curspws = vb->spectralWindows();
                  // as the vb row 0 is used for getFrequencies, the same row 0 is used here
                  Vector<Double> chanWidths = vi.subtableColumns().spectralWindow().chanWidth()(curspws[0]);  
                  // freqs are in channel center freq so add the half the width to the values to return the edge frequencies 
                  if (nfreq==1) {
		    freqEnd=max(freqEnd, localmax+fabs(chanWidths[0]/2.0));
		    freqStart=min(freqStart, localmin-fabs(chanWidths[0]/2.0));
                  }
                  else {
		    freqEnd=max(freqEnd, localmax+fabs(chanWidths[localmaxpos[0]]/2.0));
		    freqStart=min(freqStart, localmin-fabs(chanWidths[localminpos[0]]/2.0));
                  }
		   
		}
	}
    freqMin=freqStart;
    freqMax=freqEnd;
  }

  void VisBufferUtil::getFreqRangeFromRange(casacore::Double& outfreqMin, casacore::Double& outfreqMax,  const casacore::MFrequency::Types inFreqFrame, const casacore::Double infreqMin, const casacore::Double infreqMax, vi::VisibilityIterator2& vi, casacore::MFrequency::Types outFreqFrame){


    if(inFreqFrame==outFreqFrame){
      outfreqMin=infreqMin;
      outfreqMax=infreqMax;
      return;
    }

    vi.originChunks();
    vi.origin();

    outfreqMin=C::dbl_max;
    outfreqMax=0;
    vi::VisBuffer2* vb=vi.getVisBuffer();
    ROMSColumns msc(vi.ms());
  // The nominal epoch
    MEpoch ep=msc.timeMeas()(0);
    
    // The nominal position
    String observatory;
    MPosition pos;
    if (msc.observation().nrow() > 0) {
      observatory = msc.observation().telescopeName()
	(msc.observationId()(0));
    }
    if (observatory.length() == 0 ||
	!MeasTable::Observatory(pos,observatory)) {
      // unknown observatory, use first antenna
      pos=msc.antenna().positionMeas()(0);
    }

  // The nominal direction
  MDirection dir=vb->phaseCenter();
  MeasFrame mFrame(ep, pos, dir);
   // The conversion engine:
     MFrequency::Convert toNewFrame(inFreqFrame, 
				    MFrequency::Ref(outFreqFrame, mFrame));
  
      for (vi.originChunks(); vi.moreChunks();vi.nextChunk())
    	{
	  for (vi.origin(); vi.more();vi.next()){
	    //assuming time is fixed in visbuffer
	    mFrame.resetEpoch(vb->time()(0)/86400.0);

	    // Reset the direction (ASSUMES phaseCenter is constant in the VisBuffer)
	    mFrame.resetDirection(vb->phaseCenter());
	    Double temp=toNewFrame(infreqMin).getValue().getValue();
	    if(temp < outfreqMin)
	      outfreqMin = temp;
	    
	    temp=toNewFrame(infreqMax).getValue().getValue();
	    if(temp > outfreqMax)
	      outfreqMax = temp;	      
	  }
	}
      //cerr << "min " << outfreqMin << " max " << outfreqMax << endl;

  }

void VisBufferUtil::convertFrequency(Vector<Double>& outFreq, 
				     const VisBuffer& vb, 
				     const MFrequency::Types freqFrame){
   Int spw=vb.spectralWindow();
   MFrequency::Types obsMFreqType=(MFrequency::Types)(vb.msColumns().spectralWindow().measFreqRef()(spw));

   // The input frequencies (a reference)
   Vector<Double> inFreq(vb.frequency());

   // The output frequencies
   outFreq.resize(inFreq.nelements());

   MFrequency::Types newMFreqType=freqFrame;
   if (freqFrame==MFrequency::N_Types)
     // Opt out of conversion
     newMFreqType=obsMFreqType;


   // Only convert if the requested frame differs from observed frame
   if(obsMFreqType != newMFreqType){

     // Setting epoch to the first in this iteration
     //     MEpoch ep=vb.msColumns().timeMeas()(0);
     //     MEpoch ep(MVEpoch(vb.time()(0)/86400.0),MEpoch::UTC);
     //     cout << "Time = " << ep.getValue()  << endl;

     // Reset the timestamp (ASSUMES TIME is constant in the VisBuffer)
     mframe_.resetEpoch(vb.time()(0)/86400.0);

     // Reset the direction (ASSUMES phaseCenter is constant in the VisBuffer)
     mframe_.resetDirection(vb.msColumns().field().phaseDirMeasCol()(vb.fieldId())(IPosition(1,0)));

     //     cout << "Frame = " << mframe_ << endl;

     // The conversion engine:
     MFrequency::Convert toNewFrame(obsMFreqType, 
				    MFrequency::Ref(newMFreqType, mframe_));

     // Do the conversion
     for (uInt k=0; k< inFreq.nelements(); ++k)
       outFreq(k)=toNewFrame(inFreq(k)).getValue().getValue();
     
   }
   else{
     // The requested frame is the same as the observed frame
     outFreq=inFreq;
   }

 }

void VisBufferUtil::convertFrequency(Vector<Double>& outFreq, 
				     const vi::VisBuffer2& vb, 
				     const MFrequency::Types freqFrame){
  Int spw=vb.spectralWindows()(0);
  MFrequency::Types obsMFreqType=(MFrequency::Types)(ROMSColumns(vb.getVi()->ms()).spectralWindow().measFreqRef()(spw));

  

   // The input frequencies 
  Vector<Int> chanNums=vb.getChannelNumbers(0);
  
  Vector<Double> inFreq(chanNums.nelements());
  Vector<Double> spwfreqs=ROMSColumns(vb.getVi()->ms()).spectralWindow().chanFreq().get(spw);
  for (uInt k=0; k < chanNums.nelements(); ++k){

    inFreq[k]=spwfreqs[chanNums[k]];
  }

   // The output frequencies
   outFreq.resize(inFreq.nelements());

   MFrequency::Types newMFreqType=freqFrame;
   if (freqFrame==MFrequency::N_Types)
     // Opt out of conversion
     newMFreqType=obsMFreqType;


   // Only convert if the requested frame differs from observed frame
   if(obsMFreqType != newMFreqType){

     // Setting epoch to the first in this iteration
     //     MEpoch ep=vb.msColumns().timeMeas()(0);
     //     MEpoch ep(MVEpoch(vb.time()(0)/86400.0),MEpoch::UTC);
     //     cout << "Time = " << ep.getValue()  << endl;

     // Reset the timestamp (ASSUMES TIME is constant in the VisBuffer)
     mframe_.resetEpoch(vb.time()(0)/86400.0);

     // Reset the direction (ASSUMES phaseCenter is constant in the VisBuffer)
     mframe_.resetDirection(vb.phaseCenter());

     //     cout << "Frame = " << mframe_ << endl;

     // The conversion engine:
     MFrequency::Convert toNewFrame(obsMFreqType, 
				    MFrequency::Ref(newMFreqType, mframe_));

     // Do the conversion
     for (uInt k=0; k< inFreq.nelements(); ++k)
       outFreq(k)=toNewFrame(inFreq(k)).getValue().getValue();
     
   }
   else{
     // The requested frame is the same as the observed frame
     outFreq=inFreq;
   }

   //cerr << std::setprecision(9) << " infreq " << inFreq[152] << "   " << outFreq[152] << " vb freq " << vb.getFrequencies(0, freqFrame)[152] << endl;

 }

 void VisBufferUtil::toVelocity(Vector<Double>& outVel, 
				const VisBuffer& vb, 
				const MFrequency::Types freqFrame,
				const MVFrequency restFreq,
				const MDoppler::Types veldef){

   // The input frequencies (a reference)
   Vector<Double> inFreq(vb.frequency());

   // The output velocities
   outVel.resize(inFreq.nelements());

   // Reset the timestamp (ASSUMES TIME is constant in the VisBuffer)
   mframe_.resetEpoch(vb.time()(0)/86400.0);
   
   // Reset the direction (ASSUMES phaseCenter is constant in the VisBuffer)
   //mframe_.resetDirection(vb.phaseCenter());
   mframe_.resetDirection(vb.msColumns().field().phaseDirMeasCol()(vb.fieldId())(IPosition(1,0)));
 
   // The frequency conversion engine:
   Int spw=vb.spectralWindow();
   MFrequency::Types obsMFreqType=(MFrequency::Types)(vb.msColumns().spectralWindow().measFreqRef()(spw));

   MFrequency::Types newMFreqType=freqFrame;
   if (freqFrame==MFrequency::N_Types)
     // Don't convert frame
     newMFreqType=obsMFreqType;

   MFrequency::Convert toNewFrame(obsMFreqType, 
				  MFrequency::Ref(newMFreqType, mframe_));

   // The velocity conversion engine:
   MDoppler::Ref dum1(MDoppler::RELATIVISTIC);
   MDoppler::Ref dum2(veldef);
   MDoppler::Convert dopConv(dum1, dum2);

   // Cope with unspecified rest freq
   MVFrequency rf=restFreq;
   if (restFreq.getValue()<=0.0)
     rf=toNewFrame(inFreq(vb.nChannel()/2)).getValue();

   // Do the conversions
   for (uInt k=0; k< inFreq.nelements(); ++k){
	 //cerr <<"old freq " << toNewFrame(inFreq(k)).getValue().get().getValue() << endl;
     MDoppler eh = toNewFrame(inFreq(k)).toDoppler(rf);
     MDoppler eh2 = dopConv(eh);
     outVel(k)=eh2.getValue().get().getValue();
   }

 }

 void VisBufferUtil::toVelocity(Vector<Double>& outVel,
 				const vi::VisBuffer2& vb,
 				const MFrequency::Types freqFrame,
 				const MVFrequency restFreq,
 				const MDoppler::Types veldef, const Int row){

	 	 Vector<Double> inFreq;
	 	 inFreq=vb.getFrequencies(row, freqFrame);
	 	// cerr << "Freqs " << inFreq << endl;
	 	// The output velocities
	 	 outVel.resize(inFreq.nelements());
	 	 // The velocity conversion engine:
	 	 MDoppler::Ref dum1(MDoppler::RELATIVISTIC);
	 	 MDoppler::Ref dum2(veldef);
	 	 MDoppler::Convert dopConv(dum1, dum2);

	 	 // Cope with unspecified rest freq
	 	 MVFrequency rf=restFreq;
	 	 if (restFreq.getValue()<=0.0)
	 	      rf=inFreq(inFreq.nelements()/2);

	 	 // Do the conversions
	 	 for (uInt k=0; k< inFreq.nelements(); ++k){
	 		 MDoppler eh = MFrequency(Quantity(inFreq(k), "Hz"), freqFrame).toDoppler(rf);
	 		 MDoppler eh2 = dopConv(eh);
	 		 outVel(k)=eh2.getValue().get().getValue();
	 	 	}


 }
 void VisBufferUtil::toVelocity(Vector<Double>& outVel,
  				const vi::VisBuffer2& vb,
  				const vi::VisibilityIterator2& iter,
  				const MFrequency::Types freqFrame,
  				const MVFrequency restFreq,
  				const MDoppler::Types veldef, const Int row){

 	 // The input frequencies (a reference)
 	 Vector<Double> inFreq(vb.getFrequencies(row));
 	 ROMSColumns msc(iter.ms());

     MEpoch ep(Quantity(vb.time()(row)/86400.0, "d"), msc.timeMeas()(0).getRef());
 	 MDirection dir(msc.field().phaseDirMeasCol()(vb.fieldId()(row))(IPosition(1,0)));
 	 Int spw=vb.spectralWindows()(row);
 	 MFrequency::Types obsMFreqType=(MFrequency::Types)(msc.spectralWindow().measFreqRef()(spw));
 	 toVelocity(outVel, freqFrame, inFreq, obsMFreqType, ep, dir, restFreq, veldef);
  }
 void VisBufferUtil::toVelocity(Vector<Double>& outVel,
   		  const MFrequency::Types outFreqFrame,
   		  const Vector<Double>& inFreq,
   		  const MFrequency::Types inFreqFrame,
   		  const MEpoch& ep,
   		  const MDirection& dir,
   		  const MVFrequency restFreq,
   		  const MDoppler::Types veldef){



	 // The output velocities
	 outVel.resize(inFreq.nelements());

	 // Reset the timestamp
	 mframe_.resetEpoch(ep);

	 // Reset the direction
	 mframe_.resetDirection(dir);

	 // The frequency conversion engine:

	 MFrequency::Types newMFreqType=outFreqFrame;
	 if (outFreqFrame==MFrequency::N_Types)
		 // Don't convert frame
		 newMFreqType=inFreqFrame;

	 MFrequency::Convert toNewFrame(inFreqFrame,
	 				  MFrequency::Ref(newMFreqType, mframe_));

	 // The velocity conversion engine:
	 MDoppler::Ref dum1(MDoppler::RELATIVISTIC);
	 MDoppler::Ref dum2(veldef);
	 MDoppler::Convert dopConv(dum1, dum2);

	 // Cope with unspecified rest freq
	 MVFrequency rf=restFreq;
	 if (restFreq.getValue()<=0.0)
	      rf=toNewFrame(inFreq(inFreq.nelements()/2)).getValue();

	 // Do the conversions
	 for (uInt k=0; k< inFreq.nelements(); ++k){
		 MDoppler eh = toNewFrame(inFreq(k)).toDoppler(rf);
		 MDoppler eh2 = dopConv(eh);
		 outVel(k)=eh2.getValue().get().getValue();
	 	}




 }

 MDirection VisBufferUtil::getPointingDir(const VisBuffer& vb, const Int antid, const Int vbrow){
	 Timer tim;
	 tim.mark();
	 //MDirection outdir;
	 if(oldMSId_p != vb.msId()){
		 tim.mark();
		 oldMSId_p=vb.msId();
		 if(timeAntIndex_p.shape()(0) < (oldMSId_p+1)){
			 timeAntIndex_p.resize(oldMSId_p+1, true);
		 	 cachedPointingDir_p.resize(oldMSId_p+1, true);
		 }
		 if(  timeAntIndex_p[oldMSId_p].empty()){
			 Vector<Double> tOrig;
			 vb.msColumns().time().getColumn(tOrig);
			 Vector<Double> t;
			 rejectConsecutive(tOrig, t);
			 Vector<uInt>  uniqIndx;
			 uInt nTimes=GenSortIndirect<Double>::sort (uniqIndx, t, Sort::Ascending, Sort::QuickSort|Sort::NoDuplicates);
			 uInt nAnt=vb.msColumns().antenna().nrow();
			 const ROMSPointingColumns& mspc=vb.msColumns().pointing();
			 Int guessIndex=0;
			 for (uInt k=0; k <nTimes; ++k){
				 for (uInt a=0; a < nAnt; ++a){
					 std::ostringstream oss;
					 oss.precision(13);
					 oss << t[uniqIndx[k]] << "_" << a;
					 String key=oss.str();
					 //String key=String::toString(t[uniqIndx[k]])+String("_")+String::toString(a);
					 Int row=mspc.pointingIndex(a, t[uniqIndx[k]], guessIndex);
					 cerr << "String "<< key << "pointing row "<< row << endl;
					 timeAntIndex_p[oldMSId_p][key]=row > -1 ? cachedPointingDir_p[oldMSId_p].shape()[0] : -1;
					 guessIndex=row;
					 if(row >-1){
						 cachedPointingDir_p[oldMSId_p].resize(cachedPointingDir_p[oldMSId_p].nelements()+1, true);
						 cachedPointingDir_p[oldMSId_p][cachedPointingDir_p[oldMSId_p].nelements()-1]=mspc.directionMeas(row);
					 }

				 }
			 }

		 }
		 tim.show("After caching all ant pointings");
	 }

	 /////
	 //	 String index=String::toString(vb.time()(vbrow))+String("_")+String::toString(antid);
	 std::ostringstream oss;
	 oss.precision(13);
	 oss << vb.time()(vbrow) << "_" << antid  ;
	 String index=oss.str();
	 Int rowincache=timeAntIndex_p[oldMSId_p][index];
	 cerr << "key "<< index << " index " << rowincache << endl;
	 tim.show("retrieved cache");
	 if(rowincache <0)
		 return vb.phaseCenter();
	 return cachedPointingDir_p[oldMSId_p][rowincache];



 }
  MDirection VisBufferUtil::getPointingDir(const vi::VisBuffer2& vb, const Int antid, const Int vbrow){
	 Timer tim;
	 tim.mark();
	 ROMSColumns msc(vb.ms());
	 //MDirection outdir;
	 if(oldMSId_p != vb.msId()){
		 tim.mark();
		 cerr << "MSID: "<< oldMSId_p <<  "    " << vb.msId() << endl;
		 oldMSId_p=vb.msId();
		 if(timeAntIndex_p.shape()(0) < (oldMSId_p+1)){
			 timeAntIndex_p.resize(oldMSId_p+1, true);
		 	 cachedPointingDir_p.resize(oldMSId_p+1, true);
		 }
		 if(  timeAntIndex_p[oldMSId_p].empty()){
			 Vector<Double> tOrig;
			 msc.time().getColumn(tOrig);
			 Vector<Double> t;
			 rejectConsecutive(tOrig, t);
			 Vector<uInt>  uniqIndx;
			 uInt nTimes=GenSortIndirect<Double>::sort (uniqIndx, t, Sort::Ascending, Sort::QuickSort|Sort::NoDuplicates);
			 uInt nAnt=msc.antenna().nrow();
			 const ROMSPointingColumns& mspc=msc.pointing();
			 Int guessIndex=0;
			 for (uInt k=0; k <nTimes; ++k){
				 for (uInt a=0; a < nAnt; ++a){
					 std::ostringstream oss;
					 oss.precision(13);
					 oss << t[uniqIndx[k]] << "_" << a;
					 String key=oss.str();
					 //String key=String::toString(t[uniqIndx[k]])+String("_")+String::toString(a);
					 Int row=mspc.pointingIndex(a, t[uniqIndx[k]], guessIndex);
					 //cerr << "String "<< key << " pointing row "<< row << endl;
					 timeAntIndex_p[oldMSId_p][key]=row > -1 ? cachedPointingDir_p[oldMSId_p].shape()[0] : -1;
					 guessIndex=row;
					 if(row >-1){
						 cachedPointingDir_p[oldMSId_p].resize(cachedPointingDir_p[oldMSId_p].nelements()+1, true);
						 cachedPointingDir_p[oldMSId_p][cachedPointingDir_p[oldMSId_p].nelements()-1]=mspc.directionMeas(row);
					 }

				 }
			 }

		 }
		 tim.show("After caching all ant pointings");
	 }

	 /////
	 //	 String index=String::toString(vb.time()(vbrow))+String("_")+String::toString(antid);
	 std::ostringstream oss;
	 oss.precision(13);
	 oss << vb.time()(vbrow) << "_" << antid  ;
	 String index=oss.str();
	 Int rowincache=timeAntIndex_p[oldMSId_p][index];
	 //cerr << "key "<< index << " index " << rowincache << endl;
	 tim.show("retrieved cache");
	 if(rowincache <0)
		 return vb.phaseCenter();
	 return cachedPointingDir_p[oldMSId_p][rowincache];



 }
 //utility to reject consecutive similar value for sorting
 void VisBufferUtil::rejectConsecutive(const Vector<Double>& t, Vector<Double>& retval){
     uInt n=t.nelements();
     if(n >0){
       retval.resize(n);
       retval[0]=t[0];
     }
     else
       return;
     Int prev=0;
     for (uInt k=1; k < n; ++k){
       if(t[k] != retval(prev)){
    	   ++prev;

    	   retval[prev]=t[k];
       }
     }
     retval.resize(prev+1, true);

   }
// helper function to swap the y and z axes of a Cube
 void   VisBufferUtil::swapyz(Cube<Complex>& out, const Cube<Complex>& in)
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

// helper function to swap the y and z axes of a Cube
void VisBufferUtil::swapyz(Cube<Bool>& out, const Cube<Bool>& in)
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

 


} //# NAMESPACE CASA - END

