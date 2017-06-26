//# SolveDataBuffer.cc: Implementation of SolveDataBuffer.h
//# Copyright (C) 2008
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

#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisBufferComponents2.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayPartMath.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/iomanip.h>
#include <casa/Exceptions/Error.h>
#include <casa/Containers/Block.h>
#include <casa/Utilities/Assert.h>
#include <casa/Quanta/MVTime.h>
#include <casa/BasicSL/Constants.h>

using namespace casacore;
using namespace casa::vi;


using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

SolveDataBuffer::SolveDataBuffer() : 
  vb_(0),
  freqs_(0),
  feedPa_(0),
  focusChan_p(-1),
  infocusFlagCube_p(),
  infocusWtSpec_p(),
  infocusVisCube_p(),
  infocusModelVisCube_p(),
  residuals_p(),
  residFlagCube_p(),
  diffResiduals_p()
{}

SolveDataBuffer::SolveDataBuffer(const vi::VisBuffer2& vb) :
  vb_(0),
  freqs_(0),
  feedPa_(0),
  focusChan_p(-1),
  infocusFlagCube_p(),
  infocusWtSpec_p(),
  infocusVisCube_p(),
  infocusModelVisCube_p(),
  residuals_p(),
  residFlagCube_p(),
  diffResiduals_p()

{
  
  vb_= vi::VisBuffer2::factory(VbPlain,VbRekeyable);

  initFromVB(vb);
}

SolveDataBuffer::SolveDataBuffer(const SolveDataBuffer& sdb) :
  vb_(),
  freqs_(0),
  feedPa_(0),
  focusChan_p(-1),
  infocusFlagCube_p(),
  infocusWtSpec_p(),
  infocusVisCube_p(),
  infocusModelVisCube_p(),
  residuals_p(),
  residFlagCube_p(),
  diffResiduals_p()
{
  // Copy from the other's VB2
  initFromVB(*(sdb.vb_));

  // copy over freqs_, feedPa
  freqs_.assign(sdb.freqs_);
  feedPa_.assign(sdb.feedPa_);

}

SolveDataBuffer::~SolveDataBuffer()
{
  if (vb_) delete vb_;
  vb_=0;
 
}

//SolveDataBuffer& SolveDataBuffer::operator=(const VisBuffer& other)

  
Bool SolveDataBuffer::Ok() {
  // Ok if net unflagged weight is positive
  if (nfalse(flagCube())>0) {
    Float wtsum=sum(weightSpectrum()(!flagCube()));
    return (wtsum>0.0f);
  }
  else
    return false;
}

/*
// Divide corrected by model
void SolveDataBuffer::divideCorrByModel() {

  Int nCor(nCorrelations());
  Int nChan(nChannels());
  Int nRow(nRows());
  Float amp(1.0);
  Complex cor(1.0);
    
  Cube<Complex> vC; vC.reference(visCubeCorrected());
  Cube<Complex> vM; vM.reference(visCubeModel());
  Cube<Float> wS; wS.reference(weightSpectrum());

  for (Int irow=0;irow<nRow;++irow) {
    if (!flagRow()(irow)) {
      for (Int ich=0;ich<nChan;++ich) {
	for (Int icorr=0;icorr<nCor;icorr++) {
	  if (!flagCube()(icorr,ich,irow)) {
	    amp=abs(vM(icorr,ich,irow));
	    if (amp>0.0f) {
	      // Divide corr by model
	      vC(icorr,ich,irow)/=vM(icorr,ich,irow);
	      // Adjust weight by square of model amp
	      wS(icorr,ich,irow)*=square(amp);
	    }
	  } // !*fl
	  else {
	    // zero data and weight
	    vC(icorr,ich,irow)=Complex(0.0);
	    wS(icorr,ich,irow)=0.0;
	  }
	  // model always unity after division
	  vM(icorr,ich,irow)=Complex(1.0);
	} // icorr
      } // ich
    } // !*flR
  } // irow
}
*/

void SolveDataBuffer::enforceAPonData(const String& apmode)
{

  // TBD:  apply to model, too?  (never matters, since we prob /modelVis already?

  // ONLY if something to do
  if (apmode=="A" || apmode=="P") {
    Int nCor(nCorrelations());
    Int nChan(nChannels());
    Int nRow(nRows());
    Float amp(1.0);
    Complex cor(1.0);
    
    Cube<Complex> vC; vC.reference(visCubeCorrected());
    Cube<Float> wS; wS.reference(weightSpectrum());

    for (Int irow=0;irow<nRow;++irow) {
      if (!flagRow()(irow)) {
	for (Int ich=0;ich<nChan;++ich) {
	  for (Int icorr=0;icorr<nCor;icorr++) {
	    if (!flagCube()(icorr,ich,irow)) {
	      
	      amp=abs(vC(icorr,ich,irow));
	      if (amp>0.0f) {
		
		if (apmode=="P") {
		  // we will scale by amp to make data phase-only
		  cor=Complex(amp,0.0);
		  // Adjust weight by square of amp
		  wS(icorr,ich,irow)*=square(amp);
		}
		else if (apmode=="A")
		  // we will scale by "phase" to make data amp-only
		  //  no weight adjustment
		  cor=vC(icorr,ich,irow)/amp;
		
		// Apply the complex scaling
		vC(icorr,ich,irow)/=cor;
	      }
	    } // !*fl
	    else {
	      vC(icorr,ich,irow)=Complex(0.0);
	      wS(icorr,ich,irow)=0.0;
	    }
	  } // icorr
	} // ich
      } // !*flR
    } // irow

  } // phase- or amp-only

}

void SolveDataBuffer::enforceSolveWeights(const Bool phandonly)
{
  // If requested set cross-hand weights to zero (also---or only---flags?)
  if(phandonly && this->nCorrelations()>2)
      this->weightSpectrum()(Slice(1, 2, 1), Slice(), Slice()).set(0.0);

  // Set flagged weights to zero, ensuring they don't get used in accumulations
  this->weightSpectrum()(this->flagCube())=0.0f;

  // Set flagged data to zero (some solve types don't look at flags)
  Cube<Complex> vCC(this->visCubeCorrected());
  vCC(this->flagCube())=Complex(0.0f);

}



void SolveDataBuffer::setFocusChan(const Int focusChan) 
{

  // Nominally focus on the whole data array
  IPosition focusblc(3,0,0,0);
  IPosition focustrc(vb_->getShape());
  focustrc-=1;
  
  // if focusChan non-negative, select the single channel
  if (focusChan>-1) 
    focusblc(1)=focustrc(1)=focusChan;
  
  infocusFlagCube_p.reference(flagCube()(focusblc,focustrc));
  infocusWtSpec_p.reference(weightSpectrum()(focusblc,focustrc));
  infocusVisCube_p.reference(visCubeCorrected()(focusblc,focustrc));
  infocusModelVisCube_p.reference(visCubeModel()(focusblc,focustrc));

  // Remember current in-focus channel
  focusChan_p=focusChan;

}

void SolveDataBuffer::sizeResiduals(const Int& nPar,
				    const Int& nDiff)
{

  IPosition ip1(vb_->getShape());
  if (focusChan_p>-1)
    ip1(1)=1;
  residuals_p.resize(ip1);
  residuals_p.set(0.0);
  residFlagCube_p.resize(ip1);
  residFlagCube_p.set(false);

  if (nPar>0 && nDiff>0) {
    IPosition ip2(5,ip1(0),nPar,ip1(1),ip1(2),nDiff);
    diffResiduals_p.resize(ip2);
    diffResiduals_p.set(0.0);
  }

}

void SolveDataBuffer::initResidWithModel() 
{

  // Copy (literally) the in-focus model to the residual workspace
  // TBD:  weights?
  residuals_p = infocusModelVisCube_p;
  residFlagCube_p = infocusFlagCube_p;
    
  // Ensure contiguity, because CalSolver will depend on this
  AlwaysAssert(residFlagCube_p.contiguousStorage(),AipsError);
  AlwaysAssert(residuals_p.contiguousStorage(),AipsError);


}

void SolveDataBuffer::finalizeResiduals() 
{

  // Subtract in-focus obs data from residuals workspace
  residuals_p -= infocusVisCube_p;

  // TBD: zero flagged samples here?

}

void SolveDataBuffer::reportData()
{

  Slice corrs(0,2,3), sl;  // p-hands only 

  Vector<Int> a1(this->antenna1()), a2(this->antenna2());
  Cube<Complex> vCC(this->visCubeCorrected()(corrs,sl,sl));
  Cube<Complex> vCM(this->visCubeModel()(corrs,sl,sl));
  Cube<Float> wtS(this->weightSpectrum()(corrs,sl,sl));
  Cube<Bool> flg(this->flagCube()(corrs,sl,sl));

  /*
  cout << "******Making weights globally uniform!!!!!!!!!!!!!!!!" << endl;
    // this makes VI2 match VI1
  wtS.set(1.0f);
  //*/

  /*
  cout << "*******Zeroing zero-wt data (for K solutions)!!!!!!!!!!!!!" << endl;
    // this makes VI2 match VI1 solint='int' when some channels flagged
  vCC(wtS==0.0f)=Complex(0.0);
  //*/


  /*
  cout << "*******making weights uniform (for K solution testing)!!!!!!!!!!!!!" << endl;
  wtS.set(max(wtS));
  //*/

  /*
  cout << "*******zero off-wt data (for K solution testing)!!!!!!!!!!!!!" << endl;
  vCC(wtS!=max(wtS))=0.0f;
  wtS(wtS!=max(wtS))=0.0f;
  //*/

  /*
    // this makes VI2 match VI1 for  solint='inf' (and > 'int')
    cout << "*******renormalizing off-wt data (for K solution testing)!!!!!!!!!!!!!" << endl;
  vCC*=wtS;
  vCC/=max(wtS);
  //*/


  cout << "Time=" << MVTime(time()(0)/C::day).string(MVTime::YMD,7) << endl;
  cout.precision(8);
  for (Int irow=0;irow<nRows();++irow) {
    for (Int ich=0;ich<nChannels();++ich) {
      cout << std::setw(2) << a1(irow) << "-" << std::setw(2) << a2(irow) << " ";
      if (nChannels()>1) cout << "ich=" << ich << " ";
      cout << "fl=" << flg(Slice(),Slice(ich),Slice(irow)).nonDegenerate() << " ";
      cout << "wt=" << wtS(Slice(),Slice(ich),Slice(irow)).nonDegenerate() << " ";
      cout << "cM=" << vCM(Slice(),Slice(ich),Slice(irow)).nonDegenerate() << " ";
      cout << "cC=" << vCC(Slice(),Slice(ich),Slice(irow)).nonDegenerate() << " ";
      cout << endl;
    }
  }
  cout << "*****************************************************************************" << endl;
}

void SolveDataBuffer::initFromVB(const vi::VisBuffer2& vb) 
{

  // The required VB2 components
  vi::VisBufferComponents2 comps =
    vi::VisBufferComponents2::these({VisBufferComponent2::ObservationId,
 	                            VisBufferComponent2::ArrayId,
				    VisBufferComponent2::Scan,
				    VisBufferComponent2::FieldId,
				    VisBufferComponent2::DataDescriptionIds,
				    VisBufferComponent2::SpectralWindows,
				    VisBufferComponent2::Antenna1,
				    VisBufferComponent2::Antenna2,
				    VisBufferComponent2::Time,
				    VisBufferComponent2::TimeCentroid,
				    VisBufferComponent2::NCorrelations,
				    VisBufferComponent2::NChannels,
				    VisBufferComponent2::NRows,
				    VisBufferComponent2::FlagRow,
				    VisBufferComponent2::FlagCube,
	  VisBufferComponent2::WeightSpectrum,
	  VisBufferComponent2::VisibilityCubeCorrected,	  
	  VisBufferComponent2::VisibilityCubeModel});

  // Copy required components from the supplied VB2:
  vb_->copyComponents(vb,comps,True,True);   // will fetch things, if needed

  // Set weights for flagged data to zero
  Cube<Float> wtsp(vb_->weightSpectrum());
  wtsp(vb_->flagCube())=0.0;

  // Store the frequeny info
  //  TBD: also need bandwidth info....
  if (vb.isAttached())
    freqs_.assign(vb.getFrequencies(0));
  else {
    // Probably only needed in testing....  (gmoellen, 2016Aug04)
    cout << "The supplied VisBuffer2 is not attached to a ViImplementation2," << endl
	 << " which is necessary to generate accurate frequency info." << endl
	 << " This is probably just a test with a naked VisBuffer2." << endl
	 << " Spoofing freq axis with 1 MHz channels at 100 GHz." << endl;
    freqs_.resize(vb.nChannels());
    indgen(freqs_);
    freqs_*=1e6;
    freqs_+=100.0005e9; // _edge_ of first channel at 100 GHz.
  }

  // Store the feedPa info
  if (vb.isAttached())
    // Assumes vb.time() is constant!
    feedPa_.assign(vb.feedPa(vb.time()(0)));

}
void SolveDataBuffer::cleanUp() 
{

  // Zero-size all workspaces
  infocusFlagCube_p.resize();
  infocusWtSpec_p.resize();
  infocusVisCube_p.resize();
  infocusModelVisCube_p.resize();

  residuals_p.resize();
  residFlagCube_p.resize();
  diffResiduals_p.resize();

}


SDBList::SDBList() :
  nSDB_(0),
  SDB_(),
  freqs_()
{}

SDBList::~SDBList() 
{
  // Delete the SDBs
  for (Int i=0;i<nSDB_;++i)
    delete SDB_[i];
  SDB_.resize(0);
  nSDB_=0;
}

void SDBList::add(const vi::VisBuffer2& vb)
{

  // Enlarge the list, copying existing SDB pointers
  SDB_.resize(nSDB_+1,true);

  // Generate the new SolveDataBuffer
  SDB_[nSDB_] = new SolveDataBuffer(vb);

  // increment the count
  nSDB_++;

  // Clear the freqs_ info (forces recalculation)
  freqs_.resize(0);

}

SolveDataBuffer& SDBList::operator()(Int i) 
{
  if (i<nSDB_)
    return *SDB_[i];
  else
    throw(AipsError("SDBList::operator(): requests non-existent SolveDataBuffer."));

}

Int SDBList::aggregateObsId() const {
  if (nSDB_>0)
    // Obs Id from first SDB
    return SDB_[0]->observationId()(0);
  throw(AipsError("SDBList::aggregateObsId(): No SDBs in this SDBList yet."));
}

Int SDBList::aggregateScan() const {
  if (nSDB_>0)
    // Scan number from first SDB
    return SDB_[0]->scan()(0);
  throw(AipsError("SDBList::aggregateScan(): No SDBs in this SDBList yet."));
}

Int SDBList::aggregateSpw() const {
  if (nSDB_>0)
    // from first SDB
    return SDB_[0]->spectralWindow()(0);
  throw(AipsError("SDBList::aggregateSpw(): No SDBs in this SDBList yet."));
}

Int SDBList::aggregateFld() const {
  if (nSDB_>0)
    // from first SDB
    return SDB_[0]->fieldId()(0);
  throw(AipsError("SDBList::aggregateFld(): No SDBs in this SDBList yet."));
}

Double SDBList::aggregateTime() const {

  // Simple average of the mean times in each SDB
  //  (TBD: Improve with attention to flags/weights?)
  if (nSDB_>0) {
    Double aTime(0.0);
    for (Int isdb=0;isdb<nSDB_;++isdb)
      aTime+=mean(SDB_[isdb]->time());
    aTime/=nSDB_;
    return aTime;
  }
  else
    throw(AipsError("SDBList::aggregateTime(): No SDBs in this SDBList yet."));
}

Double SDBList::aggregateTimeCentroid() const {

  // Weighted average of SDBs' timeCentroids
  if (nSDB_>0) {
    Double aTime(0.0);
    Double aWt(0.0);
    Vector<Double> wtvD;
    for (Int isdb=0;isdb<nSDB_;++isdb) {
      SolveDataBuffer& sdb(*SDB_[isdb]);
      Vector<Float> wtv(partialSums(sdb.weightSpectrum(),IPosition(2,0,1)));
      wtvD.resize(wtv.nelements());
      convertArray(wtvD,wtv);
      aTime+=sum(wtvD*sdb.timeCentroid());
      aWt+=sum(wtvD);
    }
    if (aWt>0.0)
      aTime/=aWt;
    else
      // Use aggregateTime if no unflagged data
      aTime=aggregateTime();

    return aTime;
  }
  else
    throw(AipsError("SDBList::aggregateFld(): No SDBs in this SDBList yet."));
}


// How many data chans?
//   Currently, this insists on uniformity over all SDBs
//   In future, we may _sum_ the SDBs nChans, and
//    enable forming aggregate spectra (e.g., for common normalization)
//    This will require focusChan loop over SDBs...
Int SDBList::nChannels() const {

  Int nChan=SDB_[0]->nChannels();  // from first

  // Trap non-uniformity, for now
  for (Int isdb=1;isdb<nSDB_;++isdb)
    AlwaysAssert((SDB_[isdb]->nChannels()==nChan),AipsError);

  // Reach here, then ok
  return nChan;

}

const Vector<Double>& SDBList::freqs() const {

  if (nSDB_==0)
    throw(AipsError("SDBList::freqs(): No SDBs in this SDBList yet."));
  

  if (nSDB_==1) {
    // Only one SDB, just return that one's freqs
    const Vector<Double>& f(SDB_[0]->freqs());  // from first SDB
    return f;
  }

  // Reach here, more than one SDB, need to gather info
  if (freqs_.nelements()==0) {
    // Haven't accumumlated yet
    
    // How many channels in aggregate?
    //  (This will insist on uniformity over SDBs)
    Int nchan(this->nChannels());   
    
    // Will accumulate mean freqs here
    freqs_.resize(nchan);
    freqs_.set(0.0);
    
    // Average over SDBs, counting each spw exactly _once_
    // Map to keep track of unique spws
    std::set<Int> spws;
    Int nSpw(0); 
    for (Int isdb=0;isdb<nSDB_;++isdb) {
      Int ispw=SDB_[isdb]->spectralWindow()(0);
      if (spws.count(ispw)<1) {
	freqs_+=SDB_[isdb]->freqs();
	spws.insert(ispw);  // Record that we got this spw
	++nSpw;
      }
    }
    // Divide by nSpw
    freqs_/=Double(nSpw);
    
  }
  // else:  freqs_ already filled previously...

  return freqs_;
  
}

// How many correlations?
//   Currently, this insists on uniformity over all SDBs
Int SDBList::nCorrelations() const {

  Int nCorr=SDB_[0]->nCorrelations();

  // Trap non-uniformity, for now
  for (Int isdb=1;isdb<nSDB_;++isdb)
    AlwaysAssert((SDB_[isdb]->nCorrelations()==nCorr),AipsError);

  // Reach here, then ok
  return nCorr;

}





Bool SDBList::Ok() {

  for (Int i=0;i<nSDB_;++i)
    if (SDB_[i]->Ok()) return True;

  // If we get here, either no SDBs, or none have non-zero weight.
  return False;

}

/*
void SDBList::divideCorrByModel()
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->divideCorrByModel();
}
*/

void SDBList::enforceAPonData(const String& apmode)
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->enforceAPonData(apmode);
}
void SDBList::enforceSolveWeights(const Bool pHandOnly) 
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->enforceSolveWeights(pHandOnly);
}

void SDBList::sizeResiduals(const Int& nPar, const Int& nDiff)
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->sizeResiduals(nPar,nDiff);
}
void SDBList::initResidWithModel()
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->initResidWithModel();
}
void SDBList::finalizeResiduals()
{
  for (Int i=0;i<nSDB_;++i)
    SDB_[i]->finalizeResiduals();
}

void SDBList::reportData()
{
  cout << "nSDB=" << nSDB_ << endl;
  for (Int i=0;i<nSDB_;++i) {
    cout << "isdb=" << i << endl;
    SDB_[i]->reportData();
  }
}




} //# NAMESPACE CASA - END

