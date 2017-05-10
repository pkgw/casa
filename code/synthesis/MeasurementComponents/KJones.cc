//# KJones.cc: Implementation of KJones
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003,2011
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

#include <synthesis/MeasurementComponents/KJones.h>

#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisBuffAccumulator.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MSOper/MSMetaData.h>
#include <synthesis/MeasurementEquations/VisEquation.h>  // *
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>
#include <casa/Quanta/QuantumHolder.h>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/System/Aipsrc.h>

#include <casa/sstream.h>

#include <measures/Measures/MCBaseline.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasTable.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>


using namespace casa::vi;
using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN


// **********************************************************
// DelayFFT Implementations
//

// Construct from freq info and a data-like Cube<Complex>
DelayFFT::DelayFFT(Double f0, Double df, Double padBW, 
		   Cube<Complex> V) :
  f0_(f0),
  df_(df),
  padBW_(padBW),
  nCorr_(V.shape()(0)),
  nPadChan_(Int(0.5+padBW/df)),
  nElem_(V.shape()(2)),
  refant_(-1),  // ok?
  Vpad_(),
  delay_(),
  flag_()
{

  //  cout << "pad ch=" << padBW/df_ << " " << nPadChan_*df << endl;

  IPosition ip1(3,nCorr_,nPadChan_,nElem_);
  Vpad_.resize(ip1);
  Vpad_.set(0.0);

  Slicer sl1(Slice(),Slice(0,V.shape()(1),1),Slice());
  Vpad_(sl1)=V;

  //this->state();
 }


// Construct from freq info and shape, w/ initialization
DelayFFT::DelayFFT(Double f0, Double df, Double padBW, 
		   Int nCorr, Int nElem, Int refant, Complex v0) :
  f0_(f0),
  df_(df),
  padBW_(padBW),
  nCorr_(nCorr),
  nPadChan_(Int(0.5+padBW/df)),
  nElem_(nElem),
  refant_(refant), 
  Vpad_(),
  delay_(),
  flag_()
{

  //  cout << "ctor0: " << f0 << " " << df << " " << padBW << " " << nPadChan_ << endl;
  //  cout << "ctor0: " << f0_ << " " << df_ << " " << padBW_ << " " << nPadChan_ << endl;

  //  cout << "pad ch=" << padBW/df_ << " " << nPadChan_*df << endl;

  Vpad_.resize(nCorr_,nPadChan_,nElem_);
  Vpad_.set(v0);

  //  this->state();

}


DelayFFT::DelayFFT(const VisBuffer& vb,Double padBW,Int refant) :
  f0_(vb.frequency()(0)/1.e9),      // GHz
  df_(vb.frequency()(1)/1.e9-f0_),
  padBW_(padBW),
  nCorr_(0),    // set in body
  nPadChan_(Int(0.5+padBW/df_)),
  nElem_(vb.numberAnt()), // antenna-based
  refant_(refant),
  Vpad_(),
  delay_(),
  flag_()
{

  //  cout << "DelayFFT(vb)..." << endl;
  //  cout << "ctor1: " << f0_ << " " << df_ << " " << padBW_ << " " << nPadChan_ << endl;

  // VB facts
  Int nCorr=vb.nCorr();
  Int nChan=vb.nChannel();
  Int nRow=vb.nRow();

  // Discern effective shapes
  nCorr_=(nCorr>1 ? 2 : 1); // number of p-hands
  Int sC=(nCorr>2 ? 3 : 1); // step for p-hands
  
  // Shape the data
  IPosition ip1(3,nCorr_,nPadChan_,nElem_);
  Vpad_.resize(ip1);
  Vpad_.set(Complex(0.0));

  //  this->state();

  // Fill the relevant data
  Int iant=0;
  for (Int irow=0;irow<nRow;++irow) {
    Int a1(vb.antenna1()(irow)), a2(vb.antenna2()(irow));
    if (!vb.flagRow()(irow) && a1!=a2) {

      if (a1==refant)
	iant=a2;
      else if (a2==refant) 
	iant=a1;
      else
	continue;  // an irrelevant baseline

      Slicer sl0(Slice(0,nCorr_,sC),Slice(),Slice(irow,1,1)); // this visCube slice
      Slicer sl1(Slice(),Slice(0,nChan,1),Slice(iant,1,1)); // this Vpad_ slice

      Cube<Complex> vC(vb.visCube()(sl0));

      // Divide by non-zero amps
      Cube<Float> vCa(amplitude(vC));
      vCa(vCa<FLT_EPSILON)=1.0;
      vC/=vCa;

      // Zero flagged channels
      Slicer fsl0(Slice(),Slice(irow,1,1));
      Array<Bool> fl(vb.flag()(fsl0).nonDegenerate(IPosition(1,0)));
      for (Int icor=0;icor<nCorr_;++icor) 
	vC(Slice(icor,1,1),Slice(),Slice()).nonDegenerate(IPosition(1,1))(fl)=Complex(0.0);

      // TBD: apply weights
      //      Matrix<Float> wt(vb.weightMat()(Slice,Slice(irow,1,1)));

      // Acquire this baseline for solving
      Vpad_(sl1)=vC;
    }
  }
  //cout << "...end DelayFFT(vb)" << endl;
	      
}
DelayFFT::DelayFFT(SolveDataBuffer& sdb,Double padBW,Int refant,Int nElem) :
  f0_(sdb.freqs()(0)/1e9),      // GHz
  df_(sdb.freqs()(1)/1e9-f0_),
  padBW_(padBW),
  nCorr_(0),    // set in body
  nPadChan_(Int(0.5+padBW_/df_)),
  nElem_(nElem), // antenna-based
  refant_(refant),
  Vpad_(),
  delay_(),
  flag_()
{

  // VB facts
  Int nCorr=sdb.nCorrelations();
  Int nChan=sdb.nChannels();
  Int nRow=sdb.nRows();

  // Discern effective shapes
  nCorr_=(nCorr>1 ? 2 : 1); // number of p-hands
  Int sC=(nCorr>2 ? 3 : 1); // step for p-hands
  
  // Shape the data
  IPosition ip1(3,nCorr_,nPadChan_,nElem_);
  Vpad_.resize(ip1);
  Vpad_.set(Complex(0.0));

  //  this->state();

  // Fill the relevant data
  Int iant=0;
  for (Int irow=0;irow<nRow;++irow) {
    Int a1(sdb.antenna1()(irow)), a2(sdb.antenna2()(irow));
    if (!sdb.flagRow()(irow) && a1!=a2) {

      if (a1==refant)
	iant=a2;
      else if (a2==refant) 
	iant=a1;
      else
	continue;  // an irrelevant baseline

      Slicer sl0(Slice(0,nCorr_,sC),Slice(),Slice(irow,1,1)); // this visCube slice
      Slicer sl1(Slice(),Slice(0,nChan,1),Slice(iant,1,1)); // this Vpad_ slice

      Cube<Complex> vC(sdb.visCubeCorrected()(sl0));

      // Divide by non-zero amps
      Cube<Float> vCa(amplitude(vC));
      vCa(vCa<FLT_EPSILON)=1.0;
      vC/=vCa;

      // Zero flagged channels
      Cube<Bool> fl(sdb.flagCube()(sl0));
      vC(fl)=Complex(0.0);

      // TBD: apply weights
      //      Matrix<Float> wt(vb.weightMat()(Slice,Slice(irow,1,1)));

      // Acquire this baseline for solving
      Vpad_(sl1)=vC;
    }
  }
  //cout << "...end DelayFFT(sdb)" << endl;
	      
}

void DelayFFT::FFT() {

  //  cout << "DelayFFT::FFT()..." << endl;

  // We always transform only the chan axis (1)
  Vector<Bool> ax(3,false);
  ax(1)=true;

  // Transform
  //  TBD: can Vpad_ be an ArrayLattice?
  ArrayLattice<Complex> c(Vpad_);
  LatticeFFT::cfft0(c,ax,true);

  //  cout << "...end DelayFFT::FFT()" << endl;

}

void DelayFFT::shift(Double f) {

  //  cout << "DelayFFT::shift(f)..." << endl;
    
  Double shift=-(f0_-f)/df_;                    // samples  
  Vector<Double> ph(nPadChan_);  indgen(ph);   // indices
  //  ph-=Double(nPadChan_/2);                     // centered
  ph/=Double(nPadChan_);                       // cycles/sample
  ph*=shift;                                   // cycles
  ph*=(C::_2pi);                               // rad
  Vector<Double> fsh(nPadChan_*2);
  fsh(Slice(0,nPadChan_,2))=cos(ph);
  fsh(Slice(1,nPadChan_,2))=sin(ph);
  Vector<DComplex> csh(nPadChan_);
  RealToComplex(csh,fsh);
  Vector<Complex> sh(nPadChan_);
  convertArray(sh,csh);  // downcovert to Complex

  // Apply to each elem, corr
  for (Int ielem=0;ielem<nElem_;++ielem)
    for (Int icorr=0;icorr<nCorr_;++icorr) {
      Vector<Complex> v(Vpad_.xyPlane(ielem).row(icorr));
      v*=sh;
    }

  //  cout << "... end DelayFFT::shift(f)" << endl;

}



void DelayFFT::add(const DelayFFT& other) {

  //  cout << "DelayFFT::add(x)..." << endl;

  IPosition osh=other.Vpad_.shape();

  //  cout << "add: " << nPadChan_ << " " << other.nPadChan_ << endl;


  AlwaysAssert( (other.nCorr_==nCorr_), AipsError);
  AlwaysAssert( (other.nPadChan_<=nPadChan_), AipsError);
  AlwaysAssert( (other.nElem_==nElem_), AipsError);

  Int oNchan=other.nPadChan_;
  Int lo(0);
  while (lo < nPadChan_) {
    Int nch=min(oNchan,nPadChan_-lo);  // ensure we don't overrun (impossible anyway?)
    Slicer sl0(Slice(),Slice(0,nch,1),Slice());
    Slicer sl1(Slice(),Slice(lo,nch,1),Slice());
    Cube<Complex> v1(Vpad_(sl1)), v0(other.Vpad_(sl0));
    v1+=v0;
    lo+=oNchan;
  }

  //  cout << "...end DelayFFT::add(x)" << endl;

}

void DelayFFT::searchPeak() {

  //  cout << "DelayFFT::searchPeak()..." << endl;

  delay_.resize(nCorr_,nElem_);
  delay_.set(0.0);
  flag_.resize(nCorr_,nElem_);
  flag_.set(true);  // all flagged
  Vector<Float> amp;
  Int ipk;
  Float alo,amax,ahi,fpk;
  for (Int ielem=0;ielem<nElem_;++ielem) {
    for (Int icorr=0;icorr<nCorr_;++icorr) {
      amp=amplitude(Vpad_(Slice(icorr,1,1),Slice(),Slice(ielem,1,1)));
      amax=-1.0;
      ipk=0;
      for (Int ich=0;ich<nPadChan_;++ich) {
	if (amp[ich]>amax) {
	  ipk=ich;
	  amax=amp[ich];
	}
      }

      alo=amp(ipk>0 ? ipk-1 : (nPadChan_-1));
      ahi=amp(ipk<(nPadChan_-1) ? ipk+1 : 0);

      Float denom=(alo-2.*amax+ahi);
      if (amax>0.0 && abs(denom)>0.0) {
	fpk=Float(ipk)+0.5-(ahi-amax)/denom;
	Float delay=fpk/Float(nPadChan_);  // cycles/sample
	if (delay>0.5) delay-=1.0;         // fold
	delay/=df_;                        // nsec

	/*
	cout << "delay=" 
	     << Float(ipk>0 ? ipk-1 : (nPadChan_-1))/Float(nPadChan_)/df_ << "..."
	     << Float(ipk)/Float(nPadChan_)/df_ << "..."
	     << Float(ipk<(nPadChan_-1) ? ipk+1 : 0)/Float(nPadChan_)/df_ << "-->"
	     << delay << "    "
	     << alo << " " << amax << " " << ahi 
	     << endl;
	*/

	delay_(icorr,ielem)=delay;     
	flag_(icorr,ielem)=false;
      }
    }
  }
  
  // delays for ants>refant must be negated!
  if (refant_>-1 && refant_<(nElem_-1)) {  // at least 1 such ant
    Matrix<Float> d2(delay_(Slice(),Slice(refant_+1,nElem_-refant_-1,1)));
    d2*=-1.0f;
  }

  // refant (if specified) is unflagged
  if (refant_>-1)
    flag_(Slice(),Slice(refant_,1,1))=false;

  //  cout << "...end DelayFFT::searchPeak()" << endl;


}

void DelayFFT::state() {

  cout << "DelayFFT::state()..." << endl << " ";
  cout << " f0_=" << f0_ 
       << " df_=" << df_
       << " padBW_=" << padBW_ << endl << " "
       << " nCorr_=" << nCorr_
       << " nPadChan_=" << nPadChan_
       << " nElem_=" << nElem_
       << " refant_=" << refant_ << endl << " "
       << " Vpad_.shape()=" << Vpad_.shape()
       << " delay_.shape()=" << delay_.shape()
       << " flag_.shape()=" << flag_.shape()
       << endl;
    
}




// **********************************************************
//  KJones Implementations
//

KJones::KJones(VisSet& vs) :
  VisCal(vs),             // virtual base
  VisMueller(vs),         // virtual base
  GJones(vs)             // immediate parent
{
  if (prtlev()>2) cout << "K::K(vs)" << endl;

  // Extract per-spw ref Freq for phase(delay) calculation                                                                                                      
  //  TBD: these should be in the caltable!!                                                                                                                    
  MSSpectralWindow msSpw(vs.spectralWindowTableName());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,true);
  KrefFreqs_/=1.0e9;  // in GHz

}

KJones::KJones(String msname,Int MSnAnt,Int MSnSpw) :
  VisCal(msname,MSnAnt,MSnSpw),             // virtual base
  VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
  GJones(msname,MSnAnt,MSnSpw)             // immediate parent
{
  if (prtlev()>2) cout << "K::K(msname,MSnAnt,MSnSpw)" << endl;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  TBD: these should be in the caltable!!

  /* Apparently deprecated; we get it from the ct_ upon setApply...

  NEEDED IN SOLVE CONTEXT?

  MSSpectralWindow msSpw(vs.spectralWindowTableName());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,true);
  KrefFreqs_/=1.0e9;  // in GHz
  */


}

KJones::KJones(const MSMetaInfoForCal& msmc) :
  VisCal(msmc),             // virtual base
  VisMueller(msmc),         // virtual base
  GJones(msmc)             // immediate parent
{
  if (prtlev()>2) cout << "K::K(msmc)" << endl;
}

KJones::KJones(const Int& nAnt) :
  VisCal(nAnt), 
  VisMueller(nAnt),
  GJones(nAnt)
{
  if (prtlev()>2) cout << "K::K(nAnt)" << endl;
}

KJones::~KJones() {
  if (prtlev()>2) cout << "K::~K()" << endl;
}

void KJones::setApply(const Record& apply) {

  // Call parent to do conventional things
  GJones::setApply(apply);

  if (calWt()) 
    logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;

  // Enforce calWt() = false for delays
  calWt()=false;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  from the CalTable
  MSSpectralWindow msSpw(ct_->spectralWindow());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,true);
  KrefFreqs_/=1.0e9;  // in GHz

  /// Re-assign KrefFreq_ according spwmap (if any)
  if (spwMap().nelements()>0) {
    Vector<Double> tmpfreqs;
    tmpfreqs.assign(KrefFreqs_);
    for (uInt ispw=0;ispw<spwMap().nelements();++ispw)
      if (spwMap()(ispw)>-1)
	KrefFreqs_(ispw)=tmpfreqs(spwMap()(ispw));
  }

    
}
void KJones::setApply() {

  // Call parent to do conventional things
  GJones::setApply();

  // Enforce calWt() = false for delays
  calWt()=false;

  // Set the ref freqs to something usable
  KrefFreqs_.resize(nSpw());
  KrefFreqs_.set(0.0);

}

void KJones::setCallib(const Record& callib,
		       const MeasurementSet& selms) {

  // Call parent to do conventional things
  SolvableVisCal::setCallib(callib,selms);

  if (calWt()) 
    logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;

  // Enforce calWt() = false for delays
  calWt()=false;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  from the CalTable
  KrefFreqs_.assign(cpp_->refFreqIn());
  KrefFreqs_/=1.0e9;  // in GHz

  // Re-assign KrefFreq_ according spwmap (if any)
  if (spwMap().nelements()>0) {
    Vector<Double> tmpfreqs;
    tmpfreqs.assign(KrefFreqs_);
    for (uInt ispw=0;ispw<spwMap().nelements();++ispw)
      if (spwMap()(ispw)>-1)
	KrefFreqs_(ispw)=tmpfreqs(spwMap()(ispw));
  }

    
}

void KJones::setSolve(const Record& solve) {

  // Call parent to do conventional things
  GJones::setSolve(solve);

  // Trap unspecified refant:
  if (refant()<0)
    throw(AipsError("Please specify a good reference antenna (refant) explicitly."));

}


void KJones::specify(const Record& specify) {

  return SolvableVisCal::specify(specify);

}

void KJones::calcAllJones() {

  //if (prtlev()>6) cout << "       VJ::calcAllJones()" << endl;

  // Should handle OK flags in this method, and only
  //  do Jones calc if OK

  Vector<Complex> oneJones;
  Vector<Bool> oneJOK;
  Vector<Float> onePar;
  Vector<Bool> onePOK;

  ArrayIterator<Complex> Jiter(currJElem(),1);
  ArrayIterator<Bool>    JOKiter(currJElemOK(),1);
  ArrayIterator<Float>   Piter(currRPar(),1);
  ArrayIterator<Bool>    POKiter(currParOK(),1);

  Double phase(0.0);
  for (Int iant=0; iant<nAnt(); iant++) {

    for (Int ich=0; ich<nChanMat(); ich++) {
      
      oneJones.reference(Jiter.array());
      oneJOK.reference(JOKiter.array());
      onePar.reference(Piter.array());
      onePOK.reference(POKiter.array());

      for (Int ipar=0;ipar<nPar();++ipar) {
	if (onePOK(ipar)) { 
	  phase=2.0*C::pi*onePar(ipar)*(currFreq()(ich)-KrefFreqs_(currSpw()));
	  oneJones(ipar)=Complex(cos(phase),sin(phase));
	  oneJOK(ipar)=true;
	}
      }
      
      // Advance iterators
      Jiter.next();
      JOKiter.next();
      if (freqDepPar()) {
        Piter.next();
        POKiter.next();
      }

    }
    // Step to next antenns's pars if we didn't in channel loop
    if (!freqDepPar()) {
      Piter.next();
      POKiter.next();
    }
  }
}

void KJones::selfSolveOne(VisBuffGroupAcc& vbga) {

  // Forward to MBD solver if more than one VB (more than one spw, probably)
  if (vbga.nBuf()!=1) 
    //    throw(AipsError("KJones can't process multi-VB vbga."));
    this->solveOneVBmbd(vbga);

  // otherwise, call the single-VB solver with the first VB in the vbga
  else
    this->solveOneVB(vbga(0));

}


void KJones::selfSolveOne(SDBList& sdbs) {

  // Forward to MBD solver if more than one VB (more than one spw, probably)
  if (sdbs.nSDB()!=1) 
    this->solveOneSDBmbd(sdbs);

  // otherwise, call the single-VB solver with the first SDB in the SDBList
  else
    this->solveOneSDB(sdbs(0));

}



void KJones::solveOneVBmbd(VisBuffGroupAcc& vbga) {

  Int nbuf=vbga.nBuf();

  Vector<Int> nch(nbuf,0);
  Vector<Double> f0(nbuf,0.0);
  Vector<Double> df(nbuf,0.0);
  Double flo(1e15),fhi(0.0);

  for (Int ibuf=0;ibuf<nbuf;++ibuf) {
    VisBuffer& vb(vbga(ibuf));
    Vector<Double>& chf(vb.frequency());
    nch(ibuf)=vbga(ibuf).nChannel();
    f0(ibuf)=chf(0)/1.0e9;           // GHz
    df(ibuf)=(chf(1)-chf(0))/1.0e9;  // GHz
    flo=min(flo,f0[ibuf]);
    fhi=max(fhi,f0[ibuf]+nch[ibuf]*df[ibuf]);
  }
  Double tbw=fhi-flo;
  //  cout << "tbw = " << tbw << "  (" << flo << "-" << fhi << ")" << endl;

  Double ptbw=tbw*8;  // pad total bw by 8X
  // TBD:  verifty that all df are factors of tbw

  Int nCor=vbga(0).nCorr();

  DelayFFT sumfft(f0[0],min(df),ptbw,(nCor>1 ? 2 : 1),vbga(0).numberAnt(),refant(),Complex(0.0));
  for (Int ibuf=0;ibuf<nbuf;++ibuf) {
    DelayFFT delfft1(vbga(ibuf),ptbw,refant());
    delfft1.FFT();
    delfft1.shift(f0[0]);
    sumfft.add(delfft1);
    delfft1.searchPeak();
    //    cout << ibuf << " "
    //	 << delfft1.delay()(Slice(0,1,1),Slice()) << endl;
  }

  sumfft.searchPeak();

  //  cout << "sum: " << sumfft.delay()(Slice(0,1,1),Slice()) << endl;
  
  // Keep solutions
  Matrix<Float> sRP(solveRPar().nonDegenerate(1));
  sRP=sumfft.delay();  
  Matrix<Bool> sPok(solveParOK().nonDegenerate(1));
  sPok=(!sumfft.flag());

}


void KJones::solveOneSDBmbd(SDBList& sdbs) {

  Int nbuf=sdbs.nSDB();
  
  Vector<Int> nch(nbuf,0);
  Vector<Double> f0(nbuf,0.0);
  Vector<Double> df(nbuf,0.0);
  Double flo(1e15),fhi(0.0);

  for (Int ibuf=0;ibuf<nbuf;++ibuf) {
    SolveDataBuffer& sdb(sdbs(ibuf));
    Vector<Double> chf(sdb.freqs());
    nch(ibuf)=sdbs(ibuf).nChannels();
    f0(ibuf)=chf(0)/1.0e9;           // GHz
    df(ibuf)=(chf(1)-chf(0))/1.0e9;  // GHz
    flo=min(flo,f0[ibuf]);
    fhi=max(fhi,f0[ibuf]+nch[ibuf]*df[ibuf]);
  }
  Double tbw=fhi-flo;

  Double ptbw=tbw*8;  // pad total bw by 8X
  // TBD:  verifty that all df are factors of tbw

  /*
  cout << "tbw = " << tbw << "  (" << flo << "-" << fhi << ")" << " resoln=" << 1.0/tbw
       << ";   ptbw = " << ptbw << " resoln=" << 1/ptbw
       << endl;
  */

  Int nCor=sdbs(0).nCorrelations();
  
  DelayFFT sumfft(f0[0],min(df),ptbw,(nCor>1 ? 2 : 1),nAnt(),refant(),Complex(0.0));
  for (Int ibuf=0;ibuf<nbuf;++ibuf) {
    DelayFFT delfft1(sdbs(ibuf),ptbw,refant(),nAnt());
    delfft1.FFT();
    delfft1.shift(f0[0]);
    sumfft.add(delfft1);

    delfft1.searchPeak();
  }

  sumfft.searchPeak();

  // Keep solutions
  Matrix<Float> sRP(solveRPar().nonDegenerate(1));
  sRP=sumfft.delay();  
  Matrix<Bool> sPok(solveParOK().nonDegenerate(1));
  sPok=(!sumfft.flag());

}


// Do the FFTs
void KJones::solveOneVB(const VisBuffer& vb) {

  Int nChan=vb.nChannel();

  solveRPar()=0.0;
  solveParOK()=false;

  //  cout << "solveRPar().shape() = " << solveRPar().shape() << endl;
  //  cout << "vb.nCorr() = " << vb.nCorr() << endl;
  //  cout << "vb.corrType() = " << vb.corrType() << endl;

  // FFT parallel-hands only
  Int nCorr=vb.nCorr();
  Int nC= (nCorr>1 ? 2 : 1);  // number of parallel hands
  Int sC= (nCorr>2 ? 3 : 1);  // step by 3 for full pol data

  // I/O shapes
  Int fact(8);
  Int nPadChan=nChan*fact;

  IPosition ip0=vb.visCube().shape();
  IPosition ip1=ip0;
  ip1(0)=nC;    // the number of correlations to FFT 
  ip1(1)=nPadChan; // padded channel axis

  // I/O slicing
  Slicer sl0(Slice(0,nC,sC),Slice(),Slice());  
  Slicer sl1(Slice(),Slice(nChan*(fact-1)/2,nChan,1),Slice());

  // Fill the (padded) transform array
  //  TBD: only do ref baselines
  Cube<Complex> vpad(ip1);
  Cube<Complex> slvis=vb.visCube();
  vpad.set(Complex(0.0));
  vpad(sl1)=slvis(sl0);

  // We will only transform frequency axis of 3D array
  Vector<Bool> ax(3,false);
  ax(1)=true;
  
  // Do the FFT
  ArrayLattice<Complex> c(vpad);
  LatticeFFT::cfft(c,ax);        

  // Find peak in each FFT
  Int ipk=0;
  Float amax(0.0);
  Vector<Float> amp;

  //  cout << "Time=" << MVTime(refTime()/C::day).string(MVTime::YMD,7)
  //       << " Spw=" << currSpw() << ":" << endl;

  for (Int irow=0;irow<vb.nRow();++irow) {
    if (!vb.flagRow()(irow) &&
	vb.antenna1()(irow)!=vb.antenna2()(irow) &&
	(vb.antenna1()(irow)==refant() || 
	 vb.antenna2()(irow)==refant()) ) {

      for (Int icor=0;icor<ip1(0);++icor) {
	amp=amplitude(vpad(Slice(icor,1,1),Slice(),Slice(irow,1,1)));
	ipk=1;
	amax=0;
	for (Int ich=1;ich<nPadChan-1;++ich) {
	  if (amp(ich)>amax) {
	    ipk=ich;
	    amax=amp(ich);
	  }
	} // ich

	/*
	cout << vb.antenna1()(irow) << " " << vb.antenna2()(irow) << " "
	     << ntrue(vb.flagCube()(Slice(0,2,3),Slice(),Slice(irow,1,1))) << " "
	     << ntrue(vb.flag()(Slice(),Slice(irow,1,1))) 
	     << endl;
	*/

       	// Derive refined peak (fractional) channel
	// via parabolic interpolation of peak and neighbor channels

	Vector<Float> amp3(amp(IPosition(1,ipk-1),IPosition(1,ipk+1)));
	Float denom(amp3(0)-2.0*amp3(1)+amp3(2));

	if (amax>0.0 && abs(denom)>0) {

	  Float fipk=Float(ipk)+0.5-(amp3(2)-amp3(1))/denom;
	    
	  // Handle FFT offset and scale
	  Float delay=(fipk-Float(nPadChan/2))/Float(nPadChan); // cycles/sample
	  
	  // Convert to cycles/Hz and then to nsec
	  Double df=vb.frequency()(1)-vb.frequency()(0);
	  delay/=df;
	  delay/=1.0e-9;
	  
	  //	  cout << " Antenna ID=";
	  if (vb.antenna1()(irow)==refant()) {
	    //	    cout << vb.antenna2()(irow) 
	    //		 << ", pol=" << icor << " delay(nsec)="<< -delay; 
	    solveRPar()(icor,0,vb.antenna2()(irow))=-delay;
	    solveParOK()(icor,0,vb.antenna2()(irow))=true;
	  }
	  else if (vb.antenna2()(irow)==refant()) {
	    //	    cout << vb.antenna1()(irow) 
	    //		 << ", pol=" << icor << " delay(nsec)="<< delay;
	    solveRPar()(icor,0,vb.antenna1()(irow))=delay;
	    solveParOK()(icor,0,vb.antenna1()(irow))=true;
	  }
	  //	  cout << " (refant ID=" << refant() << ")" << endl;

	  /*	  
	  cout << irow << " " 
	       << vb.antenna1()(irow) << " " 
	       << vb.antenna2()(irow) << " " 
	       << icor << " "
	       << ipk << " "
	       << fipk << " "
	       << delay << " "
	       << endl;
	  */
	} // amax > 0
    /*
	else {
	  cout << "No solution found for antenna ID= ";
	  if (vb.antenna1()(irow)==refant())
	    cout << vb.antenna2()(irow);
	  else if (vb.antenna2()(irow)==refant()) 
	    cout << vb.antenna1()(irow);
	  cout << " in polarization " << icor << endl;
	}
    */
	
      } // icor
    } // !flagrRow, etc.

  } // irow

  // Ensure refant has zero delay and is NOT flagged
  solveRPar()(Slice(),Slice(),Slice(refant(),1,1)) = 0.0;
  solveParOK()(Slice(),Slice(),Slice(refant(),1,1)) = true;

  /*  
  if (nfalse(solveParOK())>0) {
    cout << "NB: No delay solutions found for antenna IDs: ";
    Int nant=solveParOK().shape()(2);
    for (Int iant=0;iant<nant;++iant)
      if (!(solveParOK()(0,0,iant)))
	cout << iant << " ";
    cout << endl;
  }
  */

}
// Do the FFTs
void KJones::solveOneSDB(SolveDataBuffer& sdb) {

  Int nChan=sdb.nChannels();

  solveRPar()=0.0;
  solveParOK()=false;

  //  cout << "solveRPar().shape()   = " << solveRPar().shape() << endl;
  //  cout << "sdb.nCorrelations()    = " << sdb.nCorrelations() << endl;
  //  cout << "sdb.correlationTypes() = " << sdb.correlationTypes() << endl;

  // FFT parallel-hands only
  Int nCorr=sdb.nCorrelations();
  Int nC= (nCorr>1 ? 2 : 1);  // number of parallel hands
  Int sC= (nCorr>2 ? 3 : 1);  // step by 3 for full pol data

  // I/O shapes
  Int fact(8);
  Int nPadChan=nChan*fact;

  IPosition ip0=sdb.visCubeCorrected().shape();
  IPosition ip1=ip0;
  ip1(0)=nC;    // the number of correlations to FFT 
  ip1(1)=nPadChan; // padded channel axis

  // I/O slicing
  Slicer sl0(Slice(0,nC,sC),Slice(),Slice());  
  Slicer sl1(Slice(),Slice(nChan*(fact-1)/2,nChan,1),Slice());

  // Fill the (padded) transform array
  //  TBD: only do ref baselines
  Cube<Complex> vpad(ip1);
  Cube<Complex> slvis=sdb.visCubeCorrected();
  vpad.set(Complex(0.0));
  vpad(sl1)=slvis(sl0);

  // We will only transform frequency axis of 3D array
  Vector<Bool> ax(3,false);
  ax(1)=true;
  
  // Do the FFT
  ArrayLattice<Complex> c(vpad);
  LatticeFFT::cfft(c,ax);        

  // Find peak in each FFT
  Int ipk=0;
  Float amax(0.0);
  Vector<Float> amp;

  //  cout << "Time=" << MVTime(refTime()/C::day).string(MVTime::YMD,7)
  //       << " Spw=" << currSpw() << ":" << endl;

  for (Int irow=0;irow<sdb.nRows();++irow) {
    if (!sdb.flagRow()(irow) &&
	sdb.antenna1()(irow)!=sdb.antenna2()(irow) &&
	(sdb.antenna1()(irow)==refant() || 
	 sdb.antenna2()(irow)==refant()) ) {

      for (Int icor=0;icor<ip1(0);++icor) {
	amp=amplitude(vpad(Slice(icor,1,1),Slice(),Slice(irow,1,1)));
	ipk=1;
	amax=0;
	for (Int ich=1;ich<nPadChan-1;++ich) {
	  if (amp(ich)>amax) {
	    ipk=ich;
	    amax=amp(ich);
	  }
	} // ich

	/*
	cout << sdb.antenna1()(irow) << " " << sdb.antenna2()(irow) << " "
	     << ntrue(sdb.flagCube()(Slice(0,2,3),Slice(),Slice(irow,1,1))) << " "
	     << ntrue(sdb.flag()(Slice(),Slice(irow,1,1))) 
	     << endl;
	*/

       	// Derive refined peak (fractional) channel
	// via parabolic interpolation of peak and neighbor channels

	Vector<Float> amp3(amp(IPosition(1,ipk-1),IPosition(1,ipk+1)));
	Float denom(amp3(0)-2.0*amp3(1)+amp3(2));

	if (amax>0.0 && abs(denom)>0) {

	  Float fipk=Float(ipk)+0.5-(amp3(2)-amp3(1))/denom;
	    
	  // Handle FFT offset and scale
	  Float delay=(fipk-Float(nPadChan/2))/Float(nPadChan); // cycles/sample
	  
	  // Convert to cycles/Hz and then to nsec
	  Double df=sdb.freqs()(1)-sdb.freqs()(0);
	  delay/=df;
	  delay/=1.0e-9;
	  
	  //	  cout << " Antenna ID=";
	  if (sdb.antenna1()(irow)==refant()) {
	    //	    cout << vb.antenna2()(irow) 
	    //		 << ", pol=" << icor << " delay(nsec)="<< -delay; 
	    solveRPar()(icor,0,sdb.antenna2()(irow))=-delay;
	    solveParOK()(icor,0,sdb.antenna2()(irow))=true;
	  }
	  else if (sdb.antenna2()(irow)==refant()) {
	    //	    cout << vb.antenna1()(irow) 
	    //		 << ", pol=" << icor << " delay(nsec)="<< delay;
	    solveRPar()(icor,0,sdb.antenna1()(irow))=delay;
	    solveParOK()(icor,0,sdb.antenna1()(irow))=true;
	  }
	  //	  cout << " (refant ID=" << refant() << ")" << endl;

	  /*	  
	  cout << irow << " " 
	       << vb.antenna1()(irow) << " " 
	       << vb.antenna2()(irow) << " " 
	       << icor << " "
	       << ipk << " "
	       << fipk << " "
	       << delay << " "
	       << endl;
	  */
	} // amax > 0
    /*
	else {
	  cout << "No solution found for antenna ID= ";
	  if (vb.antenna1()(irow)==refant())
	    cout << vb.antenna2()(irow);
	  else if (vb.antenna2()(irow)==refant()) 
	    cout << vb.antenna1()(irow);
	  cout << " in polarization " << icor << endl;
	}
    */
	
      } // icor
    } // !flagrRow, etc.

  } // irow

  // Ensure refant has zero delay and is NOT flagged
  solveRPar()(Slice(),Slice(),Slice(refant(),1,1)) = 0.0;
  solveParOK()(Slice(),Slice(),Slice(refant(),1,1)) = true;

  /*  
  if (nfalse(solveParOK())>0) {
    cout << "NB: No delay solutions found for antenna IDs: ";
    Int nant=solveParOK().shape()(2);
    for (Int iant=0;iant<nant;++iant)
      if (!(solveParOK()(0,0,iant)))
	cout << iant << " ";
    cout << endl;
  }
  */

}

// **********************************************************
//  KcrossJones Implementations
//

KcrossJones::KcrossJones(VisSet& vs) :
  VisCal(vs),             // virtual base
  VisMueller(vs),         // virtual base
  KJones(vs)             // immediate parent
{
  if (prtlev()>2) cout << "Kx::Kx(vs)" << endl;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  TBD: these should be in the caltable!!
  MSSpectralWindow msSpw(vs.spectralWindowTableName());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,true);
  KrefFreqs_/=1.0e9;  // in GHz

}

KcrossJones::KcrossJones(String msname,Int MSnAnt,Int MSnSpw) :
  VisCal(msname,MSnAnt,MSnSpw),             // virtual base
  VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
  KJones(msname,MSnAnt,MSnSpw)              // immediate parent
{
  if (prtlev()>2) cout << "Kx::Kx(msname,MSnAnt,MSnSpw)" << endl;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  TBD: these should be in the caltable!!
  /*  DEPRECATED, because we get it from ct_?
  MSSpectralWindow msSpw(vs.spectralWindowTableName());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,true);
  KrefFreqs_/=1.0e9;  // in GHz
  */
}

KcrossJones::KcrossJones(const MSMetaInfoForCal& msmc) :
  VisCal(msmc),             // virtual base
  VisMueller(msmc),         // virtual base
  KJones(msmc)              // immediate parent
{
  if (prtlev()>2) cout << "Kx::Kx(msmc)" << endl;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  TBD: these should be in the caltable!!
  /*  DEPRECATED, because we get it from ct_?
  MSSpectralWindow msSpw(vs.spectralWindowTableName());
  ROMSSpWindowColumns msCol(msSpw);
  msCol.refFrequency().getColumn(KrefFreqs_,True);
  KrefFreqs_/=1.0e9;  // in GHz
  */
}

KcrossJones::KcrossJones(const Int& nAnt) :
  VisCal(nAnt), 
  VisMueller(nAnt),
  KJones(nAnt)
{
  if (prtlev()>2) cout << "Kx::Kx(nAnt)" << endl;
}

KcrossJones::~KcrossJones() {
  if (prtlev()>2) cout << "Kx::~Kx()" << endl;
}

void KcrossJones::selfSolveOne(VisBuffGroupAcc& vbga) {

  // Trap MBD attempt (NYI)
  if (vbga.nBuf()!=1) 
    throw(AipsError("KcrossJones does not yet support MBD"));
  //    this->solveOneVBmbd(vbga);

  // otherwise, call the single-VB solver with the first VB in the vbga
  else
    this->solveOneVB(vbga(0));

}


void KcrossJones::selfSolveOne(SDBList& sdbs) {

  // Trap MBD attempt (NYI)
  if (sdbs.nSDB()!=1) 
    throw(AipsError("KcrossJones does not yet support MBD"));
  //    this->solveOneVBmbd(vbga);

  // otherwise, call the single-VB solver with the first VB in the vbga
  else
    this->solveOneSDB(sdbs(0));

}


// Do the FFTs
void KcrossJones::solveOneVB(const VisBuffer& vb) {

  solveRPar()=0.0;
  solveParOK()=false;

  Int fact(8);
  Int nChan=vb.nChannel();
  Int nPadChan=nChan*fact;

  // Collapse cross-hands over baseline
  Vector<Complex> sumvis(nPadChan);
  sumvis.set(Complex(0.0));
  Vector<Complex> slsumvis(sumvis(Slice(nChan*(fact-1)/2,nChan,1)));
  Vector<Float> sumwt(nChan);
  sumwt.set(0.0);
  for (Int irow=0;irow<vb.nRow();++irow) {
    if (!vb.flagRow()(irow) &&
	vb.antenna1()(irow)!=vb.antenna2()(irow)) {

      for (Int ich=0;ich<nChan;++ich) {

	if (!vb.flag()(ich,irow)) {
	  // 1st cross-hand
	  slsumvis(ich)+=(vb.visCube()(1,ich,irow)*vb.weightMat()(1,irow));
	  sumwt(ich)+=vb.weightMat()(1,irow);
	  // 2nd cross-hand
	  slsumvis(ich)+=conj(vb.visCube()(2,ich,irow)*vb.weightMat()(2,irow));
	  sumwt(ich)+=vb.weightMat()(2,irow);
	}
      }
    }
  }
  // Normalize the channelized sum
  for (int ich=0;ich<nChan;++ich)
    if (sumwt(ich)>0)
      slsumvis(ich)/=sumwt(ich);
    else
      slsumvis(ich)=Complex(0.0);

  /*
  cout << "slsumvis.nelements() = " << slsumvis.nelements() << endl;
  cout << "sumwt = " << sumwt/sumwt(0) << endl;
  cout << "amp   = " << amplitude(slsumvis) << endl;
  cout << "phase = " << phase(slsumvis)*180.0/C::pi << endl;
  */

  // Do the FFT
  ArrayLattice<Complex> c(sumvis);
  LatticeFFT::cfft(c,true);        
      
  // Find peak in each FFT
  Vector<Float> amp=amplitude(sumvis);

  Int ipk=0;
  Float amax(0.0);
  for (Int ich=0;ich<nPadChan;++ich) {
    if (amp(ich)>amax) {
      ipk=ich;
      amax=amp(ich);
    }
  } // ich
	
  // Derive refined peak (fractional) channel
  // via parabolic interpolation of peak and neighbor channels
  Float fipk=ipk;
  // Interpolate the peak (except at edges!)
  if (ipk>0 && ipk<(nPadChan-1)) {
    Vector<Float> amp3(amp(IPosition(1,ipk-1),IPosition(1,ipk+1)));
    fipk=Float(ipk)+0.5-(amp3(2)-amp3(1))/(amp3(0)-2.0*amp3(1)+amp3(2));
    Vector<Float> pha3=phase(sumvis(IPosition(1,ipk-1),IPosition(1,ipk+1)));
  }

  // Handle FFT offset and scale
  Float delay=(fipk-Float(nPadChan/2))/Float(nPadChan); // cycles/sample

  // Convert to cycles/Hz and then to nsec
  Double df=vb.frequency()(1)-vb.frequency()(0);
  delay/=df;
  delay/=1.0e-9;

  solveRPar()(Slice(0,1,1),Slice(),Slice())=delay;
  solveParOK()=true;

  logSink() << " Time="<< MVTime(refTime()/C::day).string(MVTime::YMD,7)
	    << " Spw=" << currSpw()
	    << " Global cross-hand delay=" << delay << " nsec"
	    << LogIO::POST;
}

// Do the FFTs
void KcrossJones::solveOneSDB(SolveDataBuffer& sdb) {

  solveRPar()=0.0;
  solveParOK()=false;

  Int fact(8);
  Int nChan=sdb.nChannels();
  Int nPadChan=nChan*fact;

  // Collapse cross-hands over baseline
  Vector<Complex> sumvis(nPadChan);
  sumvis.set(Complex(0.0));
  Vector<Complex> slsumvis(sumvis(Slice(nChan*(fact-1)/2,nChan,1)));
  Vector<Float> sumwt(nChan);
  sumwt.set(0.0);
  for (Int irow=0;irow<sdb.nRows();++irow) {
    if (!sdb.flagRow()(irow) &&
	sdb.antenna1()(irow)!=sdb.antenna2()(irow)) {

      for (Int ich=0;ich<nChan;++ich) {
	if (!sdb.flagCube()(1,ich,irow)) {
	  // 1st cross-hand
	  Float& wt(sdb.weightSpectrum()(1,ich,irow));
	  slsumvis(ich)+=(sdb.visCubeCorrected()(1,ich,irow)*wt);
	  sumwt(ich)+=wt;
	}
	if (!sdb.flagCube()(2,ich,irow)) {
	  // 2nd cross-hand
	  Float& wt(sdb.weightSpectrum()(2,ich,irow));
	  slsumvis(ich)+=conj(sdb.visCubeCorrected()(2,ich,irow)*wt);
	  sumwt(ich)+=wt;
	}
      }
    }
  }
  // Normalize the channelized sum
  for (int ich=0;ich<nChan;++ich)
    if (sumwt(ich)>0)
      slsumvis(ich)/=sumwt(ich);
    else
      slsumvis(ich)=Complex(0.0);

  /*  
  cout << "slsumvis.nelements() = " << slsumvis.nelements() << endl;
  cout << "amp   = " << amplitude(slsumvis) << endl;
  cout << "phase = " << phase(slsumvis)*180.0/C::pi << endl;
  cout << "sumwt = " << sumwt/sumwt(0) << endl;
  */

  // Do the FFT
  ArrayLattice<Complex> c(sumvis);
  LatticeFFT::cfft(c,true);        
      
  // Find peak in each FFT
  Vector<Float> amp=amplitude(sumvis);

  Int ipk=0;
  Float amax(0.0);
  for (Int ich=0;ich<nPadChan;++ich) {
    if (amp(ich)>amax) {
      ipk=ich;
      amax=amp(ich);
    }
  } // ich
	
  // Derive refined peak (fractional) channel
  // via parabolic interpolation of peak and neighbor channels
  Float fipk=ipk;
  // Interpolate the peak (except at edges!)
  if (ipk>0 && ipk<(nPadChan-1)) {
    Vector<Float> amp3(amp(IPosition(1,ipk-1),IPosition(1,ipk+1)));
    fipk=Float(ipk)+0.5-(amp3(2)-amp3(1))/(amp3(0)-2.0*amp3(1)+amp3(2));
    Vector<Float> pha3=phase(sumvis(IPosition(1,ipk-1),IPosition(1,ipk+1)));
  }

  // Handle FFT offset and scale
  Float delay=(fipk-Float(nPadChan/2))/Float(nPadChan); // cycles/sample

  // Convert to cycles/Hz and then to nsec
  Double df=sdb.freqs()(1)-sdb.freqs()(0);
  delay/=df;
  delay/=1.0e-9;

  solveRPar()(Slice(0,1,1),Slice(),Slice())=delay;
  solveParOK()=true;

  logSink() << " Time="<< MVTime(refTime()/C::day).string(MVTime::YMD,7)
	    << " Spw=" << currSpw()
	    << " Global cross-hand delay=" << delay << " nsec"
	    << LogIO::POST;
}

// **********************************************************
//  KMBDJones Implementations
//

KMBDJones::KMBDJones(VisSet& vs) :
  VisCal(vs),             // virtual base
  VisMueller(vs),         // virtual base
  KJones(vs)             // immediate parent
{
  if (prtlev()>2) cout << "Kmbd::Kmbd(vs)" << endl;

  // For MBD, the ref frequencies are zero
  KrefFreqs_.resize(nSpw());
  KrefFreqs_.set(0.0);
}

KMBDJones::KMBDJones(String msname,Int MSnAnt,Int MSnSpw) :
  VisCal(msname,MSnAnt,MSnSpw),             // virtual base
  VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
  KJones(msname,MSnAnt,MSnSpw)             // immediate parent
{
  if (prtlev()>2) cout << "Kmbd::Kmbd(msname,MSnAnt,MSnSpw)" << endl;

  // For MBD, the ref frequencies are zero
  KrefFreqs_.resize(MSnSpw);
  KrefFreqs_.set(0.0);
}

KMBDJones::KMBDJones(const MSMetaInfoForCal& msmc) :
  VisCal(msmc),             // virtual base
  VisMueller(msmc),         // virtual base
  KJones(msmc)             // immediate parent
{
  if (prtlev()>2) cout << "Kmbd::Kmbd(msmc)" << endl;

  // For MBD, the ref frequencies are zero
  KrefFreqs_.resize(nSpw());
  KrefFreqs_.set(0.0);
}

KMBDJones::KMBDJones(const Int& nAnt) :
  VisCal(nAnt), 
  VisMueller(nAnt),
  KJones(nAnt)
{

  if (prtlev()>2) cout << "Kmbd::Kmbd(nAnt)" << endl;
  // For MBD, the ref frequencies are zero
  //  TBD: these should be in the caltable!!
  KrefFreqs_.resize(nSpw());
  KrefFreqs_.set(0.0);

}

KMBDJones::~KMBDJones() {
  if (prtlev()>2) cout << "Kmbd::~Kmbd()" << endl;
}


void KMBDJones::setApply(const Record& apply) {
  if (prtlev()>2) cout << "Kmbd::setApply()" << endl;
  KJones::setApply(apply);
  KrefFreqs_.set(0.0);  // MBD is ALWAYS ref'd to zero freq
}


// **********************************************************
//  KAntPosJones Implementations
//

KAntPosJones::KAntPosJones(VisSet& vs) :
  VisCal(vs),             // virtual base
  VisMueller(vs),         // virtual base
  KJones(vs),             // immediate parent
  doTrDelCorr_(false),
  userEterm_(0.0),
  eterm_(0.0)
{
  if (prtlev()>2) cout << "Kap::Kap(vs)" << endl;

  epochref_p="UTC";

}

KAntPosJones::KAntPosJones(String msname,Int MSnAnt,Int MSnSpw) :
  VisCal(msname,MSnAnt,MSnSpw),             // virtual base
  VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
  KJones(msname,MSnAnt,MSnSpw),              // immediate parent
  doTrDelCorr_(false),
  userEterm_(0.0),
  eterm_(0.0)
{
  if (prtlev()>2) cout << "Kap::Kap(msname,MSnAnt,MSnSpw)" << endl;

  epochref_p="UTC";

}

KAntPosJones::KAntPosJones(const MSMetaInfoForCal& msmc) :
  VisCal(msmc),             // virtual base
  VisMueller(msmc),         // virtual base
  KJones(msmc),              // immediate parent
  doTrDelCorr_(false),
  userEterm_(0.0),
  eterm_(0.0)
{
  if (prtlev()>2) cout << "Kap::Kap(msmc)" << endl;

  epochref_p="UTC";

}


KAntPosJones::KAntPosJones(const Int& nAnt) :
  VisCal(nAnt), 
  VisMueller(nAnt),
  KJones(nAnt),
  doTrDelCorr_(false),
  userEterm_(0.0),
  eterm_(0.0)
{
  if (prtlev()>2) cout << "Kap::Kap(nAnt)" << endl;
}

KAntPosJones::~KAntPosJones() {
  if (prtlev()>2) cout << "Kap::~Kap()" << endl;
}

void KAntPosJones::setApply(const Record& apply) {

  //  TBD: Handle missing solutions in spw 0?

  // Force spwmap to all 0  (antpos is not spw-dep)
  //  NB: this is required before calling parents, because
  //   SVC::setApply sets up the CTPatchedInterp with spwMap()
  logSink() << " (" << this->typeName() 
	    << ": Overriding with spwmap=[0] since " << this->typeName() 
	    << " is not spw-dependent)"
	    << LogIO::POST;
  spwMap().assign(Vector<Int>(1,0));

  // Remove spwmap from record, and pass along to generic code
  Record newapply;
  newapply=apply;
  if (newapply.isDefined("spwmap"))
    newapply.removeField("spwmap");
  
  // Call parent to do conventional things
  KJones::setApply(newapply);

  // Arrange for Trop Del correction, if applicable
  //   (force check for CalTable keyword
  if (vlaTrDelCorrApplicable(true))
    initTrDelCorr();

}



void KAntPosJones::setCallib(const Record& callib,
			     const MeasurementSet& selms) 
{

  //  cout << "KAntPosJones::setCallib()" << endl;

  // Call generic to do conventional things
  SolvableVisCal::setCallib(callib,selms);

  if (calWt()) 
    logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;

  // Enforce calWt() = false for delays
  calWt()=false;


  // Force spwmap to all 0  (antpos is not spw-dep)
  //  NB: this is required before calling parents, because
  //   SVC::setApply sets up the CTPatchedInterp with spwMap()
  logSink() << " (" << this->typeName() 
	    << ": Overriding with spwmap=[0] since " << this->typeName() 
	    << " is not spw-dependent)"
	    << LogIO::POST;
  spwMap().assign(Vector<Int>(1,0));

  // Extract per-spw ref Freq for phase(delay) calculation
  //  from the CalTable
  KrefFreqs_.assign(cpp_->refFreqIn());
  KrefFreqs_/=1.0e9;  // in GHz

  // Re-assign KrefFreq_ according spwmap (if any)
  if (spwMap().nelements()>0) {
    Vector<Double> tmpfreqs;
    tmpfreqs.assign(KrefFreqs_);
    for (uInt ispw=0;ispw<spwMap().nelements();++ispw)
      if (spwMap()(ispw)>-1)
	KrefFreqs_(ispw)=tmpfreqs(spwMap()(ispw));
  }

  // Arrange for Trop Del correction, if applicable
  //   (force check for CalTable keyword
  if (vlaTrDelCorrApplicable(true))
    initTrDelCorr();

    
}

void KAntPosJones::specify(const Record& specify) {

  LogMessage message(LogOrigin("KAntPosJones","specify"));

  Vector<Int> spws;
  Vector<Int> antennas;
  Vector<Double> parameters;

  Int Nant(0);

  // Handle old VLA rotation, if necessary
  Bool doVLARot(false);
  Matrix<Double> vlaRot=Rot3D(0,0.0);
  if (specify.isDefined("caltype") ) {
    String caltype=upcase(specify.asString("caltype"));
    if (upcase(caltype).contains("VLA")) {
      doVLARot=true;
      MPosition vlaCenter;
      AlwaysAssert(MeasTable::Observatory(vlaCenter,"VLA"),AipsError);
      Double vlalong=vlaCenter.getValue().getLong();
      //      vlalong=-107.617722*C::pi/180.0;
      cout << "We will rotate specified offsets by VLA longitude = " 
	   << vlalong*180.0/C::pi << endl;
      vlaRot=Rot3D(2,vlalong);
    }
  }
  
  IPosition ip0(3,0,0,0);
  IPosition ip1(3,2,0,0);

  if (specify.isDefined("antenna")) {
    // TBD: the antennas (in order) identifying the solutions
    antennas=specify.asArrayInt("antenna");
    //    cout << "antenna indices = " << antennas << endl;
    Nant=antennas.nelements();
    if (Nant<1) {
      // Use specified values for _all_ antennas implicitly
      Nant=1;   // For the antenna loop below
      ip0(2)=0;
      ip1(2)=nAnt()-1;
    }
    else {
      // Point to first antenna
      ip0(2)=antennas(0);
      ip1(2)=ip0(2);
    }
  }

  if (specify.isDefined("parameter")) {
    // TBD: the actual cal values
    parameters=specify.asArrayDouble("parameter");
    //    cout << "parameters = ]" << parameters << "[" << endl;
  }

  Int npar=parameters.nelements();
  
  // Can't proceed of no parameters were specified
  if (npar==0)
    throw(AipsError("No antenna position corrections specified!"));

  // Must be a multiple of 3
  if (npar%3 != 0)
    throw(AipsError("For antenna position corrections, 3 parameters per antenna are required."));
  
  //  cout << "Shapes = " << parameters.nelements() << " " 
  //       << Nant*3 << endl;

  //  cout << "parameters = " << parameters << endl;

  // Always _ONLY_ spw=0 for antpos corrections
  currSpw()=0;

  // Loop over specified antennas
  Int ipar(0);
  for (Int iant=0;iant<Nant;++iant) {
    if (Nant>1)
      ip1(2)=ip0(2)=antennas(iant);

    // make sure ipar doesn't exceed specified list
    ipar=ipar%npar;
    
    // The current 3-vector of position corrections
    Vector<Double> apar(parameters(IPosition(1,ipar),IPosition(1,ipar+2)));

    // If old VLA, rotate them
    if (doVLARot) {
      cout << "id = " << antennas(iant) << " " << apar;
      apar = product(vlaRot,apar);
      cout << "--(rotation VLA to ITRF)-->" << apar << endl;
    }

    // Loop over 3 parameters, each antenna
    for (Int ipar0=0;ipar0<3;++ipar0) {
      ip1(0)=ip0(0)=ipar0;

      Array<Float> sl(solveAllRPar()(ip0,ip1));
    
      // Acccumulation is addition for ant pos corrections
      sl+=Float(apar(ipar0));
      ++ipar;
    }
  }

  // Store in the memory caltable
  //  (currSpw()=0 is the only one we need)
  keepNCT();

  // Detect if Trop Del correction applicable
  if (vlaTrDelCorrApplicable()) {
    markCalTableForTrDelCorr();
  }

}

// KAntPosJones needs ant-based phase direction and position frame
void KAntPosJones::syncMeta(const VisBuffer& vb) {

  // Call parent (sets currTime())
  KJones::syncMeta(vb);

  phasedir_p=vb.msColumns().field().phaseDirMeas(currField());
  antpos0_p=vb.msColumns().antenna().positionMeas()(0);

  if (doTrDelCorr_)
    // capture azel info
    azel_.reference(vb.azel(currTime()));

}

// KAntPosJones needs ant-based phase direction and position frame
void KAntPosJones::syncMeta2(const vi::VisBuffer2& vb) {

  // Call parent (sets currTime())
  KJones::syncMeta2(vb);

  phasedir_p=vb.subtableColumns().field().phaseDirMeas(currField());
  antpos0_p=vb.subtableColumns().antenna().positionMeas()(0);


  if (doTrDelCorr_)
    // capture azel info
    azel_.reference(vb.azel(currTime()));

}


void KAntPosJones::calcAllJones() {

  if (prtlev()>6) cout << "       Kap::calcAllJones()" << endl;

  // The relevant timestamp 
  MEpoch epoch(Quantity(currTime(),"s"));
  epoch.setRefString(epochref_p);

  // The frame in which we convert our MBaseline from earth to sky and to uvw
  MeasFrame mframe(antpos0_p,epoch,phasedir_p);

  // template MBaseline, that will be used in calculations below
  MBaseline::Ref mbearthref(MBaseline::ITRF,mframe);
  MBaseline mb;
  MVBaseline mvb;
  mb.set(mvb,mbearthref); 

  // A converter that takes the MBaseline from earth to sky frame
  MBaseline::Ref mbskyref(MBaseline::fromDirType(MDirection::castType(phasedir_p.myType())));
  MBaseline::Convert mbcverter(mb,mbskyref);

  Double phase(0.0);
  for (Int iant=0; iant<nAnt(); iant++) {

    Vector<Float> rpars(currRPar().xyPlane(iant).column(0));
    Vector<Double> dpars(rpars.nelements());
    convertArray(dpars,rpars);

    // We need the w offset (in direction of source) implied
    //  by the antenna position correction
    Double dw(0.0);
      
    // Only do complicated calculation if there 
    //   is a non-zero ant pos error
    if (max(abs(rpars))>0.0) {

      // The current antenna's error as an MBaseline (earth frame)
      mvb=MVBaseline(dpars);
      mb.set(mvb,mbearthref);
      
      // Convert to sky frame
      MBaseline mbdir = mbcverter(mb);

      // Get implied uvw
      MVuvw uvw(mbdir.getValue(),phasedir_p.getValue());

      // dw is third element
      dw=uvw.getVector()(2);  // in m

    }
    //if (iant==26) cout << iant << " " << dw << " -> ";

    // Add on the Tropo Delay correction
    if (doTrDelCorr_)
      dw+=calcTrDelError(iant); // still in m
    //if (iant==26) cout << dw;
    
    if (abs(dw)>0.0) {

      // To delay units 
      dw/=C::c;    // to sec
      dw*=1.0e9;   // to nsec

      //if (iant==26) cout << " -> " << dw << endl;

      // Form the complex corrections per chan (freq)
      for (Int ich=0; ich<nChanMat(); ++ich) {
        
	// NB: currFreq() is in GHz
	phase=2.0*C::pi*dw*currFreq()(ich);
	currJElem()(0,ich,iant)=Complex(cos(phase),sin(phase));
	currJElemOK()(0,ich,iant)=true;
	
      }
    }
    else {
      // No correction
      currJElem().xyPlane(iant)=Complex(1.0);
      currJElemOK().xyPlane(iant)=true;
    }

  } // iant


}

bool KAntPosJones::vlaTrDelCorrApplicable(bool checkCalTable) {

  // Nominally OFF
  doTrDelCorr_=false;

  Int nobs(msmc().msmd().nObservations());
  //  MUST be VLA
  Vector<String> obsnames(msmc().msmd().getObservatoryNames());
  for (Int iobs=0;iobs<nobs;++iobs)
    if (!obsnames(iobs).contains("VLA")) return false;
  
  // Reach here only if all obs are VLA

  // Parse boundary dates
  if (MJDlim_.nelements()==0) {
    MJDlim_.resize(2);
    QuantumHolder qh;
    String er;
    qh.fromString(er,MJD0);
    MJDlim_[0]=qh.asQuantity().getValue();
    qh.fromString(er,MJD1);
    MJDlim_[1]=qh.asQuantity().getValue();
  }

  // Are we in affected date range?
  double iMJD=msmc().msmd().getTimeRangesOfObservations()[0].first.getValue().get();  // days
  if (iMJD > MJDlim_[0] &&
      iMJD < MJDlim_[1]) {
    // TURN IT ON!
    doTrDelCorr_=true;
    logSink() << "NB: This EVLA dataset appears to fall within the period" << endl
	      << "  of semester 16B during which the online tropospheric" << endl
	      << "  delay model was mis-applied." << LogIO::WARN << LogIO::POST;
  }
  

  // Check table for user-specified scale, if requested
  //  (setApply context)
  if (checkCalTable) {
    const TableRecord& tr(ct_->keywordSet());
    if (tr.isDefined("VLATrDelCorr")) {
      userEterm_ =tr.asDouble("VLATrDelCorr");
      if (userEterm_==0.0) {
	// keyword value says turn it off!
	if (doTrDelCorr_)
	  logSink() << "Found VLATrDelCorr=0.0 in the antpos caltable; turning trop delay correction OFF."
		    << LogIO::WARN << LogIO::POST;
	doTrDelCorr_=false;
      }
      else {
	if (!doTrDelCorr_)
	  // user (via caltable) is insisting, even if we are out of date range
	  logSink() << "Found VLATrDelCorr keyword in the antpos caltable; turning trop delay correction ON."
		    << LogIO::WARN << LogIO::POST;
	doTrDelCorr_=true;
      }

    } 
    else {
      // Keyword not present in caltable == TURN IT OFF
      doTrDelCorr_=false;
      logSink() << "No VLATrDelCorr keyword in the antpos caltable; turning trop delay correction OFF."
		<< LogIO::WARN << LogIO::POST;
    }
  }

  if (doTrDelCorr_)
    logSink() << "A correction for the online tropospheric delay model error WILL BE APPLIED!"
	      << LogIO::WARN << LogIO::POST;

  return doTrDelCorr_;
}

void KAntPosJones::markCalTableForTrDelCorr() {

  // Only do this if turned on!
  AlwaysAssert(doTrDelCorr_,AipsError);

  logSink() << "Marking antpos caltable to turn ON the trop delay correction."
	    << LogIO::WARN << LogIO::POST;

  // Add a Table keyword to signal the correction
  TableRecord& tr(ct_->rwKeywordSet());
  userEterm_=1.0;  // non-zero is meaningful, 1.0 is nominal
  tr.define("VLATrDelCorr",userEterm_);

}


void KAntPosJones::initTrDelCorr() {

  // Must have turned TrDelCorr on!
  AlwaysAssert(doTrDelCorr_,AipsError);

  // correction scale factor
  eterm_=-1.0e-15;  // s/m   (nominal)

  // Apply user-supplied factor from keyword
  if (userEterm_!=1.0) {
    logSink() << "Applying user-supplied scalefactor="+String::toString<Double>(userEterm_)+" to the trop delay correction."
	      << LogIO::WARN << LogIO::POST;
    eterm_*=userEterm_;
  }

  logSink() << "Tropospheric delay error correction coefficient="+String::toString<Double>(eterm_/1e-12)+" (ps/m)"
	    << LogIO::WARN << LogIO::POST;

  losDist_.resize(nAnt());  // distance from center
  armAz_.resize(nAnt());

  // Local VLA geometry
  Vector<String> sta(msmc().msmd().getAntennaStations());
  Vector<Double> vlapos(msmc().msmd().getObservatoryPosition(0).get("m").getValue("m"));
  Vector<MPosition> mpos(msmc().msmd().getAntennaPositions());

  // Process each antenna
  for (int iant=0;iant<nAnt();++iant) {

    String pre("");
    if (sta(iant).startsWith("EVLA:")) pre="EVLA:";
    if (sta(iant).startsWith("VLA:")) pre="VLA:";

    // Set VLA arm Az (deg)
    if (sta(iant).contains(pre+"N")) {
      armAz_(iant) = 355.0;
    }
    else if (sta(iant).contains(pre+"E")) {
      armAz_(iant) = 115.0;
    }
    else if (sta(iant).contains(pre+"W")) {
      armAz_(iant) = 236.0;
    }

    // Set distance from center (m)
    Vector<Double> ipos(mpos(iant).get("m").getValue("m"));
    losDist_(iant)=sqrt(square(ipos(0)-vlapos(0))+
			square(ipos(1)-vlapos(1))+
			square(ipos(2)-vlapos(2)));  // m from center

  }

  // Handle units
  armAz_*=(C::pi/180.0); // to rad
  losDist_*=eterm_;   // to s (light travel time)
  losDist_*=C::c;     // to m (light travel distance)
}



Double KAntPosJones::calcTrDelError(Int iant) {

  // Time-dep part of the calculation
  Double dgeo(1.0);
  Vector<Double> iazel(azel_(iant).getAngle().getValue());
  Double az=iazel(0);
  Double el=iazel(1);
  dgeo*=cos(az-armAz_(iant));
  dgeo/=tan(el);
  dgeo/=sin(el);

  //  if (iant==26) cout << az*180/C::pi << " " << el*180/C::pi << " " << dgeo << endl;

  return dgeo*losDist_(iant);

}






} //# NAMESPACE CASA - END
