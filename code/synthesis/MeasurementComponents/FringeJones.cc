//# FringeJones.cc: Implementation of FringeJones
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

#include <synthesis/MeasurementComponents/FringeJones.h>

#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/VisBuffAccumulator.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <synthesis/MeasurementEquations/VisEquation.h>  // *
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>

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



using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN


// **********************************************************
//  CTRateAwareTimeInterp1 Implementations
//

CTRateAwareTimeInterp1::CTRateAwareTimeInterp1(NewCalTable& ct,
					       const casacore::String& timetype,
					       casacore::Array<Float>& result,
					       casacore::Array<Bool>& rflag) :
  CTTimeInterp1(ct,timetype,result,rflag)
{}

// Destructor (nothing to do locally)
CTRateAwareTimeInterp1::~CTRateAwareTimeInterp1() {}

Bool CTRateAwareTimeInterp1::interpolate(Double newtime) {
  
  // Call generic first
  if (CTTimeInterp1::interpolate(newtime)) {
    // Only if generic yields new calibration
    // NB: lastWasExact_=exact in generic
    applyPhaseRate(timeType().contains("nearest") || lastWasExact_);
    return true;
  }
  else
    // No change
    return false;

}

// Do the phase rate math
void CTRateAwareTimeInterp1::applyPhaseRate(Bool single)
{
  Int ispw=mcols_p->spwId()(0);
  MSSpectralWindow msSpw(ct_.spectralWindow());
  ROMSSpWindowColumns msCol(msSpw);
  Vector<Double> refFreqs;
  msCol.refFrequency().getColumn(refFreqs,True);

  //  cout << "time = " << (currTime_ - timeRef_) << endl;

  if (single) {
    for (Int ipol=0;ipol<2;ipol++) {
      Double dtime=(currTime_-timeRef_)-timelist_(currIdx_);
      Double phase=result_(IPosition(2,ipol*3,0));
      Double rate=result_(IPosition(2,ipol*3+2,0));
      phase+=2.0*C::pi*rate*refFreqs(ispw)*dtime;
      result_(IPosition(2,ipol*3,0))=phase;
    }
  } else {
    Vector<uInt> rows(2); indgen(rows); rows+=uInt(currIdx_);
    Cube<Float> r(mcols_p->fparamArray("",rows));

    Vector<Double> dtime(2);
    dtime(0)=(currTime_-timeRef_)-timelist_(currIdx_);
    dtime(1)=(currTime_-timeRef_)-timelist_(currIdx_+1);
    Double wt=dtime(1) / (dtime(1)-dtime(0));


    for (Int ipol=0;ipol<2;ipol++) {
      Vector<Double> phase(2), rate(2);
      phase(0)=r.xyPlane(0)(IPosition(2,ipol*3,0));
      phase(1)=r.xyPlane(1)(IPosition(2,ipol*3,0));
      rate(0)=r.xyPlane(0)(IPosition(2,ipol*3+2,0));
      rate(1)=r.xyPlane(1)(IPosition(2,ipol*3+2,0));

      phase(0)+=2.0*C::pi*rate(0)*refFreqs(ispw)*dtime(0);
      phase(1)+=2.0*C::pi*rate(1)*refFreqs(ispw)*dtime(1);

      Vector<Complex> ph(2);
      ph(0)=Complex(cos(phase(0)),sin(phase(0)));
      ph(1)=Complex(cos(phase(1)),sin(phase(1)));
      ph(0)=Float(wt)*ph(0) + Float(1.0-wt)*ph(1);
      result_(IPosition(2,ipol*3,0))=arg(ph(0));
    }
  }
}




// **********************************************************
//  FringeJones Implementations
//

FringeJones::FringeJones(VisSet& vs) :
  VisCal(vs),             // virtual base
  VisMueller(vs),         // virtual base
  SolvableVisJones(vs)       // immediate parent
{
  if (prtlev()>2) cout << "FringeJones::FringeJones(vs)" << endl;
}

FringeJones::FringeJones(String msname,Int MSnAnt,Int MSnSpw) :
  VisCal(msname,MSnAnt,MSnSpw),             // virtual base
  VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
  SolvableVisJones(msname,MSnAnt,MSnSpw)    // immediate parent
{
  if (prtlev()>2) cout << "FringeJones::FringeJones(msname,MSnAnt,MSnSpw)" << endl;
}

FringeJones::FringeJones(const MSMetaInfoForCal& msmc) :
  VisCal(msmc),             // virtual base
  VisMueller(msmc),         // virtual base
  SolvableVisJones(msmc)    // immediate parent
{
  if (prtlev()>2) cout << "FringeJones::FringeJones(msmc)" << endl;
}

FringeJones::FringeJones(const Int& nAnt) :
  VisCal(nAnt), 
  VisMueller(nAnt),
  SolvableVisJones(nAnt)
{
  if (prtlev()>2) cout << "FringeJones::FringeJones(nAnt)" << endl;
}

FringeJones::~FringeJones() {
  if (prtlev()>2) cout << "FringeJones::~FringeJones()" << endl;
}

void FringeJones::setApply(const Record& apply) {

  // Call parent to do conventional things
  SolvableVisJones::setApply(apply);

  if (calWt()) 
    logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;

  // Enforce calWt() = false for delays
  calWt()=false;

  // Extract per-spw ref Freq for phase(delay) calculation
  //  from the CalTable
  // TBD:  revise as per refFreq decisions
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

void FringeJones::setCallib(const Record& callib,
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

void FringeJones::setSolve(const Record& solve) {

  // Call parent to do conventional things
  SolvableVisJones::setSolve(solve);

  // Trap unspecified refant:
  if (refant()<0)
    throw(AipsError("Please specify a good reference antenna (refant) explicitly."));

}

void FringeJones::calcAllJones() {

  if (prtlev()>6) cout << "       FringeJones::calcAllJones()" << endl;

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

  Double phase;
  for (Int iant=0; iant<nAnt(); iant++) {

    for (Int ich=0; ich<nChanMat(); ich++) {
      
      oneJones.reference(Jiter.array());
      oneJOK.reference(JOKiter.array());
      onePar.reference(Piter.array());
      onePOK.reference(POKiter.array());

      for (Int ipar=0;ipar<nPar();ipar+=3) {
        if (onePOK(ipar)) {
          phase=onePar(ipar);
          phase+=2.0*C::pi*onePar(ipar+1)*
            (currFreq()(ich)-KrefFreqs_(currSpw()));
          phase+=2.0*C::pi*onePar(ipar+2)*KrefFreqs_(currSpw())*1e9*
            (currTime() - refTime());
          oneJones(ipar/3)=Complex(cos(phase),sin(phase));
          oneJOK(ipar/3)=True;
        }
      }
      
      // Advance iterators
      Jiter.next();
      JOKiter.next();
    }
    // Step to next antenns's pars
    Piter.next();
    POKiter.next();
  }
}

void FringeJones::selfSolveOne(VisBuffGroupAcc&) {

  // NB: Use "FringeJones::selfSolveOne(SDBList& sdbs)" instead!
  throw(AipsError("VisBuffGroupAcc is invalid for FringeJones"));

}


void FringeJones::selfSolveOne(SDBList& sdbs) {

  Int nSDB=sdbs.nSDB();
  cout << "FringeJones::selfSolveOne: nSDB = " << nSDB << endl;   // remove this

  // Implement actual solve here!!!

  //  E.g., gather data from visCubeCorrected in the SolveDataBuffers 
  //   in the SDBList and feed to solving code
  //  Typically, each SolveDataBuffer will contain a single timestamp 
  //    in a single spw
  //  E.g., see KJones::solveOneSDBmbd(SDBList&)

}

} //# NAMESPACE CASA - END
