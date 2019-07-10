//# tFringeJones_GT.cc: Tests the FringeJones
//# Copyright (C) 1995,1999,2000,2001
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
//#
//# $Id$

//#include <casa/aips.h>
//#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <synthesis/MeasurementComponents/FringeJones.h>

#include <gtest/gtest.h>

#include "VisCalTestBase_GT.h"

// <summary>
// Test program for FringeJones
// </summary>


using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;


#define FRINGEJONES_TEST_VERBOSE false

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

class DelayRateFFTTest : public ::testing::Test {

public:
    virtual Array<Complex> appdel(Int nchan, Int nt, Float f0, Float df, Float dt, Float delay, Float fringerate) {
      // f0 is unused
      Matrix<Complex> cph(IPosition(2, nt, nchan));
      for (Int i=0; i!=nt; i++) {
        for (Int j=0; j!=nchan; j++) {
          Float p = C::_2pi*(-i*dt*fringerate*f0*1e9 -j*df*delay );
          cph(IPosition(2, i, j)) = Complex(cos(p),sin(p));
       }
    }
    return cph;
  }
};

TEST_F(DelayRateFFTTest, BasicDelayRateFFTTest) {
  Int nPadfactor(4);
  Int nchan(32), nt(20), nelem(2);
  Float df(0.1), dt(2.0);;
  Float f0(9.0);
  Float delay(-2.0);
  Float rate(4.0e-13);
  SDBList s = SDBList();
  
  Array<Complex> Vobs0(IPosition(4, 1, nelem, nt, nchan));
  const Array<Complex>& V2 = Vobs0(Slicer(IPosition(4, 0, nelem-1, 0, 0),
                                          IPosition(4, 0, nelem-1, nt-1, nchan-1),
                                          IPosition(4, 1, 1, 1, 1),
                                          Slicer::endIsLast)).nonDegenerate();
  const Array<Complex>& V3 = this->appdel(nchan, nt, f0, df, dt, delay, rate);
  V2.nonDegenerate() = V3;

  Array<Double> delayWindow(IPosition(1, 2));
  Array<Double> rateWindow(IPosition(1, 2));
  delayWindow(IPosition(1, 0)) = -100.0;
  delayWindow(IPosition(1, 1)) = +100.0;
  rateWindow(IPosition(1, 0)) = -100.0;
  rateWindow(IPosition(1, 1)) = +100.0;

  IPosition ds = delayWindow.shape();
  
  DelayRateFFT drfft0(Vobs0, nPadfactor, f0, df, dt, s, delayWindow, rateWindow);
  drfft0.FFT();
  drfft0.searchPeak();
  Float rate_resn = 1.0f/(nPadfactor*nt*dt*1e9*f0);
  Float delay_resn = 1.0f/(nPadfactor*nchan*df);
  
  if (FRINGEJONES_TEST_VERBOSE) {
    cerr << boolalpha;
    cerr << "param = " << drfft0.param() << endl;
    cerr << "delay = " << drfft0.delay()(0,1) << endl;
    cerr << "delay resolution = " << delay_resn << endl;
    cerr << "scaled error in delay = " << (drfft0.delay()(0, 1) - delay)/delay_resn << endl;
    cerr << "rate = " << drfft0.rate()(0, 1) << endl;
    cerr << "rate resolution = " << rate_resn << endl;
    cerr << "scaled error in rate = " << (drfft0.rate()(0, 1) - rate)/rate_resn << endl;
  }
  ASSERT_TRUE(allNearAbs(drfft0.delay()(0,1), delay, 0.5*delay_resn));
  ASSERT_TRUE(allNearAbs(drfft0.rate()(0,1), rate, 0.5*rate_resn));
}


TEST_F(DelayRateFFTTest, BasicDelayRateFFTTest2) {
  Int nPadfactor(4);
  Int nchan(32), nt(20), nelem(2);
  Float df(0.1), dt(2.0);;
  Float f0(9.0);
  Float delay(-2.8);
  Float rate(4.0e-13);
  SDBList s = SDBList();
  
  Array<Complex> Vobs0(IPosition(4, 1, nelem, nt, nchan));
  const Array<Complex>& V2 = Vobs0(Slicer(IPosition(4, 0, nelem-1, 0, 0),
                                          IPosition(4, 0, nelem-1, nt-1, nchan-1),
                                          IPosition(4, 1, 1, 1, 1),
                                          Slicer::endIsLast)).nonDegenerate();
  const Array<Complex>& V3 = this->appdel(nchan, nt, f0, df, dt, delay, rate);
  V2.nonDegenerate() = V3;

  Array<Double> delayWindow(IPosition(1, 2));
  Array<Double> rateWindow(IPosition(1, 2));
  Float minDelay = -1.0;
  Float maxDelay = +1.0;
  delayWindow(IPosition(1, 0)) = minDelay;
  delayWindow(IPosition(1, 1)) = maxDelay;
  rateWindow(IPosition(1, 0)) = -100.0;
  rateWindow(IPosition(1, 1)) = +100.0;

  IPosition ds = delayWindow.shape();
  
  DelayRateFFT drfft0(Vobs0, nPadfactor, f0, df, dt, s, delayWindow, rateWindow);
  drfft0.FFT();
  drfft0.searchPeak();
  Float delay_out = drfft0.delay()(0, 1);
  Float rate_resn = 1.0f/(nPadfactor*nt*dt*1e9*f0);
  Float delay_resn = 1.0f/(nPadfactor*nchan*df);
  
  if (FRINGEJONES_TEST_VERBOSE) {
    cerr << boolalpha;
    cerr << "param = " << drfft0.param() << endl;
    cerr << "delay out = " << delay_out << endl;
    cerr << "delay in = " << delay << endl;
    cerr << "delay = " << delay << endl;
    cerr << "delay resolution = " << delay_resn << endl;
    cerr << "scaled error in delay = " << (delay_out - delay)/delay_resn << endl;
    cerr << "rate = " << drfft0.rate()(0, 1) << endl;
    cerr << "rate resolution = " << rate_resn << endl;
    cerr << "scaled error in rate = " << (drfft0.rate()(0, 1) - rate)/rate_resn << endl;
  }
  ASSERT_TRUE((delay_out < maxDelay) && (delay_out >= minDelay));
  ASSERT_TRUE(allNearAbs(drfft0.rate()(0, 1), rate, 0.5*rate_resn));
}



class FringeJonesTest : public VisCalTestBase {
public:
  // test values for solutions
  Cube<Float> fpar;

  FringeJonesTest() :
    VisCalTestBase(1,  // nfield
		   1,  // nscan
		   1,  // nspw
		   7,  // nant
		   4,  // ncorr
		   64, // nchan per spw
		   16, // ntime per scan
		   false),  // unpolarized
    // nPar, 1, {1 | nAntennas}
    fpar(8, 1, VisCalTestBase::nAnt, 0.0f)  // 8 pars per antenna
  {
      // Add FringeJonesTest specific init
      //  e.g., fill fpar with interesting values
    for (Int i=1; i!=VisCalTestBase::nAnt; i++) {
      // 1, 4 are delay.
      fpar(1, 0, i) = 2.3;
      fpar(5, 0, i) = -1.7;
      // Parameters 2 and 6 are rate.  But VisCal::setMeta currently
      // only supports a single double for time, and this meta data is
      // used to generate all the Jones matrices for the FringeJones
      // (or any other VisCal subclass) calibrater; I don't feel
      // comfortable enough with the calibration stack to try to
      // intervene on this, so for now we can only test zero rates for
      // the FringeJones class proper.
      fpar(2, 0, i) = 0.0;
      fpar(6, 0, i) = 0.0;
    }
    // uncomment to see data shape summary from
    //VisCalTestBase::summary("FringeJonesTest");  
  }
};


TEST_F(FringeJonesTest, FringeJonesApplyState) {

  FringeJones ff(VisCalTestBase::msmc);  
  ff.setApply();

  ASSERT_EQ(VisCalEnum::JONES,ff.matrixType());
  ASSERT_EQ(VisCal::K,ff.type());
  ASSERT_EQ(String("Fringe Jones"),ff.typeName());
  //ASSERT_EQ(6,ff.nPar());
  ASSERT_FALSE(ff.freqDepPar());
  ASSERT_TRUE(ff.freqDepMat());
  ASSERT_FALSE(ff.freqDepCalWt());
  ASSERT_EQ(True,ff.timeDepMat());
  ASSERT_TRUE(ff.isApplied());
  ASSERT_TRUE(ff.isSolvable());
  ASSERT_TRUE(ff.useGenericGatherForSolve());
  ASSERT_FALSE(ff.useGenericSolveOne());
  
}



TEST_F(FringeJonesTest, FringeJones_selfSolveOneTest) {
  // Apply-able FringeJones
  FringeJones FJapp(msmc); // "<noms>",nAnt,nSpw);
  FJapp.setApply();
  // FJapp.setPrtlev(7);
  
  // Fill FJapp with actual parameters
  for (Int ispw=0;ispw<nSpw;++ispw) {
    FJapp.setMeta(0,1,0.0,
                 ispw,ssvp.freqs(ispw),
                 nChan);
    FJapp.sizeApplyParCurrSpw(nChan);
    // Enabled now phase model implemented...
    FJapp.setApplyParCurrSpw(fpar,true,false);  // don't invert
  }
  FringeJones FJsol(VisCalTestBase::msmc);
  // FJsol.setPrtlev(7);
  Record solvePar;
  solvePar.define("table",String("test.Fringe"));  // not used
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Array<Int> refant(IPosition(1,3));
  refant(IPosition(1, 0)) = 12;
  refant(IPosition(1, 1)) = 0;
  refant(IPosition(1, 2)) = 1;
  if (FRINGEJONES_TEST_VERBOSE) {
      cerr << "Refant " << refant << endl;
  }
  solvePar.define("refant",refant);
  solvePar.define("globalsolve", true);
  solvePar.define("weightfactor", 2);
  solvePar.define("niter", 20);
  solvePar.define("zerorates", true);
  Array<Double> delayWindow(IPosition(1, 2));
  Array<Double> rateWindow(IPosition(1, 2));
  delayWindow(IPosition(1, 0)) = -100.0;
  delayWindow(IPosition(1, 1)) = +100.0;
  rateWindow(IPosition(1, 0)) = -100.0;
  rateWindow(IPosition(1, 1)) = +100.0;
  solvePar.define("delaywindow", delayWindow);
  solvePar.define("ratewindow", rateWindow);

  FJsol.setSolve(solvePar);


  SDBList sdbs;
  Double reftime;

  Bool isFirst = true;
  for (vi2.originChunks();vi2.moreChunks();vi2.nextChunk()) {
    for (vi2.origin();vi2.more();vi2.next()) {
      if (isFirst) {
          reftime = vb2->time()(0);
          isFirst = false;
      }
      Int ispw=vb2->spectralWindows()(0);
      Int obsid(vb2->observationId()(0));
      Int scan(vb2->scan()(0));
      Int fldid(vb2->fieldId()(0));
      Vector<Double> freqs(vb2->getFrequencies(0));
      Vector<Int> a1(vb2->antenna1());
      Vector<Int> a2(vb2->antenna2());

      vb2->resetWeightsUsingSigma();
      vb2->setVisCubeCorrected(vb2->visCube());
      vb2->setFlagCube(vb2->flagCube());

      // Corrupt with FJapp
      // cerr << "Corrupting..." << endl;
      FJapp.setMeta(obsid,scan,reftime,
		    ispw,freqs,
		    fldid);
      FJapp.correct2(*vb2,false,false,false); // (trial?, doWtSp?, dosync)
      
      // Add vb2 to the SDBList
      sdbs.add(*vb2);
    }
  }


  
  // Setup meta & sizes for the solve
  FJsol.setMeta(sdbs.aggregateObsId(),
		sdbs.aggregateScan(),
		sdbs.aggregateTime(),
		sdbs.aggregateSpw(),
		sdbs.freqs(),
		sdbs.aggregateFld());
  FJsol.sizeSolveParCurrSpw(sdbs.nChannels()); 
  FJsol.selfSolveOne(sdbs);

  // Add comparison tests here
  Matrix<Float> p = FJsol.solveRPar().nonDegenerate();
  Float delay1 = 2.3;
  Float delay2 = -1.7;
  Float rate1 = 0.0;
  Float rate2 = 0.0;    

  if (FRINGEJONES_TEST_VERBOSE) {
    cerr << "delay1 results " << p(1,1) << endl;
    cerr << "delay2 results " << p(4,1) << endl;
    cerr << "rate1 results " << p(2,1) << endl;
    cerr << "rate2 results " << p(5,1) << endl;
    cerr << "Parameters out: " << p << endl;
  }
  ASSERT_TRUE(allNearAbs(p(1, 1), delay1, 2e-2));
  ASSERT_TRUE(allNearAbs(p(5, 1), delay2, 2e-2));
  
  ASSERT_TRUE(allNearAbs(p(2, 1), rate1, 1e-5));
  ASSERT_TRUE(allNearAbs(p(6, 1), rate2, 1e-5));

  ASSERT_TRUE(FJsol.refant()==0);
  
}
