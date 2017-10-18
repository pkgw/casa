//# tAccorJones: test amplitude normalization
//# Copyright (C) 2013,2017
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

#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/BasicMath/Math.h>

#include <synthesis/MeasurementComponents/AccorJones.h>
#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/VisBuffer2.h>

#include <gtest/gtest.h>

#include "VisCalTestBase_GT.h"

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;

// Control verbosity
#define ACCORJONES_TEST_VERBOSE false

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

class AccorJonesTest : public VisCalTestBase {

public:

  Cube<Complex> scale;

  AccorJonesTest() :
    VisCalTestBase(1,1,1,4,4,8,1,true,"circ",true),
    scale(2,1,nAnt,Complex(1.0))
  {
    // canned scaling factors
    scale(0, 0, 1) = 0.9;
    scale(0, 0, 3) = 1.1;
    scale(1, 0, 3) = 1.15;

    if (ACCORJONES_TEST_VERBOSE)
      summary("AccorJonesTest");
  }

};

TEST_F(AccorJonesTest, Basic) {

  AccorJones Aapp(msmc);
  Aapp.setApply();

  ASSERT_EQ(VisCalEnum::JONES,Aapp.matrixType());
  ASSERT_EQ(VisCal::G,Aapp.type());
  ASSERT_EQ(String("Accor Jones"),Aapp.typeName());
  //  ASSERT_EQ(2,Aapp.nPar());
  ASSERT_FALSE(Aapp.freqDepPar());
  ASSERT_FALSE(Aapp.freqDepMat());
  ASSERT_FALSE(Aapp.freqDepCalWt());
  ASSERT_FALSE(Aapp.timeDepMat());
  ASSERT_TRUE(Aapp.isApplied());
  ASSERT_TRUE(Aapp.isSolvable());

  if (ACCORJONES_TEST_VERBOSE)
    Aapp.state();
}

TEST_F(AccorJonesTest, SolveTest) {

  AccorJones Aapp(msmc);
  Aapp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) { 
    Aapp.setMeta(0,1,0.0,ispw,ssvp.freqs(ispw),nChan);
    Aapp.sizeApplyParCurrSpw(nChan);
    Aapp.setApplyParCurrSpw(scale, true, false); // corrupt
  }

  AccorJones Asol(msmc);
  Record solvePar;
  solvePar.define("table",String("test.Accor"));
  solvePar.define("solint",String("inf"));
  Asol.setSolve(solvePar);

  SDBList sdbs;
  for (vi2.originChunks();vi2.moreChunks();vi2.nextChunk()) {
    for (vi2.origin();vi2.more();vi2.next()) {

      Int ispw=vb2->spectralWindows()(0);
      Int obsid(vb2->observationId()(0));
      Int scan(vb2->scan()(0));
      Double timestamp(vb2->time()(0));
      Int fldid(vb2->fieldId()(0));
      Vector<Double> freqs(vb2->getFrequencies(0));
      Vector<Int> a1(vb2->antenna1());
      Vector<Int> a2(vb2->antenna2());

      vb2->resetWeightsUsingSigma();
      vb2->setVisCubeCorrected(vb2->visCube());
      vb2->setFlagCube(vb2->flagCube());

      Aapp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Aapp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)

      Cube<Float> ratio =
	amplitude(vb2->visCubeCorrected()) / amplitude(vb2->visCube());

      for (Int irow=0;irow<vb2->nRows();++irow) {
	for (Int icor=0;icor<nCorr;icor+=3) {
	  Int p1 = icor == 0 ? 0 : 1;
	  Int p2 = icor == 0 ? 0 : 1;
	  Complex temp = scale(p1,0,a1(irow)) * scale(p2,0,a2(irow));
	  Float diff = ratio(icor,0,irow) - abs(temp);
	  ASSERT_TRUE(diff<1.0e-6);
	}
      }

      sdbs.add(*vb2);
    }
  }

  Asol.setMeta(sdbs.aggregateObsId(),
	      sdbs.aggregateScan(),
	      sdbs.aggregateTime(),
	      sdbs.aggregateSpw(),
	      sdbs.freqs(),
	      sdbs.aggregateFld());
  Asol.sizeSolveParCurrSpw(sdbs.nChannels()); 

  Asol.selfSolveOne(sdbs);

  Cube<Float> soldiff=amplitude(Asol.solveAllCPar()-scale);

  ASSERT_TRUE(allNearAbs(soldiff,0.0f,1e-6));
}
