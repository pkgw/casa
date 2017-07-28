//# tGJones_GT.cc: Tests the GJones viscal
//# Copyright (C) 1995,1999,2000,2001,2016
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

#include <casa/aips.h>
//#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/SimpleSimVi2.h>
//#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/OS/Timer.h>
#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <synthesis/MeasurementComponents/DJones.h>
#include <synthesis/MeasurementComponents/KJones.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>

#include <gtest/gtest.h>

#define SHOWSTATE false
using namespace casacore;
using namespace casa;
using namespace casa::vi;

class VisCalTest : public ::testing::Test {

public:
  
  VisCalTest() :
    nFld(1),
    nScan(1),
    nSpw(1),
    nAnt(5),
    nCorr(4),
    nChan(1,32),
    ss(nFld,nScan,nSpw,nAnt,nCorr,Vector<Int>(1,1),nChan),
    msmc(ss)
  {}

  Int nFld,nScan,nSpw,nAnt,nCorr;
  Vector<Int> nChan;
  SimpleSimVi2Parameters ss;
  MSMetaInfoForCal msmc;

};

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




TEST_F(VisCalTest, GJonesApplyState) {
  
  VisCal *G = new GJones(msmc);
  G->setApply();

  G->setMeta(0,0,0.0,
	     0,ss.freqs(0),
	     0);
  G->sizeApplyParCurrSpw(ss.nChan_(0));
  G->setDefApplyParCurrSpw(true,true);  // sync, w/ doInv=T

  if (SHOWSTATE)
    G->state();

  ASSERT_EQ(VisCalEnum::JONES,G->matrixType());
  ASSERT_EQ(VisCal::G,G->type());
  ASSERT_EQ(String("G Jones"),G->typeName());
  ASSERT_EQ(2,G->nPar());
  ASSERT_FALSE(G->freqDepPar());
  ASSERT_FALSE(G->freqDepMat());
  ASSERT_FALSE(G->freqDepCalWt());
  ASSERT_FALSE(G->timeDepMat());
  ASSERT_TRUE(G->isApplied());
  ASSERT_TRUE(G->isSolvable());

  /*
  IPosition sh(3,2,1,nAnt);  // nChan=1 for G
  ASSERT_TRUE(G->currCPar().shape()==sh);
  ASSERT_TRUE(G->currParOK().shape()==sh);
  ASSERT_TRUE(G->currJElem().shape()==sh);
  ASSERT_TRUE(G->currJElemOK().shape()==sh);
  ASSERT_EQ(G->currParOK().data(),G->currJElemOK().data()); // ok addr equal
  */

  delete G;
}


TEST_F(VisCalTest, GJonesSolveState) {

  //  MSMetaInfoForCal msmc(ss);
  SolvableVisCal *G = new GJones(msmc);

  Record solvePar;
  String caltablename("test.G"); solvePar.define("table",caltablename);
  String solint("int");          solvePar.define("solint",solint);
  Vector<Int> refantlist(1,3);   solvePar.define("refant",refantlist);

  G->setSolve(solvePar);

  G->setMeta(0,0,0.0,
	     0,ss.freqs(0),
	     0);
  G->sizeSolveParCurrSpw(ss.nChan_(0));
  G->setDefSolveParCurrSpw(true);

  if (SHOWSTATE)
    G->state();

  ASSERT_EQ(VisCalEnum::JONES,G->matrixType());
  ASSERT_EQ(VisCal::G,G->type());
  ASSERT_EQ(String("G Jones"),G->typeName());
  ASSERT_EQ(2,G->nPar());
  ASSERT_FALSE(G->freqDepPar());
  ASSERT_FALSE(G->freqDepMat());
  ASSERT_FALSE(G->freqDepCalWt());
  ASSERT_FALSE(G->timeDepMat());
  ASSERT_FALSE(G->isApplied());
  ASSERT_TRUE(G->isSolved());
  ASSERT_TRUE(G->isSolvable());
  
  ASSERT_EQ(caltablename,G->calTableName());
  ASSERT_EQ(solint,G->solint());
  ASSERT_EQ(refantlist[0],G->refant());
  ASSERT_TRUE(allEQ(refantlist,G->refantlist()));
  

  delete G;
}



TEST_F(VisCalTest, GJonesSpecifyTest) {
  
  SolvableVisJones *G = new GJones(msmc);

  String caltablename("GJonesSpecifyTest.G");
  Record spec;
  spec.define("caltable",caltablename);

  G->setSpecify(spec);


  if (SHOWSTATE)
    G->state();

  ASSERT_EQ(VisCalEnum::JONES,G->matrixType());
  ASSERT_EQ(VisCal::G,G->type());
  ASSERT_EQ(String("G Jones"),G->typeName());
  ASSERT_EQ(2,G->nPar());
  ASSERT_FALSE(G->freqDepPar());
  ASSERT_FALSE(G->freqDepMat());
  ASSERT_FALSE(G->freqDepCalWt());
  ASSERT_FALSE(G->timeDepMat());
  ASSERT_FALSE(G->isApplied());
  ASSERT_FALSE(G->isSolved());
  ASSERT_TRUE(G->isSolvable());

  Float dpr(180.0/C::pi);
  Int nTimes(35);

  Int nAnt=ss.nAnt_;

  Cube<Complex> g(G->solveAllCPar());
  for (Int iant=0;iant<nAnt;++iant) {
    Float polph=-iant*10.0/dpr;
    g(Slice(1,1,1),Slice(),Slice(iant,1,1))=Complex(cos(polph),sin(polph));
  }

  Int refant(2);
  G->refantlist().assign(Vector<Int>(1,refant));

  for (Int itime=0;itime<nTimes;++itime) {
    Double t(itime);
    for (Int iant=1;iant<nAnt;++iant) {
      Float ph(itime*iant/5/dpr);
      Complex dg(cos(ph),sin(ph));
      g(Slice(),Slice(),Slice(iant,1,1))=g(Slice(),Slice(),Slice(iant,1,1))*dg;
    }

    // Flag 5 intervals on refant, 2nd pol
    if (itime>24 && itime<30) 
      G->solveParOK()(Slice(1,1,1),Slice(),Slice(refant,1,1))=false;
    else
      G->solveParOK()(Slice(1,1,1),Slice(),Slice(refant,1,1))=true;
      

    G->setMeta(0,0,t,0,ss.freqs(0),0);
    G->keepNCT();
  }

  //G->calTableName()="refant_raw.G";
  //G->storeNCT();

  G->refantmode()="flex";
  G->applyRefAnt();

  G->calTableName()="refant_flex.G";
  G->storeNCT();


  Cube<Float> blph(2,nAnt*(nAnt-1)/2,nTimes,0);

  {
    //cout << endl << endl << "Testing flex refant...................." << endl;
    Table flex("refant_flex.G");
    //cout << "Table name=" << flex.tableName() << endl;

    // Check gains
    ArrayColumn<Complex> ccol(flex,"CPARAM");
    Cube<Complex> cparam;  ccol.getColumn(cparam);
    //cout << "cparam=" << phase(cparam(Slice(0),Slice(),Slice(4,nTimes,nAnt)).nonDegenerate())*dpr << endl;

    // Check cross-hand gain
    Vector<Float> RLph(phase(cparam(Slice(0,1,1),Slice(),Slice())/cparam(Slice(1,1,1),Slice(),Slice())).nonDegenerate());
    for (Int iant=0;iant<nAnt;++iant) {
      Float mRLph=mean(RLph(Slice(iant,nTimes,5)));
      Float dRLph=mRLph-(iant-refant)*10.0/dpr;
      //cout << "RLph=" << mRLph*dpr
      //<< " dph=" << dRLph
      //	   << endl;
      EXPECT_NEAR(0.0,dRLph,5e-7);
    }

    // Check flags
    ArrayColumn<bool> fcol(flex,"FLAG");
    Cube<Bool> cflag; fcol.getColumn(cflag);
    ASSERT_EQ(uInt(5),ntrue(cflag));  // 5 are flagged
    Vector<Bool> cflagged(cflag(Slice(1,1,1),Slice(),Slice(25*nAnt+refant,5,5)).nonDegenerate());
    //cout << "cflag=" << boolalpha <<  cflagged << endl;
    ASSERT_TRUE(allTrue(cflagged));

    // Check baseline gains
    for (Int itime=0;itime<nTimes;++itime) {
      Int i0=itime*nAnt;
      Int ibl=0;
      for (Int a1=0;a1<nAnt-1;++a1) {
	Vector<Complex> g1=cparam(Slice(),Slice(),Slice(i0+a1,1,1)).nonDegenerate();
	for (Int a2=a1+1;a2<nAnt;++a2) {
	  Vector<Complex> g2=cparam(Slice(),Slice(),Slice(i0+a2,1,1)).nonDegenerate();
	  Vector<Float> blph0(phase(g1/g2)*dpr);
	  blph(Slice(),Slice(ibl,1,1),Slice(itime,1,1))=blph0;
	  //if (itime==nTimes/2-1)
	  //  cout << a1 << "-" << a2 << " " 
	  //	 << phase(g1)*dpr << " " << phase(g2)*dpr 
	  //	 << " " << blph0
	  //	 << endl;
	  ++ibl;
	}
      }
    }

    ScalarColumn<Int> a2col(flex,"ANTENNA2");
    Vector<Int> a2; a2col.getColumn(a2);
    //cout << "a2=" << a2 << endl;

  }

  //cout << "dblph=" << blph << endl;



  G->refantmode()="strict";
  G->applyRefAnt();

  G->calTableName()="refant_strict.G";
  G->storeNCT();


  {
    //cout << endl << endl << "Testing strict refant...................." << endl;
    Table strict("refant_strict.G");
    //cout << "Table name=" << strict.tableName() << endl;
    ArrayColumn<Complex> ccol(strict,"CPARAM");
    Cube<Complex> cparam;  ccol.getColumn(cparam);
    //cout << "cparam=" << phase(cparam)*dpr << endl;
    Vector<Float> RLph(phase(cparam(Slice(0,1,1),Slice(),Slice())/cparam(Slice(1,1,1),Slice(),Slice())).nonDegenerate());
    for (Int iant=0;iant<nAnt;++iant) {
      Float mRLph=mean(RLph(Slice(iant,nTimes,5)));
      Float dRLph=mRLph-(iant-refant)*10.0/dpr;
      //cout << "RLph=" << mRLph*dpr
      //<< " dph=" << dRLph
      //<< endl;
      EXPECT_NEAR(0.0,dRLph,5e-7);
    }
    ArrayColumn<bool> fcol(strict,"FLAG");
    Cube<Bool> cflag; fcol.getColumn(cflag);
    ASSERT_EQ(uInt(50),ntrue(cflag));  // 50 should be flagged
    Matrix<Bool> cflagged(cflag(Slice(),Slice(),Slice(25*nAnt,5*nAnt,1)).nonDegenerate());
    //cout << "cflag=" << boolalpha <<  cflagged << endl;
    ASSERT_TRUE(allTrue(cflagged));   // which ones are flagged

    for (Int itime=0;itime<nTimes;++itime) {
      Int i0=itime*nAnt;
      Int ibl=0;
      for (Int a1=0;a1<nAnt-1;++a1) {
	Vector<Complex> g1=cparam(Slice(),Slice(),Slice(i0+a1,1,1)).nonDegenerate();
	for (Int a2=a1+1;a2<nAnt;++a2) {
	  Vector<Complex> g2=cparam(Slice(),Slice(),Slice(i0+a2,1,1)).nonDegenerate();
	  Vector<Float> blph0(phase(g1/g2)*dpr);
	  blph(Slice(),Slice(ibl,1,1),Slice(itime,1,1))=blph(Slice(),Slice(ibl,1,1),Slice(itime,1,1)).nonDegenerate()-blph0;
	  //if (itime==nTimes/2-1)
	  //  cout << a1 << "-" << a2 << " " 
	  //	 << phase(g1)*dpr << " " << phase(g2)*dpr << " " 
	  //	 << blph0 << " "
	  //	 << endl;
	  ++ibl;
	}
      }
    }

    ScalarColumn<Int> a2col(strict,"ANTENNA2");
    Vector<Int> a2; a2col.getColumn(a2);
    //cout << "a2=" << a2 << endl;

  }

  //cout << "dblph=" << blph << endl;

  ASSERT_TRUE(allNearAbs(blph,0.0f,5e-5));

  delete G;

  Table::deleteTable("refant_flex.G");
  Table::deleteTable("refant_strict.G");

}
