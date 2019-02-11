//# tDJones: test polarization terms
//# Copyright (C) 2013
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

#include <synthesis/MeasurementComponents/XJones.h>
#include <synthesis/MeasurementComponents/StandardVisCal.h>

#include <gtest/gtest.h>

#include "VisCalTestBase_GT.h"

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;

// <summary>
// Test program for XJones-related classes
// </summary>

// Control verbosity
#define XJONES_TEST_VERBOSE false


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


class XJonesTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  XJonesTest() :
    VisCalTestBase(1,5,1,4,4,64,1,true),
    x(1,nChan,nAnt,Complex(0.0))
  {
    for (Int ich=0;ich<nChan;++ich) {
      Float a((C::pi/4)*cos(ich*C::pi/90.0));
      x(Slice(),Slice(ich,1,1),Slice())=Complex(cos(a),sin(a));
    }

    //summary("XJonesTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};
  

TEST_F(XJonesTest, XfJonesTest) {

  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  // Apply-able Xf
  XfJones Xapp(msmc);
  Xapp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) { 
    Xapp.setMeta(0,1,0.0,
		 ispw,ssvp.freqs(ispw),
		 nChan);
    Xapp.sizeApplyParCurrSpw(nChan);
    
    Xapp.setApplyParCurrSpw(x,true,false);  // correct below will corrupt
  }

  XfJones Xsol(msmc);
  Record solvePar;
  solvePar.define("table",String("test.Xf"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  Xsol.setSolve(solvePar);

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

      Cube<Complex> vC(vb2->visCube());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());

      Xapp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Xapp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)

      /*
      cerr << "vCC=" << Vector<Complex>(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      */

      Papp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Papp.corrupt2(*vb2);
      //cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;

      // Add hoc divideCorrByModel (normally done by VisEquation)
      // NB: formally, should also do weights, but doesn't matter here (point source)
      Cube<Complex> vCC(vb2->visCubeCorrected());
      Cube<Complex> vCM(vb2->visCubeModel());
      vCC/=vCM;
      vCM.set(Complex(1.0));

      /*
      cerr << "vCC=" << Vector<Complex>(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0))) 
	   << Vector<Float>(phase(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0)))*180.0f/C::pi)
	   << endl;
      cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      */

      sdbs.add(*vb2);

      //cerr  << endl << endl;
    }
  }

  // Setup meta & sizes for the solve
  Xsol.setMeta(sdbs.aggregateObsId(),
	       sdbs.aggregateScan(),
	       sdbs.aggregateTime(),
	       sdbs.aggregateSpw(),
	       sdbs.freqs(),
	       sdbs.aggregateFld());
  Xsol.sizeSolveParCurrSpw(sdbs.nChannels()); 

  // Call the specialized solver
  Xsol.selfSolveOne(sdbs);
  
  Cube<Float> soldiff=amplitude(Xsol.solveAllCPar()-x);
  
  //cerr << "Xsol.solveCPar() = " << phase(Xsol.solveAllCPar())*180/C::pi << endl;
  //cerr << "Diff = " << soldiff  << endl;
  
  ASSERT_TRUE(allNearAbs(soldiff,0.0f,1e-7)); 



}



class GlinXphJonesTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  GlinXphJonesTest() :
    VisCalTestBase(1,10,1,4,4,64,1,true,"lin"),
    x(1,nChan,nAnt,Complex(0.0))
  {
    for (Int ich=0;ich<nChan;++ich) {
      Float a((C::pi/4)*cos(ich*C::pi/90.0));
      x(Slice(),Slice(ich,1,1),Slice())=Complex(cos(a),sin(a));
    }

    //summary("GlinXphJonesTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};
  

TEST_F(GlinXphJonesTest, GlinXphfJonesTest){

  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  // Apply-able Xf
  XfJones Xapp(msmc);
  Xapp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) { 
    Xapp.setMeta(0,1,0.0,
		 ispw,ssvp.freqs(ispw),
		 nChan);
    Xapp.sizeApplyParCurrSpw(nChan);
    
    Xapp.setApplyParCurrSpw(x,true,false);  // correct below will corrupt
  }

  GlinXphfJones XYsol(msmc);
  Record solvePar;
  solvePar.define("table",String("test.Xf"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  XYsol.setSolve(solvePar);

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

      Cube<Complex> vC(vb2->visCube());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());

      Xapp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Xapp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)


      //cerr << "vCC=" << Vector<Complex>(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      //cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;


      sdbs.add(*vb2);
      //cerr << "sdb.feedPa() = " << sdbs(sdbs.nSDB()-1).feedPa() << endl;
    }
  }

  // Setup meta & sizes for the solve
  XYsol.setMeta(sdbs.aggregateObsId(),
		sdbs.aggregateScan(),
		sdbs.aggregateTime(),
		sdbs.aggregateSpw(),
		sdbs.freqs(),
		sdbs.aggregateFld());
  XYsol.sizeSolveParCurrSpw(sdbs.nChannels()); 

  // Call the specialized solver
  XYsol.selfSolveOne(sdbs);

  //cerr << "srcPolPar() = " << real(XYsol.srcPolPar()) << endl;
  ASSERT_NEAR(0.04,real(XYsol.srcPolPar()(0)),1e-3);
  ASSERT_NEAR(0.03,real(XYsol.srcPolPar()(1)),1e-3);
  
  Cube<Complex> soln(XYsol.solveAllCPar()(Slice(0),Slice(),Slice()));
  //cerr << "soln = " << phase(soln.xyPlane(0))*180/C::pi << endl;

  Cube<Float> soldiff=amplitude(soln-x);
  //cerr << "Diff = " << soldiff  << endl;
  
  ASSERT_TRUE(allNearAbs(soldiff,0.0f,1e-7)); 

}




class XfparangJonesLINTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  XfparangJonesLINTest() :
    VisCalTestBase(1,10,1,4,4,64,1,true,"lin"),
    //VisCalTestBase(1,3,1,4,4,64,1,true,"lin"),
    x(1,nChan,nAnt,Complex(0.0))
  {
    for (Int ich=0;ich<nChan;++ich) {
      //Float a((C::pi/4)*cos(ich*C::pi/90.0));         // nominal ambiguity
      Float a((C::pi/4)*cos(5*ich*C::pi/90.0)+C::pi/2);   // wrong ambiguity
      x(Slice(),Slice(ich,1,1),Slice())=Complex(cos(a),sin(a));
    }

    //summary("XfparangJonesLINTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};


TEST_F(XfparangJonesLINTest, XfparangJonesTest){

  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  // Apply-able Xf
  XfJones Xapp(msmc);
  Xapp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) { 
    Xapp.setMeta(0,1,0.0,
		 ispw,ssvp.freqs(ispw),
		 nChan);
    Xapp.sizeApplyParCurrSpw(nChan);
    
    Xapp.setApplyParCurrSpw(x,true,false);  // correct below will corrupt
  }

  XfparangJones XYsol(msmc);
  Record solvePar;
  solvePar.define("table",String("XfparangJonesLINTest.Xf"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  XYsol.setSolve(solvePar);

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

      Cube<Complex> vC(vb2->visCube());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());

      Xapp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Xapp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)


      Papp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Papp.corrupt2(*vb2);


      //cerr << "time=" << timestamp-86400.0*floor(timestamp/86400.0) << endl;
      //cerr << "vCC=" << Vector<Complex>(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      //cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;


      if (scan==2) {
	Cube<Bool> fC(vb2->flagCube());
	Cube<Complex> vCC(vb2->visCubeCorrected());
	fC(1,3,2)=true;
	vCC(1,3,2)=Complex(50,1);
      }

      sdbs.add(*vb2);
      //cerr << "sdb.feedPa() = " << sdbs(sdbs.nSDB()-1).feedPa() << endl;
    }
  }

  XYsol.createMemCalTable2();


  // Setup meta & sizes for the solve
  XYsol.setMeta(sdbs.aggregateObsId(),
		sdbs.aggregateScan(),
		sdbs.aggregateTime(),
		sdbs.aggregateSpw(),
		sdbs.freqs(),
		sdbs.aggregateFld());

  XYsol.setOrVerifyCTFrequencies(sdbs.aggregateSpw());


  XYsol.sizeSolveParCurrSpw(sdbs.nChannels()); 

  // Call the specialized solver
  XYsol.selfSolveOne(sdbs);

  cerr << "srcPolPar() = " << real(XYsol.srcPolPar()) << endl;
  ASSERT_NEAR(0.04,real(XYsol.srcPolPar()(0)),1e-3);
  ASSERT_NEAR(0.03,real(XYsol.srcPolPar()(1)),1e-3);
  
  Cube<Complex> soln(XYsol.solveAllCPar()(Slice(0),Slice(),Slice()));
  //cerr << "soln = " << phase(soln.xyPlane(0))*180/C::pi << endl;
  //cerr << "x = " << phase(x)*180/C::pi << endl;

  Cube<Float> soldiff=amplitude(soln-x);
  cerr << "maxDiff = " << max(soldiff)  << endl;
  
  ASSERT_TRUE(allNearAbs(soldiff,0.0f,1e-6)); 

  XYsol.keepNCT();

  XYsol.globalPostSolveTinker();  // writes QU to header

  XYsol.storeNCT();

  cout << "solveActionRec = " << XYsol.solveActionRec() << endl;

}



class XfparangJonesCIRCTest : public VisCalTestBase {

public:

  Cube<Complex> x;

  XfparangJonesCIRCTest() :
    VisCalTestBase(1,10,1,4,4,64,1,true,"circ"),
    x(1,nChan,nAnt,Complex(0.0))
  {
    for (Int ich=0;ich<nChan;++ich) {
      //Float a((C::pi/4)*cos(ich*C::pi/90.0));         // nominal ambiguity
      Float a((C::pi/4)*cos(5*ich*C::pi/90.0)+C::pi/2);   // wrong ambiguity
      x(Slice(),Slice(ich,1,1),Slice())=Complex(cos(a),sin(a));
    }

    //summary("XfparangJonesCIRCTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};


TEST_F(XfparangJonesCIRCTest, XfparangJonesTest){

  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  // Apply-able Xf
  XfJones Xapp(msmc);
  Xapp.setApply();

  for (Int ispw=0;ispw<nSpw;++ispw) { 
    Xapp.setMeta(0,1,0.0,
		 ispw,ssvp.freqs(ispw),
		 nChan);
    Xapp.sizeApplyParCurrSpw(nChan);
    
    Xapp.setApplyParCurrSpw(x,true,false);  // correct below will corrupt
  }

  XfparangJones XYsol(msmc);
  Record solvePar;
  solvePar.define("table",String("XfparangJonesCIRCTest.Xf"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  XYsol.setSolve(solvePar);

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

      Cube<Complex> vC(vb2->visCube()), vCM(vb2->visCubeModel());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());


      Cube<Complex> vCC(vb2->visCubeCorrected());
      Cube<Complex> vC1(vCC(Slice(1),Slice(),Slice()));
      vC1+=Complex(0.1,0.2);
      Cube<Complex> vC2(vCC(Slice(2),Slice(),Slice()));
      vC2+=Complex(0.1,-0.2);



      Xapp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Xapp.correct2(*vb2,false,false,false);  // (trial?,doWtSp?,dosync?)


      //cerr << "Forcing Q=1, U=0!!!!!!!!!!!!!" << endl;
      //vCM(Slice(1,2,1),Slice(),Slice()).set(Complex(1.0));


      Papp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      Papp.corrupt2(*vb2);


      if (scan==2) {
	Cube<Bool> fC(vb2->flagCube());
	fC(1,3,2)=true;
	//fC(2,3,2)=true;
	vCC(1,3,2)=Complex(500,1);
      }


      //      cerr << "C="<< vCC(Slice(),Slice(0),Slice(0)).nonDegenerate() << " M=" << vCM(Slice(),Slice(0),Slice(0)).nonDegenerate() << endl;



      //cerr << "time=" << timestamp-86400.0*floor(timestamp/86400.0) << endl;
      //cerr << "vCC=" << Vector<Complex>(vb2->visCubeCorrected()(Slice(),Slice(nChan/2),Slice(0))) << endl;
      //cerr << "vCM=" << Vector<Complex>(vb2->visCubeModel()(Slice(),Slice(nChan/2),Slice(0))) << endl;


      sdbs.add(*vb2);
      //cerr << "sdb.feedPa() = " << sdbs(sdbs.nSDB()-1).feedPa() << endl;
    }
  }

  XYsol.createMemCalTable2();


  // Setup meta & sizes for the solve
  XYsol.setMeta(sdbs.aggregateObsId(),
		sdbs.aggregateScan(),
		sdbs.aggregateTime(),
		sdbs.aggregateSpw(),
		sdbs.freqs(),
		sdbs.aggregateFld());

  XYsol.setOrVerifyCTFrequencies(sdbs.aggregateSpw());


  XYsol.sizeSolveParCurrSpw(sdbs.nChannels()); 

  // Call the specialized solver
  XYsol.selfSolveOne(sdbs);

  cerr << "srcPolPar() = " << real(XYsol.srcPolPar()) << endl;
  ASSERT_NEAR(0.05,real(XYsol.srcPolPar()(0)),1e-5);
  ASSERT_NEAR(0.00,real(XYsol.srcPolPar()(1)),1e-5);
  
  Cube<Complex> soln(XYsol.solveAllCPar()(Slice(0),Slice(),Slice()));
  //cerr << "soln = " << phase(soln.xyPlane(0))*180/C::pi << endl;
  //cerr << "x = " << phase(x)*180/C::pi << endl;

  Cube<Float> soldiff=amplitude(soln-x);
  cerr << "maxDiff = " << max(soldiff)  << endl;
  
  ASSERT_TRUE(allNearAbs(soldiff,0.0f,1e-6)); 

  XYsol.keepNCT();
  XYsol.globalPostSolveTinker();  // writes QU to header

  XYsol.storeNCT();

  Record sAR=XYsol.solveActionRec();

  cout << "solveActionRec = " << sAR  << endl;
  cout << "sAR.nfields() = " << sAR.nfields() << endl;


}


class PosAngJonesLINTest : public VisCalTestBase {

public:

  PosAngJonesLINTest() :
    VisCalTestBase(1,5,1,4,4,64,1,true,"lin")
  {
    //summary("XJonesTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};

TEST_F(PosAngJonesLINTest, PATest1){

  {
  PosAngJones PA(msmc);
  PA.setApply();
  PA.setMeta(0,0,0.0,
	     0,ssvp.freqs(0),
	     0);
  PA.sizeApplyParCurrSpw(nChan);
  //PA.setDefApplyParCurrSpw(true,true);  // sync, w/ doInv=T
  PA.state();

  }


  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  
  PosAngJones PAsol(msmc);
  Record solvePar;
  solvePar.define("table",String("PosAngJonesLINTest.PA"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  PAsol.setSolve(solvePar);

  PAsol.createMemCalTable2();
      

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

      Cube<Complex> vC(vb2->visCube()), vCM(vb2->visCubeModel());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());
      Cube<Complex> vCC(vb2->visCubeCorrected());

      cerr << "Scan=" << scan << " Bln=" << a1(0) << "-" << a2(0) 
	   << " C=" << vCC(Slice(),Slice(0),Slice(0)).nonDegenerate()
	   << " M=" << vCM(Slice(),Slice(0),Slice(0)).nonDegenerate()
	   << endl;

      Complex oneI(0,1);
      Complex RL((vCC(0,0,0)-vCC(3,0,0))+oneI*(vCC(1,0,0)+vCC(2,0,0)) );
      Float dang(arg(RL)*90/C::pi);
      Complex RLm((vCM(0,0,0)-vCM(3,0,0))+oneI*(vCM(1,0,0)+vCM(2,0,0)) );
      Float mang(arg(RLm)*90/C::pi);
      Float diff(arg(RL/RLm)*90/C::pi);

      cerr << dang << " - " << mang << " = " << diff << endl;

      //cerr << "Forcing Q=1, U=0!!!!!!!!!!!!!" << endl;
      //vCM(Slice(1,2,1),Slice(),Slice()).set(Complex(1.0));

     
      Papp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      //Papp.corrupt2(*vb2);
      

      
      cerr << "C="<< vCC(Slice(),Slice(0),Slice(0)).nonDegenerate() << " M=" << vCM(Slice(),Slice(0),Slice(0)).nonDegenerate() << endl;

      SDBList sdbs;
      sdbs.add(*vb2);

      // Setup meta & sizes for the solve
      PAsol.setMeta(sdbs.aggregateObsId(),
		    sdbs.aggregateScan(),
		    sdbs.aggregateTime(),
		    sdbs.aggregateSpw(),
		    sdbs.freqs(),
		    sdbs.aggregateFld());


      PAsol.setOrVerifyCTFrequencies(sdbs.aggregateSpw());


      PAsol.sizeSolveParCurrSpw(sdbs.nChannels()); 


      // Call the specialized solver
      PAsol.selfSolveOne(sdbs);

      Float soln(PAsol.solveAllRPar()(0,nChan/2,0));

      cerr << "soln = " << soln*180/C::pi << endl << endl;

      PAsol.keepNCT();

    }
  }
  PAsol.storeNCT();


}


class PosAngJonesCIRCTest : public VisCalTestBase {

public:

  PosAngJonesCIRCTest() :
    VisCalTestBase(1,5,1,4,4,64,1,true,"circ")
  {
    //summary("XJonesTest");
    //cerr << "x = " << phase(x)*180/C::pi << endl;
  }


};

TEST_F(PosAngJonesCIRCTest, PATest1){

  {
  PosAngJones PA(msmc);
  PA.setApply();
  PA.setMeta(0,0,0.0,
	     0,ssvp.freqs(0),
	     0);
  PA.sizeApplyParCurrSpw(nChan);
  //PA.setDefApplyParCurrSpw(true,true);  // sync, w/ doInv=T
  PA.state();

  }


  // Apply-able parang
  PJones Papp(msmc);
  Papp.setApply();

  
  PosAngJones PAsol(msmc);
  Record solvePar;
  solvePar.define("table",String("PosAngJonesCIRCTest.PA"));
  solvePar.define("solint",String("inf"));
  solvePar.define("combine",String(""));
  Vector<Int> refant(1,0); solvePar.define("refant",refant);
  PAsol.setSolve(solvePar);

  PAsol.createMemCalTable2();
  
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

      Cube<Complex> vC(vb2->visCube()), vCM(vb2->visCubeModel());
      vb2->setVisCubeCorrected(vC);
      vb2->setFlagCube(vb2->flagCube());
      Cube<Complex> vCC(vb2->visCubeCorrected());

      cerr << "Scan=" << scan << " Bln=" << a1(0) << "-" << a2(0) 
	   << " C=" << vCC(Slice(),Slice(0),Slice(0)).nonDegenerate()
	   << " M=" << vCM(Slice(),Slice(0),Slice(0)).nonDegenerate()
	   << endl;

      Float dang(arg(vCC(1,0,0))*90/C::pi), mang(arg(vCM(1,0,0))*90/C::pi);
      Float diff(arg(vCC(1,0,0)/vCM(1,0,0))*90/C::pi);

      cerr << dang << " - " << mang << " = " << diff << endl;


      //cerr << "Forcing Q=1, U=0!!!!!!!!!!!!!" << endl;
      //vCM(Slice(1,2,1),Slice(),Slice()).set(Complex(1.0));

     
      Papp.setMeta(obsid,scan,timestamp,
		   ispw,freqs,
		   fldid);
      //Papp.corrupt2(*vb2);
      

      
      cerr << "C="<< vCC(Slice(),Slice(0),Slice(0)).nonDegenerate() << " M=" << vCM(Slice(),Slice(0),Slice(0)).nonDegenerate() << endl;

      SDBList sdbs;
      sdbs.add(*vb2);

      // Setup meta & sizes for the solve
      PAsol.setMeta(sdbs.aggregateObsId(),
		    sdbs.aggregateScan(),
		    sdbs.aggregateTime(),
		    sdbs.aggregateSpw(),
		    sdbs.freqs(),
		    sdbs.aggregateFld());


      PAsol.setOrVerifyCTFrequencies(sdbs.aggregateSpw());


      PAsol.sizeSolveParCurrSpw(sdbs.nChannels()); 


      // Call the specialized solver
      PAsol.selfSolveOne(sdbs);

      Float soln(PAsol.solveAllRPar()(0,nChan/2,0));

      cerr << "soln = " << soln*180/C::pi << endl << endl;

      PAsol.keepNCT();

    }
  }
  PAsol.storeNCT();

}
