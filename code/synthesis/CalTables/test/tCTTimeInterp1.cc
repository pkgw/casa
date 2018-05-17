//# tCTTimeInterp1.cc: Test program for CTTimeInterp1 class
//# Copyright (C) 2011
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
//# $Id: tNewCalTable.cc 15602 2011-07-14 00:03:34Z tak.tsutsumi $

#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTTimeInterp1.h>
#include <synthesis/CalTables/CTPatchedInterp.h>
#include <scimath/Mathematics/InterpolateArray1D.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/BasicMath/Math.h>
#include <casa/namespace.h>

// <summary>
// Test program for CTInterp class.
// </summary>

// Control verbosity
#define CTTIMEINTERP1TEST_VERBOSE false

using namespace casa;

void doTest1 (Bool verbose=false) {

  {  
  // The correctness tolerance
  Double tol(2e-7);

  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(1), nAnt(1), nSpw(1), nObs(1),nScan(1), nTime(100);
  Vector<Int> nChan(nSpw,1);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(60.0);
  Bool disk(verbose);
  NewCalTable tnct("tCTTimeInterp1_test1.ct","Complex",
		   nObs,nScan,nTime,
		   nAnt,nSpw,nChan,
		   nFld,
                   refTime,tint,disk,false); 

  // some sanity checks on the test NewCalTable
  if (verbose) cout << "Table::Type: " << tnct.tableType() 
                    << " (should be " << Table::Memory << ")"
                    << endl;
  AlwaysAssert( (tnct.tableType() == Table::Memory), AipsError);
  if (verbose) cout << "nrow = " << tnct.nrow() 
                    << " (should be " << nObs*nScan*nTime*nSpw*nAnt << ")"
                    << endl;
  AlwaysAssert( (tnct.nrow()==nObs*nScan*nTime*nSpw*nAnt), AipsError);

  // Make a CTInterp
  Matrix<Float> result(2,1);
  result(0,0)=1.0; result(1,0)=0.0;
  Matrix<Bool> rflag(1,1); rflag.set(false);
  String itype("nearest");
  CTTimeInterp1 sci(tnct,itype,result,rflag);

  if (verbose)
    sci.state();

  if (verbose) cout << "sci.timeType() = " << sci.timeType() << endl;
  AlwaysAssert( (sci.timeType()==itype), AipsError);

  {
    // Testing nearest in various ways
    Vector<Double> timelist(6), ntimelist(6);
    timelist(0)=refTime-1.e5*tint; ntimelist(0)=refTime;   // an early point
    timelist(1)=refTime+tint;      ntimelist(1)=refTime+tint;  // an exact point
    timelist(2)=refTime+2.2*tint;  ntimelist(2)=refTime+2*tint;  // nearest left
    timelist(3)=refTime+2.9*tint;  ntimelist(3)=refTime+3*tint;  // nearest right
    timelist(4)=refTime+3.1*tint;  ntimelist(4)=refTime+3*tint;  // new time, no change
    timelist(5)=refTime+1.e5*tint; ntimelist(5)=refTime+(nTime-1)*tint;  // a late point
    
    if (verbose) cout << "Testing 'nearest' (tol="<<tol<<")..." << endl;
    for (uInt itime=0;itime<timelist.nelements();++itime) {
      // Call time interpolation
      sci.interpolate(timelist(itime));
      Complex cfval=NewCalTable::NCTtestvalueC(0,0,0,ntimelist(itime),refTime,tint);
      RIorAPArray f(result);
      //    Matrix<Complex> c(f.c());
      Complex rc=Matrix<Complex>(f.c())(0,0);
      Float diff=abs(rc-cfval);
      //if (verbose)
	cout << "t="<<timelist(itime)<<" (dt="<<(timelist(itime)-refTime)/tint<<")"
	     << " result = " << rc << " cf=" << cfval << "  diff="<< diff << " near=" << boolalpha << nearAbs(diff,0.0f,tol) << endl;
      AlwaysAssert(nearAbs(diff,0.0f,tol),AipsError);
    }
  }

  // Reset to linear
  itype="linear";
  sci.setInterpType(itype);

  if (verbose)
    sci.state();

  if (verbose) cout << "sci.timeType() = " << sci.timeType() << endl;
  AlwaysAssert( (sci.timeType()==itype), AipsError);



  {
    // Testing linear in various ways
    Vector<Double> timelist(6), ntimelist(6);
    timelist(0)=refTime-1.e5*tint; ntimelist(0)=refTime;   // a way early point
    timelist(1)=refTime+tint;      ntimelist(1)=timelist(1);  // an exact point
    timelist(2)=refTime+2.2*tint;  ntimelist(2)=timelist(2);  // nearest to left
    timelist(3)=refTime+2.9*tint;  ntimelist(3)=timelist(3);  // nearest to right (same interval)
    timelist(4)=refTime+3.1*tint;  ntimelist(4)=timelist(4);  // new time, new interval
    timelist(5)=refTime+1.e5*tint; ntimelist(5)=refTime+(nTime-1)*tint;  // a way late point
    
    if (verbose) cout << "Testing 'linear' (tol="<<tol<<")..." << endl;
    for (uInt itime=0;itime<timelist.nelements();++itime) {
      // Call time interpolation
      sci.interpolate(timelist(itime));
      Complex cfval=NewCalTable::NCTtestvalueC(0,0,0,ntimelist(itime),refTime,tint);
      RIorAPArray f(result);
      Complex rc=Matrix<Complex>(f.c())(0,0);
      Float diff=abs(rc-cfval);
      //if (verbose)
	cout << "t="<<timelist(itime)<<" (dt="<<(timelist(itime)-refTime)/tint<<")"
	     << " result = " << rc << " cf=" << cfval << "  diff="<<diff<< " near=" << boolalpha << nearAbs(diff,0.0f,tol) << endl;
      AlwaysAssert(nearAbs(diff,0.0f,tol),AipsError);
    }
  }

  // Reset to cubic
  itype="cubic";
  sci.setInterpType(itype);

  if (verbose)
    sci.state();

  if (verbose) cout << "sci.timeType() = " << sci.timeType() << endl;
  AlwaysAssert( (sci.timeType()==itype), AipsError);



  {
    // Testing cubic in various ways
    Vector<Double> timelist(6), ntimelist(6);
    timelist(0)=refTime-1.e5*tint; ntimelist(0)=refTime;   // a way early point
    timelist(1)=refTime+tint;      ntimelist(1)=timelist(1);  // an exact point
    timelist(2)=refTime+2.2*tint;  ntimelist(2)=timelist(2);  // nearest to left
    timelist(3)=refTime+2.9*tint;  ntimelist(3)=timelist(3);  // nearest to right (same interval)
    timelist(4)=refTime+3.1*tint;  ntimelist(4)=timelist(4);  // new time, new interval
    timelist(5)=refTime+1.e5*tint; ntimelist(5)=refTime+(nTime-1)*tint;  // a way late point
    
    if (verbose) cout << "Testing 'cubic' (tol="<<tol<<")..." << endl;
    for (uInt itime=0;itime<timelist.nelements();++itime) {
      // Call time interpolation
      sci.interpolate(timelist(itime));
      Complex cfval=NewCalTable::NCTtestvalueC(0,0,0,ntimelist(itime),refTime,tint);
      RIorAPArray f(result);
      Complex rc=Matrix<Complex>(f.c())(0,0);
      Float diff=abs(rc-cfval);
      //if (verbose)
	cout << "t="<<timelist(itime)<<" (dt="<<(timelist(itime)-refTime)/tint<<")"
	     << " result = " << rc << " cf=" << cfval << "  diff="<<diff<< " near=" << boolalpha << nearAbs(diff,0.0f,tol) << endl;
      AlwaysAssert(nearAbs(diff,0.0f,tol),AipsError);
    }
  }

  }


}

/* Experimenting with pointers to factory methods

class A {
public:
  
  A() {}
  A(Int i) { cout << "A ctor w/ i=" << i << endl; }
  virtual ~A() {}

  static A* factory(Int i) { return new A(i); }
  
  virtual void name() {  cout << "I'm a A" << endl; }
  
  
};

class B: public A {
public:

  B():A() {}
  B(Int i):A() { cout << "B ctor w/ i=" << i << endl; }

  static A* factory(Int i) { return new B(i); }
  
  virtual void name() {  cout << "I'm a B" << endl; }
  
  virtual ~B() {}
  
};

typedef A* (*AFactory)(Int i);

class X {
public:
  X(): a_(0) { cout << "X ctor" << endl;}
  virtual ~X() { if (a_) delete a_; }
  virtual void init(Int i)  { cout << "X::init" << endl; a_=(*af())(i); a_->name(); }

  A* a_;

private:
  virtual AFactory af() { cout << "Using A::factory" << endl; return &A::factory; }
};


class Y: public X {
public:
  Y(): X() {cout << "Y ctor" << endl;}
  virtual ~Y() {}

private:
  virtual AFactory af() { cout << "Using B::factory" << endl; return &B::factory; }
};




void doTestFactory (AFactory f, Int i) {
  A* a=(*f)(i);
  a->name();
  delete a;
}

void doTest2() {

{
  AFactory p=&B::factory;
  A* a = (*p)(3);
  a->name();
  delete a;
}

 cout <<endl<< "***********" << endl;

 doTestFactory(&A::factory,4);
 cout << "***********" << endl;
 doTestFactory(&B::factory,5);

 cout << endl << "*********************" << endl;

 X x;
 x.init(13);
 cout << "***********" << endl;
 Y y;
 y.init(21);


}

//*/

int main ()
{
  try {

    doTest1(CTTIMEINTERP1TEST_VERBOSE);

    // doTest2();

  } catch (AipsError x) {
    cout << "Unexpected exception: " << x.getMesg() << endl;
    exit(1);
  } catch (...) {
    cout << "Unexpected unknown exception" << endl;
    exit(1);
  }
  cout << "OK" << endl;
  exit(0);
};
