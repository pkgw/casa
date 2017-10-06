//# tCTSelection.cc: Test program for selecting on a NewCalTable
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

#include <casa/namespace.h>
#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/CalTables/CTInterface.h>
#include <synthesis/CalTables/CTSelection.h>
#include <ms/MSSel/MSSelectionTools.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/BasicMath/Math.h>

// <summary>
// Test program for CTSelection class.
// </summary>

using namespace casacore;
using namespace casa;

// Control verbosity
#define CTSELECTION_VERBOSE true


void doTest1 (Bool verbose=false) {

  cout << "****----doTest1()----****" << endl;
  
  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(20), nAnt(10), nSpw(4), nObs(1), nScan(20),nTime(10);
  Vector<Int> nChan(nSpw,32);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(60.0);

  Bool disk(verbose);
  NewCalTable tnct("tCTSelection1.ct","Complex",
		   nObs,nScan,nTime,
		   nAnt,nSpw,nChan,
		   nFld,
		   refTime,tint,disk,false);

  if (verbose) {
    cout << "Wrote NewCalTable out to tCTSelection1.ct" << endl;
    cout << "Reference antenna id is 0" << endl;
  }

  // some sanity checks on the test NewCalTable
  AlwaysAssert( (tnct.tableType() == Table::Memory), AipsError);
  if (verbose) cout << "Table::Type: " << tnct.tableType() 
		    << "  (should be 1)"
		    << endl;
  AlwaysAssert( (tnct.nrow()==nObs*nScan*nTime*nSpw*nAnt), AipsError);
  if (verbose) cout << "nrow = " << tnct.nrow() 
		    << "  (should be " << nObs*nScan*nTime*nSpw*nAnt << ")"
		    << endl;

  ROCTColumns ctc(tnct);

  // Extract some field names to select
  Vector<String> fldnames=ctc.field().name().getColumn();
  Vector<Int> fldids(3);
  fldids(0)=3; fldids(1)=10; fldids(2)=17;
  String fieldsel("");
  for (uInt i=0;i<fldids.nelements();++i) {
    if (i>0) fieldsel+=",";
    fieldsel+=fldnames(fldids(i));
  }

  // Some spws to select
  Vector<Int> spwids(2);
  spwids(0)=1; spwids(1)=3;
  String spwsel("");
  for (uInt i=0;i<spwids.nelements();++i) {
    if (i>0) spwsel+=",";
    spwsel+=String::toString(spwids(i));
  }

  // Extract some antenna names to select; antenna 0 is reference antenna
  Vector<String> antnames=ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  // select reference antenna 0 to check that not all baselines/rows 
  // were selected
  antids(0)=0; antids(1)=3; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0;i<antids.nelements();++i) {
    if (i>0) antsel+=",";
    antsel+=antnames(antids(i));
  }

  if (verbose) {
    cout << "Selection strings:" << endl;
    cout << "  fieldsel = " << fieldsel << endl;
    cout << "  spwsel = " << spwsel << endl;
    cout << "  antsel = " << antsel << endl;
  }

  NewCalTable selnct(tnct);
  CTInterface cti(tnct);
  CTSelection cts;
  cts.setFieldExpr(fieldsel);
  cts.setSpwExpr(spwsel);
  cts.setAntennaExpr(antsel);

  TableExprNode ten=cts.toTableExprNode(&cti);

  if (verbose) {
    cout << "CTSelection index results (get*List): " << endl;
    cout << "  Field list: " << cts.getFieldList() << endl; // OK?
    cout << "  Antenna list: " << cts.getAntenna1List() << endl;
    cout << "  Spw list: " << cts.getSpwList() << endl;
  }
  AlwaysAssert( allEQ(cts.getFieldList(),fldids), AipsError );
  AlwaysAssert( allEQ(cts.getAntenna1List(),antids), AipsError );
  AlwaysAssert( allEQ(cts.getSpwList(),spwids), AipsError );
  if (verbose) 
    cout << "CTSelection index results (get*List): all OK " << endl;

  getSelectedTable(selnct,tnct,ten,"");

  if (verbose) {
      cout << "Selected table: nrow=" << selnct.nrow()<< " ";
      cout << "(should be " << antids.nelements()*nTime*spwids.nelements()*fldids.nelements() << ")" << endl;
  }
  AlwaysAssert( (selnct.nrow()==antids.nelements()*nTime*spwids.nelements()*fldids.nelements()), AipsError);

  ROCTMainColumns ctmc(selnct);
  Vector<Int> selfieldcol, selantcol, selant2col, selspwcol;
  ctmc.fieldId().getColumn(selfieldcol);
  ctmc.antenna1().getColumn(selantcol);
  ctmc.antenna2().getColumn(selant2col);
  ctmc.spwId().getColumn(selspwcol);

  Vector<Bool> fldok(selfieldcol.nelements(),false);
  for (uInt i=0;i<fldids.nelements();++i) {
    fldok|=(selfieldcol==fldids(i));
  }
  Vector<Bool> spwok(selspwcol.nelements(),false);
  for (uInt i=0;i<spwids.nelements();++i) {
    spwok|=(selspwcol==spwids(i));
  }
  Vector<Bool> antok(selantcol.nelements(),false);
  for (uInt i=0;i<antids.nelements();++i) {
    antok|=(selantcol==antids(i));
  }

  Bool allfldok=allEQ(fldok,true);
  Bool allspwok=allEQ(spwok,true);
  Bool allantok=allEQ(antok,true);
  if (verbose) {
    cout << boolalpha;
    cout << "allEQ(fldok,true) = " << allfldok  << endl;
    cout << "allEQ(spwok,true) = " << allspwok  << endl;
    cout << "allEQ(antok,true) = " << allantok  << endl;
  }
  AlwaysAssert( allfldok , AipsError );
  AlwaysAssert( allspwok , AipsError );
  AlwaysAssert( allantok , AipsError );

  cout << "\n**** reset() with different antenna selection ****" 
       << endl;
  antsel = "0~3;!1&0";
  antids(0)=0; antids(1)=1; antids(2)=2; antids(3)=3; antids(4)=-1; 
  if (verbose) {
    cout << "Selection strings:" << endl;
    cout << "  fieldsel = " << fieldsel << endl;
    cout << "  spwsel = " << spwsel << endl;
    cout << "  antsel = " << antsel << endl;
  }
  cts.reset(cti, casacore::MSSelection::PARSE_NOW, "", antsel,
    fieldsel, spwsel, "", "", "", "");
  if (verbose) {
    cout << "CTSelection index results (get*List): " << endl;
    cout << "  Field list: " << cts.getFieldList() << endl;
    cout << "  Antenna list: " << cts.getAntenna1List() << endl;
    cout << "  Spw list: " << cts.getSpwList() << endl;
  }
  AlwaysAssert( allEQ(cts.getAntenna1List(),antids), AipsError );
  AlwaysAssert( allEQ(cts.getFieldList(),fldids), AipsError );
  AlwaysAssert( allEQ(cts.getSpwList(),spwids), AipsError );
  if (verbose) 
    cout << "CTSelection parsing results (get*List): all OK" << endl;

  ten = cts.getTEN();
  getSelectedTable(selnct,tnct,ten,"");
  ROCTMainColumns ctmain(selnct);
  Vector<Int> selfieldcol2, selantcol2, selspwcol2;
  ctmain.fieldId().getColumn(selfieldcol2);
  ctmain.antenna1().getColumn(selantcol2);
  ctmain.spwId().getColumn(selspwcol2);
  /*
  if (verbose) {
    cout << "ctmc.fieldId().getColumn()  = " << selfieldcol2 << endl;
    cout << "ctmc.antenna1().getColumn() = " << selantcol2 << endl;
    cout << "ctmc.spwId().getColumn()    = " << selspwcol2 << endl;
  }
  */

  // antenna 1 is selected and negated
  antids.resize(3);
  antids(0)=0; antids(1)=2; antids(2)=3; 
  if (verbose)
    {
      cout << "Selected table: nrow=" << selnct.nrow()<< " ";
      cout << "(should be " << antids.size() * nTime * spwids.nelements() * 
          fldids.nelements() << ")" << endl;
    }
  AlwaysAssert( (selnct.nrow() == antids.size() * nTime * spwids.nelements() *
    fldids.nelements()), AipsError);
}

void doTest2 (Bool verbose=false) {

  cout << "\n****----doTest2()----****" << endl;
  
  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(1), nAnt(10), nSpw(4), nObs(3), nScan=6, nTime(6);
  Vector<Int> nChan(nSpw,1);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(10.0);

  Bool disk(verbose);
  NewCalTable tnct("tCTSelection2.ct","Complex",
		   nObs,nScan,nTime,
		   nAnt,nSpw,nChan,
		   nFld,
		   refTime,tint,disk,false);

  if (verbose) {
    cout << "Wrote test NewCalTable out to tCTSelection2.ct" << endl;
    cout << "Reference antenna id is 0" << endl;
  }

  // some sanity checks on the test NewCalTable
  AlwaysAssert( (tnct.tableType() == Table::Memory), AipsError);
  if (verbose) cout << "Table::Type: " << tnct.tableType() 
		    << "  (should be 1)"
		    << endl;
  AlwaysAssert( (tnct.nrow()==nObs*nScan*nTime*nSpw*nAnt), AipsError);
  if (verbose) cout << "nrow = " << tnct.nrow() 
		    << "  (should be " << nObs*nScan*nTime*nSpw*nAnt << ")"
		    << endl;

  ROCTColumns ctc(tnct);

  // And ObsID selection:
  Vector<Int> obsids(1);
  obsids(0)=1;
  String obssel=String::toString(obsids(0));

  // A scan selection:
  Vector<Int> scans(1);
  scans(0)=7;
  String scansel=String::toString(scans(0));

  // A time selection expression:
  //  this should pick 3 timestamps in obsid=1, scan=7
  String timesel("2012/01/06/12:06:15.0~2012/01/06/12:06:45");
  Vector<Double> timebounds(3);
  timebounds(0)=375.0+refTime;
  timebounds(1)=405.0+refTime;
  timebounds(2)=0.0; // integration time

  // Some spws to select
  Vector<Int> spwids(2);  
  spwids(0)=1; spwids(1)=3;
  String spwsel("");
  for (uInt i=0;i<spwids.nelements();++i) {
    if (i>0) spwsel+=",";
    spwsel+=String::toString(spwids(i));
  }

  // Extract some antenna names to select
  Vector<String> antnames=ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  antids(0)=2; antids(1)=3; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0;i<antids.nelements();++i) {
    if (i>0) antsel+=",";
    antsel+=antnames(antids(i));
  }

  if (verbose) {
    cout << "Selection strings:" << endl;
    cout << "  obssel  = " << obssel << endl;
    cout << "  scansel = " << scansel << endl;
    cout << "  timesel = " << timesel << endl;
    cout << "  spwsel  = " << spwsel << endl;
    cout << "  antsel  = " << antsel << endl;
  }

  NewCalTable selnct(tnct);
  CTInterface cti(tnct);
  CTSelection cts;
  cts.setObservationExpr(obssel);
  cts.setScanExpr(scansel);
  cts.setTimeExpr(timesel);
  cts.setSpwExpr(spwsel);
  cts.setAntennaExpr(antsel);

  TableExprNode ten = cts.toTableExprNode(&cti);

  if (verbose) {
    cout << "CTSelection parsing results (get*List): " << endl;
    cout << "  Obs list: " << cts.getObservationList() << endl;
    cout << "  Scan list: " << cts.getScanList() << endl;
    cout.precision(15);
	casacore::Vector<Double> parsedTimeBounds=Vector<Double>(cts.getTimeList()); 
    cout << "  Time bounds: " << parsedTimeBounds  
	     << " [" << parsedTimeBounds(Slice(0,2,1))-refTime << "]" << endl;
	cout << "    Expecting: " << timebounds << " ["  
		 << timebounds(Slice(0,2,1))-refTime << "]" << endl;
    cout << "    (time diff = " << parsedTimeBounds-timebounds << ")" 
		 << endl;
    cout << "  Spw list: " << cts.getSpwList() << endl;
  }
  AlwaysAssert( allEQ(cts.getObservationList(),obsids), AipsError );  
  AlwaysAssert( allEQ(cts.getScanList(),scans), AipsError );
  AlwaysAssert( allNear(Vector<Double>(cts.getTimeList()),timebounds,1.0e-9), AipsError );
  AlwaysAssert( allEQ(cts.getSpwList(),spwids), AipsError );
  AlwaysAssert( allEQ(cts.getAntenna1List(),antids), AipsError );
  AlwaysAssert( allEQ(cts.getObservationList(),obsids), AipsError );  
  if (verbose)
    cout << "CTSelection parsing results (get*List): all OK" << endl;

  getSelectedTable(selnct,tnct,ten,"");

  if (verbose)
    cout << "Selected table: nrow=" << selnct.nrow() 
	 << " (should be " << 3*antids.nelements()*spwids.nelements() << ")"
	 << endl;
  AlwaysAssert( (selnct.nrow() == 3 * antids.nelements() * spwids.nelements()),
    AipsError);

  ROCTMainColumns ctmc(selnct);
  Vector<Int> selantcol,selspwcol,selobscol,selscancol;
  Vector<Double> seltimecol;
  ctmc.obsId().getColumn(selobscol);
  ctmc.scanNo().getColumn(selscancol);
  ctmc.time().getColumn(seltimecol);
  ctmc.spwId().getColumn(selspwcol);
  ctmc.antenna1().getColumn(selantcol);

  Vector<Bool> obsok(selobscol.nelements(),false);
  for (uInt i=0;i<obsids.nelements();++i) 
    obsok|=(selobscol==obsids(i));

  Vector<Bool> scanok(selscancol.nelements(),false);
  for (uInt i=0;i<scans.nelements();++i) 
    scanok|=(selscancol==scans(i));

  Vector<Bool> timeok(seltimecol.nelements(),false);
  timeok|=(seltimecol>timebounds(0));
  timeok|=(seltimecol<timebounds(1));

  Vector<Bool> spwok(selspwcol.nelements(),false);
  for (uInt i=0;i<spwids.nelements();++i) 
    spwok|=(selspwcol==spwids(i));

  Vector<Bool> antok(selantcol.nelements(),false);
  for (uInt i=0;i<antids.nelements();++i) 
    antok|=(selantcol==antids(i));
  
  Bool allobsok=allEQ(obsok,true);
  Bool allscanok=allEQ(scanok,true);
  Bool alltimeok=allEQ(timeok,true);
  Bool allspwok=allEQ(spwok,true);
  Bool allantok=allEQ(antok,true);
  if (verbose) {
    cout << boolalpha;
    cout << "allEQ(obsok,true) = " << allobsok  << endl;
    cout << "allEQ(scanok,true) = " << allscanok  << endl;
    cout << "allEQ(timeok,true) = " << alltimeok  << endl;
    cout << "allEQ(spwok,true) = " << allspwok  << endl;
    cout << "allEQ(antok,true) = " << allantok  << endl;
  }
  AlwaysAssert( allobsok , AipsError );
  AlwaysAssert( allscanok , AipsError );
  AlwaysAssert( alltimeok , AipsError );
  AlwaysAssert( allspwok , AipsError );
  AlwaysAssert( allantok , AipsError );

  cout << "\n**** repeat test2 using ctor with selection expressions ****" 
       << endl;
  CTSelection* newcts = new CTSelection(tnct, casacore::MSSelection::PARSE_NOW,
    timesel, antsel, "", spwsel, "", scansel, "", obssel);
  AlwaysAssert( allEQ(newcts->getObservationList(),obsids), AipsError );  
  AlwaysAssert( allEQ(newcts->getScanList(),scans), AipsError );
  AlwaysAssert( allNear(Vector<Double>(newcts->getTimeList()),timebounds,1.0e-9), 
    AipsError );
  AlwaysAssert( allEQ(newcts->getSpwList(),spwids), AipsError );
  AlwaysAssert( allEQ(newcts->getAntenna1List(),antids), AipsError );
  AlwaysAssert( allEQ(newcts->getObservationList(),obsids), AipsError );  
  if (verbose) 
    cout << "CTSelection parsing results (get*List) all OK" << endl;

  ten = newcts->getTEN();
  getSelectedTable(selnct,tnct,ten,"");
  if (verbose)
    cout << "Selected table: nrow=" << selnct.nrow() 
	     << " (should be " << 3 * antids.nelements() * spwids.nelements() << ")"
	     << endl;
  AlwaysAssert( (selnct.nrow() == 3 * antids.nelements() * spwids.nelements()),
    AipsError);

  ROCTMainColumns ctmain(selnct);
  Vector<Int> selantcol2, selspwcol2, selobscol2, selscancol2;
  Vector<Double> seltimecol2;
  ctmain.antenna1().getColumn(selantcol2);
  ctmain.spwId().getColumn(selspwcol2);
  ctmain.obsId().getColumn(selobscol2);
  ctmain.scanNo().getColumn(selscancol2);
  ctmain.time().getColumn(seltimecol2);

  obsok = false;
  for (uInt i=0;i<obsids.nelements();++i) 
    obsok |= (selobscol2==obsids(i));

  scanok = false;
  for (uInt i=0;i<scans.nelements();++i) 
    scanok |= (selscancol2==scans(i));

  timeok = false;
  timeok |= (seltimecol2>timebounds(0));
  timeok |= (seltimecol2<timebounds(1));

  spwok = false;
  for (uInt i=0;i<spwids.nelements();++i)
    spwok |= (selspwcol2==spwids(i));

  antok = false;
  for (uInt i=0;i<antids.nelements();++i)
    antok |= (selantcol2==antids(i));
  
  AlwaysAssert( allEQ(obsok,true), AipsError );
  AlwaysAssert( allEQ(scanok,true), AipsError );
  AlwaysAssert( allEQ(timeok,true), AipsError );
  AlwaysAssert( allEQ(spwok,true), AipsError );
  AlwaysAssert( allEQ(antok,true), AipsError );
  if (verbose)
    cout << "CTSelection table contents results all OK" << endl;
}

int main ()
{
  try {

    doTest1(CTSELECTION_VERBOSE);
    doTest2(CTSELECTION_VERBOSE);

  } catch (AipsError x) {
    cout << "Unexpected exception: " << x.getMesg() << endl;
    exit(1);
  } catch (...) {
    cout << "Unexpected unknown exception" << endl;
    exit(1);
  }
  cout << "\nAll OK" << endl;
  exit(0);
};

