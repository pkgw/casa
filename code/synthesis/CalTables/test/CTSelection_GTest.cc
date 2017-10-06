//# CTSelection_GTest.cc: Google Test program for selecting on a NewCalTable
//# Copyright (C) 2017
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

#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/CalTables/CTInterface.h>
#include <synthesis/CalTables/CTSelection.h>
#include <ms/MSSel/MSSelectionTools.h>
#include <casa/Arrays/ArrayLogical.h>

#include <gtest/gtest.h>

// <summary>
// Google Test program for CTSelection class.
// Adapted from tCTSelection.cc legacy test.
// </summary>

using namespace casacore;
using namespace casa;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


TEST( CTSelectionTest, CTSelection1WithReset ) {
  // Select fields, spws, and antennas, then reset antenna selection

  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(20), nAnt(10), nSpw(4), nObs(1), nScan(20),nTime(10);
  Vector<Int> nChan(nSpw, 32);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(60.0);
  NewCalTable tnct("tCTSelection1.ct","Complex",
		   nObs, nScan, nTime,
		   nAnt, nSpw, nChan,
		   nFld,
		   refTime, tint);
  // some sanity checks on the test NewCalTable
  EXPECT_EQ( tnct.tableType(), Table::Memory );
  EXPECT_EQ( tnct.nrow(), nObs * nScan * nTime * nSpw * nAnt);

  ROCTColumns ctc(tnct);
  // Field selection by name
  Vector<String> fldnames = ctc.field().name().getColumn();
  Vector<Int> fldids(3);
  fldids(0)=3; fldids(1)=10; fldids(2)=17;
  String fieldsel("");
  for (uInt i=0; i<fldids.nelements(); ++i) {
    if (i>0) fieldsel += ",";
    fieldsel += fldnames(fldids(i));
  }
  // Spw selection
  Vector<Int> spwids(2);
  spwids(0)=1; spwids(1)=3;
  String spwsel("");
  for (uInt i=0; i<spwids.nelements(); ++i) {
    if (i>0) spwsel += ",";
    spwsel += String::toString(spwids(i));
  }
  // Antenna selection by name, including ref antenna 0 
  // to check that not all baselines/rows were selected
  Vector<String> antnames=ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  antids(0)=0; antids(1)=3; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0; i<antids.nelements(); ++i) {
    if (i>0) antsel += ",";
    antsel += antnames(antids(i));
  }

  // Create CTSelection and set selections
  CTSelection cts;
  cts.setFieldExpr(fieldsel);
  cts.setSpwExpr(spwsel);
  cts.setAntennaExpr(antsel);
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(cts.getFieldList(), fldids) );
  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(cts.getSpwList(), spwids) );

  // Create selected table and check nrow
  NewCalTable selnct(tnct);
  getSelectedTable(selnct,tnct,ten,"");
  ASSERT_EQ( selnct.nrow(), 
    antids.nelements() * nTime * spwids.nelements() * fldids.nelements());

  // get columns from selected table and check elements
  ROCTMainColumns ctmc(selnct);
  Vector<Int> selfieldcol, selantcol, selant2col, selspwcol;
  ctmc.fieldId().getColumn(selfieldcol);
  ctmc.antenna1().getColumn(selantcol);
  ctmc.antenna2().getColumn(selant2col);
  ctmc.spwId().getColumn(selspwcol);

  Vector<Bool> fldok(selfieldcol.nelements(), false);
  for (uInt i=0; i<fldids.nelements(); ++i) {
    fldok |= (selfieldcol==fldids(i));
  }
  Vector<Bool> spwok(selspwcol.nelements(), false);
  for (uInt i=0; i<spwids.nelements(); ++i) {
    spwok |= (selspwcol==spwids(i));
  }
  Vector<Bool> antok(selantcol.nelements(), false);
  for (uInt i=0; i<antids.nelements(); ++i) {
    antok |= (selantcol==antids(i));
  }

  EXPECT_TRUE( allEQ(fldok, true) );
  EXPECT_TRUE( allEQ(spwok, true) );
  EXPECT_TRUE( allEQ(antok, true) );

  // reset with different antenna selection, same field & spw selections
  antsel = "0~3;!1&0";
  antids(0)=0; antids(1)=1; antids(2)=2; antids(3)=3; antids(4)=-1; 
  cts.reset(cti, casacore::MSSelection::PARSE_NOW, "", antsel,
    fieldsel, spwsel, "", "", "", "");

  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(cts.getFieldList(), fldids) );
  EXPECT_TRUE( allEQ(cts.getSpwList(), spwids) );

  // Create selected table and check nrow
  ten = cts.getTEN();
  getSelectedTable(selnct, tnct, ten, "");
  // antenna 1 is selected and negated
  size_t antids_size = 3;
  ASSERT_EQ( selnct.nrow(), 
    antids_size * nTime * spwids.nelements() * fldids.nelements() );
}


TEST( CTSelectionTest, CTSelection2 ) {
  // Test another selection and CTSelection ctor with selection expressions

  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(1), nAnt(10), nSpw(4), nObs(3), nScan(6), nTime(6);
  Vector<Int> nChan(nSpw, 1);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(10.0);
  NewCalTable tnct("tCTSelection2.ct", "Complex",
		   nObs,nScan,nTime,
		   nAnt,nSpw,nChan,
		   nFld,
		   refTime,tint);
  // some sanity checks on the test NewCalTable
  EXPECT_EQ( tnct.tableType(), Table::Memory );
  EXPECT_EQ( tnct.nrow(), nObs * nScan * nTime * nSpw * nAnt );

  // ObsID selection
  Vector<Int> obsids(1);
  obsids(0) = 1;
  String obssel = String::toString(obsids(0));
  // Scan selection
  Vector<Int> scans(1);
  scans(0) = 7;
  String scansel = String::toString(scans(0));
  // Time selection expression:
  //  this should pick 3 timestamps in obsid=1, scan=7
  String timesel("2012/01/06/12:06:15.0~12:06:45");
  Vector<Double> timebounds(3);
  timebounds(0) = 375.0+refTime;
  timebounds(1) = 405.0+refTime;
  timebounds(2) = 0.0; // dT
  // Spw selection
  Vector<Int> spwids(2);  
  spwids(0)=1; spwids(1)=3;
  String spwsel("");
  for (uInt i=0; i<spwids.nelements(); ++i) {
    if (i>0) spwsel += ",";
    spwsel += String::toString(spwids(i));
  }
  // Antenna selection by name
  ROCTColumns ctc(tnct);
  Vector<String> antnames = ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  antids(0)=2; antids(1)=3; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0; i<antids.nelements(); ++i) {
    if (i>0) antsel += ",";
    antsel += antnames(antids(i));
  }

  // Create CTSelection and set selections
  CTSelection cts;
  cts.setObservationExpr(obssel);
  cts.setScanExpr(scansel);
  cts.setTimeExpr(timesel);
  cts.setSpwExpr(spwsel);
  cts.setAntennaExpr(antsel);
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(cts.getObservationList(),obsids) );  
  EXPECT_TRUE( allEQ(cts.getScanList(), scans) );
  EXPECT_TRUE( allNear(Vector<Double>(cts.getTimeList()),timebounds,1.0e-9) );
  EXPECT_TRUE( allEQ(cts.getSpwList(), spwids) );
  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(cts.getObservationList(), obsids) );  

  // Create selected table and check nrow
  NewCalTable selnct(tnct);
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ( selnct.nrow(), 3 * antids.nelements() * spwids.nelements() );

  // get columns from selected table and check elements
  ROCTMainColumns ctmc(selnct);
  Vector<Int> selantcol, selspwcol, selobscol, selscancol;
  Vector<Double> seltimecol;
  ctmc.obsId().getColumn(selobscol);
  ctmc.scanNo().getColumn(selscancol);
  ctmc.time().getColumn(seltimecol);
  ctmc.spwId().getColumn(selspwcol);
  ctmc.antenna1().getColumn(selantcol);

  Vector<Bool> obsok(selobscol.nelements(), false);
  for (uInt i=0; i<obsids.nelements(); ++i) 
    obsok |= (selobscol==obsids(i));

  Vector<Bool> scanok(selscancol.nelements(), false);
  for (uInt i=0; i<scans.nelements(); ++i) 
    scanok |= (selscancol==scans(i));

  Vector<Bool> timeok(seltimecol.nelements(), false);
  timeok |= (seltimecol>timebounds(0));
  timeok |= (seltimecol<timebounds(1));

  Vector<Bool> spwok(selspwcol.nelements(), false);
  for (uInt i=0; i<spwids.nelements(); ++i) 
    spwok |= (selspwcol==spwids(i));

  Vector<Bool> antok(selantcol.nelements(), false);
  for (uInt i=0; i<antids.nelements(); ++i) 
    antok |= (selantcol==antids(i));
  
  EXPECT_TRUE( allEQ(obsok, true) );
  EXPECT_TRUE( allEQ(scanok, true) );
  EXPECT_TRUE( allEQ(timeok, true) );
  EXPECT_TRUE( allEQ(spwok, true) );
  EXPECT_TRUE( allEQ(antok, true) );

  // Repeat test using CTSelection constructor with selection expressions
  CTSelection* newcts = new CTSelection(tnct, casacore::MSSelection::PARSE_NOW,
    timesel, antsel, "", spwsel, "", scansel, "", obssel);
  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(newcts->getObservationList(), obsids) );  
  EXPECT_TRUE( allEQ(newcts->getScanList(), scans) );
  EXPECT_TRUE( allNear(Vector<Double>(newcts->getTimeList()), timebounds, 
    1.0e-9) );
  EXPECT_TRUE( allEQ(newcts->getSpwList(), spwids) );
  EXPECT_TRUE( allEQ(newcts->getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(newcts->getObservationList(), obsids) );  

  // Create selected table and check nrow
  ten = newcts->getTEN();
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ( selnct.nrow(), 3 * antids.nelements() * spwids.nelements() );
}

