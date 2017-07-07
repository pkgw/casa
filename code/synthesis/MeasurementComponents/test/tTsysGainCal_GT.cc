//# tTsysGainCal_GT.cc: test Tsys calibration
//# Copyright (C) 1995,1999,2000,2001,2016,2017
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

#include <casa/aips.h>
#include <casa/iostream.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/OS/EnvVar.h>
#include <casa/OS/Path.h>
#include <casa/OS/Timer.h>
#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <synthesis/MeasurementComponents/TsysGainCal.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>

#include <gtest/gtest.h>

#include "VisCalTestBase_GT.h"

using namespace casacore;
using namespace casa;
using namespace casa::vi;

// Control verbosity
#define TSYS_TEST_VERBOSE false

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(TsysSpecifyTest, TsysSpectrum) {

  // Path to an MS in the DR
  String *casapath = new String[2];
  split(EnvironmentVariable::get("CASAPATH"), casapath, 2, String(" "));
  // Use of Path().absoluteName() absorbs relative stuff in casapath
  String mspath(Path(casapath[0]+"/data/regression/unittest/MSMetaData/MSMetaData.ms").absoluteName());

  MSMetaInfoForCal msm(mspath);
  StandardTsys Tsysapp(msm);
  Record specpar;

  specpar.define("uniform", true);
  Tsysapp.setSpecify(specpar);
  Tsysapp.specify(specpar);

  ASSERT_EQ(VisCalEnum::JONES,Tsysapp.matrixType());
  ASSERT_EQ(VisCal::B,Tsysapp.type());
  ASSERT_EQ(String("B TSYS"),Tsysapp.typeName());
  ASSERT_TRUE(Tsysapp.freqDepPar());
  ASSERT_TRUE(Tsysapp.freqDepMat());
  ASSERT_FALSE(Tsysapp.freqDepCalWt());
  ASSERT_FALSE(Tsysapp.timeDepMat());
  ASSERT_FALSE(Tsysapp.isApplied());
  ASSERT_TRUE(Tsysapp.isSolvable());

  if (TSYS_TEST_VERBOSE)
    Tsysapp.state();
}

TEST(TsysSpecifyTest, Tsys) {

  // Path to an MS in the DR
  String *casapath = new String[2];
  split(EnvironmentVariable::get("CASAPATH"), casapath, 2, String(" "));
  // Use of Path().absoluteName() absorbs relative stuff in casapath
  String mspath(Path(casapath[0]+"/data/regression/evn/n08c1.ms").absoluteName());

  MSMetaInfoForCal msm(mspath);
  StandardTsys Tsysapp(msm);
  Record specpar;

  specpar.define("uniform", true);
  Tsysapp.setSpecify(specpar);
  Tsysapp.specify(specpar);

  ASSERT_EQ(VisCalEnum::JONES,Tsysapp.matrixType());
  ASSERT_EQ(VisCal::B,Tsysapp.type());
  ASSERT_EQ(String("B TSYS"),Tsysapp.typeName());
  ASSERT_FALSE(Tsysapp.freqDepPar());
  ASSERT_FALSE(Tsysapp.freqDepMat());
  ASSERT_FALSE(Tsysapp.freqDepCalWt());
  ASSERT_FALSE(Tsysapp.timeDepMat());
  ASSERT_FALSE(Tsysapp.isApplied());
  ASSERT_TRUE(Tsysapp.isSolvable());

  if (TSYS_TEST_VERBOSE)
    Tsysapp.state();
}
