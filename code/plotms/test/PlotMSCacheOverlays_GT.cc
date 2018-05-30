//# PlotMSCacheOverlays_GT.cc:: GoogleTest for plotms cache overlays
//# Copyright (C) 2018
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
//#
//# $Id$

#include <gtest/gtest.h>
#include <iostream>

#include <plotms/test/tUtil.h>
#include <plotms/Data/MSCache.h>
#include <casa/Arrays/ArrayLogical.h>
//#include <casa/Arrays/ArrayMath.h>

using namespace casa;

class PlotMSCacheTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		// Using visstat2 regression dataset:
		// All datacolumns present: data, model, corrected
		dataPath = tUtil::getFullPath( "ngc5921_add_corect_model.ms", "visstat2" );
		// specific to this dataset:
		expNChunk = 60;
		expNRow = 351;  // rows in first chunk
		expNChan = 63;
		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
	Int expNChunk, expNRow;
	uInt expNChan;  // expected values
	MeasurementSet sortedMS;
	// use defaults (none!)
	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;
	MSCache* cache;
};



int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

TEST_F( PlotMSCacheTest, testOverlays) {
	// Test overlay axis options
	// ATM, TSKY

	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	std::vector<PMS::Axis> atmAxes {PMS::ATM, PMS::TSKY};
	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA, PMS::DATA, PMS::DATA};

	// check values for first chunk
	Int ichunk(0);
	for (auto yaxis : atmAxes) {
		std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::FREQUENCY, yaxis, PMS::FREQUENCY};
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		ASSERT_EQ(expNChunk, cache->nChunk());
		ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
		// Check arrays with ArrayLogical.h
		switch(yaxis) {
			case PMS::ATM: {  // percent
				Vector<Double> atm(cache->atm(ichunk));
				ASSERT_EQ(expNChan, atm.size());
				ASSERT_TRUE(allGT(atm, 99.0));
				ASSERT_TRUE(allLT(atm, 100.0));
				break;
			}
			case PMS::TSKY: { // Kelvin
				Vector<Double> tsky(cache->tsky(ichunk));
				ASSERT_EQ(expNChan, tsky.size());
				ASSERT_TRUE(allGT(tsky, 0.0));
				ASSERT_TRUE(allLT(tsky, 10.0));
				break;
			}
			default:
				break;
		}
	}
}

