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
#include <ms/MSSel/MSSelectionTools.h>

using namespace casa;

class PlotMSCacheTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		// Using mstranform regression dataset:
		// ASDM_RECEIVER table for showimage
		dataPath = tUtil::getFullPath("split_ddid_mixedpol_CAS-12283.ms", "mstransform");

		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);

		// specific to this dataset:
		expNChunk = 1; // one spw selected
		expNRow = 25;  // rows in first chunk
		expNChan = 64;
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
	Int expNChunk, expNRow;
	uInt expNChan;
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

TEST_F( PlotMSCacheTest, testAtmTskyOverlays) {
	// Test overlay axis options
	// ATM, TSKY

	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// selection, axes, datacolumns
	itsSelection.setSpw("58");
	std::vector<PMS::Axis> atmAxes {PMS::ATM, PMS::TSKY};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA, PMS::DATA, PMS::DATA};

	// check values for first chunk
	Int ichunk(0);
	for (auto yaxis : atmAxes) {
		std::vector<PMS::Axis> loadAxes{PMS::FREQUENCY, PMS::FREQUENCY, PMS::AMP, yaxis};
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		ASSERT_EQ(expNChunk, cache->nChunk());
		ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
		// Check arrays with ArrayLogical
		switch(yaxis) {
			case PMS::ATM: {  // percent
				Vector<Double> atm(cache->atm(ichunk));
				ASSERT_EQ(expNChan, atm.size());
				ASSERT_TRUE(allGT(atm, 97.0));
				ASSERT_TRUE(allLT(atm, 98.0));
				break;
			}
			case PMS::TSKY: { // Kelvin
				Vector<Double> tsky(cache->tsky(ichunk));
				ASSERT_EQ(expNChan, tsky.size());
				ASSERT_TRUE(allGT(tsky, 9.0));
				ASSERT_TRUE(allLT(tsky, 10.0));
				break;
			}
			default:
				break;
		}
	}
}

TEST_F( PlotMSCacheTest, testImageOverlay) {
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	// selection, axes, datacolumns
	itsSelection.setSpw("58");
	std::vector<PMS::Axis> loadAxes {PMS::FREQUENCY, PMS::FREQUENCY, PMS::FREQUENCY, PMS::AMP, PMS::ATM, PMS::IMAGESB};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA, PMS::DATA, PMS::DATA, PMS::DATA, PMS::DATA};

	// check values for first (only) chunk
	int ichunk(0);

	cache->load(loadAxes, loadData, dataPath, 
				itsSelection, itsAveraging,
				itsTransformations, itsCalibration,
				nullptr ); // ThreadCommunication* is nullptr 

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// Check arrays with ArrayLogical
	Vector<Double> atm(cache->atm(ichunk));
	ASSERT_EQ(expNChan, atm.size());
	ASSERT_TRUE(allGT(atm, 97.0));
	ASSERT_TRUE(allLT(atm, 98.0));
	Vector<Double> image(cache->imageSideband(ichunk));
	ASSERT_EQ(expNChan, image.size());
	ASSERT_TRUE(allGT(image, 98.0));
	ASSERT_TRUE(allLT(image, 99.0));
}
