//# PlotMSCacheCallib_GT.cc:: GoogleTest for plotms OTF calibration
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
#include <ms/MeasurementSets/MSColumns.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/OS/Directory.h>

using namespace casa;

/* test loading corrected data with OTF calibration */

class PlotMSCacheTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		// Using mstransform regression dataset with callib file
		String origPath = tUtil::getFullPath( "ngc5921_regression", "mstransform" );
		Directory origDir(origPath);
		newPath = "ngc5921_regression";
		// Copy test dir because callib paths are relative
		origDir.copy(newPath);
		dataPath = newPath + "/ngc5921.ms";

		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);

		MeasurementSet ms(dataPath);
		Block<String> sortCols(4); // Use default sort columns
		sortCols[0] = MS::columnName(MS::ARRAY_ID);
		sortCols[1] = MS::columnName(MS::FIELD_ID);
		sortCols[2] = MS::columnName(MS::DATA_DESC_ID);
		sortCols[3] = MS::columnName(MS::TIME);
		sortedMS = ms.sort(sortCols);
		// specific to this dataset:
		expNChunk = 60;
		expNRow = 351;  // rows in first chunk
		expNChan = 63;
	}

	virtual void TearDown() {
		if (cache != nullptr) delete cache;
		Directory newDir(newPath);
		newDir.removeRecursive();
	}

	String newPath, dataPath;
	Int expNChunk, expNRow, expNChan;  // expected values
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

TEST_F( PlotMSCacheTest, testCallib ) {
	// dataset only has DATA column, will request CORRECTED
	// with calibration enabled
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// main table columns : DATA
	ROMSMainColumns msmc(sortedMS);
	Array<Complex> visdata(msmc.data().getColumn());
	IPosition visshape = visdata.shape();
	visshape.setLast(IPosition(1,expNRow));
	visdata.adjustLastAxis(visshape);
	Array<Float> ampdata(amplitude(visdata));

	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::CORRECTED, PMS::DATA};
	// no corrected data, uses DATA instead
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,	itsTransformations,
		itsCalibration, nullptr );
	EXPECT_EQ(expNChunk, cache->nChunk());
	Int ichunk(0);
	EXPECT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	EXPECT_TRUE(allEQ(ampdata, cache->ampData(ichunk)));

	// no corrected data, use OTF calibration
	itsCalibration.setUseCallib(True);
	String callibPath(newPath + "/ngc5921_callib.txt");
	itsCalibration.setCalLibrary(callibPath);
	//itsCalibration.setCalLibrary("ngc5921_callib.txt");
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,	itsTransformations,
		itsCalibration, nullptr );
	EXPECT_EQ(expNChunk, cache->nChunk());
	EXPECT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// get Amp:corrected, make sure didn't load DATA again 
	EXPECT_TRUE(allNE(ampdata, cache->ampCorr(ichunk)));
}
