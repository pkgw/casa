//# PlotMSCacheInvalid_GT.cc:: GoogleTest for invalid axes (for test MS)
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

using namespace casa;

class PlotMSCacheTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		// Using visstat2 regression dataset:
		// All datacolumns present: data, model, corrected
		dataPath = tUtil::getFullPath( "ngc5921_add_corect_model.ms", "visstat2" );
		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
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

TEST_F( PlotMSCacheTest, testInvalidAxes) {
	// test cal table and ephem axes
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	std::vector<PMS::Axis> invalidAxes {PMS::GAMP, PMS::GPHASE,
		PMS::GREAL, PMS::GIMAG, PMS::DELAY, PMS::SWP, PMS::TSYS,
		PMS::OPAC, PMS::SNR, PMS::TEC, PMS::ANTPOS,
		  // Ephemeris
	    PMS::RADIAL_VELOCITY, PMS::RHO};

	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;

	for (auto axis : invalidAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};

		ASSERT_THROW(cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,	itsTransformations,
			itsCalibration, nullptr ), AipsError);
	}
}
