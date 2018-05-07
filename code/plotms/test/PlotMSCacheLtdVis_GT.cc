//# PlotMSCacheLtdVis_GT.cc:: GoogleTest for limited visibilities (one col)
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

using namespace casa;

/* test dataset visibility axes with one column: DATA or FLOAT_DATA */

int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

TEST( PlotMSCacheTest, testDataOnly ) {
	// dataset only has DATA column
	String dataPath = tUtil::getFullPath( "pm_ngc5921.ms", "plotms" );
	Int ichunk(0), expNChunk(60), expNRow(351);

	MeasurementSet ms(dataPath);
	Block<String> sortCols(4); // Use default sort columns
	sortCols[0] = MS::columnName(MS::ARRAY_ID);
	sortCols[1] = MS::columnName(MS::FIELD_ID);
	sortCols[2] = MS::columnName(MS::DATA_DESC_ID);
	sortCols[3] = MS::columnName(MS::TIME);
	MeasurementSet sortedMS = ms.sort(sortCols);
	// main table columns
	ROMSMainColumns msmc(sortedMS);
	Array<Complex> visdata(msmc.data().getColumn());
	IPosition visshape = visdata.shape();
	visshape.setLast(IPosition(1,expNRow));
	visdata.adjustLastAxis(visshape);
	Array<Float> ampdata(amplitude(visdata));

	MSCache* cache = new MSCache(nullptr);
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;

	// Datacolumn options:
	std::vector<PMS::DataColumn> visCols {PMS::DATA, PMS::CORRECTED,
		PMS::MODEL, PMS::CORRMODEL, PMS::DATAMODEL,
		PMS::DATA_DIVIDE_MODEL, PMS::CORRECTED_DIVIDE_MODEL};
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	for (auto datacol : visCols) {
		std::vector<PMS::DataColumn> loadData{datacol, PMS::DATA};
		switch (datacol) {
			case PMS::DATA:
			case PMS::CORRECTED:
			case PMS::CORRMODEL:
			case PMS::CORRECTED_DIVIDE_MODEL: {
				// no corrected data, uses data instead
				cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,	itsTransformations,
					itsCalibration, nullptr );
				EXPECT_EQ(expNChunk, cache->nChunk());
				Int nrow = cache->chunkShapes()(IPosition(2,2,ichunk));
				EXPECT_EQ(expNRow, nrow);
				EXPECT_TRUE(allEQ(ampdata, cache->ampData(ichunk)));
				break;
				}
			case PMS::MODEL:
			case PMS::DATAMODEL:
			case PMS::DATA_DIVIDE_MODEL: {
				// generates model dynamically
				ASSERT_NO_THROW(cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,	itsTransformations,
					itsCalibration, nullptr ));
				ASSERT_EQ(expNChunk, cache->nChunk());
				Int nrow = cache->chunkShapes()(IPosition(2,2,ichunk));
				ASSERT_EQ(expNRow, nrow);
				break;
				}
			case PMS::FLOAT_DATA: {
				// FLOAT_DATA throws exception for interferometry data
				ASSERT_THROW(cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,	itsTransformations,
					itsCalibration, nullptr ), AipsError);
				break;
				}
		}
	}
	delete cache;
}

TEST( PlotMSCacheTest, testFloatOnly) {
	// Test visibility axis options for singledish:
	// Axes:
	// 	amp = real => float data
	// 	phase = imag => invalid
	// Data columns:
	// 	DATA = CORRECTED = FLOAT => float data
	// 	MODEL = data/corrected residuals => invalid

	// use singledish dataset, only has FLOAT_DATA column
	String dataPath = tUtil::getFullPath( "pointing6.ms", "sdimaging" );
	MeasurementSet ms(dataPath);
	Block<String> sortCols(4); // Use default sort columns
	sortCols[0] = MS::columnName(MS::ARRAY_ID);
	sortCols[1] = MS::columnName(MS::FIELD_ID);
	sortCols[2] = MS::columnName(MS::DATA_DESC_ID);
	sortCols[3] = MS::columnName(MS::TIME);
	MeasurementSet sortedMS = ms.sort(sortCols);

	// main table columns: get first row (==first chunk)
	ROMSMainColumns msmc(sortedMS);
	Array<Float> floatData = msmc.floatData().get(0); // matrix
	Cube<Float> floatCube = floatData.addDegenerate(1); // cube
	Int expNChunk(1000), expNRow(1);  // since only 1 "baseline"

	// Visibility axis options:
	std::vector<PMS::Axis> visAxes {PMS::AMP, PMS::PHASE,
		PMS::REAL, PMS::IMAG};
	// Datacolumn options:
	std::vector<PMS::DataColumn> visCols {PMS::DATA, PMS::CORRECTED,
		PMS::MODEL, PMS::CORRMODEL, PMS::DATAMODEL,
		PMS::DATA_DIVIDE_MODEL, PMS::CORRECTED_DIVIDE_MODEL};

	MSCache* cache = new MSCache(nullptr);
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;

	Int ichunk(0);
	for (auto axis : visAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		for (auto datacol : visCols) {
			std::vector<PMS::DataColumn> loadData{datacol, PMS::DATA};
			switch (axis) {
				case PMS::AMP:
				case PMS::REAL: {
					switch (datacol) {
						case PMS::DATA:
						case PMS::CORRECTED:
						case PMS::FLOAT_DATA: {
							cache->load(loadAxes, loadData, dataPath, 
								itsSelection, itsAveraging,
								itsTransformations, itsCalibration, 
								nullptr ); // ThreadCommunication* is nullptr
							ASSERT_EQ(expNChunk, cache->nChunk());
							Int nrow = cache->chunkShapes()(IPosition(2,2,ichunk));
							ASSERT_EQ(expNRow, nrow);
							if (axis==PMS::AMP)
								ASSERT_TRUE(allEQ(floatCube, cache->ampFloat(ichunk)));
							if (axis==PMS::REAL)
								ASSERT_TRUE(allEQ(floatCube, cache->realData(ichunk)));
							break;
							}
						case PMS::MODEL:
						case PMS::DATAMODEL:
						case PMS::DATA_DIVIDE_MODEL:
						case PMS::CORRMODEL:
						case PMS::CORRECTED_DIVIDE_MODEL: {
							ASSERT_THROW(cache->load(loadAxes, loadData, dataPath, 
								itsSelection, itsAveraging,	itsTransformations,
								itsCalibration, nullptr ), AipsError);
							break;
							}
					}
					break;
				}
				case PMS::PHASE:
				case PMS::IMAG: {
					ASSERT_THROW(cache->load(loadAxes, loadData, dataPath, 
						itsSelection, itsAveraging,	itsTransformations,
						itsCalibration, nullptr ), AipsError);
					break;
				}
				default:
					break;
			}
		}
	}
	delete cache;
}

