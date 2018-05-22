//# PlotMSCacheAveraging_GT.cc:: GoogleTest for plotms cache averaging modes
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
#include <ms/MSSel/MSSelectionTools.h>
#include <casa/Arrays/ArrayLogical.h>

using namespace casa;

class PlotMSCacheTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		// Using visstat2 regression dataset:
		// All datacolumns present: data, model, corrected
		dataPath = tUtil::getFullPath( "ngc5921_add_corect_model.ms", "visstat2" );
		MeasurementSet ms(dataPath);
		Block<String> sortCols(4); // Use default sort columns
		sortCols[0] = MS::columnName(MS::ARRAY_ID);
		sortCols[1] = MS::columnName(MS::FIELD_ID);
		sortCols[2] = MS::columnName(MS::DATA_DESC_ID);
		sortCols[3] = MS::columnName(MS::TIME);
		sortedMS = ms.sort(sortCols);
		// specific to this dataset:
		expNChunk = 60;
		expNCorr = 2;
		expNChan = 63;
		expNRow = 351;  // rows in first chunk
		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
	Int expNChunk, expNCorr, expNChan, expNRow;
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

TEST_F( PlotMSCacheTest, testAvgchannel) {
	// Test visibilities with channel averaging
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	// average all channels; pipeline uses large number
	String avgchan("1e6");
	itsAveraging.setChannel(true); 
	itsAveraging.setChannelValue(String::toDouble(avgchan)); 
	expNChan = 1;

	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(sortedMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());
	// adjust cube shapes for first chunk
	IPosition visShape(corrData.shape());
	visShape.setLast(IPosition(1,expNRow));
	corrData.adjustLastAxis(visShape);
	wtsp.adjustLastAxis(visShape);

	// average channel axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int row=0; row<expNRow; ++row) {
		for (Int corr=0; corr<expNCorr; ++corr) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(), Slice(row)));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(), Slice(row)));
			dataSlice *= wtspSlice;
			avgCorrData(corr, 0, row) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNChan, cache->ampCorr(ichunk).shape()(1));
	ASSERT_TRUE(allEQ(expAmpCorr, cache->ampCorr(ichunk)));
}

TEST_F( PlotMSCacheTest, testAvgTime) {
	// Test visibilities with time averaging
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	expNChunk = 7; // since chunks are averaged together

	// average all times in chunk
	String avgtime("1e6");
	itsAveraging.setTime(true); 
	itsAveraging.setTimeValue(String::toDouble(avgtime));
	// First chunk is first scan (1) and first field (0).
    // Time averaging is by baseline; select one baseline for test
	String antennaExpr("2&3"), scanExpr("1"), fieldExpr("0");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		fieldExpr, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		scanExpr, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());
	// adjust cube shapes for first chunk
	IPosition visShape(corrData.shape());

	// average time (row) axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(chan), Slice()));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			avgCorrData(corr, chan, 0) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmpCorr, cache->ampCorr(ichunk)));
}

TEST_F( PlotMSCacheTest, testAvgScan) {
	// Test visibilities with time averaging over scan
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	expNChunk = 3; // since chunks are averaged together

	// average all times in chunk
	String avgtime("1e6");
	itsAveraging.setTime(true); 
	itsAveraging.setTimeValue(String::toDouble(avgtime));
	itsAveraging.setScan(true);
	// First chunk is first field (0) only, averaging over scan
    // Time averaging is by baseline; select one baseline for test
	String antennaExpr("2&3"), fieldExpr("0");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		fieldExpr, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		""/*scan*/, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());

	// average time (row) axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(chan), Slice()));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			avgCorrData(corr, chan, 0) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmpCorr, cache->ampCorr(ichunk)));
}

TEST_F( PlotMSCacheTest, testAvgField) {
	// Test visibilities with time averaging over field
	// NOTE: chunks are per scan, so averaging over field makes no difference
	// (compared to plain time-averaging) since fields do not change mid-scan

	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	expNChunk = 7; // since chunks are averaged together

	// average all times in chunk
	String avgtime("1e6");
	itsAveraging.setTime(true); 
	itsAveraging.setTimeValue(String::toDouble(avgtime));
	itsAveraging.setField(true);
	// First chunk is first scan (1)
    // Time averaging is by baseline; select one baseline for test
	String antennaExpr("2&3"), scanExpr("1");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		""/*field*/, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		scanExpr, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());

	// average time (row) axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(chan), Slice()));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			avgCorrData(corr, chan, 0) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmpCorr, cache->ampCorr(ichunk)));
}

TEST_F( PlotMSCacheTest, testAvgScanField) {
	// Test visibilities with time averaging over scan and field
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	expNChunk = 1; // since *all* chunks are averaged together

	// average all times in chunk
	String avgtime("1e6");
	itsAveraging.setTime(true); 
	itsAveraging.setTimeValue(String::toDouble(avgtime));
	itsAveraging.setScan(true);
	itsAveraging.setField(true);
	// First chunk is first field (0) only, averaging over scan
    // Time averaging is by baseline; select one baseline for test
	String antennaExpr("2&3");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		""/*field*/, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		""/*scan*/, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());

	// average time (row) axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(chan), Slice()));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			avgCorrData(corr, chan, 0) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allNear(expAmpCorr, cache->ampCorr(ichunk), .000001));
}

TEST_F( PlotMSCacheTest, testAvgBaseline) {
	// Select two baselines and average them together
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	Int vbNRow(2); // each chunk has one row per baseline

	// selections
	itsAveraging.setBaseline(True);
	String antennaExpr("2&3; 2&4");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		""/*field*/, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		""/*scan*/, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	// PlotMSVBAverager uses weight instead of weight spectrum
	Matrix<Float> wt(msmc.weight().getColumn());
	// vb has two rows, one per baseline
	IPosition visShape(corrData.shape());
	visShape.setLast(IPosition(1,vbNRow));
	corrData.adjustLastAxis(visShape);
	IPosition wtShape(wt.shape());
	wtShape.setLast(IPosition(1,vbNRow));
	wt.adjustLastAxis(wtShape);

	// make wtsp cube from wt: fill in same value for all chans
	Cube<Float> wtsp(visShape);
	Float weightVal;
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			for (Int row=0; row<visShape(2); ++ row) {
				// PlotMSVBAverager uses corr 0 for weight
				weightVal = wt(0,row);
				if (weightVal < FLT_MIN) weightVal = FLT_MIN;
				wtsp(corr, chan, row) = weightVal;
			}
		}
	}

	// average baselines (row) axis
	Cube<Complex> avgCorrData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Complex> dataSlice = corrData(Slicer(Slice(corr), Slice(chan), Slice()));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			avgCorrData(corr, chan, 0) = sum(dataSlice) / sum(wtspSlice);
		}
	}
	// amplitude of averaged data
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmpCorr, cache->ampCorr(ichunk)));
}

TEST_F( PlotMSCacheTest, testAvgAntenna) {
	// Select one antenna and average data for antennas in chunk
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 26; // nAnt for chunk
	Int vbNRow(25); // each chunk has one row per baseline

	// set selection and averaging
	String antennaExpr("0");
	itsAveraging.setAntenna(True);
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		""/*field*/, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		""/*scan*/, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Cube<Complex> corrData(msmc.correctedData().getColumn());
	// PlotMSVBAverager uses weight instead of weight spectrum
	Matrix<Float> wt(msmc.weight().getColumn());
	Vector<Int> ant1(msmc.antenna1().getColumn()),
		ant2(msmc.antenna2().getColumn());

	// resize for first chunk
	IPosition visShape(corrData.shape());
	visShape.setLast(IPosition(1,vbNRow));
	corrData.adjustLastAxis(visShape);
	IPosition wtShape(wt.shape());
	wtShape.setLast(IPosition(1,vbNRow));
	wt.adjustLastAxis(wtShape);
	ant1.resize(vbNRow, true);
	ant2.resize(vbNRow, true);
	// set of all antenna ids
	Vector<Int> allAntIds = concatenateArray(ant1, ant2);
	Int nAnts = GenSort<Int>::sort(allAntIds, Sort::Ascending, Sort::NoDuplicates);
	allAntIds.resize(nAnts, true);

	// make wtsp cube from wt: fill in same value for all chans
	Cube<Float> wtsp(visShape);
	Float weightVal;
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			for (Int row=0; row<visShape(2); ++ row) {
				// PlotMSVBAverager uses corr 0 for weight
				weightVal = wt(0,row);
				if (weightVal < FLT_MIN) weightVal = FLT_MIN;
				wtsp(corr, chan, row) = weightVal;
			}
		}
	}

	// accumulate weighted data and weights, one row per antenna id
	IPosition outShape(3, expNCorr, expNChan, nAnts);
	Cube<Complex> accumCorrData(outShape, 0.0);
	Cube<Float> wtPerAnt(outShape, 0.0);
	for (Int antIdx=0; antIdx<nAnts; ++antIdx) {
		Int antId(allAntIds(antIdx));
		for (uInt row=0; row<ant1.size(); ++row) {
			if (ant1(row)==antId || ant2(row)==antId) {
				Matrix<Float> rowWt = wtsp.xyPlane(row);
				Matrix<Complex> antWtData = corrData.xyPlane(row) * rowWt;
				// no += for Matrix
				Matrix<Complex> accumAntWtData = antWtData + accumCorrData.xyPlane(antIdx);
				accumCorrData.xyPlane(antIdx) = accumAntWtData;
				wtPerAnt.xyPlane(antIdx) = wtPerAnt.xyPlane(antIdx) + rowWt;
			}
		}
	}
	// divide accum data by accum weight to get expected averaged corrected data
	Cube<Complex> avgCorrData(expNCorr, expNChan, nAnts);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			for (Int ant=0; ant<nAnts; ++ant) {
				avgCorrData(corr, chan, ant) =
					accumCorrData(corr, chan, ant) / wtPerAnt(corr, chan, ant);
			}
		}
	}
	Cube<Float> expAmpCorr(amplitude(avgCorrData));

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allNear(expAmpCorr, cache->ampCorr(ichunk), .000001));
}

TEST_F( PlotMSCacheTest, testAvgSpw) {
	// average each channel over spw (spw nchans must match)
	// Use dataset with more than one spw:
	dataPath = tUtil::getFullPath( "Four_ants_3C286.ms", "flagdata" );
	MeasurementSet ms(dataPath);
	Block<String> sortCols(4); // change order of sort cols to average spws
	sortCols[0] = MS::columnName(MS::ARRAY_ID);
	sortCols[1] = MS::columnName(MS::FIELD_ID);
	sortCols[2] = MS::columnName(MS::TIME);
	sortCols[3] = MS::columnName(MS::DATA_DESC_ID);
	sortedMS = ms.sort(sortCols);

	// Select 2 spws and 2 channels; will average chan0 values and chan1 values
	// separately over the 2 spws.
	// Also select first timestamp for rows in first chunk only.
	String spwExpr("1~2:0~1"), timeExpr("14:45:08.50");
	// specific to this dataset:
	expNChunk=179; expNCorr=4; expNChan=2; expNRow=6;
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, "", timeExpr, "", "",
		spwExpr, "", "", "", "", "", "", "", "");
	ROMSMainColumns msmc(selMS);
	Cube<Complex> visData(msmc.data().getColumn());
	Cube<Bool> flags(msmc.flag().getColumn());
	Matrix<Float> wt(msmc.weight().getColumn());
	IPosition visShape(visData.shape());

	// make wtsp cube: fill in same value for all chans
	Cube<Float> wtsp(visShape);
	Float weightVal;
	for (Int corr=0; corr<visShape(0); ++corr) {
		for (Int chan=0; chan<visShape(1); ++chan) {
			for (Int row=0; row<visShape(2); ++ row) {
				weightVal = wt(0,row);
				if (weightVal < FLT_MIN) weightVal = FLT_MIN;
				wtsp(corr, chan, row) = weightVal;
			}
		}
	}

	Cube<Complex> avgVisData(expNCorr, expNChan, expNRow);
	Complex	dataSpw1, dataSpw2;
	Float wtSpw1, wtSpw2;
	Bool flagSpw1, flagSpw2;
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			for (Int row=0; row<expNRow; ++row) {
				// for spw1 
				flagSpw1 = flags(corr,chan,row);
				dataSpw1 = visData(corr, chan, row);
				wtSpw1 = wtsp(corr, chan, row);
				// for spw2 
				dataSpw2 = visData(corr, chan, row+expNRow);
				flagSpw2 = flags(corr,chan,row+expNRow);
				wtSpw2 = wtsp(corr, chan, row+expNRow);
				// average channels over both spws
				if (flagSpw1!=flagSpw2) { // use unflagged data only
					if (flagSpw1)
						avgVisData(corr,chan,row) = dataSpw2;
					else
						avgVisData(corr,chan,row) = dataSpw1;
				} else {  // average data together (both are flagged/unflagged)
					avgVisData(corr,chan, row) = (dataSpw1*wtSpw1 + dataSpw2*wtSpw2) /
						(wtSpw1 + wtSpw2);
				}
			}
		}
	}
	Cube<Float> expAmpData(amplitude(avgVisData));

	// set up selection and averaging for cache
	PlotMSSelection itsSelection;
	itsSelection.setSpw(spwExpr);
	PlotMSAveraging itsAveraging;
	itsAveraging.setSpw(true);
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;

	// check values for first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::CHANNEL};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr
	Int nchunk(cache->nChunk()), nchan(cache->chunkShapes()(IPosition(2,1,ichunk))),
	   nrow(cache->chunkShapes()(IPosition(2,2,ichunk)));
	Cube<Float> cacheAmpData(cache->ampData(ichunk));

	ASSERT_EQ(expNChunk, nchunk);
	ASSERT_EQ(expNChan, nchan);
	ASSERT_EQ(expNRow, nrow);
	ASSERT_TRUE(allNear(expAmpData, cacheAmpData, .0000001));
}

TEST_F( PlotMSCacheTest, testAvgTimeScalar) {
	// Test visibilities with time averaging
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	expNRow = 1; // all rows averaged into 1 row
	expNChunk = 7; // since chunks are averaged together

	// average all times in chunk
	String avgtime("1e6");
	itsAveraging.setTime(true); 
	itsAveraging.setTimeValue(String::toDouble(avgtime));
	itsAveraging.setScalarAve(true);
	// First chunk is first scan (1) and first field (0).
    // Time averaging is by baseline; select one baseline for test
	String antennaExpr("2&3"), scanExpr("1"), fieldExpr("0");
	itsSelection.setAntenna(antennaExpr);	

	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, antennaExpr,
		fieldExpr, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		scanExpr, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns; ignoring flags since all unflagged data
	ROMSMainColumns msmc(selMS);
	Array<Complex> corrData(msmc.correctedData().getColumn());
	Array<Float> wtsp(msmc.weightSpectrum().getColumn());
	// adjust cube shapes for first chunk
	IPosition visShape(corrData.shape());

	// average amplitudes over time (row) axis
	Cube<Float> expScalarAvgData(expNCorr, expNChan, expNRow);
	for (Int corr=0; corr<expNCorr; ++corr) {
		for (Int chan=0; chan<expNChan; ++chan) {
			Array<Float> dataSlice = amplitude(corrData(Slicer(Slice(corr), Slice(chan), Slice())));
			Array<Float> wtspSlice = wtsp(Slicer(Slice(corr), Slice(chan), Slice()));
			dataSlice *= wtspSlice;
			expScalarAvgData(corr,chan,0) = sum(dataSlice) / sum(wtspSlice);
		}
	}

	// load cache and check first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes {PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::CORRECTED, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration, 
		nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allNear(expScalarAvgData, cache->ampCorr(ichunk), .000001));
}

