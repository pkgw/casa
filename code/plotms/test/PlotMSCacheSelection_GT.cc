//# PlotMSCacheSelection_GT.cc:: GoogleTest for plotms cache visibility axes
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
#include <ms/MeasurementSets/MSMainColumns.h>
#include <ms/MSSel/MSSelectionTools.h>
//#include <casa/Arrays/ArrayLogical.h>

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

		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
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

TEST_F( PlotMSCacheTest, testFieldSelection) {
	/* select field, check nchunk, amp and field */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String fieldExpr("2");
	Int expNRow(378), expNChunk(34), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, fieldExpr, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	Int expField(msmc.fieldId().get(0));
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setField(fieldExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	// check amp for selected field in chunk 0
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	// check selected field in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_EQ(expField, cache->field(chunk));
	}
}

TEST_F( PlotMSCacheTest, testSpwSelection) {
	// select spw:chan and check chunks, shapes, amp, and spw/chan numbers
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// MS only has one spw, select channels
	String spwExpr("0:16~31");
	Int expNChunk(60), expNChan(16), expNRow(351), expSpw(0), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, spwExpr,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	// resize for nchans, nrows in first chunk
	Slicer visSlicer(Slice(), chanSlices(0)(0), Slice(0,expNRow));
	Array<Complex> selVis(visData(visSlicer));
	Array<Float> expAmp(amplitude(selVis));
	// make chan vector
	Vector<Int> expChans(expNChan);
	indgen(expChans);
	expChans += 16;  // start at channel 16

	// set PlotMSSelection and load cache
	itsSelection.setSpw(spwExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNChan, cache->chunkShapes()(IPosition(2,1,ichunk)));
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// check amp for selected channels in chunk 0
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	// check selected spw, chans in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_EQ(expSpw, cache->spw(chunk));
		ASSERT_TRUE(allEQ(expChans, cache->chan(chunk)));
	}
}

TEST_F( PlotMSCacheTest, testTimeSelection) {
	/* select time, check nchunk, amp, and first time */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String timeExpr("10:22:00~10:47:00");
	Int expNRow(378), expNChunk(23), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		timeExpr, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	Double expTime(msmc.time().get(0));
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setTimerange(timeExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	// check amp for selected field in chunk 0
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	ASSERT_EQ(expTime, cache->time(ichunk));
}

TEST_F( PlotMSCacheTest, testUVSelection) {
	/* select uvdist, check nchunk, values > min */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String uvDistExpr(">1km");
	// number of chunks reduced from 60
	Int expNChunk(16);
	Double minUVDist(1000.0);

	// set PlotMSSelection and load cache
	itsSelection.setUvrange(uvDistExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::UVDIST, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	// check all uVDist > minUVDist
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_TRUE(allGT(cache->uVDist(chunk), minUVDist));
	}
}

TEST_F( PlotMSCacheTest, testAntennaSelection) {
	/* select antenna baseline, check ant1 and ant2 */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String antennaExpr("1&2"); 
	Int expNRow(1), // one baseline per chunk
		expNChunk(60), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, antennaExpr, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Vector<Int> expAnt1(msmc.antenna1().getColumn()),
		expAnt2(msmc.antenna2().getColumn());

	// set PlotMSSelection and load cache
	itsSelection.setAntenna(antennaExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::ANTENNA1, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// check selected antenna in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_TRUE(allEQ(cache->ant1(chunk), expAnt1(chunk)));
		ASSERT_TRUE(allEQ(cache->ant2(chunk), expAnt2(chunk)));
	}
}

TEST_F( PlotMSCacheTest, testScanSelection) {
	/* select scan, check amp, scan numbers */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String scanExpr("2"); 
	Int expNRow(378), expNChunk(5), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, scanExpr,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Int expScan(msmc.scanNumber().get(0));
	Array<Complex> visData(msmc.data().getColumn());
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setScan(scanExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// check selected scan in chunks
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_EQ(cache->scan(chunk), expScan);
	}
}

TEST_F( PlotMSCacheTest, testCorrSelection) {
	// select correlation and check shape of chunks, amp, corr
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	String polnExpr("LL");
	Int expNChunk(60), expNCorr(1), expCorr(8), expNRow(351), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, polnExpr, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	// slice for corrSlices, nrows in first chunk
	Slicer visSlicer(corrSlices(0)(0), Slice(), Slice(0,expNRow));
	Array<Complex> selVis(visData(visSlicer));
	Array<Float> expAmp(amplitude(selVis));

	// set PlotMSSelection and load cache
	itsSelection.setCorr(polnExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNCorr, cache->chunkShapes()(IPosition(2,0,ichunk)));
	// check amp in chunk 0
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	// check selected corr in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_TRUE(allEQ(expCorr, cache->corr(chunk)));
	}
}

TEST_F( PlotMSCacheTest, testArraySelection) {
	/* select arrayID, check nchunk, nrow, amp
	 * (no cache getter for array id) */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String arrayExpr("0"); // only one array ID in MS, not really a selection
	Int expNRow(351), expNChunk(60), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		arrayExpr, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setArray(arrayExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr
	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
}

TEST_F( PlotMSCacheTest, testObservationSelection) {
	/* select observationID, check nchunk, nrow, amp, obsid */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	String obsExpr("0"); // only one obs ID in MS, not really a selection
	Int expNRow(351), expNChunk(60), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, obsExpr, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	Int expObs(msmc.observationId().get(0));
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setObservation(obsExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr
	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	// check selected obs id in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_TRUE(allEQ(expObs, cache->obsid(chunk)));
	}
}

TEST_F( PlotMSCacheTest, testIntentSelection) {
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	// no STATE table, STATE_ID col is -1
	// so this will throw an exception
	itsSelection.setIntent("CALIBRATE*");
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	ASSERT_THROW(cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ), AipsError); // ThreadCommunication* is nullptr
}

TEST_F( PlotMSCacheTest, testFeedSelection) {
	/* select feedID, check nchunk, nrow, amp, feed id */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// only one feed ID in MS, not really a selection
	// Use "auto-correlation" syntax to select feed1==feed2==0
	String feedExpr("0&&&");
	Int expNRow(351), expNChunk(60), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, ""/*taQLExpr*/, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, feedExpr);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Int expFeed(msmc.feed1().get(0));

	// set PlotMSSelection and load cache
	itsSelection.setFeed(feedExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::FEED1, PMS::FEED2};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr
	ASSERT_EQ(expNChunk, cache->nChunk());
	ASSERT_EQ(expNRow, cache->chunkShapes()(IPosition(2,2,ichunk)));
	// check selected feed id in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_TRUE(allEQ(expFeed, cache->feed1(chunk)));
		ASSERT_TRUE(allEQ(expFeed, cache->feed2(chunk)));
	}
}

TEST_F( PlotMSCacheTest, testTaQLSelection) {
	/* select field using taQL, check nchunk, amp and field */
	ASSERT_NE(nullptr, cache);  // make sure we have a cache
	String taQLExpr("FIELD_ID >= 1");
	Int expNRow(378), expNChunk(48), ichunk(0);

	// get selected MS
	MeasurementSet selMS;
	Vector<Vector<Slice> > chanSlices, corrSlices;
	mssSetData2(sortedMS, selMS, chanSlices, corrSlices,	""/*outMSName*/,
		""/*timeExpr*/, ""/*antennaExpr*/, ""/*fieldExpr*/, ""/*spwExpr*/,
		""/*uvDistExpr*/, taQLExpr, ""/*polnExpr*/, ""/*scanExpr*/,
		""/*arrayExpr*/, ""/*stateExpr*/, ""/*obsExpr*/, ""/*feedExpr*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn());
	Int expField(msmc.fieldId().get(0));
	// resize for nrows
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> expAmp(amplitude(visData));

	// set PlotMSSelection and load cache
	itsSelection.setMsselect(taQLExpr);	
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr

	ASSERT_EQ(expNChunk, cache->nChunk());
	// check amp for selected field in chunk 0
	ASSERT_TRUE(allEQ(expAmp, cache->ampData(ichunk)));
	ASSERT_EQ(expField, cache->field(ichunk));
	// check selected field in every chunk
	for (Int chunk=0; chunk<expNChunk; ++chunk) {
		ASSERT_GE(cache->field(chunk), expField);
	}
}

