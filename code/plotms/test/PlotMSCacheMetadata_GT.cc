//# PlotMSCacheMetadata_GT.cc:: GoogleTest for plotms cache metadata axes
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
#include <casa/Arrays/ArrayMath.h>

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
		expNRow = 351;  // rows in first chunk
		expNChan = 63;
		// Set up plotms cache object with no parent (PlotMSApp*)
		cache = new MSCache(nullptr);
	}

	virtual void TearDown() {
		delete cache;
	}

	String dataPath;
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

TEST_F( PlotMSCacheTest, testMetadata) {
	// Test metadata axis options
	// SCAN,FIELD,TIME,TIME_INTERVAL,SPW,CHANNEL,FREQUENCY,VELOCITY,CORR,
	// ANTENNA1,ANTENNA2,BASELINE,ROW,OBSERVATION,INTENT,FEED1,FEED2

	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// Load expected values from MeasurementSet main table columns
	ROMSMainColumns msmc(sortedMS);
	Int expScan(msmc.scanNumber().get(0)),
		expField(msmc.fieldId().get(0));
	Double expTime(msmc.time().get(0)),
		expInterval(msmc.interval().get(0));
	Vector<Int> expAnt1(msmc.antenna1().getColumn()),
		expAnt2(msmc.antenna2().getColumn()),
		expObsId(msmc.observationId().getColumn()),
		expIntent(msmc.stateId().getColumn()),
		expFeed1(msmc.feed1().getColumn()),
		expFeed2(msmc.feed2().getColumn()),
		expBaseline;
	Int ddID(msmc.dataDescId().get(0)); // for spw, pol

	// subtable columns
	ROMSColumns mscol(sortedMS);
	Int polID(mscol.dataDescription().polarizationId().get(ddID));
	uInt nAnt(mscol.antenna().nrow());
	Vector<Int> expChan,
		expCorr(mscol.polarization().corrType().get(polID));
	Int expSpw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> expFreq(mscol.spectralWindow().chanFreq().get(expSpw)),
		expVel; // converted from frequency
	Vector<MFrequency> freqMeas(mscol.spectralWindow().chanFreqMeas()(expSpw));
	uInt nchan = expFreq.size();
	expVel.resize(nchan);
	expChan.resize(nchan);
	// get freq, chan, velocity
	Double restfreq = expFreq(nchan/2); // before freq conversion
	expFreq /= 1.0e9; // convert to GHz
	indgen(expChan);
	for (uInt chan=0; chan<nchan; ++chan) {
		// Measures: .get(Unit) returns Quantity (value+units), then .getValue()
		expVel(chan) = casacore::MDoppler::Convert(freqMeas(chan).toDoppler(restfreq))().get("km/s").getValue();
	}

	// Resize MS columns to expected rows in first chunk
	Vector<uInt> expRow;
	expRow.resize(expNRow);
	indgen(expRow);
	expAnt1.resize(expNRow, True);
	expAnt2.resize(expNRow, True);
	expObsId.resize(expNRow, True);
	expIntent.resize(expNRow, True);
	expFeed1.resize(expNRow, True);
	expFeed2.resize(expNRow, True);
	// get baseline numbers
	expBaseline.resize(expNRow);
	for (Int i=0; i<expNRow; ++i) {
		Int ant1(expAnt1(i)), ant2(expAnt2(i));
		expBaseline(i) = (nAnt+1)*ant1-(ant1*(ant1+1))/2+ant2;
	}

	// Automatically loaded metadata axes: time, field, spw, chan, freq,
	// corr, scan, ant1, ant2, baseline, flag, obs, intent
	// Request other metadata axes:
	std::vector<PMS::Axis> metaAxes {PMS::TIME_INTERVAL, PMS::VELOCITY,
		PMS::ROW, PMS::FEED1, PMS::FEED2};
	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};

	// check values for first chunk
	Int ichunk(0);
	for (auto axis : metaAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		ASSERT_EQ(expNChunk, cache->nChunk());
		Int nrow = cache->chunkShapes()(IPosition(2,2,ichunk));
		ASSERT_EQ(expNRow, nrow);
		// Check scalars with ASSERT_EQ
		// Check arrays with allEQ in ArrayLogical.h
		switch(axis) {
			case PMS::TIME_INTERVAL:
				ASSERT_EQ(expInterval, cache->timeIntr(ichunk));
				break;
			case PMS::VELOCITY:
				ASSERT_TRUE(allNear(expVel, cache->vel(ichunk), .0000001));
				break;
			case PMS::ROW:
				ASSERT_TRUE(allEQ(expRow, cache->row(ichunk)));
				break;
			case PMS::FEED1:
				ASSERT_TRUE(allEQ(expFeed1, cache->feed1(ichunk)));
				break;
			case PMS::FEED2:
				ASSERT_TRUE(allEQ(expFeed2, cache->feed2(ichunk)));
				break;
			default:
				break;
		}
	}
	ASSERT_EQ(expScan, cache->scan(ichunk));
	ASSERT_EQ(expField, cache->field(ichunk));
	ASSERT_EQ(expTime, cache->time(ichunk));
	ASSERT_EQ(expSpw, cache->spw(ichunk));
	ASSERT_TRUE(allEQ(expChan, cache->chan(ichunk)));
	ASSERT_TRUE(allEQ(expFreq, cache->freq(ichunk)));
	ASSERT_TRUE(allEQ(expCorr, cache->corr(ichunk)));
	ASSERT_TRUE(allEQ(expAnt1, cache->ant1(ichunk)));
	ASSERT_TRUE(allEQ(expAnt2, cache->ant2(ichunk)));
	ASSERT_TRUE(allEQ(expBaseline, cache->bsln(ichunk)));
	ASSERT_TRUE(allEQ(expObsId, cache->obsid(ichunk)));
	ASSERT_TRUE(allEQ(expIntent, cache->intent(ichunk)));
}

