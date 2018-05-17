//# PlotMSCacheTransform_GT.cc:: GoogleTest for plotms transformations
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
#include <casa/Arrays/ArrayMath.h>
#include <measures/Measures/MeasTable.h>
#include <ms/MeasurementSets/MSColumns.h>

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
	}

	String dataPath;
	MeasurementSet sortedMS;
	// use defaults (none!)
	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;
};


int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

TEST_F( PlotMSCacheTest, testFreqFrame) {
	// Test frequency conversion transformations
	// Load values from MeasurementSet main table
	ROMSMainColumns msmc(sortedMS);
	Int ddID(msmc.dataDescId().get(0)); // for spw / freq
	Int fieldID(msmc.fieldId().get(0)); // for phasedir
	Int obsID(msmc.observationId().get(0)); // for phasedir
	MEpoch epoch(msmc.timeMeas()(0));
	// values from subtables for freq conversion
	ROMSColumns mscol(sortedMS);
	Int spw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> inputFreq(mscol.spectralWindow().chanFreq()(spw));
	Int refFrame(mscol.spectralWindow().measFreqRef()(spw));
	Vector<MDirection> phaseDirCol(mscol.field().phaseDirMeasCol()(fieldID));
	MDirection phaseDir(phaseDirCol(0));
	// get observatory position
	String telName(mscol.observation().telescopeName().get(obsID));
	MPosition obsPosition;
	if (!MeasTable::Observatory(obsPosition, telName))
		obsPosition = mscol.antenna().positionMeas()(0);
	// make freq converter
	MeasFrame measFrame(epoch, obsPosition, phaseDir);
	MFrequency::Ref observedFrame(refFrame, measFrame);
	Unit unit(String("Hz"));

	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::Axis> loadAxes{PMS::FREQUENCY, PMS::CHANNEL};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};
	std::vector<casacore::String> frames{"LSRK",
		"LSRD", "BARY", "GEO", "TOPO", "GALACTO", "LGROUP", "CMB"};
	Int ichunk(0);  // check values for first chunk
	uInt nchan(inputFreq.size());
	MFrequency::Convert freqNewFrame;
	Vector<Double> expFreq(nchan);
	for (auto frame : frames) {
		// transform MS freqs with converter (expected frequencies)
		MFrequency::Types newFrameType(MFrequency::typeFromString(frame));
		freqNewFrame = MFrequency::Convert(unit, observedFrame, newFrameType);
		for (uInt chan=0; chan<nchan; ++chan) 
			expFreq(chan) = freqNewFrame(inputFreq(chan)).getValue();
		expFreq /= 1.0e9;  // convert to GHz

		// get cache value
		// Set up plotms cache object each time else uses cached values for same axis
		MSCache* cache = new MSCache(nullptr);
		itsTransformations.setFrame(frame);
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		// check expected vs cache freqs
		ASSERT_TRUE(allEQ(expFreq, cache->freq(ichunk)));
		delete cache;
	}
}

TEST_F( PlotMSCacheTest, testRestFreq) {
	// Test velocity rest frequency options
	// Load values from MeasurementSet main table
	ROMSMainColumns msmc(sortedMS);
	Int ddID(msmc.dataDescId().get(0)); // for spw / freq
	// values from subtables for freq conversion
	ROMSColumns mscol(sortedMS);
	Int spw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> inputFreq(mscol.spectralWindow().chanFreq()(spw));
	Vector<MFrequency> inputFreqMeas(mscol.spectralWindow().chanFreqMeas()(spw));
	// get velocity with various rest freqs (test only, not a science case)
	uInt nchan(inputFreq.size());
	Double midrestfreq(inputFreq(nchan/2)),
		   lowrestfreq(min(inputFreq)), highrestfreq(max(inputFreq));

	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::Axis> loadAxes{PMS::VELOCITY, PMS::CHANNEL};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};
	std::vector<Double> testfreqs{lowrestfreq, midrestfreq, highrestfreq};
	Int ichunk(0);  // check values for first chunk
	Vector<Double> expVel(nchan); // converted from frequency
	for (auto restfreq : testfreqs) {
		// get expected velocities
		for (uInt chan=0; chan<nchan; ++chan) {
			// Measures: .get(Unit) returns Quantity (value+units), then .getValue()
			expVel(chan) = casacore::MDoppler::Convert(inputFreqMeas(chan).toDoppler(restfreq))().get("km/s").getValue();
		}

		// get cache value
		// Set up plotms cache object each time else uses cached values for same axis
		MSCache* cache = new MSCache(nullptr);
		itsTransformations.setRestFreq(restfreq / 1.0e6); // in MHz
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		// check expected vs cache freqs
		ASSERT_TRUE(allNear(expVel, cache->vel(ichunk), .0000001));
		delete cache;
	}
}

TEST_F( PlotMSCacheTest, testVelDef) {
	// Test velocity definition options
	// Load values from MeasurementSet main table
	ROMSMainColumns msmc(sortedMS);
	Int ddID(msmc.dataDescId().get(0)); // for spw / freq
	// values from subtables for freq conversion
	ROMSColumns mscol(sortedMS);
	Int spw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> inputFreq(mscol.spectralWindow().chanFreq()(spw));
	Vector<MFrequency> inputFreqMeas(mscol.spectralWindow().chanFreqMeas()(spw));
	// get velocity with various rest freqs (test only, not a science case)
	uInt nchan(inputFreq.size());
	Double restfreq(inputFreq(nchan/2));

	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::Axis> loadAxes{PMS::VELOCITY, PMS::CHANNEL};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};
	std::vector<casacore::String> veldefs{"RADIO", "OPTICAL", "TRUE"};
	Int ichunk(0);  // check values for first chunk
	Vector<Double> expVel(nchan); // converted from frequency
	MDoppler::Types dopptype; // converted from string
	for (auto veldef : veldefs) {
		// get expected velocities
		MDoppler::getType(dopptype, veldef);
		for (uInt chan=0; chan<nchan; ++chan) {
			// Measures: .get(Unit) returns Quantity (value+units), then .getValue()
			expVel(chan) = casacore::MDoppler::Convert(inputFreqMeas(chan).toDoppler(restfreq), dopptype)().get("km/s").getValue();
		}

		// get cache value
		// Set up plotms cache object each time else uses cached values for same axis
		MSCache* cache = new MSCache(nullptr);
		itsTransformations.setVelDef(veldef);
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		// check expected vs cache freqs
		ASSERT_TRUE(allNear(expVel, cache->vel(ichunk), .0000001));
		delete cache;
	}
}

TEST_F( PlotMSCacheTest, testPhaseShift) {
	// Test phase shift transformations
	// Load values from MeasurementSet main table
	ROMSMainColumns msmc(sortedMS);
	Matrix<Double> uvw(msmc.uvw().getColumn());
	Cube<Complex> visData(msmc.data().getColumn());
	Int ddID(msmc.dataDescId().get(0)); // for spw / freq
	// values from subtables for freq conversion
	ROMSColumns mscol(sortedMS);
	Int spw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> inputFreq(mscol.spectralWindow().chanFreq()(spw));

	// resize expected nrow in first chunk
	uInt expNRow(351);
	uvw.resize(IPosition(2,3,expNRow), True);
	Vector<Double> vecU(uvw.row(0)), vecV(uvw.row(1));
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);

	Double dx(0.5), dy(0.5);  // phase shift values
	Double toRadians = C::pi / 180.0 / 3600.0;
	Double dxrad(dx*toRadians), dyrad(dy*toRadians); // radians
	vecU *= dxrad;
	vecV *= dyrad;
	Vector<Double> phases(vecU + vecV);
	phases *= (-2.0 * C::pi / C::c); // radians/Hz 
	uInt ncorr(visData.shape()(0));
	Double phaRad;
	Complex factor;
	for (uInt row=0; row<expNRow; ++row) {
		for (uInt chan=0; chan<inputFreq.size(); ++chan) {
			phaRad = phases(row) * inputFreq(chan);
			if (phaRad != 0.0) {
				factor = Complex(cos(phaRad), sin(phaRad));
				for (uInt corr=0; corr<ncorr; ++corr) {
					visData(corr, chan, row) *= factor;
				}
			}
		}
	}
	// convert to degrees
	Cube<Float> expPhase(phase(visData) * 180.0 / C::pi);

	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::Axis> loadAxes{PMS::PHASE, PMS::TIME};
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};
	Int ichunk(0);  // check values for first chunk

	// get cache value
	itsTransformations.setXpcOffset(dx);
	itsTransformations.setYpcOffset(dy);
	MSCache* cache = new MSCache(nullptr);
	cache->load(loadAxes, loadData, dataPath, 
		itsSelection, itsAveraging,
		itsTransformations, itsCalibration,
		nullptr ); // ThreadCommunication* is nullptr 

	// check expected vs cache phase
	ASSERT_TRUE(allEQ(expPhase, cache->phaData(ichunk)));
	delete cache;
}
