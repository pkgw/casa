//# PlotMSCacheGeom_GT.cc:: GoogleTest for plotms cache observational geometry axes
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

TEST_F( PlotMSCacheTest, testObsGeom) {
	// Test observational geometry axis options:
	//  UVDIST,UVDIST_L,U,V,W,UWAVE,VWAVE,WWAVE,AZ0,EL0,HA0,PA0,
	//  ANTENNA,AZIMUTH,ELEVATION,PARANG
	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	// Load expected/needed values from MeasurementSet
	// main table columns
	ROMSMainColumns msmc(sortedMS);
	Matrix<Double> uvw(msmc.uvw().getColumn());
	Int ddID(msmc.dataDescId().get(0)),
		fieldId(msmc.fieldId().get(0));
	MEpoch epoch(msmc.timeMeas()(0));
	// subtable columns
	ROMSColumns mscol(sortedMS);
	Int spw(mscol.dataDescription().spectralWindowId().get(ddID));
	Vector<Double> chanfreqs(mscol.spectralWindow().chanFreq().get(spw));
  	uInt nAnt = mscol.antenna().nrow();
	Vector<MDirection> phasedir(mscol.field().phaseDirMeasCol()(fieldId));
	Cube<Double> recepAngle(mscol.feed().receptorAngle().getColumn());
	Matrix<Double> angles = recepAngle.xyPlane(0);

	MSDerivedValues msd;  // for azel
	msd.setAntennas(mscol.antenna());
	msd.setFieldCenter(phasedir(0));
	msd.setEpoch(epoch);

	// Resize for first chunk and do calculations
	// get u,v,w
	uvw.resize(IPosition(2,3,expNRow), True);
	Vector<Double> expU(uvw.row(0)), expV(uvw.row(1)), expW(uvw.row(2)),
		expUVDist;
	// calculate uvdist
	uInt vecsize = expU.size();
	expUVDist.resize(vecsize);
	for (uInt i=0; i<expU.size(); ++i) {
		Double u(expU(i)), v(expV(i));
		expUVDist(i) = sqrt(u*u + v*v);
	}
	// calculate uvdistl, uwave, vwave, wwave
	Matrix<Double> expUVDistL, expUWave, expVWave, expWWave;
	// c is the speed of light in km/s (casa/BasicSL/Constants.cc):
	Vector<Double> uvDistM(expUVDist/C::c);
	Vector<Double> uM(expU/C::c), vM(expV/C::c), wM(expW/C::c);
	expUVDistL.resize(expNChan, expNRow);
	expUWave.resize(expNChan, expNRow);
	expVWave.resize(expNChan, expNRow);
	expWWave.resize(expNChan, expNRow);
	for (Int row=0; row<expNRow; ++row) {
		expUVDistL.column(row) = uvDistM(row) * chanfreqs;
		expUWave.column(row) = uM(row) * chanfreqs;
		expVWave.column(row) = vM(row) * chanfreqs;
		expWWave.column(row) = wM(row) * chanfreqs;
	}
	// calculate per-antenna az and el
	Vector<Double> azel, expAz(nAnt), expEl(nAnt);
	Vector<Float> expParAng(nAnt);
	for (uInt ant=0; ant<nAnt; ++ ant) {
		msd.setAntenna(ant);
		MDirection azeldir = msd.azel();
		azel = azeldir.getAngle("deg").getValue();
		expAz(ant) = azel(0);
		expEl(ant) = azel(1);
		Float antParAng = msd.parAngle();
		antParAng += angles(0,ant);
		antParAng *= 180.0/C::pi; // degrees
		expParAng(ant) = antParAng;
	}
	msd.setAntenna(-1); // use observatory position
	MDirection azel0 = msd.azel();
	azel = azel0.getAngle("deg").getValue();
	Double expAz0(azel(0)), expEl0(azel(1));
	Double expHA0(msd.hourAngle()*12/C::pi), // convert to hours
		   expPA0(msd.parAngle()*180.0/C::pi); // convert to degrees
	if (expPA0 < 0.0) expPA0 += 360.00;  // positive degrees

	// make antenna vector
	Vector<Int> expAnt(nAnt);
	indgen(expAnt);
	
	// Request geometry axes:
	std::vector<PMS::Axis> geomAxes {PMS::UVDIST, PMS::UVDIST_L, PMS::U, PMS::V,
		PMS::W, PMS::UWAVE, PMS::VWAVE, PMS::WWAVE, PMS::AZ0, PMS::EL0, PMS::HA0,
		PMS::PA0, PMS::ANTENNA, PMS::AZIMUTH, PMS::ELEVATION, PMS::PARANG};
	// Datacolumn for non-vis axes always DATA
	std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};

	// check values for first chunk
	Int ichunk(0);
	for (auto axis : geomAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration,
					nullptr ); // ThreadCommunication* is nullptr 

		// Check arrays with allEQ in ArrayLogical.h
		switch(axis) {
			case PMS::UVDIST:
				ASSERT_TRUE(allEQ(expUVDist, cache->uVDist(ichunk)));
				break;
			case PMS::UVDIST_L:
				ASSERT_TRUE(allEQ(expUVDistL, cache->uVDistL(ichunk)));
				break;
			case PMS::U:
				ASSERT_TRUE(allEQ(expU, cache->u(ichunk)));
				break;
			case PMS::V:
				ASSERT_TRUE(allEQ(expV, cache->v(ichunk)));
				break;
			case PMS::W:
				ASSERT_TRUE(allEQ(expW, cache->w(ichunk)));
				break;
			case PMS::UWAVE:
				ASSERT_TRUE(allEQ(expUWave, cache->uWave(ichunk)));
				break;
			case PMS::VWAVE:
				ASSERT_TRUE(allEQ(expVWave, cache->vWave(ichunk)));
				break;
			case PMS::WWAVE:
				ASSERT_TRUE(allEQ(expWWave, cache->wWave(ichunk)));
				break;
			case PMS::AZ0:
				ASSERT_EQ(expAz0, cache->az0(ichunk));
				break;
			case PMS::EL0:
				ASSERT_EQ(expEl0, cache->el0(ichunk));
				break;
			case PMS::HA0:
				ASSERT_NEAR(expHA0, cache->ha0(ichunk), .00001);
				break;
			case PMS::PA0:
				ASSERT_NEAR(expPA0, cache->pa0(ichunk), .00001);
				break;
			case PMS::ANTENNA:
				ASSERT_TRUE(allEQ(expAnt, cache->ant(ichunk)));
				break;
			case PMS::AZIMUTH:
				ASSERT_TRUE(allEQ(expAz, cache->az(ichunk)));
				break;
			case PMS::ELEVATION:
				ASSERT_TRUE(allEQ(expEl, cache->el(ichunk)));
				break;
			case PMS::PARANG:
				ASSERT_TRUE(allNear(expParAng, cache->parAng(ichunk), .00001));
				break;
			default:
				break;
		}
	}
}

