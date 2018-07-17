//# PlotMSCacheVis_GT.cc:: GoogleTest for plotms cache visibility axes
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
		expNChunk = 12;
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

TEST_F( PlotMSCacheTest, testVisibilities) {
	// Test visibility and flag axis options:
	// AMP,PHASE,REAL,IMAG,WT,WTxAMP,WTSP,
	// SIGMA,SIGMASP,FLAG,FLAG_ROW
	// and datacolumn options:
	// DATA, CORRECTED, MODEL, CORRMODEL_V, corrMODEL_S,
	// DATAMODEL_V, DATAMODEL_S, CORR_DIV_MODEL_V, CORR_DIV_MODEL_S,
	// DATA_DIV_MODEL_V, DATA_DIV_MODEL_S

	ASSERT_NE(nullptr, cache);  // make sure we have a cache

	String scanExpr("1"); // add selection to reduce test time
	itsSelection.setScan(scanExpr);
	// Get selected MS:
	MeasurementSet selMS;
	mssSetData2(sortedMS, selMS, ""/*outms*/, ""/*time*/, ""/*ant*/,
		""/*field*/, ""/*spw*/, ""/*uvdist*/, ""/*taql*/, ""/*poln*/,
		scanExpr, ""/*array*/, ""/*state*/, ""/*obs*/, ""/*feed*/);
	// main table columns
	ROMSMainColumns msmc(selMS);
	Array<Complex> visData(msmc.data().getColumn()),
		modelData(msmc.modelData().getColumn()),
		corrData(msmc.correctedData().getColumn());
	Array<Bool> expFlag(msmc.flag().getColumn());
	Vector<Bool> expFlagrow(msmc.flagRow().getColumn());
	Array<Float> expWt(msmc.weight().getColumn()),
		expWtSp(msmc.weightSpectrum().getColumn()),
		expSigma(msmc.sigma().getColumn()),
		expSigmaSp;
	// adjust cube shapes for first chunk
	IPosition flagShape(expFlag.shape());
	flagShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(flagShape);
	modelData.adjustLastAxis(flagShape);
	corrData.adjustLastAxis(flagShape);
	expFlag.adjustLastAxis(flagShape);
	expWtSp.adjustLastAxis(flagShape);
	// adjust wt arrays
	IPosition wtShape(expWt.shape());
	wtShape.setLast(IPosition(1,expNRow));
	expWt.adjustLastAxis(wtShape);
	expSigma.adjustLastAxis(wtShape);
	// adjust vector
	expFlagrow.resize(expNRow, True);
	// make sigmasp (no column) by adding chan axis to sigma
	tUtil::makeSigmaSpFromSigma(expSigmaSp, expSigma, expNChan);

	// make residual data cubes
	Array<Complex> corrmodelData(corrData - modelData);
	Array<Complex> vismodelData(visData - modelData);
	Array<Complex> corrDivmodelData(corrData / modelData);
	Array<Complex> visDivmodelData(visData / modelData);

	// Visibility axis options:
	std::vector<PMS::Axis> visAxes {PMS::AMP, PMS::PHASE,
		PMS::REAL, PMS::IMAG, PMS::WTxAMP};
	// Datacolumn options:
	std::vector<PMS::DataColumn> visCols {PMS::DATA, PMS::CORRECTED,
		PMS::MODEL, PMS::CORRMODEL_V, PMS::CORRMODEL_S,
		PMS::DATAMODEL_V, PMS::DATAMODEL_S,
		PMS::CORR_DIV_MODEL_V, PMS::CORR_DIV_MODEL_S,
		PMS::DATA_DIV_MODEL_V, PMS::DATA_DIV_MODEL_S};

	// check values for first chunk: each vis axis with each datacolumn
	// Expected values:
	Array<Float> expVis;  // datacube with axis operation applied
	Int ichunk(0);
	for (auto axis : visAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		for (auto datacol : visCols) {
			std::vector<PMS::DataColumn> loadData{datacol, PMS::DATA};
			cache->load(loadAxes, loadData, dataPath, 
					itsSelection, itsAveraging,
					itsTransformations, itsCalibration, 
					nullptr ); // ThreadCommunication* is nullptr

			ASSERT_EQ(expNChunk, cache->nChunk());
			switch (axis) {
				case PMS::AMP: {
					switch (datacol) {
						case PMS::DATA: {
							expVis = amplitude(visData);
							ASSERT_TRUE(allEQ(expVis, cache->ampData(ichunk)));
							break;
						}
						case PMS::CORRECTED: {
							expVis = amplitude(corrData);
							ASSERT_TRUE(allEQ(expVis, cache->ampCorr(ichunk)));
							break;
						}
						case PMS::MODEL: {
							expVis = amplitude(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->ampModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_V: {
							expVis = amplitude(corrmodelData);
							ASSERT_TRUE(allEQ(expVis, cache->ampCorrModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_S: {
							expVis = amplitude(corrData) - amplitude(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->ampCorrModelS(ichunk)));
							break;
						}
						case PMS::DATAMODEL_V: {
							expVis = amplitude(vismodelData);
							ASSERT_TRUE(allEQ(expVis, cache->ampDataModel(ichunk)));
							break;
						}
						case PMS::DATAMODEL_S: {
							expVis = amplitude(visData) - amplitude(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->ampDataModelS(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_V: {
							expVis = amplitude(corrDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->ampCorrDivModel(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_S: {
							expVis = amplitude(corrData) / amplitude(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->ampCorrDivModelS(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_V: {
							expVis = amplitude(visDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->ampDataDivModel(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_S: {
							expVis = amplitude(visData) / amplitude(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->ampDataDivModelS(ichunk)));
							break;
						}
						case PMS::FLOAT_DATA:
							break;
					}
					break;
				}
				case PMS::PHASE: { // convert to degrees
					switch (datacol) {
						case PMS::DATA: {
							expVis = phase(visData) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaData(ichunk)));
							break;
						}
						case PMS::CORRECTED: {
							expVis = phase(corrData) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaCorr(ichunk)));
							break;
						}
						case PMS::MODEL: {
							expVis = phase(modelData) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_V: {
							expVis = phase(corrmodelData) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaCorrModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_S: {
							expVis = (phase(corrData) - phase(modelData)) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaCorrModelS(ichunk)));
							break;
						}
						case PMS::DATAMODEL_V: {
							expVis = phase(vismodelData) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaDataModel(ichunk)));
							break;
						}
						case PMS::DATAMODEL_S: {
							expVis = (phase(visData) - phase(modelData)) * 180.0 / C::pi;
							ASSERT_TRUE(allEQ(expVis, cache->phaDataModelS(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_V: {
							expVis = phase(corrDivmodelData) * 180.0 / C::pi;
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->phaCorrDivModel(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_S: {
							expVis = (phase(corrData) / phase(modelData)) * 180.0 / C::pi;
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->phaCorrDivModelS(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_V: {
							expVis = phase(visDivmodelData) * 180.0 / C::pi;
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->phaDataDivModel(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_S: {
							expVis = (phase(visData) / phase(modelData)) * 180.0 / C::pi;
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->phaDataDivModelS(ichunk)));
							break;
						}
						case PMS::FLOAT_DATA:
							break;
					}
					break;
				}
				case PMS::REAL: {
					switch (datacol) {
						case PMS::DATA: {
							expVis = real(visData);
							ASSERT_TRUE(allEQ(expVis, cache->realData(ichunk)));
							break;
						}
						case PMS::CORRECTED: {
							expVis = real(corrData);
							ASSERT_TRUE(allEQ(expVis, cache->realCorr(ichunk)));
							break;
						}
						case PMS::MODEL: {
							expVis = real(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->realModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_V: {
							expVis = real(corrmodelData);
							ASSERT_TRUE(allEQ(expVis, cache->realCorrModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_S: {
							expVis = real(corrData) - real(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->realCorrModelS(ichunk)));
							break;
						}
						case PMS::DATAMODEL_V: {
							expVis = real(vismodelData);
							ASSERT_TRUE(allEQ(expVis, cache->realDataModel(ichunk)));
							break;
						}
						case PMS::DATAMODEL_S: {
							expVis = real(visData) - real(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->realDataModelS(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_V: {
							expVis = real(corrDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->realCorrDivModel(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_S: {
							expVis = real(corrData) / real(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->realCorrDivModelS(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_V: {
							expVis = real(visDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->realDataDivModel(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_S: {
							expVis = real(visData) / real(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->realDataDivModelS(ichunk)));
							break;
						}
						case PMS::FLOAT_DATA:
							break;
					}
					break;
				}
				case PMS::IMAG: {
					switch (datacol) {
						case PMS::DATA: {
							expVis = imag(visData);
							ASSERT_TRUE(allEQ(expVis, cache->imagData(ichunk)));
							break;
						}
						case PMS::CORRECTED: {
							expVis = imag(corrData);
							ASSERT_TRUE(allEQ(expVis, cache->imagCorr(ichunk)));
							break;
						}
						case PMS::MODEL: {
							expVis = imag(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->imagModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_V: {
							expVis = imag(corrmodelData);
							ASSERT_TRUE(allEQ(expVis, cache->imagCorrModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_S: {
							expVis = imag(corrData) - imag(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->imagCorrModelS(ichunk)));
							break;
						}
						case PMS::DATAMODEL_V: {
							expVis = imag(vismodelData);
							ASSERT_TRUE(allEQ(expVis, cache->imagDataModel(ichunk)));
							break;
						}
						case PMS::DATAMODEL_S: {
							expVis = imag(visData) - imag(modelData);
							ASSERT_TRUE(allEQ(expVis, cache->imagDataModelS(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_V: {
							expVis = imag(corrDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->imagCorrDivModel(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_S: {
							expVis = imag(corrData) / imag(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->imagCorrDivModelS(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_V: {
							expVis = imag(visDivmodelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->imagDataDivModel(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_S: {
							expVis = imag(visData) / imag(modelData);
							ASSERT_TRUE(tUtil::allEQDiv(expVis, cache->imagDataDivModelS(ichunk)));
							break;
						}
						case PMS::FLOAT_DATA:
							break;
					}
					break;
				}
				case PMS::WTxAMP: {
					switch (datacol) {
						case PMS::DATA: {
							expVis = amplitude(visData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampData(ichunk)));
							break;
						}
						case PMS::CORRECTED: {
							expVis = amplitude(corrData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampCorr(ichunk)));
							break;
						}
						case PMS::MODEL: {
							expVis = amplitude(modelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_V: {
							expVis = amplitude(corrmodelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampCorrModel(ichunk)));
							break;
						}
						case PMS::CORRMODEL_S: {
							expVis = amplitude(corrData) - amplitude(modelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampCorrModelS(ichunk)));
							break;
						}
						case PMS::DATAMODEL_V: {
							expVis = amplitude(vismodelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampDataModel(ichunk)));
							break;
						}
						case PMS::DATAMODEL_S: {
							expVis = amplitude(visData) - amplitude(modelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(allEQ(expWtAmp, cache->wtxampDataModelS(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_V: {
							expVis = amplitude(corrDivmodelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(tUtil::allEQDiv(expWtAmp, cache->wtxampCorrDivModel(ichunk)));
							break;
						}
						case PMS::CORR_DIV_MODEL_S: {
							expVis = amplitude(corrData) / amplitude(modelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(tUtil::allEQDiv(expWtAmp, cache->wtxampCorrDivModelS(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_V: {
							expVis = amplitude(visDivmodelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(tUtil::allEQDiv(expWtAmp, cache->wtxampDataDivModel(ichunk)));
							break;
						}
						case PMS::DATA_DIV_MODEL_S: {
							expVis = amplitude(visData) / amplitude(modelData);
							Array<Float> expWtAmp = tUtil::getWtAmp(expWt, expVis);
							ASSERT_TRUE(tUtil::allEQDiv(expWtAmp, cache->wtxampDataDivModelS(ichunk)));
							break;
						}
						case PMS::FLOAT_DATA:
							break;
					}
					break;
				}
				default:
					break;
			}
		}
	}
	ASSERT_TRUE(allEQ(expFlag, cache->flag(ichunk)));

	// check wt, wtsp, sigma, sigmasp, flagrow 
	std::vector<PMS::Axis> wtAxes {PMS::WT, PMS::WTSP,
		PMS::SIGMA, PMS::SIGMASP, PMS::FLAG_ROW};
	for (auto axis : wtAxes) {
		std::vector<PMS::Axis> loadAxes{axis, PMS::TIME};
		std::vector<PMS::DataColumn> loadData {PMS::DATA, PMS::DATA};
		cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr 
		switch(axis) {
			case PMS::WT:
				ASSERT_TRUE(allEQ(expWt, cache->wt(ichunk)));
				break;
			case PMS::WTSP:
				ASSERT_TRUE(allEQ(expWtSp, cache->wtsp(ichunk)));
				break;
			case PMS::SIGMA:
				ASSERT_TRUE(allEQ(expSigma, cache->sigma(ichunk)));
				break;
			case PMS::SIGMASP:
				ASSERT_TRUE(allEQ(expSigmaSp, cache->sigmasp(ichunk)));
				break;
			case PMS::FLAG_ROW:
				ASSERT_TRUE(allEQ(expFlagrow, cache->flagrow(ichunk)));
				break;
			default:
				break;
		}
	}
}


