//# Copyright (C) 2008
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

#include <plotms/test/tUtil.h>
#include <plotms/Data/MSCache.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <ms/MeasurementSets/MSColumns.h>

#include <iostream>

// Tests whether simple plot (metadata axes) can be loaded
using namespace casa;

int main(int /*argc*/, char** /*argv[]*/) {
	String dataPath = tUtil::getFullPath( "ngc5921_add_corect_model.ms", "visstat2" );
	Int expNRow(351);

	MeasurementSet ms(dataPath);
	Block<String> sortCols(4); // Use default sort columns
	sortCols[0] = MS::columnName(MS::ARRAY_ID);
	sortCols[1] = MS::columnName(MS::FIELD_ID);
	sortCols[2] = MS::columnName(MS::DATA_DESC_ID);
	sortCols[3] = MS::columnName(MS::TIME);
	MeasurementSet sortedMS = ms.sort(sortCols);
	ROMSMainColumns msmc(sortedMS);
	Array<Complex> visData(msmc.data().getColumn());
	IPosition visShape(visData.shape());
	visShape.setLast(IPosition(1,expNRow));
	visData.adjustLastAxis(visShape);
	Array<Float> amp(amplitude(visData));

	PlotMSSelection itsSelection;
	PlotMSAveraging itsAveraging;
	PlotMSTransformations itsTransformations;
	PlotMSCalibration itsCalibration;
	MSCache* cache = new MSCache(nullptr);

	// check values for first chunk
	Int ichunk(0);
	std::vector<PMS::Axis> loadAxes{PMS::AMP, PMS::TIME};
	std::vector<PMS::DataColumn> loadData{PMS::DATA, PMS::DATA};
	cache->load(loadAxes, loadData, dataPath, 
			itsSelection, itsAveraging,
			itsTransformations, itsCalibration, 
			nullptr ); // ThreadCommunication* is nullptr
	cout << "chunk:" << cache->nChunk() << endl;
	cout << "nrow:" << cache->chunkShapes()(IPosition(2,2,ichunk)) << endl;
	cout << "exp=got:" << allEQ(amp, cache->ampData(ichunk)) << endl;
	bool cacheOk = (cache != nullptr);
	delete cache;
	return cacheOk;    
}

