//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA

#include <casacore/casa/Containers/ValueHolder.h>
#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/ms/MSOper/MSMetaData.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/TableProxy.h>
#include <casacore/tables/DataMan/TiledShapeStMan.h>

#include <mstransform/MSTransform/StatWt.h>
#include <mstransform/TVI/StatWtTVI.h>
#include <mstransform/TVI/StatWtTVILayerFactory.h>
#include <msvis/MSVis/ViImplementation2.h>
#include <msvis/MSVis/IteratingParameters.h>
#include <msvis/MSVis/LayeredVi2Factory.h>

using namespace casacore;

namespace casa { 

StatWt::StatWt(MeasurementSet* ms) : _ms(ms), _saf() {
    ThrowIf(! _ms, "Input MS pointer cannot be NULL");
}

StatWt::~StatWt() {}

void StatWt::setOutputMS(const casacore::String& outname) {
    _outname = outname;
}

void StatWt::setTimeBinWidth(const casacore::Quantity& binWidth) {
    ThrowIf(
        ! binWidth.isConform(Unit("s")),
        "Time bin width unit must be a unit of time"
    );
    setTimeBinWidth(binWidth.getValue("s"));
}

void StatWt::setTimeBinWidth(Double binWidth) {
    ThrowIf(
        binWidth <= 0, "time bin width must be positive"
    );
    _timeBinWidth = binWidth;
}

void StatWt::setTimeBinWidthUsingInterval(uInt n) {
    MSMetaData msmd(_ms, 100.0);
    auto stats = msmd.getIntervalStatistics();
    ThrowIf(
        stats.max/stats.median - 1 > 0.25
        || 1 - stats.min/stats.median > 0.25,
        "There is not a representative integration time in the INTERVAL column"
    );
    auto width = n*stats.median;
    _log << LogOrigin("StatWt", __func__) << LogIO::NORMAL
        << "Determined representative integration time of "
        << stats.median << "s. Setting time bin width to "
        << width << "s" << LogIO::POST;
    setTimeBinWidth(width);
}

void StatWt::setCombine(const String& combine) {
    _combine = downcase(combine);
}

void StatWt::setPreview(casacore::Bool preview) {
    _preview = preview;
}

void StatWt::setTVIConfig(const Record& config) {
    _tviConfig = config;
}

void StatWt::writeWeights() const {
    auto doWtSp = _ms->isColumn(MSMainEnums::WEIGHT_SPECTRUM);
    auto needWtSp = ! doWtSp && _tviConfig.isDefined(vi::StatWtTVI::CHANBIN);
    if (needWtSp) {
        doWtSp = True;
        // from Calibrater.cc
        // Nominal default tile shape
        IPosition dts(3, 4, 32, 1024);
        // Discern DATA's default tile shape and use it
        const auto dminfo = _ms->dataManagerInfo();
        for (uInt i=0; i<dminfo.nfields(); ++i) {
            Record col = dminfo.asRecord(i);
            if (anyEQ(col.asArrayString("COLUMNS"), String("DATA"))) {
                dts = IPosition(col.asRecord("SPEC").asArrayInt("DEFAULTTILESHAPE"));
                break;
             }
         }
         // Add the column
         String colWtSp = MS::columnName(MS::WEIGHT_SPECTRUM);
         TableDesc tdWtSp;
         tdWtSp.addColumn(ArrayColumnDesc<Float>(colWtSp, "weight spectrum", 2));
         TiledShapeStMan wtSpStMan("TiledWgtSpectrum", dts);
         _ms->addColumn(tdWtSp, wtSpStMan);
    }
    // default sort columns are from MSIter and are ARRAY_ID, FIELD_ID, DATA_DESC_ID, and TIME
    // I'm adding scan and state because, according to the statwt requirements, by default, scan
    // and state changes should mark boundaries in the weights computation
    std::vector<Int> scs;
    scs.push_back(MS::ARRAY_ID);
    if (! _combine.contains("scan")) {
        scs.push_back(MS::SCAN_NUMBER);
    }
    if (! _combine.contains("state")) {
        scs.push_back(MS::STATE_ID);
    }
    if (! _combine.contains("field")) {
        scs.push_back(MS::FIELD_ID);
    }
    scs.push_back(MS::DATA_DESC_ID);
    scs.push_back(MS::TIME);
    Block<int> sort(scs.size());
    uInt i = 0;
    for (const auto& col: scs) {
        sort[i] = col;
        ++i;
    }
    vi::SortColumns sc(sort, False);
    vi::IteratingParameters ipar(_timeBinWidth, sc);
    vi::VisIterImpl2LayerFactory data(_ms, ipar, True);
    unique_ptr<Record> config(dynamic_cast<Record*>(_tviConfig.clone()));
    vi::StatWtTVILayerFactory statWtLayerFactory(*config);
    Vector<vi::ViiLayerFactory*> facts(2);
    facts[0] = &data;
    facts[1] = &statWtLayerFactory;
    vi::VisibilityIterator2 vi(facts);
    vi::VisBuffer2 *vb = vi.getVisBuffer();
    Vector<Int> vr(1);
    for (vi.originChunks(); vi.moreChunks(); vi.nextChunk()) {
        for (vi.origin(); vi.more(); vi.next()) {
            if (_preview) {
                // just need to run the flags to accumulate
                // flagging info
                vb->flagCube();
            }
            else {
                if (needWtSp) {
                    Int nrow = vb->nRows();
                    Int nchan = vb->nChannels();
                    Int ncor = vb->nCorrelations();
                    Cube<Float> newwtsp(0, 0, 0);
                    newwtsp.resize(ncor, nchan, nrow);
                    newwtsp.set(0.0);
                    vb->initWeightSpectrum(newwtsp);
                    vb->writeChangesBack();
                }
                if (doWtSp) {
                    vb->setWeightSpectrum(vb->weightSpectrum());
                }
                vb->setWeight(vb->weight());
                vb->setFlagCube(vb->flagCube());
                vb->setFlagRow(vb->flagRow());
                vb->writeChangesBack();
            }
        }
    }
    if (_preview) {
        LogIO log(LogOrigin("StatWt", __func__));
        log << LogIO::NORMAL
            << "Running in preview mode, no data were changed"
            << LogIO::POST;
    }
    statWtLayerFactory.getTVI()->summarizeFlagging();
}

}

