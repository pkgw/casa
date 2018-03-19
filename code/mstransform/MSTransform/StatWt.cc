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
#include <casacore/casa/System/ProgressMeter.h>
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
    _timeBinWidth = vi::StatWtTVI::getTimeBinWidthInSec(binWidth);
}

void StatWt::setTimeBinWidth(Double binWidth) {
    vi::StatWtTVI::checkTimeBinWidth(binWidth);
    _timeBinWidth = binWidth;
}

void StatWt::setTimeBinWidthUsingInterval(uInt n) {
    _timeBinWidth = vi::StatWtTVI::getTimeBinWidthUsingInterval(_ms, n);
    _log << LogOrigin("StatWt", __func__) << LogIO::NORMAL
        << "Determined representative integration time of "
        << (_timeBinWidth/(Double)n) << "s. Setting time bin width to "
        << _timeBinWidth << "s" << LogIO::POST;
}

void StatWt::setCombine(const String& combine) {
    _combine = downcase(combine);
}

void StatWt::setPreview(casacore::Bool preview) {
    _preview = preview;
}

void StatWt::setTVIConfig(const Record& config) {
    static const String field = "datacolumn";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported type for field '" + field + "'"
        );
        auto val = config.asString(field);
        if (! val.empty()) {
            val.downcase();
            ThrowIf (
                ! (
                    val.startsWith("c") || val.startsWith("d")
                    || val.startsWith("residual") || val.startsWith("residual_")
                ),
                "Unsupported value for " + field + ": " + val
            );
            _possiblyWriteSigma = val.startsWith("d") || val.startsWith("residual_");
        }
    }
    _tviConfig = config;
}

Record StatWt::writeWeights() const {
    auto hasWtSp = _ms->isColumn(MSMainEnums::WEIGHT_SPECTRUM);
    auto mustWriteSigma = _possiblyWriteSigma && ! _preview;
    auto hasSigSp = _ms->isColumn(MSMainEnums::SIGMA_SPECTRUM);
    auto mustWriteSigmaSp = mustWriteSigma && hasSigSp;
    auto mustWriteWt = ! _preview
        && (
            ! mustWriteSigma
            || (
                mustWriteSigma && ! _ms->isColumn(MSMainEnums::CORRECTED_DATA)
            )
        );
    auto mustWriteWtSp = mustWriteWt
        && _tviConfig.isDefined(vi::StatWtTVI::CHANBIN);
    if (mustWriteWtSp) {
        auto type = _tviConfig.type(_tviConfig.fieldNumber(vi::StatWtTVI::CHANBIN));
        if (type == TpArrayBool) {
            // default variant type
            mustWriteWtSp = False;
        }
        else if (type == TpString) {
            auto val = _tviConfig.asString(vi::StatWtTVI::CHANBIN);
            val.downcase();
            if (val == "spw") {
                mustWriteWtSp = False;
            }
        }
    }
    auto mustInitWtSp = False;
    auto mustInitSigSp = False;
    // this conditional structure supports the
    // case of ! hasWtSp && ! mustWriteWtSp, in which case,
    // nothing need be done
    if (! hasWtSp) {
        if (mustWriteWtSp) {
            // we must create WEIGHT_SPECTRUM
            hasWtSp = True;
            mustInitWtSp = True;
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
    }
    else if (mustWriteWt) {
        // check to see if extant WEIGHT_SPECTRUM needs to be initialized
        ArrayColumn<Float> col(*_ms, MS::columnName(MS::WEIGHT_SPECTRUM));
        try {
            col.get(0);
            // its initialized, so even if we are using the full spw for
            // binning, we still need to update WEIGHT_SPECTRUM
            mustWriteWtSp = True;
        }
        catch (const AipsError& x) {
            // its not initialized, so we aren't going to write to it unless
            // chanbin has been specified to be less than the spw width
            mustInitWtSp = mustWriteWtSp;
        }
    }
    // *******
    if (! hasSigSp) {
        if (mustWriteSigmaSp) {
            // we must create SIGMA_SPECTRUM
            hasSigSp = True;
            mustInitSigSp = True;
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
            String colSigSp = MS::columnName(MS::SIGMA_SPECTRUM);
            TableDesc tdSigSp;
            tdSigSp.addColumn(ArrayColumnDesc<Float>(colSigSp, "sigma spectrum", 2));
            TiledShapeStMan sigSpStMan("TiledWgtSpectrum", dts);
            _ms->addColumn(tdSigSp, sigSpStMan);
            cout << "added column sigma spectrum" << endl;
        }
    }
    else if (mustWriteSigma) {
        // check to see if extant SIGMA_SPECTRUM needs to be initialized
        ArrayColumn<Float> col(*_ms, MS::columnName(MS::SIGMA_SPECTRUM));
        try {
            col.get(0);
            // its initialized, so even if we are using the full spw for
            // binning, we still need to update SIGMA_SPECTRUM
            mustWriteSigmaSp = True;
            cout << "sigma spec already initialized" << endl;
        }
        catch (const AipsError& x) {
            // its not initialized, so we aren't going to write to it unless
            // chanbin has been specified to be less than the spw width
            mustInitSigSp = mustWriteSigmaSp;
            cout << "sigma spec not initialized, we must do that below" << endl;
        }
    }
    LogIO log(LogOrigin("StatWt", __func__));
    if (mustWriteWt) {
        if (mustWriteSigma) {
            log << LogIO::NORMAL
                << "CORRECTED_DATA is not present. Updating the "
                << "SIGMA/SIGMA_SPECTRUM and WEIGHT/WEIGHT_SPECTRUM values "
                << "based on calculations using the DATA column."
                << LogIO::POST;
        }
        else {
            log << LogIO::NORMAL
                << "Updating the WEIGHT/WEIGHT_SPECTRUM values. SIGMA/SIGMA_SPECTRUM "
                << "values will not be recalculated as they are related to the values "
                << "in the DATA column." << LogIO::POST;
            }
    }
    else if (mustWriteSigma) {
        log << LogIO::NORMAL
        << "Updating the SIGMA/SIGMA_SPECTRUM values. WEIGHT/WEIGHT_SPECTRUM will "
        << "not be recalculated as they are related to the values in the "
        << "CORRECTED_DATA column." << LogIO::POST;
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
    ProgressMeter pm(0, _ms->nrow(), "StatWt Progress");
    uInt64 count = 0;
    for (vi.originChunks(); vi.moreChunks(); vi.nextChunk()) {
        for (vi.origin(); vi.more(); vi.next()) {
            auto nrow = vb->nRows();
            if (_preview) {
                // just need to run the flags to accumulate
                // flagging info
                vb->flagCube();
            }
            else {
                if (mustInitWtSp) {
                    auto nchan = vb->nChannels();
                    auto ncor = vb->nCorrelations();
                    Cube<Float> newwtsp(ncor, nchan, nrow, 0);
                    vb->initWeightSpectrum(newwtsp);
                    vb->writeChangesBack();
                }
                cout << "mustinitsigsp " << mustInitSigSp << endl;
                if (mustInitSigSp) {
                    auto nchan = vb->nChannels();
                    auto ncor = vb->nCorrelations();
                    Cube<Float> newsigsp(ncor, nchan, nrow, 0);
                    cout << "init sig sp" << endl;
                    vb->initSigmaSpectrum(newsigsp);
                    cout << "sig sp inited, write changes" << endl;
                    vb->writeChangesBack();
                    cout << "changes written" << endl;
                }
                if (mustWriteWtSp) {
                    vb->setWeightSpectrum(vb->weightSpectrum());
                }
                if (mustWriteSigmaSp) {
                    cout << "writing new sigma sp" << endl;
                    vb->setSigmaSpectrum(vb->sigmaSpectrum());
                }
                if (mustWriteWt) {
                    vb->setWeight(vb->weight());
                }
                if (mustWriteSigma) {
                    vb->setSigma(vb->sigma());
                }
                vb->setFlagCube(vb->flagCube());
                vb->setFlagRow(vb->flagRow());
                cout << "writing back all changes" << endl;
                vb->writeChangesBack();
                cout << "changes written" << endl;
            }
            count += nrow;
            pm.update(count);
        }
    }
    if (_preview) {
        LogIO log(LogOrigin("StatWt", __func__));
        log << LogIO::NORMAL
            << "RAN IN PREVIEW MODE. NO WEIGHTS NOR FLAGS WERE CHANGED."
            << LogIO::POST;
    }
    statWtLayerFactory.getTVI()->summarizeFlagging();
    Double mean, variance;
    statWtLayerFactory.getTVI()->summarizeStats(mean, variance);
    Record ret;
    ret.define("mean", mean);
    ret.define("variance", variance);
    return ret;
}

}

