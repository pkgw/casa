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

Record StatWt::writeWeights() {
    auto mustWriteWt = False;
    auto mustWriteWtSp = False;
    auto mustInitWtSp = False;
    auto mustWriteSig = False;
    auto mustWriteSigSp = False;
    auto mustInitSigSp = False;
    _columnInitWrite(
        mustWriteWt, mustWriteWtSp, mustInitWtSp,
        mustWriteSig, mustWriteSigSp, mustInitSigSp
    );
    shared_ptr<vi::VisibilityIterator2> vi;
    std::shared_ptr<vi::StatWtTVILayerFactory> factory;
    _constructVi(vi, factory);
    vi::VisBuffer2 *vb = vi->getVisBuffer();
    ProgressMeter pm(0, _ms->nrow(), "StatWt Progress");
    uInt64 count = 0;
    for (vi->originChunks(); vi->moreChunks(); vi->nextChunk()) {
        for (vi->origin(); vi->more(); vi->next()) {
            auto nrow = vb->nRows();
            if (_preview) {
                // just need to run the flags to accumulate
                // flagging info
                vb->flagCube();
            }
            else {
                if (mustInitWtSp || mustInitSigSp) {
                    auto nchan = vb->nChannels();
                    auto ncor = vb->nCorrelations();
                    Cube<Float> newsp(ncor, nchan, nrow, 0);
                    if (mustInitWtSp) {
                        vb->initWeightSpectrum(newsp);
                    }
                    if (mustInitSigSp) {
                        vb->initSigmaSpectrum(newsp);
                    }
                    vb->writeChangesBack();
                }
                if (mustWriteWtSp) {
                    vb->setWeightSpectrum(vb->weightSpectrum());
                }
                if (mustWriteSigSp) {
                    vb->setSigmaSpectrum(vb->sigmaSpectrum());
                }
                if (mustWriteWt) {
                    vb->setWeight(vb->weight());
                }
                if (mustWriteSig) {
                    vb->setSigma(vb->sigma());
                }
                vb->setFlagCube(vb->flagCube());
                vb->setFlagRow(vb->flagRow());
                vb->writeChangesBack();
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
    factory->getTVI()->summarizeFlagging();
    Double mean, variance;
    factory->getTVI()->summarizeStats(mean, variance);
    Record ret;
    ret.define("mean", mean);
    ret.define("variance", variance);
    return ret;
}

void StatWt::_columnInitWrite(
    casacore::Bool& mustWriteWt, casacore::Bool& mustWriteWtSp,
    casacore::Bool& mustInitWtSp, casacore::Bool& mustWriteSig,
    casacore::Bool& mustWriteSigSp, casacore::Bool& mustInitSigSp
) {
    auto hasWtSp = _ms->isColumn(MSMainEnums::WEIGHT_SPECTRUM);
    mustWriteSig = _possiblyWriteSigma && ! _preview;
    auto hasSigSp = _ms->isColumn(MSMainEnums::SIGMA_SPECTRUM);
    mustWriteSigSp = mustWriteSig && hasSigSp;
    mustWriteWt = ! _preview
        && (
            ! mustWriteSig
            || (
                mustWriteSig && ! _ms->isColumn(MSMainEnums::CORRECTED_DATA)
            )
        );
    mustWriteWtSp = mustWriteWt
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
    mustInitWtSp = False;
    mustInitSigSp = False;
    static const auto colNameWtSp = MS::columnName(MS::WEIGHT_SPECTRUM);
    static const auto descWtSp = "weight spectrum";
    _dealWithSpectrumColumn(
        hasWtSp, mustWriteWtSp, mustInitWtSp,
        mustWriteWt, colNameWtSp, descWtSp
    );
    static const auto colNameSigSp = MS::columnName(MS::SIGMA_SPECTRUM);
    static const auto descSigSp = "sigma spectrum";
    _dealWithSpectrumColumn(
        hasSigSp, mustWriteSigSp, mustInitSigSp,
        mustWriteSig, colNameSigSp, descSigSp
    );
    LogIO log(LogOrigin("StatWt", __func__));
    if (mustWriteWt) {
        if (mustWriteSig) {
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
    else if (mustWriteSig) {
        log << LogIO::NORMAL
            << "Updating the SIGMA/SIGMA_SPECTRUM values. WEIGHT/WEIGHT_SPECTRUM will "
            << "not be recalculated as they are related to the values in the "
            << "CORRECTED_DATA column." << LogIO::POST;
    }
}

void StatWt::_constructVi(
    std::shared_ptr<vi::VisibilityIterator2>& vi,
    std::shared_ptr<vi::StatWtTVILayerFactory>& factory
) const {
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
    factory.reset(new vi::StatWtTVILayerFactory(*config));
    Vector<vi::ViiLayerFactory*> facts(2);
    facts[0] = &data;
    facts[1] = factory.get();
    vi.reset(new vi::VisibilityIterator2(facts));
}

void StatWt::_dealWithSpectrumColumn(
    Bool& hasSpec, Bool& mustWriteSpec, Bool& mustInitSpec,
    Bool mustWriteNonSpec, const String& colName, const String& descName
) {
    // this conditional structure supports the
    // case of ! hasSpec && ! mustWriteSpec, in which case,
    // nothing need be done
    if (! hasSpec) {
        if (mustWriteSpec) {
            // we must create spectrum column
            hasSpec = True;
            mustInitSpec = True;
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
            String colWtSp = colName;
            TableDesc tdWtSp;
            tdWtSp.addColumn(ArrayColumnDesc<Float>(colWtSp, descName, 2));
            TiledShapeStMan wtSpStMan("TiledWgtSpectrum", dts);
            _ms->addColumn(tdWtSp, wtSpStMan);
        }
    }
    else if (mustWriteNonSpec) {
        // check to see if extant spectrum column needs to be initialized
        ArrayColumn<Float> col(*_ms, colName);
        try {
            col.get(0);
            // its initialized, so even if we are using the full spw for
            // binning, we still need to update WEIGHT_SPECTRUM
            mustWriteSpec = True;
        }
        catch (const AipsError& x) {
            // its not initialized, so we aren't going to write to it unless
            // chanbin has been specified to be less than the spw width
            mustInitSpec = mustWriteSpec;
        }
    }
}

}

