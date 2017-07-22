//# StatWtTVI.cc: This file contains the implementation of the StatWtTVI class.
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

#include <mstransform/TVI/StatWtTVI.h>

#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/ms/MSOper/MSMetaData.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#include <iomanip>

using namespace casacore;
using namespace casac;

namespace casa { 
namespace vi { 

const String StatWtTVI::CHANBIN = "stchanbin";

StatWtTVI::StatWtTVI(ViImplementation2 * inputVii, const Record &configuration)
    : TransformingVi2 (inputVii) {
	// Parse and check configuration parameters
	// Note: if a constructor finishes by throwing an exception, the memory
	// associated with the object itself is cleaned up there is no memory leak.
    ThrowIf(
        ! _parseConfiguration(configuration),
	    "Error parsing StatWtTVI configuration"
    );
    ThrowIf(
        _useCorrected && ! ms().isColumn(MSMainEnums::CORRECTED_DATA),
        "StatWtTVI requires the MS to have a "
        "CORRECTED_DATA column. This MS does not"
    );
    ThrowIf(
        ! _useCorrected && ! ms().isColumn(MSMainEnums::DATA),
        "StatWtTVI requires the MS to have a "
        "DATA column. This MS does not"
    );
	_initialize();
	// Initialize attached VisBuffer
	setVisBuffer(createAttachedVisBuffer(VbRekeyable));
}

StatWtTVI::~StatWtTVI() {}

Bool StatWtTVI::_parseConfiguration(const Record& config) {
    String field = CHANBIN;
    if (config.isDefined(field)) {
        // channel binning
        auto fieldNum = config.fieldNumber(field);
        switch (config.type(fieldNum)) {
        case DataType::TpArrayBool:
            // because this is the actual default variant type, no matter
            // what is specified in the xml
            ThrowIf(
                ! config.asArrayBool(field).empty(),
                "Unsupported data type for " + field
            );
            _setDefaultChanBinMap();
            break;
        case DataType::TpInt:
            Int binWidth;
            config.get(CHANBIN, binWidth);
            _setChanBinMap(binWidth);
            break;
        case DataType::TpString:
        {
            auto chanbin = config.asString(field);
            if (chanbin == "spw") {
                // bin using entire spws
                _setDefaultChanBinMap();
                break;
            }
            else {
                QuantumHolder qh(casaQuantity(chanbin));
                _setChanBinMap(qh.asQuantity());
            }
            break;
        }
        default:
            ThrowCc("Unsupported data type for " + field);
        }
    }
    else {
        _setDefaultChanBinMap();
    }
    field = "minsamp";
    if (config.isDefined(field)) {
        config.get(field, _minSamp);
        ThrowIf(_minSamp < 2, "Minimum size of sample must be >= 2.");
    }
    field = "combine";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported data type for combine"
        );
        _combineCorr = config.asString(field).contains("corr");
    }
    field = "wtrange";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpArrayDouble,
            "Unsupported type for field '" + field + "'"
        );
        auto myrange = config.asArrayDouble(field);
        if (! myrange.empty()) {
            ThrowIf(
                myrange.size() != 2,
                "Array specified in '" + field + "' must have exactly two values"
            );
            ThrowIf(
                casacore::anyLT(myrange, 0.0),
                "Both values specified in '" + field + "' array must be non-negative"
            );
            std::set<Double> rangeset(myrange.begin(), myrange.end());
            ThrowIf(
                rangeset.size() == 1, "Values specified in '" + field + "' array must be unique"
            );
            auto iter = rangeset.begin();
            _wtrange.reset(new std::pair<Double, Double>(*iter, *(++iter)));
        }
    }
    field = "excludechans";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported type for field '" + field + "'"
        );
        auto val = config.asString(field);
        if (! val.empty()) {
            MSSelection sel(ms());
            sel.setSpwExpr(val);
            auto chans = sel.getChanList();
            auto nrows = chans.nrow();
            //const auto& myms = ms();
            MSMetaData md(&ms(), 50);
            auto nchans = md.nChans();
            IPosition start(3, 0);
            IPosition stop(3, 0);
            IPosition step(3, 1);
            for (uInt i=0; i<nrows; ++i) {
                auto row = chans.row(i);
                const auto& spw = row[0];
                if (_chanSelFlags.find(spw) == _chanSelFlags.end()) {
                    _chanSelFlags[spw] = Cube<Bool>(1, nchans[spw], 1, False);
                }
                start[1] = row[1];
                stop[1] = row[2];
                step[1] = row[3];
                Slicer slice(start, stop, step, Slicer::endIsLast);
                _chanSelFlags[spw](slice) = True;
            }
        }
    }
    field = "datacolumn";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported type for field '" + field + "'"
        );
        auto val = config.asString(field);
        if (! val.empty()) {
            val.downcase();
            ThrowIf (
                ! (val.startsWith("c") || val.startsWith("d")),
                "Unsupported value for " + field + ": " + val
            );
            _useCorrected = val.startsWith("c");
        }
    }
    _configureStatAlg(config);
    return True;
}

void StatWtTVI::_configureStatAlg(const Record& config) {
    String field = "statalg";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported type for field '" + field + "'"
        );
        auto alg = config.asString(field);
        alg.downcase();
        if (alg.startsWith("cl")) {
            _statAlg = new ClassicalStatistics<
                Double, Array<Float>::const_iterator,
                Array<Bool>::const_iterator
            >();
        }
        else {
            StatisticsAlgorithmFactory<
                Double, Array<Float>::const_iterator,
                Array<Bool>::const_iterator
            > saf;
            if (alg.startsWith("ch")) {
                Int maxiter = -1;
                field = "maxiter";
                if (config.isDefined(field)) {
                    ThrowIf(
                        config.type(config.fieldNumber(field)) != TpInt,
                        "Unsupported type for field '" + field + "'"
                    );
                    maxiter = config.asInt(field);
                }
                Double zscore = -1;
                field = "zscore";
                if (config.isDefined(field)) {
                    ThrowIf(
                        config.type(config.fieldNumber(field)) != TpDouble,
                        "Unsupported type for field '" + field + "'"
                    );
                    zscore = config.asDouble(field);
                }
                saf.configureChauvenet(zscore, maxiter);
            }
            else if (alg.startsWith("f")) {
                FitToHalfStatisticsData::CENTER center = FitToHalfStatisticsData::CMEAN;
                field = "center";
                if (config.isDefined(field)) {
                    ThrowIf(
                        config.type(config.fieldNumber(field)) != TpString,
                        "Unsupported type for field '" + field + "'"
                    );
                    auto cs = config.asString(field);
                    cs.downcase();
                    if (cs == "mean") {
                        center = FitToHalfStatisticsData::CMEAN;
                    }
                    else if (cs == "median") {
                        center = FitToHalfStatisticsData::CMEDIAN;
                    }
                    else if (cs == "zero") {
                        center = FitToHalfStatisticsData::CVALUE;
                    }
                    else {
                        ThrowCc("Unsupported value for '" + field + "'");
                    }
                }
                field = "lside";
                FitToHalfStatisticsData::USE_DATA ud = FitToHalfStatisticsData::LE_CENTER;
                if (config.isDefined(field)) {
                    ThrowIf(
                        config.type(config.fieldNumber(field)) != TpBool,
                        "Unsupported type for field '" + field + "'"
                    );
                    ud = config.asBool(field)
                        ? FitToHalfStatisticsData::LE_CENTER
                        : FitToHalfStatisticsData::GE_CENTER;
                }
                saf.configureFitToHalf(center, ud, 0);
            }
            else if (alg.startsWith("h")) {
                Double fence = -1;
                field = "fence";
                if (config.isDefined(field)) {
                    ThrowIf(
                        config.type(config.fieldNumber(field)) != TpDouble,
                        "Unsupported type for field '" + field + "'"
                    );
                    fence = config.asDouble(field);
                }
                saf.configureHingesFences(fence);
            }
            else {
                ThrowCc("Unsupported value for 'statalg'");
            }
            _statAlg = saf.createStatsAlgorithm();
        }
    }
    else {
        _statAlg = new ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator
        >();
    }
    std::set<StatisticsData::STATS> stats;
    stats.insert(StatisticsData::VARIANCE);
    _statAlg->setStatsToCalculate(stats);
}

void StatWtTVI::_setChanBinMap(const casacore::Quantity& binWidth) {
    if (! binWidth.isConform(Unit("Hz"))) {
        ostringstream oss;
        oss << "If specified as a quantity, channel bin width must have frequency units. "
            << binWidth << " does not.";
        ThrowCc(oss.str());
    }
    ThrowIf(
        binWidth.getValue() <= 0,
        "channel bin width must be positive"
    );
    MSMetaData msmd(&ms(), 100.0);
    auto chanFreqs = msmd.getChanFreqs();
    auto nspw = chanFreqs.size();
    auto binWidthHz = binWidth.getValue("Hz");
    for (uInt i=0; i<nspw; ++i) {
        auto cfs = chanFreqs[i].getValue("Hz");
        auto citer = cfs.begin();
        auto cend = cfs.end();
        ChanBin bin;
        bin.start = 0;
        bin.end = 0;
        uInt chanNum = 0;
        auto startFreq = *citer;
        auto nchan = cfs.size();
        for (; citer!=cend; ++citer, ++chanNum) {
            if (abs(*citer - startFreq) > binWidthHz) {
                // start new bin
                _chanBins[i].push_back(bin);
                bin.start = chanNum;
                startFreq = *citer;
            }
            bin.end = chanNum;
            if (chanNum + 1 == nchan) {
                // need to add the last bin
                _chanBins[i].push_back(bin);
            }
        }
    }
    // weight spectrum must be written
    _mustComputeWtSp.reset(new Bool(True));
}

void StatWtTVI::_setChanBinMap(Int binWidth) {
    ThrowIf(binWidth < 2, "Channel bin width must >= 2");
    MSMetaData msmd(&ms(), 100.0);
    auto nchans = msmd.nChans();
    auto nspw = nchans.size();
    ChanBin bin;
    for (uInt i=0; i<nspw; ++i) {
        auto lastChan = nchans[i]-1;
        for (uInt j=0; j<nchans[i]; j += binWidth) {
            bin.start = j;
            bin.end = min(j+binWidth-1, lastChan);
            _chanBins[i].push_back(bin);
        }
    }
    // weight spectrum must be written
    _mustComputeWtSp.reset(new Bool(True));
}

void StatWtTVI::_setDefaultChanBinMap() {
    MSMetaData msmd(&ms(), 100.0);
    auto nchans = msmd.nChans();
    auto niter = nchans.begin();
    auto nend = nchans.end();
    Int i = 0;
    ChanBin bin;
    bin.start = 0;
    for (; niter!=nend; ++niter, ++i) {
        bin.end = *niter - 1;
        _chanBins[i].push_back(bin);
    }
}

void StatWtTVI::_initialize() {}

void StatWtTVI::weightSpectrum(Cube<Float> & newWtsp) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! *_mustComputeWtSp) {
        newWtsp.resize(IPosition(3, 0));
        return;
    }
    if (! _newWtSp.empty()) {
        // already calculated
        newWtsp = _newWtSp.copy();
        return;
    }
    getVii()->weightSpectrum(newWtsp);
    Vector<Int> ant1, ant2, spws;
    antenna1(ant1);
    antenna2(ant2);
    spectralWindows(spws);
    IPosition blc(3, 0);
    auto trc = newWtsp.shape() - 1;
    auto nrows = nRows();
    for (Int i=0; i<nrows; ++i) {
        blc[2] = i;
        trc[2] = i;
        BaselineChanBin blcb;
        blcb.baseline = _baseline(ant1[i], ant2[i]);
        auto spw = spws[i];
        blcb.spw = spw;
        auto bins = _chanBins.find(spw)->second;
        auto biter = bins.begin();
        auto bend = bins.end();
        for (; biter!=bend; ++biter) {
            blc[1] = biter->start;
            trc[1] = biter->end;
            blcb.chanBin = *biter;
            auto weights = _weights.find(blcb)->second;
            auto ncorr = weights.size();
            for (uInt corr=0; corr<ncorr; ++corr) {
                blc[0] = _combineCorr ? 0 : corr;
                trc[0] = _combineCorr ? newWtsp.shape()[0] - 1 : corr;
                newWtsp(blc, trc) = weights[corr];
            }
        }
    }
    // cache it
    _newWtSp = newWtsp.copy();
}

void StatWtTVI::weight(Matrix<Float> & wtmat) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! _newWt.empty()) {
        wtmat = _newWt.copy();
        return;
    }
    auto nrows = nRows();
    getVii()->weight(wtmat);
    if (*_mustComputeWtSp) {
        // always use classical algorithm to get median for weights
        ClassicalStatistics<Double, Array<Float>::const_iterator, Array<Bool>::const_iterator> cs;
        Cube<Float> newWtsp;
        Cube<Bool> flagCube;
        weightSpectrum(newWtsp);
        flag(flagCube);
        IPosition blc(3, 0);
        IPosition trc = newWtsp.shape() - 1;
        const auto ncorr = newWtsp.shape()[0];
        for (Int i=0; i<nrows; ++i) {
            blc[2] = i;
            trc[2] = i;
            if (_combineCorr) {
                auto flags = flagCube(blc, trc);
                if (allTrue(flags)) {
                    wtmat.column(i) = 0;
                }
                else {
                    auto weights = newWtsp(blc, trc);
                    auto mask = ! flags;
                    cs.setData(weights.begin(), mask.begin(), weights.size());
                    wtmat.column(i) = cs.getMedian();
                }
            }
            else {
                for (uInt corr=0; corr<ncorr; ++corr) {
                    blc[0] = corr;
                    trc[0] = corr;
                    auto weights = newWtsp(blc, trc);
                    auto flags = flagCube(blc, trc);
                    if (allTrue(flags)) {
                        wtmat(corr, i) = 0;
                    }
                    else {
                        auto mask = ! flags;
                        cs.setData(weights.begin(), mask.begin(), weights.size());
                        wtmat(corr, i) = cs.getMedian();
                    }
                }
            }
        }
    }
    else {
        // the only way this can happen is if there is a single channel bin
        // for each baseline/spw pair
        Vector<Int> ant1, ant2, spws;
        antenna1(ant1);
        antenna2(ant2);
        spectralWindows(spws);
        BaselineChanBin blcb;
        for (Int i=0; i<nrows; ++i) {
            auto bins = _chanBins.find(spws[i])->second;
            blcb.baseline = _baseline(ant1[i], ant2[i]);
            blcb.spw = spws[1];
            blcb.chanBin = bins[0];
            auto weights = _weights.find(blcb)->second;
            if (_combineCorr) {
                wtmat.column(i) = weights[0];
            }
            else {
                auto corr = 0;
                for (const auto weight: weights) {
                    wtmat(corr, i) = weight;
                    ++corr;
                }
            }
        }
    }
    _newWt = wtmat.copy();
}

void StatWtTVI::flag(Cube<Bool>& flagCube) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! _newFlag.empty()) {
        flagCube = _newFlag.copy();
        return;
    }
    getVii()->flag(flagCube);
    _nTotalPts += flagCube.size();
    auto nOrigFlagged = ntrue(flagCube);
    _nOrigFlaggedPts += nOrigFlagged;
    Vector<Int> ant1, ant2, spws;
    antenna1(ant1);
    antenna2(ant2);
    spectralWindows(spws);
    auto nrows = nRows();
    IPosition blc(3, 0);
    auto trc = flagCube.shape() - 1;
    auto ncorr = _combineCorr ? 1 : flagCube.shape()[0];
    BaselineChanBin blcb;
    auto checkFlags = False;
    for (Int i=0; i<nrows; ++i) {
        blcb.baseline = _baseline(ant1[i], ant2[i]);
        auto spw = spws[i];
        blcb.spw = spw;
        auto bins = _chanBins.find(spw)->second;
        auto biter = bins.begin();
        auto bend = bins.end();
        blc[2] = i;
        trc[2] = i;
        for (; biter!=bend; ++biter) {
            blc[1] = biter->start;
            trc[1] = biter->end;
            blcb.chanBin = *biter;
            auto weights = _weights.find(blcb)->second;
            for (uInt corr=0; corr<ncorr; ++corr) {
                blc[0] = _combineCorr ? 0 : corr;
                trc[0] = _combineCorr ? flagCube.shape()[0] - 1 : corr;
                auto wt = weights[corr];
                if (
                    wt == 0
                    || (_wtrange && (wt < _wtrange->first || wt > _wtrange->second))
                ) {
                    checkFlags = True;
                    flagCube(blc, trc) = True;
                }
            }
        }
    }
    if (checkFlags) {
        _nNewFlaggedPts += ntrue(flagCube) - nOrigFlagged;
    }
    _newFlag = flagCube.copy();
}

void StatWtTVI::flagRow (Vector<Bool>& flagRow) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! _newFlagRow.empty()) {
        flagRow = _newFlagRow.copy();
        return;
    }
    Cube<Bool> flags;
    flag(flags);
    getVii()->flagRow(flagRow);
    auto nrows = nRows();
    for (Int i=0; i<nrows; ++i) {
        flagRow[i] = allTrue(flags.xyPlane(i));
    }
    _newFlagRow = flagRow.copy();
}

void StatWtTVI::originChunks(Bool forceRewind) {
    // Drive next lower layer
    getVii()->originChunks(forceRewind);
    _weightsComputed = False;
    _gatherAndComputeWeights();
    _weightsComputed = True;
    _clearCache();
    // re-origin this chunk in next layer
    //  (ensures wider scopes see start of the this chunk)
    getVii()->origin();
}

void StatWtTVI::nextChunk() {
    // Drive next lower layer
    getVii()->nextChunk();
    _weightsComputed = False;
    _gatherAndComputeWeights();
    _weightsComputed = True;
    _clearCache();
    // re-origin this chunk next layer
    //  (ensures wider scopes see start of the this chunk)
    getVii()->origin();
}

void StatWtTVI::_clearCache() {
    _newWtSp.resize(0, 0, 0);
    _newWt.resize(0, 0);
    _newFlag.resize(0, 0, 0);
    _newFlagRow.resize(0);
}

void StatWtTVI::_gatherAndComputeWeights() const {
    // Drive NEXT LOWER layer's ViImpl to gather data into allvis:
    //  Assumes all sub-chunks in the current chunk are to be used
    //   for the variance calculation
    //  Essentially, we are sorting the incoming data into
    //   allvis, to enable a convenient variance calculation
    _weights.clear();
    ViImplementation2* vii = getVii();
    VisBuffer2* vb = vii->getVisBuffer();
    _newRowIDs.resize(vii->nRowsInChunk());
    // baseline to visibility, flag maps
    std::map<BaselineChanBin, Cube<Complex>> data;
    std::map<BaselineChanBin, Cube<Bool>> flags;
    IPosition blc(3, 0);
    auto trc = blc;
    auto doChanSelFlags = ! _chanSelFlags.empty();
    auto initChanSelFlags = doChanSelFlags;
    Cube<Bool> chanSelFlagCube;
    Cube<Bool> myChanSelFlags;
    IPosition mystart(3, 0);
    IPosition mystop(3, 0);
    Slicer sl(mystart, mystop, Slicer::endIsLast);
    auto chunkChecked = False;
    for (vii->origin(); vii->more(); vii->next()) {
        const auto& rowIDs = vb->rowIds();
        if (! chunkChecked) {
            if (_processedRowIDs.find(rowIDs[0]) == _processedRowIDs.end()) {
                // haven't processed this chunk
                _processedRowIDs.insert(rowIDs[0]);
                chunkChecked = True;
            }
            else {
                // this chunk has been processed, this can happen at the end
                // when the last chunk is processed twice
                return;
            }
        }
        if (! _mustComputeWtSp) {
            _mustComputeWtSp.reset(new Bool(vb->existsColumn(VisBufferComponent2::WeightSpectrum)));
        }
        const auto& ant1 = vb->antenna1();
        const auto& ant2 = vb->antenna2();
        // [nC,nF,nR)
        const auto& dataCube = _useCorrected
            ? vb->visCubeCorrected() : vb->visCube();
        IPosition dataCubeBLC(3, 0);
        auto dataCubeTRC = dataCube.shape() - 1;
        dataCubeTRC[2] = 0;
        const auto flagCube = vb->flagCube();
        const auto nrows = vb->nRows();
        const auto npol = dataCube.nrow();
        const auto spws = vb->spectralWindows();
        if (initChanSelFlags) {
            // this can be done just once because all the rows
            // in the chunk are guaranteed to have the same spw
            // because each subchunk is guaranteed to have a single
            // data description ID.
            auto spw = *spws.begin();
            auto chanSelFlagIter = _chanSelFlags.find(spw);
            doChanSelFlags = chanSelFlagIter != _chanSelFlags.end();
            if (doChanSelFlags) {
                chanSelFlagCube = chanSelFlagIter->second;
            }
            initChanSelFlags = False;
        }
        if (doChanSelFlags) {
            auto dataShape = dataCube.shape();
            myChanSelFlags.resize(dataShape, False);
            auto nchan = dataShape[1];
            auto ncorr = dataShape[0];
            mystop[1] = nchan-1;
            for (uInt corr=0; corr<ncorr; ++corr) {
                mystart[0] = corr;
                mystop[0] = corr;
                for (Int row=0; row<nrows; ++row) {
                    mystart[2] = row;
                    mystop[2] = row;
                    sl.setStart(mystart);
                    sl.setEnd(mystop);
                    myChanSelFlags(sl) = chanSelFlagCube;
                }
            }
        }
        for (Int row=0; row<nrows; ++row) {
            dataCubeBLC[2] = row;
            dataCubeTRC[2] = row;
            BaselineChanBin blcb;
            blcb.baseline = _baseline(ant1[row], ant2[row]);
            auto spw = spws[row];
            auto bins = _chanBins.find(spw)->second;
            blcb.spw = spw;
            auto citer = bins.begin();
            auto cend = bins.end();
            for (; citer!=cend; ++citer) {
                dataCubeBLC[1] = citer->start;
                dataCubeTRC[1] = citer->end;
                blcb.chanBin.start = citer->start;
                blcb.chanBin.end = citer->end;
                auto dataSlice = dataCube(dataCubeBLC, dataCubeTRC);
                auto flagSlice = doChanSelFlags
                    ? flagCube(dataCubeBLC, dataCubeTRC)
                        || myChanSelFlags(dataCubeBLC, dataCubeTRC)
                    : flagCube(dataCubeBLC, dataCubeTRC);
                if (data.find(blcb) == data.end()) {
                    data[blcb] = dataSlice;
                    flags[blcb] = flagSlice;
                }
                else {
                    auto myshape = data[blcb].shape();
                    auto nplane = myshape[2];
                    auto nchan = myshape[1];
                    data[blcb].resize(npol, nchan, nplane+1, True);
                    flags[blcb].resize(npol, nchan, nplane+1, True);
                    trc = myshape - 1;
                    // because we've extended the cube by one plane since
                    // myshape was determined.
                    ++trc[2];
                    blc[2] = trc[2];
                    data[blcb](blc, trc) = dataSlice;
                    flags[blcb](blc, trc) = flagSlice;
                }
            }
        }
    }
    // data has been gathered, now compute weights
    _computeWeights(data, flags);
}

void StatWtTVI::initWeightSpectrum (const casacore::Cube<casacore::Float>& wtspec) {
    // Pass to next layer down
    getVii()->initWeightSpectrum(wtspec);
}


void StatWtTVI::writeBackChanges(VisBuffer2 *vb) {
    // Pass to next layer down
    getVii()->writeBackChanges(vb);
}

StatWtTVI::Baseline StatWtTVI::_baseline(uInt ant1, uInt ant2) {
    Baseline baseline;
    if (ant1 < ant2) {
        // this may always be the case, but I'm not certain,
        baseline.first = ant1;
        baseline.second = ant2;
    }
    else {
        baseline.first = ant2;
        baseline.second = ant1;
    }
    return baseline;
}

void StatWtTVI::_computeWeights(
    const map<BaselineChanBin, Cube<Complex>>& data,
    const map<BaselineChanBin, Cube<Bool>>& flags
) const {
    auto diter = data.begin();
    auto dend = data.end();
    auto fiter = flags.begin();
    const auto nActCorr = diter->second.shape()[0];
    const auto ncorr = _combineCorr ? 1 : nActCorr;
    for (; diter!=dend; ++diter, ++fiter) {
        auto blcb = diter->first;
        auto spw = blcb.spw;
        if (_samples.find(spw) == _samples.end()) {
            _samples[spw].first = 0;
            _samples[spw].second = 0;
        }
        auto dataForBLCB = diter->second;
        auto flagsForBLCB = fiter->second;
        for (uInt corr=0; corr<ncorr; ++corr) {
            IPosition start(3, 0);
            IPosition end = dataForBLCB.shape() - 1;
            if (! _combineCorr) {
                start[0] = corr;
                end[0] = corr;
            }
            Slicer slice(start, end, Slicer::endIsLast);
            auto dataChunk = dataForBLCB(slice);
            const auto npts = dataChunk.size();
            if ((Int)npts < _minSamp) {
                // not enough points, trivial
                _weights[blcb].push_back(0);
            }
            else {
                auto flagChunk = flagsForBLCB(slice);
                if ((Int)nfalse(flagChunk) < _minSamp) {
                    // not enough points, trivial
                    _weights[blcb].push_back(0);
                }
                else {
                    // some data not flagged
                    const auto realPart = real(dataChunk);
                    const auto imagPart = imag(dataChunk);
                    const auto mask = ! flagChunk;
                    const auto riter = realPart.begin();
                    const auto iiter = imagPart.begin();
                    const auto miter = mask.begin();
                    _statAlg->setData(riter, miter, npts);
                    auto realVar = _statAlg->getStatistic(StatisticsData::VARIANCE);
                    _statAlg->setData(iiter, miter, npts);
                    auto imagVar = _statAlg->getStatistic(StatisticsData::VARIANCE);
                    auto varSum = realVar + imagVar;
                    _weights[blcb].push_back(varSum == 0 ? 0 : 2/varSum);
                    if (varSum > 0) {
                        ++_samples[spw].first;
                        if (imagVar == 0 || realVar == 0) {
                            ++_samples[spw].second;
                        }
                        else {
                            auto ratio = imagVar/realVar;
                            auto inverse = 1/ratio;
                            if (ratio > 1.5 || inverse > 1.5) {
                                ++_samples[spw].second;
                            }
                        }
                    }
                }
            }
        }
    }
}

void StatWtTVI::summarizeFlagging() const {
    auto orig = (Double)_nOrigFlaggedPts/(Double)_nTotalPts*100;
    auto stwt = (Double)_nNewFlaggedPts/(Double)_nTotalPts*100;
    auto total = orig + stwt;
    LogIO log(LogOrigin("StatWtTVI", __func__));
    log << LogIO::NORMAL << "Originally, " << orig
        << "% of the data were flagged. StatWtTVI flagged an "
        << "additional " << stwt << "%."  << LogIO::POST;
    log << LogIO::NORMAL << "TOTAL FLAGGED DATA AFTER RUNNING STATWT: "
        << total << "%" << LogIO::POST;
    log << LogIO::NORMAL << std::endl << LogIO::POST;
    String col0 = "SPECTRAL_WINDOW";
    String col1 = "SAMPLES_WITH_NON-ZERO_VARIANCE";
    String col2 = "SAMPLES_WHERE_REAL_PART_VARIANCE_DIFFERS_BY_>50%_FROM_IMAGINARY_PART";
    log << LogIO::NORMAL << col0 << " " << col1 << " " << col2 << LogIO::POST;
    auto n0 = col0.size();
    auto n1 = col1.size();
    auto n2 = col2.size();
    for (const auto& sample: _samples) {
        ostringstream oss;
        oss << std::setw(n0) << sample.first << " " << std::setw(n1)
            << sample.second.first << " " << std::setw(n2) << sample.second.second;
        log << LogIO::NORMAL << oss.str() << LogIO::POST;
    }
}

void StatWtTVI::origin() {
    // Drive underlying ViImplementation2
    getVii()->origin();
    // Synchronize own VisBuffer
    configureNewSubchunk();
    _clearCache();
}

void StatWtTVI::next() {
    // Drive underlying ViImplementation2
    getVii()->next();
    // Synchronize own VisBuffer
    configureNewSubchunk();
    _clearCache();
}

}

}
