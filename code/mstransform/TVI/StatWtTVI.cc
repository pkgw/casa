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

#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/ms/MSOper/MSMetaData.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    // FIXME when the TVI framework has methods to
    // check for metadata, like the existence of
    // columns, remove references to the original MS
    const auto& origMS = ms();
    ThrowIf(
        (_column == CORRECTED || _column == RESIDUAL)
        && ! origMS.isColumn(MSMainEnums::CORRECTED_DATA),
        "StatWtTVI requires the MS to have a CORRECTED_DATA column. This MS "
        "does not"
    );
    ThrowIf(
        (_column == DATA || _column == RESIDUAL_DATA)
        && ! origMS.isColumn(MSMainEnums::DATA),
        "StatWtTVI requires the MS to have a DATA column. This MS does not"
    );
    _mustComputeSigma = (_column == DATA || _column == RESIDUAL_DATA);
    _updateWeight = ! _mustComputeSigma 
        || (_mustComputeSigma && ! ms().isColumn(MSMainEnums::CORRECTED_DATA));
    _noModel = (_column == RESIDUAL || _column == RESIDUAL_DATA)
        && ! origMS.isColumn(MSMainEnums::MODEL_DATA)
        && ! origMS.source().isColumn(MSSourceEnums::SOURCE_MODEL);
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
                "Array specified in '" + field
                + "' must have exactly two values"
            );
            ThrowIf(
                casacore::anyLT(myrange, 0.0),
                "Both values specified in '" + field
                + "' array must be non-negative"
            );
            std::set<Double> rangeset(myrange.begin(), myrange.end());
            ThrowIf(
                rangeset.size() == 1, "Values specified in '" + field
                + "' array must be unique"
            );
            auto iter = rangeset.begin();
            _wtrange.reset(new std::pair<Double, Double>(*iter, *(++iter)));
        }
    }
    auto excludeChans = False;
    field = "excludechans";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpBool,
            "Unsupported type for field '" + field + "'"
        );
        excludeChans = config.asBool(field);
    }
    field = "fitspw";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpString,
            "Unsupported type for field '" + field + "'"
        );
        auto val = config.asString(field);
        if (! val.empty()) {
            // FIXME references underlying MS
            const auto& myms = ms();
            MSSelection sel(myms);
            sel.setSpwExpr(val);
            auto chans = sel.getChanList();
            auto nrows = chans.nrow();
            MSMetaData md(&myms, 50);
            auto nchans = md.nChans();
            IPosition start(3, 0);
            IPosition stop(3, 0);
            IPosition step(3, 1);
            for (uInt i=0; i<nrows; ++i) {
                auto row = chans.row(i);
                const auto& spw = row[0];
                if (_chanSelFlags.find(spw) == _chanSelFlags.end()) {
                    _chanSelFlags[spw]
                        = Cube<Bool>(1, nchans[spw], 1, ! excludeChans);
                }
                start[1] = row[1];
                ThrowIf(
                    start[1] < 0, "Invalid channel selection in spw "
                    + String::toString(spw))
                ;
                stop[1] = row[2];
                step[1] = row[3];
                Slicer slice(start, stop, step, Slicer::endIsLast);
                _chanSelFlags[spw](slice) = excludeChans;
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
                ! (
                    val.startsWith("c") || val.startsWith("d")
                    || val.startsWith("residual") || val.startsWith("residual_")
                ),
                "Unsupported value for " + field + ": " + val
            );
            _column = val.startsWith("c") ? CORRECTED
                : val.startsWith("d") ? DATA
                : val.startsWith("residual_") ? RESIDUAL_DATA
                : RESIDUAL;
        }
    }
    field = "slidetimebin";
    if (config.isDefined(field)) {
        ThrowIf(
            config.type(config.fieldNumber(field)) != TpBool,
            "Unsupported type for field '" + field + "'"
        );
        _timeBlockProcessing = ! config.asBool(field);
        if (! _timeBlockProcessing) {
            field = "timebin";
            ThrowIf(! config.isDefined(field), field + " not defined");
            auto mytype = config.type(config.fieldNumber(field));
            ThrowIf(
                ! (
                    mytype == TpString || mytype == TpDouble
                    || mytype == TpInt
                ),
                "Unsupported type for field '" + field + "'"
            );
            switch(mytype) {
            case TpDouble: {
                _slidingTimeWindowWidth = config.asDouble(field);
                checkTimeBinWidth(_slidingTimeWindowWidth);
                break;
            }
            case TpInt:
                _slidingTimeWindowWidth = getTimeBinWidthUsingInterval(
                    &ms(), config.asInt(field)
                );
                break;
            case TpString: {
                QuantumHolder qh(casaQuantity(config.asString(field)));
                _slidingTimeWindowWidth = getTimeBinWidthInSec(qh.asQuantity());
                break;
            }
            default:
                ThrowCc("Logic Error: Unhandled type");
            }
            LogIO log(LogOrigin("StatWtTVI", __func__));
            log << LogIO::NORMAL << "Using sliding time window of width "
                << _slidingTimeWindowWidth << " s" << LogIO::POST;

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
            casacore::StatisticsAlgorithmFactory<
                Double, Array<casacore::Float>::const_iterator,
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
                auto center = FitToHalfStatisticsData::CMEAN;
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
                auto ud = FitToHalfStatisticsData::LE_CENTER;
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
    std::set<StatisticsData::STATS> stats {StatisticsData::VARIANCE};
    _statAlg->setStatsToCalculate(stats);
    // also configure the _wtStats object here
    _wtStats.reset(
        new ClassicalStatistics<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator
        >()
    );
    stats.insert(StatisticsData::MEAN);
    _wtStats->setStatsToCalculate(stats);
    _wtStats->setCalculateAsAdded(True);
}

void StatWtTVI::_logUsedChannels() const {
    // FIXME uses underlying MS
    MSMetaData msmd(&ms(), 100.0);
    const auto nchan = msmd.nChans();
    // uInt nspw = nchan.size();
    LogIO log(LogOrigin("StatWtTVI", __func__));
    log << LogIO::NORMAL << "Weights are being computed using ";
    const auto cend = _chanSelFlags.cend();
    const auto nspw = _samples.size();
    uInt spwCount = 0;
    // for (uInt i=0; i<nspw; ++i) {
    for (const auto& kv: _samples) {
        const auto spw = kv.first;
        log << "SPW " << spw << ", channels ";
        const auto flagCube = _chanSelFlags.find(spw);
        if (flagCube == cend) {
            log << "0~" << (nchan[spw] - 1);
        }
        else {
            vector<pair<uInt, uInt>> startEnd;
            const auto flags = flagCube->second.tovector();
            bool started = false;
            std::unique_ptr<pair<uInt, uInt>> curPair;
            for (uInt j=0; j<nchan[spw]; ++j) {
                if (started) {
                    if (flags[j]) {
                        // found a bad channel, end current range
                        startEnd.push_back(*curPair);
                        started = false;
                    }
                    else {
                        // found a "good" channel, update end of current range
                        curPair->second = j;
                    }
                }
                else if (! flags[j]) {
                    // found a good channel, start new range
                    started = true;
                    curPair.reset(new pair<uInt, uInt>(j, j));
                }
            }
            if (curPair) {
                if (started) {
                    // The last pair won't get added inside the previous loop, 
                    // so add it here
                    startEnd.push_back(*curPair);
                }
                auto nPairs = startEnd.size();
                for (uInt i=0; i<nPairs; ++i) {
                    log  << startEnd[i].first << "~" << startEnd[i].second;
                    if (i < nPairs - 1) {
                        log << ", ";
                    }
                }
            }
            else {
                // if the pointer never got set, all the channels are bad
                log << "no channels";
            }
        }
        if (spwCount < (nspw - 1)) {
            log << ";";
        }
        ++spwCount;
    }
    log << LogIO::POST;
}

void StatWtTVI::_setChanBinMap(const casacore::Quantity& binWidth) {
    if (! binWidth.isConform(Unit("Hz"))) {
        ostringstream oss;
        oss << "If specified as a quantity, channel bin width must have "
            << "frequency units. " << binWidth << " does not.";
        ThrowCc(oss.str());
    }
    ThrowIf(binWidth.getValue() <= 0, "channel bin width must be positive");
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
            // both could be true, in which case both conditionals
            // must be executed
            if (abs(*citer - startFreq) > binWidthHz) {
                // add bin to list
                bin.end = chanNum - 1;
                _chanBins[i].push_back(bin);
                bin.start = chanNum;
                startFreq = *citer;
            }
            if (chanNum + 1 == nchan) {
                // add last bin
                bin.end = chanNum;
                _chanBins[i].push_back(bin);
            }
        }
    }
    // weight spectrum must be computed
    _mustComputeWtSp.reset(new Bool(True));
}

void StatWtTVI::_setChanBinMap(Int binWidth) {
    ThrowIf(binWidth < 1, "Channel bin width must be positive");
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
    // weight spectrum must be computed
    _mustComputeWtSp.reset(new Bool(True));
}

void StatWtTVI::_setDefaultChanBinMap() {
    MSMetaData msmd(&ms(), 0.0);
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

Double StatWtTVI::getTimeBinWidthInSec(const casacore::Quantity& binWidth) {
    ThrowIf(
        ! binWidth.isConform(Unit("s")),
        "Time bin width unit must be a unit of time"
    );
    auto v = binWidth.getValue("s");
    checkTimeBinWidth(v);
    return v;
}

void StatWtTVI::checkTimeBinWidth(Double binWidth) {
    ThrowIf(binWidth <= 0, "time bin width must be positive");
}

Double StatWtTVI::getTimeBinWidthUsingInterval(
    const MeasurementSet *const ms, Int n
) {
    ThrowIf(n <= 0, "number of time intervals must be positive");
    MSMetaData msmd(ms, 0.0);
    auto stats = msmd.getIntervalStatistics();
    ThrowIf(
        stats.max/stats.median - 1 > 0.25 || 1 - stats.min/stats.median > 0.25,
        "There is not a representative integration time in the INTERVAL "
        "column, likely due to different visibility integration times across "
        "different scans. Please select only parts of the MeasurementSet that "
        "have a uniform integration time for each execution of statwt. "
        "Multiple statwt executions may be needed to cover all MS rows with "
        "different integration (INTERVAL) values."
    );
    return n*stats.median;
}

void StatWtTVI::sigmaSpectrum(Cube<Float>& sigmaSp) const {
    if (_mustComputeSigma) {
        {
            Cube<Float> wtsp;
            // this computes _newWtsp, ignore wtsp
            weightSpectrum(wtsp);
        }
        sigmaSp = Float(1.0)/sqrt(_newWtSp);
        if (anyEQ(_newWtSp, Float(0))) {
            auto iter = sigmaSp.begin();
            auto end = sigmaSp.end();
            auto witer = _newWtSp.cbegin();
            for ( ; iter != end; ++iter, ++witer) {
                if (*witer == 0) {
                    *iter = -1;
                }
            }
        }
    }
    else {
        TransformingVi2::sigmaSpectrum(sigmaSp);
    }
}

void StatWtTVI::weightSpectrum(Cube<Float>& newWtsp) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! *_mustComputeWtSp) {
        newWtsp.resize(IPosition(3, 0));
        return;
    }
    if (! _newWtSp.empty()) {
        // already calculated
        if (_updateWeight) {
            newWtsp = _newWtSp.copy();
        }
        else {
            TransformingVi2::weightSpectrum(newWtsp);
        }
        return;
    }
    _computeWeightSpectrumAndFlags();
    if (_updateWeight) {
        newWtsp = _newWtSp.copy();
    }
    else {
        TransformingVi2::weightSpectrum(newWtsp);
    }
}

void StatWtTVI::_computeWeightSpectrumAndFlags() const {
    size_t nOrigFlagged;
    auto mypair = _getLowerLayerWtSpFlags(nOrigFlagged);
    auto& wtsp = mypair.first;
    auto& flagCube = mypair.second;
    if (*_mustComputeWtSp && wtsp.empty()) {
        // This can happen in preview mode if
        // WEIGHT_SPECTRUM doesn't exist or is empty
        wtsp.resize(flagCube.shape());
    }
    auto checkFlags = False;
    if (_timeBlockProcessing) {
        _weightSpectrumFlagsTimeBlockProcessing(wtsp, flagCube, checkFlags);
    }
    else {
        _weightSpectrumFlagsSlidingTimeWindow(wtsp, flagCube, checkFlags);
    }
    if (checkFlags) {
        _nNewFlaggedPts += ntrue(flagCube) - nOrigFlagged;
    }
    _newWtSp = wtsp;
    _newFlag = flagCube;
}

void StatWtTVI::_weightSpectrumFlagsSlidingTimeWindow(
    Cube<Float>& wtsp, Cube<Bool>& flagCube, Bool& checkFlags
) const {
    // fish out the rows relevant to this subchunk
    Vector<uInt> rowIDs;
    getRowIds(rowIDs);
    auto start = _rowIDInMSTorowIndexInChunk.find(*rowIDs.begin());
    ThrowIf(
        start == _rowIDInMSTorowIndexInChunk.end(),
        "Logic Error: Cannot find requested subchunk in stored chunk"
    );
    // this is the row index in the chunk
    auto chunkRowIndex = start->second;
    auto chunkRowEnd = chunkRowIndex + nRows();
    Slicer slice(IPosition(3, 0), flagCube.shape(), Slicer::endIsLength);
    auto sliceStart = slice.start();
    auto sliceEnd = slice.end();
    auto nCorrBins = _combineCorr ? 1 : flagCube.shape()[0];
    Vector<Int> spws;
    spectralWindows(spws);
    auto spw = *spws.begin();
    const auto& chanBins = _chanBins.find(spw)->second;
    auto subChunkRowIndex = 0;
    for (; chunkRowIndex < chunkRowEnd; ++chunkRowIndex, ++subChunkRowIndex) {
        sliceStart[2] = subChunkRowIndex;
        sliceEnd[2] = subChunkRowIndex;
        auto iChanBin = 0;
        for (const auto& chanBin: chanBins) {
            sliceStart[1] = chanBin.start;
            sliceEnd[1] = chanBin.end;
            auto corr = 0;
            for (; corr < nCorrBins; ++corr) {
                if (! _combineCorr) {
                    sliceStart[0] = corr;
                    sliceEnd[0] = corr;
                }
                slice.setStart(sliceStart);
                slice.setEnd(sliceEnd);
                _updateWtSpFlags(
                    wtsp, flagCube, checkFlags, slice,
                    _slidingTimeWindowWeights(corr, iChanBin, chunkRowIndex)
                );
            }
            ++iChanBin;
        }
    }
}

void StatWtTVI::_weightSpectrumFlagsTimeBlockProcessing(
    Cube<Float>& wtsp, Cube<Bool>& flagCube, Bool& checkFlags
) const {
    Vector<Int> ant1, ant2, spws;
    antenna1(ant1);
    antenna2(ant2);
    spectralWindows(spws);
    Slicer slice(IPosition(3, 0), flagCube.shape(), Slicer::endIsLength);
    auto sliceStart = slice.start();
    auto sliceEnd = slice.end();
    auto nrows = nRows();
    for (Int i=0; i<nrows; ++i) {
        sliceStart[2] = i;
        sliceEnd[2] = i;
        BaselineChanBin blcb;
        blcb.baseline = _baseline(ant1[i], ant2[i]);
        auto spw = spws[i];
        blcb.spw = spw;
        auto bins = _chanBins.find(spw)->second;
        for (const auto& bin: bins) {
            sliceStart[1] = bin.start;
            sliceEnd[1] = bin.end;
            blcb.chanBin = bin;
            auto weights = _weights.find(blcb)->second;
            auto ncorr = weights.size();
            for (uInt corr=0; corr<ncorr; ++corr) {
                if (! _combineCorr) {
                    sliceStart[0] = corr;
                    sliceEnd[0] = corr;
                }
                slice.setStart(sliceStart);
                slice.setEnd(sliceEnd);
                _updateWtSpFlags(
                    wtsp, flagCube, checkFlags, slice, weights[corr]
                );
            }
        }
    }
}

void StatWtTVI::_updateWtSpFlags(
    Cube<Float>& wtsp, Cube<Bool>& flags, Bool& checkFlags,
    const Slicer& slice, Float wt
) const {
    // writable array reference
    auto flagSlice = flags(slice);
    if (*_mustComputeWtSp) {
        // writable array reference
        auto wtSlice = wtsp(slice);
        wtSlice = wt;
        // update global stats before we potentially flag data
        auto mask = ! flagSlice;
        _wtStats->addData(wtSlice.begin(), mask.begin(), wtSlice.size());
    }
    else if (! allTrue(flagSlice)) {
        // we don't need to compute WEIGHT_SPECTRUM, and the slice isn't
        // entirely flagged, so we need to update the WEIGHT column stats
        _wtStats->addData(Array<Float>(IPosition(1, 1), wt).begin(), 1);
    }
    if (
        wt == 0
        || (_wtrange && (wt < _wtrange->first || wt > _wtrange->second))
    ) {
        checkFlags = True;
        flagSlice = True;
    }
}

std::pair<Cube<Float>, Cube<Bool>> StatWtTVI::_getLowerLayerWtSpFlags(
    size_t& nOrigFlagged
) const {
    auto mypair = std::make_pair(Cube<Float>(), Cube<Bool>());
    if (*_mustComputeWtSp) {
        getVii()->weightSpectrum(mypair.first);
    }
    getVii()->flag(mypair.second);
    _nTotalPts += mypair.second.size();
    nOrigFlagged = ntrue(mypair.second);
    _nOrigFlaggedPts += nOrigFlagged;
    return mypair;
}

void StatWtTVI::sigma(Matrix<Float>& sigmaMat) const {
    if (_mustComputeSigma) {
        if (_newWt.empty()) {
            Matrix<Float> wtmat;
            weight(wtmat);
        }
        sigmaMat = Float(1.0)/sqrt(_newWt);
        if (anyEQ(_newWt, Float(0))) {
            Matrix<Float>::iterator iter = sigmaMat.begin();
            Matrix<Float>::iterator end = sigmaMat.end();
            Matrix<Float>::iterator witer = _newWt.begin();
            for ( ; iter != end; ++iter, ++witer) {
                if (*witer == 0) {
                    *iter = -1;
                }
            }
        }
    }
    else {
        TransformingVi2::sigma(sigmaMat);
    }
}

void StatWtTVI::weight(Matrix<Float> & wtmat) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! _newWt.empty()) {
        if (_updateWeight) {
            wtmat = _newWt.copy();
        }
        else {
            TransformingVi2::weight(wtmat);
        }
        return;
    }
    auto nrows = nRows();
    getVii()->weight(wtmat);
    if (*_mustComputeWtSp) {
        // always use classical algorithm to get median for weights
        ClassicalStatistics<Double, Array<Float>::const_iterator, Array<Bool>::const_iterator> cs;
        Cube<Float> wtsp;
        Cube<Bool> flagCube;
        // this computes _newWtsP which is what we will use, so
        // just ignore wtsp
        weightSpectrum(wtsp);
        flag(flagCube);
        IPosition blc(3, 0);
        IPosition trc = _newWtSp.shape() - 1;
        const auto ncorr = _newWtSp.shape()[0];
        for (Int i=0; i<nrows; ++i) {
            blc[2] = i;
            trc[2] = i;
            if (_combineCorr) {
                auto flags = flagCube(blc, trc);
                if (allTrue(flags)) {
                    wtmat.column(i) = 0;
                }
                else {
                    auto weights = _newWtSp(blc, trc);
                    auto mask = ! flags;
                    cs.setData(weights.begin(), mask.begin(), weights.size());
                    wtmat.column(i) = cs.getMedian();
                }
            }
            else {
                for (uInt corr=0; corr<ncorr; ++corr) {
                    blc[0] = corr;
                    trc[0] = corr;
                    auto weights = _newWtSp(blc, trc);
                    auto flags = flagCube(blc, trc);
                    if (allTrue(flags)) {
                        wtmat(corr, i) = 0;
                    }
                    else {
                        auto mask = ! flags;
                        cs.setData(
                            weights.begin(), mask.begin(), weights.size()
                        );
                        wtmat(corr, i) = cs.getMedian();
                    }
                }
            }
        }
    }
    else {
        // the only way this can happen is if there is a single channel bin
        // for each baseline/spw pair
        if (_timeBlockProcessing) {
            _weightSingleChanBinBlockTimeProcessing(wtmat, nrows);
        }
        else {
            _weightSingleChanBinSlidingTimeWindow(wtmat, nrows);
        }
    }
    _newWt = wtmat.copy();
    if (! _updateWeight) {
        wtmat = Matrix<Float>(wtmat.shape()); 
        TransformingVi2::weight(wtmat);
    }
}

void StatWtTVI::_weightSingleChanBinBlockTimeProcessing(
    Matrix<Float>& wtmat, Int nrows
) const {
    Vector<Int> ant1, ant2, spws;
    antenna1(ant1);
    antenna2(ant2);
    spectralWindows(spws);
    // all spws in a given chunk are guaranteed to be the same
    auto spw = *spws.begin();
    BaselineChanBin blcb;
    blcb.spw = spw;
    for (Int i=0; i<nrows; ++i) {
        auto bins = _chanBins.find(spw)->second;
        blcb.baseline = _baseline(ant1[i], ant2[i]);
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

void StatWtTVI::_weightSingleChanBinSlidingTimeWindow(
    Matrix<Float>& wtmat, Int nrows
) const {
    Vector<uInt> rowIDs;
    getRowIds(rowIDs);
    auto start = _rowIDInMSTorowIndexInChunk.find(*rowIDs.begin());
    ThrowIf(
        start == _rowIDInMSTorowIndexInChunk.end(),
        "Logic Error: Cannot find requested subchunk in stored chunk"
    );
    // this is the row index in the chunk
    auto chunkRowIndex = start->second;
    auto ncorr = wtmat.nrow();
    for (Int i=0; i<nrows; ++i) {
        if (_combineCorr) {
            wtmat.column(i) = _slidingTimeWindowWeights(0, 0, chunkRowIndex);
        }
        else {
            for (uInt corr=0; corr<ncorr; ++corr) {
                wtmat(corr, i) = _slidingTimeWindowWeights(
                    corr, 0, chunkRowIndex
                );
            }
        }
    }
}

void StatWtTVI::flag(Cube<Bool>& flagCube) const {
    ThrowIf(! _weightsComputed, "Weights have not been computed yet");
    if (! _newFlag.empty()) {
        flagCube = _newFlag.copy();
        return;
    }
    _computeWeightSpectrumAndFlags();
    flagCube = _newFlag.copy();
}

void StatWtTVI::flagRow(Vector<Bool>& flagRow) const {
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

void StatWtTVI::_gatherAndComputeWeightsSlidingTimeWindow() const {
    ViImplementation2* vii = getVii();
    VisBuffer2* vb = vii->getVisBuffer();
    // map of rowID to rowIDs that should be used to compute
    // the weight for the key rowID
    std::vector<std::set<uInt>> rowMap;
    auto firstTime = True;
    // each subchunk is guaranteed to represent exactly one time stamp
    std::vector<Double> subChunkToTimeStamp;
    // Baseline-subchunk number pair to row index in that chunk
    std::map<std::pair<Baseline, uInt>, uInt> baselineSubChunkToIndex;
    uInt subChunkID = 0;
    auto halfWidth = _slidingTimeWindowWidth/2;
    // we build up the chunk data in the chunkData and chunkFlags cubes
    Cube<Complex> chunkData;
    Cube<Bool> chunkFlags;
    uInt subchunkStartIndex = 0;
    auto initChanSelTemplate = True;
    Cube<Bool> chanSelFlagTemplate, chanSelFlags;
    // we cannot know the spw until inside the subchunk loop
    Int spw = -1;
    _rowIDInMSTorowIndexInChunk.clear();
    Slicer sl(IPosition(3, 0), IPosition(3, 1));
    auto slStart = sl.start();
    auto slEnd = sl.end();
    for (vii->origin(); vii->more(); vii->next(), ++subChunkID) {
        if (_checkFirsSubChunk(spw, firstTime, vb)) {
            return;
        }
        if (! _mustComputeWtSp) {
            _mustComputeWtSp.reset(
                new Bool(
                    vb->existsColumn(VisBufferComponent2::WeightSpectrum)
                )
            );
        }
        _rowIDInMSTorowIndexInChunk[*vb->rowIds().begin()] = subchunkStartIndex;
        const auto& ant1 = vb->antenna1();
        const auto& ant2 = vb->antenna2();
        // [nC,nF,nR)
        const auto nrows = vb->nRows();
        // there is no guarantee a previous subchunk will be included,
        // eg if the timewidth is small enough
        Int firstIncludedSubChunk = -1;
        auto subchunkTime = vb->time()[0];
        if (subChunkID > 0) {
            auto siter = subChunkToTimeStamp.rbegin();
            auto send = subChunkToTimeStamp.rend();
            uInt idx = subChunkID - 1;
            for (; siter!=send; ++siter, --idx) {
                if (subchunkTime - *siter <= halfWidth) {
                    firstIncludedSubChunk = idx;
                }
                else {
                    // should be time sorted, so we can break
                    // out of the loop here
                    break;
                }
            }
        }
        auto rowInChunk = subchunkStartIndex;
        for (Int row=0; row<nrows; ++row, ++rowInChunk) {
            const auto baseline = _baseline(ant1[row], ant2[row]);
            std::pair<Baseline, uInt> mypair(baseline, subChunkID);
            baselineSubChunkToIndex[mypair] = rowInChunk;
            std::set<uInt> myRowNums;
            myRowNums.insert(rowInChunk);
            if (subChunkID > 0 && firstIncludedSubChunk >= 0) {
                auto subchunk = firstIncludedSubChunk;
                for (; subchunk < (Int)subChunkID; ++subchunk) {
                    auto myend = baselineSubChunkToIndex.end();
                    mypair.second = subchunk;
                    if (baselineSubChunkToIndex.find(mypair) != myend) {
                        auto existingRowNum = baselineSubChunkToIndex[mypair];
                        myRowNums.insert(existingRowNum);
                        rowMap[existingRowNum].insert(rowInChunk);
                    }
                }
            }
            rowMap.push_back(myRowNums);
        }
        const auto dataCube = _dataCube(vb);
        auto resultantFlags = _getResultantFlags(
            chanSelFlagTemplate, chanSelFlags, initChanSelTemplate,
            spw, vb->flagCube()
        );
        // build up chunkData and chunkFlags one subchunk at a time
        if (chunkData.empty()) {
            chunkData = dataCube;
            chunkFlags = resultantFlags;
        }
        else {
            auto newShape = chunkData.shape();
            newShape[2] += nrows;
            chunkData.resize(newShape, True);
            chunkFlags.resize(newShape, True);
            slStart[2] = subchunkStartIndex;
            sl.setStart(slStart);
            slEnd = newShape - 1;
            sl.setEnd(slEnd);
            chunkData(sl) = dataCube;
            chunkFlags(sl) = resultantFlags;
        }
        subChunkToTimeStamp.push_back(subchunkTime);
        subchunkStartIndex += nrows;
    }
    _computeWeightsSlidingTimeWindow(chunkData, chunkFlags, rowMap, spw);
}

const Cube<Complex> StatWtTVI::_dataCube(const VisBuffer2 *const vb) const {
    switch (_column) {
    case CORRECTED:
        return vb->visCubeCorrected();
    case DATA:
        return vb->visCube();
    case RESIDUAL:
        if (_noModel) {
            return vb->visCubeCorrected();
        }
        else {
            return vb->visCubeCorrected() - vb->visCubeModel();
        }
    case RESIDUAL_DATA:
        if(_noModel) {
            return vb->visCube();
        }
        else {
            return vb->visCube() - vb->visCubeModel();
        }
    default:
        ThrowCc("Logic error: column type not handled");
    }
}

void StatWtTVI::_computeWeightsSlidingTimeWindow(
    const Cube<Complex>& data, const Cube<Bool>& flags,
    const std::vector<std::set<uInt>>& rowMap, uInt spw
) const {
    auto chunkShape = data.shape();
    const auto nActCorr = chunkShape[0];
    const auto ncorr = _combineCorr ? 1 : nActCorr;
    const auto& chanBins = _chanBins.find(spw)->second;
    _slidingTimeWindowWeights.resize(
        IPosition(3, ncorr, chanBins.size(), chunkShape[2]),
        False
    );
    const auto nRows = rowMap.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t iRow=0; iRow<nRows; ++iRow) {
        IPosition chunkSliceStart(3, 0);
        auto chunkSliceLength = chunkShape;
        chunkSliceLength[2] = 1;
        Slicer chunkSlice(
            chunkSliceStart, chunkSliceLength, Slicer::endIsLength
        );
        auto chunkSliceEnd = chunkSlice.end();
        auto appendingSlice = chunkSlice;
        auto appendingSliceStart = appendingSlice.start();
        auto appendingSliceEnd = appendingSlice.end();
        auto intraChunkSlice = appendingSlice;
        auto intraChunkSliceStart = intraChunkSlice.start();
        auto intraChunkSliceEnd = intraChunkSlice.end();
        intraChunkSliceEnd[0] = nActCorr - 1;
        const auto& rowsToInclude = rowMap[iRow];
        auto dataShape = chunkShape;
        dataShape[2] = rowsToInclude.size();
        Cube<Complex> dataArray(dataShape);
        Cube<Bool> flagArray(dataShape);
        auto siter = rowsToInclude.begin();
        auto send = rowsToInclude.end();
        uInt n = 0;
        // create an array with only the rows that should
        // be used in the computation of weights for the
        // current row
        for (; siter!=send; ++siter, ++n) {
            appendingSliceStart[2] = n;
            appendingSlice.setStart(appendingSliceStart);
            appendingSliceEnd[2] = n;
            appendingSlice.setEnd(appendingSliceEnd);
            chunkSliceStart[2] = *siter;
            chunkSlice.setStart(chunkSliceStart);
            chunkSliceEnd[2] = *siter;
            chunkSlice.setEnd(chunkSliceEnd);
            dataArray(appendingSlice) = data(chunkSlice);
            flagArray(appendingSlice) = flags(chunkSlice);
        }
        // slice up for correlations and channel binning
        intraChunkSliceEnd[2] = dataShape[2] - 1;
        for (uInt corr=0; corr<ncorr; ++corr) {
            if (! _combineCorr) {
                intraChunkSliceStart[0] = corr;
                intraChunkSliceEnd[0] = corr;
            }
            auto citer = chanBins.begin();
            auto cend = chanBins.end();
            auto iChanBin = 0;
            for (; citer!=cend; ++citer, ++iChanBin) {
                intraChunkSliceStart[1] = citer->start;
                intraChunkSliceEnd[1] = citer->end;
                intraChunkSlice.setStart(intraChunkSliceStart);
                intraChunkSlice.setEnd(intraChunkSliceEnd);
                _slidingTimeWindowWeights(corr, iChanBin, iRow)
                    = _computeWeight(
                        dataArray(intraChunkSlice), flagArray(intraChunkSlice),
                        spw
                    );
            }
        }
    }
}

void StatWtTVI::_gatherAndComputeWeights() const {
    if (_timeBlockProcessing) {
        _gatherAndComputeWeightsTimeBlockProcessing();
    }
    else {
        _gatherAndComputeWeightsSlidingTimeWindow();
    }
}

void StatWtTVI::_gatherAndComputeWeightsTimeBlockProcessing() const {
    // Drive NEXT LOWER layer's ViImpl to gather data into allvis:
    //  Assumes all sub-chunks in the current chunk are to be used
    //   for the variance calculation
    //  Essentially, we are sorting the incoming data into
    //   allvis, to enable a convenient variance calculation
    _weights.clear();
    auto* vii = getVii();
    auto* vb = vii->getVisBuffer();
    std::map<BaselineChanBin, Cube<Complex>> data;
    std::map<BaselineChanBin, Cube<Bool>> flags;
    IPosition blc(3, 0);
    auto trc = blc;
    auto initChanSelTemplate = True;
    Cube<Bool> chanSelFlagTemplate, chanSelFlags;
    auto firstTime = True;
    // we cannot know the spw until we are in the subchunks loop
    Int spw = -1;
    for (vii->origin(); vii->more(); vii->next()) {
        if (_checkFirsSubChunk(spw, firstTime, vb)) {
            return;
        }
        if (! _mustComputeWtSp) {
            _mustComputeWtSp.reset(
                new Bool(
                    vb->existsColumn(VisBufferComponent2::WeightSpectrum)
                )
            );
        }
        const auto& ant1 = vb->antenna1();
        const auto& ant2 = vb->antenna2();
        // [nC,nF,nR)
        const auto& dataCube = _dataCube(vb);
        const auto& flagCube = vb->flagCube();
        const auto nrows = vb->nRows();
        const auto npol = dataCube.nrow();
        auto resultantFlags = _getResultantFlags(
            chanSelFlagTemplate, chanSelFlags, initChanSelTemplate,
            spw, flagCube
        );
        auto bins = _chanBins.find(spw)->second;
        BaselineChanBin blcb;
        blcb.spw = spw;
        IPosition dataCubeBLC(3, 0);
        auto dataCubeTRC = dataCube.shape() - 1;
        dataCubeTRC[2] = 0;
        for (Int row=0; row<nrows; ++row) {
            dataCubeTRC[2] = 0;
            dataCubeBLC[2] = row;
            dataCubeTRC[2] = row;
            blcb.baseline = _baseline(ant1[row], ant2[row]);
            auto citer = bins.begin();
            auto cend = bins.end();
            for (; citer!=cend; ++citer) {
                dataCubeBLC[1] = citer->start;
                dataCubeTRC[1] = citer->end;
                blcb.chanBin.start = citer->start;
                blcb.chanBin.end = citer->end;
                auto dataSlice = dataCube(dataCubeBLC, dataCubeTRC);
                auto flagSlice = resultantFlags(dataCubeBLC, dataCubeTRC);
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
    _computeWeightsTimeBlockProcessing(data, flags);
}

Cube<Bool> StatWtTVI::_getResultantFlags(
    Cube<Bool>& chanSelFlagTemplate, Cube<Bool>& chanSelFlags,
    Bool& initTemplate, Int spw, const Cube<Bool>& flagCube
) const {
    if (_chanSelFlags.find(spw) == _chanSelFlags.cend()) {
        // no selection of channels to ignore
        return flagCube;
    }
    if (initTemplate) {
        // this can be done just once per chunk because all the rows
        // in the chunk are guaranteed to have the same spw
        // because each subchunk is guaranteed to have a single
        // data description ID.
        chanSelFlagTemplate = _chanSelFlags.find(spw)->second;
        initTemplate = False;
    }
    auto dataShape = flagCube.shape();
    chanSelFlags.resize(dataShape, False);
    auto ncorr = dataShape[0];
    auto nrows = dataShape[2];
    IPosition start(3, 0);
    IPosition end = dataShape-1;
    Slicer sl(start, end, Slicer::endIsLast);
    for (uInt corr=0; corr<ncorr; ++corr) {
        start[0] = corr;
        end[0] = corr;
        for (Int row=0; row<nrows; ++row) {
            start[2] = row;
            end[2] = row;
            sl.setStart(start);
            sl.setEnd(end);
            chanSelFlags(sl) = chanSelFlagTemplate;
        }
    }
    return flagCube || chanSelFlags;
}

Bool StatWtTVI::_checkFirsSubChunk(
    Int& spw, Bool& firstTime, const VisBuffer2 * const vb
) const {
    if (! firstTime) {
        // this chunk has already been checked, it has not
        // been processed previously
        return False;
    }
    const auto& rowIDs = vb->rowIds();

    if (_processedRowIDs.find(rowIDs[0]) == _processedRowIDs.end()) {
        // haven't processed this chunk
        _processedRowIDs.insert(rowIDs[0]);
        // the spw is the same for all subchunks, so it only needs to
        // be set once
        spw = *vb->spectralWindows().begin();
        if (_samples.find(spw) == _samples.end()) {
            _samples[spw].first = 0;
            _samples[spw].second = 0;
        }
        firstTime = False;
        return False;
    }
    else {
        // this chunk has been processed, this can happen at the end
        // when the last chunk is processed twice
        return True;
    }
}

void StatWtTVI::initWeightSpectrum (const Cube<Float>& wtspec) {
    // Pass to next layer down
    getVii()->initWeightSpectrum(wtspec);
}

void StatWtTVI::initSigmaSpectrum (const Cube<Float>& sigspec) {
    // Pass to next layer down
    getVii()->initSigmaSpectrum(sigspec);
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

void StatWtTVI::_computeWeightsTimeBlockProcessing(
    const map<BaselineChanBin, Cube<Complex>>& data,
    const map<BaselineChanBin, Cube<Bool>>& flags
) const {
    auto diter = data.cbegin();
    auto dend = data.cend();
    const auto nActCorr = diter->second.shape()[0];
    const auto ncorr = _combineCorr ? 1 : nActCorr;
    // spw will be the same for all members
    const auto& spw = data.begin()->first.spw;
    std::vector<BaselineChanBin> keys(data.size());
    auto idx = 0;
    for (; diter!=dend; ++diter, ++idx) {
        const auto& blcb = diter->first;
        keys[idx] = blcb;
        _weights[blcb].resize(ncorr);
    }
    auto n = keys.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t i=0; i<n; ++i) {
        auto blcb = keys[i];
        auto dataForBLCB = data.find(blcb)->second;
        auto flagsForBLCB = flags.find(blcb)->second;
        for (uInt corr=0; corr<ncorr; ++corr) {
            IPosition start(3, 0);
            auto end = dataForBLCB.shape() - 1;
            if (! _combineCorr) {
                start[0] = corr;
                end[0] = corr;
            }
            Slicer slice(start, end, Slicer::endIsLast);
            _weights[blcb][corr] = _computeWeight(
                dataForBLCB(slice), flagsForBLCB(slice), spw
            );
        }
    }
}

casacore::Double StatWtTVI::_computeWeight(
    const Cube<Complex>& data, const Cube<Bool>& flags, uInt spw
) const {
    const auto npts = data.size();
    if ((Int)npts < _minSamp || (Int)nfalse(flags) < _minSamp) {
        // not enough points, trivial
        return 0;
    }
    // called in multi-threaded mode
    std::unique_ptr<
        StatisticsAlgorithm<
            Double, Array<Float>::const_iterator,
            Array<Bool>::const_iterator
        >
    > statAlg(_statAlg->clone());
    // some data not flagged
    const auto realPart = real(data);
    const auto imagPart = imag(data);
    const auto mask = ! flags;
    const auto riter = realPart.begin();
    const auto iiter = imagPart.begin();
    const auto miter = mask.begin();
    statAlg->setData(riter, miter, npts);
    auto realVar = statAlg->getStatistic(StatisticsData::VARIANCE);
    statAlg->setData(iiter, miter, npts);
    auto imagVar = statAlg->getStatistic(StatisticsData::VARIANCE);
    auto varSum = realVar + imagVar;
    // _samples.second can be updated in two different places, so use
    // a local (per thread) variable and update the object's private field in one
    // place
    uInt updateSecond = False;
    if (varSum > 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++_samples[spw].first;
        if (imagVar == 0 || realVar == 0) {
            updateSecond = True;
        }
        else {
            auto ratio = imagVar/realVar;
            auto inverse = 1/ratio;
            updateSecond = ratio > 1.5 || inverse > 1.5;
        }
        if (updateSecond) {
#ifdef _OPENMP
#pragma omp atomic
#endif
            ++_samples[spw].second;
        }
    }
    return varSum == 0 ? 0 : 2/varSum;
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
    if (_nOrigFlaggedPts == _nTotalPts) {
        log << LogIO::WARN << "IT APPEARS THAT ALL THE DATA IN THE INPUT "
            << "MS/SELECTION WERE FLAGGED PRIOR TO RUNNING STATWT"
            << LogIO::POST;
        log << LogIO::NORMAL << std::endl << LogIO::POST;
    }
    else if (_nOrigFlaggedPts + _nNewFlaggedPts == _nTotalPts) {
        log << LogIO::WARN << "IT APPEARS THAT STATWT FLAGGED ALL THE DATA "
            "IN THE REQUESTED SELECTION THAT WASN'T ORIGINALLY FLAGGED"
            << LogIO::POST;
        log << LogIO::NORMAL << std::endl << LogIO::POST;
    }
    String col0 = "SPECTRAL_WINDOW";
    String col1 = "SAMPLES_WITH_NON-ZERO_VARIANCE";
    String col2 = "SAMPLES_WHERE_REAL_PART_VARIANCE_DIFFERS_BY_>50%_FROM_"
        "IMAGINARY_PART";
    log << LogIO::NORMAL << col0 << " " << col1 << " " << col2 << LogIO::POST;
    auto n0 = col0.size();
    auto n1 = col1.size();
    auto n2 = col2.size();
    for (const auto& sample: _samples) {
        ostringstream oss;
        oss << std::setw(n0) << sample.first << " " << std::setw(n1)
            << sample.second.first << " " << std::setw(n2)
            << sample.second.second;
        log << LogIO::NORMAL << oss.str() << LogIO::POST;
    }
}

void StatWtTVI::summarizeStats(Double& mean, Double& variance) const {
    LogIO log(LogOrigin("StatWtTVI", __func__));
    _logUsedChannels();
    try {
        mean = _wtStats->getStatistic(StatisticsData::MEAN);
        variance = _wtStats->getStatistic(StatisticsData::VARIANCE);
        log << LogIO::NORMAL << "The mean of the computed weights is "
            << mean << LogIO::POST;
        log << LogIO::NORMAL << "The variance of the computed weights is "
            << variance << LogIO::POST;
        log << LogIO::NORMAL << "Weights which had corresponding flags of True "
            << "prior to running this application were not used to compute these "
            << "stats." << LogIO::POST;
    }
    catch (const AipsError& x) {
        log << LogIO::WARN << "There was a problem calculating the mean and "
            << "variance of the weights computed by this application. Perhaps there "
            << "was something amiss with the input MS and/or the selection criteria. "
            << "Examples of such issues are that all the data were originally flagged "
            << "or that the sample size was consistently too small for computations "
            << "of variances" << LogIO::POST;
        setNaN(mean);
        setNaN(variance);
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
