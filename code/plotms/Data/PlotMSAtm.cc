//# PlotMSAtm.cc: Implementation of PlotMSAtm.h
//# Copyright (C) 2017 
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

#include <plotms/Data/PlotMSAtm.h>

#include <iomanip>

#include <casa/OS/Path.h>
#include <casa/OS/File.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayUtil.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MeasTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <synthesis/CalTables/CTColumns.h>
#include <atmosphere/ATM/ATMSpectralGrid.h>
#include <atmosphere/ATM/ATMRefractiveIndexProfile.h>
#include <atmosphere/ATM/ATMSkyStatus.h>

using namespace casacore;

namespace casa {

PlotMSAtm::PlotMSAtm(casacore::String filename, PlotMSSelection& userSel,
        bool showatm, bool isMS, bool xAxisIsChan, PlotMSCacheBase* parent):
    filename_(filename),
    tableName_(""),
    telescopeName_(""),
    showatm_(showatm),
    isMS_(isMS),
    xIsChan_(xAxisIsChan),
    canCalculatePwv_(false),
    canCalculateWeather_(false),
    parent_(parent),
    selection_(userSel),
    ms_(nullptr),
    caltable_(nullptr),
    selectedSpw_(-1),
    selectedScan_(-1),
    pwv_(0.0),
    airmass_(0.0),
    MAX_ATM_CALC_CHAN_(512) {

    if (isMS) { // create MS
        setUpMS(filename, userSel);
    } else { // create NCT and MS if possible
        setUpCalTable(filename, userSel);
    }

    // check and log this once
    if (!tableName_.empty()) {
        casacore::Table mstab(tableName_);
        if (mstab.keywordSet().fieldNumber("ASDM_CALWVR") > -1) {
            canCalculatePwv_ = true;
        }
        if (mstab.keywordSet().fieldNumber("WEATHER") > -1) {
            canCalculateWeather_ = true;
        }
    }
    if (!canCalculatePwv_) {
        parent_->logmesg("load_cache",
            "No ASDM_CALWVR or ASDM_CALATMOSPHERE table to calculate pwv.");
        parent_->logmesg("load_cache",
            "Using default pwv for telescope " + telescopeName_ + ".");
    }
    if (!canCalculateWeather_) {
        parent_->logmesg("load_cache",
            "No WEATHER table, using default weather values for Atm/Tsky.");
    } 

    selms_ = new MeasurementSet(); // selected per chunk
    selct_ = new NewCalTable();    // selected per chunk
}

PlotMSAtm::~PlotMSAtm() {
    if (caltable_) {
        delete caltable_;
        caltable_ = nullptr;
    }
    if (ms_) {
        delete ms_;
        ms_ = nullptr;
    }
    if (selms_) {
        delete selms_;
        selms_ = nullptr;
    }
    if (selct_) {
        delete selct_;
        selct_ = nullptr;
    }
}

void PlotMSAtm::setUpMS(casacore::String filename, PlotMSSelection& userSel) {
    // create MeasurementSet, apply user selection, get telescope name
    try {
        ms_ = new MeasurementSet(filename);
    } catch (AipsError& err) {
        throw(AipsError("MeasurementSet setup failed.\n" + err.getMesg()));
        return;
    }
    tableName_ = ms_->tableName();

    // apply user-selection to ms_
    if (!userSel.isEmpty()) { 
        applyMSSelection(userSel, *ms_);
    }

    // get telescope for defaults
    ROMSColumns msCol(*ms_);
    telescopeName_ = msCol.observation().telescopeName().get(0);
    //getMSTimes(*ms_); // for weather and pwv
}

void PlotMSAtm::setUpCalTable(casacore::String filename,
        PlotMSSelection& userSel) {
    // create NewCalTable, apply user selection, get telescope name,
    caltable_ = new NewCalTable(filename); // original table
    getCalMS(); // ms_ associated with cal table

    // get all times in cal table
    getCalTimes(*caltable_);

    // apply user-selection to caltable_ and ms_
    if (!userSel.isEmpty()) {
        applyCalSelection(userSel, *caltable_);
        if (ms_ != nullptr) {
            applyMSSelection(userSel, *ms_);
        }
    }

    // get telescope for defaults
    ROCTColumns ctCol(*caltable_);
    telescopeName_ = ctCol.observation().telescopeName().get(0);
}

void PlotMSAtm::getMSTimes(MeasurementSet& ms) {
    // unique times in TIME column of ms
    ROMSMainColumns msmc(ms);
    casacore::Vector<casacore::Double> times = msmc.time().getColumn();
    getUniqueTimes(times, mstimes_);
}

void PlotMSAtm::getCalTimes(NewCalTable& ct) {
    // unique times in TIME column of cal table;
    ROCTMainColumns ctmc(ct);
    casacore::Vector<casacore::Double> times = ctmc.time().getColumn();
    getUniqueTimes(times, caltimes_);
}

void PlotMSAtm::getUniqueTimes(casacore::Vector<casacore::Double> inputTimes,
    casacore::Vector<casacore::Double>& uniqueTimes) {
    // find unique times in input times vector
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(inputTimes.data(), TpDouble);
    sorter.sort(indexVector, inputTimes.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    uniqueTimes.resize(nUnique);
    for (uInt i=0; i<nUnique; ++i) {
        uniqueTimes(i) = inputTimes.data()[indexVector(uniqueVector(i))];
    }
} 

void PlotMSAtm::getTimeRange(casacore::Double& mintime, casacore::Double& maxtime) {
    // return cal table time range if it exists, else selected ms time range
    if (!caltimes_.empty()) {
        mintime = min(caltimes_);
        maxtime = max(caltimes_);
    } else {
        mintime = min(mstimes_);
        maxtime = max(mstimes_);
    }
}

void PlotMSAtm::getMSFields(MeasurementSet& ms) {
    // unique fields in FIELD_ID column of ms
    // needed later for elevation calculation for airmass
    ROMSMainColumns ctmc(ms);
    casacore::Vector<casacore::Int> fields = ctmc.fieldId().getColumn();
    getUniqueFields(fields);
}

void PlotMSAtm::getCalFields(NewCalTable& ct) {
    // unique fields in FIELD_ID column of cal table
    // needed later for elevation calculation for airmass
    ROCTMainColumns ctmc(ct);
    casacore::Vector<casacore::Int> fields = ctmc.fieldId().getColumn();
    getUniqueFields(fields);
}

void PlotMSAtm::getUniqueFields(casacore::Vector<casacore::Int> allfields) {
    // find unique fields_
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(allfields.data(), TpInt);
    sorter.sort(indexVector, allfields.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    fields_.resize(nUnique);
    for (uInt i=0; i<nUnique; ++i) {
        fields_(i) = allfields.data()[indexVector(uniqueVector(i))];
    }
}

void PlotMSAtm::getCalMS() {
    // if caltable has keyword for ms name, sets ms_ and tableName_
    String msname("");
    if (caltable_->keywordSet().fieldNumber("MSName") > -1) { 
        msname = caltable_->keywordSet().asString("MSName");
    }
    if (!msname.empty()) {
        casacore::Path path(caltable_->tableName());
        casacore::String fullpath = path.dirName() + "/" + msname;
        casacore::File msname(fullpath);
        if (msname.exists() && msname.isDirectory()) {
            try {
                ms_ = new MeasurementSet(fullpath);
                tableName_ = ms_->tableName();
            } catch (AipsError& err) {
                parent_->logmesg("load_cache", 
                    "MeasurementSet setup failed.\n" + err.getMesg(), 
                    PlotLogger::MSG_WARN);
                parent_->logmesg("load_cache", "Proceeding without ms",
                    PlotLogger::MSG_WARN);
            }
        }
    }
}

void PlotMSAtm::applyMSSelection(PlotMSSelection& selection,
        casacore::MeasurementSet& selms) {
    Vector<Vector<Slice> > chansel, corrsel;
    selection.apply(*ms_, selms, chansel, corrsel);
}

void PlotMSAtm::applyCalSelection(PlotMSSelection& selection,
        NewCalTable& selct) {
    Vector<Vector<Slice> > chansel, corrsel;
    selection.apply(*caltable_, selct, chansel, corrsel);
}

void PlotMSAtm::calcAtmTskyCurve(casacore::Vector<casacore::Double>& curve,
    casacore::Int spw, casacore::Int scan,
    const casacore::Vector<casacore::Double>& chanFreqs) {
    // calculate atm transmission or tsky for given spw and scan
    unsigned int numChan(chanFreqs.nelements());
    curve.resize(numChan);
    curve.set(0.0);

    // cannot calculate curve for single channel
    if (numChan == 1) {
        return;
    }

    curve = calcOverlayCurve(spw, scan, chanFreqs);
    return;
}

void PlotMSAtm::calcImageCurve(casacore::Vector<casacore::Double>& curve,
    casacore::Int spw, casacore::Int scan,
    const casacore::Vector<casacore::Double>& chanFreqs) {
    // calculate image sideband transmission (atm or tsky) for given spw and scan
    unsigned int numChan(chanFreqs.nelements());
    curve.resize(numChan);
    curve.set(0.0);

    // cannot calculate curve for single channel
    if (numChan == 1) {
        return;
    }

    // must have ASDM_RECEIVER table and not split MS
    if (!canShowImageCurve()) {
        curve.set(doubleNaN()); // NaN will not be plotted
        return;
    }

    casacore::Vector<casacore::Double> imageFreqs;
    bool imageFreqsReversed = calcImageFrequencies(imageFreqs, spw, chanFreqs);

    if (allTrue(isNaN(imageFreqs))) { // get LO1 failed
        curve.set(doubleNaN()); // NaN will not be plotted
        return;
    }

    // calculate curve with image frequencies
    curve = calcOverlayCurve(spw, scan, imageFreqs);
    if (imageFreqsReversed) {
        casacore::reverseArray(curve, 0);
    }
}

casacore::Vector<casacore::Double> PlotMSAtm::calcOverlayCurve(
        casacore::Int spw, casacore::Int scan,
        const casacore::Vector<casacore::Double>& chanFreqs) {
    // Implements algorithm in CAS-9053 to get overlay curves
    // (atm, tsky, image sideband)

    PlotMSSelection pmsSel;
    pmsSel.setSpw(String::toString(spw));
    pmsSel.setScan(String::toString(scan));
    if (isMS_) {
        // apply selection to create selected ms
        applyMSSelection(pmsSel, *selms_);
        getMSTimes(*selms_);   // update times_ for airmass calc
        getMSFields(*selms_);  // update fields_ for airmass calc
    } else {
        // apply selection to create selected ct
        applyCalSelection(pmsSel, *selct_);
        getCalTimes(*selct_);  // update times_ for airmass calc
        getCalFields(*selct_); // update fields_ for airmass calc
        if (ms_ != nullptr) {
            applyMSSelection(pmsSel, *selms_);
            getMSTimes(*selms_);   // update times_ for airmass calc
        }
    }

    if ((spw != selectedSpw_) || (scan != selectedScan_)) {
        // save new selection and recalculate
        selectedSpw_ = spw;
        selectedScan_ = scan;
        getMedianPwv();
        getMeanWeather();
    }

    unsigned int numChan(chanFreqs.nelements());
    casacore::Vector<casacore::Double> curve(numChan);
    curve.set(doubleNaN());

    unsigned int midChan(numChan / 2);
    // Set the reference freq to be the middle of the middle two channels
    casacore::Double refFreq = 0.5 * (chanFreqs(midChan-1) + chanFreqs(midChan));
    // reduce number of channels to shorten calculation time
    unsigned int numCalcChan(numChan);
    if (numChan > MAX_ATM_CALC_CHAN_) { 
        while (numCalcChan > MAX_ATM_CALC_CHAN_)
            numCalcChan /= 2;
    }
    unsigned int refChan((numCalcChan - 1) / 2);
    casacore::Double chanSep = (chanFreqs(numChan-1) - chanFreqs(0)) / (numCalcChan - 1);
    if (!isfinite(refFreq) || !isfinite(chanSep)) { // if ms has nan frequencies
        return curve;
    }

    if (numCalcChan % 2 == 0) {
        // if even number of channels, the band center becomes channel boundary
        refFreq -= chanSep * 0.5;
    }

    // set atm parameters
    atm::SpectralGrid* specGrid = new atm::SpectralGrid(numCalcChan, refChan,
        atm::Frequency(refFreq, "GHz"), atm::Frequency(chanSep, "GHz"));
    atm::AtmProfile* atmProfile = getAtmProfile();
    atm::RefractiveIndexProfile* refIdxProfile =
        new atm::RefractiveIndexProfile(*specGrid, *atmProfile);
    atm::SkyStatus* skyStatus = new atm::SkyStatus(*refIdxProfile);
    skyStatus->setUserWH2O(atm::Length(pwv_, "mm"));

    // calculate opacities and airmass
    casacore::Vector<casacore::Double> dryOpacity(numCalcChan);
    casacore::Vector<casacore::Double> wetOpacity(numCalcChan);
    casacore::Vector<casacore::Double> atmTransmission, TebbSky;
    for (uInt chan=0; chan<numCalcChan; ++chan) {
        dryOpacity(chan) = refIdxProfile->getDryOpacity(0,chan).get("neper");
        wetOpacity(chan) = skyStatus->getWetOpacity(0, chan).get("mm-1");
    }
    airmass_ = computeMeanAirmass();

    // calculate atm/tsky
    atmTransmission = exp(-airmass_ * (wetOpacity + dryOpacity));
    if (!showatm_) {
        TebbSky.resize(numCalcChan);
        for (uInt chan=0; chan<numCalcChan; ++chan) {
            TebbSky(chan) = skyStatus->getTebbSky(0, chan).get("K");
        }
    }

    // final calculations for atm/tsky
    casacore::Vector<casacore::Double> calcCurve(numCalcChan);
    if (showatm_) {
        calcCurve = atmTransmission * 100.0; // percent
    } else {
        calcCurve = TebbSky *
            (1.0-atmTransmission) / (1.0-exp((-1.0*wetOpacity)-dryOpacity));
    }

    // clean up
    delete specGrid;
    delete atmProfile;
    delete refIdxProfile;
    delete skyStatus;

    if (numCalcChan == numChan) {
        curve = calcCurve;
    } else {  // fill in curve with calcCurve values, set rest to NaN
        curve.set(casacore::doubleNaN());
        unsigned int inc(numChan / numCalcChan);
        for (unsigned int i=0; i<numCalcChan; ++i) {
            curve(i*inc) = calcCurve(i);
        }
    }
    return curve;
}

atm::AtmProfile* PlotMSAtm::getAtmProfile() {
    // AtmProfile uses weather info and altitude
    casacore::Double altitude(5059.0);
    // Use at tool defaults for other values:
    // Tropospheric lapse rate, scale height, pressure step,
    // pressure step factor, top altitude
    casacore::Double TLR(-5.6), h0(2.0), dp(10.0), dPm(1.2), maxAlt(48.0);
    unsigned int atmtype(atm::midlatWinter);
    atm::AtmProfile* atmProfile = new atm::AtmProfile(
        atm::Length(altitude,"m"), 
        atm::Pressure(weather_.asDouble("pressure"), "mb"), 
        atm::Temperature(weather_.asDouble("temperature"), "K"), TLR, 
        atm::Humidity(weather_.asDouble("humidity"),"%"), atm::Length(h0,"km"), 
        atm::Pressure(dp, "mb"), dPm, atm::Length(maxAlt, "km"), atmtype);
    return atmProfile;
}

void PlotMSAtm::getMedianPwv() {
    // Get pwv (precipitable water vapor) from MS subtable in mm
    casacore::Double pwv(0.0);
    // values from ASDM_CAL subtables if they exist
    if (canCalculatePwv_) {
        casacore::String subtableName;
        casacore::Table mstab(tableName_), subtable;
        casacore::Vector<casacore::Double> waterCol, timesCol;
        try {
          if (mstab.keywordSet().fieldNumber("ASDM_CALWVR") > -1) {
            subtableName = tableName_ + "::ASDM_CALWVR";
            subtable = Table::openTable(subtableName);
            if (subtable.nrow() > 0) {
              waterCol = ScalarColumn<casacore::Double>(subtable, "water").getColumn();
              timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            }
          }
          if (waterCol.empty() && mstab.keywordSet().fieldNumber("ASDM_CALATMOSPHERE") > -1) {
              subtableName = tableName_ + "::ASDM_CALATMOSPHERE";
              subtable = Table::openTable(subtableName);
              if (subtable.nrow() > 0) {
                Array<Double> waterColArray = ArrayColumn<casacore::Double>(subtable, "water").getColumn();
                waterCol = waterColArray(Slicer(Slice(0), Slice()));
                timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
              }
          }
          mstab.closeSubTables();

          if (!waterCol.empty()) {
              casacore::Double mintime(0.0), maxtime(0.0);
              getTimeRange(mintime, maxtime);  // min/max times for selected ms/ct
              // find values in time range
              casacore::Vector<casacore::Double> water =
                  getValuesInTimeRange(waterCol, timesCol, mintime, maxtime);
              if (!water.empty()) {
                  pwv = median(water) * 1000.0; // in mm
              }
          }
        } catch (AipsError & err) {
            // openTable failed, use default pwv
        }
    }

    if (pwv == 0.0) {
        // use default value in mm
        if (telescopeName_=="ALMA") {
            pwv = 1.0;
        } else {
            pwv = 5.0;
        }
    }
    pwv_ = pwv;
}

void PlotMSAtm::getMeanWeather() {
    // Fill weather_ Record with mean values from MS WEATHER table per scan times
    // Set defaults; CAS-9053 default pressure depends on telescope
    casacore::Float humidity(20.0), temperature(273.15);
    casacore::Float pressure = (telescopeName_=="ALMA" ? 563.0 : 786.0);  // mb

    // Use values from WEATHER table if it exists
    if (canCalculateWeather_) {
        try {
            casacore::Table wtable = Table::openTable(tableName_ + "::WEATHER");
            TableColumn tempCol = TableColumn(wtable, "TEMPERATURE");
            String tempUnits = tempCol.keywordSet().asArrayString("QuantumUnits").tovector()[0];
            TableColumn pressCol = TableColumn(wtable, "PRESSURE");
            String pressUnits = pressCol.keywordSet().asArrayString("QuantumUnits").tovector()[0];
            // select valid stations and values in range
            casacore::Table selwtable = selectWeatherTable(wtable, tempUnits, pressUnits);
            
            // get columns from *selected* weather table
            casacore::Vector<casacore::Float> 
                pressureCol(ScalarColumn<casacore::Float>(selwtable, "PRESSURE").getColumn()),
                humidityCol(ScalarColumn<casacore::Float>(selwtable, "REL_HUMIDITY").getColumn()),
                temperatureCol(ScalarColumn<casacore::Float>(selwtable, "TEMPERATURE").getColumn());
            casacore::Vector<casacore::Double>
                timeCol(ScalarColumn<casacore::Double>(selwtable, "TIME").getColumn());

            // select values based on cal times
            casacore::Vector<casacore::Float> selpressure, selhumidity, seltemperature;
            casacore::Float meanP(0.0), meanH(0.0), meanT(0.0);
            casacore::Double mintime(min(mstimes_)), maxtime(max(mstimes_));
            // pressure
            if (!pressureCol.empty()) {
                selpressure = getValuesInTimeRange(pressureCol, timeCol, mintime, maxtime);
                if (!selpressure.empty()) {
                    meanP = mean(selpressure);
                }
                if ((pressUnits == "Pa") && (meanP > 1013.25)) {
                    meanP /= (float)100.0;  // Pa to hPa
                }
                if (meanP != 0.0) {
                    pressure = meanP;
                }
            }
            // humidity
            if (!humidityCol.empty()) {
                selhumidity = getValuesInTimeRange(humidityCol, timeCol, mintime, maxtime);
                if (!selhumidity.empty()) {
                    meanH = mean(selhumidity);
                }
                if (meanH != 0.0) {
                    humidity = meanH;
                }
            }

            // temperature
            if (!temperatureCol.empty()) {
                seltemperature = getValuesInTimeRange(temperatureCol, timeCol, mintime, maxtime);
                if (!seltemperature.empty()) {
                    meanT = mean(seltemperature);
                }
                if ((meanT > 0.0) && (tempUnits == "C")) {
                    meanT += (float)273.15;  // convert C to K
                }
                if (meanT != 0.0) {
                    temperature = meanT;
                }
            }
        } catch (AipsError & err) {
            std::cout << "Failed to read WEATHER table:" << err.getMesg() << std::endl;
        }
    }

    // to use in atm::AtmProfile 
    weather_.define("humidity", humidity);       // %
    weather_.define("pressure", pressure);       // mb
    weather_.define("temperature", temperature); // K
}

Table PlotMSAtm::selectWeatherTable(Table& wtable, String tempUnits, String pressureUnits) {
    // Remove invalid station, values from weather table
    TableExprNode ten;
    if (telescopeName_=="ALMA") {
        try {
            // do not use weather station "MeteoItinerant"
            casacore::Table staTable = Table::openTable(tableName_ + "::ASDM_STATION");
            casacore::Table result = staTable(staTable.col("name")=="MeteoItinerant");
            Vector<uInt> antIds = result.rowNumbers(staTable);
            if (!antIds.empty()) ten = (wtable.col("NS_WX_STATION_ID") != antIds(0));
        } catch (AipsError & err) {
            // table does not exist 
        }
    }
    // SCIREQ-1241 ranges
    // temperature in the range -30 to +30C
    Float mintemp = (tempUnits=="C" ? -30.0 : 243.15);
    Float maxtemp = (tempUnits=="C" ? 30.0 : 303.15);
    ten = ten && (wtable.col("TEMPERATURE")>=mintemp && wtable.col("TEMPERATURE")<=maxtemp);
    // pressure in the range 450 to 800mb
    // CAS-11528 units are "Pa" but might be hPa
    // need to check, for legacy datasets
    if (pressureUnits=="Pa") {
        TableColumn tabcol(wtable, "PRESSURE");
        for (uInt irow=0; irow<tabcol.nrow(); ++irow) {
            float val(tabcol.asfloat(irow));
            if (val > 0.0 && val < 1013.25) { // avg sea-level pressure mb
                pressureUnits = "hPa";
                break;
            }
        }
    }
    Float minpress = (pressureUnits=="hPa" ? 450.0 : 45000.0);
    Float maxpress = (pressureUnits=="hPa" ? 900.0 : 90000.0);
    ten = ten && (wtable.col("PRESSURE")>=minpress && wtable.col("PRESSURE")<=maxpress);
    // Humidity in the range -10 to +110%
    ten = ten && (wtable.col("REL_HUMIDITY")>=-10.0 && wtable.col("REL_HUMIDITY")<=110.0);
    casacore::Table selwtable(wtable(ten));
    return selwtable;
}

casacore::Double PlotMSAtm::computeMeanAirmass() {
    // Calculate airmass from elevation
    casacore::Double airmass(0.0);
    casacore::Double elevation = getPointingElevation();
    if (elevation == 0.0) { // could not get from POINTING table
        casacore::Vector<casacore::Double> airmasses(fields_.size());
        for (uInt i=0; i<fields_.size(); ++i) {
            elevation = getFieldElevation(fields_(i));
            if (elevation <= 3.0) {
                elevation = 45.0;
            }
            airmasses(i) = 1.0 / std::cos((90.0 - elevation) * C::pi / 180.0);
        }
        airmass = mean(airmasses);
    } else {
        airmass = 1.0 / std::cos((90.0 - elevation) * C::pi / 180.0);
    }
    return airmass;
}

casacore::Double PlotMSAtm::getPointingElevation() {
    // Determine mean elevation from POINTING table based on selected ms times
    casacore::Double elevation(0.0);
    ROMSColumns* msCol;
    if (selms_ != nullptr) {
        msCol = new ROMSColumns(*selms_);
    } else if (ms_ != nullptr) {
        msCol = new ROMSColumns(*ms_);
    } else {
        return elevation;
    }

    if (msCol != nullptr) {
        Array<Double> pointingAzEl = msCol->pointing().direction().getColumn();
        Vector<Double> pointingTime = msCol->pointing().time().getColumn();
        Vector<Int> pointingAntennaId = msCol->pointing().antennaId().getColumn();
        delete msCol;

        // find elevation measurements for the selected ms times and antenna 0
        casacore::Double mintime(min(mstimes_)), maxtime(max(mstimes_));
        Vector<Double> elevations;
        for (size_t i=0; i < pointingTime.size(); ++i) {
            if ((pointingAntennaId(i) == 0) && (pointingTime(i) >= mintime) &&
                (pointingTime(i) <= maxtime)) {
                Double dirElevationRad = pointingAzEl(IPosition(3, 1, 0, i));
                size_t elSize = elevations.size();
                elevations.resize(elSize + 1, true);
                elevations(elSize) = dirElevationRad * 180.0 / C::pi;
            }
        }
        if (!elevations.empty()) {
            elevation = mean(elevations);
        }
    }
    return elevation;
}

casacore::Double PlotMSAtm::getFieldElevation(casacore::Int fieldId) {
    // Determine elevation from FIELD table for selected fieldId
    casacore::Array<casacore::Double> raDec;
    if (isMS_) {
        ROMSColumns msCol(*ms_);
        raDec = msCol.field().delayDir().get(fieldId);
    } else {
        ROCTColumns ctCol(*caltable_);
        raDec = ctCol.field().delayDir().get(fieldId);
    }
    casacore::Double ra = raDec(IPosition(2,0,0));
    casacore::Double dec = raDec(IPosition(2,1,0));
    // convert J2000 to AzEl
    // need observatory position and timestamp for output reference
    casacore::MPosition obspos;
    casacore::MeasTable::Observatory(obspos, telescopeName_);
    // caltimes selected for scan
    casacore::Double timestamp = getMeanScantime();
    casacore::MEpoch ts(casacore::Quantity(timestamp/86400.0, "d"));
    casacore::MeasFrame frame(obspos, ts);
    casacore::MDirection::Ref inputRef(casacore::MDirection::J2000);
    casacore::MDirection inputDir(casacore::Quantity(ra, "rad"), 
    casacore::Quantity(dec, "rad"), inputRef);
    casacore::MDirection::Ref outputRef(casacore::MDirection::AZEL, frame);
    // do conversion
    casacore::MDirection::Convert j2toAzel(inputDir, outputRef);
    casacore::Vector<casacore::Double> azel = j2toAzel().getAngle("deg").getValue();
    casacore::Double elevation = azel(1);
    return elevation;
}

casacore::Double PlotMSAtm::getMeanScantime() {
    // get mean timestamp from selected table
    casacore::Vector<casacore::Double> times;
    if (isMS_) {
        ROMSMainColumns msmc(*selms_);
        times = msmc.time().getColumn();
    } else {
        ROCTMainColumns ctmc(*selct_);
        times = ctmc.time().getColumn();
    }
    casacore::Double meantime(mean(times));
    return meantime;
}

template <typename T>
casacore::Vector<T> PlotMSAtm::getValuesInTimeRange(casacore::Vector<T> inputData, 
        casacore::Vector<casacore::Double> times, casacore::Double mintime,
        casacore::Double maxtime) {
    // Use values with timestamps in times_ range
    casacore::Vector<T> outputData;
    casacore::uInt outsize(0);
    for (size_t j = 0; j < times.size(); ++j) {
        if ((times(j) >= mintime) && (times(j) <= maxtime)) {
            outputData.resize(outsize+1, True);
            outputData(outsize++) = inputData(j);
        }
    }
    if (outputData.empty()) {
        getClosestValues(outputData, times, inputData, mintime, maxtime);
    }
    return outputData;
}

template <typename T>
void PlotMSAtm::getClosestValues(casacore::Vector<T>& values,
    casacore::Vector<casacore::Double>& times, casacore::Vector<T>& data,
    double mintime, double maxtime) { 
    // return value(s) with the closest time outside range mintime~maxtime
    casacore::Double mintimediff(1.0e12);
    for (size_t i = 0; i < times.size(); ++i) {
        casacore::Double timediff = min(abs(times(i)-mintime), abs(times(i)-maxtime));
        if (timediff < mintimediff) {
            values.resize(1);
            values(0) = data(i); // save value closest to timerange
            mintimediff = timediff;
        } else if (timediff == mintimediff) {
            size_t valuesSize = values.size();
            values.resize(valuesSize + 1, true);
            values(valuesSize) = data(i);
        }
    }
}

bool PlotMSAtm::hasReceiverTable() {
    // axis can be calculated only if ms has ASDM_RECEIVER table
    bool hasTable(false);
    if (tableName_.empty()) {
        return hasTable;
    }
    // check if ms has receiver subtable
    casacore::Table mstab(tableName_);
    hasTable = (mstab.keywordSet().fieldNumber("ASDM_RECEIVER") > -1);
    return hasTable;
}

bool PlotMSAtm::canGetLOsForSpw() {
    // spw IDs must be consistent between ASDM_RECEIVER table and SPECTRAL_WINDOW table
    // (if ms was split, spws were reindexed)
    if (tableName_.empty()) {
        return false;
    }
    // spws from ASDM_RECEIVER table
    casacore::Table receiverTable = Table::openTable(tableName_ + "::ASDM_RECEIVER");
    casacore::Vector<casacore::String> receiverSpws =
        ScalarColumn<casacore::String>(receiverTable, "spectralWindowId").getColumn();
    int maxReceiverSpw(0);
    for (auto spwId : receiverSpws) {
        casacore::String spwIdNumber = spwId.after(spwId.rfind("_")); // e.g. "SpectralWindow_11" -> "11"
        try {
            int spw = std::stoi(spwIdNumber);
            maxReceiverSpw = std::max(maxReceiverSpw, spw);
        } catch (std::invalid_argument& error) {
            // ignore: cannot convert spwIdNumber to int 
        }
    }
    // spws from SPECTRAL_WINDOW table
    casacore::Table spwTable = Table::openTable(tableName_ + "::SPECTRAL_WINDOW");
    int nSpwSpws = ScalarColumn<casacore::String>(spwTable, "NAME").getColumn().size();

    return (maxReceiverSpw == nSpwSpws-1);
}

bool PlotMSAtm::calcImageFrequencies(
      casacore::Vector<casacore::Double>& imageFreqs,
      casacore::Int spw, const casacore::Vector<casacore::Double>& chanFreqs) {
    // Get LO1 for spw and calculate image frequencies for signal frequencies;
    // returns whether frequencies are reversed in image sideband
    unsigned int numChan(chanFreqs.nelements());
    imageFreqs.resize(numChan);

    double freqLO1(0.0);
    if (!getLO1FreqForSpw(freqLO1, spw)) { // spw not found
        imageFreqs.set(doubleNaN());
        return false;
    }

    // use image frequencies based on LO1 frequency
    casacore::Vector<casacore::Double> lo1Freqs(numChan,  (2.0 * freqLO1));
    imageFreqs = lo1Freqs - chanFreqs;

    // calculate number of chans to use for calculation, reference channel,
    // reference frequency, and channel separation for atm::SpectralGrid
    unsigned int midChan(numChan/2);
    // Set the reference freq to be the middle of the middle two channels
    casacore::Double refFreq = 0.5 * (imageFreqs(IPosition(2, midChan-1, 0)) +
        imageFreqs(IPosition(2, midChan, 0)));
    unsigned int refChan((numChan - 1) / 2);
    casacore::Double chanSep = (imageFreqs(IPosition(2, numChan-1, 0)) -
        imageFreqs(IPosition(2, 0, 0))) / (numChan - 1);
    if (!isfinite(refFreq) || !isfinite(chanSep)) {
        imageFreqs.set(doubleNaN());
        return false;
    }
        
    if (numChan % 2 == 0) {
        refFreq -= chanSep*0.5;
    }

    // Set up SpectralGrid and get frequencies
    atm::SpectralGrid* specGrid = new atm::SpectralGrid(numChan, refChan,
        atm::Frequency(refFreq, "GHz"), atm::Frequency(chanSep, "GHz"));
    for (unsigned int chan=0; chan<numChan; ++chan) {
        imageFreqs[chan] = specGrid->getChanFreq(0, static_cast<unsigned int>(chan)).get("GHz");
    }

    casacore::Vector<casacore::Double> adjustedImageFreqs = lo1Freqs - imageFreqs;
    return !(allNearAbs(adjustedImageFreqs, chanFreqs, 1e-7));
}

bool PlotMSAtm::getLO1FreqForSpw(double& freqGHz, int spw) {
    // For showimage, get the freqLO[0] column value from ASDM_RECEIVER table.
    // Assigns freqGHz argument, returns whether spw was found
    bool foundSpw(false);
    if (!hasReceiverTable()) {
        return foundSpw;
    }
    if (loFreqForSpw_.count(spw)) {
        freqGHz = loFreqForSpw_[spw];
        foundSpw = true;
    } else {
        // get names and freqLOs from ASDM_RECEIVER table, only keep first "WVR" freq
        casacore::Table receiverTable = Table::openTable(tableName_ + "::ASDM_RECEIVER");
        casacore::Vector<casacore::String> spwNames =
            ScalarColumn<casacore::String>(receiverTable, "spectralWindowId").getColumn();
        casacore::Vector<casacore::String> receiverNames =
            ScalarColumn<casacore::String>(receiverTable, "name").getColumn();
        ArrayColumn<casacore::Double> freqLOCol(receiverTable, "freqLO");
        casacore::Vector<casacore::Double> freqLOs;
        casacore::Vector<casacore::String> spwIds;
        bool gotWvr(false);
        for (size_t i = 0; i < receiverNames.size(); ++i) {
            // only use first value for spw id
            casacore::String spwIdStr = spwNames(i).after("SpectralWindow_");
            bool foundSpwId(false);
            for (size_t j=0; j < spwIds.size(); ++j) {
                if (spwIds(j) == spwIdStr) {
                    foundSpwId = True;
                    break;
                }
            }
            if (foundSpwId) {
               continue;
            }

            size_t idsSize = spwIds.size();
            spwIds.resize(idsSize + 1, true);
            spwIds(idsSize) = spwIdStr;

            // only use first value for WVR
            if (receiverNames(i).contains("WVR")) {
                if (gotWvr) {
                    continue;
                } else {
                    gotWvr = true;
                }
            }

            // get freqLO[0] and convert to GHz
            double freqLO = *(freqLOCol.get(i).data()); // first element of array
            freqLO *= 1e-9;
            // add this freqLO value
            size_t loSize(freqLOs.size());
            freqLOs.resize(loSize + 1, true);
            freqLOs(loSize) = freqLO;
        }

        int freqsize = freqLOs.size();
        if (freqsize > spw) {
            freqGHz = freqLOs(spw);
            // only log this when calculated for this spw
            parent_->logmesg("load_cache",
                "Image sideband LO1=" + String::toString(freqGHz) + " for spw " + String::toString(spw));
            // cache this for later chunks
            loFreqForSpw_[spw] = freqGHz;
            foundSpw = true;
        }
    }
    return foundSpw;
}

}
