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
#include <casa/Arrays/ArrayMath.h>
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
        bool isMS, PlotMSCacheBase* parent):
    isMS_(isMS),
    parent_(parent),
    ms_(nullptr),
    caltable_(nullptr),
    tableName_(""),
    telescopeName_(""),
    pwv_(0.0),
    airmass_(0.0),
    MAX_ATM_CALC_CHAN_(512) {

    // set up table and needed data (times, fields)
    if (isMS) {
        setUpMS(filename, userSel);
    } else {
        setUpCalTable(filename, userSel);
    }

    getMeanWeather();
    getMedianPwv();
    selms_ = new MeasurementSet();
    selct_ = new NewCalTable();
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
    try {
        ms_ = new MeasurementSet(filename);
    } catch (AipsError& err) {
        throw(AipsError("MeasurementSet setup failed.\n" + err.getMesg()));
        return;
    }
    tableName_ = ms_->tableName();
    if (!userSel.isEmpty()) 
        applyMSSelection(userSel, *ms_);
    ROMSColumns msCol(*ms_);
    telescopeName_ = msCol.observation().telescopeName().get(0);
    getMSTimes(); // for weather and pwv
}

void PlotMSAtm::setUpCalTable(casacore::String filename,
        PlotMSSelection& userSel) {
    caltable_ = new NewCalTable(filename); // original table
    getCalMS(); // MeasurementSet* ms_ associated with cal table
    if (!userSel.isEmpty()) {
        applyCalSelection(userSel, *caltable_); // now user-selected table
    }
    ROCTColumns ctCol(*caltable_);
    telescopeName_ = ctCol.observation().telescopeName().get(0);
    getCalTimes();
}

void PlotMSAtm::getMSTimes() {
    // unique times in TIME column of ms
    // ms could be original or selected ms
    ROMSMainColumns msmc(*ms_);
    casacore::Vector<casacore::Double> times = msmc.time().getColumn();
    getUniqueTimes(times);
}

void PlotMSAtm::getCalTimes() {
    // unique times in TIME column of cal table;
    // nct could be original or selected cal table
    ROCTMainColumns ctmc(*caltable_);
    casacore::Vector<casacore::Double> times = ctmc.time().getColumn();
    getUniqueTimes(times);
}

void PlotMSAtm::getUniqueTimes(casacore::Vector<casacore::Double> alltimes) {
    // find unique times_
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(alltimes.data(), TpDouble);
    sorter.sort(indexVector, alltimes.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    times_.resize(nUnique);
    for (uInt i=0; i<nUnique; ++i) {
        times_(i) = alltimes.data()[indexVector(uniqueVector(i))];
    }
} 

void PlotMSAtm::getMSFields() {
    // unique fields in FIELD_ID column of ms
    // needed later for elevation calculation for airmass
    ROMSMainColumns ctmc(*selms_);
    casacore::Vector<casacore::Int> fields = ctmc.fieldId().getColumn();
    getUniqueFields(fields);
}

void PlotMSAtm::getCalFields() {
    // unique fields in FIELD_ID column of cal table
    // needed later for elevation calculation for airmass
    ROCTMainColumns ctmc(*selct_);
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
    // if caltable has keyword for ms name, sets ms_
    String msname("");
    if (caltable_->keywordSet().fieldNumber("MSName") > -1) 
        msname = caltable_->keywordSet().asString("MSName");
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

casacore::Vector<casacore::Double> PlotMSAtm::calcOverlayCurve(
        casacore::Int spw, casacore::Int scan,
        const casacore::Vector<casacore::Double>& chanFreqs,
        bool atm) {
    // Implements algorithm in CAS-9053 to get overlay curves
    // (atm or tsky) per spw + scan
    unsigned int numChan(chanFreqs.nelements());
    casacore::Vector<casacore::Double> curve(numChan, 0.0);
    if (numChan == 1)
        return curve;

    PlotMSSelection pmsSel;
    pmsSel.setSpw(String::toString(spw));
    pmsSel.setScan(String::toString(scan));
    if (isMS_) {
        // apply selection to create selected ms
        applyMSSelection(pmsSel, *selms_);
        getMSFields();  // update fields for airmass calc
    } else {
        // apply selection to create selected ct
        applyCalSelection(pmsSel, *selct_);
        getCalFields();  // update fields for airmass calc
    }

    unsigned int numCalcChan(numChan);
    unsigned int midChan(numCalcChan/2);
    // Set the reference freq to be the middle of the middle two channels
    casacore::Double refFreq = 0.5 * (chanFreqs(IPosition(2, midChan-1, 0))
        + chanFreqs(IPosition(2, midChan, 0)));
    // reduce number of channels to shorten calc time
    if (numChan > MAX_ATM_CALC_CHAN_) { 
        while (numCalcChan > MAX_ATM_CALC_CHAN_)
            numCalcChan /= 2;
        midChan = numCalcChan/2;
    }
    unsigned int refChan((numCalcChan - 1) / 2);
    casacore::Double chanSep = (chanFreqs(IPosition(2, numChan-1, 0))
        - chanFreqs(IPosition(2, 0, 0))) / (numCalcChan - 1);
    if (numCalcChan % 2 == 0) refFreq -= chanSep*0.5;

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
    atmTransmission = exp(-airmass_ * (wetOpacity + dryOpacity));
    if (!atm) {
        TebbSky.resize(numCalcChan);
        for (uInt chan=0; chan<numCalcChan; ++chan) {
            TebbSky(chan) = skyStatus->getTebbSky(0, chan).get("K");
        }
    }
    // final calculations
    casacore::Vector<casacore::Double> calcCurve(numCalcChan);
    if (atm) {
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
    if (!tableName_.empty()) {
        casacore::String subname;
        casacore::Table mstab(tableName_), subtable;
        casacore::Vector<casacore::Double> waterCol, timesCol;
        try {
          if (mstab.keywordSet().fieldNumber("ASDM_CALWVR") > -1) {
            subname = tableName_ + "::ASDM_CALWVR";
            subtable = Table::openTable(subname);
            if (subtable.nrow() > 0) {
              waterCol = ScalarColumn<casacore::Double>(subtable, "water").getColumn();
              timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            }
          }
          if (waterCol.empty() && mstab.keywordSet().fieldNumber("ASDM_CALATMOSPHERE") > -1) {
            subname = tableName_ + "::ASDM_CALATMOSPHERE";
            subtable = Table::openTable(subname);
            if (subtable.nrow() > 0) {
              Array<Double> waterColArray = ArrayColumn<casacore::Double>(subtable, "water").getColumn();
              waterCol = waterColArray(Slicer(Slice(0), Slice()));
              timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            }
          }
          mstab.closeSubTables();
          if (waterCol.empty()) {
            parent_->logmesg("load_cache",
              "ASDM_CALWVR and ASDM_CALATMOSPHERE tables could not be opened or have zero rows");
          } else {
            casacore::Vector<casacore::Double> water = 
                getValuesNearTimes(waterCol, timesCol);
            if (!water.empty()) pwv = median(water) * 1000.0; // in mm
          }
        } catch (AipsError & err) {
            // openTable failed, use default pwv
        }
    }
    // else use default value in mm
    if (pwv == 0.0) {
        if (telescopeName_=="ALMA") pwv = 1.0;
        else pwv = 5.0;
        parent_->logmesg("load_cache",
            "Using default pwv " + casacore::String::toString(pwv) + " for telescope " + telescopeName_);
    }
    pwv_ = pwv;
}

void PlotMSAtm::getMeanWeather() {
    // Fill weather_ Record with mean values from MS WEATHER table
    // Set defaults; CAS-9053 default pressure depends on telescope
    casacore::Float humidity(20.0), temperature(273.15);
    casacore::Float pressure = (telescopeName_=="ALMA" ? 563.0 : 786.0);  // mb

    // Use values from WEATHER table if it exists
    bool noWeather(true);
    if (!tableName_.empty()) {
        try {
            // openTable throws exception if table doesn't exist
            casacore::Table wtable = Table::openTable(tableName_ + "::WEATHER");
            TableColumn tempCol = TableColumn(wtable, "TEMPERATURE");
            String tempUnits = tempCol.keywordSet().asArrayString("QuantumUnits").tovector()[0];
            TableColumn pressCol = TableColumn(wtable, "PRESSURE");
            String pressUnits = pressCol.keywordSet().asArrayString("QuantumUnits").tovector()[0];
            // select valid stations and values in range
            casacore::Table selwtable = selectWeatherTable(wtable, tempUnits, pressUnits);
            
            // get columns
            casacore::Vector<casacore::Float> 
                pressureCol(ScalarColumn<casacore::Float>(selwtable, "PRESSURE").getColumn()),
                humidityCol(ScalarColumn<casacore::Float>(selwtable, "REL_HUMIDITY").getColumn()),
                temperatureCol(ScalarColumn<casacore::Float>(selwtable, "TEMPERATURE").getColumn());
            casacore::Vector<casacore::Double>
                timeCol(ScalarColumn<casacore::Double>(selwtable, "TIME").getColumn());

            // select values based on cal times
            casacore::Vector<casacore::Float> selpressure, selhumidity, seltemperature;
            casacore::Float meanP(0.0), meanH(0.0), meanT(0.0);
            // pressure
            if (!pressureCol.empty()) {
                selpressure = getValuesNearTimes(pressureCol, timeCol);
                if (!selpressure.empty()) {
                    meanP = mean(selpressure);
                }
                if (meanP == 0.0) {
                    parent_->logmesg("load_cache", "WEATHER pressure is zero, using default value instead.");
                } else {
                    if ((pressUnits == "Pa") && (meanP > 1013.25)) {
                        meanP /= (float)100.0;  // Pa to hPa
                    }
                    pressure = meanP;
                }
            }
            // humidity
            if (!humidityCol.empty()) {
                selhumidity = getValuesNearTimes(humidityCol, timeCol);
                if (!selhumidity.empty()) {
                    meanH = mean(selhumidity);
                }
                if (meanH == 0.0) {
                    parent_->logmesg("load_cache", "WEATHER humidity is zero, using default value instead.");
                } else {
                    humidity = meanH;
                }
            }

            // temperature
            if (!temperatureCol.empty()) {
                seltemperature = getValuesNearTimes(temperatureCol, timeCol);
                if (!seltemperature.empty()) {
                    meanT = mean(seltemperature);
                }
                if (meanT == 0.0) {
                    parent_->logmesg("load_cache", "WEATHER temperature is zero, using default value instead.");
                } else {
                    if (tempUnits=="C") {
                        meanT += (float)273.15;  // convert C to K
                    }
                    temperature = meanT;
                }
            }
            noWeather = false;
        } catch (AipsError & err) {
            cout << err.getMesg() << endl;
        }
    }

    String msg = "Using";
    msg += (noWeather ? " default " : " ");
    msg += "weather values for Atm/Tsky:\n";
    msg += "  humidity " + casacore::String::toString(humidity) + 
        ", pressure " + casacore::String::toString(pressure) +
        ", temperature " + casacore::String::toString(temperature);
    parent_->logmesg("load_cache", msg);
    // to use in atmosphere.initAtmProfile (tool)
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
    casacore::Vector<casacore::Double> airmasses;
    airmasses.resize(fields_.size());
    casacore::Double elevation;
    for (uInt i=0; i<fields_.size(); ++i) {
        elevation = getElevation(fields_(i));
        if (elevation <= 3.0) elevation = 45.0;
        airmasses(i) = 1.0 / std::cos((90.0 - elevation) * C::pi / 180.0);
    }
    casacore::Double airmass = mean(airmasses);
    return airmass;
}

casacore::Double PlotMSAtm::getElevation(casacore::Int fieldId) {
    // Get RADec (DELAY_DIR) from FIELD table
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
    casacore::Double el = azel(IPosition(2,1,0));
    return el;
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
casacore::Vector<T> PlotMSAtm::getValuesNearTimes(
        casacore::Vector<T> inputCol, 
        casacore::Vector<casacore::Double> timesCol) {
    // Use values with timestamps close to times_ from selected table
    // First, find time difference in closest match
    casacore::Double mintimediff(1.0e12);
    for (uInt i=0; i<times_.size(); ++i) {
        for (uInt j=0; j<timesCol.size(); ++j) {
            mintimediff = min(mintimediff, abs(times_(i)-timesCol(j)));
        }
    }

    // now get index of values within 1 sec of minimum timediff
    casacore::Vector<casacore::uInt> values;
    casacore::uInt vsize(0);
    for (casacore::uInt i=0; i<times_.size(); ++i) {
        for (casacore::uInt j=0; j<timesCol.size(); ++j) {
            if ((abs(times_(i)-timesCol(j)) - mintimediff) < 1.0) {
                values.resize(vsize+1, True);
                values(vsize++) = j;
            }
        }
    }
    // unique values
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    sorter.sortKey(values.data(), TpUInt);
    sorter.sort(indexVector, values.size());
    uInt nUnique = sorter.unique(uniqueVector, indexVector);
    casacore::Vector<T> outputCol(nUnique);
    uInt index;
    for (casacore::uInt i=0; i<nUnique; ++i) {
        index = values.data()[indexVector(uniqueVector(i))];
        outputCol(i) = inputCol(index);
    }
    return outputCol;
}

bool PlotMSAtm::hasReceiverTable() {
    // axis can be calculated only if ms has ASDM_RECEIVER table
    if (tableName_.empty()) {
        return false;
    }
    // check if ms has receiver subtable
    casacore::Table mstab(tableName_);
    return mstab.keywordSet().fieldNumber("ASDM_RECEIVER") > -1;
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

casacore::Vector<casacore::Double> PlotMSAtm::calcImageSidebandCurve(
        casacore::Int /*spw*/, casacore::Int /*scan*/,
        const casacore::Vector<casacore::Double>& chanFreqs) {
    // must have >1 chan and MeasurementSet to calculate sideband curve
    unsigned int numChan(chanFreqs.nelements());
    casacore::Vector<casacore::Double> curve(numChan, 0.0);
    if (numChan == 1) {
        return curve;
    }
    if (!canShowImageCurve()) {
        curve.set(doubleNaN()); // will not be plotted
        return curve;
    }

    curve.set(80.0); // TODO: for testing
    return curve;
}

}
