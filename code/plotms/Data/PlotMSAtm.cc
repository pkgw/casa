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
    ms_(NULL),
    caltable_(NULL),
    tableName_(""),
    telescopeName_(""),
    pwv_(0.0),
    airmass_(0.0),
    MAX_ATM_CALC_CHAN_(512) {

    // set up table and needed data (times, fields)
    if (isMS) 
        setUpMS(filename, userSel);
    else
        setUpCalTable(filename, userSel);

    getMeanWeather();
    getMedianPwv();
    selms_ = new MeasurementSet();
    selct_ = new NewCalTable();
}

PlotMSAtm::~PlotMSAtm() {
    if (caltable_) {
        delete caltable_;
        caltable_ = NULL;
    }
    if (ms_) {
        delete ms_;
        ms_ = NULL;
    }
    if (selms_) {
        delete selms_;
        selms_ = NULL;
    }
    if (selct_) {
        delete selct_;
        selct_ = NULL;
    }
}

void PlotMSAtm::setUpMS(casacore::String filename, PlotMSSelection& userSel) {
    try {
        ms_ = new MeasurementSet(filename);
    } catch (AipsError& err) {
        throw(AipsError("MeasurementSet setup failed.\n" + err.getMesg()));
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
    getCalMS();
    if (!userSel.isEmpty())
        applyCalSelection(userSel, *caltable_); // now user-selected table
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
    if (numChan==1)
        return curve;

    PlotMSSelection pmsSel;
    pmsSel.setSpw(String::toString(spw));
    pmsSel.setScan(String::toString(scan));
    if (isMS_) {
        applyMSSelection(pmsSel, *selms_);
        getMSFields();  // update fields for airmass calc
    } else {
        applyCalSelection(pmsSel, *selct_);
        getCalFields();  // update fields for airmass calc
    }

    unsigned int midChan(numChan/2), numCalcChan(numChan);
    unsigned int refChan((numChan - 1) / 2);
    if (numChan > MAX_ATM_CALC_CHAN_) { 
        while (numCalcChan > MAX_ATM_CALC_CHAN_)
            numCalcChan /= 2;
    }
    casacore::Double refFreq = 0.5 * (chanFreqs(IPosition(2, midChan-1, 0))
        + chanFreqs(IPosition(2, midChan, 0)));
    casacore::Double chanSep = (chanFreqs(IPosition(2, numChan-1, 0))
        - chanFreqs(IPosition(2, 0, 0))) / (numChan - 1);
    if (numChan % 2 == 0) refFreq -= chanSep*0.5;
    // set atm parameters
    atm::SpectralGrid* specGrid = new atm::SpectralGrid(numCalcChan, refChan,
        atm::Frequency(refFreq, "GHz"), atm::Frequency(chanSep, "GHz"));
    atm::AtmProfile* atmProfile = getAtmProfile();
    atm::RefractiveIndexProfile* refIdxProfile =
        new atm::RefractiveIndexProfile(*specGrid, *atmProfile);
    atm::SkyStatus* skyStatus = new atm::SkyStatus(*refIdxProfile);
    skyStatus->setUserWH2O(atm::Length(pwv_, "mm"));
    // opacity vectors may have fewer elements (numCalcChan) than numChan
    // to save time
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
    // final calculations; calcCurve may be smaller than curve
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
            waterCol = ScalarColumn<casacore::Double>(subtable, "water").getColumn();
            timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            mstab.closeSubTables();
          } else if (mstab.keywordSet().fieldNumber("ASDM_CALATMOSPHERE") > -1) {
            subname = tableName_ + "::ASDM_CALATMOSPHERE";
            subtable = Table::openTable(subname);
            waterCol = ArrayColumn<casacore::Double>(subtable, "water").getColumn();
            timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            mstab.closeSubTables();
          }
          if (!waterCol.empty()) {
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
            "Could not open ASMD_CALWVR or ASDM_CALATMOSPHERE table; using default pwv " + casacore::String::toString(pwv) + " for telescope " + telescopeName_);
    }
    pwv_ = pwv;
}

void PlotMSAtm::getMeanWeather() {
    // Info from MS WEATHER table in weather_ Record
    // set defaults
    casacore::Float pressure, humidity(20.0), temperature(273.15);
    // NB: plotbandpass uses default pressure 563 in all cases;
    // see CAS-9053 algorithm #2
    if (telescopeName_=="ALMA") pressure = 563.0;  // mb
    else pressure = 786.0;  // mb

    // values from WEATHER table if it exists
    bool noWeather(false);
    if (!tableName_.empty()) {
        casacore::Table mstab(tableName_);
        if (mstab.keywordSet().fieldNumber("WEATHER") > -1) {
            casacore::String subname = tableName_ + "::WEATHER";
            try {
                casacore::Table subtable = Table::openTable(subname);
                // get columns and temp. units
                casacore::Vector<casacore::Float> pressureCol, humidityCol, 
                    temperatureCol;
                casacore::Vector<casacore::Double> timeCol;
                pressureCol = ScalarColumn<casacore::Float>(subtable, "PRESSURE").getColumn();
                humidityCol = ScalarColumn<casacore::Float>(subtable, "REL_HUMIDITY").getColumn();
                ScalarColumn<casacore::Float> tempColumn = 
                    ScalarColumn<casacore::Float>(subtable, "TEMPERATURE");
                temperatureCol = tempColumn.getColumn();
                casacore::String units = 
                    tempColumn.keywordSet().asArrayString("QuantumUnits").tovector()[0];
                timeCol = ScalarColumn<casacore::Double>(subtable, "TIME").getColumn();
                mstab.closeSubTables();

                // now select based on cal times
                casacore::Vector<casacore::Float> selpressure, selhumidity,
                    seltemperature;
                // pressure
                if (!pressureCol.empty()) {
                    selpressure = getValuesNearTimes(pressureCol, timeCol);
                    if (!selpressure.empty()) pressure = mean(selpressure);
                }
                // humidity
                if (!humidityCol.empty()) {
                    selhumidity = getValuesNearTimes(humidityCol, timeCol);
                    if (!selhumidity.empty()) humidity = mean(selhumidity);
                }
                // temperature
                if (!temperatureCol.empty()) {
                    seltemperature = getValuesNearTimes(temperatureCol, timeCol);
                    if (!seltemperature.empty()) {
                        temperature = mean(seltemperature);
                        if (units.compare("K")!=0) 
                            temperature += (float)273.15;  // convert C to K
                    }
                }
            } catch (AipsError & err) {
                noWeather = true;
            }
        } else {
            noWeather = true;
        }
    } else {
        noWeather = true;
    }
    if (noWeather) {
        parent_->logmesg("load_cache", "Could not open WEATHER table, using default values: humidity " + casacore::String::toString(humidity) + ", pressure " + casacore::String::toString(pressure) + ", temperature " + casacore::String::toString(temperature) + " for Atm/Tsky");
    }

    // to use in atmosphere.initAtmProfile (tool)
    weather_.define("humidity", humidity);       // %
    weather_.define("pressure", pressure);       // mb
    weather_.define("temperature", temperature); // K
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

}
