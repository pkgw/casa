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

using namespace casacore;

namespace casa {

PlotMSAtm::PlotMSAtm(casacore::String bptable):
    telescopeName_(""),
    ms_(NULL) {

    // set up cal table and get needed data
    bptable_ = new NewCalTable(bptable);
    setCalTimes();
    setTelescope();
    setFields();
    setMS();
    selection_ = PlotMSSelection();
}

PlotMSAtm::~PlotMSAtm() {
    if (bptable_) {
        delete bptable_;
        bptable_ = NULL;
    }
    if (ms_) {
        delete ms_;
        ms_ = NULL;
    }
}

void PlotMSAtm::setCalTimes() {
    // unique times in TIME column of bandpass table
    ROCTMainColumns ctmc(*bptable_);
    casacore::Vector<casacore::Double> times = ctmc.time().getColumn();
    // find unique times
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(times.data(), TpDouble);
    sorter.sort(indexVector, times.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    caltimes_.resize(nUnique);
    for (uInt i=0; i<nUnique; ++i) {
        caltimes_(i) = times.data()[indexVector(uniqueVector(i))];
    }
} 

void PlotMSAtm::setTelescope() {
    ROCTColumns ctCol(*bptable_);
    telescopeName_ = ctCol.observation().telescopeName().get(0);
}

void PlotMSAtm::setFields() {
    // unique fields in FIELD_ID column of bandpass table
    // needed later for elevation calculation for airmass
    ROCTMainColumns ctmc(*bptable_);
    casacore::Vector<casacore::Int> fields = ctmc.fieldId().getColumn();
    // find unique fields
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(fields.data(), TpInt);
    sorter.sort(indexVector, fields.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    calfields_.resize(nUnique);
    for (uInt i=0; i<nUnique; ++i) {
        calfields_(i) = fields.data()[indexVector(uniqueVector(i))];
    }
}

casacore::String PlotMSAtm::getSelectionExpr(casacore::Vector<casacore::Int> intVector) {
    // convert Int array to comma-separated selection string
    casacore::String selExpr;
    for (uInt i=0; i<intVector.size(); ++i) {
        selExpr += String::toString(intVector(i));
        if (i < intVector.size()-1) selExpr += ",";
    }
    return selExpr;
}

void PlotMSAtm::setMS() {
    // if caltable has keyword for ms name, sets ms_ 
    String msname("");
    if (bptable_->keywordSet().fieldNumber("MSName") > -1) 
        msname = bptable_->keywordSet().asString("MSName");
    if (!msname.empty()) {
        casacore::Path path(bptable_->tableName());
        casacore::String fullpath = path.dirName() + "/" + msname;
        casacore::File msname(fullpath);
        if (msname.exists() && msname.isDirectory())
            ms_ = new MeasurementSet(fullpath);
    }
}

void PlotMSAtm::setSelection(PlotMSSelection& selection) {
    selection_ = selection;
    // apply selection to caltable
    Vector<Vector<Slice> > chansel, corrsel;
    NewCalTable selct;
    selection.apply(*bptable_, selct, chansel, corrsel);
    delete bptable_;
    bptable_ = new NewCalTable(selct);
    setCalTimes();  // update based on selection
}

casacore::Vector<casacore::Int> PlotMSAtm::calcAtmTransmission() {
    // Implements algorithm in CAS-9053 to get atmospheric transmission curve
    casacore::Vector<casacore::Int> transmission;
    // Get info from MS
    casacore::Double pwv = getMedianPwv();
    casacore::Record weather = getMeanWeather();
    casacore::Double airmass = computeMeanAirmass();
    return transmission;
}

casacore::Double PlotMSAtm::getMedianPwv() {
    // Get pwv (precipitable water vapor) from MS subtable
    // or use default value for telescope
    casacore::Double pwv(0.0);

    if (ms_) {
        // values from subtables if they exist
        casacore::String msname(ms_->tableName()), subname;
        casacore::Table mstab(msname), subtable;
        casacore::Vector<casacore::Double> waterCol, timesCol;
        if (mstab.keywordSet().fieldNumber("ASDM_CALWVR") > -1) {
            subname = msname + "::ASDM_CALWVR";
            subtable = Table::openTable(subname);
            waterCol = ScalarColumn<casacore::Double>(subtable, "water").getColumn();
            timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            mstab.closeSubTables();
        } else if (mstab.keywordSet().fieldNumber("ASDM_CALATMOSPHERE") > -1) {
            subname = msname + "::ASDM_CALATMOSPHERE";
            subtable = Table::openTable(subname);
            waterCol = ArrayColumn<casacore::Double>(subtable, "water").getColumn();
            timesCol = ScalarColumn<casacore::Double>(subtable, "startValidTime").getColumn();
            mstab.closeSubTables();
        }
        if (!waterCol.empty()) {
            casacore::Vector<casacore::Double> water = 
                getValuesNearCalTimes(waterCol, timesCol);
            if (!water.empty())
                pwv = median(water) * 1000.0; // convert to mm
        }
    }

    if (pwv == 0.0) {  // no asdm tables, use default value
        if (isAlma()) pwv = 1.0;
        else pwv = 5.0;
    }
    return pwv;
}

casacore::Record PlotMSAtm::getMeanWeather() {
    // Info from MS WEATHER table
    // set defaults
    casacore::Float pressure, humidity(20.0), temperature(273.15);
    if (isAlma())
        pressure = 563.0;  // mb
    else 
        pressure = 786.0;  // mb

    // values from WEATHER table if it exists
    if (ms_) {
        casacore::String msname(ms_->tableName()), subname;
        casacore::Table mstab(msname), subtable;
        if (mstab.keywordSet().fieldNumber("WEATHER") > -1) {
            subname = msname + "::WEATHER";
            subtable = Table::openTable(subname);
            // get columns and temp. units
            casacore::Vector<casacore::Float> pressureCol, humidityCol, temperatureCol;
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
            casacore::Vector<casacore::Float> selpressure, selhumidity, seltemperature;
            // pressure
            if (!pressureCol.empty()) {
                selpressure = getValuesNearCalTimes(pressureCol, timeCol);
                if (!selpressure.empty())
                    pressure = mean(selpressure);
            }
            // humidity
            if (!humidityCol.empty()) {
                selhumidity = getValuesNearCalTimes(humidityCol, timeCol);
                if (!selhumidity.empty())
                    humidity = mean(selhumidity);
            }
            // temperature
            if (!temperatureCol.empty()) {
                seltemperature = getValuesNearCalTimes(temperatureCol, timeCol);
                if (!seltemperature.empty()) {
                    temperature = mean(seltemperature);
                    if (units.compare("K")!=0) 
                        temperature += (float)273.15;  // convert C to K
                }
            }
        }
    }

    // to use in atmosphere.initAtmProfile (tool)
    casacore::Record weatherInfo;
    weatherInfo.define("humidity", humidity);       // %
    weatherInfo.define("pressure", pressure);       // mb
    weatherInfo.define("temperature", temperature); // K
    return weatherInfo;
}

casacore::Double PlotMSAtm::computeMeanAirmass() {
    // Calculate airmass from elevation
    casacore::Vector<casacore::Double> airmasses;
    airmasses.resize(calfields_.size());
    casacore::Double elevation;
    for (uInt i=0; i<calfields_.size(); ++i) {
        elevation = getElevation(calfields_(i));
        if (elevation <= 3.0) elevation = 45.0;
        airmasses(i) = 1.0 / std::cos((90.0 - elevation) * C::pi / 180.0);
    }
    casacore::Double airmass = mean(airmasses);
    return airmass;
}

casacore::Double PlotMSAtm::getElevation(casacore::Int fieldId) {
    // Get RADec (DELAY_DIR) from FIELD table
    ROCTColumns ctCol(*bptable_);
    casacore::Array<casacore::Double> raDec = 
        ctCol.field().delayDir().get(fieldId);
    casacore::Double ra = raDec(IPosition(2,0,0));
    casacore::Double dec = raDec(IPosition(2,1,0));
    // convert J2000 to AzEl
    // need observatory position and timestamp for output reference
    casacore::MPosition obspos;
    casacore::MeasTable::Observatory(obspos, telescopeName_);
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
    // get mean timestamp for first scan in bandpass table
    casacore::Double meantime(mean(caltimes_));
    ROCTMainColumns ctmc(*bptable_);
    casacore::Vector<casacore::Int> scanCol = ctmc.scanNo().getColumn();
    // sort and unique scan numbers
    Sort sorter;
    casacore::Vector<casacore::uInt> indexVector, uniqueVector;
    uInt nUnique;
    sorter.sortKey(scanCol.data(), TpInt);
    sorter.sort(indexVector, scanCol.size());
    nUnique = sorter.unique(uniqueVector, indexVector);
    if (nUnique > 1) {
        // find times for first scan
        uInt end(uniqueVector(1));
        casacore::Vector<casacore::Double> timeCol = ctmc.time().getColumn();
        casacore::Vector<casacore::Double> scantimes;
        scantimes.resize(end);
        for (uInt i=0; i<end; ++i) {
            scantimes(i) = timeCol.data()[indexVector(i)];
        }
        meantime = mean(scantimes);
    }
    return meantime;
}

template <typename T>
casacore::Vector<T> PlotMSAtm::getValuesNearCalTimes(
        casacore::Vector<T> inputCol, 
        casacore::Vector<casacore::Double> timesCol) {
    // Use values with timestamps close to cal times
    // first, find time difference in closest match
    casacore::Double mintimediff(1.0e12);
    for (uInt i=0; i<caltimes_.size(); ++i) {
        for (uInt j=0; j<timesCol.size(); ++j) {
            mintimediff = min(mintimediff, abs(caltimes_(i)-timesCol(j)));
        }
    }
    // now get index of values within 1 sec of minimum timediff
    casacore::Vector<casacore::uInt> values;
    casacore::uInt vsize(0);
    for (casacore::uInt i=0; i<caltimes_.size(); ++i) {
        for (casacore::uInt j=0; j<timesCol.size(); ++j) {
            if ((abs(caltimes_(i)-timesCol(j)) - mintimediff) < 1.0) {
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
