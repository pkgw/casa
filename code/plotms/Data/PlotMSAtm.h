//# PlotMSAtm.h: calculate atmospheric transmission curve
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

#ifndef PLOTMSATM_H_
#define PLOTMSATM_H_

#include <plotms/Data/PlotMSCacheBase.h>
#include <plotms/PlotMS/PlotMSSelection.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <synthesis/CalTables/NewCalTable.h>
#include <atmosphere/ATM/ATMProfile.h>

namespace casa {

// <summary> 
// PlotMSAtm: plotms class for calculating atmospheric transmission curve
// for bandpass plots
// </summary>

// <reviewed reviewer="" date="" tests="" demos="">

// <prerequisite>
//   <li> <linkto class="CalCache">CalCache</linkto> module
// </prerequisite>
//
// <etymology>
// From "plotms", "atmospheric".
// </etymology>
//
// <synopsis>
// The PlotMSAtm class computes atmospheric transmission curves based on the
// algorithm in the plotbandpass task.  Needed information is obtained from
// the bandpass table and its MeasurementSet if possible, else defaults based
// on observatory are used. The calculated values are returned in a vector
// whose size is equal to the number of channels.
// </synopsis>
//
// <example>
// <srcblock>
// </srcblock>
// </example>
//
// <motivation>
// This class is used by the plotms CalCache class to obtain atmospheric
// transmission curve values for plotting the curve overlay on a bandpass
// plot.
// </motivation>
//

class PlotMSAtm {

public:

    // construct with bandpass table name
    PlotMSAtm(casacore::String filename, PlotMSSelection& userSel, bool isMS,
        PlotMSCacheBase* parent);
    ~PlotMSAtm();

    // Returns curve vector (atm if atm=true, else tsky);
    // pass in frequencies (GHz) in case of channel-averaging
    // (let VIVB2 do the work)
    casacore::Vector<casacore::Double> calcOverlayCurve(
        casacore::Int spw, casacore::Int scan, 
        const casacore::Vector<casacore::Double>& chanFreqs,
        bool atm);

    inline casacore::Double getPwv() { return pwv_; }
    inline casacore::Double getAirmass() { return airmass_; }

private:

    PlotMSAtm(const PlotMSAtm& other);
    PlotMSAtm& operator= (const PlotMSAtm& other);

    // info from MS
    void setUpMS(casacore::String filename, PlotMSSelection& userSel);
    void getMSTimes();
    void getMSFields();

    // info from cal tables
    void setUpCalTable(casacore::String filename, PlotMSSelection& userSel);
    void getCalTimes();
    void getCalFields();
    void getCalMS();

    // for user selection then each chunk's spw and scan
    void applyMSSelection(PlotMSSelection& selection,
            casacore::MeasurementSet& selms);
    void applyCalSelection(PlotMSSelection& selection,
            NewCalTable& selct);

    // calculated values
    void getMeanWeather();  // stored in weather_ Record
	casacore::Table selectWeatherTable(casacore::Table& intable,
        casacore::String tempUnits);
    void getMedianPwv();    // stored in pwv_
    casacore::Double computeMeanAirmass();
    casacore::Double getElevation(casacore::Int fieldId);
    casacore::Double getMeanScantime();

    // atmosphere tool
    atm::AtmProfile* getAtmProfile();

    // utility functions
    // Sets times_ = unique values in TIME column
    void getUniqueTimes(casacore::Vector<casacore::Double> alltimes);
    // Sets fields_ = unique values in FIELD column
    void getUniqueFields(casacore::Vector<casacore::Int> allfields);
    template <typename T>
    casacore::Vector<T> getValuesNearTimes(
        casacore::Vector<T> inputCol, 
        casacore::Vector<casacore::Double> timesCol);

    bool isMS_;
    PlotMSCacheBase* parent_; // for log messages
    casacore::MeasurementSet *ms_, *selms_; // selected MS for each spw/scan
    NewCalTable *caltable_, *selct_;  // selected CT for each spw/scan
    casacore::String tableName_, telescopeName_;
    casacore::Double pwv_, airmass_;
    casacore::Vector<casacore::Double> times_;
    casacore::Vector<casacore::Int> fields_;
    casacore::Record weather_;
	const unsigned int MAX_ATM_CALC_CHAN_;
};

}
#endif
