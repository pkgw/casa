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
    PlotMSAtm(casacore::String filename, PlotMSSelection& userSel, bool showatm,
        bool isMS, bool xAxisIsChan, PlotMSCacheBase* parent);
    ~PlotMSAtm();

    // accessors
    inline casacore::String filename() { return filename_; }
    inline PlotMSSelection selection() { return selection_; }
    inline bool showatm() { return showatm_; } // false is tsky
    inline bool xAxisIsChan() { return xIsChan_; }

    inline void setShowAtm(bool showatm) { showatm_ = showatm; }
    inline void setXAxisIsChan(bool isChan) { xIsChan_ = isChan; }

    // passes arguments through to calcOverlayCurve
    void calcAtmTskyCurve(casacore::Vector<casacore::Double>& curve,
        casacore::Int spw, casacore::Int scan,
        const casacore::Vector<casacore::Double>& chanFreqs);
    // calculates image frequencies then calls calcOverlayCurve
    void calcImageCurve(casacore::Vector<casacore::Double>& curve,
        casacore::Int spw, casacore::Int scan,
        const casacore::Vector<casacore::Double>& chanFreqs);

    inline casacore::Double getPwv() { return pwv_; }
    inline casacore::Double getAirmass() { return airmass_; }

    // image sideband curve helpers
    inline bool canShowImageCurve() { return (hasReceiverTable() && canGetLOsForSpw()); }
    bool hasReceiverTable();
    bool canGetLOsForSpw();

private:

    PlotMSAtm(const PlotMSAtm& other);
    PlotMSAtm& operator= (const PlotMSAtm& other);

    // info from MS
    void setUpMS(casacore::String filename, PlotMSSelection& userSel);
    void getMSTimes(MeasurementSet& ms);
    void getMSFields(MeasurementSet& ms);

    // info from cal tables
    void setUpCalTable(casacore::String filename, PlotMSSelection& userSel);
    void getCalTimes(NewCalTable& ct);
    void getCalFields(NewCalTable& ct);
    void getCalMS(); // uses original caltable_

    // common function for plotbandpass CalcAtmTransmission algorithm
    // Returns curve vector (atm, tsky, image sideband);
    casacore::Vector<casacore::Double> calcOverlayCurve(
        casacore::Int spw, casacore::Int scan,
        const casacore::Vector<casacore::Double>& chanFreqs);

    // for user selection then each chunk's spw and scan
    void applyMSSelection(PlotMSSelection& selection,
            casacore::MeasurementSet& selms);
    void applyCalSelection(PlotMSSelection& selection,
            NewCalTable& selct);

    // calculated values
    void getMeanWeather();  // stored in weather_ Record
    casacore::Table selectWeatherTable(casacore::Table& intable,
        casacore::String tempUnits, casacore::String pressureUnits);
    void getMedianPwv();    // stored in pwv_
    casacore::Double computeMeanAirmass();
    casacore::Double getPointingElevation();
    casacore::Double getFieldElevation(casacore::Int fieldId); // no pointing table
    casacore::Double getMeanScantime();

    // atmosphere tool
    atm::AtmProfile* getAtmProfile();

    // image sideband curve
    bool getLO1FreqForSpw(double& freq, int spw);
    bool calcImageFrequencies(casacore::Vector<casacore::Double>& imageFreqs,
        casacore::Int spw, const casacore::Vector<casacore::Double>& chanFreqs);

    // utility functions
    // Determine unique time values in input vector
    void getUniqueTimes(casacore::Vector<casacore::Double> inputTimes,
        casacore::Vector<casacore::Double>& uniqueTimes);
    // Sets fields_ = unique values in FIELD column
    void getUniqueFields(casacore::Vector<casacore::Int> allfields);
    template <typename T>
    casacore::Vector<T> getValuesInTimeRange(casacore::Vector<T> inputCol, 
        casacore::Vector<casacore::Double> timesCol, casacore::Double mintime,
        casacore::Double maxtime);
    template <typename T>
    void getClosestValues(casacore::Vector<T>& values,
        casacore::Vector<casacore::Double>& times, casacore::Vector<T>& data,
        double mintime, double maxtime);
    // use cal table times if available, else ms times
    void getTimeRange(casacore::Double& mintime, casacore::Double& maxtime);

    casacore::String filename_, tableName_, telescopeName_;
    bool showatm_;         // true=showatm, false=showtsky
    bool isMS_;            // true=MS, false=CalTable
    bool xIsChan_;         // image curve changes for chan/freq x-axis
    bool canCalculatePwv_;     // has CALWVR or CALATMOSPHERE subtable
    bool canCalculateWeather_; // has WEATHER subtable
    PlotMSCacheBase* parent_; // for log messages
    PlotMSSelection selection_;
    casacore::MeasurementSet *ms_, *selms_; // selected MS for each spw/scan
    NewCalTable *caltable_, *selct_;        // selected CT for each spw/scan

    // updated for every spw/scan selection:
    int selectedSpw_, selectedScan_;
    casacore::Double pwv_, airmass_;
    casacore::Vector<casacore::Double> mstimes_, caltimes_;
    casacore::Vector<casacore::Int> fields_;
    casacore::Record weather_;
    std::map<int, double> loFreqForSpw_;
    const unsigned int MAX_ATM_CALC_CHAN_;
};

}
#endif
