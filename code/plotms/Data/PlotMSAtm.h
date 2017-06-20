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
    PlotMSAtm(casacore::String bptable);
    ~PlotMSAtm();

    // apply selection to cal table
    void setUserSelection(PlotMSSelection& selection);

    // returns transmission vector
    casacore::Vector<casacore::Double> calcAtmTransmission(casacore::Int spw,
            casacore::Int scan);

    inline casacore::Double getPwv() { return pwv_; }
    inline casacore::Double getAirmass() { return airmass_; }

private:

    PlotMSAtm(const PlotMSAtm& other);
    PlotMSAtm& operator= (const PlotMSAtm& other);

    // info from cal tables
    void setCalTimes(NewCalTable* nct=NULL);
    void setTelescope();
    void setFields();
    void setMS();

    // for each chunk's spw and scan
    NewCalTable* applySelection(PlotMSSelection& selection);

    inline bool isAlma() { return (telescopeName_=="ALMA"); }

    // calculated values
    void getMeanWeather();  // stored in weather_ Record
    void getMedianPwv();    // stored in pwv_
    casacore::Double computeMeanAirmass();
    casacore::Double getElevation(casacore::Int fieldId);
    casacore::Double getMeanScantime();

    // atmosphere tool
    atm::AtmProfile* getAtmProfile();

    // utility functions
    template <typename T>
    casacore::Vector<T> getValuesNearCalTimes(
        casacore::Vector<T> inputCol, 
        casacore::Vector<casacore::Double> timesCol);
    casacore::String getSelectionExpr(
        casacore::Vector<casacore::Int> intVector);

    NewCalTable *bptable_, *seltable_;
    casacore::Vector<casacore::Double> caltimes_;
    casacore::String telescopeName_;
    casacore::Vector<casacore::Int> calfields_;
    casacore::MeasurementSet* ms_;
    casacore::Int scan_;  // needed for elevation
    casacore::Record weather_;
    casacore::Double pwv_, airmass_;
};

}
#endif
