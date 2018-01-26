//# CalCache.h: CalTable-specific casacore::Data cache for plotms.
//# Copyright (C) 2009
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
//# $Id: $
#ifndef CALCACHE_H_
#define CALCACHE_H_

#include <plotms/Data/PlotMSCacheBase.h>

#include <plotms/PlotMS/PlotMSAveraging.h>
#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSFlagging.h>
#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTIter.h>
#include <synthesis/CalTables/SolvableVJMCol.h>
#include <synthesis/CalTables/CalDescColumns.h>
#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/Containers/Block.h>

namespace casa {

//# Forward declarations.
class PlotMSApp;
class PlotMSIndexer;

class CalCache : public PlotMSCacheBase {
    
  // Friend class declarations.
  friend class PlotMSIndexer;

public:    
  
  // Constructor which takes parent PlotMS.
  CalCache(PlotMSApp* parent);
  
  // Destructor
  virtual ~CalCache();

  // Identify myself
  PlotMSCacheBase::Type cacheType() const { return PlotMSCacheBase::CAL; };

  // Is the underlying table complex?
  inline casacore::Bool parsAreComplex() { return parsAreComplex_; };

  // ...not yet CAL-specific... (or ever?)
  // Set up indexing for the plot
  //  void setUpIndexer(PMS::Axis iteraxis=PMS::SCAN,
  //    casacore::Bool globalXRange=false, casacore::Bool globalYRange=false);

  // Convert poln index->name and name->index
  virtual casacore::String polname(casacore::Int ipol);

  // given filename, get cal type
  void setFilename(casacore::String filename);

protected:

  // CAL-specific loadIt method
  virtual void loadIt(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData,
    ThreadCommunication* thread = NULL);

private:
    
  // Forbid copy for now
  CalCache(const CalCache&);

  // NewCalTable:
  void setUpCalIter(const casacore::String& calname,
    PlotMSSelection& selection,
    casacore::Bool readonly=true,
    casacore::Bool chanselect=true,
    casacore::Bool corrselect=true);
  void countChunks(ROCTIter& ci,
    vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData,
    ThreadCommunication* thread);
  void loadCalChunks(ROCTIter& ci, const vector<PMS::Axis> loadAxes,
    ThreadCommunication* thread);
  void loadCalAxis(ROCTIter& cti, casacore::Int chunk, PMS::Axis axis,
    casacore::String pol);

  // CalTable:
  void countChunks(casacore::Int nrowMain, vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData,
    ThreadCommunication* thread);
  void loadCalChunks(ROSolvableVisJonesMCol& mcol,
    ROCalDescColumns& dcol, const vector<PMS::Axis> loadAxes,
    ThreadCommunication* thread);
  void loadCalAxis(ROSolvableVisJonesMCol& mcol, ROCalDescColumns& dcol,
    casacore::Int chunk, PMS::Axis axis, casacore::String pol);
  // helper functions for loading CalTable chunks
  casacore::String getMSAbsPath(casacore::String msname);
  void getChanFreqsFromMS(casacore::String fullmsname);
  void getBPolyDataAxis(PMS::Axis axis,
    casacore::Cube<casacore::Complex>& viscube, casacore::Int chunk);

  // Trap attempt to use to much memory (too many points)
  //  void trapExcessVolume(map<PMS::Axis,casacore::Bool> pendingLoadAxes);

  // Check axis and slice param column appropriately
  casacore::Slice getParSlice(casacore::String axis, casacore::String polnSel);
  // Get axis string for VisCal Slice code
  casacore::String toVisCalAxis(PMS::Axis axis);

  // Check for divide-by-zero (=inf); set to 1.0 and flag it
  void checkRatioArray(casacore::Array<float>& array, int chunk);

  // Set flags in the CalTable
  virtual void flagToDisk(const PlotMSFlagging& flagging,
    casacore::Vector<casacore::Int>& chunks, 
    casacore::Vector<casacore::Int>& relids,
    casacore::Bool flag,
    PlotMSIndexer* indexer, int index);

  // cal table iterator pointers
  ROCTIter* ci_p;
  CTIter* wci_p;

  // The polarization basis
  casacore::String basis_;
  // Is parameter column complex?
  casacore::Bool parsAreComplex_;
  // Check divide-by-zero in ratio plot (checkRatioArray)
  bool divZero_;
  // get from MeasurementSet
  casacore::Array<casacore::Double> chanfreqs_;
  // for CalTables
  casacore::Int nrow_;
  // Volume meter for volume calculation
  //  PMSCacheVolMeter vm_;
 
};
typedef casacore::CountedPtr<CalCache> CalCachePtr;


}

#endif /* CALCACHE_H_ */
