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
#include <synthesis/CalTables/BJonesMCol.h>
#include <synthesis/CalTables/GJonesMCol.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>
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
  virtual void loadIt(std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData,
      ThreadCommunication* thread = nullptr);

private:
    
  // Forbid copy for now
  CalCache(const CalCache&);

  // NewCalTable:
  void loadNewCalTable(std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData, ThreadCommunication* thread = nullptr);
  void setUpCalIter(NewCalTable& selct, casacore::Bool readonly=True);
  void countChunks(ROCTIter& ci,
      std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData,
      ThreadCommunication* thread);
  void loadCalChunks(ROCTIter& ci, const std::vector<PMS::Axis> loadAxes,
      ThreadCommunication* thread);
  void loadCalAxis(ROCTIter& cti, casacore::Int chunk, PMS::Axis axis,
      casacore::String pol);
  virtual void flagToDisk(const PlotMSFlagging& flagging,
      casacore::Vector<casacore::Int>& chunks, 
      casacore::Vector<casacore::Int>& relids,
      casacore::Bool flag, PlotMSIndexer* indexer, int index);

  // CalTable:
  void countChunks(casacore::Int nrowMain, std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData, ThreadCommunication* thread);
  void setMSname(casacore::String msname); // set msname_; adds path to name
  void getNamesFromMS();    // for locate
  void setUpLoad(ThreadCommunication* thread, casacore::Slice& parSlice);
  void getCalDataAxis(PMS::Axis axis, casacore::Cube<casacore::Complex>& viscube,
      casacore::Int chunk);  // get axes derived from raw viscube

  // BPOLY CalTable:
  void loadBPoly(std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData, ThreadCommunication* thread = nullptr);
  void loadCalChunks(ROBJonesPolyMCol& mcol, ROCalDescColumns& dcol,
      casacore::Int nrow, const std::vector<PMS::Axis> loadAxes,
	  casacore::Vector<casacore::Vector<casacore::Slice> >& chansel,
	  ThreadCommunication* thread);
  void loadCalAxis(ROSolvableVisJonesMCol& mcol, ROCalDescColumns& dcol,
      casacore::Int chunk, PMS::Axis axis);
  void getChanFreqsFromMS(casacore::Vector< casacore::Vector<casacore::Double> >& mschanfreqs);
  void getSelFreqsForSpw(casacore::Vector<casacore::Slice>& chansel,
	  casacore::Vector<casacore::Double>& chanFreqs,
      casacore::Vector<casacore::Int>& chanNums);
  // cube selected by channel:
  template<class T>
  void getSelectedCube(const casacore::Cube<T>& inputCube,
	const casacore::Vector<casacore::Slice>& chanSlices,
	casacore::Cube<T>& outputCube);

  // GSPLINE CalTable:
  void loadGSpline(std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData, ThreadCommunication* thread = nullptr);
  void loadCalChunks(ROGJonesSplineMCol& mcol, ROCalDescColumns& dcol,
      casacore::Int nsample, const std::vector<PMS::Axis> loadAxes,
	  casacore::Vector<int>& selectedAnts, ThreadCommunication* thread);
  void checkAxes(const std::vector<PMS::Axis>& loadAxes);
  // cube selected by antenna1:
  template<class T>
  void getSelectedCube(casacore::Cube<T>& inputcube,
      const casacore::Vector<casacore::Int> selectedRows);

  // Utilities for all cal tables:
  // Get axis string for VisCal Slice code
  casacore::String toVisCalAxis(PMS::Axis axis);
  // Check axis and slice param column appropriately
  casacore::Slice getParSlice(casacore::String axis, casacore::String polnSel);
  // Check for divide-by-zero (=inf); set to 1.0 and flag it
  void checkRatioArray(casacore::Array<float>& array, int chunk);
  // Trap attempt to use to much memory (too many points)
  //  void trapExcessVolume(map<PMS::Axis,casacore::Bool> pendingLoadAxes);

  // All
  // Check divide-by-zero in ratio plot (checkRatioArray)
  bool divZero_;
  // Volume meter for volume calculation
  //  PMSCacheVolMeter vm_;

  // NewCalTable iterator pointers
  ROCTIter* ci_p;
  CTIter* wci_p;
  // The polarization basis
  casacore::String basis_;
  // Is parameter column complex?
  casacore::Bool parsAreComplex_;

  // CalTable (cannot plot BPOLY or GSPLINE without MS)
  casacore::String msname_;
};
typedef casacore::CountedPtr<CalCache> CalCachePtr;


}

#endif /* CALCACHE_H_ */
