//# PlotMSIndexer.h: Cache indexer for plotms.
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
#ifndef PLOTMSINDEXER_H_
#define PLOTMSINDEXER_H_
#include <map>
#include <sstream>

#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/Containers/Block.h>
#include <graphics/GenericPlotter/PlotData.h>

#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/Data/PlotMSCacheBase.h>

#include <unordered_map>

namespace casa {

//# Forward declarations.
class PlotMSApp;
class PlotMSCacheBase;
class PlotMSIndexer;  // needed for method pointer typedefs

typedef casacore::Double(PlotMSCacheBase::*CacheMemPtr)(casacore::Int,casacore::Int);
typedef    casacore::Int(PlotMSIndexer::*IndexerMethPtr)(casacore::Int,casacore::Int);
typedef   void(PlotMSIndexer::*CollapseMethPtr)(casacore::Int,casacore::Array<casacore::Bool>&);


class PlotMSIndexer : public PlotMaskedPointData, public PlotBinnedData {

public:

  // Convenient access to class name.
  static const casacore::String CLASS_NAME;
  
  // A ctor that makes an empty Indexer (for plot initialization)
  PlotMSIndexer();
    
  // Constructor which takes parent PlotMSCache, x and y axes (non-iteration)
  PlotMSIndexer(PlotMSCacheBase* plotmscache, PMS::Axis xAxis, 
    PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData,
    casacore::String xconnect, bool timeconnect, int index);

protected:
  // Constructor which supports iteration
  PlotMSIndexer(PlotMSCacheBase* plotmscache, PMS::Axis xAxis, 
    PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData,
    PMS::Axis iterAxis, casacore::Int iterValue,
    casacore::String xconnect, bool timeconnect, int index);
  friend class PlotMSIndexerFactory;

public:

  // Destructor
  virtual ~PlotMSIndexer();

  // Implemented PlotData methods.
  // <group>
  bool willDeleteData() const { return true; }
  void setDeleteData(bool del = true)   { (void)del; }
  bool isValid() const { return true;};
  // </group>

  // Implemented PlotPointData methods.
  // <group>
  unsigned int size() const;
  virtual double xAt(unsigned int i) const;
  virtual double yAt(unsigned int i) const;
  virtual void xAndYAt(unsigned int index, double& x, double& y) const;
  bool minsMaxes(double& xMin, double& xMax, double& yMin, double& yMax);
  // </group>
    
  // Implemented PlotMaskedPointData methods.
  // <group>
  bool maskedAt(unsigned int index) const;
  virtual void xyAndMaskAt(unsigned int index, double& x, double& y, bool& mask) const;
  // </group>
    
  // Unimplemented PlotMaskedPointData methods.
  // <group>
  unsigned int sizeMasked() const { return sizeMasked_; }
  unsigned int sizeUnmasked() const { return sizeUnMasked_; }
  bool maskedMinsMaxes(double& xMin, double& xMax, double& yMin,double& yMax);
  bool unmaskedMinsMaxes(double& xMin, double& xMax, double& yMin,double& yMax);
  // </group>
    
  // PlotBinnedData methods
  // <group>
  unsigned int numBins() const;
  unsigned int binAt(unsigned int i) const;
  bool isBinned() const;  // affects colorizing points
  // </group>

  // Set up indexing for the plot 
  void setUpIndexing();

  // Set global min/max flag
  void setGlobalMinMax(casacore::Bool globalX=false,casacore::Bool globalY=false);
  bool isGlobalXRange() const;
  bool isGlobalYRange() const;

  // Report per-chunk point counters
  casacore::Vector<casacore::uInt> nPoints() { return nPoints_; };
  casacore::Vector<casacore::uInt> nCumulative() { return nCumulative_; };

  // Return if the indexer is ready (setUpPlot has been run)
  inline casacore::Bool indexerReady() const { return indexerReady_; };

  // Locate datum nearest to specified x,y
  casacore::Record getPointMetaData(casacore::Int i);
  casacore::Record locateInfo(const casacore::Vector<PlotRegion>& regions,
                    casacore::Bool showUnflagged, casacore::Bool showFlagged,
                    casacore::Bool selectAll = true);
  PlotLogMessage* locateRange(const casacore::Vector<PlotRegion>& regions,
      casacore::Bool showUnflagged, casacore::Bool showFlagged);
  PlotLogMessage* flagRange(const PlotMSFlagging& flagging,
      const casacore::Vector<PlotRegion>& regions, casacore::Bool flag = true);


  // Report meta info for current value of currChunk_/irel_
  void reportMeta(casacore::Double x, casacore::Double y, casacore::Bool masked, std::stringstream& ss);

  // Set flags in the cache
  void flagInCache(const PlotMSFlagging& flagging,casacore::Bool flag);

  // Iteration label
  casacore::String iterLabel();
  casacore::String iterValue();
  casacore::String fileLabel();

  // Access to raw min/max data (no auto-global)
  bool maskedMinsMaxesRaw(double& xMin, double& xMax, double& yMin,double& yMax);
  bool unmaskedMinsMaxesRaw(double& xMin, double& xMax, double& yMin,double& yMax);


  // Directly implemented index calculators
  //  (generic index methods point to one of these depending upon axis choice)
  casacore::Int getIndex0000(casacore::Int ch,casacore::Int irel) { return 0;  (void)irel; (void)ch; };
  casacore::Int getIndex1000(casacore::Int ch,casacore::Int irel) { return irel%icorrmax_(ch);};
  casacore::Int getIndex0100(casacore::Int ch,casacore::Int irel) { return (irel/nperchan_(ch))%ichanmax_(ch);};
  casacore::Int getIndex0010(casacore::Int ch,casacore::Int irel) { return (irel/nperbsln_(ch))%ibslnmax_(ch);};
  casacore::Int getIndex0110(casacore::Int ch,casacore::Int irel) { return (irel/nperchan_(ch))%ichanbslnmax_(ch);};
  casacore::Int getIndex1010(casacore::Int ch,casacore::Int irel) { return (irel/nperbsln_(ch))*nperchan_(ch) + irel%nperchan_(ch);};
  casacore::Int getIndex1110(casacore::Int ch,casacore::Int irel) { return irel%idatamax_(ch);};
  casacore::Int getIndex0001(casacore::Int ch,casacore::Int irel) { return (irel/nperant_(ch))%iantmax_(ch);};


  // set colorize and whether binned data has coloraxis;
  // connected points are binned too
  bool colorize(bool doColorize, PMS::Axis colorizeAxis);

  bool setConnect(casacore::String xconnect, bool timeconnect);
  virtual bool reverseConnect(unsigned int index) const;
  unsigned int connectBinAt(unsigned int i) const;

  bool plotConjugates() const { return (PMS::axisIsUV(currentX_) && 
          PMS::axisIsUV(currentY_)); }

protected:
  PlotMSCacheBase* getCache() { return plotmscache_; }
  // Set currChunk_ according to a supplied index
  void setChunk(casacore::uInt i, bool ignoreReindex=false) const;
  void setReady(bool isReady=true) { indexerReady_ = isReady; };

private:
    
  // Forbid copy for now
  PlotMSIndexer(const PlotMSIndexer& mc);

  // get method for x and y axes depends on column
  void setMethod(CacheMemPtr& getmethod, PMS::Axis axis, PMS::DataColumn data);

  // index for iteraxis
  void setIndexer(IndexerMethPtr& indexmethod, PMS::Axis axis);

  // Reindex for connecting points
  void reindexForConnect();

  //  void setCollapser(CollapseMethPtr& collmethod, PMS::Axis axis);

  // Generate collapsed versions of the plmask 
  /* not needed?  (gmoellen, 2011Mar15)
  void collapseMask0000(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask1000(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask0100(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask0010(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask0110(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask1010(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask1110(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  void collapseMask0001(casacore::Int ch,casacore::Array<casacore::Bool>& collmask);
  */

  // Report the number of chunks
  casacore::Int nChunk() const; // { return (plotmscache_ ? plotmscache_->nChunk() : 0); };

  // Report the reference time for this cache (in seconds)
  inline casacore::Double refTime(); // { return plotmscache_->refTime(); };

  // Computes the X and Y limits for the currently set axes.  In the future we
  // may want to cache ALL ranges for all loaded values to avoid recomputation.
  void computeRanges();

  // Adjust the Y range for overlays to make room at top of plot
  void adjustYRange(double& yMin, double& yMax);

  // Compute baseline's length in meters between ant1 and ant2
  casacore::Double computeBaselineLength(casacore::Int ant1, casacore::Int ant2);

  // Convenience methods that call log() with the given method name and the
  // appropriate event type.
  // <group>
  void logInfo(const casacore::String& method, const casacore::String& message) {
      log(method, message, PlotLogger::MSG_INFO); }
  void logDebug(const casacore::String& method, const casacore::String& message) {
      log(method, message, PlotLogger::MSG_DEBUG); }
  void logWarn(const casacore::String& method, const casacore::String& message) {
      log(method, message, PlotLogger::MSG_WARN); }
  void logError(const casacore::String& method, const casacore::String& message) {
      log(method, message, PlotLogger::MSG_ERROR); }
  // </group>
  
  // Logs the given message from the given method name as the given event type
  // (see PlotLogger).
  void log(const casacore::String& method, const casacore::String& message, int eventType);



  // Private data
   
  // Parent plotms.
  PlotMSCacheBase* plotmscache_;

  // Pointers to methods for axis flexibility
  CacheMemPtr getXFromCache_, getYFromCache_,getColFromCache_;
  IndexerMethPtr XIndexer_, YIndexer_, ColIndexer_;
  //  CollapseMethPtr collapseXMask_, collapseYMask_;

protected:
  // The in-focus chunk and relative index offset
  mutable casacore::Int currChunk_, irel_;

private:
  mutable casacore::uInt lasti_;

  // The number of points per chunk
  casacore::Vector<casacore::uInt> nPoints_;

  // The cumulative running total of points
  casacore::Vector<casacore::uInt> nCumulative_;

  // Segment point-counting Vectors
  casacore::Int nSegment_;
  mutable casacore::Int currSeg_;
  casacore::Vector<casacore::uInt> nSegPoints_,nCumulPoints_,cacheChunk_,cacheOffset_;
  
  // Current setup/state.
  PMS::Axis currentX_, currentY_;
  PMS::DataColumn currentXdata_, currentYdata_;
  bool indexerReady_, connectReady_;

  // Indexing parameters
  casacore::Vector<casacore::Int> icorrmax_, ichanmax_, ibslnmax_, idatamax_;
  casacore::Vector<casacore::Int> nperchan_, nperbsln_, nperant_;
  casacore::Vector<casacore::Int> ichanbslnmax_;
  casacore::Vector<casacore::Int> iantmax_;

  // Nominal axes limits
  casacore::Double xmin_,ymin_,xflmin_,yflmin_,xmax_,ymax_,xflmax_,yflmax_;
  casacore::Int sizeMasked_, sizeUnMasked_;
  casacore::Bool globalXMinMax_,globalYMinMax_;

  // Iteration
  // <group>
  casacore::Bool iterate_;
  PMS::Axis iterAxis_;
  casacore::Int iterValue_;
  // </group>

  // Colorization
  // <group>
  bool itsColorize_;
  PMS::Axis itsColorizeAxis_;
  // </group>

  // Reindex and bin when connecting points
  casacore::String itsXConnect_;
  bool itsTimeConnect_;
  // map<chunk, decrease> : frequencies decrease with channel in chunk
  // (for reverseConnect)
  std::map<casacore::uInt, bool> freqsDecrease_;

  // Cope with const-ness in the get methods
protected:
  PlotMSIndexer* self;

private:

  int dataIndex_;

};

typedef casacore::CountedPtr<PlotMSIndexer> PlotMSIndexerPtr;

class PlotMSRaDecIndexer : public PlotMSIndexer {

public:

	// Convenient access to class name.
	static const casacore::String CLASS_NAME;

	using RaDecData = casacore::PtrBlock<casacore::Vector<casacore::Double>*>;
	using RaDecMap = std::map<DirectionAxisParams,RaDecData>;

	static const RaDecData EMPTY_DATA;

	// Empty Indexer
	PlotMSRaDecIndexer();

	// Constructor which takes parent PlotMSCache, x and y axes (non-iteration)
	//PlotMSRaDecIndexer(PlotMSCacheBase* plotmscache, PMS::Axis xAxis,
	//PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData, int index)
	//: PlotMSIndexer(plotmscache, xAxis,xData,yAxis,yData,index) {}

	// Constructor which supports iteration
	PlotMSRaDecIndexer(PlotMSCacheBase* plotmscache, PMS::Axis xAxis,
	PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData,
	PMS::Axis iterAxis, casacore::Int iterValue, casacore::String xconnect,
	bool timeconnect, int index);

	// Destructor
	~PlotMSRaDecIndexer() {}

	// Overridden PlotPointData methods.
	double xAt(unsigned int i) const;
	double yAt(unsigned int i) const;
	void xAndYAt(unsigned int index, double& x, double& y) const;

	void xyAndMaskAt(unsigned int index,
			double& x, double& y,
			bool& mask) const;

	inline casacore::Double getRaX(casacore::Int chnk,casacore::Int irel) const {
		return *(xRa_[chnk]->data()+irel);
	}
	inline casacore::Double getRaY(casacore::Int chnk,casacore::Int irel) const {
		return *(yRa_[chnk]->data()+irel);
	}
	inline casacore::Double getDecX(casacore::Int chnk,casacore::Int irel) const {
		return *(xDec_[chnk]->data()+irel);
	}
	inline casacore::Double getDecY(casacore::Int chnk,casacore::Int irel) const {
		return *(yDec_[chnk]->data()+irel);
	}

private:


	const bool cacheOk_;
	const bool axesOk_;
	const bool xAxisIsRa_;
	const bool xAxisIsDec_;
	const bool xAxisIsRaOrDec_;
	const bool yAxisIsRa_;
	const bool yAxisIsDec_;
	const bool yAxisIsRaOrDec_;

	const RaDecData& xRa_;
	const RaDecData& xDec_;

	const RaDecData& yRa_;
	const RaDecData& yDec_;
};

class PlotMSIndexerFactory {
public:
static PlotMSIndexer* initIndexer(PlotMSCacheBase* plotmscache, PMS::Axis xAxis,
  PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData,
  PMS::Axis iterAxis, casacore::Int iterValue, casacore::String xconnect,
  bool timeconnect, int index, bool makeRaDecIndexer=false);
};

}

#endif /* PLOTMSCACHEINDER_H_ */
