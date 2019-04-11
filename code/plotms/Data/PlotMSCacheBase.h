//# PlotMSCacheBase.h: Generic Data cache for plotms.
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
#ifndef PLOTMSCACHEBASE_H_
#define PLOTMSCACHEBASE_H_

#include <plotms/PlotMS/PlotMSAveraging.h>
#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSFlagging.h>
#include <plotms/PlotMS/PlotMSTransformations.h>
#include <plotms/PlotMS/PlotMSCalibration.h>

#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/PlotMS/PlotMSParameters.h>

#include <plotms/Data/PageHeaderCache.h>

#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/Containers/Block.h>
#include <measures/Measures/MFrequency.h>

#include <QVector>
#include <map>

namespace casa {

//# Forward declarations.
class PlotMSApp;
class PlotMSPlot;
class PlotMSIndexer;
class ThreadCommunication;
class PlotMSAtm;


class PlotMSCacheBase {
  
  // Friend class declarations.
  friend class PlotMSIndexer;

  //TBD:    friend class PlotMSData;

public:    

  // Varieties of cache
  // TBD: move to PlotMSConstants?
  enum Type {MS, CAL};
  
  static const unsigned int THREAD_SEGMENT;
  static const PMS::Axis METADATA[];
  static const unsigned int N_METADATA;
    
  static bool axisIsMetaData(PMS::Axis axis);

  
  // Constructor which takes parent PlotMS.
  PlotMSCacheBase(PlotMSApp* parent, PlotMSPlot* plot = nullptr);
  
  // Destructor
  virtual ~PlotMSCacheBase();

  // Identify myself
  // (MS or CAL)
  virtual PlotMSCacheBase::Type cacheType() const = 0;

  // Access to pol names
  virtual casacore::String polname(casacore::Int ipol)=0;

  // keep MS/CT filename (sets calType_)
  virtual void setFilename(casacore::String filename) = 0;
  casacore::String calType() const { return calType_; };
  bool polnRatio() const { return polnRatio_; };

  // Meta axes info
  int nmetadata() const {return N_METADATA;};
  PMS::Axis metadata(int i) {return METADATA[i];};

  // loaded ATM or TSKY
  bool hasOverlay();

  // Reference an indexer; returns -1 if there is no indexer
  // for the given dataIndex.
  PlotMSIndexer& indexer( int dataIndex, casacore::uInt i) {
      return (*indexer_[dataIndex][i]);
  };
  PlotMSIndexer& indexer0() {
      return *indexer0_;
  };
  void resizeIndexer( int size );
  int getDataCount() const {
      return currentX_.size();
  }
  casacore::Int nIter( int dataIndex ) const;

  PMS::Axis getIterAxis() const;

  // Report the number of chunks
  casacore::Int nChunk() const { return nChunk_; };

  // Returns whether cache is filled
  bool cacheReady() const { return dataLoaded_; }

  // Returns whether user canceled during loading chunks
  bool wasCanceled() const { return userCanceled_; }

  // Report the data shapes
  inline casacore::Matrix<casacore::Int>& chunkShapes() {return chshapes_;};

  // A chunk is good (T) if it contains data
  //  (when averaging, some chunks may have nrows=0)
  inline casacore::Bool goodChunk(casacore::Int ichunk) {return goodChunk_(ichunk); };

  // Is there a reference value for the specified axis?
  // TBD: actually make this axis-dep?
  inline bool hasReferenceValue(PMS::Axis axis) { return (axis==PMS::TIME && cacheReady()); };
  inline double referenceValue(PMS::Axis axis) { return (hasReferenceValue(axis) ? refTime() : 0.0); };
  
  // Report the reference time for this cache (in seconds)
  inline casacore::Double refTime() const { return refTime_p; };

  // Frequency frame in original casacore::MS or requested by user
  inline casacore::MFrequency::Types getFreqFrame() const { return freqFrame_; };

  // Returns which axes have been loaded into the cache, including metadata.
  // Also includes the size (number of points) for each axis (which will
  // eventually be used for a cache manager to let the user know the
  // relative memory use of each axis).
  std::vector<PMS::Axis> loadedAxes() const;

  // Returns true if RA/DEC axes data
  // - for the given parameters - have been loaded
  bool areRaDecAxesLoaded(const DirectionAxisParams &params) const;

  // Access to averaging state in the cache:
  PlotMSAveraging& averaging() { return averaging_; }

  // Access to transformations state in the cache
  PlotMSTransformations& transformations() { return transformations_; }

  // Loads the cache for the given axes and data
  // columns.  IMPORTANT: this method assumes that any currently loaded data is
  // valid for the given VisIter; i.e., if the meta-information or either of
  // the axes are already loaded, then they don't need to be reloaded.  If this
  // is not the case, then clear() should be called BEFORE append().  If a
  // PlotMSCacheThreadHelper object is given, it will be used to report
  // progress information.
  virtual void load(const std::vector<PMS::Axis>& axes,
            const std::vector<PMS::DataColumn>& data,
            const casacore::String& filename,
            const PlotMSSelection& selection,
            const PlotMSAveraging& averaging,
            const PlotMSTransformations& transformations,
            const PlotMSCalibration& calibration,
            /*PlotMSCacheThread**/ThreadCommunication* thread = NULL);
  
  // Clears the cache of all stored values.  This should be called when the
  // underlying casacore::MS or casacore::MS selection is changed, thus invalidating stored data.
  void clear();
  void clearRanges();
  
  // Releases the given axes from the cache.
  void release(const std::vector<PMS::Axis>& axes);
  
  // Set up indexing for the plot
  bool isIndexerInitialized( PMS::Axis iteraxis, casacore::Bool globalXRange,
    casacore::Bool globalYRange, int dataIndex ) const;
  void setUpIndexer(PMS::Axis iteraxis=PMS::SCAN,
    casacore::Bool globalXRange=false, casacore::Bool globalYRange=false, 
    const casacore::String& xconnect = "none", bool timeconnect=false,
    int dataIndex = 0);

  // Some metadata axes not loaded for certain calibration tables
  inline bool hasChan() { return ( chan_[0]->size() > 0); };
  inline bool hasAnt2() { return ( antenna2_[0]->size() > 0); };

  /*** Axis-specific generic gets, per chunk and relative index (from indexer) ***/
  inline casacore::Double getScan(casacore::Int chnk,casacore::Int irel)     { return scan_(chnk);   (void)irel; };
  inline casacore::Double getField(casacore::Int chnk,casacore::Int irel)    { return field_(chnk);  (void)irel; };
  inline casacore::Double getTime(casacore::Int chnk,casacore::Int irel)     { return time_(chnk);  (void)irel; };
  inline casacore::Double getTimeIntr(casacore::Int chnk,casacore::Int irel) { return timeIntr_(chnk);  (void)irel; };
  inline casacore::Double getSpw(casacore::Int chnk,casacore::Int irel)      { return spw_(chnk);  (void)irel; };
  inline casacore::Double getFreq(casacore::Int chnk,casacore::Int irel) { return *(freq_[chnk]->data()+irel); };
  inline casacore::Double getVel(casacore::Int chnk,casacore::Int irel)  { return *(vel_[chnk]->data()+irel); };
  inline casacore::Double getChan(casacore::Int chnk,casacore::Int irel) { return *(chan_[chnk]->data()+irel); };
  inline casacore::Double getCorr(casacore::Int chnk,casacore::Int irel) { return *(corr_[chnk]->data()+irel); };
  inline casacore::Double getAnt1(casacore::Int chnk,casacore::Int irel) { return *(antenna1_[chnk]->data()+irel); };
  inline casacore::Double getAnt2(casacore::Int chnk,casacore::Int irel) { return *(antenna2_[chnk]->data()+irel); };
  inline casacore::Double getBsln(casacore::Int chnk,casacore::Int irel) { return *(baseline_[chnk]->data()+irel); };
  inline casacore::Double getRow(casacore::Int chnk,casacore::Int irel) { return *(row_[chnk]->data()+irel); };
  inline casacore::Double getObsid(casacore::Int chnk,casacore::Int irel) { return *(obsid_[chnk]->data()+irel); };
  // this metadata axis is "loaded" for cal tables so check for empty array
  inline casacore::Double getIntent(casacore::Int chnk,casacore::Int irel) { return ((*intent_[chnk]).empty() ? -1 : *(intent_[chnk]->data()+irel)); };
  inline casacore::Double getFeed1(casacore::Int chnk,casacore::Int irel) { return *(feed1_[chnk]->data()+irel); };
  inline casacore::Double getFeed2(casacore::Int chnk,casacore::Int irel) { return *(feed2_[chnk]->data()+irel); };

  inline casacore::Double getAmp(casacore::Int chnk,casacore::Int irel)  { return *(amp_[chnk]->data()+irel); };
  inline casacore::Double getAmpCorr(casacore::Int chnk,casacore::Int irel)  { return *(ampCorr_[chnk]->data()+irel); };
  inline casacore::Double getAmpModel(casacore::Int chnk,casacore::Int irel)  { return *(ampModel_[chnk]->data()+irel); };
  inline casacore::Double getAmpCorrMod(casacore::Int chnk,casacore::Int irel)  { return *(ampCorrModel_[chnk]->data()+irel); };
  inline casacore::Double getAmpCorrModS(casacore::Int chnk,casacore::Int irel)  { return *(ampCorrModelS_[chnk]->data()+irel); };
  inline casacore::Double getAmpDataMod(casacore::Int chnk,casacore::Int irel)  { return *(ampDataModel_[chnk]->data()+irel); };
  inline casacore::Double getAmpDataModS(casacore::Int chnk,casacore::Int irel)  { return *(ampDataModelS_[chnk]->data()+irel); };
  inline casacore::Double getAmpDataDivMod(casacore::Int chnk,casacore::Int irel)  { return *(ampDataDivModel_[chnk]->data()+irel); };
  inline casacore::Double getAmpDataDivModS(casacore::Int chnk,casacore::Int irel)  { return *(ampDataDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getAmpCorrDivMod(casacore::Int chnk,casacore::Int irel)  { return *(ampCorrDivModel_[chnk]->data()+irel); };
  inline casacore::Double getAmpCorrDivModS(casacore::Int chnk,casacore::Int irel)  { return *(ampCorrDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getAmpFloat(casacore::Int chnk,casacore::Int irel)  { return *(ampFloat_[chnk]->data()+irel); };

  inline casacore::Double getPha(casacore::Int chnk,casacore::Int irel)  { return *(pha_[chnk]->data()+irel); };
  inline casacore::Double getPhaCorr(casacore::Int chnk,casacore::Int irel)  { return *(phaCorr_[chnk]->data()+irel); };
  inline casacore::Double getPhaModel(casacore::Int chnk,casacore::Int irel)  { return *(phaModel_[chnk]->data()+irel); };
  inline casacore::Double getPhaCorrMod(casacore::Int chnk,casacore::Int irel)  { return *(phaCorrModel_[chnk]->data()+irel); };
  inline casacore::Double getPhaCorrModS(casacore::Int chnk,casacore::Int irel)  { return *(phaCorrModelS_[chnk]->data()+irel); };
  inline casacore::Double getPhaDataMod(casacore::Int chnk,casacore::Int irel)  { return *(phaDataModel_[chnk]->data()+irel); };
  inline casacore::Double getPhaDataModS(casacore::Int chnk,casacore::Int irel)  { return *(phaDataModelS_[chnk]->data()+irel); };
  inline casacore::Double getPhaDataDivMod(casacore::Int chnk,casacore::Int irel)  { return *(phaDataDivModel_[chnk]->data()+irel); };
  inline casacore::Double getPhaDataDivModS(casacore::Int chnk,casacore::Int irel)  { return *(phaDataDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getPhaCorrDivMod(casacore::Int chnk,casacore::Int irel)  { return *(phaCorrDivModel_[chnk]->data()+irel); };
  inline casacore::Double getPhaCorrDivModS(casacore::Int chnk,casacore::Int irel)  { return *(phaCorrDivModelS_[chnk]->data()+irel); };

  inline casacore::Double getReal(casacore::Int chnk,casacore::Int irel) { return *(real_[chnk]->data()+irel); };
  inline casacore::Double getRealCorr(casacore::Int chnk,casacore::Int irel)  { return *(realCorr_[chnk]->data()+irel); };
  inline casacore::Double getRealModel(casacore::Int chnk,casacore::Int irel)  { return *(realModel_[chnk]->data()+irel); };
  inline casacore::Double getRealCorrMod(casacore::Int chnk,casacore::Int irel)  { return *(realCorrModel_[chnk]->data()+irel); };
  inline casacore::Double getRealCorrModS(casacore::Int chnk,casacore::Int irel)  { return *(realCorrModelS_[chnk]->data()+irel); };
  inline casacore::Double getRealDataMod(casacore::Int chnk,casacore::Int irel)  { return *(realDataModel_[chnk]->data()+irel); };
  inline casacore::Double getRealDataModS(casacore::Int chnk,casacore::Int irel)  { return *(realDataModelS_[chnk]->data()+irel); };
  inline casacore::Double getRealDataDivMod(casacore::Int chnk,casacore::Int irel)  { return *(realDataDivModel_[chnk]->data()+irel); };
  inline casacore::Double getRealDataDivModS(casacore::Int chnk,casacore::Int irel)  { return *(realDataDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getRealCorrDivMod(casacore::Int chnk,casacore::Int irel)  { return *(realCorrDivModel_[chnk]->data()+irel); };
  inline casacore::Double getRealCorrDivModS(casacore::Int chnk,casacore::Int irel)  { return *(realCorrDivModelS_[chnk]->data()+irel); };

  inline casacore::Double getImag(casacore::Int chnk,casacore::Int irel) { return *(imag_[chnk]->data()+irel); };
  inline casacore::Double getImagCorr(casacore::Int chnk,casacore::Int irel)  { return *(imagCorr_[chnk]->data()+irel); };
  inline casacore::Double getImagModel(casacore::Int chnk,casacore::Int irel)  { return *(imagModel_[chnk]->data()+irel); };
  inline casacore::Double getImagCorrMod(casacore::Int chnk,casacore::Int irel)  { return *(imagCorrModel_[chnk]->data()+irel); };
  inline casacore::Double getImagCorrModS(casacore::Int chnk,casacore::Int irel)  { return *(imagCorrModelS_[chnk]->data()+irel); };
  inline casacore::Double getImagDataMod(casacore::Int chnk,casacore::Int irel)  { return *(imagDataModel_[chnk]->data()+irel); };
  inline casacore::Double getImagDataModS(casacore::Int chnk,casacore::Int irel)  { return *(imagDataModelS_[chnk]->data()+irel); };
  inline casacore::Double getImagDataDivMod(casacore::Int chnk,casacore::Int irel)  { return *(imagDataDivModel_[chnk]->data()+irel); };
  inline casacore::Double getImagDataDivModS(casacore::Int chnk,casacore::Int irel)  { return *(imagDataDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getImagCorrDivMod(casacore::Int chnk,casacore::Int irel)  { return *(imagCorrDivModel_[chnk]->data()+irel); };
  inline casacore::Double getImagCorrDivModS(casacore::Int chnk,casacore::Int irel)  { return *(imagCorrDivModelS_[chnk]->data()+irel); };

  inline casacore::Double getWtxAmp(casacore::Int chnk, casacore::Int irel) { return *(wtxamp_[chnk]->data()+irel); }
  inline casacore::Double getWtxAmpCorr(casacore::Int chnk,casacore::Int irel)  { return *(wtxampCorr_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpModel(casacore::Int chnk,casacore::Int irel)  { return *(wtxampModel_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpCorrMod(casacore::Int chnk,casacore::Int irel)  { return *(wtxampCorrModel_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpCorrModS(casacore::Int chnk,casacore::Int irel)  { return *(wtxampCorrModelS_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpDataMod(casacore::Int chnk,casacore::Int irel)  { return *(wtxampDataModel_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpDataModS(casacore::Int chnk,casacore::Int irel)  { return *(wtxampDataModelS_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpDataDivMod(casacore::Int chnk,casacore::Int irel)  { return *(wtxampDataDivModel_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpDataDivModS(casacore::Int chnk,casacore::Int irel)  { return *(wtxampDataDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpCorrDivMod(casacore::Int chnk,casacore::Int irel)  { return *(wtxampCorrDivModel_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpCorrDivModS(casacore::Int chnk,casacore::Int irel)  { return *(wtxampCorrDivModelS_[chnk]->data()+irel); };
  inline casacore::Double getWtxAmpFloat(casacore::Int chnk,casacore::Int irel)  { return *(wtxampFloat_[chnk]->data()+irel); };

  inline casacore::Double getFlag(casacore::Int chnk,casacore::Int irel) { return *(flag_[chnk]->data()+irel); };
  inline casacore::Double getFlagRow(casacore::Int chnk,casacore::Int irel) { return *(flagrow_[chnk]->data()+irel); };

  inline casacore::Double getWt(casacore::Int chnk,casacore::Int irel) { return *(wt_[chnk]->data()+irel); };
  inline casacore::Double getWtSp(casacore::Int chnk,casacore::Int irel) { return *(wtsp_[chnk]->data()+irel); };
  inline casacore::Double getSigma(casacore::Int chnk,casacore::Int irel) { return *(sigma_[chnk]->data()+irel); };
  inline casacore::Double getSigmaSp(casacore::Int chnk,casacore::Int irel) { return *(sigmasp_[chnk]->data()+irel); };

  inline casacore::Double getUVDist(casacore::Int chnk,casacore::Int irel) { return *(uvdist_[chnk]->data()+irel); };
  inline casacore::Double getUVDistL(casacore::Int chnk,casacore::Int irel) { return *(uvdistL_[chnk]->data()+irel); };
  inline casacore::Double getU(casacore::Int chnk,casacore::Int irel) { return *(u_[chnk]->data()+irel); };
  inline casacore::Double getV(casacore::Int chnk,casacore::Int irel) { return *(v_[chnk]->data()+irel); };
  inline casacore::Double getW(casacore::Int chnk,casacore::Int irel) { return *(w_[chnk]->data()+irel); };
  inline casacore::Double getUwave(casacore::Int chnk,casacore::Int irel) { return *(uwave_[chnk]->data()+irel); };
  inline casacore::Double getVwave(casacore::Int chnk,casacore::Int irel) { return *(vwave_[chnk]->data()+irel); };
  inline casacore::Double getWwave(casacore::Int chnk,casacore::Int irel) { return *(wwave_[chnk]->data()+irel); };

  // These are array-global (one value per chunk)
  inline casacore::Double getAz0(casacore::Int chnk,casacore::Int irel) { return az0_(chnk);  (void)irel; };
  inline casacore::Double getEl0(casacore::Int chnk,casacore::Int irel) { return el0_(chnk);  (void)irel; };
  inline casacore::Double getRadialVelocity0(casacore::Int chnk, casacore::Int irel){ return radialVelocity_(chnk); (void)irel;};
  inline casacore::Double getRHO0(casacore::Int chnk, casacore::Int irel){return rho_(chnk); (void)irel; };
  inline casacore::Double getHA0(casacore::Int chnk,casacore::Int irel) { return ha0_(chnk);  (void)irel; };
  inline casacore::Double getPA0(casacore::Int chnk,casacore::Int irel) { return pa0_(chnk);  (void)irel; };

  // These are antenna-based
  inline casacore::Double getAntenna(casacore::Int chnk,casacore::Int irel) { return *(antenna_[chnk]->data()+irel); };
  inline casacore::Double getAz(casacore::Int chnk,casacore::Int irel)      { return *(az_[chnk]->data()+irel); };
  inline casacore::Double getEl(casacore::Int chnk,casacore::Int irel)             { return *(el_[chnk]->data()+irel); };
  casacore::Double getRa(casacore::Int chnk,casacore::Int irel) {
	  throw AipsError("PlotMS internal error. PlotMsCacheBase::getRa() was called");
	  return chnk + irel; /* return *(ra_[chnk]->data()+irel); */
  };
  inline casacore::Double getDec(casacore::Int chnk,casacore::Int irel) {
		throw AipsError("PlotMS internal error. PlotMsCacheBase::getDec() was called");
	  return chnk + irel; /*(dec_[chnk]->data()+irel); */
  };
  inline casacore::Double getParAng(casacore::Int chnk,casacore::Int irel)  { return *(parang_[chnk]->data()+irel); };

  // These support generic non-complex calibration
  inline casacore::Double getPar(casacore::Int chnk,casacore::Int irel)  { return *(par_[chnk]->data()+irel); };
  inline casacore::Double getSnr(casacore::Int chnk,casacore::Int irel)  { return *(snr_[chnk]->data()+irel); };
  inline casacore::Double getAntPos(casacore::Int chnk,casacore::Int irel)  { return *(antpos_[chnk]->data()+irel); };

  // Curve overlays
  inline casacore::Double getAtm(casacore::Int chnk,casacore::Int irel) { return *(atm_[chnk]->data()+irel); };
  inline casacore::Double getTsky(casacore::Int chnk,casacore::Int irel) { return *(tsky_[chnk]->data()+irel); };

  /* -----------------------------------------------------------------------*/
  /*** Axis-specific generic gets, per chunk (for unit tests) ***/
  // metadata axes
  inline casacore::Int scan(casacore::Int chnk)     { return scan_(chnk);};
  inline casacore::Int field(casacore::Int chnk)    { return field_(chnk); };
  inline casacore::Double time(casacore::Int chnk)     { return time_(chnk); };
  inline casacore::Double timeIntr(casacore::Int chnk) { return timeIntr_(chnk); };
  inline casacore::Int spw(casacore::Int chnk)      { return spw_(chnk); };
  inline casacore::Vector<casacore::Int>& chan(casacore::Int chnk) { return *(chan_[chnk]); };
  inline casacore::Vector<casacore::Double>& freq(casacore::Int chnk) { return *(freq_[chnk]); };
  inline casacore::Vector<casacore::Double>& vel(casacore::Int chnk)  { return *(vel_[chnk]); };
  inline casacore::Vector<casacore::Int>& corr(casacore::Int chnk) { return *(corr_[chnk]); };
  inline casacore::Vector<casacore::Int>& ant1(casacore::Int chnk) { return *(antenna1_[chnk]); };
  inline casacore::Vector<casacore::Int>& ant2(casacore::Int chnk) { return *(antenna2_[chnk]); };
  inline casacore::Vector<casacore::Int>& bsln(casacore::Int chnk) { return *(baseline_[chnk]); };
  inline casacore::Vector<casacore::uInt>& row(casacore::Int chnk) { return *(row_[chnk]); };
  inline casacore::Vector<casacore::Int>& obsid(casacore::Int chnk) { return *(obsid_[chnk]); };
  inline casacore::Vector<casacore::Int>& intent(casacore::Int chnk) { return *(intent_[chnk]); };
  inline casacore::Vector<casacore::Int>& feed1(casacore::Int chnk) { return *(feed1_[chnk]); };
  inline casacore::Vector<casacore::Int>& feed2(casacore::Int chnk) { return *(feed2_[chnk]); };

  // visibility and flag axes (S is for scalar residuals)
  inline casacore::Array<casacore::Float>& ampData(casacore::Int chnk)  { return *(amp_[chnk]); };
  inline casacore::Array<casacore::Float>& ampCorr(casacore::Int chnk)  { return *(ampCorr_[chnk]); };
  inline casacore::Array<casacore::Float>& ampModel(casacore::Int chnk)  { return *(ampModel_[chnk]); };
  inline casacore::Array<casacore::Float>& ampCorrModel(casacore::Int chnk)  { return *(ampCorrModel_[chnk]); };
  inline casacore::Array<casacore::Float>& ampCorrModelS(casacore::Int chnk)  { return *(ampCorrModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& ampDataModel(casacore::Int chnk)  { return *(ampDataModel_[chnk]); };
  inline casacore::Array<casacore::Float>& ampDataModelS(casacore::Int chnk)  { return *(ampDataModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& ampDataDivModel(casacore::Int chnk)  { return *(ampDataDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& ampDataDivModelS(casacore::Int chnk)  { return *(ampDataDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& ampCorrDivModel(casacore::Int chnk)  { return *(ampCorrDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& ampCorrDivModelS(casacore::Int chnk)  { return *(ampCorrDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& ampFloat(casacore::Int chnk)  { return *(ampFloat_[chnk]); };
  inline casacore::Array<casacore::Float>& phaData(casacore::Int chnk)  { return *(pha_[chnk]); };
  inline casacore::Array<casacore::Float>& phaCorr(casacore::Int chnk)  { return *(phaCorr_[chnk]); };
  inline casacore::Array<casacore::Float>& phaModel(casacore::Int chnk)  { return *(phaModel_[chnk]); };
  inline casacore::Array<casacore::Float>& phaCorrModel(casacore::Int chnk)  { return *(phaCorrModel_[chnk]); };
  inline casacore::Array<casacore::Float>& phaCorrModelS(casacore::Int chnk)  { return *(phaCorrModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& phaDataModel(casacore::Int chnk)  { return *(phaDataModel_[chnk]); };
  inline casacore::Array<casacore::Float>& phaDataModelS(casacore::Int chnk)  { return *(phaDataModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& phaDataDivModel(casacore::Int chnk)  { return *(phaDataDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& phaDataDivModelS(casacore::Int chnk)  { return *(phaDataDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& phaCorrDivModel(casacore::Int chnk)  { return *(phaCorrDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& phaCorrDivModelS(casacore::Int chnk)  { return *(phaCorrDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& realData(casacore::Int chnk)  { return *(real_[chnk]); };
  inline casacore::Array<casacore::Float>& realCorr(casacore::Int chnk)  { return *(realCorr_[chnk]); };
  inline casacore::Array<casacore::Float>& realModel(casacore::Int chnk)  { return *(realModel_[chnk]); };
  inline casacore::Array<casacore::Float>& realCorrModel(casacore::Int chnk)  { return *(realCorrModel_[chnk]); };
  inline casacore::Array<casacore::Float>& realCorrModelS(casacore::Int chnk)  { return *(realCorrModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& realDataModel(casacore::Int chnk)  { return *(realDataModel_[chnk]); };
  inline casacore::Array<casacore::Float>& realDataModelS(casacore::Int chnk)  { return *(realDataModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& realDataDivModel(casacore::Int chnk)  { return *(realDataDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& realDataDivModelS(casacore::Int chnk)  { return *(realDataDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& realCorrDivModel(casacore::Int chnk)  { return *(realCorrDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& realCorrDivModelS(casacore::Int chnk)  { return *(realCorrDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& imagData(casacore::Int chnk)  { return *(imag_[chnk]); };
  inline casacore::Array<casacore::Float>& imagCorr(casacore::Int chnk)  { return *(imagCorr_[chnk]); };
  inline casacore::Array<casacore::Float>& imagModel(casacore::Int chnk)  { return *(imagModel_[chnk]); };
  inline casacore::Array<casacore::Float>& imagCorrModel(casacore::Int chnk)  { return *(imagCorrModel_[chnk]); };
  inline casacore::Array<casacore::Float>& imagCorrModelS(casacore::Int chnk)  { return *(imagCorrModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& imagDataModel(casacore::Int chnk)  { return *(imagDataModel_[chnk]); };
  inline casacore::Array<casacore::Float>& imagDataModelS(casacore::Int chnk)  { return *(imagDataModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& imagDataDivModel(casacore::Int chnk)  { return *(imagDataDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& imagDataDivModelS(casacore::Int chnk)  { return *(imagDataDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& imagCorrDivModel(casacore::Int chnk)  { return *(imagCorrDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& imagCorrDivModelS(casacore::Int chnk)  { return *(imagCorrDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampData(casacore::Int chnk)  { return *(wtxamp_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampCorr(casacore::Int chnk)  { return *(wtxampCorr_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampModel(casacore::Int chnk)  { return *(wtxampModel_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampCorrModel(casacore::Int chnk)  { return *(wtxampCorrModel_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampCorrModelS(casacore::Int chnk)  { return *(wtxampCorrModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampDataModel(casacore::Int chnk)  { return *(wtxampDataModel_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampDataModelS(casacore::Int chnk)  { return *(wtxampDataModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampDataDivModel(casacore::Int chnk)  { return *(wtxampDataDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampDataDivModelS(casacore::Int chnk)  { return *(wtxampDataDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampCorrDivModel(casacore::Int chnk)  { return *(wtxampCorrDivModel_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampCorrDivModelS(casacore::Int chnk)  { return *(wtxampCorrDivModelS_[chnk]); };
  inline casacore::Array<casacore::Float>& wtxampFloat(casacore::Int chnk)  { return *(wtxampFloat_[chnk]); };
  inline casacore::Array<casacore::Bool>& flag(casacore::Int chunk) { return *flag_[chunk]; };
  inline casacore::Vector<casacore::Bool>& flagrow(casacore::Int chunk) { return *flagrow_[chunk]; };

  // weight axes
  inline casacore::Array<casacore::Float>& wt(casacore::Int chnk)  { return *(wt_[chnk]); };
  inline casacore::Array<casacore::Float>& wtsp(casacore::Int chnk)  { return *(wtsp_[chnk]); };
  inline casacore::Array<casacore::Float>& sigma(casacore::Int chnk)  { return *(sigma_[chnk]); };
  inline casacore::Array<casacore::Float>& sigmasp(casacore::Int chnk)  { return *(sigmasp_[chnk]); };

  // observational geometry axes
  inline casacore::Vector<casacore::Double>& uVDist(casacore::Int chnk)  { return *(uvdist_[chnk]); };
  inline casacore::Matrix<casacore::Double>& uVDistL(casacore::Int chnk)  { return *(uvdistL_[chnk]); };
  inline casacore::Vector<casacore::Double>& u(casacore::Int chnk)  { return *(u_[chnk]); };
  inline casacore::Vector<casacore::Double>& v(casacore::Int chnk)  { return *(v_[chnk]); };
  inline casacore::Vector<casacore::Double>& w(casacore::Int chnk)  { return *(w_[chnk]); };
  inline casacore::Matrix<casacore::Double>& uWave(casacore::Int chnk)  { return *(uwave_[chnk]); };
  inline casacore::Matrix<casacore::Double>& vWave(casacore::Int chnk)  { return *(vwave_[chnk]); };
  inline casacore::Matrix<casacore::Double>& wWave(casacore::Int chnk)  { return *(wwave_[chnk]); };
  inline casacore::Double az0(casacore::Int chnk)  { return az0_[chnk]; };
  inline casacore::Double el0(casacore::Int chnk)  { return el0_[chnk]; };
  inline casacore::Double ha0(casacore::Int chnk)  { return ha0_[chnk]; };
  inline casacore::Double pa0(casacore::Int chnk)  { return pa0_[chnk]; };
  inline casacore::Vector<casacore::Int>& ant(casacore::Int chnk)  { return *(antenna_[chnk]); };
  inline casacore::Vector<casacore::Double>& az(casacore::Int chnk)  { return *(az_[chnk]); };
  inline casacore::Vector<casacore::Double>& el(casacore::Int chnk)  { return *(el_[chnk]); };
  inline casacore::Vector<casacore::Float>& parAng(casacore::Int chnk)  { return *(parang_[chnk]); };

  // ephemeris axes
  inline casacore::Double radvel(casacore::Int chnk)  { return radialVelocity_[chnk]; };
  inline casacore::Double rho(casacore::Int chnk)  { return rho_[chnk]; };

  // calibration axes
  inline casacore::Array<casacore::Float>& par(casacore::Int chnk)  { return *(par_[chnk]); };
  inline casacore::Array<casacore::Float>& snr(casacore::Int chnk)  { return *(snr_[chnk]); };
  inline casacore::Array<casacore::Float>& antpos(casacore::Int chnk)  { return *(antpos_[chnk]); };

  // curve overlay axes
  inline casacore::Vector<casacore::Double>& atm(casacore::Int chnk)  { return *(atm_[chnk]); };
  inline casacore::Vector<casacore::Double>& tsky(casacore::Int chnk)  { return *(tsky_[chnk]); };

  /* -----------------------------------------------------------------------*/

  // Returns a list of channel numbers that were averaged together in that chunk
  inline casacore::Vector<casacore::Int> getChansPerBin(casacore::Int chnk,casacore::Int irel) { return (*chansPerBin_[chnk])[irel]; };

  casacore::Record locateInfo(int plotIterIndex, 
    const casacore::Vector<PlotRegion>& regions, bool showUnflagged, 
    bool showFlagged, bool selectAll );

  PlotLogMessage* locateRange( int plotIterIndex, 
    const casacore::Vector<PlotRegion> & regions, bool showUnflagged, 
    bool showFlagged);

  PlotLogMessage* flagRange( int plotIterIndex, casa::PlotMSFlagging& flagging,
     const casacore::Vector<PlotRegion>& regions, bool showFlagged);

  //Return a formatted string for time iteration plots giving the time range.
  casacore::String getTimeBounds( int iterValue );
  // Return the time as doubles 
  pair<casacore::Double,casacore::Double> getTimeBounds() const;
  // Return the axes ranges
  pair<casacore::Double,casacore::Double> getXAxisBounds(int index) const;
  pair<casacore::Double,casacore::Double> getYAxisBounds(int index) const;

  inline PMS::DataColumn getXDataColumn() { return currentXData_[0]; };
  inline PMS::DataColumn getYDataColumn(int index) { return currentYData_[index]; };

  void setPlot(PlotMSPlot *plot);

  using RaDecData = casacore::PtrBlock<casacore::Vector<casacore::Double>*>;
  using RaDecMap = std::map<DirectionAxisParams,RaDecData>;

  bool isValidRaDecIndex(int index) const;
  const RaDecData & getRaDataX(int index) const;
  const RaDecData & getRaDataY(int index) const;
  const RaDecData & getDecDataX(int index) const;
  const RaDecData & getDecDataY(int index) const;

  // public log method
  inline void logmesg(const casacore::String& method, 
    const casacore::String& message, int type=PlotLogger::MSG_INFO) 
      { log(method, message, type); };

protected:
    
  // Forbid copy for now
  PlotMSCacheBase(const PlotMSCacheBase&);

  // Resize storage for the number of chunks
  // increaseCache parameter:
  //   false to initialize with empty Arrays before loading cache
  //   true  to increase: copy values then add empty Arrays while loading cache
  void setCache(casacore::Int newnChunk, const std::vector<PMS::Axis>& loadAxes,
    const std::vector<PMS::DataColumn>& loadData, bool increaseCache=false);
  template<typename T> void addArrays(
    casacore::PtrBlock<casacore::Array<T>*>& input, bool increaseCache=false);
  template<typename T> void addMatrices(
    casacore::PtrBlock<casacore::Matrix<T>*>& input, bool increaseCache=false);
  template<typename T> void addVectors(
    casacore::PtrBlock<casacore::Vector<T>*>& input, bool increaseCache=false);

  // Specialized method for loading the cache
  //  (pure virtual: implemented specifically in child classes)
  virtual void loadIt(std::vector<PMS::Axis>& loadAxes,
      std::vector<PMS::DataColumn>& loadData,
      /*PlotMSCacheThread**/ThreadCommunication* thread = NULL)=0;

  virtual void flagToDisk(const PlotMSFlagging& flagging,
    casacore::Vector<casacore::Int>& chunks, 
    casacore::Vector<casacore::Int>& relids,
    casacore::Bool flag,
    PlotMSIndexer* indexer, int dataIndex)=0;
  
  // Clean up the PtrBlocks
  void deleteCache();
  void deleteIndexer();

  // helpers for atm/tsky overlays
  void deleteAtm();
  void printAtmStats(casacore::Int scan);

  virtual bool isEphemeris() {return false;};
  bool isEphemerisAxis( PMS::Axis axis ) const;
  // Set the net axes mask (defines how to collapse flags for the chosen plot axes)
  void setAxesMask(PMS::Axis axis,casacore::Vector<casacore::Bool>& axismask);

  // Return the net axes mask for the currently set plot axes
  casacore::Vector<casacore::Bool> netAxesMask(PMS::Axis xaxis,PMS::Axis yaxis);

  // Derive the plot mask by appropriately collapsing the flags
  void setPlotMask( casacore::Int dataIndex);           // all chunks
  void setPlotMask(casacore::Int dataIndex, casacore::Int chunk);  // per chunk

  // Delete the whole plot mask
  void deletePlotMask();

  // Returns the number of points loaded for the given axis or 0 if not loaded.
  // Only used for PlotMSCacheTab, no longer used in plotms
  //unsigned int nPointsForAxis(PMS::Axis axis) const;
  
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
  
  void logLoad(const casacore::String& message) {
      log(PMS::LOG_ORIGIN_LOAD_CACHE, message, PMS::LOG_EVENT_LOAD_CACHE); }
  void logFlag(const casacore::String& message) {
      log(PMS::LOG_ORIGIN_FLAG, message, PMS::LOG_EVENT_FLAG); }
  // </group>
  
  // Logs the given message from the given method name as the given event type
  // (see PlotLogger).
  void log(const casacore::String& method, const casacore::String& message, int eventType);

  //Return the color lookup index for the chunk.
  int findColorIndex( int chunk, bool initialize );

  // Check access to vector
  template<typename T>
  T checkIndex(int index, const std::vector<T>& v, const std::string &vname) const;

  // Private data
  
  // Parent plotms.
  //  (used only for access to logger, so far)
  PlotMSApp* plotms_;

  // Parent PlotMSPlot, if any
  PlotMSPlot* itsPlot_;

  // An empty indexer (its an empty PlotData object used for initialization)
  PlotMSIndexer* indexer0_;

  // The indexer into the cache
  std::vector<casacore::PtrBlock<PlotMSIndexer*> > indexer_;
  
  // The number of chunks in the cache
  casacore::Int nChunk_;

  // The reference time for this cache, in seconds
  casacore::Double refTime_p;

  // The number of antennas
  casacore::Int nAnt_;

  // Set frame from VI if not specified by user
  // (for VI2::getFrequencies and axis label)
  casacore::MFrequency::Types freqFrame_;

  // Global min/max
  casacore::Double minX_,maxX_,minY_,maxY_;

  // The fundamental meta-data cache
  casacore::Matrix<casacore::Int> chshapes_;
  casacore::Vector<casacore::Bool> goodChunk_;
  casacore::Vector<casacore::Double> time_, timeIntr_;
  casacore::Vector<casacore::Int> field_, spw_, scan_;
  casacore::PtrBlock<casacore::Vector<casacore::uInt>*> row_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> antenna1_, antenna2_, baseline_;
  casacore::PtrBlock<casacore::Vector<casacore::Double>*> uvdist_, u_, v_, w_;
  casacore::PtrBlock<casacore::Matrix<casacore::Double>*> uvdistL_, uwave_, vwave_, wwave_;
  casacore::PtrBlock<casacore::Vector<casacore::Double>*> freq_, vel_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> chan_;
  casacore::PtrBlock<casacore::Array<casacore::Int>*> chansPerBin_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> corr_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> obsid_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> intent_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> feed1_, feed2_;

  // Optional parts of the cache
  casacore::PtrBlock<casacore::Vector<casacore::Float>*> pa_;

  // casacore::Data (the heavy part)
  casacore::PtrBlock<casacore::Array<casacore::Float>*> amp_, 
      ampCorr_, ampModel_, ampCorrModel_, ampCorrModelS_, ampDataModel_, 
      ampDataModelS_, ampDataDivModel_, ampDataDivModelS_, ampCorrDivModel_,
      ampCorrDivModelS_, ampFloat_;
  casacore::PtrBlock<casacore::Array<casacore::Float>*> pha_, 
      phaCorr_, phaModel_, phaCorrModel_, phaCorrModelS_, phaDataModel_, 
      phaDataModelS_, phaDataDivModel_, phaDataDivModelS_, phaCorrDivModel_,
      phaCorrDivModelS_;  // no phase for FLOAT_DATA
  casacore::PtrBlock<casacore::Array<casacore::Float>*> real_, 
      realCorr_, realModel_, realCorrModel_, realCorrModelS_, realDataModel_,
      realDataModelS_, realDataDivModel_, realDataDivModelS_, realCorrDivModel_,
      realCorrDivModelS_;  // use real_ for FLOAT_DATA
  casacore::PtrBlock<casacore::Array<casacore::Float>*> imag_,
      imagCorr_, imagModel_, imagCorrModel_, imagCorrModelS_, imagDataModel_,
      imagDataModelS_, imagDataDivModel_, imagDataDivModelS_, imagCorrDivModel_,
      imagCorrDivModelS_;  // no imag for FLOAT_DATA
  casacore::PtrBlock<casacore::Array<casacore::Float>*> wtxamp_,
      wtxampCorr_, wtxampModel_, wtxampCorrModel_, wtxampCorrModelS_,
      wtxampDataModel_, wtxampDataModelS_, wtxampDataDivModel_, 
      wtxampDataDivModelS_, wtxampCorrDivModel_, wtxampCorrDivModelS_,
      wtxampFloat_;

  casacore::PtrBlock<casacore::Array<casacore::Bool>*> flag_;
  casacore::PtrBlock<casacore::Vector<casacore::Bool>*> flagrow_;
  
  casacore::PtrBlock<casacore::Array<casacore::Float>*> wt_,wtsp_;
  casacore::PtrBlock<casacore::Array<casacore::Float>*> sigma_,sigmasp_;

  casacore::PtrBlock<casacore::Vector<casacore::Float>*> parang_;
  casacore::PtrBlock<casacore::Vector<casacore::Int>*> antenna_;
  casacore::PtrBlock<casacore::Vector<casacore::Double>*> az_,el_;
  casacore::PtrBlock<casacore::Vector<casacore::Double>*> ra_,dec_;
  std::map<DirectionAxisParams,RaDecData> raMap_,decMap_;

  casacore::Vector<casacore::Double> radialVelocity_, rho_;
  casacore::Vector<casacore::Double> az0_,el0_,ha0_,pa0_;

  casacore::PtrBlock<casacore::Vector<casacore::Double>*> atm_, tsky_;

  // for cal tables
  casacore::PtrBlock<casacore::Array<casacore::Float>*> par_, snr_;
  casacore::PtrBlock<casacore::Array<casacore::Float>*> antpos_; 

  // Current setup/state.
  bool dataLoaded_;
  bool userCanceled_;
  std::vector<PMS::Axis> currentX_;
  std::vector<PMS::Axis> currentY_;
  std::vector<PMS::DataColumn> currentXData_;
  std::vector<PMS::DataColumn> currentYData_;
  std::vector<PMS::CoordSystem> currentXFrame_;
  std::vector<PMS::CoordSystem> currentYFrame_;
  std::vector<PMS::InterpMethod> currentXInterp_;
  std::vector<PMS::InterpMethod> currentYInterp_;
  std::vector<PMS::CoordSystem> xyFrame_;
  std::vector<PMS::InterpMethod> xyInterp_;
  std::vector<PMS::CoordSystem> loadXYFrame_;
  std::vector<PMS::InterpMethod> loadXYInterp_;
  decltype(raMap_)::mapped_type *loadRa_;
  decltype(decMap_)::mapped_type *loadDec_;
  map<PMS::Axis, bool> loadedAxes_;
  map<PMS::Axis, casacore::Record> loadedAxesData_;
  //map<PMS::Axis, std::set<PMS::DataColumn>> loadedAxesData_;
  map<PMS::Axis, bool> pendingLoadAxes_;

  // Global ranges (unflagged and flagged, per indexer)
  casacore::Vector<casacore::Double> xminG_, xmaxG_, yminG_, ymaxG_,
      xflminG_, xflmaxG_, yflminG_, yflmaxG_;

  // A copy of the casacore::Data parameters 
  casacore::String filename_;
  PlotMSSelection selection_;
  PlotMSAveraging averaging_;
  PlotMSTransformations transformations_;
  PlotMSCalibration calibration_;

  // Axes mask
  std::vector<casacore::Vector<casacore::Bool> > netAxesMask_;

  // collapsed flag mask for plotting
  std::vector<casacore::PtrBlock<casacore::Array<casacore::Bool>* > > plmask_;

  // meta info for locate output
  casacore::Vector<casacore::String> antnames_;
  casacore::Vector<casacore::String> stanames_;
  casacore::Vector<casacore::String> antstanames_;
  casacore::Vector<casacore::String> fldnames_;
  casacore::Vector<casacore::String> intentnames_;
  casacore::Array<casacore::Double> positions_;

  PMS::Axis iterAxis;
  bool ephemerisInitialized;
  ::QVector<double> uniqueTimes;

  // The calibration type (casacore::Table subType)
  casacore::String calType_;
  // polarization selection is ratio ("/")
  bool polnRatio_;

  // Page header items
  PageHeaderCache pageHeaderCache_;

  // For atm/tsky overlays
  PlotMSAtm* plotmsAtm_;


private:
  void _updateAntennaMask( casacore::Int a, casacore::Vector<casacore::Bool>& antMask, const casacore::Vector<casacore::Int> selectedAntennas );
  bool axisIsValid(PMS::Axis axis, const PlotMSAveraging& averaging);

};
typedef casacore::CountedPtr<PlotMSCacheBase> PlotMSCacheBasePtr;


}

#endif /* PLOTMSCACHEBASE_H_ */
