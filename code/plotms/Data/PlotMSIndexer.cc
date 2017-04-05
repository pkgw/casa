//# PlotMSIndexer.cc: Cache indexer for plotms.
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
#include <plotms/Data/PlotMSIndexer.h>

#include <casa/Quanta/MVTime.h>
#include <casa/Utilities/Sort.h>
#include <casa/OS/Timer.h>
#include <plotms/PlotMS/PlotMS.h>
#include <tables/Tables/Table.h>
#include <measures/Measures/Stokes.h>
//#include <QtCore/qmath.h>
#include <QDebug>

using namespace casacore;
namespace casa {

PlotMSIndexer::PlotMSIndexer():
		  plotmscache_(NULL),
		  currChunk_(0),
		  irel_(0),
		  lasti_(-1),
		  nPoints_(),
		  nCumulative_(),
		  nSegment_(0),currSeg_(0),
		  nSegPoints_(),nCumulPoints_(),cacheChunk_(),cacheOffset_(),
		  currentX_(PMS::SCAN),
		  currentY_(PMS::SCAN),
		  indexerReady_(false),
		  icorrmax_(),
		  ichanmax_(),
		  ibslnmax_(),
		  idatamax_(),
		  nperchan_(),
		  nperbsln_(),
		  nperant_(),
		  ichanbslnmax_(),
		  iantmax_(),
		  xmin_(DBL_MAX),
		  ymin_(DBL_MAX),
		  xflmin_(DBL_MAX),
		  yflmin_(DBL_MAX),
		  xmax_(-DBL_MAX),
		  ymax_(-DBL_MAX),
		  xflmax_(-DBL_MAX),
		  yflmax_(-DBL_MAX),
		  sizeMasked_(0),
		  sizeUnMasked_(0),
		  globalXMinMax_(false),
		  globalYMinMax_(false),
		  iterate_(false),
		  iterAxis_(PMS::NONE),
		  iterValue_(-999),
		  itsColorize_(false),
		  itsColorizeAxis_(PMS::DEFAULT_COLOR_AXIS),
		  self(const_cast<PlotMSIndexer*>(this))
{
	dataIndex = 0;
	}

PlotMSIndexer::PlotMSIndexer(PlotMSCacheBase* parent, PMS::Axis xAxis,
        PMS::DataColumn xData, PMS::Axis yAxis, PMS::DataColumn yData,
        int index ):
		  plotmscache_(parent),
		  currChunk_(0),
		  irel_(0),
		  lasti_(-1),
		  nPoints_(),
		  nCumulative_(),
		  nSegment_(0),currSeg_(0),
		  nSegPoints_(),nCumulPoints_(),cacheChunk_(),cacheOffset_(),
		  currentX_(xAxis),
		  currentY_(yAxis),
		  currentXdata_(xData),
		  currentYdata_(yData),
		  indexerReady_(false),
		  icorrmax_(),
		  ichanmax_(),
		  ibslnmax_(),
		  idatamax_(),
		  nperchan_(),
		  nperbsln_(),
		  nperant_(),
		  ichanbslnmax_(),
		  iantmax_(),
		  xmin_(DBL_MAX),
		  ymin_(DBL_MAX),
		  xflmin_(DBL_MAX),
		  yflmin_(DBL_MAX),
		  xmax_(-DBL_MAX),
		  ymax_(-DBL_MAX),
		  xflmax_(-DBL_MAX),
		  yflmax_(-DBL_MAX),
		  sizeMasked_(0),
		  sizeUnMasked_(0),
		  globalXMinMax_(false),
		  globalYMinMax_(false),
		  iterate_(false),
		  iterAxis_(PMS::SCAN),
		  iterValue_(-999),
		  itsColorize_(false),
		  itsColorizeAxis_(PMS::DEFAULT_COLOR_AXIS),
		  self(const_cast<PlotMSIndexer*>(this))
{
	dataIndex = index;
	setUpIndexing();
}

PlotMSIndexer::PlotMSIndexer(PlotMSCacheBase* parent,
        PMS::Axis xAxis, PMS::DataColumn xDataColumn, 
        PMS::Axis yAxis, PMS::DataColumn yDataColumn, 
        PMS::Axis iterAxis, Int iterValue, int index ):
		plotmscache_(parent),
		currChunk_(0),
		irel_(0),
		lasti_(-1),
		nPoints_(),
		nCumulative_(),
		nSegment_(0),currSeg_(0),
		nSegPoints_(),nCumulPoints_(),cacheChunk_(),cacheOffset_(),
		currentX_(xAxis),
		currentY_(yAxis),
        currentXdata_(xDataColumn),
        currentYdata_(yDataColumn),
		indexerReady_(false),
		icorrmax_(),
		ichanmax_(),
		ibslnmax_(),
		idatamax_(),
		nperchan_(),
		nperbsln_(),
		nperant_(),
		ichanbslnmax_(),
		iantmax_(),
		xmin_(DBL_MAX),
		ymin_(DBL_MAX),
		xflmin_(DBL_MAX),
		yflmin_(DBL_MAX),
		xmax_(-DBL_MAX),
		ymax_(-DBL_MAX),
		xflmax_(-DBL_MAX),
		yflmax_(-DBL_MAX),
		sizeMasked_(0),
		sizeUnMasked_(0),
		globalXMinMax_(false),
		globalYMinMax_(false),
		iterate_(iterAxis!=PMS::NONE),
		iterAxis_(iterAxis),
		iterValue_(iterValue),
		itsColorize_(false),
		itsColorizeAxis_(PMS::DEFAULT_COLOR_AXIS),
		self(const_cast<PlotMSIndexer*>(this))
{ 
	dataIndex = index;
	setUpIndexing();
}

PlotMSIndexer::~PlotMSIndexer() {} // anything?

unsigned int PlotMSIndexer::size() const { 
	//  return (nChunk()>0 ? nCumulative_(nChunk()-1) : 0);
	return (nSegment_>0 ? nCumulPoints_(nSegment_-1) : 0);
}

double PlotMSIndexer::xAt(unsigned int i) const {
	setChunk(i);  // sets chunk and relative index in chunk

	double x= (plotmscache_->*getXFromCache_)(currChunk_,
			(self->*XIndexer_)(currChunk_,irel_));
	return x;
}
double PlotMSIndexer::yAt(unsigned int i) const {
	setChunk(i);  // sets chunk and relative index in chunk
	return (plotmscache_->*getYFromCache_)(currChunk_,
			(self->*YIndexer_)(currChunk_,irel_));
}
void PlotMSIndexer::xAndYAt(unsigned int index, 
		double& x, double& y) const {
	setChunk(index);  // sets chunk and relative index in chunk
	x=(plotmscache_->*getXFromCache_)(currChunk_,
			(self->*XIndexer_)(currChunk_,irel_));
	y=(plotmscache_->*getYFromCache_)(currChunk_,
			(self->*YIndexer_)(currChunk_,irel_));
}

bool PlotMSIndexer::minsMaxes(double& xMin, double& xMax, 
		double& yMin, double& yMax) {

	if (this->size()<1) return false;

	// return the collective (flagged/unflagged) min/max
	// X:
	if (globalXMinMax_ || (sizeMasked()==0 && sizeUnmasked()==0)) {
		// calculate global
		xMin=min(plotmscache_->xminG_,plotmscache_->xflminG_);
		xMax=max(plotmscache_->xmaxG_,plotmscache_->xflmaxG_);
	}
	else {
		xMin=min(xmin_,xflmin_);
		xMax=max(xmax_,xflmax_);
	}

	// Y:
	if (globalYMinMax_ || (sizeMasked()==0 && sizeUnmasked()==0)) {
		// calculate global
		yMin=min(plotmscache_->yminG_,plotmscache_->yflminG_);
		yMax=max(plotmscache_->ymaxG_,plotmscache_->yflmaxG_);
	}
	else {
		yMin=min(ymin_,yflmin_);
		yMax=max(ymax_,yflmax_);
	}
	return true;
}

bool PlotMSIndexer::maskedAt( unsigned int index) const {
	setChunk(index);
	return !(*(plotmscache_->plmask_[dataIndex][currChunk_]->data()+irel_));
}
void PlotMSIndexer::xyAndMaskAt(unsigned int index,
		double& x, double& y,
		bool& mask) const {
	setChunk(index);
	x=(plotmscache_->*getXFromCache_)(currChunk_,
			(self->*XIndexer_)(currChunk_,irel_));
	y=(plotmscache_->*getYFromCache_)(currChunk_,
			(self->*YIndexer_)(currChunk_,irel_));
	mask=!(*(plotmscache_->plmask_[dataIndex][currChunk_]->data()+irel_));
}

bool PlotMSIndexer::maskedMinsMaxes(double& xMin, double& xMax, 
		double& yMin, double& yMax) {

	if (this->size()<1) return false;

	// return the collective (flagged) min/max
	// X:
	if (globalXMinMax_ || sizeMasked()==0) {
		// Use globals from the cache
		xMin=plotmscache_->xflminG_;
		xMax=plotmscache_->xflmaxG_;
	}
	else {
		// get local ones
		xMin=xflmin_;
		xMax=xflmax_;
	}
	// Y:
	if (globalYMinMax_ || sizeMasked()==0) {
		// Use globals from the cache
		yMin=plotmscache_->yflminG_;
		yMax=plotmscache_->yflmaxG_;
	}
	else {
		// use local ones
		yMin=yflmin_;
		yMax=yflmax_;
	}
	return true;
}



bool PlotMSIndexer::maskedMinsMaxesRaw(double& xMin, double& xMax, 
		double& yMin, double& yMax) {

	if (this->size()<1) return false;

	xMin=xflmin_;
	xMax=xflmax_;
	yMin=yflmin_;
	yMax=yflmax_;
	return true;
}

bool PlotMSIndexer::unmaskedMinsMaxes(double& xMin, double& xMax, 
		double& yMin, double& yMax) {

	if (this->size()<1) return false;

	// return the collective (unflagged) min/max
	// X:
	if (globalXMinMax_ || sizeUnmasked()==0 ) {
		// Use globals from the cache
		xMin=plotmscache_->xminG_;
		xMax=plotmscache_->xmaxG_;
	}
	else {
		// get local ones
		xMin=xmin_;
		xMax=xmax_;
	}

	// Y:
	if (globalYMinMax_ || sizeUnmasked()==0 ) {
		// Use globals from the cache
		yMin=plotmscache_->yminG_;
		yMax=plotmscache_->ymaxG_;
	}
	else {
		// get local ones
		yMin=ymin_;
		yMax=ymax_;
	}
	return true;
}


bool PlotMSIndexer::unmaskedMinsMaxesRaw(double& xMin, double& xMax, 
		double& yMin, double& yMax) {

	if (this->size()<1) return false;

	xMin=xmin_;
	xMax=xmax_;
	yMin=ymin_;
	yMax=ymax_;
	return true;
}



unsigned int PlotMSIndexer::numBins() const {
	// TODO
	return PMS::COLORS_LIST().size();
}

unsigned int PlotMSIndexer::binAt(unsigned int i) const {
	unsigned int binValue = 0;
	if(itsColorize_) {
		setChunk(i);
		unsigned int val = (unsigned int)(plotmscache_->*getColFromCache_)(currChunk_,
				(self->*ColIndexer_)(currChunk_,irel_));

		if ( itsColorizeAxis_ != PMS::TIME ){
			binValue = val % numBins();
		}
		else {
			if ( plotmscache_->averaging_.time() ){
				double timeInterval= plotmscache_->averaging_.timeValue();
				double baseTime = plotmscache_->getTime( 0, 0 );
				Double time = plotmscache_->getTime(currChunk_, 0);
				Double timeDiff = time - baseTime;
				int timeIndex = static_cast<int>( timeDiff / timeInterval );
				binValue = timeIndex % numBins();
			}
			else {
				int timeIndex = 0;
				if ( currChunk_ == 0 ){
					timeIndex = plotmscache_->findColorIndex( currChunk_, false );
				}
				else {
					timeIndex = plotmscache_->findColorIndex( currChunk_, false );
				}
				binValue = timeIndex % numBins();

			}
		}
	}
	return binValue;
}

bool PlotMSIndexer::isBinned() const {
	// TODO
	return itsColorize_;
}

bool PlotMSIndexer::colorize(bool doColorize, PMS::Axis colorizeAxis) {
	// TODO
	bool changed = (doColorize != itsColorize_) ||
			(doColorize && colorizeAxis != itsColorizeAxis_);
	itsColorize_ = doColorize;
	itsColorizeAxis_ = colorizeAxis;

	if (itsColorize_) {
        // None of the coloraxis options have a datacolumn so doesn't matter
		setMethod(getColFromCache_,itsColorizeAxis_, PMS::DEFAULT_DATACOLUMN);
		setIndexer(ColIndexer_,itsColorizeAxis_);
	}

	//cout << "COLORIZE!! " << boolalpha << itsColorize_ << " " << PMS::axis(itsColorizeAxis_) << endl;

	return changed;
}


void PlotMSIndexer::setUpIndexing() {

	// Forbid antenna-based/baseline-based combination plots, for now
	//  (e.g., data vs. _antenna-based_ elevation)
	if (plotmscache_->netAxesMask_[dataIndex](2)&&plotmscache_->netAxesMask_[dataIndex](3))
		throw(AipsError("Cannot yet support antenna-based and baseline-based data in same plot."));

	// Refer to the chunk shape matrix in the cache
	Matrix<Int>& chsh(plotmscache_->chunkShapes());

	icorrmax_.reference(chsh.row(0));
	ichanmax_.reference(chsh.row(1));
	ibslnmax_.reference(chsh.row(2));
	iantmax_.reference(chsh.row(3));

	idatamax_.resize(nChunk());
	idatamax_ = chsh.row(0);
	idatamax_ *= chsh.row(1);
	idatamax_ *= chsh.row(2);

	ichanbslnmax_.resize(nChunk());
	ichanbslnmax_ = chsh.row(1);
	ichanbslnmax_ *= chsh.row(2);

	nperchan_.resize(nChunk());
	nperchan_.set(1);
	if (plotmscache_->netAxesMask_[dataIndex](0)) nperchan_ *= chsh.row(0);

	nperbsln_.resize(nChunk());
	nperbsln_.set(1);
	if (plotmscache_->netAxesMask_[dataIndex](0)) nperbsln_ *= chsh.row(0);
	if (plotmscache_->netAxesMask_[dataIndex](1)) nperbsln_ *= chsh.row(1);

	nperant_.reference(nperbsln_);


  /*cout << "ichanbslnmax_ = " << ichanbslnmax_ << endl;
  cout << "nperchan_     = " << nperchan_ << endl;
  cout << "nperbsln_     = " << nperbsln_ << endl;
  cout << "nperant_      = " << nperant_ << endl;

  cout << "...done." << endl << "Set methods..." << flush;
*/

	// Set up method pointers for the chosen axes
	setMethod(getXFromCache_, currentX_, currentXdata_);
	setMethod(getYFromCache_, currentY_, currentYdata_);

	// And the indexers
	setIndexer(XIndexer_, currentX_);
	setIndexer(YIndexer_, currentY_);

	//  cout << "done." << endl;

	// And the mask collapsers
	//  setCollapser(collapseXMask_,currentX_);
	//  setCollapser(collapseYMask_,currentY_);

	// Count up the total number of points we will plot
	//   (keep a cumualtive running total)

	//  cout << "*>*>*>*>*>*> iterate_=" << boolalpha << iterate_ << " Axis: " << PMS::axis(iterAxis_) << " iterValue_ = " << iterValue_ << endl;

	//  cout << "Count points..." << flush;

	// Count data segments in this iteration
	switch (iterAxis_) {
	case PMS::ANTENNA: {
		nSegment_=0;
		for (Int ich=0; ich<nChunk(); ++ich)
			// only check for non-empty chunks
			if (plotmscache_->goodChunk(ich)) {
				for (Int ibl=0; ibl<chsh(2,ich); ++ibl)
					if ( (*(plotmscache_->antenna1_[ich]->data()+ibl) == iterValue_) ||
							(*(plotmscache_->antenna2_[ich]->data()+ibl) == iterValue_) )
						++nSegment_;
			}
		break;
	}
	case PMS::TIME : {
		nSegment_ = 0;
		for ( Int ich=0; ich<nChunk(); ++ich ){
			if (plotmscache_->goodChunk(ich)){
				if (plotmscache_->getTime(ich, 0) == plotmscache_->getTime(iterValue_,0)){
					++nSegment_;
				}
			}
		}
		break;
	}
    case PMS::CORR: {
		nSegment_=0;
		for (Int ich=0; ich<nChunk(); ++ich)
			// only check for non-empty chunks
			if (plotmscache_->goodChunk(ich)) {
				for (Int icorr=0; icorr<chsh(0,ich); ++icorr)
					if (*(plotmscache_->corr_[ich]->data()+icorr) == iterValue_) 
						nSegment_ += ichanbslnmax_[ich];
			}
		break;
	}
	default:
		// most iteration axis have nChunk segments
		nSegment_ = nChunk();
	}

	// Size up and initialize the indexing helpers and counters
	nSegPoints_.resize(nSegment_,0);
	nSegPoints_.set(0);
	nCumulPoints_.resize(nSegment_,0);
	nCumulPoints_.set(0);
	cacheChunk_.resize(nSegment_,0);
	cacheChunk_.set(0);
	cacheOffset_.resize(nSegment_,0);
	cacheOffset_.set(0);

	// Refer to simple axes, for convenience below
	Vector<Int> iterAxisRef;
	switch (iterAxis_) {
	case PMS::SCAN:
		iterAxisRef.reference(plotmscache_->scan_);
		break;
	case PMS::SPW:
		iterAxisRef.reference(plotmscache_->spw_);
		break;
	case PMS::FIELD:
		iterAxisRef.reference(plotmscache_->field_);
		break;
	default:
		break;
	}

	// Count per segment
	Int iseg(-1);
	Vector<Bool>& nAM(plotmscache_->netAxesMask_[dataIndex]);
	double timeInterval = 1;
	bool averagingTime = plotmscache_->averaging_.time();
	if ( averagingTime ){
		timeInterval = plotmscache_->averaging_.timeValue();
	}
	double iterTime = plotmscache_->time_[iterValue_];

	for (Int ic=0; ic<nChunk(); ++ic) {

		// skip this chunk if empty
		if (!plotmscache_->goodChunk(ic)){
			continue;
		}

		switch (iterAxis_) {
		case PMS::NONE:
		case PMS::FIELD:
		case PMS::SPW:
		case PMS::SCAN: {
			if (iterAxis_==PMS::NONE || iterAxisRef(ic)==iterValue_) {
				++iseg;
				cacheChunk_(iseg)=ic;
				cacheOffset_(iseg)=0;
				nSegPoints_(iseg)=((nAM(0) ? icorrmax_(ic) : 1)*
						(nAM(1) ? ichanmax_(ic) : 1)*
						(nAM(2) ? ibslnmax_(ic) : 1)*
						(nAM(3) ? iantmax_(ic) : 1));
			}
			break;
		}
		case PMS::TIME: {

			bool assignValues = false;
			if ( averagingTime ){
				double icTime = plotmscache_->time_[ic];
				if ( icTime >= iterTime && icTime < iterTime + timeInterval ){
					assignValues = true;
				}
			}
			else {
				if ( iterTime  == plotmscache_->time_[ic] ){
					assignValues = true;
				}
			}

			if ( assignValues ){
				++iseg;
				cacheChunk_(iseg) = ic;
				cacheOffset_(iseg) = 0;
				nSegPoints_(iseg) = ((nAM(0) ? icorrmax_(ic) : 1)*
						(nAM(1) ? ichanmax_(ic) : 1)*
						(nAM(2) ? ibslnmax_(ic) : 1)*
						(nAM(3) ? iantmax_(ic) : 1));
			}
			break;
		}
		case PMS::BASELINE: {
			Int ibsln=0, nBsln=chsh(2,ic);
			while (ibsln<nBsln && (*(plotmscache_->baseline_[ic]->data()+ibsln)!=iterValue_))
				++ibsln;
			if (ibsln<nBsln) {
				// Found the baseline in this chunk
				++iseg;
				cacheChunk_(iseg)=ic;
				cacheOffset_(iseg)=ibsln*Int(nAM(2))*nperbsln_(ic); // non-zero only if there is a baseline axis
				nSegPoints_(iseg)=nperbsln_(ic);
			}
			break;
		}
		case PMS::ANTENNA: {
			Int nBsln=chsh(2,ic);
			for (Int ibsln=0;ibsln<nBsln;++ibsln) {
				if (*(plotmscache_->antenna1_[ic]->data()+ibsln)==iterValue_ ||
						*(plotmscache_->antenna2_[ic]->data()+ibsln)==iterValue_) {
					// found antenna for this iteration
					++iseg;
					cacheChunk_(iseg)=ic;
					cacheOffset_(iseg)=ibsln*Int(nAM(2))*nperbsln_(ic);
					nSegPoints_(iseg)=nperbsln_(ic);
				}
			}
			break;
		}
        case PMS::CORR: {
            Int nCorr = chsh(0,ic);
            for (Int icorr=0; icorr<nCorr; ++icorr) {
                if (*(plotmscache_->corr_[ic]->data()+icorr) == iterValue_) {
                    for (Int nPoints=0; nPoints<ichanbslnmax_[ic]; ++nPoints) {
                        ++iseg;
                        cacheChunk_(iseg)  = ic;
                        cacheOffset_(iseg) = icorr + nPoints*nCorr;
					    nSegPoints_(iseg)  = 1;
                    }
                }
            }
            break;
        }
        default:
			// shouldn't reach here...
			throw(AipsError("Unsupported iteration axis: "+PMS::axis(iterAxis_)));
			break;
		}
	}

	// Contract nSegment_ if we aren't using them all
	if (iseg+1 < nSegment_) {
		nSegment_ = iseg+1;

		// Cope with no segments found
		//  (this happens when all data is flagged, time-averaging is on,
		//    and iteration is off, in v3.2)
		if (nSegment_==0) nSegment_=1;

		// (w/ copy)
		nSegPoints_.resize(nSegment_,true);
		nCumulPoints_.resize(nSegment_,true);
		cacheChunk_.resize(nSegment_,true);
		cacheOffset_.resize(nSegment_,true);
	}

	// Fill cumulative counter
	nCumulPoints_(0) = nSegPoints_(0);
	for (Int iseg=1; iseg<nSegment_; ++iseg)
		nCumulPoints_(iseg) = nCumulPoints_(iseg-1) + nSegPoints_(iseg);

	nPoints_.reference(nSegPoints_);
	nCumulative_.reference(nCumulPoints_);

	//  cout << "done." << endl;

	// Compute the nominal plot ranges
	computeRanges();

	// The indexer is now ready for plotting
	indexerReady_ = true;

}

bool PlotMSIndexer::isGlobalXRange() const {
	return globalXMinMax_;
}

bool PlotMSIndexer::isGlobalYRange() const {
	return globalYMinMax_;
}

void PlotMSIndexer::setGlobalMinMax(Bool globalX, Bool globalY ) {
    globalXMinMax_=globalX;
    globalYMinMax_=globalY;
};

void PlotMSIndexer::setChunk(uInt i) const {

	// NB: this method assumes that i>=lasti, for now

	if (i==lasti_)
		// already found this one on previous call (e.g., for mask
		//   or the other axis), so change nothing
		return;

	// reset to the first chunk if very first or earlier point requested
	if (i==0 || i<lasti_) currSeg_=0;

	// Bump at segment boundaries
	while (i > (nCumulPoints_(currSeg_)-1)) ++currSeg_;

	// TBD:  Back up to a previous chunk?
	//  while (i < (nCumulative_(currChunk_)) && currChunk_>0) --currChunk_;


	// Found current Indexer _segment_, so set current _cache_ chunk:
	currChunk_=cacheChunk_(currSeg_);

	// Calculate the offset into the current chunk
	if (currSeg_>0)
		irel_=Int(i-nCumulPoints_(currSeg_-1));
	else
		irel_=Int(i);

	// Offset into the cache (e.g., non-zero for baseline iteration)
	irel_+=cacheOffset_(currSeg_);

	// Remember this i next time around
	lasti_=i;

}

void PlotMSIndexer::setMethod(CacheMemPtr& getmethod,PMS::Axis axis,
        PMS::DataColumn datacol) {

	// Set axis-specific get methods
	switch(axis) {
    // Metadata
	case PMS::SCAN:
		getmethod = &PlotMSCacheBase::getScan;
		break;
	case PMS::FIELD:
		getmethod = &PlotMSCacheBase::getField;
		break;
	case PMS::TIME:
		getmethod = &PlotMSCacheBase::getTime;
		break;
	case PMS::TIME_INTERVAL:
		getmethod = &PlotMSCacheBase::getTimeIntr;
		break;
	case PMS::SPW:
		getmethod = &PlotMSCacheBase::getSpw;
		break;
	case PMS::CHANNEL:
		getmethod = &PlotMSCacheBase::getChan;
		break;
	case PMS::FREQUENCY:
		getmethod = &PlotMSCacheBase::getFreq;
		break;
	case PMS::VELOCITY:
		getmethod = &PlotMSCacheBase::getVel;
		break;
	case PMS::CORR:
		getmethod = &PlotMSCacheBase::getCorr;
		break;
	case PMS::ANTENNA1:
		getmethod = &PlotMSCacheBase::getAnt1;
		break;
	case PMS::ANTENNA2:
		getmethod = &PlotMSCacheBase::getAnt2;
		break;
	case PMS::BASELINE:
		getmethod = &PlotMSCacheBase::getBsln;
		break;
	case PMS::ROW:
		getmethod = &PlotMSCacheBase::getRow;
		break;
	case PMS::OBSERVATION:
		getmethod = &PlotMSCacheBase::getObsid;
		break;
	case PMS::INTENT:
		getmethod = &PlotMSCacheBase::getIntent;
		break;
	case PMS::FEED1:
		getmethod = &PlotMSCacheBase::getFeed1;
		break;
	case PMS::FEED2:
		getmethod = &PlotMSCacheBase::getFeed2;
		break;

    // Data
	case PMS::AMP: {
        switch(datacol) {
            case PMS::DATA:
                getmethod = &PlotMSCacheBase::getAmp;
                break;
            case PMS::CORRECTED:
                getmethod = &PlotMSCacheBase::getAmpCorr;
                break;
            case PMS::MODEL:
                getmethod = &PlotMSCacheBase::getAmpModel;
                break;
            case PMS::CORRMODEL:
                getmethod = &PlotMSCacheBase::getAmpCorrMod;
                break;
            case PMS::DATAMODEL:
                getmethod = &PlotMSCacheBase::getAmpDataMod;
                break;
            case PMS::DATA_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getAmpDataDivMod;
                break;
            case PMS::CORRECTED_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getAmpCorrDivMod;
                break;
            case PMS::FLOAT_DATA:
                getmethod = &PlotMSCacheBase::getAmpFloat;
                break;
            }
        }
		break;
	case PMS::PHASE: {
        switch(datacol) {
            case PMS::DATA:
		        getmethod = &PlotMSCacheBase::getPha;
                break;
            case PMS::CORRECTED:
                getmethod = &PlotMSCacheBase::getPhaCorr;
                break;
            case PMS::MODEL:
                getmethod = &PlotMSCacheBase::getPhaModel;
                break;
            case PMS::CORRMODEL:
                getmethod = &PlotMSCacheBase::getPhaCorrMod;
                break;
            case PMS::DATAMODEL:
                getmethod = &PlotMSCacheBase::getPhaDataMod;
                break;
            case PMS::DATA_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getPhaDataDivMod;
                break;
            case PMS::CORRECTED_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getPhaCorrDivMod;
                break;
            case PMS::FLOAT_DATA:
                break;
            }
        }
		break;
	case PMS::REAL: {
        switch(datacol) {
            case PMS::DATA:
		        getmethod = &PlotMSCacheBase::getReal;
                break;
            case PMS::CORRECTED:
                getmethod = &PlotMSCacheBase::getRealCorr;
                break;
            case PMS::MODEL:
                getmethod = &PlotMSCacheBase::getRealModel;
                break;
            case PMS::CORRMODEL:
                getmethod = &PlotMSCacheBase::getRealCorrMod;
                break;
            case PMS::DATAMODEL:
                getmethod = &PlotMSCacheBase::getRealDataMod;
                break;
            case PMS::DATA_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getRealDataDivMod;
                break;
            case PMS::CORRECTED_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getRealCorrDivMod;
                break;
            case PMS::FLOAT_DATA:
                getmethod = &PlotMSCacheBase::getReal;
                break;
            }
        }
		break;
	case PMS::IMAG: {
        switch(datacol) {
            case PMS::DATA:
		        getmethod = &PlotMSCacheBase::getImag;
                break;
            case PMS::CORRECTED:
                getmethod = &PlotMSCacheBase::getImagCorr;
                break;
            case PMS::MODEL:
                getmethod = &PlotMSCacheBase::getImagModel;
                break;
            case PMS::CORRMODEL:
                getmethod = &PlotMSCacheBase::getImagCorrMod;
                break;
            case PMS::DATAMODEL:
                getmethod = &PlotMSCacheBase::getImagDataMod;
                break;
            case PMS::DATA_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getImagDataDivMod;
                break;
            case PMS::CORRECTED_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getImagCorrDivMod;
                break;
            case PMS::FLOAT_DATA:
                break;
            }
        }
		break;
	case PMS::WTxAMP: {
        switch(datacol) {
            case PMS::DATA:
		        getmethod = &PlotMSCacheBase::getWtxAmp;
                break;
            case PMS::CORRECTED:
                getmethod = &PlotMSCacheBase::getWtxAmpCorr;
                break;
            case PMS::MODEL:
                getmethod = &PlotMSCacheBase::getWtxAmpModel;
                break;
            case PMS::CORRMODEL:
                getmethod = &PlotMSCacheBase::getWtxAmpCorrMod;
                break;
            case PMS::DATAMODEL:
                getmethod = &PlotMSCacheBase::getWtxAmpDataMod;
                break;
            case PMS::DATA_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getWtxAmpDataDivMod;
                break;
            case PMS::CORRECTED_DIVIDE_MODEL:
                getmethod = &PlotMSCacheBase::getWtxAmpCorrDivMod;
                break;
            case PMS::FLOAT_DATA:
                getmethod = &PlotMSCacheBase::getWtxAmpFloat;
                break;
            }
        }
		break;

	case PMS::WT:
		getmethod = &PlotMSCacheBase::getWt;
		break;
	case PMS::WTSP:
		getmethod = &PlotMSCacheBase::getWtSp;
		break;
	case PMS::SIGMA:
		getmethod = &PlotMSCacheBase::getSigma;
		break;
	case PMS::SIGMASP:
		getmethod = &PlotMSCacheBase::getSigmaSp;
		break;

	case PMS::FLAG:
		getmethod = &PlotMSCacheBase::getFlag;
		break;
	case PMS::FLAG_ROW:
		getmethod = &PlotMSCacheBase::getFlagRow;
		break;

	case PMS::UVDIST:
		getmethod = &PlotMSCacheBase::getUVDist;
		break;
	case PMS::UVDIST_L:
		getmethod = &PlotMSCacheBase::getUVDistL;
		break;
	case PMS::U:
		getmethod = &PlotMSCacheBase::getU;
		break;
	case PMS::V:
		getmethod = &PlotMSCacheBase::getV;
		break;
	case PMS::W:
		getmethod = &PlotMSCacheBase::getW;
		break;
	case PMS::UWAVE:
		getmethod = &PlotMSCacheBase::getUwave;
		break;
	case PMS::VWAVE:
		getmethod = &PlotMSCacheBase::getVwave;
		break;
	case PMS::WWAVE:
		getmethod = &PlotMSCacheBase::getWwave;
		break;

	case PMS::AZ0:
		getmethod = &PlotMSCacheBase::getAz0;
		break;
	case PMS::EL0:
		getmethod = &PlotMSCacheBase::getEl0;
		break;
	case PMS::HA0:
		getmethod = &PlotMSCacheBase::getHA0;
		break;
	case PMS::PA0:
		getmethod = &PlotMSCacheBase::getPA0;
		break;
	case PMS::ANTENNA:
		getmethod = &PlotMSCacheBase::getAntenna;
		break;
	case PMS::AZIMUTH:
		getmethod = &PlotMSCacheBase::getAz;
		break;
	case PMS::ELEVATION:
		getmethod = &PlotMSCacheBase::getEl;
		break;
	case PMS::PARANG:
		getmethod = &PlotMSCacheBase::getParAng;
		break;

    // Calibration tables
	case PMS::GAMP:
		getmethod = &PlotMSCacheBase::getAmp;
		break;
	case PMS::GPHASE:
		getmethod = &PlotMSCacheBase::getPha;
		break;
	case PMS::GREAL:
		getmethod = &PlotMSCacheBase::getReal;
		break;
	case PMS::GIMAG:
		getmethod = &PlotMSCacheBase::getImag;
		break;
	case PMS::DELAY:
	case PMS::SWP:
	case PMS::TSYS:
	case PMS::OPAC:
	case PMS::TEC:
		getmethod = &PlotMSCacheBase::getPar;
		break;
	case PMS::SNR:
		getmethod = &PlotMSCacheBase::getSnr;
		break;
	case PMS::RADIAL_VELOCITY:
		getmethod = &PlotMSCacheBase::getRadialVelocity0;
		break;
	case PMS::RHO:
		getmethod = &PlotMSCacheBase::getRHO0;
		break;
	default:
		throw(AipsError("Can't find get method for "+PMS::axis(axis)+"."));
		break;
	}

}

void PlotMSIndexer::setIndexer(IndexerMethPtr& indexmethod,PMS::Axis axis) {

	// Set axis-specific indexing method
	switch(axis) {

	// Degenerate axes (no corr-,chan-,bsln,antenna-dependence)
	case PMS::SCAN:
	case PMS::FIELD:
	case PMS::TIME:
	case PMS::TIME_INTERVAL:
	case PMS::SPW:
	case PMS::OBSERVATION:
	case PMS::INTENT:
	case PMS::FEED1:
	case PMS::FEED2:
		indexmethod = &PlotMSIndexer::getIndex0000;
		break;

		// corr-dep
	case PMS::CORR:
		indexmethod = &PlotMSIndexer::getIndex1000;
		break;

		// corr-,bsln-dep
	case PMS::WT:
	case PMS::SIGMA:
		indexmethod = &PlotMSIndexer::getIndex1010;
		break;

		// corr-,chan-,bsln-dep
	case PMS::AMP:
	case PMS::PHASE:
	case PMS::REAL:
	case PMS::IMAG:
	case PMS::FLAG:
	case PMS::GAMP:
	case PMS::GPHASE:
	case PMS::GREAL:
	case PMS::GIMAG:
	case PMS::DELAY:
	case PMS::SWP:
	case PMS::TSYS:
	case PMS::OPAC:
	case PMS::SNR:
	case PMS::TEC:
	case PMS::WTxAMP:
	case PMS::WTSP:
	case PMS::SIGMASP:
		indexmethod = &PlotMSIndexer::getIndex1110;
		break;

		// chan-dep
	case PMS::FREQUENCY:
	case PMS::VELOCITY:
	case PMS::CHANNEL:
		indexmethod = &PlotMSIndexer::getIndex0100;
		break;

		// chan-,bsln-dep
	case PMS::UVDIST_L:
	case PMS::UWAVE:
	case PMS::VWAVE:
	case PMS::WWAVE:
		indexmethod = &PlotMSIndexer::getIndex0110;
		break;

		// bsln-dep
	case PMS::ROW:
	case PMS::ANTENNA1:
	case PMS::ANTENNA2:
	case PMS::BASELINE:
	case PMS::UVDIST:
	case PMS::U:
	case PMS::V:
	case PMS::W:
	case PMS::FLAG_ROW:
		indexmethod = &PlotMSIndexer::getIndex0010;
		break;

		// chunk-dep geometry
	case PMS::AZ0:
	case PMS::EL0:
	case PMS::RADIAL_VELOCITY:
	case PMS::RHO:
	case PMS::HA0:
	case PMS::PA0:
		indexmethod = &PlotMSIndexer::getIndex0000;
		break;

		// antenna-dep
	case PMS::ANTENNA:
	case PMS::AZIMUTH:
	case PMS::ELEVATION:
	case PMS::PARANG:
		indexmethod = &PlotMSIndexer::getIndex0001;
		break;

	default:
		throw(AipsError("Help! No index method available!"));
	}

}

/* not needed?  (gmoellen 2011Mar15)
void PlotMSIndexer::setCollapser(CollapseMethPtr& collmethod,PMS::Axis axis) {

  // Set axis-specific mask collapsing method
  switch(axis) {

    // Degenerate axes (no corr-,chan-,bsln,antenna-dependence)
  case PMS::SCAN:
  case PMS::FIELD:
  case PMS::TIME:
  case PMS::TIME_INTERVAL:
  case PMS::SPW: {
    collmethod = &PlotMSIndexer::collapseMask0000;
    break;
  }
    // corr-dep
  case PMS::CORR:
    collmethod = &PlotMSIndexer::collapseMask1000;
    break;

    // corr-,bsln-dep
  case PMS::WT:
    collmethod = &PlotMSIndexer::collapseMask1010;
    break;

    // corr-,chan-,bsln-dep
  case PMS::AMP:
  case PMS::PHASE:
  case PMS::REAL:
  case PMS::IMAG:
  case PMS::FLAG:
  case PMS::WTSP:
    collmethod = &PlotMSIndexer::collapseMask1110;
    break;

    // chan-dep
  case PMS::FREQUENCY:
  case PMS::CHANNEL:
    collmethod = &PlotMSIndexer::collapseMask0100;
    break;

    // chan-,bsln-dep
  case PMS::UVDIST_L:
  case PMS::UWAVE:
  case PMS::VWAVE:
  case PMS::WWAVE:
    collmethod = &PlotMSIndexer::collapseMask0110;
    break;

    // bsln-dep
  case PMS::ROW:
  case PMS::ANTENNA1:
  case PMS::ANTENNA2:
  case PMS::BASELINE:
  case PMS::UVDIST:
  case PMS::U:
  case PMS::V:
  case PMS::W:
  case PMS::FLAG_ROW:
    collmethod = &PlotMSIndexer::collapseMask0010;
    break;

    // antenna-dep
  case PMS::ANTENNA:
  case PMS::AZIMUTH:
  case PMS::ELEVATION:
  case PMS::PARANG:
    collmethod = &PlotMSIndexer::collapseMask0001;
    break;

  default:
    throw(AipsError("Help! No index method available!"));
  }

}
 */

Record PlotMSIndexer::getPointMetaData(Int i) {
	Double thisx, thisy;
	xAndYAt(i, thisx, thisy);
	// Collect meta data
	Int ichan = getIndex0100(currChunk_, irel_);
	Int chan = Int(plotmscache_->getChan(currChunk_,ichan));
	Int scan = Int(plotmscache_->getScan(currChunk_,0));
	Int field = Int(plotmscache_->getField(currChunk_,0));
	Int ant1 = Int(plotmscache_->getAnt1(currChunk_,
			getIndex0010(currChunk_,irel_)));
	Int ant2 = Int(plotmscache_->getAnt2(currChunk_,
			getIndex0010(currChunk_,irel_)));
	String ant1name;
	String ant2name;
	if(ant1 == -1) {
		ant1name = "*";
	} else {
		ant1name = plotmscache_->antstanames_(ant1);
	}
	if(ant2 == -1) {
		ant2name = "*";
	} else {
		ant2name = plotmscache_->antstanames_(ant2);
	}
	Double time = plotmscache_->getTime(currChunk_, 0);
	Int spw = Int(plotmscache_->getSpw(currChunk_, 0));
	Double freq = plotmscache_->getFreq(currChunk_, ichan);
	Int icorr = Int(plotmscache_->getCorr(currChunk_,getIndex1000(currChunk_,irel_)));
	String corr = plotmscache_->polname(icorr);
	Int obsId = Int(plotmscache_->getObsid(currChunk_, 0));

	Int offset = (currChunk_ > 0 ? (nCumulative_(currChunk_-1)+irel_) : irel_);
	// Collate meta data
	Record r;
	r.define("chan", chan);
	r.define("scan", scan);
	r.define("field", field);
	r.define("time", time);
	r.define("ant1", ant1);
	r.define("ant2", ant2);
	r.define("ant1name", ant1name);
	r.define("ant2name", ant2name);
	r.define("time", time);
	r.define("spw", spw);
	r.define("freq", freq);
	r.define("corr", corr);
	r.define("x", thisx);
	r.define("y", thisy);
	r.define("obsid", obsId);
	r.define("offset", offset);
	r.define("currchunk", currChunk_);
	r.define("irel", irel_);
	return r;
}

Record PlotMSIndexer::locateInfo(const Vector<PlotRegion>& regions,
		Bool showUnflagged, Bool showFlagged,
		Bool selectAll) {
	int nFound = 0, n = size();
	Int nFoundMasked = 0, nFoundUnmasked = 0;
	Bool m = false;
	Record result;
	result.define("xaxis", PMS::axis(currentX_));
	result.define("yaxis", PMS::axis(currentY_));
	for(Int i = 0; i < n; ++i) {
		m = maskedAt(i);
		// Skip point if it is not displayed
		if((!m && !showUnflagged) || (m && !showFlagged)) {
			continue;
		}
		for(uInt j = 0; j < regions.size(); ++j) {
			Double thisx, thisy;
			xAndYAt(i, thisx, thisy);
			// If a point falls inside a bounding region...
			if(thisx > regions[j].left() && thisx < regions[j].right() &&
					thisy > regions[j].bottom() && thisy < regions[j].top()) {
				// Stuff it in the result Record
				ostringstream ss;
				ss << nFound;
				result.defineRecord(String(ss.str()), getPointMetaData(i));
				// Increment point counts
				++nFound;
				if(m) { ++nFoundMasked; }
				else  { ++nFoundUnmasked; }
				break;
			}
		}
		// If there are no selections, grab all points that are displayed
		if(selectAll && (regions.size() == 0) &&
				((m && showFlagged) || (!m && showUnflagged))) {
			ostringstream ss;
			ss << nFound;
			result.defineRecord(String(ss.str()), getPointMetaData(i));
			++nFound;
			if(m) { ++nFoundMasked; }
			else  { ++nFoundUnmasked; }
		}
	}
	return result;
}


PlotLogMessage* PlotMSIndexer::locateRange(const Vector<PlotRegion>& regions,
		Bool showUnflagged, Bool showFlagged) {

	Timer locatetimer;
	locatetimer.mark();

	Double thisx, thisy;
	stringstream ss;
	Int nFound = 0, n = size();
	Int nFoundMasked(0),nFoundUnmasked(0);

	Bool m(false);
	for(Int i = 0; i < n; i++) {

		m=maskedAt(i);

		// Only locate if point is visible
		if ( (!m && showUnflagged) || (m && showFlagged) ) {

			xAndYAt(i, thisx, thisy);

			for(uInt j = 0; j < regions.size(); j++) {
				if (thisx > regions[j].left() && thisx < regions[j].right() &&
						thisy > regions[j].bottom() && thisy < regions[j].top()) {
					nFound++;
					(m ? ++nFoundMasked : ++nFoundUnmasked);
					// only report first 1000, so logger isn't overloaded
					if (nFound<1001) {
						reportMeta(thisx, thisy, m, ss);
						ss << '\n';
					}
					break;
				}
			}
		}
	}

	if (nFound>1000)
		ss << "NB: Only first 1000 points reported above." << '\n';

	ss << "Found " << nFound << " points (";
	if (showUnflagged) ss << nFoundUnmasked << " unflagged" << (showFlagged ? ", " : "");
	if (showFlagged) ss << nFoundMasked << " flagged";
	ss << ") among " << n << " in "
			<< locatetimer.all_usec()/1.0e6 << "s.";
    if (plotConjugates())
        ss << "\nNote: only points in the original MS can be located, not the plotted conjugate points.";

	return new PlotLogMessage(PMS::LOG_ORIGIN,PMS::LOG_ORIGIN_LOCATE,ss.str(),PMS::LOG_EVENT_LOCATE);
}


void PlotMSIndexer::reportMeta(Double x, Double y, Bool masked,stringstream& ss) {

	// This method assumes currChunk_ and irel_ already set correctly!
    // Note some values not set (-1) in some cal tables
	Bool showindices(true);

	ss << "Scan=";
    Int scan= Int(plotmscache_->getScan(currChunk_,0));
    if (scan < 0)
        ss << "* ";
    else
        ss << scan << " ";

	ss << "Field=";
	Int fld=Int(plotmscache_->getField(currChunk_,0));
	if (fld<0)
		ss << "* ";
	else {
		ss << plotmscache_->fldnames_(fld);
		if (showindices) ss << "[" << fld << "]";
        ss << " ";
	}

    ss << "Time=";
    Double timeVal = plotmscache_->getTime(currChunk_,0);
    if (timeVal == 0.0)
        ss << "* ";
    else
	    ss << MVTime(timeVal/C::day).string(MVTime::YMD,10) << " ";

	Int ant1=Int( plotmscache_->getAnt1(currChunk_,getIndex0010(currChunk_,irel_)) );
	Int ant2=Int( plotmscache_->getAnt2(currChunk_,getIndex0010(currChunk_,irel_)) );
    if (ant2 < 0)
	    ss << "ANT1=";
    else
	    ss << "BL=";
	// Antenna Names
	if (!plotmscache_->netAxesMask_[dataIndex](2) || ant1<0)
		ss << "*";
	else
		ss << plotmscache_->antstanames_(ant1);
	if (!plotmscache_->netAxesMask_[dataIndex](2))
		ss << " & * ";
    else if (ant1==ant2)
		ss << " && " << plotmscache_->antstanames_(ant2);
	else if (ant2>=0)
		ss << " & " << plotmscache_->antstanames_(ant2);
	// Antenna indices
	if (showindices) {
		ss << " [";
		if (!plotmscache_->netAxesMask_[dataIndex](2) || ant1<0)
			ss << "*";
		else
			ss << ant1;
		if (!plotmscache_->netAxesMask_[dataIndex](2))
			ss << "&*";
        else if (ant1==ant2)
			ss << "&&" << ant2;
		else if (ant2>=0)
			ss << "&" << ant2;
		ss << "]";
	}
	ss << " ";

	ss << "Spw=";
	Int spw=Int(plotmscache_->getSpw(currChunk_,0));
	if (spw<0)
		ss << "* ";
	else
		ss << spw << " ";

	ss << "Chan=";
	Int ichan=getIndex0100(currChunk_,irel_);
    PlotMSAveraging& pmsave(plotmscache_->averaging());
	if (plotmscache_->netAxesMask_[dataIndex](1)) {
		if (pmsave.channel() && pmsave.channelValue()>1) {
            Vector<Int> chansPerBin = plotmscache_->getChansPerBin(currChunk_, ichan);
			ss << "<" << chansPerBin[0] << "~" << chansPerBin[chansPerBin.size()-1] << ">";
		}
		else
			ss << Int(plotmscache_->getChan(currChunk_,ichan));
	}
	else
		ss << "*";
	ss << " ";

	if (plotmscache_->netAxesMask_[dataIndex](1)) {
		if (pmsave.channel() && pmsave.channelValue()>1) {
	        ss << "Avg Freq=";
        } else {
	        ss << "Freq=";
        }
		ss << plotmscache_->getFreq(currChunk_,ichan) << " ";
    } else {
		ss << "Freq=*        ";
    }

    if (plotmscache_->cacheType()==PlotMSCacheBase::CAL)
	    ss << "Poln=";
    else
	    ss << "Corr=";
	if (plotmscache_->netAxesMask_[dataIndex](0))
		ss << plotmscache_->polname(Int(plotmscache_->getCorr(currChunk_,getIndex1000(currChunk_,irel_))));
	else
		ss << "*";
	ss << " ";

	ss << "X=" << x << " ";
	ss << "Y="  << y;
	ss << ( masked ? " F " : " ");
	ss << "Observation=" << plotmscache_->getObsid(currChunk_,0) << " ";
}

PlotLogMessage* PlotMSIndexer::flagRange(const PlotMSFlagging& flagging,
		const Vector<PlotRegion>& regions, Bool flag) {
	Timer flagtimer;
	flagtimer.mark();

	// List of flags
	Vector<Int> flagchunk(1000,-1),flagindex(1000,-1);

	Double thisx, thisy;
	stringstream ss;
	Int nFound = 0, n = size(), flsz;

	for(Int i = 0; i < n; i++) {

		// The following sets currChunk_ and irel_ (as needed below)
		if ((!maskedAt(i) && flag) ||    // not yet flagged and we are flagging
				(maskedAt(i) && !flag) ) {   // already flagged and we are unflagging

			xAndYAt(i, thisx, thisy);

			for(uInt j = 0; j < regions.size(); j++) {
				if(thisx > regions[j].left() && thisx < regions[j].right() &&
						thisy > regions[j].bottom() && thisy < regions[j].top()) {
					nFound++;

					// The following assumes currChunk_ and irel_ are properly set...
					flagInCache(flagging, flag);

					// Record this flags indices so we can apply to MS (VisSet) below
					flsz = flagchunk.nelements();
					if(flsz < nFound) {
						// Add 50% more space (w/ copy!)
						flagchunk.resize(Int(flsz * 1.5), true);
						flagindex.resize(Int(flsz * 1.5), true);
					}
					flagchunk(nFound - 1) = currChunk_;
					flagindex(nFound - 1) = irel_;

					//  	reportMeta(thisx, thisy, ss);
					//  ss << '\n';
				}
			}
		}
	}

	//cout << "Found " << nFound << " points to " << (flag ? "flag." : "unflag.") << endl;

	// Apply (un)flags only if some found
	if (nFound > 0) {
		// Refresh the plot mask to reflect newly flagged data
		//  TBD: only do chunks that need it!
		plotmscache_->setPlotMask(dataIndex);

		//    cout << "Finished in-memory flagging." << endl;

		// shrink flag list to correct size
		if(flagchunk.nelements() > uInt(nFound)) {
			flagchunk.resize(nFound, true);
			flagindex.resize(nFound, true);
		}

		//    cout << "flagchunk = " << flagchunk << endl;
		//    cout << "flagindex = " << flagindex << endl;

		// Set the flags in the MS
		plotmscache_->flagToDisk(flagging, flagchunk, flagindex, flag, this, dataIndex);


		// Recompute ranges
		computeRanges();

	}


	ss << (flag ? "FLAGGED " : "UNFLAGGED ") << nFound
			<< " points among " << n << " in "
			<< flagtimer.all_usec()/1.0e6 << "s.";
    if (plotConjugates())
        ss << "\nNote: only points in the original MS can be flagged, not the plotted conjugate points.";

	return new PlotLogMessage(PMS::LOG_ORIGIN,
			flag ? PMS::LOG_ORIGIN_FLAG : PMS::LOG_ORIGIN_UNFLAG,
					ss.str(),
					flag ? PMS::LOG_EVENT_FLAG : PMS::LOG_EVENT_UNFLAG);

}

String PlotMSIndexer::iterValue() {

	String itername(PMS::axis(iterAxis_));

	switch (iterAxis_) {
	case PMS::SCAN:
	case PMS::SPW:
		return String::toString(iterValue_);
		break;
	case PMS::FIELD:
		return plotmscache_->fldnames_(iterValue_);
		break;
	case PMS::TIME:{
		return plotmscache_->getTimeBounds( iterValue_);
		break;
	}
	case PMS::BASELINE: {
		maskedAt(0);  // sets currChunk_ and irel_ for first point so we can get ant indices
		Int ant1=Int(plotmscache_->getAnt1(currChunk_,getIndex0010(currChunk_,irel_)));
		Int ant2=Int(plotmscache_->getAnt2(currChunk_,getIndex0010(currChunk_,irel_)));
		String label;
		label += (ant1>-1 ? plotmscache_->antstanames_(ant1) : "*")+" & ";
		label += (ant2>-1 ? plotmscache_->antstanames_(ant2) : "*");
        // CAS-4239 add baseline length to plot title
        String bsnLen = ((ant1>-1 && ant2>-1) ? String::format("_%.0fm", computeBaselineLength(ant1, ant2)) : "_*m");
        label += bsnLen;
		return label;
		break;
	}
	case PMS::ANTENNA:
		return plotmscache_->antstanames_(iterValue_);
		break;
	case PMS::CORR:
        return plotmscache_->polname(iterValue_);
        break;
	default:
		return String("");
		//    throw(AipsError("Unsupported iteration axis: "+PMS::axis(iterAxis_)));
		break;
	}

	return String("");
}

Double PlotMSIndexer::computeBaselineLength(Int ant1, Int ant2) {
    Vector<Double> ant1pos = plotmscache_->positions_[ant1];
    Vector<Double> ant2pos = plotmscache_->positions_[ant2];
    Double length = sqrt(pow((ant1pos[0] - ant2pos[0]), 2) + 
                         pow((ant1pos[1] - ant2pos[1]), 2) +
                         pow((ant1pos[2] - ant2pos[2]), 2));
    return length;
}

String PlotMSIndexer::iterLabel() {
	String itername(PMS::axis(iterAxis_));
    if ((plotmscache_->cacheType()==PlotMSCacheBase::CAL) &&
        (itername=="Corr"))
        itername="Poln";
	String iterVal = iterValue();
	String iterLabel = itername + ": " +iterVal;
	if (iterAxis_ == PMS::TIME ) {
		iterLabel = ": "+iterVal;
	}
	return iterLabel;
}

String PlotMSIndexer::fileLabel(){
	String itername(PMS::axis(iterAxis_));
    if ((plotmscache_->cacheType()==PlotMSCacheBase::CAL) &&
        (itername=="Corr"))
        itername="Poln";
	String iterVal = iterValue();
	String iterLabel = itername +iterVal;
	if (iterAxis_ == PMS::TIME ) {
		iterLabel = "Time" + iterVal;
	}
	return iterLabel;
}


void PlotMSIndexer::flagInCache(const PlotMSFlagging& flagging,Bool flag) {

	Slice corr,chan,bsln;

	// Set flag range on correlation axis:
	Int icorr(0);
	if (plotmscache_->netAxesMask_[dataIndex](0) && !flagging.corrAll()) {
		// specific correlation
		icorr=getIndex1000(currChunk_,irel_); // (irel_%icorrmax_(currChunk_));
		corr=Slice(icorr,1,1);
	}
	else
		// All correlations
		corr=Slice(0,plotmscache_->chunkShapes()(0,currChunk_),1);

	// Set Flag range on channel axis:
	Int ichan(-1);
	if (plotmscache_->netAxesMask_[dataIndex](1) && !flagging.channel()) {
		// specific channel
		ichan=getIndex0100(currChunk_,irel_); // (irel_%icorrmax_(currChunk_));  //Int(getChan());
		/* ....old way require convert from chan value to channel index...
    if (averaging_.channel()) {
      // correct to the _in-cache_ channel _index_ 
      Int dch=Int(averaging_.channelValue());
      if (dch>1) {
	ichan/=dch;
	ichan-=((*chan_[currChunk_])(0)/dch);
      }
    }
    else
      // when not averaging, maybe first channel isn't 0th
      ichan-=(*chan_[currChunk_])(0);  
		 */

		chan=Slice(ichan,1,1);
	}
	else
		// All channels
		chan=Slice(0,plotmscache_->chunkShapes()(1,currChunk_),1);


	// Set Flag range on baseline axis:
	Int ibsln(-1);
	if (plotmscache_->netAxesMask_[dataIndex](2)) {
		// specific correlation
		ibsln=getIndex0010(currChunk_,irel_);   //(irel_/nperbsln_(currChunk_))%ibslnmax_(currChunk_);
		bsln=Slice(ibsln,1,1);
	}
	else
		// All baselines
		bsln=Slice(0,plotmscache_->chunkShapes()(2,currChunk_),1);

	//  cout << "Bsln, chan, corr: " << ibsln << " " << ichan << " " << icorr << endl;

	// Set the sliced flag
	Cube<Bool> flagcube(plotmscache_->flag(currChunk_));
	flagcube(corr,chan,bsln)=flag;

	// unset flagrow when unflagging (if present in cache)
	if (!flag) {
		Vector<Bool> flagrow(plotmscache_->flagrow(currChunk_));
		if (flagrow.nelements()>0)
			flagrow(bsln)=false;
	}
}

/* These may not ever be needed? (gmoellen 2011March15)

void PlotMSIndexer::collapseMask0000(Int ch,Array<Bool>& collmask) {
  collmask.resize(IPosition(1,1));
  collmask(IPosition(1,0)) = (ntrue( *(plotmscache_->plmask_[ch]) ) > uInt(0));
}

void PlotMSIndexer::collapseMask1000(Int ch,Array<Bool>& collmask) {
  collmask.resize();
  // collapse on chan, row
  collmask = operator>(partialNTrue(*(plotmscache_->plmask_[ch]),IPosition(2,1,2)),uInt(0));
}

void PlotMSIndexer::collapseMask0100(Int ch,Array<Bool>& collmask) {
  collmask.resize();
  // collapse on corr, row
  collmask = operator>(partialNTrue(*(plotmscache_->plmask_[ch]),IPosition(2,0,2)),uInt(0));
}

void PlotMSIndexer::collapseMask0010(Int ch,Array<Bool>& collmask) {
  collmask.resize();
  // collapse on corr, chan
  collmask = operator>(partialNTrue(*(plotmscache_->plmask_[ch]),IPosition(2,0,1)),uInt(0));
}

void PlotMSIndexer::collapseMask0110(Int ch,Array<Bool>& collmask) {
  collmask.resize();
  // collapse on corr
  collmask = operator>(partialNTrue(*(plotmscache_->plmask_[ch]),IPosition(1,0)),uInt(0));
}

void PlotMSIndexer::collapseMask1010(Int ch,Array<Bool>& collmask) {
  collmask.resize();
  // collapse on chan
  collmask = operator>(partialNTrue(*(plotmscache_->plmask_[ch]),IPosition(1,1)),uInt(0));
}

void PlotMSIndexer::collapseMask1110(Int ch,Array<Bool>& collmask) {
  // just reference plmask directly
  collmask.reference(plotmscache_->*plmask_[ch]);
}

void PlotMSIndexer::collapseMask0001(Int ch,Array<Bool>& collmask) {
  //  TBD: generate antenna-based mask from baseline-based flags
  IPosition nsh(3,1,1,iantmax_(ch));
  collmask.resize(nsh);
  collmask.set(true);   
}
 */

void PlotMSIndexer::computeRanges() {

	// Initialize limits
	xmin_=ymin_=xflmin_=yflmin_=DBL_MAX;
	xmax_=ymax_=xflmax_=yflmax_=-DBL_MAX;

	// We will count up flagged/unflagged here
	sizeMasked_=sizeUnMasked_=0;

	// Loop over all points to detected min/max
	for (uInt i=0;i<size();++i) {

		Double x,y;
		Bool m;
		xyAndMaskAt(i,x,y,m);

		// CAS-8019 nan>ymax_ caused error in autorange
		if ( !m ) {
			++sizeUnMasked_;
			if (!isNaN(x)) {
			    xmin_ = min(xmin_,x);
			    xmax_ = max(xmax_,x);
			}
			if (!isNaN(y)) {
			    ymin_ = min(ymin_,y);
			    ymax_ = max(ymax_,y);
			}
		}
		else {
			++sizeMasked_;
			if (!isNaN(x)) {
			    xflmin_ = min(xflmin_,x);
			    xflmax_ = max(xflmax_,x);
			}
			if (!isNaN(y)) {
			    yflmin_ = min(yflmin_,y);
			    yflmax_ = max(yflmax_,y);
			}
		}
	}
}

void PlotMSIndexer::log(const String& method, const String& message,
		int /*eventType*/) {

	cout << method << ": " << message << endl;

	//  plotms_->getLogger()->postMessage(PMS::LOG_ORIGIN, method, message, eventType);}

}

}

