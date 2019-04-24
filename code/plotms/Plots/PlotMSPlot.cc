//# PlotMSPlot.cc: High level plot concept across potentially multiple objects.
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
#include <plotms/Plots/PlotMSPlot.h>

#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Data/CacheFactory.h>
#include <plotms/Data/PlotMSCacheBase.h>
#include <plotms/Data/MSCache.h>
#include <plotms/Data/CalCache.h>
#include <QDebug>

using namespace casacore;
namespace casa {

////////////////////////////
// PLOTMSPLOT DEFINITIONS //
////////////////////////////

// Static //

PlotMSPlotParameters PlotMSPlot::makeParameters(PlotMSApp* plotms) {
	PlotMSPlotParameters p(plotms->getPlotFactory());
	makeParameters(p, plotms);
	return p;
}

void PlotMSPlot::makeParameters(PlotMSPlotParameters& params, PlotMSApp* /*plotms*/) {
	// Add data parameters if needed.
	if(params.typedGroup<PMS_PP_MSData>() == NULL)
		params.setGroup<PMS_PP_MSData>();

	// Add cache parameters if needed
	if(params.typedGroup<PMS_PP_Cache>() == NULL)
		params.setGroup<PMS_PP_Cache>();

	// Add axes parameters if needed.
	if(params.typedGroup<PMS_PP_Axes>() == NULL)
		params.setGroup<PMS_PP_Axes>();

	// Add canvas parameters if needed.
	if(params.typedGroup<PMS_PP_Canvas>() == NULL)
		params.setGroup<PMS_PP_Canvas>();

	// Add display parameters if needed.
	if(params.typedGroup<PMS_PP_Display>() == NULL)
		params.setGroup<PMS_PP_Display>();

	// Add page header parameters if needed.
	if(params.typedGroup<PMS_PP_PageHeader>() == NULL)
		params.setGroup<PMS_PP_PageHeader>();

	// Add iteration parameters if needed.
	if(params.typedGroup<PMS_PP_Iteration>() == NULL)
		params.setGroup<PMS_PP_Iteration>();
}

const uInt PlotMSPlot::PIXEL_THRESHOLD = 1000000;
const uInt PlotMSPlot::MEDIUM_THRESHOLD = 10000;
const uInt PlotMSPlot::LARGE_THRESHOLD = 1000;
const uInt PlotMSPlot::XLARGE_THRESHOLD = 50;

// Constructors/Destructors //

PlotMSPlot::PlotMSPlot(PlotMSApp* parent) : 
	itsParent_(parent),
	itsFactory_(parent->getPlotFactory()),
	itsParams_(itsFactory_),
	itsCache_(NULL),
	iter_(0),
	iterStep_(1) {
  
	itsCache_ = new MSCache(itsParent_);
	cacheUpdating = false;
	constructorSetup();
}

PlotMSPlot::~PlotMSPlot() {
	if (itsCache_)
		delete itsCache_;
	itsCache_ = NULL;
}

void PlotMSPlot::customizeAutoSymbol( const PlotSymbolPtr& baseSymbol, uInt dataSize ){
	if( baseSymbol->symbol() == PlotSymbol::AUTOSCALING) {
		if( dataSize > PIXEL_THRESHOLD ) {
			baseSymbol->setSymbol(PlotSymbol::PIXEL);
			baseSymbol->setSize(1,1);
		}
		else if( dataSize > MEDIUM_THRESHOLD ) {
			baseSymbol->setSymbol( PlotSymbol::CIRCLE);
			baseSymbol->setSize(2,2);
		}
		else if( dataSize > LARGE_THRESHOLD ) {
			baseSymbol->setSymbol( PlotSymbol::CIRCLE );
			baseSymbol->setSize(4,4);
		}
		else {
			baseSymbol->setSymbol( PlotSymbol::CIRCLE );
			baseSymbol->setSize(6,6);
		}
	}
}

void PlotMSPlot::customizeOverlaySymbol( const PlotSymbolPtr& baseSymbol, uInt dataSize ){
	if( dataSize > MEDIUM_THRESHOLD ) {
		baseSymbol->setSize(2,2);
	} else if( dataSize > LARGE_THRESHOLD ) {
		baseSymbol->setSize(3,3);
	} else if( dataSize > XLARGE_THRESHOLD ) {
		baseSymbol->setSize(4,4);
	} else {
		baseSymbol->setSize(6,6);
	}
}

// Public Methods //

void PlotMSPlot::resize(PlotMSPages &pages, uInt rows, uInt cols) {
	// Resize canvases and plots
	int plotCanvasRowCount = 1;
	int plotCanvasColCount = 1;

	if ( isIteration()  ){
		plotCanvasRowCount = rows;
		plotCanvasColCount = cols;
	}

	//Number of canvases is based on the grid.
	itsCanvases_.resize( plotCanvasRowCount );
	for( int r = 0; r < plotCanvasRowCount; ++r) {
		itsCanvases_[r].resize(plotCanvasColCount);
	}

	//Number of plots is based on how many overplots we are supporting (dataCount)
	//and on the iteration count over the data.
	Int plotRows = 1;
	Int plotCols = 1;
	getPlotSize( plotRows, plotCols );
	resizePlots( plotRows, plotCols );


	// Resize pages
	iterStep_ = plotCanvasRowCount * plotCanvasColCount;
	int ownedCanvasStart = this->getPageIterationCount( pages[0]);
	int pageSize = static_cast<int>(ceil( (plotCols*1.0f+ownedCanvasStart) / iterStep_));
	pages.resize( pageSize );
	for(size_t i = 0; i < pages.totalPages(); ++i) {
		pages[i].resize(rows, cols);
	}
}

String PlotMSPlot::name() const {
	const PMS_PP_MSData *data = itsParams_.typedGroup<PMS_PP_MSData>();
	const PMS_PP_Cache *cache = itsParams_.typedGroup<PMS_PP_Cache>();
	const PMS_PP_Display *display = itsParams_.typedGroup<PMS_PP_Display>();

	if(data == NULL || cache == NULL || display == NULL || !data->isSet())
		return "Over Plot";
	return display->titleFormat().getLabel(cache->xAxis(), cache->yAxis());
}

vector<MaskedScatterPlotPtr> PlotMSPlot::plots() const {
	if((itsPlots_.size() == 0) || (itsPlots_[0].size() == 0))
		return vector<MaskedScatterPlotPtr>();
	int index = 0;

	int dataCount = itsPlots_.size();
	int nIter = itsPlots_[0].size();

	int plotCount = nIter * dataCount;
	vector<MaskedScatterPlotPtr> v( plotCount );
	for(unsigned int i = 0; i < itsPlots_.size(); i++) {
		for(unsigned int j = 0; j < itsPlots_[i].size(); j++) {
			if(index >= plotCount)
				break;
			v[index] = itsPlots_[i][j];
			++index;
		}
	}
	return v;
}

vector<PlotCanvasPtr> PlotMSPlot::canvases() const {

	if(( itsCanvases_.size() == 0) || (itsCanvases_[0].size() == 0))
		return vector<PlotCanvasPtr>();
	uInt index = 0;
	uInt nIter = itsCache_->nIter(0);
	int canvasCount = std::min(nIter, uInt(itsCanvases_.size() * itsCanvases_[0].size()));
	vector<PlotCanvasPtr> v( canvasCount );
	for(uInt i = 0; i < itsCanvases_.size(); i++) {
		for(uInt j = 0; j < itsCanvases_[i].size(); j++) {
			if(index >= nIter) break;
			v[index] = itsCanvases_[i][j];
			++index;
		}
	}
	return v;
}

void PlotMSPlot::attachToCanvases() {
	Int nIter = itsCache_->nIter(0);
	if ( nIter <= 0 ){
		nIter = 1;
	}
	int canvasRows = itsCanvases_.size();
	for( int r = 0; r < canvasRows; ++r) {
		int canvasCols = itsCanvases_[r].size();
		for( int c = 0; c < canvasCols; ++c) {
			if(!itsCanvases_[r][c].null()) {
				if ( ! isIteration() ){
					//There is just one canvas for this plot,
					//but we may be adding several sets of data to it.
					int dataRowCount = itsPlots_.size();
					for ( int i = 0; i < dataRowCount; i++ ){
						int dataColCount = itsPlots_[i].size();
						for ( int j = 0; j < dataColCount; j++ )
							itsCanvases_[r][c]->plotItem( itsPlots_[i][j]);
					}
				} else {
					QList<PlotMSPlot*> canvasPlots = itsParent_->getPlotManager().getCanvasPlots(r,c);
					if ( canvasPlots.contains( this )){
						//For an iteration plot, there is one canvas per iteration.
						//In the case of overplotting with an iteration, we may be 
						//adding several sets of data to each canvas.
						PlotMSPage page = itsParent_->getPlotManager().itsPages_.currentPage();
						int iterationIndex = getIterationIndex (r, c, page );
						//int iterationIndex = r * canvasCols + c + iter;
						if ( iterationIndex < nIter ){
							int dataRowCount = itsPlots_.size();
							for ( int i = 0; i < dataRowCount; i++ ){
								if(!itsPlots_[i][iterationIndex].null())
									itsCanvases_[r][c]->plotItem(itsPlots_[i][iterationIndex]);
							}
						}
					}
				}
				((&*itsCanvases_[r][c]))->show();
				((&*itsCanvases_[r][c]))->setMinimumSize(5,5);
			}
		}
	}
}

void PlotMSPlot::detachFromCanvases() {
	for(uInt r = 0; r < itsCanvases_.size(); ++r) {
		for(uInt c = 0; c < itsCanvases_[r].size(); ++c) {
			if(!itsCanvases_[r][c].null()) {
				if(itsCanvases_[r][c]->numPlotItems() > 0) {
					int dataRowCount = itsPlots_.size();
					for ( int i = 0; i < dataRowCount; i++ ){
						int dataColCount = itsPlots_[i].size();
						for ( int j = 0; j < dataColCount; j++ ){
							itsCanvases_[r][c]->removePlotItem(itsPlots_[i][j]);
						}
					}
				}
				//This is necessary in scripting mode so that we don't see detached canvases.
				if (!itsParent_->guiShown() ||
						itsCanvases_[r][c]->numPlotItems() ==0 ){
					((&*itsCanvases_[r][c]))->hide();
				}
			}
		}
	}
}

void PlotMSPlot::dataMissing(){
	cacheUpdating = false;
	detachFromCanvases();
	initializePlot();
	releaseDrawing();
	itsCache_->clear();
}

vector<PMS::DataColumn> PlotMSPlot::getCachedData(){
	PMS_PP_Cache* cache = itsParams_.typedGroup<PMS_PP_Cache>();
	int xAxisCount = cache->numXAxes();
	int yAxisCount = cache->numYAxes();
	int count = xAxisCount + yAxisCount;
	vector<PMS::DataColumn> cdata( count );
	for( int i = 0; i < xAxisCount; ++i)
		cdata[i] = cache->xDataColumn(i);
	for( int i = xAxisCount; i < count; ++i)
		cdata[i] = cache->yDataColumn(i - xAxisCount);
	return cdata;
}

vector<PMS::Axis> PlotMSPlot::getCachedAxes() {
    PMS_PP_Cache* c = itsParams_.typedGroup<PMS_PP_Cache>();
    // get default axes if not given by user
    for(uInt i=0; i<c->numXAxes(); i++){
        if (c->xAxis(i) == PMS::NONE) 
            c->setXAxis(getDefaultXAxis(), i);
    }
    for(uInt i=0; i<c->numYAxes(); i++){
        if (c->yAxis(i) == PMS::NONE) {
            if (itsCache_->calType().startsWith("Xf")) {
                c->setYAxis(PMS::GPHASE, i);
            } else if (itsCache_->calType() == "GSPLINE") {
                PMS_PP_MSData* d = itsParams_.typedGroup<PMS_PP_MSData>();
                c->setYAxis(getGsplineAxis(d->filename()), i);
            } else {
                c->setYAxis(PMS::DEFAULT_YAXIS, i);
            }
        }
    }

    // add ATM/TSKY yaxis "under the hood" if valid xaxis
    if (c->showAtm() || c->showTsky()) {
        PMS::Axis xaxis = c->xAxis();
        bool validXAxis = (xaxis==PMS::CHANNEL || xaxis==PMS::FREQUENCY );
        if (!validXAxis) {
            c->setShowAtm(false);
            c->setShowTsky(false);
            itsParent_->showWarning("Overlays are valid only when xaxis is Channel or Frequency");
        } else {
            // add here for script client
            bool found(false);
            const vector<PMS::Axis> yAxes = c->yAxes();
            PMS::Axis atmAxis = (c->showAtm() ? PMS::ATM : PMS::TSKY);
            for (uInt i=0; i<yAxes.size(); ++i) {
                if (yAxes[i] == atmAxis) {
                    found=True;
                    break;
                }
            }
            if (!found) {
                // add ATM/TSKY to Cache axes
                int index = c->numXAxes();
                c->setAxes(xaxis, atmAxis, c->xDataColumn(0), 
                        PMS::DEFAULT_DATACOLUMN, index);
                // set Axes positions
                PMS_PP_Axes* a = itsParams_.typedGroup<PMS_PP_Axes>();
                a->resize(index+1, true);  // copy values
                a->setAxes(a->xAxis(index-1), Y_RIGHT, index);
                // keep same xaxis range
                a->setXRange(a->xRangeSet(index-1), a->xRange(index-1), index);
                // set Display symbol color
                PMS_PP_Display* disp = itsParams_.typedGroup<PMS_PP_Display>();
                PlotSymbolPtr atmSymbol = disp->unflaggedSymbol(index);
                atmSymbol->setSymbol("circle");
                atmSymbol->setSize(2,2);
                atmSymbol->setColor("#FF00FF");
                disp->setUnflaggedSymbol(atmSymbol, index);
                PlotSymbolPtr flaggedSymbol = disp->flaggedSymbol();
                disp->setFlaggedSymbol(flaggedSymbol, index);
            }
        }
    }
	vector<PMS::Axis> axes;
	for(uInt i=0; i<c->numXAxes(); i++)
		axes.push_back(c->xAxis(i));
	for(uInt i=0; i<c->numYAxes(); i++)
		axes.push_back(c->yAxis(i));
	return axes;
}

PMS::Axis PlotMSPlot::getDefaultXAxis() {
	PMS::Axis xaxis = PMS::TIME;
	if (itsCache_->cacheType() == PlotMSCacheBase::CAL) {
		String caltype = itsCache_->calType();
		if (caltype.contains("BPOLY"))
			xaxis = PMS::FREQUENCY;
		else if (caltype.contains("TSYS") || caltype[0]=='B' || caltype.contains("Mf") || caltype[0]=='X' )
			xaxis = PMS::CHANNEL;
		else if (caltype[0]=='D' || caltype[0]=='K')
			xaxis = PMS::ANTENNA1;
	}
	return xaxis;
}

PMS::Axis PlotMSPlot::getGsplineAxis(const String filename) {
	// When not set, set axis based on 'mode' column.
	// Modes are "AMP", "PHAS", and "A&P"
	PMS::Axis y(PMS::GAMP);
	Table tab(filename);
	TableColumn polymode(tab, "POLY_MODE");
	String mode = polymode.asString(0);
	if (mode == "PHAS") y = PMS::GPHASE;
	return y;
}

const PlotMSPlotParameters& PlotMSPlot::parameters() const{ return itsParams_;}

PlotMSPlotParameters& PlotMSPlot::parameters() { return itsParams_; }

vector<PlotCanvasPtr> PlotMSPlot::visibleCanvases() const {
	vector<PlotCanvasPtr> v;
	vector<PlotCanvasPtr> canv = canvases();
	for(unsigned int i = 0; i < canv.size(); i++) {
		bool visible = itsParent_->isVisible( canv[i] );
		if ( visible ) v.push_back( canv[i] );
	}
	return v;
}

Record PlotMSPlot::locateInfo(int plotIterIndex, const Vector<PlotRegion>& regions, bool showUnflagged,
		bool showFlagged, bool selectAll ) const {
	Record resultRecord = itsCache_->locateInfo(plotIterIndex, regions, showUnflagged, showFlagged, selectAll);
	return resultRecord;
}

PlotLogMessage* PlotMSPlot::locateRange( int canvasIndex, const Vector<PlotRegion> & regions, bool showUnflagged,
		bool showFlagged){
	int iterIndex = iter() + canvasIndex;
	PlotLogMessage* m = itsCache_->locateRange(iterIndex, regions, showUnflagged, showFlagged);
	return m;
}

PlotLogMessage* PlotMSPlot::flagRange( int canvasIndex, casa::PlotMSFlagging& flagging,
		const Vector<PlotRegion>& regions, bool showFlagged){
	int iterIndex = iter() + canvasIndex;
	PlotLogMessage* m = itsCache_->flagRange(iterIndex, flagging, regions, showFlagged);
	return m;
}

PlotMSRegions PlotMSPlot::selectedRegions() const {
	return selectedRegions(canvases()); }

PlotMSRegions PlotMSPlot::visibleSelectedRegions() const {
	return selectedRegions(visibleCanvases()); }

bool PlotMSPlot::initializePlot(PlotMSPages& pages) {
	bool hold = allDrawingHeld();
	if(!hold)
		holdDrawing();

	// Initialize plot objects and assign canvases.
	if(!assignCanvases(pages) || !initializePlot()) {
		if(!hold) releaseDrawing();
		return false;
	}

	// Set up page.
	pages.setupCurrentPage();

	// Attach plot objects to their assigned canvases.
	attachToCanvases();
    
	// Update objects with set parameters.
	parameters().releaseNotification();
	parametersHaveChanged(parameters(),
		PlotMSWatchedParameters::ALL_UPDATE_FLAGS());
    
	// Draw if necessary.
	if(!hold) releaseDrawing();
    
	return true;
}

bool PlotMSPlot::updateCache() {
	PMS_PP_MSData* data = itsParams_.typedGroup<PMS_PP_MSData>();
	PMS_PP_Cache* cache = itsParams_.typedGroup<PMS_PP_Cache>();
	PMS_PP_Iteration* iter = itsParams_.typedGroup<PMS_PP_Iteration>();
	PMS_PP_Display* disp = itsParams_.typedGroup<PMS_PP_Display>();
	if (data==NULL || cache==NULL || iter==NULL || disp==NULL ){
		return false;
	}

	// Don't load if data isn't set or there was an error during data opening.
	if (!data->isSet())
		return false;

	// Trap bad averaging/iteration combo
	if (data->averaging().baseline() && iter->iterationAxis()==PMS::ANTENNA) {
		logMessage( "Cannot iterate on Antenna if averaging over baseline, so turning off iteration.");
		iter->setIterationAxis(PMS::NONE);
	}

	// Notify the plots that the data will change
	updatePlots();

	// Set up cache loading parameters
	if (cache->numXAxes() != cache->numYAxes())
		return false;

	itsParent_->getLogger()->markMeasurement(PMS::LOG_ORIGIN,
			PMS::LOG_ORIGIN_LOAD_CACHE,
			PMS::LOG_EVENT_LOAD_CACHE);
	itsTCLParams_.endCacheLog = true;

	// Delete existing cache if it doesn't match
	String filename = data->filename();
	cacheUpdating = true;
	if (CacheFactory::needNewCache(itsCache_, filename)) {
		if(itsCache_) {
			clearPlotData(); //plot has ptr to indexer about to be deleted
			delete itsCache_;
			itsCache_ = NULL;
		}
		itsCache_ = CacheFactory::getCache(filename, itsParent_);
		if(itsCache_ == NULL) {
			throw AipsError("Failed to create a new Cache object!");
		} else {
			itsCache_->setFilename(filename);
			data->setType(itsCache_->cacheType());
		}
	}
	// Trap bad caltable/iteration and caltable/coloraxis combo
	String caltype = itsCache_->calType();
	checkColoraxis(caltype, disp);
	checkIteraxis(caltype, iter);

	bool result = true;
	try {
		result = itsParent_->updateCachePlot( this,
			PlotMSPlot::cacheLoaded, true );
			
	}
	catch( AipsError& error ){
		cacheUpdating = false;
		logMessage( error.getMesg().c_str() );
		result = false;
	}
	return result;
}

void PlotMSPlot::checkColoraxis(String caltype, PMS_PP_Display* display) {
	// check coloraxis for cal tables, xconnect for MS
	if (caltype.empty()) {
		if (display->xConnect() != "none")
			logMessage("WARNING: Connecting points is implemented for calibration tables only");
		return;
	}
	PMS::Axis coloraxis = display->colorizeAxis();
	if (coloraxis==PMS::INTENT) {
		logMessage("WARNING: Cannot colorize Intent for cal tables.\nTurning off colorize.");
		display->setColorize(false, PMS::NONE);
		return;
	}
	bool isBPoly(caltype=="BPOLY"), isGSpline(caltype=="GSPLINE");
	if ((isBPoly || isGSpline) && (coloraxis==PMS::BASELINE || coloraxis==PMS::ANTENNA2)) {
		logMessage("WARNING: Cannot colorize Baseline or Antenna2 for this cal table type.\nTurning off colorize.");
		display->setColorize(false, PMS::NONE);
	}
	if (isGSpline && (coloraxis==PMS::CHANNEL || coloraxis==PMS::FREQUENCY)) {
		logMessage("WARNING: Cannot colorize Channel or Frequency for GSPLINE cal table.\nTurning off colorize.");
		display->setColorize(false, PMS::NONE);
	}
}

void PlotMSPlot::checkIteraxis(String caltype, PMS_PP_Iteration* iter) {
	// check iteraxis for cal tables
	if (caltype.empty())
		return;
	bool isBPoly(caltype=="BPOLY"), isGSpline(caltype=="GSPLINE");
	PMS::Axis iteraxis = iter->iterationAxis();
	if ((isBPoly || isGSpline) && (iteraxis==PMS::BASELINE)) {
		logMessage("WARNING: Cannot iterate on Baseline for this cal table type.\nTurning off iteration.");
		iter->setIterationAxis(PMS::NONE);
		return;
	}
	if (isGSpline && (iteraxis==PMS::CHANNEL || iteraxis==PMS::FREQUENCY)) {
		logMessage("WARNING: Cannot iterate on Channel or Frequency for GSPLINE cal table.\nTurning off iteration.");
		iter->setIterationAxis(PMS::NONE);
	}
}

bool PlotMSPlot::updateCanvas() {

	PMS_PP_Axes* axes = itsParams_.typedGroup<PMS_PP_Axes>();
	PMS_PP_Cache* cache = itsParams_.typedGroup<PMS_PP_Cache>();
	PMS_PP_Canvas* canv = itsParams_.typedGroup<PMS_PP_Canvas>();
	PMS_PP_Iteration* iter = itsParams_.typedGroup<PMS_PP_Iteration>();
	PMS_PP_MSData* data = itsParams_.typedGroup<PMS_PP_MSData>();
	PMS_PP_Display* display = itsParams_.typedGroup<PMS_PP_Display>();
	if(axes==NULL || cache==NULL || canv==NULL || iter==NULL || data==NULL || display==NULL) {
		return false;
	}

	uInt nIter = itsCache_->nIter(0);
	uInt rows = itsCanvases_.size();
	for(uInt r = 0; r < rows; ++r) {
		uInt cols = itsCanvases_[r].size();
		for(uInt c = 0; c < cols; ++c) {
			PlotMSPage page = itsParent_->getPlotManager().itsPages_.currentPage();
			uInt iteration = getIterationIndex( r, c, page );
			if(iteration >= nIter) {
				clearCanvasProperties( r, c );
			} else {
				int numplots = rows*cols;
				setCanvasProperties(r, c, numplots, iteration, 
					axes, cache, canv, iter, data, display);
			}
		}
	}
	return true;
}

bool PlotMSPlot::updateDisplay() {
	try {
		PMS_PP_Cache *cache = itsParams_.typedGroup<PMS_PP_Cache>();
		PMS_PP_Axes *axes = itsParams_.typedGroup<PMS_PP_Axes>();
		PMS_PP_Display *display = itsParams_.typedGroup<PMS_PP_Display>();
		if(cache == NULL || axes == NULL || display == NULL)
			return false;

		MaskedScatterPlotPtr plot;
		int nIter = itsCache_->nIter(0);
		if ( nIter <= 0 )
			nIter = 1;
		uInt rows = itsPlots_.size();
		for(uInt row = 0; row < rows; ++row) {
			PMS::Axis x = cache->xAxis(row);
			PMS::Axis y = cache->yAxis(row);
			uInt cols = itsPlots_[row].size();
			for(uInt col = 0; col < cols; ++col) {
				// Set symbols.
				PlotSymbolPtr unflaggedSym = display->unflaggedSymbol(row);
				PlotSymbolPtr symbolUnmasked = itsParent_->createSymbol(unflaggedSym);
				uInt dataSize = itsCache_->indexer(row,col).sizeUnmasked();
				if (y==PMS::ATM || y==PMS::TSKY) 
					customizeOverlaySymbol( symbolUnmasked, dataSize );
				else 
					customizeAutoSymbol( symbolUnmasked, dataSize );

				PlotSymbolPtr flaggedSym = display->flaggedSymbol(row);
				PlotSymbolPtr symbolMasked = itsParent_->createSymbol(flaggedSym);
				dataSize = itsCache_->indexer(row,col).sizeMasked();
				if (y==PMS::ATM || y==PMS::TSKY)
					customizeOverlaySymbol( symbolMasked, dataSize );
				else
					customizeAutoSymbol( symbolMasked, dataSize );

				plot = itsPlots_[row][col];
				if (plot.null()) continue;

				plot->setSymbol(symbolUnmasked);
				plot->setMaskedSymbol(symbolMasked);
				// Colorize and set data changed, if redraw is needed
				String caltype = itsCache_->calType();
				checkColoraxis(caltype, display);
				bool colorizeChanged = itsCache_->indexer(row,col).colorize(display->colorizeFlag(), display->colorizeAxis());

				// Set xconnector in indexer and plot;
				// time connector changes indexer only
				String xconnector = display->xConnect();
				bool connectorChanged(false);
				if (itsCache_->cacheType()==PlotMSCacheBase::CAL) {
					bool timeconnector = display->timeConnect();
					connectorChanged = itsCache_->indexer(row,col).setConnect(xconnector, timeconnector);
					if (xconnector == "none") {
						plot->setLinesShown(false);
						plot->setMaskedLinesShown(false);
					} else {
						// only connect symbols being plotted
						if (symbolUnmasked->symbol() != PlotSymbol::NOSYMBOL) {
							plot->setLinesShown(true);
							plot->setLine(symbolUnmasked->getColor());
						}
						if (symbolMasked->symbol() != PlotSymbol::NOSYMBOL) {
							plot->setMaskedLinesShown(true);
							plot->setMaskedLine(symbolMasked->getColor());
						}
						bool step = (xconnector == "step");
						plot->setLinesStep(step);
						plot->setMaskedLinesStep(step);
					}
				}

				if ((nIter > 0) && (colorizeChanged || connectorChanged))
					plot->dataChanged();

				// Set item axes
				plot->setAxes(axes->xAxis(row), axes->yAxis(row));

				// Set plot title for legend; convert axes for cal table
				if (itsCache_->cacheType()==PlotMSCacheBase::CAL) {
					String caltype = itsCache_->calType();
					x = getCalAxis(caltype, x);
					y = getCalAxis(caltype, y);
				}
				vector<PMS::Axis> yAxes(1, y);
				vector<bool> yRefs(1, itsCache_->hasReferenceValue(y));
				vector<double> yRefValues(1, itsCache_->referenceValue(y));
				casacore::String title(display->titleFormat().getLabel(x, yAxes,
						itsCache_->hasReferenceValue(x),
						itsCache_->referenceValue(x),
						yRefs, yRefValues ));
				if (itsCache_->cacheType()==PlotMSCacheBase::CAL)
					title.gsub("Corr", "Poln");
				plot->setTitle(title);
			}
		}
	} catch(AipsError &err) {
		String errorMsg = "Could not update plot: " + err.getMesg();
		qDebug() << errorMsg.c_str();
		itsParent_->showError( errorMsg );
		return false;
	} catch(...) {
		String errorMsg = "Could not update plot, for unknown reasons!";
		qDebug() << errorMsg.c_str();
		itsParent_->showError( errorMsg );
		return false;
	}
	return true;
}

void PlotMSPlot::setColors() {
	uInt nIter = itsCache_->nIter(0);
	uInt rows = itsPlots_.size();
	itsColoredPlots_.resize(rows);
	for(uInt row = 0; row < rows; ++row) {
		uInt cols = itsPlots_[row].size();
		itsColoredPlots_[row].resize(cols);
		for(uInt col = 0; col < cols; ++col) {
			uInt iteration = row * cols + col;
			if(iteration >= nIter) break;
			itsColoredPlots_[row][col] = ColoredPlotPtr(dynamic_cast<ColoredPlot*>(&*itsPlots_[row][col]), false);
			if (!itsColoredPlots_[row][col].null()) {
				const vector<String> &colors = PMS::COLORS_LIST();
				for(uInt i = 0; i < colors.size(); ++i) 
					itsColoredPlots_[row][col]->setColorForBin(i, itsFactory_->color(colors[i]));
			} else {
				itsParent_->showError("Could not convert a plot in a ColoredPlot");
			}
		}
	}
}

bool PlotMSPlot::updateData() {
	itsCache_->clear();
	return true;
};

void PlotMSPlot::clearCanvases() {
	int rowCount = itsCanvases_.size();
	for ( int i = 0; i < rowCount; i++ ){
		int colCount = itsCanvases_[i].size();
		for ( int j = 0; j < colCount; j++ )
			clearCanvasProperties( i, j );
	}
}

bool PlotMSPlot::isCacheUpdating() const {
	return cacheUpdating;
}

void PlotMSPlot::setCacheUpdating( bool updating ){
	cacheUpdating = updating;
}

void PlotMSPlot::updatePlots() {
	for(uInt row = 0; row < itsPlots_.size(); ++row) {
		for(uInt col = 0; col < itsPlots_[row].size(); ++col) {
			bool plottable = itsParent_->getPlotManager().isPlottable(this);
			if(!itsPlots_[row][col].null() && plottable )
				itsPlots_[row][col]->dataChanged();
		}
	}
}

void PlotMSPlot::clearPlotData() {
	for(uInt row = 0; row < itsPlots_.size(); ++row) {
		for(uInt col = 0; col < itsPlots_[row].size(); ++col) {
			bool plottable = itsParent_->getPlotManager().isPlottable(this);
			if(!itsPlots_[row][col].null() && plottable )
				itsPlots_[row][col]->clearData();
		}
	}
}
bool PlotMSPlot::updateIndexing() {
	PMS_PP_Iteration *iter = itsParams_.typedGroup<PMS_PP_Iteration>();
	PMS_PP_Axes* axes = itsParams_.typedGroup<PMS_PP_Axes>();
	PMS_PP_Display* disp = itsParams_.typedGroup<PMS_PP_Display>();
	bool globalX = iter->isGlobalScaleX();
	bool globalY = iter->isGlobalScaleY();
	PMS::Axis iterAxis = iter->iterationAxis();

	String caltype = itsCache_->calType();
	checkIteraxis(caltype, iter);

	int dataCount = axes->numYAxes();
	//Only update if we need to.
	bool requiredUpdate = false;

	for ( int i = 0; i < dataCount; i++ ){
		bool iterationInitialized = itsCache_->isIndexerInitialized(iterAxis, globalX, globalY, i);
		if ( !iterationInitialized ){
			requiredUpdate = true;
			break;
		}
	}

	if ( requiredUpdate ){
		itsCache_->clearRanges();
		//Set up the indexer.
		for ( int i = 0; i < dataCount; i++ ) {
	        String xconnect(disp->xConnect());
	        bool timeconnect(disp->timeConnect());
			itsCache_->setUpIndexer(iterAxis, globalX, globalY, xconnect, timeconnect, i);
        }
	}
	return true;
}

void PlotMSPlot::logPoints() {
	PMS_PP_Display *display = itsParams_.typedGroup<PMS_PP_Display>();
	bool showUnflagged = display->unflaggedSymbol()->symbol() != PlotSymbol::NOSYMBOL;
	bool showFlagged = display->flaggedSymbol()->symbol() != PlotSymbol::NOSYMBOL;
	bool allFlagged = false;

	stringstream ss;
	ss << "Plotting ";
	if(showUnflagged) {
		if ( itsCache_->nIter(0) > iter_ ){
			uInt nUnflaggedPoints = itsCache_->indexer(0,iter_).sizeUnmasked();
			ss << nUnflaggedPoints << " unflagged" << (showFlagged ? ", " : "");
			if (nUnflaggedPoints==0) allFlagged = true;
		} else {
			ss << "0 unflagged" <<(showFlagged ? ", " : "");
		}
	}
	if(showFlagged) {
		if ( itsCache_->nIter(0) > iter_ ){
			ss << itsCache_->indexer(0,iter_).sizeMasked() << " flagged";
		} else {
			ss << "0 flagged";
		}
	}
	ss << " points.";

	itsParent_->getLogger()->postMessage(PMS::LOG_ORIGIN,
			PMS::LOG_ORIGIN_PLOT,
			ss.str(),
			PMS::LOG_EVENT_PLOT);
	if (allFlagged) {
		itsParent_->showWarning("All selected data are flagged.");
	} else { //clear warning
		itsParent_->clearMessage();
	}
}

void PlotMSPlot::logIter(Int iter, Int nIter) {
	if(nIter > 1) {
		stringstream ss;
		ss << "Stepping to iteration = " << iter+1
				<< " (of " << nIter << "): "
				<< itsCache_->indexer(0,iter).iterLabel();
		itsParent_->getLogger()->postMessage(PMS::LOG_ORIGIN,
				PMS::LOG_ORIGIN_PLOT,
				ss.str(),
				PMS::LOG_EVENT_PLOT);
	}
}

void PlotMSPlot::parametersHaveChanged(const PlotMSWatchedParameters& p,
        int updateFlag ) {
    if ( isCacheUpdating() )
		return;
    // Make sure it's this plot's parameters.
    if( &p != &parameters() )
		return;

    //A plot not to be shown.
    bool plottable = itsParent_->getPlotManager().isPlottable( this );
    if ( ! plottable ){
        //Clear the plot
        detachFromCanvases();
        return;
	}

    vector<String> updates = PlotMSWatchedParameters::UPDATE_FLAG_NAMES(updateFlag);
    if(updates.size() == 0)
		return;
    
    // Log what we're going to be updating.
    stringstream ss;
    ss << "Updating: ";
    for(unsigned int i = 0; i < updates.size(); i++) {
        if(i > 0) ss << ", ";
        ss << updates[i];
    }
    ss << ".";
    itsParent_->getLogger()->postMessage(PMS::LOG_ORIGIN,
            PMS::LOG_ORIGIN_PARAMS_CHANGED, ss.str(),
            PMS::LOG_EVENT_PARAMS_CHANGED);
    int updateRedraw= updateFlag & PMS_PP::UPDATE_REDRAW;
    bool releaseWhenDone = !allDrawingHeld() && updateRedraw;
    if(releaseWhenDone)
        holdDrawing();
    
    // Update MS as needed.
    const PMS_PP_MSData* d = parameters().typedGroup<PMS_PP_MSData>();
    bool dataSuccess = d->isSet();
    if(dataSuccess && (updateFlag & PMS_PP::UPDATE_MSDATA)){
        bool fileExists = true;
        String fileName = d->filename();
        ifstream ifile( fileName.c_str() );
        if ( !ifile ) {
            fileExists = false;
        } else {
            ifile.close();
        }

        if ( fileExists ){
            dataSuccess = updateData();
        } else {
            String errorMessage( "Please check that the file ");
            errorMessage.append( fileName );
            errorMessage.append( " is valid.");
            itsParent_->getLogger()->postMessage(PMS::LOG_ORIGIN,
                PMS::LOG_ORIGIN_PLOT, errorMessage, PlotLogger::MSG_WARN);
            dataSuccess = false;
        }
    }

    // If something went wrong, clear the cache and plots.
    if (!dataSuccess) {
        itsCache_->clear();
        plotDataChanged();
    }
    
    // Let the child handle the rest of the parameter changes, and release
    // drawing if needed.
    bool result = parametersHaveChanged_(p,updateFlag,releaseWhenDone);
    if( result && releaseWhenDone){
		//Note::this was put in because when reload was checked from the gui
		//We were getting a segfault because the plot was redrawing before the
		//cache was loaded from a thread.  There seems to be a mechanism in
		//place to release the drawing later after the cache is loaded.
		if ( ! itsParent_->guiShown() )
			releaseDrawing();
    } 
}

void PlotMSPlot::plotDataChanged() {
    bool hold = allDrawingHeld();
    if(!hold) holdDrawing();
    
    vector<MaskedScatterPlotPtr> p = plots();
    for(unsigned int i = 0; i < p.size(); i++){
        if(!p[i].null())
            p[i]->dataChanged();
    }
    if(!hold)
        releaseDrawing();
}

bool PlotMSPlot::isIteration() const {
	const PMS_PP_Iteration *iter = itsParams_.typedGroup<PMS_PP_Iteration>();
	bool iterationPlot = false;
	if ( iter != NULL )
		iterationPlot = iter->isIteration();
	return iterationPlot;
}



bool PlotMSPlot::exportToFormat(const PlotExportFormat& format) {
    vector<PlotCanvasPtr> canv = canvases();
    bool exportSuccess = true;

    //Determine how many pages we need to print.
    int pageCount = 1;
    //Store the current page.
    Int currentIter = iter();

    PlotMSExportParam& exportParams = itsParent_->getExportParameters();
    PMS::ExportRange range = exportParams.getExportRange();
    if ( range == PMS::PAGE_ALL ){
    	int iterationCount = itsCache_->nIter( 0 );
    	float divResult = (iterationCount * 1.0f) / canv.size();
    	pageCount = static_cast<int>(ceil( divResult ));
    	//If we are an iteration plot and we don't own the first few plots on the
    	//page we may need to bump the page count up by one.
    	if ( isIteration() ){
    		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
    		PlotMSPage firstPage = pages.getFirstPage();
    		int firstPagePlotCount = getPageIterationCount( firstPage );
    		if ( firstPagePlotCount < static_cast<int>(canv.size()) ){
    			int notOwnedCount = canv.size() - firstPagePlotCount;
    			int excessSpace = (pageCount * canv.size()) - (notOwnedCount + iterationCount );
    			if ( excessSpace < 0 )
    				pageCount = pageCount + 1;
    		}
    	}
    	firstIter();
    }

    PlotExportFormat exportFormat( format );
    String baseFileName = format.location;
    String suffix = "";
    int periodIndex = baseFileName.find_last_of( ".");
    // Remove the last '.' from the storage location.
    if ( periodIndex != static_cast<int>(String::npos) ){
        suffix = baseFileName.substr( periodIndex, baseFileName.size() - periodIndex);
        baseFileName = baseFileName.substr(0, periodIndex );
    }

    // Loop over all the iterations, exporting them
    waitOnCanvases();
    PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
    const String sep( "_");
    bool shortenName = false;
    for ( int i = 0; i < pageCount; i++ ){
        String pageStr, itersInclude;
        if ( i > 0 )
            pageStr = String::toString( i+1 );
        size_t amp;
        if (isIteration()){
            int iterStart = this->iter_;
            int iterEnd = getPageIterationCount( pages[i]);
            int lastIndex = std::min(this->nIter()-iterStart, iterEnd) + iterStart;
            int index = iterStart;

            while ( index < lastIndex ){
                String iterId;
                if ( index == iterStart ){
                    iterId = itsCache_->indexer(0,index).fileLabel();
                    amp = iterId.find(" & ");
                    if (amp != string::npos)
                        iterId.replace(amp, 3, "_with_");
                } else {
                    iterId = itsCache_->indexer(0,index).iterValue();
                    amp = iterId.find(" & ");
                    if (amp != string::npos)
                        iterId.replace(amp, 3, "_with_");
                }
                if ( index < lastIndex - 1 )
                    iterId = iterId + ",";
                itersInclude = itersInclude + iterId;
                index++;
            }
        }

        String fileId;
        if ( itersInclude.size() > 0 )
            fileId = sep + itersInclude;
        if ( pageStr.size() > 0 )
            fileId = fileId + sep + pageStr;

        std::string::size_type filepos = baseFileName.rfind('/');
        // returned by os.getcwd() in task_plotms.py
        String path = baseFileName.substr(0, filepos+1);
        // user's plotfile, minus suffix
        String filename = baseFileName.substr(filepos+1, baseFileName.size());
        // This could add iteration names, plus suffix
        String exportFileName = filename + fileId + suffix;
        // CAS-7777 - file(s) will not export if filename > 255
        // Need cushion for sep, pageStr, additional digits, etc.
        // Hopefully first filename is representative of the rest!
        if ((exportFileName.length() > 256) || shortenName) {
            if (pageStr.size() > 0) {
                // shorten to 'basename_#.ext' e.g. 'test_2.jpg'
                exportFileName = filename + sep + pageStr + suffix;
            } else {
                // no iteration, just use 'basename.ext'
                exportFileName = filename + suffix;
            }
            // if shorten one name, shorten them all
            shortenName = true;
        }

        // check if shortening filename didn't work
        if (exportFileName.length() > 255) {
            logMessage("ERROR: Export filename exceeds length limit (256).  Export failed.");
            exportSuccess = false;
        } else {
            exportFormat.location = path + exportFileName;
            exportSuccess = itsParent_->exportToFormat( exportFormat );
            if (exportSuccess) {
                // let user know exported filename (if added iteration label)
                String msg = "Exported " + exportFormat.location;
                logMessage(msg.c_str());
            }
            waitOnCanvases();
            if ( i < pageCount - 1 )
                nextIter();
            waitOnCanvases();
        }
    }

    // Warn user if shortened plotfile name
    if (exportSuccess && shortenName)
        logMessage("Export filenames do not include iteration labels so that plotfile names do not exceed length limit (256).");

    //Restore the current page
    setIter( currentIter );
    return exportSuccess;
}

void PlotMSPlot::exportToFormatCancel(){
	vector<PlotCanvasPtr> canv = canvases();
	PlotOperationPtr op;
	for(unsigned int i = 0; i < canv.size(); i++) {
		if(canv[i].null()) continue;
		op = canv[i]->operationExport();
		if(op.null()) continue;
		op->setCancelRequested(true);
	}
}

void PlotMSPlot::cacheLoaded_(bool wasCanceled) {
    // Ensure we fail gracefully if cache loading yielded nothing or was cancelled

    if ( itsCache_ == NULL ){
        return;
    }
    if (!itsCache_->cacheReady() || wasCanceled) {
        dataMissing();
        return;
    }
    // Report we are done
    if(itsTCLParams_.endCacheLog)
        itsParent_->getLogger()->releaseMeasurement();

    // Make this more specific than canvas-triggered
    if (itsTCLParams_.updateCanvas || itsTCLParams_.updateIteration ){
        updateIndexing();
    }

    // Reset the iterator (if data are new)
    bool iterRecalculated = resetIter();

    // These are called in recalculateIteration, so don't call again unless necessary
    if (!iterRecalculated) {
        // Let the plot know that the data has been changed as needed,
        // unless the thread was canceled.
        updatePlots();

        // Update display as needed.  Put this before update canvas so
        // that the legend item keys will have the correct color.
        if(itsTCLParams_.updateDisplay)
            updateDisplay();

        // Update canvas as needed.
        if(itsTCLParams_.updateCanvas)
            updateCanvas();
    }

    // Release drawing if needed.
    if(itsTCLParams_.releaseWhenDone && !isCacheUpdating() )
        releaseDrawing();
}

void PlotMSPlot::setRelease( bool b ){
    itsTCLParams_.releaseWhenDone = b;
}

void PlotMSPlot::canvasWasDisowned(PlotCanvasPtr canvas) {
    if(canvas.null()) return;

    vector<MaskedScatterPlotPtr> p = plots();
    for(unsigned int i = 0; i < p.size(); i++)
        if(!p[i].null()) canvas->removePlotItem(p[i]);
}


// Protected Methods //

bool PlotMSPlot::initializePlot() {
	Int rows = 1;
	Int cols = 1;
	getPlotSize( rows, cols );
	resizePlots( rows, cols );
	setColors();
	return true;
}

bool PlotMSPlot::parametersHaveChanged_(const PlotMSWatchedParameters &p,
		int updateFlag, bool releaseWhenDone) {

	if(&p != &itsParams_)
		return false;

	const PMS_PP_MSData *data = itsParams_.typedGroup<PMS_PP_MSData>();
	const PMS_PP_Iteration *iter = itsParams_.typedGroup<PMS_PP_Iteration>();
	const PMS_PP_Axes *axes = itsParams_.typedGroup<PMS_PP_Axes>();
	const PMS_PP_Display *display = itsParams_.typedGroup<PMS_PP_Display>();
	if(data == NULL || iter == NULL || axes == NULL || display == NULL)
		return true;

	bool isConnected(false);
	std::vector<String> connects = display->xConnects();
	for (size_t i=0; i<connects.size(); ++i) {
		if (connects[i]=="line" || connects[i]=="step")
			isConnected = true;
    }

	itsTCLParams_.releaseWhenDone = releaseWhenDone;
	itsTCLParams_.updateCanvas = (updateFlag & PMS_PP::UPDATE_AXES) ||
			(updateFlag & PMS_PP::UPDATE_CACHE) ||
			(updateFlag & PMS_PP::UPDATE_CANVAS) ||
			(updateFlag & PMS_PP::UPDATE_ITERATION) ||
			(updateFlag & PMS_PP::UPDATE_MSDATA) ||
			isConnected || !data->isSet();

	itsTCLParams_.updateDisplay = updateFlag & PMS_PP::UPDATE_DISPLAY;
	itsTCLParams_.endCacheLog = false;

	// Clear selection if axes change
	// UPDATE_CACHE should be close enough for now (I hope)
	int updateCacheFlag = updateFlag & PMS_PP::UPDATE_CACHE;
	if( updateCacheFlag ) {
		for(size_t r = 0; r < itsCanvases_.size(); ++r) {
			for(size_t c = 0; c < itsCanvases_[r].size(); ++c) {
				PlotCanvasPtr plotCanvas = itsCanvases_[r][c];
				if ( ! plotCanvas.null() ){
					plotCanvas->standardMouseTools()->selectTool()->clearSelectedRects();
					plotCanvas->clearAnnotations();
					plotCanvas->clearShapes();
				}
			}
		}
	}

	//See if the iteration parameters have changed.
	bool commonAxisX = iter->isCommonAxisX();
	bool commonAxisY = iter->isCommonAxisY();
	Int rows(0), cols(0);
	getPlotSize( rows, cols );
	PlotAxis locationAxisX = axes->xAxis();
	PlotAxis locationAxisY = axes->yAxis();
	int displayRow = iter->getGridRow();
	int displayCol = iter->getGridCol();
	int plotRows = itsPlots_.size();
	int plotCols(0);
	if ( plotRows > 0 )
		plotCols = itsPlots_[0].size();
	bool locationChange(false);

	if ( (gridRow != displayRow || gridCol != displayCol) && gridRow != -1 ){
		locationChange = true;
		//This removes the title and axes from previous plot location.
		QList<PlotMSPlot*> canvasPlots = itsParent_->getPlotManager().getCanvasPlots( gridRow, gridCol);
		if ( canvasPlots.size() == 1 ){
			//We are the sole occupant of the old spot (no overplotting)
			//so we erase all evidence of there being a plot
			itsParent_->getPlotManager().clearCanvas(gridRow, gridCol);
		} else if ( canvasPlots.size() > 1 ){
			//Just erase ourselves from the canvas.
			itsParent_->getPlotManager().itsPages_.disown( gridRow, gridCol, this );
			detachFromCanvases();
			//Tell the other plots to redraw
			for ( int i = 0; i < canvasPlots.size(); i++ ){
				if ( canvasPlots[i] != this )
					canvasPlots[i]->parametersHaveChanged( canvasPlots[i]->parameters(),PMS_PP::UPDATE_REDRAW );
			}
		}
	}

	bool updateIter = updateFlag & PMS_PP::UPDATE_ITERATION;
	itsTCLParams_.updateIteration = ( updateIter ||
		((plotRows != rows) || (plotCols != cols)) ||
		(itsParent_->isCommonAxisX() != commonAxisX) ||
		(itsParent_->isCommonAxisY() != commonAxisY) ||
		(itsParent_->getAxisLocationX() != locationAxisX) ||
		(itsParent_->getAxisLocationY() != locationAxisY) ||
		locationChange );
	itsParent_->setCommonAxes( commonAxisX, commonAxisY);
	itsParent_->setAxisLocation( locationAxisX, locationAxisY);
	gridRow = displayRow;
	gridCol = displayCol;

	//We are not plotting this particular plot so just clear it and return.
	if ( displayRow == -1 || displayCol == -1 ){
		clearCanvases();
		return true;
	}

	bool dataSet = data->isSet();
	bool updateData = (updateFlag & PMS_PP::UPDATE_MSDATA) || (updateFlag & PMS_PP::UPDATE_CACHE);

	//If the iteration count has changed, ie from an iteration to a
	//non-iteration or just a change in the iteration axis, we may need
	//to clear the cache and update it.
	if ( updateIter ){
		PMS::Axis newAxis = iter->iterationAxis();
		PMS::Axis cacheIterationAxis = itsCache_->getIterAxis();
		if ( newAxis != cacheIterationAxis ){
			updateData = true;
		}
	}

	// Update cache if needed
	bool handled = true;
	if( dataSet && updateData ) {
		try {
			handled = updateCache();
		}
		catch( AipsError& error ){
			cerr << "Could not update cache: "<<error.getMesg().c_str()<<endl;
			cacheLoaded_(false);
			handled = false;
		}
	}
	return handled;
}

void PlotMSPlot::constructorSetup() {
	itsCache_->setPlot(this);
	PlotMSPlotParameters& params = parameters();
	params.addWatcher(this);
	// hold notification until initializePlot is called
	params.holdNotification(this);
	gridRow = -1;
	gridCol = -1;
	makeParameters(params, itsParent_);
}

bool PlotMSPlot::allDrawingHeld() {
	vector<PlotCanvasPtr> canv = canvases();
	bool allDrawingHeld = true;
	int canvasCount = canv.size();
	for(int i = 0; i < canvasCount; i++){
		if(!canv[i].null()){
			bool canvasDrawingHeld = canv[i]->drawingIsHeld();
			if ( !canvasDrawingHeld ){
				allDrawingHeld = false;
				break;
			}
		}
	}
	return allDrawingHeld;
}

void PlotMSPlot::holdDrawing() {
	vector<PlotCanvasPtr> canv = canvases();
	for(unsigned int i = 0; i < canv.size(); i++){
		if ( !canv[i].null() ){
			bool canvasDrawing = canv[i]->isDrawing();
			if ( canvasDrawing ){
				waitOnCanvas( canv[i]);
			}
			canv[i]->holdDrawing();
		}
	}
}

void PlotMSPlot::releaseDrawing() {
	vector<PlotCanvasPtr> canv = canvases();
	for(unsigned int i = 0; i < canv.size(); i++){
		if(!canv[i].null()){
			if ( canv[i]->drawingIsHeld()){
				canv[i]->releaseDrawing();
			}
		}
	}
}

void PlotMSPlot::waitOnCanvas( const PlotCanvasPtr& canvas ){
	if ( !canvas.null()){
		int callIndex = 0;
		int maxCalls =  60;

		bool scriptClient = !itsParent_->guiShown();
		if ( scriptClient )
			return;

		bool canvasDrawing = canvas->isDrawing( );
		while(  canvasDrawing && callIndex < maxCalls ){
			usleep(1000000);
			callIndex++;
			canvasDrawing = canvas->isDrawing();
		}
	}
}

void PlotMSPlot::waitOnCanvases(){
	vector<PlotCanvasPtr> canv = canvases();
		for (unsigned int i = 0; i < canv.size(); i++ ){
			if ( !canv[i].null()){
				waitOnCanvas( canv[i]);
			}
		}
}

void PlotMSPlot::waitForDrawing( bool holdDrawing ){
	vector<PlotCanvasPtr> canv = canvases();
	for (unsigned int i = 0; i < canv.size(); i++ ){
		if ( !canv[i].null()){
			waitOnCanvas( canv[i]);
			if ( holdDrawing ){
				canv[i]->holdDrawing();
			}
		}
	}
	detachFromCanvases();
}


bool PlotMSPlot::firstIter() {
	Int nIter = itsCache_->nIter(0);
	if( (nIter > 1) && (iter_ != 0) ) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		pages.firstPage();
		iter_ = 0;
		recalculateIteration();
		return true;
	}
	return false;
}

bool PlotMSPlot::prevIter() {
	Int nIter = itsCache_->nIter(0);
	if( nIter > 1 && iter_ > 0 ) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		iter_ -= getPageIterationCount( pages.currentPage() );
		if (iter_ < 0) iter_=0;  // just in case
		pages.previousPage();
		recalculateIteration();
		return true;
	}
	return false;
}

bool PlotMSPlot::nextIter() {
	Int nIter = itsCache_->nIter(0);
	if( nIter > 1) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		int pageIterCount = getPageIterationCount(pages.currentPage());
		if((iter_+pageIterCount) < nIter ) {
			iter_ += pageIterCount;
			pages.nextPage();
			recalculateIteration();
			return true;
		}
	}
	return false;
}

bool PlotMSPlot::lastIter() {
	Int nIter = itsCache_->nIter(0);
	if((nIter > 0) && (iter_ < (nIter - iterStep_))) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		int firstPageIterCount = getPageIterationCount( pages.getFirstPage() );
		iter_ = int(double(nIter-1) / iterStep_) * iterStep_;
		if(iterStep_ == 1){
			iter_ = nIter - 1;
		} else {
			if ( firstPageIterCount < iterStep_ )
				iter_ = iter_ - (iterStep_ - firstPageIterCount );
		}
		pages.lastPage();
		recalculateIteration();
		return true;
	}
	return false;
}

bool PlotMSPlot::setIter( int index ){
	Int nIter = itsCache_->nIter(0);
	bool successful = false;
	if( nIter > 1 && index < nIter && index >= 0) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		pages.setCurrentPageNum( index );
		iter_ = index;
		recalculateIteration();
		successful = true;
	}
	return successful;
}

bool PlotMSPlot::resetIter() {
	Int nIter = itsCache_->nIter(0);
	if(nIter > 0 ) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		pages.firstPage();
		iter_ = 0;
		recalculateIteration();
		return true;
	}
	return false;
}

void PlotMSPlot::recalculateIteration( ) {
	bool drawingHeld = allDrawingHeld();
	if ( !drawingHeld )
		this->holdDrawing();

	int nIter = itsCache_->nIter(0);
	if ( nIter <= 0 )
		nIter = 1;

	detachFromCanvases();
	if(itsTCLParams_.updateIteration  || isIteration()) {
		PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
		assignCanvases(pages);
	}

	//Put the data into the plot
	uInt rows = itsPlots_.size();
	for(uInt r = 0; r < rows; ++r) {
		uInt cols = itsPlots_[r].size();
		for(uInt c = 0; c < cols; ++c) {
			int iterationIndex = c;
			if(iterationIndex >= nIter ) break;
			logIter(iterationIndex, nIter);
			PlotMaskedPointDataPtr data(&(itsCache_->indexer(r,c)), false);
			itsPlots_[r][c] = itsFactory_->maskedPlot(data);
		}
	}

	setColors();
	itsTCLParams_.updateDisplay = true;

	//Update display should come before update canvas so that the
	//legend items get the correct color.
	updateDisplay();
	updateCanvas();
	attachToCanvases();
	updatePlots();
	if ( !isCacheUpdating() && !drawingHeld ){
		releaseDrawing();
	}

	//Add for CAS-6928/CAS-7014.  Tools were not being reset when
	//the iteration plot page changed.
	// CAS-3125 Need to do this after plots are set for axes stack
	itsParent_->resetTools();
	logPoints();
}

Int PlotMSPlot::nIter() {
	Int iterationCount = 0;
	if ( itsCache_ != NULL )
		iterationCount = itsCache_->nIter(0);
	return iterationCount;
}

int PlotMSPlot::getPageIterationCount( const PlotMSPage& page ) {
	int rows = itsCanvases_.size();
	int cols = 0;
	if ( rows > 0 )
		cols = itsCanvases_[0].size();
	int iterationCanvasCount = getIterationIndex(rows,cols,page);
	iterationCanvasCount = iterationCanvasCount - iter_;
	return iterationCanvasCount;
}

void PlotMSPlot::updateLocation(){

	PlotMSPages &pages = itsParent_->getPlotManager().itsPages_;
	//Initializes the canvases for this plot
	assignCanvases(pages);
	//Put the plot data on the canvas.
	attachToCanvases();
	//For scripting mode, we get plots without axes if the call is not preset.
	if ( !itsParent_->guiShown()  ){
		//Put the plot axis on the canvas.
		updateCanvas();
	}
}

PlotMSRegions PlotMSPlot::selectedRegions(const vector<PlotCanvasPtr>& canvases) const {
	PlotMSRegions r;
	PMS::Axis x = (PMS::Axis)PMS_PP_RETCALL(itsParams_, PMS_PP_Cache, xAxis, 0);
	PMS::Axis y = (PMS::Axis)PMS_PP_RETCALL(itsParams_, PMS_PP_Cache, yAxis, 0);

	for(uInt i = 0; i < canvases.size(); ++i) 
		r.addRegions(x, y, canvases[i]);
	return r;
}

bool PlotMSPlot::assignCanvases(PlotMSPages &pages) {
	if(pages.totalPages() == 0) {
		pages.insertPage();
		pages.firstPage();
	}

	//Resize based on the row and column count
	PlotMSParameters params = itsParent_->getParameters();
	int rows = params.getRowCount();
	int cols = params.getColCount();
	resize( pages, rows, cols );
	int currentPage = pages.currentPageNumber();
	PlotMSPage& page = pages[currentPage];

	const PMS_PP_Iteration* iterParams = itsParams_.typedGroup<PMS_PP_Iteration>();
	int rowIndex = 0;
	int colIndex = 0;
	if ( iterParams != NULL ){
		rowIndex = iterParams->getGridRow();
		colIndex = iterParams->getGridCol();
	}

	page.disown( this );
	if ( rowIndex >= 0 && colIndex >= 0 ){
		//Find a canvas for this plot.
		for(int r = 0; r < rows; ++r) {
			bool assigned = false;
			for(int c = 0; c < cols; ++c) {
				if ( isIteration() ){
					if ( (r > rowIndex) || (r == rowIndex && c>=colIndex) || currentPage > 0 ){
						if( !page.isOwned(r, c)) {
							page.setOwner(r, c, this);
							itsCanvases_[r][c] = page.canvas(r, c);
						}
					}
				} else {
					//If it is not an iteration plot, there is just one canvas for this plot.
					if ( rowIndex == r && colIndex == c){
						page.setOwner(r, c, this);
						itsCanvases_[0][0] = page.canvas(r,c);
						assigned = true;
						break;
					}
				}
			}
			if ( assigned )
				break;
		}
	}

	page.setupPage();
	return true;
}

void PlotMSPlot::resizePlots( int rows, int cols ){
	itsPlots_.resize( rows );
	for ( int r = 0; r < rows; ++r) {
		itsPlots_[r].resize( cols );
		for ( int c = 0; c < cols; ++c) {
			//Put empty data into the plot.
			PlotMaskedPointDataPtr data(&(itsCache_->indexer0()), false);
			itsPlots_[r][c] = itsFactory_->maskedPlot(data);
		}
	}
}

void PlotMSPlot::getPlotSize( Int& rows, Int& cols ){
	rows = 1;
	cols = 1;
	//Number of plots is based on how many overplots we are supporting (dataCount)
	//and on the iteration count over the data.
	const PMS_PP_Axes *axes = itsParams_.typedGroup<PMS_PP_Axes>();
	if ( axes != NULL )
		rows = axes->numYAxes();

	int iterationCount = itsCache_->nIter(0);
	if ( iterationCount > 0 )
		cols = iterationCount;
}

int PlotMSPlot::getIterationIndex( int r, int c, const PlotMSPage& page ){
	int iterationIndex = iter_;
	bool found = false;
	int rows = page.canvasRows();
	int cols = page.canvasCols();
	for ( int i = 0; i < rows; i++ ){
		for ( int j = 0; j <= cols; j++ ){
			if ( i == r && j == c ){
				found =true;
				break;
			} else {
				bool ownsCanvas = page.isOwner(i,j, this);
				if ( ownsCanvas )
					iterationIndex++;
			}
		}
		if ( found )
			break;
	}
	return iterationIndex;
}

void PlotMSPlot::logMessage( const QString& msg ) const {
	if ( itsParent_ != NULL ){
		stringstream ss;
		ss << msg.toStdString().c_str();
		itsParent_->getLogger()->postMessage(PMS::LOG_ORIGIN,
					PMS::LOG_ORIGIN_PLOT,
					ss.str(),
					PMS::LOG_EVENT_PLOT);
	}
}

void PlotMSPlot::clearCanvasProperties( int row, int col){
	PlotCanvasPtr canvas = itsCanvases_[row][col];
	if(canvas.null())
		return;
	canvas->showAllAxes( false );
	canvas->setTitle( "" );
	canvas->setCommonAxes( false, false );
}

void PlotMSPlot::setCanvasProperties (int row, int col, int numplots, uInt iteration,
		PMS_PP_Axes* axesParams, PMS_PP_Cache* cacheParams, PMS_PP_Canvas *canvParams,
		PMS_PP_Iteration *iterParams, PMS_PP_MSData* dataParams, PMS_PP_Display* displayParams) {

	PlotCanvasPtr canvas = itsCanvases_[row][col];
	if(canvas.null())
		return;
	canvas->showAllAxes(false);
	canvas->clearAxesLabels();

	// used throughout
	int yAxisCount = axesParams->numYAxes();
	bool set = dataParams->isSet();

	// Whether to share common axes for iterated plots on grid
	bool commonX = iterParams->isCommonAxisX();
	bool commonY = iterParams->isCommonAxisY();
	canvas->setCommonAxes( commonX, commonY );
	// showX and showY determine whether axes are visible at all.
	bool showX = set && canvParams->xAxisShown();
	bool showY = set && canvParams->yAxisShown();
	PlotAxis cx = axesParams->xAxis();
	canvas->showAxis(cx, showX);
	for ( int i = 0; i < yAxisCount; i++ ) {
		PlotAxis cy = axesParams->yAxis( i );
		canvas->showAxis(cy, showY);
	}

	// xaxis, scale (TIME/NORMAL), ref value
	PMS::Axis x = cacheParams->xAxis();
	if (x==PMS::NONE) {
		x = getDefaultXAxis();
		cacheParams->setXAxis(x);
	}
	canvas->setAxisScale(cx, PMS::axisScale(x));
	bool xref = itsCache_->hasReferenceValue(x);
	double xrefval = itsCache_->referenceValue(x);
	canvas->setAxisReferenceValue(cx, xref, xrefval);

	// yaxis scale(s), ref value
	for ( int i = 0; i < yAxisCount; i++ ){
		PMS::Axis y = cacheParams->yAxis( i );
		if (y==PMS::NONE) {
			String caltype = itsCache_->calType();
			if (caltype.startsWith("Xf")) {
				y = PMS::GPHASE;
			} else if (caltype == "GSPLINE") {
				y = getGsplineAxis(dataParams->filename());
			} else {
				y = PMS::DEFAULT_YAXIS;
			}
			cacheParams->setYAxis(y, i);
		}
		// yaxis scale
		PlotAxis cy = axesParams->yAxis( i );
		canvas->setAxisScale(cy, PMS::axisScale(y));
		// yaxis ref value
		bool yref = itsCache_->hasReferenceValue(y);
		double yrefval = itsCache_->referenceValue(y);
		canvas->setAxisReferenceValue(cy, yref, yrefval);
	}

	// if shown, set axis fonts
	casacore::Int pointsize;
	if (set && showX) {
		PlotFontPtr xFont = canvas->axisFont(cx);
		pointsize = (canvParams->xFontSet()) ? canvParams->xAxisFont(): std::max(12.-numplots+1., 8.);
		xFont->setPointSize(pointsize);
		canvas->setAxisFont(cx, xFont);
	}
	if (set && showY) {
		pointsize = (canvParams->yFontSet()) ? canvParams->yAxisFont(): std::max(12.-numplots+1., 8.);
		PlotFontPtr yFont = canvas->axisFont(Y_LEFT);
		yFont->setPointSize(pointsize);
		canvas->setAxisFont(Y_LEFT, yFont);
		yFont = canvas->axisFont(Y_RIGHT);
		yFont->setPointSize(pointsize);
		canvas->setAxisFont(Y_RIGHT, yFont);
	}

	// x and y axis ranges
	canvas->setAxesAutoRescale(true);
	bool makeSquare(false), waveplot(false);  // true if uv/uvwave plot
	if (set) {
		double xmin, xmax, ymin, ymax, xymax = 0;

		bool displayUnflagged = (displayParams->unflaggedSymbol()->symbol() != PlotSymbol::NOSYMBOL);
		bool displayFlagged = (displayParams->flaggedSymbol()->symbol() != PlotSymbol::NOSYMBOL);
		if (displayUnflagged && !displayFlagged) {        // get range of unflagged data only
			itsCache_->indexer(0,iteration).unmaskedMinsMaxesRaw(xmin, xmax, ymin, ymax);
		} else if (displayFlagged && !displayUnflagged) { // get range of flagged data only
			itsCache_->indexer(0,iteration).maskedMinsMaxesRaw(xmin, xmax, ymin, ymax);
		} else {                                          // get range of all data
			itsCache_->indexer(0,iteration).minsMaxes(xmin, xmax, ymin, ymax);
		}
		bool xPtsToPlot(xmin != DBL_MAX), yPtsToPlot(ymin != DBL_MAX);

		// x range
		bool xIsUV(false), xIsUVwave(false);
		if ( axesParams->xRangeSet() ){
			// Custom axes ranges set by user
			canvas->setAxisRange(cx, axesParams->xRange());
		} else if (xPtsToPlot) {
			setAxisRange(x, cx, xmin, xmax, canvas);
			if (PMS::axisIsUV(x)) {
				xIsUV = true;
				if (x==PMS::UWAVE || x==PMS::VWAVE)
					xIsUVwave = true;
				xymax = canvas->axisRange(cx).second;  // should be equal
			}
		}
		// y range
		for ( int i = 0; i < yAxisCount; i++ ){
			PlotAxis cy = axesParams->yAxis( i );
			if ( axesParams->yRangeSet(i) ){
				// Custom axes ranges set by user
				canvas->setAxisRange(cy, axesParams->yRange(i));
			} else if (yPtsToPlot) {
				PMS::Axis y = cacheParams->yAxis(i);
				// add margin if showAtm so overlay doesn't overlap plot
				if ((cacheParams->showAtm() && y!=PMS::ATM) ||
					(cacheParams->showTsky() && y!=PMS::TSKY)) {
					ymax += (ymax-ymin)*0.5;
					pair<double, double> ybounds = make_pair(ymin, ymax);
					canvas->setAxisRange(cy, ybounds);
				}
				setAxisRange(y, cy, ymin, ymax, canvas);
				if (PMS::axisIsUV(y) && xIsUV) {
					// set x and y ranges equally
					double ymax = canvas->axisRange(cy).first;
					xymax = max(xymax, ymax);
					pair<double, double> xybounds = make_pair(-xymax, xymax);
					canvas->setAxisRange(cx, xybounds);
					canvas->setAxisRange(cy, xybounds);
					makeSquare = true;
					if (xIsUVwave && (y==PMS::UWAVE || y==PMS::VWAVE))
						waveplot=true;
                } else if (y==PMS::ATM || y==PMS::TSKY) {
                    itsCache_->indexer(1,iteration).minsMaxes(xmin, xmax, ymin, ymax);
                    pair<double,double> atmrange;
                    if (y==PMS::ATM) atmrange = make_pair(0, min(ymax+1.0, 100.0));
                    else atmrange = make_pair(0, ymax+0.1);
                    canvas->setAxisRange(cy, atmrange);
                }
			}
		}
	}
	itsParent_->getPlotter()->makeSquarePlot(makeSquare, waveplot);

	// For title and axis labels, need all plots on this canvas
	int gridRow(iterParams->getGridRow());
	int gridCol(iterParams->getGridCol());
	QList<PlotMSPlot*> canvasPlots = itsParent_->getPlotManager().getCanvasPlots(gridRow, gridCol);
	int canvasPlotCount = canvasPlots.size();

	// determine which are MS and which are CalTable cache types
	// Needed for cal axes, axis labels, title
	casacore::Vector<casacore::Int> cacheTypes(canvasPlotCount, 0);   // default MS
	casacore::Vector<casacore::String> calTypes(canvasPlotCount, ""); // default no caltype
	casacore::Int calTableType = PlotMSCacheBase::CAL;
	for (int i=0; i<canvasPlotCount; ++i) {
		PlotMSPlotParameters plotParams = canvasPlots[i]->parameters();
		PMS_PP_MSData* dataParams = plotParams.typedGroup<PMS_PP_MSData>();
		PMS_PP_Cache* cacheParams = plotParams.typedGroup<PMS_PP_Cache>();
		if (dataParams==NULL || cacheParams==NULL)
			continue;
		cacheTypes(i) = dataParams->type();
		if (cacheTypes(i)==calTableType) {
			casacore::String filename = dataParams->filename();
			NewCalTable ct(NewCalTable::createCT(filename, Table::Old, Table::Plain));
			calTypes(i) = ct.tableInfo().subType();
		}
	}

	// needed for title and axis labels
	bool allCalTables = allEQ(cacheTypes, calTableType);
	bool polnRatio = itsCache_->polnRatio();
	PlotMSAveraging averaging = dataParams->averaging();

	if (set) {
		PMS::DataColumn xDataColumn(itsCache_->getXDataColumn());
		vector<PMS::Axis> yAxes;
		vector<bool> yRefs;
		vector<double> yRefVals;
		vector<PMS::DataColumn> yDatas;
		// x-axis label 
		if (showX) {
			if (allCalTables && PMS::axisIsData(x)) // convert xaxis to cal axis depending on type (e.g. "Amp"->"Tsys")
				x = getCalAxis(calTypes(0), x);
			// data col may have been changed during loading if no col
			casacore::String xLabelSingle = canvParams->xLabelFormat().getLabel(x, xref, xrefval, xDataColumn, polnRatio);
			if (x==PMS::TIME && xLabelSingle.contains("1858")) // xrefval==0
				xLabelSingle.gsub("(from 1858/11/17)", "");
			if (x == PMS::FREQUENCY)
				xLabelSingle = addFreqFrame(xLabelSingle);
			if (axisIsAveraged(x, averaging) && !allCalTables)
				xLabelSingle = "Average " + xLabelSingle;
			if (allCalTables && xLabelSingle.contains("Corr")) 
				xLabelSingle.gsub("Corr", "Poln");
			canvas->setAxisLabel(cx, xLabelSingle);
		}
		// y-axis label(s)
		if(showY) {
			casacore::String yLabelLeft(""), yLabelRight("");
			for ( int i=0; i<canvasPlotCount; i++ ){
				PlotMSPlotParameters plotParams = canvasPlots[i]->parameters();
				PMS_PP_Cache *plotCacheParams = plotParams.typedGroup<PMS_PP_Cache>();
				PMS_PP_Axes * plotAxisParams = plotParams.typedGroup<PMS_PP_Axes>();
				if ( plotCacheParams == NULL || plotAxisParams == NULL )
					continue;
				PlotMSCacheBase& plotCacheBase = canvasPlots[i]->cache();
				bool isCalTable(cacheTypes(i)==calTableType);
				int plotYAxisCount = plotAxisParams->numYAxes();
				for ( int j=0; j<plotYAxisCount; j++ ){
					PMS::Axis y = plotCacheParams->yAxis( j );
					if (isCalTable && PMS::axisIsData(y))
						y = getCalAxis(calTypes(i), y);
					yAxes.push_back(y);  // save for title
					PlotAxis cy = plotAxisParams->yAxis( j );
					bool yref = plotCacheBase.hasReferenceValue(y);
					yRefs.push_back(yref); // save for title
					double yrefval = plotCacheBase.referenceValue(y);
					yRefVals.push_back(yrefval); // save for title
					// data col may have been changed during loading if no col
					PMS::DataColumn yDataColumn = plotCacheBase.getYDataColumn(j);
					yDatas.push_back(yDataColumn); // save for title
					casacore::String yLabelSingle = canvParams->yLabelFormat( ).getLabel(y, yref, yrefval, yDataColumn, polnRatio);
					if (y==PMS::TIME && yLabelSingle.contains("1858")) // yrefval==0
						yLabelSingle.gsub("(from 1858/11/17)", "");
					if (y == PMS::FREQUENCY)
						yLabelSingle = addFreqFrame(yLabelSingle);
					if (axisIsAveraged(y, averaging) && !isCalTable)
						yLabelSingle = "Average " + yLabelSingle;
					if (isCalTable && yLabelSingle.contains("Corr"))
						yLabelSingle.gsub("Corr", "Poln");
					if ( cy == Y_LEFT ){
						if ( yLabelLeft.size() > 0 )
							yLabelLeft.append( ", ");
						yLabelLeft.append( yLabelSingle );
					} else {
						if ( yLabelRight.size() > 0 )
							yLabelRight.append( ", ");
						yLabelRight.append( yLabelSingle );
					}
				}
			}
			canvas->setAxisLabel(Y_LEFT, yLabelLeft);
			canvas->setAxisLabel(Y_RIGHT, yLabelRight);
		}

		// Title font
		PlotFontPtr font = canvas->titleFont();
		pointsize = (canvParams->titleFontSet()) ? canvParams->titleFont() :
			std::max(16.-numplots+1., 8.);
		font->setPointSize(pointsize);
		font->setBold(true);
		canvas->setTitleFont(font);
		// Title text
		casacore::String iterTxt("");
		if((iterParams->iterationAxis()!=PMS::NONE) && itsCache_->nIter(0) > 0) 
			iterTxt = itsCache_->indexer(0,iteration).iterLabel();
		casacore::String title = canvParams->titleFormat().getLabel(x, yAxes, xref, xrefval,
				yRefs, yRefVals, xDataColumn, yDatas, polnRatio) + " " + iterTxt;
		// change "Corr" ->"Poln" for cal tables:
		if (title.contains("Corr") && anyEQ(cacheTypes, calTableType)) {
			if (allCalTables) {
				title.gsub("Corr", "Poln");
			} else {  // mixed MS/CT
				// change Corr->Pol for CT yaxis only
				if (title.startsWith("Corr") && cacheTypes(0)==calTableType)
					title.replace(0, 4, "Poln");
				else if (title.contains(", Corr") && cacheTypes(1)==calTableType)
					title.gsub(", Corr", ", Poln");
			}
		}
		canvas->setTitle(title);
	}

	// Legend
	canvas->showLegend(set && canvParams->legendShown(), canvParams->legendPosition());

	// Grid lines
	canvas->showGrid(canvParams->gridMajorShown(), canvParams->gridMinorShown(),
		canvParams->gridMajorShown(), canvParams->gridMinorShown());
	// major
	PlotLinePtr major_line = itsFactory_->line(canvParams->gridMajorLine());
	if (!canvParams->gridMajorShown())
		major_line->setStyle(PlotLine::NOLINE);
	canvas->setGridMajorLine(major_line);
	// minor
	PlotLinePtr minor_line = itsFactory_->line(canvParams->gridMinorLine());
	if (!canvParams->gridMinorShown()) 
		minor_line->setStyle(PlotLine::NOLINE);
	canvas->setGridMinorLine(minor_line);
}

void PlotMSPlot::setAxisRange(PMS::Axis axis, PlotAxis paxis, 
		double minval, double maxval, PlotCanvasPtr& canvas) {
	pair<double, double> bounds;

	// don't override larger axis range 
	bool rangeSet(canvas->numPlots() > 0);  // already a plot!
	double canvRangeMin, canvRangeMax;
	if (rangeSet) { // get range and check for default
		canvRangeMin = canvas->axisRange(paxis).first;
	    canvRangeMax = canvas->axisRange(paxis).second;
		rangeSet &= ((canvRangeMin != 0.0) && (canvRangeMax != 1000.0));
	}
	if (rangeSet) {  // possibly by user on first plot
		minval = min(minval, canvas->axisRange(paxis).first);
		maxval = max(maxval, canvas->axisRange(paxis).second);
	}

	// CAS-3263 points near zero are not plotted, so add lower margin
	if ((minval > -0.5) && (minval < 1.0) && (maxval > 10.0)) {
		if (maxval > 100.0) minval -= 1.0; // add larger margin for larger range
		else minval -= 0.1;
		bounds = make_pair(minval, maxval);
		canvas->setAxisRange(paxis, bounds);
	}
	
	if (axis==PMS::TIME) {
		// explicitly set range so can set time scale 
		double diff = maxval - minval;
		if (diff>120.0) {  // seconds (2 minutes)
			bounds = make_pair(minval, maxval);
			canvas->setAxisRange(paxis, bounds);
		} else if (diff==0.0) {
			// override autoscale which sets crazy tick marks;
			// add 2-sec margins
			bounds = make_pair(minval-2.0, maxval+2.0);
			canvas->setAxisRange(paxis, bounds);
		}
	} else if (PMS::axisIsUV(axis)) {
		// make range symmetrical for uv plot
		if ((minval != DBL_MAX) && (maxval != -DBL_MAX)) {
			double maximum = round(max(abs(minval),maxval)) + 10.0;
			minval = -maximum;
			maxval = maximum;
			bounds = make_pair(minval, maxval);
			canvas->setAxisRange(paxis, bounds);
		}
	}
}

bool PlotMSPlot::axisIsAveraged(PMS::Axis axis, PlotMSAveraging averaging) {
    bool avgAxis = false;
    switch (axis) {
        case PMS::TIME:
            if (averaging.time()) avgAxis = true;
            break;
        case PMS::CHANNEL:
            if (averaging.channel()) avgAxis = true;
            break;
        case PMS::BASELINE:
            if (averaging.baseline()) avgAxis = true;
            break;
        default:
            break;
    }
    return avgAxis;
}

String PlotMSPlot::addFreqFrame(String freqLabel) {
    if (itsCache_->cacheType() == PlotMSCacheBase::MS) {
        String freqType = MFrequency::showType(itsCache_->getFreqFrame());
        return freqLabel + " " + freqType;
    } else {
        return freqLabel;
    }
}

PMS::Axis PlotMSPlot::getCalAxis(String calType, PMS::Axis axis) {
    if (axis==PMS::AMP) {
        if (calType.contains("TSYS")) return PMS::TSYS;
        if (calType.contains("SWPOW")) return PMS::SWP;
        if (calType.contains("Opac")) return PMS::OPAC;
        if (calType.contains("SD")) return PMS::GREAL;
        if (calType[0]=='F') return PMS::TEC;
		if (calType.startsWith("KAntPos")) return PMS::ANTPOS;
        if (calType[0]=='K') return PMS::DELAY;
        return PMS::GAMP;
    }
    if (axis==PMS::PHASE) return PMS::GPHASE;
    if (axis==PMS::REAL) return PMS::GREAL;
    if (axis==PMS::IMAG) return PMS::GIMAG;
    return axis;
}

}
