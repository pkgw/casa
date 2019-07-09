//# PlotMS.cc: Main controller for plotms.
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: $
#include <plotms/PlotMS/PlotMS.h>

#include <plotms/Gui/PlotMSPlotter.qo.h>
#if ! defined(WITHOUT_DBUS)
#include <plotms/PlotMS/PlotMSDBusApp.h>
#else
#include <plotms/PlotMS/grpcPlotMSAdaptor.qo.h>
#include <casagrpc/protos/registrar.grpc.pb.h>

using casatools::rpc::Registrar;

#include <fcntl.h>
#include <unistd.h>
#endif
#include <plotms/Client/ClientFactory.h>
#include <plotms/Actions/ActionFactory.h>
#include <plotms/Actions/ActionCacheLoad.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/PlotMS/PlotMSParameters.h>
#include <plotms/PlotMS/PlotMSExportParam.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

#include <QDebug>

using namespace casacore;
namespace casa {

// TODO PlotMSAction: iteration, release cache.  action for new plots.  update
//      hold/release action/button text.
// TODO PlotMSCache: multi-region locate
// TODO PlotMSLogger: log source (std out, text widget, casapy logger), better
//      locate message, log message for parameters updated and action
//      execution, log event flag for flag/unflag
// TODO PlotMSParameters: canvas background, fonts, grids, spacing, cartesian
//      axes, limit zoom/pan cache size
// TODO PlotMSPlot: different colors within one plot, different types, shared
//      caches
// TODO PlotMSPlotter: range padding, customize toolbars/tabs
// TODO PlotMSThread: background, pause/resume
// TODO PlotMSWidgets: label creator
// TODO PlotTool: set tracker font

////////////////////////
// PLOTMS DEFINITIONS //
////////////////////////

// Constructors/Destructors //

PlotMSApp::PlotMSApp(bool connectToDBus, bool userGui
#if defined(WITHOUT_DBUS)
                        , const casacore::String &casapy_address
#endif
                     ) :
		itsLastPlotter_(NULL), isGUI_(userGui), allow_popups(true), 
		itsExportFormat( PlotExportFormat::JPG, "")
#if ! defined(WITHOUT_DBUS)
		, itsDBus_(NULL)
#endif
{
	initialize(connectToDBus, userGui
#if defined(WITHOUT_DBUS)
                        , casapy_address
#endif
               );
}

PlotMSApp::PlotMSApp(const PlotMSParameters& params, bool connectToDBus, bool userGui
#if defined(WITHOUT_DBUS)
                        , const casacore::String &casapy_address
#endif
                     ) :
		itsPlotter_(NULL), itsLastPlotter_(NULL), isGUI_(userGui), 
		allow_popups(true), itsParameters_(params),
		itsExportFormat( PlotExportFormat::JPG, "")
#if ! defined(WITHOUT_DBUS)
        , itsDBus_(NULL)
#endif
{
	initialize(connectToDBus, userGui
#if defined(WITHOUT_DBUS)
                        , casapy_address
#endif
               );
}

PlotMSApp::~PlotMSApp() {
#if ! defined(WITHOUT_DBUS)
    if(itsDBus_ != NULL) delete itsDBus_;
#endif
}


// Public Methods //

//PlotMSPlotter* PlotMSApp::getPlotter() { return itsPlotter_; }
bool PlotMSApp::guiShown() const { return itsPlotter_->guiShown(); }
int PlotMSApp::execLoop() { return itsPlotter_->execLoop(); }
int PlotMSApp::showAndExec(bool show) { return itsPlotter_->showAndExec(show); }
void PlotMSApp::close() { itsPlotter_->close(); }

void PlotMSApp::showError(const String& message, const String& title,
        bool isWarning) { itsPlotter_->showError(message, title, isWarning); }
void PlotMSApp::showWarning(const String& message, const String& title) {
    itsPlotter_->showError(message, title, true); }
void PlotMSApp::showMessage(const String& message, const String& title) {
    itsPlotter_->showMessage(message, title);
}
void PlotMSApp::clearMessage() { itsPlotter_->clearMessage(); }

void PlotMSApp::showGUI(bool show) {
    if (show != isGUI_) {
        close();
        isGUI_ = show; 
        initialize(false, show);  // DBus has already been taken care of 
    } 
    itsPlotter_->showGUI(show);
}

PlotMSParameters& PlotMSApp::getParameters() {
	return itsParameters_;
}

void PlotMSApp::setParameters(const PlotMSParameters& params) {
    itsParameters_ = params;
    parametersHaveChanged(itsParameters_,
                PlotMSWatchedParameters::ALL_UPDATE_FLAGS());
}

PlotMSExportParam& PlotMSApp::getExportParameters() {
	return itsExportParameters_;
}
void PlotMSApp::setExportParameters(const PlotMSExportParam& params) {
    itsExportParameters_ = params;
}

PlotExportFormat PlotMSApp::getExportFormat(){
	return itsExportFormat;
}

void PlotMSApp::setExportFormat( const PlotExportFormat format ){
	itsExportFormat = format;
}

void PlotMSApp::parametersHaveChanged(const PlotMSWatchedParameters& params,
            int updateFlag) {
	(void)updateFlag;
    // We only care about PlotMS's parameters.
    if(&params == &itsParameters_) {
        itsLogger_->setSinkLocation(itsParameters_.logFilename());
        itsLogger_->setFilterEventFlags(itsParameters_.logEvents());
        itsLogger_->setFilterMinPriority(itsParameters_.logPriority());
        
        int rowCount = itsParameters_.getRowCount();
        int colCount = itsParameters_.getColCount();
        bool gridChanged = itsPlotManager_.pageGridChanged( rowCount, colCount, false );
        if ( itsPlotter_ != NULL ){
        	pair<int, int> cis = itsParameters_.cachedImageSize();
        	itsPlotter_->setCanvasCachedAxesStackImageSize( cis.first, cis.second );
        	if ( gridChanged ){
        		itsPlotter_->gridSizeChanged( rowCount, colCount );
        		itsPlotter_->plot();
        	}
        }
    }
}

PlotSymbolPtr PlotMSApp::createSymbol( const PlotSymbolPtr& copy ){
	PlotSymbolPtr symbol;
	if ( itsPlotter_ != NULL ){
		symbol = itsPlotter_->createSymbol( copy );
	}
	return symbol;
}

PlotSymbolPtr PlotMSApp::createSymbol (const String& descriptor,
    		Int size, const String& color,
        	const String& fillPattern, bool outline ){
	PlotSymbolPtr symbol;
	if ( itsPlotter_ != NULL ){
		symbol = itsPlotter_->createSymbol( descriptor, size, color, fillPattern, outline );
	}
	return symbol;
}

PlotLoggerPtr PlotMSApp::getLogger() { return itsLogger_; }
PlotMSPlotManager& PlotMSApp::getPlotManager() { return itsPlotManager_; }


PlotMSPlot* PlotMSApp::addOverPlot(const PlotMSPlotParameters* p) {
    return itsPlotManager_.addOverPlot(p);
}

void PlotMSApp::clearPlots(){
	return itsPlotManager_.clearPlotsAndCanvases();
}

bool PlotMSApp::isDrawing() const {
	return itsPlotter_->isDrawing();
}

void PlotMSApp::setAnnotationModeActive( PlotMSAction::Type type, bool active ){
	itsPlotter_->setAnnotationModeActive( type, active );
}

bool PlotMSApp::isClosed() const {
	return itsPlotter_ == NULL ||
               itsPlotter_->isClosed();
}

bool PlotMSApp::save(const PlotExportFormat& format) {
	return itsPlotter_->exportPlot(format, false);
}

PlotFactoryPtr PlotMSApp::getPlotFactory(){
	return itsPlotter_->getPlotFactory();
}

void PlotMSApp::quitApplication(){
	CountedPtr<PlotMSAction> quitAction = ActionFactory::getAction( PlotMSAction::QUIT, NULL );
	quitAction->doAction( this );
}

PlotMSFlagging PlotMSApp::getFlagging() const {
	return itsPlotter_->getFlagging();
}

void PlotMSApp::setFlagging(PlotMSFlagging flag){
	itsPlotter_->setFlagging( flag );
}

void PlotMSApp::canvasAdded( PlotCanvasPtr canvas ){
	itsPlotter_->canvasAdded( canvas );
}

bool PlotMSApp::exportToFormat(const PlotExportFormat& format){
	bool exportOK = itsPlotter_->exportToFormat( format );
    if (!exportOK && ((format.width > 0) || format.height > 0))
        showError("Export failed! Try decreasing width/height settings.", "Export");
    return exportOK;
}

bool PlotMSApp::isVisible(PlotCanvasPtr& canvas ) {
	return itsPlotter_->isVisible( canvas );
}

Record PlotMSApp::locateInfo( Bool& success, String& errorMessage ){
	CountedPtr<PlotMSAction> action = ActionFactory::getAction(PlotMSAction::SEL_INFO, NULL );
	Record retval;
	success = action->doActionWithResponse( this, retval );
	if(! success ) {
		errorMessage = action->doActionResult();
	}
	return retval;
}

PlotterPtr PlotMSApp::getPlotter(){
	return itsPlotter_->getPlotter();
}

bool PlotMSApp::updateCachePlot( PlotMSPlot* plot,
		void (*f)(void*, bool), bool setupPlot){
	vector<PlotMSPlot*> updatePlots;
	if ( plot != NULL ){
		updatePlots.push_back( plot );
	}
	 ActionCacheLoad loadCacheAction( itsPlotter_, updatePlots, f);
	 bool threading = itsPlotter_->isInteractive();
	 loadCacheAction.setUseThreading( threading );
	 loadCacheAction.setSetupPlot( setupPlot );
	 bool result = loadCacheAction.doAction( this );
	 if ( result ){
		 //So that in fixed mode the plot gets redrawn.
		 plot->setRelease( true );
	 }

	 return result;
}

void PlotMSApp::setCommonAxes(bool commonX, bool commonY ){
	itsPlotter_->setCommonAxes( commonX, commonY );
}

bool PlotMSApp::isCommonAxisX() const {
	return itsPlotter_->isCommonAxisX();
}

void PlotMSApp::resetTools(){
	if ( guiShown()){
		// Reset enabled Tool actions across iteration pages
		PlotMSAction::Type actionType;
		bool enabled;
		CountedPtr<PlotMSAction> toolAction;
		for(int actionInt = PlotMSAction::TOOL_MARK_REGIONS; 
		    actionInt <=  PlotMSAction::TRACKER_ENABLE_DISPLAY; 
		    actionInt++) {
			actionType = static_cast<PlotMSAction::Type>(actionInt); 
			enabled = itsPlotter_->isActionEnabled(actionType);
			if (enabled) {
				toolAction = ActionFactory::getAction( actionType, itsPlotter_ );
				toolAction->doAction( this );
			}
		}
	}
}

bool PlotMSApp::isCommonAxisY() const {
	return itsPlotter_->isCommonAxisY();
}

void PlotMSApp::setAxisLocation( PlotAxis locationX, PlotAxis locationY ){
	itsPlotter_->setAxisLocation( locationX, locationY );
}

PlotAxis PlotMSApp::getAxisLocationX() const {
	return itsPlotter_->getAxisLocationX();
}

PlotAxis PlotMSApp::getAxisLocationY() const {
	return itsPlotter_->getAxisLocationY();
}

bool PlotMSApp::isOperationCompleted() const {
	return operationCompleted;
}

void PlotMSApp::setOperationCompleted( bool completed ){
	operationCompleted = completed;
}


vector<String> PlotMSApp::getFiles() const {
	vector<String> files = itsPlotManager_.getFiles();
	return files;
}

// Private Methods //

void PlotMSApp::initialize(bool connectToDBus, bool userGui
#if defined(WITHOUT_DBUS)
                        , const casacore::String &casapy_address
#endif
                           ) {

	char *server_string = 0;
	std::function<void( )> do_notify = []( ) { };
	operationCompleted = true;
    itsParameters_.addWatcher(this);

    if ((itsLastPlotter_ != NULL) && (itsLastPlotter_->isInteractive() == userGui)) {
        Client* tmpPlotter = itsLastPlotter_;
        itsLastPlotter_ = itsPlotter_;
        itsPlotter_ = tmpPlotter;
        itsPlotter_->clearMessage(); // clear status bar
    } else {
        itsLastPlotter_ = itsPlotter_;
        ClientFactory::ClientType clientType = ClientFactory::GUI;
        if ( !userGui ){
    	    clientType = ClientFactory::SCRIPT;
        }
        itsPlotter_ = ClientFactory::makeClient( clientType, this );
    }
    itsLogger_ = itsPlotter_->getLogger();

    // do this after setting up plotter
    itsPlotManager_.setParent(this);
    
    // Update internal state to reflect parameters.
    parametersHaveChanged(itsParameters_,
            PlotMSWatchedParameters::ALL_UPDATE_FLAGS());
    if(connectToDBus) {
#if ! defined(WITHOUT_DBUS)
        itsDBus_ = new PlotMSDBusApp(*this);
        itsDBus_->connectToDBus();
#else
        auto qplotobj = dynamic_cast<PlotMSPlotter*>(itsPlotter_);
        if ( ! qplotobj ) {
            fprintf(stderr, "cannot start plot server because no Qt plotter is available");
            exit(1);
        }

        static const auto debug = getenv("GRPC_DEBUG");
        grpcPlotMSState *state = new grpcPlotMSState(this);

        //***
        //*** set up a default address (grpc picks port) and address buffers
        //***
        char address_buf[100];
        constexpr char address_template[] = "0.0.0.0:%d";
        snprintf(address_buf,sizeof(address_buf),address_template,0);
        std::string server_address(address_buf);
        int selected_port = 0;

        //***
        //*** build grpc service
        //***
        grpc::ServerBuilder builder;
        // Listen on the given address without any authentication mechanism.
        builder.AddListeningPort(server_address, grpc::InsecureServerCredentials(), &selected_port);
        // Register "service" as the instance through which we'll communicate with
        // clients. In this case it corresponds to an *synchronous* service.
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // plotms service (and currently interactive clean service though this needs
        // to eventually move to a seperate grpc service description which could e.g. be
        // shared with carta
        auto plotms_svc = state->plotms_service.get( );
        builder.RegisterService(plotms_svc);
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // shutdown service is used by casatools etc. to notify gui services
        // when the system is shutting down...
        auto shutdown_svc = state->shutdown_service.get( );
        builder.RegisterService(shutdown_svc);
        QObject::connect( shutdown_svc, SIGNAL(exit_now( )),
                          qplotobj, SLOT(grpc_exit_now( )) );
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // all gui operations must happen in the "gui thread" because Qt is not
        // thread-safe... connect qt events from services to PlotMSPlotter so that
        // grpc events/commands can be executed in the Qt GUI thread...
        plotms_svc->set_plotter(qplotobj);
        QObject::connect( plotms_svc, SIGNAL(new_op( )),
                          qplotobj, SLOT(grpc_handle_op( )) );

        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // Launch server...
        state->server = builder.BuildAndStart( );
        if ( selected_port > 0 ) {
				// if an available port can be found, selected_port is set to a value greater than zero
				snprintf(address_buf,sizeof(address_buf),address_template,selected_port);
				state->uri = address_buf;
				if ( debug ) {
					std::cout << "plotms service available at " << state->uri << std::endl;
					fflush(stdout);
				}

				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				// complete startup
				grpc_.reset(state);

				server_string = casapy_address.size( ) > 0 ? strdup(casapy_address.c_str( )) : 0;
				if ( server_string ) {
					if ( access( server_string, F_OK ) == -1 ) {
						do_notify = [=]( ) {
							// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
							// create the connection to the registrar for service registration using the uri
							// provided on the command line...
							std::unique_ptr<Registrar::Stub> proxy =
								Registrar::NewStub(grpc::CreateChannel(server_string, grpc::InsecureChannelCredentials( )));
							// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
							// register our "shutdown, "image-view" and "interactive-clean" services with
							// the registrar...
							casatools::rpc::ServiceId sid;
							sid.set_id("casaplotms");
							sid.set_uri(state->uri);
							sid.add_types("shutdown");
							sid.add_types("plotms");
							grpc::ClientContext context;
							casatools::rpc::ServiceId accepted_sid;
							if ( debug ) {
								std::cout << "registering services with registrar (at " << server_string << ")" << std::endl;
								fflush(stdout);
							}
							::grpc::Status status = proxy->add(&context,sid,&accepted_sid);
							if ( ! status.ok( ) ) {
								// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
								// if registration was not successful, we exit...
								std::cerr << "registration failed, exiting..." << std::endl;
								fflush(stderr);
								state->server->Shutdown( );
								QCoreApplication::exit(1);
								exit(1);
							}
							if ( debug ) {
								std::cout << "accepted service id: ( " << accepted_sid.id( ) << ", " << accepted_sid.uri( ) << ", ";
								for ( auto i=accepted_sid.types( ).begin( ); i != accepted_sid.types( ).end( ); ++i )
									std::cout << "'" << (*i) << "' ";
								std::cout << ")" << std::endl;
								fflush(stdout);
							}
						};
					} else {
						do_notify = [=]( ) {
							int fd;
							if ( (fd = open(server_string, O_WRONLY)) != -1 ) {
								char uri[strlen(state->uri.c_str( ))+2];
								sprintf(uri,"%s\n",state->uri.c_str( ));
								if ( (size_t) write( fd, uri, strlen(uri)) != strlen(uri) ) {
									qWarning("server failed to write gRPC URI to named pipe...");
									qFatal("exiting...");
									exit(1);
								}
								::close(fd);
							} else {
								qWarning("server failed to open gRPC URI named pipe...");
								qFatal("exiting...");
								exit(1);
							}
						};
					}
				} else {
					grpc_.reset( );
					if ( debug ) {
						std::cout << "no registrar provided... skipped registering." << std::endl;
						fflush(stdout);
					}
				}
			} else grpc_.reset( );

#endif
	}
	showGUI( userGui );
	do_notify( );
	free(server_string);
}



}

