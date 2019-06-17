//# QtViewer.cc: Qt implementation of main viewer supervisory object
//#		 -- Gui level.
//# Copyright (C) 2005
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
//# $Id$
#if ! defined(WITHOUT_DBUS)
#include <display/QtViewer/QtDBusViewerAdaptor.qo.h>
#else
#include <display/QtViewer/grpcViewerAdaptor.qo.h>
#include <grpc++/grpc++.h>
#include <casagrpc/protos/registrar.grpc.pb.h>
#include <casagrpc/protos/img.grpc.pb.h>

using casatools::rpc::Registrar;

#include <unistd.h>
#endif
#include <display/QtViewer/QtViewer.qo.h>
#include <display/QtViewer/QtDisplayPanelGui.qo.h>
#include <display/QtViewer/QtCleanPanelGui.qo.h>
#include <display/QtViewer/QtCleanPanelGui2.qo.h>

extern int qInitResources_QtViewer();

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

#if defined(WITHOUT_DBUS)
	inline std::string to_string(const QString &other) { return std::string((const char*) other.toLatin1().data()); }
#endif
	QString QtViewer::name_;

	const QString &QtViewer::name( ) {
		return name_;
	}

	QtViewer::QtViewer( const std::list<std::string> &args, bool is_server, const char *server_string ) :
		QtViewerBase(is_server),
#if ! defined(WITHOUT_DBUS)
		dbus_(NULL),
#endif
		args_(args), is_server_(is_server) {

		name_ = (is_server_ ? "view_server" : "viewer");
		server_string_ = (server_string ? server_string : 0);

		qInitResources_QtViewer();
		// Makes QtViewer icons, etc. available via Qt resource system.
		//
		// You would normally use this macro for the purpose instead:
		//
		//   Q_INIT_RESOURCE(QtViewer);
		//
		// It translates as:
		//
		//   extern int qInitResources_QtViewer();
		//   qInitResources_QtViewer();
		//
		// It doesn't work here because it makes the linker looks for
		//   casa::qInitResources_QtViewer()     :-)   dk

		if ( is_server_ ) {
#if ! defined(WITHOUT_DBUS)
			dbus_ = new QtDBusViewerAdaptor(this);
			dbus_->connectToDBus(server_string_);
		} else {
			dbus_ = 0;
#else
			static const auto debug = getenv("GRPC_DEBUG");

			grpcViewerState *state = new grpcViewerState(this);

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
			// image viewer service (and currently interactive clean service though this needs
			// to eventually move to a seperate grpc service description which could e.g. be
			// shared with carta
			auto viewer_svc = state->viewer_service.get( );
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// all gui operations must happen in the "gui thread" because Qt is not
			// thread-safe... so we need create & result signals and slots
			connect( viewer_svc, SIGNAL(new_op( )),
					 this, SLOT(grpc_handle_op( )) );
			connect( viewer_svc, SIGNAL(exit_now( )),
					 this, SLOT(quit( )) );
			builder.RegisterService(viewer_svc);
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// shutdown service is used by casatools etc. to notify gui services
			// when the system is shutting down...
			auto shutdown_svc = state->shutdown_service.get( );
			builder.RegisterService(shutdown_svc);
			connect( shutdown_svc, SIGNAL(exit_now( )),
					 this, SLOT(quit( )) );


			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// Launch server...
			state->server = builder.BuildAndStart( );
			if ( selected_port > 0 ) {
				// if an available port can be found, selected_port is set to a value greater than zero
				snprintf(address_buf,sizeof(address_buf),address_template,selected_port);
				state->uri = address_buf;
				if ( debug ) {
					std::cout << "viewer service available at " << state->uri << std::endl;
					fflush(stdout);
				}
				if ( server_string ) {
					if ( access( server_string, F_OK ) == -1 ) {
						// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
						// create the connection to the registrar for service registration using the uri
						// provided on the command line...
						std::unique_ptr<Registrar::Stub> proxy =
							Registrar::NewStub(grpc::CreateChannel(server_string, grpc::InsecureChannelCredentials( )));
						// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
						// register our "shutdown, "image-view" and "interactive-clean" services with
						// the registrar...
						casatools::rpc::ServiceId sid;
						sid.set_id("casaviewer");
						sid.set_uri(state->uri);
						sid.add_types("shutdown");
						sid.add_types("image-view");
						sid.add_types("interactive-clean");
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
					} else {
						int fd;
						if ( (fd = open(server_string, O_WRONLY)) != -1 ) {
							char uri[strlen(state->uri.c_str( ))+2];
							sprintf(uri,"%s\n",state->uri.c_str( ));
							if ( write( fd, uri, strlen(uri)) != strlen(uri) ) {
								qWarning("server failed to write gRPC URI to named pipe...");
								qFatal("exiting...");
								exit(1);
							}
							close(fd);
						} else {
							qWarning("server failed to open gRPC URI named pipe...");
							qFatal("exiting...");
							exit(1);
						}
					}

					// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
					// complete startup
					grpc_.reset(state);
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

	}

	QtViewer::~QtViewer() {
		// wonder if we need to delete dbus adaptor...
	}

	QtDisplayPanelGui *QtViewer::createDPG() {
		// Create a main display panel Gui.
		//
		QtDisplayPanelGui* dpg = new QtDisplayPanelGui(this,0,"dpg",args_);
		// don't need to worry about QtCleanPanelGui objs because they
		// do not need to support cursor tracking and the P/V tool...
		// although, there is nothing preventing QtCleanPanelGui objs
		// being added to this list...
		panels.push_back(dpg);

		// Previously casaviewer.cc created a QtDisplayPanelGui directly,
		// now it uses this function to ensure consistent behavior with
		// creation of QtDisplayPanelGui objects...
		if ( is_server_ ) {

			dpg->setAttribute(Qt::WA_DeleteOnClose);
			// Deletes this window (only) when user closes it.
			// (Note that 'closed' Qt windows are not deleted by default --
			// not unless this attribute is explicitly set).
			//
			// Noone will hold the dpg pointer (! -- at least in the current
			// revision).  Essentially, the gui user manages dpg's storage....
			//
			// Note, however: DPGs do not have to be created via this routine.
			// QtClean, e.g., constructs its own dpg directly, which does
			// not delete on close, and which is re-opened on successive
			// restarts of the Qt event loop.  In that case, QtClean manages
			// its own dpg storage.

			dpg->show();
		}

		return dpg;
	}

	QtDisplayPanelGui *QtViewer::createDPG(const QtDisplayPanelGui *other) {
		// Create a main display panel Gui.
		//
		QtDisplayPanelGui* dpg = new QtDisplayPanelGui(other,0,"dpg",args_);
		// don't need to worry about QtCleanPanelGui objs because they
		// do not need to support cursor tracking and the P/V tool...
		// although, there is nothing preventing QtCleanPanelGui objs
		// being added to this list...
		panels.push_back(dpg);

		// Previously casaviewer.cc created a QtDisplayPanelGui directly,
		// now it uses this function to ensure consistent behavior with
		// creation of QtDisplayPanelGui objects...
		if ( is_server_ ) {

			dpg->setAttribute(Qt::WA_DeleteOnClose);
			// Deletes this window (only) when user closes it.
			// (Note that 'closed' Qt windows are not deleted by default --
			// not unless this attribute is explicitly set).
			//
			// Noone will hold the dpg pointer (! -- at least in the current
			// revision).  Essentially, the gui user manages dpg's storage....
			//
			// Note, however: DPGs do not have to be created via this routine.
			// QtClean, e.g., constructs its own dpg directly, which does
			// not delete on close, and which is re-opened on successive
			// restarts of the Qt event loop.  In that case, QtClean manages
			// its own dpg storage.

			dpg->show();
		}

		return dpg;
	}

	void QtViewer::dpgDeleted( QtDisplayPanelGui *dpg ) {
		panel_list_t::iterator iter = std::find(panels.begin( ), panels.end( ), dpg );
		if ( iter != panels.end( ) ) panels.erase(iter);
	}

	QtCleanPanelGui *QtViewer::createInteractiveCleanGui( ) {
		QtCleanPanelGui* cpg = new QtCleanPanelGui(this,0,args_);
		return cpg;
	}

	QtCleanPanelGui2 *QtViewer::createInteractiveCleanGui2( ) {
		QtCleanPanelGui2* cpg = new QtCleanPanelGui2(this,0,args_);
		return cpg;
	}

	void QtViewer::activate( bool state ) {
		for ( panel_list_t::iterator iter = panels.begin( ); iter != panels.end( ); ++iter )
			(*iter)->activate( state );
	}

	void QtViewer::quit() {
#if defined(WITHOUT_DBUS)
		static const auto debug = getenv("GRPC_DEBUG");
		if ( grpc_ && grpc_->server ) {
			if ( debug ) {
				std::cout << "entering QtViewer::quit( )..." << std::endl;
				std::cout << "		  ...shutting down grpc server..." << std::endl;
				fflush(stdout);
			}
			grpc_->server->Shutdown( );
		}
		
		if ( debug && grpc_->server ) {
			std::cout << "		  ...shutting down qt..." << std::endl;
			fflush(stdout);
		}

#endif
		QtViewerBase::quit(); 
#if defined(WITHOUT_DBUS)
		if ( grpc_ && grpc_->server ) {
			QCoreApplication::exit( );
			// calling the system exit( ) here causes immediate
			// shutdown, but does not allow global cleanup...
		}
#endif
	}

#if defined(WITHOUT_DBUS)

	void QtViewer::grpc_handle_op( ) {
		std::lock_guard<std::mutex> exc(grpc_queue_mutex);
		if ( ! grpc_queue.empty( ) ) {
			std::function<void()> f = grpc_queue.front( );
			grpc_queue.pop( );
			f( );
		}
	}

#endif


} //# NAMESPACE CASA - END
