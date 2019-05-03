//# grpcViewerAdaptor.h: provides viewer services via grpc
//# Copyright (C) 2019
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

#include <future>
#include <display/QtViewer/QtViewer.qo.h>
#include <display/QtViewer/QtCleanPanelGui.qo.h>
#include <display/QtViewer/QtCleanPanelGui2.qo.h>
#include <display/QtViewer/grpcViewerAdaptor.qo.h>

namespace casa { //# NAMESPACE CASA - BEGIN

    grpcShutdown::grpcShutdown( QtViewer *viewer ) : viewer_(viewer) {	}

    ::grpc::Status grpcShutdown::now( ::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty* ) {
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received shutdown notification..." << std::endl;
            fflush(stdout);
        }
        static auto bye_bye = std::async( std::launch::async, [&]( ) { sleep(2); emit exitnow( ); } );
        return grpc::Status::OK;
    }

	QtDisplayPanelGui *grpcImageViewer::findpanel( int key /*, bool create*/ ) {

		if ( key == 0 ) key = INT_MAX;

		if ( managed_panels.find( key ) != managed_panels.end( ) )
			return managed_panels.find( key )->second->panel( );
        /* reenabling this will require generating a signal to cause
         * QtViewer object to create a panel, and wait for the result
         * signal to be recieved, so that the created panel pointer
         * can be returned here...
		if ( key == INT_MAX && create ) {
			QtDisplayPanelGui *dpg = create_panel( );
			managed_panels.insert(panelmap::value_type(INT_MAX, new panel_desc(dpg)));
			return dpg;
		}
        */
		return 0;
	}

    int grpcImageViewer::get_id( ) {
		int index = QtId::get_id( );
		managed_panels.insert(panelmap::value_type(index, new panel_desc( )));
		return index;
    }

    void grpcImageViewer::panel_result( QtDisplayPanelGui *panel, int id ) {
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received qt panel_result( " << id << " ) signal..." << std::endl;
            fflush(stdout);
        }
        std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);
        managed_panels[id]->panel( ) = panel;
    }

	void grpcImageViewer::handle_destroyed_panel( QObject *panel ) {
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received qt handle_destroyed_panel( ) signal..." << std::endl;
            fflush(stdout);
        }
        std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);
        std::lock_guard<std::recursive_mutex> d_guard(managed_datas_mutex);
		for ( panelmap::iterator iter = managed_panels.begin();
		        iter != managed_panels.end(); ++iter ) {
			if ( iter->second->panel( ) == panel ) {
                if (getenv("GRPC_DEBUG")) {
                    std::cout << "found destroyed panel..." << std::endl;
                    fflush(stdout);
                }
				for ( std::list<int>::iterator diter = iter->second->data().begin();
				        diter != iter->second->data().end(); ++diter ) {
					datamap::iterator kill = managed_datas.find(*diter);
					if ( kill != managed_datas.end( ) ) {
						delete kill->second;
						managed_datas.erase(kill);
					}
				}
				delete iter->second;
				managed_panels.erase(iter);
				break;
			}
		}
	}


	grpcImageViewer::grpcImageViewer(QtViewer  *viewer) : viewer_(viewer) {	}

    ::grpc::Status grpcImageViewer::panel(::grpc::ServerContext *context, const ::rpc::img::NewPanel *req, ::rpc::img::PanelId *reply) {
        int result = get_id( );
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received grpc panel( " << req->type( ) << " ) event..." << std::endl;
            std::cout << "sending qt panel( " << result << " ) signal..." << std::endl;
            fflush(stdout);
        }
        emit panel( QString(req->type( ).c_str( )), req->hidden( ), result );
        reply->set_id(result);
        return grpc::Status::OK;
    }
            
}
