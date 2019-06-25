//# grpcPlotMSAdaptor.h: provides plotms services via grpc
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

#include <thread>
#include <future>
#include <unistd.h>
#include <plotms/PlotMS/grpcPlotMSAdaptor.qo.h>
#include <plotms/PlotMS/PlotMSParameters.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

namespace casa {

    constexpr char grpcPlotMS::APP_SERVER_SWITCH[];

    ::grpc::Status grpcPMSShutdown::now(::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty*) {
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received shutdown notification..." << std::endl;
            fflush(stdout);
        }
        if ( itsPlotms_ ) {
            itsPlotms_ = NULL;
            static auto bye_bye = std::async( std::launch::async, [=]( ) { sleep(2); emit exit_now( ); } );
        }
        return grpc::Status::OK;
    }

    grpcPlotMS::grpcPlotMS( PlotEngine *pe ) : itsPlotms_(pe), plotter_(0) { }
    void grpcPlotMS::set_plotter( PlotMSPlotter *plt ) { plotter_ = plt; }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotMS::update_parameters( int index ) {
        PlotMSPlotParameters* sp = itsPlotms_->getPlotManager().plotParameters(index);
        if ( sp == NULL ) {
            fprintf(stderr, "internal inconsistency, plot parameter %d not available\n", index);
            return;
        }
        sp->holdNotification( );
        for ( auto iter = param_groups.begin( ); iter != param_groups.end( ); ++iter ) {
            if ( iter->first.first != index ) continue;
            switch ( iter->first.second ) {
                case T_MSDATA:
                    *sys_ppdata(sp) = *((PMS_PP_MSData*)iter->second);
            }
        }
        if ( sp != NULL) sp->releaseNotification();
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    PMS_PP_MSData *grpcPlotMS::sys_ppdata(PlotMSPlotParameters* sp) {
        PMS_PP_MSData *plotdata  = sp->typedGroup<PMS_PP_MSData>( );
        if ( plotdata == NULL) {
            sp->setGroup<PMS_PP_MSData>( );
            plotdata  = sp->typedGroup<PMS_PP_MSData>( );
        }
        return plotdata;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    PMS_PP_MSData *grpcPlotMS::ppdata(int index) {
        std::pair<int,int> idx(index,T_MSDATA);
        auto ptr = (PMS_PP_MSData*) param_groups[idx];
        if ( ptr == NULL ) {
            param_groups[idx] = ptr = new PMS_PP_MSData(itsPlotms_->getPlotFactory( ));
            PlotMSPlotParameters *sp = itsPlotms_->getPlotManager().plotParameters(index);
            sp->holdNotification( );
            *ptr = *sys_ppdata(sp);
            sp->releaseNotification();
        }
        return (PMS_PP_MSData*) param_groups[idx];
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotMS::qtGO( std::function<void()> func ) {
        std::lock_guard<std::mutex> exc(plotter_->grpc_queue_mutex);
        plotter_->grpc_queue.push(func);
        emit new_op( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::getPlotMSPid( ::grpc::ServerContext *context,
                                             const ::google::protobuf::Empty*,
                                             ::rpc::plotms::Pid *reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cout << "received getPlotMSPid( ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        //
        // the DBus version of plotms does some dance with getting and *setting* the pid...
        // I doubt that's necessary with grpc...
        reply->set_id(getpid( ));
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::setShowGui( ::grpc::ServerContext *context,
                                           const ::rpc::plotms::Toggle *req,
                                           ::google::protobuf::Empty* ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cout << "received setShowGui( " << req->state( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        qtGO( [=]( ){ itsPlotms_->showGUI(req->state( )); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::setGridSize( ::grpc::ServerContext *context,
                                            const ::rpc::plotms::GridSize *req,
                                            ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cout << "received setGridSize( " << req->rows( ) << ", " << req->cols( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        auto rows = req->rows( );
        auto cols = req->cols( );
		qtGO( [=]( ){ itsPlotms_->getParameters().setGridSize( rows, cols ); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::isDrawing( ::grpc::ServerContext *context,
                                          const ::google::protobuf::Empty*,
                                          ::rpc::plotms::Toggle *reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cout << "received isDrawing( ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        //
        // the DBus version of plotms does some dance with getting and *setting* the pid...
        // I doubt that's necessary with grpc...
        bool isdrawing = itsPlotms_->isDrawing( );
        if (debug) {
            std::cout << "\t\tresult: " << isdrawing << std::endl;
            fflush(stdout);
        }
        reply->set_state(isdrawing);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::clearPlots( ::grpc::ServerContext *context,
                                           const ::google::protobuf::Empty*,
                                           ::google::protobuf::Empty* ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cout << "received clearPlots( ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
		qtGO( [=]( ){ itsPlotms_->clearPlots( ); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotMS::setPlotMSFilename( ::grpc::ServerContext *context,
                                                  const ::rpc::plotms::SetVis *req,
                                                  ::google::protobuf::Empty* ) {
        int index = req->index( );
        if ( invalid_index(index) ) return grpc::Status(grpc::StatusCode::OUT_OF_RANGE, "index out of range");

        bool update = req->update( );
        std::string name = req->name( );

        ppdata(index)->setFilename(name);
        if ( update ) qtGO( [=]( ) { update_parameters(index); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    bool grpcPlotMS::invalid_index(int index) {
        return index < 0 || index >= itsPlotms_->getPlotManager().numPlots( );
    }

}
