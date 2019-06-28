//# grpcPlotMSAdaptor.h.h: provides viewer services via grpc
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

#ifndef GRPCPLOTMSADAPTOR_H_
#define GRPCPLOTMSADAPTOR_H_
#include <grpc++/grpc++.h>
#include <plotms/PlotMS/PlotEngine.h>
#include <casagrpc/protos/plotms.grpc.pb.h>
#include <casagrpc/protos/shutdown.grpc.pb.h>
#include <plotms/Gui/PlotMSPlotter.qo.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <QObject>
#include <map>

namespace casa { //# NAMESPACE CASA - BEGIN

    class grpcPMSShutdown : public QObject, public ::casatools::rpc::Shutdown::Service {

        Q_OBJECT    //# Allows slot/signal definition.  Must only occur in
                    //# implement/.../*.h files
    public:
        grpcPMSShutdown( PlotEngine *qtv ) : itsPlotms_(qtv) { }
        ::grpc::Status now(::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty*);

    signals:
        void exit_now( );

    private:
        PlotEngine *itsPlotms_;
    };

    class grpcPlotMS : public QObject, public ::rpc::plotms::app::Service {

        Q_OBJECT    //# Allows slot/signal definition.  Must only occur in
                    //# implement/.../*.h files
    public:
        constexpr static char APP_SERVER_SWITCH[] = "--server";
        grpcPlotMS( PlotEngine *pe );
        void set_plotter( PlotMSPlotter *plt );

        ::grpc::Status getPlotMSPid( ::grpc::ServerContext *context,
                                     const ::google::protobuf::Empty*,
                                     ::rpc::plotms::Pid *reply );
        ::grpc::Status setShowGui( ::grpc::ServerContext *context,
                                   const ::rpc::plotms::Toggle *req,
                                   ::google::protobuf::Empty* );
        ::grpc::Status setGridSize( ::grpc::ServerContext *context,
                                    const ::rpc::plotms::GridSize *req,
                                    ::google::protobuf::Empty* );
        ::grpc::Status isDrawing( ::grpc::ServerContext *context,
                                  const ::google::protobuf::Empty*,
                                  ::rpc::plotms::Toggle *reply );
        ::grpc::Status clearPlots( ::grpc::ServerContext *context,
                                   const ::google::protobuf::Empty*,
                                   ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSFilename( ::grpc::ServerContext *context,
                                          const ::rpc::plotms::SetVis *req,
                                          ::google::protobuf::Empty* );
        ::grpc::Status setPlotAxes( ::grpc::ServerContext *context,
                                    const ::rpc::plotms::SetAxes *req,
                                    ::google::protobuf::Empty* );
        ::grpc::Status setShowAtm( ::grpc::ServerContext *context,
                                   const ::rpc::plotms::SetToggle *req,
                                   ::google::protobuf::Empty* );
        ::grpc::Status setShowTsky( ::grpc::ServerContext *context,
                                    const ::rpc::plotms::SetToggle *req,
                                    ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSSelection( ::grpc::ServerContext *context,
                                           const ::rpc::plotms::SetSelection *req,
                                           ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSAveraging( ::grpc::ServerContext *context,
                                           const ::rpc::plotms::SetAveraging *req,
                                           ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSTransformations( ::grpc::ServerContext *context,
                                                 const ::rpc::plotms::SetTransform *req,
                                                 ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSCalibration( ::grpc::ServerContext *context,
                                             const ::rpc::plotms::SetCalibration *req,
                                             ::google::protobuf::Empty* );
        ::grpc::Status setFlagExtension( ::grpc::ServerContext *context,
                                         const ::rpc::plotms::SetFlagExtension *req,
                                         ::google::protobuf::Empty* );
        ::grpc::Status setExportRange( ::grpc::ServerContext *context,
                                       const ::rpc::plotms::ExportRange *req,
                                       ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSIterate( ::grpc::ServerContext *context,
                                         const ::rpc::plotms::SetIterate *req,
                                         ::google::protobuf::Empty* );
        ::grpc::Status setColorAxis( ::grpc::ServerContext *context,
                                     const ::rpc::plotms::SetString *req,
                                     ::google::protobuf::Empty* );
        ::grpc::Status setSymbol( ::grpc::ServerContext *context,
                                  const ::rpc::plotms::SetSymbol *req,
                                  ::google::protobuf::Empty* );
        ::grpc::Status setFlaggedSymbol( ::grpc::ServerContext *context,
                                         const ::rpc::plotms::SetSymbol *req,
                                         ::google::protobuf::Empty* );
        ::grpc::Status setConnect( ::grpc::ServerContext *context,
                                   const ::rpc::plotms::SetConnect *req,
                                   ::google::protobuf::Empty* );
        ::grpc::Status setLegend( ::grpc::ServerContext *context,
                                   const ::rpc::plotms::SetLegend *req,
                                   ::google::protobuf::Empty* );
        ::grpc::Status setTitle( ::grpc::ServerContext *context,
                                 const ::rpc::plotms::SetString *req,
                                 ::google::protobuf::Empty* );
        ::grpc::Status setTitleFont( ::grpc::ServerContext *context,
                                     const ::rpc::plotms::SetInt *req,
                                     ::google::protobuf::Empty* );
        ::grpc::Status setXAxisFont( ::grpc::ServerContext *context,
                                     const ::rpc::plotms::SetInt *req,
                                     ::google::protobuf::Empty* );
        ::grpc::Status setYAxisFont( ::grpc::ServerContext *context,
                                     const ::rpc::plotms::SetInt *req,
                                     ::google::protobuf::Empty* );
        ::grpc::Status setXAxisLabel( ::grpc::ServerContext *context,
                                      const ::rpc::plotms::SetString *req,
                                      ::google::protobuf::Empty* );
        ::grpc::Status setYAxisLabel( ::grpc::ServerContext *context,
                                      const ::rpc::plotms::SetString *req,
                                      ::google::protobuf::Empty* );
        ::grpc::Status setGridParams( ::grpc::ServerContext *context,
                                      const ::rpc::plotms::SetGrid *req,
                                      ::google::protobuf::Empty* );
        ::grpc::Status setXRange( ::grpc::ServerContext *context,
                                  const ::rpc::plotms::SetRange *req,
                                  ::google::protobuf::Empty* );
        ::grpc::Status setYRange( ::grpc::ServerContext *context,
                                  const ::rpc::plotms::SetRange *req,
                                  ::google::protobuf::Empty* );
        ::grpc::Status setPlotMSPageHeaderItems( ::grpc::ServerContext *context,
                                                 const ::rpc::plotms::SetString *req,
                                                 ::google::protobuf::Empty* );
        ::grpc::Status save( ::grpc::ServerContext *context,
                             const ::rpc::plotms::Save *req,
                             ::google::protobuf::Empty* );
        ::grpc::Status update( ::grpc::ServerContext *context,
                               const ::google::protobuf::Empty *req,
                               ::google::protobuf::Empty* );

    protected:
        // update parameters
        void update_parameters( );
        void update_parameters(int index);
        void populate_selection( const ::rpc::plotms::SetSelection &req, PlotMSSelection &sel );
        enum group_tags { T_MSDATA, T_CACHE, T_AXES, T_ITER, T_DISP, T_CAN, T_HEAD };
        std::map<std::pair<int,int>,PlotMSPlotParameters::Group*> param_groups;
        PMS_PP_MSData *sys_ppdata(PlotMSPlotParameters* sp);
        PMS_PP_Cache *sys_ppcache(PlotMSPlotParameters* sp);
        PMS_PP_Axes *sys_ppaxes(PlotMSPlotParameters* sp);
        PMS_PP_Iteration *sys_ppiter(PlotMSPlotParameters* sp);
        PMS_PP_Display *sys_ppdisp(PlotMSPlotParameters* sp);
        PMS_PP_Canvas *sys_ppcan(PlotMSPlotParameters* sp);
        PMS_PP_PageHeader *sys_pphead(PlotMSPlotParameters* sp);
        PMS_PP_MSData *ppdata(int index);
        PMS_PP_Cache *ppcache(int index);
        PMS_PP_Axes *ppaxes(int index);
        PMS_PP_Iteration *ppiter(int index);
        PMS_PP_Display *ppdisp(int index);
        PMS_PP_Canvas *ppcan(int index);
        PMS_PP_PageHeader *pphead(int index);
        void qtGO( std::function<void()> );

    signals:
        void new_op( );

    private:
        bool invalid_index( int );
        PlotEngine *itsPlotms_;
        PlotMSPlotter *plotter_;
    };

    class grpcPlotMSState {
    public:
        grpcPlotMSState(PlotEngine *v) :
            plotms_service(new grpcPlotMS(v)),
            shutdown_service(new grpcPMSShutdown(v)) { }

        std::string uri;

        std::unique_ptr<grpc::Server> server;

        std::unique_ptr<grpcPlotMS> plotms_service;
        std::unique_ptr<grpcPMSShutdown> shutdown_service;

        ~grpcPlotMSState( ) {
            if (getenv("GRPC_DEBUG")) {
                fprintf(stdout, "stopping grpc server...\n");
                fflush(stdout);
            }
            if ( server ) server->Shutdown( );
        }
    };
}
#endif
