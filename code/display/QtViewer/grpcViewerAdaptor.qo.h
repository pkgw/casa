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

#ifndef GRPCVIEWERADAPTOR_H_
#define GRPCVIEWERADAPTOR_H_

#include <map>
#include <list>
#include <memory>
#include <functional>
#include <display/QtViewer/QtDisplayPanelGui.qo.h>
#include <casagrpc/protos/img.grpc.pb.h>
#include <casagrpc/protos/shutdown.grpc.pb.h>
#include <grpc++/grpc++.h>

namespace casa { //# NAMESPACE CASA - BEGIN

    class grpcShutdown : public QObject, public ::casatools::rpc::Shutdown::Service {

        Q_OBJECT    //# Allows slot/signal definition.  Must only occur in
                    //# implement/.../*.h files
    public:
        grpcShutdown( QtViewer *qtv );
        ::grpc::Status now(::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty*);

    signals:
        void exit_now( );

    private:
        QtViewer *viewer_;
        
    };

    class grpcImageViewer : public QObject, public ::rpc::img::view::Service {

        Q_OBJECT    //# Allows slot/signal definition.  Must only occur in
                    //# implement/.../*.h files
    protected:
        void qtGO( std::function<void()> );

    private:

        class data_desc {
        public:
            data_desc( int idx, const std::string &pathx, const std::string &typex,
                       QtDisplayData *ddx, QtDisplayPanelGui *dpx ) :
                id_(idx), path_(pathx), type_(typex), dd_(ddx), dp_(dpx) { }

            data_desc( int idx ) : id_(idx), dd_(0) { }
            data_desc( ) : id_(0), dd_(0) { }

            int &id( ) {
                return id_;
            }
            int id( ) const {
                return id_;
            }
            std::string &path( ) {
                return path_;
            }
            const std::string &path( ) const {
                return path_;
            }
            std::string &type( ) {
                return type_;
            }
            const std::string &type( ) const {
                return type_;
            }
            QtDisplayData *&data( ) {
                return dd_;
            }
            const QtDisplayData *data( ) const {
                return dd_;
            }
            QtDisplayPanelGui *&panel( ) {
                return dp_;
            }
            const QtDisplayPanelGui *panel( ) const {
                return dp_;
            }

        private:
            int id_;
            std::string path_;
            std::string type_;
            QtDisplayData *dd_;
            QtDisplayPanelGui *dp_;

            // QtDisplayData does not have a copy constructor...
            // wonder if we'll need to copy our descriptor...
            data_desc( const data_desc &other);
            data_desc &operator=( const data_desc &);
        };


        class panel_desc {
        public:

            panel_desc( ) : panel_(0) { }

            panel_desc(QtDisplayPanelGui*p) : panel_(p) { }

            std::list<int> &data( ) {
                return data_;
            }
            const std::list<int> &data( ) const {
                return data_;
            }
            QtDisplayPanelGui *&panel( ) {
                return panel_;
            }
            const QtDisplayPanelGui *panel( ) const {
                return panel_;
            }

        private:
            std::list<int> data_;
            QtDisplayPanelGui *panel_;
        };

    public:

        // Constructor which takes the application.
        grpcImageViewer( QtViewer * );

        ::grpc::Status panel( ::grpc::ServerContext *context,
                              const ::rpc::img::NewPanel *req,
                              ::rpc::img::Id *reply );
        ::grpc::Status hide( ::grpc::ServerContext *context,
                             const ::rpc::img::Id *req,
                             ::google::protobuf::Empty* );
        ::grpc::Status show( ::grpc::ServerContext *context,
                             const ::rpc::img::Id *req,
                             ::google::protobuf::Empty* );
        ::grpc::Status freeze( ::grpc::ServerContext *context,
                               const ::rpc::img::Id *req,
                               ::google::protobuf::Empty* );
        ::grpc::Status unfreeze( ::grpc::ServerContext *context,
                                 const ::rpc::img::Id *req,
                                 ::google::protobuf::Empty* );
        ::grpc::Status close( ::grpc::ServerContext *context,
                              const ::rpc::img::Id *req,
                              ::google::protobuf::Empty* );
        ::grpc::Status release( ::grpc::ServerContext *context,
                                const ::rpc::img::Id *req,
                                ::google::protobuf::Empty* );

        ::grpc::Status axes( ::grpc::ServerContext *context,
                             const ::rpc::img::Axes *req,
                             ::google::protobuf::Empty* );

        ::grpc::Status popup( ::grpc::ServerContext *context,
                              const ::rpc::img::PopUp *req,
                              ::google::protobuf::Empty* );

        ::grpc::Status load( ::grpc::ServerContext *context,
                             const ::rpc::img::NewData *req,
                             ::rpc::img::Id *reply );
        ::grpc::Status reload( ::grpc::ServerContext *context,
                             const ::rpc::img::Id *req,
                             ::google::protobuf::Empty* );
        ::grpc::Status unload( ::grpc::ServerContext *context,
                               const ::rpc::img::Id *req,
                               ::google::protobuf::Empty* );

        ::grpc::Status restore( ::grpc::ServerContext *context,
                                const ::rpc::img::Restore *req,
                                ::rpc::img::Id *reply );

        ::grpc::Status colormap( ::grpc::ServerContext *context,
                                  const ::rpc::img::ColorMap *req,
                                  ::google::protobuf::Empty* );

        ::grpc::Status colorwedge( ::grpc::ServerContext *context,
                                   const ::rpc::img::Toggle *req,
                                   ::google::protobuf::Empty* );

        ::grpc::Status datarange( ::grpc::ServerContext *context,
                                  const ::rpc::img::DataRange *req,
                                  ::google::protobuf::Empty* );

        ::grpc::Status contourlevels( ::grpc::ServerContext *context,
                                      const ::rpc::img::ContourLevels *req,
                                      ::google::protobuf::Empty* );
        ::grpc::Status contourthickness( ::grpc::ServerContext *context,
                                         const ::rpc::img::ContourThickness *req,
                                         ::google::protobuf::Empty* );
        ::grpc::Status contourcolor( ::grpc::ServerContext *context,
                                     const ::rpc::img::ContourColor *req,
                                     ::google::protobuf::Empty* );

        ::grpc::Status channel( ::grpc::ServerContext *context,
                                const ::rpc::img::SetChannel *req,
                                ::google::protobuf::Empty* );
        ::grpc::Status zoomlevel( ::grpc::ServerContext *context,
                                  const ::rpc::img::SetZoomLevel *req,
                                  ::google::protobuf::Empty* );
        ::grpc::Status zoombox( ::grpc::ServerContext *context,
                                const ::rpc::img::SetZoomBox *req,
                                ::google::protobuf::Empty* );

        ::grpc::Status output( ::grpc::ServerContext *context,
                               const ::rpc::img::Output *req,
                               ::google::protobuf::Empty* );


        ::grpc::Status fileinfo( ::grpc::ServerContext *context,
                                 const ::rpc::img::Path *req,
                                 ::rpc::img::FileInfo *reply );
        ::grpc::Status keyinfo( ::grpc::ServerContext *context,
                                 const ::rpc::img::Id *req,
                                 ::rpc::img::KeyInfo *reply );
        // if path is null string, return current path
        // if path is not null (and represents a valid dir),
        //            cd to new path && return new path
        ::grpc::Status cwd( ::grpc::ServerContext *context,
                            const ::rpc::img::Path *req,
                            ::rpc::img::Path *reply );
        // duplicates the shutdown::now rpc call...
        // duplicates grpcShutdown::now(...)
        ::grpc::Status done( ::grpc::ServerContext*,
                             const ::google::protobuf::Empty*,
                             ::google::protobuf::Empty* );

        int get_id( QtDisplayPanelGui* );
        int get_id( QtDisplayPanelGui *, QtDisplayData *, const std::string &path, const std::string &type );

        QtDisplayData *finddata( int key );
        QtDisplayPanelGui *findpanel( int key /*, bool create=true*/ );
        void erase_panel( int panel );
        void erase_data( int );

        bool load_data( QtDisplayPanelGui*, int );
        bool unload_data( QtDisplayPanelGui*, int, bool = true );

        bool printraster( QtDisplayPanel *panel, const QString &type,
                          const QString &file, double scale );
        bool printps( QtDisplayPanel *panel, const QString &type, const QString &file,
                      int dpi, const QString &orientation, const QString &media );
        void  adjusteps( const char *from, const char *to,
                         const QSize &wcmax, const QRect &viewport );

    signals:

        void exit_now( );
        void new_op( );

    public slots:

        void handle_destroyed_panel(QObject*);

    private:
        typedef std::map<int,panel_desc*> panelmap;
        std::recursive_mutex managed_panels_mutex;
        panelmap managed_panels;

        typedef std::map<int,data_desc*> datamap;
        std::recursive_mutex managed_datas_mutex;
        datamap managed_datas;

        QtViewer *viewer_;
        
    };

    class grpcViewerState {
    public:
        grpcViewerState(QtViewer *v) : 
            viewer_service(new grpcImageViewer(v)),
            shutdown_service(new grpcShutdown(v)) { }

        std::string uri;

        std::unique_ptr<grpc::Server> server;

        std::unique_ptr<grpcImageViewer> viewer_service;
        std::unique_ptr<grpcShutdown> shutdown_service;

        ~grpcViewerState( ) {
            if (getenv("GRPC_DEBUG")) {
                fprintf(stdout, "stopping grpc server...\n");
                fflush(stdout);
            }
            if ( server ) server->Shutdown( );
        }
    };

}

#endif
