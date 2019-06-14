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
#include <climits>
#include <algorithm>
#include <display/QtViewer/QtViewer.qo.h>
#include <display/QtViewer/QtCleanPanelGui.qo.h>
#include <display/QtViewer/QtCleanPanelGui2.qo.h>
#include <display/QtViewer/QtDisplayData.qo.h>
#include <display/QtViewer/grpcViewerAdaptor.qo.h>
#include <display/Display/DisplayState.h>

namespace casa { //# NAMESPACE CASA - BEGIN

    bool ends_with( const std::string& str, const std::string& ending ) {
        return ( str.size( ) >= ending.size( ) ) && equal( ending.rbegin( ), ending.rend( ), str.rbegin( ) );
    };

    /**********************************************************************************************************
     **********************************************************************************************************
     *****  grpcShutdown                                                                                  *****
     **********************************************************************************************************
     **********************************************************************************************************/
    grpcShutdown::grpcShutdown( QtViewer *viewer ) : viewer_(viewer) {  }

    ::grpc::Status grpcShutdown::now( ::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty* ) {
        if (getenv("GRPC_DEBUG")) {
            std::cout << "received shutdown notification..." << std::endl;
            fflush(stdout);
        }
        static auto bye_bye = std::async( std::launch::async, [&]( ) { sleep(2); emit exit_now( ); } );
        return grpc::Status::OK;
    }

    /**********************************************************************************************************
     **********************************************************************************************************
     *****  grpcImageViewer                                                                               *****
     **********************************************************************************************************
     **********************************************************************************************************/
    int grpcImageViewer::get_id( QtDisplayPanelGui *panel ) {
        std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);

        for ( panelmap::iterator iter = managed_panels.begin(); iter != managed_panels.end(); ++iter ) {
            if ( iter->second->panel() == panel )
                return iter->first;
        }

        int index = QtId::get_id( );
        managed_panels.insert(panelmap::value_type(index, new panel_desc(panel)));
        return index;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    int grpcImageViewer::get_id( QtDisplayPanelGui *panel, QtDisplayData *data,
                                 const std::string &path, const std::string &type ) {

        int index = 0;
        {
            std::lock_guard<std::recursive_mutex> p_guard(managed_datas_mutex);

            for ( datamap::iterator iter = managed_datas.begin(); iter != managed_datas.end(); ++iter ) {
                if ( iter->second->data() == data )
                    return iter->second->id();
            }

            index = QtId::get_id( );
            data_desc *dd = new data_desc(index, path, type, data, panel );
            managed_datas.insert(datamap::value_type(index, dd));
        }

        if ( index != 0 ) {
            std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);

            for ( panelmap::iterator dpiter = managed_panels.begin(); dpiter != managed_panels.end(); ++dpiter ) {
                if ( dpiter->second->panel() == panel ) {
                    dpiter->second->data().push_back(index);
                    break;
                }
            }
        }

        return index;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    QtDisplayData *grpcImageViewer::finddata( int key ) {

        datamap::iterator iter = managed_datas.find( key );
        if ( iter != managed_datas.end( ) )
            return iter->second->data( );

        return 0;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    QtDisplayPanelGui *grpcImageViewer::findpanel( int key /*, bool create*/ ) {
        std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);

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

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcImageViewer::erase_panel( int panel ) {
        std::lock_guard<std::recursive_mutex> p_guard(managed_panels_mutex);

        if ( panel == 0 ) panel = INT_MAX;

        QtDisplayPanelGui *win = NULL;
        panelmap::iterator dpiter = managed_panels.find( panel );
        if ( dpiter != managed_panels.end( ) ) {
            std::list<int> &data = dpiter->second->data();
            for ( std::list<int>::iterator diter = data.begin( ); diter != data.end(); ++diter ) {
                erase_data(*diter);
            }
            win = dpiter->second->panel();
            delete dpiter->second;
            managed_panels.erase(dpiter);
        }

    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcImageViewer::erase_data( int index ) {
        std::lock_guard<std::recursive_mutex> p_guard(managed_datas_mutex);

        datamap::iterator dditer = managed_datas.find(index);
        if ( dditer != managed_datas.end( ) ) {
            delete dditer->second;
            dditer->second = 0;
            managed_datas.erase(dditer);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
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


    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    grpcImageViewer::grpcImageViewer(QtViewer  *viewer) : viewer_(viewer) { }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcImageViewer::qtGO( std::function<void()> func ) {
        std::lock_guard<std::mutex> exc(viewer_->grpc_queue_mutex);
        viewer_->grpc_queue.push(func);
        emit new_op( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::panel( ::grpc::ServerContext *context,
                                           const ::rpc::img::NewPanel *req,
                                           ::rpc::img::Id *reply) {
        static const auto debug = getenv("GRPC_DEBUG");

        if (debug) {
            std::cout << "received grpc panel( " << req->type( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        std::string type = req->type( );
        bool hidden = req->hidden( );
        QtDisplayPanelGui *result = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {
                  if (debug) {
                      std::cout << "creating qt panel ( " << req->type( ) << " )... (thread " <<
                          std::this_thread::get_id() << ")" << std::endl;
                      fflush(stdout);
                  }

                  if ( type == "clean" ) {

                      // <drs> need to ensure that this is not leaked...
                      result = viewer_->createInteractiveCleanGui( );
                      if ( hidden ) result->hide( );

//*grpc*              connect(result, SIGNAL(interact(QVariant)), this, SLOT(handle_interact(QVariant)));

                  } else if ( type == "clean2" ) {

                      // <drs> need to ensure that this is not leaked...
                      result = viewer_->createInteractiveCleanGui2( );

                      if ( hidden ) result->hide( );

//*grpc*              connect(cpg_, SIGNAL(interact(QVariant)), this, SLOT(handle_interact(QVariant)));

                  } else {

                      result = viewer_->createDPG( );
//*grpc*              connect( result, SIGNAL(destroyed(QObject*)), SLOT(handle_destroyed_panel(QObject*)) );

                      if ( ends_with(type, ".rstr") ) {
                          struct stat buf;
                          if ( stat( type.c_str( ), &buf ) == 0 ) {
                              result->restorePanelState(type);
                          }
                      }

                      if ( hidden ) result->hide( );
                  }

                  prom.set_value(get_id(result));
            } );


        if (debug) {
            std::cout << "waiting for grpc panel( " << req->type( ) << " ) result... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        auto fut = prom.get_future( );
        fut.wait( );
        auto ret = fut.get( );
        if (debug) {
            std::cout << "returning grpc panel( " << req->type( ) << " ) result: " <<
                ret << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }
        reply->set_id(ret);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
#define PANEL_OP(PANELFUNC,NAME,ACTIVE_DESC,PAST_DESC,ACTION)                                          \
    ::grpc::Status grpcImageViewer::NAME( ::grpc::ServerContext *context,                              \
                                          const ::rpc::img::Id *req,                                   \
                                          ::google::protobuf::Empty* ) {                               \
                                                                                                       \
        static const auto debug = getenv("GRPC_DEBUG");                                                \
                                                                                                       \
        if (debug) {                                                                                   \
            std::cout << "received grpc " #NAME "( " << req->id( ) << " ) event... (thread " <<        \
                std::this_thread::get_id() << ")" << std::endl;                                        \
            fflush(stdout);                                                                            \
        }                                                                                              \
                                                                                                       \
        int id = req->id( );                                                                           \
        auto panel = findpanel(id);                                                                    \
                                                                                                       \
        if ( panel == 0 ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such panel id");        \
                                                                                                       \
        ACTION                                                                                         \
                                                                                                       \
        std::promise<bool> prom;                                                                       \
                                                                                                       \
        qtGO( [&]( ) {                                                                                 \
                  if (debug) {                                                                         \
                      std::cout << #ACTIVE_DESC " qt panel ( " << id << " )... (thread " <<            \
                          std::this_thread::get_id() << ")" << std::endl;                              \
                      fflush(stdout);                                                                  \
                  }                                                                                    \
                  panel->PANELFUNC( );                                                                 \
                  prom.set_value(true);                                                                \
            } );                                                                                       \
                                                                                                       \
        if (debug) {                                                                                   \
            std::cout << "waiting for panel to be " #PAST_DESC " ( " << id << " )... (thread " <<      \
                std::this_thread::get_id() << ")" << std::endl;                                        \
            fflush(stdout);                                                                            \
        }                                                                                              \
        auto fut = prom.get_future( );                                                                 \
        fut.wait( );                                                                                   \
        auto ret = fut.get( );                                                                         \
        if (debug) {                                                                                   \
            std::cout << "completed " #ACTIVE_DESC " ( " << ret << " )" <<                             \
                " (thread " << std::this_thread::get_id() << ")" << std::endl;                         \
            fflush(stdout);                                                                            \
        }                                                                                              \
        return grpc::Status::OK;                                                                       \
    }

    PANEL_OP(show,show,showing,shown,;)
    PANEL_OP(hide,hide,hiding,hid,;)
    PANEL_OP(hold,freeze,freezing,frozen,;)
    PANEL_OP(release,unfreeze,unfreezing,unfrozen,;)
    PANEL_OP(closeMainPanel,close,closing,closed,erase_panel( id );)
    PANEL_OP(releaseMainPanel,release,releasing,released,erase_panel( id );)

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::axes( ::grpc::ServerContext *context,
                                          const ::rpc::img::Axes *req,
                                          ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto panel = req->panel( ).id( );
        auto x = req->x( );
        auto y = req->y( );
        auto z = req->z( );

        if (debug) {
            std::cout << "received grpc axes( " << panel << ", " << x << ", " << y << ", " << z <<
                req->panel( ).id( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        casacore::Record rec;
        bool update_axes = false;

        // -- should not be necessary, but when "setOptions( )" is called  --
        // -- from a script the axes values being changed are not updated  --
        // -- in the data options panel...                                 --
        casacore::Record optionsChangedRec;

        if ( x.size( ) > 0 ) {
            rec.define( "xaxis", x );
            update_axes = true;
            casacore::Record xaxis;
            xaxis.define( "value", x );
            optionsChangedRec.defineRecord("xaxis",xaxis);
        }

        if ( y.size( ) > 0 ) {
            rec.define( "yaxis", y );
            update_axes = true;
            casacore::Record yaxis;
            yaxis.define( "value", y );
            optionsChangedRec.defineRecord("yaxis",yaxis);
        }

        if ( z.size( ) > 0 ) {
            rec.define( "zaxis", z );
            update_axes = true;
            casacore::Record zaxis;
            zaxis.define( "value", z );
            optionsChangedRec.defineRecord("zaxis",zaxis);
        }

        if ( update_axes == false ) {
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes specified");
        }

        panelmap::iterator dpiter = managed_panels.find( panel == 0 ? INT_MAX : panel );

        if ( dpiter == managed_panels.end( ) ) {
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot find panel id");
        }

        bool set_option = false;
        std::list<int> &data = dpiter->second->data( );
        for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
            datamap::iterator dditer = managed_datas.find(*diter);
            if ( dditer != managed_datas.end( ) ) {
                QtDisplayData *dd = dditer->second->data( );
                dd->setOptions(rec,true);
                set_option = true;
                qtGO( [&]( ) {
                          // -- it seems like a better idea to signal change here (when they    --
                          // -- fail to update due to scripting ops ) instead of littering the  --
                          // -- code with unnecessary updates for the GUI case                  --
                          dd->emitOptionsChanged( optionsChangedRec );
                      } );
            }
        }

        if ( set_option ) return grpc::Status::OK;
        else return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes made");

    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::popup( ::grpc::ServerContext*,
                                           const ::rpc::img::PopUp *req,
                                           ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        if (debug) {
            std::cout << "received grpc popup( " << req->name( ) << ", " <<
                req->panel( ).id( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        int id  = req->panel( ).id( );
        auto name = req->name( );
        auto panel = findpanel(id);

        if ( panel == 0 ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such panel id");

        std::promise<bool> prom;

        if ( name == "open" ) {

            qtGO( [&]( ) {
                  if (debug) {
                      std::cout << " opening qt data options for panel ( " << id << " )... (thread " <<
                          std::this_thread::get_id() << ")" << std::endl;
                      fflush(stdout);
                  }
                  panel->showDataManager( );
                  prom.set_value(true);
                } );

        } else return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "unknown popup-widow name, expected 'open'");

        // must wait otherwise our stack state (upon which the lambda depends) will disappear...
        auto fut = prom.get_future( );
        fut.wait( );
        auto ret = fut.get( );

        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::load( ::grpc::ServerContext *context,
                                          const ::rpc::img::NewData *req,
                                          ::rpc::img::Id *reply ) {

        static const auto debug = getenv("GRPC_DEBUG");

        if (debug) {
            std::cout << "received grpc load( " << req->path( ) << ", " <<
                req->type( ) << ", " << req->type( ) << ", " << req->scale( ) <<
                req->panel( ).id( ) << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        auto id = req->panel( ).id( );
        auto panel = findpanel(id);
        auto path = req->path( );
        auto displaytype = req->type( );
        auto scaling = req->scale( );

        if ( panel == 0 ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such panel id");

        struct stat buf;
        if ( stat( path.c_str( ), &buf ) < 0 ) {
            // file (or dir) does not exist
            return grpc::Status(grpc::StatusCode::NOT_FOUND, "image path does not exist");
        }

        std::string datatype = viewer_->filetype(path);

        if ( datatype != "image" && datatype != "ms" )
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "invalid/unknown/unloadable image file type");
        if ( displaytype != "raster" && displaytype != "contour" &&
             displaytype != "vector" && displaytype == "marker" )
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "invalid/unknown display raster type");

        std::promise<int> prom;

        qtGO( [&]( ) {
                if (debug) {
                    std::cout << " opening qt data ( " << path << ", " << id << " )... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stdout);
                }

                viewer::DisplayDataOptions ddo;
                auto data =  panel->createDD( path, datatype, displaytype, true,
                                              -1, false, false, false, ddo );

                if ( scaling != 0.0 ) data->setRasterPowerScaling(scaling);

                panel->addedData( displaytype.c_str( ), data );

                prom.set_value( get_id( panel, data, path, displaytype ) );
            } );

        // must wait otherwise our stack state (upon which the lambda depends) will disappear...
        auto fut = prom.get_future( );
        fut.wait( );
        auto ret = fut.get( );

        reply->set_id(ret);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    bool grpcImageViewer::load_data( QtDisplayPanelGui *panel, int index ) {

        static const auto debug = getenv("GRPC_DEBUG");

        datamap::iterator iter = managed_datas.find( index );
        if ( iter == managed_datas.end( ) ) return false;

        struct stat buf;
        const std::string &path = iter->second->path( );
        const std::string &displaytype = iter->second->type( );
        if ( stat(path.c_str( ),&buf) < 0 ) return false;

        std::string datatype = viewer_->filetype(path);
        if ( datatype == "image" ) {
            if ( displaytype == "raster" || displaytype == "contour" ||
                 displaytype == "vector" || displaytype == "marker" ) {

                std::promise<bool> prom;

                qtGO( [&]( ) {
                        if (debug) {
                            std::cout << " qt data load ( " << path << ", " << displaytype << " )... (thread " <<
                                std::this_thread::get_id() << ")" << std::endl;
                            fflush(stdout);
                        }

                        QtDisplayData *dp = 0;
                        if ( panel != 0 ) {
                            if ( iter->second->data( ) != 0 ) {
                                panel->removeDD( iter->second->data( ) );
                                iter->second->data( ) = 0;
                            }

                            panel->autoDDOptionsShow = false;
                            dp = panel->createDD( path, datatype, displaytype, false);
                            panel->displayPanel()->registerDD(dp);
                            panel->autoDDOptionsShow = true;
                            panel->addedData( displaytype.c_str( ), dp );
                            iter->second->data( ) = dp;
                            prom.set_value(true);
                        } else prom.set_value(false);
                    } );

                auto fut = prom.get_future( );
                fut.wait( );
                return fut.get( );
            } else return false;
        } else return false;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    bool grpcImageViewer::unload_data( QtDisplayPanelGui */*panel*/, int index, bool erase ) {
        datamap::iterator iter = managed_datas.find( index );
        if ( iter == managed_datas.end( ) ) return false;
        if ( iter->second->data( ) != 0 ) {
            std::promise<bool> prom;
            qtGO( [&]( ) {
                    iter->second->panel()->removeDD( iter->second->data( ) );
//***               fails to notify the wrench that things have changed...
//                  panel->unregisterDD( iter->second->data( ) );
                    iter->second->data( ) = 0;
                    if ( erase ) managed_datas.erase(iter);
                    prom.set_value(true);
                } );

            auto fut = prom.get_future( );
            fut.wait( );
            return fut.get( );
        } else return false;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::reload( ::grpc::ServerContext *context,
                                            const ::rpc::img::Id *req,
                                            ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto id = req->id( );

        if (debug) {
            std::cout << "received grpc reload( " << id << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        std::promise<bool> prom;

        managed_panels_mutex.lock( );
        panelmap::iterator dpiter = managed_panels.find( id );
        if ( dpiter != managed_panels.end( ) ) {

            qtGO( [&]( ) {
                    if (debug) {
                        std::cout << "qt reload of data ( " << id << " )... (thread " <<
                            std::this_thread::get_id() << ")" << std::endl;
                        fflush(stdout);
                    }
                    QtDisplayPanel::panel_state state = dpiter->second->panel( )->displayPanel( )->getPanelState( );
                    dpiter->second->panel( )->displayPanel()->hold( );

                    std::list<int> &data = dpiter->second->data( );
                    for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                        unload_data( dpiter->second->panel( ), *diter, false );
                    }

                    for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                        load_data( dpiter->second->panel( ), *diter );
                    }
                    dpiter->second->panel( )->displayPanel()->setPanelState( state );
                    dpiter->second->panel( )->displayPanel()->release( );
                    prom.set_value(true);
                } );

            managed_panels_mutex.unlock( );

        } else {

            managed_panels_mutex.unlock( );
            std::lock_guard<std::recursive_mutex> d_guard(managed_datas_mutex);

            datamap::iterator dmiter = managed_datas.find( id );
            if ( dmiter != managed_datas.end( ) ) {
                load_data( dmiter->second->panel(), dmiter->first );
                prom.set_value(true);
            } else {
                prom.set_value(false);
            }
        }

        // must wait otherwise our stack state (upon which the lambda depends) will disappear...
        auto fut = prom.get_future( );
        fut.wait( );
        auto ret = fut.get( );

        return ret ? grpc::Status::OK : grpc::Status(grpc::StatusCode::NOT_FOUND, "no such panel id");
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::unload( ::grpc::ServerContext *context,
                                            const ::rpc::img::Id *req,
                                            ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto id = req->id( );

        if (debug) {
            std::cout << "received grpc unload( " << id << " ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        std::lock_guard<std::recursive_mutex> d_guard(managed_datas_mutex);

        datamap::iterator dmiter = managed_datas.find( id );
        if ( dmiter == managed_datas.end( ) )
            return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such data id");

        if ( unload_data( dmiter->second->panel( ), id ) )
            return grpc::Status::OK;
        else
            return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such data id");
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::restore( ::grpc::ServerContext *context,
                                             const ::rpc::img::Restore *req,
                                             ::rpc::img::Id *reply ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto id = req->panel( ).id( );
        auto path = req->path( );

        if (debug) {
            std::cout << "received grpc restore( " << path << ", " << id <<
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        struct stat buf;
        if ( stat(path.c_str( ),&buf) < 0 || ! S_ISREG(buf.st_mode) ) {
            // file (or dir) does not exist
            return grpc::Status(grpc::StatusCode::NOT_FOUND, "restore file not found");
        }

        auto dpg = findpanel( id );
        if ( ! dpg ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "panel id not found");

        std::promise<bool> prom;
        qtGO( [&]( ) { prom.set_value(dpg->displayPanel()->restorePanelState(path)); } );

        auto fut = prom.get_future( );
        fut.wait( );
        if ( fut.get( ) ) {
            reply->set_id(get_id( dpg ));
            return grpc::Status::OK;
        } else {
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "unable to load restore file");
        }

    }


    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::colormap( ::grpc::ServerContext *context,
                                              const ::rpc::img::ColorMap *req,
                                              ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto panel_or_data = req->id( ).id( );
        auto map = req->map( );

        if (debug) {
            std::cout << "received grpc colormap( " << panel_or_data << ", " << map <<
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        QtDisplayData *dd = finddata(panel_or_data);

        if ( dd == 0 ) {

            panelmap::iterator dpiter = managed_panels.find( panel_or_data == 0 ? INT_MAX : panel_or_data );

            if ( dpiter == managed_panels.end( ) ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot find panel/data id");
            }

            std::list<int> &data = dpiter->second->data( );
            // first verify that the colormap name is valid for all qt display datas
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->hasColormap( ) && dd->isValidColormap( QString(map.c_str( )) ) == false ) {
                        return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "invalid colormap for one (or more) display data(s)");
                    }
                }
            }
            // next replace the colormap
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->hasColormap( ) ) {
                        qtGO( [&]( ) { dd->setColormap( map ); } );
                    }
                }
            }
            return grpc::Status::OK;
        }

        if ( dd->hasColormap( ) ) {
            if ( dd->isValidColormap( QString(map.c_str()) ) == false ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "invalid colormap");
            }
            qtGO( [&]( ) { dd->setColormap(map); } );
        }

        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::colorwedge( ::grpc::ServerContext *context,
                                                const ::rpc::img::Toggle *req,
                                                ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto panel_or_data = req->id( ).id( );
        auto show = req->state( );

        if (debug) {
            std::cout << "received grpc colorwedge( " << panel_or_data << ", " << show <<
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        QtDisplayData *dd = finddata(panel_or_data);

        casacore::Record cw;
        cw.define( "wedge", show ? "Yes" : "No" );

        if ( dd == 0 ) {

            panelmap::iterator dpiter = managed_panels.find( panel_or_data == 0 ? INT_MAX : panel_or_data );

            if ( dpiter == managed_panels.end( ) ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot find panel/data id");
            }

            std::list<int> &data = dpiter->second->data( );
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->hasColorBar( ) ) {
                        qtGO( [=]( ) { dd->setOptions( cw, true); } );
                    }
                }
            }
            return grpc::Status::OK;
        }

        if ( dd->hasColorBar( ) ) {
            qtGO( [=]( ) { dd->setOptions( cw, true); } );
        }

        return grpc::Status::OK;

    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::datarange( ::grpc::ServerContext *context,
                                               const ::rpc::img::DataRange *req,
                                               ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto data = req->data( ).id( );

        casacore::Vector<float> value;
        value.resize(2);
        value(0) = req->min( );
        value(1) = req->max( );

        if (debug) {
            std::cout << "received grpc datarange( " << value(0) << ", " << value(1) << ", " << data <<
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        casacore::Record rec;
        casacore::Record minmax;
        minmax.define( "value", value );
        rec.defineRecord( "minmaxhist", minmax );

        QtDisplayData *dd = finddata(data);

        if ( dd == 0 ) {

            if ( debug ) {
                std::cout << "searching for data id..." << std::endl;
                fflush(stdout);
            }

            // if we have a "id-less" panel (INT_MAX), see if we can
            // set the range on rasters there...
            if ( data == 0 ) {
                panelmap::iterator dpiter = managed_panels.find( INT_MAX );

                if ( dpiter != managed_panels.end( ) ) {
                    std::list<int> &data = dpiter->second->data( );
                    bool found = false;
                    for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                        datamap::iterator dditer = managed_datas.find(*diter);
                        if ( dditer != managed_datas.end( ) ) {
                            QtDisplayData *dd = dditer->second->data( );
                            if ( dd->displayType( ) == "raster" ) {
                                qtGO( [=]( ) { dd->setOptions(rec,true); } );
                                found = true;
                            }
                        }
                    }

                    if ( found ) return grpc::Status::OK;
                }
            }

            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "data id not found");
        }

        if ( debug ) {
            std::cout << "found data id (" << dd << ")..." << std::endl;
            fflush(stdout);
        }

        qtGO( [=]( ){ dd->setOptions(rec,true); } );

        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::contourlevels( ::grpc::ServerContext *context,
                                                   const ::rpc::img::ContourLevels *req,
                                                   ::google::protobuf::Empty* ) {


        static const auto debug = getenv("GRPC_DEBUG");

        auto panel_or_data = req->id( ).id( );

        if (debug) {
            std::cout << "received grpc contourlevels( " << panel_or_data << ", ..." << 
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        casacore::Record rec;
        bool update_contours = false;

        float baselevel = req->baselevel( );
        float unitlevel = req->unitlevel( );

        if ( req->levels( ).size( ) > 0 ) {
            casacore::Vector<float> value;
            value.resize( req->levels( ).size( ) );
            int i = 0;
            for ( auto iter = req->levels( ).begin( ); iter != req->levels( ).end( ); ++iter ) {
                value(i++) = *iter;
            }
            rec.define( "rellevels", value );
            update_contours = true;
        }

        if ( baselevel != 2147483648.0 ) {
            rec.define( "basecontour", baselevel );
            update_contours = true;
        }

        if ( unitlevel != 2147483648.0 ) {
            rec.define( "unitcontour", unitlevel );
            update_contours = true;
        }

        if ( update_contours == false ) {
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes specified");
        }

        QtDisplayData *dd = finddata(panel_or_data);

        if ( dd == 0 ) {

            panelmap::iterator dpiter = managed_panels.find( panel_or_data == 0 ? INT_MAX : panel_or_data );

            if ( dpiter == managed_panels.end( ) ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "could not find data/panel based on id");
            }

            bool set_contour = false;
            std::list<int> &data = dpiter->second->data( );
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->displayType( ) == "contour" ) {
                        qtGO( [=]( ) { dd->setOptions(rec,true); } );
                        set_contour = true;
                    }
                }
            }
            if ( set_contour ) return grpc::Status::OK;
            else return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes made");
        }

        qtGO( [=]( ) { dd->setOptions(rec,true); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::contourthickness( ::grpc::ServerContext *context,
                                                      const ::rpc::img::ContourThickness *req,
                                                      ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto panel_or_data = req->id( ).id( );
        auto thickness = req->thickness( );

        if (debug) {
            std::cout << "received grpc contourthickness( " << panel_or_data << ", " << thickness << 
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        if ( thickness <= 0 || thickness > 5 ) {
            return grpc::Status(grpc::StatusCode::OUT_OF_RANGE, "thickness should be greater than 0 and less than or equal to 5");
        }

        casacore::Record rec;
        rec.define( "line", thickness );

        QtDisplayData *dd = finddata(panel_or_data);

        if ( dd == 0 ) {

            panelmap::iterator dpiter = managed_panels.find( panel_or_data == 0 ? INT_MAX : panel_or_data );

            if ( dpiter == managed_panels.end( ) ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot find panel/data id");
            }

            bool set_thickness = false;
            std::list<int> &data = dpiter->second->data( );
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->displayType( ) == "contour" ) {
                        qtGO( [=]( ) { dd->setOptions(rec,true); } );
                        set_thickness = true;
                    }
                }
            }
            
            if ( set_thickness ) return grpc::Status::OK;
            else return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes made");

        }

        qtGO( [=]( ) { dd->setOptions(rec,true); } );
        return grpc::Status::OK;

    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::contourcolor( ::grpc::ServerContext *context,
                                                  const ::rpc::img::ContourColor *req,
                                                  ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto panel_or_data = req->id( ).id( );
        auto color = req->color( );

        if (debug) {
            std::cout << "received grpc contourthickness( " << panel_or_data << ", " << color << 
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        casacore::Record rec;
        std::list<std::string> valid { "foreground", "background", "black", "white", "red",
                                       "green", "blue", "cyan", "magenta", "yellow", "gray" };

        if ( std::find( valid.begin( ), valid.end(), color ) == valid.end( ) ) {
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "invalid color");
        }

        rec.define( "color", color );

        QtDisplayData *dd = finddata(panel_or_data);

        if ( dd == 0 ) {

            panelmap::iterator dpiter = managed_panels.find( panel_or_data == 0 ? INT_MAX : panel_or_data );

            if ( dpiter == managed_panels.end( ) ) {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot find panel/data id");
            }

            bool set_color = false;
            std::list<int> &data = dpiter->second->data( );
            for ( std::list<int>::iterator diter = data.begin(); diter != data.end(); ++diter ) {
                datamap::iterator dditer = managed_datas.find(*diter);
                if ( dditer != managed_datas.end( ) ) {
                    QtDisplayData *dd = dditer->second->data( );
                    if ( dd->displayType( ) == "contour" ) {
                        qtGO( [=]( ) { dd->setOptions(rec,true); } );
                        set_color = true;
                    }
                }
            }

            if ( set_color ) return grpc::Status::OK;
            else return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "no changes made");

        }

        qtGO( [=]( ) { dd->setOptions(rec,true); } );
        return grpc::Status::OK;

    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::channel( ::grpc::ServerContext *context,
                                             const ::rpc::img::SetChannel *req,
                                             ::google::protobuf::Empty* ) {

        auto id = req->panel( ).id( );
        auto chan = req->number( );

        QtDisplayPanelGui *dpg = findpanel( id );
        if ( ! dpg ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "panel id not found");
        if ( chan < 0 ) return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "channel number must be >= 0");

        std::promise<bool> prom;
        qtGO( [&]( ) {
                dpg->displayPanel()->goTo(chan);
                prom.set_value(true);
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::zoomlevel( ::grpc::ServerContext *context,
                                               const ::rpc::img::SetZoomLevel *req,
                                               ::google::protobuf::Empty* ) {

        auto id = req->panel( ).id( );
        auto level = req->level( );

        QtDisplayPanelGui *dpg = findpanel( id );

        if ( ! dpg ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "panel id not found");

        std::promise<bool> prom;
        qtGO( [&]( ) {
                if ( level == 0 )
                    dpg->displayPanel()->unzoom( );
                else if ( level < 0 )
                    dpg->displayPanel()->zoomOut( abs(level) );
                else
                    dpg->displayPanel()->zoomIn( level );
                prom.set_value(true);
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::zoombox( ::grpc::ServerContext *context,
                                             const ::rpc::img::SetZoomBox *req,
                                             ::google::protobuf::Empty* ) {

        auto id = req->panel( ).id( );
        QtDisplayPanelGui *dpg = findpanel( id );

        if ( ! dpg ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "panel id not found");

        auto coordinates = req->coord_type( );

        if ( coordinates != "pixel" && coordinates != "world" )
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "coordinates must be either 'world' or 'pixel'");

        std::vector<float> blc { req->blc( ).x( ), req->blc( ).y( ) };
        std::vector<float> trc { req->trc( ).x( ), req->trc( ).y( ) };

        casacore::Vector<double> v_blc((unsigned int)2);
        casacore::Vector<double> v_trc((unsigned int)2);
        for ( int i = 0; i < 2; ++i ) {
            v_blc[i] = blc[i];
            v_trc[i] = trc[i];
        }

        std::promise<bool> prom;
        if ( coordinates == "pixel" ) {

            qtGO( [&]( ) {
                    dpg->displayPanel()->zoom( v_blc, v_trc );
                    prom.set_value(true);
                } );

        } else if ( coordinates == "world" ) {

            casacore::Vector<double> vp_blc((unsigned int)2);
            casacore::Vector<double> vp_trc((unsigned int)2);
            bool OK = ( dpg->displayPanel( )->worldToLin( vp_blc, v_blc ) &&
                        dpg->displayPanel( )->worldToLin( vp_trc, v_trc ) );
            if ( OK == false )
                return grpc::Status( grpc::StatusCode::INVALID_ARGUMENT,
                                     "conversion from world to pixel coordinates failed" );

            qtGO( [&]( ) {
                    dpg->displayPanel()->zoom( vp_blc, vp_trc );
                    prom.set_value(true);
                } );

        }

        auto fut = prom.get_future( );
        fut.wait( );
        return grpc::Status::OK;
    }

    //  bool QtDBusViewerAdaptor::output( const QString &device, const QString &devicetype, int panel, double scale,
    //                                    int dpi, const QString &format, const QString &orientation, const QString &media ) {
    ::grpc::Status grpcImageViewer::output( ::grpc::ServerContext *context,
                                            const ::rpc::img::Output *req,
                                            ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        auto id = req->panel( ).id( );

        if (debug) {
            std::cout << "received grpc output( " << id <<
                " ) event... (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        auto dpg = findpanel( id );
        if ( ! dpg ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "panel id not found");

        auto device = QString(req->device( ).c_str( ));
        auto devicetype = QString(req->devicetype( ).c_str( ));
        auto format = QString(req->format( ).c_str( ));
        auto orientation = QString(req->orientation( ).c_str( ));
        auto media = QString(req->media( ).c_str( ));
        auto scale = req->scale( );
        auto dpi = req->dpi( );

        QString base;
        QString path;
        QString suffix;

        char printer_file[80];
        char printer_base[80];

        if ( devicetype == "printer" || devicetype == "ghostview" )
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "printer and ghostview devicetype no longer supported");
        else {
            QFileInfo file(device);
            base = file.completeBaseName( );
            path = file.absolutePath();
            suffix = file.suffix( ).toLower( );

            if ( suffix != "jpg" && suffix != "pdf" && suffix != "eps" && suffix != "ps" &&
                    suffix != "png" && suffix != "xbm" && suffix != "xpm" && suffix != "ppm" ) {
                suffix = format.toLower( );
                if ( suffix != "jpg" && suffix != "pdf" && suffix != "eps" && suffix != "ps" &&
                        suffix != "png" && suffix != "xbm" && suffix != "xpm" && suffix != "ppm" ) {
                    suffix = "jpg";
                }
            }
        }

        if ( suffix == "pdf" || suffix == "ps" || suffix == "eps" ) {
            printps( dpg->displayPanel(), suffix, path + "/" + base + "." + suffix, dpi, orientation, media );
        } else if ( suffix == "jpg" || suffix == "png" ||
                    suffix == "xbm" || suffix == "xpm" || suffix == "ppm" ) {
            printraster( dpg->displayPanel(), suffix, path + "/" + base + "." + suffix, scale );
        } else {
            return grpc::Status(grpc::StatusCode::NOT_FOUND, "unsupported output type");
        }

        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::fileinfo( ::grpc::ServerContext *context,
                                              const ::rpc::img::Path *req,
                                              ::rpc::img::FileInfo *reply ) {
        auto path = req->path( );

        struct stat buf;
        if ( stat(path.c_str( ),&buf) < 0 ) {
            // file (or dir) does not exist
            reply->set_type("nonexistent");
        } else {
            reply->set_type(viewer_->filetype(path));
        }
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::keyinfo( ::grpc::ServerContext *context,
                                             const ::rpc::img::Id *req,
                                             ::rpc::img::KeyInfo *reply ) {

        auto id = req->id( );

        if ( managed_datas.find( id ) != managed_datas.end( ) ) {
            reply->set_type("data");
        } else if ( managed_panels.find( id ) != managed_panels.end( ) ) {
            reply->set_type("panel");
        } else {
            reply->set_type("not found");
        }

        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::cwd( ::grpc::ServerContext *context,
                                         const ::rpc::img::Path *req,
                                         ::rpc::img::Path *reply ) {
        auto new_path = req->path( );

        char p[MAXPATHLEN+1];
        if ( new_path != "" ) {
            struct stat buf;
            sprintf( p, "%s", new_path.c_str( ) );
            if ( stat(p,&buf) == 0 && S_ISDIR(buf.st_mode) ) {
                if ( chdir(p) != 0 )
                    return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot change to directory");
            } else if ( stat(p,&buf) != 0 )
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "directory does not exist");
            else if ( ! S_ISDIR(buf.st_mode) )
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "path is not a directory");
            else
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "cannot change to directory");
        }
        getcwd(p,MAXPATHLEN+1);
        reply->set_path((const char*)p);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcImageViewer::done( ::grpc::ServerContext*,
                                          const ::google::protobuf::Empty*,
                                          ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        if (debug) {
            std::cout << "received grpc done( ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stdout);
        }

        static auto bye_bye = std::async( std::launch::async, [&]( ) { sleep(2); emit exit_now( ); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    bool grpcImageViewer::printraster( QtDisplayPanel *panel, const QString &type,
                                           const QString &file, double scale ) {
        QSize s = panel->canvasSize();
        int width = s.width();
        int height = s.height();

        width =  (int)(((double) width  * scale) + 0.5);
        height = (int)(((double) height * scale) + 0.5);

        // ensure that bitmap is redrawn to avoid displaying labeling within the plot
        if ( s.width() == width && s.height() == height ) {
            width += 1;
            height += 1;
        }

        std::promise<bool> prom;

        qtGO( [&]( ) {
                // resized
                QApplication::setOverrideCursor(Qt::WaitCursor);
                display::state::instance().beginFileOutputMode( );

                QSize oldSize = panel->size();
                QSize scaledSize = s;
                int dw = oldSize.width() - scaledSize.width(),
                    dh = oldSize.height() - scaledSize.height();
                scaledSize.scale(width, height, Qt::KeepAspectRatio);

                panel->setUpdateAllowed(false);
                // (Prevent display widget flashing during temporary resize.)
                panel->resize(scaledSize.width() + dw, scaledSize.height() + dh);
                QPixmap* mp = panel->contents();
                display::state::instance().endFileOutputMode( );
                panel->setUpdateAllowed(true);
                panel->resize(oldSize);
                QCoreApplication::processEvents();

                QApplication::restoreOverrideCursor();

                char *ctype = strdup(type.toLatin1().constData( ));
                if ( ! mp->save(file, ctype ) ) {
                    free( ctype );
                    delete mp;
                    prom.set_value(false);
                } else {
                    free( ctype );
                    delete mp;
                    prom.set_value(true);
                }
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        return fut.get( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    bool grpcImageViewer::printps( QtDisplayPanel *panel, const QString &type, const QString &file,
                                   int dpi, const QString &orientation, const QString &media ) {
        std::promise<std::pair<QSize,QRect>> prom;

        char eps_file_name[40];

        qtGO( [&]( ) {
                QPrinter *printer = new QPrinter;

#if QT_VERSION < 0x050000
                if ( type == "ps" || type == "eps" ) {
                    printer->setOutputFormat(QPrinter::PostScriptFormat);
                } else
#endif
                    if ( type == "pdf" ) {
                        printer->setOutputFormat(QPrinter::PdfFormat);
                    }

                printer->setResolution(dpi);
                if ( orientation == "landscape" ) {
                    printer->setOrientation(QPrinter::Landscape);
                } else {
                    printer->setOrientation(QPrinter::Portrait);
                }

                if ( media == "a4" ) {
                    printer->setPageSize(QPrinter::A4);
                } else {
                    printer->setPageSize(QPrinter::Letter);
                }

                if ( type == "eps" ) {
                    pid_t pid = getpid( );
                    sprintf( eps_file_name, "/tmp/eps-out.%06d", pid );
                    printer->setFullPage(true);
                    printer->setOutputFileName(eps_file_name);
                } else {
                    printer->setOutputFileName(file);
                }

                QPainter painter(printer);
                QRect viewport = painter.viewport();
                QRect rect = viewport;
                rect.adjust(72, 72, -72, -72);

                QSize sz = panel->size();

                int xs = sz.width();
                int ys = sz.height();
                QSize pSz(xs, ys);

                double ratio = 1;
                if ( orientation == "landscape" ) {
                    double rx = (double)rect.height() / xs;
                    double ry = (double)rect.width() / ys;
                    ratio = std::min(rx, ry);
                } else {
                    double rx = (double)rect.width() / xs;
                    double ry = (double)rect.height() / ys;
                    ratio = std::min(rx, ry);
                }

                pSz.setWidth(int(xs * ratio));
                pSz.setHeight(int(ys * ratio));

                QPixmap pmp(pSz);
                QApplication::setOverrideCursor(Qt::WaitCursor);
                panel->setAllowBackToFront(false);
                casacore::String backColor = "white";
                casacore::String foreColor = "black";
                panel->setBackgroundPS(backColor, foreColor);
                panel->beginLabelAndAxisCaching( );
                panel->resize(pSz);
                pmp = panel->getBackBuffer()->copy();
                painter.drawPixmap(0,0,pmp);
                panel->endLabelAndAxisCaching( painter );
                painter.end();
                panel->setBackgroundPS(backColor, foreColor);
                panel->setAllowBackToFront(true);
                panel->resize(sz);
                QApplication::restoreOverrideCursor();

                delete printer;
                prom.set_value(std::pair<QSize,QRect>(pmp.size( ),viewport));
            } );

        auto fut = prom.get_future( );
        fut.wait( );

        auto dim = fut.get( );
        if ( type == "eps" ) {
            char path[MAXPATHLEN+1];
            sprintf( path, "%s", file.toLatin1( ).constData( ) );
            adjusteps( eps_file_name, path, dim.first, dim.second );
            remove( eps_file_name );
        }

        return true;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    // adjust EPS bounding box...
    void  grpcImageViewer::adjusteps( const char *from, const char *to,
                                      const QSize &wcmax, const QRect &viewport ) {

        FILE *in = fopen( from, "r" );
        FILE *out = fopen( to, "w" );
        bool found = false;
        char buf[2049];
        while ( ! feof(in) ) {
            char *line = fgets( buf, 2049, in );
            if ( line ) {
                if ( ! found && ! strncmp( "%%BoundingBox: ", line, 15 ) ) {

                    float xmin, xmax, ymin, ymax;
                    int g = sscanf( line, "%%%%BoundingBox: %f %f %f %f", &xmin, &ymin, &xmax, &ymax );

                    if ( g != 4 ) {
                        fputs( line, out );
                    } else {
                        float ratio_y = ymax / float(viewport.height());
                        float ratio_x = xmax / float(viewport.width());

                        fprintf( out, "%%%%BoundingBox: 0 %d %d %d\n",
                                 int((ymax - (wcmax.height() * ratio_y)) + 1),
                                 int((wcmax.width() * ratio_x) + 1),
                                 int(ymax + 1) );
                    }

                    found = true;

                } else {
                    fputs( line, out );
                }
            }
        }
        fclose( out );
        fclose( in );
    }

}
