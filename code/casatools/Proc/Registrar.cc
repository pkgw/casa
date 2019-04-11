//# Registrar.cc: maintain registry of services
//# Copyright (C) 2017
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
#include <time.h>
#include <stdlib.h>
#include <map>
#include <algorithm>
#include <casatools/Proc/Registrar.h>
#ifdef USE_GRPC
#include <grpc++/grpc++.h>
#include "registrar.grpc.pb.h"
using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::Status;
#endif

using std::find_if;

namespace casatools {   /** namespace for CASAtools classes within "CASA code" **/

#ifdef USE_GRPC

    // generated stubs/base-classes are in casatools::rpc namespace...
    class grpcRegistrar final : public rpc::Registrar::Service {

        Status add( ServerContext* context, const rpc::ServiceId* req, rpc::ServiceId* reply) override {
            std::lock_guard<std::mutex> guard(registry_mutex);
            if ( registry == 0 || req->id( ).size( ) == 0 ||
                 req->type( ).size( ) == 0 || req->uri( ).size( ) == 0 )
                return Status::CANCELLED;
            ServiceId actual = registry->add(ServiceId(req->id( ),req->type( ),req->uri( )));
            reply->set_id(actual.id( ));
            reply->set_type(actual.type( ));
            reply->set_uri(actual.uri( ));
            return Status::OK;
        }

        Status remove( ServerContext* context, const rpc::ServiceId* req, ::google::protobuf::BoolValue* reply) override {
            std::lock_guard<std::mutex> guard(registry_mutex);
            if ( registry == 0 || req->id( ).size( ) == 0 )
                return Status::CANCELLED;
            reply->set_value(registry->remove(req->id( )));
            return Status::OK;
        }

        Status services( ServerContext* context, const ::google::protobuf::Empty* req, rpc::ServiceIds* reply) override {
            std::list<ServiceId> services;

            {
                std::lock_guard<std::mutex> guard(registry_mutex);
                if ( registry == 0 )
                    return Status::CANCELLED;
                services = registry->services( );
            }

            for ( std::list<ServiceId>::const_iterator ci = services.begin( ); ci != services.end( ); ++ci ) {
                rpc::ServiceId *serv = reply->add_service( );
                serv->set_id(ci->id( ));
                serv->set_type(ci->type( ));
                serv->set_uri(ci->uri( ));
            }
            
            return Status::OK;
        }

        std::mutex registry_mutex;
        Registrar *registry;

    public:
        grpcRegistrar( Registrar *r ) : registry(r) { }

    };

    struct grpcState {
        std::unique_ptr<Server> server;
        std::unique_ptr<grpcRegistrar> reg;
        ~grpcState( ) {
            if (getenv("GRPC_DEBUG")) fprintf(stdout, "stopping registry\n");
            fflush(stdout);
            if ( server ) server->Shutdown( );
        }
    };
    
#endif



    Registrar::Registrar( ) {
        srand(time(0));
#ifdef USE_GRPC
        grpc_state = 0;
        grpcState *state = new grpcState;
        state->reg.reset(new grpcRegistrar(this));
        ServerBuilder builder;
        char address_buf[100];
        constexpr char address_template[] = "0.0.0.0:%d";
        snprintf(address_buf,sizeof(address_buf),address_template,0);
        std::string server_address(address_buf);
        int selected_port = 0;

        // Listen on the given address without any authentication mechanism.
        builder.AddListeningPort(server_address, grpc::InsecureServerCredentials(), &selected_port);
        // Register "service" as the instance through which we'll communicate with
        // clients. In this case it corresponds to an *synchronous* service.
        builder.RegisterService(state->reg.get( ));
        // Finally assemble the server.
        state->server = builder.BuildAndStart( );
        if ( selected_port > 0 ) {
            // if an available port can be found, selected_port is set to a value greater than zero
            snprintf(address_buf,sizeof(address_buf),address_template,selected_port);
            uri_ = address_buf;
            if (getenv("GRPC_DEBUG")) std::cout << "registry available at " << uri_ << std::endl;
            grpc_state = (void*) state;
        } else delete state;

#endif
    }

    Registrar::~Registrar( ) {
#ifdef USE_GRPC
        if ( grpc_state ) delete (grpcState*) grpc_state;
#endif
    }
    
    bool Registrar::remove( std::string id ) {
        std::lock_guard<std::mutex> guard(service_list_mutex);
        auto search = find_if(service_list.begin(), service_list.end(), [=](const ServiceId &sid) { return sid == id; });
        if ( search == service_list.end( ) )
            return false;
        service_list.erase(search);
        return true;
    }

    ServiceId Registrar::add( const ServiceId &proposed ) {
        const int maxindex = 65535;
        char buf[8];
        sprintf( buf, ":%04x", rand( ) % maxindex + 1 );
        std::lock_guard<std::mutex> guard(service_list_mutex);
        auto search = find_if(service_list.begin(), service_list.end(), [=](const ServiceId &id) { return id == (proposed.id() + buf); });
        while ( search != service_list.end( ) ) {
            sprintf( buf, ":%04x", rand( ) % maxindex + 1 );
            search = find_if(service_list.begin(), service_list.end(), [=](const ServiceId &id) { return id == (proposed.id() + buf); });
        }
        ServiceId result(proposed.id() + buf, proposed.type(), proposed.uri());
        service_list.push_back(result);
        return result;
    }

}
