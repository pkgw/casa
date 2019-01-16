//# Registrar.h: maintain registry of services
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
#ifndef CASATOOLS_PROC_REGISTRAR_H
#define CASATOOLS_PROC_REGISTRAR_H
#include <list>
#include <mutex>
#include <string>

namespace casatools {   /** namespace for CASAtools classes within "CASA code" **/

    class ServiceId {
      public:
        ServiceId( std::string id_p, std::string type_p, std::string uri_p, unsigned int priority_p=0 ) : id_(id_p), type_(type_p), uri_(uri_p), priority_(priority_p) { }
        ServiceId( const ServiceId &other ) : id_(other.id_), type_(other.type_), uri_(other.uri_), priority_(other.priority_) { }
        ~ServiceId( ) { }

        std::string id( ) const { return id_; }
        std::string type( ) const { return type_; }
        std::string uri( ) const { return uri_; }
        int priority( ) const { return priority_; }

        operator std::string( ) const { return id_; }

        bool operator<( const ServiceId &other ) const { return id_ < other.id_; }
        bool operator>( const ServiceId &other ) const { return id_ > other.id_; }
        bool operator==( const ServiceId &other ) const { return id_ == other.id_; }
        bool operator!=( const ServiceId &other ) const { return id_ != other.id_; }
        bool operator<=( const ServiceId &other ) const { return id_ <= other.id_; }
        bool operator>=( const ServiceId &other ) const { return id_ >= other.id_; }

        bool operator<( const std::string &id_o ) const { return id_ < id_o; }
        bool operator>( const std::string &id_o ) const { return id_ > id_o; }
        bool operator==( const std::string &id_o ) const { return id_ == id_o; }
        bool operator!=( const std::string &id_o ) const { return id_ != id_o; }
        bool operator<=( const std::string &id_o ) const { return id_ <= id_o; }
        bool operator>=( const std::string &id_o ) const { return id_ >= id_o; }

      private:
        std::string id_;
        std::string type_;
        std::string uri_;
        int priority_;
    };
    
    class Registrar {
      public:

        Registrar( );

        virtual ~Registrar( );
        
        // get map of registrations
        std::list<ServiceId> services( ) {
            std::lock_guard<std::mutex> guard(service_list_mutex);
            return service_list;
        }

        // returns true if a registration for 'id' was found
        bool remove( std::string id );

        // returns assigned identifier (likely based upon the proposed_id)
        ServiceId add( const ServiceId &proposed );

        std::string uri( ) {
            std::lock_guard<std::mutex> guard(uri_mutex);
            return uri_;
        }
        
      private:
        std::mutex service_list_mutex;
        std::list<ServiceId> service_list;
        std::mutex uri_mutex;
        std::string uri_;
#ifdef USE_GRPC
        void *grpc_state;
#endif
    };

}

#endif
