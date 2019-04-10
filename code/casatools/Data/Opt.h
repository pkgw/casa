//# Opt.h: simple C++11 alternative to C++17 std::optional
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
#ifndef CASATOOLS_DATA_OPT_H_
#define CASATOOLS_DATA_OPT_H_
#include <stdexcept>

namespace casatools {

    template <typename T> class OptionValue {
    public:
        const T &get( ) {
            if ( ! has_value_ )
                throw std::runtime_error("attempt to retrieve value from none (option)");
            return value_;
        }
        bool has_value( ) { return has_value_; }
    private:
        template<typename U> friend class opt;
        OptionValue( ) : has_value_(false) { }
        OptionValue( const T &val ): has_value_(true), value_(val) { }
        bool has_value_;
        T value_;
    };

    template <typename T> class opt {
    public:
        static OptionValue<T> some(const T& val) { return OptionValue<T>(val); }
        static OptionValue<T> none( ) { return OptionValue<T>( ); }
    };

}

#endif
