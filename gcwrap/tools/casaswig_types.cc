#include "casaswig_types.h"

namespace casac {

    //***
    //*** In file included from /opt/casa/02/include/python2.7/unicodeobject.h:4:0,
    //*** from /opt/casa/02/include/python2.7/Python.h:85,
    //***    from coordsysPYTHON_wrap.cxx:174:
    //*** ../../../../src/gcwrap/tools/casaswig_types.h: In function ‘std::vector<T> casac::initialize_vector(int, T, ...) [with T = bool]’:
    //*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: warning: ‘bool’ is promoted to ‘int’ when passed through ‘...’
    //***    T val = va_arg(ap,T);
    //***                      ^
    //*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: note: (so you should pass ‘int’ not ‘bool’ to ‘va_arg’)
    //*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: note: if this code is reached, the program will abort
    //***
    template<>
    std::vector<bool> initialize_vector<bool>(int count, bool v1, ...) {
        va_list ap;
        va_start(ap, v1);
        std::vector<bool> result(count);
        result[0] = v1;
        for ( int i=1; i < count; ++i ) {
            bool val = va_arg(ap,int);
            result[i] = val;
        }
        return result;
    }

}
