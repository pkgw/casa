#define NO_INITIALIZE_STATICS 1
#include <imagemetadata_cmpt.h>
#undef NO_INITIALIZE_STATICS
#include <stdcasa/StdCasa/CasacSupport.h>

#include <memory>

namespace casacore{
	class LogIO;
}

namespace casa 
{
    template<class T> class ImageMetaDataRW;
}
