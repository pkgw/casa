#define NO_INITIALIZE_STATICS 1
#include <coordsys_cmpt.h>
#undef NO_INITIALIZE_STATICS
#include <stdcasa/StdCasa/CasacSupport.h>
#include <casa/Utilities/CountedPtr.h>
#include <imageanalysis/ImageTypedefs.h>
#include <imageanalysis/ImageAnalysis/ImageDecimatorData.h>
#include <measures/Measures/Stokes.h>
#include <components/ComponentModels/ComponentType.h>

#include <set>

namespace casacore{

	class GaussianBeam;
	class ImageBeamSet;
	class ImageRegion;
	class LatticeBase;
	class LogIO;
	template<class T> class ImageStatistics;
	template<class T> class PtrHolder;
	template<class T> class SubImage;
	class LatticeExprNode;
    class String;
	class DirectionCoordinate;
}

namespace casa {
	template<class T> class ImageHistograms;
	template<class T> class ImageMetaData;
	template<class T> class ImageRegridderBase;
	template<class T> class ImageStatsCalculator;
	class SkyComponent;
}
