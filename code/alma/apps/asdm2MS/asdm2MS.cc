// #ifdef _OPENMP
// #include <omp.h>
// #endif 
//#else
//  #define omp_get_num_threads() 0
//  #define omp_get_thread_num() 0
//#endif

#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

#define DDPRIORITY 1
//#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <regex>
#include <algorithm>
#include <vector>
#include <iomanip>

#include <stdcasa/optionparser.h>
#include <alma/Options/AlmaArg.h>
using namespace alma;

#include <alma/Enumerations/CWindowFunction.h>
#include <alma/ASDM/ASDMAll.h>
#include <alma/ASDM/Misc.h>

#include <alma/ASDMBinaries/SDMBinData.h>
using namespace sdmbin;

#include <exception>
using namespace asdm;
#include <alma/ASDM/IllegalAccessException.h>

#include <alma/apps/asdm2MS/UvwCoords.h>
#include <alma/apps/asdm2MS/ASDM2MSFiller.h>

#include <measures/Measures/Stokes.h>
#include <measures/Measures/MFrequency.h>
using namespace casacore;
#include <tables/Tables/Table.h>
#include <tables/Tables/PlainTable.h>
#include <tables/Tables/TableCopy.h>
#include <tables/Tables/TableInfo.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/OS/Path.h>
#include <casa/Quanta/Quantum.h>
#include <casa/BasicMath/Math.h>
#include <alma/Enumerations/CBasebandName.h>
#include <alma/Enumerations/CCalibrationDevice.h>
#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>
using namespace CalibrationDeviceMod;
#include <alma/Enumerations/CFrequencyReferenceCode.h>
#include <alma/Enumerations/CPolarizationType.h>
#include <alma/Enumerations/CProcessorSubType.h>
#include <alma/Enumerations/CProcessorType.h>
#include <alma/Enumerations/CScanIntent.h>
#include <alma/Enumerations/CSubscanIntent.h>
using namespace SubscanIntentMod;
#include <alma/Enumerations/CStokesParameter.h>

#include <alma/apps/asdm2MS/Name2Table.h>
#include <alma/apps/asdm2MS/ASDMVerbatimFiller.h>

using namespace asdmbinaries;
#include <alma/ASDMBinaries/SDMDataObjectReader.h>
#include <alma/ASDMBinaries/SDMDataObject.h>

#include <alma/ASDM/TableStreamReader.h>
#include <alma/apps/asdm2MS/asdm2MSGeneric.h>

#include <asdmstman/AsdmStMan.h>
#include <alma/apps/asdm2MS/BDF2AsdmStManIndex.h>

#include <alma/apps/asdm2MS/ScansParser.h>

#include <alma/apps/asdm2MS/ASDM2MSException.h>

#include	<time.h>
#if	defined(__sysv__)
#include	<sys/time.h>
#elif	defined(__bsd__)
#define		ftime	FTIME
#include	<sys/timeb.h>			/* from system */
#undef		ftime
extern	void	ftime( struct timeb * );	/* this is the funtion */
#else
#endif

#include <iostream>
#include <sstream>

using namespace std;

// The macro below was defined at the time when the code was changed to make the order of the baselines produced by the lazy filler  
// the same as the ones produced by the hard-working filler. 1 means yes same order, 0 means transposed order. The utilization of
// this macro should be withdrawn when the new version is validated.  
#define TRANSPOSE_BL_NUM 0

void	myTimer( double *cpu_time ,		/* cpu timer */
		 double *real_time ,		/* real timer */
		 int  *mode ) {			/* the mode */
  clock_t	tc;				/* clock time */
  double	ct;				/* cpu time in seconds */
  double	rt = 0.0 ;           		/* real time in seconds */

#if	defined(__sysv__)
  struct timeval 	Tp;
  struct timezone	Tzp;
#elif	defined(__bsd__)
  struct timeb tr;				/* struct from ftime */
#else
#endif
  tc = clock( );				/* get clock time */
  ct = (double)(tc) / (double)CLOCKS_PER_SEC;	/* to seconds */
#if	defined(__sysv__)
  gettimeofday( &Tp, &Tzp );			/* get timeofday */
  rt = (double) Tp.tv_sec + 0.000001 * (double) Tp.tv_usec;
#elif	defined(__bsd__)
  ftime( &tr );				/* get real time */
  rt = (double) tr.time + 0.001 * (double) tr.millitm;	/* in seconds */
#else
#endif
  if (*mode) {					/* calculate difference */
    (*cpu_time)  = ct - (*cpu_time);		/* cpu time */
    (*real_time) = rt - (*real_time);		/* real time */
  } else {
    (*cpu_time)  = ct;			/* set cpu time */
    (*real_time) = rt;			/* set real time */
  }
}

ASDM2MSFiller* msFiller;
string appName;
bool verbose = true;
bool isEVLA = false;
bool lazy = false;

double au_m = Quantity(1.0, "AU").getValue("m");  // AU in meters

void info (const string& message) {
  
  if (!verbose){
    return;
  }

  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void warning (const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void error(const string& message, int status=1) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
  //os << LogIO::POST;
  // cout << message << endl;
  exit(status);
}

/*
** A simplistic tracing toolbox.
*/
bool debug = (getenv("ASDM_DEBUG") != NULL);
vector<char> logIndent;
// #define LOGENTER(name) if (debug) {for_each(logIndent.begin(), logIndent.end(), cout << _1); logIndent.push_back('\t'); cout << #name ": entering" << endl;}
// #define LOGEXIT(name)  if (debug) {logIndent.pop_back(); for_each(logIndent.begin(), logIndent.end(), cout << _1); cout << #name ": exiting" << endl;}
// #define LOG(msg) if (debug) {for_each(logIndent.begin(), logIndent.end(), cout << _1); cout << msg << endl;}


ostringstream errstream;
ostringstream infostream;

template<typename T> string TO_STRING(const T &v) {
  char buffer[128];
  sprintf(buffer,"%d",v);
  return string(buffer);
}

string TO_STRING(const long unsigned int &v) {
    char buffer[128];
    sprintf(buffer,"%ld", v);
    return string(buffer);
}

string TO_STRING(const long int &v) {
    char buffer[128];
    sprintf(buffer,"%ld", v);
    return string(buffer);
}

string TO_STRING(const long long &v) {
    char buffer[128];
    sprintf(buffer,"%lld", v);
    return string(buffer);
}

string TO_STRING(const float &v) {
  const float shift = 100000000;
  const float rnd = .000000005;
  char buffer[128];
  int whole = int(v);
  int fraction = int((v-whole+rnd)*shift);
  sprintf(buffer,"%d.%08d",whole,fraction);
  return string(buffer);
}

string TO_STRING(const double &d) {
  char buffer[128];
  sprintf(buffer,"%g",d);
  return string(buffer);
}

// A facility to get rid of blanks at start and end of a string.
// 
string lrtrim(std::string& s,const std::string& drop = " ")
{
  std::string r=s.erase(s.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}

// These classes provide mappings from some ALMA Enumerations to their CASA counterparts.
class StokesMapper {
private :
  Stokes::StokesTypes* sa;

public :
  StokesMapper();
  ~StokesMapper();

  static Stokes::StokesTypes value(StokesParameterMod::StokesParameter s);
  Stokes::StokesTypes* to1DArray(const vector<StokesParameterMod::StokesParameter>& v);
  static vector<Stokes::StokesTypes> toVectorST(const vector<StokesParameterMod::StokesParameter>& v);
  static vector<int> toVectorI(const vector<StokesParameterMod::StokesParameter>& v);
}
  ;

StokesMapper::StokesMapper() {
  sa = 0;
}

StokesMapper::~StokesMapper() {
  if (sa)  delete[] sa;
}

Stokes::StokesTypes StokesMapper::value(StokesParameterMod::StokesParameter s) {
  switch (s) {
  case StokesParameterMod::I : return Stokes::I; 
  case StokesParameterMod::Q : return Stokes::Q; 
  case StokesParameterMod::U : return Stokes::U; 
  case StokesParameterMod::V : return Stokes::V; 
  case StokesParameterMod::RR : return Stokes::RR; 
  case StokesParameterMod::RL : return Stokes::RL; 
  case StokesParameterMod::LR : return Stokes::LR; 
  case StokesParameterMod::LL : return Stokes::LL; 
  case StokesParameterMod::XX : return Stokes::XX; 
  case StokesParameterMod::XY : return Stokes::XY; 
  case StokesParameterMod::YX : return Stokes::YX; 
  case StokesParameterMod::YY : return Stokes::YY; 
  case StokesParameterMod::RX : return Stokes::RX; 
  case StokesParameterMod::RY : return Stokes::RY; 
  case StokesParameterMod::LX : return Stokes::LX; 
  case StokesParameterMod::LY : return Stokes::LY; 
  case StokesParameterMod::XR : return Stokes::XR; 
  case StokesParameterMod::XL : return Stokes::XL; 
  case StokesParameterMod::YR : return Stokes::YR; 
  case StokesParameterMod::YL : return Stokes::YL; 
  case StokesParameterMod::PP : return Stokes::PP; 
  case StokesParameterMod::PQ : return Stokes::PQ; 
  case StokesParameterMod::QP : return Stokes::QP; 
  case StokesParameterMod::QQ : return Stokes::QQ; 
  case StokesParameterMod::RCIRCULAR : return Stokes::RCircular; 
  case StokesParameterMod::LCIRCULAR : return Stokes::LCircular; 
  case StokesParameterMod::LINEAR : return Stokes::Linear; 
  case StokesParameterMod::PTOTAL : return Stokes::Ptotal; 
  case StokesParameterMod::PLINEAR : return Stokes::Plinear; 
  case StokesParameterMod::PFTOTAL : return Stokes::PFtotal; 
  case StokesParameterMod::PFLINEAR : return Stokes::PFlinear; 
  case StokesParameterMod::PANGLE : return Stokes::Pangle;  
  }
  return Stokes::Undefined;
}

Stokes::StokesTypes* StokesMapper::to1DArray(const vector<StokesParameterMod::StokesParameter>& v) {
  if (v.size() == 0) return 0;

  if (sa) {
    delete[] sa;
    sa = 0;
  }

  sa = new Stokes::StokesTypes[v.size()];
  for (unsigned int i = 0; i < v.size(); i++) 
    sa[i] = value(v.at(i));

  return sa;
}

vector<int> StokesMapper::toVectorI(const vector<StokesParameterMod::StokesParameter>& v) {
  vector<int> result;

  for (unsigned int i = 0; i < v.size(); i++)
    result.push_back(value(v[i]));

  return result;
}

vector<Stokes::StokesTypes> StokesMapper::toVectorST(const vector<StokesParameterMod::StokesParameter>& v) {
  vector<Stokes::StokesTypes> result;

  for (unsigned int i = 0; i < v.size(); i++)
    result.push_back(value(v[i]));

  return result;
}

class FrequencyReferenceMapper {

public :
  FrequencyReferenceMapper();
  ~FrequencyReferenceMapper();

  static MFrequency::Types value(FrequencyReferenceCodeMod::FrequencyReferenceCode frc);
}
  ;

FrequencyReferenceMapper::FrequencyReferenceMapper() { }

FrequencyReferenceMapper::~FrequencyReferenceMapper() { }

MFrequency::Types FrequencyReferenceMapper::value(FrequencyReferenceCodeMod::FrequencyReferenceCode frc) {
  switch (frc) {
  case FrequencyReferenceCodeMod::LABREST : 
    errstream.str("");
    errstream << " Can't map FrequencyReferenceCode::LABREST to an MFrequency::Types" << endl;
    error(errstream.str());
    break;
    
  case FrequencyReferenceCodeMod::LSRD : return MFrequency::LSRD; 
  case FrequencyReferenceCodeMod::LSRK : return MFrequency::LSRK; 
  case FrequencyReferenceCodeMod::BARY : return MFrequency::BARY;
  case FrequencyReferenceCodeMod::REST : return MFrequency::REST;
  case FrequencyReferenceCodeMod::GEO  : return MFrequency::GEO; 
  case FrequencyReferenceCodeMod::GALACTO : return MFrequency::GALACTO;
  case FrequencyReferenceCodeMod::TOPO : return MFrequency::TOPO; 
  }
  return MFrequency::TOPO; // Never happens.
}

class PolTypeMapper {
private :
  vector<string> polType;

public :
  PolTypeMapper();
  ~PolTypeMapper();

  static char value(PolarizationTypeMod::PolarizationType p);
  vector<string> toStringVector(const vector<PolarizationTypeMod::PolarizationType>& v);
};


PolTypeMapper::PolTypeMapper() {
  ;
}

PolTypeMapper::~PolTypeMapper() {
  ;
}

char PolTypeMapper::value(PolarizationTypeMod::PolarizationType p) {
  return (CPolarizationType::name(p)).at(0);
}

vector<string> PolTypeMapper::toStringVector(const vector<PolarizationTypeMod::PolarizationType>& v) {
  polType.clear();
  for (unsigned int i = 0; i < v.size(); i++) {
    polType.push_back((CPolarizationType::name(v.at(i))));
  }
  return polType;
}


// These classes provide methods to convert from vectors of "anything" into vectors of basic types.
class FConverter {
public:

  static vector<float> toVectorF(const vector<vector<float> >& vvF, bool tranpose=false);

    template<class T> static vector<float> toVectorF(const vector<vector<T> >& vvT, bool transpose=false) {
    vector<float> result;

    if (transpose == false) {
      //
      // Simply linearize the vector of vectors.
      //
      for (unsigned int i = 0; i < vvT.size(); i++)
	for (unsigned int j = 0; j < vvT.at(i).size(); j++)
	  result.push_back(vvT.at(i).at(j).get());
    } else {
      //
      // We want to transpose.
      // Let's consider the very general case where our vector of vectors does not represent
      // a rectangular matrix, i.e. all the elements of this vector are vectors with possibly
      // different size. 
      //
      unsigned int maxsize = 0;
      unsigned int cursize = 0;
      for (unsigned int i = 0; i < vvT.size(); i++) {
	cursize = vvT.at(i).size();
	if (cursize > maxsize) maxsize = cursize;
      }
      
      for (unsigned int i = 0; i < maxsize; i++)
	for (unsigned int j = 0; j < vvT.size(); j++)
	  if (i < vvT.at(j).size())
	    result.push_back(vvT.at(j).at(i).get());     
    }
    return result;
  }
};

float d2f(double d) { return (float) d; }

vector<float> FConverter::toVectorF(const vector< vector <float> > & vvF, bool transpose) {
  vector<float> result;
  if (transpose == false) {
    //
    // Simply linearize the vector of vectors.
    //
    for (unsigned int i = 0; i < vvF.size(); i++)
      for (unsigned int j = 0; j < vvF.at(i).size(); j++)
	result.push_back(vvF.at(i).at(j));
  } else {
    //
    // We want to transpose.
    // Let's consider the very general case where our vector of vectors does not represent
    // a rectangular matrix, i.e. all the elements of this vector are vectors with possibly
    // different size. 
    //
    unsigned int maxsize = 0;
    unsigned int cursize = 0;
    for (unsigned int i = 0; i < vvF.size(); i++) {
      cursize = vvF.at(i).size();
      if (cursize > maxsize) maxsize = cursize;
    }
    
    for (unsigned int i = 0; i < maxsize; i++)
      for (unsigned int j = 0; j < vvF.size(); j++)
	if (i < vvF.at(j).size())
	  result.push_back(vvF.at(j).at(i));     
  }
  
  return result;
}   

class DConverter {
  
public :

  static vector<double> toVectorD(const vector<vector<double> >& vv);
  template<class T> static vector<double> toVectorD(const vector<T>& v) {
    vector<double> result;
    for (typename vector<T>::const_iterator iter = v.begin(); iter != v.end(); ++iter)
      result.push_back(iter->get());
    return result;
  }

  template<class T> static vector<double> toVectorD(const vector<vector<T> >& vv) {
    vector<double> result;
    for (typename vector<vector<T> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter)
      for (typename vector<T>::const_iterator iiter = iter->begin(); iiter != iter->end(); ++iiter)
	result.push_back(iiter->get());
    return result;
  }

  template<class T> static vector<vector<double> > toMatrixD(const vector<vector<T> >& vv) {
    vector<vector<double> > result;
    vector<double> vD;
    for( vector<T> v: vv ) {
      vD.clear();
      for( T x: v ) {
	vD.push_back(x.get());
      }
      result.push_back(vD);
    }
    return result;
  }
};

vector<double> DConverter::toVectorD(const vector<vector<double> >& vv) {
  vector<double> result;
  
  for (vector<vector<double> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter)
    result.insert(result.end(), iter->begin(), iter->end());

  return result;
}

class IConverter {
  
public :
  vector<int> static toVectorI(vector<Tag>& v);
};

vector<int> IConverter::toVectorI(vector<Tag>& v) {
  vector<int> result(v.size());
  vector<Tag>::const_iterator iiter = v.begin();
  vector<int>::iterator oiter = result.begin();
  for (; iiter != v.end(); ++iiter, ++oiter)
    *oiter = iiter->getTagValue();

  return result;
}

class SConverter {
public :
  template<class Enum, class EnumHelper>
  static vector<string> toVectorS(const vector<Enum>& vEnum) {
    vector<string> result(vEnum.size());
    typename vector<Enum>::const_iterator iiter = vEnum.begin();
    vector<string>::iterator     oiter = result.begin();

    for (; iiter != vEnum.begin(); ++iiter, ++oiter)
      *oiter = EnumHelper::name(*iiter);

    return result;
  }
};

class  CConverter {
public :
  static vector<std::complex<float> > toVectorCF(const vector<vector<asdm::Complex> >& vv);
};

vector<std::complex<float> > CConverter::toVectorCF(const vector<vector<asdm::Complex> >& vv) {
  vector<std::complex<float> > result;

  for (vector<vector<asdm::Complex> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter) 
    result.insert(result.end(), iter->begin(), iter->end());

  return result;
}

class ComplexDataFilter {

public:
  ComplexDataFilter();
  virtual ~ComplexDataFilter();
  virtual float* to4Pol(int numPol, int numChan, float* cdata);

private:
  vector<float *> storage; 
  vector<float> storage_v;
};

ComplexDataFilter::ComplexDataFilter() { }

ComplexDataFilter::~ComplexDataFilter() {
  for (unsigned int i = 0; i < storage.size(); i++) 
    delete[] storage.at(i);
}

float *ComplexDataFilter::to4Pol(int numCorr, int numChan, float* cdata) {
  // Do nothing if numCorr != 3
  if (numCorr != 3) return cdata;

  // Allocate storage for 2 * numCorr * numChan with numCorr changed from 3 to 4.
  float* filtered = new float[ 2 * 4 * numChan ];

  storage.push_back(filtered);

  for (int i = 0; i < numChan; i++) {
    // The 1st row goes to the first row.
    filtered[ 8 * i ]     = cdata[ 6 * i ] ;
    filtered[ 8 * i + 1 ] = cdata[ 6 * i + 1 ] ;

    // The second row goes the second row.
    filtered [ 8 * i + 2 ] = cdata[ 6 * i + 2 ] ;
    filtered [ 8 * i + 3 ] = cdata[ 6 * i + 3 ] ;

    // The second row's conjugate goes to the third row.
    filtered [ 8 * i + 4 ] = cdata[ 6 * i + 2 ] ;
    filtered [ 8 * i + 5 ] = - cdata[ 6 * i + 3 ] ;

    // The third row goes to the third row.
    filtered [ 8 * i + 6 ] = cdata[ 6 * i + 4 ] ;
    filtered [ 8 * i + 7 ] = cdata[ 6 * i + 5 ] ;    
  }
  return filtered;
}

//
//  A collection of functions to combine target and offset
//  into a pointing direction in the Pointing table.
//
/**
 *
 * From a Fortran 90 subroutine, written by Didier Despois (1980)
 * and revised by Michel Perault (1984).
 *
 * Computes a coodinates conversion matrix for a rotation defined
 * by Euler angles psi, the and phi (values in radians).
 *
 * The product of the euclidean coordinates of a vector in the original
 * base by this matrix gives the coordinates in the rotated base.
 * 
 */
void eulmat(double psi, double the, double phi, vector<vector<double> >& mat) {
  double cpsi, spsi, cthe, sthe, cphi, sphi;
  double x1, x2;

  cpsi = cos(psi);
  spsi = sin(psi);
  cthe = cos(the);
  sthe = sin(the);
  cphi = cos(phi);
  sphi = sin(phi);

  x1 = spsi*cthe;
  x2 = cpsi * cthe; 
  mat.at(0).at(0) =  cpsi * cphi - x1 * sphi;
  mat.at(0).at(1) =  spsi * cphi + x2 * sphi;
  mat.at(0).at(2) =  sthe * sphi;

  mat.at(1).at(0) = -cpsi * sphi - x1 * cphi;
  mat.at(1).at(1) = -spsi * sphi + x2 * cphi;
  mat.at(1).at(2) = sthe * cphi;

  mat.at(2).at(0) =  spsi * sthe;
  mat.at(2).at(1) =  -cpsi * sthe;
  mat.at(2).at(2) =  cthe;
}

/**
 * Performs a product Matrix x Vector.
 *
 * Attention ! vFactor , mFactor and vResult must be correctly
 * sized ! And also vFactor and vResult must be different objects.
 */
void matvec (const vector<vector<double> >& mFactor, const vector<double>& vFactor, vector<double>& vResult) {
  vResult.at(0) = vFactor.at(0)*mFactor.at(0).at(0) + vFactor.at(1)*mFactor.at(0).at(1) + vFactor.at(2)*mFactor.at(0).at(2);
  vResult.at(1) = vFactor.at(0)*mFactor.at(1).at(0) + vFactor.at(1)*mFactor.at(1).at(1) + vFactor.at(2)*mFactor.at(1).at(2);
  vResult.at(2) = vFactor.at(0)*mFactor.at(2).at(0) + vFactor.at(1)*mFactor.at(2).at(1) + vFactor.at(2)*mFactor.at(2).at(2);
}

/**
 * Spherical to cartesian conversion.
 * Radius assumed to be equal to 1.
 *
 * Attention ! a and x must be correctly sized !
 */
void rect(const vector<double>& s, vector<double>& x) {
  x.at(0) = cos(s.at(0)) * cos(s.at(1));
  x.at(1) = sin(s.at(0)) * cos(s.at(1));
  x.at(2) = sin(s.at(1));
}

/**
 * Cartesian to spherical conversion.
 * Radius assumed to be equal to 1.
 *
 * Attention ! a and x must be correctly sized !
 */
void spher(const vector<double>& x, vector<double>& s) {
  s.at(0) = atan2(x.at(1) , x.at(0));
  s.at(1) = atan2(x.at(2), sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1)));
}


/** 
 * a rotation matrix from local (topocentric) coordinates to ITRF geocentric coordinates
 * lambda - longitude, phi - latitude 
 */
void topo2geomat(double lambda, double phi, vector<vector<double> >& mat) {
 
  double clam, slam, cphi, sphi;

  clam = cos(lambda);
  slam = sin(lambda);
  cphi = cos(phi);
  sphi = sin(phi);

  mat.at(0).at(0) =  -slam;
  mat.at(0).at(1) =  -sphi * clam;
  mat.at(0).at(2) = cphi * clam; 

  mat.at(1).at(0) = clam;
  mat.at(1).at(1) = -sphi * slam;
  mat.at(1).at(2) = cphi * slam;

  mat.at(2).at(0) =  0; 
  mat.at(2).at(1) =  cphi;
  mat.at(2).at(2) =  sphi;
}

/**
 * Reorder the collection of spectral windows ids.
 *
 * Given a dataset 'ds', the statement 'ds.spectralWindowTable.get()' returns
 * a vector of pointers on instances of 'SpectralWindowRow' and 
 * defines implicitely an ordered collection of Tag (of type TagType::SpectralWindow) instances
 * obtained by using 'getSpectralWindowId()' on each instance.
 *
 * This method partitions this ordered collection into two collections
 * based on the criterium : does the Tag appear in at least one row
 * of the DataDescription table in the SpectralWindowId attribute or not.
 *
 * The methods returns a vector of (SpactralWindow) Tag  obtained by appending the collections
 * of (SpectralWindow) Tag which appear in the DataDescription table
 * followed by the collection of (SpectralWindow) Tag which does not appear
 * in the DataDescription table. The int value which associated with a Tag
 * is the order number of the insertion of the pair (Tag, int) in the map.
 * 
 */

struct TagCmp {
  bool operator() (const Tag& t1, const Tag& t2) const {
    return t1.getTagValue() < t2.getTagValue();
  }
};

vector<Tag> reorderSwIds(const ASDM& ds) {
  vector<SpectralWindowRow *> swRs = ds.getSpectralWindow().get();
  vector<DataDescriptionRow *> ddRs = ds.getDataDescription().get();
  map<Tag, bool, TagCmp> isInDD;

  for (vector<SpectralWindowRow *>::size_type  i = 0; i < swRs.size(); i++) isInDD[swRs[i]->getSpectralWindowId()] = false;
  for (vector<DataDescriptionRow *>::size_type i = 0; i < ddRs.size(); i++) isInDD[ddRs[i]->getSpectralWindowId()] = true;

  vector<Tag> swIdsDD, swIdsNoDD;
  for (map<Tag, bool, TagCmp>::iterator iter = isInDD.begin(); iter != isInDD.end(); ++iter)
    if (iter->second) swIdsDD.push_back(iter->first);
    else swIdsNoDD.push_back(iter->first);

  vector<Tag> result (swIdsDD.begin(), swIdsDD.end());
  for (vector<Tag>::size_type i = 0; i < swIdsNoDD.size(); i++)
    result.push_back(swIdsNoDD[i]);

  return result;
}

template<class T> void checkVectorSize(const string& vectorAttrName, const vector<T>& vectorAttr,
                                       const string& sizeAttrName, unsigned int sizeAttr,
                                       const string& tableName, unsigned int rowNumber) {
  if (vectorAttr.size() != sizeAttr) {
    errstream.str("");
    errstream << "In the '"
	      << tableName 
	      << " table, at row #"
	      << rowNumber
	      << ", I found '"
	      << vectorAttrName
	      << "' with a size of '"
	      << vectorAttr.size() 
	      << "', I was expecting it to be equal to the size of '"
	      << sizeAttrName
	      <<"' which is '"
	      <<"'"
	      <<sizeAttr
	      <<"'. I can't go further."
	      << endl;
    error(errstream.str());
  }
}

EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection> apcLiterals(const ASDM& ds) {
    EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection> result;

  vector<MainRow *> mRs = ds.getMain().get();
  
  for (unsigned int i = 0; i < mRs.size(); i++) {
    ConfigDescriptionRow * configDescriptionRow = mRs.at(i)->getConfigDescriptionUsingConfigDescriptionId();
    vector<AtmPhaseCorrectionMod::AtmPhaseCorrection> apc = configDescriptionRow -> getAtmPhaseCorrection();
    for (unsigned int i = 0; i < apc.size(); i++)
      result.set(apc.at(i));
  }
  return result;
}

bool hasCorrectedData(const EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection>& es) {
    return es[AtmPhaseCorrectionMod::AP_CORRECTED];
}

bool hasUncorrectedData(const EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection>& es) {
    return es[AtmPhaseCorrectionMod::AP_UNCORRECTED];
}
		    
//
// A number of EnumSet to encode the different selection criteria.
//
EnumSet<CorrelationModeMod::CorrelationMode>         es_cm;
EnumSet<SpectralResolutionTypeMod::SpectralResolutionType>  es_srt;
EnumSet<TimeSamplingMod::TimeSampling>            es_ts;
Enum<CorrelationModeMod::CorrelationMode>            e_query_cm; 
EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection>      es_query_apc;    

// An EnumSet to store the different values of AtmPhaseCorrection present
// in the binary data (apc in datastruct).
//
EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection>      es_apc;

//
// By default the resulting MS will not contain compressed columns
// unless the 'compress' option has been given.
// 
bool                             withCompression = false;

//
// A function to determine if overTheTop is present in a given row of the Pointing table.
//
bool overTheTopExists(PointingRow* row) { return row->isOverTheTopExists(); }


map<int, int> swIdx2Idx ;                       // A map which associates old and new index of Spectral Windows before/after reordering.

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Some functions defined for the processing of the SysPower table 
// with 'functional programming' techniques.
//
int sysPowerAntennaId(SysPowerRow* row) {
  return row->getAntennaId().getTagValue();
}

int sysPowerSpectralWindowId(const SysPowerRow* row) {
  return swIdx2Idx[row->getSpectralWindowId().getTagValue()];
}

int sysPowerFeedId(const SysPowerRow* row) {
  return row->getFeedId();
}

double sysPowerMidTimeInSeconds(const SysPowerRow* row) {
  //if (isEVLA) {
  //  return row->getTimeInterval().getStartInMJD()*86400 ; 
  //}
  return row->getTimeInterval().getStartInMJD()*86400 + ((double) row->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond / 2.;
}

double sysPowerIntervalInSeconds(const SysPowerRow* row) {
  return ((double) row->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
}

int sysPowerNumReceptor(const SysPowerRow* row) {
  return row->getNumReceptor();
}


struct sysPowerSwitchedPowerDifference {
private:
  vector<float>::iterator iter;
  
public:
  sysPowerSwitchedPowerDifference(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getSwitchedPowerDifference();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};

struct sysPowerSwitchedPowerSum {
private:
  vector<float>::iterator iter;
  
public:
  sysPowerSwitchedPowerSum(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getSwitchedPowerSum();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};


struct sysPowerRequantizerGain {
private:
  vector<float>::iterator iter;

public:
  sysPowerRequantizerGain(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getRequantizerGain();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};

/** 
 * This function tries to allocate the number of integrations into slices
 * based on the average size of one integration and the requested maximum size of one slice.
 * The sum of all of the values in the returned vector is always equal to nIntegrations.
 * @parameter nIntegrations the number of integrations in the BDF
 * @parameter bdfSize the size of the BDF, in bytes
 * @parameter approxSizeInMemory the approximate size that one wants for one slice, in bytes.
 */
vector<unsigned int> getIntegrationSlices(int nIntegrations, uint64_t bdfSize, uint64_t approxSizeInMemory) {
  if (debug) cout << "getIntegrationSlices: entering" << endl;
  vector<unsigned int> result;
  if (nIntegrations > 0) {
    if (bdfSize < approxSizeInMemory) {
      result.push_back(nIntegrations);
    } else {
      uint64_t avgIntSize = bdfSize/nIntegrations;
      unsigned int nIntPerSlice = approxSizeInMemory/avgIntSize;
      unsigned int nSlice = nIntegrations/nIntPerSlice;
      result.resize(nSlice,nIntPerSlice);
      unsigned int residualInts = nIntegrations % nIntPerSlice;
      if (residualInts > (nSlice*nIntPerSlice/5)) {
	// if the number of left over integrations is more than 1/5 of what's in the other slices
	// then this makes an acceptable last slice - no need to redistribute
	result.push_back(residualInts);
      } else {
	// final slice would too small, redistribute it to the other slices
	while(residualInts > 0) {
	  for (unsigned int i = 0; residualInts > 0 && i < result.size(); i++) {
	    result[i]++; 
	    residualInts--;
	  }
	}
      }
    }
  }
  if (debug) cout << "getItegrationSlices: exiting" << endl;
  return result;
}

vector<map<AtmPhaseCorrectionMod::AtmPhaseCorrection,ASDM2MSFiller*> >  msFillers_v;

map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*> msFillers; // There will be one filler per value of the axis APC.

vector<int>	dataDescriptionIdx2Idx;
int		ddIdx;
map<MainRow*, int>     stateIdx2Idx;

set<int> SwIdUsed;

double radian2degree(double radian) {
  return radian / M_PI * 180.0;
}

double m2au(double m) {
  return m / au_m;
}

double mpers2auperd(double mpers) {
  return m2au(mpers) * 24. * 3600.;
}

#define LOG_EPHEM(message) if (getenv("FILLER_LOG_EPHEM")) cout << message;

void linearInterpCoeff(uint32_t npoints, const vector<double>& time_v, const vector<double>& k_v, vector<vector<double> >& coeff_vv) {
  LOGENTER("linearInterpCoeff");
  coeff_vv.clear();
  coeff_vv.resize(npoints-1);

  LOG_EPHEM("linearInterpCoeff for npoints = " + TO_STRING(npoints));
  
  for (uint32_t i = 0; i < npoints-1; i++) {
    vector<double> coeff_v (2);
    coeff_v[0] = k_v[i];
    coeff_v[1] = (k_v[i+1] - k_v[i]) / (time_v[i+1] - time_v[i]);
    coeff_vv[i] = coeff_v;

    LOG_EPHEM(TO_STRING(i) + " : " + TO_STRING(k_v[i]) + ", " + TO_STRING(k_v[i+1]) + ", " + TO_STRING(time_v[i]) + ", " + TO_STRING(time_v[i+1]) + ": " + TO_STRING(coeff_v[0]) + ", " + TO_STRING(coeff_v[1]) + "\n");
  }
  LOGEXIT("linearInterpCoeff");
}

double evalPoly (unsigned int numCoeff, const vector<double>& coeff, double timeOrigin, double time) {
  LOGENTER("evalPoly");
  LOG( "numCoeff=" + TO_STRING(numCoeff) + ", size of coeff=" + TO_STRING(coeff.size())) ;
  LOG( "time=" + TO_STRING(time) + ", timeOrigin=" + TO_STRING(timeOrigin));
  //
  // Let's use the Horner schema to evaluate the polynomial.
  double result = coeff[numCoeff-1];
  for (int i = numCoeff - 2; i >= 0; i--) 
    result = coeff[i] + result*(time-timeOrigin);
  LOG (" result= " + TO_STRING(result));
  LOGEXIT("evalPoly");
  return result;
}

void deleteEphemeris(map<AtmPhaseCorrectionMod::AtmPhaseCorrection, Table*>& apc2EphemTable_m) {
    for ( map<AtmPhaseCorrectionMod::AtmPhaseCorrection, Table*>::iterator iter = apc2EphemTable_m.begin(); iter!=apc2EphemTable_m.end(); ++iter)
    if (apc2EphemTable_m[iter->first] != NULL) delete apc2EphemTable_m[iter->first];
}

std::map<int, double>ephemStartTime_m;
void fillEphemeris(ASDM* ds_p, uint64_t timeStepInNanoSecond, bool interpolate_ephemeris, string telescopeName) {
  LOGENTER("fillEphemeris");

  // division by timeStepInNanoSecond below causes FPE - ensure it's set to something non-zero
  // default to 0.001s = 1e6 ns
  timeStepInNanoSecond = (timeStepInNanoSecond == 0 ? 1e6 : timeStepInNanoSecond);
  
  try {
    // Retrieve the Ephemeris table's content.
    EphemerisTable & eT = ds_p->getEphemeris();
    const vector<EphemerisRow *>&  eR_v = eT.get();
  
    infostream.str("");
    infostream << "The dataset has " << eT.size() << " ephemeris row(s)...";
    info(infostream.str());

    // Let's partition the vector of Ephemeris rows based on the value of the field ephemerisId.
    vector<EphemerisRow *> empty_v;
    map<int, vector<EphemerisRow *> > i2e_m; // A map which associates a value of ephemerisId to the vector of Ephemeris rows 
    // having this value in their field ephemerisId.
    set<int> ephemerisId_s;
    for( EphemerisRow * eR_p: eR_v ) {
      int ephemerisId = eR_p -> getEphemerisId();
      if (i2e_m.find(ephemerisId) == i2e_m.end()) {
	ephemerisId_s.insert(ephemerisId);
	i2e_m[ephemerisId] = empty_v;
      }
      i2e_m[ephemerisId].push_back(eR_p);
    }

    // Let's create and fill the MS ephemeris tables.

    for(int ephemerisId: ephemerisId_s) {
      /**
       * Check if there is at least one ASDM::Field row refering to this ephemerisId.
       */
      const vector<FieldRow *> fR_v = ds_p->getField().get();
      vector<FieldRow *> relatedField_v;    // This vector elements will contain pointers to all the Fields refering to this ephemerisId.
      for(const FieldRow* fR_p: fR_v) {
	if (fR_p->isEphemerisIdExists() && (fR_p->getEphemerisId() == ephemerisId))
	  relatedField_v.push_back(const_cast<FieldRow *>(fR_p));
      }
    
      if (relatedField_v.size() == 0) {
	infostream.str("");
	infostream << "No ASDM Field found with 'ephemerisId=='"<< ephemerisId << "'. An ephemeris table will be created in the MS though, but it will not be connected a Field.";
	info(infostream.str());
      }

      // Let's use the value of 'fieldName' in first element (if any) of relatedField_v to build the name of the MS Ephemeris table to
      // be created.
      //
      string fieldName;
      if (relatedField_v.size() > 0)
	fieldName = relatedField_v[0]->getFieldName();

      /**
       * For each MS ephemeris table we need :
       *
       * 1) MJD0
       * 2) DMJD0
       * 3) NAME
       * 4) GeoLong in deg.
       * 5) GeoLat  in deg.
       * 6) GeoDist in km.
       * 7) origin (a textual summary of the origin of this ephemeris)
       *
       * We derive these values from the informations found in the first element of i2e_m[ephemerisId]
       */
      vector<EphemerisRow *>&	ephRow_v	 = i2e_m[ephemerisId];    
      vector<double>		observerLocation = ephRow_v[0]->getObserverLocation();
      double			geoLong		 = radian2degree(observerLocation[0]); // in order to get degrees.
      double			geoLat		 = radian2degree(observerLocation[1]); // in order to get degrees.
      double                    geoDist          = observerLocation[2] / 1000.0;             // in order to get km (supposedly above the reference ellipsoid)

      int64_t	t0ASDM = ephRow_v[0]->getTimeInterval().getStart().get();	// The first time recorded for this ephemerisId.
      int64_t	q      = t0ASDM / timeStepInNanoSecond;
      int64_t	r      = t0ASDM % timeStepInNanoSecond;
      int64_t	t0MS   = t0ASDM;
      if ( r != 0 ) {  
	q = q + 1;
	t0MS = q * timeStepInNanoSecond;
      }

      double mjd0 = ArrayTime(t0MS).getMJD();
     
      double dmjd = interpolate_ephemeris ? 0.001 : ephRow_v[0]->getTimeInterval().getDuration().get() / 1000000000LL / 86400.0;
      // Grid time step == 0.001 if ephemeris interpolation requested
      // otherwise == the interval of time of the first element of ephemeris converted in days.
      // *SUPPOSEDLY* constant over all the ephemeris. 
 
      // determine the position reference system
      double equator =  ephRow_v[0]->getEquinoxEquator();
      string posref = "unknown";
      if (equator == 2000.) { // the Ephemeris table presently only stores the equator
	posref = "ICRF/ICRS";
      }

      // Prepare the table keywords with the values computed above.
      TableDesc tableDesc;
    
      tableDesc.comment() = ephRow_v[0]->getOrigin();
      time_t now = time(0);
      struct tm tstruct;
      char buf[80];

      tstruct = *gmtime(&now);

      strftime(buf, sizeof(buf), "%Y/%m/%d/%H:%M", &tstruct);
      string creationDate(buf);
      
      tableDesc.rwKeywordSet().define("VS_CREATE", creationDate);
      tableDesc.rwKeywordSet().define("VS_DATE", creationDate);
      tableDesc.rwKeywordSet().define("VS_VERSION", "0001.0001");
      tableDesc.rwKeywordSet().define("VS_TYPE", "Table of comet/planetary positions");
      tableDesc.rwKeywordSet().define("MJD0", casacore::Double(mjd0));
                // Actually this value is temporary, it'll be updated when the table will be populated.
      tableDesc.rwKeywordSet().define("dMJD", casacore::Double(dmjd));
      tableDesc.rwKeywordSet().define("NAME", fieldName);
      tableDesc.rwKeywordSet().define("GeoLong", casacore::Double(geoLong));
      tableDesc.rwKeywordSet().define("GeoLat", casacore::Double(geoLat));
      tableDesc.rwKeywordSet().define("GeoDist", casacore::Double(geoDist));
      if (geoDist == 0.0) 
	tableDesc.rwKeywordSet().define("obsloc", "GEOCENTRIC");
      else
	tableDesc.rwKeywordSet().define("obsloc", telescopeName);
    
      tableDesc.rwKeywordSet().define("posrefsys", posref);

      // Then the fields definitions and keywords.
      ScalarColumnDesc<casacore::Double> mjdColumn("MJD");
      mjdColumn.rwKeywordSet().define("UNIT", "d");
      tableDesc.addColumn(mjdColumn);
    
      ScalarColumnDesc<casacore::Double> raColumn("RA");
      raColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(raColumn);
    
      ScalarColumnDesc<casacore::Double> decColumn("DEC");
      decColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(decColumn);
    
      ScalarColumnDesc<casacore::Double> rhoColumn("Rho");
      rhoColumn.rwKeywordSet().define("UNIT", "AU");
      tableDesc.addColumn(rhoColumn);
    
      ScalarColumnDesc<casacore::Double> radVelColumn("RadVel");
      radVelColumn.rwKeywordSet().define("UNIT", "AU/d");
      tableDesc.addColumn(radVelColumn);
    
      ScalarColumnDesc<casacore::Double> diskLongColumn("diskLong");
      diskLongColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(diskLongColumn);
    
      ScalarColumnDesc<casacore::Double> diskLatColumn("diskLat");
      diskLatColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(diskLatColumn);
  
      string tableName = "EPHEM"
	+ TO_STRING(ephemerisId)
	+ "_"
	+ fieldName
	+ "_"
	+ TO_STRING(mjd0)
	+ ".tab";

      // The following characters are replaced with "_" : (){}[]/ and " " (space)
      std::regex e("[\\[\\]\\(\\)\\{\\}\\/ ]");
      tableName = std::regex_replace(tableName, e, std::string("_"));
      
      map<AtmPhaseCorrectionMod::AtmPhaseCorrection, Table*> apc2EphemTable_m;
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	string tablePath = (string) iter->second->ms()->tableName() + string("/FIELD/") + tableName;  
	SetupNewTable tableSetup(tablePath, tableDesc, Table::New);
    
	Table * table_p = new Table(tableSetup, TableLock(TableLock::PermanentLockingWait));
	TableInfo& info = table_p->tableInfo();
	info.setType("IERS");
	info.setSubType("Solar System");
	info.readmeAddLine("generated by asdm2ms");

	AlwaysAssert(table_p, AipsError);
	(const_cast<casacore::MeasurementSet*>(iter->second->ms()))->rwKeywordSet().defineTable(tableName, *table_p);
	table_p->flush();
	apc2EphemTable_m[iter->first] = table_p;
	LOG("Empty ephemeris table '" + tablePath + "' created");
      }

      // 
      // Now it's time to fill the EPHEM table(s).
      //
      // Below the vectors which will be used to populate slices of Ephemeris tables in the MSes.
      //
      vector<double> mjdMS_v;
      vector<double> raMS_v;
      vector<double> decMS_v;
      vector<double> distanceMS_v;
      vector<double> radVelMS_v;

      bool	numPolyDirIsOne	   = ephRow_v[0]->getNumPolyDir() == 1;
      bool	numPolyDistIsOne   = ephRow_v[0]->getNumPolyDist() == 1;
      bool	radVelExists	   = ephRow_v[0]->isRadVelExists() && ephRow_v[0]->isNumPolyRadVelExists();
      bool	numPolyRadVelIsOne = radVelExists ? ephRow_v[0]->getNumPolyRadVel() == 1 : false; 
      bool      anyNumPolyIsOne = numPolyDirIsOne || numPolyDistIsOne || (radVelExists && numPolyRadVelIsOne);
      bool      allNumPolyIsOne = numPolyDirIsOne && numPolyDistIsOne && (!radVelExists || numPolyRadVelIsOne);

      LOG ("numPolyDirIsOne = " + TO_STRING(numPolyDirIsOne));
      LOG ("numPolyDistIsOne = " + TO_STRING(numPolyDistIsOne));
      LOG ("radVelExists = " + TO_STRING(radVelExists));
      LOG ("numPolyRadVelIsOne = " + TO_STRING(numPolyRadVelIsOne));
      LOG ("anyNumPolyIsOne = " + TO_STRING(anyNumPolyIsOne));
      LOG ("allNumPolyIsOne = " + TO_STRING(allNumPolyIsOne));

      // the single row case
      if (ephRow_v.size()==1) {

	infostream.str("");
	infostream << "The MS Ephemeris table for ephemerisId = '" << ephemerisId
		   << "' will be produced by tabulating the polynomials found in the single row for 'dir', 'distance' and optionally 'radVel' with a timestep of '"
		   << timeStepInNanoSecond / 1.e9 / 86400. << "' days"; 
	info(infostream.str());
	
	dmjd =  timeStepInNanoSecond / 1.e9 / 86400.; // to be written in the keyword DMJD
	
	//
	// Calculate the grid of times where the polynomials will be tabulated.
	// This grid contains all the times which :
	// - are multiple of the tabulation time steps,
	// - contained in the arraytime interval of validity of the current ASDM Ephemeris row.
	// 
	// At this point times are expressed in nanoseconds.
	//

	// this case has just one element in ephRow_v

	int64_t	tstartASDM = ephRow_v[0]->getTimeInterval().getStart().get();
	int64_t	tendASDM   = tstartASDM + ephRow_v[0]->getTimeInterval().getDuration().get();
	int64_t	q	   = tstartASDM / timeStepInNanoSecond;
	int64_t	r	   = tstartASDM % timeStepInNanoSecond;
	int64_t	tstartMS   = tstartASDM;
	  
	if ( r!= 0 ) {
	  q = q + 1;
	  tstartMS = q * timeStepInNanoSecond;
	}
	  
	vector<int64_t> tabulation_time_v;
	int64_t t = tstartMS - timeStepInNanoSecond;  // One extra timestep before the beginning of the interval of validity.
	do {
	  tabulation_time_v.push_back(t);
	  t += timeStepInNanoSecond;
	} while (t <= tendASDM);
	tabulation_time_v.push_back(t); // One extra timestep after the end of the interval of validity. 
	  
	//
	// Tabulate the MS Ephemeris columns for each tabulation time, in s. 
	//
	double	timeOrigin = ephRow_v[0]->getTimeOrigin().get() * 1.0e-09 / 86400. ; // to "days"
	  
	// Convert from `radians to degrees.
	const vector<vector<double> >& dir_v =  ephRow_v[0]->getDir();
	vector<double>	ra_coeff_v;
	vector<double>	dec_coeff_v;
	for (unsigned int idir = 0; idir < dir_v.size(); idir++) {
	  ra_coeff_v.push_back(dir_v[idir][0]);
	  dec_coeff_v.push_back(dir_v[idir][1]);
	}
				 
	LOG(" There are " + TO_STRING(ra_coeff_v.size()) + " RA coeffs and " + TO_STRING(dec_coeff_v.size()) + " DEC coeffs");
	  
	// Convert radian to degree.
	std::transform(ra_coeff_v.begin(), ra_coeff_v.end(), ra_coeff_v.begin(), radian2degree);
	std::transform(dec_coeff_v.begin(), dec_coeff_v.end(), dec_coeff_v.begin(), radian2degree);
	  
	// Convert from m to AU.
	vector<double> distance_coeff_v = ephRow_v[0]->getDistance();
	std::transform(distance_coeff_v.begin(), distance_coeff_v.end(), distance_coeff_v.begin(), m2au);
	  
	vector<double> radvel_coeff_v;
	if (ephRow_v[0]->isRadVelExists()) {
	  radvel_coeff_v = ephRow_v[0]->getRadVel();
	  // Convert from m per s to AU per day.
	  std::transform(radvel_coeff_v.begin(), radvel_coeff_v.end(), radvel_coeff_v.begin(), mpers2auperd);
	}
	
	// And proceed...
	LOG ("There will be " + TO_STRING(tabulation_time_v.size()) + " time steps used to tabulate the polynomials.");
	for (unsigned int itab = 0; itab < tabulation_time_v.size(); itab++) {
	  double tabulation_time = tabulation_time_v[itab] * 1.0e-09 / 86400.0 ;  // It appeared that times should be expressed in "day" !!
	      
	  // MJD
	  mjdMS_v.push_back(ArrayTime(tabulation_time_v[itab]).getMJD());
	    
	  // RA / DEC
	  raMS_v.push_back(evalPoly(ra_coeff_v.size(), ra_coeff_v, timeOrigin, tabulation_time));

	  decMS_v.push_back(evalPoly(dec_coeff_v.size(), dec_coeff_v, timeOrigin, tabulation_time));
	    
	  // DISTANCE
	  distanceMS_v.push_back(evalPoly(distance_coeff_v.size(), distance_coeff_v, timeOrigin, tabulation_time));
	    
	  // RADVEL
	  if (radVelExists) {
	    radVelMS_v.push_back(evalPoly(radvel_coeff_v.size(),
					  radvel_coeff_v,
					  timeOrigin,
					  tabulation_time));
	  }
	}
      } else {
	// both of these cases require that the numPoly* is either all 1 or always greater than 1. Check here.
	errstream.str("");
	for (unsigned int i = 1; i < ephRow_v.size(); i++) {
	  if (numPolyDirIsOne != (ephRow_v[i]->getNumPolyDir() == 1)) {
	    errstream << "In the table Ephemeris the value of the field 'numPolyDir' is expected to be always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	    error(errstream.str());
	  }
	  
	  if (numPolyDistIsOne != (ephRow_v[i]->getNumPolyDist() == 1)) {
	    errstream << "In the table Ephemeris the value of the field 'numPolyDist' is expected to be always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	    error(errstream.str());
	  }
	  
	  if (radVelExists != (ephRow_v[i]->isRadVelExists() && ephRow_v[i]->isNumPolyRadVelExists())) {
	    errstream << "In the table Ephemeris the fields 'radVel' and 'numPolyRadVel' are expected to be always absent or always present. This rule is violated at line #" << i <<".";
	    error(errstream.str());
	  }
	  
	  if (radVelExists) {
	    if (numPolyRadVelIsOne != (ephRow_v[i]->getNumPolyRadVel() == 1)) {
	      errstream << "In the table Ephemeris the value of the field 'numPolyRadVel' is expected to be always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	      error(errstream.str());
	    }	 
	  }
	}
	if (!interpolate_ephemeris && allNumPolyIsOne) {
	  // interpolation is NOT requested and all possible polynomial columns are simple scalars, numPoly==1
	  // Just copy ephemeris without any interpolation. Just adapt the units.
	  infostream.str("");
	  infostream << "The MS Ephemeris table for ephemerisId = '" << ephemerisId
		     << "' will be produced by copying the values found in the ASDM with no interpolation";
	  info(infostream.str());
	  for(const EphemerisRow *eR_p: ephRow_v) {
	    mjdMS_v.push_back(eR_p->getTimeInterval().getMidPoint().getMJD()); // MJD
	    vector<vector<double> > dir = eR_p->getDir();
	    raMS_v.push_back(radian2degree(dir[0][0]));  // deg
	    decMS_v.push_back(radian2degree(dir[0][1])); // deg
	    distanceMS_v.push_back(m2au(eR_p->getDistance()[0])); // AU
	    if (radVelExists) {
	      radVelMS_v.push_back(mpers2auperd(eR_p->getRadVel()[0])); // AU/d
	    }
	  }
	} else {
	  // the general case, polynomials are evaluated and scalars are interpolated to the same timesteps.
	  infostream.str("");
	  infostream << "The MS Ephemeris table for ephemerisId = '" << ephemerisId
		     << "' will be produced by tabulating the polynomials found in the 'dir', 'distance' and optionally 'radVel' columns with a timestep of '"
		     << timeStepInNanoSecond / 1.e9 / 86400. << "' days"; 
	  info(infostream.str());
	  if (anyNumPolyIsOne) {
	    infostream.str("");
	    infostream << "the polynomials of order 0 will be interpolated between mid-interval values using the same timestep.";
	    info(infostream.str());
	  }

	  // the times that will be used depends on whether there are any numPoly*IsOne cases.  In that case,
	  // times must run from the first interval mid-point through the last interval mid-point.  Otherwise
	  // it's all evaluating the existing polynomials and the times can run from the first start time to the end
	  // of the last interval.

	  // Check that for each polynomial column the degree is always null or never null and that the optional fields are always present or always absent.
	  // And also verify that there is no "hole" in the time range covered by the sequence of ArrayTime intervals when the degree is == 0.

	  // this vector is only used in the linear interpolation case
	  vector<double> time_v;      //  "      "
	  
	  if (anyNumPolyIsOne) {
	    // need to recalculate t0MS using the midPoint
	    t0ASDM = ephRow_v[0]->getTimeInterval().getMidPoint().get();
	    q = t0ASDM / timeStepInNanoSecond;
	    r = t0ASDM % timeStepInNanoSecond;
	    t0MS = t0ASDM;
	    if (r != 0) {
	      q = q+1;
	      t0MS = q * timeStepInNanoSecond;
	    }

	    // and populate time_v and look for gaps
	    time_v.push_back(1.0e-09*t0ASDM);  

	    for (unsigned int i = 1; i < ephRow_v.size(); i++) {
	      int64_t midPoint_i = ephRow_v[i]->getTimeInterval().getMidPoint().get() ;
	      int64_t midPoint_i_1 = ephRow_v[i-1]->getTimeInterval().getMidPoint().get();
	      int64_t duration_i_1 = ephRow_v[i-1]->getTimeInterval().getDuration().get();
	      // look for gaps
	      // comparing midPoint times should be equivalent to start times for this purpose, and we need the mid points for time_v
	      if (midPoint_i != (midPoint_i_1 + duration_i_1)) {
		infostream.str("");
		infostream << "The value of 'timeInterval' at row #" << i-1 << " does not cover the time range up to the start time of the next row. The polynomial will be evaluated despite the presence of this 'hole'";
		info(infostream.str());
	      }
	      time_v.push_back(1.0e-09*midPoint_i); 
	    }
	  }
    
	  // Determine the timely ordered sequence of indexes in ephRow_v which will be used to tabulate the ephemeris data to be put into the MS table.
	  LOG("Prepare the time ordered sequence of indexes used to tabulate the ephemeris data to be written in the MS table.");

	  // The polynomials change at the interval boundaries, the linear interpolations change at the interval mid-points
	  // It may be necessary to have both sets of (index,time) pairs.

	  typedef pair<uint32_t, int64_t> atiIdxMStime_pair;
	  vector<atiIdxMStime_pair>  atiIdxMStimePoly_v;
	  vector<atiIdxMStime_pair>  atiIdxMStimeLine_v;

	  // BUT - the time values used in each pair MUST be the same, starting from t0MS through either
	  // the end of the final interval OR through the midpoint of the final interval
	  // all times here the integer times in nano seconds
	  vector<int64_t> intMStimes_v;
	  int64_t tMS = t0MS;
	  int32_t lastIntIndx = ephRow_v.size()-1;
	  int64_t end = ephRow_v[lastIntIndx]->getTimeInterval().getStart().get() + ephRow_v[lastIntIndx]->getTimeInterval().getDuration().get();
	  // there must be a vector way to do this - revisit this part later, and even a for loop would be better
	  do {
	    intMStimes_v.push_back(tMS);
	    tMS += timeStepInNanoSecond;
	  } while (tMS < end);
	  
	  // populate the pairs as needed
	  if (!allNumPolyIsOne) {
	    // some polynomials exist
	    uint32_t index = 0;  
	    int64_t  start =  ephRow_v[index]->getTimeInterval().getStart().get();
	    int64_t  end   =  start + ephRow_v[index]->getTimeInterval().getDuration().get();

	    // position index so that t0MS is < end
	    while ((t0MS >= end) && (index < ephRow_v.size())) {
	      index++;
	      end = ephRow_v[index]->getTimeInterval().getStart().get() + ephRow_v[index]->getTimeInterval().getDuration().get();
	    }

	    if (index >= ephRow_v.size()) {
	      // this might happen if the user chose a time step that's really bad for the ephemeris table being interpolated
	      // can not continue
	      errstream.str("");
	      errstream << "The ephemeris table can't be evaluated or interpolated using the chosen timestep : "
			<< timeStepInNanoSecond / 1.e9 / 86400. << "' days.  Choose a smaller value.";
	      error(errstream.str());
	    }

	    for (int64_t tMS : intMStimes_v) {
	      // make sure this index is appropriate for tMS
	      // there may be a gap so it's important to keep advancing index until there's nothing left
	      while ((tMS >= end) && (index < ephRow_v.size())) {
		index++;
		end = ephRow_v[index]->getTimeInterval().getStart().get() + ephRow_v[index]->getTimeInterval().getDuration().get();
	      }
	      if (index < ephRow_v.size()) {
		atiIdxMStimePoly_v.push_back(atiIdxMStime_pair(index, tMS));
		LOG ("size of atiIdxMStimePoly_v="+TO_STRING(atiIdxMStimePoly_v.size())+", index = "+TO_STRING(index)+", tMS = "+TO_STRING(tMS));
	      }
	    }
	  }
	  if (anyNumPolyIsOne) {
	    // some linear interpolation is needed - breaks at the interval midpoint
	    uint32_t index = 0;  
	    int64_t  start =  ephRow_v[index]->getTimeInterval().getMidPoint().get();
	    int64_t  end   =  start + ephRow_v[index]->getTimeInterval().getDuration().get();

	    // the linear interpolation can't ever start with the final index number, that would required interpolating beyond the end
	    // position index so that t0MS is < end
	    while ((t0MS >= end) && (index < (ephRow_v.size()-1))) {
	      index++;
	      end = ephRow_v[index]->getTimeInterval().getMidPoint().get() + ephRow_v[index]->getTimeInterval().getDuration().get();
	    }

	    if (index >= (ephRow_v.size()-1)) {
	      // this might happen if the user chose a time step that's really bad for the ephemeris table being interpolated
	      // can not continue
	      errstream.str("");
	      errstream << "The ephemeris table can't be evaluated or interpolated using the chosen timestep : "
			<< timeStepInNanoSecond / 1.e9 / 86400. << "' days.  Choose a smaller value.";
	      error(errstream.str());
	    }

	    for (int64_t tMS : intMStimes_v) {
	      // make sure this index is appropriate for tMS
	      // there may be a gap, so keep advancing until there's nothing left
	      while (tMS >= end && (index < (ephRow_v.size()-1))) {
		index++;
		end = ephRow_v[index]->getTimeInterval().getMidPoint().get() + ephRow_v[index]->getTimeInterval().getDuration().get();
	      }
	      if (index < (ephRow_v.size()-1)) {
		atiIdxMStimeLine_v.push_back(atiIdxMStime_pair(index, tMS));
		LOG ("size of atiIdxMStimeLine_v="+TO_STRING(atiIdxMStimeLine_v.size())+", index = "+TO_STRING(index)+", tMS = "+TO_STRING(tMS));
	      }
	    }
	  }
    
	  LOG("atiIdxMStimePoly_v has " + TO_STRING(atiIdxMStimePoly_v.size()) + " elements.");
	  LOG("atiIdxMStimeLine_v has " + TO_STRING(atiIdxMStimeLine_v.size()) + " elements.");

	  // Prepare the coefficients which will be used for the tabulation.
 
	  LOG("Prepare the coefficients which will be used for the tabulations.");
	  vector<vector<double> >  raASDM_vv;
	  vector<double>           raASDM_v;
	  vector<vector<double> >  decASDM_vv;
	  vector<double>           decASDM_v;
	  vector<vector<double> >  distanceASDM_vv;
	  vector<double>           distanceASDM_v;
	  vector<vector<double> >  radVelASDM_vv;
	  vector<double>           radVelASDM_v;
	  vector<double>           empty_v;
	  vector<double>           temp_v;

	  cout.precision(10);
	  for (unsigned int i = 0; i < ephRow_v.size(); i++) {
	    LOG_EPHEM("original " + TO_STRING (ArrayTime(ephRow_v[i]->getTimeInterval().getStart().get()).getMJD()));
	    vector<vector<double> > temp_vv = ephRow_v[i]->getDir();
	    if (numPolyDirIsOne) {
	      raASDM_v.push_back(radian2degree(temp_vv[0][0]));
	      decASDM_v.push_back(radian2degree(temp_vv[0][1]));
	      LOG_EPHEM (" " + TO_STRING(raASDM_v.back()) + " " + TO_STRING(decASDM_v.back()));
	    } else {
	      raASDM_vv.push_back(empty_v);
	      decASDM_vv.push_back(empty_v);
	      for (int j = 0; j < ephRow_v[i]->getNumPolyDir(); j++) {
		raASDM_vv.back().push_back(radian2degree(temp_vv[j][0]));
		decASDM_vv.back().push_back(radian2degree(temp_vv[j][1]));          ;
	      }
	    }
	    
	    temp_v = ephRow_v[i]->getDistance();      
	    if (numPolyDistIsOne) {
	      distanceASDM_v.push_back(m2au(temp_v[0]));           // AU
	      LOG_EPHEM (" " + TO_STRING(distanceASDM_v.back()));
	    } else {
	      distanceASDM_vv.push_back(empty_v);
	      for (int j = 0; j < ephRow_v[i]->getNumPolyDist(); j++)
		distanceASDM_vv.back().push_back(m2au(temp_v[j])); // AU
	    }
	    
	    if (radVelExists) {
	      temp_v = ephRow_v[i]->getRadVel();
	      if (numPolyRadVelIsOne) { 
		radVelASDM_v.push_back(mpers2auperd(temp_v[0]));      // AU/d
		LOG_EPHEM(" " + TO_STRING(radVelASDM_v.back()));
	      } else {
		radVelASDM_vv.push_back(empty_v);
		for (int j = 0; j < ephRow_v[i]->getNumPolyRadVel(); j++)
		  radVelASDM_vv.back().push_back(mpers2auperd(temp_v[j]));   // AU/d
	      }	
	    }
	    LOG_EPHEM("\n");
	  }
	  
	  // Preparing the coefficients of piecewise polynomial of degree 1.
	  if (numPolyDirIsOne) {
	    LOG("Compute the linear interpolation coefficients for RAD");
	    linearInterpCoeff(ephRow_v.size(), time_v, raASDM_v, raASDM_vv);
	      
	    LOG("Compute the linear interpolation coefficients for DEC");
	    linearInterpCoeff(ephRow_v.size(), time_v, decASDM_v, decASDM_vv);
	  }
	    
	  if (numPolyDistIsOne)  {
	    LOG("Compute the linear interolation coefficients for Dist");
	    linearInterpCoeff(ephRow_v.size(), time_v, distanceASDM_v, distanceASDM_vv);
	  }
	    
	  if (radVelExists && numPolyRadVelIsOne) {
	    LOG("Compute the linear interpolation coefficients for RadVel");
	    linearInterpCoeff(ephRow_v.size(), time_v, radVelASDM_v, radVelASDM_vv);
	  }     

	  // loops over the two vectors of pairs
	  // the polynomials
	  for (atiIdxMStime_pair atiIdxMStime: atiIdxMStimePoly_v) {
	    //
	    // MJD
	    // mjdMS_v is always unfilled at this point, this is always appropriate here
	    mjdMS_v.push_back(ArrayTime(atiIdxMStime.second).getMJD());
	    LOG_EPHEM(TO_STRING(mjdMS_v.size()) + ": ");
	    LOG_EPHEM( "resampled " + TO_STRING(mjdMS_v.back()));
	    LOG("mjdMS_v -> "+TO_STRING(mjdMS_v.back()));
	    
	    double timeOrigin = 1.0e-09 * ephRow_v[atiIdxMStime.first]->getTimeOrigin().get();
	    double time       = 1.0e-09 * atiIdxMStime.second;
	    
	    LOG("timeOrigin="+TO_STRING(timeOrigin)+", time="+TO_STRING(time));
	    
	    // RA / DEC
	    if (!numPolyDirIsOne) {
	      LOG("Eval poly for RA");
	      LOG("atiIdxMStime.first = " + TO_STRING(atiIdxMStime.first));
	      raMS_v.push_back(evalPoly(raASDM_vv[atiIdxMStime.first].size(),
					raASDM_vv[atiIdxMStime.first], timeOrigin, time));
	      LOG_EPHEM(" " + TO_STRING(raMS_v.back()));
	      LOG("raMS_v -> "+TO_STRING(raMS_v.back()));
	    
	      LOG("Eval poly for DEC");
	      decMS_v.push_back(evalPoly(decASDM_vv[atiIdxMStime.first].size(),
					 decASDM_vv[atiIdxMStime.first], timeOrigin, time));
	      LOG_EPHEM(" " + TO_STRING(decMS_v.back()));
	    }
	    
	    // Distance
	    if (!numPolyDistIsOne) {
	      LOG("Eval poly for distance");
	      distanceMS_v.push_back(evalPoly(distanceASDM_vv[atiIdxMStime.first].size(),
					      distanceASDM_vv[atiIdxMStime.first], timeOrigin, time));
	      LOG_EPHEM(" " + TO_STRING(distanceMS_v.back()));

	      // Radvel
	      if (radVelExists && !numPolyRadVelIsOne) { 
		LOG("Eval poly for radvel");
		radVelMS_v.push_back(evalPoly(radVelASDM_vv[atiIdxMStime.first].size(),
					      radVelASDM_vv[atiIdxMStime.first], timeOrigin, time));
		LOG_EPHEM(" " + TO_STRING(radVelMS_v.back()));
	      }
	      LOG_EPHEM("\n");
	    }
	  } // end of polynomial evaluations
	  
	  // the linear interpolations
	  for (atiIdxMStime_pair atiIdxMStime: atiIdxMStimeLine_v) {
	    //
	    // MJD
	    // mjdMS_v has already been filled if the poly pair vector has a non-zero size
	    if (atiIdxMStimePoly_v.size() == 0) {
	      mjdMS_v.push_back(ArrayTime(atiIdxMStime.second).getMJD());
	      LOG_EPHEM(TO_STRING(mjdMS_v.size()) + ": ");
	      LOG_EPHEM( "resampled " + TO_STRING(mjdMS_v.back()));
	      LOG("mjdMS_v -> "+TO_STRING(mjdMS_v.back()));
	    }

	    // linear interpolations are evaluated relative to the mid point
	    double timeMidPoint = 1.0e-09 * ephRow_v[atiIdxMStime.first]->getTimeInterval().getMidPoint().get();
	    double time       = 1.0e-09 * atiIdxMStime.second;
	    
	    LOG("timeMidPoint="+TO_STRING(timeMidPoint)+", time="+TO_STRING(time));
	    
	    // RA / DEC
	    if (numPolyDirIsOne) {
	      LOG("Eval line for RA");
	      LOG("atiIdxMStime.first = " + TO_STRING(atiIdxMStime.first));
	      raMS_v.push_back(evalPoly(raASDM_vv[atiIdxMStime.first].size(),
					raASDM_vv[atiIdxMStime.first], timeMidPoint, time));
	      LOG_EPHEM(" " + TO_STRING(raMS_v.back()));
	      LOG("raMS_v -> "+TO_STRING(raMS_v.back()));
	    
	      LOG("Eval line for DEC");
	      decMS_v.push_back(evalPoly(decASDM_vv[atiIdxMStime.first].size(),
					 decASDM_vv[atiIdxMStime.first], timeMidPoint, time));
	      LOG_EPHEM(" " + TO_STRING(decMS_v.back()));
	    }
	    
	    // Distance
	    if (numPolyDistIsOne) {
	      LOG("Eval line for distance");
	      distanceMS_v.push_back(evalPoly(distanceASDM_vv[atiIdxMStime.first].size(),
					      distanceASDM_vv[atiIdxMStime.first], timeMidPoint, time));
	      LOG_EPHEM(" " + TO_STRING(distanceMS_v.back()));

	      // Radvel
	      if (radVelExists && numPolyRadVelIsOne) { 
		LOG("Eval line for radvel");
		radVelMS_v.push_back(evalPoly(radVelASDM_vv[atiIdxMStime.first].size(),
					      radVelASDM_vv[atiIdxMStime.first], timeMidPoint, time));
		LOG_EPHEM(" " + TO_STRING(radVelMS_v.back()));
	      }
	      LOG_EPHEM("\n");
	    } // end of linear evaluations
	  }
	} // end of the general interpolation case
      } // end of the non-single row case, everything is ready to push to the MS now

      // 
      // Record the starting time + one time step so that Field can use it if needed.
      //
      ephemStartTime_m[ephemerisId] = ArrayTime(mjdMS_v[1]).get()*1.e-09;

      // Now the data are ready to be written to the MS Ephemeris table.
      // Let's proceed, using Slicers.
    
      unsigned int numRows = raMS_v.size();
      Slicer slicer(IPosition(1, 0), IPosition(1, numRows-1), Slicer::endIsLast);

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	Table * table_p = apc2EphemTable_m[iter->first];

	// Update the MJD0 Keyword to the correct value, i.e. mjdMS_v[0] - dmjd
	table_p->rwKeywordSet().define("MJD0", casacore::Double(mjdMS_v[0] - dmjd));
	
	// Update the dMJD keyword.
	table_p->rwKeywordSet().define("dMJD", casacore::Double(dmjd));

	// And fill the table
	table_p->addRow(numRows);
	LOG ("Added "+TO_STRING(numRows)+" rows to table "+((string)table_p->tableName()));
      
	LOG("Filling column MJD");
	Vector<casacore::Double> MJD_V(IPosition(1, numRows), &mjdMS_v[0], SHARE);
	ScalarColumn<casacore::Double> MJD(*table_p, "MJD");
	MJD.putColumnRange(slicer, MJD_V);
      
	LOG("Filling column RA");
	Vector<casacore::Double> RA_V(IPosition(1, numRows), &raMS_v[0], SHARE);
	ScalarColumn<casacore::Double> RA(*table_p,  "RA");
	RA.putColumnRange(slicer, RA_V);
      
	LOG("Filling column DEC");
	Vector<casacore::Double> DEC_V(IPosition(1, numRows), &decMS_v[0], SHARE);
	ScalarColumn<casacore::Double> DEC(*table_p, "DEC");
	DEC.putColumnRange(slicer, DEC_V);
      
	LOG ("Filling column Rho");
	Vector<casacore::Double> Rho_V(IPosition(1, numRows), &distanceMS_v[0], SHARE);
	ScalarColumn<casacore::Double> Rho(*table_p, "Rho");
	Rho.putColumnRange(slicer, Rho_V);
      
	if (radVelExists) {
	  LOG ("Filling column RadVel");
	  Vector<casacore::Double> RadVel_V(IPosition(1, numRows), &radVelMS_v[0], SHARE);
	  ScalarColumn<casacore::Double> RadVel(*table_p, "RadVel");
	  RadVel.putColumnRange(slicer, RadVel_V);
	}
      
	infostream.str("");
	infostream << "converted in " << table_p->nrow() << " ephemeris rows in the table '" << table_p->tableName() << "'.";
	info(infostream.str());
      }
      deleteEphemeris(apc2EphemTable_m); // Get rid of the ephemeris tables now that we have finished with them.
    }
  } catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillEphemeris");
  
}

/** 
 * This function fills the MS Field table.
 * given :
 * @parameter ds_p a pointer to the ASDM dataset.
 * @parameter considerEphemeris take into account the reference to Ephemeris table(s).
 */
void fillField(ASDM* ds_p, bool considerEphemeris) {
  LOGENTER("fillField");
  vector<pair<int, int> > idxEphemerisId_v;

  try {
    FieldTable& fieldT = ds_p->getField();
    FieldRow* r = 0;
    int nField = fieldT.size();
    infostream.str("");
    infostream << "The dataset has " << nField << " field(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nField; i++) {
      if ((r=fieldT.getRowByKey(Tag(i, TagType::Field))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Field table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      /*
       * The field name should not contain any "&" because this character has a specific meaning in the MS SELECTION syntax
       * Any occurrence of character "&" in a field name is replaced by a character "#" and a warning is emitted in the log
       * about this replacement.
       */
      string fieldName = r->getFieldName();
      if (fieldName.front()=='&') {
	replace(fieldName.begin(), fieldName.end(), '&', '#');
	infostream.str("");
	infostream << "ATTENTION !!! In row #" << i << " of the Field table, the character '&' has been replaced by the character '#' in the field name." << endl;
	info(infostream.str());
      }

      string code = r->getCode();
      DirectionReferenceCodeMod::DirectionReferenceCode dirRefCode = DirectionReferenceCodeMod::J2000;
      if(r->isDirectionCodeExists()){
	dirRefCode = r->getDirectionCode();
	//cout << "found directionCode for field " << fieldName << ": ";
      }
      //       else{
      // 	cout << "No directionCode in input table. Assuming ";
      //       }
      string directionCode = CDirectionReferenceCode::name(dirRefCode);
      //cout << directionCode << endl;
      
      vector<vector<double> > delayDir     = DConverter::toMatrixD<Angle>(r->getDelayDir());
      vector<vector<double> > phaseDir     = DConverter::toMatrixD<Angle>(r->getPhaseDir());
      vector<vector<double> > referenceDir = DConverter::toMatrixD<Angle>(r->getReferenceDir());
      
      if (r->isEphemerisIdExists()) idxEphemerisId_v.push_back(pair<int, int>(i, r->getEphemerisId()));

      /*
       * The value of time requires some attention.
       * if :
       * 1) there is an ephemeris associated to the field and
       * 2) time is not present in the ASDM Field table or if time is present but equal to 0 
       * then:
       *  write a time in the MS Field table equal to the start of the first arrayTime interval of validity in the Ephemeris.
       *
       * else :
       * write the time found in the ASDM Field relevant columns.
       *
       */
      vector<EphemerisRow *> eR_v;
      
      if (r->isEphemerisIdExists())
	eR_v = *(r->getTable().getContainer().getEphemeris().getByContext(r->getEphemerisId()));
	
      bool getTimeFromEphemeris = (eR_v.size() > 0) && (!r->isTimeExists() || (r->isTimeExists() && (r->getTime().get() == 0)));
      
      double time = 0.0;
      
      if (getTimeFromEphemeris)
	time = ephemStartTime_m[r->getEphemerisId()]; // It's already in second.
      else if (r->isTimeExists())
	time = ((double) r->getTime().get()) * 1.e-09 ;
      else
	time = 0.0;

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addField( fieldName, code, time, r->getNumPoly(), delayDir, phaseDir, referenceDir,
				directionCode, r->isSourceIdExists()?r->getSourceId():0 );
      }
    }
    
    if (considerEphemeris && (idxEphemerisId_v.size() > 0)) 
        for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->updateEphemerisIdInField(idxEphemerisId_v);
      }

    if (nField) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->field().nrow() << "  field(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  LOGEXIT("fillField");
}

/**
 * Utility function to sort out the numBin value given appropriate for a given row in the SpectralWindow table
 */
int getNumBin(SpectralWindowRow *spwRow, const string &telescopeName) {
  // assumes that spRow is already a non-null pointer
  if (spwRow->isNumBinExists()) {
    // use any numBin value that already exists
    return spwRow->getNumBin();
  } else {
    // 1 will be returned if numBin <=0 at the end
    int numBin = 0;

    // infer it based on the telescope name
    // both methods require that these 2 scalar values exist
    if (spwRow->isChanWidthExists() && spwRow->isResolutionExists()) {
      if (telescopeName == string("EVLA")) {
	numBin = std::nearbyint(spwRow->getChanWidth().get()/spwRow->getResolution().get());
      } else if (telescopeName == string("ALMA")) {
	// only relevant for HANNING
	if (spwRow->getWindowFunction()==WindowFunctionMod::HANNING) {
	  // also requires that effectiveBw exists
	  if (spwRow->isEffectiveBwExists()) {
	    double chanWidth = spwRow->getChanWidth().get();
	    if (spwRow->getResolution().get() != chanWidth) {
	      // effectiveBw and resolution have been set according to the Hills memo, get implied numBin
	      // the expected values of the ratio of effectiveBw/chanWidth are from the Hills memo
	      // second table on p. 6
	      double ratio = spwRow->getEffectiveBw().get()/std::abs(chanWidth);
	      // no other way to do this than one value at a time, starting from numBin=1
	      // the Hills memo values are before decimation, ratio above is after decimation

	      // tolerance here is 5.0e-4, i.e. precision of values in Hills memo
	      double tolerance=5.0e-4;
	      if (std::abs(2.667-ratio)<=tolerance) {
		numBin = 1;
	      } else if (std::abs(3.200-ratio*2.0)<=tolerance) {
		numBin = 2;
	      } else if (std::abs(4.923-ratio*4.0)<=tolerance) {
		numBin = 4;
	      } else if (std::abs(8.828-ratio*8.0)<=tolerance) {
		numBin = 8;
	      } else if (std::abs(16.787-ratio*16.0)<tolerance) {
		numBin=16;
	      }
	    }
	  }
	}
      }
    }
    if (numBin <= 0) {
      // was not set to a valid value
      // none of these cases generate an error or warning
      //   scalar columns don't exist - that's expected and so not an error
      //   unrecognzed telescope - that's also probably not an error, no rule for inferring numBin
      //   bad value for the VLA from that ratio, that might be an reportable error
      //   other windowFunction for ALMA - that's expected and so no an error, but might report if that's not UNIFORM
      //   no match within tolerance for ALMA - that might be a reportable error
      numBin = 1;
    }
    return numBin;
  }
  // there's no way to get to here
  throw ASDM2MSException("Impossible code line reached in getNumBin");
}

/** 
 * This function fills the MS Spectral Window table.
 * given :
 * @parameter ds_p a pointer to the ASDM dataset.
 */
void fillSpectralWindow(ASDM* ds_p, map<unsigned int, double>& effectiveBwPerSpwId_m, const string &telescopeName) {
  LOGENTER("fillSpectralWindow");

  effectiveBwPerSpwId_m.clear();

  try {
    SpectralWindowTable& spwT = ds_p->getSpectralWindow();      
    vector<Tag> reorderedSwIds = reorderSwIds(*ds_p); // The vector of Spectral Window Tags in the order they will be inserted in the MS.
    //for (vector<Tag>::size_type i = 0; i != reorderedSwIds.size() ; i++) cerr<<" reorderedSwIds["<<i<<"]="<<reorderedSwIds[i].getTagValue()<<endl;
 
    for (vector<Tag>::size_type i = 0; i != reorderedSwIds.size() ; i++) swIdx2Idx[reorderedSwIds[i].getTagValue()] = i;

    SpectralWindowRow* r = 0;
    int nSpectralWindow = spwT.size();
    
    infostream.str("");
    infostream << "The dataset has " << nSpectralWindow << " spectral window(s)..."; 
    info(infostream.str());
    
    for (vector<Tag>::size_type i = 0; i < reorderedSwIds.size(); i++) {
      if ((r = spwT.getRowByKey(reorderedSwIds[i])) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the SpectralWindow table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
            
      /* Processing the chanFreqXXX
       *
       * Either (chanFreqStart, chanFreqStep) or (chanFreqArray) must be present
       * with a priority to chanFreqArray if both are present. If chanFreqArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      bool chanStartAndStep = r->isChanFreqStartExists() && r->isChanFreqStepExists();
      bool chanArray = r->isChanFreqArrayExists();
      if (!chanStartAndStep && !chanArray) {
	errstream.str("");
	errstream << "Did not find (chanFreqStart, chanFreqStep) nor chanFreqArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* chanFreq1D = (double *) 0;
      //DConverter chanFreqConverter;
      vector<double> chanFreq1D;
      vector<Frequency> chanFreqArray;
      Frequency chanFreqStart, chanFreqStep;
      if (chanArray) { // Frequency channels are specified by an array.
	chanFreqArray = r->getChanFreqArray();
	if (chanFreqArray.size() != (unsigned int)r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of chanFreqArray ('"
		    << chanFreqArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      } else { // Frequency channels are specified by a (start, step) pair.
	chanFreqStart = r->getChanFreqStart();
	chanFreqStep  = r->getChanFreqStep();
	for (int i = 0; i < r->getNumChan(); i++)
	  chanFreqArray.push_back(chanFreqStart + i * chanFreqStep);
      }
      //chanFreq1D = chanFreqConverter.to1DArray(chanFreqArray);
      chanFreq1D = DConverter::toVectorD<Frequency>(chanFreqArray);
      
      /* Processing the chanWidthXXX
       *
       * Either chanWidth or chanWidthArray must be present
       * with a priority to chanWidthArray if both are present. If chanWidthArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isChanWidthExists() && !r->isChanWidthArrayExists()) {
	errstream.str("");
	errstream << "Did not find chanWidth nor chanWidthArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* chanWidth1D = (double *) 0;
      //DConverter chanWidthConverter;
      vector<double> chanWidth1D;
      vector<Frequency> chanWidthArray;
      if (r->isChanWidthArrayExists()) { // Frequency channels widths are specified by an array.
	chanWidthArray = r->getChanWidthArray();
	if (chanWidthArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of chanWidthArray ('"
		    << chanWidthArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      } else { // Frequency channels widths are specified by a constant value.
	chanWidthArray.resize(r->getNumChan());
	chanWidthArray.assign(chanWidthArray.size(), r->getChanWidth());
      }
      //chanWidth1D = chanWidthConverter.to1DArray(chanWidthArray);
      chanWidth1D = DConverter::toVectorD<Frequency>(chanWidthArray);

      // To answer JIRA ticket CAS-3265 / Dirk Petry
      if (chanFreq1D[chanFreq1D.size() - 1] < chanFreq1D[0] )
	transform(chanWidth1D.begin(), chanWidth1D.end(), chanWidth1D.begin(), negateFunctor<double>());
      
      /* Processing the effectiveBwXXX
       *
       * Either effectiveBw or effectiveBwArray must be present
       * with a priority to effectiveBwArray if both are present. If effectiveBwArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isEffectiveBwExists() && !r->isEffectiveBwArrayExists()) {
	errstream.str("");
	errstream << "Did not find effectiveBw nor effectiveBwArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* effectiveBw1D = (double *) 0;
      vector<double> effectiveBw1D;
      //DConverter effectiveBwConverter;
      vector<Frequency> effectiveBwArray;
      if (r->isEffectiveBwArrayExists()) { // Effective BWs are specified by an array.
	effectiveBwArray = r->getEffectiveBwArray();
	if (effectiveBwArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of effectiveBwArray ('"
		    << effectiveBwArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further." ;
	  error(errstream.str());
	}
      } else { // Effective BWs are specified by a constant value.
	effectiveBwArray.resize(r->getNumChan());
	effectiveBwArray.assign(effectiveBwArray.size(), r->getEffectiveBw());
      }
      //effectiveBw1D = effectiveBwConverter.to1DArray(effectiveBwArray);
      effectiveBw1D = DConverter::toVectorD<Frequency>(effectiveBwArray);

      effectiveBwPerSpwId_m[r->getSpectralWindowId().getTagValue()] = effectiveBw1D.at(0);
      
      
      /* Processing the resolutionXXX
       *
       * Either resolution or resolutionArray must be present
       * with a priority to resolutionArray if both are present. If resolutionArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isResolutionExists() && !r->isResolutionArrayExists()) {
	errstream.str("");
	errstream << "Did not find resolution nor resolutionArray. Can't go further";
	error(errstream.str());
      }
      
      //double* resolution1D = (double *) 0;
      vector<double> resolution1D;
      //DConverter resolutionConverter;
      vector<Frequency> resolutionArray;
      if (r->isResolutionArrayExists()) { // Resolutions are specified by an array.
	resolutionArray = r->getResolutionArray();
	if (resolutionArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of resolutionArray ('"
		    << resolutionArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      } else { // Resolutions are specified by a constant value.
	resolutionArray.resize(r->getNumChan());
	resolutionArray.assign(resolutionArray.size(), r->getResolution());
      }
      //resolution1D = resolutionConverter.to1DArray(resolutionArray);
      resolution1D = DConverter::toVectorD<Frequency>(resolutionArray);

      /*
       * associated spectral windows and and natures.
       */
      
      unsigned int numAssocValues = 0;
      if (r->isNumAssocValuesExists()) numAssocValues = r->getNumAssocValues();
      
      // Test the simultaneous presence or absence of assocNature and assoSpectralWindowId
      if ((r->isAssocNatureExists() && !r->isAssocSpectralWindowIdExists()) ||
	  (!r->isAssocNatureExists() && r->isAssocSpectralWindowIdExists())) {
	errstream.str("");
	errstream << "Only one of the attributes assocSpectralWindowId and assocNature is present. Can't go further."
		  << endl;
	error(errstream.str());
      }
      
      vector<int> assocSpectralWindowId_;
      vector<string> assocNature_ ;

      if (r->isAssocSpectralWindowIdExists()) { // it exists then the assocNature exists also, given the test
	                                        // which is done before.
	vector<Tag> assocSpectralWindowId = r->getAssocSpectralWindowId();
	if (numAssocValues != assocSpectralWindowId.size()) {
	  infostream.str("");
	  infostream << "The size of assocSpectralWindowId ('"
		     << assocSpectralWindowId.size()
		     << "') is not equal to the value announced in numAssocValues ('"
		     << numAssocValues
		     << "'). Ignoring the difference and sending the full vector assocSpectralWindowId to the filler";
	  info(infostream.str());
	}
	numAssocValues = assocSpectralWindowId.size();
	assocSpectralWindowId_ = IConverter::toVectorI(assocSpectralWindowId);

	// Take into account the re ordering of the spectral window indices.
	for (unsigned int iAssocSw = 0; iAssocSw < numAssocValues; iAssocSw++)
	  assocSpectralWindowId_[iAssocSw] =  swIdx2Idx[assocSpectralWindowId_[iAssocSw]];

	vector<SpectralResolutionTypeMod::SpectralResolutionType> assocNature = r->getAssocNature();

	if (assocNature.size() != assocSpectralWindowId_.size()) {
	  infostream.str("");
	  infostream << "The size of assocNature ('"
		     << assocNature.size() 
		     << "') is not equal to the size of assocSpectralWindowId ('"
		     << assocSpectralWindowId.size()
		     << "'). Ignoring the difference and sending the full assocNature vector to the filler.";
	  info(infostream.str());
	}
	assocNature_ = SConverter::toVectorS<SpectralResolutionTypeMod::SpectralResolutionType, CSpectralResolutionType>(r->getAssocNature());
      }
      
      int numChan           = r->getNumChan();
      string name           = r->isNameExists()?r->getName():"";
      double refFreq        = r->getRefFreq().get();
      int measFreqRef       = r->isMeasFreqRefExists()?FrequencyReferenceMapper::value(r->getMeasFreqRef()):MFrequency::TOPO;
      double totalBandwidth = r->getTotBandwidth().get();
      int netSideband       = r->getNetSideband();
      int bbcNo             = r->getBasebandName();
      int ifConvChain       = 0;
      int freqGroup         = r->isFreqGroupExists()?r->getFreqGroup():0;
      string freqGroupName  = r->isFreqGroupNameExists()?r->getFreqGroupName().c_str():"";
      // windowFunction is a required field
      std::string windowFunction = CWindowFunction::name(r->getWindowFunction());
      int numBin = getNumBin(r, telescopeName);
      if (telescopeName == "EVLA" && (numBin>1) && (!r->isNumBinExists())) {
	// numBin has been inferred for EVLA data, adjust resolution 
	resolution1D = chanWidth1D;
      }

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addSpectralWindow(numChan,
					name,
					refFreq,
					chanFreq1D,  
					chanWidth1D, 
					measFreqRef,
					effectiveBw1D, 
					resolution1D,
					totalBandwidth,
					netSideband,
					bbcNo,
					ifConvChain,
					freqGroup,
					freqGroupName,
					numAssocValues,
					assocSpectralWindowId_,
					assocNature_,
					windowFunction,
					numBin);
      }      
    }
    if (nSpectralWindow) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->spectralWindow().nrow() << " spectral window(s) in the measurement set(s).";
      info(infostream.str());
    }
  } catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillSpectralWindow");
}  

/**
 * This function fills the MS State table.
 * given :
 * @parameter r_p a pointer on the row of the ASDM Main table being processed.
 *
 */

void fillState(MainRow* r_p) {
  LOGENTER("fillState");

  ASDM&			ds	   = r_p -> getTable() . getContainer();
  ScanRow*		scanR_p	   = ds.getScan().getRowByKey(r_p -> getExecBlockId(),	r_p -> getScanNumber());
  vector<ScanIntentMod::ScanIntent>	scanIntent = scanR_p -> getScanIntent();
  SubscanRow*		sscanR_p   = ds.getSubscan().getRowByKey(r_p -> getExecBlockId(),
								 r_p -> getScanNumber(),
								 r_p -> getSubscanNumber());
  if (sscanR_p == 0) {
    errstream.str("");
    errstream << "Could not find a row in the Subscan table for the following key value (execBlockId=" << r_p->getExecBlockId().toString()
	      <<", scanNumber="<< r_p->getScanNumber()
	      <<", subscanNum=" << r_p->getSubscanNumber() << ").";
    throw ASDM2MSException(errstream.str());
  }	  

  SubscanIntent subscanIntent = sscanR_p->getSubscanIntent();
  string obs_mode;
  if (scanIntent.size() > 0) {
    obs_mode = CScanIntent::name(scanIntent.at(0))+"#"+CSubscanIntent::name(subscanIntent);
   
    for (unsigned int iScanIntent = 1; iScanIntent < scanIntent.size(); iScanIntent++) {
      obs_mode += ",";
      obs_mode +=  CScanIntent::name(scanIntent.at(iScanIntent))+"#"+CSubscanIntent::name(subscanIntent);
    }
  }

  const vector<StateRow *>& sRs =  ds.getState().get() ;
  for (unsigned int iState = 0; iState < sRs.size(); iState++) {							     	    
    bool pushed = false;
    
    for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	 iter != msFillers.end();
	 ++iter) {
      int retId = iter->second->addUniqueState(sRs[iState]->getSig(),
					       sRs[iState]->getRef(),
					       0.0, 
					       0.0, 
					       r_p->getSubscanNumber(),
					       obs_mode, 
					       false);
      if (!pushed) {
	stateIdx2Idx[r_p] = retId;
	pushed = true;
      }
    }	    
  }
  LOGEXIT("fillState");
}

template<typename T>
void v2oss(std::vector<T> v,
	   ostringstream& oss,
	   const std::string& oChar,
	   const std::string& cChar,
	   const std::string& sepChar) {
  oss << oChar;
  if (v.size() > 0) {
    oss << v[0];
    if (v.size() > 1) {
      for (unsigned int i = 1; i < v.size(); i++) 
	oss << sepChar << v[i];
    }
  }
  oss << cChar;
}

typedef struct MainRowCUStruct {
  MainRow* mR_p;             // A pointer on a row of the Main ASDM table,
  string   bdfName;          // The path to the BDF which contains the data,
  int32_t  index;            // The (0-based) index of the row in the Main ASDM table, 
  bool     uncorrected;       // true if this row contains uncorrected data,
  bool     corrected;        // true if this row contains corrected data,
} MainRowCUStruct;

void fillMainLazily( const string& dsName, ASDM* ds_p,
                     std::map<int, std::set<int> >& selected_eb_scan_m,
                     std::map<unsigned int , double>& effectiveBwPerDD_m,
                     Enum<CorrelationModeMod::CorrelationMode> e_query_cm, bool checkdupints) {

  LOGENTER("fillMainLazily");

  ostringstream oss;

  MainTable&			mainT	= ds_p->getMain();
  ConfigDescriptionTable&	cfgDscT = ds_p->getConfigDescription();
  ProcessorTable&		procT   = ds_p->getProcessor();

  MainRowCUStruct mRCU_s;
  vector<MainRowCUStruct> mRCU_s_v;

  vector<int32_t>	mRIndexUncorrected_v;
  vector<string>	bdfNamesUncorrected_v;

  vector<int32_t>	mRIndexCorrected_v;
  vector<string>	bdfNamesCorrected_v;

  bool	produceUncorrected = msFillers.find(AtmPhaseCorrectionMod::AP_UNCORRECTED) != msFillers.end();
  bool	produceCorrected   = msFillers.find(AtmPhaseCorrectionMod::AP_CORRECTED) != msFillers.end();

  const vector<MainRow *>& temp = mainT.get();

  for ( vector<MainRow *>::const_iterator iter_v = temp.begin(); iter_v != temp.end(); iter_v++) {
    map<int, set<int> >::iterator iter_m = selected_eb_scan_m.find((*iter_v)->getExecBlockId().getTagValue());
    if ( iter_m != selected_eb_scan_m.end() && iter_m->second.find((*iter_v)->getScanNumber()) != iter_m->second.end() ) {
      string dataUID = (*iter_v)->getDataUID().getEntityId().toString();
      replace(dataUID.begin(),dataUID.end(),':','_');
      replace(dataUID.begin(),dataUID.end(),'/','_');
      string abspath = Path(dsName + "/ASDMBinary/" + dataUID).absoluteName();

      // Are these data radiometric , if yes consider them both for corrected and uncorrected ms?
      // ProcessorType processorType = procT.getRowByKey(cfgDscT.getRowByKey((*iter_v)->getConfigDescriptionId())->getProcessorId())->getProcessorType();
      ConfigDescriptionRow *configDescRow = cfgDscT.getRowByKey((*iter_v)->getConfigDescriptionId());
      ProcessorTypeMod::ProcessorType processorType = procT.getRowByKey(configDescRow->getProcessorId())->getProcessorType();
      CorrelationModeMod::CorrelationMode corrMode = configDescRow->getCorrelationMode();

      // some of these can be skipped immediately
      if ( (e_query_cm == CorrelationModeMod::AUTO_ONLY && corrMode == CorrelationModeMod::CROSS_ONLY) ||
	   (e_query_cm == CorrelationModeMod::CROSS_ONLY && corrMode == CorrelationModeMod::AUTO_ONLY) ) {
	infostream.str("");
	infostream << "Skipping file " << abspath << " due to correlation mode selection.";
	info(infostream.str());
	continue;
      }

      mRCU_s.mR_p = *iter_v; 
      mRCU_s.bdfName = abspath;
      mRCU_s.index = iter_v - temp.begin();

      if (processorType == ProcessorTypeMod::RADIOMETER) {
	// one considers that radiometric are uncorrected and corrected data.
	mRCU_s.uncorrected = true; mRCU_s.corrected = true; 
      }
      else if (processorType == ProcessorTypeMod::CORRELATOR) {
	// We are in front of CORRELATOR data. what's their status regarding AP correction ?
          vector<AtmPhaseCorrectionMod::AtmPhaseCorrection> apc_v =  cfgDscT.getRowByKey((*iter_v)->getConfigDescriptionId())->getAtmPhaseCorrection();
	mRCU_s.uncorrected = find(apc_v.begin(), apc_v.end(), AtmPhaseCorrectionMod::AP_UNCORRECTED) != apc_v.end();
	mRCU_s.corrected   = find(apc_v.begin(), apc_v.end(), AtmPhaseCorrectionMod::AP_CORRECTED) != apc_v.end(); 
      }
      
      if ( mRCU_s.uncorrected && produceUncorrected) { 
	bdfNamesUncorrected_v.push_back(abspath);
      }

      if ( mRCU_s.corrected && produceCorrected ) {
	bdfNamesCorrected_v.push_back(abspath);
      }

      if (mRCU_s.uncorrected or mRCU_s.corrected) {
	mRCU_s_v.push_back(mRCU_s);  // We will consider only the rows which have data of at least one of the two kinds
                                     // (practically we know that the rows will have at least uncorrected data.
      }
    }
  }

  // Let's determine the byte order of the binary parts.
  // We make here the realistic but strong assumption that *all* binary parts will have the same byte order.
  bool isBigEndian;
  SDMDataObjectStreamReader sdosr;
  sdosr.open(mRCU_s_v[0].bdfName);
  isBigEndian = sdosr.byteOrder() == asdmbinaries::ByteOrder::Big_Endian;
  sdosr.close();

  // Let's have instance(s) of BDF2AsdmStMainIndex to create the asdmindex(es) which will be used by the asdmstman.
  BDF2AsdmStManIndex bdf2AsdmStManIndexU;
  BDF2AsdmStManIndex bdf2AsdmStManIndexC;

  if ( bdfNamesUncorrected_v.size() and produceUncorrected ) {
    const casacore::MeasurementSet* ms_p = msFillers.find(AtmPhaseCorrectionMod::AP_UNCORRECTED)->second->ms();
    oss.str("");
    if (ms_p->tableDesc().isColumn("DATA")) {
      oss << RODataManAccessor(*ms_p, "DATA", true).dataManagerSeqNr();
    } else {
      oss << RODataManAccessor(*ms_p, "FLOAT_DATA", true).dataManagerSeqNr();
    }
    bdf2AsdmStManIndexU.init(bdfNamesUncorrected_v, isBigEndian, ms_p->tableName() + "/table.f" + String(oss.str()));
  }
  
  if ( bdfNamesCorrected_v.size() and produceCorrected ) {
    const casacore::MeasurementSet* ms_p = msFillers.find(AtmPhaseCorrectionMod::AP_CORRECTED)->second->ms();
    oss.str("");
    if (ms_p->tableDesc().isColumn("DATA")) {
      oss << RODataManAccessor(*ms_p, "DATA", true).dataManagerSeqNr();
    } else {
      oss << RODataManAccessor(*ms_p, "FLOAT_DATA", true).dataManagerSeqNr();
    }
    bdf2AsdmStManIndexC.init(bdfNamesCorrected_v, isBigEndian, ms_p->tableName() + "/table.f" + String(oss.str()));
  } 

  // Initialize an UVW coordinates engine.
  UvwCoords uvwCoords(ds_p);  

  //
  // Some informations
  // 
  infostream.str("");
  infostream << "The dataset has " << mainT.size() << " main(s)...";
  if ( bdfNamesUncorrected_v.size() and produceUncorrected ) 
    infostream << bdfNamesUncorrected_v.size() << " of them in the selected exec blocks / scans for the uncorrected data." << endl;
  if ( bdfNamesCorrected_v.size() and produceCorrected )
    infostream << bdfNamesCorrected_v.size() << " of them in the selected exec blocks / scans for the corrected data." << endl;
  info(infostream.str());

  // Now traverse the BDFs : 
  //   * to write the indexes for asdmstman
  //   * to populate all the columns other than the DATA's one in the non lazy way.
  //

  uInt		lastMSNUrows = 0;
  uInt		lastMSNCrows = 0;

  // used in checking for duplicate integrations in the WVR (Radiometer) case
  // This holds the most recent last integration time for each configDescriptionId - but only for Radiometer data.
  map<Tag, double> lastTimeMap;

  try {
    unsigned int mainRowIndex;
    vector<MainRowCUStruct>::iterator iter;
    for (iter=mRCU_s_v.begin(), mainRowIndex=0; iter!=mRCU_s_v.end(); iter++, mainRowIndex++) {
      MainRow* mR_p = iter->mR_p;

      SDMDataObjectStreamReader sdosr;
      sdosr.open(iter->bdfName);
      LOG("Processing " + iter->bdfName);
      unsigned int numberOfAntennas = sdosr.numAntenna();
      unsigned int numberOfBaselines = numberOfAntennas * (numberOfAntennas - 1) / 2 ;
      ProcessorTypeMod::ProcessorType processorType = sdosr.processorType();
      CorrelationModeMod::CorrelationMode correlationMode = sdosr.correlationMode();
      const SDMDataObject::DataStruct& dataStruct = sdosr.dataStruct();

      infostream.str("");
      infostream << "ASDM Main row #" << iter->index
		 << " (scan #" << mR_p->getScanNumber()
		 <<", subscan #" <<  mR_p->getSubscanNumber()
		 <<", " << mR_p->getConfigDescriptionId().toString() << ") "
		 << " contains data produced by a '" << CProcessorType::name(processorType) << "'." ;
      info(infostream.str());

      // skip rows that don't match the selected mode
      if ( (correlationMode == CorrelationModeMod::CROSS_ONLY && e_query_cm == CorrelationModeMod::AUTO_ONLY) || 
	   (correlationMode == CorrelationModeMod::AUTO_ONLY && e_query_cm == CorrelationModeMod::CROSS_ONLY) ) {
	infostream.str("");
	infostream << "The main row # " << iter->index << " is ignored because the correlationMode is excluded by the selected mode to fill.";
	info(infostream.str());
	continue;
      }

      /**
       * Take care of the MS State table prior to the Main.
       */
      try {
	fillState(mR_p);
      }
      catch (ASDM2MSException e) {
	// The State table could not be filled, then let's forget this Main table row.
	infostream.str("");
	infostream << e.getMessage() << ". The main row # " << iter->index << " is ignored.";
	info(infostream.str());
	continue;
      }
      
      /**
       * And then work on the MS Main rows
       */
      ConfigDescriptionRow*	cdR		   = ds_p->getConfigDescription().getRowByKey(mR_p->getConfigDescriptionId());
      vector<Tag>		antennaIds	   = cdR->getAntennaId();
      vector<Tag>		dataDescriptionIds = cdR->getDataDescriptionId();
      vector<int>		feedIds		   = cdR->getFeedId();
      int			fieldId		   = mR_p->getFieldId().getTagValue();
      int			observationId	   = mR_p->getExecBlockId().getTagValue();
      int			processorId	   = cdR->getProcessorId().getTagValue();
      int			scanNumber	   = mR_p->getScanNumber();
      int			arrayId		   = 0;
            
      if (iter->uncorrected and produceUncorrected)
	bdf2AsdmStManIndexU.setNumberOfDataDescriptions(dataDescriptionIds.size());

      if (iter->corrected and produceCorrected)
	bdf2AsdmStManIndexC.setNumberOfDataDescriptions(dataDescriptionIds.size());
      
      unsigned int numberOfSpectralWindows = 0;
      for (const SDMDataObject::Baseband& bb: dataStruct.basebands()) {
	numberOfSpectralWindows += bb.spectralWindows().size();
      }

      if (debug) {
	oss.str("");
	oss << "There are " << numberOfSpectralWindows << " spectral windows." << endl;
	LOG(oss.str());
	oss.str("");
	oss << "There are " << dataDescriptionIds.size() << " data descriptions." << endl;
	LOG(oss.str());
      }

      vector<unsigned int> numberOfChannels_v;
      vector<unsigned int> numberOfSDPolarizations_v;
      vector<unsigned int> numberOfCrossPolarizations_v;
      for (const SDMDataObject::Baseband& bb: dataStruct.basebands()) {
	for (const SDMDataObject::SpectralWindow& spw: bb.spectralWindows()) {
	  numberOfChannels_v.push_back(spw.numSpectralPoint());
	  if (correlationMode != CorrelationModeMod::AUTO_ONLY)
	    numberOfCrossPolarizations_v.push_back(spw.crossPolProducts().size());
	  else
	    numberOfCrossPolarizations_v.push_back(0);

	  if (correlationMode != CorrelationModeMod::CROSS_ONLY)
	    numberOfSDPolarizations_v.push_back(spw.sdPolProducts().size());
	  else
	    numberOfSDPolarizations_v.push_back(0);
	}
      }
      if (debug) {
	oss.str("");
	oss << "numbers of Channels : " ;
	v2oss(numberOfChannels_v, oss, "{", "}", ", "); 
	LOG(oss.str());
	oss.str("");
	oss << "numbers of SD Polarizations : ";
	v2oss(numberOfSDPolarizations_v, oss, "{", "}", ", "); 
	oss << "numbers of Cross Polarizations : ";
	v2oss(numberOfCrossPolarizations_v, oss, "{", "}", ", "); 
	LOG(oss.str());
      }

      // Prepare vectors of scale factors
      vector<double> crossScaleFactors;
      vector<double> autoScaleFactors;
      
      // The cross data scale factors exist.
      if (correlationMode != CorrelationModeMod::AUTO_ONLY) { 
	for (const SDMDataObject::Baseband& bb: dataStruct.basebands()) {
	  for (const SDMDataObject::SpectralWindow& spw: bb.spectralWindows()) {
	    crossScaleFactors.push_back(spw.scaleFactor());
	  }
	}
	if (debug) {
	  oss.str("");
	  oss << "crossScaleFactors : " ;
	  v2oss(crossScaleFactors, oss, "{", "}", ", "); 
	  LOG(oss.str());
	}
      }
      
      // The auto data scale factors are fake.
      if (correlationMode != CorrelationModeMod::CROSS_ONLY) {
	for (unsigned int i = 0; i < numberOfSpectralWindows; i++)
	  autoScaleFactors.push_back(1.0);
	if (debug) {
	  oss.str("");
	  oss << "autoScaleFactors : " ;
	  v2oss(autoScaleFactors, oss, "{", "}", ", "); 
	  LOG(oss.str());
	}
      }
            
      // 
      // The number of values between to consecutive baselines (or antennas), stepBl, is :
      //
      unsigned int	stepSDBl     = 0; 
      unsigned int	stepCrossBl  = 0;

      int factor = (iter->uncorrected and iter->corrected) ? 2 : 1; 
      for (unsigned int i = 0; i < numberOfSpectralWindows; i++) {
	stepSDBl		    += numberOfChannels_v[i]*(numberOfSDPolarizations_v[i]==3?4:numberOfSDPolarizations_v[i]);
	stepCrossBl		    += factor * numberOfChannels_v[i]*numberOfCrossPolarizations_v[i];
      }

      if (debug) {
	oss.str("");

	oss << "stepSDBl : " << stepSDBl << "\n" << "stepCrossBl : " << stepCrossBl; 
	LOG(oss.str());
      }

      //
      // The offsets to the beginning of the i-th spectral window, spwOffset_v, is:
      std::vector<uint32_t>	spwSDOffset_v(numberOfSpectralWindows);
      //std::vector<uint32_t>	spwCrossOffset_v(numberOfSpectralWindows);
      std::vector<uint32_t>	spwCrossOffsetU_v(numberOfSpectralWindows);
      std::vector<uint32_t>	spwCrossOffsetC_v(numberOfSpectralWindows);

      spwSDOffset_v[0]	   = 0;
      spwCrossOffsetU_v[0] = 0;
      spwCrossOffsetC_v[0] = numberOfChannels_v[0] * numberOfCrossPolarizations_v[0];

      bool onlyUncorrected = iter->uncorrected and !iter->corrected;
      bool uncorrectedANDcorrected = iter->uncorrected and iter->corrected;
      if (!onlyUncorrected and !uncorrectedANDcorrected) {
	errstream.str("");
	errstream.str("I don't know how to process data with uncorrected = " +
		      TO_STRING(iter->uncorrected) +
		      " and corrected = " +
		      TO_STRING(iter->corrected) );
	error(errstream.str());
      }

      for (uint32_t i = 1; i < numberOfSpectralWindows; i++) {
	spwSDOffset_v[i] = spwSDOffset_v[i-1] +
	  numberOfChannels_v[i-1] * (numberOfSDPolarizations_v[i-1]==3?4:numberOfSDPolarizations_v[i-1]);
	if (onlyUncorrected)
	  spwCrossOffsetU_v[i] = spwCrossOffsetU_v[i-1] +
	    numberOfChannels_v[i-1] * numberOfCrossPolarizations_v[i-1];
	else if (uncorrectedANDcorrected){
	  spwCrossOffsetU_v[i] = spwCrossOffsetU_v[i-1] +
	    2 * numberOfChannels_v[i-1] * numberOfCrossPolarizations_v[i-1];
	  spwCrossOffsetC_v[i] = spwCrossOffsetC_v[i-1] +
	    numberOfChannels_v[i-1] * numberOfCrossPolarizations_v[i-1]+
	    numberOfChannels_v[i] * numberOfCrossPolarizations_v[i];
	}	  
      }

      if (debug) {
	oss.str("");
	//oss << "spwOffset_v : " ;
	//v2oss(spwOffset_v, oss, "{", "}", ", "); 
	oss << "spwSDOffset_v : " ;
	v2oss(spwSDOffset_v, oss, "{", "}", ", ");
	oss << "spwCrossOffsetU_v : " ;
	v2oss(spwCrossOffsetU_v, oss, "{", "}", ", "); 
	oss << "spwCrossOffsetC_v : " ;
	v2oss(spwCrossOffsetC_v, oss, "{", "}", ", "); 
	LOG(oss.str());
      }

      //
      // Now delegate to bdf2AsdmStManIndex the creation of the AsmdIndex 'es.
      // 
      if (processorType == ProcessorTypeMod::RADIOMETER && sdosr.hasPackedData()) {

	//
	// Declare some containers required to populate the columns of the MS MAIN table in a non lazy way.
	vector<vector<int> >           antenna1_vv(dataDescriptionIds.size());	// Column ANTENNA1
	vector<vector<int> >           antenna2_vv(dataDescriptionIds.size());	// Column ANTENNA2
	vector<vector<int> >           dataDescId_vv(dataDescriptionIds.size());	// Column DATA_DESC_ID
	vector<vector<double> >        exposure_vv(dataDescriptionIds.size());	// Column EXPOSURE
	vector<vector<double> >        interval_vv(dataDescriptionIds.size());	// Column INTERVAL
	vector<vector<double> >        time_vv(dataDescriptionIds.size());	// Column TIME    
	vector<vector<int> >           feed1_vv(dataDescriptionIds.size());	// Column FEED1
	vector<vector<int> >           feed2_vv(dataDescriptionIds.size());	// Column FEED2
	vector<vector<bool> >          flagRow_vv(dataDescriptionIds.size());	// Column FLAG_ROW
	vector<vector<int> >           stateId_vv(dataDescriptionIds.size());	// Column STATE_ID
	vector<vector<double> >        timeCentroid_vv(dataDescriptionIds.size());	// Column TIME_CENTROID
	vector<vector<pair<int, int> > >    nChanNPol_vv(dataDescriptionIds.size());  // numChan , numPol information 
	vector<vector<double> >        uvw_vv(dataDescriptionIds.size());       // Column UVW
	vector<vector<double> >        weight_vv(dataDescriptionIds.size());    // Column WEIGHT
	vector<vector<double> >        sigma_vv(dataDescriptionIds.size());     // Column SIGMA

	//
	// Everything is contained in *one* SDMDataSubset.
	//
	const SDMDataSubset& sdmDataSubset = sdosr.getSubset();

	int64_t  deltaTime = sdmDataSubset.interval() / sdosr.numTime();
	int64_t startTime = (int64_t)sdmDataSubset.time() -  (int64_t)sdmDataSubset.interval()/2LL + deltaTime/2LL;
	double   interval = deltaTime / 1000000000.0;
	
	// should the first integration be skipped? Any actual skipping happens later.
	bool skipFirstIntegration = checkdupints && lastTimeMap[mR_p->getConfigDescriptionId()] == ArrayTime(startTime).getMJD();
	if (debug && skipFirstIntegration) {
	  cout << "Duplicate time seen in Row : " << mainRowIndex
	       << " cdId : " << mR_p->getConfigDescriptionId()
	       << " " << mR_p->getDataUID().getEntityId().toString()
	       << " numTime : " << sdosr.numTime()
	       << " num MS rows : " << dataDescriptionIds.size()*sdosr.numTime()*antennaIds.size()
	       << endl;
	}
	lastTimeMap[mR_p->getConfigDescriptionId()] = ArrayTime(startTime+(sdosr.numTime()-1)*deltaTime).getMJD();
	
	for (unsigned int iDD = 0; iDD < dataDescriptionIds.size(); iDD++) {
	  //
	  // Prepare a pair<int, int> to transport the shape of some cells
	  //
	  pair<int,int> nChanNPol = make_pair<int, int>(numberOfChannels_v[iDD],
							numberOfSDPolarizations_v[iDD]);

	  //
	  // Compute weight and sigma which depend on the data description id and on the interval
	  //
	  double	weight = 1.0 * effectiveBwPerDD_m[dataDescriptionIdx2Idx[dataDescriptionIds[iDD].getTagValue()]] * interval;
	  weight	       = (weight == 0.0) ? 1.0 : weight;
	  double	sigma  = 1./sqrt(weight);

	  for (unsigned int itime = 0; itime < sdosr.numTime(); itime++) {
	    if (skipFirstIntegration && itime==0) continue;
	    for (unsigned int iA = 0; iA < antennaIds.size(); iA++) {
	      antenna1_vv[iDD].push_back(antennaIds[iA].getTagValue());
	      antenna2_vv[iDD].push_back(antennaIds[iA].getTagValue());
	      dataDescId_vv[iDD].push_back(dataDescriptionIdx2Idx[dataDescriptionIds[iDD].getTagValue()]);
	      exposure_vv[iDD].push_back(interval);
	      interval_vv[iDD].push_back(interval);
	      time_vv[iDD].push_back(ArrayTime(startTime + itime * deltaTime).getMJD() * 86400.0);
	      feed1_vv[iDD].push_back(feedIds[iA]);
	      feed2_vv[iDD].push_back(feedIds[iA]);
	      flagRow_vv[iDD].push_back(false);
	      stateId_vv[iDD].push_back(stateIdx2Idx[mR_p]);
	      timeCentroid_vv[iDD].push_back(time_vv[iDD].back());
	      nChanNPol_vv[iDD].push_back(nChanNPol);
	      uvw_vv[iDD].push_back(0.0);uvw_vv[iDD].push_back(0.0);uvw_vv[iDD].push_back(0.0);
	      weight_vv[iDD].push_back(weight);
	      sigma_vv[iDD].push_back(sigma);
	    }
	    // If we have uncorrected data (which is in practice always the case I think) and want to output those then populate the asdmindex of uncorrected data MS.
	    if (iter->uncorrected and produceUncorrected) {
	      bdf2AsdmStManIndexU.appendWVRIndex(iDD,
						 iter->bdfName,
						 numberOfAntennas,
						 numberOfSpectralWindows,
						 numberOfChannels_v[iDD],
						 numberOfSDPolarizations_v[iDD],
						 stepSDBl, //numberOfSpectralWindows * numberOfChannels * numberOfPolarizations,
						 iDD, // this will be used as an index in the seq of windows in the BDFs
						 autoScaleFactors,
						 sdmDataSubset.autoDataPosition() + itime * numberOfAntennas * stepSDBl * sizeof(AUTODATATYPE),
						 spwSDOffset_v[iDD]);
	    }
	    // If we have corrected data  and want to output those  then populate the asdmindex of corrected data MS.
	    if (iter->corrected and produceCorrected) {
	      bdf2AsdmStManIndexC.appendWVRIndex(iDD,
						 iter->bdfName,
						 numberOfAntennas,
						 numberOfSpectralWindows,
						 numberOfChannels_v[iDD],
						 numberOfSDPolarizations_v[iDD],
						 stepSDBl, //numberOfSpectralWindows * numberOfChannels * numberOfPolarizations,
						 iDD, // this will be used as an index in the seq of windows in the BDFs
						 autoScaleFactors,
						 sdmDataSubset.autoDataPosition() + itime * numberOfAntennas * stepSDBl * sizeof(AUTODATATYPE),
						 spwSDOffset_v[iDD]);
	    }  
	  }
	}

	// If we have uncorrected data (which is in practice always the case I think) and want to output those then populate the asdmindex of uncorrected data MS.
	if (iter->uncorrected and produceUncorrected) 
	  bdf2AsdmStManIndexU.dumpAutoCross();

	// If we have corrected data  and want to output those  then populate the asdmindex of corrected data MS.
	if (iter->corrected and produceCorrected)
	  bdf2AsdmStManIndexC.dumpAutoCross();

	//
	// It's now time to populate the columns of the MAIN table but the DATA's one.
	for (unsigned int iDD = 0; iDD < dataDescriptionIds.size(); iDD++) {
            for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator msfIter = msFillers.begin();
	       msfIter != msFillers.end();
	       ++msfIter) {
	    if (time_vv[iDD].size() > 0) {
	      msfIter->second->addData(true,             // Yes ! these are complex data.
				       time_vv[iDD],
				       antenna1_vv[iDD],
				       antenna2_vv[iDD],
				       feed1_vv[iDD],
				       feed2_vv[iDD],
				       dataDescId_vv[iDD],
				       processorId,
				       fieldId,
				       interval_vv[iDD],
				       exposure_vv[iDD],
				       timeCentroid_vv[iDD],
				       scanNumber, 
				       arrayId,
				       observationId,
				       stateId_vv[iDD],
				       nChanNPol_vv[iDD],
				       uvw_vv[iDD],
				       weight_vv[iDD],
				       sigma_vv[iDD]);
	    }
	  }
	}
      }

      else if (!sdosr.hasPackedData() && (processorType == ProcessorTypeMod::CORRELATOR || processorType == ProcessorTypeMod::RADIOMETER)) {

	// duplicate time skipping is not supported here

	//
	// Declare some containers required to populate the columns of the MS MAIN table in a non lazy way.
	//
	// We use vectors of vectors in order to be able to build separate vectors for different data description
	// and then output these vectors in the appropriate order.
	//
	// The cross correlation chapter.
	//
	vector<vector<int> >     cross_antenna1_vv(dataDescriptionIds.size());      // Column ANTENNA1 per Data Description
	vector<vector<int> >     cross_antenna2_vv(dataDescriptionIds.size());      // Column ANTENNA2 per Data Description
	vector<vector<int> >     cross_dataDescId_vv(dataDescriptionIds.size());    // Column DATA_DESC_ID per Data Description
	vector<vector<double> >  cross_exposure_vv(dataDescriptionIds.size());      // Column EXPOSURE per Data Description
	vector<vector<double> >  cross_interval_vv(dataDescriptionIds.size());      // Column INTERVAL per Data Description
	vector<vector<double> >  cross_time_vv(dataDescriptionIds.size());          // Column TIME per Data Description    
	vector<vector<int> >     cross_feed1_vv(dataDescriptionIds.size());         // Column FEED1 per Data Description
	vector<vector<int> >     cross_feed2_vv(dataDescriptionIds.size());         // Column FEED2 per Data Description
	vector<vector<bool> >    cross_flagRow_vv(dataDescriptionIds.size());       // Column FLAG_ROW per Data Description
	vector<vector<int> >     cross_stateId_vv(dataDescriptionIds.size());       // Column STATE_ID per Data Description
	vector<vector<double> >  cross_timeCentroid_vv(dataDescriptionIds.size());  // Column TIME_CENTROID per Data Description
	vector<vector<pair<int, int> > >    cross_nChanNPol_vv(dataDescriptionIds.size());  // numChan , numPol information 
	vector<vector<double> >  cross_uvw_vv(dataDescriptionIds.size());           // Column UVW
	vector<vector<double> >  cross_weight_vv(dataDescriptionIds.size());        // Column WEIGHT
	vector<vector<double> >  cross_sigma_vv(dataDescriptionIds.size());         // Column SIGMA

	// The auto correlation chapter.
	vector<vector<int> >     auto_antenna1_vv(dataDescriptionIds.size());      // Column ANTENNA1 per Data Description
	vector<vector<int> >     auto_antenna2_vv(dataDescriptionIds.size());      // Column ANTENNA2 per Data Description
	vector<vector<int> >     auto_dataDescId_vv(dataDescriptionIds.size());    // Column DATA_DESC_ID per Data Description
	vector<vector<double> >  auto_exposure_vv(dataDescriptionIds.size());      // Column EXPOSURE per Data Description
	vector<vector<double> >  auto_interval_vv(dataDescriptionIds.size());      // Column INTERVAL per Data Description
	vector<vector<double> >  auto_time_vv(dataDescriptionIds.size());          // Column TIME per Data Description    
	vector<vector<int> >     auto_feed1_vv(dataDescriptionIds.size());         // Column FEED1 per Data Description
	vector<vector<int> >     auto_feed2_vv(dataDescriptionIds.size());         // Column FEED2 per Data Description
	vector<vector<bool> >    auto_flagRow_vv(dataDescriptionIds.size());       // Column FLAG_ROW per Data Description
	vector<vector<int> >     auto_stateId_vv(dataDescriptionIds.size());       // Column STATE_ID per Data Description
	vector<vector<double> >  auto_timeCentroid_vv(dataDescriptionIds.size());  // Column TIME_CENTROID per Data Description
	vector<vector<pair<int, int> > >    auto_nChanNPol_vv(dataDescriptionIds.size());  // numChan , numPol information 
	vector<vector<double> >  auto_uvw_vv(dataDescriptionIds.size());           // Column UVW
	vector<vector<double> >  auto_weight_vv(dataDescriptionIds.size());           // Column WEIGHT
	vector<vector<double> >  auto_sigma_vv(dataDescriptionIds.size());            // Column SIGMA
	
	// Ignore cross data if AUTO_ONLY is selected.
	bool hasCrossData = ((correlationMode == CorrelationModeMod::CROSS_AND_AUTO || correlationMode == CorrelationModeMod::CROSS_ONLY) && e_query_cm != CorrelationModeMod::AUTO_ONLY);

	// Ignore auto data if CROSS_ONLY has been selected.
	bool hasAutoData = ((correlationMode == CorrelationModeMod::CROSS_AND_AUTO || correlationMode == CorrelationModeMod::AUTO_ONLY) && e_query_cm != CorrelationModeMod::CROSS_ONLY);

	//
	// Traverse all the integrations.
	//

	while (sdosr.hasSubset()) {

	  const SDMDataSubset& sdmDataSubset = sdosr.getSubset();
	  
	  string time_s = ArrayTime((int64_t) sdmDataSubset.time()).toFITS();
	  double time = ArrayTime((int64_t) sdmDataSubset.time()).getMJD() * 86400.0;
	  double interval =  sdmDataSubset.interval() / 1000000000.0;
#if TRANSPOSE_BL_NUM
	  pair<bool, bool> dataOrder(true, false);  // 1st: reverse bls YES, 2nd: autotrailing NO
#else
	  pair<bool, bool> dataOrder(false, false);  // 1st: reverse bls NO, 2nd: autotrailing NO
#endif
	  vector<Vector<casacore::Double> > vv_uvw;
	  vector<double> time_v(dataDescriptionIds.size() * (numberOfBaselines + numberOfAntennas),
				time);

	  if ( correlationMode != CorrelationModeMod::AUTO_ONLY ) {
	    uvwCoords.uvw_bl(mR_p,
			     time_v, 
			     correlationMode,
			     dataOrder,
			     vv_uvw);
	  }
	  
	  //
	  // If we have autocorrelations and cross correlations , ignore the numberOfAntennas * dataDescriptionIds.size()
	  // first element of vv_uvw
	  // 
	  unsigned int uvwIndexBase = 0;
	  if (correlationMode == CorrelationModeMod::CROSS_AND_AUTO) {
	    uvwIndexBase += numberOfAntennas * dataDescriptionIds.size();
	  }

	  if (hasCrossData) {
	    for (unsigned int iDD = 0; iDD < dataDescriptionIds.size(); iDD++) {
	      unsigned int uvwIndex = uvwIndexBase + iDD;
	      unsigned int ddIndex = dataDescriptionIdx2Idx[dataDescriptionIds[iDD].getTagValue()];

	      //
	      // Prepare a pair<int, int> to transport the shape of some cells
	      //
	      pair<int,int> nChanNPol = make_pair<int, int>(numberOfChannels_v[iDD],
							    numberOfCrossPolarizations_v[iDD]);
	      //
	      // Compute weight and sigma which depend on interval and iDD.
	      //
	      double weight = 2.0*interval*effectiveBwPerDD_m[dataDescriptionIdx2Idx[dataDescriptionIds[iDD].getTagValue()]];
	      weight = (weight == 0.0) ? 1.0 : weight;
	      double sigma = 1.0 / sqrt (weight);

#if !TRANSPOSE_BL_NUM
	      for (unsigned int iA1 = 0; iA1 < antennaIds.size(); iA1++)
		for (unsigned int iA2 = iA1+1; iA2 < antennaIds.size(); iA2++) {
		  cross_antenna1_vv[iDD].push_back(antennaIds[iA1].getTagValue());
		  cross_antenna2_vv[iDD].push_back(antennaIds[iA2].getTagValue());
		  cross_dataDescId_vv[iDD].push_back(ddIndex);
		  cross_exposure_vv[iDD].push_back(interval);
		  cross_interval_vv[iDD].push_back(interval);
		  cross_time_vv[iDD].push_back(time);
		  cross_feed1_vv[iDD].push_back(feedIds[iA1]);
		  cross_feed2_vv[iDD].push_back(feedIds[iA2]);
		  cross_flagRow_vv[iDD].push_back(false);
		  cross_stateId_vv[iDD].push_back(stateIdx2Idx[mR_p]);
		  cross_timeCentroid_vv[iDD].push_back(time);
		  cross_nChanNPol_vv[iDD].push_back(nChanNPol);
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](0));
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](1));
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](2));
		  uvwIndex += dataDescriptionIds.size();
		  cross_weight_vv[iDD].push_back(weight);
		  cross_sigma_vv[iDD].push_back(sigma);
		}
#else
	      for (unsigned int iA2 = 1; iA2 < antennaIds.size(); iA2++)
		for (unsigned int iA1 = 0; iA1 < iA2; iA1++) {
		  cross_antenna1_vv[iDD].push_back(antennaIds[iA1].getTagValue());
		  cross_antenna2_vv[iDD].push_back(antennaIds[iA2].getTagValue());
		  cross_dataDescId_vv[iDD].push_back(ddIndex);
		  cross_exposure_vv[iDD].push_back(interval);
		  cross_interval_vv[iDD].push_back(interval);
		  cross_time_vv[iDD].push_back(time);
		  cross_feed1_vv[iDD].push_back(feedIds[iA1]);
		  cross_feed2_vv[iDD].push_back(feedIds[iA2]);
		  cross_flagRow_vv[iDD].push_back(false);
		  cross_stateId_vv[iDD].push_back(stateIdx2Idx[mR_p]);
		  cross_timeCentroid_vv[iDD].push_back(time);
		  cross_nChanNPol_vv[iDD].push_back(nChanNPol);
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](0));
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](1));
		  cross_uvw_vv[iDD].push_back(vv_uvw[uvwIndex](2));
		  uvwIndex += dataDescriptionIds.size();
		  cross_weight_vv[iDD].push_back(weight);
		  cross_sigma_vv[iDD].push_back(sigma);
		}
#endif	     
	      // If we have uncorrected data (which is in practice always the case I think) and want to output those then populate the asdmindex of uncorrected data MS.
	      if (iter->uncorrected and produceUncorrected) 
		bdf2AsdmStManIndexU.appendCrossIndex(iDD,
						     iter->bdfName,
						     numberOfBaselines,
						     numberOfSpectralWindows,
						     numberOfChannels_v[iDD],
						     numberOfCrossPolarizations_v[iDD],
						     stepCrossBl, 
						     iDD, // this will be used as an index in the seq of windows in the BDFs
						     crossScaleFactors,
						     sdmDataSubset.crossDataPosition(),
						     spwCrossOffsetU_v[iDD],
						     sdmDataSubset.crossDataType());

	      // If we have corrected data  and want to output those  then populate the asdmindex of corrected data MS.
	      if (iter->corrected and produceCorrected) 
		bdf2AsdmStManIndexC.appendCrossIndex(iDD,
						     iter->bdfName,
						     numberOfBaselines,
						     numberOfSpectralWindows,
						     numberOfChannels_v[iDD],
						     numberOfCrossPolarizations_v[iDD],
						     stepCrossBl, 
						     iDD, // this will be used as an index in the seq of windows in the BDFs
						     crossScaleFactors,
						     sdmDataSubset.crossDataPosition(),
						     spwCrossOffsetC_v[iDD],
						     sdmDataSubset.crossDataType());
	    }
	  }
	  
	  if (hasAutoData) {
	    for (unsigned int iDD = 0; iDD < dataDescriptionIds.size(); iDD++) {
	      unsigned int ddIndex = dataDescriptionIdx2Idx[dataDescriptionIds[iDD].getTagValue()];
	      //
	      // Prepare a pair<int, int> to transport the shape of some cells
	      //
	      pair<int,int> nChanNPol = make_pair<int, int>(numberOfChannels_v[iDD],
							    numberOfSDPolarizations_v[iDD] == 3 ? 4 : numberOfSDPolarizations_v[iDD] );

	      //
	      // Compute weight and sigma which depend on interval and iDD.
	      //
	      double weight = 1.0*interval*effectiveBwPerDD_m[ddIndex];
	      weight = (weight == 0.0) ? 1.0 : weight;
	      double sigma = 1.0 / sqrt (weight);


	      for (unsigned int iA = 0; iA < antennaIds.size(); iA++) {
		auto_antenna1_vv[iDD].push_back(antennaIds[iA].getTagValue());
		auto_antenna2_vv[iDD].push_back(antennaIds[iA].getTagValue());
		auto_dataDescId_vv[iDD].push_back(ddIndex);
		auto_exposure_vv[iDD].push_back(interval);
		auto_interval_vv[iDD].push_back(interval);
		auto_time_vv[iDD].push_back(time);
		auto_feed1_vv[iDD].push_back(feedIds[iA]);
		auto_feed2_vv[iDD].push_back(feedIds[iA]);
		auto_flagRow_vv[iDD].push_back(false);
		auto_stateId_vv[iDD].push_back(stateIdx2Idx[mR_p]);
		auto_timeCentroid_vv[iDD].push_back(time);
		auto_nChanNPol_vv[iDD].push_back(nChanNPol);
		auto_uvw_vv[iDD].push_back(0.0);auto_uvw_vv[iDD].push_back(0.0);auto_uvw_vv[iDD].push_back(0.0);
		auto_weight_vv[iDD].push_back(weight);
		auto_sigma_vv[iDD].push_back(sigma);		
	      }

	      // If we have uncorrected data (which is in practice always the case I think) and want to output those then populate the asdmindex of uncorrected data MS.
	      if (iter->uncorrected and produceUncorrected)
		bdf2AsdmStManIndexU.appendAutoIndex(iDD,
						    iter->bdfName,
						    numberOfAntennas,
						    numberOfSpectralWindows,
						    numberOfChannels_v[iDD],
						    numberOfSDPolarizations_v[iDD],
						    stepSDBl, 
						    iDD,
						    autoScaleFactors,
						    sdmDataSubset.autoDataPosition(),
						    spwSDOffset_v[iDD]);
	      // If we have corrected data  and want to output those  then populate the asdmindex of corrected data MS.
	      if (iter->corrected and produceCorrected)
		bdf2AsdmStManIndexC.appendAutoIndex(iDD,
						    iter->bdfName,
						    numberOfAntennas,
						    numberOfSpectralWindows,
						    numberOfChannels_v[iDD],
						    numberOfSDPolarizations_v[iDD],
						    stepSDBl, //numberOfSpectralWindows * numberOfChannels * numberOfSDPolarizations,
						    iDD,
						    autoScaleFactors,
						    sdmDataSubset.autoDataPosition(),
						    spwSDOffset_v[iDD]);
	    }	      
	  }
	}
	// If we have uncorrected data (which is in practice always the case I think) and want to output those then populate the asdmindex of uncorrected data MS.
	if (iter->uncorrected and produceUncorrected)
	  bdf2AsdmStManIndexU.dumpAutoCross();
	// If we have corrected data  and want to output those  then populate the asdmindex of corrected data MS.
	if (iter->corrected and produceCorrected)
	  bdf2AsdmStManIndexC.dumpAutoCross();


	//
	// It's now time to populate the columns of the MAIN table but the DATA's one.
	// This is done with data descriptions varying the more slowly.
	//
	for (unsigned int iDD = 0; iDD < dataDescriptionIds.size(); iDD++) {
	  if (hasAutoData)
              for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator msfIter = msFillers.begin();
		 msfIter != msFillers.end();
		 ++msfIter)
	      if ((msfIter->first == AtmPhaseCorrectionMod::AP_UNCORRECTED and iter->uncorrected and produceUncorrected) ||
		  (msfIter->first == AtmPhaseCorrectionMod::AP_CORRECTED and iter->corrected and produceCorrected)) {
		msfIter->second->addData(true,             // Yes ! these are complex data.
					 auto_time_vv[iDD],
					 auto_antenna1_vv[iDD],
					 auto_antenna2_vv[iDD],
					 auto_feed1_vv[iDD],
					 auto_feed2_vv[iDD],
					 auto_dataDescId_vv[iDD],
					 processorId,
					 fieldId,
					 auto_interval_vv[iDD],
					 auto_exposure_vv[iDD],
					 auto_timeCentroid_vv[iDD],
					 scanNumber, 
					 arrayId,
					 observationId,
					 auto_stateId_vv[iDD],
					 auto_nChanNPol_vv[iDD],
					 auto_uvw_vv[iDD],
					 auto_weight_vv[iDD],
					 auto_sigma_vv[iDD]);
	      }
	  if (hasCrossData) 
              for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator msfIter = msFillers.begin();
		 msfIter != msFillers.end();
		 ++msfIter)
	      if ((msfIter->first == AtmPhaseCorrectionMod::AP_UNCORRECTED and iter->uncorrected and produceUncorrected) ||
		  (msfIter->first == AtmPhaseCorrectionMod::AP_CORRECTED and iter->corrected and produceCorrected)) {
		msfIter->second->addData(true,             // Yes ! these are complex data.
					 cross_time_vv[iDD],
					 cross_antenna1_vv[iDD],
					 cross_antenna2_vv[iDD],
					 cross_feed1_vv[iDD],
					 cross_feed2_vv[iDD],
					 cross_dataDescId_vv[iDD],
					 processorId,
					 fieldId,
					 cross_interval_vv[iDD],
					 cross_exposure_vv[iDD],
					 cross_timeCentroid_vv[iDD],
					 scanNumber, 
					 arrayId,
					 observationId,
					 cross_stateId_vv[iDD],
					 cross_nChanNPol_vv[iDD],
					 cross_uvw_vv[iDD],
					 cross_weight_vv[iDD],
					 cross_sigma_vv[iDD]);
	      }
	}
      }
      else 
	cout << "Processor not supported in lazy mode." << endl;
      
      sdosr.close();

      if (iter->uncorrected and produceUncorrected) {
	infostream.str("");
	infostream << "ASDM Main row #" << iter->index << " produced a total of " << msFillers[AtmPhaseCorrectionMod::AP_UNCORRECTED]->ms()->nrow() - lastMSNUrows << " MS Main rows with uncorrected data." << endl;
	info(infostream.str());
	lastMSNUrows = msFillers[AtmPhaseCorrectionMod::AP_UNCORRECTED]->ms()->nrow();	
      }

      if (iter->corrected and produceCorrected) {
	infostream.str("");
	infostream << "ASDM Main row #" << iter->index << " produced a total of " << msFillers[AtmPhaseCorrectionMod::AP_CORRECTED]->ms()->nrow() - lastMSNCrows << " MS Main rows with corrected data." << endl;
	info(infostream.str());
	lastMSNCrows = msFillers[AtmPhaseCorrectionMod::AP_CORRECTED]->ms()->nrow();	
      }
    }
    infostream.str("");
    if (produceUncorrected) 
      infostream << "The MS main table for wvr uncorrected data contains " << msFillers[AtmPhaseCorrectionMod::AP_UNCORRECTED]->ms()->nrow() << " rows." << endl;

    if (produceCorrected) 
      infostream << "The MS main table for wvr corrected data contains " << msFillers[AtmPhaseCorrectionMod::AP_CORRECTED]->ms()->nrow() << " rows." << endl;


    info(infostream.str());
  }
  catch (SDMDataObjectStreamReaderException e) {
    cout << e.getMessage() << endl;
  }
  catch (SDMDataObjectException e) {
    cout << e.getMessage() << endl;
  }
  bdf2AsdmStManIndexU.done();
  bdf2AsdmStManIndexC.done();
  LOGEXIT("fillMainLazily");

}

template<class T> vector<T> reorder(const vector<T>& v, vector<int> index) {
  vector<T> result(v.size());
  for (unsigned int i = 0; i < result.size(); i++)
    result.at(i) = v.at(index.at(i));
  return result;
}

/**
 * This function fills the MS Main table from an ASDM Main table which refers to correlator data.
 *
 * given:
 * @parameter r_p a pointer to the MainRow being processed.
 * @parameter sdmBinData a reference to the SDMBinData containing a lot of information about the binary data being processed. Useful to know the requested ordering of data.
 * @parameter uvwCoords a reference to the UVW calculator.
 * @parameter complexData a bool which says if the DATA is going to be filled (true) or if it will be the FLOAT_DATA (false).
 * @parameter mute if the value of this parameter is false then nothing is written in the MS .
 * @parameter skipFirstTime if the value of this parameter is true, then all rows with time equal to the first time seen will be skipped (not filled).
 *
 * !!!!! One must be carefull to the fact that fillState must have been called before fillMain. During the execution of fillState , the global vector<int> msStateID
 * is filled and will be used by fillMain.
 */ 
void fillMain( MainRow* r_p, SDMBinData& sdmBinData, const VMSData* vmsData_p, UvwCoords& uvwCoords,
               std::map<unsigned int, double>& effectiveBwPerDD_m, bool complexData, bool mute,
               bool ac_xc_per_timestamp, bool skipFirstTime=false) {
  
  if (debug) cout << "fillMain : entering" << endl;

  // ASDM & ds = r_p -> getTable() . getContainer();

  // Then populate the Main table.
  ComplexDataFilter filter; // To process the case numCorr == 3
  
  if (vmsData_p->v_antennaId1.size() == 0) {
    infostream.str("");
    infostream << "No MS data produced for the current row." << endl;
    info(infostream.str());
    return;
  }

  // msRowReIndex_v maps from location in the vmsData_p vectors to the output MS row
  // this may be reorderd and the first time integration may be skipped (in which case msRowReIndex_v[i] is -1

  vector <int> msRowReIndex_v(vmsData_p->v_antennaId1.size());

  // cout << "skipFirstTime ; " << skipFirstTime << endl;

  CorrelationModeMod::CorrelationMode cm = r_p->getTable().getContainer().getConfigDescription().getRowByKey(r_p->getConfigDescriptionId())->getCorrelationMode();
  if (cm == CorrelationModeMod::CROSS_AND_AUTO and ac_xc_per_timestamp) {
    // this is the case where the rows may be reordered

    unsigned int numDD = r_p->getTable().getContainer().getConfigDescription().getRowByKey(r_p->getConfigDescriptionId())->getDataDescriptionId().size();
    unsigned int numOfMSRowsPerIntegration = 0;
    unsigned int numAntenna = r_p->getNumAntenna();
    unsigned int numBl =  r_p->getNumAntenna() * (r_p->getNumAntenna() - 1) / 2;
    if (cm != CorrelationModeMod::CROSS_ONLY) numOfMSRowsPerIntegration += numAntenna;
    if (cm != CorrelationModeMod::AUTO_ONLY) numOfMSRowsPerIntegration += numBl;
    unsigned int numIntegrations =  vmsData_p->v_antennaId1.size() / numOfMSRowsPerIntegration / numDD;
    if (vmsData_p->v_antennaId1.size() % numOfMSRowsPerIntegration != 0) {
      errstream.str("");
      errstream << "The total number of rows to be writtren into the MS (" << vmsData_p->v_antennaId1.size()
		<<") is not a multiple of the number of MS rows per integration (" << numOfMSRowsPerIntegration << ").";
      errstream << "This is not normal. Aborting !";
      throw ASDM2MSException(errstream.str());
    }

    unsigned int i = 0;
    for (unsigned int iDD = 0; iDD < numDD; iDD++) {
      unsigned int ddOffset = iDD * numIntegrations * (numAntenna + numBl);
      for (unsigned int integration = 0; integration < numIntegrations; integration++) {
	// First the auto correlations.
	for (unsigned int iAnt = 0; iAnt < numAntenna; iAnt++) {
	  if (skipFirstTime && integration == 0) {
	    msRowReIndex_v[i++] = -1;
	  } else {
	    msRowReIndex_v[i++] = ddOffset + numAntenna * integration + iAnt;
	  }
	}
	// Then the cross correlations.
	for (unsigned iBl = 0; iBl < numBl; iBl++)
	  if (skipFirstTime && integration == 0) {
	    msRowReIndex_v[i++] = -1;
	  } else {
	    msRowReIndex_v[i++] = ddOffset + numIntegrations * numAntenna + integration * numBl + iBl;
	  }	  
      }
    }
  } else {
    // set them in order, skipping times as appropriate
    // it's not obvious that the vmsData_p rows are time sorted, but the first time there should be the time to be skipped
    double timeToSkip = vmsData_p->v_time[0];
    for (unsigned int i = 0; i < msRowReIndex_v.size(); i++) {
      if (skipFirstTime && (vmsData_p->v_time[i] == timeToSkip)) {
	msRowReIndex_v[i] = -1;
      } else {
	msRowReIndex_v[i] = i;
      }
    }
  }
  
  vector<vector<unsigned int> > filteredShape_vv = vmsData_p->vv_dataShape;
  for (unsigned int ipart = 0; ipart < filteredShape_vv.size(); ipart++) {
    if (filteredShape_vv.at(ipart).at(0) == 3) filteredShape_vv.at(ipart).at(0) = 4;
  }

  vector<int> filteredDD;
  for (unsigned int idd = 0; idd < vmsData_p->v_dataDescId.size(); idd++){
    filteredDD.push_back(dataDescriptionIdx2Idx.at(vmsData_p->v_dataDescId.at(idd)));
  }

  vector<float *> uncorrectedData_v;
  vector<float *> correctedData_v;

  // these sizes of these are not known immediately if there skipFirstTime is true
  vector<double> uvw_v;              // when set here, this is reordered and can be used as is when writing to the uncorrected MS
  vector<double> weight_v;           // these are the weights, in the order seen. They must be reordered to use.
  vector<double> correctedWeight_v;  // to be put into the MS, can only be set later when other corrected column values are set
  vector<double> uncorrectedWeight_v;  // to be put into the MS. May be set here when skipFirstTime is not true.
  vector<double> sigma_v;            // the sigmas, in the order seen. They must be reordered to use.
  vector<double> correctedSigma_v;   // to be put into the MS, can only be set later when other corrected column values are set
  vector<double> uncorrectedSigma_v; // to be put into the MS. May be set here when skipFirstTime is not true.

  /* compute the UVW - do this for all times, time skipping happens later */
  vector<casacore::Vector<casacore::Double> > vv_uvw(vmsData_p->v_time.size());
#if DDPRIORITY
  uvwCoords.uvw_bl(r_p, sdmBinData.timeSequence(), e_query_cm, 
		   sdmbin::SDMBinData::dataOrder(),
		   vv_uvw);
#else
  uvwCoords.uvw_bl(r_p, vmsData_p->v_timeCentroid, e_query_cm, 
		   sdmbin::SDMBinData::dataOrder(),
		   vv_uvw);
#endif

  // set sigma_v and weight_v first, in order.
  weight_v.resize(vmsData_p->v_time.size());
  sigma_v.resize(vmsData_p->v_time.size());
  for (unsigned int i = 0; i < weight_v.size(); i++) {
    weight_v[i] = vmsData_p->v_exposure[i] * effectiveBwPerDD_m[filteredDD[i]];
    if (vmsData_p->v_antennaId1[i] != vmsData_p->v_antennaId2[i])
      weight_v[i] *= 2.0;

    if (weight_v[i] == 0.0) weight_v[i] = 1.0;

    sigma_v[i] = 1.0 / sqrt(weight_v[i]);

  }
   
  if (!skipFirstTime) {
    // several values can be fully filled out in this case. 
    // it's more efficient to just do that, most of the time this block will be executed

    /*
    ** Let's apply the reindexing on the UVW coordinates !!!
    */
    uvw_v.resize(3*vmsData_p->v_time.size());
    int k = 0;
    for (unsigned int iUvw = 0; iUvw < vv_uvw.size(); iUvw++) {
      uvw_v[k++] = vv_uvw[msRowReIndex_v[iUvw]](0); 
      uvw_v[k++] = vv_uvw[msRowReIndex_v[iUvw]](1);
      uvw_v[k++] = vv_uvw[msRowReIndex_v[iUvw]](2);
    } 

    /*
    ** Let's apply the reindexing on weight and sigma.
    */
    uncorrectedWeight_v.resize(weight_v.size());
    uncorrectedSigma_v.resize(weight_v.size());
    for (unsigned int i = 0; i < weight_v.size(); i++) {
      uncorrectedWeight_v[i] = weight_v.at(msRowReIndex_v[i]);
      uncorrectedSigma_v[i] = sigma_v.at(msRowReIndex_v[i]);
    }
  }

  ComplexDataFilter cdf;
  map<AtmPhaseCorrectionMod::AtmPhaseCorrection, float*>::const_iterator iter;

  vector<double>	uncorrectedTime_v;
  vector<int>		uncorrectedAntennaId1_v;
  vector<int>		uncorrectedAntennaId2_v;
  vector<int>		uncorrectedFeedId1_v;
  vector<int>		uncorrectedFeedId2_v;
  vector<int>		uncorrectedFieldId_v;
  vector<int>		uncorrectedFilteredDD_v;
  vector<vector< unsigned int> > uncorrectedFilteredShape_vv;
  vector<double>	uncorrectedInterval_v;
  vector<double>	uncorrectedExposure_v;
  vector<double>	uncorrectedTimeCentroid_v;
  vector<unsigned int>	uncorrectedFlag_v;

  vector<double>	correctedTime_v;
  vector<int>		correctedAntennaId1_v;
  vector<int>		correctedAntennaId2_v;
  vector<int>		correctedFeedId1_v;
  vector<int>		correctedFeedId2_v;
  vector<int>		correctedFieldId_v;
  vector<int>           correctedFilteredDD_v;
  vector<vector< unsigned int> > correctedFilteredShape_vv;
  vector<double>	correctedInterval_v;
  vector<double>	correctedExposure_v;
  vector<double>	correctedTimeCentroid_v;
  vector<double>        correctedUvw_v;
  vector<unsigned int>	correctedFlag_v;

  Tag configDescriptionId = r_p -> getConfigDescriptionId();
  ConfigDescriptionTable & cfgDescT = r_p -> getTable() . getContainer() . getConfigDescription();
  ConfigDescriptionRow * cfgDescR_p = cfgDescT.getRowByKey(configDescriptionId);
  const vector<AtmPhaseCorrectionMod::AtmPhaseCorrection >& apc_v = cfgDescR_p->getAtmPhaseCorrection();
  bool subscanHasCorrectedData = std::find(apc_v.begin(), apc_v.end(), AtmPhaseCorrectionMod::AP_CORRECTED)!=apc_v.end();

  // Do we have to fill an MS with uncorrected data + radiometric data (radiometric data are considered as uncorrected data)  ?
  // Apply here the redindexing on uncorrected data !!!!
  //
  for (unsigned int iData = 0; iData < vmsData_p->v_m_data.size(); iData++) {
    if (skipFirstTime && msRowReIndex_v[iData] < 0) {
      continue;
    }

    if ((iter=vmsData_p->v_m_data.at(msRowReIndex_v[iData]).find(AtmPhaseCorrectionMod::AP_UNCORRECTED)) != vmsData_p->v_m_data.at(msRowReIndex_v[iData]).end()){
      uncorrectedData_v.push_back(cdf.to4Pol(vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(0),
					     vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(1),
					     iter->second));
      if (skipFirstTime) {
	// uncorrected values must also be explicitly added here since they can not be put in as entire vectors
	uncorrectedTime_v.push_back(vmsData_p->v_time.at(msRowReIndex_v[iData]));
	uncorrectedAntennaId1_v.push_back(vmsData_p->v_antennaId1.at(msRowReIndex_v[iData]));
	uncorrectedAntennaId2_v.push_back(vmsData_p->v_antennaId2.at(msRowReIndex_v[iData]));
	uncorrectedFeedId1_v.push_back(vmsData_p->v_feedId1.at(msRowReIndex_v[iData]));
	uncorrectedFeedId2_v.push_back(vmsData_p->v_feedId2.at(msRowReIndex_v[iData]));
	uncorrectedFilteredDD_v.push_back(filteredDD.at(msRowReIndex_v[iData]));
	uncorrectedFilteredShape_vv.push_back(filteredShape_vv.at(msRowReIndex_v[iData]));
	uncorrectedFieldId_v.push_back(vmsData_p->v_fieldId.at(msRowReIndex_v[iData]));
	uncorrectedInterval_v.push_back(vmsData_p->v_interval.at(msRowReIndex_v[iData]));
	uncorrectedExposure_v.push_back(vmsData_p->v_exposure.at(msRowReIndex_v[iData]));
	uncorrectedTimeCentroid_v.push_back(vmsData_p->v_timeCentroid.at(msRowReIndex_v[iData]));
	uvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](0));
	uvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](1));
	uvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](2));
	uncorrectedFlag_v.push_back(vmsData_p->v_flag.at(msRowReIndex_v[iData]));
	uncorrectedWeight_v.push_back(weight_v.at(msRowReIndex_v[iData]));
	uncorrectedSigma_v.push_back(sigma_v.at(msRowReIndex_v[iData]));
      }
    }

    // Have we asked to write an MS with corrected data + radiometric data ?
    
    // Are we with radiometric data ? Then we assume that the data are labelled AP_UNCORRECTED.
    if (sdmBinData.processorType(r_p) == ProcessorTypeMod::RADIOMETER) {
      if ((iter=vmsData_p->v_m_data.at(msRowReIndex_v[iData]).find(AtmPhaseCorrectionMod::AP_UNCORRECTED)) != vmsData_p->v_m_data.at(msRowReIndex_v[iData]).end()){
	
	correctedTime_v.push_back(vmsData_p->v_time.at(msRowReIndex_v[iData]));
	correctedAntennaId1_v.push_back(vmsData_p->v_antennaId1.at(msRowReIndex_v[iData]));
	correctedAntennaId2_v.push_back(vmsData_p->v_antennaId2.at(msRowReIndex_v[iData]));
	correctedFeedId1_v.push_back(vmsData_p->v_feedId1.at(msRowReIndex_v[iData]));
	correctedFeedId2_v.push_back(vmsData_p->v_feedId2.at(msRowReIndex_v[iData]));
	correctedFilteredDD_v.push_back(filteredDD.at(msRowReIndex_v[iData]));
	correctedFilteredShape_vv.push_back(filteredShape_vv.at(msRowReIndex_v[iData]));
	correctedFieldId_v.push_back(vmsData_p->v_fieldId.at(msRowReIndex_v[iData]));
	correctedInterval_v.push_back(vmsData_p->v_interval.at(msRowReIndex_v[iData]));
	correctedExposure_v.push_back(vmsData_p->v_exposure.at(msRowReIndex_v[iData]));
	correctedTimeCentroid_v.push_back(vmsData_p->v_timeCentroid.at(msRowReIndex_v[iData]));
	correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](0));
	correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](1));
	correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](2));
	// this is probably the most recent uncorrectedData, but it might not be - else why the if blocks here
	correctedData_v.push_back(cdf.to4Pol(vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(0),
					     vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(1),
					     iter->second));
	correctedFlag_v.push_back(vmsData_p->v_flag.at(msRowReIndex_v[iData]));
	correctedWeight_v.push_back(weight_v.at(msRowReIndex_v[iData]));
	correctedSigma_v.push_back(sigma_v.at(msRowReIndex_v[iData]));
      }
    }
    else {  // We assume that we are in front of CORRELATOR data, but do we have corrected data on that specific subscan ?
      if (subscanHasCorrectedData) {
	// Then we know that we have AP_CORRECTED data.
	if  (vmsData_p->v_antennaId1.at(msRowReIndex_v[iData]) == vmsData_p->v_antennaId2.at(msRowReIndex_v[iData]) ) {
	  /*
	  ** do not forget to prepend the autodata copied from the uncorrected data, because the lower layers of the software do not put the (uncorrected) autodata in the
	  ** corrected data.
	  */
	  correctedTime_v.push_back(vmsData_p->v_time.at(msRowReIndex_v[iData]));
	  correctedAntennaId1_v.push_back(vmsData_p->v_antennaId1.at(msRowReIndex_v[iData]));
	  correctedAntennaId2_v.push_back(vmsData_p->v_antennaId2.at(msRowReIndex_v[iData]));
	  correctedFeedId1_v.push_back(vmsData_p->v_feedId1.at(msRowReIndex_v[iData]));
	  correctedFeedId2_v.push_back(vmsData_p->v_feedId2.at(msRowReIndex_v[iData]));
	  correctedFilteredDD_v.push_back(filteredDD.at(msRowReIndex_v[iData]));
	  correctedFieldId_v.push_back(vmsData_p->v_fieldId.at(msRowReIndex_v[iData]));
	  correctedInterval_v.push_back(vmsData_p->v_interval.at(msRowReIndex_v[iData]));
	  correctedExposure_v.push_back(vmsData_p->v_exposure.at(msRowReIndex_v[iData]));
	  correctedTimeCentroid_v.push_back(vmsData_p->v_timeCentroid.at(msRowReIndex_v[iData]));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](0));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](1));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](2));	  
	  // this is probably the most recent uncorrectedData, but it might not be - else why the if blocks here
	  correctedData_v.push_back(cdf.to4Pol(vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(0),
					       vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(1),
					       iter->second));
	  correctedFlag_v.push_back(vmsData_p->v_flag.at(msRowReIndex_v[iData]));
	  correctedFilteredShape_vv.push_back(filteredShape_vv.at(msRowReIndex_v[iData]));
	  correctedWeight_v.push_back(weight_v.at(msRowReIndex_v[iData]));
	  correctedSigma_v.push_back(sigma_v.at(msRowReIndex_v[iData]));

	
	}
	else {
	  /*
	  ** And now finally the correlation corrected data.
	  */
	  correctedTime_v.push_back(vmsData_p->v_time.at(msRowReIndex_v[iData]));
	  correctedAntennaId1_v.push_back(vmsData_p->v_antennaId1.at(msRowReIndex_v[iData]));
	  correctedAntennaId2_v.push_back(vmsData_p->v_antennaId2.at(msRowReIndex_v[iData]));
	  correctedFeedId1_v.push_back(vmsData_p->v_feedId1.at(msRowReIndex_v[iData]));
	  correctedFeedId2_v.push_back(vmsData_p->v_feedId2.at(msRowReIndex_v[iData]));
	  correctedFilteredDD_v.push_back(filteredDD.at(msRowReIndex_v[iData]));
	  correctedFieldId_v.push_back(vmsData_p->v_fieldId.at(msRowReIndex_v[iData]));
	  correctedInterval_v.push_back(vmsData_p->v_interval.at(msRowReIndex_v[iData]));
	  correctedExposure_v.push_back(vmsData_p->v_exposure.at(msRowReIndex_v[iData]));
	  correctedTimeCentroid_v.push_back(vmsData_p->v_timeCentroid.at(msRowReIndex_v[iData]));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](0));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](1));
	  correctedUvw_v.push_back(vv_uvw[msRowReIndex_v[iData]](2));
	  iter=vmsData_p->v_m_data.at(msRowReIndex_v[iData]).find(AtmPhaseCorrectionMod::AP_CORRECTED);
	  float* theData = cdf.to4Pol(vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(0),
				      vmsData_p->vv_dataShape.at(msRowReIndex_v[iData]).at(1),
				      iter->second);
	  correctedData_v.push_back(theData);
	  correctedFlag_v.push_back(vmsData_p->v_flag.at(msRowReIndex_v[iData]));
	  correctedFilteredShape_vv.push_back(filteredShape_vv.at(msRowReIndex_v[iData]));
	  correctedWeight_v.push_back(weight_v.at(msRowReIndex_v[iData]));
	  correctedSigma_v.push_back(sigma_v.at(msRowReIndex_v[iData]));
	}
      }
    }
  }
 
  if (uncorrectedData_v.size() > 0 && (msFillers.find(AtmPhaseCorrectionMod::AP_UNCORRECTED) != msFillers.end())) {
    if (! mute) { // Here we make the assumption that we have always uncorrected data. This realistic even if not totally rigorous.

      if (!skipFirstTime) {
	// need to set these now, but they can be set as vectors that need to be reordered
	uncorrectedTime_v		  = reorder<double>( vmsData_p->v_time, msRowReIndex_v);
	uncorrectedAntennaId1_v	  = reorder<int> (vmsData_p->v_antennaId1, msRowReIndex_v ) ;
	uncorrectedAntennaId2_v	  = reorder<int>(vmsData_p->v_antennaId2, msRowReIndex_v );
	uncorrectedFeedId1_v	  = reorder<int>(vmsData_p->v_feedId1, msRowReIndex_v );
	uncorrectedFeedId2_v	  = reorder<int>(vmsData_p->v_feedId2, msRowReIndex_v );
	uncorrectedFieldId_v	  = reorder<int>(vmsData_p->v_fieldId, msRowReIndex_v );
	uncorrectedFilteredDD_v	  = reorder<int>(filteredDD, msRowReIndex_v );
	uncorrectedFilteredShape_vv = reorder<vector< unsigned int> >(filteredShape_vv, msRowReIndex_v );
	uncorrectedInterval_v	  = reorder<double>(vmsData_p->v_interval, msRowReIndex_v );
	uncorrectedExposure_v	  = reorder<double>(vmsData_p->v_exposure, msRowReIndex_v );
	uncorrectedTimeCentroid_v	  = reorder<double>(vmsData_p->v_timeCentroid, msRowReIndex_v );
	uncorrectedFlag_v		  = reorder<unsigned int>(vmsData_p->v_flag, msRowReIndex_v );
	// uvw_v, uncorrectedSigma_v, and uncorrectedWeight_v have all been previously set
      }

      // Here we make the assumption that the State is the same for all the antennas and let's use the first State found in the vector stateId contained in the ASDM Main Row
      // state must have already been filled so that the stateIdx2IDx map is available to extract the state ID for this main row pointer (r_p).
      vector<int> msStateId_v(uncorrectedTime_v.size(), stateIdx2Idx[r_p]);

      msFillers[AtmPhaseCorrectionMod::AP_UNCORRECTED]->addData(complexData,
					 uncorrectedTime_v	, // this is already time midpoint
					 uncorrectedAntennaId1_v,
					 uncorrectedAntennaId2_v,
					 uncorrectedFeedId1_v,
					 uncorrectedFeedId2_v,
					 uncorrectedFilteredDD_v,
					 (int) vmsData_p->processorId,
					 uncorrectedFieldId_v,
					 uncorrectedInterval_v,
					 uncorrectedExposure_v,
					 uncorrectedTimeCentroid_v,
					 (int) r_p->getScanNumber(), 
					 0,                                               // Array Id
					 (int) r_p->getExecBlockId().getTagValue(), // Observation Id
					 msStateId_v,
					 uvw_v,               
					 uncorrectedFilteredShape_vv, // vmsData_p->vv_dataShape after filtering the case numCorr == 3
					 uncorrectedData_v,
					 uncorrectedFlag_v,
					 uncorrectedWeight_v,
					 uncorrectedSigma_v);
    }
  }


  if (correctedData_v.size() > 0 && (msFillers.find(AtmPhaseCorrectionMod::AP_CORRECTED) != msFillers.end())) {
    if (! mute) {
      // Here we make the assumption that the State is the same for all the antennas and let's use the first State found in the vector stateId contained in the ASDM Main Row
      // state must have already been filled so that the stateIdx2IDx map is available to extract the state ID for this main row pointer (r_p).
      vector<int>  correctedMsStateId_v(correctedTime_v.size(), stateIdx2Idx[r_p]);
      msFillers[AtmPhaseCorrectionMod::AP_CORRECTED]->addData(complexData,
				       correctedTime_v, // this is already time midpoint
				       correctedAntennaId1_v, 
				       correctedAntennaId2_v,
				       correctedFeedId1_v,
				       correctedFeedId2_v,
				       correctedFilteredDD_v,
				       vmsData_p->processorId,
				       correctedFieldId_v,
				       correctedInterval_v,
				       correctedExposure_v,
				       correctedTimeCentroid_v,
				       (int) r_p->getScanNumber(), 
				       0,                                               // Array Id
				       (int) r_p->getExecBlockId().getTagValue(), // Observation Id
				       correctedMsStateId_v,
				       correctedUvw_v,
				       correctedFilteredShape_vv, // vmsData_p->vv_dataShape after filtering the case numCorr == 3
				       correctedData_v,
				       correctedFlag_v,
				       correctedWeight_v,
				       correctedSigma_v);
    }
  }
  if (debug) cout << "fillMain : exiting" << endl;
}

void fillSysPower_aux (const vector<SysPowerRow *>& sysPowers, map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>& msFillers_m) {
   LOGENTER("fillSysPower_aux");

  vector<int>		antennaId;
  vector<int>		spectralWindowId;
  vector<int>		feedId;
  vector<double>	time;
  vector<double>	interval;
  vector<int>		numReceptor;
  vector<float>		switchedPowerDifference;
  vector<float>		switchedPowerSum;
  vector<float>		requantizerGain;

  LOG("fillSysPower_aux : resizing the arrays (" + TO_STRING(sysPowers.size()) + ") to populate the columns of the MS SYSPOWER table.");

  antennaId.resize(sysPowers.size());
  spectralWindowId.resize(sysPowers.size());
  feedId.resize(sysPowers.size());
  time.resize(sysPowers.size());
  interval.resize(sysPowers.size());
  
  /*
   * Prepare the mandatory attributes.
   */
  LOG("fillSysPower_aux : filling the arrays to populate the columns of the MS SYSPOWER table.");

  transform(sysPowers.begin(), sysPowers.end(), antennaId.begin(), sysPowerAntennaId);
  transform(sysPowers.begin(), sysPowers.end(), spectralWindowId.begin(), sysPowerSpectralWindowId);
  transform(sysPowers.begin(), sysPowers.end(), feedId.begin(), sysPowerFeedId);
  transform(sysPowers.begin(), sysPowers.end(), time.begin(), sysPowerMidTimeInSeconds);
  transform(sysPowers.begin(), sysPowers.end(), interval.begin(), sysPowerIntervalInSeconds);
  
  /*
   * Prepare the optional attributes.
   */
  LOG("fillSysPower_aux : working on the optional attributes.");

  unsigned int numReceptor0 = sysPowers[0]->getNumReceptor();
  LOG("fillSysPower_aux : numReceptor = " + TO_STRING(numReceptor0));
 
  bool switchedPowerDifferenceExists0 = sysPowers[0]->isSwitchedPowerDifferenceExists();
  if (switchedPowerDifferenceExists0) {
    switchedPowerDifference.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerSwitchedPowerDifference(switchedPowerDifference.begin()));
  }

  bool switchedPowerSumExists0 = sysPowers[0]->isSwitchedPowerSumExists();
  if (switchedPowerSumExists0) {
    switchedPowerSum.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerSwitchedPowerSum(switchedPowerSum.begin()));
  }
  
  bool requantizerGainExists0 = sysPowers[0]->isRequantizerGainExists();
  if (requantizerGainExists0) {
    requantizerGain.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerRequantizerGain(requantizerGain.begin()));  
  }

  LOG("fillSysPower_aux : about to append a slice to the MS SYSPOWER table.");
 
  for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator msIter = msFillers_m.begin();
       msIter != msFillers_m.end();
       ++msIter) {
    msIter->second->addSysPowerSlice(antennaId.size(),
				     antennaId,
				     spectralWindowId,
				     feedId,
				     time,
				     interval,
				     (unsigned int) numReceptor0,
				     switchedPowerDifference,
				     switchedPowerSum,
				     requantizerGain);
  }
  infostream << "Appended " << sysPowers.size() << " rows to the MS SYSPOWER table." << endl;
  LOGEXIT("fillSysPower_aux");
}

/**
 * This function fills the MS SysPower table from an ASDM SysPower table.
 *
 * @param ds the ASDM dataset the ASDM SysPower table belongs to.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
 * @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
 * @param msFillers_m a map of ASDM2MSFillers depending on AtmosphericPhaseCorrection.
 *
 */
void fillSysPower(const string asdmDirectory, ASDM* ds_p, bool ignoreTime, const vector<ScanRow *>& selectedScanRow_v, map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>& msFillers_m) {
  LOGENTER("fillSysPower");

  const SysPowerTable& sysPowerT = ds_p->getSysPower();

  infostream.str("");
  infostream << "The dataset has " << sysPowerT.size() << " syspower(s).";
  info(infostream.str()); 
  
  if (sysPowerT.size() > 0 ) {
    try {
      // Prepare a row filter based on the time intervals of the selected scans.
      rowsInAScanbyTimeIntervalFunctor<SysPowerRow> selector(selectedScanRow_v);

      //
      // We can assume that there is an SysPower table , but we don't know yet if it's stored in a binary or an XML file.
      // 
      if (file_exists(uniqSlashes(asdmDirectory + "/SysPower.bin"))) {

	LOG("fillSysPower : working with SysPower.bin by successive slices.");

	TableStreamReader<SysPowerTable, SysPowerRow> tsrSysPower;
	tsrSysPower.open(asdmDirectory);

	// We can process the SysPower table by slice when it's stored in a binary file so let's do it.
	while (tsrSysPower.hasRows()) {
	  const vector<SysPowerRow*>&	sysPowerRows = tsrSysPower.untilNBytes(50000000);
	  infostream.str("");
	  infostream << "(considering the next " << sysPowerRows.size() << " rows of the SysPower table. ";

	  LOG("fillSysPower : determining which rows are in the selected scans.");
	  const vector<SysPowerRow *>& sysPowers = selector(sysPowerRows, ignoreTime);
  
	  if (!ignoreTime) 
	    infostream << sysPowers.size() << " of them are in the selected exec blocks / scans";
  
	  infostream << ")";	     
	  info(infostream.str());
	  
	  infostream.str("");
	  errstream.str("");
  
	  if (sysPowers.size() > 0)
	    fillSysPower_aux(sysPowerRows, msFillers_m);
	}
	tsrSysPower.close();
      }
      
      else if (file_exists(uniqSlashes(asdmDirectory + "/SysPower.xml"))) {

	LOG("fillSysPower : working with SysPower.xml read with a TableSAXReader");

	//
	// Instantiate a TableSAXReader functor with T==SysPowerTable, R==SysPowerRow 
	// and RFilter==rowsInAScanbyTimeIntervalFunctor<SysPowerRow>.
	//
	TableSAXReader<SysPowerTable, SysPowerRow, rowsInAScanbyTimeIntervalFunctor<SysPowerRow> >
	  tableSAXReader(verbose,
			 selector,
			 &fillSysPower_aux,
			 msFillers_m);
	
	// Execute the functor
	tableSAXReader(asdmDirectory, ignoreTime);

      }
      else 
	throw ConversionException ("fillSysPower: no file found for SysPower", "SysPower");

      unsigned int numMSSysPowers =  (const_cast<casacore::MeasurementSet*>(msFillers_m.begin()->second->ms()))->rwKeywordSet().asTable("SYSPOWER").nrow();
      if (numMSSysPowers > 0) {
	infostream.str("");
	infostream << "converted in " << numMSSysPowers << " syspower(s) in the measurement set.";
	info(infostream.str());
      }
      
    } // end of filling SysPower by slice.
    
    catch (ConversionException e) {
      errstream.str("");
      errstream << e.getMessage();
      error(errstream.str());
    }
    catch ( std::exception & e) {
      errstream.str("");
      errstream << e.what();
      error(errstream.str());      
    }
    
  }
  LOGEXIT("fillSysPower");
}

void fillPointingRows(const vector<PointingRow *> &vpr, bool withPointingCorrection, map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller *> &msFillers_m, map<Tag, double> &lastTime_m) {
    LOGENTER("fillPointingRows");

    int nPointing = vpr.size();

    // Check some assertions.
    //
    // All rows of ASDM-Pointing must have their attribute usePolynomials equal to false
    // and their numTerm attribute equal to 1. Use the opportunity of this check
    // to compute the number of rows to be created in the MS-Pointing by summing
    // all the numSample attributes values. Watch for duplicate times due and avoid.
    //
    int numMSPointingRows = 0;

    // set this to true for any rows where the first element should be skipped because the time is a duplicate
    vector<bool> vSkipFirst(vpr.size(),false);

    for (unsigned int i = 0; i < vpr.size(); i++) {
        if (vpr[i]->getUsePolynomials()) {
            errstream.str("");
            errstream << "Found usePolynomials equal to true at row #" << i <<". Can't go further.";
            error(errstream.str());
        }
        
        // look for duplicate rows - time at the start of row is near the time at the end of the previous row for that antennaId
        int numSample = vpr[i]->getNumSample();
        Tag antId = vpr[i]->getAntennaId();
        double tfirst = 0.0;
        double tlast = 0.0;
        if (vpr[i]->isSampledTimeIntervalExists()) {
            tfirst = ((double) (vpr[i]->getSampledTimeInterval().at(0).getStart().get())) / ArrayTime::unitsInASecond;
            tlast = ((double) (vpr[i]->getSampledTimeInterval().at(numSample-1).getStart().get())) / ArrayTime::unitsInASecond;
        } else {
            double tstart = ((double) vpr[i]->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond;
            double tint = ((double) vpr[i]->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond / numSample;
            tfirst = tstart + tint/2.0;
            tlast = tfirst + tint*(numSample-1);
        }
        if (casacore::near(tfirst,lastTime_m[antId])) {
            infostream.str("");
            infostream << "First time in Pointing row " << i << " for antenna " << antId << " is near the previous last time for that antenna, skipping that first time in the output MS POINTING table" << endl;
            info(infostream.str());
            numSample--;
            lastTime_m[antId] = tlast;
            vSkipFirst[i] = true;
        }
        lastTime_m[antId] = tlast;
        numMSPointingRows += numSample;
    }

    //
    // Ok now we have verified the assertions and we know the number of rows
    // to be created into the MS-Pointing, we can proceed.

    PointingRow* r = 0;

    vector<int>	antenna_id_(numMSPointingRows, 0);
    vector<double>	time_(numMSPointingRows, 0.0);
    vector<double>	interval_(numMSPointingRows, 0.0);
    vector<double>	direction_(2 * numMSPointingRows, 0.0);
    vector<double>	target_(2 * numMSPointingRows, 0.0);
    vector<double>	pointing_offset_(2 * numMSPointingRows, 0.0);
    vector<double>	encoder_(2 * numMSPointingRows, 0.0);
    vector<bool>	tracking_(numMSPointingRows, false);

    //
    // Let's check if the optional attribute overTheTop is present somewhere in the table.
    //
    unsigned int numOverTheTop = count_if(vpr.begin(), vpr.end(), overTheTopExists);
    bool overTheTopExists4All = vpr.size() == numOverTheTop;

    vector<bool> v_overTheTop_ ;

    vector<s_overTheTop> v_s_overTheTop_;
    
    if (overTheTopExists4All) 
        v_overTheTop_.resize(numMSPointingRows);
    else if (numOverTheTop > 0) 
        v_overTheTop_.resize(numOverTheTop);

    int iMSPointingRow = 0;
    for (int i = 0; i < nPointing; i++) {     // Each row in the ASDM-Pointing ...
        r = vpr.at(i);

        // Let's prepare some values.
        int antennaId = r->getAntennaId().getTagValue();
      
        double time = 0.0, interval = 0.0;
        if (!r->isSampledTimeIntervalExists()) { // If no sampledTimeInterval then
            // then compute the first value of MS TIME and INTERVAL.
            interval   = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond / r->getNumSample();
            // if (isEVLA) {
            //   time = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond;
            // }
            // else {
            time = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0;
            //}
        }

        //
        // The size of each vector below 
        // should be checked against numSample !!!
        //
        int numSample = r->getNumSample();
        const vector<vector<Angle> > encoder = r->getEncoder();
        checkVectorSize<vector<Angle> >("encoder", encoder, "numSample", (unsigned int) numSample, "Pointing", (unsigned int)i);
        
        const vector<vector<Angle> > pointingDirection = r->getPointingDirection();
        checkVectorSize<vector<Angle> >("pointingDirection", pointingDirection, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);
        
        const vector<vector<Angle> > target = r->getTarget();
        checkVectorSize<vector<Angle> >("target", target, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);

        const vector<vector<Angle> > offset = r->getOffset();
        checkVectorSize<vector<Angle> >("offset", offset, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);

        bool   pointingTracking = r->getPointingTracking();
 
        //
        // Prepare some data structures and values required to compute the
        // (MS) direction.
        vector<double> cartesian1(3, 0.0);
        vector<double> cartesian2(3, 0.0);
        vector<double> spherical1(2, 0.0);
        vector<double> spherical2(2, 0.0);
        vector<vector<double> > matrix3x3;
        for (unsigned int ii = 0; ii < 3; ii++) {
            matrix3x3.push_back(cartesian1); // cartesian1 is used here just as a way to get a 1D vector of size 3.
        }
        double PSI = M_PI_2;
        double THETA;
        double PHI;
      
        vector<ArrayTimeInterval> timeInterval ;
        if (r->isSampledTimeIntervalExists()) timeInterval = r->getSampledTimeInterval();

        // and now insert the values into the MS POINTING table
        // watch for the skipped initial sample
        int numSampleUsed = numSample;
        if (vSkipFirst.at(i)) numSampleUsed--;
        
        // Use 'fill' from algorithm for the cases where values remain constant.
        // ANTENNA_ID
        fill(antenna_id_.begin()+iMSPointingRow, antenna_id_.begin()+iMSPointingRow+numSampleUsed, antennaId);

        // TRACKING 
        fill(tracking_.begin()+iMSPointingRow, tracking_.begin()+iMSPointingRow+numSampleUsed, pointingTracking);

        // OVER_THE_TOP 
        if (overTheTopExists4All)
            // it's present everywhere
            fill(v_overTheTop_.begin()+iMSPointingRow, v_overTheTop_.begin()+iMSPointingRow+numSampleUsed,
                 r->getOverTheTop());
        else if (r->isOverTheTopExists()) {
            // it's present only in some rows.
            s_overTheTop saux ;
            saux.start = iMSPointingRow; saux.len = numSampleUsed; saux.value = r->getOverTheTop();
            v_s_overTheTop_.push_back(saux);
        }
       
        // Use an explicit loop for the other values.
        for (int j = 0 ; j < numSample; j++) { // ... must be expanded in numSample MS-Pointing rows.
            if (j == 0 && vSkipFirst.at(i)) continue;   // the first element is to be skipped
                
            // TIME and INTERVAL
            if (r->isSampledTimeIntervalExists()) { //if sampledTimeInterval is present use its values.	           
                // Here the size of timeInterval will have to be checked against numSample !!
                interval_[iMSPointingRow] = ((double) timeInterval.at(j).getDuration().get()) / ArrayTime::unitsInASecond ;
                time_[iMSPointingRow] = ((double) timeInterval.at(j).getStart().get()) / ArrayTime::unitsInASecond
                    + interval_[iMSPointingRow]/2;	  
            }
            else {                                     // otherwise compute TIMEs and INTERVALs from the first values.
                interval_[iMSPointingRow]            = interval;
                time_[iMSPointingRow]                = time + j*interval;
            }

            // DIRECTION
            THETA			      = target.at(j).at(1).get();
            PHI				      = -M_PI_2 - target.at(j).at(0).get();
            spherical1[0]		      = offset.at(j).at(0).get();
            spherical1[1]		      = offset.at(j).at(1).get();
            rect(spherical1, cartesian1);
            eulmat(PSI, THETA, PHI, matrix3x3);
            matvec(matrix3x3, cartesian1, cartesian2);
            spher(cartesian2, spherical2);
            direction_[2*iMSPointingRow]      = spherical2[0] ;
            direction_[2*iMSPointingRow+1]    = spherical2[1] ;
            if (withPointingCorrection) { // Cf CSV-2878 and ICT-1532
                direction_[2*iMSPointingRow]   += encoder.at(j).at(0).get() - pointingDirection.at(j).at(0).get();
                direction_[2*iMSPointingRow+1] += encoder.at(j).at(1).get() - pointingDirection.at(j).at(1).get() ;
            }

            // TARGET
            target_[2*iMSPointingRow]     = target.at(j).at(0).get();
            target_[2*iMSPointingRow+1]   = target.at(j).at(1).get();
                
            // POINTING_OFFSET
            pointing_offset_[2*iMSPointingRow]   = offset.at(j).at(0).get();
            pointing_offset_[2*iMSPointingRow+1] = offset.at(j).at(1).get();

            // ENCODER
            encoder_[2*iMSPointingRow]           = encoder.at(j).at(0).get();
            encoder_[2*iMSPointingRow+1]         = encoder.at(j).at(1).get();

            // increment the row number in MS Pointing.
            iMSPointingRow++;	
        }
    }
    
    
    for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers_m.begin();
         iter != msFillers_m.end();
         ++iter) {
        iter->second->addPointingSlice(numMSPointingRows,
                                       antenna_id_,
                                       time_,
                                       interval_,
                                       direction_,
                                       target_,
                                       pointing_offset_,
                                       encoder_,
                                       tracking_,
                                       overTheTopExists4All,
                                       v_overTheTop_,
                                       v_s_overTheTop_);
    }
    LOGEXIT("fillPointingRows");
}

/**
 * This function fills the MS POINTING table from an ASDM Pointing table. For a binary pointing table it streams the table to try and minimize memory use. 
 * 
 * @param ds : the ASDM dataset that the Pointing is found in.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all of the table is going to be processed.
 * @param withPointingCorrection add (ASDM::Pointing::encoder - ASDM::POINTING::pointingDirection) to the value going to MS::POINTING::DIRECTION
 * @param selectedScanRow_v is a vector of pointers on ScanRow used to determine which rows of Pointing are going to be processed.
 * @param msFillers_m a map of ASDM2MSFiller depending on AtmosphericPhaseCorrection
 * 
 */
void fillPointing(const string asdmDirectory, ASDM* ds_p, bool ignoreTime, bool withPointingCorrection, const vector<ScanRow *> &selectedScanRow_v, map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller *> &msFillers_m) {
    LOGENTER("fillPointing");
    
    const PointingTable& pointingT = ds_p->getPointing();
    infostream.str("");
    infostream << "The dataset has " << pointingT.size() << " pointing(s)...";
    info(infostream.str());

    if (pointingT.size() > 0) {

	// initialize the lastTime to 0.0 for all antennaIds
	map<Tag, double> lastTime;
	const AntennaTable& antennaT = ds_p->getAntenna();
	const vector<AntennaRow *>& vAntRow = antennaT.get();
	for (unsigned int i=0; i < vAntRow.size(); i++) {
            lastTime[vAntRow[i]->getAntennaId()] = 0.0;
	}
        // Prepare a row filter based on the time intervals of the selected scans.
        rowsInAScanbyTimeIntervalFunctor<PointingRow> selector(selectedScanRow_v);

        if (file_exists(uniqSlashes(asdmDirectory + "/Pointing.bin"))) {
            LOG("fillPointing : working with Pointing.bin by successive slices.");

            TableStreamReader<PointingTable, PointingRow> tsrPointing;
            tsrPointing.open(asdmDirectory);

            // we can process the Pointing table by slice when it's stored in a binary file so let's do it.
            while (tsrPointing.hasRows()) {
                // arbitrarily uses the same limit here that's used in SysPower
                const vector<PointingRow *>& pointingRowsSlice = tsrPointing.untilNBytes(50000000);
                infostream.str("");
                infostream << "(considering the next " << pointingRowsSlice.size() << " rows of the Pointing table. ";
                LOG ("fillPointing : determining which rows are in the selected scans.");
                const vector<PointingRow *>& pointingRows = selector(pointingRowsSlice, ignoreTime);
                if (!ignoreTime)
                    infostream << pointingRows.size() << " of them are in the selected exec blocks / scans";

                infostream << ")";
                info(infostream.str());

                infostream.str("");
                errstream.str("");

                if (pointingRows.size() > 0) 
                    fillPointingRows(pointingRows, withPointingCorrection, msFillers_m, lastTime);
            }
            tsrPointing.close();
        } else {
            // process the full Pointing table, presumably this is in Pointing.xml
            LOG("fillPointing : determining which rows are in the selected scans.");
            const vector<PointingRow *> &pointingRows = selector(pointingT.get(), ignoreTime);
            if (!ignoreTime)  {
                infostream.str("");
                infostream << pointingRows.size() << " of them in the selected exec blocks / scans ... ";
                info(infostream.str());
            }
            fillPointingRows(pointingRows, withPointingCorrection, msFillers_m, lastTime);
        }
    }
    
    unsigned int numMSPointings = msFillers_m.begin()->second->ms()->pointing().nrow();
    if (numMSPointings) {
        infostream.str("");
        infostream << "converted in " << numMSPointings << " pointing(s) in the measurement set." ;
        info(infostream.str()); 
    }
    LOGEXIT("fillPointing");
}

    

class MSMainRowsInSubscanChecker {
public:
  MSMainRowsInSubscanChecker();
  virtual ~MSMainRowsInSubscanChecker();
  void check(const VMSData* vmsData_p, MainRow* mainRow_p, unsigned int mainRowIndex, const string& BDFName);
  const vector<string>& report() const;
  void reset();

private:
  vector<string> report_v;
};

MSMainRowsInSubscanChecker::MSMainRowsInSubscanChecker() {;}
MSMainRowsInSubscanChecker::~MSMainRowsInSubscanChecker() {;}
void MSMainRowsInSubscanChecker::reset() {
  LOGENTER("MSMainRowsInSubscanChecker::reset");
  report_v.clear();
  LOGEXIT("MSMainRowsInSubscanChecker::reset");
}

void MSMainRowsInSubscanChecker::check( const VMSData* vmsData_p,
					MainRow* mainRow_p,
					unsigned int mainRowIndex,
					const string& BDFName ) {
  LOGENTER("MSMainRowsInSubscanChecker::check");
  SubscanTable & subscanTable = mainRow_p->getTable().getContainer().getSubscan();

  SubscanRow* subscanRow_p = subscanTable.getRowByKey(mainRow_p->getExecBlockId(),
						      mainRow_p->getScanNumber(),
						      mainRow_p->getSubscanNumber());
  if (subscanRow_p == NULL) {
    infostream.str("");
    infostream << "Could not find a row in the subscan table with the key 'execBlockId = "<< mainRow_p->getExecBlockId()
	       << ", scanNumber = " << mainRow_p->getScanNumber()
	       << ", subscanNumber = " << mainRow_p->getSubscanNumber()
	       << "'. I can't check if the BDF contents is in the subscan's time range.";
    info(infostream.str());
    LOGEXIT("MSMainRowsInSubscanChecker::check");
    return;
  }

  // We make the assumption that the content pointed by vmsData_p is ordered by time.
  double subscanStartTime = subscanRow_p->getStartTime().getMJD()*86400.0;
  double subscanEndTime   = subscanRow_p->getEndTime().getMJD()*86400.0;

  //
  // Now detect one of two abnormal situations : the 1st data time is anterior to the subscan start time or the last data time
  // is posterior to the subscan end time. 
  if ( (vmsData_p->v_time[0] < subscanStartTime) || (subscanEndTime < vmsData_p->v_time[vmsData_p->v_time.size() - 1])) {
    ostringstream oss;
    oss << "Main row #" << mainRowIndex
	<< " - The BDF '" << BDFName << "' contained data not in the time range of scan=" << mainRow_p->getScanNumber()
	<< ", subscan=" << mainRow_p->getSubscanNumber() << ".";
    string s = oss.str();
    if (!(report_v.size() > 0 && s == report_v.back()))
      report_v.push_back(s);
  }
  LOGEXIT("MSMainRowsInSubscanChecker::check");
}

const vector<string>& MSMainRowsInSubscanChecker::report() const {
  return report_v;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * The main function.
 */
int main(int argc, char *argv[]) {

  string dsName;
  string msNamePrefix;
  string msNameExtension;

  appName = string(argv[0]);

  ofstream ofs;

  //LogSinkInterface& lsif = LogSink::globalSink();
  static_cast<void>(LogSink::globalSink());

  uint64_t bdfSliceSizeInMb = 0; // The default size of the BDF slice hold in memory.

  bool mute = false;

  bool ac_xc_per_timestamp = false; // for the time being the option is 'preserve the old order'

  bool		interpolate_ephemeris	       = false; 
  bool		tabulate_ephemeris_polynomials = false;
  double	polyephem_tabtimestep	       = 0.001;
  bool          checkRowUniqueness = false; 
  string        scansOptionInfo;
  string        asisOption;
  bool	ignoreTime = false;
  bool	processSysPower = true;
  bool	processCalDevice = true;
  bool  processPointing	= true;
  bool  withPointingCorrection = false;
  bool  processEphemeris = true;
  bool checkdupints = true;
  
  //   Process command line options and parameters.

  // all of the non-positional options need to be enumerated here
  // note that asdm-directory and ms-directory-prefix can be specified as both
  // positional arguments and as named arguments. The positional argument takes
  // precedence.

  enum optionIndex { UNKNOWN, HELP, ICM, ISRT, ITS, OCM, COMPRESSION, 
		     ASIS, WVRCORRDATA, SCANS, LOGFILE, VERBOSE, REVISION, 
		     DRYRUN, IGNORETIME, NOCALDEV, NOEPHEMERIS, NOSYSPOWER,
		     NOPOINTING, CHECKROWUNIQ, BDFSLICESIZE, LAZY,
		     WITHPCORR, ACXCPERTIME, POLYEPHTSTEP, INTEPHEM,
                     ASDMDIR, MSDIRPREFIX, CHECKDUPINTS };


  try {

    // remove the program name
    argc--;
    argv++;

    string usageIntro = 
      "Converts an ASDM dataset into a CASA measurement set.\n"
      "Usage : " + appName + " [options] asdm-directory [ms-directory-prefix]\n\n"
      "Command parameters: \n";

    // Descriptor elements are: OptionIndex, OptionType, shortopt, longopt, check_arg, help
    option::Descriptor usage[] = {
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  usageIntro.c_str()},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  " \tasdm-directory :  \tthe pathname to the ASDM dataset to be converted"},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  
	" \tms-directory-prefix :  \tthe prefix of the pathname(s) of the measurement "
	"set(s) to be created. This prefix is completed by a suffix to form the "
	"name(s) of the resulting measurement set(s). "
	"this suffix depends on the selected options (see options compression and wvr-corrected-data)."},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  "\nAllowed options:\n"},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  0 }, // helps with formatting

      // these are the actual options
      { HELP, 0, "", "help", AlmaArg::None, " --help  \tproduces this help message."},
      { ICM, 0, "", "icm",  AlmaArg::Required, 
	" --icm arg (=all) \tspecifies the correlation mode to be considered on input. "
	"A quoted string containing a sequence of 'ao' 'co' 'ac' 'all' separated by whitespaces is expected"},
      { ISRT, 0, "", "isrt", AlmaArg::Required, 
	" --isrt arg (=all) \tspecifies the spectral resolution type to be considered on input. "
	"A quoted string containing a sequence of 'fr' 'ca' 'bw' 'all' separated by whitespaces is expected"},
      { ITS, 0, "", "its",  AlmaArg::Required, 
	" --its arg (=all) \tspecifies the time sampling (INTEGRATION and/or SUBINTEGRATION)  to be considered on input. "
	"A quoted string containing a sequence of 'i' 'si' 'all' separated by whitespaces is expected"},  
      { OCM, 0, "", "ocm", AlmaArg::Required,
	" --ocm arg (=ca) \toutput data for correlation mode AUTO_ONLY (ao) or CROSS_ONLY (co) or CROSS_AND_AUTO (ca)"},
      { COMPRESSION, 0, "c", "compression", AlmaArg::None,
	" --c [--compression]  \tproduces compressed columns in the resulting measurement set "
	"(not set by default). When this option is selected the string '-compressed' is inserted "
	"in the pathname of the resulting measurement set."},
      { ASIS, 0, "", "asis", AlmaArg::Required, 
	" --asis arg \tcreates verbatim copies of the ASDM tables in the output measurement set. "
	"The value given to this option must be a quoted string containing a list of table names "
	"separated by space characters; the wildcard character '*' is allowed in table names."},
      { WVRCORRDATA, 0, "", "wvr-corrected-data", AlmaArg::Required,  
	" --wvr-corrected-data arg (=no) \tspecifies wich values are considered in the ASDM binary data "
	"to fill the DATA column in the MAIN table of the MS. Expected values for this option are "
	"'no' for the uncorrected data (this is the default), 'yes' for the corrected data and "
	"'both' for corrected and uncorrected data. In the latter case, two measurement sets are "
	"created, one containing the uncorrected data and the other one, whose name is suffixed "
	"by '-wvr-corrected', containing the corrected data."},
      { SCANS, 0, "s", "scans", AlmaArg::Required,
	" --s [--scans] arg \tprocesses only the scans specified in the option's value. "
	"This value is a semicolon separated list of scan specifications. A scan specification "
	"consists of an exec bock index followed by the character ':' followed by a comma separated "
	"list of scan indexes or scan index ranges. A scan index is relative to the exec block it "
	"belongs to. Scan indexes are 1-based while exec blocks's are 0-based. "
	"\"0:1\" or \"2:2~6\" or \"0:1,1:2~6,8;2:,3:24~30\" \"1,2\" are valid values for the option. "
	"\"3:\" alone will be interpreted as 'all the scans of the exec block#3'. A scan index or a "
	"scan index range not preceded by an exec block index will be interpreted as 'all the scans with "
	"such indexes in all the exec blocks'.  By default all the scans are considered."},
      { LOGFILE, 0, "l", "logfile", AlmaArg::Required, 
	" -l [--logfile] arg \tspecifies the log filename. If the option is not used then the "
	"logged informations are written to the standard error stream."},
      { VERBOSE, 0, "v", "verbose", AlmaArg::None, " -v [--verbose]  \tlogs numerous informations as the filler is working."},
      { REVISION, 0, "r", "revision", AlmaArg::None, " -r [--revision]  \tlogs information about the revision of this application."},
      { DRYRUN, 0, "m", "dry-run", AlmaArg::None, " -m [--dry-run]  \tdoes not fill the MS MAIN table."},
      { IGNORETIME, 0, "t", "ignore-time", AlmaArg::None, 
	" -t [--ignore-time]  \tall the rows of the tables Feed, History, Pointing, Source, "
	"SysCal, CalDevice, SysPower and Weather are processed independently of the time range of the "
	"selected exec block / scan."},
      { NOCALDEV, 0, "", "no-caldevice", AlmaArg::None, " --no-caldevice  \tThe CalDevice table will be ignored."},
      { NOEPHEMERIS, 0, "", "no-ephemeris", AlmaArg::None, " --no-ephemeris  \tThe ephemeris table will be ignored."},
      { NOSYSPOWER, 0, "", "no-syspower", AlmaArg::None, " --no-syspower  \tThe SysPower table will be  ignored."},
      { NOPOINTING, 0, "", "no-pointing", AlmaArg::None, " --no-pointing  \tThe Pointing table will be ignored."},
      { CHECKROWUNIQ, 0, "", "check-row-uniqueness", AlmaArg::None, " --check-row-uniqueness  \tThe row uniqueness constraint will be checked in the tables where it's defined"},
      { BDFSLICESIZE, 0, "", "bdf-slice-size", AlmaArg::Long,  
	" --bdf-slice-size arg (=500) \tThe maximum amount of memory expressed as an integer "
	"in units of megabytes (1024*1024) allocated for BDF data. The default is 500 (megabytes)"},
      { LAZY, 0, "", "lazy", AlmaArg::None, " --lazy  \tdefers the production of the observational data in the MS Main table (DATA column) - Purely experimental, don't use in production !"},
      { WITHPCORR, 0, "", "with-pointing-correction", AlmaArg::None, 
	" --with-pointing-correction  \tadd (ASDM::Pointing::encoder - ASDM::Pointing::pointingDirection) "
	"to the value to be written in MS::Pointing::direction - (related with JIRA tickets CSV-2878 and ICT-1532))"},
      { ACXCPERTIME, 0, "", "ac-xc-per-timestamp", AlmaArg::Required, 
	" --ac-xc-per-timestamp arg (=no) \tif set to yes, then the filler writes in that order autocorrelations "
	"and cross correlations rows for one given data description and timestamp. Otherwise auto "
	"correlations data are grouped for a sequence of time stamps and then come the cross correlations "
	"data for the same sequence of timestamps."},
      { POLYEPHTSTEP, 0, "", "polyephem-tabtimestep", AlmaArg::Float, 
	" --polyephem-tabtimestep arg (=0.001) \tDefines the time step used to tabulate the polynomials "
	"found in the columns 'dir', 'distance' and optionally 'radVel' of the ASDM Ephemeris table. "
	"The unit to express the time step is the day and the default value is 0.001. If 'radvel' "
	"is not present then the radial velocity will be obtained by tabulating the derivative of "
	"the polynomial found in 'distance'."},
      { INTEPHEM, 0, "", "interpolate-ephemeris", AlmaArg::Required, 
	" --interpolate-ephemeris arg (=no) \tif set to 'yes' then the filler will resample the sequence "
	"of times found in the ASDM Ephemeris table into an evenly spaced sequence of times on which "
	"the ephemeris paarameters will obtained by an interpolation of degree 1. Otherwise (!= 'yes') "
	"the ephemeris parameters will be copies of what's in the ASDM Ephemeris table on the same "
	"sequence of times"},
      // these can be set by the command line, but are not shown the user. A config file may set these to be used as a parameter.
      { ASDMDIR, 0, "", "asdm-directory", AlmaArg::Required, 0 },
      { MSDIRPREFIX, 0, "", "ms-directory-prefix", AlmaArg::Required, 0},
      // checkdupints is used by the unit tests to turn off checks for duplicate integration times in RADIOMETER data, not intended for normal use
      { CHECKDUPINTS, 0, "", "checkdupints", AlmaArg::Bool, 0},
      { 0, 0, 0, 0, 0, 0 } };

    // Defaults are set by parsing an argv-like set of options where the values are the defaults
    const char *defaults[] = { "--icm=all",
			       "--isrt=all",
			       "--its=all",
			       "--ocm=ca",
			       "--wvr-corrected-data=no",
			       "--bdf-slice-size=500",
			       "--ac-xc-per-timestamp=no",
			       "--polyephem-tabtimestep=0.001",
			       "--interpolate-ephemeris=no",
			       "--checkdupints=true",
			       (const char *)-1};   // unambiguously signal the end

    // count defaults, more robust than setting a value here that must be changed when defaults changes
    int defaultCount = 0;
    while (defaults[defaultCount] != (const char *) -1) ++defaultCount;

    // parse defaults, argv
    // establish sizes
    option::Stats stats;
    // true here turns on re-ordering of args so that positional argument are always seen last
    stats.add(true, usage, defaultCount, defaults);
    stats.add(true, usage, argc, argv);

    // buffers to hold the parsed options
    // options has one element per optionIndex, last value is the last time it was set
    // buffer has one element for each option encountered, in order. Not used here.
    option::Option *options = new option::Option[stats.options_max];
    option::Option *buffer = new option::Option[stats.buffer_max];
    option::Parser parse;
    // parse the defaults first, then argv. User set options always come last
    // true here has same meaning as in stats above. This may not be necessary here, I think
    // the stats usage above has already reorderded argv in place.
    parse.parse(true, usage, defaultCount, defaults, options, buffer);
    parse.parse(true, usage, argc, argv, options, buffer);
    
    if (parse.error()) {
      errstream.str("");
      errstream << "Problem parsing the command line arguments";
      error(errstream.str());
    }
    
    // User-specified logfile?
    if (options[LOGFILE] != NULL && options[LOGFILE].last()->arg != NULL) {
#if 0     
      // Replaced the change at the LogSink level ...
      ofs.open(options[LOGFILE].last()->arg, ios_base::app);
      LogSinkInterface *theSink = new casacore::StreamLogSink(&ofs);
      LogSink::globalSink(theSink);
#else
      // ... with a change at the cerr (stderr) level since by default global logs are going to cerr (stderr).
      freopen(options[LOGFILE].last()->arg, "a", stderr);
#endif
    }

    // Help ? displays help's content and don't go further.

    if (options[HELP] || (argc==0) ) {
      errstream.str("");
      option::printUsage(errstream, usage,80);
      error(errstream.str());
    }

    // too many positional arguments?
    if (parse.nonOptionsCount() > 2) {
      errstream.str("");
      errstream << "Too many positional options" << endl;
      error(errstream.str());
    }

    // Verbose or quiet ?
    verbose = options[VERBOSE] != NULL;
   
    // Revision ? displays revision's info and don't go further if there is no dataset
    // to process otherwise proceed....
    string revision = "$Id: asdm2MS.cpp,v 1.84 2011/10/25 14:56:48 mcaillat Exp $\n";
    if (options[REVISION]) {
      if (options[ASDMDIR] || parse.nonOptionsCount() > 0) {
	infostream.str("");
	infostream << revision ;
	info(infostream.str());
      } else {
	errstream.str("");
	errstream << revision ;
	error(errstream.str());
      }
    }

    // set the non-string options with required values - always available because of defaults
    // just make sure the last one provided is always the one used
    {
      stringstream str(options[BDFSLICESIZE].last()->arg);
      str >> bdfSliceSizeInMb;
      if (!str) {
	// unlikely given that any value was already checked by AlmaArg::Long
	errstream.str("");
	errstream << "There was an error converting the bdf-slice-size value to an integer.";
	error(errstream.str());
      }
    }
    {
      stringstream str(options[POLYEPHTSTEP].last()->arg);
      str >> polyephem_tabtimestep;
      if (!str) {
	// unlikely given that any value was already checked by AlmaArg::Float
	errstream.str("");
	errstream << "There was an error converting the polyephem-tabtimestep value to a double";
	error(errstream.str());
      }
    }

    // non-default input data selection is incompatible with the lazy mode
    bool lazyModeOK = true;

    // this always has a value because of defaults
    string checkdupintsOpt(options[CHECKDUPINTS].last()->arg);
    trim(checkdupintsOpt);
    checkdupintsOpt = str_tolower(checkdupintsOpt);
    // AlmaArg::Bool already guarantees this is either "true" or "false"
    checkdupints = (checkdupintsOpt == "true");

    // Selection of correlation mode of data to be considered on input.
    istringstream iss;
    string token;

    // this always has a value because of defaults
    string icm_opt = string(options[ICM].last()->arg);
    iss.clear();
    iss.str(icm_opt);
    
    while (iss >> token) {
      if (token.compare("co") == 0) {
	es_cm.fromString("CROSS_ONLY", false);
	lazyModeOK = false;
      } else if (token.compare("ao") == 0) {
	es_cm.fromString("AUTO_ONLY", false);
	lazyModeOK = false;
      } else if (token.compare("ac") == 0) {
	es_cm.fromString("CROSS_AND_AUTO", false);
	lazyModeOK = false;
      } else if (token.compare("all") == 0) {
	es_cm.fromString("CROSS_ONLY AUTO_ONLY CROSS_AND_AUTO", false);
      } else {
	errstream.str("");
	errstream << "Token '" << token << "' invalid for --icm option." << endl;
	option::printUsage(errstream, usage);
	error(errstream.str());
      }
    }

    // Selection of spectral resolution type of data to be considered.

    // this always has a value because of defaults
    string isrt_opt = string(options[ISRT].last()->arg);
    iss.clear();
    iss.str(isrt_opt);

    while (iss >> token) {
      if (token.compare("fr") == 0) {
	lazyModeOK = false;
	es_srt.fromString("FULL_RESOLUTION", false);
      } else if (token.compare("ca") == 0) {
	lazyModeOK = false;
	es_srt.fromString("CHANNEL_AVERAGE", false);
      } else if (token.compare("bw") == 0) {
	lazyModeOK = false;
	es_srt.fromString("BASEBAND_WIDE", false);
      } else if (token.compare("all") == 0) {
	es_srt.fromString("FULL_RESOLUTION CHANNEL_AVERAGE BASEBAND_WIDE", false);
      } else { 
	errstream.str("");
	errstream << "Token '" << token << "' invalid for --isrt option." << endl;
	option::printUsage(errstream, usage);
	error(errstream.str());
      }
    }


    // Selection of the time sampling of data to be considered (integration and/or subintegration)
    // this always has a value because of defaults
    string its_opt = string(options[ITS].last()->arg);
    iss.clear();
    iss.str(its_opt);

    while ( iss >> token ) {
      if (token.compare("i") == 0) {
	lazyModeOK = false;
	es_ts.fromString("INTEGRATION",false);
      } else if (token.compare("si") == 0) {
	lazyModeOK = false;
	es_ts.fromString("SUBINTEGRATION", false);
      } else if (token.compare("all") == 0) {
	es_ts.fromString("INTEGRATION SUBINTEGRATION", false);
      } else {
	errstream.str("");
	errstream << "Token '" << token << "' invalid for its option." << endl;
	option::printUsage(errstream, usage);
	error(errstream.str());
      }
    }

    // Selection of the correlation mode of data to be produced in the measurement set.
    // this always has a value because of defaults
    string ocm_opt = string(options[OCM].last()->arg);
    if ( ocm_opt.compare("co") == 0 )
        e_query_cm = CorrelationModeMod::CROSS_ONLY;
    else if ( ocm_opt.compare("ao") == 0 )
        e_query_cm = CorrelationModeMod::AUTO_ONLY;
    else if ( ocm_opt.compare("ca") == 0 )
        e_query_cm = CorrelationModeMod::CROSS_AND_AUTO;
    else {
      errstream.str("");
      errstream << "Token '" << ocm_opt << "' invalid for ocm option." << endl;
      option::printUsage(errstream, usage);
      error(errstream.str());
    }

    if (parse.nonOptionsCount() > 0 || options[ASDMDIR]) {
      string dummy;
      if (parse.nonOptionsCount() > 0) {
	dummy = string(parse.nonOption(0));
      } else {
	// ASDMDIR must have been set
	dummy = string(options[ASDMDIR].last()->arg);
      }
      dsName = lrtrim(dummy) ;
      if (dsName.back()=='/') dsName.erase(dsName.size()-1);
    } else {
      errstream.str("");
      option::printUsage(errstream, usage);
      error(errstream.str());
    }
    
    if (parse.nonOptionsCount() > 1 || options[MSDIRPREFIX]) {
      string dummyMSName;
      if (parse.nonOptionsCount() > 1) {
	dummyMSName = string(parse.nonOption(1));
      } else {
	// MSDIRPREFIX must be set
	dummyMSName = string(options[MSDIRPREFIX].last()->arg);
      }
      dummyMSName = lrtrim(dummyMSName);
      if (dummyMSName.back()=='/') dummyMSName.erase(dummyMSName.size()-1);
      Path msPath(dummyMSName);
      string msDirectory = msPath.dirName();
      if (msDirectory.size() == 0) msDirectory = ".";
      // extract the prefix and extension. Prefix is everything before any final "."
      // extension is everything after any final ".", including that final dot.
      string msPathBasename = msPath.baseName();
      size_t rdot = msPathBasename.find_last_of('.');
      if (rdot != std::string::npos) {
	msNameExtension = msPathBasename.substr(rdot,msPathBasename.size()-rdot);
	msNamePrefix = msPathBasename.substr(0,rdot);
      } else {
	msNameExtension = "";
	msNamePrefix = msPathBasename;
      } 
      msNamePrefix = msDirectory + "/" + msNamePrefix;
    } else {
      msNamePrefix = dsName;
      msNameExtension = ".ms";
    }
    
    // Does the user want compressed columns in the resulting MS ?
    if ((withCompression = (options[COMPRESSION] != NULL))) {
      infostream.str("");
      infostream << "Compressed columns in the resulting MS(s) : Yes" ;
      info(infostream.str());
    } else {
      infostream.str("");
      infostream << "Compressed columns in the resulting MS(s) : No" ;
      info(infostream.str());
    }

    // WVR uncorrected and|or corrected data required ?
    // always available because of defaults
    string wvr_corrected_data = string(options[WVRCORRDATA].last()->arg);
    if (wvr_corrected_data.compare("no") == 0)
      es_query_apc.fromString("AP_UNCORRECTED");
    else if (wvr_corrected_data.compare("yes") == 0)
      es_query_apc.fromString("AP_CORRECTED");
    else if (wvr_corrected_data.compare("both") == 0)
      es_query_apc.fromString("AP_CORRECTED AP_UNCORRECTED");
    else {
      errstream.str("");
      errstream << "Token '" << wvr_corrected_data << "' invalid for wvr-corrected-data." << endl;
      option::printUsage(errstream, usage);
      error(errstream.str());
    }
    
    // Do we want an MS Main table to be filled or not ?
    mute = (options[DRYRUN] != NULL);
    if (mute) {
      infostream.str("");
      infostream << "option dry-run is used, the MS Main table will not be filled" << endl;
      info(infostream.str());
    }

    // What is the amount of memory allocated to the BDF slices.
    infostream.str("");
    infostream << "the BDF slice size is set to " << bdfSliceSizeInMb << " megabytes." << endl;
    info(infostream.str());

    lazy = options[LAZY] != NULL;

    if (lazy && !lazyModeOK) {
      lazy = false;
      infostream.str("");
      infostream << "The lazy filler can not be used with any input data selection (icm, isrt, its)." << endl;
      infostream << "The non-lazy version of the filler will be used." << endl;
      // an option is being ignored, warn the user even if not verbose
      warning(infostream.str());
    }

    if (debug) {
      cout << "checkdupints : " << checkdupints << endl;
    }

    // Do we consider another order than ac_xc_per_timestamp ?
    // always available because of defaults
    string acXcOpt = string(options[ACXCPERTIME].last()->arg);
    ac_xc_per_timestamp = str_tolower(acXcOpt) == "yes";

    // Do we want to tabulate polynomial present in the ephemeris table ?
    // This is inferred by seeing if the user specifically set this.
    // But it's set at least once because of defaults, so, look for there being
    // more than one of these in options.
    tabulate_ephemeris_polynomials = (options[POLYEPHTSTEP].count() > 1);
    
    if (tabulate_ephemeris_polynomials) {
      // If we tabluate then ignore all the other options about ephemeris    
      //  infostream.str();
      //  infostream << "The MS Ephemeris table(s) will be produced by tabulating the polynomials found in the columns 'dir', 'distance' and optionally 'radVel' with a timestep of '"
      //		 << polyephem_tabtimestep << "' day, i.e. '"<< ((uint64_t) (polyephem_tabtimestep * 86400 * 1.e09)) <<"' nanoseconds."; 
      // info(infostream.str());
    } else {
      // Do we want interpolate the values found in the ASDM Ephemeris table or not ?
      // always available because of defaults
      string intEphOpt = string(options[INTEPHEM].last()->arg);
      interpolate_ephemeris = str_tolower(intEphOpt) == "yes";
      infostream.str("");
      if (interpolate_ephemeris) {
	infostream << "the MS Ephemeris table(s) will be produced by interpolation of the values present in the ASDM Ephemeris table on a resampled time grid."; 
      } else {
	infostream << "the MS Ephemeris tables(s) will be produced by simple copies of the values found in the ASDM Ephemeris table with just units conversion.";
      }
      info(infostream.str());
    }

    checkRowUniqueness = options[CHECKROWUNIQ] != NULL;
    if (options[SCANS]) {
      scansOptionInfo = string(options[SCANS].last()->arg);
    }
    if (options[ASIS]) {
      asisOption = string(options[ASIS].last()->arg);
    }
    ignoreTime = options[IGNORETIME] != NULL;
    processSysPower = options[NOSYSPOWER] == NULL;
    processCalDevice = options[NOCALDEV] == NULL;
    processPointing = options[NOPOINTING] == NULL;
    withPointingCorrection = options[WITHPCORR] != NULL;
    processEphemeris = options[NOEPHEMERIS] == NULL;
  } catch (std::exception& e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }

  //
  // Try to open an ASDM dataset whose name has been passed as a parameter on the command line
  //
  if ( (dsName.size() > 0) && dsName.at(dsName.size()-1) == '/' ) dsName.erase(dsName.size()-1);

  double cpu_time_parse_xml  = 0.0;
  double real_time_parse_xml = 0.0;
  int mode;
  mode = 0; myTimer(&cpu_time_parse_xml, &real_time_parse_xml, &mode);

  ASDM* ds = new ASDM();

  infostream.str("");
  if (checkRowUniqueness) 
    infostream << "Row uniqueness constraint will be applied." << endl;
  else
    infostream << "Row uniqueness constraint will be ignored." << endl;

  info(infostream.str());

  try {
    infostream.str("");
    infostream << "Input ASDM dataset : " << dsName << endl;
    info(infostream.str());
    
    ASDMParseOptions parse = ASDMParseOptions().loadTablesOnDemand(true).checkRowUniqueness(checkRowUniqueness);
    ds->setFromFile(dsName, parse);
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (std::exception e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }
  catch (...) {
    errstream.str("");
    errstream << "Uncaught exception !" << endl;
    error(errstream.str());
  }
  
  mode = 1; myTimer(&cpu_time_parse_xml, &real_time_parse_xml, &mode);
  infostream.str("");
  infostream << "Time spent parsing the ASDM medata : " << cpu_time_parse_xml << " s.";
  info(infostream.str());
  
  //
  // What are the apc literals present in the binary data.
  //
  try {
    es_apc = apcLiterals(*ds);
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  
  //
  // Determine what kind of data complex (DATA column) or float (FLOAT_DATA) will be
  // stored in the measurement set by using the method isDataComplex on the first row of 
  // the ASDM Dataset. 
  // This method called on all remaining rows of the ASDM should return the same result
  // otherwise the filling process will stop.
  //
  SDMBinData sdmBinData(ds, dsName);

  // Define the SDM Main table subset of rows to be accepted
  sdmBinData.select( es_cm, es_srt, es_ts);   
  
  // From now we decide to extract the data with all atmospheric phase corrections.
  // The selection will be done on output. Michel Caillat Thur 18 Sept 2014 - CAS-6935
  EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection> allAPCs;
  allAPCs.fromString("AP_CORRECTED AP_UNCORRECTED");
  sdmBinData.selectDataSubset(e_query_cm, allAPCs);
  
  //
  // Selection of the scans to consider.
  //
  vector<ScanRow *>	scanRow_v;
  try {
    scanRow_v   = ds->getScan().get();
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }

  map<int, set<int> > all_eb_scan_m;
  for (vector<ScanRow *>::size_type i = 0; i < scanRow_v.size(); i++)
    all_eb_scan_m[scanRow_v[i]->getExecBlockId().getTagValue()].insert(scanRow_v[i]->getScanNumber());
  
  vector<ScanRow *>	selectedScanRow_v;
  map<int, set<int> >   selected_eb_scan_m;
  
  if (scansOptionInfo.size()>0) {
    map<int, set<int> > eb_scan_m;
    int status = scansParser(scansOptionInfo, eb_scan_m);
        
    if (status == 0) {
      errstream.str("");
      errstream << "'" << scansOptionInfo << "' is an invalid scans selection." << endl;
      error(errstream.str());
    }

    vector<ScanRow *> scanRow_v = ds->getScan().get();
    map<int, set<int> >::iterator iter_m = eb_scan_m.find(-1);

    if (iter_m != eb_scan_m.end())
      for ( map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin();
            iterr_m != all_eb_scan_m.end(); iterr_m++ ) {
          if ((iter_m->second).empty())
              selected_eb_scan_m[iterr_m->first] = iterr_m->second;
          else
              selected_eb_scan_m[iterr_m->first] = SetAndSet<int>(iter_m->second, iterr_m->second);
      }

    for ( map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin();
          iterr_m != all_eb_scan_m.end(); iterr_m++)
        if ((iter_m=eb_scan_m.find(iterr_m->first)) != eb_scan_m.end()) {
            if ((iter_m->second).empty())
                selected_eb_scan_m[iterr_m->first].insert((iterr_m->second).begin(), (iterr_m->second).end());
            else {
                set<int> s = SetAndSet<int>(iter_m->second, iterr_m->second);
                selected_eb_scan_m[iterr_m->first].insert(s.begin(), s.end());
            }
        }
    
    ostringstream	oss;
    oss << "The following scans will be processed : " << endl;
    for ( map<int, set<int> >::const_iterator iter_m = selected_eb_scan_m.begin();
          iter_m != selected_eb_scan_m.end(); iter_m++ ) {
      oss << "eb#" << iter_m->first << " -> " << displaySet<int>(iter_m->second) << endl;
      Tag execBlockTag  = Tag(iter_m->first, TagType::ExecBlock);
      for ( set<int>::const_iterator iter_s = iter_m->second.begin();
	   iter_s != iter_m->second.end();
	   iter_s++ )
          selectedScanRow_v.push_back(ds->getScan().getRowByKey(execBlockTag, *iter_s));

    }

    scansOptionInfo = oss.str();
  } else {
    selectedScanRow_v = ds->getScan().get();
    selected_eb_scan_m = all_eb_scan_m;
    scansOptionInfo = "All scans of all exec blocks will be processed \n";
  }

  //
  // Report the selection's parameters.
  //
  infostream.str("");
  infostream << "Correlation modes requested : " << e_query_cm.str() << endl;
  infostream << "Spectral resolution types requested : " << es_srt.str() << endl;
  infostream << "Time sampling requested : " << es_ts.str() << endl;
  infostream << "WVR uncorrected and|or corrected data requested : " << es_query_apc.str() << endl;
  if (selectedScanRow_v.size() == 0) { 
    errstream.str("");
    errstream << "No scan number corresponding to your request. Can't go further.";
    error(errstream.str());
  }

  infostream << scansOptionInfo;

  if (ignoreTime)
    infostream << "All rows of the tables depending on time intervals will be processed independently of the selected exec block / scan.";
  info(infostream.str());

  infostream.str("");
  if (!processSysPower)   infostream << "The SysPower table will not be processed." << endl;
  if (!processCalDevice)  infostream << "The CalDevice table will not be processed." << endl;
  if (!processPointing)   infostream << "The Pointing table will not be processed." << endl;
  if (processPointing && withPointingCorrection ) infostream << "The correction (encoder - pointingDirection) will be applied" << endl;
  if (!processEphemeris)  infostream << "The Ephemeris table will not be processed." << endl;
  if (ac_xc_per_timestamp)
      infostream << "For each data description for each timestamp auto correlations followed by cross correlations will be written in the Main table" << endl;
  else 
    infostream << "For each data description auto correlations for a sequence of timestamps followed by cross correlations for the same sequence will be written in the Main table" << endl;

  info(infostream.str());
  //
  // Shall we have Complex or Float data ?
  //
  bool complexData = true;
  try {
    complexData =  sdmBinData.isComplexData();
  } catch (Error & e) {
    errstream.str("");
    errstream << e.getErrorMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }


  //
  // Prepare a map AtmPhaseCorrection -> name of measurement set.
  // Three cases are possible :
  // only AP_CORRECTED -> MS name is suffixed with "-wvr-corrected",
  // only AP_UNCORRECTED -> MS name has no particular suffix,
  // AP_CORRECTED and AP_UNCORRECTED -> 2 MSs whith names defined by the two above rules.
  //
  map<AtmPhaseCorrectionMod::AtmPhaseCorrection, string> msNames;
  if (hasCorrectedData(es_apc) && es_query_apc[AtmPhaseCorrectionMod::AP_CORRECTED]) {
    msNames[AtmPhaseCorrectionMod::AP_CORRECTED] = msNamePrefix + "-wvr-corrected";
  }
  if (hasUncorrectedData(es_apc) && es_query_apc[AtmPhaseCorrectionMod::AP_UNCORRECTED]) {
    msNames[AtmPhaseCorrectionMod::AP_UNCORRECTED] = msNamePrefix;
  }
      
  if (msNames.size() == 0) {
    //
    // no MS can be produced due to the selection parameters values.
    // 
    infostream.str("");
    infostream << "No measurement set can be produced with  your selection criteria on '" << dsName << "'" << endl;
    info(infostream.str());
    delete ds;
    exit(1);
  } else {
    //
    // OK, we are going to produce at least one MS.
    // If the '--compression' option has been used then append a suffix ".compressed" to
    // the MS name.
    // And eventually always suffix with ".ms".
    //
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, string>::iterator iter=msNames.begin(); iter != msNames.end(); ++iter) {
      if (withCompression)
	iter->second = iter->second + ".compressed";
      iter->second +=  msNameExtension;
    }
  }

  infostream.str("");
  infostream << "The resulting measurement set will contain a '" << ((complexData) ? "DATA" : "FLOAT_DATA") << "' column" << endl; 

#if DDPRIORITY
  // The data returned by getDataCols will be ordered in priority by Data Descriptions.
  sdmBinData.setPriorityDataDescription();
#endif

  //get numCorr, numChan, telescope name for setupMS
  
  ExecBlockTable& temp_execBlockT = ds->getExecBlock();
  //take first row of the table (assuming telescope name is all the same)
  ExecBlockRow* temp_ebtrow = temp_execBlockT.get()[0];
  string telName  = temp_ebtrow->getTelescopeName();
  //cout<<"telName="<<telName<<endl;

  int maxNumCorr =1;
  PolarizationTable& temp_polT = ds->getPolarization();
  PolarizationRow* temp_poltrow;
  for (unsigned int i=0; i<temp_polT.size(); i++) {
    temp_poltrow = temp_polT.get()[i];
    maxNumCorr=max(maxNumCorr, temp_poltrow->getNumCorr());
  }

  //need to add analysis of max NumChan
  int maxNumChan=1;
  SpectralWindowTable& temp_spwT = ds->getSpectralWindow();
  SpectralWindowRow* temp_spwtrow;

  vector<int> SwIds;
  try {
    for (unsigned int i=0; i<temp_spwT.size(); i++) {
      temp_spwtrow = temp_spwT.get()[i];
      maxNumChan=max(maxNumChan, temp_spwtrow->getNumChan());
      SwIds.push_back(temp_spwtrow->getSpectralWindowId().getTagValue());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }

  // Create the measurement set(s). 
  if (!false) {
    try {
      if (lazy)  casa::AsdmStMan::registerClass();
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, string>::iterator iter = msNames.begin(); iter != msNames.end(); ++iter) {
	info("About to create a filler for the measurement set '" + msNames[iter->first] + "'");
	msFillers[iter->first] = new ASDM2MSFiller(msNames[iter->first], 0.0, complexData,
						   withCompression, telName, maxNumCorr, maxNumChan,
						   false, lazy);
      }
    } catch(AipsError & e) {
      errstream.str("");
      errstream << e.getMesg();
      error(errstream.str());
    } catch (std::exception & e) {
      errstream.str("");
      errstream << e.what();
      error(errstream.str());
    }

    msFiller = msFillers.begin()->second;    
  }

  //
  // Firstly convert the basic tables.
  //
  // For purpose of simplicity we assume that in all ASDM basic tables having a Tag identifier
  // these Tag values are forming a sequence of integer 0 -> table.size() - 1. If that's not the case, the
  // program aborts.
  //

  //
  // Process the Antenna table.
  // 
  // (This part needs the Station table)
  //
  // At the same time, we populate a map ASDM Station Tag -> MS ANTENNA ID
  // which will be useful when the Weather table will be converted.
  //
  //unsigned int numTrueAntenna
  
  try { 
    AntennaTable& antennaT = ds->getAntenna();
    AntennaRow*   r   = 0;

    int nAntenna = antennaT.size();
    infostream.str("");
    infostream << "The dataset has " << nAntenna << " antenna(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nAntenna; i++) {
      if ((r = antennaT.getRowByKey(Tag(i, TagType::Antenna))) == 0){
	errstream.str("");
	errstream << "Problem while reading the Antenna table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl;
	error(errstream.str());
      }

      // The MS Antenna position is defined as the sum of the ASDM station position and
      // of the ASDM Antenna position after applying to it a coordinate system transformation.
      // Since the ASDM Antenna position is 0,0,0 for now, we only use the ASDM station position.
      // Update - 2012-03-22
      // Now the ASDM Antenna position contains non-zeros  so need to take account
      // for this now as shown below. For EVLA, this is still 0,0,0.       
      vector<Length> position = r->getStationUsingStationId()->getPosition();
      double xStation = position.at(0).get();
      double yStation = position.at(1).get();
      double zStation = position.at(2).get();
      
      // ---- transform antenna position to geocentric coordinates ----
      // Method 1 - assume z axis of antenna position lines up with station vector
      //            and transformation is done by a rotation matrix based on 
      //            geocntric longitude and latitude. Good enough for current 
      //            antenna position measurement accuracy.
      //
      // geocentric longitude  and latitude
      double glat = atan2(zStation,sqrt(xStation*xStation + yStation*yStation)); 
      double glon = atan2(yStation,xStation);

      // get ASDM Antenna position vector
      vector<Length> antPosition = r->getPosition();
      
      vector<double> cartesianAnt1(3, 0.0);
      vector<double> cartesianAnt2(3, 0.0);
      vector<vector<double> > matrixAnt3x3;
      for (unsigned int ii = 0; ii < 3; ii++) {
        matrixAnt3x3.push_back(cartesianAnt1);
      }
      cartesianAnt1[0] = antPosition.at(0).get();
      cartesianAnt1[1] = antPosition.at(1).get();
      cartesianAnt1[2] = antPosition.at(2).get();

      if (cartesianAnt1[0]!=0.0 || cartesianAnt1[1]!=0 || cartesianAnt1[2]!=0.0) {
        topo2geomat(glon,glat,matrixAnt3x3);
        matvec(matrixAnt3x3,cartesianAnt1,cartesianAnt2); 
      }
      
      /*** 
       // Method 2 - use Measures and let Measure figure out
       //            transoformation for local geodetic coordinates (with
       //            earth's oblateness taken account) to geocentric
       //            coordinates.
       //            Use AZELGEO as a coordinate ref to be more precise
       casacore::Vector<casacore::Quantity> vq; vq.resize(3);
       vq[0] = casacore::Quantity(antPosition.at(0).get(),"m");
       vq[1] = casacore::Quantity(antPosition.at(1).get(),"m");
       vq[2] = casacore::Quantity(antPosition.at(2).get(),"m");
       casacore::MVPosition mvp(vq);
       casacore::MVBaseline mvb(mvp);

       // setup conversion template
       double anttime =  ((double) r->getTime().get()) / ArrayTime::unitsInASecond ;
       casacore::MEpoch ep(casacore::Quantity(anttime,"s"), casacore::MEpoch::UTC);
       casacore::Vector<casacore::Quantity> rvq; rvq.resize(3);
       //station vector in ITRF
       rvq[0] = casacore::Quantity(xStation,"m");
       rvq[1] = casacore::Quantity(yStation,"m");
       rvq[2] = casacore::Quantity(zStation,"m");
       casacore::MVPosition rmvp(rvq);
       casacore::MPosition rmp(rmvp,casacore::MPosition::ITRF);
       // set the direction to the pole
       casacore::MVDirection mvd = casacore::MVDirection();
       // this approximate the antenna position to be in topocentric, z is parallel
       // to the station vector. 
       casacore::MDirection mdir(mvd,casacore::MDirection::AZEL);
       // to be precise set the ref to geodetic local coordinates  
       //casacore::MDirection mdir(mvd,casacore::MDirection::AZELGEO);
       casacore::MeasFrame mFrame(rmp,ep,mdir);

       casacore::MBaseline baseMeas;
       casacore::MVBaseline mantv;
       casacore::MBaseline::Ref baseref(MBaseline::AZEL, mFrame);
       // geodetic local coordinates case
       //casacore::MBaseline::Ref baseref(MBaseline::AZELGEO, mFrame);
       baseMeas.set(mantv, baseref);
       baseMeas.getRefPtr()->set(mFrame);

       casacore::MBaseline mb(mvb,baseref);
       casacore::MBaseline::Convert antvconv(baseMeas, MBaseline::Ref(MBaseline::ITRF));
       casacore::MBaseline mbantp = antvconv(mb);
      
       //compare transformed antenna positions in the two methods
       //for the measure frame with AZEL those two values should be identical
       cerr<<"rotated measAnt(0)="<<mbantp.getValue()(0)<<" measAnt(1)="<<mbantp.getValue()(1)<<" measAnt(2)="<<mbantp.getValue()(2)<<endl;
       // end of Method2
       ***/

      //cerr<<"rotated cartAnt(0)="<<cartesianAnt2.at(0)<<" cartAnt(1)="<<cartesianAnt2.at(1)<<" cartAnt(2)="<<cartesianAnt2.at(2)<<endl;

      //add antenna position: for now use ones obtained by the topo2geomat() 
      double xPosition = position.at(0).get() + cartesianAnt2.at(0);
      double yPosition = position.at(1).get() + cartesianAnt2.at(1);
      double zPosition = position.at(2).get() + cartesianAnt2.at(2);

      vector<Length> offset = r->getOffset();
      double xOffset = offset.at(0).get();
      double yOffset = offset.at(1).get();
      double zOffset = offset.at(2).get();

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	/*
	** "&" are not recommanded in antenna names in order to avoid bumps with MS Selection syntax.
	*/ 
	string aName = r->getName();
	if (aName.find_first_of('&')) replace(aName.begin(), aName.end(), '&', '#');
	
	static_cast<void>(iter->second->addAntenna(aName, r->getStationUsingStationId()->getName(),
						   xPosition, yPosition, zPosition, xOffset, yOffset,
						   zOffset, (float)r->getDishDiameter().get()));
      }
    }

    int numTrueAntenna = msFillers.begin()->second->ms()->antenna().nrow();
    if (numTrueAntenna) {
      infostream.str("");
      infostream << "converted in " << numTrueAntenna << " antenna(s)  in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  //
  // Process the SpectralWindow table.
  //
  map<unsigned int, double> effectiveBwPerSpwId_m;
  fillSpectralWindow(ds, effectiveBwPerSpwId_m, telName);

  //
  // Process the Polarization table
  //
  Stokes::StokesTypes linearCorr[] = { Stokes::XX, Stokes::XY, Stokes::YX, Stokes::YY };
  Stokes::StokesTypes circularCorr[] = { Stokes::RR, Stokes::RL, Stokes::LR, Stokes::LL };
  int corrProduct1[] = { 0, 0 };
  int corrProduct2[] = { 0, 0, 1, 1};
  int corrProduct4[] = { 0, 0, 0, 1, 1, 0, 1, 1 };
			 
  vector<int> polarizationIdx2Idx;
  int pIdx;

  try {
    PolarizationTable& polT = ds->getPolarization();  
    PolarizationRow* r = 0;
    int nPolarization = polT.size();
    infostream.str("");
    infostream << "The dataset has " << nPolarization << " polarization(s)..."; 
    info(infostream.str());

    for (int i = 0; i < nPolarization; i++) {
      if ((r=polT.getRowByKey(Tag(i, TagType::Polarization))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Polarization table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      int numCorr = r->getNumCorr();
      if (numCorr < 1 || numCorr > 4) {
	ostringstream oss ;
	oss << "a polarization row cannot be processed due to  'numCorr = " << numCorr << "'.";
	throw ASDM2MSException(oss.str());
      }

      //Stokes::StokesTypes * corrType;
      //vector<Stokes::StokesTypes> corrType;
      vector<int> corrType;
      StokesMapper stokesMapper;
      if (numCorr != 3) {
	corrType = StokesMapper::toVectorI(r->getCorrType());
      } else {
	numCorr  = 4;
	StokesParameterMod::StokesParameter sp = r->getCorrType()[0];
	if ((sp == StokesParameterMod::RR) ||
	    (sp == StokesParameterMod::LL) ||
	    (sp == StokesParameterMod::RL) ||
	    (sp == StokesParameterMod::LR)) {
	  corrType.resize(4);
	  copy (circularCorr, circularCorr+4, corrType.begin());
	} else if ((sp == StokesParameterMod::XX) ||
		 (sp == StokesParameterMod::XY) ||
		 (sp == StokesParameterMod::YX) ||
		 (sp == StokesParameterMod::YY)) {
	  corrType.resize(4);
	  copy (linearCorr, linearCorr+4, corrType.begin());
	} else {
	  errstream.str("");
	  errstream << " I don't know what to do with the given Stokes parameters for autocorrelation data" << endl;
	  error(errstream.str());
	}
	  
      }
      
      
      /*int* corrProduct = 0;*/
      vector<int> corrProduct; 
      switch (numCorr) {
      case 1: corrProduct.resize(2); copy(corrProduct1, corrProduct1+2, corrProduct.begin()); break;
      case 2: corrProduct.resize(4); copy(corrProduct2, corrProduct2+4, corrProduct.begin()); break;
      case 4: corrProduct.resize(8); copy(corrProduct4, corrProduct4+8, corrProduct.begin()); break;
      }


      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	pIdx = iter->second->addUniquePolarization( numCorr, corrType, corrProduct );
      }
      polarizationIdx2Idx.push_back(pIdx);
    }
    if (nPolarization) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->polarization().nrow() << " polarization(s)." ;
      info(infostream.str());
    }
  } catch (ASDM2MSException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
   
  //
  // Process the DataDescription table.
  //

  std::map<unsigned int, double> effectiveBwPerDD_m;
  try {
    DataDescriptionTable& ddT = ds->getDataDescription();
    DataDescriptionRow* r = 0;
    int nDataDescription = ddT.size();
    infostream.str("");
    infostream << "The dataset has " << nDataDescription << " data description(s)...";
    info(infostream.str());

    for (int i = 0; i < nDataDescription; i++) {
      if ((r=ddT.getRowByKey(Tag(i, TagType::DataDescription))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the DataDescription table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	ddIdx = iter->second->addUniqueDataDescription(swIdx2Idx[r->getSpectralWindowId().getTagValue()],
						       polarizationIdx2Idx.at(r->getPolOrHoloId().getTagValue()));
      }
      dataDescriptionIdx2Idx.push_back(ddIdx); 

      effectiveBwPerDD_m[ddIdx] = effectiveBwPerSpwId_m[r->getSpectralWindowId().getTagValue()];
    }

    if (nDataDescription) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->dataDescription().nrow() << " data description(s)  in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  // 
  // Process the ExecBlock table,
  // in order to build the MS Observation table.
  // 
  string telescopeName; 
  try {
    const ExecBlockTable& execBlockT = ds->getExecBlock(); 
    ExecBlockRow* r = 0;
    int nExecBlock = execBlockT.size();
    infostream.str("");
    infostream << "The dataset has " << nExecBlock << " execBlock(s) ...";

    vector<ExecBlockRow *> temp_v = execBlockT.get();
    vector<ExecBlockRow *> v;
    for (vector<ExecBlockRow *>::iterator iter_v = temp_v.begin(); iter_v != temp_v.end(); iter_v++)
      if ( selected_eb_scan_m.find((*iter_v)->getExecBlockId().getTagValue()) != selected_eb_scan_m.end() )
	v.push_back(*iter_v);
    
    vector<string> schedule; schedule.resize(2);

    infostream << v.size() << " of them in the selected exec blocks / scans ... ";
    info(infostream.str());

    for (unsigned int i = 0; i < v.size(); i++) {
      r = v.at(i);
      
      telescopeName	     = r->getTelescopeName();
      double	startTime    = r->getStartTime().getMJD()*86400;
      double	endTime      = r->getEndTime().getMJD()*86400;
      string	observerName = r->getObserverName();

      vector<string> observingLog = r->getObservingLog();

      string scheduleType(r->getTelescopeName());
      schedule[0] = "SchedulingBlock " + ds->getSBSummary().getRowByKey(r->getSBSummaryId())->getSbSummaryUID().getEntityId().toString();
      schedule[1] = "ExecBlock " + r->getExecBlockUID().getEntityId().toString();
      string project(r->getProjectUID().getEntityId().toString());
      double releaseDate = r->isReleaseDateExists() ? r->getReleaseDate().getMJD():0.0;

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addObservation( telescopeName, startTime, endTime, observerName, observingLog,
                                      scheduleType, schedule, project, releaseDate );
      }
      if (i==0) { // assume same telescope for all execBlocks 
        if (telescopeName.find("EVLA")!=string::npos) {
          isEVLA=true;
        } else if (telescopeName.find("OSF")!=string::npos || 
                 telescopeName.find("AOS")!=string::npos || 
                 telescopeName.find("ATF")!=string::npos) {
          isEVLA=false;
        }     
        string telname = (isEVLA ? "EVLA" : "ALMA");
        infostream.str("");
        infostream << "Telescope Name:" <<telescopeName << ", process as "<<telname<<" data." ; 
        info(infostream.str());
      }
    }
    if (nExecBlock) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->observation().nrow() << " observation(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  //
  // Process the Feed table
  // Issues :
  //    - time (epoch) : at the moment it takes directly the time as it is stored in the ASDM.
  //    - focusLength (in AIPS++) is no defined.
  try {
    const FeedTable& feedT = ds->getFeed();
    FeedRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << feedT.size() << " feed(s)...";
    rowsInAScanbyTimeIntervalFunctor<FeedRow> selector(selectedScanRow_v);
    
    const vector<FeedRow *>& v = selector(feedT.get(), ignoreTime);
    if (!ignoreTime)
      infostream << v.size() << " of them in the exec blocks / selected scans ... ";
    
    info(infostream.str());
    int nFeed = v.size();
    for (int i = 0; i < nFeed; i++) {
      r = v.at(i);
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval/2.0;
      //}

      vector<double> beam_offset_ =  DConverter::toVectorD(r->getBeamOffset());
      vector<std::string> polarization_type_ = PolTypeMapper().toStringVector(r->getPolarizationTypes());
      vector<complex<float> > polarization_response_ = CConverter::toVectorCF(r->getPolResponse());
      vector<double> xyzPosition (3, 0.0);
      if (r->isPositionExists()) {
	vector<Length> position = r->getPosition();
	if (position.size() != 3) {
	  errstream.str("");
	  errstream << "The size of attribute position ('" 
		    << position.size()
		    << "') is not equal to 3. Can't go further."
		    << endl;
	  error(errstream.str());
	}
	
	xyzPosition = DConverter::toVectorD<Length>(position);
      }
      vector<double> receptor_angle_ = DConverter::toVectorD<Angle>(r->getReceptorAngle());
      
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addFeed( (int) r->getAntennaId().getTagValue(),r->getFeedId(),
                               swIdx2Idx[r->getSpectralWindowId().getTagValue()],
                               time, interval, r->getNumReceptor(), -1,  // We ignore the beamId array
                               beam_offset_, polarization_type_,
                               polarization_response_, xyzPosition, receptor_angle_);
      } 
    }
    if (nFeed) {
      infostream.str("");
      infostream <<  "converted in " << msFillers.begin()->second->ms()->feed().nrow() << " feed(s) in the measurement set." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  // Process the Ephemeris table.
  //
  // Create and fill the MS ephemeris table(s) with a time interpolation time step set to 86400000000 nanoseconds ( 1/1000 day).
  if (processEphemeris) {
    uint64_t timeStepInNanoSeconds = polyephem_tabtimestep * 86400 * 1.e09;
    fillEphemeris(ds, timeStepInNanoSeconds, interpolate_ephemeris, telescopeName);
  }

  // Process the Field table.
  // Now it respects the degree of the polynomials but it ignores the ephemerisId.
  // The ephemerisId will be processed during the call to fillEphemeris.
  //
  fillField(ds, processEphemeris);
   
  // Process the FlagCmd table.
  //
  try {
    const FlagCmdTable& flagCmdT  = ds->getFlagCmd();
    FlagCmdRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << flagCmdT.size() << " FlagCmd(s)...";
    rowsInAScanbyTimeIntervalFunctor<FlagCmdRow> selector(selectedScanRow_v);

    const vector<FlagCmdRow *>& v = selector(flagCmdT.get(), ignoreTime);
    if (!ignoreTime)
      infostream << v.size() << " of them in the exec blocks / selected scans ... ";

    info(infostream.str());
    int nFlagCmd = v.size();
    for (int i = 0; i < nFlagCmd; i++) {
      r = v.at(i);
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond ;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval/2.0;
      //}
      string type = r->getType();
      string reason = r->getReason();
      string command = r->getCommand();
   
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addFlagCmd( time, interval, type, reason, r->getLevel(), r->getSeverity(),
                                  r->getApplied() ? 1 : 0, command);
      }
    }
    if (nFlagCmd) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->flagCmd().nrow() << " in the measurement set." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  // Process the History table.
  // Issues :
  // - use executeBlockId for observationId ...to be discussed with Francois.
  // - objectId : not taken into account (it's a string while the MS expects an int).
  try {
    const HistoryTable& historyT = ds->getHistory();
    HistoryRow* r = 0;
    int nHistory = historyT.size();
    infostream.str("");
    infostream << "The dataset has " << nHistory << " history(s)...";
    rowsInAScanbyTimeFunctor<HistoryRow> selector(selectedScanRow_v);

    const vector<HistoryRow *>& v = selector(historyT.get(), ignoreTime);;
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected exec blocks / scans ... ";

    info(infostream.str()); 

    for (int i = 0; i < nHistory; i++) {
      r = v.at(i);
      double time =  ((double) r->getTime().get()) / ArrayTime::unitsInASecond ;
      string message     = r->getMessage();
      string priority    = r->getPriority();
      string origin      = r->getOrigin();
      string application = r->getApplication();
      string cliCommand  = r->getCliCommand();
      string appParams   = r->getAppParms();
 
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addHistory( time, r->getExecBlockId().getTagValue(), message,
                                  priority, origin, -1, application, cliCommand,
                                  appParams);
      }
    }
    if (nHistory) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->history().nrow() << " history(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  //
  // Process the Pointing table.
  // Issues :
  // - pointingModelId , phaseTracking, sourceOffset and overTheTop not taken into account.

  if (processPointing) 
    try {
        fillPointing(dsName, ds, ignoreTime, withPointingCorrection, selectedScanRow_v, msFillers);
    }
    catch (ConversionException e) {
      errstream.str("");
      errstream << e.getMessage();
      error(errstream.str());
    } catch ( std::exception & e) {
      errstream.str("");
      errstream << e.what();
      error(errstream.str());      
    }

    
  // Process the processor table.
  //

  try {
    ProcessorTable& processorT = ds->getProcessor();
    ProcessorRow* r = 0;
    int nProcessor = processorT.size();

    infostream.str("");
    infostream << "The dataset has " << nProcessor << " processor(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nProcessor; i++) {
      if ((r=processorT.getRowByKey(Tag(i, TagType::Processor))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Processor table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      string processorType    = CProcessorType::name(r->getProcessorType());
      string processorSubType = CProcessorSubType::name(r->getProcessorSubType());
      
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addProcessor( processorType, processorSubType,
                                    -1,    // Since there is no typeId in the ASDM.
                                    r->getModeId().getTagValue());
      }  
    }
    if (nProcessor) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->processor().nrow() << " processor(s) in the measurement set." ;
      info(infostream.str());
    } 
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  // Process the Source table.
  //
  const SourceTable& sourceT = ds->getSource();
  try {
    SourceRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << sourceT.size() << " sources(s)...";
    rowsInAScanbyTimeIntervalFunctor<SourceRow> selector(selectedScanRow_v);
    
    const vector<SourceRow *>& v = selector(sourceT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";

    info(infostream.str());
    int nSource = v.size();

    for (int i = 0; i < nSource; i++) {
      r = v.at(i);
      //
      // Check some assertions. 
      // For each row of the Source table, if any of the optional attributes which is an array
      // and depend on numLines for the size of one of its dimensions then the (optional) \
      // attribute numLines must be present and the dimensions depending on on numLines must be 
      // consistent with the value of numLines.
      //
      int numLines = r->isNumLinesExists() ? r->getNumLines() : 0;
      
      if (r->isTransitionExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'transition' exists but the attribute 'numLines' which serves to define its shape is missing. Can't go further.";
	  error(errstream.str());
	}

	int transitionSize = r->getTransition().size();
	if (numLines != transitionSize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'transition' (" << transitionSize << "). Can't go further.";
	  error(errstream.str());
	}
      }

      if (r->isRestFrequencyExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'restFrequency' exists but the attribute 'numLines' which serves to define its shape is missing. Cant' go further.";
	  error(errstream.str());
	}
	
	int restFrequencySize = r->getRestFrequency().size();
	if (numLines != restFrequencySize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'restFrequency' (" << restFrequencySize << "). Can't go further.";
	  error(errstream.str());
	}
      }

      if (r->isSysVelExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'sysVel' exists but the attribute 'numLines' which serves to define its shape is missing. Cant' go further.";
	  error(errstream.str());
	}
	
	int sysVelSize = r->getSysVel().size();
	if (numLines != sysVelSize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'sysVel' (" << sysVelSize << "). Can't go further.";
	  error(errstream.str());
	}
      }          

      int sourceId = r->getSourceId();
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  r->getTimeInterval().getStartInMJD()*86400 ;
      // }
      // else {
      time =  r->getTimeInterval().getStartInMJD()*86400 + interval / 2.0 ;
      //}

      int spectralWindowId = swIdx2Idx[r->getSpectralWindowId().getTagValue()];

      string sourceName = r->getSourceName();

      int calibrationGroup = r->isCalibrationGroupExists() ? r->getCalibrationGroup() : 0;
      
      string code = r->getCode();

      vector<double> direction = DConverter::toVectorD(r->getDirection());
 
      DirectionReferenceCodeMod::DirectionReferenceCode dirRefCode = DirectionReferenceCodeMod::J2000;
      if(r->isDirectionCodeExists()){
	dirRefCode = r->getDirectionCode();
	//cout << "found directionCode for source " << sourceName << ": ";
      }
      //else{
      //  cout << "No directionCode in input table. Assuming ";
      //}
      string directionCode = CDirectionReferenceCode::name(dirRefCode);
      //cout << directionCode << endl;

      vector<double> position ;
      if (r->isPositionExists()){
	position = DConverter::toVectorD<Length>(r->getPosition());
      } 
				
      vector<double> properMotion = DConverter::toVectorD(r->getProperMotion());
  
      vector<string> transition;
      if (r->isTransitionExists()) {
	transition = r->getTransition();
      }

      vector<double> restFrequency;
      if (r->isRestFrequencyExists()) {
	restFrequency = DConverter::toVectorD<Frequency>(r->getRestFrequency());
      }

      vector<double> sysVel;
      if (r->isSysVelExists()) {
	sysVel = DConverter::toVectorD<Speed>(r->getSysVel());
      }
   
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addSource( sourceId, time, interval, spectralWindowId, numLines,
                                 sourceName, calibrationGroup, code, direction,
                                 directionCode, position, properMotion, transition,
                                 restFrequency, sysVel);
      }
    }
    
    if (nSource) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->source().nrow() <<" source(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (IllegalAccessException& e) {
    errstream.str("");
    error(errstream.str());
  } catch (ConversionException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }

  //
  // Process the SysCal table.
  //
  const SysCalTable& sysCalT = ds->getSysCal();
  try {
    SysCalRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << sysCalT.size() << " sysCal(s)...";
    rowsInAScanbyTimeIntervalFunctor<SysCalRow> selector(selectedScanRow_v);

    const vector<SysCalRow *>& v = selector(sysCalT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";

    info(infostream.str());
    int nSysCal = v.size();

    for (int i = 0; i < nSysCal; i++) {
      r = v.at(i);
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond ;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0 ;
      //}

      pair<bool, bool> tcal_flag_pair;
      tcal_flag_pair.first   = r->isTcalFlagExists();
      tcal_flag_pair.second  = r->isTcalFlagExists() ? r->getTcalFlag() : false;

      pair<bool, vector<float> > tcal_spectrum_pair;
      tcal_spectrum_pair.first  =  r->isTcalSpectrumExists() ;
      if (tcal_spectrum_pair.first)
	tcal_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTcalSpectrum(), true);

      pair<bool, bool> trx_flag_pair;
      trx_flag_pair.first   = r->isTrxFlagExists();
      trx_flag_pair.second  = r->isTrxFlagExists() ? r->getTrxFlag() : false;

      pair<bool, vector<float> > trx_spectrum_pair;
      trx_spectrum_pair.first  =  r->isTrxSpectrumExists() ;
      if (trx_spectrum_pair.first)
	trx_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTrxSpectrum(), true);

      pair<bool, bool> tsky_flag_pair;
      tsky_flag_pair.first   = r->isTskyFlagExists();
      tsky_flag_pair.second  = r->isTskyFlagExists() ? r->getTskyFlag() : false;

      pair<bool, vector<float> > tsky_spectrum_pair;
      tsky_spectrum_pair.first  =  r->isTskySpectrumExists() ;
      if (tsky_spectrum_pair.first)
	tsky_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTskySpectrum(), true);

      pair<bool, bool> tsys_flag_pair;
      tsys_flag_pair.first   = r->isTsysFlagExists();
      tsys_flag_pair.second  = r->isTsysFlagExists() ? r->getTsysFlag() : false;

      pair<bool, vector<float> > tsys_spectrum_pair;
      tsys_spectrum_pair.first  =  r->isTsysSpectrumExists() ;
      if (tsys_spectrum_pair.first)
	tsys_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTsysSpectrum(), true);

      pair<bool, bool> tant_flag_pair;
      tant_flag_pair.first   = r->isTantFlagExists();
      tant_flag_pair.second  = r->isTantFlagExists() ? r->getTantFlag() : false;

      pair<bool, vector<float> > tant_spectrum_pair;
      tant_spectrum_pair.first  =  r->isTantSpectrumExists() ;
      if (tant_spectrum_pair.first)
	tant_spectrum_pair.second = FConverter::toVectorF(r->getTantSpectrum(), true);

      pair<bool, bool> tant_tsys_flag_pair;
      tant_tsys_flag_pair.first   = r->isTantTsysFlagExists();
      tant_tsys_flag_pair.second  = r->isTantTsysFlagExists() ? r->getTantTsysFlag() : false;

      pair<bool, vector<float> > tant_tsys_spectrum_pair;
      tant_tsys_spectrum_pair.first  =  r->isTantTsysSpectrumExists() ;
      if (tant_tsys_spectrum_pair.first)
	tant_tsys_spectrum_pair.second = FConverter::toVectorF(r->getTantTsysSpectrum(), true);

      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	iter->second->addSysCal( (int) r->getAntennaId().getTagValue(), (int) r->getFeedId(),
				 (int) swIdx2Idx[r->getSpectralWindowId().getTagValue()],
                                 time, interval, r->getNumReceptor(), r->getNumChan(),
                                 tcal_spectrum_pair, tcal_flag_pair, trx_spectrum_pair,
                                 trx_flag_pair, tsky_spectrum_pair, tsky_flag_pair,
                                 tsys_spectrum_pair, tsys_flag_pair, tant_spectrum_pair,
                                 tant_flag_pair, tant_tsys_spectrum_pair, tant_tsys_flag_pair);				
      }
    }
    if (nSysCal) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->sysCal().nrow() <<" sysCal(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  
  //
  // Process the CalDevice table.
  try {
    const CalDeviceTable& calDeviceT = ds->getCalDevice();
    infostream.str("");
    infostream << "The dataset has " << calDeviceT.size() << " calDevice(s)...";

    if (processCalDevice && calDeviceT.size() > 0) {
      rowsInAScanbyTimeIntervalFunctor<CalDeviceRow> selector(selectedScanRow_v);    
      const vector<CalDeviceRow *>& calDevices = selector(calDeviceT.get(), ignoreTime);
      if (!ignoreTime) 
	infostream << calDevices.size() << " of them in the selected exec blocks / scans ... ";
      
      info(infostream.str());
      
      for (vector<CalDeviceRow*>::const_iterator iter = calDevices.begin(); iter != calDevices.end(); iter++) {
	bool ignoreThisRow = false;
	unsigned int numCalload = 0;
	unsigned int numReceptor = 0;
	
	//
	// Let's make some checks on the attributes.
	errstream.str("");
	infostream.str("");

	//
	// Is numCalload > 0 ?
	if ((numCalload = (*iter)->getNumCalload()) <= 0) { 
	  errstream << "In the table CalDevice, the attribute 'numCalload' in row #"
		    << (unsigned int) (iter - calDevices.begin())
		    << " has an invalid value '("
		    << numCalload << "'), a strictly positive value is expected."
		    << endl; 
	  ignoreThisRow = true;
	}
	
	//
	// Do we have enough elements in calLoadNames ?
	vector<CalibrationDeviceMod::CalibrationDevice> temp = (*iter)->getCalLoadNames();
	vector<string> calLoadNames;
	if (temp.size() < numCalload) { 
	  errstream  << "In the table CalDevice, the size of the attribute 'calLoadNames' in row #"
		     << (unsigned int) (iter - calDevices.begin())
		     << " is too small. It should be greater than or equal to the value of the atttribute 'numCalload' ("
		     << numCalload
		     <<")."
		     << endl;
	  ignoreThisRow = true;
	} else {
	  calLoadNames.resize(temp.size());
	  transform(temp.begin(), temp.end(), calLoadNames.begin(), stringValue<CalibrationDeviceMod::CalibrationDevice, CCalibrationDevice>);	  
	}
	
	//
	// Do we have numReceptor ?
	if ((*iter)->isNumReceptorExists()) {
	  numReceptor = (*iter)->getNumReceptor();
	  if (numReceptor == 0) {
	    errstream << "In the table CalDevice, the value of the attribute 'numReceptor' in row #"
		      << (unsigned int) (iter - calDevices.begin())
		      << " is invalid (" 
		      << numReceptor 
		      << "). It is expected to be strictly positive."
		      << endl;
	    ignoreThisRow = true;
	  }
	}
	
	//
	// Do we have calEff ?
	vector<vector<float> > calEff;
	if ((*iter)->isCalEffExists()) {
	  //
	  // Do we take it into account ?
	  if (numReceptor == 0) {
	    infostream << "In the table CalDevice, the attribute 'calEff' is present in row #"
		       << (unsigned int) (iter - calDevices.begin())
		       << " but it will be ignored due to the fact that the attribute 'numReceptor' is null."
		       << endl;
	  } else {
	    calEff = (*iter)->getCalEff();
	  
	    //
	    // Let's check the sizes of its content.
	    if (calEff.size() < numReceptor) {
	      errstream << "In the table CalDevice, the size of the attribute 'calEff' in row #"
			<< (unsigned int) (iter - calDevices.begin())
			<< " is too small. It should be greater than or equal to the value of the attribute 'numReceptor' ("
			<< numReceptor
			<<")."
			<< endl;
	      ignoreThisRow = true;
	    } else {
	      if (find_if(calEff.begin(), calEff.end(), size_lt<float>(numCalload)) != calEff.end()) {
		errstream << "In the table CalDevice, the attribute 'calEff' in row #"
			  << (unsigned int) (iter - calDevices.begin())
			  << " has at least one element whose size is too small. All its elements should have their size"
			  << " greater then  or equal to the value of the attribute 'numCalload' ("
			  << numCalload
			  << ")."
			  << endl;
		ignoreThisRow = true;
	      }
	    }
	  }	
	}
      
	//
	// In priority let's see if we have coupledNoiseCal ?
	vector<vector<float> > coupledNoiseCal;
	if ((*iter)->isCoupledNoiseCalExists()) {
	  //
	  // Do we take it into account ?
	  if (numReceptor == 0) {
	    infostream << "In the table CalDevice, the attribute 'coupledNoiseCal' is present in row #"
		       << (unsigned int) (iter - calDevices.begin())
		       << " but it will be ignored due to the fact that the attribute 'numReceptor' is null."
		       << endl;
	  } else {
	    coupledNoiseCal = (*iter)->getCoupledNoiseCal();
	  
	    //
	    // Let's check the sizes of its content.
	    if (coupledNoiseCal.size() < numReceptor) {
	      errstream << "In the table CalDevice, the size of the attribute 'coupledNoiseCal' in row #"
			<< (unsigned int) (iter - calDevices.begin())
			<< " is too small. It should be greater than or equal to the value of the attribute 'numReceptor' ("
			<< numReceptor
			<<")."
			<< endl;
	      ignoreThisRow = true;
	    } else {
	      if (find_if(coupledNoiseCal.begin(), coupledNoiseCal.end(), size_lt<float>(numCalload)) != coupledNoiseCal.end()) {
		errstream << "In the table CalDevice, the attribute 'coupledNoiseCal' in row #"
			  << (unsigned int) (iter - calDevices.begin())
			  << " has at least one element whose size is too small. All its elements should have their size"
			  << " greater than or equal to the value of the attribute 'numCalload' (=="
			  << numCalload
			  << ")."
			  << endl;
		ignoreThisRow = true;
	      }
	    }
	  }	
	}
	// Ok we don't have coupledNoiseCal , but maybe we have noiseCal ?
	else if ((*iter)->isNoiseCalExists()) {
	  //
	  // Do we take it into account ?
	  vector<double> noiseCal = (*iter)->getNoiseCal();
	
	  if (noiseCal.size() < numCalload) {
	    infostream << "In the table CalDevice, the size of the attribute 'noiseCal' in row #"
		       << (unsigned int) (iter - calDevices.begin())
		       << " is too small. It should be greater than or equal to the value of the attribute 'numCalload' ("
		       << numCalload
		       << ")."
		       << endl;
	    ignoreThisRow = true;
	  } else {
	    // So yes we have a noiseCal attribute, then pretend we have coupledNoiseCal. 
	    // Artificially force numReceptor to 2 and fill coupledNoiseCal by replicating what we have in noiseCal :
	    // coupledNoiseCal[0] = noiseCal
	    // coupledNoiseCal[1] = noiseCal
	    //
	    // infostream << "In the table CalDevice  there is no attribute 'coupledNoiseCal' but there an attribute 'noiseCal' in row #"
	    // 	       << (unsigned int) (iter - calDevices.begin())
	    // 	       << " which we are going to use to fill the MS NOISE_CAL by replicating its values."
	    // 	       << endl;
	  
	    numReceptor = 2;
	    coupledNoiseCal.resize(numReceptor);
	    for (unsigned int iReceptor = 0; iReceptor < numReceptor; iReceptor++) {
	      coupledNoiseCal[iReceptor].resize(numCalload);
	      transform(noiseCal.begin(), noiseCal.begin()+numCalload, coupledNoiseCal[iReceptor].begin(), d2f);
	    } 
	  }
	}
      
	//
	// Do we have temperatureLoad ?
	vector<double> temperatureLoad;
	if ((*iter)->isTemperatureLoadExists()) {
	  vector<Temperature> temp = (*iter)->getTemperatureLoad();
	  if (temp.size() < numCalload) {
	    errstream  << "In the table CalDevice, the size of the attribute 'temperatureLoad' in row #"
		       << (unsigned int) (iter - calDevices.begin())
		       << " is too small. It should be greater than or equal to the value of the atttribute 'numCalload' ("
		       << numCalload
		       <<")."
		       << endl;
	    ignoreThisRow = true;
	  } else {
	    temperatureLoad.resize(temp.size());
	    transform(temp.begin(), temp.end(), temperatureLoad.begin(), basicTypeValue<Temperature, float>);	  
	  }
	}
      
	if (errstream.str().size() > 0) 
	  error(errstream.str());
      
	if (infostream.str().size() > 0)
	  info(infostream.str());

	if (ignoreThisRow) {
	  infostream.str("");
	  infostream << "This row will be ignored." << endl;
	  info(infostream.str());
	  continue;
	}

	//
	// And finally we can add a new row to the MS CALDEVICE table.
	double interval = ((double) (*iter)->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
	double time;
	// if (isEVLA) {
	//   time =  (*iter)->getTimeInterval().getStartInMJD()*86400 ;
	// }
	// else {
	time =  (*iter)->getTimeInterval().getStartInMJD()*86400 + interval / 2.0 ;
	//}
      
	
	for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator msIter = msFillers.begin();
	     msIter != msFillers.end(); ++msIter) {
	  msIter->second->addCalDevice( (*iter)->getAntennaId().getTagValue(), (*iter)->getFeedId(),
                                        swIdx2Idx[(*iter)->getSpectralWindowId().getTagValue()],
                                        time, interval, numCalload, calLoadNames, numReceptor,
                                        calEff, coupledNoiseCal, temperatureLoad);
	}      
      }

      unsigned int numMSCalDevices = (const_cast<casacore::MeasurementSet*>(msFillers.begin()->second->ms()))->rwKeywordSet().asTable("CALDEVICE").nrow();
      if (numMSCalDevices > 0) {
	infostream.str("");
	infostream << "converted in " << numMSCalDevices << " caldevice(s) in the measurement set.";
	info(infostream.str());
      }
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
 
  //
  // Process the SysPower table.
  if ( processSysPower )
    fillSysPower(dsName, ds, ignoreTime, selectedScanRow_v, msFillers);

  //
  // Load the weather table
  const WeatherTable& weatherT = ds->getWeather();

  try {
    WeatherRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << weatherT.size() << " weather(s)...";
    rowsInAScanbyTimeIntervalFunctor<WeatherRow> selector(selectedScanRow_v);
    
    const vector<WeatherRow *>& v = selector(weatherT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";

    info(infostream.str());
    int nWeather = v.size();

    infostream.str("");
    infostream << "The dataset has " << nWeather << " weather(s)...";
    info(infostream.str());
    
    pair<bool, float>
      pressureOpt,
      relHumidityOpt,
      temperatureOpt,
      windDirectionOpt,
      windSpeedOpt,
      dewPointOpt;
    
#define OPT_ATTR_PAIR( rowPtr, AttributeName ) rowPtr -> is ## AttributeName ## Exists() ? make_pair ( true, rowPtr -> get ## AttributeName ().get()) : make_pair( false, 0.)        
    
    for (int i = 0; i < nWeather; i++) {
      r			 = v.at(i);      
      double	interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double	time	 = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0;
      
      pressureOpt			   = OPT_ATTR_PAIR(r, Pressure);
      pressureOpt.second		  /= 100. ;	// We consider that ASDM stores Pascals & MS expects hectoPascals
      relHumidityOpt			   = OPT_ATTR_PAIR(r, RelHumidity);
      temperatureOpt			   = OPT_ATTR_PAIR(r, Temperature);
      windDirectionOpt			   = OPT_ATTR_PAIR(r, WindDirection);
      windSpeedOpt			   = OPT_ATTR_PAIR(r, WindSpeed);
      dewPointOpt			   = OPT_ATTR_PAIR(r, DewPoint);
      int		wxStationId        = r->getStationId().getTagValue();
      vector<double>	wxStationPosition  = DConverter::toVectorD(r->getStationUsingStationId()->getPosition());
    
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
          iter->second->addWeather( -1, time, interval, pressureOpt, relHumidityOpt, temperatureOpt,
                                    windDirectionOpt, windSpeedOpt, dewPointOpt, wxStationId, wxStationPosition);
      }
    }

    if (nWeather) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->weather().nrow() <<" weather(s) in the measurement set." ;
      info(infostream.str());
    }
  } catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  } catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  // And then finally process the state and the main table.
  //
  if (lazy) {
    fillMainLazily(dsName, ds, selected_eb_scan_m,effectiveBwPerDD_m,e_query_cm,checkdupints);
  } else {

    const MainTable&		mainT  = ds->getMain();
    const StateTable&		stateT = ds->getState();

    vector<MainRow*>			v;
    vector<int32_t>			mainRowIndex; 
    //
    //
    // Consider only the Main rows whose execBlockId and scanNumber attributes correspond to the selection.
    // (execBlockId, scanNumber, wvr-corrected-data option)
    //
    vector<AtmPhaseCorrectionMod::AtmPhaseCorrection>		queriedAPC_v						  = es_apc.toEnumType();
    const vector<MainRow *>&		temp							  = mainT.get();
    for ( vector<MainRow *>::const_iterator iter_v = temp.begin(); iter_v			 != temp.end(); iter_v++) {
      map<int, set<int> >::iterator iter_m = selected_eb_scan_m.find((*iter_v)->getExecBlockId().getTagValue());
      if ( iter_m != selected_eb_scan_m.end() && iter_m->second.find((*iter_v)->getScanNumber()) != iter_m->second.end() ) {
	mainRowIndex.push_back(iter_v - temp.begin());
	v.push_back(*iter_v);
      }
    }
      
    infostream.str("");
    infostream << "The dataset has " << mainT.size() << " main(s)...";
    infostream << v.size() << " of them in the selected exec blocks / scans." << endl;
    info(infostream.str());

    MSMainRowsInSubscanChecker msMainRowsInSubscanChecker;

    unsigned int  nMain = v.size();      
    const VMSData *vmsDataPtr = 0;

    // Initialize an UVW coordinates engine.
    UvwCoords uvwCoords(ds);
      
    ostringstream oss;
    EnumSet<AtmPhaseCorrectionMod::AtmPhaseCorrection> es_query_ap_uncorrected;
    es_query_ap_uncorrected.fromString("AP_UNCORRECTED");
    
    // used in checking for duplicate integrations in the WVR (Radiometer) case
    // This holds the most recent last integration time for each configDescriptionId - but only for Radiometer data.
    map<Tag, double> lastTimeMap;
    // for debugging - to report on what uid the lastTime being compared came from
    // map<Tag, string> lastTimeUIDMap;

    for (unsigned int i = 0; i < nMain; i++) {
      try {
	// What's the processor for this Main row ?
	Tag cdId = v[i]->getConfigDescriptionId();
	ConfigDescriptionTable& cT = ds->getConfigDescription();
	ConfigDescriptionRow* cR = cT.getRowByKey(cdId);
	Tag pId = cR->getProcessorId();


        ProcessorTypeMod::ProcessorType processorType = ds->getProcessor().getRowByKey(pId)->getProcessorType();
	infostream.str("");
	infostream << "ASDM Main row #" << mainRowIndex[i] << " contains data produced by a '" << CProcessorType::name(processorType) << "'." ;
	info(infostream.str());

	string dataUID = v[i]->getDataUID().getEntityId().toString();
	replace(dataUID.begin(),dataUID.end(),':','_');
	replace(dataUID.begin(),dataUID.end(),'/','_');
	string absBDFpath = Path(dsName + "/ASDMBinary/" + dataUID).absoluteName();
	infostream.str("");
	infostream << "ASDM Main row #" << mainRowIndex[i]
		   << " (scan #" << v[i]->getScanNumber()
		   <<", subscan #" <<  v[i]->getSubscanNumber()
		   <<", " << v[i]->getConfigDescriptionId().toString() << ")"
		   << " - BDF file '" << absBDFpath << "' - Size is " << v[i]->getDataSize() << " bytes for " << v[i]->getNumIntegration() << " integrations." << endl;
	info(infostream.str());

        if(v[i]->getNumIntegration()==0 ||v[i]->getDataSize()==0) {
	  infostream.str("");
          infostream << "No valid data in this BDF. The main row # "<<  mainRowIndex[i] << " is ignored." ;
          info(infostream.str());
          continue;
        } 

	// Populate the State table.
	try {
	  fillState(v[i]);
	} catch (ASDM2MSException e) {
	  // The State table could not be filled, then let's forget this Main table row.
	  infostream.str("");
	  infostream << e.getMessage() << ". The main row # " <<  mainRowIndex[i] << " is ignored.";
	  info(infostream.str());
	  continue;	  
	}

	if (processorType == ProcessorTypeMod::RADIOMETER) {
	  if (!sdmBinData.acceptMainRow(v[i])) {
	    infostream.str("");
	    infostream <<"No data retrieved in the Main row #" << mainRowIndex[i] << " (" << sdmBinData.reasonToReject(v[i]) <<")" << endl;
	    info(infostream.str());
	    continue;
	  }
	  vmsDataPtr = sdmBinData.getDataCols();
	  bool skipFirstIntegration = checkdupints && lastTimeMap[cdId] == vmsDataPtr->v_time[0];
	  unsigned int skipValues = 0;
	  // useful just for debugging
	  if (skipFirstIntegration) {
	    // work out the number of elements in this first row
	    // possibly this is available reliably elsewhere, but thats' not obvious
	    // this should be rare and only a few elements should tested, so relatively quick
	    bool firstIntegration = True;
	    while (firstIntegration && skipValues < vmsDataPtr->v_time.size()) {
	      firstIntegration = lastTimeMap[cdId] == vmsDataPtr->v_time[skipValues];
	      if (firstIntegration) skipValues++;
	    }
	    if (debug) {
	      cout << "Duplicate time seen in Row : " << i
		   <<  " cdId : " << cdId 
		   << " " << v[i]->getDataUID().getEntityId().toString()
		//   << " duplicates time at end of " << lastTimeUIDMap[cdId]
		   << " time size " << vmsDataPtr->v_time.size()
		   << " antennaId1 size " << vmsDataPtr->v_antennaId1.size()
		   << " skipValues : " << skipValues
		   << endl;
	    }
	  }

	  lastTimeMap[cdId] = vmsDataPtr->v_time[vmsDataPtr->v_time.size()-1];
	  // lastTimeUIDMap[cdId] = v[i]->getDataUID().getEntityId().toString();
	  
	  fillMain( v[i], sdmBinData, vmsDataPtr, uvwCoords, effectiveBwPerDD_m,
                    complexData, mute, ac_xc_per_timestamp, skipValues);
          
	  infostream.str("");
	  infostream << "ASDM Main row #" << mainRowIndex[i] << " produced a total of " << (vmsDataPtr->v_antennaId1.size()-skipValues) << " MS Main rows." << endl;
	  info(infostream.str());
	} else { // Assume we are in front of a Correlator.
	  // Open its associate BDF.

	  bool rowOK = sdmBinData.openMainRow(v[i]);
	  if (!rowOK) {
	    infostream.str("");
	    infostream << "No data retrieved in the Main row #" << mainRowIndex[i] << " (" << sdmBinData.reasonToReject(v[i]) << ")" << endl;
	    info(infostream.str());
	    continue;
	  }
	  
	  int N = v[i]->getNumIntegration();
	  uint64_t  bdfSize = v[i]->getDataSize();
	  vector<unsigned int> integrationSlices(getIntegrationSlices(N, bdfSize, bdfSliceSizeInMb*1024*1024));
	  int32_t			numberOfMSMainRows	 = 0;
	  int32_t			numberOfIntegrations	 = 0;
	  int32_t			numberOfReadIntegrations = 0;

	  // For each slice of the BDF with a size approx equal to the required size
	  for (unsigned int j = 0; j < integrationSlices.size(); j++) {
	    numberOfIntegrations = integrationSlices[j];
	    if (numberOfIntegrations) {
	      infostream.str("");
	      infostream << "ASDM Main row #" << mainRowIndex[i] << " - " << numberOfReadIntegrations  << " integrations done so far - the next " << numberOfIntegrations << " integrations produced " ;
	      vmsDataPtr = sdmBinData.getNextMSMainCols(numberOfIntegrations);

	      msMainRowsInSubscanChecker.check(vmsDataPtr, v[i], mainRowIndex[i], absBDFpath);
	      numberOfReadIntegrations += numberOfIntegrations;
	      numberOfMSMainRows += vmsDataPtr->v_antennaId1.size();
	      fillMain(v[i], sdmBinData, vmsDataPtr, uvwCoords, effectiveBwPerDD_m, complexData,  mute, ac_xc_per_timestamp);
	      infostream << vmsDataPtr->v_antennaId1.size()  << " MS Main rows." << endl;
	      info(infostream.str());
	    }
	  }
	  
	  // this should no longer be necessary, but keep this block in place just in case.
	  uint32_t numberOfRemainingIntegrations = N - numberOfReadIntegrations;
	  if (numberOfRemainingIntegrations) { 
	    vmsDataPtr = sdmBinData.getNextMSMainCols(numberOfRemainingIntegrations);

	    if (vmsDataPtr != NULL && vmsDataPtr->v_antennaId1.size() > 0) {
	      infostream.str("");
	      infostream << "ASDM Main row #" << mainRowIndex[i] << " - " << numberOfReadIntegrations  << " integrations done so far - the next " << numberOfRemainingIntegrations << " integrations produced " ;
	      
	      msMainRowsInSubscanChecker.check(vmsDataPtr, v[i], mainRowIndex[i], absBDFpath);
	      fillMain(v[i], sdmBinData, vmsDataPtr, uvwCoords, effectiveBwPerDD_m, complexData, mute, ac_xc_per_timestamp);
	      
	      infostream << vmsDataPtr->v_antennaId1.size()  << " MS Main rows." << endl;
	      info(infostream.str());
	      numberOfMSMainRows += vmsDataPtr->v_antennaId1.size();
	      infostream.str("");
	      infostream << "ASDM Main row #" << mainRowIndex[i] << "produced a total of " << numberOfMSMainRows << " MS Main rows." << endl;
	    }
	  }
	}
      } catch ( ConversionException& e) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      } catch ( IllegalAccessException& e) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      } catch ( SDMDataObjectParserException& e) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      } catch ( SDMDataObjectStreamReaderException& e ) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      } catch ( SDMDataObjectReaderException& e ) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      } catch (ASDM2MSException& e) {
	infostream.str("");
	infostream << e.getMessage();
	info(infostream.str());
      }
      /*
	catch ( std::exception & e) {
	infostream.str("");
	infostream << e.what();
	info(infostream.str());      
	}
      */
      catch (Error & e) {
	infostream.str("");
	infostream << e.getErrorMessage();
	info(infostream.str());
      }
    }
  
    
    // Did we have problem with BDF with data not falling in the time range of their scan/subscan pair ?
    const vector<string>& report = msMainRowsInSubscanChecker.report();
    for_each(report.begin(), report.end(), warning ); 
    
    infostream.str("");
    infostream << "The dataset has "  << stateT.size() << " state(s)..." ;
    info(infostream.str());
    
    if (stateT.size()) {
      infostream.str("");
      infostream << "converted in " << msFiller->ms()->state().nrow() << " state(s) in the measurement set.";
      info(infostream.str());
    }
    
    infostream.str("");
    infostream << "The dataset has " << mainT.size() << " main(s)...";
    info(infostream.str());
    
    if (mainT.size()) {
        for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end(); ++iter) {
	string kindOfData = (iter->first == AtmPhaseCorrectionMod::AP_UNCORRECTED) ? "wvr uncorrected" : "wvr corrected";
	infostream.str("");
	infostream << "converted in " << iter->second->ms()->nrow() << " main(s) rows in the measurement set containing the " << kindOfData << " data.";
	info(infostream.str());
      }
    }
  }
    
  // Do we also want to store the verbatim copies of some tables of the ASDM dataset ?
  if (asisOption.size() > 0) {
    try {
      istringstream iss;
      iss.str(asisOption);
      string word;
      vector<string> tablenames;
      while (iss>>word)
	tablenames.push_back(word);
      for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin(); iter != msFillers.end(); ++iter){   
	ASDMVerbatimFiller avf(const_cast<casacore::MS*>(iter->second->ms()), Name2Table::find(tablenames, verbose));
	avf.fill(*ds);
      }
    } catch (ConversionException e) {
      infostream.str("");
      infostream << e.getMessage();
      info(infostream.str());
    }
  }
  
  for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
       iter != msFillers.end(); ++iter)
    iter->second->end();
  
  
  for (map<AtmPhaseCorrectionMod::AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
       iter != msFillers.end(); ++iter)
    delete iter->second;
  
  delete ds;
  return 0;
}
