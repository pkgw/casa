/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version January 2010.
   Maintained by ESO since 2013. 
   
   This file is part of LibAIR and is licensed under GNU Public
   License Version 2
   
   \file wvrgcal.cpp

   This is a command line program that reads WVR data from a CASA
   measurement set, computes the predicted complex gain of each
   antenna as a function of time from these WVR data and then writes
   these solutions out to a CASA gain table. 
   
*/


#include <iostream>
#include <numeric>
#include <string>

#include <stdcasa/optionparser.h>
#include <ms/MeasurementSets/MeasurementSet.h>

#include "../casawvr/mswvrdata.hpp"
#include "../casawvr/msgaintable.hpp"
#include "../casawvr/msutils.hpp"
#include "../casawvr/msspec.hpp"
#include "../casawvr/msantdata.hpp"
#include "../src/apps/arraydata.hpp"
#include "../src/apps/arraygains.hpp"
#include "../src/apps/almaabs.hpp"
#include "../src/apps/dtdlcoeffs.hpp"
#include "../src/apps/almaresults.hpp"
#include "../src/apps/segmentation.hpp"
#include "../src/libair_main.hpp"

#include "wvrgcalerrors.hpp"
#include "wvrgcalfeedback.hpp"
#include "wvrgcalargs.hpp"

using LibAIR2::fatalMsg;
using LibAIR2::errorMsg;
using LibAIR2::warnMsg;


/// Check the options and parameters supplied by the user for internal
/// consistency 
bool checkPars( const option::Option vm[] )
{
  if (vm[wvr::arg::index::MS].count( ) <1)
  {
    fatalMsg("No input measurement sets given -- aborting ");
    return true;
  }

  if (vm[wvr::arg::index::OUTPUT].count( ) <1)
  {
    fatalMsg("No output file give -- aborting ");
    return true;
  }

  if (vm[wvr::arg::index::SEGFIELD].count( ) )
  {
    fatalMsg("The --segfield option has been removed because it does not handle mosaic observations well.");
    fatalMsg("Please use the --segsource option instead");
    return true;
  }

  if (vm[wvr::arg::index::SEGSOURCE].count( ) && vm[wvr::arg::index::CONT].count( ))
  {
    fatalMsg("Multiple retrievals using continuum estimation not yet supported");
    fatalMsg("Use only one of  --segfield OR --cont");
    return true;
  }

  if (vm[wvr::arg::index::NSOL].count( ) && vm[wvr::arg::index::CONT].count( ))
  {
    try {
        if ( std::stoi(vm[wvr::arg::index::NSOL].last( )->arg) > 1) {
            fatalMsg("Multiple retrievals using continuum estimation not yet supported");
            fatalMsg("You can not use the nsol parameter yet with cont");
            return true;
        }
    } catch (std::exception const &e) {
        fatalMsg("You can not use the nsol parameter yet with cont");
        return true;
    }
    
  }

  if (vm[wvr::arg::index::SOURCEFLAG].count( ) && !vm[wvr::arg::index::SEGSOURCE].count( ))
  {
    fatalMsg("Can only flag a source using --sourceflag if the --segsource option is also used");
    fatalMsg("Please either remove the --sourceflag option or also specify the --segsource option");
    return true;
  }

  if (vm[wvr::arg::index::TIE].count( ) && !vm[wvr::arg::index::SEGSOURCE].count( ))
  {
    fatalMsg("Can only tie sources together if the --segsource option is also used");
    fatalMsg("Please either remove the --tie option or also specify the --segsource option");
    return true;
  }

  if (vm[wvr::arg::index::REVERSE].count( ) and vm[wvr::arg::index::REVERSESPW].count( ))
  {
    warnMsg("You are specifying both the reverse and reversespw options;"
	    "the latter will be ignored and all spectral windows will be reversed");
  }

  if (vm[wvr::arg::index::SMOOTH].count( ))
  {
      try {
          if ( std::stod(vm[wvr::arg::index::SMOOTH].last( )->arg) < 1 ) {
              warnMsg("Smooth parameter must be 1 or greater");
              return true;
          }
      } catch (std::exception const &e) { }
  }

  if ( vm[wvr::arg::index::MAXDISTM].count( ) )
  {
      try {
          if ( std::stod(vm[wvr::arg::index::MAXDISTM].last( )->arg) < 1 ) {
              warnMsg("maxdistm parameter must be 0. or greater");
              return true;
          }
      } catch (std::exception const &e) { }
  }

  if (vm[wvr::arg::index::OFFSETS].count( ))
  {
      std::string offsetstable=vm[wvr::arg::index::OFFSETS].last( )->arg;
      try{
          //casacore::Table f(offsetstable);
      } catch (std::exception const &e) {
//    catch(const casacore::AipsError rE){
          //fatalMsg(rE.getMesg());
          fatalMsg("The --offsets option needs to point to an existing CASA Table containing temperature offsets.");
          return true;
      }
  }
  return false;
  
}

void checkWarnPars(const option::Option vm[] )
{

  if (vm[wvr::arg::index::STATFIELD].count( ))
  {
    warnMsg("The use of \"statfield\" is not recommended as mosaiced"
	    "observations are recorded as many separate fields. Recomended option"
	    "to use is \"statsource\". ");
  }
}

/* Checks on parameter that can only be done once the measurement set
   is opened.  Do not put computationally intensive checks otherwise
   feedback to the user will be too slow.
*/
void checkMSandPars(const casacore::MeasurementSet &ms, const option::Option vm[])
{
  if (vm[wvr::arg::index::STATSOURCE].count( ))
  {
    std::string srcname=vm[wvr::arg::index::STATSOURCE].last( )->arg;
    std::set<size_t> fselect=LibAIR2::getSrcFields(ms,
						  srcname);
    if (fselect.size() == 0)
    {
      std::cout<<"WARNING: No Source table entries appear to be identified with source "<< srcname <<std::endl
	       <<"         that you supplied to the statsource option"<<std::endl
	       <<"         Statistics will be corrupted and wvrgcal may fail"<<std::endl
	       <<std::endl;
    }
  }

}

/** \brief Take a parameter than can be specified as a sequence of
    antenna numbers or names and always return as a sequence of
    antenna numbers.
 */


						   
#if __cplusplus < 201103L
struct hack01 {
	bool operator()(bool acc, const std::map<size_t, std::string >::value_type &p){
		bool cmp = p.second == ele;
		if ( cmp ) match = p.first;
		return acc || cmp;
	}
	hack01( size_t &m, const std::string &e ) : match(m), ele(e) { }
	size_t &match;
	const std::string &ele;
};

struct hack02 {
	bool operator()(bool acc, const std::map<size_t, std::string >::value_type &p) { return acc || (p.first == n); }
	hack02( int x ) : n(x) { }
	int n;
};
#endif


LibAIR2::AntSet getAntPars(const std::string &s,
			   const option::Option *par,
			   const casacore::MeasurementSet &ms)
{
  using namespace LibAIR2;
  LibAIR2::AntSet res;
  aname_t anames=getAName(ms);
  for (auto p=par; p; p=p->next( ))
  {
	  size_t match;
	  if (std::accumulate( anames.begin( ),
						   anames.end( ), false,
#if __cplusplus >= 201103L
						   [&](bool acc, const aname_t::value_type &p){
							   bool cmp = p.second == par->arg;
							   if ( cmp ) match = p.first;
							   return acc || cmp;
						   }
#else
						   hack01(match,par)
#endif
))
    {
      res.insert(match);
    }
    else
    {
      // should be an antenna number
      try {
    int n=std::stoi(par->arg);
    fprintf( stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %d\n", n );
    if ( std::accumulate( anames.begin(),
                          anames.end(), false,
#if __cplusplus >= 201103L
                          [=](bool acc, const aname_t::value_type &p) {
                              return acc || ((int)p.first == n);
                          }
#else
						  hack02(n)
#endif
) == false )
    {
	  throw AntIDError(n,
			   anames);
	}
	res.insert(n);
      }
      catch (...)
      {
	throw AntIDError(par->arg,
			 anames);

      }

    }
  }


  return res;
}


/** Simple function to turn the command line parameters into a single
    string.

    Return by value as only needs to be called once in program
*/
static std::string buildCmdLine(int argc,
				char* argv[])
{
  std::string cmdline;
  for (size_t i=0; i< (size_t)argc; ++i)
  {
    cmdline+=std::string(argv[i]);
    cmdline+=" ";
  }
  return cmdline;
}


/* Determine the nearest n antennas not accepting ones which
   are flagged or have a distance > maxdist_m
*/
 
LibAIR2::AntSetWeight limitedNearestAnt(const LibAIR2::antpos_t &pos,
				       size_t i,
				       const LibAIR2::AntSet &flag,
				       size_t n,
				       double maxdist_m)
{
  LibAIR2::AntSetD dist=LibAIR2::antsDist(pos, i, flag);
  LibAIR2::AntSetWeight res;
    
  double total=0;
  size_t limitedn=0;
  LibAIR2::AntSetD::const_iterator s=dist.begin();
  for (size_t j=0; j<n; ++j)
  {
    if(s!=dist.end() and s->first <= maxdist_m)
    {
      total+=s->first;
      ++s;
      ++limitedn;
    }
  }

  s=dist.begin();
  for (size_t j=0; j<limitedn; ++j)
  {
    res.insert(std::make_pair(s->first/total, s->second));
    ++s;
  }

  return res;

}


/** \brief Flag and interpolate WVR data
 */
void flagInterp(const casacore::MeasurementSet &ms,
		const LibAIR2::AntSet &wvrflag,
		LibAIR2::InterpArrayData &d,
		const double maxdist_m,
		const int minnumants,
		LibAIR2::AntSet &interpImpossibleAnts)
{

  LibAIR2::antpos_t apos;
  LibAIR2::getAntPos(ms, apos);
  LibAIR2::AntSet wvrflag_s(wvrflag.begin(), 
			   wvrflag.end());

  for(LibAIR2::AntSet::const_iterator i=wvrflag.begin();
      i!=wvrflag.end(); 
      ++i)
  {

    LibAIR2::AntSetWeight near=limitedNearestAnt(apos, 
						 *i, 
						 wvrflag_s, 
						 3,
						 maxdist_m);
    if(near.size()>= static_cast<unsigned int>(minnumants)){
      //LibAIR2::interpBadAntW(d, *i, near);
      const LibAIR2::InterpArrayData::wvrdata_t &data(d.g_wvrdata());
      for(size_t ii=0; ii<d.g_time().size(); ++ii){
	for(size_t k=0; k < 4; ++k){
	  double p=0;
	  for(LibAIR2::AntSetWeight::const_iterator j=near.begin(); j!=near.end(); ++j){
        double thisData = data(ii,j->second,k);
	    if(thisData>0){
	      p+=thisData*j->first;
	    }
	    else{ // no good data; set solution to zero => will be flagged later
	      p=0.;
	      break;
	    }
	  }
	  d.set(ii, *i, k, p);
	}
      }
    }
    else
    { 
      std::ostringstream oss;
      oss << "Antenna " << *i 
	  << " has bad or no WVR and only " << near.size() << " near antennas ("
	  << maxdist_m << " m max. distance) to interpolate from. Required are " 
	  << minnumants << "." << std::endl;
      std::cout << std::endl << "*** " << oss.str();
      std::cerr << oss.str();
      for(size_t j=0; j<d.g_time().size(); ++j)
      {
        for(size_t k=0; k < 4; ++k)
        {
          d.set(j, *i, k, 0.); // set all WVR data of this antenna to zero => will be flagged later
        }
      }
      interpImpossibleAnts.insert(*i);
    }
  }
  
}


/// Work out which spectral windows might need to be reversed
std::set<size_t> reversedSPWs(const LibAIR2::MSSpec &sp,
                              const option::Option vm[])
{
  std::set<size_t> reverse;
  if (vm[wvr::arg::index::REVERSE].count( ))
  {
    for (size_t spw =0; spw<sp.spws.size(); ++spw)
      reverse.insert(spw);
  }    
  if (vm[wvr::arg::index::REVERSESPW].count( ))
  {
    auto revspws = vm[wvr::arg::index::REVERSESPW];
    for(option::Option *cur = revspws; cur; cur = cur->next( ))
    {
      int curi = std::stoi(cur->arg);
      if (curi < 0 or curi >= (int)sp.spws.size())
      {
	throw LibAIR2::SPWIDError(curi,sp.spws.size());
      }
      reverse.insert(curi);
    }
  }    
  return reverse;
}

void printExpectedPerf(const LibAIR2::ArrayGains &g,
		       const LibAIR2::dTdLCoeffsBase &coeffs,
		       const std::vector<std::pair<double, double> > &tmask)
{

  std::cout<<"  Expected performance "<<std::endl
	   <<"------------------------------------------------------------------"<<std::endl;
  
  std::vector<double> cr, err;
  coeffs.repr(cr, err);
  std::cout<<"* Estimated WVR thermal contribution to path fluctuations (micron per antenna): "
	   <<LibAIR2::thermal_error(cr)/1e-6
	   <<std::endl;
  const double grmsbl=g.greatestRMSBl(tmask);
  std::cout<<"* Greatest Estimated path fluctuation is (micron on a baseline): "
	   <<grmsbl/1e-6
	   <<std::endl;
  if (cr[1]>0.){
    std::cout<<"* Rough estimate path error due to coefficient error (micron on a baseline): "
	     <<grmsbl* (err[1]/cr[1])/1e-6
	     <<std::endl
	     <<std::endl;
  }
  else if(err[1]==0.) {
    std::cout<<"* Rough estimate path error due to coefficient error could not be calculated (is nominally zero)."<<std::endl;
  }
  else {
    std::cout<<"* Rough estimate path error due to coefficient error could not be calculated."<<std::endl;
    std::cerr<<"* Rough estimate path error due to coefficient error could not be calculated."<<std::endl;
  }
}

/** Compute the time intervals over which the statistics should be
    computed
 */
void statTimeMask(const casacore::MeasurementSet &ms,
		  const option::Option vm[],
		  std::vector<std::pair<double, double> > &tmask,
		  const std::vector<size_t> &sortedI,
		  const std::vector<int> &wvrspws)
{
  std::vector<int> flds;
  std::vector<double> time;
  std::vector<int> src;
  LibAIR2::fieldIDs(ms,
		   time,
		   flds,
		   src,
		   sortedI);
  std::vector<size_t> spws;
  LibAIR2::dataSPWs(ms, spws, sortedI);

  if (vm[wvr::arg::index::STATFIELD].count( ) == 0 && vm[wvr::arg::index::STATSOURCE].count( ) == 0)
  {
    tmask.resize(0);
    tmask.push_back(std::pair<double, double>(time[0], time[time.size()-1]));
  }
  else if ( vm[wvr::arg::index::STATSOURCE].count( ) > 0)
  {
    std::set<size_t> fselect=LibAIR2::getSrcFields(ms, vm[wvr::arg::index::STATSOURCE].last( )->arg);

    LibAIR2::fieldTimes( time, flds, spws, fselect, (size_t) wvrspws[0], tmask );
  }
  else
  {
    std::vector<std::string> fields;
    const option::Option *args = vm[wvr::arg::index::STATFIELD];
    for (auto arg = args; arg; arg = arg->next( )) fields.push_back(arg->arg);

    LibAIR2::field_t fnames=LibAIR2::getFieldNames(ms);

    std::set<size_t> fselect;
	size_t val;
    if (std::accumulate( fnames.begin( ),
						 fnames.end( ), false,
#if __cplusplus >= 201103L
						 [&](bool acc, const LibAIR2::field_t::value_type &p) {
							 bool cmp = p.second == fields[0];
							 if ( cmp ) val = p.first;
							 return acc || cmp;
						 }
#else
						 hack01(val,fields[0])
#endif
))  // User supplied  field *name*
    {
      fselect.insert(val);
    }
    else
    {
      try {
          size_t n=std::stoi(fields[0]);
          fselect.insert(n);
      } catch (...) {
          std::cout<<"Warning: Could not understand statfield argument. Will use zeroth field."
                   <<std::endl;
      }
    }
    LibAIR2::fieldTimes(time,
		       flds,
		       spws,
		       fselect,
		       (size_t) wvrspws[0],
		       tmask);
  }
  LibAIR2::printStatTimes(std::cout,
			 time,
			 tmask);
}


/// Compute the discrepance in path estimate between channels 1 and 3
void computePathDisc(const LibAIR2::InterpArrayData &d,
		     const std::vector<std::pair<double, double> > &tmask,
		     LibAIR2::dTdLCoeffsBase  &coeffs,
		     std::vector<double> &res)
{
  LibAIR2::ArrayGains g1(d.g_time(), 
			d.g_el(),
			d.g_state(),
			d.g_field(),
			d.g_source(),
			d.nAnts);
  std::array<double, 4> c1mask = {{0, 1, 0,0}};
  std::array<double, 4> c3mask = {{0, 0, 0,1}};
  std::array<double, 4> callmask = {{1, 1, 1,1}};
    
  coeffs.chmask=c1mask;
  g1.calc(d,
	  coeffs);    
  
  LibAIR2::ArrayGains g3(d.g_time(), 
			 d.g_el(),
			 d.g_state(),
			 d.g_field(),
			 d.g_source(),
			 d.nAnts);
  coeffs.chmask=c3mask;
  g3.calc(d,
	  coeffs);    

  g1.pathDiscAnt(g3, 
		 tmask,
		 res);

  coeffs.chmask=callmask;
  
}

void printFieldSegments(const std::vector<std::pair<double, double> >  &fb,
			double tbase)
{
  for (size_t i=0; i<fb.size(); ++i)
    {
      std::cout<<fb[i].first-tbase<<","<<fb[i].second-tbase<<std::endl;
    }

}

std::vector<std::set<std::string> > getTied(const option::Option vm[])
{
  std::vector<std::set<std::string> > res;

  //only run if --tie option given on command line
  if (vm[wvr::arg::index::TIE].count( ))
    {
    const option::Option *pars = vm[wvr::arg::index::TIE];
    for (auto par = pars; par; par = par->next( ))
      {
          std::set<std::string> cs;
          constexpr char sep[] = ",";
          char *buf = strdup(par->arg);
          for( char *arg = strtok(buf,sep); arg; arg = strtok(NULL, sep) ) cs.insert(arg);
          res.push_back(cs);
          free(buf);
      }
    }
    return res;
}

// Convert tied source names to tied source IDs
std::vector<std::set<size_t> >  tiedIDs(const std::vector<std::set<std::string> > &tied,
					const casacore::MeasurementSet &ms)
{
  std::map<size_t, std::string > srcmap=LibAIR2::getSourceNames(ms);
  std::vector<std::set<size_t> > res;
  for (size_t i=0; i<tied.size(); ++i)
  {
    std::set<size_t> cs;
    for(std::set<std::string>::const_iterator j=tied[i].begin();
	j!=tied[i].end();
	++j)
    {	
      try
      {
    int srcid=std::stoi(*j);
	std::map<size_t, std::string>::const_iterator it = srcmap.find(srcid);
	if(it == srcmap.end()) { // id does not exist
	  std::cerr << "Parameter 'tie': The source id " << *j << " is an integer but not a valid numerical Source ID. Will try to interpret it as a name ..." << std::endl;
	  throw std::exception();
	}
	cs.insert(srcid);
      }
      catch (const std::exception& x)
      {
        size_t match;
        if ( std::accumulate( srcmap.begin( ),
                              srcmap.end( ), false,
#if __cplusplus >= 201103L
                              [&](bool acc, const std::map<size_t, std::string >::value_type &p) {
                                  bool cmp = p.second == *j;
                                  if ( cmp ) match = p.first;
                                  return acc || cmp;
                              }
#else
							  hack01(match,*j)
#endif
)) {
          cs.insert(match);
        } else {
          std::ostringstream oss;
          oss << "Parameter 'tie': The field " << *j << " is not recognised. Please check for typos." << std::endl;
          throw LibAIR2::WVRUserError(oss.str());
		}
      }
    } // end for
    res.push_back(cs);
  }
  return res;
}

void printTied(const std::vector<std::set<std::string> > &tied,
	       const std::vector<std::set<size_t> > &tiedi)
{
  for(size_t i=0; i<tied.size(); ++i)
  {
    std::set<std::string>::const_iterator it=tied[i].begin();
    std::cout<<"Tying: " << *it;
    for(it++; it!=tied[i].end(); it++){
      std::cout<<" and "<<*it;
    }
    std::cout<<std::endl;
  }
  if (tied.size())
    std::cout<<"Tied sets as numerical source IDs:"<<std::endl;
  for(size_t i=0; i<tiedi.size(); ++i)
  {
    std::set<size_t>::const_iterator it=tiedi[i].begin();
    std::cout<<"Tying: " << *it;
    for(it++; it!=tiedi[i].end(); it++){
      std::cout<<" and "<<*it;
    }
    std::cout<<std::endl;
  }

}

/** Compute the set of source_ids corresponding to a vector of source
    names.
 */
std::set<size_t> sourceSet(const std::vector<std::string> &sources,
			   const casacore::MeasurementSet &ms)
{
  std::map<size_t, std::string > snames=LibAIR2::getSourceNames(ms);
  std::set<size_t> sset;
  for(size_t i=0; i<sources.size(); ++i) {
	size_t match;
	if (std::accumulate( snames.begin( ),
						 snames.end( ), false,
#if __cplusplus >= 201103L
						 [&](bool acc, const std::map<size_t, std::string >::value_type &p) {
							 bool cmp = p.second == sources[i];
							 if ( cmp ) match = p.first;
							 return acc || cmp;
						 }
#else
						 hack01(match,sources[i])
#endif
)) {
		sset.insert(match);
	}
  }
  return sset;
}
  

/** Filter the set of input WVR measurements to retrieve the
    coefficients from to exclude flagged sources

    This function takes and returns two containers: the first is the
    list of WVR values to be analysed; the second is a vector of time
    ranges to use.
    
    These time ranges must be filtered together with the inputs since
    they are not referenced directly to the inputs are simply assumed
    to be "row-synchronous".
    
 */
std::pair<LibAIR2::ALMAAbsInpL,  std::vector<std::pair<double, double> > >
filterInp(const LibAIR2::ALMAAbsInpL &inp,
	  const std::vector<std::pair<double, double> > &fb,
	  const std::vector<std::string> &sourceflag,
	  const casacore::MeasurementSet &ms)
{

  std::set<size_t> flagset=sourceSet(sourceflag, ms);

  LibAIR2::ALMAAbsInpL res;
  std::vector<std::pair<double, double> > rfb;
  size_t j=0;
  for(LibAIR2::ALMAAbsInpL::const_iterator i=inp.begin();
      i!=inp.end();
      ++i, ++j)
  {
    if (flagset.count(i->source) ==0)
    {
      res.push_back(*i);
      rfb.push_back(fb[j]);
    }
  }
  return std::make_pair(res, rfb);
}

/** Filter the set of input WVR measurements to retrieve the
    coefficients from to exclude flagged data points (zero Tobs)
    
 */
std::pair<LibAIR2::ALMAAbsInpL,  std::vector<std::pair<double, double> > >
filterFlaggedInp(const LibAIR2::ALMAAbsInpL &inp,
		 const std::vector<std::pair<double, double> > &fb)
{

  LibAIR2::ALMAAbsInpL res;
  std::vector<std::pair<double, double> > rfb;
  size_t j=0;
  bool fbFilled = (fb.size()>0);
  for(LibAIR2::ALMAAbsInpL::const_iterator i=inp.begin();
      i!=inp.end();
      ++i, ++j)
  {
    if(i->TObs[0]>0.) // flagged ALMAAbsInp would have TObs==0
    {
      res.push_back(*i);
      if(fbFilled)
      {
	rfb.push_back(fb[j]);
      }
    }
  }
  return std::make_pair(res, rfb);
}

/** Return the set of antenna IDs that do not have a WVR
 */
LibAIR2::AntSet NoWVRAnts(const LibAIR2::aname_t &an)
{
  LibAIR2::AntSet res;
  for(LibAIR2::aname_t::const_iterator i=an.begin();
      i!= an.end();
      ++i)
  {
    if (i->second[0]=='C' and i->second[1]=='M')
      res.insert(i->first);
  }
  return res;
}

int main(int argc,  char* argv[])
{

  int rval = -2;

  option::Descriptor usage[] = {
      { wvr::arg::index::UNKNOWN,     0, "", "",            wvr::arg::check::Unknown,
        "\nAllowed options:\n" },
      { wvr::arg::index::HELP,        0, "", "help",        wvr::arg::check::None,
        " --help      \tPrint information about program usage" },
      { wvr::arg::index::MS,          0, "", "ms",          wvr::arg::check::Required,
        " --ms        \tInput measurement set" },
      { wvr::arg::index::OUTPUT,      0, "", "output",      wvr::arg::check::Required,
        " --output    \tName of the output file" },
      { wvr::arg::index::TOFFSET,     0, "", "toffset",     wvr::arg::check::Float,
        " --toffset   \tTime offset (in seconds) between interferometric and WVR data" },
      { wvr::arg::index::NSOL,        0, "", "nsol",        wvr::arg::check::Int,
        " --nsol      \tNumber of solutions for phase correction coefficients to make during this observation" },
      { wvr::arg::index::SEGFIELD,    0, "", "segfield",    wvr::arg::check::Deprecated,
        " --segfield  \tDo a new coefficient calculation for each field (this option is disabled, see segsource)" },
      { wvr::arg::index::SEGSOURCE,   0, "", "segsource",   wvr::arg::check::None,
        " --segsource \tDo a new coefficient calculation for each source" },
      { wvr::arg::index::REVERSE,     0, "", "reverse",     wvr::arg::check::None,
        " --reverse   \tReverse the sign of correction in all SPW (e.g. due to AIV-1740)" },
      { wvr::arg::index::REVERSESPW,  0, "", "reversespw",  wvr::arg::check::Int,
        " --reversespw \tReverse the sign correction for this spw" },
      { wvr::arg::index::DISPERSE,    0, "", "disperse" ,   wvr::arg::check::None,
        " --disperse   \tApply correction for dispersion" },
      { wvr::arg::index::WVRFLAG,     0, "", "wvrflag",     wvr::arg::check::Required,
        " --wvrflag    \tRegard this WVR (labelled with either antenna number or antenna name) as bad, and use interpolated values instead" },
      { wvr::arg::index::SOURCEFLAG,  0, "", "sourceflag",  wvr::arg::check::Required,
        " --sourceflag \tFlag the WVR data for this source and do not produce any phase corrections on it" },
      { wvr::arg::index::STATFIELD,   0, "", "statfield",   wvr::arg::check::Required,
        " --statfield  \tCompute the statistics (Phase RMS, Disc) on this field only" },
      { wvr::arg::index::STATSOURCE,  0, "", "statsource",  wvr::arg::check::Required,
        " --statsource  \tCompute the statistics (Phase RMS, Disc) on this source only" },
      { wvr::arg::index::TIE,         0, "", "tie",         wvr::arg::check::Required,
        " --tie  \tPrioritise tieing the phase of these sources as well as possible" },
      { wvr::arg::index::SMOOTH,      0, "", "smooth",      wvr::arg::check::Int,
        " --smooth  \tSmooth WVR data by this many samples before applying the correction" },
      { wvr::arg::index::SCALE,       0, "", "scale",       wvr::arg::check::Float,
        " --scale  \tScale the entire phase correction by this factor" },
      { wvr::arg::index::MAXDISTM,    0, "", "maxdistm",    wvr::arg::check::Float,
        " --maxdistm  \tmaximum distance (m) an antenna may have to be considered for being part of the <=3 antenna set for interpolation of a solution for a flagged antenna" },
      { wvr::arg::index::MINNUMANTS,  0, "", "minnumants",  wvr::arg::check::Int,
        " --minnumants  \tminimum number of near antennas (up to 3) required for interpolation" },
      { wvr::arg::index::MINGOODFRAC, 0, "", "mingoodfrac", wvr::arg::check::Float,
        " --mingoodfrac  \tIf the fraction of unflagged data for an antenna is below this value (0. to 1.), the antenna is flagged" },
      { wvr::arg::index::USEFIELDTAB, 0, "", "usefieldtab", wvr::arg::check::None,
        " --usefieldtab  \tDerive the antenna pointing information from the FIELD table instead of the POINTING table" },
      { wvr::arg::index::SPW,         0, "", "spw",         wvr::arg::check::Int,
        " --spw  \tOnly write out corrections for these SPWs" },
      { wvr::arg::index::WVRSPW,      0, "", "wvrspw",      wvr::arg::check::Int,
        " --wvrspw  \tOnly use data from these WVR SPWs" },
      { wvr::arg::index::REFANT,      0, "", "refant",      wvr::arg::check::Required,
        " --refant  \tUse the WVR data from this antenna for calculating the dT/dL parameters" },
      { wvr::arg::index::OFFSETS,     0, "", "offsets",     wvr::arg::check::Required,
        " --offsets  \tName of the optional input table containing the temperature offsets, e.g. generated by remove_cloud"},
      { 0, 0, 0, 0, 0, 0 }
  };

  // Defaults are set by parsing an argv-like set of options where the values are the defaults
  const char *defaults[] = { "--toffset=0",
                             "--nsol=1",
                             "--scale=1.0",
                             "--maxdistm=500.0",
                             "--minnumants=2",
                             "--mingoodfrac=0.8",
                             (const char *) -1 };   // unambiguously signal the end

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
  option::Option options[stats.options_max], buffer[stats.buffer_max];
  option::Parser parse;

  // parse the defaults first, then argv. User set options always come last
  // true here has same meaning as in stats above. This may not be necessary here, I think
  // the stats usage above has already reorderded argv in place.
  parse.parse(true, usage, defaultCount, defaults, options, buffer);
  parse.parse(true, usage, argc, argv, options, buffer);

    LibAIR2::printBanner(std::cout);

    if (options[wvr::arg::index::HELP] || (argc==0) ) {
        std::cout<<"Write out a gain table based on WVR data"
                 <<std::endl
                 <<std::endl
                 <<"GPL license -- you have the right to the source code. See COPYING"
                 <<std::endl
                 <<std::endl;
        option::printUsage(std::cout,usage);
        return 0;
    }

    if (checkPars(options)) return -1;

    checkWarnPars(options);

    std::vector<std::set<std::string> > tied=getTied(options);

    std::string msname(options[wvr::arg::index::MS].last( )->arg);
    casacore::MeasurementSet ms(msname);

    checkMSandPars(ms, options);

    std::vector<int> wvrspws;
    {
        LibAIR2::SPWSet thewvrspws=LibAIR2::WVRSPWIDs(ms);
        fprintf( stderr, "WVRSPW: %d\n", options[wvr::arg::index::WVRSPW].count( ) );
        if (options[wvr::arg::index::WVRSPW].count( )) {
            option::Option *args=options[wvr::arg::index::WVRSPW];
            for(auto arg = args; arg; arg = arg->next( )) {
                int spw = std::stod(arg->arg);
                if(thewvrspws.count(spw)==0) {
                    std::cout << "ERROR: SPW " << spw << " is not a WVR SPW or invalid." <<std::endl;
                    std::cerr << "ERROR: SPW " << spw << " is not a WVR SPW or invalid." <<std::endl;
                    return -10;
                }
                wvrspws.push_back(spw);
            }

            std::cout<<"Will use the following WVR SPWs:"<<std::endl;
            for(size_t i=0; i<wvrspws.size();i++){
                std::cout<< " " << wvrspws[i];
            }
            std::cout <<std::endl;
        }

        if (wvrspws.size()==0) {
            std::cout<<"Will use all WVR SPWs:"<<std::endl;
            for(LibAIR2::SPWSet::const_iterator si=thewvrspws.begin(); si!=thewvrspws.end();++si) {
                wvrspws.push_back(*si);
                std::cout<< " " <<  *si;
            }
            std::cout <<std::endl;
        }
    }

    std::vector<int> sciencespws;
    LibAIR2::SPWSet thewvrspws=LibAIR2::WVRSPWIDs(ms);

    if (options[wvr::arg::index::SPW].count( )) {
        option::Option *args=options[wvr::arg::index::SPW];
        for(auto arg = args; arg; arg = arg->next( )) {
            int spw = std::stod(arg->arg);
            if(thewvrspws.count(spw) != 0) {
                std::cout<<"WARNING: SPW "<< spw << " is a WVR SPW, not a science SPW." <<std::endl;
                std::cerr<<"WARNING: SPW "<< spw << " is a WVR SPW, not a science SPW." <<std::endl;
            }
            int nspw=LibAIR2::numSPWs(ms);
            if(spw < 0 || spw >= nspw) {
                std::cout<<"ERROR: Invalid SPW "<< spw <<std::endl;
                std::cerr<<"ERROR: Invalid SPW "<< spw <<std::endl;
                return -11;
            }
            sciencespws.push_back(spw);
        }
            
        if (sciencespws.size() > 0){
            std::cout<<"Will produce solutions for the following SPWs:"<<std::endl;
            for(size_t i=0; i<sciencespws.size();i++){
                std::cout<< " " << sciencespws[i];
            }
            std::cout <<std::endl;
        }
    }

    if(sciencespws.size()==0){
        std::cout<<"Will produce solutions for all SPWs:"<<std::endl;
        for(size_t i=0; i<LibAIR2::numSPWs(ms);i++){
            sciencespws.push_back(i);
            std::cout<< " " << i;
        }
        std::cout <<std::endl;
    }

    std::string fnameout=options[wvr::arg::index::OUTPUT].last( )->arg;

    std::string offsetstable="";
    if(options[wvr::arg::index::OFFSETS].count( )){
        offsetstable = options[wvr::arg::index::OFFSETS].last( )->arg;
    }

    std::set<size_t> useID=LibAIR2::skyStateIDs(ms);

    LibAIR2::AntSet wvrflag;
    // Prepare flagging and interpolation
    if (options[wvr::arg::index::WVRFLAG].count( ))
        {
            wvrflag=getAntPars("wvrflag", options[wvr::arg::index::WVRFLAG], ms);
        }

    LibAIR2::aname_t anames=LibAIR2::getAName(ms);
    LibAIR2::AntSet nowvr=NoWVRAnts(anames);
  
    LibAIR2::AntSet interpwvrs(wvrflag); // the antennas to interpolate solutions for
    interpwvrs.insert(nowvr.begin(), nowvr.end());

    LibAIR2::AntSet flaggedants; // the antennas flagged in the ANTENNA table are not to be interpolated
    LibAIR2::WVRAddFlaggedAnts(ms, flaggedants);

    wvrflag.insert(flaggedants.begin(),flaggedants.end());

    if(interpwvrs.size()+flaggedants.size()==ms.antenna().nrow()){
        std::cout << "No good antennas with WVR data found." << std::endl;
        std::cerr << "No good antennas with WVR data found." << std::endl;
        return -1;
    }

    int iterations = 0;

    while(rval<0 && iterations<2){

        iterations++;

        std::vector<size_t> sortedI; // to be filled with the time-sorted row number index
        std::set<int> flaggedantsInMain; // the antennas totally flagged in the MS main table
        std::shared_ptr<LibAIR2::InterpArrayData>
            d (LibAIR2::loadWVRData( ms, wvrspws, sortedI, flaggedantsInMain,
                                     std::stod(options[wvr::arg::index::MINGOODFRAC].last( )->arg),
                                     options[wvr::arg::index::USEFIELDTAB].count( )==0,
                                     offsetstable)
               );

     // For debug purposes, print the loaded WVR data: 
     // for(size_t j=0; j<d->g_time().size(); ++j)
     // {
     // 	for(size_t i=0; i < ms.antenna().nrow(); ++i)
     // 	  {
     // 	    std::cout << "row ant data " << j << " " << i << " "
     // 		      << d->g_wvrdata()[j][i][0] << " "
     // 		      << d->g_wvrdata()[j][i][1] << " " 
     // 		      << d->g_wvrdata()[j][i][2] << " "
     // 		      << d->g_wvrdata()[j][i][3] << std::endl; 
     // 	  }
     // }
     

     interpwvrs.insert(flaggedantsInMain.begin(),flaggedantsInMain.end()); // for flagInterp()
     wvrflag.insert(flaggedantsInMain.begin(),flaggedantsInMain.end());

     d->offsetTime(std::stod(options[wvr::arg::index::TOFFSET].last( )->arg));
     
     if (options[wvr::arg::index::SMOOTH].count( ))
         smoothWVR(*d, std::stoi(options[wvr::arg::index::SMOOTH].last( )->arg));
     
     d.reset(LibAIR2::filterState(*d, useID));

     LibAIR2::AntSet interpImpossibleAnts;

     // Flag and interpolate
     flagInterp( ms, interpwvrs, *d,
                 std::stod(options[wvr::arg::index::MAXDISTM].last( )->arg),
                 std::stoi(options[wvr::arg::index::MINNUMANTS].last( )->arg),
                 interpImpossibleAnts );

     // Determine the reference antenna for dTdL calculation
     int refant = -1; 

     if (options[wvr::arg::index::REFANT].count( )) {
         LibAIR2::AntSet refants=getAntPars("refant", options[wvr::arg::index::REFANT], ms);    
         for(LibAIR2::AntSet::iterator it=refants.begin(); it != refants.end(); it++){
             if(interpImpossibleAnts.count(*it)==0){
                 refant = *it; // use the first of the given list of possible ref antennas which was OK or which could be interpolated to
                 break;
             }
             else{
                 std::cout << "Given reference antenna " << *it << "==" << anames.at(*it) 
                           << " is flagged and cannot be interpolated." << std::endl;
             }	   
         }
         if(refant<0){
             std::cout << "None of the given reference antennas is usable." << std::endl;
             std::cerr << "None of the given reference antennas is usable." << std::endl;
             return -1;
         }
     } else {
         LibAIR2::AntSet wvrants=LibAIR2::WVRAntennas(ms, wvrspws);
         for(LibAIR2::AntSet::iterator it=wvrants.begin(); it != wvrants.end(); it++){
             if(interpImpossibleAnts.count(*it)==0) {
                 refant = *it; // use the first antenna which was OK or which could be interpolated to
                 break;
             }
         }
         if(refant<0) {
             std::cout << "No antennas with sufficient WVR data found." << std::endl;
             std::cerr << "No antennas with sufficient WVR data found." << std::endl;
             return -1;
         }
     }

     std::cout << "Choosing";
     if(interpwvrs.count(refant)>0){
         std::cout << " (interpolated)";
     }
     std::cout << " antenna " << refant  << " == " << anames.at(refant)
               << " as reference antenna for dTdL calculations." << std::endl;


     LibAIR2::ArrayGains g(d->g_time(), 
			   d->g_el(),
			   d->g_state(),
			   d->g_field(),
			   d->g_source(),
			   d->nAnts);
     
     std::shared_ptr<LibAIR2::dTdLCoeffsBase>  coeffs;
     
     // These are the segments on which coefficients are re-calculated
     std::vector<std::pair<double, double> >  fb;
     
     if ( options[wvr::arg::index::CONT].count( ) )
     {
         std::cout<<"[Output from \"cont\" option has not yet been updated]"
                  <<std::endl;
         coeffs.reset(LibAIR2::SimpleSingleCont(*d, refant));
     } else {
         LibAIR2::ALMAAbsInpL inp;
         if (options[wvr::arg::index::SEGSOURCE].count( ))	{
             std::vector<int> flds;
             std::vector<double> time;
             std::vector<int> src;
             LibAIR2::fieldIDs( ms, time, flds, src, sortedI );
             try {
                 std::vector<std::set<size_t> >  tiedi=tiedIDs(tied, ms);
	     
                 printTied(tied, tiedi);
                 LibAIR2::fieldSegmentsTied( time, src, tiedi, fb );
             } catch(LibAIR2::WVRUserError& x){
                 std::cout << x.what() << std::endl;
                 std::cerr << x.what() << std::endl;
                 return -1;
             }

	   //printFieldSegments(fb, time[0]);
	   
//       { // debugging output
// 	std::vector<double> tt(d->g_time());
// 	std::vector<double> te(d->g_el());
// 	std::vector<size_t> ts(d->g_state());
// 	std::vector<size_t> tf(d->g_field());
// 	std::vector<size_t> tsou(d->g_source());
// 	for(uint i=0; i<tt.size(); i++){
// 	  std::cerr << "i time el state field source " << i << " " << tt[i] << " ";
// 	  std::cerr << te[i] << " ";
// 	  std::cerr << ts[i] << " ";
// 	  std::cerr << tf[i] << " ";
// 	  std::cerr << tsou[i] << std::endl;
// 	}
// 	std::cerr << "nAnts " << d->nAnts << std::endl;
// 	for(uint i=0; i<time.size(); i++){
// 	  std::cerr << "i time  " << i << " " << time[i] << std::endl;
// 	}
// 	for(uint i=0; i<fb.size(); i++){
// 	  std::cerr << "i fb  " << i << " " << fb[i].first << " " << fb[i].second  << std::endl;
// 	}
// 	for(std::set<size_t>::iterator it = useID.begin(); it != useID.end(); it++){
// 	  std::cerr << "useID " << *it << std::endl;
// 	}
//       }

             inp=FieldMidPointI( *d, fb, useID, refant );

         } else	{
             const size_t n=std::stoi(options[wvr::arg::index::NSOL].last( )->arg);
             inp=LibAIR2::MultipleUniformI( *d, n, useID, refant );
         }

	
         if (options[wvr::arg::index::SOURCEFLAG].count( )) {
             std::vector<std::string> flags;
             option::Option *args = options[wvr::arg::index::SOURCEFLAG];
             for (auto arg = args; arg; arg = arg->next( )) flags.push_back(arg->arg);
             std::tie(inp,fb) = filterInp( inp, fb, flags, ms);
         }
	
         std::tie(inp,fb) = filterFlaggedInp( inp, fb );

         std::cerr << "Calculating the coefficients now ... " << std::endl;
         std::list<std::shared_ptr<LibAIR2::ALMAResBase> > rlist;
         LibAIR2::AntSet problemAnts;

         rval = 0;

         try {
             fprintf(stderr, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<LibAIR2::doALMAAbsRet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
             rlist=LibAIR2::doALMAAbsRet( inp, fb, problemAnts );
             fprintf(stderr, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<LibAIR2::doALMAAbsRet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
         } catch(const std::runtime_error rE) {
             rval = 1;
             std::cerr << std::endl << "WARNING: problem while calculating coefficients:"
                       << std::endl << "         LibAIR2::doALMAAbsRet: " << rE.what() << std::endl;
             std::cout << std::endl << "WARNING: problem while calculating coefficients:"
                       << std::endl << "         LibAIR2::doALMAAbsRet: " << rE.what() << std::endl;
         }
	
         if( problemAnts.size()>0 ) {
	   
             rval = -2;

             if(iterations<2){
                 for(LibAIR2::AntSet::const_iterator it=problemAnts.begin(); it!=problemAnts.end(); it++){
                     if(interpwvrs.count(*it)==0){
                         std::cerr	<< "Flagging antenna " << *it << " == " << anames.at(*it) << std::endl;
                         std::cout	<< "Flagging antenna " << *it << " == " << anames.at(*it) << std::endl;
                         interpwvrs.insert(*it); // for flagInterp()
                         wvrflag.insert(*it); // for later log output
                     }
                 }
                 std::cerr	<< "Reiterating ..." << std::endl;
                 std::cout	<< "Reiterating ..." << std::endl;
                 continue;
             } else {
                 std::cerr << "Number of remaining antennas with problematic WVR measurements: " << problemAnts.size() << std::endl;
                 std::cout << "Number of remaining antennas with problematic WVR measurements: " << problemAnts.size() << std::endl;
                 std::cerr << "Will continue without further iterations ..." << std::endl;
                 std::cout << "Will continue without further iterations ..." << std::endl;
             }	      
         }	   
	
         std::cerr<<"done!"
                  <<std::endl;
	
	
         std::cout<<"       Retrieved parameters      "<<std::endl
                  <<"----------------------------------------------------------------"<<std::endl
                  <<rlist<<std::endl;
	
         if ( options[wvr::arg::index::SEGSOURCE].count( ) ) {
             //std::vector<int> flds;
             //std::vector<double> time;
             //std::vector<int> src;
             //LibAIR2::fieldIDs(ms, 
             //		    time,
             //		    flds,
             //		    src,
             //		    sortedI);
	   
             coeffs.reset(LibAIR2::SimpleMultiple(fb,rlist));   
         } else {
             coeffs.reset(LibAIR2::ALMAAbsProcessor(inp, rlist));
         }  
	
     }
    
     try{
         g.calc(*d,*coeffs);    
     } catch(const std::runtime_error& x) {
         std::cout << "Problem while calculating gains: " << x.what() << std::endl;
         std::cerr << "Problem while calculating gains: " << x.what() << std::endl;
         return 1;
     }

     if (options[wvr::arg::index::SOURCEFLAG].count( )) {
         std::vector<std::string> flags;
         option::Option *args = options[wvr::arg::index::SOURCEFLAG];
         for (auto arg = args; arg; arg = arg->next( )) flags.push_back(arg->arg);
         std::set<size_t> flagset=sourceSet( flags, ms );
         g.blankSources(flagset);
     }
     
     std::vector<std::pair<double, double> > tmask;
     statTimeMask(ms, options, tmask, sortedI, wvrspws);
     
     std::vector<double> pathRMS;
     g.pathRMSAnt(tmask, pathRMS);
     
     
     std::vector<double> pathDisc;
     try{
       computePathDisc(*d, tmask, *coeffs, pathDisc);
     
       std::cout<<LibAIR2::AntITable(anames, wvrflag, nowvr, pathRMS, pathDisc, interpImpossibleAnts);
       
       printExpectedPerf(g, *coeffs, tmask);
       
     } catch(const std::runtime_error& x) {
         std::cout << "Problem while calculating path RMS discrepancy: " << x.what() << std::endl;
         std::cerr << "Problem while calculating path RMS discrepancy: " << x.what() << std::endl;
         return 1;
     }
     
     if (options[wvr::arg::index::SCALE].count( )) {
         g.scale(std::stod(options[wvr::arg::index::SCALE].last( )->arg));
     }
     
     LibAIR2::MSSpec sp;
     loadSpec(ms, sciencespws, sp);
     std::set<size_t> reverse=reversedSPWs(sp, options);
     
     std::cout << "Writing gain table ..." << std::endl;

     // Write new table, including history
     LibAIR2::writeNewGainTbl( g, fnameout.c_str(), sp, reverse,
                               options[wvr::arg::index::DISPERSE].count( )>0,
                               msname, buildCmdLine(argc,argv), interpImpossibleAnts);

#ifdef BUILD_HD5
     LibAIR2::writeAntPath(g,fnameout+".hd5");
#endif

  } // end while


  return rval;
}
