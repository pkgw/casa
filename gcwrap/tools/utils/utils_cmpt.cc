
/***
 * Framework independent implementation file for utils...
 *
 * Implement the utils component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <fstream>
#include <stdcasa/record.h>
#include <stdcasa/version.h>
#include <utils_cmpt.h>
#include <tools/utils/stdBaseInterface.h>
#include <tools/xerces/stdcasaXMLUtil.h>
#include <casa/Logging/LogIO.h>
#include <casa/BasicSL/String.h>
#include <casa/OS/File.h>
#include <casa/OS/DOos.h>
#include <tables/Tables/Table.h>
#include <casa/System/Aipsrc.h>
#include <casa/OS/HostInfo.h>
#ifndef NO_CRASH_REPORTER
#include <stdcasa/StdCasa/CrashReporter.h>
#endif
#include <signal.h>
#include <string>
#include <vector>
#include <cstdlib>


using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

utils::utils()
{
  myConstraints = 0;
  itsLog = new casacore::LogIO;
}

utils::~utils()
{
  if(myConstraints)
     delete myConstraints;
  delete itsLog;
}

bool
utils::verify(const ::casac::record& input, const ::casac::variant& xmldescriptor, bool throwexcept)
{

   bool rstat(true);
   record *constraints(0);
   *itsLog << LogOrigin("utils", "verify") << LogIO::NORMAL3 << "Verifying arguments....";
   switch(xmldescriptor.type()){
      case variant::STRING :
	 constraints = torecord(xmldescriptor.getString());
	  //std::cerr << "constraints record: ";
	  //dumpRecord(std::cerr, *constraints);
         break;
      case variant::RECORD :
         constraints = new record(xmldescriptor.getRecord());
         break;
      default :
         rstat = false;
         break;
   }
   if(rstat){
	   rstat = stdBaseInterface::verify(const_cast<record &>(input), *constraints, *itsLog);
           if(constraints)
	      delete constraints;
	   if(rstat){
		   *itsLog << LogOrigin("utils", "verify") << LogIO::NORMAL3 << "verified." << LogIO::POST;
	   }else{
                if(throwexcept){
                   throw(AipsError("Parameter verification failed"));
                } else {
		   *itsLog <<  LogIO::POST;
		   *itsLog << LogOrigin("utils", "verify") << LogIO::WARN << "Some arguments failed to verify!" << LogIO::POST;
		}
	   }
   }
   //std::cerr << "return from verify is " << rstat << std::endl;
   return rstat;
}

bool
utils::setconstraints(const ::casac::variant& xmldescriptor)
{
   bool rstat(true);
   if(myConstraints)
	   delete myConstraints;
   *itsLog << LogOrigin("utils", "setconstraints") << LogIO::NORMAL3 << "Setting constraints ...";
   switch(xmldescriptor.type()){
      case variant::STRING :
	 myConstraints = torecord(xmldescriptor.getString());
	 // std::cerr << "constraints record: ";
	 // dumpRecord(std::cerr, *constraints);
         break;
      case variant::RECORD :
         myConstraints = new record(xmldescriptor.getRecord());
         break;
      default :
	 rstat = false;
   }
   *itsLog << LogIO::NORMAL3 << "Constraints set." << LogIO::POST;
	 //std::cerr << "constraints record: ";
	 //dumpRecord(std::cerr, *myConstraints);
   return rstat;
}

bool
utils::verifyparam(const ::casac::record& param)
{
   bool rstat(true);
   if(myConstraints && !param.empty()){
      //dumpRecord(std::cerr, param);
      rec_map::iterator iter = myConstraints->begin(); // We need the underlying record...
      /*
      cerr << "Constraints Record " << endl;
      dumpRecord(std::cerr, *myConstraints);
      cerr << "Param  Record " << endl;
      dumpRecord(std::cerr, param);
      */
      rstat = stdBaseInterface::verifyOne(const_cast<record &>(param), (*iter).second.asRecord(), *itsLog);
   } else {
	 if(param.empty()){
            *itsLog << LogOrigin("utils", "verifyparam") << LogIO::WARN
		 << "parameter not set, unable to verify parameter" << LogIO::POST;
	 }else{
            *itsLog << LogOrigin("utils", "verifyparam") << LogIO::WARN
		 << "Constraints record not set, unable to verify parameter" << LogIO::POST;
	 }
   }
   return rstat;
}

::casac::variant*
utils::expandparam(const std::string& name , const ::casac::variant& value )
{
   ::casac::variant *rstat(0);

   if(myConstraints){
       rec_map::iterator iter = myConstraints->begin(); // We need the underlying record...
       //dumpRecord(std::cerr, (*iter).second.asRecord()["parameters"].asRecord()[name].asRecord());
       if((*iter).second.asRecord()["parameters"].asRecord().count(name) &&
          (*iter).second.asRecord()["parameters"].asRecord()[name].asRecord().count("allowed")){
          rstat = stdBaseInterface::expandEnum((*iter).second.asRecord()["parameters"].asRecord()[name].asRecord()["allowed"], value, *itsLog);
       }
       else{
	  rstat = new variant(value);
       }
   } else {
       rstat = new variant(casac::initialize_variant(""));;
       *itsLog << LogOrigin("utils", "expandparam") << LogIO::WARN
		 << "Constraints record not set, unable to expand parameter" << LogIO::POST;
   }
   return rstat;
}

::casac::record*
utils::torecord(const std::string& input)
{
   stdcasaXMLUtil xmlUtils;
   casac::record *rstat = new casac::record;
   if(!input.find("<?xml version")){
      xmlUtils.toCasaRecord(*rstat, input);
   }else{
      if(!input.find("file:///")){
         Bool ok = xmlUtils.readXMLFile(*rstat, input.substr(7));
	 if(!ok){
            *itsLog << LogIO::SEVERE << "Unable to read XML file " << input << ", unable to verify input" << LogIO::POST;
	 }
      } else {
         *itsLog << LogIO::SEVERE << "Defaults specified are not an XML string, unable to verify input" << LogIO::POST;
      }
   }
   return rstat;
}

std::string
utils::toxml(const ::casac::record& input, const bool asfile, const std::string& filename)
{   string rstat;

   stdcasaXMLUtil xmlUtils;
   if(asfile){
      std::ofstream xmlout(filename.c_str(), ios::out);
      xmlUtils.fromCasaRecord(xmlout, input);
      rstat = filename;
   } else {
      ostringstream xmlout;
      xmlUtils.fromCasaRecord(xmlout, input);
      rstat = xmlout.str();
   }
   return rstat;
}

std::string
utils::getrc(const std::string& rcvar)
{
  String rstat1;
  if(!rcvar.length()){
	  rstat1 = Aipsrc::aipsRoot();
  } else {
	  if(!Aipsrc::find(rstat1, rcvar))
		  rstat1 = "Unknown value";
  }
  string rstat(rstat1.c_str());
  return rstat;
}

bool
utils::removetable(const std::vector<std::string> &tablenames)
{
  bool rstat(true);
  try {
     *itsLog << LogOrigin("utils", "removetable");
     for(vector<std::string>::const_iterator iter = tablenames.begin();
		     iter != tablenames.end(); iter++){
       String fileName(*iter);
       if (fileName.empty()) {
          *itsLog << LogIO::WARN << "Empty filename" << LogIO::POST;
          rstat = false;
       }
       File f(fileName);
       if (! f.exists()) {
           *itsLog << LogIO::WARN << fileName << " does not exist." << LogIO::POST;
          rstat = false;
       }

// Now try and blow it away.  If it's open, tabledelete won't delete it.
       String message;
       if(rstat && Table::isReadable(fileName)){
          if (Table::canDeleteTable(message, fileName, true)) {
             Table::deleteTable(fileName, true);
          } else {
             *itsLog << LogIO::WARN << "Cannot delete file " << fileName
             << " because " << message << LogIO::POST;
          }
       } else {
           *itsLog << LogIO::WARN << "Cannot delete file " << fileName
           << " because it's not a table." << LogIO::POST;
       }
  }
    } catch (AipsError x) {
       *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
       << LogIO::POST;
       RETHROW(x);
  }
  return rstat;
}

::casac::record *utils::tableinfo(const std::string &tablename) {
    Vector<Int> info = casacore::DOos::lockInfo(tablename);
    ::casac::record *result = new record( );

    switch( info[0] ) {
    case 3:
	result->insert("lockstatus", "write");
	break;
    case 2:
	result->insert("lockstatus", "read");
	break;
    case 1:
	result->insert("lockstatus", "open");
	break;
    case 0:
	result->insert("lockstatus", "not in use");
	break;
    default:
	result->insert("lockstatus", "unknown");
    }

    result->insert("lockpid", info[1]);
    result->insert("lockperm", info[2] ? true : false);
    return result;
}

std::vector<std::string> utils::lockedtables( ) {
    Vector<String> locks = Table::getLockedTables( );
    std::vector<std::string> result;
    for (unsigned int x = 0; x < locks.nelements(); ++x ) {
	result.push_back(locks[x]);
    }
    return result;
}


typedef int SIZETCAST;
::casac::record *utils::hostinfo( ) {
    ::casac::record *result = new record( );

    ::casac::record *swap = new record( );
    swap->insert( "total", (SIZETCAST) HostInfo::swapTotal( ) );
    swap->insert( "used", (SIZETCAST) HostInfo::swapUsed( ) );
    swap->insert( "free", (SIZETCAST) HostInfo::swapFree( ) );
    result->insert( "swap", swap );

    ::casac::record *memory = new record( );
    memory->insert( "total", (SIZETCAST) HostInfo::memoryTotal( ) );
    memory->insert( "available", (SIZETCAST) HostInfo::memoryTotal(true) );
    memory->insert( "used", (SIZETCAST) HostInfo::memoryUsed( ) );
    memory->insert( "free", (SIZETCAST) HostInfo::memoryFree( ) );
    result->insert( "memory", memory );

    ::casac::record *cpus = new record( );
    cpus->insert( "total", HostInfo::numCPUs( ) );
    cpus->insert( "available", HostInfo::numCPUs(true) );
    result->insert( "cpus", cpus );

    result->insert( "endian", HostInfo::bigEndian( ) ? "big" : "little" );
    result->insert( "hostname", HostInfo::hostName( ) );
    result->insert( "pid", HostInfo::processID( ) );

    result->insert( "seconds", HostInfo::secondsFrom1970( ) );
    
    return result;
}

std::string
utils::c_exception ()
{
  String lastMessage, lastStackTrace;
  AipsError::getLastInfo (lastMessage, lastStackTrace);

  String result = lastMessage + "\n" + lastStackTrace;

  return result;
}

void
utils::c_exception_clear ()
{
  AipsError::clearLastInfo ();
}

void bogusHandler (int, siginfo_t *, void *)
{
    // Do nothing
}

string
utils::_crash_reporter_initialize (const string & crashDirectory,
                                   const string & crashPosterApplication,
                                   const string & crashPostingUrl,
				   const string & logFile)
{
#ifndef NO_CRASH_REPORTER
    // *NOTE*: Not intended for casual use!

    string status = casa::CrashReporter::initialize(crashDirectory, crashPosterApplication,
                                                    crashPostingUrl, logFile);

    return status;
#else
    return "no-op";
#endif
}

bool
utils::_trigger_segfault (int faultType)
{
    // *NOTE*: Not intended for casual use!

    switch (faultType) {

    case 0:{
	bool * p;
        long zero = 0;
	p = (bool *) zero;
	return * p;
	break;
    }

    default:
    case 1:{
	throw exception();
	break;
    }

    }

    return false;
}


std::vector<int>
utils::version( ) {
    std::vector<int> result = {
        VersionInfo::major( ),
        VersionInfo::minor( ),
        VersionInfo::patch( ),
        VersionInfo::feature( )
    };
    return result;
}

std::string
utils::version_desc( ) { return VersionInfo::desc( ); }

std::string
utils::version_info( ) { return VersionInfo::info( ); }

std::string
utils::version_string( ) { return VersionInfo::str( ); }

bool
 utils::compare_version(const  string& comparitor,  const std::vector<int>& vec) {
  vector<int> current_version = version( );
  for ( unsigned int i=0; i < vec.size( ); ++i )
    if ( vec[i] < 0 ) throw(AipsError("negative values not allowed in version numbers"));

  unsigned int limit = min(current_version.size( ),vec.size( ));
  if ( comparitor == ">" ) {
    for ( unsigned int i=0; i < limit; ++i ) {
      if ( current_version[i] > vec[i] ) return true;
      else if ( current_version[i] < vec[i] ) return false;
    }
    for ( unsigned int i=limit; i < current_version.size( ); ++i )
      if ( current_version[i] > 0 ) return true;
    return false;
  } else if ( comparitor == "<" ) {
    for ( unsigned int i=0; i < limit; ++i ) {
      if ( current_version[i] > vec[i] ) return false;
      else if ( current_version[i] < vec[i] ) return true;
    }
    return false;
  } else if ( comparitor == ">=" ) {
    for ( unsigned int i=0; i < limit; ++i ) {
      if ( current_version[i] > vec[i] ) return true;
      else if ( current_version[i] < vec[i] ) return false;
    }
    return true;
  } else if ( comparitor == "<=" ) {
    for ( unsigned int i=0; i < limit; ++i ) {
      if ( current_version[i] > vec[i] ) return false;
      else if ( current_version[i] < vec[i] ) return true;
    }
    for ( unsigned int i=limit; i < current_version.size( ); ++i )
      if ( current_version[i] > 0 ) return false;
    return true;
  } else if ( comparitor == "=" || comparitor == "==" ) {
    for ( unsigned int i=0; i < limit; ++i ) {
      if ( current_version[i] != vec[i] ) return false;
    }
    for ( unsigned int i=limit; i < current_version.size( ); ++i )
      if ( current_version[i] > 0 ) return false;
    return true;
  } else if ( comparitor == "!=" ) {
    return ! compare_version("=",vec);
  } else {
    throw(AipsError("unknown comparator"));
  }
  return false;
}

} // casac namespace

