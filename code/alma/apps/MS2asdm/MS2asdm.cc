//# MS2asdm.cc 
//#
//#  ALMA - Atacama Large Millimeter Array
//#  (c) European Southern Observatory, 2002
//#  (c) Associated Universities Inc., 2002
//#  Copyright by ESO (in the framework of the ALMA collaboration),
//#  Copyright by AUI (in the framework of the ALMA collaboration),
//#  All rights reserved.
//#  
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#  
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#  
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <string>

#include <stdcasa/optionparser.h>
#include <alma/Options/AlmaArg.h>
using namespace alma;

#include <alma/ASDM/ASDMAll.h>

#include <alma/MS2ASDM/MS2ASDM.h>

#include <exception>
#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>

string appName;
bool verbose = true;
using namespace casacore;
using namespace casa;
using namespace std;

void info(const string& message) {  
  if(!verbose){
    return;
  }
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void error(const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::SEVERE));
  exit(1);
}

// A facility to get rid of blanks at start and end of a string.
// 
string lrtrim(std::string& s,const std::string& drop = " ")
{
  std::string r=s.erase(s.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}


#include <iostream>
#include <sstream>

/**
 * The main function.
 */
int main(int argc, char *argv[]) {

  ostringstream errstream;

  String asdmfile = "";
  String msfile = "";
  String datacolumn = "DATA";
  String archiveid = "S0";
  String rangeid = "X1"; 
  bool verbose = false;
  bool showversion = false;
  // these two defaults must also appear in the defaults array below
  double subscanduration = 24.*3600.; // default is one day
  double schedblockduration = 2700.; // default is 45 minutes
  bool apcorrected = true;

  appName = string(argv[0]);
  ofstream ofs;

  //   Process command line options and parameters.

  enum optionIndex { UNKNOWN, HELP, DATACOLUMN, ARCHIVEID, RANGEID, SUBSCANDUR, 
		     SCHEDBLOCKDUR, LOGFILE, APUNCORR, VERBOSE, REVISION, MSDIR, ASDMDIR };
 
  // remove the program name
  argc--;
  argv++;

  string usageIntro = 
    "Converts an ASDM dataset into a CASA measurement set.\n"
    "Usage : " + appName +" [options] ms-directory asdm-directory\n\n"
    "Command parameters: \n";

  // Descriptor elements are: OptionIndex, OptionType, shortopt, longopt, check_arg, help
  option::Descriptor usage[] = {
    { UNKNOWN, 0, "", "", AlmaArg::Unknown, usageIntro.c_str()},
    { UNKNOWN, 0, "", "", AlmaArg::Unknown, " \tms-directory : \tms directory"},
    { UNKNOWN, 0, "", "", AlmaArg::Unknown, " \tasdm-directory : \tasdm directory"},
    { UNKNOWN, 0, "", "", AlmaArg::Unknown, "\nAllowed options:\n"},
    { UNKNOWN, 0, "", "", AlmaArg::Unknown, 0}, // helps with formatting
 
    // these are the non-positional options
    { HELP, 0, "", "help", AlmaArg::None, " --help  \tproduces this help message."},
    { DATACOLUMN, 0, "d", "datacolumn", AlmaArg::Required, " -d [--datacolumn] arg (=DATA) \tspecifies the datacolumn."},
    { ARCHIVEID, 0, "a", "archiveid", AlmaArg::Required, " -a [--archiveid] arg (=S0) \tspecifies the archive ID."},
    { RANGEID, 0, "g", "rangeid", AlmaArg::Required, " -g [--rangeid] arg (=X1) \tspecifies the range ID."},
    { SUBSCANDUR, 0, "s", "subscanduration", AlmaArg::Float, " -s [--subscanduration] arg (=86400)  \tspecifies the maximum duration of a subscan in the output ASDM (seconds). Default: 86400"},
    { SCHEDBLOCKDUR, 0, "", "schedblockduration", AlmaArg::Float, " --schedblockduration  arg (=2700) \tspecifies the maximum duration of a scheduling block in the output ASDM (seconds). Default: 2700"},
    { LOGFILE, 0, "l", "logfile", AlmaArg::Required, " -l [--logfile] arg \tspecifies the log filename. If the option is not used then the logged informations are written to the standard error stream."},
    { APUNCORR, 0, "u", "apuncorrected", AlmaArg::None, " -u [--apuncorrected] \tthe data given by datacolumn should be regarded as not having an atmospheric phase correction. Default: data is AP corrected."},
    { VERBOSE, 0, "v", "verbose", AlmaArg::None, " -v [--verbose]  \tlogs numerous informations as the filler is working."},
    { REVISION, 0, "r", "revision", AlmaArg::None, " -r [--revision]  \tlogs information about the revision of this application."},
    
    // hidden options, allowed on the command line but do not show up in help.
    // the positional options take precendence 
    { MSDIR, 0, "", "ms-directory", AlmaArg::Required, 0},
    { ASDMDIR, 0, "", "asdm-directory", AlmaArg::Required, 0},
    { 0, 0, 0, 0, 0, 0} };
  
  // Defaults are set by parsing an argv-like set of options where the values are the defaults
  const char *defaults[] = { "--subscanduration=86400",
			     "--schedblockduration=2700",
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
  
  // Where do the log messages go ?
  if (options[LOGFILE] != NULL && options[LOGFILE].last()->arg != NULL) {
    ofs.open(options[LOGFILE].last()->arg, ios_base::app);
    LogSinkInterface *theSink = new casacore::StreamLogSink(&ofs);
    LogSink::globalSink(theSink);
  }

  if (options[HELP] || (argc==0)) {
    errstream.str("");
    option::printUsage(errstream, usage, 80);
    error(errstream.str());
  }
  
  // AP corrected will be true unless APUNCORR has been set
  apcorrected = (options[APUNCORR]==NULL);
  // Verbose or quiet ?
  verbose = (options[VERBOSE]!=NULL);
  // showversion ?
  showversion = (options[REVISION]!=NULL);

  if (options[DATACOLUMN] != NULL && options[DATACOLUMN].last()->arg != NULL) {
    datacolumn = String(options[DATACOLUMN].last()->arg);
  }

  if (options[ARCHIVEID] != NULL && options[ARCHIVEID].last()->arg != NULL) {
    archiveid = String(options[ARCHIVEID].last()->arg);
  }

  if (options[RANGEID] != NULL && options[RANGEID].last()->arg != NULL) {
    msfile = String(options[RANGEID].last()->arg);
  }

  // this always has a value because of defaults
  {
    stringstream str(options[SUBSCANDUR].last()->arg);
    str >> subscanduration;
    if (!str) {
      // unlikely given that any value was already checked by AlmaArg::Float
      errstream.str("");
      errstream << "There was an error converting the subscanduration value to a double";
      error(errstream.str());
    }
  }

  // this always has a value because of defaults
  {
    stringstream str(options[SCHEDBLOCKDUR].last()->arg);
    str >> schedblockduration;
    if (!str) {
      // unlikely given that any value was already checked by AlmaArg::Float
      errstream.str("");
      errstream << "There was an error converting the schedblockduration value to a double";
      error(errstream.str());
    }
  }

  if (parse.nonOptionsCount() > 0 || options[MSDIR]) {
    if (parse.nonOptionsCount() > 0) {
      msfile = String(parse.nonOption(0));
    } else {
      // MSDIR must have been set
      msfile = String(options[MSDIR].last()->arg);
    }
  } else if (!showversion) {
    error("Error: Need to provide name of input Measurement Set."); 
  }

  if (parse.nonOptionsCount() > 1 || options[ASDMDIR]) {
    if (parse.nonOptionsCount() > 1) {
      asdmfile = String(parse.nonOption(1));
    } else {
      // ASDMDIR must have been set
      asdmfile = String(options[ASDMDIR].last()->arg);
    }
  } else if(!showversion){
    error("Error: Need to provide name of output ASDM."); 
  }

  MeasurementSet* itsMS=0;
   
  int rstat = 0; // return value 0 means "everything OK"
  MS2ASDM* m2a=0;
  try{
    if(showversion && msfile==""){
      MeasurementSet ms;
      m2a = new MS2ASDM(ms);
      error("Using ASDM version "+ m2a->showversion());
      delete m2a;
    }
    else{
      itsMS = new MeasurementSet(msfile);
      m2a = new MS2ASDM(*itsMS);
      info("Using ASDM version " + m2a->showversion());
      if (!m2a->writeASDM(asdmfile, datacolumn, archiveid, rangeid, verbose,
			  subscanduration, schedblockduration, apcorrected)) {
	delete m2a;
	delete itsMS;
	error("Conversion to ASDM failed.");
      }
      delete m2a;
      delete itsMS;
    }
  } catch (AipsError x) {
    if(m2a){
      delete m2a;
    }
    if(itsMS){
      delete itsMS;
    }
    Table::relinquishAutoLocks();
    error("Exception Reported: " + x.getMesg());
  }
  Table::relinquishAutoLocks();
  exit(rstat); // note: "return rstat" does not give the correct return value when used inside Python os.system()
  
}


