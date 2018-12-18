#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <assert.h>
#include <cmath>
#include <complex>
#include <string>
#include <regex>

#include <stdcasa/optionparser.h>
#include <alma/Options/AlmaArg.h>
using namespace alma;

#include <alma/ASDM/ASDMAll.h>
#include <alma/ASDM/Misc.h>

#include <alma/ASDMBinaries/SDMBinData.h>
using namespace sdmbin;

#include <exception>
using namespace asdm;
#include <alma/ASDM/IllegalAccessException.h>

#include <alma/ASDMBinaries/SDMDataObjectReader.h>
#include <alma/ASDMBinaries/SDMDataObject.h>

#include <asdmstman/AsdmStMan.h>

#include <casa/OS/Path.h>


void scanBDF(const string & pathToBDF) {
  SDMDataObjectStreamReader sdosr;
  
  sdosr.open(pathToBDF);
  ProcessorType processorType = sdosr.processorType();
  if (processorType == RADIOMETER) {
    // returned SDMDataSubset value is not used, ok to ignore
    sdosr.getSubset();
  }
  else if (processorType == CORRELATOR) {
    while (sdosr.hasSubset()) {
      // returned SDMDataSubset value is not used, ok to ignore
      sdosr.getSubset();
    }
  }
  else 
    cout << "Processor not supported in lazy mode." << endl;
  
  sdosr.close(); 
}


int main ( int argc, char * argv[] ) {
  string dsName;
  string appName = string(argv[0]);

  enum optionIndex { UNKNOWN, HELP, ASDMDIR};

  try {
    // remove the program name
    argc--;
    argv++;

    string usageIntro = 
      "Read sequentially all the BDFs of  ASDM dataset without any processing. "
      "It's just an application to measure the time reading the BDFs by using the class SDMDataObjectStreamReader. \n"
      "Usage : " + appName +" asdm-directory \n\n"
      "Command parameters: \n";

    // Descriptor elements are: OptionIndex, OptionType, shortopt, longopt, check_arg, help
    option::Descriptor usage[] = {
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  usageIntro.c_str()},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, " \tasdm-directory :  \tthe pathname to the ASDM dataset containing the BFDs to be read."},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, "\nAllowed options:\n"},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown,  0 }, // helps with formatting
      // help is the only visible option
      { HELP, 0, "", "help", AlmaArg::None, " --help  \tproduces help message."},
      // Hidden options, will be allowed both on command line and
      // in config file, but will not be shown to the user.
      { ASDMDIR, 0, "", "asdm-directory", AlmaArg::Required, 0},
      { 0, 0, 0, 0, 0, 0} };

    // there are no defaults

    // parse the command line
    // establish sizes
    option::Stats stats;
    stats.add(true, usage, argc, argv);

    // size the buffers to hold the parsed options
    option::Option options[stats.options_max], buffer[stats.buffer_max];

    // and parse things
    option::Parser parse;
    parse.parse(true, usage, argc, argv, options, buffer);

    if (parse.error()) {
      cerr << "Problem parsing the command line arguments" << endl;
      exit(1);
    }

    if (options[HELP] || (argc==0)) {
      option::printUsage(cout, usage, 80);
      exit(1);
    }

    if (parse.nonOptionsCount() > 0 || options[ASDMDIR]) {
      if (parse.nonOptionsCount() > 0) {
	dsName = string(parse.nonOption(0));
      } else {
	dsName = string(options[ASDMDIR].last()->arg);
      }
      cout << "dsName = " << dsName << endl;
      trim(dsName);
      if (dsName.back()=='/') dsName.pop_back();
    } else {
      option::printUsage(cout, usage, 80);
      exit (1);
    }
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
  }

  ASDM* ds_p = new ASDM();
  try {
    cout << "Input ASDM dataset : " << dsName << endl;
    
    ds_p->setFromFile(dsName, ASDMParseOptions().loadTablesOnDemand(true));
  }
  catch (ConversionException e) {
    cout << e.getMessage() << endl;
  }
  catch (std::exception e) {
    cout << e.what() << endl;
  }
  catch (...) {
    cout << "Uncaught exception !" << endl;
  } 

  const MainTable& mainT = ds_p->getMain();
  const vector<MainRow *>& v = mainT.get();

  vector<string> bdfNames;
  for ( vector<MainRow *>::const_iterator iter_v = v.begin(); iter_v != v.end(); iter_v++) {
    string dataUID = (*iter_v)->getDataUID().getEntityId().toString();
    replace(dataUID.begin(),dataUID.end(),':','_');
    replace(dataUID.begin(),dataUID.end(),'/','_');
    string abspath = casacore::Path(dsName + "/ASDMBinary/" + dataUID).absoluteName();
    cout << abspath << endl;
    bdfNames.push_back(abspath);
  }

  for (unsigned int iBDF = 0; iBDF < bdfNames.size(); iBDF++) {
    scanBDF(bdfNames[iBDF]);
    cout << iBDF+1 << "/" << bdfNames.size() << endl;
  }
  exit (0);
}
