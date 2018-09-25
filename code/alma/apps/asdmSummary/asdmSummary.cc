#include <alma/ASDM/ASDM.h>
using namespace asdm;

#include <alma/ASDM/AntennaRow.h>
#include <alma/ASDM/AntennaTable.h>
#include <alma/ASDM/ConfigDescriptionRow.h>
#include <alma/ASDM/ConfigDescriptionTable.h>
#include <alma/ASDM/DataDescriptionRow.h>
#include <alma/ASDM/DataDescriptionTable.h>
#include <alma/ASDM/ExecBlockRow.h>
#include <alma/ASDM/ExecBlockTable.h>
#include <alma/ASDM/MainRow.h>
#include <alma/ASDM/MainTable.h>
#include <alma/ASDM/PolarizationRow.h>
#include <alma/ASDM/PolarizationTable.h>
#include <alma/ASDM/ScanRow.h>
#include <alma/ASDM/ScanTable.h>
#include <alma/ASDM/SpectralWindowRow.h>
#include <alma/ASDM/SpectralWindowTable.h>
#include <alma/ASDM/StationRow.h>
#include <alma/ASDM/StationTable.h>
#include <alma/ASDM/SubscanRow.h>
#include <alma/ASDM/SubscanTable.h>

#include <alma/Enumerations/CAntennaMake.h>
using namespace AntennaMakeMod;

#include <alma/Enumerations/CAtmPhaseCorrection.h>
using namespace AtmPhaseCorrectionMod;

#include <alma/Enumerations/CCorrelationMode.h>
using namespace CorrelationModeMod;

#include <alma/Enumerations/CStokesParameter.h>
using namespace StokesParameterMod;

#include <alma/Enumerations/CFrequencyReferenceCode.h>
using namespace FrequencyReferenceCodeMod;

#include <alma/Enumerations/CScanIntent.h>
using namespace ScanIntentMod;

#include <alma/Enumerations/CSpectralResolutionType.h>
using namespace SpectralResolutionTypeMod;

#include <alma/Enumerations/CSubscanIntent.h>
using namespace SubscanIntentMod;

#include <alma/Enumerations/CTimeSampling.h>
using namespace TimeSamplingMod;

#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>
using namespace casacore;

using namespace std;

#include <stdcasa/optionparser.h>
#include <alma/Options/AlmaArg.h>
using namespace alma;

string appName;

// A facility to get rid of blanks at start and end of a string.
// 
string lrtrim(std::string& s,const std::string& drop = " ")
{
  std::string r=s.erase(s.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}

void info (const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void error(const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
  //os << LogIO::POST;
  exit(1);
}

#include <iostream>
#include <sstream>

ostringstream errstream;
ostringstream infostream;
using namespace std;

void antennaSummary(const ExecBlockRow* eb_p) {
  ASDM& ds = eb_p->getTable().getContainer();
  AntennaTable& aT = ds.getAntenna();
  StationTable& sT = ds.getStation();

  infostream.str("");
  const vector<Tag> antennaIds = eb_p->getAntennaId();
  AntennaRow * antenna_p = NULL ;
  StationRow * station_p = NULL ;
  infostream << endl;
  infostream << antennaIds.size() << " antennas have been used in this exec block." << endl;
  infostream << "        Id     Name         Make Station    Diameter         X              Y             Z" << endl;
  for (unsigned int i = 0; i < antennaIds.size(); i++) {
    antenna_p = aT.getRowByKey(antennaIds[i]);
    station_p = sT.getRowByKey(antenna_p->getStationId());
    vector<Length> position = station_p->getPosition();
    //infostream.fill(''); 
    infostream.width(12);infostream << antenna_p->getAntennaId() ;
    infostream.width(6); infostream.setf(ios::right); infostream   << antenna_p->getName() ;
    infostream.width(13); infostream  << CAntennaMake::name(antenna_p->getAntennaMake()) ;
    infostream.width(6); infostream   << station_p->getName() ;
    infostream.width(10); infostream.precision(10); infostream << antenna_p->getDishDiameter() ;
    infostream.width(15); infostream.precision(10); infostream << position[0] ;
    infostream.width(15); infostream.precision(10); infostream << position[1] ;
    infostream.width(15); infostream.precision(10); infostream << position[2] << endl;
  }
    info(infostream.str());
}

template<typename Enum, typename EnumHelper> 
void output1 (typename  vector<Enum>::iterator begin, typename vector<Enum>::iterator end, ostringstream & oss) {
  if (begin == end) return;
  oss << ',' << EnumHelper::name(*begin);
  output1<Enum, EnumHelper>(begin+1, end, oss);    
} 

template<typename Enum, typename EnumHelper> 
void output (typename std::vector<Enum>::iterator begin, typename std::vector<Enum>::iterator end, std::ostringstream & oss) {
  if (begin == end) return;
  oss << EnumHelper::name(*begin);
  output1<Enum, EnumHelper>(begin+1, end, oss);    
} 


template<typename T>
void output1 (typename vector<T>::iterator begin, typename vector<T>::iterator end, ostringstream& oss) {
  if (begin == end) return;
  oss << "," << *begin;
  output1<T> (begin+1, end, oss);
}

template<typename T>
void output (typename vector<T>::iterator begin, typename vector<T>::iterator end, ostringstream& oss) {
  if (begin == end) return;
  oss << *begin;
  output1<T> (begin+1, end, oss);
}

typedef struct SpectralWindowSummaryStruct {
  int		numChan;
  string	measFreqRef;
  Frequency     firstChan;
  Frequency	chanWidth;
  Frequency	refFreq;
} SpectralWindowSummary;

SpectralWindowSummary spectralWindowSummary(SpectralWindowRow * spw_p) {
  SpectralWindowSummary result;
  
  result.numChan = spw_p->getNumChan();
  
  if (spw_p->isChanFreqStartExists()) 
    result.firstChan = spw_p->getChanFreqStart();
  else 
    if (spw_p->isChanFreqArrayExists())
      result.firstChan = spw_p->getChanFreqArray()[0];
    else
      result.firstChan = Frequency(0.0);

  if (spw_p->isChanWidthArrayExists())
    result.chanWidth = spw_p->getChanWidthArray()[0];
  else 
    if (spw_p->isChanWidthExists())
      result.chanWidth = spw_p->getChanWidth();
    else
      result.chanWidth = Frequency(0.0);

  if (spw_p->isMeasFreqRefExists())
    result.measFreqRef = CFrequencyReferenceCode::name(spw_p->getMeasFreqRef());
  else
    result.measFreqRef = "TOPO";

  result.refFreq = spw_p->getRefFreq();
  
  return result;
}

bool notNull(int n) { return n != 0 ; }

void mainSummary(ExecBlockRow* eb_p, int scanNumber, int subscanNumber) {

  ASDM& ds = eb_p->getTable().getContainer();
  
  Tag ebId = eb_p->getExecBlockId();

  const vector<MainRow *>& mains = ds.getMain().get();
  vector<MainRow *> eb_mains;

  for(MainRow* main_p: mains) {
    if ( main_p->getExecBlockId() == ebId && main_p->getScanNumber() == scanNumber && main_p->getSubscanNumber() == subscanNumber )
      eb_mains.push_back(main_p);
  }

  DataDescriptionTable& ddT = ds.getDataDescription();
  PolarizationTable& polT = ds.getPolarization();
  SpectralWindowTable& spwT = ds.getSpectralWindow();
  ConfigDescriptionTable& cfgDescT = ds.getConfigDescription();

  for ( MainRow* main_p: eb_mains ) {
    infostream.str("");
    infostream << endl;
    infostream << "\t\t Binary data in " << main_p->getDataUID().getEntityId() << endl;
    infostream << "\t\t Number of integrations : " << main_p->getNumIntegration() << endl;
    infostream << "\t\t Time sampling : " << CTimeSampling::name(main_p->getTimeSampling()) << endl;
    ConfigDescriptionRow* cfgDesc_p = cfgDescT.getRowByKey(main_p->getConfigDescriptionId());
    infostream << "\t\t Correlation Mode : " << CCorrelationMode::name(cfgDesc_p->getCorrelationMode()) << endl;
    infostream << "\t\t Spectral resolution type : " << CSpectralResolutionType::name(cfgDesc_p->getSpectralType()) << endl;
    infostream << "\t\t Atmospheric phase correction : " ;
    vector<AtmPhaseCorrection> apcs = cfgDesc_p->getAtmPhaseCorrection();
    output<AtmPhaseCorrection, CAtmPhaseCorrection>(apcs.begin(), apcs.end(), infostream);
    infostream << endl;
    info(infostream.str());

    vector<Tag> ddIds = cfgDesc_p->getDataDescriptionId();
    for ( Tag ddId: ddIds ) {
      DataDescriptionRow * dd_p = ddT.getRowByKey(ddId);
      SpectralWindowRow * spw_p = spwT.getRowByKey(dd_p->getSpectralWindowId());
      PolarizationRow * p_p = polT.getRowByKey(dd_p->getPolOrHoloId());
      infostream.str("");
      SpectralWindowSummary spwSummary = spectralWindowSummary(spw_p);
      infostream << "\t\t " << spw_p->getSpectralWindowId() << " : numChan = " << spwSummary.numChan
		 << ", frame = " << spwSummary.measFreqRef
		 << ", firstChan = " << spwSummary.firstChan
		 << ", chandWidth = " << spwSummary.chanWidth
		 << " x " 
		 << p_p->getPolarizationId() << " : corr = " ; 
      vector<StokesParameter> corrType = p_p->getCorrType();
      output<StokesParameter, CStokesParameter>(corrType.begin(), corrType.end(), infostream);
      infostream << endl;
      info(infostream.str());
    }
  }
}

void subscanSummary(ExecBlockRow* eb_p, int scanNumber) {
      
  ASDM& ds = eb_p->getTable().getContainer();
  Tag ebId = eb_p->getExecBlockId();

  const vector<SubscanRow *>& subscans = ds.getSubscan().get();
  vector<SubscanRow *> eb_subscans;
  for (SubscanRow * sscan_p: subscans) {
    if (sscan_p->getExecBlockId() == ebId && sscan_p->getScanNumber() == scanNumber) 
      eb_subscans.push_back(sscan_p);
  }

  for (SubscanRow* sscan_p: eb_subscans) {
    infostream.str("");
    infostream << "\tSubscan #" << sscan_p->getSubscanNumber()
	       << " from " << sscan_p->getStartTime().toFITS()
	       << " to " << sscan_p->getEndTime().toFITS()
	       << endl;
    infostream << "\t\tIntent : " << CSubscanIntent::name(sscan_p->getSubscanIntent()) << endl;
    infostream << "\t\tNumber of integrations : " << sscan_p->getNumIntegration() << endl;
    vector<int> numSubintegration = sscan_p->getNumSubintegration();
    if (find_if(numSubintegration.begin(), numSubintegration.end(), notNull) != numSubintegration.end()) {
      infostream << "\t\tNumber of subintegrations per integration : ";
      output<int>(numSubintegration.begin(), numSubintegration.end(), infostream);
      infostream << endl;
    }
    info(infostream.str()); 

    mainSummary(eb_p, scanNumber, sscan_p->getSubscanNumber());
  }

}


void scanSummary(ExecBlockRow* eb_p) {

  ASDM& ds = eb_p->getTable().getContainer();  
  Tag ebId = eb_p->getExecBlockId();

  const vector<MainRow *>& mains = ds.getMain().get();
  vector<MainRow *> eb_mains;

  for(MainRow* main: mains) {
    if ( main->getExecBlockId() == ebId) eb_mains.push_back(main);
  }
  
  const vector<ScanRow*>& scans = ds.getScan().get();
  vector<ScanRow *> eb_scans;
  for(ScanRow* scan: scans) {
    if ( scan->getExecBlockId() == ebId) eb_scans.push_back(scan);
  }

  infostream.str("");
  infostream << endl;
  infostream << "Number of scans in this exec Block : " << eb_scans.size() << endl;
  info(infostream.str());
  if (eb_scans.size() > 0) {
    for (ScanRow* scan_p: eb_scans) {
      infostream.str("");
      infostream << endl;
      infostream << "scan #" << scan_p->getScanNumber()
		 << " from " << scan_p->getStartTime().toFITS()
		 << " to " <<  scan_p->getEndTime().toFITS()
		 << endl;
      
      vector<ScanIntent> scis = scan_p->getScanIntent();
      if (scis.size() > 0) {
	infostream << "\tIntents : ";
	output<ScanIntent, CScanIntent>(scis.begin(), scis.end(), infostream);
	infostream << endl;
      }
      
      if ( scan_p->isFieldNameExists() ) {
	vector<string> fields = scan_p->getFieldName();
	if (fields.size() > 0) {
	  infostream << "\tFields : ";
	  output<string>(fields.begin(), fields.end(), infostream);
	  infostream << endl;
	}
      }
      
      if ( scan_p->isSourceNameExists() ) {
	infostream << "\tSources : " << scan_p->getSourceName() << endl; 
      }
      info(infostream.str());
      subscanSummary(eb_p, scan_p->getScanNumber());
    }
  }
}

void execBlockSummary(const ASDM& ds) {
  infostream.str("");

  const vector<ExecBlockRow*>& ebs = ds.getExecBlock().get();
  for (unsigned int i = 0; i < ebs.size(); i++) {
    ExecBlockRow* eb_p = ebs[i];
    infostream << "\n";
    infostream << "Exec Block : " << eb_p->getExecBlockId() << endl;
    infostream << "Telescope : " << eb_p->getTelescopeName() << endl;
    infostream << "Configuration name : " << eb_p->getConfigName() << endl;
    infostream << "Observer name : " << eb_p->getObserverName() << endl;
    infostream << "The exec block started on " << eb_p->getStartTime().toFITS() << " and ended on " << eb_p->getEndTime().toFITS() << endl;
    if (eb_p->getAborted())
      infostream << "It was aborted." << endl;
    info(infostream.str());

    antennaSummary(eb_p);
    scanSummary(eb_p);
  }
}


void summary(const ASDM& ds, const string& dsPath) {
  infostream.str("");
  infostream << "========================================================================================" << endl;
  infostream << "ASDM dataset :" << dsPath << endl;
  infostream << "========================================================================================" << endl;
  info(infostream.str());

  execBlockSummary(ds);

}

int main (int argc, char* argv[]) {
  string dsName;
  string appName = string(argv[0]);
  ofstream ofs;

  // process command line otpions and parameters
  // asdm-directory may be given as a named parameter=value, but that usage is hidden from the user
  // the positional version of asdm-directory takes precedence over any named parameter use when both are present

  enum optionIndex { UNKNOWN, HELP, LOGFILE, ASDMDIR };

  try {
    // remove the program name
    argc--;
    argv++;

    string usageIntro = 
      "Displays a summary of the content of an ASDM dataset .\n"
      "Usage : " + appName +" + [optiopns] asdm-directory \n\n"
      "Command parameters: \n";

    // Descriptor elements are : OptionIndex, OptionType, shortopt, longopt, check_arg, help
    option::Descriptor usage[] = {
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, usageIntro.c_str()},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, " \tasdm-directory : \tthe pathname to the ASDM dataset to be reported on."},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, "\nAllowed options:\n"},
      { UNKNOWN, 0, "", "", AlmaArg::Unknown, 0}, // helps with formatting

      // these are the non-positional options
      { HELP, 0, "", "help", AlmaArg::None, " --help  \tproduces this help message."},
      { LOGFILE, 0, "l", "logfile", AlmaArg::Required, 
	" -l [--logfile] arg \tspecifies the log filename. "
	"If the option is not used then the logged informations are written to the standard error stream."},

      // hidden non-positional option, the positional one takes precedence when both are set
      { ASDMDIR, 0, "", "asdm-directory", AlmaArg::Required, 0},
      { 0, 0, 0, 0, 0, 0} };
 
    // there are no defaults to set here

    // parse argv
    // establish sizes
    option::Stats stats;
    // true here turns on re-ordering of args so that positional arguments are always seen last
    stats.add(true, usage, argc, argv);

    // buffers to hold the parsed options
    // options has one element per optionIndex, last value is the last time that option was set
    // buffer has one element for each option encountered, in order. Not used here.
    option::Option *options = new option::Option[stats.options_max];
    option::Option *buffer = new option::Option[stats.buffer_max];

    option::Parser parse;
    // true has the same meaning as in stats above. This may not be necessary here.
    parse.parse(true, usage, argc, argv, options, buffer);

    if (parse.error()) {
      errstream.str("");
      errstream << "Problem parsing the command line arguments";
      error(errstream.str());
    }

    // User-specified logfile?
    cout << "checking on LOGFILE" << endl;
    if (options[LOGFILE] != NULL && options[LOGFILE].last()->arg != NULL) {
      cout << " arg = " << options[LOGFILE].last()->arg << endl;
      //LogSinkInterface *theSink;
      ofs.open(options[LOGFILE].last()->arg, ios_base::app);
      LogSinkInterface *theSink = new casacore::StreamLogSink(&ofs);
      LogSink::globalSink(theSink);
    }

    // Help ? display help's content and don't go further
    if (options[HELP] || (argc==0)) {
      errstream.str("");
      option::printUsage(errstream, usage, 80);
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
      dsName = trim_copy(dummy) ;
      if (dsName.back()=='/') dsName.erase(dsName.size()-1);
    }
    else {
      errstream.str("");
      option::printUsage(errstream, usage);
      error(errstream.str());
    }
  }
  catch (std::exception& e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }

  try {
    infostream.str("");
    infostream << "Input ASDM dataset : " << dsName << endl;
    info(infostream.str());
    
    ASDM ds ;
    ds.setFromFile(dsName, ASDMParseOptions().loadTablesOnDemand(true).checkRowUniqueness(false));
    summary(ds, dsName);
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
  return 0;
}

  
