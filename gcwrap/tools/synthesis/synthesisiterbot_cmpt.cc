/***
 * Framework independent implementation file for imager...
 *
 * Implement the imager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 
 * @author Wes Young
 * @version 
 ***/

#include <iostream>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <ms/MeasurementSets.h>
#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <casa/Logging/LogIO.h>

#include <synthesis/ImagerObjects/SynthesisIterBot.h>

#include <synthesisiterbot_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

  synthesisiterbot::synthesisiterbot():
    itsIterBot(NULL)
{
  //  itsIterBot = new SynthesisIterBot() ;
  itsIterBot = new SynthesisIterBotWithOldGUI() ;
}

synthesisiterbot::~synthesisiterbot()
{
  if(itsIterBot)
    {
      delete itsIterBot;
    }
}


casac::record* synthesisiterbot::setupiteration(const casac::record& iterpars)
{

  try 
    {
      casacore::Record recpars = *toRecord( iterpars );
      itsIterBot->setupIteration( recpars );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return getiterationdetails();
}


casac::record* synthesisiterbot::getiterationdetails()
{
  casac::record* rstat(0);

  try 
    {
      rstat=fromRecord(itsIterBot->getIterationDetails());
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}

casac::record* synthesisiterbot::pauseforinteraction()
{
  casac::record* rstat;

  try 
    {
      rstat = fromRecord( itsIterBot->pauseForUserInteractionOld() );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}
/*
std::vector<int> synthesisiterbot::pauseforinteraction()
{
  std::vector<int> rstat;

  try 
    {
      Vector<Int> retcode = itsIterBot->pauseForUserInteraction() ;
      rstat.resize( retcode.nelements() );
      for (uInt i=0;i<retcode.nelements();i++) rstat[i]=retcode[i];
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}
*/
casac::record* synthesisiterbot::getiterationsummary()
{
  casac::record* rstat(0);

  try 
    {
      rstat=fromRecord(itsIterBot->getIterationSummary());
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}



  int synthesisiterbot::cleanComplete(const bool lastcyclecheck)
{
  Int rstat=0;

  try 
    {
      rstat = itsIterBot->cleanComplete( lastcyclecheck );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}



bool synthesisiterbot::endmajorcycle()
{
  Bool rstat(false);
  
  try 
    {
      itsIterBot->endMajorCycle();
     } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

bool synthesisiterbot::resetminorcycleinfo()
{
  Bool rstat(false);
  
  try 
    {
      itsIterBot->resetMinorCycleInfo();
     } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}



  casac::record* synthesisiterbot::getminorcyclecontrols()
{
  casac::record* rstat(0);
  try {
    rstat=fromRecord(itsIterBot->getSubIterBot());
  } catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}  

bool synthesisiterbot::mergeinitrecord(const casac::record& initrecord)
{
  Bool rstat(false);
  
  try 
    {
      casacore::Record recpars = *toRecord( initrecord );
      itsIterBot->startMinorCycle( recpars );
     } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

bool synthesisiterbot::mergeexecrecord(const casac::record& execrecord)
{
  Bool rstat(false);
  
  try 
    {
      casacore::Record recpars = *toRecord( execrecord );
      itsIterBot->endMinorCycle(recpars);
     } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  bool synthesisiterbot::changestopflag(const bool stopflag)
{
  Bool rstat(true);
  
  try 
    {
      itsIterBot->changeStopFlag(stopflag);
     } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}


bool
synthesisiterbot::done()
{
  Bool rstat(false);

  try 
    {
      if (itsIterBot)
	{
	  delete itsIterBot;
	  itsIterBot=NULL;
	}
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}



} // casac namespace
