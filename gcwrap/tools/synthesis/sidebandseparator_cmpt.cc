/*
 * sidebandseparator_cmpt.cc
 *
 *  Created on: 2017/08/01
 *      Author: kana
 */

#include <iostream>
#include <vector>
#include <sidebandseparator_cmpt.h>
#include <synthesis/MeasurementEquations/SideBandSeparator.h>

using namespace std;
using namespace casacore;
using namespace casa;

namespace casac {

sidebandseparator::sidebandseparator()
{

  itsSep=0;
  itsLog = new LogIO();

}

sidebandseparator::~sidebandseparator()
{

  if(itsSep !=0)
    delete itsSep;

  delete itsLog;

}

bool
sidebandseparator::open(const vector<string>& imagename)
{
  Bool rstat(false);
  try {

    // In case already open, close it!
    close();

    itsSep=new SideBandSeparatorII(imagename);
    if(itsSep !=0)
      rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}


bool
sidebandseparator::close()
{
  Bool rstat(false);

  try {

    if(itsSep !=0)
      delete itsSep;
    itsSep=0;
    rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;


}

bool
sidebandseparator::done()
{
	return close();
}

bool
sidebandseparator::setshift(const vector<double> &shift, const bool signal)
{
  Bool rstat(false);
  try {
	  itsSep->setShift(shift, signal);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
sidebandseparator::setlimit(const double limit)
{
  Bool rstat(false);
  try {
	  itsSep->setThreshold(limit);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
sidebandseparator::setboth(const bool getbothside)
{
  Bool rstat(false);
  try {
	  itsSep->solveBoth(getbothside);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
sidebandseparator::set_imageband_frequency(const double refpix, const Quantity& refval)
{
  Bool rstat(false);
  try {
	  itsSep->setImageBandFrequency(refpix, casaQuantity(refval));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
sidebandseparator::setsolveother(const bool subtract_from_other)
{
  Bool rstat(false);
  try {
	  itsSep->solvefromOther(subtract_from_other);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
sidebandseparator::separate(const std::string& outfile, const bool overwrite)
{
  Bool rstat(false);
  try {
	  itsSep->separate(outfile, overwrite);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

} // casac namespace

