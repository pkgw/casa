/***
 * Framework independent implementation file for imager...
 *
 * Implement the imager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 
 ***/

#include <iostream>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <casa/Logging/LogIO.h>

#include <casa/OS/Directory.h>
#include <images/Images/PagedImage.h>
#include <coordinates/Coordinates/CoordinateUtil.h>

#include <synthesis/ImagerObjects/SDMaskHandler.h>

#include <synthesismaskhandler_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

synthesismaskhandler::synthesismaskhandler()
{
  itsLog = new LogIO();
  itsMaskHandler = new SDMaskHandler();
}

synthesismaskhandler::~synthesismaskhandler()
{
  done();
}

  casac::record* synthesismaskhandler::pruneregions(const string& inmaskname, double prunesize, const std::vector<bool>& chanflag, const string& outmaskname)
{

  casac::record* rstat(0);
  *itsLog << casacore::LogOrigin("synthesismaskhandler", __func__);
  try 
    {
      PagedImage<Float> inmask(inmaskname);
      CoordinateSystem inCsys = inmask.coordinates();
      IPosition inShape = inmask.shape();
      
      String outmasknameMod;
      if (outmaskname=="") {
        outmasknameMod = String(inmaskname)+String(".pruned");
      }
      else {
       outmasknameMod = outmaskname;
      }
        
      PagedImage<Float> outmask(TiledShape(inShape), inCsys, outmasknameMod);
      Vector<Bool> chanFlag(chanflag);
      if (chanFlag.nelements()==0) {
        *itsLog<<LogIO::SEVERE<<"Please specify chanflag."<<LogIO::POST;
        return(rstat);
      }
      Int inSpecAxis = CoordinateUtil::findSpectralAxis(inCsys); 
      if (inSpecAxis != -1 && chanFlag.nelements()!=uInt(inShape(inSpecAxis))) {
        *itsLog<<LogIO::SEVERE<<"The number elements in chanflag does not match with input mask shape."<<LogIO::POST;
        return(rstat);
      }
      Vector<Bool> allpruned;
      Vector<uInt> nreg;
      Vector<uInt> npruned;
      std::shared_ptr<ImageInterface<Float> > tempIm_ptr = itsMaskHandler->YAPruneRegions(inmask, chanFlag, allpruned, nreg, npruned, prunesize);
      *itsLog<<"nreg="<<nreg<<" npruned="<<npruned<<" prunesize="<<prunesize<<LogIO::POST; 
      outmask.copyData(*(tempIm_ptr.get()));
      casacore::Record outinfo;
      Vector<Int> nregForReport(nreg.nelements());
      Vector<Int> nprunedForReport(npruned.nelements());
      for (uInt i = 0; i < chanFlag.nelements();i++) {
        if (chanFlag[i]) {
          nregForReport[i] = -1;
          nprunedForReport[i] = -1;
        }
        else {
          nregForReport[i] = (Int) nreg[i];
          nprunedForReport[i] = (Int) npruned[i];
        }
      }
      outinfo.define("N_reg",nregForReport);
      outinfo.define("N_reg_pruned",nprunedForReport);
      outinfo.define("prunesize",prunesize);
 
      rstat = fromRecord(outinfo);

    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  bool synthesismaskhandler::done()
{
  Bool rstat(false);

  try
    {
      if (itsMaskHandler)
        {
          delete itsMaskHandler;
          itsMaskHandler=NULL;
        }
      if (itsLog)
        {
          delete itsLog;
          itsLog=NULL;
        }
      rstat=true;
    }
  catch  (AipsError x)
    {
      RETHROW(x);
    }

  return rstat;
}



} // casac namespace
