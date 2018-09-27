
/***
 * Framework independent implementation file for imagepol...
 *
 * Implement the imagepol component here.
 * 
 * // imagepol_cmpt.cc: defines imagepol class which implements functionity
 * // for image polarimetry.
 * @author
 * @version 
 ***/

#include <iostream>
#include <imagepol_cmpt.h>
#include <casa/Containers/Record.h>
#include <casa/Exceptions/Error.h>
#include <casa/Logging/LogIO.h>
#include <images/Images/ImageUtilities.h>
#include <images/Images/ImageExpr.h>

#include <imageanalysis/ImageAnalysis/ImageFactory.h>
#include <imageanalysis/ImageAnalysis/ImagePolProxy.h>
#include <imageanalysis/ImageAnalysis/ImageTotalPolarization.h>
#include <imageanalysis/ImageAnalysis/LinearPolarizationAngleCalculator.h>
#include <imageanalysis/ImageAnalysis/LinearPolarizationCalculator.h>

#include <casa/namespace.h>

using namespace std;

using namespace casacore;
using namespace casa;

namespace casac {

imagepol::imagepol()
{
  itsImPol=0;
  itsLog=0;

  itsLog= new LogIO();
}

imagepol::~imagepol()
{

}

bool 
imagepol::open(const variant& image){
	try {
		*itsLog << LogOrigin("imagepol", "open");
		if(itsLog==0) {
			itsLog=new LogIO();
		}
		if(itsImPol) delete itsImPol;
		if(image.type()==variant::RECORD){
			variant localvar(image);
			std::auto_ptr<Record> tmp(
				toRecord(localvar.asRecord())
			);
			TempImage<Float> tmpim;
			String err;
			if (! tmpim.fromRecord(err, *tmp)) {
				*itsLog << LogIO::SEVERE << "Could not convert image record"
						<< LogIO::EXCEPTION;
			}
			itsImPol= new ImagePol(tmpim);
		}
		else if(image.type()== variant::STRING) {
			//PtrHolder<ImageInterface<Float> > im;
			//ImageUtilities::openImage(im, toCasaString(image));
            SPIIF imageF;
            std::tie(imageF, std::ignore, std::ignore, std::ignore)
                = ImageFactory::fromFile(toCasaString(image));
            ThrowIf(! imageF, "Can only open float valued images");
			itsImPol= new ImagePol(*imageF);
		}
		else {
			*itsLog << "Unsupported type for image input" << LogIO::EXCEPTION;
		}
		return true;
	}
	catch (AipsError x) {
		*itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
		RETHROW(x);
	}
}

bool
imagepol::imagepoltestimage(const std::string& outfile,
			    const std::vector<double>& rm,
			    const double pa0, const double sigma,
			    const int nx, const int ny, const int nf,
			    const double f0, const double bw) {

  bool rstat(false);
  bool rmdefault(false);

  try {
    *itsLog << LogOrigin("imagepol", "imagepoltestimage");
    if(itsImPol) delete itsImPol;
    itsImPol = new ImagePol();
    if(itsLog==0) itsLog=new LogIO();

    if (rm.size() == 1 and rm[0]==0.0) {
      rmdefault = true;
    }
    rstat = itsImPol->imagepoltestimage(outfile, rm, rmdefault, pa0, sigma,
    					nx, ny, nf, f0, bw);
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: "
		    << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool 
imagepol::close() {
  bool rstat(false);
  try {
    *itsLog << LogOrigin("imagepol", "close");

    if (itsImPol != 0) {
      delete itsImPol;
      itsImPol = 0;
    }
    rstat = true;

  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
imagepol::complexlinpol(const std::string& outfile)
{
  bool rstat(false);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rstat;
    }
    rstat = itsImPol->complexlinpol(String(outfile));
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
imagepol::complexfraclinpol(const std::string& outfile)
{

  bool rstat(false);
  try{
    *itsLog << LogOrigin("imagepol", "complexfraclinpol");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rstat;
    }
    itsImPol->complexFractionalLinearPolarization(String(outfile));
    rstat = true;
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

image *
imagepol::depolratio(const std::string& infile, const bool debias, const double clip, const double sigma, const std::string& outfile)
{

  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->depolratio(out, infile, debias, clip, sigma, outfile);
    if (rstat) {
        std::shared_ptr<ImageInterface<Float> > x(out);
        return new image(x);
    }
    else {
      throw(AipsError("could not attach depolratio image"));
    }
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

bool
imagepol::done()
{
  bool rstat(false);
  try {
    *itsLog << LogOrigin("imagepol", "done");

    if (itsImPol != 0) {
      delete itsImPol;
      itsImPol = 0;
    }
    rstat = true;

  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
imagepol::fourierrotationmeasure(
    const std::string& complex, const std::string& amp,
    const std::string& pa, const std::string& real,
    const std::string& imag, const bool zerolag0
) {
  bool rstat(false);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rstat;
    }
    itsImPol->fourierRotationMeasure(String(complex), String(amp),
				     String(pa), String(real),
				     String(imag), zerolag0);
    rstat = true;
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

image *
imagepol::fraclinpol(const bool debias, const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "fraclinpol");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->fracLinPol(out,debias,Float(clip),
				 Float(sigma),String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach fraclinpol image"));
    }
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::fractotpol(const bool debias, const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->fracTotPol(out,debias,Float(clip),
				 Float(sigma),String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach fractotpol image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image * imagepol::linpolint(
    const bool debias, const double clip, const double sigma, const std::string& outfile,
    const variant& region, const std::string& mask, bool stretch
) {
    try{
        *itsLog << LogOrigin("imagepol", __func__);
        if(! itsImPol){
            *itsLog << LogIO::SEVERE <<"No attached image, please use open "
                << LogIO::POST;
            return nullptr;
        }
        auto myreg = _getRegion(region, False);
        auto subImage = SubImageFactory<Float>::createSubImageRO(
            *itsImPol->getImage(), *myreg, mask,
            itsLog, AxesSpecifier(), stretch, true
        );
        LinearPolarizationCalculator lpc(subImage, outfile, False);
        lpc.setClip(clip);
        lpc.setSigma(sigma);
        lpc.setDebias(debias);
        return new image(lpc.compute());
    }
    catch (const AipsError& x) {
        *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
}

image * imagepol::linpolposang(
    const std::string& outfile, const variant& region,
    const std::string& mask, bool stretch
) {
    try{
        *itsLog << LogOrigin("imagepol", __func__);
        if(! itsImPol){
            *itsLog << LogIO::SEVERE <<"No attached image, please use open "
                << LogIO::POST;
            return nullptr;
        }
        auto myreg = _getRegion(region, False);
        auto subImage = SubImageFactory<Float>::createSubImageRO(
            *itsImPol->getImage(), *myreg, mask,
            itsLog, AxesSpecifier(), stretch, true
        );
        LinearPolarizationAngleCalculator lpac(subImage, outfile, False);
        return new image(lpac.compute(False));
    }
    catch (const AipsError& x) {
        *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
}

bool
imagepol::makecomplex(const std::string& complex, const std::string& real, const std::string& imag, const std::string& amp, const std::string& phase)
{
  bool rstat(false);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rstat;
    }
    itsImPol->makeComplex(String(complex), String(real),
			  String(imag), String(amp),
			  String(phase));
    rstat = true;
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

image *
imagepol::pol(const std::string& which, const bool debias, const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    String type(which);
    type.upcase();
    Bool rstat(false);
    if (type=="LPI") {
        return linpolint(debias, clip, sigma, outfile);
        /*
      rstat = itsImPol->linPolInt(out,debias,Float(clip),
				  Float(sigma),String(outfile));
				  */
    } else if (type=="TPI") {
        return totpolint(debias, clip, sigma, outfile);
    } else if (type=="LPPA") {
        return linpolposang(outfile);
      //rstat = itsImPol->linPolPosAng(out,String(outfile));
    } else if (type=="FLP") {
      rstat = itsImPol->fracLinPol(out,debias,Float(clip),
				   Float(sigma),String(outfile));
    } else if (type=="FTP") {
      rstat = itsImPol->fracTotPol(out,debias,Float(clip),
				   Float(sigma),String(outfile));
    } else {
      throw(AipsError("Code "+String(which)+" is unrecognized"));
    }
      if (rstat) {
          return new image(out);
      }
      else {
      throw(AipsError("could not attach pol image"));
      }    
  } catch (const AipsError& x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

bool imagepol::rotationmeasure(
    const std::string& rm, const std::string& rmerr,
    const std::string& pa0, const std::string& pa0err,
    const std::string& nturns, const std::string& chisq,
    const double sigma, const double rmfq, const double rmmax,
    const double maxpaerr, const std::string&, const int, const int
) {
    try {
        *itsLog << LogOrigin("imagepol", __func__);
        if(! itsImPol) {
            *itsLog << LogIO::SEVERE <<"No attached image, please use open "
                << LogIO::POST;
            return False;
        }
        itsImPol->rotationMeasure(
            rm, rmerr, pa0, pa0err, nturns,
            chisq, -1, sigma, rmfq, rmmax, maxpaerr
        );
        return True;
    }
    catch (const AipsError& x) {
        *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
    return True;
}

double
imagepol::sigma(const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigma(Float(clip));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

image *
imagepol::sigmadepolratio(const std::string& infile, const bool debias, const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->sigmaDepolarizationRatio(out,String(infile), debias,
					       Float(clip), Float(sigma),
					       String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach sigmadepolratio image"));
    }
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

}

image *
imagepol::sigmafraclinpol(const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->sigmaFracLinPol(out, Float(clip), Float(sigma),
				      String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach sigmafraclinpol image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::sigmafractotpol(const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->sigmaFracTotPol(out, Float(clip),
				      Float(sigma), String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach sigmafractotpol image"));
    }
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

}

double
imagepol::sigmalinpolint(const double clip, const double sigma, const std::string& )
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaLinPolInt(clip, sigma);
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

image *
imagepol::sigmalinpolposang(const double clip, const double sigma, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->sigmaLinPolPosAng(out,Float(clip),
					Float(sigma),String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
    throw(AipsError("could not attach sigmalinpolposang image"));
    }
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

}

double
imagepol::sigmastokes(const std::string& which, const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    String type(which);
    type.upcase();
    if (type=="I") {
      rvalue = itsImPol->sigmaStokesI(Float(clip));
    } else if (type=="Q") {
      rvalue = itsImPol->sigmaStokesQ(Float(clip));
    } else if (type=="U") {
      rvalue = itsImPol->sigmaStokesU(Float(clip));
    } else if (type=="V") {
      rvalue = itsImPol->sigmaStokesV(Float(clip));
    } else {
      throw(AipsError("Stokes "+String(which)+" is unrecognized"));
    }
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

double
imagepol::sigmastokesi(const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", "sigmastokesi");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaStokesI(Float(clip));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

double
imagepol::sigmastokesq(const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", "sigmastokesq");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaStokesQ(Float(clip));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

double
imagepol::sigmastokesu(const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", "sigmastokesu");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaStokesU(Float(clip));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

double
imagepol::sigmastokesv(const double clip)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", "sigmastokesv");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaStokesV(Float(clip));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

double
imagepol::sigmatotpolint(const double clip, const double sigma)
{
  double rvalue(-1.0);
  try{
    *itsLog << LogOrigin("imagepol", __func__);
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rvalue;
    }
    rvalue = itsImPol->sigmaTotPolInt(Float(clip),Float(sigma));
  } catch (AipsError x) {
     *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rvalue;
}

image *
imagepol::stokes(const std::string& which, const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "stokes");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    String type(which);
    type.upcase();
    if (type=="I") {
      rstat = itsImPol->stokesI(out,String(outfile));
    } else if (type=="Q") {
      rstat = itsImPol->stokesQ(out,String(outfile));
    } else if (type=="U") {
      rstat = itsImPol->stokesU(out,String(outfile));
    } else if (type=="V") {
      rstat = itsImPol->stokesV(out,String(outfile));
    } else {
      throw(AipsError("Stokes "+String(which)+" is unrecognized"));
    }
    if (rstat) {
        return new image(out);
    } 
    else {
      throw(AipsError("could not attach sigmadepolratio image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::stokesi(const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "stokesi");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->stokesI(out, String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
        throw(AipsError("could not attach stokesi image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::stokesq(const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "stokesq");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->stokesQ(out, String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
    throw(AipsError("could not attach stokesq image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::stokesu(const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "stokesu");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->stokesU(out, String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
      throw(AipsError("could not attach stokesu image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
}

image *
imagepol::stokesv(const std::string& outfile)
{
  try{
    *itsLog << LogOrigin("imagepol", "stokesv");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return 0;
    }
    ImageInterface<Float> *out;
    Bool rstat(false);
    rstat = itsImPol->stokesV(out, String(outfile));
    if (rstat) {
        return new image(out);
    }
    else {
    throw(AipsError("could not attach stokesv image"));
    }
    } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

}

bool
imagepol::summary()
{
  bool rstat(false);
  try{
    *itsLog << LogOrigin("imagepol", "summary");
    if(itsImPol==0){
      *itsLog << LogIO::SEVERE <<"No attached image, please use open " 
	      << LogIO::POST;
      return rstat;
    }
    itsImPol->summary();
    rstat = true;
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

image* imagepol::totpolint(
    bool debias, double clip, double sigma, const std::string& outfile,
    const variant& region, const std::string& mask, bool stretch
) {
    try{
        *itsLog << LogOrigin("imagepol", __func__);
        if(! itsImPol){
            *itsLog << LogIO::SEVERE <<"No attached image, please use open "
                << LogIO::POST;
            return nullptr;
        }
        auto myreg = _getRegion(region, False);
        auto subImage = SubImageFactory<Float>::createSubImageRO(
            *itsImPol->getImage(), *myreg, mask,
            itsLog, AxesSpecifier(), stretch, true
        );
        ImageTotalPolarization itp(subImage, outfile, False);
        itp.setClip(clip);
        itp.setSigma(sigma);
        itp.setDebias(debias);
        return new image(itp.compute());
    }
    catch (const AipsError& x) {
        *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
}

std::shared_ptr<Record> imagepol::_getRegion(
    const variant& region, bool nullIfEmpty
) const {
    switch (region.type()) {
    case variant::BOOLVEC:
        return std::shared_ptr<Record>(nullIfEmpty ? nullptr : new Record());
    case variant::STRING: {
        if (region.toString().empty()) {
            return std::shared_ptr<Record>(nullIfEmpty ? 0 : new Record());
        }
        auto image = itsImPol->getImage();
        ThrowIf(
            ! image, "No attached image. Cannot use a string value for region"
        );
        return std::shared_ptr<Record>(
            new Record(
                CasacRegionManager::regionFromString(
                    image->coordinates(), region.toString(),
                    image->name(False), image->shape(), True
                )
            )
        );
    }
    case variant::RECORD: {
        return std::shared_ptr<Record>(
            nullIfEmpty && region.size() == 0
            ? 0
            : toRecord(
                std::shared_ptr<variant>(region.clone())->asRecord()
            )
        );
    }
    default:
        ThrowCc("Unsupported type for region " + region.typeString());
    }
}

} // casac namespace
