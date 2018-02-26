//# ImagePolarimetry.h: Polarimetric analysis of images
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: ImagePolarimetry.h 19817 2006-12-22 05:24:00Z gvandiep $

#ifndef IMAGES_IMAGEPOLARIMETRY_H
#define IMAGES_IMAGEPOLARIMETRY_H


//# Includes
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Block.h> 
#include <measures/Measures/Stokes.h>
#include <casa/BasicSL/Complex.h>
#include <images/Images/ImageInterface.h>
#include <scimath/Fitting/LinearFitSVD.h>

#include <imageanalysis/ImageTypedefs.h>

namespace casacore{

template <class T> class SubImage;
template <class T> class ImageExpr;
template <class T> class Quantum;
template <class T> class LatticeStatistics;
class CoordinateSystem;
class IPosition;
class LatticeExprNode;
class LCBox;
class LogIO;
}

namespace casa { //# NAMESPACE CASA - BEGIN

//# Forward Declarations
//

// <summary>
// Polarimetric analysis of images
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> <linkto class=casacore::ImageExpr>ImageExpr</linkto>
//   <li> <linkto class=casacore::ImageInterface>ImageInterface</linkto>
// </prerequisite>

// <etymology>
//  Polarimetric analysis of Images
// </etymology>

// <synopsis>
// This class provides polarimetric image analysis capability.
// It takes an image with a casacore::Stokes axis (some combination
// of IQUV is needed) as its input.
//
// Many functions return casacore::ImageExpr objects.  These are
// read-only images.
//
// Sometimes the standard deviation of the noise is needed.
// This is for debiasing polarized intensity images or working out
// error images.  By default it is worked out for you with a
// clipped mean algorithm.  However, you can provide sigma if you
// know it accurately.   It should be the standard deviation of the noise in
// the absence of signal.  You won't measure that very well from
// casacore::Stokes I if it is dynamic range limited.  Better to get it from 
// V or Q or U.  When this class needs the standard deviation of
// the noise, it will try and get it from V or Q and U and finally I.
//
// However, note that the functions sigmaStokes{I,Q,U,V} DO return the standard
// deviation of the noise for that specific casacore::Stokes type.
//
// The casacore::ImageExpr objects returned have the brightness units and ImageInfo
// set.  The MiscInfo (a permanent record) and logSink are not set.
// 
// </synopsis>
//
// <motivation>
// Basic image analysis capability
// </motivation>

// <todo asof="1999/11/01">
//   <li> plotting for function rotationMeasure 
//   <li> some assessment of the curvature of pa-l**2
// </todo>


class ImagePolarimetry 
{
public:

// casacore::Stokes types
   enum StokesTypes {I, Q, U, V}; 

// Constructor.  The input image must have a Stokes
// axis with some subset of I,Q,U, and V
   ImagePolarimetry (const casacore::ImageInterface<casacore::Float>& image);

// Copy constructor (reference semantics)
   ImagePolarimetry(const ImagePolarimetry& other);

// Destructor
   virtual ~ImagePolarimetry ();
   
// Assignment operator (reference semantics)
   ImagePolarimetry& operator=(const ImagePolarimetry& other);

// Summary.  Just invokes the casacore::ImageSummary list function
// to summarize the header of the construction image
   void summary(casacore::LogIO& os) const;

// Get the casacore::ImageInterface pointer of the construction image
// Don't delete it !
   SPCIIF imageInterface() const {return itsInImagePtr;};

// Get the casacore::CoordinateSystem of the construction image
   casacore::CoordinateSystem coordinates() const {return itsInImagePtr->coordinates();};

// Get the shape of the construction image
   casacore::IPosition shape() const {return itsInImagePtr->shape();};

// Is the construction image masked ?
   casacore::Bool isMasked() const {return itsInImagePtr->isMasked();};

// Get the shape and casacore::CoordinateSystem of an image for a single casacore::Stokes pixel
// Thus, if the construction image shape was [10,10,4,20] where
// axis 2 (shape 4) is the casacore::Stokes axis, this function
// would return [10,10,1,20]    Specify the type of casacore::Stokes pixel
// you want.  
   casacore::IPosition singleStokesShape(casacore::CoordinateSystem& cSys, casacore::Stokes::StokesTypes type) const;

// casacore::Complex linear polarization
   casacore::ImageExpr<casacore::Complex> complexLinearPolarization ();

// casacore::Complex fractional linear polarization
   casacore::ImageExpr<casacore::Complex> complexFractionalLinearPolarization ();

// Get the casacore::Stokes I image and the standard deviation of the
// I image.  This  is worked out by first clipping 
// outliers from the mean at the specified level.
// <group>
   casacore::ImageExpr<casacore::Float> stokesI() const;
   casacore::Float sigmaStokesI (casacore::Float clip=10.0);
// </group>

// Get the casacore::Stokes Q image and the standard deviation 
// of the Q image.  This  is worked out by first clipping 
// outliers from the mean at the specified level.
// <group>
   casacore::ImageExpr<casacore::Float> stokesQ() const;
   casacore::Float sigmaStokesQ (casacore::Float clip=10.0);
// </group>

// Get the casacore::Stokes U image and the standard deviation 
// of the U image.  This  is worked out by first clipping 
// outliers from the mean at the specified level.
// <group>
   casacore::ImageExpr<casacore::Float> stokesU() const;
   casacore::Float sigmaStokesU (casacore::Float clip=10.0);
// </group>

// Get the casacore::Stokes V image and the standard deviation 
// of the V image.  This  is worked out by first clipping 
// outliers from the mean at the specified level.
// <group>
   casacore::ImageExpr<casacore::Float> stokesV() const;
   casacore::Float sigmaStokesV (casacore::Float clip=10.0);
// </group>

// Get the specified casacore::Stokes image and the standard deviation 
// of the image.  This  is worked out by first clipping 
// outliers from the mean at the specified level.
// <group>
   casacore::ImageExpr<casacore::Float> stokes(ImagePolarimetry::StokesTypes index) const;
   casacore::Float sigmaStokes (ImagePolarimetry::StokesTypes index, casacore::Float clip=10.0);
// </group>

// Get the best estimate of the statistical noise. This gives you
// the standard deviation with outliers from the mean
// clipped first. The idea is to not be confused by source or dynamic range issues.
// Generally casacore::Stokes V is empty of sources (not always), then Q and U are generally
// less bright than I.  So this function first tries V, then Q and U 
// and lastly I to give you its noise estimate
   casacore::Float sigma (casacore::Float clip=10.0);

// Get the linearly polarized intensity image and its
// standard deviation.  If wish to debias the image, you
// can either provide <src>sigma</src> (the standard 
// deviation of the termal noise ) or if <src>sigma</src> is non-positive, 
// it will  be worked out for you with a clipped mean algorithm.
// <group>
   casacore::Float sigmaLinPolInt (casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

// Get the total polarized intensity (from whatever combination
// of Q, U, and V the construction image has) image and its error 
// (standard deviation).  If wish to debias the image, you
// can either provide <src>sigma</src> (the standard  deviation 
// of the thermal noise) or if <src>sigma</src> is 
// non-positive, it will be worked out for you with a 
// clipped mean algorithm.
// <group>
   casacore::ImageExpr<casacore::Float> totPolInt(casacore::Bool debias, casacore::Float clip=10.0, casacore::Float sigma=-1.0);
   casacore::Float sigmaTotPolInt (casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

// Get linearly polarized position angle (degrees or radians) image
// and error (standard deviation).   If you provide 
// <src>sigma</src> it is the  standard deviation of 
// the termal noise.  If <src>sigma</src> is non-positive, it will be 
// worked out for you with a  clipped mean algorithm.
// <group>
   casacore::ImageExpr<casacore::Float> linPolPosAng(casacore::Bool radians) const;
   casacore::ImageExpr<casacore::Float> sigmaLinPolPosAng (casacore::Bool radians, casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

// Get fractional linear polarization image 
// and error (standard deviation).   If wish to debias the image, you
// can either provide <src>sigma</src> (the standard 
// deviation of the termal noise) or if <src>sigma</src> is non-positive, 
// it will  be worked out for you with a clipped mean algorithm.
// <group>
   casacore::ImageExpr<casacore::Float> fracLinPol(casacore::Bool debias, casacore::Float clip=10.0, casacore::Float sigma=-1.0);
   casacore::ImageExpr<casacore::Float> sigmaFracLinPol (casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

// Get Fractional total polarization and error (standard deviation)
// <src>var</src> is the standard deviation  of the thermal noise.
// If <src>sigma</src> is non-positive, 
// it will  be worked out for you with a clipped mean algorithm.
// <group>
   casacore::ImageExpr<casacore::Float> fracTotPol(casacore::Bool debias, casacore::Float clip=10.0, casacore::Float sigma=-1.0);
   casacore::ImageExpr<casacore::Float> sigmaFracTotPol (casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

// Fourier Rotation Measure.  The output image is the complex polarization
// (Q + iU) with the spectral axis replaced by a RotationMeasure axis.
// The appropriate shape and casacore::CoordinateSystem must be obtained with
// function singleStokesShape (ask for type STokes::Plinear).  
// Howeverm this function will replace the casacore::SpectralCoordinate 
// by a casacore::LinearCoordinate describing  the Rotation Measure.  
// casacore::ImageInfo, and Units are copied to the output.  MiscInfo and
// history are not.  If the output has a mask,
// and the input is masked, the mask is copied.  If the output
// has a mask, it should already have been initialized to true
   void fourierRotationMeasure(casacore::ImageInterface<casacore::Complex>& pol,
                               casacore::Bool zeroZeroLag);

// This function is used in concert with the rotationMeasure function.
// It tells you what the shape of the output RM image should be, and
// gives you its CoordinateSystem.  Because the ImagePolarimetry 
// construction image may house the frequencies coordinate description
// in a Spectral, Tabular or Linear coordinate, you may explicitly
// specify which axis is the Spectral axis (spectralAxis).  By default,
// it tries to find the SpectralCoordinate.  If there is none, it will
// look for Tabular or Linear coordinates with a "FREQ" label.
// It returns to you the frequencyAxis (i.e. the one it is concluded
// houses the frequency spectrum) and the stokesAxis that it
// finds.
   casacore::IPosition rotationMeasureShape(casacore::CoordinateSystem& cSys,
                                  casacore::Int& frequencyAxis, casacore::Int& stokesAxis, 
                                  casacore::LogIO& os, casacore::Int spectralAxis=-1) const;

// This function is used in concert with the rotationMeasure function.
// It tells you what the shape of the output Position Angle image should be, and
// gives you its CoordinateSystem.  Because the ImagePolarimetry
// construction image may house the frequencies coordinate description
// in a Spectral, Tabular or Linear coordinate, you may explicitly
// specify which axis is the Spectral axis (spectralAxis).  By default,
// it tries to find the SpectralCoordinate.  If there is none, it will
// look for Tabular or Linear coordinates with a "FREQ" label.
   casacore::IPosition positionAngleShape(casacore::CoordinateSystem& cSys,
                                casacore::Int& frequencyAxis, casacore::Int& stokesAxis, 
                                casacore::LogIO& os, casacore::Int spectralAxis=-1) const;

// This function applies a traditional (i.e. non-Fourier) Rotation casacore::Measure 
// algorithm (Leahy et al, A&A, 156, 234) approach.     For the RM images
// you must get the shape and CoordinateSYstem from function
// rotationMeasureShape.  For the position angle images, use function
// singleStokesShape.  Null pointer casacore::ImageInterface objects are 
// not accessed so you can select which output images you want.
// Any output images not masked will be given a mask.
// The position angles are all in degrees. The RM images in rad/m/m.
// casacore::ImageInfo and Units, are copied to the output.  MiscInfo and history are not.
// You specify which axis houses the frequencies, the noise level of
// Q and U  if you  know it (by default it is worked out for you) 
// for error images, the value of the foreground RM if you know it
// (helps for unwinding ambiguity), the absolute maximum RM it should 
// solve for, and the maximum error in the position angle that should
// be allowed.  The state of the plotter should be set by the caller
// (e.g. character size, number of plots in x and y etc).
   void rotationMeasure(casacore::ImageInterface<casacore::Float>*& rmPtr,  casacore::ImageInterface<casacore::Float>*& rmErrPtr, 
                        casacore::ImageInterface<casacore::Float>*& pa0Ptr, casacore::ImageInterface<casacore::Float>*& pa0ErrPtr,
                        casacore::ImageInterface<casacore::Float>*& nTurns, casacore::ImageInterface<casacore::Float>*& rChiSqPtr,
                        casacore::Int spectralAxis,  casacore::Float rmMax, casacore::Float maxPaErr=1.0e30,
                        casacore::Float sigma=-1.0, casacore::Float rmFg=0.0, casacore::Bool showProgress=false);

// Depolarization ratio image and error.  Requires two images hence static
// functions.
// <group>
   static casacore::ImageExpr<casacore::Float> depolarizationRatio (const casacore::ImageInterface<casacore::Float>& im1,
                                                const casacore::ImageInterface<casacore::Float>& im2,
                                                casacore::Bool debias, casacore::Float clip=10.0, casacore::Float sigma=-1.0);

   static casacore::ImageExpr<casacore::Float> sigmaDepolarizationRatio (const casacore::ImageInterface<casacore::Float>& im1,
                                                     const casacore::ImageInterface<casacore::Float>& im2,
                                                     casacore::Bool debias, casacore::Float clip=10.0, casacore::Float sigma=-1.0);
// </group>

private:
   SPCIIF itsInImagePtr;
   // const casacore::ImageInterface<casacore::Float>* itsInImagePtr;
   casacore::LinearFitSVD<casacore::Float>* itsFitterPtr;
   casacore::Float itsOldClip;

// These blocks are always size 4, with IQUV in slots 0,1,2,3
// If your image is IV only, they still use slots 0 and 3

   casacore::PtrBlock<casacore::ImageInterface<casacore::Float>* > itsStokesPtr;
   casacore::PtrBlock<casacore::LatticeStatistics<casacore::Float>* > itsStokesStatsPtr;

   casacore::Matrix<casacore::Bool> _beamsEqMat;


// Delete all private pointers
   void cleanup();

// Copy data and mask
   void copyDataAndMask(casacore::ImageInterface<casacore::Float>& out,
                        casacore::ImageInterface<casacore::Float>& in) const;

// For traiditional RM approach, give output a mask if possible
   casacore::Bool dealWithMask (casacore::Lattice<casacore::Bool>*& pMask, casacore::ImageInterface<casacore::Float>*& pIm, casacore::LogIO& os, const casacore::String& type) const;

// Change the casacore::Stokes casacore::Coordinate for the given float image to be of the specified casacore::Stokes type
   void fiddleStokesCoordinate(casacore::ImageInterface<casacore::Float>& ie, casacore::Stokes::StokesTypes type) const;
   void fiddleStokesCoordinate(casacore::CoordinateSystem& cSys, casacore::Stokes::StokesTypes type) const;

// Change the casacore::Stokes casacore::Coordinate for the given complex image to be of the specified casacore::Stokes type
   void fiddleStokesCoordinate(casacore::ImageInterface<casacore::Complex>& ie, casacore::Stokes::StokesTypes type) const;

// Change the time coordinate to be rotation measure
   void fiddleTimeCoordinate(casacore::ImageInterface<casacore::Complex>& ie, const casacore::Quantum<casacore::Double>& f, casacore::Int coord) const;

// Find the central frequency from the given spectral coordinate
   casacore::Quantum<casacore::Double> findCentralFrequency(const casacore::Coordinate& coord, casacore::Int shape) const;

// Fit the spectrum of position angles to find the rotation measure via Leahy algorithm
   casacore::Bool findRotationMeasure (casacore::Float& rmFitted, casacore::Float& rmErrFitted,
                             casacore::Float& pa0Fitted, casacore::Float& pa0ErrFitted, casacore::Float& rChiSqFitted, 
                             casacore::Float& nTurns,
                             const casacore::Vector<casacore::uInt>& sortidx, const casacore::Vector<casacore::Float>& wsq, 
                             const casacore::Vector<casacore::Float>& pa, 
                             const casacore::Array<casacore::Bool>& paMask, 
                             const casacore::Array<casacore::Float>& paerr, 
                             casacore::Float rmfg, casacore::Float rmmax, casacore::Float paErrMax,
                             const casacore::String& posString);

// Find the casacore::Stokes in the construction image and assign pointers
   void findStokes();

// Find the spectral coordinate. 
   casacore::Int findSpectralCoordinate(const casacore::CoordinateSystem& cSys, casacore::LogIO& os, casacore::Bool fail) const;

// FInd frequency axis
   void findFrequencyAxis (casacore::Int& spectralCoord, casacore::Int& fAxis,
                           const casacore::CoordinateSystem& cSys, casacore::Int spectralAxis) const;
// So we have Q and U ?  Excpetion if not
   void hasQU () const;

// Make a LEN for the give types of polarized intensity
   casacore::LatticeExprNode makePolIntNode(casacore::LogIO& os, casacore::Bool debias, casacore::Float clip, casacore::Float sigma,
                                  casacore::Bool doLin, casacore::Bool doCirc);

// Make an IE for the specified Stokes
   casacore::ImageExpr<casacore::Float> makeStokesExpr(casacore::ImageInterface<casacore::Float>* imPtr,
                                   casacore::Stokes::StokesTypes type, const casacore::String& name) const;

// Make a casacore::SubImage from the construction image for the specified pixel
// along the specified pixel axis
   casacore::ImageInterface<casacore::Float>* makeSubImage (casacore::IPosition& blc, casacore::IPosition& trc,
                                        casacore::Int axis, casacore::Int pix) const;

// Least squares fit to find RM from position angles
   casacore::Bool rmLsqFit (casacore::Vector<casacore::Float>& pars, const casacore::Vector<casacore::Float>& wsq, 
                  const casacore::Vector<casacore::Float> pa, const casacore::Vector<casacore::Float>& paerr) const;

// Fit the spectrum of position angles to find the rotation measure via Leahy algorithm
// for primary (n>2) points
   casacore::Bool rmPrimaryFit (casacore::Float& nTurns, casacore::Float& rmFitted, casacore::Float& rmErrFitted,
                      casacore::Float& pa0Fitted, casacore::Float& pa0ErrFitted,
                      casacore::Float& rChiSqFitted, const casacore::Vector<casacore::Float>& wsq, 
                      const casacore::Vector<casacore::Float>& pa, const casacore::Vector<casacore::Float>& paerr, 
                      casacore::Float rmmax, const casacore::String& posString);

// Fit the spectrum of position angles to find the rotation measure via Leahy algorithm
// for supplementary (n==2) points
   casacore::Bool rmSupplementaryFit (casacore::Float& nTurns, casacore::Float& rmFitted, casacore::Float& rmErrFitted,
                            casacore::Float& pa0Fitted, casacore::Float& pa0ErrFitted,
                            casacore::Float& rChiSqFitted, const casacore::Vector<casacore::Float>& wsq, 
                            const casacore::Vector<casacore::Float>& pa, const casacore::Vector<casacore::Float>& paerr);

// Return I, Q, U or V for specified integer index (0-3)
   casacore::String stokesName (ImagePolarimetry::StokesTypes index) const;

// Return I, Q, U or V for specified integer index (0-3)
   casacore::Stokes::StokesTypes stokesType (ImagePolarimetry::StokesTypes index) const;

// Find the standard deviation for the casacore::Stokes image specified by the integer index
   casacore::Float sigma (ImagePolarimetry::StokesTypes index, casacore::Float clip);

// Subtract profile mean from image
   void subtractProfileMean (casacore::ImageInterface<casacore::Float>& im, casacore::uInt axis) const;

   void _createBeamsEqMat();

   casacore::Bool _checkBeams(
		   const casacore::Vector<StokesTypes>& stokes,
		   const casacore::Bool requireChannelEquality,
		   const casacore::Bool throws=true
   ) const;

   casacore::Bool _checkIQUBeams(
		   const casacore::Bool requireChannelEquality,
		   const casacore::Bool throws=true
   ) const;

   casacore::Bool _checkIVBeams(
		   const casacore::Bool requireChannelEquality,
		   const casacore::Bool throws=true
   ) const;

   casacore::Bool _checkQUBeams(
		   const casacore::Bool requireChannelEquality,
		   const casacore::Bool throws=true
   ) const;

   static void _checkBeams(
		   const ImagePolarimetry& im1, const ImagePolarimetry& im2,
		   const casacore::Vector<StokesTypes>& stokes
   );

   void _setInfo(casacore::ImageInterface<casacore::Complex>& im, const StokesTypes stokes) const;

   void _setInfo(casacore::ImageInterface<casacore::Float>& im, const StokesTypes stokes) const;

   void _setDoLinDoCirc(casacore::Bool& doLin, casacore::Bool& doCirc) const;

};


} //# NAMESPACE CASA - END

#endif
