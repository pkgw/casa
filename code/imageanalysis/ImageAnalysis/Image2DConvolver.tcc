//# Image2DConvolver.cc:  convolution of an image by given Array
//# Copyright (C) 1995,1996,1997,1998,1999,2000,2001,2002
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
//   
#include <imageanalysis/ImageAnalysis/Image2DConvolver.h>

#include <casa/aips.h>
#include <casa/Arrays/IPosition.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Exceptions/Error.h>
#include <components/ComponentModels/GaussianDeconvolver.h>
#include <components/ComponentModels/GaussianShape.h>
#include <components/ComponentModels/SkyComponentFactory.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <lattices/LatticeMath/Fit2D.h>
#include <scimath/Functionals/Gaussian2D.h>
#include <imageanalysis/ImageAnalysis/ImageConvolver.h>
#include <imageanalysis/ImageAnalysis/ImageMetaData.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageUtilities.h>
#include <images/Images/SubImage.h>
#include <casa/Logging/LogIO.h>
#include <scimath/Mathematics/Convolver.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/Unit.h>
#include <casa/Quanta/QLogical.h>
#include <casa/iostream.h>

#include <memory>

namespace casa {

template <class T> const casacore::String Image2DConvolver<T>::CLASS_NAME = "Image2DConvolver";

template <class T> Image2DConvolver<T>::Image2DConvolver(
    const SPCIIT image, const casacore::Record *const &region,
    const casacore::String& mask, const casacore::String& outname, const casacore::Bool overwrite
) : ImageTask<T>(image, "", region, "", "", "", mask, outname, overwrite),
    _type(casacore::VectorKernel::GAUSSIAN),  _scale(0), _major(), _minor(),
    _pa(), _axes(image->coordinates().directionAxesNumbers()), _targetres(false),
    _suppressWarnings(false) {
    this->_construct(true);
}

// TODO use GaussianBeams rather than casacore::Vector<casacore::Quantity>s, this method
// can probably be eliminated.
template <class T>
std::vector<casacore::Quantity> Image2DConvolver<T>::_getConvolvingBeamForTargetResolution(
    const std::vector<casacore::Quantity>& targetBeamParms,
    const casacore::GaussianBeam& inputBeam
) const {
    casacore::GaussianBeam convolvingBeam;

    casacore::GaussianBeam targetBeam(
        targetBeamParms[0], targetBeamParms[1],
        targetBeamParms[2]
    );
    try {
        if(GaussianDeconvolver::deconvolve(convolvingBeam, targetBeam, inputBeam)) {
            // point source, or convolvingBeam nonsensical
            throw casacore::AipsError();
        }
    }
    catch (const casacore::AipsError& x) {
        ostringstream os;
        os << "Unable to reach target resolution of "
            << targetBeam << " Input image beam "
            << inputBeam << " is (nearly) identical "
            << "to or larger than the output beam size";
        ThrowCc(os.str());
    }
    std::vector<casacore::Quantity> kernelParms {
        convolvingBeam.getMajor(),
        convolvingBeam.getMinor(),
        convolvingBeam.getPA(true)
    };
    return kernelParms;
}

template <class T> void Image2DConvolver<T>::setAxes(
    const std::pair<casacore::uInt, casacore::uInt>& axes
) {
    casacore::uInt ndim = this->_getImage()->ndim();
    ThrowIf(axes.first == axes.second, "Axes must be different");
    ThrowIf(
        axes.first >= ndim || axes.second >= ndim,
        "Axis value must be less than number of axes in image"
    );
    if (_axes.size() != 2) {
        _axes.resize(2, false);
    }
    _axes[0] = axes.first;
    _axes[1] = axes.second;
}

template <class T> void Image2DConvolver<T>::setKernel(
    const casacore::String& type, const casacore::Quantity& major, const casacore::Quantity& minor,
    const casacore::Quantity& pa
) {
    ThrowIf (major < minor, "Major axis is less than minor axis");
    _type = casacore::VectorKernel::toKernelType(type);
    _major = major;
    _minor = minor;
    _pa = pa;
}

template <class T> SPIIT Image2DConvolver<T>::convolve() {
    ThrowIf(
        _axes.nelements() != 2,
        "You must give two pixel axes to convolve"
    );
    casacore::Vector<casacore::Double> inc = this->_getImage()->coordinates().increment();
    casacore::Vector<casacore::String> units = this->_getImage()->coordinates().worldAxisUnits();
    ThrowIf(
        ! near (
            casacore::Quantity(fabs(inc[_axes[0]]), units[_axes[0]]),
            casacore::Quantity(fabs(inc[_axes[1]]), units[_axes[1]])
        ),
        "Pixels must be square, please regrid your image so that they are"
    );
    auto subImage = SubImageFactory<casacore::Float>::createImage(
        *this->_getImage(), "", *this->_getRegion(), this->_getMask(),
        this->_getDropDegen(), false, false, this->_getStretch()
    );
    const casacore::Int nDim = subImage->ndim();
    ThrowIf(
        _axes(0) < 0 || _axes(0) >= nDim
        || _axes(1) < 0 || _axes(1) >= nDim,
        "The pixel axes " + casacore::String::toString(_axes) + " are illegal"
    );
    ThrowIf(
        nDim < 2,
        "The image axes must have at least 2 pixel axes"
    );
    SPIIT outImage(new casacore::TempImage<casacore::Float> (subImage->shape(), subImage->coordinates()));
    Image2DConvolver<casacore::Float>::_convolve(
        *this->_getLog(), outImage, *subImage, _type
    );
    return this->_prepareOutputImage(*outImage);
}

template <class T> void Image2DConvolver<T>::_convolve(
    casacore::LogIO& os, SPIIT imageOut,
    const casacore::ImageInterface<T>& imageIn, casacore::VectorKernel::KernelTypes kernelType
) const {
    const casacore::IPosition& inShape = imageIn.shape();
    const casacore::IPosition& outShape = imageOut->shape();
    ThrowIf(
        ! inShape.isEqual(outShape),
        "Input and output images must have the same shape"
    );
    // Generate Kernel casacore::Array (height unity)
    ThrowIf(
        _targetres && kernelType != casacore::VectorKernel::GAUSSIAN,
        "targetres can only be true for a Gaussian convolving kernel"
    );
    casacore::Array<T> kernel;
    // initialize to avoid compiler warning, kernelVolume will always be set to something
    // reasonable below before it is used.
    T kernelVolume = -1;
    std::vector<casacore::Quantity> originalParms{_major, _minor, _pa};
    if (! _targetres) {
        kernelVolume = _makeKernel(
            kernel, kernelType, originalParms, imageIn
        );
    }
    const casacore::CoordinateSystem& cSys = imageIn.coordinates();
    if (_major.getUnit().startsWith("pix")) {
        auto inc = cSys.increment()[_axes[0]];
        auto unit = cSys.worldAxisUnits()[_axes[0]];
        originalParms[0] = _major.getValue() * casacore::Quantity(abs(inc), unit);
    }
    if (_minor.getUnit().startsWith("pix")) {
        auto inc = cSys.increment()[_axes[1]];
        auto unit = cSys.worldAxisUnits()[_axes[1]];
        originalParms[1] = _minor.getValue() * casacore::Quantity(abs(inc), unit);
    }

    std::vector<casacore::Quantity> kernelParms = originalParms;

    // Figure out output image restoring beam (if any), output units and scale
    // factor for convolution kernel array

    casacore::GaussianBeam beamOut;
    const auto& imageInfo = imageIn.imageInfo();
    const auto& brightnessUnit = imageIn.units();
    casacore::String brightnessUnitOut;
    casacore::ImageInfo iiOut = imageOut->imageInfo();
    casacore::Bool logFactors = false;
    casacore::Double factor1 = -1;
    double pixelArea = 0;
    auto autoScale = _scale <= 0;
    if (autoScale) {
        auto bunitUp = brightnessUnit.getName();
        bunitUp.upcase();
        logFactors = bunitUp.contains("/BEAM");
        if (logFactors) {
            pixelArea = cSys.directionCoordinate().getPixelArea().getValue("arcsec*arcsec");
            if (! _targetres) {
                casacore::GaussianBeam kernelBeam(kernelParms);
                factor1 = pixelArea/kernelBeam.getArea("arcsec*arcsec");
            }
        }
    }
    if (imageInfo.hasMultipleBeams()) {
        _doMultipleBeams(
            os, iiOut, kernelVolume, imageOut, brightnessUnitOut,
            beamOut, factor1, imageIn, originalParms, kernelParms,
            kernel, kernelType, logFactors, pixelArea
        );
    }
    else {
        _doSingleBeam(
            os, iiOut, kernelVolume, kernelParms, kernel,
            brightnessUnitOut, beamOut, imageOut, imageIn,
            originalParms, kernelType, logFactors, factor1, 
            pixelArea
        );
        /*
        casacore::GaussianBeam inputBeam = imageInfo.restoringBeam();
        if (_targetres) {
            kernelParms = _getConvolvingBeamForTargetResolution(
                originalParms, inputBeam
            );
            os << casacore::LogIO::NORMAL << "Convolving image that has a beam of "
                << inputBeam << " with a Gaussian of "
                << casacore::GaussianBeam(kernelParms) << " to reach a target resolution of "
                << casacore::GaussianBeam(originalParms) << casacore::LogIO::POST;
            kernelVolume = _makeKernel(
                kernel, kernelType, kernelParms, imageIn
            );
        }
        T scaleFactor = _dealWithRestoringBeam(
            os, brightnessUnitOut, beamOut, kernel, kernelVolume,
            kernelType, kernelParms, cSys, inputBeam,
            brightnessUnit, true
        );
        os << casacore::LogIO::NORMAL << "Scaling pixel values by ";
        if (logFactors) {
            if (_targetres) {
                casacore::GaussianBeam kernelBeam(kernelParms);
                factor1 = pixelArea/kernelBeam.getArea("arcsec*arcsec");
            }
            casacore::Double factor2 = beamOut.getArea("arcsec*arcsec")/inputBeam.getArea("arcsec*arcsec");
            os << "inverse of area of convolution kernel in pixels (" << factor1
                << ") times the ratio of the beam areas (" << factor2 << ") = ";
        }
        os << scaleFactor << casacore::LogIO::POST;
        if (_targetres && near(beamOut.getMajor(), beamOut.getMinor(), 1e-7)) {
            // circular beam should have same PA as given by user if
            // targetres
            beamOut.setPA(originalParms[2]);
        }
        // Convolve.  We have already scaled the convolution kernel (with some
        // trickery cleverer than what ImageConvolver can do) so no more scaling
        ImageConvolver<T> aic;
        aic.convolve(
            os, *imageOut, imageIn, scaleFactor*kernel, ImageConvolver<T>::NONE,
            1.0, true
        );
        // Overwrite some bits and pieces in the output image to do with the
        // restoring beam  and image units
        casacore::Bool holdsOneSkyAxis;
        casacore::Bool hasSky = casacore::CoordinateUtil::holdsSky (holdsOneSkyAxis, cSys, _axes.asVector());
        if (hasSky && ! beamOut.isNull()) {
            iiOut.setRestoringBeam(beamOut);
        }
        else {
            // If one of the axes is in the sky plane, we must
            // delete the restoring beam as it is no longer meaningful
            if (holdsOneSkyAxis) {
                os << casacore::LogIO::WARN << "Because you convolved just one of the sky axes" << endl;
                os << "The output image does not have a valid spatial restoring beam" << casacore::LogIO::POST;
                iiOut.removeRestoringBeam();
            }
        }
        */
    }
    imageOut->setUnits(brightnessUnitOut);
    imageOut->setImageInfo(iiOut);
}

template <class T> void Image2DConvolver<T>::_doSingleBeam(
    LogIO& os, ImageInfo& iiOut, T& kernelVolume, vector<Quantity>& kernelParms, Array<T>& kernel,
    String& brightnessUnitOut, GaussianBeam& beamOut, SPIIT imageOut,
    const ImageInterface<T>& imageIn, const vector<Quantity>& originalParms,
    VectorKernel::KernelTypes kernelType, Bool logFactors, Double factor1, Double pixelArea
) const {
    casacore::GaussianBeam inputBeam = imageIn.imageInfo().restoringBeam();
    if (_targetres) {
        kernelParms = _getConvolvingBeamForTargetResolution(
            originalParms, inputBeam
        );
        os << casacore::LogIO::NORMAL << "Convolving image that has a beam of "
            << inputBeam << " with a Gaussian of "
            << casacore::GaussianBeam(kernelParms) << " to reach a target resolution of "
            << casacore::GaussianBeam(originalParms) << casacore::LogIO::POST;
        kernelVolume = _makeKernel(
            kernel, kernelType, kernelParms, imageIn
        );
    }
    const CoordinateSystem& cSys = imageIn.coordinates();
    T scaleFactor = _dealWithRestoringBeam(
        os, brightnessUnitOut, beamOut, kernel, kernelVolume,
        kernelType, kernelParms, cSys, inputBeam,
        imageIn.units(), true
    );
    os << casacore::LogIO::NORMAL << "Scaling pixel values by ";
    if (logFactors) {
        if (_targetres) {
            casacore::GaussianBeam kernelBeam(kernelParms);
            factor1 = pixelArea/kernelBeam.getArea("arcsec*arcsec");
        }
        casacore::Double factor2 = beamOut.getArea("arcsec*arcsec")/inputBeam.getArea("arcsec*arcsec");
        os << "inverse of area of convolution kernel in pixels (" << factor1
            << ") times the ratio of the beam areas (" << factor2 << ") = ";
    }
    os << scaleFactor << casacore::LogIO::POST;
    if (_targetres && near(beamOut.getMajor(), beamOut.getMinor(), 1e-7)) {
        // circular beam should have same PA as given by user if
        // targetres
        beamOut.setPA(originalParms[2]);
    }
    // Convolve.  We have already scaled the convolution kernel (with some
    // trickery cleverer than what ImageConvolver can do) so no more scaling
    ImageConvolver<T> aic;
    aic.convolve(
        os, *imageOut, imageIn, scaleFactor*kernel, ImageConvolver<T>::NONE,
        1.0, true
    );
    // Overwrite some bits and pieces in the output image to do with the
    // restoring beam  and image units
    casacore::Bool holdsOneSkyAxis;
    casacore::Bool hasSky = casacore::CoordinateUtil::holdsSky (holdsOneSkyAxis, cSys, _axes.asVector());
    if (hasSky && ! beamOut.isNull()) {
        iiOut.setRestoringBeam(beamOut);
    }
    else {
        // If one of the axes is in the sky plane, we must
        // delete the restoring beam as it is no longer meaningful
        if (holdsOneSkyAxis) {
            os << casacore::LogIO::WARN << "Because you convolved just one of the sky axes" << endl;
            os << "The output image does not have a valid spatial restoring beam" << casacore::LogIO::POST;
            iiOut.removeRestoringBeam();
        }
    }
}

template <class T> void Image2DConvolver<T>::_doMultipleBeams(
    LogIO& os, ImageInfo& iiOut, T& kernelVolume, SPIIT imageOut, String& brightnessUnitOut,
    GaussianBeam& beamOut, Double factor1, const ImageInterface<T>& imageIn, const vector<Quantity>& originalParms,
    vector<Quantity>& kernelParms, Array<T>& kernel,
    VectorKernel::KernelTypes kernelType, Bool logFactors, Double pixelArea
) const {
    ImageMetaData md(imageOut);
    casacore::uInt nChan = md.nChannels();
    casacore::uInt nPol = md.nStokes();
    // initialize all beams to be null
    iiOut.setAllBeams(nChan, nPol, casacore::GaussianBeam());
    const CoordinateSystem& cSys = imageIn.coordinates();
    casacore::Int specAxis = cSys.spectralAxisNumber();
    casacore::Int polAxis = cSys.polarizationAxisNumber();
    casacore::IPosition start(imageIn.ndim(), 0);
    casacore::IPosition end = imageIn.shape();
    if (nChan > 0) {
        end[specAxis] = 1;
    }
    if (nPol > 0) {
        end[polAxis] = 1;
    }
    casacore::Int channel = -1;
    casacore::Int polarization = -1;
    if (_targetres) {
        iiOut.removeRestoringBeam();
        iiOut.setRestoringBeam(casacore::GaussianBeam(kernelParms));
    }
    casacore::uInt count = (nChan > 0 && nPol > 0)
        ? nChan * nPol
        : nChan > 0
          ? nChan
          : nPol;
    for (casacore::uInt i=0; i<count; i++) {
        if (nChan > 0) {
            channel = i % nChan;
            start[specAxis] = channel;
        }
        if (nPol > 0) {
            polarization = nChan > 1
                ? (i - channel) % nChan
                : i;
            start[polAxis] = polarization;
        }
        casacore::Slicer slice(start, end);
        casacore::SubImage<T> subImage(imageIn, slice);
        casacore::CoordinateSystem subCsys = subImage.coordinates();
        if (subCsys.hasSpectralAxis()) {
            casacore::Vector<casacore::Double> subRefPix = subCsys.referencePixel();
            subRefPix[specAxis] = 0;
            subCsys.setReferencePixel(subRefPix);
        }
        auto inputBeam = imageIn.imageInfo().restoringBeam(channel, polarization);
        casacore::Bool doConvolve = true;
        if (_targetres) {
            os << casacore::LogIO::NORMAL;
            if (channel >= 0) {
                os << "Channel " << channel << " of " << nChan;
                if (polarization >= 0) {
                    os << ", ";
                }
            }
            if (polarization >= 0) {
                os << "Polarization " << polarization << " of " << nPol;
            }
            os << " ";
            if (near(inputBeam, casacore::GaussianBeam(originalParms), 1e-5, casacore::Quantity(1e-2, "arcsec"))) {
                doConvolve = false;
                os << casacore::LogIO::NORMAL << " Input beam is already near target resolution so this "
                    << "plane will not be convolved" << casacore::LogIO::POST;
            }
            else {
                kernelParms = _getConvolvingBeamForTargetResolution(
                    originalParms, inputBeam
                );
                kernelVolume = _makeKernel(
                    kernel, kernelType, kernelParms, imageIn
                );
                os << ": Convolving image which has a beam of " << inputBeam
                    << " with a Gaussian of "
                    << casacore::GaussianBeam(kernelParms) << " to reach a target resolution of "
                    << casacore::GaussianBeam(originalParms) << casacore::LogIO::POST;
            }
        }
        casacore::TempImage<casacore::Float> subImageOut(
            subImage.shape(), subImage.coordinates()
        );
        if (doConvolve) {
            T scaleFactor = _dealWithRestoringBeam(
                os, brightnessUnitOut, beamOut, kernel, kernelVolume,
                kernelType, kernelParms, subCsys, inputBeam,
                imageIn.units(), i == 0
            );
            {
                os << casacore::LogIO::NORMAL << "Scaling pixel values by ";
                if (logFactors) {
                    if (_targetres) {
                        casacore::GaussianBeam kernelBeam(kernelParms);
                        factor1 = pixelArea/kernelBeam.getArea("arcsec*arcsec");
                    }
                    casacore::Double factor2 = beamOut.getArea("arcsec*arcsec")/inputBeam.getArea("arcsec*arcsec");
                    os << "inverse of area of convolution kernel in pixels (" << factor1
                        << ") times the ratio of the beam areas (" << factor2 << ") = ";
                }
                os << scaleFactor
                    << " for ";
                if (channel >= 0) {
                    os << "channel number " << channel;
                    if (polarization >= 0) {
                        os << " and ";
                    }
                }
                if (polarization >= 0) {
                    os << "polarization number " << polarization;
                }
                os << casacore::LogIO::POST;
            }
            if (_targetres && near(beamOut.getMajor(), beamOut.getMinor(), 1e-7)) {
                // circular beam should have same PA as given by user if
                // targetres
                beamOut.setPA(originalParms[2]);
            }
            ImageConvolver<T> aic;
            aic.convolve(
                os, subImageOut, subImage, scaleFactor*kernel,
                ImageConvolver<T>::NONE, 1.0, true
            );
        }
        else {
            brightnessUnitOut = imageIn.units().getName();
            beamOut = inputBeam;
            subImageOut.put(subImage.get());
        }
        {
            casacore::Bool doMask = imageOut->isMasked() && imageOut->hasPixelMask();
            casacore::Lattice<casacore::Bool>* pMaskOut = 0;
            if (doMask) {
                pMaskOut = &imageOut->pixelMask();
                if (! pMaskOut->isWritable()) {
                    doMask = false;
                }
            }
            casacore::IPosition cursorShape = subImageOut.niceCursorShape();
            casacore::IPosition outPos = start;
            casacore::LatticeStepper stepper(
                subImageOut.shape(), cursorShape, casacore::LatticeStepper::RESIZE
            );
            casacore::RO_MaskedLatticeIterator<T> iter(subImageOut, stepper);
            for (iter.reset(); !iter.atEnd(); iter++) {
                casacore::IPosition cursorShape = iter.cursorShape();
                imageOut->putSlice(iter.cursor(), outPos);
                if (doMask) {
                    pMaskOut->putSlice(iter.getMask(), outPos);
                }
                outPos = outPos + cursorShape;
            }
        }
        if (! _targetres) {
            iiOut.setBeam(
                channel, polarization, beamOut
            );
        }
    }
}

template <class T> T Image2DConvolver<T>::_makeKernel(
    casacore::Array<T>& kernelArray,
    casacore::VectorKernel::KernelTypes kernelType,
    const std::vector<casacore::Quantity>& parameters,
    const casacore::ImageInterface<T>& imageIn
) const {

// Check number of parameters

   _checkKernelParameters(kernelType, parameters);

// Convert kernel widths to pixels from world.  Demands major and minor
// both in pixels or both in world, else exception

   casacore::Vector<casacore::Double> dParameters;
   const casacore::CoordinateSystem cSys = imageIn.coordinates();

// Use the reference value for the shape conversion direction

   casacore::Vector<casacore::Quantity> wParameters(5);
   for (casacore::uInt i=0; i<3; i++) {
      wParameters(i+2) = parameters[i];
   }
//
   const casacore::Vector<casacore::Double> refVal = cSys.referenceValue();
   const casacore::Vector<casacore::String> units = cSys.worldAxisUnits();
   casacore::Int wAxis = cSys.pixelAxisToWorldAxis(_axes(0));
   wParameters(0) = casacore::Quantity(refVal(wAxis), units(wAxis));
   wAxis = cSys.pixelAxisToWorldAxis(_axes(1));
   wParameters(1) = casacore::Quantity(refVal(wAxis), units(wAxis));
   SkyComponentFactory::worldWidthsToPixel (dParameters, wParameters, cSys, _axes, false);

// Create n-Dim kernel array shape

   casacore::IPosition kernelShape = _shapeOfKernel (kernelType, dParameters, imageIn.ndim());

// Create kernel array. We will fill the n-Dim array (shape non-unity
// only for pixelAxes) through its 2D casacore::Matrix incarnation. Aren't we clever.
   kernelArray = 0;
   kernelArray.resize(kernelShape);
   casacore::Array<T> kernelArray2 = kernelArray.nonDegenerate (_axes);
   casacore::Matrix<T> kernelMatrix = static_cast<casacore::Matrix<T> >(kernelArray2);

// Fill kernel casacore::Matrix with functional (height unity)

   return _fillKernel (kernelMatrix, kernelType, kernelShape, dParameters);
}

template <class T> T Image2DConvolver<T>::_dealWithRestoringBeam(
    casacore::LogIO& os, casacore::String& brightnessUnitOut,
    casacore::GaussianBeam& beamOut, const casacore::Array<T>& kernelArray,
    const T kernelVolume, const casacore::VectorKernel::KernelTypes,
    const casacore::Vector<casacore::Quantity>& parameters,
    const casacore::CoordinateSystem& cSys,
    const casacore::GaussianBeam& beamIn, const casacore::Unit& brightnessUnitIn,
    casacore::Bool emitMessage
) const {
    os << casacore::LogOrigin(CLASS_NAME, __func__);
    // Find out if convolution axes hold the sky.  Scaling from
    // Jy/beam and Jy/pixel only really makes sense if this is true
    casacore::Bool holdsOneSkyAxis;
    casacore::Bool hasSky = casacore::CoordinateUtil::holdsSky (holdsOneSkyAxis, cSys, _axes.asVector());
    if (hasSky) {
        const casacore::DirectionCoordinate dc = cSys.directionCoordinate();
        auto inc = dc.increment();
        auto unit = dc.worldAxisUnits();
        casacore::Quantity x(inc[0], unit[0]);
        casacore::Quantity y(inc[1], unit[1]);
        auto diag = sqrt(x*x + y*y);
        auto minAx = parameters[1];
        if (minAx.getUnit().startsWith("pix")) {
            minAx.setValue(minAx.getValue()*x.getValue());
            minAx.setUnit(x.getUnit());
        }
        if (minAx < diag) {
            diag.convert(minAx.getFullUnit());
            if (! _suppressWarnings) {
                os << casacore::LogIO::WARN << "Convolving kernel has minor axis "
                    << minAx << " which is less than the pixel diagonal "
                    << "length of " << diag << ". Thus, the kernel is poorly sampled, "
                    << "and so the output of this application may not be what you expect. "
                    << "You should consider increasing the kernel size or regridding "
                    << "the image to a smaller pixel size" << casacore::LogIO::POST;
            }
        }
        else if (beamIn.getMinor() < diag && beamIn != casacore::GaussianBeam::NULL_BEAM) {
            diag.convert(beamIn.getMinor().getFullUnit());
            if (! _suppressWarnings) {
                os << casacore::LogIO::WARN << "Input beam has minor axis "
                    << beamIn.getMinor() << " which is less than the pixel diagonal "
                    << "length of " << diag << ". Thus, the beam is poorly sampled, "
                    << "and so the output of this application may not be what you expect. "
                    << "You should consider regridding "
                    << "the image to a smaller pixel size." << casacore::LogIO::POST;
            }
        }
    }
    if (emitMessage) {
        os << "You are " << (hasSky ? "" : " not ") << "convolving the sky" << casacore::LogIO::POST;
    }
    beamOut = casacore::GaussianBeam();
    auto bUnitIn = upcase(brightnessUnitIn.getName());
    const auto& refPix = cSys.referencePixel();
    T scaleFactor = 1;
    brightnessUnitOut = brightnessUnitIn.getName();
    auto autoScale = _scale <= 0;
    if (hasSky && bUnitIn.contains("/PIXEL")) {
        // Easy case.  Peak of convolution kernel must be unity
        // and output units are Jy/beam.  All other cases require
        // numerical convolution of beams
        brightnessUnitOut = "Jy/beam";
        // Exception already generated if only one of major and minor in pixel units
        auto majAx = parameters(0);
        auto minAx = parameters(1);
        if (majAx.getFullUnit().getName() == "pix") {
            casacore::Vector<casacore::Double> pixelParameters(5);
            pixelParameters(0) = refPix(_axes(0));
            pixelParameters(1) = refPix(_axes(1));
            pixelParameters(2) = parameters(0).getValue();
            pixelParameters(3) = parameters(1).getValue();
            pixelParameters(4) = parameters(2).getValue(casacore::Unit("rad"));
            casacore::GaussianBeam worldParameters;
            SkyComponentFactory::pixelWidthsToWorld(
                worldParameters, pixelParameters,
                cSys, _axes, false
            );
            majAx = worldParameters.getMajor();
            minAx = worldParameters.getMinor();
        }
        beamOut = casacore::GaussianBeam(majAx, minAx, parameters(2));
        // casacore::Input p.a. is positive N->E
        if (! autoScale) {
            scaleFactor = static_cast<T>(_scale);
            os << casacore::LogIO::WARN << "Autoscaling is recommended for Jy/pixel convolution"
                << casacore::LogIO::POST;
        }
    }
    else {
        // Is there an input restoring beam and are we convolving the sky to which it
        // pertains ?  If not, all we can do is use user scaling or normalize the convolution
        // kernel to unit volume.  There is no point to convolving the input beam either as it pertains
        // only to the sky
        if (hasSky && ! beamIn.isNull()) {
            // Convert restoring beam parameters to pixels.  Output pa is pos +x -> +y in pixel frame.
            casacore::Vector<casacore::Quantity> wParameters(5);
            const casacore::Vector<casacore::Double> refVal = cSys.referenceValue();
            const casacore::Vector<casacore::String> units = cSys.worldAxisUnits();
            auto wAxis = cSys.pixelAxisToWorldAxis(_axes(0));
            wParameters(0) = casacore::Quantity(refVal(wAxis), units(wAxis));
            wAxis = cSys.pixelAxisToWorldAxis(_axes(1));
            wParameters(1) = casacore::Quantity(refVal(wAxis), units(wAxis));
            wParameters(2) = beamIn.getMajor();
            wParameters(3) = beamIn.getMinor();
            wParameters(4) = beamIn.getPA(true);
            casacore::Vector<casacore::Double> dParameters;
            SkyComponentFactory::worldWidthsToPixel(
                dParameters, wParameters, cSys, _axes, false
            );
            // Create 2-D beam array shape
            // casacore::IPosition dummyAxes(2, 0, 1);
            auto beamShape = _shapeOfKernel(
                casacore::VectorKernel::GAUSSIAN,
                dParameters, 2 /*, dummyAxes */
            );

            // Create beam casacore::Matrix and fill with height unity
   
            casacore::Matrix<T> beamMatrixIn(beamShape(0), beamShape(1));
            _fillKernel(
                beamMatrixIn, casacore::VectorKernel::GAUSSIAN, beamShape,
                /*dummyAxes,*/ dParameters
            );

            auto shape = beamMatrixIn.shape();

            // Get 2-D version of convolution kenrel
            auto kernelArray2 = kernelArray.nonDegenerate(_axes);
            auto kernelMatrix = static_cast<casacore::Matrix<T> >(kernelArray2);
            // Convolve input restoring beam array by convolution kernel array
            casacore::Matrix<T> beamMatrixOut;

            casacore::Convolver<T> conv(beamMatrixIn, kernelMatrix.shape());
            conv.linearConv(beamMatrixOut, kernelMatrix);

            // Scale kernel
            T maxValOut = max(beamMatrixOut);

            scaleFactor = autoScale ? 1/maxValOut : (T)_scale;
            // Fit output beam matrix with a Gaussian, for better or worse
            // casacore::Fit2D is not templated.  So all our templating is useless
            // other than for casacore::Float until I template Fit2D
            casacore::Fit2D fitter(os);
            const casacore::uInt n = beamMatrixOut.shape()(0);
            auto bParameters = fitter.estimate(casacore::Fit2D::GAUSSIAN, beamMatrixOut);
            casacore::Vector<casacore::Bool> bParameterMask(bParameters.nelements(), true);
            bParameters(1) = (n-1)/2;          // x centre
            bParameters(2) = bParameters(1);    // y centre
            // Set range so we don't include too many pixels in fit which will make it very slow
            fitter.addModel (casacore::Fit2D::GAUSSIAN, bParameters, bParameterMask);
            casacore::Array<casacore::Float> sigma;
            fitter.setIncludeRange(maxValOut/10.0, maxValOut+0.1);
            auto error = fitter.fit(beamMatrixOut, sigma);
            ThrowIf(
                error == casacore::Fit2D::NOCONVERGE || error == casacore::Fit2D::FAILED
                || error == casacore::Fit2D::NOGOOD,
                "Failed to fit the output beam"
            );
            auto bSolution = fitter.availableSolution();
            // Convert to world units.
            casacore::Vector<casacore::Double> pixelParameters(5);
            pixelParameters(0) = refPix(_axes(0));
            pixelParameters(1) = refPix(_axes(1));
            pixelParameters(2) = bSolution(3);
            pixelParameters(3) = bSolution(4);
            pixelParameters(4) = bSolution(5);
            SkyComponentFactory::pixelWidthsToWorld(
                beamOut, pixelParameters, cSys, _axes, false
            );
            if (! brightnessUnitIn.getName().contains(casacore::Regex(Regex::makeCaseInsensitive("beam")))) {
                scaleFactor *= beamIn.getArea("arcsec2")/beamOut.getArea("arcsec2");
            }
        }
        else {
            if (autoScale) {
                // Conserve flux is the best we can do
                scaleFactor = 1/kernelVolume;
            }
            else {
                scaleFactor = (T)_scale;
            }
        }
    }
    // Put beam position angle into range +/- 180 in case it has eluded us so far
    if (! beamOut.isNull()) {
        casacore::MVAngle pa(beamOut.getPA(true).getValue(casacore::Unit("rad")));
        pa();
        beamOut = casacore::GaussianBeam(
            beamOut.getMajor(), beamOut.getMinor(),
            casacore::Quantity(pa.degree(), casacore::Unit("deg"))
        );
    }
    return scaleFactor;
}

template <class T> void Image2DConvolver<T>::_checkKernelParameters(
    casacore::VectorKernel::KernelTypes kernelType,
    const casacore::Vector<casacore::Quantity >& parameters
) const {
    if (kernelType==casacore::VectorKernel::BOXCAR) {
        ThrowCc("Boxcar kernel not yet implemented");
        ThrowIf(
            parameters.nelements() != 3,
            "Boxcar kernels require 3 parameters"
        );
    }
    else if (kernelType==casacore::VectorKernel::GAUSSIAN) {
        ThrowIf(
            parameters.nelements() != 3,
            "Gaussian kernels require exactly 3 parameters"
        );
    }
    else {
        ThrowCc(
            "The kernel type " + casacore::VectorKernel::fromKernelType(kernelType) + " is not supported"
        );
    }
}

template <class T> casacore::IPosition Image2DConvolver<T>::_shapeOfKernel(
    const casacore::VectorKernel::KernelTypes kernelType,
    const casacore::Vector<casacore::Double>& parameters,
    const casacore::uInt ndim
) const {
//
// Work out how big the array holding the kernel should be.
// Simplest algorithm possible. Shape is presently square.
//

// Find 2D shape

   casacore::uInt n;
   if (kernelType==casacore::VectorKernel::GAUSSIAN) {
      casacore::uInt n1 = _sizeOfGaussian (parameters(0), 5.0);
      casacore::uInt n2 = _sizeOfGaussian (parameters(1), 5.0);
      n = max(n1,n2);
      if (n%2==0) n++;                                     // Make shape odd so centres well
   } else if (kernelType==casacore::VectorKernel::BOXCAR) {
      n = 2 * casacore::Int(max(parameters(0), parameters(1))+0.5);
      if (n%2==0) n++;                                     // Make shape odd so centres well
   } else {
     throw(casacore::AipsError("Unrecognized kernel type"));        // Earlier checking should prevent this
   }

// Now find the shape for the image and slot the 2D shape in
// in the correct axis locations

   casacore::IPosition shape(ndim,1);
   shape(_axes(0)) = n;
   shape(_axes(1)) = n;
   return shape;
}
   
template <class T>
uInt Image2DConvolver<T>::_sizeOfGaussian(
    const casacore::Double width, const casacore::Double nSigma
) const {
// +/- 5sigma is a volume error of less than 6e-5%

   casacore::Double sigma = width / sqrt(casacore::Double(8.0) * C::ln2);
   return  (casacore::Int(nSigma*sigma + 0.5) + 1) * 2;
}


template <class T> T Image2DConvolver<T>::_fillKernel(
    casacore::Matrix<T>& kernelMatrix,
    casacore::VectorKernel::KernelTypes kernelType,
    const casacore::IPosition& kernelShape,
    const casacore::Vector<casacore::Double>& parameters
) const {

// Centre functional in array (shape is odd)
// Need to think about these T castes for casacore::Complex images

   T xCentre = static_cast<T>((kernelShape(_axes(0)) - 1) / 2.0);
   T yCentre = static_cast<T>((kernelShape(_axes(1)) - 1) / 2.0);
   T height = static_cast<T>(1.0);

// Create functional.  We only have gaussian2d functionals
// at this point.  Later the filling code can be moved out
// of the if statement

   T maxValKernel;
   T volumeKernel = T(0);  
   T pa = static_cast<T>(parameters(2));
   T ratio = static_cast<T>(parameters(1) / parameters(0));
   T major = static_cast<T>(parameters(0));
   if (kernelType==casacore::VectorKernel::GAUSSIAN) {
       _fillGaussian (maxValKernel, volumeKernel, kernelMatrix, height,
                     xCentre, yCentre, major, ratio, pa);
   } else if (kernelType==casacore::VectorKernel::BOXCAR) {
/*
      fillBoxcar (maxValKernel, volumeKernel, kernelMatrix, height,
                  xCentre, yCentre, major, ratio, pa);
*/
   } else {
     throw(casacore::AipsError("Unrecognized kernel type"));        // Earlier checking should prevent this
   }
   return volumeKernel;
}         

template <class T> void Image2DConvolver<T>::_fillGaussian(
    T& maxVal, T& volume, casacore::Matrix<T>& pixels, T height,
    T xCentre, T yCentre, T majorAxis, T ratio,
    T positionAngle
) const {
// 
// pa positive in +x ->+y pixel coordinate frame
//
   casacore::uInt n1 = pixels.shape()(0);
   casacore::uInt n2 = pixels.shape()(1);
   AlwaysAssert(n1==n2,casacore::AipsError);
   positionAngle += C::pi_2;        // +y -> -x
   casacore::Gaussian2D<T> g2d(height, xCentre, yCentre, majorAxis,
               ratio, positionAngle);
   maxVal = -1.0e30;
   volume = 0.0;
   casacore::Vector<T> pos(2);
   for (casacore::uInt j=0; j<n1; j++) {
      pos(1) = static_cast<T>(j);
      for (casacore::uInt i=0; i<n1; i++) {
         pos(0) = static_cast<T>(i);
         T val = g2d(pos);
         pixels(i,j) = val;
//
         maxVal = max(val, maxVal);
         volume += val;
      }
   } 
}

}

