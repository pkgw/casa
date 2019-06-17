#include <image_cmpt.h>

#include <iostream>
#include <sys/wait.h>
#include <casa/Arrays/ArrayIO.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/BasicMath/Random.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Exceptions/Error.h>
#include <casa/fstream.h>
#include <casa/BasicSL/STLIO.h>
#include <casa/Logging/LogFilter.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogOrigin.h>
#include <casa/OS/Directory.h>
#include <casa/OS/EnvVar.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/RegularFile.h>
#include <casa/OS/SymLink.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Utilities/Assert.h>

#include <images/Images/ImageExpr.h>
#include <images/Images/ImageExprParse.h>
#include <images/Images/ImageFITSConverter.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/ImageStatistics.h>
#include <images/Images/ImageSummary.h>
#include <images/Images/ImageUtilities.h>
#include <images/Images/LELImageCoord.h>
#include <images/Images/PagedImage.h>
#include <images/Images/RebinImage.h>
#include <images/Images/TempImage.h>
#include <images/Regions/ImageRegion.h>
#include <images/Regions/WCLELMask.h>
#include <lattices/LatticeMath/Fit2D.h>
#include <lattices/LatticeMath/LatticeFit.h>
#include <lattices/LEL/LatticeExprNode.h>
#include <lattices/Lattices/LatticeIterator.h>
#include <lattices/LRegions/LatticeRegion.h>
#include <lattices/LatticeMath/LatticeSlice1D.h>
#include <lattices/Lattices/LatticeUtilities.h>
#include <lattices/LRegions/LCBox.h>
#include <lattices/LRegions/LCSlicer.h>
#include <lattices/Lattices/MaskedLatticeIterator.h>
#include <lattices/Lattices/PixelCurve1D.h>
#include <lattices/LRegions/RegionType.h>
#include <lattices/Lattices/TiledLineStepper.h>
#include <measures/Measures/Stokes.h>
#include <measures/Measures/MeasIERS.h>
#include <scimath/Fitting/LinearFitSVD.h>
#include <scimath/Functionals/Polynomial.h>
#include <scimath/Mathematics/VectorKernel.h>
#include <tables/LogTables/NewFile.h>

#include <components/ComponentModels/GaussianDeconvolver.h>
#include <components/ComponentModels/SkyCompRep.h>
#include <components/SpectralComponents/SpectralListFactory.h>

#include <imageanalysis/ImageAnalysis/BeamManipulator.h>
#include <imageanalysis/ImageAnalysis/CasaImageBeamSet.h>
#include <imageanalysis/ImageAnalysis/ComponentImager.h>
#include <imageanalysis/ImageAnalysis/ComponentListDeconvolver.h>
#include <imageanalysis/ImageAnalysis/Image2DConvolver.h>
#include <imageanalysis/ImageAnalysis/ImageBoxcarSmoother.h>
#include <imageanalysis/ImageAnalysis/ImageCollapser.h>
#include <imageanalysis/ImageAnalysis/ImageConcatenator.h>
#include <imageanalysis/ImageAnalysis/ImageConvolverTask.h>
#include <imageanalysis/ImageAnalysis/ImageCropper.h>
#include <imageanalysis/ImageAnalysis/ImageDecimator.h>
#include <imageanalysis/ImageAnalysis/ImageDecomposerTask.h>
#include <imageanalysis/ImageAnalysis/ImageExprCalculator.h>
#include <imageanalysis/ImageAnalysis/ImageFactory.h>
#include <imageanalysis/ImageAnalysis/ImageFFTer.h>
#include <imageanalysis/ImageAnalysis/ImageFitter.h>
#include <imageanalysis/ImageAnalysis/ImageHanningSmoother.h>
#include <imageanalysis/ImageAnalysis/ImageHistogramsCalculator.h>
#include <imageanalysis/ImageAnalysis/ImageHistory.h>
#include <imageanalysis/ImageAnalysis/ImageMaskedPixelReplacer.h>
#include <imageanalysis/ImageAnalysis/ImageMaskHandler.h>
#include <imageanalysis/ImageAnalysis/ImageMaxFitter.h>
#include <imageanalysis/ImageAnalysis/ImageMetaDataRW.h>
#include <imageanalysis/ImageAnalysis/ImageMomentsTask.h>
#include <imageanalysis/ImageAnalysis/ImagePadder.h>
#include <imageanalysis/ImageAnalysis/ImageProfileFitter.h>
#include <imageanalysis/ImageAnalysis/ImagePrimaryBeamCorrector.h>
#include <imageanalysis/ImageAnalysis/ImageRebinner.h>
#include <imageanalysis/ImageAnalysis/ImageRegridder.h>
#include <imageanalysis/ImageAnalysis/ImageRotator.h>
#include <imageanalysis/ImageAnalysis/ImageStatsCalculator.h>
#include <imageanalysis/ImageAnalysis/ImageTransposer.h>
#include <imageanalysis/ImageAnalysis/PeakIntensityFluxDensityConverter.h>
#include <imageanalysis/ImageAnalysis/PixelValueManipulator.h>
#include <imageanalysis/ImageAnalysis/PVGenerator.h>
#include <imageanalysis/ImageAnalysis/SepImageConvolverTask.h>
#include <imageanalysis/ImageAnalysis/StatImageCreator.h>
#include <imageanalysis/ImageAnalysis/SubImageFactory.h>
#include <imageanalysis/ImageAnalysis/TwoPointCorrelator.h>

#include <componentlist_cmpt.h>

#include <stdcasa/version.h>

#include <casa/namespace.h>

#include <memory>

using namespace std;

#define _ORIGIN LogOrigin(_class, __func__, WHERE)

using namespace casacore;
using namespace casa;

namespace casac {

const String image::_class = "image";

image::image() : _log(), _imageF(), _imageC() {}

// private ImageInterface constructor for on the fly components
// The constructed object will take over management of the provided pointer
// using a shared_ptr

image::image(casacore::ImageInterface<casacore::Float> *inImage) :
    _log(), _imageF(inImage), _imageC(), _imageD(), _imageDC() {
}

image::image(ImageInterface<Complex> *inImage) :
    _log(), _imageF(), _imageC(inImage), _imageD(), _imageDC() {
}

image::image(casacore::ImageInterface<casacore::Double> *inImage) :
    _log(), _imageF(), _imageC(), _imageD(inImage), _imageDC() {
}

image::image(ImageInterface<DComplex> *inImage) :
    _log(), _imageF(), _imageC(), _imageD(), _imageDC(inImage) {
}

image::image(casa::SPIIF inImage) :
    _log(), _imageF(inImage), _imageC(), _imageD(), _imageDC() {
}

image::image(casa::SPIIC inImage) :
    _log(), _imageF(), _imageC(inImage), _imageD(), _imageDC() {
}

image::image(casa::SPIID inImage)
    : _imageD(inImage) {}

image::image(casa::SPIIDC inImage)
    : _imageDC(inImage) {}

image::image(ITUPLE mytuple) {
    _setImage(mytuple);
}


image::~image() {}

image* image::adddegaxes(
    const std::string& outfile, bool direction,
    bool spectral, const std::string& stokes, bool linear,
    bool tabular, bool overwrite, bool silent
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        if (_imageF) {
            casa::SPCIIF myim = _imageF;
            return _adddegaxes(
                myim, outfile, direction, spectral, stokes,
                linear, tabular, overwrite, silent
            );
        }
        else if (_imageD) {
            casa::SPCIID myim = _imageD;
            return _adddegaxes(
                myim, outfile, direction, spectral, stokes,
                linear, tabular, overwrite, silent
            );
        }
        else if (_imageC) {
            casa::SPCIIC myim = _imageC;
            return _adddegaxes(
                myim, outfile, direction, spectral, stokes,
                linear, tabular, overwrite, silent
            );
        }
        else if (_imageDC) {
            casa::SPCIIC myim = _imageC;
            return _adddegaxes(
                myim, outfile, direction, spectral, stokes,
                linear, tabular, overwrite, silent
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> image* image::_adddegaxes(
    SPCIIT inImage,
    const std::string& outfile, bool direction,
        bool spectral, const std::string& stokes, bool linear,
        bool tabular, bool overwrite, bool silent
) {
    _log << _ORIGIN;
    PtrHolder<ImageInterface<T> > outimage;
    ImageUtilities::addDegenerateAxes(
        _log, outimage, *inImage, outfile,
        direction, spectral, stokes, linear, tabular, overwrite, silent
    );
    auto *outPtr = outimage.ptr();
    outimage.clear(false);
    SPIIT z(outPtr);
    vector<String> names {
        "outfile", "direction", "spectral", "stokes",
        "linear", "tabular", "overwrite", "silent"
    };
    vector<variant> values {
        outfile, direction, spectral, stokes,
        linear, tabular, overwrite, silent
    };
    _addHistory(z, "adddegaxes", names, values);
    return new image(z);
}

bool image::addnoise(
    const std::string& type, const std::vector<double>& pars,
    const variant& region, bool zeroIt, const vector<int>& seeds
) {
    try {
        _log << LogOrigin("image", __func__);
        if (_detached()) {
            return false;
        }
        auto pRegion = _getRegion(region, false);
        std::shared_ptr<std::pair<Int, Int> > seedPair(
            new std::pair<Int, Int>(0, 0)
        );
        if (seeds.size() >= 2) {
            seedPair->first = seeds[0];
            seedPair->second = seeds[1];
        }
        else {
            Time now;
            Double seedBase = 1e7*now.modifiedJulianDay();
            seedPair->second = (Int)((uInt)seedBase);
            seedPair->first = seeds.size() == 1
                ? seeds[0]
                : (Int)((uInt)(1e7*(seedBase - seedPair->second)));
        }
        if (_imageF) {
            PixelValueManipulator<Float>::addNoise(
                _imageF, type, *pRegion,
                pars, zeroIt, seedPair.get()
            );
        }
        else if (_imageC) {
            PixelValueManipulator<Complex>::addNoise(
                _imageC, type, *pRegion,
                pars, zeroIt, seedPair.get()
            );
        }
        else if (_imageD) {
            PixelValueManipulator<Double>::addNoise(
                _imageD, type, *pRegion,
                pars, zeroIt, seedPair.get()
            );
        }
        else if (_imageDC) {
            PixelValueManipulator<DComplex>::addNoise(
                _imageDC, type, *pRegion,
                pars, zeroIt, seedPair.get()
            );
        }
        else {
            ThrowCc("Logic error");
        }
        vector<String> names { "type", "pars", "region", "zeroit", "seeds" };
        vector<variant> values { type, pars, region, zeroIt, seeds };
        _addHistory(__func__, names, values);
        _statsF.reset();
        _statsD.reset();
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

record* image::beamarea(int channel, int polarization) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        auto dc = _imageF
            ? _imageF->coordinates().directionCoordinate()
            : _imageC->coordinates().directionCoordinate();
        auto pixelArea = dc.getPixelArea();
        auto beamInPixels = _imageF
            ? _imageF->imageInfo().getBeamAreaInPixels(
                channel, polarization, dc
            )
            : _imageC->imageInfo().getBeamAreaInPixels(
                    channel, polarization, dc
            );
        auto arcsec2 = beamInPixels*pixelArea;
        record *rec = new record();
        rec->insert("pixels", beamInPixels);
        rec->insert("arcsec2", arcsec2.getValue("arcsec2"));
        return rec;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::beamforconvolvedsize(
    const variant& source, const variant& convolved
) {
    try {
        _log << _ORIGIN;
        Vector<casacore::Quantity> sourceParam, convolvedParam;
        if (
            ! toCasaVectorQuantity(source, sourceParam)
            || sourceParam.size() != 3
        ) {
            throw(AipsError("Cannot understand source values"));
        }
        if (
            ! toCasaVectorQuantity(convolved, convolvedParam)
            && convolvedParam.size() != 3
        ) {
            throw(AipsError("Cannot understand target values"));
        }
        Angular2DGaussian mySource(sourceParam[0], sourceParam[1], sourceParam[2]);
        GaussianBeam myConvolved(convolvedParam[0], convolvedParam[1], convolvedParam[2]);
        GaussianBeam neededBeam;
        try {
            if (GaussianDeconvolver::deconvolve(neededBeam, myConvolved, mySource)) {
                // throw without a message here, it will be caught
                // in the associated catch block and a new error will
                // be thrown with the appropriate message.
                throw AipsError();
            }
        }
        catch (const AipsError& x) {
            ostringstream os;
            os << "Unable to reach target resolution of "
                << myConvolved << " Input source "
                << mySource << " is probably too large.";
            throw AipsError(os.str());
        }
        Record ret;
        QuantumHolder qh(neededBeam.getMajor());
        ret.defineRecord("major", qh.toRecord());
        qh = QuantumHolder(neededBeam.getMinor());
        ret.defineRecord("minor", qh.toRecord());
        qh = QuantumHolder(neededBeam.getPA());
        ret.defineRecord("pa", qh.toRecord());
        return fromRecord(ret);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::boundingbox(const variant& region) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        if (_imageF) {
            return _boundingbox(_imageF, region);
        }
        else if (_imageC) {
            return _boundingbox(_imageC, region);
        }
        else {
            ThrowCc("Logic error");
        }
   }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_boundingbox(
    SPIIT image, const variant& region
) const {
    auto myreg = _getRegion(region, false);
    ImageMetaData<T> md(image);
    return fromRecord(*md.getBoundingBox(*myreg));
}

image* image::boxcar(
    const string& outfile, const variant& region,
    const variant& vmask, int axis, int width, bool drop,
    const string& dmethod,
    bool overwrite, bool stretch
) {
    LogOrigin lor(_class, __func__);
    _log << lor;
    if (_detached()) {
        throw AipsError("Unable to create image");
    }
    try {
        _notSupported(__func__);
        if (axis < 0) {
            const auto& csys = _imageF
                ? _imageF->coordinates()
                : _imageC->coordinates();
            ThrowIf(
                ! csys.hasSpectralAxis(),
                "Axis not specified and image has no spectral coordinate"
            );
            axis = csys.spectralAxisNumber(false);
        }
        if (_imageF) {
            SPCIIF image = _imageF;
            return _boxcar(
                image, region, vmask, outfile,
                overwrite, stretch, axis, width, drop,
                dmethod, lor
            );
        }
        else {
            SPCIIC image = _imageC;
            return _boxcar(
                image, region, vmask, outfile,
                overwrite, stretch, axis, width, drop,
                dmethod, lor
            );
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> image* image::_boxcar(
    SPCIIT myimage, const variant& region,
    const variant& mask, const string& outfile, bool overwrite,
    bool stretch, int axis, int width, bool drop,
    const string& dmethod, const LogOrigin& lor
) {
    ImageBoxcarSmoother<T> smoother(
        myimage, _getRegion(region, true).get(),
        _getMask(mask), outfile, overwrite
    );
    smoother.setAxis(axis);
    smoother.setDecimate(drop);
    smoother.setStretch(stretch);
    smoother.setWidth(width);
    ImageDecimatorData::Function dFunction = ImageDecimatorData::NFUNCS;
    if (drop) {
        String mymethod = dmethod;
        mymethod.downcase();
        if (mymethod.startsWith("m")) {
            dFunction = ImageDecimatorData::MEAN;
        }
        else if (mymethod.startsWith("c")) {
            dFunction = ImageDecimatorData::COPY;
        }
        else {
            ThrowCc(
                "Value of dmethod must be "
                    "either 'm'(ean) or 'c'(opy)"
            );
        }
        smoother.setDecimationFunction(dFunction);
    }
    vector<String> names {
        "outfile", "region", "mask", "axis", "width",
        "drop", "dmethod", "overwrite", "stretch"
    };
    vector<variant> values {
        outfile, region, mask, axis, width,
        drop, dmethod, overwrite, stretch
    };
    if (_doHistory) {
        auto msgs = _newHistory("boxcar", names, values);
        smoother.addHistory(lor, msgs);
    }
    return new image(smoother.smooth());
}

std::string image::brightnessunit() {
    if (_detached()) {
        return "";
    }
    try {
        if (_imageF) {
            return _imageF->units().getName();
        }
        else if (_imageC) {
            return _imageC->units().getName();
        }
        else if (_imageD) {
            return _imageD->units().getName();
        }
        else if (_imageDC) {
            return _imageDC->units().getName();
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return "";
}

bool image::calc(const std::string& expr, bool verbose) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        if (_imageF) {
            ImageExprCalculator<Float>::compute2(_imageF, expr, verbose);
        }
        else {
            ImageExprCalculator<Complex>::compute2(_imageC, expr, verbose);
        }
        vector<String> names = {"expr", "verbose"};
        vector<variant> values = {expr, verbose};
        _addHistory(__func__, names, values);
        _statsF.reset();
        _statsD.reset();
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::calcmask(
    const string& mask, const string& maskName, bool makeDefault
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        Record region;
        if (_imageF) {
            ImageMaskHandler<Float> imh(_imageF);
            imh.calcmask(mask, region, maskName, makeDefault);
        }
        else if (_imageC){
            ImageMaskHandler<Complex> imh(_imageC);
            imh.calcmask(mask, region, maskName, makeDefault);
        }
        else if (_imageD) {
            ImageMaskHandler<Double> imh(_imageD);
            imh.calcmask(mask, region, maskName, makeDefault);
        }
        else if (_imageDC) {
            ImageMaskHandler<DComplex> imh(_imageDC);
            imh.calcmask(mask, region, maskName, makeDefault);
        }
        else {
            ThrowCc("Logic error");
        }
        vector<String> names {"mask", "name", "asdefault"};
        vector<variant> values {mask, maskName, makeDefault};
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::close() {
    try {
        _log << _ORIGIN;
        _reset();
        MeasIERS::closeTables();
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

image* image::collapse(
    const string& function, const variant& axes,
    const string& outfile, const variant& region, const string& box,
    const string& chans, const string& stokes, const string& mask,
    const bool overwrite, const bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        return 0;
    }
    try {
        _notSupported(__func__);
        IPosition myAxes;
        auto axesType = axes.type();
        ThrowIf(axesType == variant::BOOLVEC, "axes must be specified");
        if (axesType == variant::INT) {
            myAxes = IPosition(1, axes.toInt());
        }
        else if (axesType == variant::INTVEC) {
            myAxes = IPosition(axes.getIntVec());
        }
        else if (
            axesType == variant::STRINGVEC
            || axesType == variant::STRING
        ) {
            Vector<String> axVec = (axes.type() == variant::STRING)
                ? Vector<String> (1, axes.getString())
                : toVectorString(axes.toStringVec());
            myAxes = IPosition(
                _imageF
                ? _imageF->coordinates().getWorldAxesOrder(
                    axVec, false
                )
                : _imageC->coordinates().getWorldAxesOrder(
                    axVec, false
                )
            );
            for (
                IPosition::iterator iter = myAxes.begin();
                iter != myAxes.end(); iter++
            ) {
                ThrowIf(
                    *iter < 0,
                    "At least one specified axis does not exist"
                );
            }
        }
        else {
            ThrowCc("Unsupported type for parameter axes");
        }
        std::shared_ptr<Record> regRec = _getRegion(region, true);
        String aggString = function;
        aggString.trim();
        aggString.downcase();
        ThrowIf(
            aggString == "avdev",
            "avdev currently not supported. Let us know if you have a need for it"
        );
        vector<String> names {
            "function", "axes", "outfile", "region", "box",
            "chans", "stokes", "mask", "overwrite", "stretch"
        };
        vector<variant> values {
            function, axes, outfile, region, box,
            chans, stokes, mask, overwrite, stretch
        };

        String myStokes = stokes;
        if (_imageF) {
            CasacRegionManager rm(_imageF->coordinates());
            String diagnostics;
            uInt nSelectedChannels;
            Record myreg = rm.fromBCS(
                diagnostics, nSelectedChannels, myStokes, regRec.get(),
                "", chans, CasacRegionManager::USE_ALL_STOKES, box,
                _imageF->shape(), mask, false
            );
            ImageCollapser<Float> collapser(
                aggString, _imageF, &myreg,
                mask, myAxes, false, outfile, overwrite
            );
            collapser.setStretch(stretch);
            if (_doHistory) {
                collapser.addHistory(_ORIGIN, "ia." + String(__func__), names, values);
            }
            return new image(collapser.collapse());
        }
        else {
            CasacRegionManager rm(_imageC->coordinates());
            String diagnostics;
            uInt nSelectedChannels;
            Record myreg = rm.fromBCS(
                diagnostics, nSelectedChannels, myStokes, regRec.get(),
                "", chans, CasacRegionManager::USE_ALL_STOKES, box,
                _imageC->shape(), mask, false
            );
            ImageCollapser<Complex> collapser(
                aggString, _imageC, &myreg,
                mask, myAxes, false, outfile, overwrite
            );
            collapser.setStretch(stretch);
            if (_doHistory) {
                collapser.addHistory(_ORIGIN, "ia." + String(__func__), names, values);
            }
            return new image(collapser.collapse());
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::commonbeam() {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        ImageInfo myInfo = _imageF ? _imageF->imageInfo() : _imageC->imageInfo();
        ThrowIf(
            ! myInfo.hasBeam(),
            "This image has no beam(s)."
        );
        GaussianBeam beam;
        if (myInfo.hasSingleBeam()) {
            _log << LogIO::WARN
                << "This image only has one beam, so just returning that"
                << LogIO::POST;
            beam = myInfo.restoringBeam();
        }
        else {
            // multiple beams in this image
            beam = CasaImageBeamSet(myInfo.getBeamSet()).getCommonBeam();
        }
        beam.setPA(casacore::Quantity(beam.getPA("deg", true), "deg"));
        Record x = beam.toRecord();
        x.defineRecord("pa", x.asRecord("positionangle"));
        x.removeField("positionangle");
        return fromRecord(x);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::continuumsub(
    const string& outline, const string& outcont,
    const variant& region, const vector<int>& channels,
    const string& pol, const int in_fitorder, const bool overwrite
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        ThrowIf(in_fitorder < 0, "Polynomial order cannot be negative");
        if (! pol.empty()) {
            _log << LogIO::NORMAL << "The pol parameter is no longer "
                << "supported and will be removed in the near future. "
                << "Please set the region parameter appropriately "
                << "to select the polarization in which you are interested."
                << LogIO::POST;
        }
        std::shared_ptr<Record> leRegion = _getRegion(region, false);
        vector<Int> planes = channels;
        if (planes.size() == 1 && planes[0] == -1) {
            planes.resize(0);
        }
        Int spectralAxis = _imageF->coordinates().spectralAxisNumber();
        ThrowIf(spectralAxis < 0, "This image has no spectral axis");
        ImageProfileFitter fitter(
        _imageF, "", leRegion.get(),
            "", "", "", "", spectralAxis,
            0, overwrite
        );
        fitter.setDoMultiFit(true);
        fitter.setPolyOrder(in_fitorder);
        fitter.setModel(outcont);
        fitter.setResidual(outline);
        fitter.setStretch(false);
        fitter.setLogResults(false);
        if (! planes.empty()) {
            std::set<int> myplanes(planes.begin(), planes.end());
            ThrowIf(*myplanes.begin() < 0, "All planes must be nonnegative");
            fitter.setGoodPlanes(std::set<uInt>(myplanes.begin(), myplanes.end()));
        }
        fitter.createResidualImage(true);
        vector<String> names {
            "outline", "outcont", "region", "channels",
            "pol", "fitorder", "overwrite"
        };
        vector<variant> values {
            outline, outcont, region, channels,
            pol, in_fitorder, overwrite
        };
        if (_doHistory) {
            auto msgs = _newHistory(__func__, names, values);
            fitter.addHistory(_ORIGIN, msgs);
        }
        fitter.fit(false);
        return new image(fitter.getResidual());
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::convertflux(
    const variant& qvalue, const variant& major,
    const variant& minor,  const string& /*type*/,
    const bool toPeak,
    const int channel, const int polarization
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return 0;
        }
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        casacore::Quantity value = casaQuantity(qvalue);
        casacore::Quantity majorAxis = casaQuantity(major);
        casacore::Quantity minorAxis = casaQuantity(minor);
        Bool noBeam = false;
        PeakIntensityFluxDensityConverter<Float> converter(_imageF);
        converter.setSize(
            Angular2DGaussian(majorAxis, minorAxis, casacore::Quantity(0, "deg"))
        );
        converter.setBeam(channel, polarization);
        return recordFromQuantity(
            toPeak
            ? converter.fluxDensityToPeakIntensity(noBeam, value)
            : converter.peakIntensityToFluxDensity(noBeam, value)
        );
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: "
            << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::convolve(
    const string& outfile, const variant& kernel,
    double scale, const variant& region,
    const variant& vmask, bool overwrite,
    bool stretch
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! (_imageF || _imageD),
            "This method only supports real-valued images"
        );
        if (_imageF) {
            return _convolve(
                _imageF, outfile, kernel, scale,
                region, vmask, overwrite, stretch
            );
        }
        /*
        else if (_imageC) {
            return _convolve(
                _imageC, outfile, kernel, scale,
                region, vmask, overwrite, stretch
            );
        }
        */
        else if (_imageD) {
            return _convolve(
                _imageD, outfile, kernel, scale,
                region, vmask, overwrite, stretch
            );
        }
        /*
        else if (_imageDC) {
            return _convolve(
                _imageDC, outfile, kernel, scale,
                region, vmask, overwrite, stretch
            );
        }
        */
        else {
            ThrowCc("Logic Error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> image* image::_convolve(
    SPIIT myImage, const string& outfile, const variant& kernel, double scale,
    const variant& region, const variant& vmask, bool overwrite, bool stretch
) {
    Array<T> fkernelArray;
    String kernelFileName = "";
    if (kernel.type() == variant::DOUBLEVEC) {
        const auto kernelVector = kernel.toDoubleVec();
        const auto shape = kernel.arrayshape();
        fkernelArray.resize(IPosition(shape));
        Vector<Double> localkern(kernelVector);
        convertArray(fkernelArray, localkern.reform(IPosition(shape)));
    }
    else if (kernel.type() == variant::INTVEC) {
        const auto kernelVector = kernel.toIntVec();
        const auto shape = kernel.arrayshape();
        fkernelArray.resize(IPosition(shape));
        Vector<Int> localkern(kernelVector);
        convertArray(fkernelArray, localkern.reform(IPosition(shape)));
    }
    else if (
        kernel.type() == variant::STRING
        || kernel.type() == variant::STRINGVEC
    ) {
        kernelFileName = kernel.toString();
        fkernelArray = PagedImage<T>(kernelFileName).get();
    }
    else {
        ThrowCc("kernel is not understood, try using an array or an image");
    }
    auto theMask = _getMask(vmask);
    std::shared_ptr<Record> myregion = _getRegion(region, false);
    ImageConvolverTask<T> ic(
        myImage, myregion.get(), theMask, outfile, overwrite
    );
    ic.setScale(scale);
    ic.setStretch(stretch);
    ic.setKernel(fkernelArray);
    if (_doHistory) {
        vector<String> names {
        "outfile", "kernel", "scale", "region", "vmask", "overwrite", "stretch"
        };
        vector<variant> values {
            outfile, kernel, scale, region,
            vmask, overwrite, stretch
        };
        auto msgs = _newHistory("convolve", names, values);
        ic.addHistory(LogOrigin(_class, "convolve"), msgs);
    }
    auto z = ic.convolve();
    return new image(z);
}

image* image::convolve2d(
    const string& outFile, const vector<int>& axes,
    const string& type, const variant& major, const variant& minor,
    const variant& pa, double in_scale, const variant& region,
    const variant& vmask, bool overwrite, bool stretch,
    bool targetres, const record& beam
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! (_imageF || _imageD),
            String(__func__) + " only supports real-valued images"
        );
        if (_imageF) {
            return _convolve2d(
                _imageF, outFile, axes, type, major, minor, pa, in_scale,
                region, vmask, overwrite, stretch, targetres, beam
            );
        }
        /*
        else if (_imageC) {
            return _convolve2d(
                _imageC, outFile, axes, type, major, minor, pa, in_scale,
                region, vmask, overwrite, stretch, targetres, beam
            );
        }
        */
        else if (_imageD) {
            return _convolve2d(
                _imageD, outFile, axes, type, major, minor, pa, in_scale,
                region, vmask, overwrite, stretch, targetres, beam
            );
        }
        /*
        else if (_imageDC) {
            return _convolve2d(
                _imageDC, outFile, axes, type, major, minor, pa, in_scale,
                region, vmask, overwrite, stretch, targetres, beam
            );
        }
        */
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: "
            << x.getMesg() << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> image* image::_convolve2d(
    SPIIT myImage, const string& outFile, const vector<int>& axes,
    const string& type, const variant& major, const variant& minor,
    const variant& pa, double in_scale, const variant& region,
    const variant& vmask, bool overwrite, bool stretch,
    bool targetres, const record& beam
) {
    UnitMap::putUser("pix", UnitVal(1.0), "pixel units");
    std::shared_ptr<Record> Region(_getRegion(region, false));
    auto mask = _getMask(vmask);
    String kernel(type);
    casacore::Quantity majorKernel;
    casacore::Quantity minorKernel;
    casacore::Quantity paKernel;
    _log << _ORIGIN;
    if (! beam.empty()) {
        ThrowIf(
            ! String(type).startsWith("g") && ! String(type).startsWith("G"),
            "beam can only be given with a gaussian kernel"
        );
        ThrowIf(
            ! major.toString(false).empty()
            || ! minor.toString(false).empty()
            || ! pa.toString(false).empty(),
            "major, minor, and/or pa may not be specified if beam is specified"
        );
        ThrowIf(
            beam.size() != 3,
            "If given, beam must have exactly three fields"
        );
        ThrowIf(
            beam.find("major") == beam.end(),
            "Beam must have a 'major' field"
        );
        ThrowIf(
            beam.find("minor") == beam.end(),
            "Beam must have a 'minor' field"
        );
        ThrowIf(
            beam.find("positionangle") == beam.end()
            && beam.find("pa") == beam.end(),
            "Beam must have a 'positionangle' or 'pa' field"
        );
        std::unique_ptr<Record> nbeam(toRecord(beam));
        for (uInt i=0; i<3; ++i) {
            String key = i == 0
                ? "major"
                : i == 1
                    ? "minor"
                    : beam.find("pa") == beam.end()
                        ? "positionangle"
                        : "pa";
            casacore::Quantity x;
            auto type = nbeam->dataType(nbeam->fieldNumber(key));
            String err;
            QuantumHolder z;
            Bool success;
            if (type == TpString) {
                success = z.fromString(err, nbeam->asString(key));
            }
            else if (type == TpRecord) {
                success = z.fromRecord(err, nbeam->asRecord(key));
            }
            else {
                ThrowCc("Unsupported data type for beam");
            }
            if (! success) {
                ThrowCc("Error converting beam to Quantity");
            }
            if (key == "major") {
                majorKernel = z.asQuantity();
            }
            else if (key == "minor") {
                minorKernel = z.asQuantity();
            }
            else {
                paKernel = z.asQuantity();
            }
        }
    }
    else {
        majorKernel = _casaQuantityFromVar(major);
        minorKernel = _casaQuantityFromVar(minor);
        paKernel = _casaQuantityFromVar(pa);
    }
    _log << _ORIGIN;
    Vector<Int> Axes(axes);
    if (Axes.size() == 0) {
        Axes.resize(2);
        Axes[0] = 0;
        Axes[1] = 1;
    }
    else {
        ThrowIf(
            axes.size() != 2,
            "Number of axes to convolve must be exactly 2"
        );
    }
    Image2DConvolver<T> convolver(
        myImage, Region.get(), mask, outFile, overwrite
    );
    convolver.setAxes(std::make_pair(Axes[0], Axes[1]));
    convolver.setKernel(type, majorKernel, minorKernel, paKernel);
    convolver.setScale(in_scale);
    convolver.setStretch(stretch);
    convolver.setTargetRes(targetres);
    if (_doHistory) {
        vector<String> names = {
            "outfile", "axes", "type", "major",
            "minor", "pa", "scale", "region",
            "mask", "overwrite", "stretch",
            "targetres", "beam"
        };
        vector<variant> values = {
            outFile, axes, type, major, minor,
            pa, in_scale, region, vmask,
            overwrite, stretch, targetres, beam
        };
        auto msgs = _newHistory("convolve2d", names, values);
        convolver.addHistory(_ORIGIN, msgs);
    }
    return new image(convolver.convolve());
}

record* image::coordmeasures(
    const std::vector<double>&pixel, const string& dframe,
    const string& sframe
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        casacore::Record theDir;
        casacore::Record theFreq;
        casacore::Record theVel;
        Vector<Double> vpixel;
        if (!(pixel.size() == 1 && pixel[0] == -1)) {
            vpixel = pixel;
        }
        unique_ptr<Record> retval;
        casacore::String error;
        Record R;
        if (_imageF) {
            casacore::Quantum<Float> intensity;
            retval.reset(
                PixelValueManipulator<Float>::coordMeasures(
                    intensity, theDir, theFreq, theVel,
                    _imageF, vpixel, dframe, sframe
                )
            );
            ThrowIf(
                ! QuantumHolder(intensity).toRecord(error, R),
                "Could not convert intensity to record. "
                + error
            );
        }
        else {
            casacore::Quantum<Complex> intensity;
            retval.reset(
                PixelValueManipulator<Complex>::coordMeasures(
                    intensity, theDir, theFreq, theVel,
                    _imageC, vpixel, dframe, sframe
                )
            );
            ThrowIf(
                ! QuantumHolder(intensity).toRecord(error, R),
                "Could not convert intensity to record. "
                + error
            );
        }
        retval->defineRecord(RecordFieldId("intensity"), R);
        return fromRecord(*retval);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

coordsys* image::coordsys(const std::vector<int>& pixelAxes) {
    _log << _ORIGIN;
    try {
        if (_detached()) {
            return nullptr;
        }
        if (_imageF) {
            return _coordsys(_imageF, pixelAxes);
        }
        else if (_imageC) {
            return _coordsys(_imageC, pixelAxes);
        }
        else if (_imageD) {
            return _coordsys(_imageD, pixelAxes);
        }
        else if (_imageDC) {
            return _coordsys(_imageDC, pixelAxes);
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> coordsys* image::_coordsys(
    SPIIT image, const std::vector<int>& pixelAxes
) {
    vector<Int> myAxes = pixelAxes;
    if (pixelAxes.size() == 1 && pixelAxes[0] == -1) {
        myAxes.clear();
    }
    ImageMetaData<T> imd(image);
    auto csys =  imd.coordsys(myAxes);
    std::unique_ptr<casac::coordsys> rstat(new ::casac::coordsys());
    rstat->setcoordsys(csys);
    return rstat.release();
}

image* image::crop(
    const string& outfile, const vector<int>& axes,
    bool overwrite, const variant& region, const string& box,
    const string& chans, const string& stokes, const string& mask,
    bool  stretch, bool wantreturn
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        if (axes.size() > 0) {
            std::set<int> saxes(axes.begin(), axes.end());
            if (*saxes.begin() < 0) {
                _log << "All axes values must be >= 0" << LogIO::EXCEPTION;
            }
        }
        std::set<uInt> saxes(axes.begin(), axes.end());
        std::shared_ptr<Record> regionPtr = _getRegion(region, true);
        ImageCropper<Float> cropper(
            _imageF, regionPtr.get(), box,
            chans, stokes, mask, outfile, overwrite
        );
        cropper.setStretch(stretch);
        cropper.setAxes(saxes);
        if (_doHistory) {
            vector<String> names {
                "outfile", "axes", "overwrite",
                "region", "box", "chans", "stokes",
                "mask", "stretch",  "wantreturn"
            };
            vector<variant> values {
                outfile, axes, overwrite,
                region, box, chans, stokes,
                mask, stretch,  wantreturn
            };
            auto msgs = _newHistory(__func__, names, values);
            cropper.addHistory(_ORIGIN, msgs);
        }
        auto out(cropper.crop(wantreturn));
        if (wantreturn) {
            return new image(out);
        }
        return nullptr;

    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::decimate(
    const string& outfile, int axis, int factor, const string& method,
    const variant& region, const string& mask, bool overwrite, bool stretch
) {
    try {
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        ThrowIf(
            axis < 0,
            "The value of axis cannot be negative"
        );
        ThrowIf(
            factor < 0,
            "The value of factor cannot be negative"
        );
        String mymethod = method;
        mymethod.downcase();
        ImageDecimatorData::Function f;
        if (mymethod.startsWith("c")) {
            f = ImageDecimatorData::COPY;
        }
        else if (mymethod.startsWith("m")) {
            f = ImageDecimatorData::MEAN;
        }
        else {
            ThrowCc("Unsupported decimation method " + method);
        }
        std::shared_ptr<Record> regPtr(_getRegion(region, true));
        vector<String> msgs;
        if (_doHistory) {
            vector<String> names {
                "outfile", "axis", "factor", "method",
                "region", "mask", "overwrite", "stretch"
            };
            vector<variant> values {
                outfile, axis, factor, method,
                region, mask, overwrite, stretch
            };
            msgs = _newHistory(__func__, names, values);
        }
        if (_imageF) {
            SPCIIF myim = _imageF;
            return _decimate(
                myim, outfile, axis, factor, f,
                regPtr, mask, overwrite, stretch, msgs
            );
        }
        else {
            SPCIIC myim = _imageC;
            return _decimate(
                myim, outfile, axis, factor, f,
                regPtr, mask, overwrite, stretch, msgs
            );
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> image* image::_decimate(
    const SPCIIT myimage,
    const string& outfile, int axis, int factor,
    ImageDecimatorData::Function f,
    const std::shared_ptr<Record> region,
    const string& mask, bool overwrite, bool stretch,
    const vector<String>& msgs
) const {
    ImageDecimator<T> decimator(
        myimage, region.get(),
        mask, outfile, overwrite
    );
    decimator.setFunction(f);
    decimator.setAxis(axis);
    decimator.setFactor(factor);
    decimator.setStretch(stretch);
    decimator.addHistory(_ORIGIN, msgs);
    SPIIT out = decimator.decimate();
    return new image(out);
}

record* image::decompose(
    const variant& region, const ::casac::variant& vmask,
    bool simple, double threshold, int ncontour, int minrange,
    int naxis, bool fit, double maxrms, int maxretry, int maxiter,
    double convcriteria, bool stretch
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(! _imageF, "This application supports only real-valued images");
        ThrowIf(
            threshold < 0,
            "Threshold = " + String::toString(threshold)
            + ". You must specify a nonnegative threshold"
        );
        auto Region = _getRegion(region, false);
        String mask = _getMask(vmask);
        Matrix<Int> blcs;
        Matrix<Int> trcs;
        casacore::Record outrec1;
        ImageDecomposerTask<Float> idt(_imageF, Region.get(), mask);
        idt.setSimple(simple);
        idt.setDeblendOptions(threshold, ncontour, minrange, naxis);
        idt.setFit(fit);
        idt.setFitOptions(maxrms, maxretry, maxiter, convcriteria);
        idt.setStretch(stretch);
        outrec1.define("components", idt.decompose(blcs, trcs));
        outrec1.define("blc", blcs);
        outrec1.define("trc", trcs);
        return fromRecord(outrec1);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::deconvolvecomponentlist(
    const record& complist, int channel, int polarization
) {
    _log << _ORIGIN;
    if (_detached()) {
        return nullptr;
    }
    try {
        _notSupported(__func__);
        std::unique_ptr<Record> compList(toRecord(complist));
        ComponentList cl, clOut;
        casacore::String err;
        ThrowIf(
            ! cl.fromRecord(err, *compList),
            "Input dictionary is not a valid component list: " + err
        );
        if (_imageF) {
            ComponentListDeconvolver<Float> cld(_imageF);
            clOut = cld.deconvolve(cl, channel, polarization);
        }
        else {
            ComponentListDeconvolver<Complex> cld(_imageC);
            clOut = cld.deconvolve(cl, channel, polarization);
        }
        Record rec;
        ThrowIf(
            ! clOut.toRecord(err, rec),
            "Cannot convert resulting component list to record: " + err
        );
        return fromRecord(rec);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::deconvolvefrombeam(
    const variant& source, const variant& beam
) {
    try {
        _log << _ORIGIN;
        Vector<casacore::Quantity> sourceParam, beamParam;
        Angular2DGaussian mySource;
        if (
            ! toCasaVectorQuantity(source, sourceParam)
            || (sourceParam.nelements() == 0)
            || sourceParam.nelements() > 3
        ) {
            throw(AipsError("Cannot understand source values"));
        }
        else {
            if (sourceParam.nelements() == 1) {
                sourceParam.resize(3, true);
                sourceParam[1] = sourceParam[0];
                sourceParam[2] = casacore::Quantity(0, "deg");
            }
            else if (sourceParam.nelements() == 2) {
                sourceParam.resize(3, true);
                sourceParam[2] = casacore::Quantity(0, "deg");
            }
            mySource = Angular2DGaussian(
                sourceParam[0], sourceParam[1], sourceParam[2]
            );
        }
        if (
            ! toCasaVectorQuantity(beam, beamParam)
            || (beamParam.nelements() == 0)) {
            throw(AipsError("Cannot understand beam values"));
        }
        else {
            if (beamParam.nelements() == 1) {
                beamParam.resize(3, true);
                beamParam[1] = beamParam[0];
                beamParam[2] = casacore::Quantity(0.0, "deg");
            }
            if (beamParam.nelements() == 2) {
                beamParam.resize(3, true);
                beamParam[2] = casacore::Quantity(0.0, "deg");
            }
        }
        GaussianBeam myBeam(beamParam[0], beamParam[1], beamParam[2]);
        Bool success = false;
        Angular2DGaussian decon;
        Bool retval = false;
        try {
            retval = GaussianDeconvolver::deconvolve(decon, mySource, myBeam);
            success = true;
        }
        catch (const AipsError& x) {
            retval = false;
            success = false;
        }
        Record deconval = decon.toRecord();
        deconval.defineRecord("pa", deconval.asRecord("positionangle"));
        deconval.removeField("positionangle");
        deconval.define("success", success);
        Record outrec1;
        outrec1.define("return", retval);
        outrec1.defineRecord("fit", deconval);
        return fromRecord(outrec1);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::_detached() const {
    if ( ! (_imageF || _imageC || _imageD || _imageDC)) {
        _log <<  _ORIGIN;
        _log << LogIO::SEVERE
            << "Image is detached - cannot perform operation." << endl
            << "Call image.open('filename') to reattach." << LogIO::POST;
        return true;
    }
    return false;
}

bool image::dohistory(bool enable) {
    _doHistory = enable;
    return True;
}

bool image::done(bool remove, bool verbose) {
    try {
        _log << _ORIGIN;
        MeasIERS::closeTables();
        if (remove && !_detached()) {
            // object _reset happens in _remove()
            _remove(verbose);
        }
        else {
            _reset();
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

record* image::findsources(
    int nMax, double cutoff, const variant& region,
    const variant& vmask, bool point, int width, bool absFind
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(! _imageF, "This application supports only float-valued images");
        std::shared_ptr<Record> Region(_getRegion(region, false));
        auto mask = _getMask(vmask);
        ImageSourceFinder<Float> sf(_imageF, Region.get(), mask);
        sf.setCutoff(cutoff);
        sf.setDoPoint(point);
        sf.setWidth(width);
        sf.setAbsFind(absFind);
        auto cl = sf.findSources(nMax);
        Record rec;
        casacore::String error;
        ThrowIf (
            ! cl.toRecord(error, rec),
            "Failed to convert component list to record: " + error
        );
        return fromRecord(rec);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::fft(
    const string& realOut, const string& imagOut, const string& ampOut,
    const string& phaseOut, const std::vector<int>& axes, const variant& region,
    const variant& vmask, bool stretch, const string& complexOut
) {
    try {
        _log << LogOrigin(_class, __func__);
        if (_detached()) {
            return false;
        }
        if (_imageF) {
            return _fft(
                _imageF, realOut, imagOut, ampOut, phaseOut,
                axes, region, vmask, stretch, complexOut
            );
        }
        else if (_imageC) {
            return _fft(
                _imageC, realOut, imagOut, ampOut, phaseOut,
                axes, region, vmask, stretch, complexOut
            );
        }
        else if (_imageD) {
            return _fft(
                _imageD, realOut, imagOut, ampOut, phaseOut,
                axes, region, vmask, stretch, complexOut
            );
        }
        else if (_imageDC) {
            return _fft(
                _imageDC, realOut, imagOut, ampOut, phaseOut,
                axes, region, vmask, stretch, complexOut
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

template<class T> bool image::_fft(
    SPIIT myImage, const string& realOut, const string& imagOut,
    const string& ampOut, const string& phaseOut, const std::vector<int>& axes,
    const variant& region, const variant& vmask, bool stretch,
    const string& complexOut
) {
    std::shared_ptr<Record> myregion(_getRegion(region, false));
    String mask = vmask.toString();
    if (mask == "[]") {
        mask = "";
    }
    Vector<uInt> leAxes(0);
    if (axes.size() > 1 || (axes.size() == 1 && axes[0] >= 0)) {
        leAxes.resize(axes.size());
        for (uInt i=0; i<axes.size(); i++) {
            ThrowIf(
                axes[i] < 0,
                "None of the elements of axes may be less than zero"
            );
            leAxes[i] = axes[i];
        }
    }
    vector<String> msgs;
    if (_doHistory) {
        vector<String> names = {
            "real", "imag", "amp", "phase", "axes",
            "region", "mask", "stretch", "complex"
        };
        vector<variant> values = {
            realOut, imagOut, ampOut, phaseOut, axes,
            region, vmask, stretch, complexOut
        };
        msgs = _newHistory("fft", names, values);
    }
    ImageFFTer<T> ffter(myImage, myregion.get(), mask, leAxes);
    ffter.setStretch(stretch);
    ffter.setReal(realOut);
    ffter.setImag(imagOut);
    ffter.setAmp(ampOut);
    ffter.setPhase(phaseOut);
    ffter.setComplex(complexOut);
    if (_doHistory) {
        ffter.addHistory(_ORIGIN, msgs);
    }
    ffter.fft();
    return true;
}

record* image::fitcomponents(
    const string& box, const variant& region, const variant& chans,
    const string& stokes, const variant& vmask,
    const vector<double>& in_includepix, const vector<double>& in_excludepix,
    const string& residual, const string& model, const string& estimates,
    const string& logfile, bool append, const string& newestimates,
    const string& complist, bool overwrite, bool dooff, double offset,
    bool fixoffset, bool stretch, const variant& rms, const variant& noisefwhm,
    const string& summary
) {
    if (_detached()) {
        return nullptr;
    }
    _log << _ORIGIN;
    try {
        ThrowIf(
            ! (_imageF || _imageD),
            "This method only supports real valued images"
        );
        if (_imageF) {
            return _fitcomponents(
                _imageF, box, region, chans, stokes, vmask, in_includepix,
                in_excludepix, residual, model, estimates, logfile, append,
                newestimates, complist, overwrite, dooff, offset, fixoffset,
                stretch, rms, noisefwhm, summary
            );
        }
        else if (_imageD) {
            return _fitcomponents(
                _imageD, box, region, chans, stokes, vmask, in_includepix,
                in_excludepix, residual, model, estimates, logfile, append,
                newestimates, complist, overwrite, dooff, offset, fixoffset,
                stretch, rms, noisefwhm, summary
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << "Exception Reported: " << x.getMesg()
            << LogIO::EXCEPTION;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_fitcomponents(
    SPIIT myImage, const string& box, const variant& region,
    const variant& chans, const string& stokes, const variant& vmask,
    const vector<double>& in_includepix, const vector<double>& in_excludepix,
    const string& residual, const string& model, const string& estimates,
    const string& logfile, const bool append, const string& newestimates,
    const string& complist, bool overwrite, bool dooff, double offset,
    bool fixoffset, bool stretch, const variant& rms, const variant& noisefwhm,
    const string& summary
) {
    auto num = in_includepix.size();
    Vector<Float> includepix(num);
    num = in_excludepix.size();
    Vector<Float> excludepix(num);
    convertArray(includepix, Vector<Double> (in_includepix));
    convertArray(excludepix, Vector<Double> (in_excludepix));
    if (includepix.size() == 1 && includepix[0] == -1) {
        includepix.resize();
    }
    if (excludepix.size() == 1 && excludepix[0] == -1) {
        excludepix.resize();
    }
    auto mask = _getMask(vmask);
    auto writeControl = complist.empty()
        ? ImageFitterResults<T>::NO_WRITE
        : overwrite
          ? ImageFitterResults<T>::OVERWRITE
          : ImageFitterResults<T>::WRITE_NO_REPLACE;
    String sChans;
    if (chans.type() == variant::BOOLVEC) {
        // for some reason which eludes me, the default variant type is boolvec
        sChans = "";
    }
    else if (chans.type() == variant::STRING) {
        sChans = chans.toString();
      }
    else if (chans.type() == variant::INT) {
        sChans = String::toString(chans.toInt());
    }
    else {
        ThrowCc(
            "Unsupported type for chans. chans must "
            "be either an integer or a string"
        );
    }
    auto regionRecord = _getRegion(region, true);
    auto doImages = ! residual.empty() || ! model.empty();
    ImageFitter<T> fitter(
         myImage, "", regionRecord.get(), box, sChans,
         stokes, mask, estimates, newestimates, complist
    );
    if (includepix.size() == 1) {
        fitter.setIncludePixelRange(
            std::make_pair(includepix[0],includepix[0])
        );
    }
    else if (includepix.size() == 2) {
        fitter.setIncludePixelRange(
            std::make_pair(includepix[0],includepix[1])
        );
    }
    if (excludepix.size() == 1) {
        fitter.setExcludePixelRange(
            std::make_pair(excludepix[0],excludepix[0])
        );
    }
    else if (excludepix.size() == 2) {
        fitter.setExcludePixelRange(
            std::make_pair(excludepix[0],excludepix[1])
        );
    }
    fitter.setWriteControl(writeControl);
    fitter.setStretch(stretch);
    fitter.setModel(model);
    fitter.setResidual(residual);
    if (! logfile.empty()) {
        fitter.setLogfile(logfile);
        fitter.setLogfileAppend(append);
    }
    if (dooff) {
        fitter.setZeroLevelEstimate(offset, fixoffset);
    }
    auto myrms = (rms.type() == variant::DOUBLE || rms.type() == variant::INT)
        ? casacore::Quantity(rms.toDouble(), brightnessunit())
        : _casaQuantityFromVar(rms);
    if (myrms.getValue() > 0) {
        fitter.setRMS(myrms);
    }
    auto noiseType = noisefwhm.type();
    if (noiseType == variant::DOUBLE || noiseType == variant::INT) {
        fitter.setNoiseFWHM(noisefwhm.toDouble());
    }
    else if (noiseType == variant::BOOLVEC) {
        fitter.clearNoiseFWHM();
    }
    else if (
        noiseType == variant::STRING || noiseType == variant::RECORD
    ) {
        if (noiseType == variant::STRING && noisefwhm.toString().empty()) {
            fitter.clearNoiseFWHM();
        }
        else {
            fitter.setNoiseFWHM(_casaQuantityFromVar(noisefwhm));
        }
    }
    else {
        ThrowCc(
            "Unsupported data type for noisefwhm: " + noisefwhm.typeString()
        );
    }
    if (doImages && _doHistory) {
        vector<casacore::String> names {
            "box", "region", "chans", "stokes", "mask", "includepix",
            "excludepix", "residual", "model", "estimates", "logfile",
            "append", "newestimates", "complist", "dooff", "offset",
            "fixoffset", "stretch", "rms", "noisefwhm"
        };
        vector<variant> values {
            box, region, chans, stokes, vmask, in_includepix,
            in_excludepix, residual, model, estimates, logfile,
            append, newestimates, complist, dooff, offset, fixoffset,
            stretch, rms, noisefwhm
        };
        auto msgs = _newHistory("fitcomponents", names, values);
        fitter.addHistory(_ORIGIN, msgs);
    }
    fitter.setSummaryFile(summary);
    auto compLists = fitter.fit();
    return fromRecord(fitter.getOutputRecord());
}

record* image::fitprofile(const string& box, const variant& region,
    const string& chans, const string& stokes, int axis,
    const variant& vmask, int ngauss, int poly,
    const string& estimates, int minpts, bool multifit,
    const string& model, const string& residual, const string& amp,
    const string& amperr, const string& center, const string& centererr,
    const string& fwhm, const string& fwhmerr, const string& integral,
    const string& integralerr, bool stretch, bool logResults,
    const variant& pampest, const variant& pcenterest,
    const variant& pfwhmest, const variant& pfix, const variant& gmncomps,
    const variant& gmampcon, const variant& gmcentercon,
    const variant& gmfwhmcon, const vector<double>& gmampest,
    const vector<double>& gmcenterest, const vector<double>& gmfwhmest,
    const variant& gmfix, const string& spxtype, const vector<double>& spxest,
    const vector<bool>& spxfix, const variant& div, const string& spxsol,
    const string& spxerr, const string& logfile, bool append,
    const variant& pfunc, const vector<double>& goodamprange,
    const vector<double>& goodcenterrange,
    const vector<double>& goodfwhmrange, const variant& sigma,
    const string& outsigma, const vector<int>& planes
) {
    _log << LogOrigin(_class, __func__);
    if (_detached()) {
        return 0;
    }
    try {
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        String regionName;
        std::shared_ptr<Record> regionPtr = _getRegion(region, true);
        if (ngauss < 0) {
            _log << LogIO::WARN
                << "ngauss < 0 is meaningless. Setting ngauss = 0 "
                << LogIO::POST;
            ngauss = 0;
        }
        vector<double> mygoodamps = toVectorDouble(goodamprange, "goodamprange");
        if (mygoodamps.size() > 2) {
            _log << "Too many elements in goodamprange" << LogIO::EXCEPTION;
        }
        vector<double> mygoodcenters = toVectorDouble(goodcenterrange, "goodcenterrange");
        if (mygoodcenters.size() > 2) {
            _log << "Too many elements in goodcenterrange" << LogIO::EXCEPTION;
        }
        vector<double> mygoodfwhms = toVectorDouble(goodfwhmrange, "goodcenterrange");
        if (mygoodfwhms.size() > 2) {
            _log << "Too many elements in goodfwhmrange" << LogIO::EXCEPTION;
        }
        String mask = _getMask(vmask);
        std::unique_ptr<Array<Float> > sigmaArray;
        std::unique_ptr<PagedImage<Float> > sigmaImage;
        if (sigma.type() == variant::STRING) {
            String sigmaName = sigma.toString();
            if (! sigmaName.empty()) {
                sigmaImage.reset(new PagedImage<Float>(sigmaName));
            }
        }
        else if (
            sigma.type() == variant::DOUBLEVEC
            || sigma.type() == variant::INTVEC
        ) {
            sigmaArray.reset(new Array<Float>());
            vector<double> sigmaVector = sigma.getDoubleVec();
            Vector<Int> shape = sigma.arrayshape();
            sigmaArray->resize(IPosition(shape));
            convertArray(
                *sigmaArray,
                Vector<Double>(sigmaVector).reform(IPosition(shape))
            );
        }
        else if (sigma.type() == variant::BOOLVEC) {
            // nothing to do
        }
        else {
            _log << LogIO::SEVERE
                << "Unrecognized type for sigma. Use either a string (image name) or a numpy array"
                << LogIO::POST;
            return 0;
        }
        String myspxtype;
        vector<double> plpest, ltpest;
        vector<bool> plpfix, ltpfix;
        if (! spxtype.empty()) {
            myspxtype = String(spxtype);
            myspxtype.downcase();
            if (myspxtype == "plp") {
                plpest = spxest;
                plpfix = spxfix;
            }
            else if (myspxtype == "ltp") {
                ltpest = spxest;
                ltpfix = spxfix;
            }
            else {
                ThrowCc("Unsupported value for spxtype");
            }
        }
        SpectralList spectralList = SpectralListFactory::create(
            _log, pampest, pcenterest, pfwhmest, pfix, gmncomps,
            gmampcon, gmcentercon, gmfwhmcon, gmampest,
            gmcenterest, gmfwhmest, gmfix, pfunc, plpest, plpfix,
            ltpest, ltpfix
        );
        ThrowIf(
            ! estimates.empty() && spectralList.nelements() > 0,
            "You cannot specify both an "
            "estimates file and set estimates "
            "directly. You may only do one or "
            "the either (or neither in which "
            "case you must specify ngauss and/or poly)"
        );
        std::shared_ptr<ImageProfileFitter> fitter;
        if (spectralList.nelements() > 0) {
            fitter.reset(new ImageProfileFitter(
                _imageF, regionName, regionPtr.get(),
                box, chans, stokes, mask, axis,
                spectralList
            ));
        }
        else if (! estimates.empty()) {
            fitter.reset(new ImageProfileFitter(
                _imageF, regionName, regionPtr.get(),
                box, chans, stokes, mask, axis,
                estimates
            ));
        }
        else {
            fitter.reset(new ImageProfileFitter(
                _imageF, regionName, regionPtr.get(),
                box, chans, stokes, mask, axis,
                ngauss
            ));
        }
        fitter->setDoMultiFit(multifit);
        if (poly >= 0) {
            fitter->setPolyOrder(poly);
        }
        fitter->setModel(model);
        fitter->setResidual(residual);
        fitter->setAmpName(amp);
        fitter->setAmpErrName(amperr);
        fitter->setCenterName(center);
        fitter->setCenterErrName(centererr);
        fitter->setFWHMName(fwhm);
        fitter->setFWHMErrName(fwhmerr);
        fitter->setIntegralName(integral);
        fitter->setIntegralErrName(integralerr);
        fitter->setMinGoodPoints(minpts > 0 ? minpts : 0);
        fitter->setStretch(stretch);
        fitter->setLogResults(logResults);
        if (! planes.empty()) {
            std::set<int> myplanes(planes.begin(), planes.end());
            ThrowIf(*myplanes.begin() < 0, "All planes must be nonnegative");
            fitter->setGoodPlanes(std::set<uInt>(myplanes.begin(), myplanes.end()));
        }
        if (! logfile.empty()) {
            fitter->setLogfile(logfile);
            fitter->setLogfileAppend(append);
        }
        if (mygoodamps.size() == 2) {
            fitter->setGoodAmpRange(mygoodamps[0], mygoodamps[1]);
        }
        if (mygoodcenters.size() == 2) {
            fitter->setGoodCenterRange(mygoodcenters[0], mygoodcenters[1]);
        }
        if (mygoodfwhms.size() == 2) {
            fitter->setGoodFWHMRange(mygoodfwhms[0], mygoodfwhms[1]);
        }
        if (sigmaImage.get()) {
            fitter->setSigma(sigmaImage.get());
        }
        else if (sigmaArray.get()) {
            fitter->setSigma(*sigmaArray);
        }
        if (! outsigma.empty()) {
            if (sigmaImage.get() || sigmaArray.get()) {
                fitter->setOutputSigmaImage(outsigma);
            }
            else {
                _log << LogIO::WARN
                    << "outsigma specified but no sigma image "
                    << "or array specified. outsigma will be ignored"
                    << LogIO::POST;
            }
        }
        if (plpest.size() > 0 || ltpest.size() > 0) {
            variant::TYPE t = div.type();
            if (div.type() == variant::BOOLVEC) {
                fitter->setAbscissaDivisor(0);
            }
            else if (t == variant::INT || t == variant::DOUBLE) {
                fitter->setAbscissaDivisor(div.toDouble());
            }
            else if (t == variant::STRING || t == variant::RECORD) {
                fitter->setAbscissaDivisor(casaQuantity(div));
            }
            else {
                throw AipsError("Unsupported type " + div.typeString() + " for div");
            }
            if (! myspxtype.empty()) {
                if (myspxtype == "plp") {
                    fitter->setPLPName(spxsol);
                    fitter->setPLPErrName(spxerr);
                }
                else if (myspxtype == "ltp") {
                    fitter->setLTPName(spxsol);
                    fitter->setLTPErrName(spxerr);
                }
            }
        }
        return fromRecord(fitter->fit());
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::fromarray(
    const std::string& outfile, const variant& pixels,
    const record& csys, bool linear, bool overwrite,
    bool log, const std::string& type
) {
    try {
        _reset();
        auto mytuple = _fromarray(
            outfile, pixels, csys, linear,
            overwrite, log, type
        );
        _setImage(mytuple);
        vector<String> names {
            "pixels", "csys", "linear",
            "overwrite", "log", "type"
        };
        variant k("[...]");
        const auto* mpixels = pixels.size() <= 100 ? &pixels : &k;
        vector<variant> values {
            *mpixels, csys, linear, overwrite, log, type
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

ITUPLE image::_fromarray(
    const string& outfile,
    const ::casac::variant& pixels, const record& csys,
    bool linear, bool overwrite, bool log,
    const string& type
) {
    String mytype = type;
    mytype.downcase();
    ThrowIf(
        ! (mytype == "d" || mytype == "f"),
        "Unsupported value for type: \"" + type + "\". "
        "Supported values are \"d\" and \"f\""
    );
    auto doFloat = mytype == "f";
    ::casacore::IPosition shape = pixels.arrayshape();
    ThrowIf(
        shape.empty(), "The pixels array cannot be empty"
    );
    ::casacore::Array<::casacore::Float> floatArray;
    ::casacore::Array<::casacore::Double> doubleArray;
    ::casacore::Array<::casacore::Complex> complexArray;
    ::casacore::Array<::casacore::DComplex> dcomplexArray;
    auto pixType = pixels.type();
    if (
        pixType == ::casac::variant::DOUBLEVEC
        || pixType == ::casac::variant::INTVEC
        || pixType == ::casac::variant::LONGVEC
        || pixType == ::casac::variant::UINTVEC
    ) {
        ::casacore::Array<::casacore::Double> dv;
        if (pixType == ::casac::variant::DOUBLEVEC) {
            dv = ::casacore::Vector<::casacore::Double>(
                pixels.getDoubleVec()
            ).reform(shape);
        }
        else if (pixType == ::casac::variant::INTVEC) {
            dv = ::casacore::Vector<::casacore::Double>(
                pixels.getIntVec()
            ).reform(shape);
        }
        else if (pixType == ::casac::variant::LONGVEC) {
            dv = ::casacore::Vector<::casacore::Double>(
                pixels.getLongVec()
            ).reform(shape);
        }
        else if (pixType == ::casac::variant::UINTVEC) {
            dv = ::casacore::Vector<::casacore::Double>(
                pixels.getuIntVec()
            ).reform(shape);
        }
        else {
            ThrowCc("Logic error");
        }
        if (doFloat) {
            floatArray.resize(shape);
            ::casacore::convertArray(floatArray, dv);
        }
        else {
            doubleArray = dv;
        }
    }
    else if (pixels.type() == ::casac::variant::COMPLEXVEC) {
        auto localpix = ::casacore::Vector<::casacore::DComplex>(
            pixels.getComplexVec()
        ).reform(shape);
        if(doFloat) {
            complexArray.resize(shape);
            ::casacore::convertArray(complexArray, localpix);
        }
        else {
            dcomplexArray = localpix;
        }
    }
    else {
        ThrowCc("pixels is not understood, try using an array");
    }
    casacore::LogOrigin lor("image", __func__);
    _log << lor;
    std::unique_ptr<Record> coordinates(toRecord(csys));
    SPIIF f;
    SPIIC c;
    SPIID d;
    SPIIDC dc;
    if (! floatArray.empty()) {
        f = ImageFactory::imageFromArray(
            outfile, floatArray, *coordinates,
            linear, overwrite, log
        );
    }
    else if (! doubleArray.empty()) {
        d = ImageFactory::imageFromArray(
            outfile, doubleArray, *coordinates,
            linear, overwrite, log
        );
    }
    else if (! complexArray.empty()) {
        c = ImageFactory::imageFromArray(
            outfile, complexArray, *coordinates,
            linear, overwrite, log
        );
    }
    else {
        dc = ImageFactory::imageFromArray(
            outfile, dcomplexArray, *coordinates,
            linear, overwrite, log
        );
    }
    return ::casa::ITUPLE(f, c, d, dc);
}

bool image::fromascii(
    const string& outfile, const string& infile,
    const vector<int>& shape, const string& sep, const record& csys,
    bool linear, bool overwrite
) {
    try {
        _log << _ORIGIN;
        _log << LogIO::WARN << __func__ << "() IS DEPRECATED AND WILL BE "
            << "REMOVED IN A NEAR-FUTURE VERSION OF CASA. YOU SHOULD USE "
            << "ANOTHER SET OF IMAGE EXPORT AND IMPORT METHODS SUCH AS "
            << "tofits()/fromfits() TO EXPORT AND IMPORT CASA IMAGES. IF YOU "
            << "SIMPLY WISH TO MODIFY PIXEL VALUES, USE getchunk()/putchunk() "
            << "OR getregion()/putregion() FOR THAT" << LogIO::POST;
        ThrowIf(infile.empty(), "infile must be specified");
        ThrowIf(
            shape.size() == 1 && shape[0] == -1,
            "Image shape must be specified"
        );
        std::unique_ptr<Record> coordsys(toRecord(csys));
        _reset();
        _imageF = ImageFactory::fromASCII(
            outfile, infile, IPosition(Vector<Int>(shape)),
            sep, *coordsys, linear, overwrite
        );
        vector<String> names {
            "outfile", "infile", "shape", "sep",
            "csys", "linear",  "overwrite"
        };
        vector<variant> values {
            outfile, infile, shape, sep,
            csys, linear,  overwrite
        };
        this->_addHistory(__func__, names, values);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return true;
}

bool image::fromcomplist(
    const string& outfile, const vector<int>& shape, const variant& cl,
    const record& csys, bool overwrite, bool log, bool cache
) {
    try {
        _log << _ORIGIN;
        _reset();
        std::unique_ptr<Record> coordinates(toRecord(csys));
        auto myType = cl.type();
        std::unique_ptr<Record> mycl;
        if (myType == variant::RECORD) {
            std::unique_ptr<variant> clone(cl.clone());
            mycl.reset(toRecord(clone->asRecord()));
        }
        else if (myType == variant::STRING) {
            auto myname = cl.toString();
            ThrowIf(myname.empty(), "Component list table name cannot be empty");
            componentlist cltool;
            cltool.open(myname, True);
            std::unique_ptr<record> myrec(cltool.torecord());
            mycl.reset(toRecord(*myrec));
            cltool.done();
        }
        else {
            ThrowCc("Unsupported type for parameter cl");
        }
        _imageF = ImageFactory::createComponentListImage(
            outfile, *mycl, shape, *coordinates, overwrite, log, cache
        );
        vector<String> names {
            "outfile", "shape", "cl",
            "csys", "overwrite", "log", "cache"
        };
        vector<variant> values {
            outfile, shape, cl, csys, overwrite, log, cache
        };
        _addHistory(__func__, names, values);
        if (! outfile.empty()) {
            // force a flush to disk and reopen
            done();
            open(outfile, cache);
        }
        return True;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::fromfits(
    const string& outfile, const string& fitsfile,
    int whichrep, int whichhdu, bool zeroBlanks,
    bool overwrite
) {
    try {
        _log << _ORIGIN;
        auto im = ImageFactory::fromFITS(
            outfile, fitsfile, whichrep, whichhdu,
            zeroBlanks, overwrite
        );
        if (im) {
            _reset();
            _imageF = im;
            vector<String> names {
                "outfile", "fitsfile", "whichrep",
                "whichhdu", "zeroBlanks", "overwrite"
            };
            vector<variant> values {
                outfile, fitsfile, whichrep,
                whichhdu, zeroBlanks, overwrite
            };
            _addHistory(__func__, names, values);
            return true;
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::fromimage(
    const string& outfile, const string& infile,
    const variant& region, const variant& mask,
    bool dropdeg, bool overwrite
) {
    try {
        _log << _ORIGIN;
        String theMask = _getMask(mask);
        std::shared_ptr<Record> regionPtr(_getRegion(region, false));
        auto imagePtrs = ImageFactory::fromImage(
            outfile, infile, *regionPtr, theMask,
            dropdeg, overwrite
        );
        _setImage(imagePtrs);
        vector<String> names {
            "outfile", "infile", "region",
            "mask", "dropdeg", "overwrite"
        };
        vector<variant> values {
            outfile, infile, region, mask,
            dropdeg, overwrite
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::fromrecord(const record& imrecord, const string& outfile) {
    try {
        _log << _ORIGIN;
        std::unique_ptr<casacore::Record> tmpRecord(toRecord(imrecord));
        _reset();
        auto imagePair = ImageFactory::fromRecord(*tmpRecord, outfile);
        vector<String> names { "record", "outfile" };
        vector<variant> values { imrecord, outfile };
        auto msgs = _newHistory(__func__, names, values);
        if (imagePair.first) {
            _imageF = imagePair.first;
        }
        else {
            _imageC = imagePair.second;
        }
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        RETHROW(x);
    }
    return false;
}

bool image::fromshape(
    const string& outfile, const vector<int>& shape, const record& csys,
    const bool linear, const bool overwrite, const bool log, const string& type
) {
    try {
        LogOrigin lor("image", __func__);
        _log << lor;
        _reset();
        std::unique_ptr<Record> coordinates(toRecord(csys));
        String mytype = type;
        mytype.downcase();
        ThrowIf(
            ! (
                mytype == "f" || mytype == "c"
                || mytype == "d" || mytype == "cd"
            ),
            "Input parm type must be either 'f', 'c', 'd', or 'cd'"
        );
        if (mytype == "f") {
            _imageF = ImageFactory::floatImageFromShape(
                outfile, shape, *coordinates,
                linear, overwrite, log
            );
        }
        else if (mytype == "c") {
            _imageC = ImageFactory::complexImageFromShape(
                outfile, shape, *coordinates,
                linear, overwrite, log
            );
        }
        else if (mytype == "d") {
            _imageD = ImageFactory::doubleImageFromShape(
                outfile, shape, *coordinates,
                linear, overwrite, log
            );
        }
        else if (mytype == "cd") {
            _imageDC = ImageFactory::complexDoubleImageFromShape(
                outfile, shape, *coordinates,
                linear, overwrite, log
            );
        }
        vector<String> names {
            "outfile", "shape", "csys", "linear",
            "overwrite", "log", "type"
        };
        vector<variant> values {
            outfile, shape, csys, linear,
            overwrite, log, type
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

variant* image::getchunk(
    const std::vector<int>& blc, const std::vector<int>& trc,
    const std::vector<int>& inc, const std::vector<int>& axes,
    bool list, bool dropdeg, bool getmask
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        casacore::Record ret;
        if (_imageF) {
            ret = _getchunk<Float>(
                _imageF, blc, trc, inc,
                axes, list, dropdeg
            );
            if (! getmask) {
                Array<Float> vals = ret.asArrayFloat("values");
                vector<double> v(vals.begin(), vals.end());
                return new variant(v, vals.shape().asStdVector());
            }
        }
        else if (_imageC) {
            ret = _getchunk<Complex> (
                _imageC, blc, trc, inc,
                axes, list, dropdeg
            );
            if (! getmask) {
                Array<Complex> vals = ret.asArrayComplex("values");
                vector<std::complex<double> > v(vals.begin(), vals.end());
                return new variant(v, vals.shape().asStdVector());
            }
        }
        else if (_imageD) {
            ret = _getchunk<Double> (
                _imageD, blc, trc, inc,
                axes, list, dropdeg
            );
            if (! getmask) {
                Array<Double> vals = ret.asArrayDouble("values");
                vector<Double> v(vals.begin(), vals.end());
                return new variant(v, vals.shape().asStdVector());
            }
        }
        else if (_imageDC) {
            ret = _getchunk<DComplex> (
                _imageDC, blc, trc, inc,
                axes, list, dropdeg
            );
            if (! getmask) {
                Array<DComplex> vals = ret.asArrayDComplex("values");
                vector<DComplex> v(vals.begin(), vals.end());
                return new variant(v, vals.shape().asStdVector());
            }
        }
        else {
            ThrowCc("Logic Error");
        }
        if (getmask) {
            Array<Bool> pixelMask = ret.asArrayBool("mask");
            std::vector<bool> s_pixelmask(pixelMask.begin(), pixelMask.end());
            return new variant(s_pixelmask, pixelMask.shape().asStdVector());
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    // eliminate compiler warning, execution should never get here
    return nullptr;
}

template<class T> Record image::_getchunk(
    SPCIIT myimage,
    const vector<int>& blc, const vector<int>& trc,
    const vector<int>& inc, const vector<int>& axes,
    bool list, bool dropdeg
) {
    Array<T> pixels;
    Array<Bool> pixelMask;
    Vector<Int> iaxes(axes);
    // if default value change it to empty vector
    if (iaxes.size() == 1 && iaxes[0] < 0) {
        iaxes.resize();
    }
    uInt ndim = myimage->ndim();

    if (iaxes.size() == 1 && iaxes[0] < 0) {
        iaxes.resize();
    }
    // We have to support handling of sloppy inputs for backwards
    // compatibility. Ugh.
    vector<int> mblc(ndim);
    vector<int> mtrc(ndim);
    vector<int> minc(ndim);
    if (blc.size() == 1 && blc[0] < 0) {
        IPosition x(ndim, 0);
        mblc = x.asStdVector();
    }
    else {
        for (uInt i=0; i<ndim; i++) {
            mblc[i] = i < blc.size() ? blc[i] : 0;
        }
    }
    IPosition shape = myimage->shape();
    if (trc.size() == 1 && trc[0] < 0) {
        mtrc = (shape - 1).asStdVector();
    }
    else {
        for (uInt i=0; i<ndim; i++) {
            mtrc[i] = i < trc.size() ? trc[i] : shape[i] - 1;
        }
    }
    if (inc.size() == 1 && inc[0] == 1) {
        IPosition x(ndim, 1);
        minc = x.asStdVector();
    }
    else {
        for (uInt i=0; i<ndim; i++) {
            minc[i] = i < inc.size() ? inc[i] : 1;
        }
    }
    for (uInt i=0; i<ndim; i++) {
        if (mblc[i] < 0 || mblc[i] > shape[i] - 1) {
            mblc[i] = 0;
        }
        if (mtrc[i] < 0 || mtrc[i] > shape[i] - 1) {
            mtrc[i] = shape[i] - 1;
        }
        if (mblc[i] > mtrc[i]) {
            mblc[i] = 0;
            mtrc[i] = shape[i] - 1;
        }
        if (inc[i] > shape[i]) {
            minc[i] = 1;
        }
    }
    Vector<Double> vblc(mblc);
    Vector<Double> vtrc(mtrc);
    Vector<Double> vinc(minc);
    LCSlicer slicer(vblc, vtrc, vinc);
    Record rec;
    rec.assign(slicer.toRecord(""));
    PixelValueManipulator<T> pvm(myimage, &rec, "");
    if (axes.size() != 1 || axes[0] >= 0) {
        pvm.setAxes(IPosition(axes));
    }
    pvm.setVerbosity(
        list ? ImageTask<T>::DEAFENING : ImageTask<T>::QUIET
    );
    pvm.setDropDegen(dropdeg);
    return pvm.get();
}

record* image::getprofile(
    int axis, const string& function, const variant& region,
    const string& mask, const string& unit, bool stretch,
    const string& spectype, const variant& restfreq,
    const string& frame, const string& logfile
) {
    try {
        _log << _ORIGIN;
        ThrowIf(
            _detached(), "No image attached to tool"
        );
        _notSupported(__func__);
        ThrowIf(axis<0, "Axis must be greater than 0");
        std::shared_ptr<Record> myregion(_getRegion(region, false));
        std::shared_ptr<casacore::Quantity> rfreq;
        if (restfreq.type() != variant::BOOLVEC) {
            String rf = restfreq.toString();
            rf.trim();
            if (! rf.empty()) {
                rfreq.reset(
                    new casacore::Quantity(_casaQuantityFromVar(variant(restfreq)))
                );
            }
        }
        String regionName = region.type() == variant::STRING
            ? region.toString() : "";
        String myframe = frame;
        myframe.trim();
        if (_imageF) {
            SPCIIF myimage = _imageF;
            return fromRecord(
                _getprofile(
                    myimage, axis, function, unit,
                    *myregion, mask, stretch,
                    spectype, rfreq.get(), myframe,
                    logfile, regionName
                )
            );
        }
        else {
            SPCIIC myimage = _imageC;
            return fromRecord(
                _getprofile(
                    myimage, axis, function, unit,
                    *myregion, mask, stretch,
                    spectype, rfreq.get(), myframe,
                    logfile, regionName
                )
            );
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> Record image::_getprofile(
    SPCIIT myimage, int axis, const String& function,
    const String& unit, const Record& region, const String& mask,
    bool stretch, const String& spectype,
    const casacore::Quantity* const &restfreq, const String& frame,
    const String& logfile, const String& regionName
) {
    PixelValueManipulatorData::SpectralType type = PixelValueManipulatorData::spectralType(spectype);
    PixelValueManipulator<T> pvm(myimage, &region, mask);
    pvm.setLogfile(logfile);
    pvm.setRegionName(regionName);
    pvm.setStretch(stretch);
    Record x = pvm.getProfile(axis, function, unit, type, restfreq, frame);
    return x;
}

variant* image::getregion(
    const variant& region, const std::vector<int>& axes,
    const variant& mask, bool list, bool dropdeg,
    bool getmask, bool stretch
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        if (_imageF) {
            return _getregion2(
                _imageF, region, axes, mask, list,
                dropdeg, getmask, stretch
            );
        }
        else if (_imageC) {
            return _getregion2(
                _imageC, region, axes, mask, list,
                dropdeg, getmask, stretch
            );
        }
        else if (_imageD) {
            return _getregion2(
                _imageD, region, axes, mask, list,
                dropdeg, getmask, stretch
            );
        }
        else if (_imageDC) {
            return _getregion2(
                _imageDC, region, axes, mask, list,
                dropdeg, getmask, stretch
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> variant* image::_getregion2(
    SPIIT image, const variant& region,
    const std::vector<int>& axes, const variant& mask,
    bool list, bool dropdeg, bool getmask, bool stretch
) {
    auto Region = _getRegion(region, false);
    auto Mask = _getMask(mask);
    Vector<Int> iaxes(axes);
    // if default value change it to empty vector
    if (iaxes.size() == 1 && iaxes[0] < 0) {
        iaxes.resize();
    }
    PixelValueManipulator<T> pvm(
        image, Region.get(), Mask
    );
    pvm.setAxes(IPosition(iaxes));
    pvm.setVerbosity(
        list ? ImageTask<T>::DEAFENING : ImageTask<T>::QUIET
    );
    pvm.setDropDegen(dropdeg);
    pvm.setStretch(stretch);
    auto ret = pvm.get();
    auto pixelmask = ret.asArrayBool("mask");
    auto s_shape = pixelmask.shape().asStdVector();
    if (getmask) {
        pixelmask.shape().asVector().tovector(s_shape);
        std::vector<bool> s_pixelmask(pixelmask.begin(), pixelmask.end());
        return new ::casac::variant(s_pixelmask, s_shape);
    }
    if (_imageF || _imageD) {
        std::vector<Double> d_pixels;
        if (_imageF) {
            auto pixels = ret.asArrayFloat("values");
            d_pixels = std::vector<Double>(pixels.begin(), pixels.end());
        }
        else {
            d_pixels = ret.asArrayDouble("values").tovector();
        }
        return new ::casac::variant(d_pixels, s_shape);
    }
    else if (_imageC || _imageDC) {
        std::vector<std::complex<double> > d_pixels;
        if (_imageC) {
            auto pixels = ret.asArrayComplex("values");
            d_pixels = std::vector<std::complex<double>>(pixels.begin(), pixels.end());
        }
        else {
            d_pixels = ret.asArrayDComplex("values").tovector();
        }
        return new ::casac::variant(d_pixels, s_shape);
    }
    else {
        ThrowCc("Logic error");
    }
}

record* image::getslice(
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<int>& axes, const std::vector<int>& coord,
    int npts, const std::string& method
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        Vector<Int> ncoord(coord);
        if (ncoord.size() == 1 && ncoord[0] == -1) {
            ncoord.resize(shape().size());
            ncoord.set(0);
        }
        unique_ptr<Record> outRec;
        if (_imageF) {
            outRec.reset(
                PixelValueManipulator<Float>::getSlice(
                    _imageF, Vector<Double>(x), Vector<Double>(y),
                    Vector<Int>(axes), ncoord, npts, method
                )
            );
        }
        else if (_imageC) {
            outRec.reset(
                PixelValueManipulator<Complex>::getSlice(
                    _imageC, Vector<Double>(x), Vector<Double>(y),
                    Vector<Int>(axes), ncoord, npts, method
                )
            );
        }
        else if (_imageD) {
            outRec.reset(
                PixelValueManipulator<Double>::getSlice(
                    _imageD, Vector<Double>(x), Vector<Double>(y),
                    Vector<Int>(axes), ncoord, npts, method
                )
            );
        }
        else if (_imageDC) {
            outRec.reset(
                PixelValueManipulator<DComplex>::getSlice(
                    _imageDC, Vector<Double>(x), Vector<Double>(y),
                    Vector<Int>(axes), ncoord, npts, method
                )
            );
        }
        else {
            ThrowCc("Logic error");
        }
        return fromRecord(*outRec);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::hanning(
    const string& outfile, const variant& region,
    const variant& vmask, int axis, bool drop,
    bool overwrite, bool stretch,
    const string& dmethod
) {
    LogOrigin lor(_class, __func__);
    _log << lor;
    if (_detached()) {
        throw AipsError("Unable to create image");
    }
    try {
        _notSupported(__func__);
        auto myregion = _getRegion(
            region, true
        );
        auto mask = _getMask(vmask);
        if (axis < 0) {
            const CoordinateSystem csys = _imageF
                ? _imageF->coordinates()
                : _imageC->coordinates();
            ThrowIf(
                ! csys.hasSpectralAxis(),
                "Axis not specified and image has no spectral coordinate"
            );
            axis = csys.spectralAxisNumber(false);
        }
        ImageDecimatorData::Function dFunction = ImageDecimatorData::NFUNCS;
        if (drop) {
            String mymethod = dmethod;
            mymethod.downcase();
            if (mymethod.startsWith("m")) {
                dFunction = ImageDecimatorData::MEAN;
            }
            else if (mymethod.startsWith("c")) {
                dFunction = ImageDecimatorData::COPY;
            }
            else {
                ThrowCc(
                    "Value of dmethod must be "
                    "either 'm'(ean) or 'c'(opy)"
                );
            }
        }
        vector<variant> values { outfile, region, vmask, axis, drop, overwrite, stretch, dmethod };
        if (_imageF) {
            SPCIIF image = _imageF;
            return _hanning(
                image, myregion, mask, outfile,
                overwrite, stretch, axis, drop,
                dFunction, values
            );
        }
        else {
            SPCIIC image = _imageC;
            return _hanning(
                image, myregion, mask, outfile,
                overwrite, stretch, axis, drop,
                dFunction, values
            );
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
}

template <class T> image* image::_hanning(
    SPCIIT myimage, std::shared_ptr<const Record> region,
    const String& mask, const string& outfile, bool overwrite,
    bool stretch, int axis, bool drop,
    ImageDecimatorData::Function dFunction,
    const vector<variant> values
) const {
    ImageHanningSmoother<T> smoother(
        myimage, region.get(), mask, outfile, overwrite
    );
    smoother.setAxis(axis);
    smoother.setDecimate(drop);
    smoother.setStretch(stretch);
    if (drop) {
        smoother.setDecimationFunction(dFunction);
    }
    if (_doHistory) {
        vector<String> names {
            "outfile", "region", "mask", "axis",
            "drop", "overwrite", "stretch", "dmethod"
        };
        auto msgs = _newHistory("hanning", names, values);
        smoother.addHistory(_ORIGIN, msgs);
    }
    return new image(smoother.smooth());
}

vector<bool> image::haslock() {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return vector<bool>();
        }
        _notSupported(__func__);
        if (_imageF) {
            return vector<bool> {
                _imageF->hasLock(FileLocker::Read),
                _imageF->hasLock(FileLocker::Write)
            };
        }
        else {
            return vector<bool> {
                _imageC->hasLock(FileLocker::Read),
                _imageC->hasLock(FileLocker::Write)
            };
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return vector<bool>();
}

record* image::histograms(
    const vector<int>& axes, const variant& region, const variant& mask,
    int nbins, const vector<double>& includepix, bool cumu, bool log,
    bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        return nullptr;
    }
    try {
        ThrowIf(! (
            _imageF || _imageD),
            "This method only supports real-valued images"
        );
        if (_imageF) {
            return _histograms(
                _imageF, axes, region, mask, nbins,
                includepix, cumu, log, stretch
            );
        }
        else if (_imageD) {
            return _histograms(
                _imageD, axes, region, mask, nbins,
                includepix, cumu, log, stretch
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_histograms(
    SPIIT myImage, const vector<int>& axes, const variant& region,
    const variant& mask, int nbins, const vector<double>& includepix, bool cumu,
    bool log, bool stretch
) {
    vector<uInt> myaxes;
    if (axes.size() != 1 || axes[0] != -1) {
        ThrowIf(
            *min_element(axes.begin(), axes.end()) < 0,
            "All axes must be nonnegative"
        );
        myaxes.insert(begin(myaxes), begin(axes), end(axes));
    }
    auto regionRec = _getRegion(region, false);
    String Mask = _getMask(mask);
    vector<Double> myIncludePix;
    if (!(includepix.size() == 1 && includepix[0] == -1)) {
        myIncludePix = includepix;
    }
    ImageHistogramsCalculator<T> ihc(
        myImage, regionRec.get(), Mask
    );
    if (! myaxes.empty()) {
        ihc.setAxes(myaxes);
    }
    ihc.setNBins(nbins);
    if (! myIncludePix.empty()) {
        ihc.setIncludeRange(myIncludePix);
    }
    ihc.setCumulative(cumu);
    ihc.setDoLog10(log);
    ihc.setStretch(stretch);
    return fromRecord(ihc.compute());
}

std::vector<std::string> image::history(bool list) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return vector<string>();
        }
        if (_imageF) {
            ImageHistory<Float> hist(_imageF);
            return fromVectorString(hist.get(list));
        }
        else if (_imageC) {
            ImageHistory<Complex> hist(_imageC);
            return fromVectorString(hist.get(list));
        }
        else if (_imageD) {
            ImageHistory<Double> hist(_imageD);
            return fromVectorString(hist.get(list));
        }
        else if (_imageDC) {
            ImageHistory<DComplex> hist(_imageDC);
            return fromVectorString(hist.get(list));
        }
        else {
            ThrowCc("Logic error");
        }
    } catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return vector<string>();
}

image* image::imagecalc(
    const string& outfile, const string& pixels,
    bool overwrite, const string& imagemd, const string& prec 
) {
    try {
        ThrowIf(
            pixels.empty(),
            "You must provide an expression using the pixels parameter"
        );
        String myPrec = prec;
        myPrec.downcase();
        auto asFloat = myPrec.startsWith("f");
        ThrowIf(
            ! myPrec.startsWith("d") && ! asFloat,
            "Unsupported value for type, it must be 'float' or 'double'"
        );
        if (isReal(ImageExprParse::command(pixels).dataType())) {
            return asFloat
                ? new image(
                    _imagecalc<Float>(outfile, pixels, overwrite, imagemd)
                )
                : new image(
                    _imagecalc<Double>(outfile, pixels, overwrite, imagemd)
                );
        }
        else {
            return asFloat
                ? new image(
                    _imagecalc<Complex>(outfile, pixels, overwrite, imagemd)
                )
                : new image(
                    _imagecalc<DComplex>(outfile, pixels, overwrite, imagemd)
                );
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> SPIIT image::_imagecalc(
    const string& outfile, const string& pixels,
    bool overwrite, const string& imagemd
) {
    ImageExprCalculator<T> calculator(pixels, outfile, overwrite);
    calculator.setCopyMetaDataFromImage(imagemd);
    auto out = calculator.compute();
    if (_doHistory) {
        vector<String> names {"outfile", "pixels", "overwrite", "imagemd"};
        vector<variant> values {outfile, pixels, overwrite, imagemd};
        _addHistory(out, "imagecalc", names, values);
    }
    return out;
}

image* image::imageconcat(
    const string& outfile, const variant& infiles, int axis,
    bool relax, bool tempclose, bool overwrite, bool reorder
) {
    try {
        Vector<String> inFiles;
        if (infiles.type() == variant::BOOLVEC) {
            inFiles.resize(0); // unset
        }
        else if (infiles.type() == variant::STRING) {
            sepCommaEmptyToVectorStrings(inFiles, infiles.toString());
        }
        else if (infiles.type() == variant::STRINGVEC) {
            inFiles = toVectorString(infiles.toStringVec());
        }
        else {
            ThrowCc("Unrecognized infiles datatype");
        }
        auto imageNames = Directory::shellExpand(inFiles, false).tovector();
        ThrowIf(
            imageNames.size() < 2,
            "You must provide at least two images to concatentate"
        );
        auto first = imageNames[0];
        imageNames.erase(imageNames.begin());
        std::shared_ptr<LatticeBase> latt(ImageOpener::openImage(first));
        ThrowIf (! latt, "Unable to open image " + first);
        auto dataType = latt->dataType();
        if (dataType == TpFloat) {
            return new image(
                _concat<Float>(
                    latt, outfile, infiles, axis, relax, tempclose,
                    overwrite, reorder, imageNames
                )
            );
        }
        else if (dataType == TpComplex) {
            return new image(
                _concat<Complex>(
                    latt, outfile, infiles, axis, relax, tempclose,
                    overwrite, reorder, imageNames
                )
            );
        }
        else if (dataType == TpDouble) {
            return new image(
                _concat<Double>(
                    latt, outfile, infiles, axis, relax, tempclose,
                    overwrite, reorder, imageNames
                )
            );
        }
        else if (dataType == TpDComplex) {
            return new image(
                _concat<DComplex>(
                    latt, outfile, infiles, axis, relax, tempclose,
                    overwrite, reorder, imageNames
                )
            );
        }
        else {
            ostringstream x;
            x << dataType;
            ThrowCc("Unsupported data type " + x.str());
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> SPIIT image::_concat(
    std::shared_ptr<LatticeBase> latt, const string& outfile,
    const variant& infiles, int axis, bool relax, bool tempclose,
    bool overwrite, bool reorder, const vector<String>& imageNames
) {
    SPIIT im = std::dynamic_pointer_cast<ImageInterface<T>>(latt);
    ThrowIf(! im, "dynamic cast failed");
    ImageConcatenator<T> concat(im, outfile, overwrite);
    concat.setAxis(axis);
    concat.setRelax(relax);
    concat.setReorder(reorder);
    concat.setTempClose(tempclose);
    if (_doHistory) {
        vector<String> names {
            "outfile", "infiles", "axis", "relax", "tempclose",
            "overwrite", "reorder"
        };
        vector<variant> values {
            outfile, infiles, axis,  relax, tempclose, overwrite, reorder
        };
        concat.addHistory(_ORIGIN, "ia.imageconcat", names, values);
    }
    return concat.concatenate(imageNames);
}

bool image::insert(
    const std::string& infile, const variant& region,
    const std::vector<double>& locate, bool verbose
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        Vector<Double> locatePixel(locate);
        if (locatePixel.size() == 1 && locatePixel[0] < 0) {
            locatePixel.resize(0);
        }
        auto Region = _getRegion(region, false);
        SPCIIF imageF;
        SPCIIC imageC;
        std::tie(imageF, imageC, std::ignore, std::ignore)
            = ImageFactory::fromFile(infile);
        ThrowIf(! (imageF || imageC), "Unsupported image data type");
        if (imageF && _imageF) {
            PixelValueManipulator<Float>::insert(
                *_imageF, *imageF, *Region,
                locatePixel, verbose
            );
        }
        else if (imageC && _imageC){
            PixelValueManipulator<Complex>::insert(
                *_imageC, *imageC, *Region,
                locatePixel, verbose
            );
        }
        else {
            ThrowCc("Attached image pixel data type differs from that of " + infile);
        }
        vector<String> names = {
           "infile", "region", "locate", "verbose"
        };
        vector<variant> values = {
            infile, region, locate, verbose
        };
        _addHistory(__func__, names, values);
        _statsF.reset();
        _statsD.reset();
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::isopen() {
    try {
        _log << _ORIGIN;
        return _imageF || _imageC || _imageD || _imageDC;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::ispersistent() {
    try {
        _log << LogOrigin("image", "ispersistent");
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        if (_imageF) {
            return _imageF->isPersistent();
        }
        else {
            return _imageC->isPersistent();
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::lock(bool writelock, int nattempts) {
    try {
        _log << LogOrigin("image", __func__);
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        FileLocker::LockType locker = FileLocker::Read;
        if (writelock) {
            locker = FileLocker::Write;
        }
        uInt n = max(0, nattempts);
        if (_imageF) {
            return _imageF->lock(locker, n);
        }
        else {
            return _imageC->lock(locker, n);
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::makecomplex(
    const std::string& outfile, const std::string& imagFile,
    const variant& region, bool overwrite
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        ThrowIf(
            ! (_imageF || _imageD), "The attached image must be float valued"
        );
        std::shared_ptr<Record> Region(_getRegion(region, false));
        auto imagePtrs = ImageFactory::fromFile(imagFile);
        auto imageF = std::get<0>(imagePtrs);
        auto imageD = std::get<2>(imagePtrs);
        ThrowIf(
            ! (imageF || imageD),
            imagFile + " does not have supported real valued pixels"
        );
        ThrowIf(
            (_imageF && imageD) || (_imageD && imageF),
            "Real and imaginary images do not have the same precision"
        );
        SPIIC cImage;
        SPIIDC dcImage;
        if (_imageF) {
            cImage = ImageFactory::makeComplex<Float>(
                _imageF, imageF, outfile, *Region, overwrite
            );
        }
        else if (_imageD) {
            dcImage = ImageFactory::makeComplex<Double>(
                _imageD, imageD, outfile, *Region, overwrite
            );
        }
        else {
            ThrowCc("Logic error");
        }
        vector<String> names = {
            "outfile", "imag", "region", "overwrite"
        };
        vector<variant> values = {
            outfile, imagFile, region, overwrite
        };
        if (cImage) {
            _addHistory(cImage, __func__, names, values);
        }
        else if (dcImage) {
            _addHistory(dcImage, __func__, names, values);
        }
        else {
            ThrowCc("Logic Error");
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

image* image::deviation(
    const std::string& outfile, const variant& region,
    const string& mask, bool overwrite, bool stretch,
    const vector<int>& grid, const variant& anchor,
    const variant& xlength, const variant& ylength,
    const string& interp, const string& stattype,
    const string& statalg, double zscore, int maxiter
) {
    _log << _ORIGIN;
    try {
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        ThrowIf(
            grid.size() != 2,
            "grid must have exactly two positive integer values"
        );
        auto useRef = False;
        switch (anchor.type()) {
        case variant::INTVEC:
            ThrowIf(
                anchor.toIntVec().size() != 2,
                "anchor must have exactly two integer values"
            );
            useRef = False;
            break;
        case variant::STRING:
            ThrowIf(
                anchor.toString() != "ref",
                "Unsupported value for anchor: " + anchor.toString()
            );
            useRef = True;
            break;
        case variant::BOOLVEC:
            // because the interface always passes in a boolvec by default for a variant,
            // even if specified differently in the XML
            useRef = True;
            break;
        default:
            ThrowCc("Unsupported type for anchor");
        }
        ThrowIf(
            grid[0] <= 0 || grid[1] <= 0,
            "Both grid value(s) must be positive"
        );
        String mystatalg = statalg;
        auto myreg = _getRegion(region, False);
        auto  myxlen = xlength.type() == variant::INT
            ? casacore::String::toString(xlength.toInt()) + "pix"
            : xlength.toString();
        auto ytype = ylength.type();
        auto myylen = ytype == variant::BOOLVEC
            ? "" : ytype == variant::INT
            ? casacore::String::toString(ylength.toInt()) + "pix"
            : ylength.toString();
        String err;
        QuantumHolder qh;
        casacore::Quantity qxl, qyl;
        ThrowIf(
            ! qh.fromString(err, myxlen),
            "xlength is not a valid quantity: " + err
        );
        qxl = qh.asQuantity();
        if (myylen.empty()) {
            // circle, so we need the radius, not the diameter
            auto z = qh.asQuantity();
            qxl = z/2;
        }
        else {
            ThrowIf(
                ! qh.fromString(err, myylen),
                "ylength is not a valid quantity: " + err
            );
            qyl = qh.asQuantity();
        }
        String myinterp = interp;
        myinterp.downcase();
        Interpolate2D::Method interpAlg;
        if (myinterp.startsWith("c")) {
            interpAlg = Interpolate2D::CUBIC;
        }
        else if (myinterp.startsWith("la")) {
            interpAlg = Interpolate2D::LANCZOS;
        }
        else if (myinterp.startsWith("li")) {
            interpAlg = Interpolate2D::LINEAR;
        }
        else if (myinterp.startsWith("n")) {
            interpAlg = Interpolate2D::NEAREST;
        }
        else {
            ThrowCc("Interpolation algorithm " + interp + " is not supported.");
        }
        StatImageCreator sic(_imageF, myreg.get(), mask, outfile, overwrite);
        mystatalg.downcase();
        if (mystatalg.startsWith("cl")) {
            sic.configureClassical(ImageStatsData::AUTO);
        }
        else if (mystatalg.startsWith("ch")) {
            sic.configureChauvenet(zscore, maxiter);
        }
        else {
            ThrowCc("Unsupported stats algorithm " + statalg);
        }
        if (useRef) {
            sic.useReferencePixelAsAnchor();
        }
        else {
            auto myan = anchor.toIntVec();
            sic.setAnchorPosition(myan[0], myan[1]);
        }
        sic.setGridSpacing(grid[0], grid[1]);
        sic.setStretch(stretch);
        sic.setStatType(stattype);
        if (myylen.empty()) {
            sic.setRadius(qxl);
        }
        else {
            sic.setRectangle(qxl, qyl);
        }
        sic.setInterpAlgorithm(interpAlg);
        vector<String> names {
            "outfile", "region", "mask", "overwrite", "stretch",
            "grid", "anchor", "xlength", "ylength", "interp",
            "stattype", "statalg", "zscore", "maxiter"
        };
        vector<variant> values {
            outfile, region, mask, overwrite, stretch,
            grid, anchor, xlength, ylength, interp,
            stattype, statalg, zscore, maxiter
        };
        if (_doHistory) {
            auto msgs = _newHistory(__func__,names, values);
            sic.addHistory(_ORIGIN, msgs);
        }
        return new image(sic.compute());
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
             << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::maketestimage(
    const string& outfile, bool overwrite
) {
    try {
        _reset();
        _log << _ORIGIN;
        _imageF = ImageFactory::testImage(
            outfile, overwrite
        );
        vector<String> names = {
            "outfile", "overwrite"
        };
        vector<variant> values { outfile, overwrite };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

vector<string> image::maskhandler(
    const string& op, const vector<string>& name
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return vector<string>(0);
        }
        String oper = op;
        oper.upcase();
        vector<string> res;
        if (_imageF) {
            res = _handleMask(_imageF, oper, name);
        }
        else if (_imageC) {
            res = _handleMask(_imageC, oper, name);
        }
        else if (_imageD) {
            res = _handleMask(_imageD, oper, name);
        }
        else if (_imageDC) {
            res = _handleMask(_imageDC, oper, name);
        }
        else {
            ThrowCc("Logic error");
        }
        if (res.empty()) {
            res = vector<string>(1, "T");
        }
        if (
            oper.startsWith("SET") || oper.startsWith("DEL")
            || oper.startsWith("REN") || oper.startsWith("COPY")
        ) {
            vector<String> names {"op", "name"};
            vector<variant> values {op, name};
            _addHistory(__func__, names, values);
            _statsF.reset();
            _statsD.reset();
        }
        return res;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return vector<string>();
}

template<class T> vector<string>  image::_handleMask(
    SPIIT myimage, const String& oper,
    const vector<string>& name
) {
    ImageMaskHandler<T> imh(myimage);
    if (oper.startsWith("SET")) {
        auto myname = name.empty() ? "" : name[0];
        imh.set(myname);
        return vector<string>();
    }
    else if (oper.startsWith("DEF")) {
        return vector<string>(1, imh.defaultMask());
    }
    else if (oper.startsWith("DEL")) {
        imh.deleteMasks(std::set<casacore::String>(name.begin(), name.end()));
        return vector<string>();
    }
    else if (oper.startsWith("REN")) {
        ThrowIf(
            name.size() != 2,
            "name must be an array of size exactly two. "
            + String::toString(name.size()) + " values were "
            "given"
        );
        imh.rename(name[0], name[1]);
        return vector<string>();
    }
    else if (oper.startsWith("GET")) {
        return fromVectorString(imh.get());
    }
    else if (oper.startsWith("COP")) {
        imh.copy(name[0], name[1]);
        return vector<string>();
    }
    else {
        ThrowCc("Unknown operation " + oper);
    }
}

record* image::maxfit(
    const variant& region, bool doPoint,
    int width, bool absFind, bool list
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! _imageF,
            "This method only supports float-valued images"
        );
        auto Region = _getRegion(region, false);
        ImageMaxFitter<Float> imf(_imageF, Region.get());
        return fromRecord(imf.fit(doPoint, width, absFind, list));
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::miscinfo() {
    try {
        _log << LogOrigin("image", "miscinfo");
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        if (_imageF) {
            return fromRecord(_imageF->miscInfo());
        }
        else {
            return fromRecord(_imageC->miscInfo());
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::modify(
    const record& model, const variant& region,
    const variant& vmask, bool subtract, bool list,
    bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        return false;
    }
    try {
        ThrowIf(
            ! _imageF,
            "This method only supports float valued images"
        );
        String error;
        std::unique_ptr<Record> mymodel(toRecord(model));
        ComponentList cl;
        ThrowIf(
            ! cl.fromRecord(error, *mymodel),
            "model is an invalid componentlist record"
        );
        std::shared_ptr<Record> Region = _getRegion(region, false);
        String mask = _getMask(vmask);
        ComponentImager ci(
            _imageF, Region.get(), mask
        );
        ci.setComponentList(cl);
        ci.setSubtract(subtract);
        ci.setStretch(stretch);
        ci.modify(list);
        _statsF.reset();
        _statsD.reset();
        vector<String> names {
            "model", "region", "mask",
            "subtract", "list", "stretch"
        };
        vector<variant> values {
            model, region, vmask,
            subtract, list, stretch
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

image* image::moments(
    const vector<int>& moments, int axis,
    const variant& region, const variant& vmask,
    const vector<string>& in_method,
    const vector<int>& smoothaxes,
    const variant& smoothtypes,
    const vector<double>& smoothwidths,
    const vector<double>& d_includepix,
    const vector<double>& d_excludepix, double peaksnr,
    double stddev, const string& velocityType,
    const string& out, const string& smoothout,
    bool overwrite, bool removeAxis,
    bool stretch
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        UnitMap::putUser("pix", UnitVal(1.0), "pixel units");
        Vector<Int> whichmoments(moments);
        std::shared_ptr<Record> Region(_getRegion(region, false));
        auto mask = _getMask(vmask);
        Vector<String> kernels;
        if (smoothtypes.type() == ::casac::variant::BOOLVEC) {
            kernels.resize(0); // unset
        }
        else if (smoothtypes.type() == ::casac::variant::STRING) {
            sepCommaEmptyToVectorStrings(kernels, smoothtypes.toString());
        }
        else if (smoothtypes.type() == ::casac::variant::STRINGVEC) {
            kernels = toVectorString(smoothtypes.toStringVec());
        }
        else {
            _log << LogIO::WARN << "Unrecognized smoothtypes datatype"
                << LogIO::POST;
        }
        int num = kernels.size();
        vector<casacore::Quantity> kernelwidths(num);
        Unit u("pix");
        for (int i = 0; i < num; ++i) {
            kernelwidths[i] = casacore::Quantity(smoothwidths[i], u);
        }
        std::vector<Double> includepix;
        num = d_includepix.size();
        if (!(num == 1 && d_includepix[0] == -1)) {
            includepix = d_includepix;;
        }
        std::vector<Double> excludepix;
        num = d_excludepix.size();
        if (!(num == 1 && d_excludepix[0] == -1)) {
            excludepix = d_excludepix;
        }
        ThrowIf(
            ! includepix.empty() && ! excludepix.empty(),
            "Only one of includepix or excludepix may be specified, not both"
        );
        ImageMomentsTask<Float> momentsTask(
            _imageF, Region.get(), mask,
            smoothout, overwrite
        );
        momentsTask.setMoments(whichmoments);
        momentsTask.setAxis(axis);
        auto methods = toVectorString(in_method).tovector();
        momentsTask.setMethods(methods);
        if (
            ! smoothaxes.empty()
            && ! (smoothaxes.size() == 1 && smoothaxes[0] == -1)
        ) {
            ThrowIf (
                *std::min_element(smoothaxes.begin(), smoothaxes.end()) < 0,
                "All smoothaxes must be nonnegative"
            );
            std::vector<uInt> sa;
            for (const auto s : smoothaxes) {
                sa.push_back(s);
            }
            momentsTask.setSmoothAxes(sa);
        }
        if (! kernels.empty()) {
            momentsTask.setKernels(kernels.tovector());
        }
        if (! kernelwidths.empty()) {
            momentsTask.setKernelWidths(kernelwidths);
        }
        if (! includepix.empty() || ! excludepix.empty()) {
            auto vrange = ! includepix.empty() ? includepix : excludepix;
            std::vector<Float> range;
            for (const auto v : vrange) {
                range.push_back(v);
            }
            auto isInclude = ! includepix.empty();
            momentsTask.setIncludeExcludeRange(range, isInclude);
        }
        momentsTask.setSNR(peaksnr);
        momentsTask.setStdDev(stddev);
        momentsTask.setVelocityType(velocityType);
        momentsTask.setMomentImageName(out);
        momentsTask.setRemoveAxis(removeAxis);
        momentsTask.setStretch(stretch);
        vector<String> names {
            "moments", "axis", "region", "mask",
            "method", "smoothaxes",
            "smoothtypes", "smoothwidths",
            "includepix", "excludepix",
            "peaksnr", "stddev", "doppler",
            "outfile", "smoothout", "overwrite",
            "drop", "stretch"
        };
        vector<variant> values {
            moments, axis, region, vmask,
            in_method, smoothaxes,
            smoothtypes, smoothwidths,
            d_includepix, d_excludepix,
            peaksnr, stddev, velocityType,
            out, smoothout, overwrite,
            removeAxis, stretch
        };
        if (_doHistory) {
            auto msgs = _newHistory(__func__,names, values);
            momentsTask.addHistory(_ORIGIN, msgs);
        }
        return new image(momentsTask.makeMoments());
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                        << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

string image::name(bool strippath) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return "none";
        }
        _notSupported(__func__);
        return _name(strippath);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return "";
}

String image::_name(bool strippath) const {
    if (_imageF) {
        return _imageF->name(strippath);
    }
    else if (_imageC) {
        return _imageC->name(strippath);
    }
    else if (_imageD) {
        return _imageD->name(strippath);
    }
    else if (_imageDC) {
        return _imageDC->name(strippath);
    }
    else {
        ThrowCc("Logic error");
    }
}

image* image::newimage(const string& fileName) {
    try {
        _log << _ORIGIN;
        auto *rstat = newimagefromfile(fileName);
        ThrowIf(! rstat, "Unable to create image");
        return rstat;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::newimagefromarray(
    const string& outfile, const variant& pixels,
    const record& csys, bool linear,
    bool overwrite, bool log, const string& type
) {
    try {
        auto mytuple = _fromarray(
            outfile, pixels, csys,
            linear, overwrite, log, type
        );
        auto* res = new image(mytuple);
        vector<String> names = {
            "outfile", "pixels", "csys",
            "linear", "overwrite", "log"
        };
        variant k("[...]");
        const auto* mpixels = pixels.size() <= 100 ? &pixels : &k;
        vector<variant> values = {
            outfile, *mpixels, csys,
            linear, overwrite, log
        };
        if (_doHistory) {
            res->_addHistory(__func__, names, values);
        }
        return res;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::newimagefromfile(const string& fileName) {
    try {
        _log << _ORIGIN;
        // not adding history because all this method does is open
        // the image, it doesn't change it
        return new image(ImageFactory::fromFile(fileName));
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::newimagefromimage(
    const string& infile, const string& outfile,
    const variant& region, const variant& vmask,
    bool dropdeg, bool overwrite
) {
    try {
        _log << _ORIGIN;
        auto mask = this->_getMask(vmask);
        auto regionPtr = _getRegion(region, false, infile);
        auto ret = ImageFactory::fromImage(
            outfile, infile, *regionPtr, mask,
            dropdeg, overwrite
        );
        vector<String> names = {
            "infile", "outfile", "region",
            "vmask", "dropdeg", "overwrite"
        };
        vector<variant> values = {
            infile, outfile, region,
            vmask, dropdeg, overwrite
        };
        auto f = std::get<0>(ret);
        auto c = std::get<1>(ret);
        auto d = std::get<2>(ret);
        auto dc = std::get<3>(ret);
        if (f) {
            _addHistory(f, __func__, names, values);
        }
        else if (c) {
            _addHistory(c, __func__, names, values);
        }
        else if (d) {
            _addHistory(d, __func__, names, values);
        }
        else if (dc) {
            _addHistory(dc, __func__, names, values);
        }
        return new image(ret);
        ThrowCc("Error creating image");
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::newimagefromshape(
    const string& outfile, const vector<int>& shape, const record& csys,
    bool linear, bool overwrite, bool log, const string& type
) {
    try {
        _log << _ORIGIN;
        unique_ptr<image> ret(new image());
        ret->dohistory(False);
        ThrowIf(
            ! ret->fromshape(
                outfile, shape, csys, linear, overwrite, log, type
            ), "Failed to create image from shape"
        );
        vector<String> names {
            "outfile", "shape", "csys", "linear",
            "overwrite", "log", "type"
        };
        vector<variant> values {
            outfile, shape, csys, linear,
            overwrite, log, type
        };
        ret->dohistory(True);
        if (_doHistory) {
            ret->_addHistory(__func__, names, values);
        }
        return ret.release();
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::open(const std::string& infile, bool cache) {
    try {
        _log << _ORIGIN;
        if (_imageF || _imageC || _imageD || _imageDC) {
            _log << LogIO::WARN << "Another image is already open, closing first"
                << LogIO::POST;
        }
        _reset();
        std::tie(_imageF, _imageC, _imageD, _imageDC)
            = ImageFactory::fromFile(infile, cache);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

image* image::pad(
    const string& outfile, int npixels, double value, bool padmask,
    bool overwrite, const variant& region, const string& box,
    const string& chans, const string& stokes, const string& mask,
    bool  stretch, bool wantreturn
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        if (npixels <= 0) {
            _log << "Value of npixels must be greater than zero" << LogIO::EXCEPTION;
        }
        auto regionPtr = _getRegion(region, true);
        ImagePadder padder(
            _imageF, regionPtr.get(), box,
            chans, stokes, mask, outfile, overwrite
        );
        padder.setStretch(stretch);
        padder.setPaddingPixels(npixels, value, padmask);
        vector<String> names {
            "outfile", "npixels", "value", "padmask",
            "overwrite", "region", "box", "chans",
            "stokes", "mask", "stretch", "wantreturn"
        };
        vector<variant> values {
            outfile, npixels, value, padmask,
            overwrite, region, box, chans,
            stokes, mask, stretch, wantreturn
        };
        if (_doHistory) {
            auto msgs = this->_newHistory(__func__, names, values);
            padder.addHistory(_ORIGIN, msgs);
        }
        auto out = padder.pad(wantreturn);
        if (wantreturn) {
            return new image(out);
        }
        return nullptr;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::pbcor(
    const variant& pbimage, const string& outfile, bool overwrite,
    const string& box, const variant& region, const string& chans,
    const string& stokes, const string& mask, const string& mode, float cutoff,
    bool stretch
) {
    if (_detached()) {
        ThrowCc("Unable to create image");
    }
    try {
        _log << _ORIGIN;
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        Array<Float> pbPixels;
        SPCIIF pb_ptr;
        if (pbimage.type() == variant::DOUBLEVEC) {
            Vector<Int> shape = pbimage.arrayshape();
            pbPixels.resize(IPosition(shape));
            Vector<Double> localpix(pbimage.getDoubleVec());
            casacore::convertArray(pbPixels, localpix.reform(IPosition(shape)));
        }
        else if (pbimage.type() == variant::STRING) {
            ImageInterface<Float>* pb;
            ImageUtilities::openImage(pb, pbimage.getString());
            ThrowIf(
                ! pb, "Unable to open primary beam image " + pbimage.getString()
            );
            pb_ptr.reset(pb);
        }
        else {
            ThrowCc(
                "Unsupported type " + pbimage.typeString() + " for pbimage"
            );
        }
        auto myRegion = _getRegion(region, true);
        String modecopy = mode;
        modecopy.downcase();
        modecopy.trim();
        if (! modecopy.startsWith("d") && ! modecopy.startsWith("m")) {
            throw AipsError("Unknown mode " + mode);
        }
        ImagePrimaryBeamCorrector::Mode myMode = modecopy.startsWith("d")
            ? ImagePrimaryBeamCorrector::DIVIDE
            : ImagePrimaryBeamCorrector::MULTIPLY;
        Bool useCutoff = cutoff >= 0.0;
        SPCIIF shImage = _imageF;
        std::unique_ptr<ImagePrimaryBeamCorrector> pbcor(
            (!pb_ptr)
            ? new ImagePrimaryBeamCorrector(
                shImage, pbPixels, myRegion.get(), "", box, chans, stokes, mask,
                outfile, overwrite, cutoff, useCutoff, myMode
            )
            : new ImagePrimaryBeamCorrector(
                shImage, pb_ptr, myRegion.get(), "", box, chans, stokes, mask,
                outfile, overwrite, cutoff, useCutoff, myMode
            )
        );
        pbcor->setStretch(stretch);
        if (_doHistory) {
            vector<String> names = {
                "pbimage", "outfile", "overwrite", "box",
                "region", "chans", "stokes", "mask", "mode",
                "cutoff", "stretch"
            };
            vector<variant> values = {
                pbimage.type() == variant::DOUBLEVEC && pbimage.size() > 100
                    ? "(...)" : pbimage,
                    outfile, overwrite, box, region, chans,
                    stokes, mask, mode, cutoff, stretch
            };
            auto msgs = _newHistory(__func__, names, values);
            pbcor->addHistory(_ORIGIN, msgs);
        }
        return new image(pbcor->correct(true));
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

std::string image::pixeltype() {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return "";
        }
        if (_imageF) {
            return "float";
        }
        else if (_imageC) {
            return "complex";
        }
        else if (_imageD) {
            return "double";
        }
        else if (_imageDC) {
            return "dcomplex";
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
             << LogIO::POST;
        RETHROW(x);
    }
    return "";
}

record* image::pixelvalue(const vector<int>& pixel) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        _notSupported(__func__);
        if(_imageF) {
            PixelValueManipulator<Float> pvm(
                _imageF, nullptr, "", false
            );
            return fromRecord(pvm.pixelValue(Vector<Int> (pixel)));
        }
        else {
            PixelValueManipulator<Complex> pvm(
                _imageC, nullptr, "", false
            );
            return fromRecord(pvm.pixelValue(Vector<Int> (pixel)));
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::putchunk(
    const variant& pixels,
    const vector<int>& blc, const vector<int>& inc,
    bool list, bool locking, bool replicate
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        if (_imageF) {
            _putchunk(
                _imageF, pixels, blc, inc, list, locking, replicate
            );
        }
        else if(_imageD) {
            _putchunk(
                _imageD, pixels, blc, inc, list, locking, replicate
            );
        }
        else {
            if (
                pixels.type() == variant::COMPLEXVEC
            ) {
                Vector<DComplex> pixelVector(pixels.getComplexVec());
                auto shape = pixels.arrayshape();
                auto reshapedArray = pixelVector.reform(IPosition(shape));
                if (_imageC) {
                    Array<Complex> pixelsArray;
                    pixelsArray.resize(IPosition(shape));
                    // Vector<DComplex> localpix(pixelVector);
                    casacore::convertArray(pixelsArray, reshapedArray);
                    PixelValueManipulator<Complex>::put(
                        _imageC, pixelsArray, Vector<Int>(blc),
                        Vector<Int>(inc), list, locking,
                        replicate
                    );
                }
                else if (_imageDC) {
                    PixelValueManipulator<DComplex>::put(
                        _imageDC, reshapedArray, Vector<Int>(blc),
                        Vector<Int>(inc), list, locking,
                        replicate
                    );
                }
            }
            else if (_imageC) {
                _putchunk(
                    _imageC, pixels, blc, inc,
                    list, locking, replicate
                );
            }
            else if (_imageDC) {
                _putchunk(
                    _imageDC, pixels, blc, inc,
                    list, locking, replicate
                );
            }
        }
        vector<String> names {
            "pixels", "blc", "inc",
            "list", "locking", "replicate"
        };
        vector<variant> values {
            pixels.size() > 100 ? "[...]" : pixels, blc, inc,
            list, locking, replicate
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

template<class T> void image::_putchunk(
    SPIIT myimage, const variant& pixels,
    const vector<int>& blc, const vector<int>& inc,
    const bool list, const bool locking, const bool replicate
) {
    Array<T> pixelsArray;
    Vector<Int> shape = pixels.arrayshape();
    pixelsArray.resize(IPosition(shape));
    if (pixels.type() == variant::DOUBLEVEC) {
        std::vector<double> pixelVector = pixels.getDoubleVec();
        Vector<Double> localpix(pixelVector);
        casacore::convertArray(pixelsArray, localpix.reform(IPosition(shape)));
    }
    else if (pixels.type() == variant::INTVEC) {
        std::vector<int> pixelVector = pixels.getIntVec();
        Vector<Int> localpix(pixelVector);
        casacore::convertArray(pixelsArray, localpix.reform(IPosition(shape)));
    }
    else {
        String types = myimage->dataType() == TpFloat
            ? "doubles or ints"
            : "complexes, doubles, or ints";
        ThrowCc(
            "Unsupported type for pixels parameter. It "
            "must be either a vector of " + types
        );
    }
    PixelValueManipulator<T>::put(
        myimage, pixelsArray, Vector<Int>(blc),
        Vector<Int>(inc), list, locking,
        replicate
    );
}

bool image::putregion(
    const variant& v_pixels, const variant& v_pixelmask,
    const variant& region, bool list, bool usemask,
    bool, bool replicateArray
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        Bool ret;
        if (_imageF) {
            ret = _putregionReal(
                _imageF, v_pixels, v_pixelmask, region,
                list, usemask, replicateArray
            );
        }
        else if (_imageC) {
            ret = _putregionComplex(
                _imageC, v_pixels, v_pixelmask, region,
                list, usemask, replicateArray
            );
        }
        else if (_imageD) {
            ret = _putregionReal(
                _imageD, v_pixels, v_pixelmask, region,
                list, usemask, replicateArray
            );
        }
        else if (_imageDC) {
            ret = _putregionComplex(
                _imageDC, v_pixels, v_pixelmask, region,
                list, usemask, replicateArray
            );
        }
        else {
            ThrowCc("Logic error")
        }
        if (ret) {
            _statsF.reset();
            _statsD.reset();
            vector<String> names {
                "pixels", "pixelmask", "region",
                "list", "usemask", "replicate"
            };
            vector<variant> values {
                v_pixels.size() > 100 ? "[...]" : v_pixels,
                v_pixelmask.size() > 100 ? "[...]" : v_pixelmask,
                region, list, usemask,
                replicateArray
            };
            _addHistory(__func__, names, values);
        }
        else {
            ThrowCc("Error putting region.");
        }
        return ret;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}


template<class T> bool image::_putregionComplex(
    SPIIT image, const variant& v_pixels, const variant& v_pixelmask,
    const variant& region, bool list, bool usemask,
    bool replicateArray
) {
    Array<T> pixels;
    if (! _isUnset(v_pixels)) {
        IPosition shape(v_pixels.arrayshape());
        auto mytype = v_pixels.type();
        ThrowIf(
            mytype == variant::DOUBLEVEC || mytype == variant::INTVEC,
            "Real values cannot be put in images with complex valued pixels"
        );
        if (mytype == variant::COMPLEXVEC) {
            _convertArray(
                pixels, Vector<DComplex>(v_pixels.getComplexVec()), shape
            );
        }
        else {
            ThrowCc("pixels is not understood, try using an array of real values");
        }
    }
    return _putregion2(
        image, pixels, v_pixelmask, region, list,
        usemask, replicateArray
    );
}

template<class T> bool image::_putregionReal(
    SPIIT image, const variant& v_pixels, const variant& v_pixelmask,
    const variant& region, bool list, bool usemask,
    bool replicateArray
) {
    Array<T> pixels;
    if (! _isUnset(v_pixels)) {
        IPosition shape(v_pixels.arrayshape());
        auto mytype = v_pixels.type();
        ThrowIf(
            mytype == variant::COMPLEXVEC,
            "Complex values cannot be put in images with real valued pixels"
        );
        if (mytype == variant::DOUBLEVEC) {
            _convertArray(
                pixels, Vector<Double>(v_pixels.getDoubleVec()), shape
            );
        }
        else if (mytype == variant::INTVEC) {
            _convertArray(
                pixels, Vector<Int>(v_pixels.getIntVec()), shape
            );
        }
        else {
            ThrowCc("pixels is not understood, try using an array of real values");
        }
    }
    return _putregion2(
        image, pixels, v_pixelmask, region, list,
        usemask, replicateArray
    );
}

template<class T> bool image::_putregion2(
    SPIIT image, const Array<T>& pixels, const variant& v_pixelmask,
    const variant& region, bool list, bool usemask,
    bool replicateArray
) {
    Array<Bool> mask;
    if (! _isUnset(v_pixelmask)) {
        IPosition shape(v_pixelmask.arrayshape());
        auto mytype = v_pixelmask.type();
        if (mytype == variant::DOUBLEVEC) {
            _convertArray(
                mask, Vector<Double>(v_pixelmask.getDoubleVec()), shape
            );
        }
        else if (mytype == variant::INTVEC) {
            _convertArray(
                mask, Vector<Int>(v_pixelmask.getIntVec()), shape
            );
        }
        else if (mytype == ::casac::variant::BOOLVEC) {
            _convertArray(
                mask, Vector<Bool>(v_pixelmask.getBoolVec()), shape
            );
        }
        else {
            ThrowCc("mask is not understood, try using an array");
        }
    }
    if (pixels.empty() && mask.empty()) {
        ThrowCc(
            "You must specify at least either the pixels or the mask"
        );
    }
    auto theRegion = _getRegion(region, false);
    return PixelValueManipulator<T>::putRegion(
        image, pixels, mask, *theRegion,
        list, usemask, replicateArray
    );
}

template<class T, class U>
void image::_convertArray(
    casacore::Array<T>& out, const casacore::Vector<U>& in,
    const IPosition& shape
) {
    out.resize(shape);
    convertArray(out, in.reform(shape));
}

image* image::pv(
    const string& outfile, const variant& start,
    const variant& end, const variant& center, const variant& length,
    const variant& pa, const variant& width, const string& unit,
    bool overwrite, const variant& region, const string& chans,
    const string& stokes, const string& mask, bool stretch,
    bool wantreturn
) {
    if (_detached()) {
        return nullptr;
    }
    try {
        _log << _ORIGIN;
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        std::shared_ptr<casacore::MDirection> startMD, endMD, centerMD;
        Vector<Double> startPix, endPix, centerPix;
        std::shared_ptr<casacore::Quantity> lengthQ;
        Double lengthD = 0;
        if (! start.empty() && ! end.empty()) {
            ThrowIf(
                ! center.empty() || ! length.empty()
                || ! pa.empty(),
                "None of center, length, nor pa may be specified if start and end are specified"
            );
            ThrowIf(
                start.type() != end.type(),
                "start and end must be the same data type"
            );
            casacore::MDirection dir;
            _processDirection(startPix, dir, start, String("start"));
            if (startPix.size() == 0) {
                startMD.reset(new casacore::MDirection(dir));
            }
            _processDirection(endPix, dir, end, "end");
            if (endPix.size() == 0) {
                endMD.reset(new casacore::MDirection(dir));
            }
        }
        else if (
            ! center.empty() && ! length.empty()
            && ! pa.empty()
        ) {
            ThrowIf(
                ! start.empty() || ! end.empty(),
                "Neither start nor end may be specified "
                "if center, length, and pa are specified"
            );
            casacore::MDirection dir;
            _processDirection(centerPix, dir, center, "center");
            if (centerPix.size() == 0) {
                centerMD.reset(new casacore::MDirection(dir));
            }
            if (length.type() == variant::INT || length.type() == variant::DOUBLE) {
                lengthD = length.toDouble();
            }
            else {
                lengthQ.reset(
                    new casacore::Quantity(_casaQuantityFromVar(length))
                );
            }
        }
        else {
            ThrowCc(
                "Either both of start and end or all three of "
                "center, width, and pa must be specified"
            );
        }
        _log << _ORIGIN;

        uInt intWidth = 0;
        casacore::Quantity qWidth;
        if (width.type() == variant::INT) {
            intWidth = width.toInt();
            ThrowIf(
                intWidth % 2 == 0 || intWidth < 1,
                "width must be an odd integer >= 1"
            );
        }
        else if (
            width.type() == variant::STRING
            || width.type() == variant::RECORD
        ) {
            qWidth = _casaQuantityFromVar(width);
        }
        else if (width.type() == variant::BOOLVEC && width.empty()) {
            intWidth = 1;
        }
        else {
            ThrowCc("Unsupported data type for width " + width.typeString());
        }
        if (outfile.empty() && ! wantreturn) {
            _log << LogIO::WARN << "outfile was not specified and wantreturn is false. "
                << "The resulting image will be inaccessible" << LogIO::POST;
        }
        std::shared_ptr<Record> regionPtr = _getRegion(region, true);
        PVGenerator pv(
            _imageF, regionPtr.get(),
            chans, stokes, mask, outfile, overwrite
        );
        if (startPix.size() == 2) {
            pv.setEndpoints(
                make_pair(startPix[0], startPix[1]),
                make_pair(endPix[0], endPix[1])
            );
        }
        else if (startMD) {
            pv.setEndpoints(*startMD, *endMD);
        }
        else if (centerMD) {
            if (lengthQ) {
                pv.setEndpoints(
                    *centerMD, *lengthQ, _casaQuantityFromVar(pa)
                );
            }
            else {
                pv.setEndpoints(
                    *centerMD, lengthD, _casaQuantityFromVar(pa)
                );
            }
        }
        else {
            if (lengthQ) {
                pv.setEndpoints(
                    make_pair(centerPix[0], centerPix[1]),
                    *lengthQ, _casaQuantityFromVar(variant(pa))
                );
            }
            else {
                pv.setEndpoints(
                        make_pair(centerPix[0], centerPix[1]),
                    lengthD, _casaQuantityFromVar(variant(pa))
                );
            }
        }
        if (intWidth == 0) {
            pv.setWidth(qWidth);
        }
        else {
            pv.setWidth(intWidth);
        }
        pv.setStretch(stretch);
        pv.setOffsetUnit(unit);
        _log << _ORIGIN;
        vector<String> names {
            "outfile", "start", "end", "center", "length",
            "pa", "width", "unit", "overwrite", "region", "chans",
            "stokes", "mask", "stretch", "wantreturn"
        };
        vector<variant> values {
            outfile, start, end, center, length,
            pa, width, unit, overwrite, region, chans,
            stokes, mask, stretch, wantreturn
        };
        if (_doHistory) {
            auto msgs = _newHistory(__func__, names, values);
            pv.addHistory(_ORIGIN, msgs);
        }
        auto out = pv.generate();
        image *ret = wantreturn ? new image(out) : nullptr;
        return ret;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

void image::_processDirection(
    Vector<Double>& pixel, casacore::MDirection& dir,
    const variant& inputDirection, const String& paramName
) {
    variant::TYPE myType = inputDirection.type();
    ThrowIf(
        (
            myType == variant::INTVEC
            || myType == variant::DOUBLEVEC
            || myType == variant::STRINGVEC
        ) &&
        inputDirection.size() != 2,
        "If specified as an array, " + paramName
        + " must have exactly two elements"
    );
    pixel.resize(0);
    if (myType == variant::INTVEC || myType == variant::DOUBLEVEC) {
        pixel = Vector<Double>(_toDoubleVec(inputDirection));
    }
    else if(myType == variant::STRINGVEC) {
        vector<string> x = inputDirection.toStringVec();
        casacore::Quantity q0 = _casaQuantityFromVar(variant(x[0]));
        casacore::Quantity q1 = _casaQuantityFromVar(variant(x[1]));
        dir = casacore::MDirection(q0, q1);
    }
    else if (myType == variant::STRING) {
        string parts[3];
        split(inputDirection.toString(), parts, 3, Regex("[, \n\t\r\v\f]+"));
        casacore::MDirection::Types frame;
        casacore::MDirection::getType(frame, parts[0]);
        dir = casacore::MDirection::getType(frame, parts[0])
            ? casacore::MDirection(
                _casaQuantityFromVar(parts[1]),
                _casaQuantityFromVar(parts[2]), frame
            )
            : casacore::MDirection(
                _casaQuantityFromVar(parts[0]),
                _casaQuantityFromVar(parts[1])
            );
    }
    else {
        ThrowCc("Unsupported type for " + paramName);
    }
}

image* image::rebin(
    const string& outfile, const vector<int>& bin,
    const variant& region, const variant& vmask,
    bool dropdeg, bool overwrite,
    bool stretch, bool crop
) {
    LogOrigin lor(_class, __func__);
    _log << lor;
    ThrowIf(
        _detached(), "Unable to create image"
    );
    Vector<Int> mybin(bin);
    ThrowIf(
        anyTrue(mybin <= 0),
        "All binning factors must be positive."
    );
    try {
        _notSupported(__func__);
        vector<String> msgs;
        if (_doHistory) {
            vector<String> names {
                "outfile", "bin", "region", "mask",
                "dropdeg", "overwrite", "stretch", "crop"
            };
            vector<variant> values {
                outfile, bin, region, vmask,
                dropdeg, overwrite, stretch, crop
            };
            msgs = _newHistory(__func__, names, values);
        }
        auto mask = _getMask(vmask);
        if (_imageF) {
            SPIIF myfloat = _imageF;
            ImageRebinner<Float> rebinner(
                myfloat, _getRegion(region, true).get(),
                mask, outfile, overwrite
            );
            rebinner.setFactors(mybin);
            rebinner.setStretch(stretch);
            rebinner.setDropDegen(dropdeg);
            if (_doHistory) {
                rebinner.addHistory(lor, msgs);
            }
            rebinner.setCrop(crop);
            return new image(rebinner.rebin());
        }
        else {
            SPIIC myComplex = _imageC;
            ImageRebinner<Complex> rebinner(
                myComplex, _getRegion(region, true).get(),
                mask, outfile, overwrite
            );
            rebinner.setFactors(mybin);
            rebinner.setStretch(stretch);
            rebinner.setDropDegen(dropdeg);
            if (_doHistory) {
                rebinner.addHistory(lor, msgs);
            }
            rebinner.setCrop(crop);
            return new image(rebinner.rebin());
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::regrid(
    const string& outfile, const vector<int>& inshape, const record& csys,
    const vector<int>& inaxes, const variant& region, const variant& vmask,
    const string& method, int decimate, bool replicate, bool doRefChange,
    bool dropDegenerateAxes, bool overwrite, bool forceRegrid,
    bool specAsVelocity, bool stretch
) {
    try {
        LogOrigin lor(_class, __func__);
        _log << lor;
        if (_detached()) {
            ThrowCc("Unable to create image");
        }
        unique_ptr<Record> csysRec(toRecord(csys));
        unique_ptr<CoordinateSystem> coordinates(
            CoordinateSystem::restore(*csysRec, "")
        );
        ThrowIf (
            ! coordinates.get(), "Invalid specified coordinate system record."
        );
        auto regionPtr = _getRegion(region, true);
        String mask = _getMask(vmask);
        Vector<Int> axes;
        if (!((inaxes.size() == 1) && (inaxes[0] == -1))) {
            axes = inaxes;
        }
        vector<String> msgs;
        if (_doHistory) {
            vector<String> names {
                "outfile", "shape", "csys", "axes", "region", "mask", "method",
                "decimate", "replicate", "doref", "dropdegen", "overwrite",
                "force", "asvelocity", "stretch"
            };
            vector<variant> values {
                outfile, inshape, csys, inaxes, region, vmask, method, decimate,
                replicate, doRefChange, dropDegenerateAxes, overwrite,
                forceRegrid, specAsVelocity, stretch
            };
            msgs = _newHistory(__func__, names, values);
        }
        if (_imageF) {
            ImageRegridder<Float> regridder(
                _imageF, regionPtr.get(), mask, outfile, overwrite,
                *coordinates, IPosition(axes), IPosition(inshape)
            );
            return _regrid(
                regridder, method, decimate, replicate, doRefChange,
                forceRegrid, specAsVelocity, stretch, dropDegenerateAxes, lor,
                msgs
            );
        }
        else if (_imageC) {
            ImageRegridder<Complex> regridder(
                _imageC, regionPtr.get(), mask, outfile, overwrite,
                *coordinates, IPosition(axes), IPosition(inshape)
            );
            return _regrid(
                regridder, method, decimate, replicate, doRefChange,
                forceRegrid, specAsVelocity, stretch, dropDegenerateAxes, lor,
                msgs
            );
        }
        if (_imageD) {
            ImageRegridder<Double> regridder(
                _imageD, regionPtr.get(), mask, outfile, overwrite,
                *coordinates, IPosition(axes), IPosition(inshape)
            );
            return _regrid(
                regridder, method, decimate, replicate, doRefChange,
                forceRegrid, specAsVelocity, stretch, dropDegenerateAxes, lor,
                msgs
            );
        }
        else if (_imageDC) {
            ImageRegridder<DComplex> regridder(
                _imageDC, regionPtr.get(), mask, outfile, overwrite,
                *coordinates, IPosition(axes), IPosition(inshape)
            );
            return _regrid(
                regridder, method, decimate, replicate, doRefChange,
                forceRegrid, specAsVelocity, stretch, dropDegenerateAxes, lor,
                msgs
            );
        }
        else {
            ThrowCc("Logic Error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> image* image::_regrid(
    ImageRegridderBase<T>& regridder,
    const string& method, int decimate, bool replicate,
    bool doRefChange, bool forceRegrid,
    bool specAsVelocity, bool stretch,
    bool dropDegenerateAxes, const LogOrigin& lor,
    const vector<String>& msgs
) const {
    regridder.setMethod(method);
    regridder.setDecimate(decimate);
    regridder.setReplicate(replicate);
    regridder.setDoRefChange(doRefChange);
    regridder.setForceRegrid(forceRegrid);
    regridder.setSpecAsVelocity(specAsVelocity);
    regridder.setStretch(stretch);
    regridder.setDropDegen(dropDegenerateAxes);
    if (_doHistory) {
        regridder.addHistory(lor, msgs);
    }
    return new image(regridder.regrid());
}

bool image::remove(const bool finished, const bool verbose) {
    try {
        _log << _ORIGIN;

        if (_detached()) {
            return false;
        }
        _remove(verbose);
        if (finished) {
            done();
        }
        return true;
    }
    catch (const AipsError &x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

void image::_remove(bool verbose) {
    auto imageF = _imageF;
    auto imageC = _imageC;
    auto imageD = _imageD;
    auto imageDC = _imageDC;
    _reset();
    if (imageF) {
        ImageFactory::remove(imageF, verbose);
    }
    else if (imageC) {
        ImageFactory::remove(imageC, verbose);
    }
    else if (imageD) {
        ImageFactory::remove(imageD, verbose);
    }
    else if (imageDC) {
        ImageFactory::remove(imageDC, verbose);
    }
    else {
        ThrowCc("Logic error");
    }
}

bool image::removefile(const std::string& filename) {
    bool rstat(false);
    try {
        _log << LogOrigin("image", "removefile");

        String fileName(filename);
        if (fileName.empty()) {
            _log << LogIO::WARN << "Empty filename" << LogIO::POST;
            return rstat;
        }
        File f(fileName);
        if (!f.exists()) {
            _log << LogIO::WARN << fileName << " does not exist."
                    << LogIO::POST;
            return rstat;
        }

        // Now try and blow it away.  If it's open, tabledelete won't delete it.
        String message;
        if (Table::canDeleteTable(message, fileName, true)) {
            Table::deleteTable(fileName, true);
            rstat = true;
        } else {
            _log << LogIO::WARN << "Cannot delete file " << fileName
                    << " because " << message << LogIO::POST;
        }
    } catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return rstat;
}

bool image::rename(const string& name, bool overwrite) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        _statsF.reset();
        _statsD.reset();
        auto oldName = this->name(False);
        if (_imageF) {
            auto myimage = _imageF;
            _imageF.reset();
            ImageFactory::rename(myimage, name, overwrite);
            _imageF = myimage;
        }
        else {
            auto myimage = _imageC;
            _imageC.reset();
            ImageFactory::rename(myimage, name, overwrite);
            _imageC = myimage;
        }
        vector<String> names = { "name", "overwrite" };
        vector<variant> values = { name, overwrite };
        vector<String> appendMsgs = {"original name was " + oldName};
        _addHistory(__func__, names, values, appendMsgs);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::replacemaskedpixels(
    const variant& pixels, const variant& region,
    const variant& vmask, bool updateMask, bool list,
    bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        return false;
    }
    try {
        _notSupported(__func__);
        auto regionPtr = _getRegion(region, true);
        auto mask = _getMask(vmask);
        if (_imageF) {
            auto myfloat = _imageF;
            ImageMaskedPixelReplacer<Float> impr(
                myfloat, regionPtr.get(), mask
            );
            impr.setStretch(stretch);
            impr.replace(pixels.toString(), updateMask, list);
            _imageF = myfloat;
        }
        else {
            auto mycomplex = _imageC;
            ImageMaskedPixelReplacer<Complex> impr(
                mycomplex, regionPtr.get(), mask
            );
            impr.setStretch(stretch);
            impr.replace(pixels.toString(), updateMask, list);
            _imageC = mycomplex;
        }
        _statsF.reset();
        _statsD.reset();
        vector<String> names = {
            "pixels", "region", "mask", "update",
            "list", "stretch"
        };
        vector<variant> values = {
            pixels, region, vmask, updateMask,
            list, stretch
        };
        _addHistory(__func__, names, values);
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

record* image::restoringbeam(int channel, int polarization) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return nullptr;
        }
        if (_imageF) {
            return fromRecord(
                _imageF->imageInfo().beamToRecord(channel, polarization)
            );
        }
        else if (_imageC) {
            return fromRecord(
                _imageC->imageInfo().beamToRecord(channel, polarization)
            );
        }
        else if (_imageD) {
            return fromRecord(
                _imageD->imageInfo().beamToRecord(channel, polarization)
            );
        }
        else if (_imageDC) {
            return fromRecord(
                _imageDC->imageInfo().beamToRecord(channel, polarization)
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

image* image::rotate(
    const string& outfile, const vector<int>& inshape, const variant& inpa,
    const variant& region, const variant& vmask, const string& method,
    int decimate, bool replicate, bool dropdeg, bool overwrite, bool stretch
) {
    try {
        _log << _ORIGIN;
        ThrowIf(_detached(), "Unable to create image");
        if (_imageF) {
            return _rotate(
                _imageF, outfile, inshape, inpa, region, vmask, method,
                decimate, replicate, dropdeg, overwrite, stretch
            );
        }
        else if (_imageC) {
            return _rotate(
                _imageC, outfile, inshape, inpa, region, vmask, method,
                decimate, replicate, dropdeg, overwrite, stretch
            );
        }
        else if (_imageD) {
            return _rotate(
                _imageD, outfile, inshape, inpa, region, vmask, method,
                decimate, replicate, dropdeg, overwrite, stretch
            );
        }
        else if (_imageDC) {
            return _rotate(
                _imageDC, outfile, inshape, inpa, region, vmask, method,
                decimate, replicate, dropdeg, overwrite, stretch
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> image* image::_rotate(
    SPIIT myImage, const string& outfile, const vector<int>& inshape,
    const variant& inpa, const variant& region, const variant& vmask,
    const string& method, int decimate, bool replicate, bool dropdeg,
    bool overwrite, bool stretch
) {
    Vector<Int> shape(inshape);
    if (shape.size() == 1 && shape[0] == -1) {
        shape.resize(0);
    }
    auto pa = _casaQuantityFromVar(inpa);
    auto Region = _getRegion(region, false);
    auto mask = _getMask(vmask);
    ImageRotator<T> rotator(myImage, Region.get(), mask, outfile, overwrite);
    rotator.setShape(IPosition(shape));
    rotator.setAngle(pa);
    rotator.setInterpolationMethod(method);
    rotator.setDecimate(decimate);
    rotator.setReplicate(replicate);
    rotator.setDropDegen(dropdeg);
    rotator.setStretch(stretch);
    vector<String> names = {
        "outfile", "shape", "pa", "region", "mask", "method", "decimate",
        "replicate", "dropdeg", "overwrite", "stretch"
    };
    vector<variant> values = {
        outfile, inshape, inpa, region, vmask, method,
        decimate, replicate, dropdeg, overwrite, stretch
    };
    if (_doHistory) {
        auto msgs = _newHistory("rotate", names, values);
        rotator.addHistory(_ORIGIN, msgs);
    }
    auto x = rotator.rotate();
    _log << LogIO::NORMAL << "Using position angle rotation "
        << inpa.toString() << LogIO::POST;
    return new image(x);
}

bool image::rotatebeam(const variant& angle) {
    try {
        _log << _ORIGIN;
        if(_detached()) {
            return false;
        }
        _notSupported(__func__);
        Quantum<Double> pa(_casaQuantityFromVar(angle));
        vector<String> msgs;
        if (_doHistory) {
            vector<String> names { "angle" };
            vector<variant> values { angle };
            msgs = _newHistory(__func__, names, values);
        }
        if (_imageF) {
            BeamManipulator<Float> bManip(_imageF);
            bManip.rotate(pa, msgs);
        }
        else {
            BeamManipulator<Complex> bManip(_imageC);
            bManip.rotate(pa, msgs);
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

image* image::sepconvolve(
    const string& outfile, const vector<int>& axes,
    const vector<std::string>& types,
    const variant& widths,
    double scale, const variant& region,
    const variant& vmask, bool overwrite, bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        ThrowCc("Unable to create image");
    }
    try {
        _notSupported(__func__);
        UnitMap::putUser("pix", UnitVal(1.0), "pixel units");
        auto pRegion = _getRegion(region, false);
        auto mask = _getMask(vmask);
        Vector<Int> smoothaxes(axes);
        auto kernels = toVectorString(types);

        Int num = 0;
        Vector<Quantum<Double> > kernelwidths;
        if (widths.type() == ::casac::variant::INTVEC) {
            vector<int> widthsIVec = widths.toIntVec();
            num = widthsIVec.size();
            vector<double> widthsVec(num);
            for (int i = 0; i < num; ++i) {
                widthsVec[i] = widthsIVec[i];
            }
            kernelwidths.resize(num);
            Unit u("pix");
            for (int i = 0; i < num; ++i) {
                kernelwidths[i] = casacore::Quantity(widthsVec[i], u);
            }
        }
        else if (widths.type() == variant::DOUBLEVEC) {
            std::vector<double> widthsVec = widths.toDoubleVec();
            num = widthsVec.size();
            kernelwidths.resize(num);
            Unit u("pix");
            for (int i = 0; i < num; ++i) {
                kernelwidths[i] = casacore::Quantity(widthsVec[i], u);
            }
        }
        else if (
            widths.type() == variant::STRING
            || widths.type() == variant::STRINGVEC
        ) {
            toCasaVectorQuantity(widths, kernelwidths);
            num = kernelwidths.size();
        }
        else {
            _log << LogIO::WARN << "Unrecognized kernelwidth datatype"
                    << LogIO::POST;
            return nullptr;
        }
        if (kernels.size() == 1 && kernels[0] == "") {
            kernels.resize(num);
            for (int i = 0; i < num; ++i) {
                kernels[i] = "gauss";
            }
        }
        if (
            smoothaxes.size() == 0 || ((smoothaxes.size() == 1)
            && (smoothaxes[0] = -1))
        ) {
            smoothaxes.resize(num);
            for (int i = 0; i < num; ++i) {
                smoothaxes[i] = i;
            }
        }
        SepImageConvolverTask<Float> task(
            _imageF, pRegion.get(), mask,
            outfile, overwrite
        );
        task.setScale(scale);
        task.setSmoothAxes(smoothaxes);
        task.setKernels(kernels);
        task.setKernelWidths(kernelwidths);
        task.setStretch(stretch);
        if (_doHistory) {
            vector<String> names {
                "outfile", "axes", "types", "widths", "scale",
                "region", "mask", "overwrite", "stretch"
            };
            vector<variant> values {
                outfile, axes, types, widths, scale,
                region, vmask, overwrite, stretch
            };
            auto msgs = _newHistory(__func__, names, values);
            task.addHistory(_ORIGIN, msgs);
        }
        return new image(task.convolve());
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

bool image::set(
    const variant& vpixels, int pixelmask,
    const variant& region, bool list
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        auto pixels = vpixels.toString();
        if (pixels == "[]") {
            pixels = "";
        }
        auto pRegion = _getRegion(region, false);
        if (pixels == "" && pixelmask == -1) {
            _log << LogIO::WARN
                << "You must specify at least either the pixels or the mask to set"
                << LogIO::POST;
            return false;
        }
        if (
            PixelValueManipulator<Float>::set(
                _imageF, pixels, pixelmask, *pRegion, list
            )
        ) {
            _statsF.reset();
            _statsD.reset();
            if (_doHistory) {
                vector<String> names = {
                    "pixels", "pixelmask", "region", "list"
                };
                vector<variant> values = {
                    vpixels, pixelmask, region, list
                };
                _addHistory(__func__, names, values);
            }
            return true;
        }
        ThrowCc("Error setting pixel values.");
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::setbrightnessunit(const std::string& unit) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        Bool res = _imageF
            ?  _imageF->setUnits(Unit(unit))
            : _imageC->setUnits(Unit(unit));
        ThrowIf(! res, "Unable to set brightness unit");
        if (_doHistory) {
            vector<String> names = {"unit"};
            vector<variant> values = {unit};
            _addHistory(__func__, names, values);
        }
        _statsF.reset();
        _statsD.reset();
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::setcoordsys(const record& csys) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        unique_ptr<Record> coordinates(toRecord(csys));
        if (_imageF) {
            ImageMetaDataRW<Float> md(_imageF);
            md.setCsys(*coordinates);
        }
        else if (_imageC) {
            ImageMetaDataRW<Complex> md(_imageC);
            md.setCsys(*coordinates);
        }
        else if (_imageD) {
            ImageMetaDataRW<Double> md(_imageD);
            md.setCsys(*coordinates);
        }
        else if (_imageDC) {
            ImageMetaDataRW<DComplex> md(_imageDC);
                md.setCsys(*coordinates);
        }
        else {
            ThrowCc("Logic error");
        }
        if (_doHistory) {
            vector<String> names = {"csys"};
            vector<variant> values = {csys};
            _addHistory(__func__, names, values);
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return true;
}

bool image::sethistory(
    const string& origin, const vector<string>& history
) {
    try {
        if (_detached()) {
            return false;
        }
        if ((history.size() == 1) && (history[0].size() == 0)) {
            LogOrigin lor("image", "sethistory");
            _log << lor << "history string is empty" << LogIO::POST;
        }
        else {
            if(_imageF) {
                ImageHistory<Float> hist(_imageF);
                hist.addHistory(origin, history);
            }
            else if (_imageD) {
                ImageHistory<Double> hist(_imageD);
                hist.addHistory(origin, history);
            }
            else if (_imageC) {
                ImageHistory<Complex> hist(_imageC);
                hist.addHistory(origin, history);
            }
            else if (_imageDC) {
                ImageHistory<DComplex> hist(_imageDC);
                hist.addHistory(origin, history);
            }
            else {
                ThrowCc("Logic error");
            }
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << _ORIGIN << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::setmiscinfo(const record& info) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        std::unique_ptr<Record> tmp(toRecord(info));
        auto res = False;
        if (_imageF) {
            res = _imageF->setMiscInfo(*tmp);
        }
        else if (_imageC) {
            res = _imageC->setMiscInfo(*tmp);
        }
        else if (_imageD) {
            res = _imageD->setMiscInfo(*tmp);
        }
        else if (_imageDC) {
            res = _imageDC->setMiscInfo(*tmp);
        }
        else {
            ThrowCc("Logic error");
        }
        if (res && _doHistory) {
            vector<String> names {"info"};
            vector<variant> values {info};
            _addHistory(__func__, names, values);
        }
        return res;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::setrestoringbeam(
    const variant& major, const variant& minor, const variant& pa,
    const record& beam, bool remove, bool log, int channel, int polarization,
    const string& imagename
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return false;
        }
        std::unique_ptr<Record> rec(toRecord(beam));
        ImageBeamSet bs;
        if (! imagename.empty()) {
            ThrowIf(
                ! major.empty() || ! minor.empty() || ! pa.empty(),
                "Cannot specify both imagename and major, minor, and/or pa"
            );
            ThrowIf(
                remove, "remove cannot be true if imagename is specified"
            );
            ThrowIf(
                ! beam.empty(),
                "beam must be empty if imagename specified"
            );
            ThrowIf(
                channel >= 0 || polarization >= 0,
                "Neither channel nor polarization can be non-negative if "
                "imagename is specified"
            );
            PtrHolder<ImageInterface<Float> > k;
            ImageUtilities::openImage(k, imagename);
            if (k.ptr() == 0) {
                PtrHolder<ImageInterface<Float> > c;
                ImageUtilities::openImage(c, imagename);
                ThrowIf(
                    c.ptr() == 0,
                    "Unable to open " + imagename
                );
                bs = c->imageInfo().getBeamSet();
            }
            else {
                bs = k->imageInfo().getBeamSet();
            }
            ThrowIf(
                bs.empty(),
                "Image " + imagename + " has no beam"
            );
        }
        if (_imageF) {
            _setrestoringbeam(
                _imageF, major, minor, pa, remove, log,
                channel, polarization, *rec, bs
            );
        }
        else if (_imageC) {
            _setrestoringbeam(
                _imageC, major, minor, pa, remove, log,
                channel, polarization, *rec, bs
            );
        }
        else if (_imageD) {
            _setrestoringbeam(
                _imageD, major, minor, pa, remove, log,
                channel, polarization, *rec, bs
            );
        }
        else if (_imageDC) {
            _setrestoringbeam(
                _imageDC, major, minor, pa, remove, log,
                channel, polarization, *rec, bs
            );
        }
        else {
            ThrowCc("Logic error");
        }
        variant myb = beam;
        std::set<String> dontQuote;
        if (rec && ! rec->empty()) {
            std::ostringstream oss;
            auto paKey = rec->isDefined("pa") ? "pa" : "positionangle";
            oss << "{'major': " << _quantityRecToString(rec->asRecord("major"))
                << ", 'minor': " << _quantityRecToString(rec->asRecord("minor"))
                << ", 'positionangle': "
                << _quantityRecToString(rec->asRecord(paKey)) << "}";
            myb = oss.str();
            dontQuote.insert("beam");
        }
        if (_doHistory) {
            vector<String> names {
                "major", "minor", "pa", "beam", "remove", "log",
                "channel", "polarization", "imagename"
            };
            vector<variant> values {
                major, minor, pa, myb, remove, log,
                channel, polarization, imagename
            };
            _addHistory(
                __func__, names, values, vector<casacore::String>(), dontQuote
            );
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

template<class T> void image::_setrestoringbeam(
    SPIIT image, const variant& major, const variant& minor, const variant& pa,
    bool remove, bool log, int channel, int polarization,
    const Record& rec, const ImageBeamSet& bs
) {
    BeamManipulator<T> bManip(image);
    bManip.setVerbose(log);
    if (remove) {
        bManip.remove();
    }
    else if (! bs.empty()) {
        bManip.set(bs);
    }
    else {
        bManip.set(
            major.empty() ? casacore::Quantity() : _casaQuantityFromVar(major),
            minor.empty() ? casacore::Quantity() : _casaQuantityFromVar(minor),
            pa.empty() ? casacore::Quantity() : _casaQuantityFromVar(pa),
            rec, channel, polarization
        );
    }
}

String image::_quantityRecToString(const Record& q) {
    ostringstream oss;
    oss << "{'value': " << q.asDouble("value")
        << ", 'unit': '" << q.asString("unit") << "'}";
    return oss.str();
}

vector<int> image::shape() {
    _log << _ORIGIN;
    if (_detached()) {
        return vector<int>();
    }
    try {
        vector<int> rstat = _imageF
            ? _imageF->shape().asVector().tovector()
            : _imageC ? _imageC->shape().asVector().tovector()
            : _imageD ? _imageD->shape().asVector().tovector()
            : _imageDC ? _imageDC->shape().asVector().tovector() : vector<int>();
        return rstat;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return vector<int>();
}

record* image::statistics(
    const vector<int>& axes, const variant& region, const variant& mask,
    const vector<double>& includepix, const vector<double>& excludepix,
    bool list, bool force, bool disk, bool robust, bool verbose, bool stretch,
    const string& logfile, bool append, const string& algorithm, double fence,
    const string& center, bool lside, double zscore, int maxiter,
    const string& clmethod, int niter
) {
    _log << _ORIGIN;
    if (_detached()) {
        _log << "Image not attached" << LogIO::POST;
        return nullptr;
    }
    try {
        ThrowIf(
            ! (_imageF || _imageD),
            "This method only supports real valued images"
        );
        if (_imageF) {
            return _statistics(
                _statsF, _imageF, axes, region, mask, includepix,
                excludepix, list, force, disk, robust, verbose, stretch,
                logfile, append, algorithm, fence, center, lside,
                zscore, maxiter, clmethod, niter
            );
        }
        else if (_imageD) {
            return _statistics(
                _statsD, _imageD, axes, region, mask, includepix,
                excludepix, list, force, disk, robust, verbose, stretch,
                logfile, append, algorithm, fence, center, lside,
                zscore, maxiter, clmethod, niter
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_statistics(
    std::unique_ptr<casa::ImageStatsCalculator<T>>& stats, SPIIT myImage,
    const vector<int>& axes, const variant& region,
    const variant& mask, const vector<double>& includepix,
    const vector<double>& excludepix, bool list, bool force, bool disk,
    bool robust, bool verbose, bool stretch, const string& logfile, bool append,
    const string& algorithm, double fence, const string& center, bool lside,
    double zscore, int maxiter, const string& clmethod, int niter
) {
    auto regionRec = _getRegion(region, true);
    auto mtmp = mask.toString();
    if (mtmp == "false" || mtmp == "[]") {
        mtmp = "";
    }
    Vector<Int> tmpaxes(axes);
    if (tmpaxes.size() == 1 && tmpaxes[0] == -1) {
        tmpaxes.resize(0);
    }
    Vector<T> tmpinclude, tmpexclude;
    if (! (includepix.size() == 1 && includepix[0] == -1)) {
        tmpinclude.resize(includepix.size());
        for (uInt i=0; i<includepix.size(); ++i) {
            tmpinclude[i] = includepix[i];
        }
    }
    if (!(excludepix.size() == 1 && excludepix[0] == -1)) {
        tmpexclude.resize(excludepix.size());
        for (uInt i = 0; i < excludepix.size(); ++i) {
            tmpexclude[i] = excludepix[i];
        }
    }
    if (verbose) {
        _log << LogIO::NORMAL << "Determining stats for image "
            << _name(true) << LogIO::POST;
    }
    Record ret;
    if (force || ! stats.get()) {
        stats.reset(
            new ImageStatsCalculator<T>(myImage, regionRec.get(), mtmp, verbose)
        );
    }
    else {
        stats->setMask(mtmp);
        stats->setRegion(regionRec ? *regionRec : Record());
    }
    String myalg = algorithm;
    myalg.downcase();
    if (myalg.startsWith("b")) {
        stats->configureBiweight(niter);
        if (robust) {
            _log << LogIO::WARN << "The biweight algorithm does not support "
                << "computation of quantile-related (median, MADM, first/third "
                << "quartile, IQR) statistics (robust=True). Proceeding "
                << "without calculating those stats." << LogIO::POST;
            robust = False;
        }
    }
    else if (myalg.startsWith("ch")) {
        stats->configureChauvenet(zscore, maxiter);
    }
    else if (myalg.startsWith("cl")) {
        String mymethod = clmethod;
        mymethod.downcase();
        ImageStatsData::PreferredClassicalAlgorithm method;
        if (mymethod.startsWith("a")) {
            method = ImageStatsData::AUTO;
        }
        else if (mymethod.startsWith("t")) {
            method = ImageStatsData::TILED_APPLY;
        }
        else if (mymethod.startsWith("f")) {
            method = ImageStatsData::STATS_FRAMEWORK;
        }
        else {
            ThrowCc("Unsupported classical method " + clmethod);
        }
        stats->configureClassical(method);
    }
    else if (myalg.startsWith("f")) {
        String mycenter = center;
        mycenter.downcase();
        FitToHalfStatisticsData::CENTER centerType;
        if (mycenter.startsWith("mea")) {
            centerType = FitToHalfStatisticsData::CMEAN;
        }
        else if (mycenter.startsWith("med")) {
            centerType = FitToHalfStatisticsData::CMEDIAN;
        }
        else if (mycenter.startsWith("z")) {
            centerType = FitToHalfStatisticsData::CVALUE;
        }
        else {
            ThrowCc("Unsupported center value " + center);
        }
        FitToHalfStatisticsData::USE_DATA useData = lside
            ? FitToHalfStatisticsData::LE_CENTER
            : FitToHalfStatisticsData::GE_CENTER;
        stats->configureFitToHalf(centerType, useData, 0.0);
    }
    else if (myalg.startsWith("h")) {
        stats->configureHingesFences(fence);
    }
    else {
        ThrowCc("Unsupported algorithm " + algorithm);
    }
    stats->setAxes(tmpaxes);
    stats->setIncludePix(tmpinclude);
    stats->setExcludePix(tmpexclude);
    stats->setList(list);
    if (force) {
        stats->forceNewStorage();
    }
    stats->setDisk(disk);
    stats->setRobust(robust);
    stats->setVerbose(verbose);
    stats->setStretch(stretch);
    if (! logfile.empty()) {
        stats->setLogfile(logfile);
    }
    stats->setLogfileAppend(append);
    return fromRecord(stats->calculate());
}

image* image::subimage(
    const string& outfile, const variant& region, const variant& vmask,
    bool dropDegenerateAxes, bool overwrite, bool list, bool stretch,
    bool wantreturn, const vector<int>& keepaxes
) {
    try {
        _log << _ORIGIN;
        ThrowIf(_detached(), "Unable to create image");
        if (outfile.empty() && ! wantreturn) {
            _log << LogIO::WARN << "outfile was not specified and wantreturn "
                << "is false. The resulting image will be inaccessible"
                << LogIO::POST;
        }
        if (_imageF) {
            return _subimage<Float>(
                _imageF, outfile, region, vmask, dropDegenerateAxes, overwrite,
                list, stretch, keepaxes, wantreturn
            );
        }
        else if (_imageC) {
            return _subimage<Complex>(
                _imageC, outfile, region, vmask, dropDegenerateAxes, overwrite,
                list, stretch, keepaxes, wantreturn
            );
        }
        else if (_imageD) {
            return _subimage<Double>(
                _imageD, outfile, region, vmask, dropDegenerateAxes, overwrite,
                list, stretch, keepaxes, wantreturn
            );
        }
        else if (_imageDC) {
            return _subimage<DComplex>(
                _imageDC, outfile, region, vmask, dropDegenerateAxes, overwrite,
                list, stretch, keepaxes, wantreturn
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template<class T> image* image::_subimage(
    SPIIT clone, const String& outfile, const variant& region,
    const variant& vmask, bool dropDegenerateAxes, bool overwrite, bool list,
    bool stretch, const vector<int>& keepaxes, bool wantreturn
) {
    SPIIT im;
    auto regionRec = _getRegion(region, false);
    auto mask = _getMask(vmask);
    if (! dropDegenerateAxes || keepaxes.empty()) {
        im = SubImageFactory<T>::createImage(
            *clone, outfile, *regionRec, mask, dropDegenerateAxes,
            overwrite, list, stretch
        );
    }
    else {
        im = SubImageFactory<T>::createImage(
            *clone, outfile, *regionRec, mask,
            AxesSpecifier(IPosition(Vector<Int>(keepaxes))),
            overwrite, list, stretch
        );
    }
    if (_doHistory) {
        vector<String> names {
            "outfile", "region", "mask", "dropdeg", "overwrite",
            "list", "stretch", "wantreturn", "keepaxes"
        };
        vector<variant> values {
            outfile, region, vmask, dropDegenerateAxes,
            overwrite, list, stretch, wantreturn, keepaxes
        };
        _addHistory(im, "subimage", names, values);
    }
    return wantreturn ? new image(im) : nullptr;
}

record* image::summary(
    const string& doppler, bool list,
    bool pixelorder, bool verbose
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            return 0;
        }
        if (_imageF) {
            return _summary(
                _imageF, doppler, list, pixelorder, verbose
            );
        }
        else if (_imageC) {
            return _summary(
                _imageC, doppler, list, pixelorder, verbose
            );
        }
        else if (_imageD) {
            return _summary(
                _imageD, doppler, list, pixelorder, verbose
            );
        }
        else if (_imageDC) {
            return _summary(
                _imageDC, doppler, list, pixelorder, verbose
            );
        }
        else {
            ThrowCc("Logic error");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_summary(
    SPIIT image, 
    const string& doppler, bool list,
    bool pixelorder, bool verbose
) {
    ImageMetaData<T> md(image);
    return fromRecord(
        md.summary(doppler, list, pixelorder, verbose)
    );
}

bool image::toASCII(
    const string& outfile, const variant& region,
    const variant& mask, const string& sep,
    const string& format, double maskvalue, bool overwrite,
    bool stretch
) {
    // sep is hard-wired as ' ' which is what imagefromascii expects
    _log << _ORIGIN;
    _log << LogIO::WARN << __func__ << "() IS DEPRECATED AND WILL BE REMOVED "
        << "IN A NEAR-FUTURE VERSION OF CASA. YOU SHOULD USE ANOTHER IMAGE "
        << "EXPORT METHOD SUCH AS tofits() TO EXPORT CASA IMAGES. IF YOU "
        << "SIMPLY WISH TO MODIFY PIXEL VALUES, USE getchunk()/putchunk() OR "
        << "getregion()/putregion() FOR THAT" << LogIO::POST;
    if (_detached()) {
        return false;
    }
    try {
        _notSupported(__func__);
        String Mask;
        if (mask.type() == variant::BOOLVEC) {
            Mask = "";
        }
        else if (
            mask.type() == variant::STRING
            || mask.type() == variant::STRINGVEC
        ) {
            Mask = mask.toString();
        }
        std::shared_ptr<Record> pRegion(_getRegion(region, false));
        ImageFactory::toASCII(
            _imageF, outfile, *pRegion, Mask,
            sep, format, maskvalue, overwrite, stretch
        );
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

bool image::tofits(
    const string& fitsfile, bool velocity,
    bool optical, int bitpix, double minpix,
    double maxpix, const variant& region,
    const variant& vmask, bool overwrite,
    bool dropdeg, bool deglast, bool dropstokes,
    bool stokeslast, bool wavelength, bool airwavelength,
    bool stretch, bool history
) {
    _log << _ORIGIN;
    if (_detached()) {
        return false;
    }
    try {
        _notSupported(__func__);
        ThrowIf(
            fitsfile.empty(),
            "fitsfile must be specified"
        );
        ThrowIf(
            fitsfile == "." || fitsfile == "..",
            "Invalid fitsfile name " + fitsfile
        );
        std::shared_ptr<Record> pRegion(_getRegion(region, false));
        String mask = vmask.toString();
        if (mask == "[]") {
            mask = "";
        }
        String origin;
        {
            ostringstream buffer;
            buffer << "CASA ";
            VersionInfo::report(buffer);
            origin = String(buffer);

	    // sanitize: replace CR and LF by SPACE
	    const Char *cOrigin = origin.chars();
	    for(String::size_type i=0; i<origin.length(); i++){
	      if(cOrigin[i]==10 || cOrigin[i]==13){
		origin.at(i,(String::size_type)1) = " ";
	      }
	    }
	    origin.rtrim(' ');
	    
        }
        ThrowIf(
            ! _imageF,
            "Only writing float-valued images to FITS is supported"
        );
        ImageFactory::toFITS(
            _imageF, fitsfile, velocity, optical,
            bitpix, minpix, maxpix, *pRegion, mask, overwrite,
            dropdeg, deglast, dropstokes, stokeslast, wavelength,
            airwavelength, origin, stretch, history
        );
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

record* image::topixel(const variant& value) {
    try {
        _log << LogOrigin("image", __func__);
        if (_detached()) {
            return nullptr;
        }
        //NOT using _image->toworld as most of the math is in casac namespace
        //in coordsys...should revisit this when casac::coordsys is cleaned
        return coordsys()->topixel(value);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

record* image::torecord() {
    _log << LogOrigin("image", __func__);
    if (_detached()) {
        return new record();
    }
    try {
        Record rec;
        String err;
        Bool ret = false;
        if (_imageF) {
            ret = _imageF->toRecord(err, rec);
        }
        else if (_imageC) {
            ret = _imageC->toRecord(err, rec);
        }
        else if (_imageD) {
            ret = _imageD->toRecord(err, rec);
        }
        else if (_imageDC) {
            ret = _imageDC->toRecord(err, rec);
        }
        ThrowIf (! ret, "Could not convert to record: " + err);
        return fromRecord(rec);
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
            << LogIO::POST;
        RETHROW(x);
    }
    return new record();
}

record* image::toworld(
    const variant& value, const string& format, bool dovelocity
) {
    try {
        _log << LogOrigin("image", __func__);
        if (_detached()) {
            return nullptr;
        }
        if (_imageF) {
            return _toworld(_imageF, value, format, dovelocity);
        }
        else if (_imageC) {
            return _toworld(_imageC, value, format, dovelocity);
        }
        else if (_imageD) {
            return _toworld(_imageD, value, format, dovelocity);
        }
        else if (_imageDC) {
            return _toworld(_imageDC, value, format, dovelocity);
        }
        ThrowCc("Logic error");
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
             << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

template <class T> record* image::_toworld(
    SPIIT image, const variant& value,
    const string& format, bool dovelocity
) {
    Vector<Double> pixel;
    if (_isUnset(value)) {
        pixel.resize(0);
    }
    else if (value.type() == variant::DOUBLEVEC) {
        pixel = value.getDoubleVec();
    }
    else if (value.type() == variant::INTVEC) {
        variant vcopy = value;
        Vector<Int> ipixel = vcopy.asIntVec();
        Int n = ipixel.size();
        pixel.resize(n);
        for (int i = 0; i < n; ++i) {
            pixel[i] = ipixel[i];
        }
    }
    else if (value.type() == ::casac::variant::RECORD) {
        variant localvar(value);
        unique_ptr<Record> tmp(toRecord(localvar.asRecord()));
        if (tmp->isDefined("numeric")) {
            pixel = tmp->asArrayDouble("numeric");
        }
        else {
            ThrowCc("Unsupported record type for value");
        }
    }
    else {
        ThrowCc("Unsupported data type for value");
    }
    ImageMetaData<T> md(image);
    return fromRecord(md.toWorld(pixel, format, dovelocity));
}

image* image::transpose(
    const std::string& outfile, const variant& order
) {
    try {
        _log << _ORIGIN;
        if (_detached()) {
            throw AipsError("No image specified to transpose");
        }
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        std::unique_ptr<ImageTransposer> transposer;
        switch(order.type()) {
        case variant::INT:
            transposer.reset(
                new ImageTransposer(
                    _imageF,
                    order.toInt(), outfile
                )
            );
            break;
        case variant::STRING:
            transposer.reset(
                new ImageTransposer(
                    _imageF,
                    order.toString(), outfile
                )
            );
            break;
        case variant::STRINGVEC:
            {
                Vector<String> orderVec = toVectorString(order.toStringVec());
                transposer.reset(
                    new ImageTransposer(
                        _imageF, orderVec,
                        outfile
                    )
                );
            }
            break;
        default:
            ThrowCc(
                "Unsupported type for order parameter " + order.typeString()
                + ". Supported types are a non-negative integer, a single "
                "string containing all digits or a list of strings which "
                "unambiguously match the image axis names."
            );
            break;
        }
        if (_doHistory) {
            vector<String> names = {"outfile", "order"};
            vector<variant> values = {outfile, order};
            auto msgs = _newHistory(__func__, names, values);
            transposer->addHistory(_ORIGIN, msgs);
        }
        return new image(
            transposer->transpose()
        );
    }
    catch (const AipsError& x) {
        _log << "Exception Reported: " << x.getMesg()
           << LogIO::EXCEPTION;
        RETHROW(x);
    }
    return nullptr;
}

bool image::twopointcorrelation(
    const string& outfile,
    const variant& region, const variant& vmask,
    const vector<int>& axes, const string& method,
    bool overwrite, bool stretch
) {
    _log << _ORIGIN;
    if (_detached()) {
        return false;
    }
    try {
        _notSupported(__func__);
        String outFile(outfile);
        std::shared_ptr<Record> Region(_getRegion(region, false));
        String mask = vmask.toString();
        if (mask == "[]") {
            mask = "";
        }
        Vector<Int> iAxes;
        if (!(axes.size() == 1 && axes[0] == -1)) {
            iAxes = axes;
        }
        vector<String> msgs;
        if (_doHistory) {
            vector<String> names {
                "outfile", "region", "mask", "axes",
                "method", "overwrite", "stretch"
            };
            vector<variant> values {
                outfile, region, vmask, axes,
                method, overwrite, stretch
            };
            msgs = _newHistory(__func__, names, values);
        }
        if (_imageF) {
            auto im = _twopointcorrelation(
                _imageF, outfile, Region, mask,
                IPosition(iAxes), method, overwrite, stretch,
                _ORIGIN, msgs
            );
        }
        else {
            auto im = _twopointcorrelation(
                _imageC, outfile, Region, mask,
                IPosition(iAxes), method, overwrite, stretch,
                _ORIGIN, msgs
            );
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

template <class T> SPIIT image::_twopointcorrelation(
    SPIIT myimage, const std::string& outfile,
    std::shared_ptr<Record> region, const casacore::String& mask,
    const IPosition& axes, const std::string& method,
    bool overwrite, bool stretch, const LogOrigin& origin,
    const vector<String>& msgs
) const {
    TwoPointCorrelator<T> tpc(
        myimage, region.get(), mask, outfile, overwrite
    );
    tpc.setAxes(axes);
    tpc.setMethod(method);
    tpc.setStretch(stretch);
    if (_doHistory) {
        tpc.addHistory(origin, msgs);
    }
    return tpc.correlate();
}

string image::type() {
    return "image";
}

bool image::unlock() {
    try {
        _log << LogOrigin("image", __func__);
        if (_detached()) {
            return false;
        }
        _notSupported(__func__);
        if (_imageF) {
            _imageF->unlock();
        }
        else {
            _imageC->unlock();
        }
        return true;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return false;
}

template <class T> void image::_addHistory(
    SPIIT image, const casacore::String& method, const vector<String>& names,
    const std::vector<variant>& values,
    const std::vector<casacore::String>& appendMsgs,
    const std::set<casacore::String>& dontQuote
) {
    if (! _doHistory) {
        // history writing disabled
        return;
    }
    auto msgs = _newHistory(method, names, values, dontQuote);
    for (const auto& m: appendMsgs) {
        msgs.push_back(m);
    }
    ImageHistory<T> ih(image);
    ih.addHistory("image::" + method, msgs);
}

void image::_addHistory(
    const casacore::String& method, const vector<casacore::String>& names,
    const std::vector<casac::variant>& values,
    const vector<casacore::String>& appendMsgs,
    const std::set<casacore::String>& dontQuote
) {
    if (! _doHistory) {
        // history writing disabled
        return;
    }
    if (_imageF) {
        _addHistory(
            _imageF, method, names, values, appendMsgs, dontQuote
        );
    }
    else if (_imageC) {
        _addHistory(
            _imageC, method, names, values, appendMsgs, dontQuote
        );
    }
    else if (_imageD) {
        _addHistory(
            _imageD, method, names, values, appendMsgs, dontQuote
        );
    }
    else if (_imageDC) {
        _addHistory(
            _imageDC, method, names, values, appendMsgs, dontQuote
        );
    }
    else {
        ThrowCc("Logic error");
    }
}

String image::_getMask(const variant& mask) {
    if (mask.type() == variant::BOOLVEC) {
       return "";
    }
    else if (mask.type() == variant::STRING) {
        return mask.toString();
    }
    else {
        ThrowCc("Mask is not understood, try a valid LEL string");
    }
}

String image::_inputsString(
    const vector<std::pair<String, variant> >& inputs,
    const std::set<String>& dontQuote
) {
    String out = "(";
    String quote;
    vector<pair<String, variant> >::const_iterator begin = inputs.begin();
    vector<pair<String, variant> >::const_iterator iter = begin;
    vector<pair<String, variant> >::const_iterator end = inputs.end();
    while (iter != end) {
        if (iter != begin) {
            out += ", ";
        }
        quote = iter->second.type() == variant::STRING
            && std::find(dontQuote.begin(), dontQuote.end(), iter->first) == dontQuote.end()
            ? "\"" : "";
        out += iter->first + "=" + quote;
        auto asString = iter->second.toString();
        if (asString.size() > 300) {
            asString = "...";
        }
        out += asString;
        out += quote;
        ++iter;
    }
    out += ")";
    return out;
}

Bool image::_isUnset(const variant &theVar) {
    return theVar.type() == variant::BOOLVEC
        && theVar.size() == 0;
}

vector<String> image::_newHistory(
    const string& method, const vector<String>& names,
    const vector<variant>& values, const std::set<String>& dontQuote
) {
    AlwaysAssert(names.size() == values.size(), AipsError);
    vector<String> msgs;
    ostringstream os;
    os << "Ran ia." << method;
    msgs.push_back(os.str());
    vector<std::pair<String, variant> > inputs;
    for (uInt i=0; i<names.size(); ++i) {
        inputs.push_back(make_pair(names[i], values[i]));
    }
    os.str("");
    os << "ia." << method << _inputsString(inputs, dontQuote);
    msgs.push_back(os.str());
    return msgs;
}

void image::_reset() {
    _imageF.reset();
    _imageC.reset();
    _imageD.reset();
    _imageDC.reset();
    _statsF.reset();
    _statsD.reset();
}

image* image::newimagefromfits(
    const string& outfile, const string& fitsfile,
    const int whichrep, const int whichhdu,
    const bool zeroBlanks, const bool overwrite
) {
    try {
        auto im = ImageFactory::fromFITS(
            outfile, fitsfile, whichrep, whichhdu,
            zeroBlanks, overwrite
        );
        if (im) {
            return new image(im);
        }
        else {
            ThrowCc("Unable to create image");
        }
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return nullptr;
}

variant* image::makearray(
    double v, const vector<int>& shape
) {
    auto ndim = shape.size();
    int nelem = 1;
    for (uInt i = 0; i < ndim; ++i) {
        nelem *= shape[i];
    }
    std::vector<double> element(nelem);
    for (int i = 0; i < nelem; ++i) {
        element[i] = v;
    }
    return new variant(element, shape);
}

void image::_notSupported(const std::string& method) const {
    ThrowIf(
        _imageD,
        method + " does not support images with double precision pixel values"
    );
    ThrowIf(
        _imageDC,
        method + " does not support images with complex<double> precision pixel values"
    );
}

casac::record*
image::recordFromQuantity(const casacore::Quantity q) {
    ::casac::record *r = 0;
    try {
        _log << LogOrigin("image", "recordFromQuantity");
        String error;
        casacore::Record R;
        if (QuantumHolder(q).toRecord(error, R)) {
            r = fromRecord(R);
        }
        else {
            _log << LogIO::SEVERE << "Could not convert quantity to record."
                    << LogIO::POST;
        }
    } catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return r;
}

::casac::record*
image::recordFromQuantity(const Quantum<Vector<Double> >& q) {
    ::casac::record *r = 0;
    try {
        _log << LogOrigin("image", "recordFromQuantity");
        String error;
        casacore::Record R;
        if (QuantumHolder(q).toRecord(error, R)) {
            r = fromRecord(R);
        } else {
            _log << LogIO::SEVERE << "Could not convert quantity to record."
                    << LogIO::POST;
        }
    } catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
    return r;
}

casacore::Quantity image::_casaQuantityFromVar(const ::casac::variant& theVar) {
    try {
        _log << _ORIGIN;
        casacore::QuantumHolder qh;
        String error;
        if (
            theVar.type() == ::casac::variant::STRING
            || theVar.type() == ::casac::variant::STRINGVEC
        ) {
            ThrowIf(
                !qh.fromString(error, theVar.toString()),
                "Error " + error + " in converting quantity "
            );
        }
        else if (theVar.type() == ::casac::variant::RECORD) {
            //NOW the record has to be compatible with QuantumHolder::toRecord
            ::casac::variant localvar(theVar); //cause its const
            unique_ptr<Record> ptrRec(toRecord(localvar.asRecord()));
            ThrowIf(
                !qh.fromRecord(error, *ptrRec),
                "Error " + error + " in converting quantity "
            );
        }
        else if (theVar.type() == variant::BOOLVEC) {
            return casacore::Quantity();
        }
        return qh.asQuantity();
    }
    catch (const AipsError& x) {
        _log << LogOrigin("image", __func__);
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
}

bool image::isconform(const string& other) {
    _log << _ORIGIN;

    if (_detached()) {
        return false;
    }
    try {
        ThrowIf(
            ! _imageF,
            "This method only supports Float valued images"
        );
        SPCIIF imageF;
        SPCIIC imageC;
        std::tie(imageF, imageC, std::ignore, std::ignore)
            = ImageFactory::fromFile(other);
        ThrowIf(! (imageF || imageC), "Unsupported image pixel data type");
        auto shape = imageF ? imageF->shape()
            : imageC->shape();
        const auto& csys = imageF ? imageF->coordinates()
                : imageC->coordinates();
        std::shared_ptr<const ImageInterface<Float> > mine = _imageF;
        if (
            _imageF->shape().isEqual(shape)
            && _imageF->coordinates().near(csys)
        ) {
            auto axisnames = csys.worldAxisNames();
            Vector<String> mc = _imageF->coordinates().worldAxisNames();
            if (mc.size() != axisnames.size()) {
                return false;
            }
            for (uInt i = 0; i < mc.size(); ++i) {
                if (mc[i] != axisnames[i]) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
    catch (const AipsError& x) {
        _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                << LogIO::POST;
        RETHROW(x);
    }
}

std::shared_ptr<Record> image::_getRegion(
    const ::casac::variant& region, const bool nullIfEmpty,
    const string& otherImageName
) const {
    switch (region.type()) {
    case variant::BOOLVEC:
        return std::shared_ptr<Record>(nullIfEmpty ? 0 : new Record());
    case variant::STRING: {
        IPosition shape;
        CoordinateSystem csys;
        if (otherImageName.empty()) {
            ThrowIf(
                ! (_imageF || _imageC || _imageD || _imageDC),
                "No attached image. Cannot use a string value for region"
            );
            if (_imageF) {
                shape = _imageF->shape();
                csys = _imageF->coordinates();
            }
            else if (_imageC) {
                shape = _imageC->shape();
                csys = _imageC->coordinates();
            }
            else if (_imageD) {
                shape = _imageD->shape();
                csys = _imageD->coordinates();
            }
            else if (_imageDC) {
                shape = _imageDC->shape();
                csys = _imageDC->coordinates();
            }
            else {
                ThrowCc("Logic Error");
            }
        }
        else {
            PtrHolder<ImageInterface<Float> > image;
            ImageUtilities::openImage(image, otherImageName);
            if (image.ptr()) {
                shape = image->shape();
                csys = image->coordinates();
            }
            else {
                PtrHolder<ImageInterface<Complex> > imagec;
                ImageUtilities::openImage(imagec, otherImageName);
                ThrowIf(
                    ! imagec.ptr(),
                    "Unable to open image " + otherImageName
                );
                shape = imagec->shape();
                csys = imagec->coordinates();
            }
        }
        return std::shared_ptr<Record>(
            region.toString().empty()
                ? nullIfEmpty ? 0 : new Record()
                : new Record(
                    CasacRegionManager::regionFromString(
                        csys, region.toString(),
                        _name(false), shape, True
                    )
                )
        );
    }
    case variant::RECORD:
        {
            std::shared_ptr<variant> clon(region.clone());
            return std::shared_ptr<casacore::Record>(
                nullIfEmpty && region.size() == 0
                    ? nullptr
                    : toRecord(
                        std::shared_ptr<::casac::variant>(region.clone())->asRecord()
                    )
            );
        }
    default:
        ThrowCc("Unsupported type for region " + region.typeString());
    }
}

void image::_setImage(casa::ITUPLE mytuple) {
    auto imageF = std::get<0>(mytuple);
    auto imageC = std::get<1>(mytuple);
    auto imageD = std::get<2>(mytuple);
    auto imageDC = std::get<3>(mytuple);
    uInt n = 0;
    if (imageF) {
        ++n;
    }
    if (imageC) {
        ++n;
    }
    if (imageD) {
        ++n;
    }
    if (imageDC) {
        ++n;
    }
    ThrowIf(n == 0, "No image defined");
    ThrowIf(
        n > 1,
        "Multiple images (" + String::toString(n)
        + ") defined"
    );
    _reset();
    _imageF = imageF;
    _imageC = imageC;
    _imageD = imageD;
    _imageDC = imageDC;
}


vector<double> image::_toDoubleVec(const variant& v) {
    variant::TYPE type = v.type();
    ThrowIf(
        type != variant::INTVEC && type != variant::LONGVEC
        && type != variant::DOUBLEVEC,
        "variant is not a numeric array"
    );
    vector<double> output;
    if (type == variant::INTVEC || type == variant::LONGVEC) {
        Vector<Int> x = v.toIntVec();
        std::copy(x.begin(), x.end(), std::back_inserter(output));
    }
    if (type == variant::DOUBLEVEC) {
        Vector<Double> x = v.toDoubleVec();
        std::copy(x.begin(), x.end(), std::back_inserter(output));
    }
    return output;
}

}
