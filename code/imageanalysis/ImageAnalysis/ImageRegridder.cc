//# Copyright (C) 1998,1999,2000,2001,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: $

#include <imageanalysis/ImageAnalysis/ImageRegridder.h>

#include <imageanalysis/ImageAnalysis/SubImageFactory.h>
#include <images/Images/ImageConcat.h>
#include <images/Images/ImageRegrid.h>
#include <images/Regions/WCBox.h>

#include <imageanalysis/ImageAnalysis/ImageMetaData.h>
#include <imageanalysis/Annotations/AnnRectBox.h>

#include <casa/Arrays/ArrayLogical.h>

#include <memory>

namespace casa {

const String ImageRegridder::_class = "ImageRegridder";

ImageRegridder::ImageRegridder(
	const ImageTask::shCImFloat image,
	const Record *const regionRec,
	const String& maskInp, const String& outname, Bool overwrite,
	const CoordinateSystem& csysTo, const IPosition& axes,
	const IPosition& shape, Bool dropdeg
) : ImageTask(
		image, "", regionRec, "", "", "",
		maskInp, outname, overwrite
	),
	_csysTo(csysTo), _axes(axes), _shape(shape), _dropdeg(dropdeg),
	_specAsVelocity(False), _doRefChange(True), _replicate(False),
	_forceRegrid(False), _debug(0), _decimate(10),
	_method(Interpolate2D::LINEAR), _outputStokes(0) {
	_construct();
	_finishConstruction();
}

ImageRegridder::ImageRegridder(
	const ImageTask::shCImFloat image, const String& outname,
	const ImageTask::shCImFloat templateIm, const IPosition& axes,
	const Record *const regionRec, const String& maskInp,
	Bool overwrite, Bool dropdeg, const IPosition& shape
)  : ImageTask(
		image, "", regionRec, "", "", "",
		maskInp, outname, overwrite
	),
	_csysTo(templateIm->coordinates()), _axes(axes), _shape(shape),
	_dropdeg(dropdeg), _specAsVelocity(False), _doRefChange(True),
	_replicate(False), _forceRegrid(False), _debug(0), _decimate(10),
	_method(Interpolate2D::LINEAR),_outputStokes(0) {
	_construct();
	_finishConstruction();
}

ImageRegridder::~ImageRegridder() {}

std::tr1::shared_ptr<ImageInterface<Float> > ImageRegridder::regrid(
	Bool wantReturn
) const {
	Bool regridByVel = False;
	if (
		_specAsVelocity && _getImage()->coordinates().hasSpectralAxis()
		&& _csysTo.hasSpectralAxis()
	) {
		if (_axes.size() == 0) {
			regridByVel = True;
		}
		else {
			Int specAxis = _getImage()->coordinates().spectralAxisNumber();
			for (uInt i=0; i<_axes.size(); i++) {
				if (_axes[i] == specAxis) {
					regridByVel = True;
					break;
				}
			}
		}
	}
	std::tr1::shared_ptr<ImageInterface<Float> > workIm;
	if (regridByVel) {
		workIm = _regridByVelocity();
	}
	else {
		workIm = _regrid();
	}
	std::tr1::shared_ptr<ImageInterface<Float> > outImage = _prepareOutputImage(*workIm);
	if (! wantReturn) {
		outImage.reset();
	}
	return outImage;
}

std::tr1::shared_ptr<ImageInterface<Float> > ImageRegridder::_regrid() const {
	std::auto_ptr<ImageInterface<Float> > clone(_getImage()->cloneII());
	SubImage<Float> subImage = SubImageFactory<Float>::createSubImage(
		*clone, *_getRegion(), _getMask(), _getLog().get(),
		False, AxesSpecifier(! _dropdeg), _getStretch()
	);
	*_getLog() << LogOrigin(_class, __FUNCTION__);
	if (! anyTrue(subImage.getMask())) {
		*_getLog() << "All selected pixels are masked" << LogIO::EXCEPTION;
	}
	clone.reset(0);
	const CoordinateSystem csysFrom = subImage.coordinates();
	CoordinateSystem csysTo = _csysTo;
	csysTo.setObsInfo(csysFrom.obsInfo());
	std::set<Coordinate::Type> coordsToRegrid;
	CoordinateSystem csys = ImageRegrid<Float>::makeCoordinateSystem(
		*_getLog(), coordsToRegrid, csysTo, csysFrom, _axes, subImage.shape()
	);
	if (csys.nPixelAxes() != _shape.nelements()) {
		*_getLog()
			<< "The number of pixel axes in the output shape and Coordinate System must be the same"
			<< LogIO::EXCEPTION;
	}
	std::tr1::shared_ptr<ImageInterface<Float> > workIm(
		new TempImage<Float>(_kludgedShape, csys)
	);
	workIm->set(0.0);
	ImageUtilities::copyMiscellaneous(*workIm, subImage);
	String maskName("");
	ImageMaskAttacher<Float>::makeMask(*workIm, maskName, True, True, *_getLog(), True);
	if (csysFrom.hasDirectionCoordinate() && csys.hasDirectionCoordinate()) {
		ThrowIf (
			! _doImagesOverlap(subImage, *workIm),
			"There is no overlap between the (region chosen in) the input image"
			" and the output image with respect to the axes being regridded."
		);
	}
	ImageRegrid<Float> ir;
	ir.showDebugInfo(_debug);
	ir.disableReferenceConversions(!_doRefChange);
	ir.regrid(
		*workIm, _method, _axes, subImage,
		_replicate, _decimate, True,
		_forceRegrid
	);
	if (! _outputStokes.empty()) {
		ImageMetaData<Float> md(workIm.get());
		if (_outputStokes.size() < md.nStokes()) {
			CasacRegionManager rm(workIm->coordinates());
			String diagnostics;
			uInt nSelectedChannels = 0;
			if (_outputStokes.size() == 1) {
				/*
				String& diagnostics, uInt& nSelectedChannels, String& stokes,
					const String& chans,
					const StokesControl stokesControl, const String& box,
					const IPosition& imShape
					*/
				String stokes = _outputStokes[0];
				Record region = rm.fromBCS(
					diagnostics, nSelectedChannels, stokes,
					"", CasacRegionManager::USE_FIRST_STOKES,
					"", workIm->shape()
				).toRecord("");
				workIm = SubImageFactory<Float>::createImage(
					*workIm, "", region, "", False, False, False, False
				);
			}
			else {
				// Only include the wanted stokes
				std::tr1::shared_ptr<ImageConcat<Float> > concat(
					new ImageConcat<Float>(
						workIm->coordinates().polarizationAxisNumber(False)
					)
				);
				for (uInt i=0; i<_outputStokes.size(); i++) {
					String stokes = _outputStokes[i];
					Record region = rm.fromBCS(
						diagnostics, nSelectedChannels, stokes,
						"", CasacRegionManager::USE_FIRST_STOKES,
						"", workIm->shape()
					).toRecord("");
					concat->setImage(
						*SubImageFactory<Float>::createImage(
							*workIm, "", region, "", False, False, False, False
						), True
					);
				}
				workIm = concat;
			}
		}
	}
	ThrowIf(
		workIm->hasPixelMask() && ! anyTrue(workIm->pixelMask().get()),
		"All output pixels are masked."
	);
	return workIm;
}

std::tr1::shared_ptr<ImageInterface<Float> > ImageRegridder::_regridByVelocity() const {
	if (
		_csysTo.spectralCoordinate().frequencySystem(True)
		!= _getImage()->coordinates().spectralCoordinate().frequencySystem(True)
	) {
		*_getLog() << "Image to be regridded has different frequency system from template."
			<< LogIO::EXCEPTION;
	}
	std::auto_ptr<CoordinateSystem> csys(
		dynamic_cast<CoordinateSystem *>(_csysTo.clone())
	);
	SpectralCoordinate templateSpecCoord = csys->spectralCoordinate();
	std::auto_ptr<ImageInterface<Float> > clone(_getImage()->cloneII());
 	std::tr1::shared_ptr<SubImage<Float> > maskedClone(
		new SubImage<Float>(
			SubImageFactory<Float>::createSubImage(
				*clone, *_getRegion(), _getMask(), 0, False,
				AxesSpecifier(), _getStretch()
			)
		)
	);
 	clone.reset(0);
	std::auto_ptr<CoordinateSystem> coordClone(
		dynamic_cast<CoordinateSystem *>(maskedClone->coordinates().clone())
	);

	SpectralCoordinate newSpecCoord = coordClone->spectralCoordinate();

	Double newVelRefVal = 0;
	Double newVelInc = 0;
	for (uInt i=0; i<2; i++) {
		CoordinateSystem *cs = i == 0 ? csys.get() : coordClone.get();
		// create and replace the coordinate system's spectral coordinate with
		// a linear coordinate which describes the velocity axis. In this way
		// we can regrid by velocity.
		Int specCoordNum = cs->spectralCoordinateNumber();
		SpectralCoordinate specCoord = cs->spectralCoordinate();
		if (
			specCoord.frequencySystem(False) != specCoord.frequencySystem(True)
		) {
			// the underlying conversion system is different from the overlying one, so this
			// is pretty confusing. We want the underlying one also be the overlying one before
			// we regrid.
			Vector<Double> newRefVal;
			Double newRefPix = specCoord.referencePixel()[0];
			specCoord.toWorld(newRefVal, Vector<Double>(1, newRefPix));
			Vector<Double> newVal;
			specCoord.toWorld(newVal, Vector<Double>(1, newRefPix+1));

			specCoord = SpectralCoordinate(
				specCoord.frequencySystem(True), newRefVal[0],
				newVal[0] - newRefVal[0], newRefPix, specCoord.restFrequency()
			);
			if (cs == coordClone.get()) {
				newSpecCoord = specCoord;
			}
		}
		Double freqRefVal = specCoord.referenceValue()[0];
		Double velRefVal;
		if (! specCoord.frequencyToVelocity(velRefVal, freqRefVal)) {
			*_getLog() << "Unable to determine reference velocity" << LogIO::EXCEPTION;
		}
		Double vel0, vel1;
		if (
			! specCoord.pixelToVelocity(vel0, 0.0)
			|| ! specCoord.pixelToVelocity(vel1, 1.0)
		) {
			*_getLog() << "Unable to determine velocity increment" << LogIO::EXCEPTION;
		}
		Matrix<Double> pc(1, 1, 0);
		pc.diagonal() = 1.0;
		LinearCoordinate lin(
			Vector<String>(1, "velocity"),
			specCoord.worldAxisUnits(),
			Vector<Double>(1, velRefVal),
			Vector<Double>(1, vel1 - vel0),
			pc, specCoord.referencePixel()
		);
		if (
			! cs->replaceCoordinate(lin, specCoordNum)
			&& ! lin.near(cs->linearCoordinate(specCoordNum))) {
				*_getLog() << "Unable to replace spectral with linear coordinate"
				<< LogIO::EXCEPTION;
		}
		if (cs == csys.get()) {
			newVelRefVal = velRefVal;
			newVelInc = vel1 - vel0;
		}
		else {
			maskedClone->setCoordinateInfo(*cs);
		}
	}
	// do not pass the region or the mask, the maskedClone has already had the region and
	// mask applied
	ImageRegridder regridder(
		maskedClone, 0, "", _getOutname(),
		_getOverwrite(), *csys, _axes, _shape, _dropdeg
	);
	regridder.setMethod(_method);
	regridder.setDecimate(_decimate);
	regridder.setReplicate(_replicate);
	regridder.setDoRefChange(_doRefChange);
	regridder.setForceRegrid(_forceRegrid);
	regridder.setStretch(_getStretch());
	regridder.setSpecAsVelocity(False);

	std::tr1::shared_ptr<ImageInterface<Float> > outImage = regridder._regrid();

	// replace the temporary linear coordinate with the saved spectral coordinate
	std::auto_ptr<CoordinateSystem> newCoords(
		dynamic_cast<CoordinateSystem *>(outImage->coordinates().clone())
	);

	// make frequencies correct
	Double newRefFreq;
	if (
		! newSpecCoord.velocityToFrequency(
			newRefFreq, newVelRefVal
		)
	) {
		*_getLog() << "Unable to determine new reference frequency"
			<< LogIO::EXCEPTION;
	}
	// get the new frequency increment
	Double newFreq;
	if (
		! newSpecCoord.velocityToFrequency(
			newFreq, newVelRefVal + newVelInc
		)
	) {
		*_getLog() << "Unable to determine new frequency increment" << LogIO::EXCEPTION;
	}
	if (! newSpecCoord.setReferenceValue(Vector<Double>(1, newRefFreq))) {
		*_getLog() << "Unable to set new reference frequency" << LogIO::EXCEPTION;
	}
	if (! newSpecCoord.setIncrement((Vector<Double>(1, newFreq - newRefFreq)))) {
		*_getLog() << "Unable to set new frequency increment" << LogIO::EXCEPTION;
	}
	if (
		! newSpecCoord.setReferencePixel(
			templateSpecCoord.referencePixel()
		)
	) {
		*_getLog() << "Unable to set new reference pixel" << LogIO::EXCEPTION;
	}
	if (
		! newCoords->replaceCoordinate(
			newSpecCoord,
			maskedClone->coordinates().linearCoordinateNumber())
		&& ! newSpecCoord.near(newCoords->spectralCoordinate())
	) {
		*_getLog() << "Unable to replace coordinate for velocity regridding"
			<< LogIO::EXCEPTION;
	}
	outImage->setCoordinateInfo(*newCoords);
	return outImage;
}

void ImageRegridder::_finishConstruction() {
	Bool shapeSpecified = ! _shape.empty() && _shape[0] >= 0;
	if (! shapeSpecified) {
		IPosition imShape = _getImage()->shape();
		_shape.resize(imShape.size());
		_shape = imShape;
		if (_dropdeg) {
			IPosition tmp = _shape.nonDegenerate();
			_shape.resize(tmp.size());
			_shape = tmp;
		}
	}
	_kludgedShape = _shape;
	const CoordinateSystem csysFrom = _getImage()->coordinates();
	// enforce stokes rules CAS-4960
	if (csysFrom.hasPolarizationCoordinate() && _csysTo.hasPolarizationCoordinate()) {
		Vector<Int> templateStokes = _csysTo.stokesCoordinate().stokes();
		Vector<Int> inputStokes = csysFrom.stokesCoordinate().stokes();
		Int inputPolAxisNumber = csysFrom.polarizationAxisNumber();
		if (
			(
				_axes.empty()
				|| inputPolAxisNumber < (Int)_axes.size()
			) && templateStokes.size() > 1
		) {
			if (
				(
					_axes.empty() && inputStokes.size() > 1
				)
				|| _axes[inputPolAxisNumber] > 0
			) {
				StokesCoordinate stokesFrom = csysFrom.stokesCoordinate();
				StokesCoordinate stokesTo = _csysTo.stokesCoordinate();
				Stokes::StokesTypes valFrom, valTo;
				for (uInt i=0; i<inputStokes.size(); i++) {
					stokesFrom.toWorld(valFrom, i);
					for (uInt j=0; j<templateStokes.size(); j++) {
						stokesTo.toWorld(valTo, j);
						if (valFrom == valTo) {
							_outputStokes.push_back(Stokes::name(valFrom));
							break;
						}
					}
				}
				ThrowIf(
					_outputStokes.empty(),
					"Input image and template coordinate system have no common stokes."
				);
				ThrowIf(
					shapeSpecified && ((Int)_outputStokes.size() != _shape[inputPolAxisNumber]),
					"Specified output stokes axis length (" + String::toString(_shape[inputPolAxisNumber])
					+ ") does not match the number of common stokes ("
					+ String::toString(_outputStokes.size())
					+ ") in the input image and template coordinate system."
				);
				// This is a kludge to fool the underlying ImageRegrid constructor that the shape
				// is acceptable to it. We copy just the stokes we from the output of ImageRegrid.
				ImageMetaData<Float> md(_getImage().get());
				_kludgedShape[csysFrom.polarizationAxisNumber(False)] = md.nStokes();
			}
		}
	}
}

Bool ImageRegridder::_doImagesOverlap(
	const ImageInterface<Float>& image0,
	const ImageInterface<Float>& image1
) {
	const CoordinateSystem csys0 = image0.coordinates();
	const CoordinateSystem csys1 = image0.coordinates();
	IPosition shape0 = image0.shape();
	IPosition shape1 = image1.shape();
	if (
		csys0.hasDirectionCoordinate()
		&& csys1.hasDirectionCoordinate()
	) {
		const DirectionCoordinate dc0 = csys0.directionCoordinate();
		DirectionCoordinate dc1 = csys1.directionCoordinate();
		Bool sameFrame = dc0.directionType(True) == dc1.directionType(True);
		if (!sameFrame) {
			dc1.setReferenceConversion(dc0.directionType(True));
		}
		ImageMetaData<Float> md0(&image0);
		ImageMetaData<Float> md1(&image1);
		Vector<Int> dirShape0 = md0.directionShape();
		Vector<Int> dirShape1 = md1.directionShape();
		Vector<std::pair<Double, Double> > corners0 = _getDirectionCorners(
			dc0, dirShape0
		);
		Vector<std::pair<Double, Double> > corners1 = _getDirectionCorners(
			dc1, dirShape1
		);
		if (sameFrame) {
			Double minx0 = corners0[0].first;
			Double maxx0 = minx0;
			Double miny0 = corners0[0].second;
			Double maxy0 = miny0;
			Double minx1 = corners1[0].first;
			Double maxx1 = minx1;
			Double miny1 = corners1[0].second;
			Double maxy1 = miny1;

			for (uInt i=1; i<4; i++) {
				minx0 = min(minx0, corners0[i].first);
				maxx0 = max(maxx0, corners0[i].first);
				miny0 = min(miny0, corners0[i].second);
				maxy0 = max(maxy0, corners0[i].second);

				minx1 = min(minx1, corners1[i].first);
				maxx1 = max(maxx1, corners1[i].first);
				miny1 = min(miny1, corners1[i].second);
				maxy1 = max(maxy1, corners1[i].second);
			}
			if (
				minx0 > maxx1 || maxx0 < minx1
				|| miny0 > maxy1 || maxy0 < miny1
			) {
				return False;
			}
		}
		for (uInt i=0; i<4; i++) {
			Vector<Double> start0(2, corners0[i].first);
			start0[1] = corners0[i].second;
			Vector<Double> end0(
				2,
				i == 3 ? corners0[0].first
					: corners0[i+1].first
			);
			end0[1] = i == 3 ? corners0[0].second : corners0[i+1].second;

			for (uInt j=0; j<4; j++) {
				Vector<Double> start1(2, corners1[j].first);
				start1[1] = corners1[j].second;
				Vector<Double> end1(
					2,
					j == 3 ? corners1[0].first
						: corners1[j+1].first
				);
				end1[1] = j == 3 ? corners1[0].second : corners1[j+1].second;
				if (
					_doLineSegmentsIntersect(
						start0, end0, start1, end1
					)
				) {
					return True;
				}
			}
		}
	}
	return False;
}

Vector<std::pair<Double, Double> > ImageRegridder::_getDirectionCorners(
	const DirectionCoordinate& dc,
	const IPosition& directionShape
) {
	Vector<Double> world;
	Vector<Double> pixel(2, 0);
	Vector<String> units = dc.worldAxisUnits();
	dc.toWorld(world, pixel);
	Vector<std::pair<Double, Double> > corners(4);
	// blcx, blcy
	corners[0].first = Quantity(world[0], units[0]).getValue("rad");
	corners[0].second = Quantity(world[1], units[1]).getValue("rad");
	// trcx, blcy
	pixel[0] = directionShape[0];
	dc.toWorld(world, pixel);
	corners[1].first = Quantity(world[0], units[0]).getValue("rad");
	corners[1].second = Quantity(world[1], units[1]).getValue("rad");
	// trcx, trcy
	pixel[1] = directionShape[1];
	dc.toWorld(world, pixel);
	corners[2].first = Quantity(world[0], units[0]).getValue("rad");
	corners[2].second = Quantity(world[1], units[1]).getValue("rad");
	// blcx, trcy
	pixel[0] = 0;
	dc.toWorld(world, pixel);
	corners[3].first = Quantity(world[0], units[0]).getValue("rad");
	corners[3].second = Quantity(world[1], units[1]).getValue("rad");
	return corners;
}

Bool ImageRegridder::_doLineSegmentsIntersect(
	const Vector<Double>& line0point0,
	const Vector<Double>& line0point1,
	const Vector<Double>& line1point0,
	const Vector<Double>& line1point1
) {
	// algorithm from
	// http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
	Vector<Double> p = line0point0;
	Vector<Double> r = line0point1 - line0point0;
	Vector<Double> q = line1point0;
	Vector<Double> s = line1point1 - line1point0;
	Double rCrossS = _2DCrossProduct(r, s);
	Vector<Double> diffQP = q-p;

	if (rCrossS == 0) {
		if (_2DCrossProduct(diffQP, r) == 0) {
			// lines are coincident
			return True;
		}
		else {
			// lines are parallel
			return False;
		}
	}
	Double t = _2DCrossProduct(diffQP, s)/rCrossS;
	Double u = _2DCrossProduct(diffQP, r)/rCrossS;
	return t >= 0 && t <= 1 && u >= 0 && u <= 1;
}

Double ImageRegridder::_2DCrossProduct(
	const Vector<Double> v0, const Vector<Double> v1
) {
	return v0[0]*v1[1] - v0[1]*v1[0];
}

}


