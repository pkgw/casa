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

#include <imageanalysis/ImageAnalysis/ComplexImageRegridder.h>

#include <imageanalysis/ImageAnalysis/ImageFactory.h>
#include <imageanalysis/ImageAnalysis/ImageRegridder.h>

namespace casa {

template<class T>
const String  ComplexImageRegridder<T>::_class = "ComplexImageRegridder";

template<class T> ComplexImageRegridder<T>::ComplexImageRegridder(
	SPCIIT image, const casacore::Record *const regionRec,
	const casacore::String& maskInp, const casacore::String& outname,
	casacore::Bool overwrite, const casacore::CoordinateSystem& csysTo,
	const casacore::IPosition& axes, const casacore::IPosition& shape
) : ImageRegridderBase<T>(
		image, regionRec,
		maskInp, outname, overwrite,
		csysTo, axes, shape
	) {}

template<class T> template<class U>
ComplexImageRegridder<T>::ComplexImageRegridder(
	SPCIIT image, const casacore::String& outname, SPCIIU templateIm,
	const casacore::IPosition& axes, const casacore::Record *const regionRec,
	const casacore::String& maskInp, casacore::Bool overwrite,
	const casacore::IPosition& shape
)  : ImageRegridderBase<T>(
		image, regionRec, maskInp, outname, overwrite,
		templateIm->coordinates(), axes, shape
) {}

template<class T> ComplexImageRegridder<T>::~ComplexImageRegridder() {}

template<class T> SPIIT ComplexImageRegridder<T>::regrid() const {
	auto myimage = this->_getImage();
	auto realPart = ImageFactory::floatFromComplex(myimage, ImageFactory::REAL);
	ImageRegridder<typename casacore::NumericTraits<T>::BaseType> rgReal(
	    realPart, this->_getRegion(), this->_getMask(), "",
		false, this->_getTemplateCoords(), this->_getAxes(), this->_getShape()
	);
	rgReal.setConfiguration(*this);
	SHARED_PTR<
	    const casacore::ImageInterface<
	        typename casacore::NumericTraits<T>::BaseType
	    >
	> outReal = rgReal.regrid();
	auto imagPart = ImageFactory::floatFromComplex(
		myimage, ImageFactory::IMAG
	);
	ImageRegridder<typename casacore::NumericTraits<T>::BaseType> rgImag(
	    imagPart, this->_getRegion(), this->_getMask(), "",
		false, this->_getTemplateCoords(), this->_getAxes(), this->_getShape()
	);
	rgImag.setConfiguration(*this);
	SHARED_PTR<
	    const casacore::ImageInterface<
	        typename casacore::NumericTraits<T>::BaseType
	    >
	> outImag = rgImag.regrid();
	auto outImage = ImageFactory::makeComplexImage(outReal, outImag);
	return this->_prepareOutputImage(*outImage);
}

}
