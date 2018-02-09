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

#include <imageanalysis/ImageAnalysis/ImageCollapser.h>

namespace casa {

template<> void ImageCollapser<Float>::_doHighPerf(
    SPCIIF image, casacore::TempImage<Float>& outImage
) const {
    auto doMedian = _aggType == ImageCollapserData::MEDIAN;
    auto doMADM = _aggType == ImageCollapserData::MADM
        || _aggType == ImageCollapserData::XMADM;
    ThrowIf(
        ! doMedian && ! doMADM,
        "Logic error, unsupported aggregate type "
        + String(ImageCollapserData::funcNameMap()->at((uInt)_aggType)) + " for method "
        + String(__func__)
    );
    IPosition cursorShape(image->ndim(), 1);
    for (uInt i = 0; i < cursorShape.size(); ++i) {
        for (uInt j = 0; j < _axes.size(); ++j) {
            if (_axes[j] == i) {
                cursorShape[i] = image->shape()[i];
                break;
            }
        }
    }
    LatticeStepper stepper(image->shape(), cursorShape);
    std::unique_ptr<Array<Bool>> outMask;
    // accumtype being the same precision as the input data type is ok here,
    // since we are only computing the median/madm and not actually accumulating
    ClassicalStatistics<
        Double, Array<Float>::const_iterator, Array<Bool>::const_iterator
    > stats;
    auto hasMaskedPixels = ! ImageMask::isAllMaskTrue(*image);
    for (stepper.reset(); !stepper.atEnd(); stepper++) {
        Slicer slicer(
            stepper.position(), stepper.endPosition(), casacore::Slicer::endIsLast
        );
        auto data = image->getSlice(slicer);
        Bool isMasked = False;
        Array<Bool> maskSlice;
        if (hasMaskedPixels) {
            maskSlice = image->getMaskSlice(slicer);
            isMasked = ! allTrue(maskSlice);
        }
        if (isMasked) {
            if (! anyTrue(maskSlice)) {
                if (! outMask) {
                    outMask.reset(new Array<Bool>(outImage.shape(), true));
                }
                (*outMask)(stepper.position()) = false;
                outImage.putAt(0, stepper.position());
            }
            else if (! allTrue(maskSlice)) {
                stats.setData(data.begin(), maskSlice.begin(), data.size());
                if (doMedian) {
                    outImage.putAt(stats.getMedian(), stepper.position());
                }
                else if (doMADM) {
                    auto x = stats.getMedianAbsDevMed();
                    if (_aggType == ImageCollapserData::XMADM) {
                        x *= C::PROBIT_3_4;
                    }
                    outImage.putAt(x, stepper.position());
                }
            }
        }
        else {
            stats.setData(data.begin(), data.size());
            if (doMedian) {
                outImage.putAt(stats.getMedian(), stepper.position());
            }
            else if (doMADM) {
                auto x = stats.getMedianAbsDevMed();
                if (_aggType == ImageCollapserData::XMADM) {
                    x *= C::PROBIT_3_4;
                }
                outImage.putAt(x, stepper.position());
            }
        }
    }
    if (outMask) {
        outImage.attachMask(ArrayLattice<Bool>(*outMask));
    }
}
        
template<> void ImageCollapser<std::complex<float>>::_doHighPerf(
    SHARED_PTR<const ImageInterface<std::complex<float>>>, casacore::TempImage<std::complex<float>>&
) const {
    ThrowCc("Logic error: This version of the method should never be called");
}

}
