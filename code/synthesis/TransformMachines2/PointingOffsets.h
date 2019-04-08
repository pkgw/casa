// -*- C++ -*-
//# PointingOffsets.h: Definition of the PointingOffsets class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//# $Id$
//
#ifndef SYNTHESIS_TRANSFORM2_POINTINGOFFSETS_H
#define SYNTHESIS_TRANSFORM2_POINTINGOFFSETS_H

#include <synthesis/MeasurementComponents/SolvableVisCal.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <images/Images/ImageInterface.h>
#include <msvis/MSVis/VisBufferUtil.h>
#include <images/Images/TempImage.h>
#include <msvis/MSVis/VisBuffer2.h>

class SolvableVisJones;
namespace casa { //# NAMESPACE CASA - BEGIN
  namespace refim{
  //
  //-------------------------------------------------------------------------------------------
  //
  class PointingOffsets 
  {
  public:
    PointingOffsets(const int& convOversampling):epJ_p(), doPointing_p(false),vbUtils_p() {convOversampling_p = convOversampling;}
    ~PointingOffsets() {};
    PointingOffsets& operator=(const PointingOffsets& other);

    void setEPJones(SolvableVisJones* epJ) {epJ_p = epJ;};
    virtual casacore::Vector<casacore::Vector<casacore::Double> >findMosaicPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
									const vi::VisBuffer2& vb, const casacore::Bool& doPointing=false);

    virtual casacore::Vector<casacore::Vector<casacore::Double> >findAntennaPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
									 const vi::VisBuffer2& vb, const casacore::Bool& doPointing=true);

    virtual casacore::Vector<casacore::Vector<casacore::Double> >findPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
								  const vi::VisBuffer2& vb, const casacore::Bool doPointing=false);

    casacore::Vector<double> gradPerPixel(const casacore::Vector<double>& p);
    casacore::Vector<casacore::Double>& toPix(const vi::VisBuffer2& vb, 
					      const casacore::MDirection& dir1, const casacore::MDirection& dir2);
    void storeImageParams(const casacore::ImageInterface<casacore::Complex>& iimage, const vi::VisBuffer2& vb);

    void setDoPointing(const bool& dop=false) {doPointing_p = dop;}

  private:

    casacore::Vector<double> thePix_p, pixFieldGrad_p;
    casacore::DirectionCoordinate imageDC_p;
    casacore::ObsInfo imageObsInfo_p;
    int nx_p, ny_p, nchan_p, npol_p, convOversampling_p, directionIndex_p;
    casacore::CoordinateSystem csys_p;
    casacore::DirectionCoordinate dc_p;
    casacore::MDirection::Convert pointToPix_p;
    casacore::MeasFrame pointFrame_p;
    casacore::MEpoch::Types timeMType_p;
    casacore::Unit timeUnit_p;
    casacore::MDirection direction1_p;
    casacore::MDirection direction2_p;
    casacore::CountedPtr<SolvableVisJones> epJ_p;
    bool doPointing_p;
    VisBufferUtil vbUtils_p;
  };
  //
  //-------------------------------------------------------------------------------------------
  //
};
};
#endif
