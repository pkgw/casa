// -*- C++ -*-
//# AWConvFuncEPJones.cc: Implementation of the AWConvFuncEPJones class
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
#include <synthesis/TransformMachines2/AWConvFuncEPJones.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <images/Images/ImageInterface.h>

using namespace casacore;
namespace casa{
  using namespace vi;
  using namespace refim;
  //
  //----------------------------------------------------------------------
  //
  AWConvFuncEPJones& AWConvFuncEPJones::operator=(const AWConvFuncEPJones& other)
  {
    if(this!=&other) 
      {
	AWConvFunc::operator=(other);
	// imageDC_p = other.imageDC_p;
	// imageObsInfo_p = other.imageObsInfo_p;
      }
    return *this;
  }
  //
  //----------------------------------------------------------------------
  // Find the offset between the VB and the image phase center
  //
  Vector<Double> AWConvFuncEPJones::findPointingOffset(const ImageInterface<Complex>& image,
						       const VisBuffer2& vb)
  {
    //return po_p->findMosaicPointingOffset(image,vb);
    return po_p->findAntennaPointingOffset(image,vb);
  }
  //
  //----------------------------------------------------------------------
  //
  void AWConvFuncEPJones::makeConvFunction(const ImageInterface<Complex>& image,
					   const VisBuffer2& vb,
					   const Int wConvSize,
					   const CountedPtr<PolOuterProduct>& pop,
					   const Float pa,
					   const Float dpa,
					   const Vector<Double>& uvScale, const Vector<Double>& uvOffset,
					   const Matrix<Double>& spwFreqSel,
					   CFStore2& cfs,
					   CFStore2& cfwts,
					   Bool fillCF)
  {
    findPointingOffset(image,vb);
    AWConvFunc::makeConvFunction(image,vb,wConvSize,pop,pa,dpa,uvScale,uvOffset,spwFreqSel,cfs,cfwts,fillCF);
  }
}
