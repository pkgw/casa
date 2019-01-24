// -*- C++ -*-
//# PhaseGrad.cc: Implementation of the PhaseGrad class
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

#include <synthesis/TransformMachines2/PhaseGrad.h>
#include <synthesis/TransformMachines/SynthesisMath.h>
// #include <casa/Logging/LogIO.h>
// #include <casa/Logging/LogSink.h>
// #include <casa/Logging/LogOrigin.h>

using namespace casacore;
namespace casa{
  using namespace refim;
  //
  //----------------------------------------------------------------------
  //
  PhaseGrad& PhaseGrad::operator=(const PhaseGrad& other)
  {
    if(this!=&other) 
      {
	field_phaseGrad_p = other.field_phaseGrad_p;
	antenna_phaseGrad_p = other.antenna_phaseGrad_p;
	cached_FieldOffset_p[0] = other.cached_FieldOffset_p[0];
	cached_FieldOffset_p[1] = other.cached_FieldOffset_p[1];
	cached_AntennaOffset_p[0] = other.cached_AntennaOffset_p[0];
	cached_AntennaOffset_p[1] = other.cached_AntennaOffset_p[1];
      }
    return *this;
  }
  //
  //----------------------------------------------------------------------
  //
  void PhaseGrad::ComputeFieldPointingGrad(const Vector<double>& pointingOffset,
					   const Vector<int>&cfShape,
					   const Vector<int>& convOrigin,
					   const double& /*cfRefFreq*/,
					   const double& /*imRefFreq*/,
					   const int& spwID, const int& fieldId)
    {
      if (
	  ((fabs(pointingOffset[0]-cached_FieldOffset_p[0])) > 1e-6) ||
	  ((fabs(pointingOffset[0]-cached_FieldOffset_p[0])) > 1e-6) ||
	  (field_phaseGrad_p.shape()[0] < cfShape[0])              ||
	  (field_phaseGrad_p.shape()[1] < cfShape[1])
	  )
	{
	  // LogIO log_l(LogOrigin("AWVisResampler","cachePhaseGrad[R&D]"));
	  // log_l << "Computing phase gradiant for pointing offset " 
	  //       << pointingOffset << cfShape << " " << field_phaseGrad_p.shape() 
	  //       << "(SPW: " << spwID << " Field: " << fieldId << ")"
	  //       << LogIO::DEBUGGING
	  //       << LogIO::POST;
	  int nx=cfShape(0), ny=cfShape(1);
	  double grad;
	  Complex phx,phy;

	  field_phaseGrad_p.resize(nx,ny);
	  cached_FieldOffset_p[0] = pointingOffset[0];
	  cached_FieldOffset_p[1] = pointingOffset[1];
	
	  for(int ix=0;ix<nx;ix++)
	    {
	      grad = (ix-convOrigin[0])*pointingOffset[0];
	      double sx,cx;
	      SINCOS(grad,sx,cx);
	      phx = Complex(cx,sx);
	      for(int iy=0;iy<ny;iy++)
		{
		  grad = (iy-convOrigin[1])*pointingOffset[1];
		  Double sy,cy;
		  SINCOS(grad,sy,cy);
		  phy = Complex(cy,sy);
		  field_phaseGrad_p(ix,iy)=phx*phy;
		}
	    }
	}
    }
}
