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
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogOrigin.h>

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
	cached_FieldOffset_p = other.cached_FieldOffset_p;
      }
    return *this;
  }
  //
  //----------------------------------------------------------------------
  //
  bool PhaseGrad::needsNewPhaseGrad(const Vector<Vector<double> >& pointingOffset,
				    const VisBuffer2& vb,
				    const int& row)
  {
    unsigned int nRow=vb.nRows();
    if (cached_FieldOffset_p.nelements() < nRow) cached_FieldOffset_p.resize(nRow,true);

    return (
	    ((fabs(pointingOffset[row][0]-cached_FieldOffset_p[row](0))) > 1e-12) ||
	    ((fabs(pointingOffset[row][1]-cached_FieldOffset_p[row](1))) > 1e-12) ||
	    (field_phaseGrad_p.shape()[0] < maxCFShape_p[0])           ||
	    (field_phaseGrad_p.shape()[1] < maxCFShape_p[1])
	    );
  }
  //
  //----------------------------------------------------------------------
  //
  // bool PhaseGrad::ComputeFieldPointingGrad(const Vector<double>& pointingOffset,
  // 					   const CountedPtr<CFBuffer>& cfb,
  // 					   const Vector<int>&cfShape,
  // 					   const Vector<int>& convOrigin,
  // 					   const double& /*cfRefFreq*/,
  // 					   const double& /*imRefFreq*/,
  // 					   const int& spwID, const int& fieldId)
  bool PhaseGrad::ComputeFieldPointingGrad(const Vector<Vector<double> >& pointingOffset,
					   const CountedPtr<CFBuffer>& cfb,
					   const VisBuffer2& vb,
					   const int& row
					   )

    {
      //
      // Re-find the max. CF size if the CFB changed.
      //
      CFBuffer *thisCFB = cfb.get();
      if (thisCFB != cachedCFBPtr_p)
	{
	  maxCFShape_p[0] = maxCFShape_p[1] = cfb->getMaxCFSize();
	  {
	    // LogIO log_l(LogOrigin("PhaseGrad","computeFieldPointingGrad"));
	    //cerr << "CFB changed: "<< thisCFB << " " << cachedCFBPtr_p << " " << vb.spectralWindows()(0) << " " << vb.fieldId()(0) << " " << maxCFShape_p << endl;
	  }
	  cachedCFBPtr_p = thisCFB;
	}
      //
      // If the pointing or the max. CF size changed, recompute the phase gradient.
      //
      if (needsNewPhaseGrad(pointingOffset, vb, row))
	{
	  LogIO log_l(LogOrigin("PhaseGrad","computeFieldPointingGrad"));
	  log_l << "Computing Phase Grad: " << row << " " << pointingOffset[row][0] << " " << pointingOffset[row][1] << " " << cached_FieldOffset_p[row](0) << " " 
		<< cached_FieldOffset_p[row](1) << " " << field_phaseGrad_p.shape() << " " << maxCFShape_p[0]  << " "
		<< cached_FieldOffset_p.nelements() 
		<< LogIO::POST;

	  int nx=maxCFShape_p(0), ny=maxCFShape_p(1);
	  double grad;
	  Complex phx,phy;
	  Vector<int> convOrigin = maxCFShape_p/2;
	  
	  field_phaseGrad_p.resize(nx,ny);
	  cached_FieldOffset_p[row] = pointingOffset[row];
	  // cached_FieldOffset_p[row](0) = pointingOffset[row][0];
	  // cached_FieldOffset_p[row](1) = pointingOffset[row][1];
	  
	  for(int ix=0;ix<nx;ix++)
	    {
	      grad = (ix-convOrigin[0])*pointingOffset[row][0];
	      double sx,cx;
	      SINCOS(grad,sx,cx);
	      phx = Complex(cx,sx);
	      for(int iy=0;iy<ny;iy++)
		{
		  grad = (iy-convOrigin[1])*pointingOffset[row][1];
		  Double sy,cy;
		  SINCOS(grad,sy,cy);
		  phy = Complex(cy,sy);
		  field_phaseGrad_p(ix,iy)=phx*phy;
		}
	    }
	  return true; // New phase gradient was computed
	}
      return false;
    }
}
