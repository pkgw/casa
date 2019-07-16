//-*- C++ -*-
//# VB2CFBMap.h: Definition of the VB2CFBMap class
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
#include <synthesis/TransformMachines/SynthesisMath.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogOrigin.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Array.h>
#include <casa/Utilities/CountedPtr.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <synthesis/TransformMachines2/ConvolutionFunction.h>
#include <synthesis/TransformMachines2/CFBuffer.h>
#include <synthesis/TransformMachines2/PhaseGrad.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/TransformMachines2/VB2CFBMap.h>
namespace casa{
using namespace vi;
  namespace refim{
  Int mapAntIDToAntType(const casacore::Int& /*ant*/) {return 0;};

    VB2CFBMap::VB2CFBMap(): vb2CFBMap_p(), cfPhaseGrad_p(), phaseGradCalculator_p(),doPointing_p(false)
    {
      phaseGradCalculator_p = new PhaseGrad();
      newPhaseGradComputed_p = false;
    };

    VB2CFBMap& VB2CFBMap::operator=(const VB2CFBMap& other)
    {
      if(this!=&other) 
	{
	  phaseGradCalculator_p = other.phaseGradCalculator_p;
	  cfPhaseGrad_p.assign(other.cfPhaseGrad_p);
	  vb2CFBMap_p.assign(other.vb2CFBMap_p);
	  doPointing_p = other.doPointing_p;
	}
      return *this;
    };

    void VB2CFBMap::setPhaseGradPerRow(const casacore::Vector<casacore::Vector<double> >& pointingOffset,
				       const casacore::CountedPtr<CFBuffer>& cfb,
				       const vi::VisBuffer2& vb,
				       const int& row)
    {
      // if (doPointing_p)
      // 	{
      // 	  if (phaseGradCalculator_p->needsNewPhaseGrad(pointingOffset, vb, 0))
      // 	    {
      // 	      phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
      // 	      newPhaseGradComputed_p=true;
      // 	    }
      // 	}
      // else
	{
	  phaseGradCalculator_p->ComputeFieldPointingGrad(pointingOffset,cfb,vb, 0);
	}
	{
	  //cfPhaseGrad_p(row).assign(phaseGradCalculator_p->getFieldPointingGrad());
	  cfPhaseGrad_p(row).reference(phaseGradCalculator_p->field_phaseGrad_p);
	}
    }
    
    Int VB2CFBMap::makeVBRow2CFBMap(CFStore2& cfs,
				    const VisBuffer2& vb, 
				    const Quantity& dPA,
				    const Vector<Int>& /*dataChan2ImChanMap*/,
				    const Vector<Int>& /*dataPol2ImPolMap*/,
				    const Vector<Vector<Double> >& pointingOffset)
				    //const Bool& /*doPointing*/)
    {
      //    VBRow2CFMapType& vbRow2CFMap_p,
      const Int nRow=vb.nRows(); 
      //UNUSED: nChan=dataChan2ImChanMap.nelements(), 
      //UNUSED: nPol=dataPol2ImPolMap.nelements();
      //    vbRow2CFMap_p.resize(nPol, nChan, nRow);
      vb2CFBMap_p.resize(nRow);
      cfPhaseGrad_p.resize(nRow);

      Quantity pa(getPA(vb),"rad");
      //PolOuterProduct outerProduct;
      Int statusCode=CFDefs::MEMCACHE;

      for (Int irow=0;irow<nRow;irow++)
	{
	  //
	  // Translate antenna ID to antenna type
	  //
	  Int ant1Type = mapAntIDToAntType(vb.antenna1()(irow)),
	    ant2Type = mapAntIDToAntType(vb.antenna2()(irow));
	  //
	  // Get the CFBuffer for the given PA and baseline catagorized
	  // by the two antenna types.  For homgeneous arrays, all
	  // baselines will map to a single antenna-type pair.
	  //
	  
	  CountedPtr<CFBuffer> cfb_l;
	  try
	    {
	      cfb_l = cfs.getCFBuffer(pa, dPA, ant1Type, ant2Type);
	      //cfb_l->show("From VRB: ");
	    }
	  catch (CFNotCached& x)
	    {
	      LogIO log_l(LogOrigin("VB2CFBMap", "makeVBRow2CFBMap[R&D]"));

	      log_l << "CFs not cached for " << pa.getValue("deg") 
		    << " deg, dPA = " << dPA.getValue("deg") 
		    << " Field ID = " << vb.fieldId()(0);
	      log_l << " Ant1Type, Ant2Type = " << ant1Type << "," << ant2Type << LogIO::POST;
	      statusCode=CFDefs::NOTCACHED;
	    }
	  
	  if (statusCode==CFDefs::NOTCACHED)
	    {
	      break;
	    }
	  else
	    {
	      // Set the phase grad for the CF per VB row
	      setPhaseGradPerRow(pointingOffset, cfb_l, vb, irow);

	      // Set the CFB per VB row
	      cfb_l->setPointingOffset(pointingOffset);
	      vb2CFBMap_p(irow) = cfb_l;
	    }
	}
      // {
      // 	double n=0;
      // 	for (int i=0;i<cfPhaseGrad_p.nelements();i++)
      // 	  n+=cfPhaseGrad_p[i].shape().product()*sizeof(casacore::Complex);
      // 	log_l << "Size of VB2CFBMap::cfPhaseGrad_p = " << n << " bytes" << LogIO::POST;
      // }
      return statusCode;
    }
}
}
