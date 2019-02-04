//-*- C++ -*-
//# PhaseGrad.h: Definition of the PhaseGrad class
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

  Int VB2CFBMap::makeVBRow2CFBMap(CFStore2& cfs,
					       ConvolutionFunction& cf,
					       const VisBuffer2& vbs, 
					       const Quantity& dPA,
					       const Vector<Int>& /*dataChan2ImChanMap*/,
					       const Vector<Int>& /*dataPol2ImPolMap*/,
					       const Vector<Double>& pointingOffset)
  {
    //    VBRow2CFMapType& vbRow2CFMap_p,
    const Int nRow=vbs.nRows(); 
    //UNUSED: nChan=dataChan2ImChanMap.nelements(), 
    //UNUSED: nPol=dataPol2ImPolMap.nelements();
    //    vbRow2CFMap_p.resize(nPol, nChan, nRow);
    vbRow2CFBMap_p.resize(nRow);
    Quantity pa(getPA(vbs),"rad");
    PolOuterProduct outerProduct;
    Int statusCode=CFDefs::MEMCACHE;
    for (Int irow=0;irow<nRow;irow++)
      {
	//
	// Translate antenna ID to antenna type
	//
	Int ant1Type = cf.mapAntIDToAntType(vbs.antenna1()(irow)),
	  ant2Type = cf.mapAntIDToAntType(vbs.antenna2()(irow));
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
	    LogIO log_l(LogOrigin("VisibilityResamplerBase2", "makeVBRow2CFMap"));
	    log_l << "CFs not cached for " << pa.getValue("deg") 
		  << " deg, dPA = " << dPA.getValue("deg") 
		  << " Field ID = " << vbs.fieldId()(0);
	      //		  << " TimeStamps(0-10) = " << vbs.feedPa1(getCurrentTimeStamp(vbs)).nelements() << " ";
	      //		  << " TimeStamps(0-10) = " << vbs.feedPa1().nelements() << " ";
	    // for (Int i=0;i<10;i++) 
	    //   {
	    // 	//		log_l << MVTime(vbs.time()(i)).string(MVTime::TIME) << " ";
	    // 	log_l << vbs.time()(i)/1.0e8 << " ";
	    // 	log_l << "(" << (vbs.feed_pa(getCurrentTimeStamp(vbs))(i))*57.2956 << ") ";
	    //   }
	    log_l << " Ant1Type, Ant2Type = " << ant1Type << "," << ant2Type << LogIO::POST;
	    statusCode=CFDefs::NOTCACHED;
	  }

	if (statusCode==CFDefs::NOTCACHED)
	  {
	    break;
	  }
	else
	  {
	    cfb_l->setPointingOffset(pointingOffset);
	    vbRow2CFBMap_p(irow) = cfb_l;
	  }

	/*
	//
	// Now do the in-row mappings.
	// 
	// Convert the VB polarizations to MuelllerElements.  
	for (Int ichan=0;ichan<nChan;ichan++)
	  {
	    //	    Double freq = vb.freq_p(ichan), wVal=vbs.uvw_p(irow,2);
	    Double freq = vbs.frequency()(ichan), wVal=abs(vbs.uvw()(irow)(2));
	    wVal *= freq/C::c;
	    for (Int ipol=0;ipol<nPol;ipol++)
	      {
		Vector<Int> vbPol = vbs.corrType();
		if (dataPol2ImPolMap(ipol) >= 0)
		  {
		    // The translate global functions comes from
		    // PolOuterProduct.{cc,h}.
		    //
		    // The code below translates, e.g.,
		    // Stokes::RR-->PolCrossProduct::RR-->MuellerElement.
		    MuellerElementType muellerElement;// = outerProduct.getMuellerElement(translateStokesToCrossPol(vbPol(ipol)));
		    Bool found=false;
		    Double f,w;
		    f=cfb_l->nearestFreq(found,freq);
		    w=cfb_l->nearestWVal(found,wVal);
		    if (!found) log_l << "Nearest freq. or w value not found " 
				      << freq << " " << wVal << LogIO::EXCEPTION;

		    vbRow2CFMap_p(ipol,ichan,irow) = cfb_l->getCFCellPtr(f, w, muellerElement);

		    // Bool Dummy;
		    // if (irow == 1)
		    //   {
		    // 	cerr << "#### " << ipol << ", " << ichan << ", " << irow << " " 
		    // 	     << vbRow2CFMap_p(ipol,ichan,irow)->getStorage()->getStorage(Dummy) << endl;
		    //   }
		  }
	      }
	  }
	*/
      }
    return statusCode;
  }
}
}
