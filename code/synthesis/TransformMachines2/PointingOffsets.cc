// -*- C++ -*-
//# PointingOffsets.cc: Implementation of the PointingOffsets class
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

#include <msvis/MSVis/VisBufferUtil.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <measures/Measures/MeasTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <synthesis/TransformMachines2/PointingOffsets.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogOrigin.h>

using namespace casacore;
namespace casa{
  using namespace vi;
  using namespace refim;
  //
  //----------------------------------------------------------------------
  //
  PointingOffsets& PointingOffsets::operator=(const PointingOffsets& other)
  {
    if(this!=&other) 
      {
	imageDC_p = other.imageDC_p;
	imageObsInfo_p = other.imageObsInfo_p;
	doPointing_p = other.doPointing_p;
      }
    return *this;
  }
  //
  //----------------------------------------------------------------------
  //
  casacore::Vector<double> PointingOffsets::findMosaicPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
								     const VisBuffer2& vb, const casacore::Bool& doPointing)
  {
    VisBufferUtil vbUtils;
    //cerr << "#######: Using Mosaic Field pointing offsets" <<endl;
    storeImageParams(image,vb);
    //where in the image in pixels is this pointing
    pixFieldGrad_p.resize(2);
    int antId,numRow;
    antId = 0;
    numRow = 0;
    MDirection dir = vbUtils.getPointingDir(vb,antId,numRow,doPointing);
    //cerr << "VBDir: " << dir<< endl;
    thePix_p = toPix(vb, dir, dir);

    pixFieldGrad_p = gradPerPixel(thePix_p);
    //    MDirection fieldDir=direction1_p;
    //shift from center
    //pixFieldGrad_p(0) = pixFieldGrad_p(0) - double(nx_p / 2);
    //pixFieldGrad_p(1) = pixFieldGrad_p(1) - double(ny_p / 2);

    //Int convSampling=getOversampling(*psTerm_p,*wTerm_p,*aTerm_p);

    //phase gradient per pixel to apply
    //pixFieldGrad_p(0) = -pixFieldGrad_p(0)*2.0*C::pi/double(nx_p)/double(convOversampling_p);
    //pixFieldGrad_p(1) = -pixFieldGrad_p(1)*2.0*C::pi/double(ny_p)/double(convOversampling_p);

    return pixFieldGrad_p;
  };
  //
  //----------------------------------------------------------------------
  //
  casacore::Vector<double> PointingOffsets::findAntennaPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
								      const vi::VisBuffer2& vb, const casacore::Bool& doPointing)
  {
    casacore::Vector<double> antOffsets;
    storeImageParams(image,vb);
    VisBufferUtil vbUtils;
    int antId,numRow;
    antId = 0;
    numRow = 0;
    
    MDirection antDir =vbUtils.getPointingDir(vb, antId, numRow, doPointing); 
    //cerr << "AntDir="<<antDir << endl;
    thePix_p = toPix(vb, antDir, antDir);
    antOffsets = gradPerPixel(thePix_p);

    // // if (epJ_p.isNull())
    {
      //cerr << "#######: Using POINTING subtable to get antenna pointing offsets" << endl;
    
      //	VisBufferUtil vbUtils;

      //int vbrow=0;
      //int nant = vb.subtableColumns().antenna().nrow();
      //antOffsets.resize(nant);
      //for (int antid=0;antid<nant; antid++)
     	  {
	    //    MVDirection antdir=vbUtils.getPointingDir(vb, antid, vbrow).getValue();
     	    //    MVDirection vbdir=vb.direction1()(0).getValue();
     	    //antOffsets[antid]=antdir.separation(vbdir);
     	  }
      }
      
      return antOffsets;
  }
  //
  //----------------------------------------------------------------------
  //
  Vector<double> PointingOffsets::gradPerPixel(const Vector<double>& p)
  {
    Vector<double> gPP(2);
    gPP(0) = p(0) - double(nx_p/2);
    gPP(1) = p(1) - double(ny_p/2);

    gPP(0) = -gPP(0)*2.0*C::pi/double(nx_p)/double(convOversampling_p);
    gPP(1) = -gPP(1)*2.0*C::pi/double(ny_p)/double(convOversampling_p);

    return gPP;
  }
  //
  //----------------------------------------------------------------------
  //
  Vector<double>& PointingOffsets::toPix(const VisBuffer2& vb, 
					 const MDirection& dir1, 
					 const MDirection& dir2) 
  {
    thePix_p.resize(2);

     //    if(dc_p.directionType() !=  MDirection::castType(vb.direction1()(0).getRef().getType())){
    if(dc_p.directionType() !=  MDirection::castType(dir1.getRef().getType()))
      {
      //pointToPix_p.setModel(theDir);
      
      MEpoch timenow(Quantity(vb.time()(0), timeUnit_p), timeMType_p);
      //cout << "Ref " << vb.direction1()(0).getRefString() << " ep "
      //<< timenow.getRefString() << " time " <<
      //MVTime(timenow.getValue().getTime()).string(MVTime::YMD) <<
      //endl;
      pointFrame_p.resetEpoch(timenow);
      //////////////////////////
      //pointToPix holds pointFrame_p by reference...
      //thus good to go for conversion
      direction1_p=pointToPix_p(dir1);
      direction2_p=pointToPix_p(dir2);
      dc_p.toPixel(thePix_p, direction1_p);

     }
    else
      {
      direction1_p=dir1;
      direction2_p=dir2;
      dc_p.toPixel(thePix_p, dir1);
    }
    // Return the internal variable, just to make code more readable
    // at the place of the call.
    return thePix_p;
  };

  //
  //----------------------------------------------------------------------
  //
Vector<Double> PointingOffsets::findPointingOffset(const ImageInterface<Complex>& image,
						   const VisBuffer2& vb, const Bool doPointing)
  {

    if (!doPointing) 
      { 
	
	return findMosaicPointingOffset(image,vb,doPointing);
			
      }
    else 
      {
	return findAntennaPointingOffset(image,vb,doPointing);
	
      }
  }

  void PointingOffsets::storeImageParams(const casacore::ImageInterface<casacore::Complex>& iimage,
					 const VisBuffer2& vb) 
  {
    //image signature changed...rather simplistic for now
    if((iimage.shape().product() != nx_p*ny_p*nchan_p*npol_p) || nchan_p < 1)
      {
	csys_p=iimage.coordinates();
	int coordIndex=csys_p.findCoordinate(Coordinate::DIRECTION);
	AlwaysAssert(coordIndex>=0, AipsError);
	directionIndex_p=coordIndex;
	dc_p=csys_p.directionCoordinate(directionIndex_p);
	ObsInfo imInfo=csys_p.obsInfo();
	String tel= imInfo.telescope();
	MPosition pos;
	ROMSColumns mscol(vb.ms());
	if (vb.subtableColumns().observation().nrow() > 0) 
	  {
	    tel = vb.subtableColumns().observation().telescopeName()
	      (mscol.observationId()(0));
	  }
	if (tel.length() == 0 || !tel.contains("VLA") ||  
	    !MeasTable::Observatory(pos,tel)) 
	  {
	    // unknown observatory, use first antenna
	    int ant1 = vb.antenna1()(0);
	    pos=vb.subtableColumns().antenna().positionMeas()(ant1);
	  }
	//cout << "TELESCOPE " << tel << endl;
	//Store this to build epochs via the time access of visbuffer later

	timeMType_p=MEpoch::castType(mscol.timeMeas()(0).getRef().getType());
	timeUnit_p=Unit(mscol.timeMeas().measDesc().getUnits()(0).getName());
	// timeUnit_p=Unit("s");
	//cout << "UNIT " << timeUnit_p.getValue() << " name " << timeUnit_p.getName()  << endl;
	pointFrame_p=MeasFrame(imInfo.obsDate(), pos);
	MDirection::Ref elRef(dc_p.directionType(), pointFrame_p);
	//For now we set the conversion from this direction 
	pointToPix_p=MDirection::Convert( MDirection(), elRef);
	nx_p=iimage.shape()(coordIndex);
	ny_p=iimage.shape()(coordIndex+1);
	coordIndex=csys_p.findCoordinate(Coordinate::SPECTRAL);
	int pixAxis=csys_p.pixelAxes(coordIndex)[0];
	nchan_p=iimage.shape()(pixAxis);
	coordIndex=csys_p.findCoordinate(Coordinate::STOKES);
	pixAxis=csys_p.pixelAxes(coordIndex)[0];
	npol_p=iimage.shape()(pixAxis);
      }
  }
  };
