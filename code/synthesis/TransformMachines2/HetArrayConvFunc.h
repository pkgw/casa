//# HetArrayConvFunc.h: Definition for HetArrayConvFunc
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  General Public
//# License for more details.
//#
//# You should have received a copy of the GNU  General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$
#ifndef SYNTHESIS_TRANSFORM2_HETARRAYCONVFUNC_H
#define SYNTHESIS_TRANSFORM2_HETARRAYCONVFUNC_H

#include <synthesis/TransformMachines2/SimplePBConvFunc.h>

namespace casacore{

  template<class T> class ImageInterface;
  template<class T> class Matrix;
}

namespace casa {

  // <summary>  A class to support FTMachines get their convolution casacore::Function </summary>
  
  // <use visibility=export>
  // <prerequisite>
  //   <li> <linkto class=VisBuffer>VisBuffer</linkto> module
  // </prerequisite>
  // <etymology>
  // "HetArray" for Heterogeneous casacore::Array => different dish sizes
  // "ConvFunc" for Convolution Functions 
  //  appropriate convfunctions for each pair of antenna generated and cached
  // </etymology>
  //
  // <synopsis> 
  // FTMachines like WProjection and MosaicFT need convolution functions to 
  // deal with directional dependent issues...
  // this class and related ones provide and cache  such functions for re-use 
  //</synopsis>
  //Forward declarations
 namespace vi{class VisBuffer2;}
 namespace refim{
  
  class HetArrayConvFunc : public SimplePBConvFunc

  {
  public:
    HetArrayConvFunc();
    HetArrayConvFunc(const PBMathInterface::PBClass 
		     typeToUse, const casacore::String vpTable="");
    //Constructor from record
    //if for prediction only no need to recover fluxscale
    HetArrayConvFunc(const casacore::RecordInterface& rec, casacore::Bool calcFluxscale);
    virtual ~HetArrayConvFunc();

    //Returns the convfunctions in the Arrays...the rowMap maps the vb.row 
    //to the  plane of the convfunc appropriate...chanMap and polMap similarly 

    virtual void findConvFunction(const casacore::ImageInterface<casacore::Complex>& iimage, 
				  const vi::VisBuffer2& vb,
				    const casacore::Int& convSampling,
				  const casacore::Vector<casacore::Double>& visFreq,
				    casacore::Array<casacore::Complex>& convFunc,
				    casacore::Array<casacore::Complex>& weightConvFunc,
				    casacore::Vector<casacore::Int>& convsize,
				    casacore::Vector<casacore::Int>& convSupport,
				    casacore::Vector<casacore::Int>& polMap, casacore::Vector<casacore::Int>& chanMap, casacore::Vector<casacore::Int>& rowMap,
				  const casacore::Bool getConjConvFuncs=false,
				  const casacore::MVDirection& extraShift=casacore::MVDirection(0.0), const casacore::Bool useExtraShift=casacore::False
 								);
    virtual void findConvFunction2(const casacore::ImageInterface<casacore::Complex>& iimage, 
				  const vi::VisBuffer2& vb,
				    const casacore::Int& convSampling,
				  const casacore::Vector<casacore::Double>& visFreq,
				    casacore::Array<casacore::Complex>& convFunc,
				    casacore::Array<casacore::Complex>& weightConvFunc,
				    casacore::Vector<casacore::Int>& convsize,
				    casacore::Vector<casacore::Int>& convSupport,
				    casacore::Vector<casacore::Int>& polMap, casacore::Vector<casacore::Int>& chanMap, casacore::Vector<casacore::Int>& rowMap,
				  const casacore::Bool getConjConvFuncs=false,
				  const casacore::MVDirection& extraShift=casacore::MVDirection(0.0), const casacore::Bool useExtraShift=casacore::False
 								);
    virtual casacore::ImageInterface<casacore::Float>&  getFluxScaleImage();
    // slice flux scale images 
    virtual void sliceFluxScale(const casacore::Int npol);
    //Serialization
   virtual casacore::Bool toRecord(casacore::RecordInterface& rec);
   virtual casacore::Bool fromRecord(casacore::String& err, const casacore::RecordInterface& rec, casacore::Bool calcFluxscale=false);
   virtual void reset();
   virtual casacore::String name() {return casacore::String("HetArrayConvFunc");}
    //----------------------------------------------

    private:
   void applyGradientToYLine(const casacore::Int iy, casacore::Complex*& convFunctions, 
			     casacore::Complex*& convWeights, const casacore::Double pixXdir, const casacore::Double pixYdir, 
			     casacore::Int convSize, const casacore::Int ndishpair, const casacore::Int nchan, const casacore::Int nPol);
   void fillConjConvFunc(const casacore::Vector<casacore::Double>& beamFreqs);
   void fillConjConvFunc(casacore::Array<casacore::Complex>& out, const casacore::Array<casacore::Complex>& in , const casacore::Double& beamFreq);
   casacore::Int conjSupport(const casacore::Vector<casacore::Double>& beamFreqs);
      casacore::Int factorial(casacore::Int n);
      // the return value are -1 or false for not in cache yet but pointing direction 
      //seems to be inside image
      // 1 if value is cached..we have stopped caching..so it should not return this value
      // 2 pointing is off image ...thus valid but not useful
      casacore::Int checkPBOfField(const vi::VisBuffer2& vb, casacore::Vector<casacore::Int>& rowMap, const casacore::MVDirection& extraShift=casacore::MVDirection(0.0), const casacore::Bool useExtraShift=casacore::False);
      void findAntennaSizes(const vi::VisBuffer2& vb);
      void supportAndNormalize(casacore::Int plane, casacore::Int convSampling);
      void supportAndNormalizeLatt(casacore::Int plane, casacore::Int convSampling, casacore::TempLattice<casacore::Complex>& convFuncLat,
				   casacore::TempLattice<casacore::Complex>& weightConvFuncLat);
      casacore::Bool supportAndNormalize(casacore::Array<casacore::Complex>& convFunc,
					 casacore::Array<casacore::Complex>& weightConvFunc, casacore::Vector<casacore::Int>& supp, const casacore::Double factor=1.0);
      void init(const PBMathInterface::PBClass typeToUse);
      void makerowmap(const vi::VisBuffer2& vb, casacore::Vector<casacore::Int>& rowMap);
      casacore::Float interpLanczos( const casacore::Double& x , const casacore::Double& y, const casacore::Double& nx, const casacore::Double& ny,   const casacore::Float* data, const casacore::Float a=3);
      void diffPointingCorrection(const vi::VisBuffer2& vb, const casacore::Int origSupportSize, casacore::Vector<casacore::Int>& rowmap, const casacore::Int convSamp, const casacore::Vector<casacore::Double>& freqs, const casacore::Bool doConj=false);
      void rephaseBeams(casacore::Array<casacore::Complex>& pb,
			casacore::Array<casacore::Complex>& pb2,
			const vi::VisBuffer2& vb, const casacore::Int row,
			const casacore::Int origSize, const casacore::Double pixshift_x,
			const casacore::Double pixshift_y);
      //apply a phase gradient corresponding to (xshift, yshift) in pixels in the image domain
      void applyPhaseGradient(casacore::Array<casacore::Complex>& arr, const casacore::Double& xshift, const casacore::Double& yshift, const casacore::Int xsize, const casacore::Int ysize, const casacore::Int convsamp=1);
      casacore::Float sinc(const casacore::Float x) ;
      casacore::Array<casacore::Complex> resample(const casacore::Array<casacore::Complex>& inarray, const casacore::Double factor);
      casacore::Matrix<casacore::Complex> resample2(const casacore::Matrix<casacore::Complex>& inarray, const casacore::Double factor);
      void multiplySelfConjugate(casacore::ImageInterface<casacore::Complex>& im, const casacore::Double freq);
      void initSincCache();
      void roundShiftInPointing(casacore::Double& xshift, casacore::Double& yshift, const casacore::Int origSupport, const casacore::MDirection& dir1, const casacore::MDirection& dir2, const casacore::DirectionCoordinate dc);
      PBMathInterface::PBClass pbClass_p;
      //casacore::SimpleOrderedMap <casacore::String, casacore::Int> convFunctionMap_p;
      casacore::Vector<casacore::Int64> convFunctionMap_p;
      casacore::Int64 nDefined_p;
      casacore::SimpleOrderedMap <casacore::String, casacore::Int> antDiam2IndexMap_p;
      casacore::Vector<casacore::Int> antIndexToDiamIndex_p;
      ///The key to this map is a   tuple of 5 elements, 
      /// first 3 elements  are (ANT1_VP_ID,ANT2_VP_ID, MS_SPW_ID)
      /// last 2 elements  (x_shift, y_shift) in pixels of ANT2 pointing w.r.t ANT1 pointing
      std::map<std::tuple<casacore::Int, casacore::Int, casacore::Int, casacore::Int, casacore::Int>, casacore::Int> baselinePointingDiffIndex_p;
      casacore::Block<casacore::CountedPtr<PBMathInterface> > antMath_p;
      casacore::Int msId_p;
      casacore::Int actualConvIndex_p;
      casacore::Array<casacore::Complex> convFunc_p;
      casacore::Array<casacore::Complex> weightConvFunc_p;
      casacore::Array<casacore::Complex> convSave_p;
      casacore::Array<casacore::Complex> weightSave_p;
      casacore::Int convSize_p; 
      casacore::String vpTable_p;
      casacore::Vector<casacore::Int> convSupport_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > vpFunctions_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > pbFunctions_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > baslPBFunctions_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > baslWeightFunctions_p;
      casacore::Block<casacore::CountedPtr<casacore::Vector<casacore::Int> > > baslSupport_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > convFunctions_p;      
      casacore::Block < casacore::CountedPtr<casacore::Array<casacore::Complex> > > convFunctionsConjFreq_p;
      casacore::Block <casacore::CountedPtr<casacore::Array<casacore::Complex> > > convWeights_p;
      casacore::Block<casacore::CountedPtr<casacore::Vector<casacore::Int> > > convSizes_p;
      casacore::Block <casacore::CountedPtr<casacore::Vector<casacore::Int> > > convSupportBlock_p;
      casacore::Vector<casacore::Int> origConvSize_p;
      std::array<casacore::Float,8000> sincCache_p;
      casacore::Float *sincCachePtr_p;
      double timer1_p, timer2_p, timer3_p;
    };
}; //end of namespace refim
} // end namespace casa
#endif
