//# SIImageStoreMultiTerm.cc: Implementation of Imager.h
//# Copyright (C) 1997-2008
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
//# $Id$

#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/sstream.h>

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>

#include <casa/Logging.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>

#include <casa/OS/DirectoryIterator.h>
#include <casa/OS/File.h>
#include <casa/OS/Path.h>

#include <casa/OS/HostInfo.h>
#include <images/Images/TempImage.h>
#include <images/Images/PagedImage.h>
#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <images/Images/TempImage.h>
#include <images/Images/SubImage.h>
#include <images/Regions/ImageRegion.h>

#include <synthesis/ImagerObjects/SIImageStoreMultiTerm.h>

#include <casa/Arrays/MatrixMath.h>
#include <scimath/Mathematics/MatrixMathLA.h>

#include <sys/types.h>
#include <unistd.h>
using namespace std;

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  SIImageStoreMultiTerm::SIImageStoreMultiTerm():SIImageStore()
  {
    itsPsfs.resize(0);
    itsModels.resize(0);
    itsResiduals.resize(0);
    itsWeights.resize(0);
    itsImages.resize(0);
    itsSumWts.resize(0);
    itsPBs.resize(0);
    itsImagePBcors.resize(0);
    itsMask.reset( );
    itsGridWt.reset( );
    
    itsForwardGrids.resize(0);
    itsBackwardGrids.resize(0);

    itsNTerms=0;

    itsNFacets=1;
    itsFacetId=0;
    itsNChanChunks = 1;
    itsChanId = 0;
    itsNPolChunks = 1;
    itsPolId = 0;

    itsUseWeight=false;

    itsImageShape=IPosition(4,0,0,0,0);
    itsImageName=String("");
    itsCoordSys=CoordinateSystem();
    itsMiscInfo=Record();

    //    itsValidity = false;

    init();

    validate();

  }
  
  SIImageStoreMultiTerm::SIImageStoreMultiTerm(String imagename, 
					       CoordinateSystem &imcoordsys, 
					       IPosition imshape, 
					       const Int /*nfacets*/,
					       const Bool /*overwrite*/, 
					       uInt ntaylorterms,
					       Bool useweightimage)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","Open new Images",WHERE) );

    itsNTerms = ntaylorterms;

    itsPsfs.resize(2 * itsNTerms - 1);
    itsModels.resize(itsNTerms);
    itsResiduals.resize(itsNTerms);
    itsWeights.resize(2 * itsNTerms - 1);
    itsImages.resize(itsNTerms);
    itsSumWts.resize(2 * itsNTerms - 1);
    itsPBs.resize(itsNTerms);
    itsImagePBcors.resize(itsNTerms);

    itsMask.reset( );
    itsGridWt.reset( );

    itsForwardGrids.resize(itsNTerms);
    itsBackwardGrids.resize(2 * itsNTerms - 1);

    //cout << "Input imshape : " << imshape << endl;

    itsImageName = imagename;
    itsImageShape = imshape;
    itsCoordSys = imcoordsys;

    //    itsNFacets = nfacets;  // So that sumwt shape happens properly, via checkValidity
    //    itsFacetId = -1;

    itsNFacets=1;
    itsFacetId=0;
    itsNChanChunks = 1;
    itsChanId = 0;
    itsNPolChunks = 1;
    itsPolId = 0;


    itsUseWeight = useweightimage;

    itsMiscInfo=Record();

    init();

    validate();

  }

  SIImageStoreMultiTerm::SIImageStoreMultiTerm(String imagename, uInt ntaylorterms,const Bool ignorefacets) 
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","Open existing Images",WHERE) );

    itsNTerms = ntaylorterms;

    itsPsfs.resize(2 * itsNTerms - 1);
    itsModels.resize(itsNTerms);
    itsResiduals.resize(itsNTerms);
    itsWeights.resize(2 * itsNTerms - 1);
    itsImages.resize(itsNTerms);
    itsPBs.resize(itsNTerms);
    itsSumWts.resize(2 * itsNTerms - 1);
    itsMask.reset( );
    itsGridWt.reset( );
    itsImagePBcors.resize(itsNTerms);

    itsMiscInfo=Record();

    itsForwardGrids.resize(itsNTerms);
    itsBackwardGrids.resize(2 * itsNTerms - 1);

    itsImageName = imagename;

    itsNFacets=1;
    itsFacetId=0;
    itsNChanChunks = 1;
    itsChanId = 0;
    itsNPolChunks = 1;
    itsPolId = 0;

    Bool exists=true;
    Bool sumwtexists=true;
    for(uInt tix=0;tix<2*itsNTerms-1;tix++) 
      {
	if( tix<itsNTerms ) {
	    exists &= ( doesImageExist( itsImageName+String(".residual.tt")+String::toString(tix) ) ||
			doesImageExist( itsImageName+String(".psf.tt")+String::toString(tix) ) );
	  }else {
	    exists &= ( doesImageExist( itsImageName+String(".psf.tt")+String::toString(tix) ) );
	    sumwtexists &= ( doesImageExist( itsImageName+String(".sumwt.tt")+String::toString(tix) ) );
	  }
      }

    // The PSF or Residual images must exist. ( or the gridwt image)
    //  All this is just for the shape and coordinate system.
    if( exists || doesImageExist(itsImageName+String(".gridwt")) )
      {
	SHARED_PTR<ImageInterface<Float> > imptr;
	if( doesImageExist(itsImageName+String(".psf.tt0")) )
	  imptr.reset( new PagedImage<Float> (itsImageName+String(".psf.tt0")) );
	else if( doesImageExist(itsImageName+String(".residual.tt0")) )
	  imptr.reset( new PagedImage<Float> (itsImageName+String(".residual.tt0")) );
	else
	  imptr.reset( new PagedImage<Float> (itsImageName+String(".gridwt")) );
	  
	itsImageShape = imptr->shape();
	itsCoordSys = imptr->coordinates();
      }
    else
      {
	throw( AipsError( "Multi-term PSF or Residual Images do not exist. Please create one of them." ) );
      }

    if( doesImageExist(itsImageName+String(".residual.tt0")) || 
	doesImageExist(itsImageName+String(".psf.tt0")) )
      {
    if( sumwtexists )
      {
	SHARED_PTR<ImageInterface<Float> > imptr;
	imptr.reset( new PagedImage<Float> (itsImageName+String(".sumwt.tt0")) );
	itsNFacets = imptr->shape()[0];
	itsFacetId = 0;
	itsUseWeight = getUseWeightImage( *imptr );
	if( itsUseWeight && ! doesImageExist(itsImageName+String(".weight.tt0")) )
	  {
	    throw(AipsError("Internal error : MultiTerm Sumwt has a useweightimage=true but the weight image does not exist."));
	  }
      }
    else
      {
	throw( AipsError( "Multi-term SumWt does not exist. Please create PSFs or Residuals." ) );
      }
      }// if psf0 or res0 exist

    if( ignorefacets==true ) itsNFacets=1;

    init();
    validate();

  }

  /*
  /////////////Constructor with pointers already created else where but taken over here
  SIImageStoreMultiTerm::SIImageStoreMultiTerm(Block<SHARED_PTR<ImageInterface<Float> > > modelims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >residims,
					       Block<SHARED_PTR<ImageInterface<Float> > >psfims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >weightims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >restoredims,
					       Block<SHARED_PTR<ImageInterface<Float> > >sumwtims, 
					       SHARED_PTR<ImageInterface<Float> > newmask,
					       SHARED_PTR<ImageInterface<Float> > newalpha,
					       SHARED_PTR<ImageInterface<Float> > newbeta)
  {
    
    itsPsfs=psfims;
    itsModels=modelims;
    itsResiduals=residims;
    itsWeights=weightims;
    itsImages=restoredims;
    itsSumWts=sumwtims;
    itsMask = newmask;
    itsAlpha = newalpha;
    itsBeta = newbeta;

    itsNTerms = itsResiduals.nelements();
    itsMiscInfo=Record();

    AlwaysAssert( itsPsfs.nelements() == 2*itsNTerms-1 , AipsError ); 
    AlwaysAssert( itsPsfs.nelements()>0 && itsPsfs[0] , AipsError );
    AlwaysAssert( itsSumWts.nelements()>0 && itsSumWts[0] , AipsError );

    itsForwardGrids.resize( itsNTerms );
    itsBackwardGrids.resize( 2 * itsNTerms - 1 );

    itsImageShape=psfims[0]->shape();
    itsCoordSys = psfims[0]->coordinates();
    itsMiscInfo = psfims[0]->miscInfo();

    itsNFacets=sumwtims[0]->shape()[0];
    itsUseWeight=getUseWeightImage( *(sumwtims[0]) );

    itsImageName = String("reference");  // This is what the access functions use to guard against allocs...

    init();
    validate();
	
  }
  */

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////Constructor with pointers already created else where but taken over here
  SIImageStoreMultiTerm::SIImageStoreMultiTerm(Block<SHARED_PTR<ImageInterface<Float> > > modelims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >residims,
					       Block<SHARED_PTR<ImageInterface<Float> > >psfims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >weightims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >restoredims,
					       Block<SHARED_PTR<ImageInterface<Float> > >sumwtims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >pbims, 
					       Block<SHARED_PTR<ImageInterface<Float> > >restoredpbcorims, 
					       SHARED_PTR<ImageInterface<Float> > newmask,
					       SHARED_PTR<ImageInterface<Float> > newalpha,
					       SHARED_PTR<ImageInterface<Float> > newbeta,
					       SHARED_PTR<ImageInterface<Float> > newalphaerror,
					       SHARED_PTR<ImageInterface<Float> > newalphapbcor,
					       SHARED_PTR<ImageInterface<Float> > newbetapbcor,
					       CoordinateSystem& csys,
					       IPosition imshape,
					       String imagename,
					       const Int facet, const Int nfacets,
					       const Int chan, const Int nchanchunks,
					       const Int pol, const Int npolchunks)
  {
    
    itsPsfs=psfims;
    itsModels=modelims;
    itsResiduals=residims;
    itsWeights=weightims;
    itsImages=restoredims;
    itsSumWts=sumwtims;
    itsMask = newmask;
    itsAlpha = newalpha;
    itsBeta = newbeta;
    itsAlphaError = newalphaerror;
    itsAlphaPBcor = newalphapbcor;
    itsBetaPBcor = newbetapbcor;

    itsPBs=pbims;
    itsImagePBcors=restoredpbcorims;

    itsNTerms = itsResiduals.nelements();
    itsMiscInfo=Record();

    AlwaysAssert( itsPsfs.nelements() == 2*itsNTerms-1 , AipsError ); 
    //    AlwaysAssert( itsPsfs.nelements()>0 && itsPsfs[0] , AipsError );
    //    AlwaysAssert( itsSumWts.nelements()>0 && itsSumWts[0] , AipsError );

    itsForwardGrids.resize( itsNTerms );
    itsBackwardGrids.resize( 2 * itsNTerms - 1 );

    itsNFacets = nfacets;
    itsFacetId = facet;
    itsNChanChunks = nchanchunks;
    itsChanId = chan;
    itsNPolChunks = npolchunks;
    itsPolId = pol;

    itsParentImageShape = imshape; 
    itsImageShape = imshape;
    /////    itsImageShape = IPosition(4,0,0,0,0);

    itsCoordSys = csys; // Hopefully this doesn't change for a reference image
    itsImageName = imagename;

	
    //-----------------------
    init(); // Connect parent pointers to the images.
    //-----------------------

    // Set these to null, to be set later upon first access.
    // Setting to null will hopefully set all elements of each array, to NULL.
    itsPsfs=SHARED_PTR<ImageInterface<Float> >();  
    itsModels=SHARED_PTR<ImageInterface<Float> >();
    itsResiduals=SHARED_PTR<ImageInterface<Float> >();
    itsWeights=SHARED_PTR<ImageInterface<Float> >();
    itsImages=SHARED_PTR<ImageInterface<Float> >();
    itsSumWts=SHARED_PTR<ImageInterface<Float> >();
    itsPBs=SHARED_PTR<ImageInterface<Float> >();

    itsMask.reset( );

     validate();

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  uInt SIImageStoreMultiTerm::getNTaylorTerms(Bool dopsf)
  {
    return dopsf ? (2*itsNTerms-1) : itsNTerms;
   }

  // Check if images that are asked-for are ready and all have the same shape.
  /*
  Bool SIImageStoreMultiTerm::checkValidity(const Bool ipsf, const Bool iresidual, 
					    const Bool iweight, const Bool imodel, const Bool irestored, 
					    const Bool imask,const Bool isumwt,
					    const Bool ialpha, const Bool ibeta)
  {

    //    cout << "In MT::checkValidity imask is " << imask << endl;

    Bool valid = true;

    for(uInt tix=0; tix<2*itsNTerms-1; tix++)
      {
	
	if(ipsf==true)
	  { psf(tix); 
	    valid = valid & ( itsPsfs[tix] && itsPsfs[tix]->shape()==itsImageShape ); }
	if(iweight==true)
	  { weight(tix);  
	    valid = valid & ( itsWeights[tix] && itsWeights[tix]->shape()==itsImageShape ); }

	if(isumwt==true) {
	    IPosition useShape(itsImageShape);
	    useShape[0]=itsNFacets; useShape[1]=itsNFacets;
	    sumwt(tix);  
	    valid = valid & ( itsSumWts[tix] && itsSumWts[tix]->shape()==useShape ); 
	  }
	
	if( tix< itsNTerms )
	  {
	    if(iresidual==true)
	      { residual(tix);  
		valid = valid & ( itsResiduals[tix] && itsResiduals[tix]->shape()==itsImageShape ); }
	    if(imodel==true)
	      { model(tix);
		valid = valid & ( itsModels[tix] && itsModels[tix]->shape()==itsImageShape); }
	    if(irestored==true)
	      { image(tix);
		valid = valid & ( itsImages[tix] && itsImages[tix]->shape()==itsImageShape); }
	  }
      }
    
    if(imask==true)
      { mask(); valid = valid & ( itsMask && itsMask->shape()==itsImageShape); 
	//	cout << " Mask null ? " << (bool) itsMask << endl;
      }
    if(ialpha==true)
      { alpha();  valid = valid & ( itsAlpha && itsAlpha->shape()==itsImageShape ); }
    if(ibeta==true)
      { beta();  valid = valid & ( itsBeta && itsBeta->shape()==itsImageShape ); }

    return valid;
    
  }
  */

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  SIImageStoreMultiTerm::~SIImageStoreMultiTerm() 
  {
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  Bool SIImageStoreMultiTerm::releaseLocks() 
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","releaseLocks",WHERE) );

    for(uInt tix=0; tix<2*itsNTerms-1; tix++)
      {
	if( itsPsfs[tix] ) releaseImage( itsPsfs[tix] );
	if( itsWeights[tix] ) releaseImage( itsWeights[tix] );
	if( itsSumWts[tix] ) releaseImage( itsSumWts[tix] );
	if( tix < itsNTerms )
	  {
	    if( itsModels[tix] ) releaseImage( itsModels[tix] );
	    if( itsResiduals[tix] ) releaseImage( itsResiduals[tix] );
	    if( itsImages[tix] ) releaseImage( itsImages[tix] );
	    if( itsPBs[tix] ) releaseImage( itsPBs[tix] );
	    if( itsImagePBcors[tix] ) releaseImage( itsImagePBcors[tix] );
	  }
      }
    if( itsMask ) releaseImage( itsMask );
    if( itsAlpha ) releaseImage( itsAlpha );
    if( itsBeta ) releaseImage( itsBeta );
    if( itsAlphaError ) releaseImage( itsAlphaError );
    if( itsAlphaPBcor ) releaseImage( itsAlphaPBcor );
    if( itsBetaPBcor ) releaseImage( itsBetaPBcor );
    if( itsGridWt ) releaseImage( itsGridWt );
    
    return true; // do something more intelligent here.
  }

  Bool SIImageStoreMultiTerm::releaseComplexGrids() 
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","releaseComplexGrids",WHERE) );

    for(uInt tix=0; tix<2*itsNTerms-1; tix++)
      {
	if( itsBackwardGrids[tix] ) releaseImage( itsBackwardGrids[tix] );
	if( tix < itsNTerms )
	  {
	    if( itsForwardGrids[tix] ) releaseImage( itsForwardGrids[tix] );
	  }
      }
    return True; // do something more intelligent here.
  }


  Double SIImageStoreMultiTerm::getReferenceFrequency()
  {
    Double theRefFreq;

    Vector<Double> refpix = (itsCoordSys.spectralCoordinate()).referencePixel();
    AlwaysAssert( refpix.nelements()>0, AipsError );
    (itsCoordSys.spectralCoordinate()).toWorld( theRefFreq, refpix[0] );
    //    cout << "Reading ref freq as : " << theRefFreq << endl;
    return theRefFreq;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  Vector<String> SIImageStoreMultiTerm::getModelImageName()
  {
    Vector<String> mods(itsNTerms);
    for(uInt tix=0;tix<itsNTerms;tix++)
      {mods[tix]=itsImageName + imageExts(MODEL)+".tt"+String::toString(tix);}
    return mods;
  }

  void SIImageStoreMultiTerm::setModelImage( Vector<String> modelnames )
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","setModelImage",WHERE) );

    if( modelnames.nelements() > itsNTerms ) 
      { throw(AipsError("We currently cannot support more than nterms images as the starting model. "));
      }

    if( modelnames.nelements() > 0 && modelnames.nelements() <= itsNTerms )
      {
	for(uInt tix=0;tix<modelnames.nelements();tix++)
	  {
	    setModelImageOne( modelnames[tix], tix );
	  }
      }

  }

  /*
  void SIImageStoreMultiTerm::setModelImage( String modelname )
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","setModelImage",WHERE) );

    for(uInt tix=0;tix<itsNTerms;tix++)
      {
	
	Directory immodel( modelname+String(".model.tt")+String::toString(tix) );
	if( !immodel.exists() ) 
	  {
	    os << "Starting model image does not exist for term : " << tix << LogIO::POST;
	  }
	else
	  {
	    SHARED_PTR<PagedImage<Float> > newmodel( new PagedImage<Float>( modelname+String(".model.tt")+String::toString(tix) ) );
	    // Check shapes, coordsys with those of other images.  If different, try to re-grid here.
	    
	    if( newmodel->shape() != model(tix)->shape() )
	      {
		os << "Regridding input model to target coordinate system for term " << tix << LogIO::POST;
		regridToModelImage( *newmodel , tix);
		// For now, throw an exception.
		//throw( AipsError( "Input model image "+modelname+".model.tt"+String::toString(tix)+" is not the same shape as that defined for output in "+ itsImageName + ".model" ) );
	      }
	    else
	      {
		os << "Setting " << modelname << " as model for term " << tix << LogIO::POST;
		// Then, add its contents to itsModel.
		//itsModel->put( itsModel->get() + model->get() );
		/////( model(tix) )->put( newmodel->get() );
		model(tix)->copyData( LatticeExpr<Float> (*newmodel) );
	      }
	  }
      }//nterms
  }
  */ 

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::psf(uInt term)
  {
    AlwaysAssert( itsPsfs.nelements() > term, AipsError );
    accessImage( itsPsfs[term], itsParentPsfs[term], imageExts(PSF)+".tt"+String::toString(term) );
    return itsPsfs[term];
  }
  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::residual(uInt term)
  {
    accessImage( itsResiduals[term], itsParentResiduals[term], imageExts(RESIDUAL)+".tt"+String::toString(term) );
    return itsResiduals[term];
  }
  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::weight(uInt term)
  {
    accessImage( itsWeights[term], itsParentWeights[term], imageExts(WEIGHT)+".tt"+String::toString(term) );
    return itsWeights[term];
  }
  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::sumwt(uInt term)
  {
    accessImage( itsSumWts[term], itsParentSumWts[term], imageExts(SUMWT)+".tt"+String::toString(term) );

    
    if( itsNFacets>1 || itsNChanChunks>1 || itsNPolChunks>1 ) 
      {itsUseWeight = getUseWeightImage( *itsParentSumWts[0] );}
      setUseWeightImage( *(itsSumWts[term]) , itsUseWeight); // Sets a flag in the SumWt image. 

    return itsSumWts[term];
  }
  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::model(uInt term)
  {

    accessImage( itsModels[term], itsParentModels[term], imageExts(MODEL)+".tt"+String::toString(term) );

    // Set up header info the first time.
    ImageInfo info = itsModels[term]->imageInfo();
    String objectName("");
    if( itsMiscInfo.isDefined("OBJECT") ){ itsMiscInfo.get("OBJECT", objectName); }
    info.setObjectName(objectName);
    itsModels[term]->setImageInfo( info );
    itsModels[term]->setMiscInfo( itsMiscInfo );
    itsModels[term]->setUnits("Jy/pixel");

    return itsModels[term];
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::image(uInt term)
  {

    accessImage( itsImages[term], itsParentImages[term], imageExts(IMAGE)+".tt"+String::toString(term) );
    itsImages[term]->setUnits("Jy/beam");
    return itsImages[term];
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::pb(uInt term)
  {

    accessImage( itsPBs[term], itsParentPBs[term], imageExts(PB)+".tt"+String::toString(term) );
    return itsPBs[term];
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::imagepbcor(uInt term)
  {

    accessImage( itsImagePBcors[term], itsParentImagePBcors[term], imageExts(IMAGE)+".tt"+String::toString(term)+ ".pbcor" );
    itsImagePBcors[term]->setUnits("Jy/beam");
    return itsImagePBcors[term];
  }

    SHARED_PTR<ImageInterface<Complex> > SIImageStoreMultiTerm::forwardGrid(uInt term){
    if( itsForwardGrids[term] )// && (itsForwardGrids[term]->shape() == itsImageShape))
      return itsForwardGrids[term];
    Vector<Int> whichStokes(0);
    IPosition cimageShape;
    cimageShape=itsImageShape;
    CoordinateSystem cimageCoord = StokesImageUtil::CStokesCoord( itsCoordSys,
								  whichStokes, itsDataPolRep);
    cimageShape(2)=whichStokes.nelements();
    itsForwardGrids[term].reset(new TempImage<Complex>(TiledShape(cimageShape, tileShape()), cimageCoord, memoryBeforeLattice()));
    return itsForwardGrids[term];
  }
  SHARED_PTR<ImageInterface<Complex> > SIImageStoreMultiTerm::backwardGrid(uInt term){
  	  if( itsBackwardGrids[term] && (itsBackwardGrids[term]->shape() == itsImageShape))
  		  return itsBackwardGrids[term];
	  //	  cout << "MT : Making backward grid of shape : " << itsImageShape << endl;
    Vector<Int> whichStokes(0);
    IPosition cimageShape;
    cimageShape=itsImageShape;
    CoordinateSystem cimageCoord = StokesImageUtil::CStokesCoord( itsCoordSys,
								  whichStokes, itsDataPolRep);
    cimageShape(2)=whichStokes.nelements();
    itsBackwardGrids[term].reset(new TempImage<Complex>(TiledShape(cimageShape, tileShape()), cimageCoord, memoryBeforeLattice()));
    return itsBackwardGrids[term];
    }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::alpha()
  {
    if( itsAlpha && itsAlpha->shape() == itsImageShape ) { return itsAlpha; }
    //    checkRef( itsAlpha , "alpha" );
    itsAlpha = openImage( itsImageName+String(".alpha"), false );
    //    itsAlpha->setUnits("Alpha");
    return itsAlpha;
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::beta()
  {
    if( itsBeta && itsBeta->shape() == itsImageShape ) { return itsBeta; }
    //    checkRef( itsBeta , "beta" );
    itsBeta = openImage( itsImageName+String(".beta"), false );
    //    itsBeta->setUnits("Beta");
    return itsBeta;
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::alphaerror()
  {
    if( itsAlphaError && itsAlphaError->shape() == itsImageShape ) { return itsAlphaError; }
    //    checkRef( itsAlpha , "alpha" );
    itsAlphaError = openImage( itsImageName+String(".alpha.error"), false );
    //    itsAlpha->setUnits("Alpha");
    return itsAlphaError;
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::alphapbcor()
  {
    if( itsAlphaPBcor && itsAlphaPBcor->shape() == itsImageShape ) { return itsAlphaPBcor; }
    //    checkRef( itsAlphaPBcor , "alpha" );
    itsAlphaPBcor = openImage( itsImageName+String(".alpha.pbcor"), False );
    //    itsAlphaPBcor->setUnits("Alpha");
    return itsAlphaPBcor;
  }

  SHARED_PTR<ImageInterface<Float> > SIImageStoreMultiTerm::betapbcor()
  {
    if( itsBetaPBcor && itsBetaPBcor->shape() == itsImageShape ) { return itsBetaPBcor; }
    //    checkRef( itsBetaPBcor , "beta" );
    itsBetaPBcor = openImage( itsImageName+String(".beta.pbcor"), False );
    //    itsBetaPBcor->setUnits("Beta");
    return itsBetaPBcor;
  }



    // TODO : Move to an image-wrapper class ? Same function exists in SynthesisDeconvolver.
  Bool SIImageStoreMultiTerm::doesImageExist(String imagename)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","doesImageExist",WHERE) );
    Directory image( imagename );
    return image.exists();
  }


  void SIImageStoreMultiTerm::resetImages( Bool resetpsf, Bool resetresidual, Bool resetweight )
  {
    for(uInt tix=0;tix<2*itsNTerms-1;tix++)
      {
	if( resetpsf ) psf(tix)->set(0.0);
	if( resetweight && itsWeights[tix] ) weight(tix)->set(0.0);

	if( tix < itsNTerms ) {
	  if( resetresidual ) {
	    //removeMask( residual(tix) );
	    residual(tix)->set(0.0);
	  } 
	}
      }//nterms
  }
  
  void SIImageStoreMultiTerm::addImages( SHARED_PTR<SIImageStore> imagestoadd,
					 Bool addpsf, Bool addresidual, Bool addweight, Bool adddensity)
  {

    for(uInt tix=0;tix<2*itsNTerms-1;tix++)
      {
	
	if(addpsf)
	  {
	    LatticeExpr<Float> adderPsf( *(psf(tix)) + *(imagestoadd->psf(tix)) ); 
	    psf(tix)->copyData(adderPsf);	
	  }
	if(addweight)
	  {

	    if(getUseWeightImage( *(imagestoadd->sumwt(tix)) ) ) // Access and add weight only if it is needed.
	      {
		LatticeExpr<Float> adderWeight( *(weight(tix)) + *(imagestoadd->weight(tix)) ); 
		weight(tix)->copyData(adderWeight);
	      }

	    LatticeExpr<Float> adderSumWt( *(sumwt(tix)) + *(imagestoadd->sumwt(tix)) ); 
	    sumwt(tix)->copyData(adderSumWt);

	  }

	if(tix < itsNTerms && addresidual)
	  {
	    LatticeExpr<Float> adderRes( *(residual(tix)) + *(imagestoadd->residual(tix)) ); 
	    residual(tix)->copyData(adderRes);
	  }

	if( tix==0 && adddensity )
	  {
	    LatticeExpr<Float> adderDensity( *(gridwt()) + *(imagestoadd->gridwt()) ); 
	    gridwt()->copyData(adderDensity);
	  }

      }
  }

  void SIImageStoreMultiTerm::dividePSFByWeight(const Float /*pblimit*/)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","dividePSFByWeight",WHERE) );

    ////    for(uInt tix=0;tix<2*itsNTerms-1;tix++)
    for(Int tix=2*itsNTerms-1-1;tix>-1;tix--) // AAH go backwards so that zeroth term is normalized last..... sigh sigh sigh.
      {

	//cout << "npsfs : " << itsPsfs.nelements() << " and tix : " << tix << endl;

	normPSF(tix);
	
	if ( itsUseWeight) {
	  divideImageByWeightVal( *weight(tix) ); 
	}
	
      }     
   }

 void SIImageStoreMultiTerm::normalizePrimaryBeam(const Float pblimit)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","normalizePrimaryBeam",WHERE) );
    if ( itsUseWeight) {
      
      makePBFromWeight(pblimit);
    }
    else { makePBImage(pblimit); }
    //    calcSensitivity();
  }


  // Make another for the PSF too.
  void SIImageStoreMultiTerm::divideResidualByWeight(Float pblimit, const String normtype)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","divideResidualByWeight",WHERE) );

    if( itsUseWeight )  
    {
	itsPBScaleFactor = getPbMax();
      }

    for(uInt tix=0;tix<itsNTerms;tix++)
      {

	divideImageByWeightVal( *residual(tix) );

	//	if(doesImageExist(itsImageName+String(".weight.tt0"))  )
	if( itsUseWeight )
	{
	    
	    LatticeExpr<Float> deno;
	    if( normtype=="flatnoise"){
	      deno = LatticeExpr<Float> ( sqrt( abs(*(weight(0)) ) ) * itsPBScaleFactor );
	      os << LogIO::NORMAL1 << "Dividing " << itsImageName+String(".residual.tt")+String::toString(tix) ;
	      os << " by [ sqrt(weightimage) * " << itsPBScaleFactor ;
	      os << " ] to get flat noise with unit pb peak."<< LogIO::POST;
	      
	    }
	    if( normtype=="flatsky") {
	      deno = LatticeExpr<Float> ( *(weight(0)) );
	      os << LogIO::NORMAL1 << "Dividing " << itsImageName+String(".residual.tt")+String::toString(tix) ;
	      os << " by [ weight ] to get flat sky"<< LogIO::POST;
	    }
	    
	    Float scalepb = fabs(pblimit) * itsPBScaleFactor * itsPBScaleFactor ;
	    LatticeExpr<Float> mask( iif( (deno) > scalepb , 1.0, 0.0 ) );
	    LatticeExpr<Float> maskinv( iif( (deno) > scalepb , 0.0, 1.0 ) );
	    LatticeExpr<Float> ratio( ( (*(residual(tix))) * mask ) / ( deno + maskinv ) );

	    residual(tix)->copyData(ratio);
	  }

	if( (residual(tix)->getDefaultMask()=="") && hasPB()  && pblimit >= 0.0 )
	  {copyMask(pb(),residual(tix));}

	if( pblimit <0.0 && (residual(tix)->getDefaultMask()).matches("mask0") ) removeMask( residual(tix) );

      }
  }

  void SIImageStoreMultiTerm::divideModelByWeight(Float pblimit, const String normtype)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","divideModelByWeight",WHERE) );

        if( 	itsUseWeight // only when needed
	&& weight() )// i.e. only when possible. For an initial starting model, don't need wt anyway.
      {

	if( normtype=="flatsky") {
	  Array<Float> arrmod;
	  model(0)->get( arrmod, true );

	  os << "Model is already flat sky with peak flux : " << max(arrmod);
	  os << ". No need to divide before prediction" << LogIO::POST;
	  
	  return;
	  }

	  itsPBScaleFactor = getPbMax();

	for(uInt tix=0;tix<itsNTerms;tix++)
	  {
	    os << LogIO::NORMAL1 << "Dividing " << itsImageName+String(".model")+String::toString(tix);
	    os << " by [ sqrt(weight) / " << itsPBScaleFactor ;
	    os <<" ] to get to flat sky model before prediction" << LogIO::POST;
	    
	    LatticeExpr<Float> deno( sqrt( abs(*(weight(0))) ) / itsPBScaleFactor );
	    
	    LatticeExpr<Float> mask( iif( (deno) > fabs(pblimit) , 1.0, 0.0 ) );
	    LatticeExpr<Float> maskinv( iif( (deno) > fabs(pblimit) , 0.0, 1.0 ) );
	    LatticeExpr<Float> ratio( ( (*(model(tix))) * mask ) / ( deno + maskinv ) );

	    itsModels[tix]->copyData(ratio);
	  }    
      }
  }


  void SIImageStoreMultiTerm::multiplyModelByWeight(Float pblimit, const String normtype)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","multiplyModelByWeight",WHERE) );

    
    if(        itsUseWeight // only when needed
	&& weight() )// i.e. only when possible. For an initial starting model, don't need wt anyway.
      {

	if( normtype=="flatsky") {
	  os << "Model is already flat sky. No need to multiply back after prediction" << LogIO::POST;
	  return;
	  }

	  itsPBScaleFactor = getPbMax();

	for(uInt tix=0;tix<itsNTerms;tix++)
	  {

	    os << LogIO::NORMAL1 << "Multiplying " << itsImageName+String(".model")+String::toString(tix);
	  os << " by [ sqrt(weight) / " << itsPBScaleFactor;
	  os <<  " ] to take model back to flat noise with unit pb peak." << LogIO::POST;
	  
	  LatticeExpr<Float> deno( sqrt( abs(*(weight(0)) ) ) / itsPBScaleFactor );

	  LatticeExpr<Float> mask( iif( (deno) > fabs(pblimit) , 1.0, 0.0 ) );
	  LatticeExpr<Float> maskinv( iif( (deno) > fabs(pblimit) , 0.0, 1.0 ) );
	  LatticeExpr<Float> ratio( ( (*(model(tix))) * mask ) * ( deno + maskinv ) );

	    itsModels[tix]->copyData(ratio);
	  }    
      }
  }


  void SIImageStoreMultiTerm::restore(GaussianBeam& rbeam, String& usebeam, uInt /*term*/)
  {

    LogIO os( LogOrigin("SIImageStoreMultiTerm","restore",WHERE) );

    for(uInt tix=0; tix<itsNTerms; tix++)
      {
	SIImageStore::restore(rbeam, usebeam, tix);
      }	
   
    calculateAlphaBeta("image");
 
  }// restore

  void SIImageStoreMultiTerm::calculateAlphaBeta(String imtype)
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","calculateAlphaBeta",WHERE) );
    try
      {

	// Check if this is Stokes I.
	Bool isOnlyStokesI=True;

	Vector<Int> stokes (  (itsParentCoordSys.stokesCoordinate()).stokes() );
	AlwaysAssert( stokes.nelements()>0 , AipsError);
	if( stokes.nelements()>1 || stokes[0]!=1 ) isOnlyStokesI=False;

	if( ! isOnlyStokesI )
	  { os << "Alpha and Beta images will not be computed for images that contain anything other than stokes I in them." << LogIO::POST; }

	//Calculate alpha, beta
	    if( itsNTerms > 1 && isOnlyStokesI)
	      {
		// Calculate alpha and beta
		LatticeExprNode leMaxRes = max( *( residual(0) ) );
		Float maxres = leMaxRes.getFloat();
		Float specthreshold = maxres/10.0;  //////////// do something better here..... 
		
		os << "Calculating spectral parameters for Intensity > peakresidual/10 = " << specthreshold << " Jy/beam" << LogIO::POST;

		/////////////////////////////////////////////////////////
		/////// Calculate alpha
		if(imtype=="image")
		  {
		    LatticeExpr<Float> mask1(iif(((*(image(0))))>(specthreshold),1.0,0.0));
		    LatticeExpr<Float> mask0(iif(((*(image(0))))>(specthreshold),0.0,1.0));
		    LatticeExpr<Float> alphacalc( (((*(image(1))))*mask1)/(((*(image(0))))+(mask0)) );
		    alpha()->copyData(alphacalc);
		    
		    ImageInfo ii = image(0)->imageInfo();
		    // Set the restoring beam for alpha
		    alpha()->setImageInfo(ii);
		    //alpha()->table().unmarkForDelete();
		    
		    // Make a mask for the alpha image
		    LatticeExpr<Bool> lemask(iif(((*(image(0))) > specthreshold) , True, False));
		    removeMask( alpha() );
		    createMask( lemask, (alpha()) );

		    /////// Calculate alpha error
		    alphaerror()->set(0.0);

		    LatticeExpr<Float> alphacalcerror( abs(alphacalc) * sqrt( ( (*residual(0)*mask1)/(*image(0)+mask0) )*( (*residual(0)*mask1)/(*image(0)+mask0) ) + ( (*residual(1)*mask1)/(*image(1)+mask0) )*( (*residual(1)*mask1)/(*image(1)+mask0) )  ) );
		    alphaerror()->copyData(alphacalcerror);
		    alphaerror()->setImageInfo(ii);
		    removeMask(alphaerror());
		    createMask(lemask, alphaerror());
		    //		    alphaerror()->table().unmarkForDelete();      
		    os << "Written Spectral Index Error Image : " << alphaerror()->name() << LogIO::POST;

		if(itsNTerms>2) // calculate beta too.
		  {
			beta()->set(0.0);
			LatticeExpr<Float> betacalc( (*image(2)*mask1)/((*image(0))+(mask0))-0.5*(*alpha())*((*alpha())-1.0) );
			beta()->copyData(betacalc);
			beta()->setImageInfo(ii);
			//imbeta.setUnits(Unit("Spectral Curvature"));
			removeMask(beta());
			createMask(lemask, beta());
			os << "Written Spectral Curvature Image : " << beta()->name() << LogIO::POST;
		  }

		  }
		else
		  {
		    LatticeExpr<Float> mask1(iif(((*(imagepbcor(0))))>(specthreshold),1.0,0.0));
		    LatticeExpr<Float> mask0(iif(((*(imagepbcor(0))))>(specthreshold),0.0,1.0));
		    LatticeExpr<Float> alphacalc( (((*(imagepbcor(1))))*mask1)/(((*(imagepbcor(0))))+(mask0)) );
		    alphapbcor()->copyData(alphacalc);
		    
		    ImageInfo ii = image(0)->imageInfo();
		    // Set the restoring beam for alpha
		    alphapbcor()->setImageInfo(ii);
		    //alpha()->table().unmarkForDelete();
		    
		    // Make a mask for the alpha image
		    LatticeExpr<Bool> lemask(iif(((*(imagepbcor(0))) > specthreshold) , True, False));
		    removeMask( alphapbcor() );
		    createMask( lemask, (alphapbcor()) );

		/////////////////////////////////////////////////////////
		if(itsNTerms>2) // calculate beta too.
		  {
			betapbcor()->set(0.0);
			LatticeExpr<Float> betacalc( (*imagepbcor(2)*mask1)/((*imagepbcor(0))+(mask0))-0.5*(*alphapbcor())*((*alphapbcor())-1.0) );
			betapbcor()->copyData(betacalc);
			betapbcor()->setImageInfo(ii);
			//imbeta.setUnits(Unit("Spectral Curvature"));
			removeMask(betapbcor());
			createMask(lemask, betapbcor());
			os << "Written Spectral Curvature Image : " << betapbcor()->name() << LogIO::POST;
		  }
		
		  }// pbcor

	      }//if nterms>1
      }
    catch(AipsError &x)
      {
	throw( AipsError("Error in computing Alpha (and Beta) : " + x.getMesg() ) );
      }

  }// calculateAlphaBeta

  void SIImageStoreMultiTerm::pbcor()
  {

    LogIO os( LogOrigin("SIImageStoreMultiTerm","pbcorPlane",WHERE) );

    /// Temp Code to prevent this approximate PBCOR from happening for EVLA data
    if(1)
      {
	String telescope = itsCoordSys.obsInfo().telescope();
	if ( telescope != "ALMA" )
	  {
	    os << LogIO::WARN << "Wideband (multi-term) PB correction is not yet available via tclean in the 4.7 release. Please use the widebandpbcor task instead. "<< LogIO::POST;
	    return;
	  }
	else
	  {
	    os << LogIO::WARN << "Wideband (multi-term) PB Correction is currently only an approximation. It assumes no PB frequency dependence. This code has been added for the 4.7 release to support the current ALMA pipeline, which does not apply corrections for the frequency dependence of the primary beam across small fractional bandwidths. Please look at the help for the 'pbcor' parameter and use the widebandpbcor task if needed. " <<LogIO::POST;
	  }

	
      }


    // message saying that it's only stokes I for now...

    for(uInt tix=0; tix<itsNTerms; tix++)
      {
	SIImageStore::pbcor(tix);
      }	

    calculateAlphaBeta("pbcor");

  }

  void SIImageStoreMultiTerm::calcSensitivity()
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","calcSensitivity",WHERE) );

    // Construct Hessian.
    
    Matrix<Float> hess(IPosition(2,itsNTerms,itsNTerms));
    for(uInt tay1=0;tay1<itsNTerms;tay1++)
      for(uInt tay2=0;tay2<itsNTerms;tay2++)
	{
	  //uInt taymin = (tay1<=tay2)? tay1 : tay2;
	  //uInt taymax = (tay1>=tay2)? tay1 : tay2;
	  //uInt ind = (taymax*(taymax+1)/2)+taymin;

	  uInt ind = tay1+tay2;
	  AlwaysAssert( ind < 2*itsNTerms-1, AipsError );

	  Array<Float> lsumwt;
	  sumwt( ind )->get( lsumwt, false );
	  //	  cout << "lsumwt shape : " << lsumwt.shape() << endl;
	  AlwaysAssert( lsumwt.shape().nelements()==4, AipsError );
	  AlwaysAssert( lsumwt.shape()[0]>0, AipsError );

	  //	  hess(tay1,tay2) = lsumwt(IPosition(1,0)); //Always pick sumwt from first facet only.
	  hess(tay1,tay2) = lsumwt(IPosition(4,0,0,0,0)); //Always pick sumwt from first facet only.
	}

    os << "Multi-Term Hessian Matrix : " << hess << LogIO::POST;

    // Invert Hessian. 
    try
      {
	Float deter=0.0;
	Matrix<Float> invhess;
	invertSymPosDef(invhess,deter,hess);
	os << "Multi-Term Covariance Matrix : " << invhess << LogIO::POST;

	// Just print the sqrt of the diagonal elements. 
	
	for(uInt tix=0;tix<itsNTerms;tix++)
	  {
	    os << "[" << itsImageName << "][Taylor"<< tix << "] Theoretical sensitivity (Jy/bm):" ;
	    if( invhess(tix,tix) > 0.0 ) { os << sqrt(invhess(tix,tix)) << LogIO::POST; }
	    else { os << "none" << LogIO::POST; }
	  }
      }
    catch(AipsError &x)
      {
	os << LogIO::WARN << "Cannot invert Hessian Matrix : " << x.getMesg()  << " || Calculating approximate sensitivity " << LogIO::POST;
	
	// Approximate : 1/h^2
	for(uInt tix=0;tix<itsNTerms;tix++)
	  {
	    Array<Float> lsumwt;
	    AlwaysAssert( 2*tix < 2*itsNTerms-1, AipsError );
	    sumwt(2*tix)->get( lsumwt , false ); 
	    
	    IPosition shp( lsumwt.shape() );
	    //cout << "Sumwt shape : " << shp << " : " << lsumwt << endl;
	    //AlwaysAssert( shp.nelements()==4 , AipsError );
	    
	    os << "[" << itsImageName << "][Taylor"<< tix << "] Approx Theoretical sensitivity (Jy/bm):" ;
	    
	    IPosition it(4,0,0,0,0);
	    for ( it[0]=0; it[0]<shp[0]; it[0]++)
	      for ( it[1]=0; it[1]<shp[1]; it[1]++)
		for ( it[2]=0; it[2]<shp[2]; it[2]++)
		  for ( it[3]=0; it[3]<shp[3]; it[3]++)
		    {
		      if( shp[0]>1 ){os << "f"<< it[0]+(it[1]*shp[0]) << ":" ;}
		      if( lsumwt( it )  > 1e-07 ) 
			{ 
			  os << 1.0/(sqrt(lsumwt(it))) << " " ;
			}
		      else
			{
			  os << "none" << " ";
			}
		    }
	    
	    os << LogIO::POST;
	  }
      }

    calcFractionalBandwidth();
  }
 
  Double SIImageStoreMultiTerm::calcFractionalBandwidth()
  {

    LogIO os( LogOrigin("SIImageStoreMultiTerm","calcFractionalBandwidth",WHERE) );

    Double fbw=1.0;

    for(uInt i=0; i<itsCoordSys.nCoordinates(); i++)
    {
      if( itsCoordSys.type(i) == Coordinate::SPECTRAL )
	{
	  SpectralCoordinate speccoord(itsCoordSys.spectralCoordinate(i));
	  Double startfreq=0.0,startpixel=-0.5;
	  Double endfreq=0.0,endpixel=+0.5;
	  speccoord.toWorld(startfreq,startpixel);
	  speccoord.toWorld(endfreq,endpixel);
	  Double midfreq = (endfreq+startfreq)/2.0;
	  fbw = ((endfreq - startfreq)/midfreq) * 100.0;
	  os << "MFS frequency range : " << startfreq/1.0e+9 << " GHz -> " << endfreq/1.0e+9 << "GHz."; 
	  os << "Fractional Bandwidth : " << fbw << " %.";
	  os << "Reference Frequency for Taylor Expansion : "<< getReferenceFrequency()/1.0e+9 << "GHz." << LogIO::POST;
	}
    }
    return fbw;
  }

 
  // Check for non-zero model (this is different from getting model flux, for derived SIIMMT)
  Bool SIImageStoreMultiTerm::isModelEmpty()
  {
    /// There MUST be a more efficient way to do this !!!!!  I hope. 
    /// Maybe save this info and change it anytime model() is accessed.... 
    Bool emptymodel=true;
    for(uInt tix=0;tix<itsNTerms;tix++)
      {
	//if( fabs( getModelFlux(tix) ) > 1e-08  ) emptymodel=false;
	if( doesImageExist(itsImageName+String(".model.tt")+String::toString(tix)) && 
	    fabs( getModelFlux(tix) ) > 1e-08  ) emptymodel=false;
      } 
    return  emptymodel;
  }

  void SIImageStoreMultiTerm::printImageStats()
  {
    LogIO os( LogOrigin("SIImageStoreMultiTerm","printImageStats",WHERE) );
    Float minresmask, maxresmask, minres, maxres;
    //    findMinMax( residual()->get(), mask()->get(), minres, maxres, minresmask, maxresmask );

    if (hasMask())
      {
	findMinMaxLattice(*residual(), *mask() , maxres,maxresmask, minres, minresmask);
      }
    else
      {
	LatticeExprNode pres( max( *residual() ) );
	maxres = pres.getFloat();
	LatticeExprNode pres2( min( *residual() ) );
	minres = pres2.getFloat();
      }

    os << "[" << itsImageName << "]" ;
    os << " Peak residual (max,min) " ;
    if( minresmask!=0.0 || maxresmask!=0.0 )
      { os << "within mask : (" << maxresmask << "," << minresmask << ") "; }
    os << "over full image : (" << maxres << "," << minres << ")" << LogIO::POST;

    os << "[" << itsImageName << "] Total Model Flux : " ;
    for(uInt tix=0;tix<itsNTerms;tix++)
      {os << getModelFlux(tix) << "(tt" << String::toString(tix) << ")"; }
    os<<LogIO::POST;

  }
  
  SHARED_PTR<SIImageStore> SIImageStoreMultiTerm::getSubImageStore(const Int facet, const Int nfacets, 
							  const Int chan, const Int nchanchunks, 
							  const Int pol, const Int npolchunks)
  {
    return SHARED_PTR<SIImageStore>(new SIImageStoreMultiTerm(itsModels, itsResiduals, itsPsfs, itsWeights, itsImages, itsSumWts, itsPBs, itsImagePBcors, itsMask, itsAlpha, itsBeta, itsAlphaError,itsAlphaPBcor, itsBetaPBcor,  itsCoordSys,itsParentImageShape, itsImageName, facet, nfacets,chan,nchanchunks,pol,npolchunks));
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////

//
  //-------------------------------------------------------------
  // Initialize the internals of the object.  Perhaps other such
  // initializations in the constructors can be moved here too.
  //
  void SIImageStoreMultiTerm::init()
  {
    imageExts.resize(MAX_IMAGE_IDS);
    
    imageExts(MASK)=".mask";
    imageExts(PSF)=".psf";
    imageExts(MODEL)=".model";
    imageExts(RESIDUAL)=".residual";
    imageExts(WEIGHT)=".weight";
    imageExts(IMAGE)=".image";
    imageExts(SUMWT)=".sumwt";
    imageExts(GRIDWT)=".gridwt";
    imageExts(PB)=".pb";
    imageExts(FORWARDGRID)=".forward";
    imageExts(BACKWARDGRID)=".backward";
    //    imageExts(IMAGEPBCOR)=".image.pbcor";

    itsParentPsfs.resize(itsPsfs.nelements());
    itsParentModels.resize(itsModels.nelements());
    itsParentResiduals.resize(itsResiduals.nelements());
    itsParentWeights.resize(itsWeights.nelements());
    itsParentImages.resize(itsImages.nelements());
    itsParentSumWts.resize(itsSumWts.nelements());
    itsParentImagePBcors.resize(itsImagePBcors.nelements());
    itsParentPBs.resize(itsPBs.nelements());
   
    itsParentPsfs = itsPsfs;
    itsParentModels=itsModels;
    itsParentResiduals=itsResiduals;
    itsParentWeights=itsWeights;
    itsParentImages=itsImages;
    itsParentSumWts=itsSumWts;
    itsParentImagePBcors=itsImagePBcors;
    itsParentPBs=itsPBs;

    itsParentMask=itsMask;

    itsParentImageShape = itsImageShape;
    itsParentCoordSys = itsCoordSys;
    if( itsNFacets>1 || itsNChanChunks>1 || itsNPolChunks>1 ) { itsImageShape=IPosition(4,0,0,0,0); }

    itsOpened=0;

  }

  /*
  Bool SIImageStoreMultiTerm::getUseWeightImage()
  {  if( itsParentSumWts.nelements()==0 || ! itsParentSumWts[0] ) 
      {return false;} 
    else
      {
	Bool ret = SIImageStore::getUseWeightImage( *(itsParentSumWts[0]) );
	return ret;
      }
  }
  */
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

} //# NAMESPACE CASA - END

