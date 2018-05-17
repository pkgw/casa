// -*- C++ -*-
//
// File to include to have the inner loops for de-gridding
// (prediction) without a function call.
//
//--------------------------------------------------------------------------------
{
  Vector<Int> iloc(4,0),iCFPos(4,0);
  Bool Dummy;
  Complex wt;

  const Int * __restrict__ gridInc_p_ptr = gridInc_p.getStorage(Dummy);
  const Complex * __restrict__ gridStore = grid.getStorage(Dummy);
  const Int* scaledSupport_ptr    = support.getStorage(Dummy);
  const Float *scaledSampling_ptr = sampling.getStorage(Dummy);
  const Double *off_ptr           = off.getStorage(Dummy);
  const Int *loc_ptr              = loc.getStorage(Dummy);
  const Int* convOrigin_ptr       = convOrigin.getStorage(Dummy);

  Int * __restrict__       iGrdPos_ptr   = igrdpos.getStorage(Dummy);
  Int *iCFPos_ptr                 = iCFPos.getStorage(Dummy);
  Int *iLoc_ptr                   = iloc.getStorage(Dummy);
  Int *cfInc_p_ptr                = cfInc_p.getStorage(Dummy);
  Int phaseGradOrigin_l[2]; 

  phaseGradOrigin_l[0] = cached_phaseGrad_p.shape()(0)/2;
  phaseGradOrigin_l[1] = cached_phaseGrad_p.shape()(1)/2;

  for(Int iy=-scaledSupport_ptr[1]; iy <= scaledSupport_ptr[1]; iy++) 
    {
      iLoc_ptr[1]    = SynthesisUtils::nint(scaledSampling_ptr[1]*iy+off_ptr[1]);
      iCFPos_ptr[1]  = iLoc_ptr[1] + convOrigin_ptr[1];
      iGrdPos_ptr[1] = loc_ptr[1]+iy;

      for(Int ix=-scaledSupport_ptr[0]; ix <= scaledSupport_ptr[0]; ix++) 
	{
	  iLoc_ptr[0]    = SynthesisUtils::nint(scaledSampling_ptr[0]*ix+off_ptr[0]);
	  iCFPos_ptr[0]  = iLoc_ptr[0] + convOrigin_ptr[0];
	  iGrdPos_ptr[0] = loc_ptr[0]+ix;
	  {
	    wt=getFrom4DArray((const Complex * __restrict__ &)convFuncV,
			      iCFPos_ptr,cfInc_p_ptr);
	    if (dataWVal <= 0.0) wt = conj(wt);
	    norm(apol)+=(wt);
	    if (finitePointingOffset) 
	      {
		wt *= (cached_phaseGrad_p(iloc[0]+phaseGradOrigin_l[0],
					  iloc[1]+phaseGradOrigin_l[1]));
	      }
	    // nvalue+=wt*grid(grdpos);
	    // The following uses raw index on the 4D grid
	    // nvalue+=wt*getFrom4DArray(gridStore,iPosPtr,gridInc);
	    nvalue+=wt*getFrom4DArray(gridStore,iGrdPos_ptr,gridInc_p_ptr);
	  }
	}
    }
}
//--------------------------------------------------------------------------------
