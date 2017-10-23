################################################
#                                                                              
#           Regression/Benchmarking Script for wideband imaging with mosaicing               #
#                                                                             
################################################
#                                        
# Data : Simulated wideband, two pointing mosaic, with PB rotation and squint.
# Sky : Single flat-spectrum point source
#
# (1) CS with wideband A-projection
# (2) MS-MFS with nterms=2 and wideband A-projection
#
#                                                                            
################################################

import time
import os

# Data : reg_mawproject_apr13.ms   in  regression/wideband/
pathname=os.environ.get('CASAPATH').split()[0]
tempdir="./"

#tempdir="widebandmosaic_regression_data/"
#reuserepos=True
#if (reuserepos==False):
#   os.system("rm -rf "+tempdir)
#   os.mkdir(tempdir)
#   pathname=os.environ.get('CASAPATH').split()[0] + "/data/regression/wideband/wbmos_awproject.ms"
#   os.system("cp -r " + pathname + " " + tempdir)

# Initialize status flag
regstate = True;

# Start timers
startTime=time.time()
startProc=time.clock()

# Mark time
copyTime=time.time()

imname1=tempdir+'reg_widebandmosaic_cs.wbT.atT.mtT.psF'
imname2= tempdir+'reg_widebandmosaic_mtmfs.wbT.atT.mtT.psF'

msname = tempdir+'reg_mawproject_apr13.ms'


def locclean( vis='', imagename='', field='', spw='', phasecenter='', nterms=1, reffreq='1.5GHz', gridmode='advancedaprojection', wbawp=True, aterm=True, mterm=True, psterm=False, conjbeams=True, painc=360, rotpainc=360.0, cfcache='', imsize=512, cell='10.0arcsec', niter=10, gain=0.5, minpb=1e-04, usescratch=True ):
   
   im.open(msname,usescratch=usescratch);
   im.selectvis(spw=spw,field=field);

   models=[];
   restoreds=[];
   residuals=[];
   if(nterms>1):
      for tt in range(0,nterms):
	   models.append(imagename+'.model.tt'+str(tt));
	   restoreds.append(imagename+'.image.tt'+str(tt));
	   residuals.append(imagename+'.residual.tt'+str(tt));
	   im.defineimage(nx=imsize,ny=imsize,
			  cellx=cell,celly=cell,
                          phasecenter = [phasecenter],
			  nchan=1,stokes='I',mode='mfs')
#                          start='1.0GHz',step='2.0GHz');
	   im.make(models[tt]);
      algo = 'msmfs'
   else:
      models = [imagename+'.model'];
      restoreds = [imagename+'.image'];
      residuals = [imagename+'.residual'];
      im.defineimage(nx=imsize,ny=imsize,
			  cellx=cell,celly=cell,
			  nchan=1, 
                          phasecenter = [phasecenter],
                          stokes='IV',mode='mfs')
      im.make(models[0]);
      algo = 'hogbom'
   im.weight(type='natural');
   im.settaylorterms(ntaylorterms=nterms,
		     reffreq=(qa.convert(qa.unit(reffreq),"Hz"))['value']);
   im.setscales(scalemethod='uservector',uservector=[0]);
   if(gridmode=='advancedaprojection'):
     im.setoptions(ftmachine='awproject',
		   applypointingoffsets=False, 
		   dopbgriddingcorrections=True, 
		   cfcachedirname=cfcache, 
		   pastep=painc, rotpastep=rotpainc,pblimit=minpb, 
                   psterm=psterm, aterm=aterm, wbawp=wbawp,mterm=mterm,conjbeams=conjbeams);
   else:
     im.setoptions(ftmachine="ft");
#     im.setvp(dovp=True,usedefaultvp=True,telescope='EVLA');
#     im.makeimage(type='pb',image=imagename+'.pb');

   print('clean without mask');
   im.clean(model=models,image=restoreds,residual=residuals,algorithm=algo,threshold='0.0Jy',niter=niter,interactive=False,npercycle=niter,gain=gain,mask=''); 

   im.done();

###############################


if(regstate):
   npix=512
   phasecenter = 0 #'J2000 19h58m28.5s +40d44m01.50s'
   # Test (1) : CS clean with wideband A-Projection with freq-conjugate beams.
   print('-- CS Clean with Wideband AProjection and Mosaicing --')
   #default('clean')
   locclean( vis=msname, imagename=imname1, field='0,1', spw='', phasecenter=phasecenter, nterms=1, reffreq='1.5GHz', gridmode='advancedaprojection', wbawp=True, aterm=True, mterm=True, psterm=False, conjbeams=True, painc=360, rotpainc=360.0, cfcache=imname1+'.cfcache.dir', imsize=npix, cell='10.0arcsec', niter=10, gain=0.5, minpb=1e-04, usescratch=True)

   # Test (2) : MTMFS clean with wideband A-Projection with freq-conjugate beams.
   print('-- MTMFS Clean with Wideband AProjection and Mosaicing --')
   #default('clean')
   locclean( vis=msname, imagename=imname2, field='0,1', spw='', phasecenter=phasecenter, nterms=2, reffreq='1.5GHz', gridmode='advancedaprojection', wbawp=True, aterm=True, mterm=True, psterm=False, conjbeams=True, painc=360, rotpainc=360.0, cfcache=imname2+'.cfcache.dir', imsize=npix, cell='10.0arcsec', niter=10, gain=0.5, minpb=1e-04, usescratch=True)

   # Test (3) Mosaic with wbawp=False, and post-deconvolution PB-correction
   #print '-- MTMFS Clean with Reference-Frequency AProjection and Mosaicing + Post Deconv PB-Correction --'

# Stop timers
endProc=time.clock()
endTime=time.time()

# Start printing info.
import datetime
datestring=datetime.datetime.isoformat(datetime.datetime.today())
outfile='reg_widebandmosaic_'+datestring+'.log'
logfile=open(outfile,'w')

# Data summary
print('**************** Regression-test for wide-band mosaicing *************************', file=logfile)
# Perform the checks
print('', file=logfile)

if(not regstate):
   print('* Data file reg_mawproject_apr13.ms cannot be found', file=logfile);
else:

#####################################################
#####################################################
   # Truth for v1.0 : 20Nov2012, r22088 : 512x512 image, aliasing errors exist, no PS term.
   # Test 1
   # 17Dec2012 (SB): 
   #   Changeing correct_cs_intensity from 0.721199 to 0.80432087183.  code rev.  22329.
   #   AWProjectWBFT::pbFunc() returns 1.0 since division by avgPB in AWProjectFT::init2Vis()
   #   is not necessary (the normalization is done via the vectorized version of initi2Vis())
   #correct_cs_intensity = 0.804321
   #correct_cs_avgpb = 0.803946
   # Test 2
   #correct_mtmfs_intensity = 0.908331
   #correct_mtmfs_alpha = 0.00341234
   #correct_mtmfs_coeffpb_0 = 0.910124
   #correct_mtmfs_coeffpb_1 = 0.321995

   # 11 Apr 2011 (UR). Changed the simulated MS, with both pointings off-source.
   # These are pixel values at the center of the source, pixel 256, 315
   # Test 1
   #correct_cs_intensity = 0.692243
   #correct_cs_avgpb = 0.705166
   # Test 2
   #correct_mtmfs_intensity = 0.908868
   #correct_mtmfs_alpha = -0.0037804
   #correct_mtmfs_coeffpb_0 = 0.999975
   #correct_mtmfs_coeffpb_1 = 0.705602

   # 10 Sept 2013 (UR). 
   ## Removed CF rotation (changed from 5.0 to 360.0), and using PB=sqrt(sum_pbsq)
   ## Removed test for coeffPB_1, since this is no longer correct 
   # These are pixel values at the center of the source, pixel 256, 315
   # Test 1
   #correct_cs_intensity = 0.67087
   #correct_cs_avgpb = 0.695073
   # Test 2
   #correct_mtmfs_intensity = 0.679821
   #correct_mtmfs_alpha = 0.036
   #correct_mtmfs_avgPB.tt0 = 0.695073  # avgPB.tt0
##   correct_mtmfs_coeffpb_1 = 0.705602   ## Make this avgPB.tt1

   # 24 Oct 2013 (UR). -- This is with an oversampling of 50.
   ## Switched to using PBSQ. Search for 'PBSQWeight' in the code. (0.69*0.69 = 0.476)
   # These are pixel values at the center of the source, pixel 256, 315
   # Test 1
   #correct_cs_intensity = 0.66764
   #correct_cs_avgpb = 0.48338
   # Test 2
   #correct_mtmfs_intensity = 0.68578
   #correct_mtmfs_alpha = 0.002019
   #correct_mtmfs_avgPB.tt0 = 0.48338  # avgPB.tt0

   # 25 Oct 2013 (UR). -- This is with DEFAULT oversampling (20)
   # These are pixel values at the center of the source, pixel 256, 315
   # Test 1
   #correct_cs_intensity = 0.671525
   #correct_cs_avgpb = 0.4836
   # Test 2
   #correct_mtmfs_intensity = 0.6818
   #correct_mtmfs_alpha = 0.0188
   #correct_mtmfs_avgPB.tt0 = 0.4836

   # 28 Oct 2013 (UR). After fixing more PBSQ and sqrt errors.
   # These are pixel values at the center of the source, pixel 256, 315

   # 29 Jan 2014 (SB): Added Stokes-V testing 
   # Test 1
   correct_cs_intensity = 0.6710
   #correct_cs_intensity_v = -0.0007073
   correct_cs_intensity_v = -0.002090328
   # 03June15 (SB): Changed the v-intensity after switch to full-pol. CF comp. (PJ's) code
   # This can be un-done in VLACalcICF.cc by changing BeamCalc->calculateAperture(...,inStokes)
   # and setting tempmueller=-1
   correct_cs_intensity_v = -0.001594
   correct_cs_avgpb = 0.4830
   # Test 2
   correct_mtmfs_intensity = 0.6803
   correct_mtmfs_alpha = 0.0360
   # 03June15 (SB): Changed the alpha after switch to full-pol. CF comp. (PJ's) code
   # This can be un-done in VLACalcICF.cc by changing BeamCalc->calculateAperture(...,inStokes)
   # and setting tempmueller=-1
   correct_mtmfs_alpha = 0.0324
   correct_mtmfs_avgPB_tt0 = 0.4830


###################################################
###################################################
   print('*********************** Comparison of results **********************************', file=logfile)
   print('**                                                                            **', file=logfile)
   print('**      (1) CS with WideBand A-Projection         **', file=logfile)
   print('**                                                                            **', file=logfile)
   print('********************************************************************************', file=logfile)

   # Test 1 : Intensity
   if(os.path.exists(imname1+'.image')):
      ia.open(imname1+'.image');
      midpix = ia.pixelvalue([npix/2,npix/2])
      midpix = ia.pixelvalue([256,315,0,0])
      midpix_v = ia.pixelvalue([256,315,1,0])
      ia.close();
      diff_cs_intensity = abs( midpix['value']['value'] - correct_cs_intensity )/ abs(correct_cs_intensity);
      if(diff_cs_intensity<0.02): 
         print('* Passed Test 1 : peak cs_intensity_I test ', file=logfile);
      else: 
         print('* FAILED Test 1 : peak cs_intensity_I test at the 2-percent level ', file=logfile)
	 regstate = False;

      diff_cs_intensity_v = abs( midpix_v['value']['value'] - correct_cs_intensity_v )/ abs(correct_cs_intensity_v);
      if(diff_cs_intensity_v<0.02): 
         print('* Passed Test 1 : peak cs_intensity_V test ', file=logfile);
      else: 
         print('* FAILED Test 1 : peak cs_intensity_V test at the 2-percent level ', file=logfile)
	 regstate = False;

      print('-- Test 1 : peak cs_intensity_I : ' + str(midpix['value']['value']) + ' (' + str(correct_cs_intensity) + ')', file=logfile);
      print('-- Test 1 : peak cs_intensity_V : ' + str(midpix_v['value']['value']) + ' (' + str(correct_cs_intensity_v) + ')', file=logfile);
   else:
      print('-- FAILED Test 1 : No cs_intensity map generated', file=logfile);
      regstate = False;

   print('\n', file=logfile)

   # Test 1 : Primary Beam from inside the cfcache
   if(os.path.exists(imname1+'.cfcache.dir/avgPB')):
      ia.open(imname1+'.cfcache.dir/avgPB');
      midpix = ia.pixelvalue([npix/2,npix/2])
      midpix = ia.pixelvalue([256,315])
      ia.close();
      diff_cs_avgpb = abs( midpix['value']['value'] - correct_cs_avgpb )/ abs(correct_cs_avgpb);
      if(diff_cs_avgpb<0.02): 
         print('* Passed Test 1 : peak cs_avgpb test ', file=logfile);
      else: 
         print('* FAILED Test 1 : peak cs_avgpb test at the 2-percent level ', file=logfile)
	 regstate = False;
      print('-- Test 1 : peak cs_avgpb : ' + str(midpix['value']['value']) + ' (' + str(correct_cs_avgpb) + ')', file=logfile);
   else:
      print('-- FAILED Test 1 : No cs_avgpb map generated', file=logfile);
      regstate = False;

   print('\n', file=logfile)

###################################################
   print('*********************** Comparison of results **********************************', file=logfile)
   print('**                                                                            **', file=logfile)
   print('**      (2) MTMFS with WideBand A-Projection         **', file=logfile)
   print('**                                                                            **', file=logfile)
   print('********************************************************************************', file=logfile)

   # Test 2 : Intensity
   if(os.path.exists(imname2+'.image.tt0')):
      ia.open(imname2+'.image.tt0');
      midpix = ia.pixelvalue([npix/2,npix/2])
      midpix = ia.pixelvalue([256,315])
      ia.close();
      diff_mtmfs_intensity = abs( midpix['value']['value'] - correct_mtmfs_intensity )/ abs(correct_mtmfs_intensity);
      if(diff_mtmfs_intensity<0.02): 
         print('* Passed Test 2 : peak mtmfs_intensity test ', file=logfile);
      else: 
         print('* FAILED Test 2 : peak mtmfs_intensity test at the 2-percent level ', file=logfile)
	 regstate = False;
      print('-- Test 2 : peak mtmfs_intensity : ' + str(midpix['value']['value']) + ' (' + str(correct_mtmfs_intensity) + ')', file=logfile);
   else:
      print('-- FAILED Test 2 : No mtmfs_intensity map generated', file=logfile);
      regstate = False;

   print('\n', file=logfile)

   # Test 2 : Spectral Index
   if(os.path.exists(imname2+'.image.alpha')):
      ia.open(imname2+'.image.alpha');
      midpix = ia.pixelvalue([npix/2,npix/2])
      midpix = ia.pixelvalue([256,315])
      ia.close();
      diff_mtmfs_alpha = abs( midpix['value']['value'] - correct_mtmfs_alpha )/ abs(correct_mtmfs_alpha);
      if(diff_mtmfs_alpha<0.02): 
         print('* Passed Test 2 : peak mtmfs_alpha test ', file=logfile);
      else: 
         print('* FAILED Test 2 : peak mtmfs_alpha test at the 2-percent level ', file=logfile)
	 regstate = False;
      print('-- Test 2 : peak mtmfs_alpha : ' + str(midpix['value']['value']) + ' (' + str(correct_mtmfs_alpha) + ')', file=logfile);
   else:
      print('-- FAILED Test 2 : No mtmfs_alpha map generated', file=logfile);
      regstate = False;

   print('\n', file=logfile)

   # Test 2 : PB 0
   if(os.path.exists(imname2+'.cfcache.dir/avgPB.tt0')):
      ia.open(imname2+'.cfcache.dir/avgPB.tt0');
      midpix = ia.pixelvalue([npix/2,npix/2])
      midpix = ia.pixelvalue([256,315])
      ia.close();
      diff_mtmfs_avgPB_tt0 = abs( midpix['value']['value'] - correct_mtmfs_avgPB_tt0 )/ abs(correct_mtmfs_avgPB_tt0);
      if(diff_mtmfs_avgPB_tt0<0.02): 
         print('* Passed Test 2 : peak mtmfs_avgPB.tt0 test ', file=logfile);
      else: 
         print('* FAILED Test 2 : peak mtmfs_avgPB.tt0 test at the 2-percent level ', file=logfile)
	 regstate = False;
      print('-- Test 2 : peak mtmfs_avgPB.tt0 : ' + str(midpix['value']['value']) + ' (' + str(correct_mtmfs_avgPB_tt0) + ')', file=logfile);
   else:
      print('-- FAILED Test 2 : No mtmfs_avgPB.tt0 map generated', file=logfile);
      regstate = False;

   print('\n', file=logfile)

#   # Test 2 : PB 1
#   if(os.path.exists(imname2+'.cfcache.dir/coeffPB_1')):
#      ia.open(imname2+'.cfcache.dir/coeffPB_1');
#      midpix = ia.pixelvalue([npix/2,npix/2])
#      midpix = ia.pixelvalue([256,315])
#      ia.close();
#      diff_mtmfs_coeffpb_1 = abs( midpix['value']['value'] - correct_mtmfs_coeffpb_1 )/ abs(correct_mtmfs_coeffpb_1);
#      if(diff_mtmfs_coeffpb_1<0.02): 
#         print >>logfile,'* Passed Test 2 : peak mtmfs_coeffpb_1 test ';
#      else: 
#         print >>logfile,'* FAILED Test 2 : peak mtmfs_coeffpb_1 test at the 2-percent level '
#	 regstate = False;
#      print >>logfile,'-- Test 2 : peak mtmfs_coeffpb_1 : ' + str(midpix['value']['value']) + ' (' + str(correct_mtmfs_coeffpb_1) + ')';
#   else:
#      print >>logfile,'-- FAILED Test 2 : No mtmfs_coeffpb_1 map generated';
#      regstate = False;
#
#   print >>logfile, '\n'


###################################################

# Final verdict
if(regstate):
   print('PASSED regression test for wideband-mosaic-imaging.', file=logfile)
   print('')
   print('Regression PASSED')
   print('')
else:
   print('FAILED regression test for wideband-mosaic-imaging.', file=logfile)
   print('')
   print('Regression FAILED')
   print('')

print('', file=logfile)

# Print timing info
print('********************************************************************************', file=logfile)
print('**                         Benchmarking                                       **', file=logfile)
print('********************************************************************************', file=logfile)
print('Total wall clock time was: '+str(endTime - startTime), file=logfile)
print('Total CPU        time was: '+str(endProc - startProc), file=logfile)
print('Processing rate MB/s  was: '+str(278./(endTime - startTime)), file=logfile)
print('* Breakdown:                                                                   *', file=logfile)
print('*   copy         time was: '+str(copyTime-startTime), file=logfile)
print('*   imaging      time was: '+str(endTime-copyTime), file=logfile)
print('*                                                                              *', file=logfile)
print('********************************************************************************', file=logfile)

logfile.close()

