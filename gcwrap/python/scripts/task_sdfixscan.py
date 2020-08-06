# sd task for image processing (fft_mask or model)
import os
import time
import numpy
import numpy.fft as npfft

from taskinit import casalog, gentools, utilstool, qatool

import sdutil

def create_4d_image(infile, outfile):
    (ia,) = gentools(['ia'])
    ia.open(infile)
    image_shape = ia.shape()
    try:
        if len(image_shape) < 4:
            # add degenerate axes

            cs = ia.coordsys()
            axistypes = cs.axiscoordinatetypes()
            no_stokes = 'Stokes' not in axistypes
            no_spectral = 'Spectral' not in axistypes
            stokes_axis = 'I' if no_stokes else ''
            outimage = ia.adddegaxes(outfile=outfile, spectral=no_spectral,
                                           stokes=stokes_axis)
        else:
            # generage complete copy of input image using subimage
            outimage = ia.subimage(outfile=outfile)
    finally:
        if len(image_shape) < 4: cs.done()
        ia.close()
        
    return outimage


@sdutil.sdtask_decorator
def sdfixscan(infiles, mode, numpoly, beamsize, smoothsize, direction, maskwidth, tmax, tmin, outfile, overwrite):
    with sdutil.sdtask_manager(sdfixscan_worker, locals()) as worker:
        worker.initialize()
        worker.execute()
        worker.finalize()
    
class sdfixscan_worker(sdutil.sdtask_interface):
    def __init__(self, **kwargs):
        super(sdfixscan_worker,self).__init__(**kwargs)

    def __del__(self, base=sdutil.sdtask_interface):
        # cleanup method must be called when the instance is
        # deleted
        self.cleanup()
        super(sdfixscan_worker,self).__del__()

    def initialize(self):
        self.parameter_check()
        
        # temporary filename
        tmpstr = time.ctime().replace( ' ', '_' ).replace( ':', '_' )
        self.tmpmskname = 'masked.'+tmpstr+'.im'
        self.tmpconvname = 'convolve2d.'+tmpstr+'.im'
        self.tmppolyname = 'polyfit.'+tmpstr+'.im'
        # set tempolary filename
        self.tmprealname = []
        self.tmpimagname = []
        self.image = None
        self.convimage = None
        self.polyimage = None
        self.imageorg = None
        self.realimage = None
        self.imagimage = None
        if type(self.infiles) == str:
            self.tmprealname.append( 'fft.'+tmpstr+'.real..0.im' )
            self.tmpimagname.append( 'fft.'+tmpstr+'.imag.0.im' )
        else:
            for i in range(len(self.infiles)):
                self.tmprealname.append( 'fft.%s.%s.real.im' % (tmpstr,i) )
                self.tmpimagname.append( 'fft.%s.%s.imag.im' % (tmpstr,i) )

        # default output filename
        if self.outfile == '':
            self.outfile = 'sdfixscan.out.im'
        casalog.post( 'outfile=%s' % self.outfile )

        # threshold
        self.nolimit = 'nolimit'
        if self.tmin == 0.0 and self.tmax == 0.0:
            self.thresh = []
        elif self.tmin > self.tmax:
            casalog.post('tmin > tmax. Swapped.' )
            self.thresh = [self.tmax, self.tmin]
        elif self.tmin == self.tmax:
            if self.tmin > 0.0:
                casalog.post( 'tmin == tmax. Use tmin as minumum threshold.' )
                self.thresh = [self.tmin, self.nolimit]
            else:
                casalog.post( 'tmin == tmax. Use tmax as maximum threshold.' )
                self.thresh = [self.nolimit, self.tmin]
        else:
            self.thresh = [self.tmin, self.tmax]

    def parameter_check(self):
        if self.mode.lower() == 'model':
            # Pressed-out method
            # check input file
            if type(self.infiles) == list:
                if len(self.infiles) != 1:
                    raise Exception("infiles allows only one input file for pressed-out method.") 
                else:
                    self.infiles = self.infiles[0]
            # check direction
            if type(self.direction) == list:
                if len(self.direction) != 1:
                    raise Exception("direction allows only one direction for pressed-out method.")
                else:
                    self.direction = self.direction[0]
        elif self.mode.lower() == 'fft_mask':
            # FFT-based basket-weaving method
            # check input file
            if type(self.infiles) == str or \
                   (type(self.infiles) == list and len(self.infiles) < 2):
                raise Exception("infiles should be a list of input images for Basket-Weaving.")

            # check direction
            if type(self.direction) == float:
                raise Exception('direction must have at least two different direction.')
            else:
                if len(self.direction) < 2:
                    raise Exception('direction must have at least two different direction.')
        else:
            raise Exception('Unsupported processing mode: %s'%(self.mode))

    def execute(self):
        if self.mode.lower() == 'model':
           self.__execute_press()
        elif self.mode.lower() == 'fft_mask':
            self.__execute_basket_weaving()

    def __execute_press(self):
        ###
        # Pressed-out method (Sofue & Reich 1979)
        ###
        casalog.post( 'Apply Pressed-out method' )

        # CAS-5410 Use private tools inside task scripts
        ia = gentools(['ia'])[0]

        # mask
        self.image = ia.newimagefromimage(infile=self.infiles,outfile=self.tmpmskname)
        # back-up original mask name
        is_initial_mask = (self.image.maskhandler('default')[0] != '')
        temp_maskname = "temporal"
        imshape = self.image.shape()
        ndim = len(imshape)
        nx = imshape[0]
        ny = imshape[1]
        if len(self.thresh) == 0:
            casalog.post( 'Use whole region' )
        else:
            # mask pixels beyond thresholds
            maskstr = ("mask('%s')" % self.tmpmskname)
            if self.thresh[0] != self.nolimit:
                maskstr += (" && '%s'>=%f" % (self.tmpmskname, self.thresh[0]))
            if self.thresh[1] != self.nolimit:
                maskstr += (" && '%s'<=%f" % (self.tmpmskname, self.thresh[1]))
            # Need to flush to image once to calcmask ... sigh
            self.image.done()
            self.image = ia.newimage( self.tmpmskname )
            self.image.calcmask(mask=maskstr, name=temp_maskname, asdefault=True)

        # smoothing
        #bmajor = 0.0
        #bminor = 0.0
        # CAS-5410 Use private tools inside task scripts
        qa = qatool()
        if type(self.beamsize) == str:
            qbeamsize = qa.quantity(self.beamsize)
        else:
            qbeamsize = qa.quantity(self.beamsize,'arcsec')
        if type(self.smoothsize) == str:
            #bmajor = smoothsize
            #bminor = smoothsize
            qsmoothsize = qa.quantity(self.smoothsize)
        else:
            #bmajor = '%sarcsec' % (beamsize*smoothsize)
            #bminor = '%sarcsec' % (beamsize*smoothsize)
            qsmoothsize = qa.mul(qbeamsize,self.smoothsize)
        bmajor = qsmoothsize
        bminor = qsmoothsize
        pa = qa.quantity(0.0, 'deg')
        # masked channels are replaced by zero and convolved here.
        self.convimage = self.image.convolve2d( outfile=self.tmppolyname, major=bmajor, minor=bminor, pa=pa, 
                                                overwrite=True )
        self.convimage.done()

        # get dTij (original - smoothed)
        self.convimage = ia.imagecalc(outfile=self.tmpconvname, 
                                      pixels='"{org}" - "{conv}"'.format(org=self.tmpmskname,
                                                                         conv=self.tmppolyname),
                                      overwrite=True)

        # polynomial fit
        fitaxis = 0
        if self.direction == 0.0:
            fitaxis = 0
        elif self.direction == 90.0:
            fitaxis = 1
        else:
            raise Exception("Sorry, the task don't support inclined scan with respect to horizontal or vertical axis, right now.")
        # Replace duplicated method ia.fitpolynomial with
        # ia.fitprofile 
        #polyimage = convimage.fitpolynomial( fitfile=tmppolyname, axis=fitaxis, order=numpoly, overwrite=True )
        #polyimage.done()
        if os.path.exists( self.tmppolyname ):
            # CAS-5410 Use private tools inside task scripts
            cu = utilstool()
            cu.removetable([self.tmppolyname])
        self.convimage.setbrightnessunit('K')
        # Unfortunately, ia.fitprofile is very fragile.
        # Using numpy instead for fitting with masked pixels (KS, 2014/07/02)
        #resultdic = self.convimage.fitprofile( model=self.tmppolyname, axis=fitaxis, poly=self.numpoly, ngauss=0, multifit=True, gmncomps=0 )
        self.__polynomial_fit_model(image=self.tmpmskname, model=self.tmppolyname,axis=fitaxis, order=self.numpoly)
        polyimage = ia.newimage( self.tmppolyname )
        # set back defalut mask (need to get from self.image)
        avail_mask = polyimage.maskhandler('get')
        if is_initial_mask:
            casalog.post("copying mask from %s"%(self.infiles))
            polyimage.calcmask("mask('%s')"%self.infiles,asdefault=True)
        else: #no mask in the original image
            polyimage.calcmask('T', asdefault=True)
        if temp_maskname in avail_mask:
            polyimage.maskhandler('delete', name=temp_maskname)

        # subtract fitted image from original map
        subtracted = ia.imagecalc(outfile=self.outfile, 
                                  pixels='"{org}" - "{fit}"'.format(org=self.infiles,
                                                                    fit=self.tmppolyname),
                                  overwrite=self.overwrite)
        subtracted.done()

        # finalization
        polyimage.done(remove=True)
        self.convimage.done(remove=True)
        self.image.done()

    def __polynomial_fit_model(self, image=None, model=None, axis=0, order=2):
        if not image or not os.path.exists(image):
            raise RuntimeError("No image found to fit.")
        if os.path.exists( model ):
            # CAS-5410 Use private tools inside task scripts
            cu = utilstool()
            cu.removetable([model])
        tmpia = gentools(['ia'])[0]
        modelimg = tmpia.newimagefromimage(infile=image,outfile=model)
        try:
            if tmpia.isopen(): tmpia.close()
            imshape = modelimg.shape()
            # the axis order of [ra, dec, chan(, pol)] is assumed throughout the task.
            ndim = len(imshape)
            nx = imshape[0]
            ny = imshape[1]
            # an xy-plane can be fit simultaneously (doing per plane to save memory)
            if ndim == 3:
                get_blc = lambda i, j: [0, 0, i]
                get_trc = lambda i, j: [nx-1, ny-1, i]
                imshape2 = imshape[2]
                imshape3 = 1
            elif ndim == 4:
                get_blc = lambda i, j: [0, 0, i, j]
                get_trc = lambda i, j: [nx-1, ny-1, i, j]
                imshape2 = imshape[2]
                imshape3 = imshape[3]
            else: # ndim == 2 
                get_blc = lambda i, j: [0,0]
                get_trc = lambda i, j: [nx-1, ny-1]
                imshape2 = 1
                imshape3 = 1
                
            for i3 in range(imshape3):
                for i2 in range(imshape2):
                    blc = get_blc(i2, i3)
                    trc = get_trc(i2, i3)
                    dslice = modelimg.getchunk(blc, trc)
                    mslice = modelimg.getchunk(blc, trc, getmask=True)
                    model = self._get_polyfit_model_array(dslice.reshape(nx,ny),
                                                          mslice.reshape(nx,ny),
                                                          axis, order)
                    modelimg.putchunk(model, blc)
                
            # the fit model image itself is free from invalid pixels
            modelimg.calcmask('T', asdefault=True)
        except: raise
        finally: modelimg.close()

    def _get_polyfit_model_array(self, data, mask, axis, order):
        if axis ==1:
            tmp  = data.transpose()
            data = tmp
            tmp  = mask.transpose()
            mask = tmp
            del tmp
        nx = data.shape[0]
        ny = data.shape[1]
        x = list(range(nx))
        flag = mask ^ True #invert mask for masked array
        mdata = numpy.ma.masked_array(data,flag)
        retc = numpy.ma.polyfit(x,mdata,order)
        del flag
        coeffs = retc.transpose()
        tmpmodel = numpy.zeros([nx, ny])
        for iy in range(ny):
            tmpmodel[:,iy] = numpy.poly1d(coeffs[iy])(x)
        if axis == 1: return tmpmodel.transpose()
        return tmpmodel

    def __execute_basket_weaving(self):
        ###
        # Basket-Weaving (Emerson & Grave 1988)
        ###
        casalog.post( 'Apply Basket-Weaving' )

        # CAS-5410 Use private tools inside task scripts
        ia = gentools(['ia'])[0]

        # initial setup
        outimage = ia.newimagefromimage( infile=self.infiles[0], outfile=self.outfile, overwrite=self.overwrite )
        imshape_out = outimage.shape()
        ndim_out = len(imshape_out)
        coordsys = outimage.coordsys()
        axis_types = coordsys.axiscoordinatetypes()
        # direction axis should always exist
        try:
            direction_axis0 = axis_types.index('Direction')
            direction_axis1 = axis_types[direction_axis0+1:].index('Direction') + direction_axis0 + 1
        except IndexError:
            raise RuntimeError('Direction axes don\'t exist.')
        finally:
            coordsys.done()
        nx = imshape_out[direction_axis0]
        ny = imshape_out[direction_axis1]
        tmp=[]
        nfile = len(self.infiles)
        for i in range(nfile):
            tmp.append(numpy.zeros(imshape_out,dtype=float))
        maskedpixel=numpy.array(tmp)
        del tmp

        # direction
        dirs = []
        if len(self.direction) == nfile:
            dirs = self.direction
        else:
            casalog.post( 'direction information is extrapolated.' )
            for i in range(nfile):
                dirs.append(self.direction[i%len(self.direction)])

        # maskwidth
        masks = []
        if isinstance(self.maskwidth, int) or isinstance(self.maskwidth, float):
            for i in range(nfile):
                masks.append(self.maskwidth)
        elif isinstance(self.maskwidth, list):#  and nfile != len(self.maskwidth):
            for i in range(nfile):
                masks.append(self.maskwidth[i%len(self.maskwidth)])
        for i in range(len(masks)):
            masks[i] = 0.01 * masks[i]
        
        # mask
        for i in range(nfile):
            self.realimage = create_4d_image(self.infiles[i], self.tmprealname[i])
            self.imagimage = self.realimage.subimage(outfile=self.tmpimagname[i])
            
            # replace masked pixels with 0.0
            if not self.realimage.getchunk(getmask=True).all():
                casalog.post("Replacing masked pixels with 0.0 in %d-th image" % (i))
                self.realimage.replacemaskedpixels(0.0)
            self.realimage.close()
            self.imagimage.close()
          
        # Below working images are all 4D regardless of dimension of input images  
        # image shape for temporary images (always 4D)
        ia.open(self.tmprealname[0])
        imshape = ia.shape()
        ndim = len(imshape)
        ia.close()

        if len(self.thresh) == 0:
            casalog.post( 'Use whole region' )
        else:
            for i in range(nfile):
                self.realimage = ia.newimage( self.tmprealname[i] )
                for iaxis2 in range(imshape[2]):
                    for iaxis3 in range(imshape[3]):
                        pixmsk = self.realimage.getchunk( [0,0,iaxis2,iaxis3], [nx-1,ny-1,iaxis2,iaxis3])
                        for ix in range(pixmsk.shape[0]):
                            for iy in range(pixmsk.shape[1]):
                                if self.thresh[0] == self.nolimit:
                                    if pixmsk[ix][iy] > self.thresh[1]:
                                        maskedpixel[i][ix][iy][iaxis2][iaxis3]=pixmsk[ix][iy]
                                        pixmsk[ix][iy] = 0.0
                                elif self.thresh[1] == self.nolimit:
                                    if pixmsk[ix][iy] < self.thresh[0]:
                                        maskedpixel[i][ix][iy][iaxis2][iaxis3]=pixmsk[ix][iy]
                                        pixmsk[ix][iy] = 0.0
                                else:
                                    if pixmsk[ix][iy] < self.thresh[0] or pixmsk[ix][iy] > self.thresh[1]:
                                        maskedpixel[i][ix][iy][iaxis2][iaxis3]=pixmsk[ix][iy]
                                        pixmsk[ix][iy] = 0.0
                        self.realimage.putchunk( pixmsk, [0,0,iaxis2,iaxis3] )
                self.realimage.close()
        maskedvalue=None
        if any(maskedpixel.flatten()!=0.0):
                maskedvalue=maskedpixel.mean(axis=0)
        del maskedpixel

        # set weight factor
        weights = numpy.ones(shape=(nfile,nx,ny), dtype=float)
        eps = 1.0e-5
        dtor = numpy.pi / 180.0
        for i in range(nfile):
            scan_direction = ''
            if abs(numpy.sin(dirs[i]*dtor)) < eps: # direction is around 0 deg
                maskw = 0.5 * nx * masks[i]
                scan_direction = 'horizontal'
            elif abs(numpy.cos(dirs[i]*dtor)) < eps: # direction is around 90 deg
                maskw = 0.5 * ny * masks[i]
                scan_direction = 'vertical'
            else:
                maskw = 0.5 * numpy.sqrt(nx*ny) * masks[i]
            for ix in range(nx):
                halfwx = (nx-1)/2
                for iy in range(ny):
                    halfwy = (ny-1)/2
                    if scan_direction == 'horizontal':
                        #dd = abs(float(ix) - 0.5*(nx-1))
                        dd = abs(float(ix) - halfwx) # for CAS-9434
                    elif scan_direction == 'vertical':
                        #dd = abs(float(iy) - 0.5*(ny-1))
                        dd = abs(float(iy) - halfwy) # for CAS-9434
                    else:
                        tand = numpy.tan((dirs[i]-90.0)*dtor)
                        #dd = abs((float(ix) - 0.5*(nx-1)) * tand - (float(iy) - 0.5*(ny-1)))
                        dd = abs((float(ix) - halfwx) * tand - (float(iy) - halfwy)) # for CAS-9434
                        dd = dd / numpy.sqrt(1.0 + tand*tand)
                    if dd < maskw:
                        cosd = numpy.cos(0.5*numpy.pi*dd/maskw)
                        weights[i][ix][iy] = 1.0 - cosd * cosd
                    if weights[i][ix][iy] == 0.0:
                        weights[i][ix][iy] += eps*0.01
            """
            if abs(numpy.sin(dirs[i]*dtor)) < eps:
                # direction is around 0 deg
                maskw = 0.5 * nx * masks[i] 
                for ix in range(nx):
                    for iy in range(ny):
                        dd = abs( float(ix) - 0.5 * (nx-1) )
                        if dd < maskw:
                            cosd = numpy.cos(0.5*numpy.pi*dd/maskw)
                            weights[i][ix][iy] = 1.0 - cosd * cosd
                        if weights[i][ix][iy] == 0.0:
                            weights[i][ix][iy] += eps*0.01
            elif abs(numpy.cos(dirs[i]*dtor)) < eps:
                # direction is around 90 deg
                maskw = 0.5 * ny * masks[i]
                for ix in range(nx):
                    for iy in range(ny):
                        dd = abs( float(iy) - 0.5 * (ny-1) )
                        if dd < maskw:
                            cosd = numpy.cos(0.5*numpy.pi*dd/maskw)
                            weights[i][ix][iy] = 1.0 - cosd * cosd
                        if weights[i][ix][iy] == 0.0:
                            weights[i][ix][iy] += eps*0.01
            else:
                maskw = 0.5 * numpy.sqrt( nx * ny ) * masks[i]
                for ix in range(nx):
                    for iy in range(ny):
                        tand = numpy.tan((dirs[i]-90.0)*dtor)
                        dd = abs( ix * tand - iy - 0.5 * (nx-1) * tand + 0.5 * (ny-1) )
                        dd = dd / numpy.sqrt( 1.0 + tand * tand )
                        if dd < maskw:
                            cosd = numpy.cos(0.5*numpy.pi*dd/maskw)
                            weights[i][ix][iy] = 1.0 - cosd * cosd
                        if weights[i][ix][iy] == 0.0:
                            weights[i][ix][iy] += eps*0.01 
            """
            # shift
            xshift = -((ny-1)/2)
            yshift = -((nx-1)/2)
            for ix in range(xshift,0,1):
                tmp = weights[i,:,0].copy()
                weights[i,:,0:ny-1] = weights[i,:,1:ny].copy()
                weights[i,:,ny-1] = tmp
            for iy in range(yshift,0,1):
                tmp = weights[i,0:1].copy()
                weights[i,0:nx-1] = weights[i,1:nx].copy()
                weights[i,nx-1:nx] = tmp

        # FFT
        for i in range(nfile):
            self.realimage = ia.newimage( self.tmprealname[i] )
            self.imagimage = ia.newimage( self.tmpimagname[i] )
            for iaxis2 in range(imshape[2]):
                for iaxis3 in range(imshape[3]):
                    pixval = self.realimage.getchunk( [0,0,iaxis2,iaxis3], [nx-1,ny-1,iaxis2,iaxis3] )
                    pixval = pixval.reshape((nx,ny))
                    pixfft = npfft.fft2( pixval )
                    pixfft = pixfft.reshape((nx,ny,1,1))
                    self.realimage.putchunk( pixfft.real, [0,0,iaxis2,iaxis3] )
                    self.imagimage.putchunk( pixfft.imag, [0,0,iaxis2,iaxis3] )
                    del pixval, pixfft
            self.realimage.close()
            self.imagimage.close()

        # weighted mean
        for ichan in range(imshape[2]):
            for iaxis3 in range(imshape[3]):
                pixout = numpy.zeros( shape=(nx,ny), dtype=complex )
                denom = numpy.zeros( shape=(nx,ny), dtype=float )
                for i in range(nfile):
                    self.realimage = ia.newimage( self.tmprealname[i] )
                    self.imagimage = ia.newimage( self.tmpimagname[i] )
                    pixval = self.realimage.getchunk( [0,0,ichan,iaxis3], [nx-1,ny-1,ichan,iaxis3] ) \
                        + self.imagimage.getchunk( [0,0,ichan,iaxis3], [nx-1,ny-1,ichan,iaxis3] ) * 1.0j
                    pixval = pixval.reshape((nx,ny))
                    pixout = pixout + pixval * weights[i]
                    denom = denom + weights[i]
                    self.realimage.close()
                    self.imagimage.close()
                pixout = pixout / denom
                pixout = pixout.reshape((nx,ny,1,1))
                self.realimage = ia.newimage( self.tmprealname[0] )
                self.imagimage = ia.newimage( self.tmpimagname[0] )
                self.realimage.putchunk( pixout.real, [0,0,ichan,iaxis3] )
                self.imagimage.putchunk( pixout.imag, [0,0,ichan,iaxis3] )
                self.realimage.close()
                self.imagimage.close()

        # inverse FFT
        self.realimage = ia.newimage( self.tmprealname[0] )
        self.imagimage = ia.newimage( self.tmpimagname[0] )
        for ichan in range(imshape[2]):
            for iaxis3 in range(imshape[3]):
                pixval = self.realimage.getchunk( [0,0,ichan,iaxis3], [nx-1,ny-1,ichan,iaxis3] ) \
                    + self.imagimage.getchunk( [0,0,ichan,iaxis3], [nx-1,ny-1,ichan,iaxis3] ) * 1.0j
                pixval = pixval.reshape((nx,ny))
                pixifft = npfft.ifft2( pixval )
                pixifft = pixifft.reshape((nx,ny,1,1))
                self.realimage.putchunk(pixifft.real, blc=[0,0,ichan,iaxis3])
                del pixval, pixifft
        if maskedvalue is not None:
            self.realimage.putchunk(self.realimage.getchunk()+maskedvalue)
            
        # put result into outimage
        chunk = self.realimage.getchunk()
        outimage.putchunk(chunk.reshape(imshape_out))
        # handling of output image mask
        maskstr = ""
        for name in self.infiles:
            if len(maskstr) > 0: maskstr += " || "
            maskstr += ("mask('%s')" % (name))
        outimage.calcmask(maskstr,name="basketweaving",asdefault=True)
        
        self.realimage.close()
        self.imagimage.close()
        outimage.close()

    def finalize(self):
        pass

    def cleanup(self):
        # finalize image analysis tool
        if hasattr(self,'image') and self.image is not None:
            if self.image.isopen(): self.image.done()
        tools = ['convimage', 'imageorg', 'realimage', 'imagimage']
        for t in tools:
            if hasattr(self,t):
                v = getattr(self,t)
                if v and v.isopen():
                    v.done(remove=True)

        # remove tempolary files
        filelist = ['tmpmskname', 'tmpconvname', 'tmppolyname',
                    'tmprealname', 'tmpimagname']
        existing_files = []
        for s in filelist:
            if hasattr(self, s):
                f = getattr(self, s)
                if isinstance(f,list):
                    for g in f:
                        if os.path.exists(g):
                            existing_files.append(g)
                else:
                    if os.path.exists(f):
                        existing_files.append(f)
        # CAS-5410 Use private tools inside task scripts
        if len(existing_files) > 0:
            cu = utilstool()
            cu.removetable(existing_files)
    
