# sd task for imaging
import os
import re
import numpy
import shutil
import contextlib
import functools

from taskinit import casalog, gentools, qatool

# import sdutil
import sdbeamutil
from cleanhelper import cleanhelper

## (1) Import the python application layer
from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

image_suffix = '.residual'
weight_suffix = '.weight'
associate_suffixes = ['.psf', '.sumwt', weight_suffix]

def _configure_spectral_axis(mode, nchan, start, width, restfreq):
    # TODO: implement the function
    imnchan = nchan
    imstart = start
    imwidth = width
    return imnchan, imstart, imwidth

def _handle_grid_defaults(value):
    ret = ''
    if isinstance(value, int) or isinstance(value, float):
        ret = str(value)
    elif isinstance(value, str):
        ret = value
    return ret
   
def _get_param(ms_index, param):
    if isinstance(param, str):
        return param
    elif hasattr(param, '__iter__'):
        if len(param) == 1:
            return param[0]
        else:
            return param[ms_index]
    else:
        raise RuntimeError('Invalid parameter')

def _remove_image(imagename):
    if os.path.exists(imagename):
        if os.path.isdir(imagename):
            shutil.rmtree(imagename)
        elif os.path.isfile(imagename):
            os.remove(imagename)
        else:
            # could be a symlink
            os.remove(imagename)
            
def _get_restfreq_if_empty(vislist, spw, restfreq):
    qa = qatool()
    rf = None
    # if restfreq is nonzero float value, return it
    if isinstance(restfreq, float):
        if restfreq != 0.0:
            rf = restfreq
    # if restfreq is valid frequency string, return it
    # numeric string is interpreted as a value in the unit of Hz
    elif isinstance(restfreq, str):
        q = qa.convert(qa.quantity(restfreq), 'Hz')
        if q['unit'] == 'Hz' and q['value'] > 0.0:
            rf = restfreq
    # if restfreq is valid quantity, return it
    elif isinstance(restfreq, dict):
        q = qa.convert(restfreq, 'Hz')
        if q['unit'] == 'Hz' and q['value'] > 0.0:
            rf = restfreq
            
    if isinstance(vislist, str):
        vis = vislist
    elif hasattr(vislist, '__iter__'):
        vis = vislist[0]
    else:
        raise RuntimeError('Internal Error: invalid vislist \'{0}\''.format(vislist))
    
    if isinstance(spw, str):
        spwsel = spw
    elif hasattr(spw, '__iter__'):
        spwsel = spw[0]
    else:
        raise RuntimeError('Internal Error: invalid spw selection \'{0}\''.format(spw))
    
    with open_ms(vis) as ms:
        ms.msselect({'spw': spwsel})
        ndx = ms.msselectedindices()
        spwid = ndx['spw'][0]
        
    if rf is None:
        # if restfrequency is defined in SOURCE table, return it
        with open_table(os.path.join(vis, 'SOURCE')) as tb:
            if 'REST_FREQUENCY' in tb.colnames():
                tsel = tb.query('SPECTRAL_WINDOW_ID == {0}'.format(spwid))
                try:
                    nrow = tsel.nrows()
                    if nrow > 0:
                        for irow in xrange(nrow):
                            if tsel.iscelldefined('REST_FREQUENCY', irow):
                                rf = tsel.getcell('REST_FREQUENCY', irow)[0]
                                break
                finally:
                    tsel.close()
                            
    if rf is None:
        # otherwise, return mean frequency of given spectral window
        with open_table(os.path.join(vis, 'SPECTRAL_WINDOW')) as tb:
            cf = tb.getcell('CHAN_FREQ', spwid)
            restfreq = cf.mean()    

    assert rf is not None
    
    return rf

class OldImagerBasedTools(object):
    def __init__(self):
        self.imager = gentools(['im'])[0]
        
    @contextlib.contextmanager    
    def open_old_imager(self, vis):
        try:
            self.imager.open(vis)
            yield self.imager
        finally:
            self.imager.close()
        
    def test(self, vis):
        with self.open_old_imager(vis) as im:
            print 'test'
            raise RuntimeError('ERROR!')
        
    def get_pointing_sampling_params(self, vis, field, spw, baseline, scan, intent, outref, movingsource, pointingcolumntouse, antenna_name):
        with self.open_old_imager(vis) as im:
            im.selectvis(field=field,
                        spw=spw,
                        nchan=-1,
                        start=0,
                        step=1,
                        baseline=baseline,
                        scan=scan,
                        intent=intent)
            sampling_params = im.pointingsampling(pattern='raster',
                                                ref=outref,
                                                movingsource=movingsource,
                                                pointingcolumntouse=pointingcolumntouse,
                                                antenna='{0}&&&'.format(antenna_name))
        return sampling_params

@contextlib.contextmanager
def open_ia(imagename):
    (ia,) = gentools(['ia'])
    ia.open(imagename)
    try:
        yield ia
    finally:
        ia.close()

@contextlib.contextmanager
def open_ms(vis):
    (ms,) = gentools(['ms'])
    ms.open(vis)
    try:
        yield ms
    finally:
        ms.close()

@contextlib.contextmanager
def open_table(vis):
    (tb,) = gentools(['tb'])
    tb.open(vis)
    try:
        yield tb
    finally:
        tb.close()
        
def set_beam_size(vis, imagename,
                  field, spw, baseline, scan, intent,
                  ephemsrcname, pointingcolumntouse, antenna_name, antenna_diameter,
                  restfreq, gridfunction, convsupport, truncate, gwidth, jwidth):
    """
    Set estimated beam size to the image.
    """
    is_alma = antenna_name[0:2] in ['PM', 'DV', 'DA', 'CM']
    blockage = '0.75m' if is_alma else '0.0m' 
    
    with open_ia(imagename) as ia:
        csys = ia.coordsys()
        outref = csys.referencecode('direction')[0]
        cell = list(csys.increment(type='direction', format='s')['string'])

    old_tool = OldImagerBasedTools()
    sampling_params = old_tool.get_pointing_sampling_params(vis, field, spw, baseline, scan, intent,
                                                        outref=outref,
                                                        movingsource=ephemsrcname,
                                                        pointingcolumntouse=pointingcolumntouse,
                                                        antenna_name=antenna_name)
    qa = qatool()
    casalog.post('sampling_params={0}'.format(sampling_params))
    xsampling, ysampling = qa.getvalue(qa.convert(sampling_params['sampling'], 'arcsec'))
    angle = qa.getvalue(qa.convert(sampling_params['angle'], 'deg'))[0]
    
    casalog.post('Detected raster sampling = [{0:f}, {1:f}] arcsec'.format(xsampling, ysampling))
    
    # handling of failed sampling detection
    valid_sampling = True
    # TODO: copy from sdimaging implementation
    sampling = [xsampling, ysampling]
    if abs(xsampling) < 2.2e-3 or not numpy.isfinite(xsampling):
        casalog.post("Invalid sampling=%s arcsec. Using the value of orthogonal direction=%s arcsec" % (xSampling, ySampling), priority="WARN")
        sampling = [ ysampling ]
        angle = 0.0
        valid_sampling = False
    if abs(ysampling) < 1.0e-3 or not numpy.isfinite(ysampling):
        if valid_sampling:
            casalog.post("Invalid sampling=%s arcsec. Using the value of orthogonal direction=%s arcsec" % (ySampling, xSampling), priority="WARN")
            sampling = [ xsampling ]
            angle = 0.0
            valid_sampling = True
    # reduce sampling and cell if it's possible
    if len(sampling) > 1 and abs(sampling[0] - sampling[1]) <= 0.01 * abs(sampling[0]):
        sampling = [sampling[0]]
        angle = 0.0
        if cell[0] == cell[1]: cell = [cell[0]]
    if valid_sampling:
        # actual calculation of beam size
        bu = sdbeamutil.TheoreticalBeam()
        bu.set_antenna(antenna_diameter, blockage)
        bu.set_sampling(sampling, "%fdeg" % angle)
        bu.set_image_param(cell, restfreq, gridfunction,
                           convsupport, truncate, gwidth,
                           jwidth, is_alma)
        bu.summary()
        imbeam_dict = bu.get_beamsize_image()
        casalog.post("Setting image beam: major=%s, minor=%s, pa=%s" % 
                     (imbeam_dict['major'], imbeam_dict['minor'],
                      imbeam_dict['pa'],))
        # set beam size to image
        with open_ia(imagename) as ia:
            ia.setrestoringbeam(**imbeam_dict)
    else:
        # BOTH sampling was invalid
        casalog.post("Could not detect valid raster sampling. Exitting without setting beam size to image", priority='WARN')
    
def do_weight_mask(imagename, weightimage, minweight):
    # Mask image pixels whose weight are smaller than minweight.
    # Weight image should have 0 weight for pixels below < minweight
    casalog.post("Start masking the map using minweight = %f" % \
                 minweight, "INFO")
    with open_ia(weightimage) as ia:
        try:
            stat=ia.statistics(mask="'"+weightimage+"' > 0.0", robust=True)
            valid_pixels=stat['npts']
        except RuntimeError, e:
            if e.message.find('No valid data found.') >= 0:
                valid_pixels = [0]
            else:
                raise e
            
    if len(valid_pixels) == 0 or valid_pixels[0] == 0:
        casalog.post("All pixels weight zero. This indicates no data in MS is in image area. Mask will not be set. Please check your image parameters.","WARN")
        return
    median_weight = stat['median'][0]
    weight_threshold = median_weight * minweight
    casalog.post("Median of weight in the map is %f" % median_weight, \
                 "INFO")
    casalog.post("Pixels in map with weight <= median(weight)*minweight = %f will be masked." % \
                 (weight_threshold),"INFO")
    ###Leaving the original logic to calculate the number of masked pixels via
    ###product of median of and min_weight (which i don't understand the logic)
    ### if one wanted to find how many pixel were masked one could easily count the
    ### number of pixels set to false 
    ### e.g  after masking self.outfile below one could just do this 
    ### nmasked_pixels=tb.calc('[select from "'+self.outfile+'"/mask0'+'"  giving [nfalse(PagedArray )]]')
    my_tb = gentools(['tb'])[0]
    nmask_pixels=0
    nchan=stat['trc'][3]+1
    casalog.filter('ERROR') ### hide the useless message of tb.calc

   
    ### doing it by channel to make sure it does not go out of memory
    ####tab.calc try to load the whole chunk in ram 
    for k in range(nchan):
        nmask_pixels += my_tb.calc('[select from "'+weightimage+'"  giving [ntrue(map[,,,'+str(k)+'] <='+str(median_weight*minweight)+')]]')['0'][0]
    casalog.filter()  ####set logging back to normal
    
    # Modify default mask
    with open_ia(imagename) as ia:
        imsize=numpy.product(ia.shape())
        ia.calcmask("'%s'>%f" % (weightimage, weight_threshold), asdefault=True)

    masked_fraction = 100.*(1. - (imsize - nmask_pixels) / float(valid_pixels[0]) )
    casalog.post("This amounts to %5.1f %% of the area with nonzero weight." % \
                ( masked_fraction ),"INFO")
    casalog.post("The weight image '%s' is returned by this task, if the user wishes to assess the results in detail." \
                 % (weightimage), "INFO")
        

def tsdimaging(infiles, outfile, overwrite, field, spw, antenna, scan, intent, mode, nchan, start, width, veltype, outframe,
               gridfunction, convsupport, truncate, gwidth, jwidth, imsize, cell, phasecenter, projection, ephemsrcname,
               pointingcolumn, restfreq, stokes, minweight, brightnessunit, clipminmax):
    
    origin = 'tsdimaging'
    imager = None
 
    casalog.post('spw = \'{0}\' type {1}'.format(spw, type(spw)))
 
    try: 
        
        # handle overwrite parameter
        _outfile = outfile.rstrip('/')
        presumed_imagename = _outfile + image_suffix
        if os.path.exists(presumed_imagename):
            if overwrite == False:
                raise RuntimeError('Output file \'{0}\' exists.'.format(presumed_imagename))
            else:
                # delete existing images
                casalog.post('Removing \'{0}\''.format(presumed_imagename))
                _remove_image(presumed_imagename)
                assert not os.path.exists(presumed_imagename)
                for _suffix in associate_suffixes:
                    casalog.post('Removing \'{0}\''.format(_outfile + _suffix))
                    _remove_image(_outfile + _suffix)
                    assert not os.path.exists(_outfile + _suffix)
    
        # parse parameter for spectral axis 
        imnchan, imstart, imwidth = _configure_spectral_axis(mode, nchan, start, width, restfreq)
        _restfreq = _get_restfreq_if_empty(infiles, spw, restfreq)
        
        # translate some default values into the ones that are consistent with the current framework
        gtruncate = _handle_grid_defaults(truncate)
        ggwidth = _handle_grid_defaults(gwidth)
        gjwidth = _handle_grid_defaults(jwidth)
        
        ## (2) Set up Input Parameters 
        ##       - List all parameters that you need here
        ##       - Defaults will be assumed for unspecified parameters
        ##       - Nearly all parameters are identical to that in the task. Please look at the
        ##         list of parameters under __init__ using  " help ImagerParameters " )
        casalog.post('*** Creating paramList ***', origin=origin)
        paramList = ImagerParameters(
            # input file name
            msname =infiles,#'sdimaging.ms',
            # data selection
            field=field,#'',
            spw=spw,#'0',
            antenna=antenna,
            scan=scan,
            state=intent,
            # image parameters
            imagename=_outfile,#'try2',
            nchan=imnchan,#1024,
            start=imstart,#'0',
            width=imwidth,#'1',
            outframe=outframe,
            veltype=veltype,
            restfreq=_restfreq,
            phasecenter=phasecenter,#'J2000 17:18:29 +59.31.23',
            imsize=imsize,#[75,75], 
            cell=cell,#['3arcmin', '3arcmin'], 
            projection=projection,
            stokes=stokes,
            # fix specmode to 'cubedata'
            # output spectral coordinate will be determined based on mode, start, and width 
            specmode='cubedata',
            gridder='singledish',
            # single dish specific parameters
            gridfunction=gridfunction,
            convsupport=convsupport,
            truncate=gtruncate,
            gwidth=ggwidth,
            jwidth=gjwidth,
            pointingcolumntouse=pointingcolumn,
            minweight=minweight,
            clipminmax=clipminmax,
            # normalizer
            normtype='flatsky',
            pblimit=1e-16
        )
        
        # TODO: hadnle ephemsrcname
        # TODO: handle brightnessunit
        # TODO: handle overwrite
        # TODO: output image name
            
        ## (3) Construct the PySynthesisImager object, with all input parameters
        
        casalog.post('*** Creating imager object ***', origin=origin)
        imager = PySynthesisImager(params=paramList)
            
        ## (4) Initialize various modules. 
        ##       - Pick only the modules you will need later on. For example, to only make
        ##         the PSF, there is no need for the deconvolver or iteration control modules.
        
        ## Initialize modules major cycle modules
        
        casalog.post('*** Initializing imagers ***', origin=origin)
        imager.initializeImagers()
        casalog.post('*** Initializing normalizers ***', origin=origin)
        imager.initializeNormalizers()
        #imager.setWeighting()
        
        ## (5) Make the initial images 
        
        #imager.makePSF()
        casalog.post('*** Executing runMajorCycle ***', origin=origin)
        casalog.post('NF = {0}'.format(imager.NF), origin=origin)
        #imager.runMajorCycle()  # Make initial dirty / residual image
        imager.makeSdPSF()
        imager.makeSdImage()
        
    except Exception as e:
        #print 'Exception : ' + str(e)
        casalog.post('Exception from task_tsdimaging : ' + str(e), "SEVERE", origin=origin)
#         if imager != None:
#             imager.deleteTools() 

        larg = list(e.args)
        larg[0] = 'Exception from task_tsdimaging : ' + str(larg[0])
        e.args = tuple(larg)
        raise
    
    finally:
        ## (8) Close tools.
        
        casalog.post('*** Cleaning up tools ***', origin=origin)
        if imager is not None:
            imager.deleteTools()
        

    # set beam size
    # TODO: re-define related functions in the new tool framework (sdms?)
    imagename = outfile + image_suffix
    ms_index = 0
    rep_ms = _get_param(0, infiles)
    rep_field = _get_param(0, field)
    rep_spw = _get_param(0, spw)
    rep_antenna = _get_param(0, antenna)
    rep_scan = _get_param(0, scan)
    rep_intent = _get_param(0, intent)
    if len(rep_antenna) > 0:
        baseline = '{0}&&&'.format(rep_antenna)
    else:
        baseline = '*&&&'
    with open_ms(rep_ms) as ms:
        ms.msselect({'baseline': baseline})
        ndx = ms.msselectedindices()
        antenna_index = ndx['antenna1']
    with open_table(os.path.join(rep_ms, 'ANTENNA')) as tb:
        antenna_name = tb.getcell('NAME', antenna_index)
        antenna_diameter = tb.getcell('DISH_DIAMETER', antenna_index)
    set_beam_size(rep_ms, imagename,
                  rep_field, rep_spw, baseline, rep_scan, rep_intent,
                  ephemsrcname, pointingcolumn, antenna_name, antenna_diameter,
                  _restfreq, gridfunction, convsupport, truncate, gwidth, jwidth)
    
    # mask low weight pixels 
    #weightimage = outfile + weight_suffix
    #do_weight_mask(imagename, weightimage, minweight)
