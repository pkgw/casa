# sd task for imaging
import os
import re
import numpy
import shutil

from taskinit import casalog, gentools, qatool

#import sdutil
#import sdbeamutil
from cleanhelper import cleanhelper

## (1) Import the python application layer
from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

image_suffix = '.residual'
associate_suffixes = ['.psf', '.sumwt', '.weight']

def _configure_spectral_axis(mode, nchan, start, width, restfreq):
    #TODO: implement the function
    imnchan = nchan
    imstart = start
    imwidth = width
    return imnchan, imstart, imwidth

def _handle_grid_defaults(value):
    ret = ''
    if isinstance(value, int) or isinstance(value, float):
        if value != -1:
            ret = string(value)
    elif isinstance(value, str):
        ret = value
    return ret

def _remove_image(imagename):
    if os.path.exists(imagename):
        if os.path.isdir(imagename):
            shutil.rmtree(imagename)
        elif os.path.isfile(imagename):
            os.remove(imagename)
        else:
            # could be a symlink
            os.remove(imagename)

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
            restfreq=restfreq,
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
            normtype='flatsky'
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
        

