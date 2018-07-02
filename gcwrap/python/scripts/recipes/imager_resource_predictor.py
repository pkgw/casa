from taskinit import casalog
from imagerhelpers.input_parameters import ImagerParameters


def predict_tclean_mem_use(imager_params):
    """
    Takes an ImagerParameters object as input and estimates/predict memory use in tclean.

    :param imager_params: tclean parameters as ImagerParameters object

    :returns: a dictionary with memory use estimations
    """
    if not isinstance(imager_params, ImagerParameters):
        raise ValueError('Error: the input parameter must be a valid ImagerParameters'
                         'object, but its type is: {0}'.format(type(imager_params)))

    if not imager_params.checkParameters():
        raise ValueError('Error: this ImagerParameter objects has not been built correctly. '
                         'This should not have happened if the constructor of the object '
                         'finished without errors.')

    im_pars = imager_params.getImagePars()['0']
    imsize = im_pars['imsize']
    one_plane = _mem_per_plane(imsize)
    casalog.post('One image takes: {0:5.2f} MB'.format(one_plane))

    dec_pars = imager_params.getDecPars()['0']
    nterms = im_pars['nterms']
    decpeak = _predict_minor_cycle(dec_pars, one_plane, imsize, nterms)
    casalog.post('Estimated peak memory in minor cycle: {0:5.2f} MB ({1:5.1f} images)'.
                 format(decpeak, decpeak/one_plane))

    grid_pars = imager_params.getGridPars()['0']
    ftmpeak = _predict_major_cycle(grid_pars, one_plane, im_pars['nchan'])

    if dec_pars['deconvolver'] == 'mtmfs':
        casalog.post('Estimated peak memory in major cycle: {0:5.2f} MB ({1:5.1f} '
                     'images) for the PSF and {2:5.2f} MB ({3:5.1f} images) thereafter'.
                     format(ftmpeak*(2*nterms-1), ftmpeak*(2*nterms-1)/one_plane,
                            ftmpeak*(nterms), ftmpeak*(nterms)/one_plane))
        # (2*nterms-1)  times for the PSF calc, 2 times otherwise
        # It's less than a factor of 2. Maybe the complex image grids disappear after each
        # FTM
        ftmpeak = ftmpeak * nterms
    else:
        casalog.post('Estimated peak memory in major cycle: {0:5.2f} MB ({1:5.1f} images)'.
                     format(ftmpeak, ftmpeak/one_plane))

    result = {'peak': ftmpeak+decpeak,
              'peak_in_minor_cycle': decpeak,
              'peak_in_major_cycle': ftmpeak,
              'mem_per_image': one_plane}
    casalog.post('Result from memory predictor: {0}'.format(result))
    return result


def _mem_per_plane(imsize):
    """
    Get MB per image/plane.

    :param imsize: image/plane size in pixels

    :returns: memory per image/plane in MB
    """
    size_float = 4  # bytes
    bytes_per_mb = 1024*1024.0
    # backgroundmem = 150 # ~150MB mem use at casa startup

    # Peak memory should be a multiple of a single image plane, depending on algorithm
    one_plane = (imsize * imsize * 1 * size_float) / bytes_per_mb
    return one_plane


def _predict_minor_cycle(dec_pars, one_plane, imsize, nterms):
    """
    Estimate memory peak for the minor cycles (deconvolvers). This adds up on top of the
    major cycle peaks (see _predict_major_cycle).

    :param dec_pars: dict of deconvolver parameters for tclean, as you would get from
                     ImagerParameters.getDecPars()
    :param one_plane: memory for one plane/image, in MB
    :param imsize: image/plane size in pixels, from image parameters, needed for mtmfs
    :param nterms: nterms from image parameters, needed for mtmfs

    :returns: memory use for the deconvolver / minor cycles, in MB
    """
    # these need to come from outside later
    nscales = 1
    ntotal = 0

    if dec_pars['deconvolver'] == 'hogbom':
        npsf = 1
        nresidual = 1
        nmodel = 1
        nmask = 1
        ndeltamodel = 0
        misc = 0

        ntotal = npsf + nresidual + nmodel + nmask + ndeltamodel + misc

    if dec_pars['deconvolver'] == 'multiscale':
        # itsPsfConvScales , #copy of input PSF
        npsf = nscales + 1
        # itsDirty, itsXfr(complex), itsDirtyConvScales, copy of input dirty
        nresidual = 1 + 2 + nscales + 1
        nmodel = 1
        nmask = 1 + 1 + nscales  # itsMask, itsScaleMasks, copy of input mask
        misc = nscales*2  # itsScales, itsScalesXfrs ,
        # Wow.... MatrixCleaner COPIES all of the images it gets ?!!!!

        nfftserver = 1 + 2*2  # Created for each FFT - so lots of spikes, but not present in accumulated mem use.

        ntotal = npsf + nresidual + nmodel + nmask + misc

    if dec_pars['deconvolver'] == 'mtmfs':
        # psf patches in the Hessian
        n4d = (nscales * (nscales+1) / 2.0) * (nterms * (nterms+1) / 2)

        nsupport = ((100.0/imsize)**2 * (n4d + nscales))

        nfull = 2 + 2 + 3 * nscales + 3 * nterms + (2 * nterms - 1) + 2 * nterms * nscales
        nfftserver = 1 + 2*2  # 1 float and 2 complex

        mystery = 1 + 1  # TODO

        ntotal = nsupport + nfull + nfftserver + mystery

    if dec_pars['deconvolver'] == 'clark':
        npsf = 1
        nresidual = 2
        nmodel = 1
        nmask = 1
        ndeltamodel = 1

        ntotal = npsf + nresidual + nmodel + nmask + ndeltamodel

    # Deconvolver total
    decpeak = one_plane * ntotal

    return decpeak


def _predict_major_cycle(grid_pars, one_plane, nchan):
    """
    Estimate memory peak for the major cycles (gridder).

    :param grid_pars: dict of gridder parameters for tclean, as you would get from
                      ImagerParameters.getGridPars()
    :param one_plane: memory for one plane/image, in MB
    :param nchan: number of channels, from the image parameters

    :returns: memory use for the gridder / major cycles, per channel, in MB
    """
    ntotalgrid = 0

    if grid_pars['gridder'] == 'gridft':
        # Gridder peak memory use is in getImage when copying double to float grids
        ngrid = 6  # griddedData(float complex)=2  + griddedData2(double complex)=4
        # De-gridder uses only float grid.
        ndegrid = 2  # griddedData(float complex)=2

        # TODO: no idea. needed to make the prediction match measurements.
        # Needs tracking down.
        nmystery = 1

        ntotalgrid = ngrid + ndegrid + nmystery

    if grid_pars['gridder'] == 'wproject':
        ngrid = 6
        ndegrid = 2
        ntotalgrid = ngrid + ndegrid

    if grid_pars['gridder'] == 'mosaic':
        ngrid = 6
        ndegrid = 2
        ntotalgrid = ngrid + ndegrid

    if grid_pars['gridder'] == 'awproject':
        ngrid = 6
        ndegrid = 2
        ntotalgrid = ngrid + ndegrid

    ftmpeak = ntotalgrid * one_plane * nchan
    return ftmpeak
