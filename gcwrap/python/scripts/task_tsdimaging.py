# sd task for imaging
import os
import re
import numpy
import shutil
import contextlib
import functools

from taskinit import casalog, gentools, qatool

import sdutil
import sdbeamutil
from cleanhelper import cleanhelper

## (1) Import the python application layer
from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

image_suffix = '.image'
residual_suffix = '.residual'
weight_suffix = '.weight'
associate_suffixes = ['.psf', '.sumwt', weight_suffix, residual_suffix]

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
        
class SelectionHandler(object):
    def __init__(self, sel):
        self.sel = sel
        if isinstance(self.sel, str):
            self.selector = self._select0
        elif len(self.sel) == 1:
            self.selector = self._select1
        else:
            self.selector = self._select2
    
    def __call__(self, i):
        return self.selector(i)
        
    def _select0(self, i):
        return self.sel

    def _select1(self, i):
        return self.sel[0]
    
    def _select2(self, i):
        return self.sel[i]
    
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
            
    @contextlib.contextmanager
    def open_and_select_old_imager(self, vislist, field, spw, antenna, scan, intent):
        if isinstance(vislist, str):
            with self.open_old_imager(vislist) as im:
                im.selectvis(field=field,
                             spw=spw,
                             nchan=-1,
                             start=0,
                             step=1,
                             baseline=antenna,
                             scan=scan,
                             intent=intent)
                yield im
        else:
            fieldsel = SelectionHandler(field)
            spwsel = SelectionHandler(spw)
            antennasel = SelectionHandler(antenna)
            scansel = SelectionHandler(scan)
            intentsel = SelectionHandler(intent)
            try:
                for i in range(len(vislist)):
                    vis = vislist[i]
                    _field = fieldsel(i)
                    _spw = spwsel(i)
                    _antenna = antennasel(i)
                    _scan = scansel(i)
                    _intent = intentsel(i)
                    if len(_antenna) == 0:
                        _baseline = _antenna
                    elif len(_antenna) < 4 or _antenna[:-3] != '&&&':
                        _baseline = _antenna + '&&&'
                    else:
                        _baseline = _antenna
                    self.imager.selectvis(vis, field=_field, spw=_spw, nchan=-1, start=0, step=1,
                                          baseline=_baseline, scan=_scan, intent=_intent)
                yield self.imager
            finally:
                self.imager.close()
        
    def test(self, vis):
        with self.open_old_imager(vis) as im:
            print('test')
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
    
    def get_map_extent(self, vislist, field, spw, antenna, scan, intent, 
                       ref, movingsource, pointingcolumntouse):
        
        with self.open_and_select_old_imager(vislist=vislist, field=field,
                                             spw=spw, antenna=antenna, scan=scan, 
                                             intent=intent) as im:
            map_param = im.mapextent(ref=ref, movingsource=movingsource, 
                                     pointingcolumntouse=pointingcolumntouse)
        return map_param
    
    def sort_vis(self, vislist, spw, mode, width, field, antenna, scan, intent):
        if isinstance(vislist, str) or len(vislist) == 1:
            return vislist, field, spw, antenna, scan, intent
        imhelper = cleanhelper(imtool=self.imager, vis=vislist, casalog=casalog)
        imhelper.sortvislist(spw=spw, mode=mode, width=width)
        sorted_idx = list(imhelper.sortedvisindx)
        # reverse the order
        sorted_idx.reverse()
        sorted_vislist = [vislist[i] for i in sorted_idx]
        fieldsel = SelectionHandler(field)
        sorted_field = [fieldsel(i) for i in sorted_idx]
        spwsel = SelectionHandler(spw)
        sorted_spw = [spwsel(i) for i in sorted_idx]
        antennasel = SelectionHandler(antenna)
        sorted_antenna = [antennasel(i) for i in sorted_idx]
        scansel = SelectionHandler(scan)
        sorted_scan = [scansel(i) for i in sorted_idx]
        intentsel = SelectionHandler(intent)
        sorted_intent = [intentsel(i) for i in sorted_idx]
        return sorted_vislist, sorted_field, sorted_spw, sorted_antenna, sorted_scan, sorted_intent

def _configure_spectral_axis(mode, nchan, start, width, restfreq):
    # fix default
    if mode == 'channel':
        if start == '': start = 0
        if width == '': width = 1
    else:
        if start == 0: start = ''
        if width == 1: width = ''
    # fix unit
    if mode == 'frequency':
        myunit = 'Hz'
    elif mode == 'velocity':
        myunit = 'km/s'
    else: # channel
        myunit = ''

    tmp_start = _format_quantum_unit(start, myunit)
    if tmp_start == None:
        raise ValueError("Invalid unit for %s in mode %s: %s" % ('start', mode, start))
    start = tmp_start
    if mode == 'channel':
        start = int(start)
    tmp_width = _format_quantum_unit(width, myunit)
    if tmp_width == None:
        raise ValueError("Invalid unit for %s in mode %s: %s" % ('width', mode, width))
    width = tmp_width
    if mode == 'channel':
        width = int(width)

    #TODO: work for nchan
    imnchan = nchan
    imstart = start
    imwidth = width
    return imnchan, imstart, imwidth

def _format_quantum_unit(data, unit):
    """
    Returns False if data has an unit which in not a variation of
    input unit.
    Otherwise, returns input data as a quantum string. The input
    unit is added to the return value if no unit is in data.
    """
    my_qa = qatool()
    if data == '' or my_qa.compare(data, unit):
        return data
    if my_qa.getunit(data) == '':
        casalog.post("No unit specified. Using '%s'" % unit)
        return '%f%s' % (data, unit)
    return None

def _handle_grid_defaults(value):
    ret = ''
    if isinstance(value, int) or isinstance(value, float):
        ret = str(value)
    elif isinstance(value, str):
        ret = value
    return ret

def _calc_PB(vis, antenna_id, restfreq):
    """
    Calculate the primary beam size of antenna, using dish diamenter
    and rest frequency
    Average antenna diamter and reference frequency are adopted for
    calculation.
    The input argument should be a list of antenna IDs.
    """
    casalog.post("Calculating Pirimary beam size:")
    # CAS-5410 Use private tools inside task scripts
    my_qa = qatool()
    
    pb_factor = 1.175
    # Reference frequency
    ref_freq = restfreq
    if type(ref_freq) in [float, numpy.float64]:
        ref_freq = my_qa.tos(my_qa.quantity(ref_freq, 'Hz'))
    if not my_qa.compare(ref_freq, 'Hz'):
        msg = "Could not get the reference frequency. " + \
              "Your data does not seem to have valid one in selected field.\n" + \
              "PB is not calculated.\n" + \
              "Please set restreq or cell manually to generate an image."
        raise Exception(msg)
    # Antenna diameter
    with open_table(os.path.join(vis, 'ANTENNA')) as tb:
        antdiam_ave = tb.getcell('DISH_DIAMETER', antenna_id)
    #antdiam_ave = self._get_average_antenna_diameter(antenna)
    # Calculate PB
    wave_length = 0.2997924 / my_qa.convert(my_qa.quantity(ref_freq),'GHz')['value']
    D_m = my_qa.convert(antdiam_ave, 'm')['value']
    lambda_D = wave_length / D_m * 3600. * 180 / numpy.pi
    PB = my_qa.quantity(pb_factor*lambda_D, 'arcsec')
    # Summary
    casalog.post("- Antenna diameter: %s m" % D_m)
    casalog.post("- Reference Frequency: %s" % ref_freq)
    casalog.post("PB size = %5.3f * lambda/D = %s" % (pb_factor, my_qa.tos(PB)))
    return PB

def _get_imsize(width, height, dx, dy):
    casalog.post("Calculating pixel size.")
    # CAS-5410 Use private tools inside task scripts
    my_qa = qatool()
    ny = numpy.ceil( ( my_qa.convert(height, my_qa.getunit(dy))['value'] /  \
                       my_qa.getvalue(dy) ) )
    nx = numpy.ceil( ( my_qa.convert(width, my_qa.getunit(dx))['value'] /  \
                       my_qa.getvalue(dx) ) )
    casalog.post("- Map extent: [%s, %s]" % (my_qa.tos(width), my_qa.tos(height)))
    casalog.post("- Cell size: [%s, %s]" % (my_qa.tos(dx), my_qa.tos(dy)))
    casalog.post("Image pixel numbers to cover the extent: [%d, %d] (projected)" % \
                 (nx+1, ny+1))
    return [int(nx+1), int(ny+1)]

def _get_pointing_extent(phasecenter, vislist, field, spw, antenna, scan, intent, 
                         pointingcolumntouse, ephemsrcname):
    ### MS selection is ignored. This is not quite right.
    casalog.post("Calculating map extent from pointings.")
    # CAS-5410 Use private tools inside task scripts
    my_qa = qatool()
    ret_dict = {}
    
    if isinstance(vislist, str):
        vis = vislist
    else:
        vis = vislist[0]
    
    colname = pointingcolumntouse.upper()

    if phasecenter == "":
        # defaut is J2000
        base_mref = 'J2000'
    elif isinstance(phasecenter, int) or phasecenter.isdigit():
        # may be field id
        with open_table(os.path.join(vis, 'FIELD')) as tb:
            base_mref = tb.getcolkeyword('PHASE_DIR', 'MEASINFO')['Ref']
    else:
        # may be phasecenter is explicitly specified
        pattern = '^([\-\+]?[0-9.]+([eE]?-?[0-9])?)|([\-\+]?[0-9][:h][0-9][:m][0-9.]s?)|([\-\+]?[0-9][.d][0-9][.d][0-9.]s?)$'
        items = phasecenter.split()
        base_mref = 'J2000'
        for i in items:
            s = i.strip()
            if re.match(pattern, s) is None:
                base_mref = s
                break

    t = OldImagerBasedTools()
    mapextent = t.get_map_extent(vislist, field, spw, antenna, scan, intent, 
                                 ref=base_mref, movingsource=ephemsrcname, 
                                 pointingcolumntouse=pointingcolumntouse)
    #mapextent = self.imager.mapextent(ref=base_mref, movingsource=ephemsrcname, 
    #                                  pointingcolumntouse=colname)
    if mapextent['status'] is True:
        qheight = my_qa.quantity(mapextent['extent'][1], 'rad')
        qwidth = my_qa.quantity(mapextent['extent'][0], 'rad')
        qcent0 = my_qa.quantity(mapextent['center'][0], 'rad')
        qcent1 = my_qa.quantity(mapextent['center'][1], 'rad')
        scenter = '%s %s %s'%(base_mref, my_qa.formxxx(qcent0, 'hms'), 
                              my_qa.formxxx(qcent1, 'dms'))

        casalog.post("- Pointing center: %s" % scenter)
        casalog.post("- Pointing extent: [%s, %s] (projected)" % (my_qa.tos(qwidth), \
                                                              my_qa.tos(qheight)))
        ret_dict['center'] = scenter
        ret_dict['width'] = qwidth
        ret_dict['height'] = qheight
    else:
        casalog.post('Failed to derive map extent from the MSs registered to the imager probably due to mising valid data.', priority='SEVERE')
        ret_dict['center'] = ''
        ret_dict['width'] = my_qa.quantity(0.0, 'rad')
        ret_dict['height'] = my_qa.quantity(0.0, 'rad')
    return ret_dict

def _handle_image_params(imsize, cell, phasecenter, 
                         vislist, field, spw, antenna, scan, intent, 
                         restfreq, pointingcolumntouse, ephemsrcname):
    # round-up imsize
    _imsize = sdutil._to_list(imsize, int) or sdutil._to_list(imsize, numpy.integer)
    if _imsize is None:
        _imsize = imsize if hasattr(imsize, '__iter__') else [ imsize ]
        _imsize = [ int(numpy.ceil(v)) for v in _imsize ]
        casalog.post("imsize is not integers. force converting to integer pixel numbers.", priority="WARN")
        casalog.post("rounded-up imsize: %s --> %s" % (str(imsize), str(_imsize)))

    # calculate cell based on PB if it is not given
    _cell = cell
    if _cell == '' or _cell[0] == '':
        # calc PB
        if isinstance(vislist, str):
            vis = vislist
        else:
            vis = vislist[0]
        if isinstance(antenna, str):
            antsel = antenna
        else:
            antsel = antenna[0]
        if antsel == '':
            antenna_id = 0
        else:
            if len(antsel) > 3 and antsel[:-3] == '&&&':
                baseline = antsel
            else:
                baseline = antsel + '&&&'
            with open_ms(vis) as ms:
                ms.msselect({'baseline': baseline})
                ndx = ms.msselectedindices()
            antenna_id = ndx['antenna1'][0]
        grid_factor = 3.
        casalog.post("The cell size will be calculated using PB size of antennas in the first MS")
        qpb = _calc_PB(vis, antenna_id, restfreq)
        _cell = '%f%s' % (qpb['value']/grid_factor, qpb['unit'])
        casalog.post("Using cell size = PB/%4.2F = %s" % (grid_factor, _cell))
        
    # Calculate Pointing center and extent (if necessary)
    _phasecenter = phasecenter
    if _phasecenter == '' or len(_imsize) == 0 or _imsize[0] < 1:
        # return a dictionary with keys 'center', 'width', 'height'
        map_param = _get_pointing_extent(_phasecenter, vislist, field, spw, antenna, scan, intent, 
                                         pointingcolumntouse, ephemsrcname)
        # imsize
        (cellx,celly) = sdutil.get_cellx_celly(_cell, unit='arcmin')
        if len(_imsize) == 0 or _imsize[0] < 1:
            _imsize = _get_imsize(map_param['width'], map_param['height'], cellx, celly)
            if _phasecenter != "":
                casalog.post("You defined phasecenter but not imsize. The image will cover as wide area as pointing in MS extends, but be centered at phasecenter. This could result in a strange image if your phasecenter is a part from the center of pointings", priority='WARN')
            if _imsize[0] > 1024 or _imsize[1] > 1024:
                casalog.post("The calculated image pixel number is larger than 1024. It could take time to generate the image depending on your computer resource. Please wait...", priority='WARN')

        # phasecenter
        # if empty, it should be determined here...
        if _phasecenter == "":
            _phasecenter = map_param['center']
    
    return _imsize, _cell, _phasecenter

def _calc_pblimit(minweight):
    if minweight == 0.0:
        # set tiny value
        pblimit = 1e-16
    else:
        pblimit = minweight
        
    # disable pixel mask by pblimit
    pblimit = 1e-16 
    return pblimit
   
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
            
def _get_restfreq_if_empty(vislist, spw, field, restfreq):
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
    
    if isinstance(field, str):
        fieldsel = field
    elif hasattr(field, '__iter__'):
        fieldsel = field[0]
    else:
        raise RuntimeError('Internal Error: invalid field selection \'{0}\''.format(field))
       
    
    with open_ms(vis) as ms:
        ms.msselect({'spw': spwsel, 'field': fieldsel})
        ndx = ms.msselectedindices()
        if len(ndx['spw']) > 0:
            spwid = ndx['spw'][0]
        else:
            spwid = None
        if len(ndx['field']) > 0:
            fieldid = ndx['field'][0]
        else:
            fieldid = None
    sourceid = None
    if fieldid is not None:
        with open_table(os.path.join(vis, 'FIELD')) as tb:
            sourceid = tb.getcell('SOURCE_ID', fieldid)
        if sourceid < 0:
            sourceid = None
    if rf is None:
        # if restfrequency is defined in SOURCE table, return it
        with open_table(os.path.join(vis, 'SOURCE')) as tb:
            if 'REST_FREQUENCY' in tb.colnames():
                tsel = None
                taql = ''
                if spwid is not None:
                    taql = 'SPECTRAL_WINDOW_ID == {0}'.format(spwid)
                if sourceid is not None:
                    delimiter = '&&' if len(taql) > 0 else ''
                    taql += '{0}SOURCE_ID == {1}'.format(delimiter, sourceid)
                if len(taql) > 0:
                    tsel = tb.query(taql)
                    t = tsel
                else:
                    t = tb
                try:
                    nrow = t.nrows()
                    if nrow > 0:
                        for irow in range(nrow):
                            if t.iscelldefined('REST_FREQUENCY', irow):
                                rfs = t.getcell('REST_FREQUENCY', irow)
                                if len(rfs) > 0:
                                    rf = rfs[0]
                                    break
                finally:
                    if tsel is not None:
                        tsel.close()
                            
    if rf is None:
        if spwid is None:
            spwid = 0
        # otherwise, return mean frequency of given spectral window
        with open_table(os.path.join(vis, 'SPECTRAL_WINDOW')) as tb:
            cf = tb.getcell('CHAN_FREQ', spwid)
            rf = cf.mean()
            
    assert rf is not None
    
    return rf

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
        csys.done()

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
        casalog.post("Invalid sampling=%s arcsec. Using the value of orthogonal direction=%s arcsec" % (xsampling, ysampling), priority="WARN")
        sampling = [ ysampling ]
        angle = 0.0
        valid_sampling = False
    if abs(ysampling) < 1.0e-3 or not numpy.isfinite(ysampling):
        if valid_sampling:
            casalog.post("Invalid sampling=%s arcsec. Using the value of orthogonal direction=%s arcsec" % (ysampling, xsampling), priority="WARN")
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
        except RuntimeError as e:
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
        
def get_ms_column_unit(tb, colname):
    col_unit = ''
    if colname in tb.colnames():
        cdkw = tb.getcoldesc(colname)['keywords']
        if 'QuantumUnits' in cdkw:
            u = cdkw['QuantumUnits']
            if isinstance(u, str):
                col_unit = u.strip()
            elif isinstance(u, list):
                col_unit = u[0].strip()
    return col_unit

def get_brightness_unit_from_ms(msname):
    image_unit = ''
    with open_table(msname) as tb:
        image_unit = get_ms_column_unit(tb, 'DATA')
        if image_unit == '': image_unit = get_ms_column_unit(tb, 'FLOAT_DATA')
    if image_unit.upper() == 'K':
        image_unit = 'K'
    else:
        image_unit = 'Jy/beam'

    return image_unit



def tsdimaging(infiles, outfile, overwrite, field, spw, antenna, scan, intent, mode, nchan, start, width, veltype, outframe,
               gridfunction, convsupport, truncate, gwidth, jwidth, imsize, cell, phasecenter, projection, 
               pointingcolumn, restfreq, stokes, minweight, brightnessunit, clipminmax):
    
    origin = 'tsdimaging'
    imager = None
  
    try: 
        # if spw starts with ':', add '*' at the beginning
        if isinstance(spw, str):
            _spw = '*' + spw if spw.startswith(':') else spw
        else:
            _spw = ['*' + v if v.startswith(':') else v for v in spw]
            
        # if antenna doesn't contain '&&&', append it
        def antenna_to_baseline(s):
            if len(s) == 0:
                return s
            elif len(s) > 3 and s.endswith('&&&'):
                return s
            else:
                return '{0}&&&'.format(s)
        if isinstance(antenna, str):
            baseline = antenna_to_baseline(antenna)
        else:
            baseline = [antenna_to_baseline(a) for a in antenna]
            
        
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
        _restfreq = _get_restfreq_if_empty(infiles, _spw, field, restfreq)
        
        # translate some default values into the ones that are consistent with the current framework
        gtruncate = _handle_grid_defaults(truncate)
        ggwidth = _handle_grid_defaults(gwidth)
        gjwidth = _handle_grid_defaults(jwidth)
        
        _ephemsrcname = ''
        if isinstance(phasecenter, str) and phasecenter.strip().upper() in ['MERCURY', 'VENUS', 'MARS', 'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE', 'PLUTO', 'SUN', 'MOON', 'TRACKFIELD']:
            _ephemsrcname = phasecenter

        # handle image parameters
        if isinstance(infiles, str) or len(infiles) == 1:
            _imsize, _cell, _phasecenter = _handle_image_params(imsize, cell, phasecenter, infiles, 
                                                                field, _spw, antenna, scan, intent,
                                                                _restfreq, pointingcolumn, _ephemsrcname)
        else:
            # sort input data using cleanhelper function to get consistent result with older sdimaging
            o = OldImagerBasedTools()
            _sorted = o.sort_vis(infiles, _spw, mode, imwidth, field, antenna, scan, intent)
            sorted_vis, sorted_field, sorted_spw, sorted_antenna, sorted_scan, sorted_intent = _sorted
            _imsize, _cell, _phasecenter = _handle_image_params(imsize, cell, phasecenter, sorted_vis, 
                                                                sorted_field, sorted_spw, sorted_antenna, 
                                                                sorted_scan, sorted_intent,
                                                                _restfreq, pointingcolumn, _ephemsrcname)

        # calculate pblimit from minweight
        pblimit = _calc_pblimit(minweight)
        
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
            spw=_spw,#'0',
            antenna=baseline,
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
            phasecenter=_phasecenter,#'J2000 17:18:29 +59.31.23',
            imsize=_imsize,#[75,75], 
            cell=_cell,#['3arcmin', '3arcmin'], 
            projection=projection,
            stokes=stokes,
            # fix specmode to 'cubedata'
            # output spectral coordinate will be determined based on mode, start, and width 
            specmode='cube',
            #specmode='cubesource',
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
            pblimit=pblimit
        )
        
        # handle brightnessunit (CAS-11503)
        image_unit = ''
        if len(brightnessunit) > 0:
            if brightnessunit.lower() == 'k':
                image_unit = 'K'
            elif brightnessunit.lower() == 'jy/beam':
                image_unit = 'Jy/beam'
            else:
                raise ValueError("Invalid brightness unit, %s" % brightnessunit)
        
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
        casalog.post('*** makeSdPSF... ***', origin=origin)
        imager.makeSdPSF()
        casalog.post('*** makeSdImage... ***', origin=origin)
        imager.makeSdImage()
        
    except Exception as e:
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
            
        # change image suffix from .residual to .image
        if os.path.exists(outfile + residual_suffix):
            os.rename(outfile + residual_suffix, outfile + image_suffix)
        

    # set beam size
    # TODO: re-define related functions in the new tool framework (sdms?)
    imagename = outfile + image_suffix
    ms_index = 0
    rep_ms = _get_param(0, infiles)
    rep_field = _get_param(0, field)
    rep_spw = _get_param(0, _spw)
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
        antenna_index = ndx['antenna1'][0]
    with open_table(os.path.join(rep_ms, 'ANTENNA')) as tb:
        antenna_name = tb.getcell('NAME', antenna_index)
        antenna_diameter = tb.getcell('DISH_DIAMETER', antenna_index)
    set_beam_size(rep_ms, imagename,
                  rep_field, rep_spw, baseline, rep_scan, rep_intent,
                  _ephemsrcname, pointingcolumn, antenna_name, antenna_diameter,
                  _restfreq, gridfunction, convsupport, truncate, gwidth, jwidth)
    
    # set brightness unit (CAS-11503)
    if len(image_unit) == 0:
        image_unit = get_brightness_unit_from_ms(rep_ms)
    if len(image_unit) > 0:
        with open_ia(imagename) as ia:
            casalog.post("Setting brightness unit '%s' to image." % image_unit)
            ia.setbrightnessunit(image_unit)

    # mask low weight pixels 
    weightimage = outfile + weight_suffix
    do_weight_mask(imagename, weightimage, minweight)
