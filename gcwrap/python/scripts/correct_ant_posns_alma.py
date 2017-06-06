from taskinit import *

(tb,)=gentools(['tb'])

SRV_WSDL_URL = 'http://asa.alma.cl/axis2/services/TMCDBAntennaPadService?wsdl'
ALMA_CONFIG_NAME = 'CURRENT.AOS'

def query_tmcdb_antennas(ant_names, timestamp):
    """
    Retrieves information for a list of antennas, at a given timestamp,
    by querying the ALMA TMC DB AntennaPad web service.
    See service documentation:
    https://ictwiki.alma.cl/twiki/bin/view/SoftOps/TMCDBAntennaPadService

    :param ant_names: list of antenna names as strings
    :param timestamp: timestamp specification string (in
                     '2017-01-01T16:53:54' format)

    :returns: response as retrieved from the web service
    :raises urllib2.URLError: if network or server issues
    """

    from suds.client import Client

    # build a string of space-separated antenna names
    ant_names_str = ' '.join(ant for ant in ant_names)
    casalog.post('Querying positions for these antennas: {0}, '
                 'at time {1}'.format(ant_names_str, timestamp), 'INFO')

    srv_wsdl_url = SRV_WSDL_URL
    wsclient = Client(srv_wsdl_url)

    config_name = ALMA_CONFIG_NAME
    resp = wsclient.service.getAntennaPositions(config_name,
                                                ant_names_str, timestamp)
    return resp

def process_tmcdb_response_for_gencal(resp):
    """
    Takes a response from the TMC DB AntennaPad service and produces a
    list of antenna names and a list of their position correction
    parameters

    :param resp: response object from web service

    :returns: a tuple with the second and third elements returned by
    correct_ant_posns_alma() (and correct_ant_posns()). The tuple has
    as first element a string of comma-separated antenna names, and as
    second element a vector of parameters (positions as Bx, By, Bz for
    every antenna).
    """

    def check_pos(pos):
        """
        Basic sanity checks on a position object

        :param pos: position object as returned for one antenna by the
        TMC DB AntennaPad service.
        """
        if not pos:
            casalog.post('Empty position object received for antenna with index {}'.
                         format(ant_idx), 'ERROR')

        try:
            pos.completion
        except AttributeError, exc:
            raise RuntimeError('Got response but with only a single completion '
                               'object: {0}. errorList: {1}'.
                               format(pos, pos.errorList))
    def get_ant_params(pos):
        """
        Get the position parameters for an antenna from a position object

        :param pos: position object as returned for one antenna by the
        TMC DB AntennaPad service.

        :returns: tuple: position as a vector [x,y,z], position found
        """
        ant_params = []
        pos_found = False
        # when the antenna is not found
        if not pos.completion.status == 'true' or not pos.position:
            casalog.post('Did not find position parameters for antenna {0}.'
                         'Error description: {1}'.
                         format(pos.name, pos.completion.errorList), 'WARN')
            ant_params = [0, 0, 0]
        else:
            ant_params = [pos.position.x, pos.position.y, pos.position.z]
            pos_found = True

        return (pos_found, ant_params)

    if not resp:
        raise RuntimeError('No response received: {}'.format(resp))
    elif not resp.antennaPositionList or 0 == len(resp.antennaPositionList):
        raise RuntimeError('Response with no list: {}'.format(resp))

    casalog.post('Got this response from ALMA TMC DB: {}'.
                 format(resp), 'DEBUG')

    accum_ant_names = []
    accum_params = []
    cnt_params_found = 0
    for _ant_idx, pos in enumerate(resp.antennaPositionList):
        check_pos(pos)

        accum_ant_names.append(pos.name)
        (found, ant_params) = get_ant_params(pos)
        if found:
            cnt_params_found += 1
        accum_params.extend(ant_params)

    casalog.post('Found position parameters for {} antennas out of {} '
                 'queried in total '.
                 format(cnt_params_found, len(accum_ant_names)), 'INFO')

    # build a string of comma-separated antenna names
    ant_names_str = ','.join(str(ant) for ant in accum_ant_names)

    return (ant_names_str, accum_params)

def correct_ant_posns_alma(vis_name, print_offsets=False):
    """
    Implements the interface of correct_ant_posns with specific logic
    for ALMA.

    This produces a return error code as returned by correct_ant_posns()
    by capturing some common exceptions.
    """

    def get_ant_pad_names(vis_name):
        """
        :param vis_name: name of an MS

        :returns: a tuple: 1) list of antenna names from the MS, 2) list
        of pad/station names from the MS
        """
        tb.open(vis_name + '/ANTENNA')
        ant_names = tb.getcol('NAME')
        pad_names = tb.getcol('STATION')
        tb.close()
        return (ant_names, pad_names)

    def get_time_range_from_obs(vis_name):
        observation = tb.open(vis_name + '/OBSERVATION')
        if 'ASDM_ANTENNA' not in tb.keywordnames():
            casalog.post('The ANTENNA tables from the ASDM were not '
                         'transferred to the Measurement Set', 'WARN')
        if 'ASDM_STATION' not in tb.keywordnames():
            casalog.post('The ANTENNA and/or STATION tables from the ASDM '
                         'were not transferred to the Measurement Set',
                         'WARN')
        time_range = tb.getcol('TIME_RANGE')
        tb.close()
        return time_range

    def build_obs_time(time_range):
        """
        :returns: time specification string formatted as for example
        2011-11-13T04:10:12
        """
        import time
        obs_time = (time_range[0] + time_range[1]) / 2.0
        obs_time = ((obs_time / 86400.0) + 2400000.5 - 2440587.5) * 86400.0
        obs_time = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(obs_time))
        return obs_time

    time_range = get_time_range_from_obs(vis_name)
    ant_names, _pad_names = get_ant_pad_names(vis_name)
    obs_time = build_obs_time(time_range)

    # the three element to return
    ret_code = 1
    ant_names_str = ''
    ant_params = []
    try:
        import urllib2
        response = query_tmcdb_antennas(ant_names, obs_time)
        (ant_names_str, ant_params) = process_tmcdb_response_for_gencal(response)
        ret_code = 0
    except urllib2.URLError as exc:
        casalog.post('Network or server issue found while querying ALMA TMC DB '
                     'service. Details: {}'.format(exc), 'ERROR')
        ret_code = 2
    except RuntimeError as exc:
        casalog.post('Issue found while querying ALMA TMC DB service. Details: {}'.
                     format(exc), 'ERROR')
        ret_code = 1

    return [ret_code, ant_names_str, ant_params]
