import numpy as np

from taskinit import *

(tb_tool,) = gentools(['tb'])

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
    config_name = ALMA_CONFIG_NAME
    casalog.post('Querying positions for configuration: {0}, for these '
                 'antennas: {1}, at time: {2}'.
                 format(config_name, ant_names_str, timestamp), 'INFO')

    wscli = Client(SRV_WSDL_URL)
    resp = None
    try:
        resp = wscli.service.getAntennaPositions(configurationName=config_name,
                                                 antennaNames=ant_names_str,
                                                 timestamp=timestamp)
    except AttributeError as exc:
        import traceback
        traceback.print_exc(file=sys.stdout)
        raise RuntimeError('Got a response back from the server but it is '
                           'wrong or could not be parsed as a valid SOAP '
                           'response. This occurs for example when a wrong '
                           'timestamp is given. Response object: {0}. '
                           'Error description:'.
                           format(resp, exc))

    return resp


def process_tmcdb_response_for_gencal(resp):
    """
    Takes a response from the TMC DB AntennaPad service and produces a
    list of antenna names and a list of their position correction
    parameters.

    :param resp: response object from web service

    :returns: a tuple with three elements. The tuple has
    as first element a list of antenna names as strings, and as
    second element a vector of positions (as Bx, By, Bz) for
    every antenna. The third element holds the positions of the
    corresponding pads.
    """

    def check_pos(pos, ant_idx):
        """
        Basic sanity checks on a position object

        :param pos: position object as returned for one antenna by the
        TMCDB AntennaPad service.
        :param ant_idx: index of the antenna
        """
        if not pos:
            casalog.post('Empty position object received for antenna with '
                         'index {}'.format(ant_idx), 'ERROR')

        try:
            pos.completion
        except AttributeError as exc:
            raise RuntimeError('Got response but with only a single '
                               'completion object: {0}. errorList: {1}. '
                               'Antenna index: {2}, exception: {3}'.
                               format(pos, pos.errorList, ant_idx, exc))

    def get_ant_pad_position(pos):
        """
        Get the position parameters for an antenna from a position object

        :param pos: position object as returned for one antenna by the
        TMCDB AntennaPad service.

        :returns: tuple: position found, position as a list [x,y,z], pad
        position as a list [x, y, z]
        """

        pos_found = False
        # when the antenna is not found
        if not pos.completion.status == 'true' or not pos.position:
            casalog.post('Did not find position parameters for antenna {0}. '
                         'Error description: {1}'.
                         format(pos.name, pos.completion.errorList), 'WARN')
            ant_position = np.array([0, 0, 0], dtype=float)
            pad_position = np.array([0, 0, 0], dtype=float)
        else:
            ant_position = np.asarray([pos.position.x, pos.position.y,
                                       pos.position.z], dtype=float)
            pad_position = np.asarray([pos.pad.position.x, pos.pad.position.y,
                                       pos.pad.position.z], dtype=float)
            pos_found = True

        return (pos_found, ant_position, pad_position)

    if not resp:
        raise RuntimeError('No response received: {}'.format(resp))
    elif not resp.antennaPositionList or 0 == len(resp.antennaPositionList):
        raise RuntimeError('Response with no list: {}'.format(resp))

    casalog.post('Got this response from ALMA TMCDB: {}'.
                 format(resp), 'DEBUG')

    accum_ant_names = []
    accum_positions = []
    accum_pad_positions = []
    cnt_pos_found = 0
    for ant_idx, pos in enumerate(resp.antennaPositionList):
        check_pos(pos, ant_idx)

        accum_ant_names.append(pos.name)
        (found, ant_position, pad_position) = get_ant_pad_position(pos)
        if found:
            cnt_pos_found += 1
        accum_positions.append(ant_position)
        accum_pad_positions.append(pad_position)

    casalog.post('Found position parameters for {} antennas out of {} '
                 'queried in total '.
                 format(cnt_pos_found, len(accum_ant_names)), 'INFO')

    return (accum_ant_names, accum_positions, accum_pad_positions)


def correct_ant_posns_alma(vis_name, print_offsets=False):
    """
    Implements the interface of correct_ant_posns with specific logic
    for ALMA. The parameters (third returned variable) are calculated
    from the the difference between the pad+antenna position found from
    the web service - the pad+antenna position found in the MS (ASDM).

    This produces a return error code as returned by correct_ant_posns()
    by capturing some common exceptions.
    """

    def get_ant_pad_names_posns(vis_name):
        """
        Find positions recorded in the MS (ASDM).

        :param vis_name: name of an MS

        :returns: a tuple: 1) list of antenna names from the MS, 2) list
        of pad/station names from the MS, 3) positions [x, y, z] of the
        antennas as recorded in the MS.
        """
        tb_tool.open(vis_name + '/ANTENNA')
        ant_names = tb_tool.getcol('NAME')
        pad_names = tb_tool.getcol('STATION')
        tb_tool.close()

        try:
            tb_tool.open(vis_name)
            if 'ASDM_ANTENNA' not in tb_tool.keywordnames():
                raise RuntimeError('The ANTENNA tables from the ASDM were not '
                                   'transferred to the Measurement Set')
            if 'ASDM_STATION' not in tb_tool.keywordnames():
                raise RuntimeError('The ANTENNA and/or STATION tables from '
                                   'the ASDM were not transferred to the '
                                   'Measurement Set')
            tb_tool.close()

            tb_tool.open(vis_name + '/ASDM_ANTENNA')
            asdm_ant_pos = tb_tool.getcol('position')
            tb_tool.close()

            tb_tool.open(vis_name + '/ASDM_STATION')
            asdm_pad_pos = tb_tool.getcol('position')
            tb_tool.close()

            asdm_pad_pos = asdm_pad_pos[:, 0:asdm_ant_pos.shape[1]]
            ant_posns = asdm_pad_pos + asdm_ant_pos
        except RuntimeError as exc:
            casalog.post('Could not find pad and antenna position information '
                         'from the ASDM_ANTENNA and ASDM_STATION subtables. '
                         'Defaulting to the POSITION column of the ANTENNA '
                         'subtable. This means that the parameter calculates '
                         'are most likely inaccurate. Error description: {0}'.
                         format(exc), 'WARN')
            tb_tool.open(vis_name + '/ANTENNA')
            ant_posns = tb_tool.getcol('POSITION')
            tb_tool.close()

        ant_posns_lists = map(list, zip(ant_posns[0], ant_posns[1],
                                        ant_posns[2]))
        return (ant_names, pad_names, ant_posns_lists)

    def get_time_range_from_obs(vis_name):
        tb_tool.open(vis_name + '/OBSERVATION')
        time_range = tb_tool.getcol('TIME_RANGE')
        tb_tool.close()
        return time_range

    def build_obs_time(time_range):
        """
        Produces an observation time string, picking the middle point in time
        between start and end of observation.

        :returns: time specification string formatted as for example
        2011-11-13T04:10:12
        """
        import time
        obs_time = (time_range[0] + time_range[1]) / 2.0
        obs_time = ((obs_time / 86400.0) + 2400000.5 - 2440587.5) * 86400.0
        obs_time = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(obs_time))
        return obs_time

    def calc_ant_params_from_positions(ant_names, ant_corr_posns, pad_posns,
                                       ant_posns_ms):
        """
        Produces the antenna position correction offset parameters required by
        gencal. Uses positions from the TMCDB and the MS.

        :param ant_names: list of names of the antennas
        :param ant_corr_posns: antenna (corrected) positions from the TMCDB
        :param pad_posns: pad positions from the TMCDB
        :param ant_posns_ms: antenna positions recorded in the MS

        :returns: position correction parameters for the antennas listed in the
                  first parameter
        """

        def calc_param_diff(ant_corr_pos, pad_pos, ant_pos_ms):
            """
            Subtracts A (found from the TMCDB) - B (found in the MS), where
            A is: pad position + antenna position correction/vector
            B is: antenna position in MS
            to produce the correctin parameters expected by gencal

            :param ant_corr_pos: antenna position from the TMCDB
            :param pad_pos: pad position from the TMCDB
            :param ant_pos_ms: antenna position recorded in the MS

            :returns: position correction parameters for one antenna
                      (Bx, By, Bz)
            """
            import math

            lat = (math.asin(pad_pos[2] /
                             math.sqrt(pad_pos[0]**2
                                       + pad_pos[1]**2 + pad_pos[2]**2)))
            lon = math.atan2(pad_pos[1], pad_pos[0])

            itrf_diff = pad_pos + ant_corr_pos - ant_pos_ms
            param = np.array([0, 0, 0], dtype=float)
            param[0] = (-math.sin(lon) * itrf_diff[0]
                        - math.cos(lon) * math.sin(lat) * itrf_diff[1]
                        + math.cos(lon) * math.cos(lat) * itrf_diff[2])
            param[1] = (math.cos(lon) * itrf_diff[0]
                        - math.sin(lon) * math.sin(lat) * itrf_diff[1]
                        + math.sin(lon) * math.cos(lat) * itrf_diff[2])
            param[2] = (math.cos(lat) * itrf_diff[1]
                        + math.sin(lat) * itrf_diff[2])

            return param

        ant_params = []
        casalog.post('Antenna positions found in TMCDB: {0}, \nPositions '
                     'found in MS: {1}'.format(ant_corr_posns, ant_posns_ms),
                     'DEBUG')

        for idx, _name in enumerate(ant_names):
            if np.all(ant_corr_posns[idx] == 0) and np.all(pad_posns[idx] == 0):
                param = np.array([0, 0, 0], dtype=float)
            else:
                param = calc_param_diff(ant_corr_posns[idx], pad_posns[idx],
                                        ant_posns_ms[idx])

            ant_params.extend(param)

        # build a string of comma-separated antenna names
        ant_names_str = ','.join(str(ant) for ant in ant_names)

        return (ant_names_str, ant_params)

    time_range = get_time_range_from_obs(vis_name)
    ant_names, _pad_names, ant_posns_ms = get_ant_pad_names_posns(vis_name)
    obs_time = build_obs_time(time_range)

    # the three element to return
    ret_code = 1
    ant_names_str = ''
    # Get corrected positions by querying the TMC database via the TMCDB
    # AntennaPad service
    try:
        import urllib2
        response = query_tmcdb_antennas(ant_names, obs_time)
        (ant_names_db, ant_corr_posns, pad_posns) = (
            process_tmcdb_response_for_gencal(response))
        if (ant_names_db != ant_names).any():
            raise RuntimeError('The antenna names found in the MS (which were '
                               'used to query the TMCDB) do not match the '
                               'names returned by the database.\nFound in MS: '
                               '{0}.\nFound in TMCDB: {1}'.
                               format(ant_names_db, ant_names))
        ret_code = 0
    except urllib2.URLError as exc:
        casalog.post('Network or server issue found while querying ALMA TMCDB '
                     'AntennaPad service. Details: {}'.format(exc), 'ERROR')
        ret_code = 2
    except RuntimeError as exc:
        casalog.post('Issue found while querying ALMA TMCDB AntennaPad '
                     'service. Details: {}'.format(exc), 'ERROR')
        ret_code = 1

    ant_names_str, ant_params = (
        calc_ant_params_from_positions(ant_names, ant_corr_posns,
                                       pad_posns, ant_posns_ms))

    casalog.post('Parameter values produced for gencal: {0}'.
                 format(ant_params), 'DEBUG')

    return [ret_code, ant_names_str, ant_params]
