from taskinit import *
from correct_ant_posns_alma import correct_ant_posns_alma
from correct_ant_posns_evla import correct_ant_posns_evla

# for getting a single tool in gentools
(tb,) = gentools(['tb'])


def correct_ant_posns(vis_name, print_offsets=False):
    """
    Given an input visibility MS name (vis_name), find the antenna
    position offsets that should be applied.  This application should
    be via the gencal task, using caltype='antpos'.

    If the print_offsets parameter is True, will print out each of
    the found offsets (or indicate that none were found), otherwise
    runs silently.

    A list is returned where the first element is the returned error
    code, the second element is a string of the antennas, and the
    third element is a list of antenna Bx,By,Bz offsets. An example
    return list might look like:
    [ 0, 'ea01,ea19', [0.0184, -0.0065, 0.005, 0.0365, -0.0435, 0.0543] ]

    The second and third elements of the list returned are in the format
    expected by the calibrater tool method cb.specifycal() for the
    parameters antenna and parameter, respectively.

    Usage examples:

       CASA <1>: antenna_offsets = correct_ant_posns('test.ms')
       CASA <2>: if (antenna_offsets[0] == 0):
       CASA <3>:     gencal(vis='test.ms', caltable='cal.G', \
                     caltype='antpos', antenna=antenna_offsets[1], \
                     parameter=antenna_offsets[2])

    For specific details for the EVLA see correct_ant_posns_evla.
    For specific details for ALMA see correct_ant_posns_alma.
    """
    tb.open(vis_name+'/OBSERVATION')
    # specific code for different telescopes
    tel_name = tb.getcol('TELESCOPE_NAME')
    tb.close()
    if tel_name == 'EVLA' or tel_name == 'VLA':
        return correct_ant_posns_evla(vis_name, print_offsets)
    elif tel_name == 'ALMA':
        return correct_ant_posns_alma(vis_name, print_offsets)
    else:
        msg = 'Currently only work for EVLA or ALMA observations'
        if (print_offsets):
            print(msg)
        else:
            # send to casalogger
            casalog.post(msg, "WARN")
        return [1, '', []]
