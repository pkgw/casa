import shutil
import unittest
import os
import numpy
import math
import sys
import exceptions
import filecmp
import glob
from tasks import nrobeamaverage
from taskinit import mstool, tbtool
from __main__ import default
import testhelper as th
from sdutil import tbmanager, toolmanager, table_selector

# Define the root for the data files
datapath = os.environ.get('CASAPATH').split()[0] + "/data/regression/unittest/nrobeamaverage/"

def check_eq(val, expval, tol=None):
    """Checks that val matches expval within tol."""
    if type(val) == dict:
        for k in val:
            check_eq(val[k], expval[k], tol)
    else:
        try:
            if tol and hasattr(val, '__rsub__'):
                are_eq = abs(val - expval) < tol #absolute check
                if not are_eq:
                    are_eq = abs(val - expval)/abs(expval) < tol #relative check
            else:
                are_eq = val == expval
            if hasattr(are_eq, 'all'):
                are_eq = are_eq.all()
            if not are_eq:
                raise ValueError('!=')
        except ValueError:
            errmsg = "%r != %r" % (val, expval)
            if (len(errmsg) > 66): # 66 = 78 - len('ValueError: ')
                errmsg = "\n%r\n!=\n%r" % (val, expval)
            raise ValueError(errmsg)
        except Exception as e:
            print("Error comparing", val, "to", expval)
            raise e

class test_nrobeamaverage(unittest.TestCase):
    def setUp(self):
        default(nrobeamaverage)

        self.i_ms = "onon.ms"
        os.system('cp -RL '+ datapath + self.i_ms +' '+ self.i_ms)
        self.o_ms = "bave.ms"
        self.args = {'infile': self.i_ms, 'outfile': self.o_ms}

        self.antid = self._get_antid()
        self.min_antid = 0
        self.st_onsrc = self._get_onsource_stateid()
        self.interval = 21 # in seconds
        self.tol = 1e-5

    def tearDown(self):
        os.system('rm -rf ' + self.i_ms)
        os.system('rm -rf ' + self.o_ms)

    def _get_antid(self):
        with tbmanager(self.i_ms + '/ANTENNA') as tb:
            acol = tb.getcol('NAME')
        return range(len(acol))

    def _get_onsource_stateid(self):
        with tbmanager(self.i_ms + '/STATE') as tb:
            ocol = tb.getcol('OBS_MODE')
        res = None
        for i in range(len(ocol)):
            if ocol[i] == 'OBSERVE_TARGET#ON_SOURCE':
                res = i
                break
        if res is None: raise Exception('State ID for on_source data not found.')
        return res

    def run_task(self, aux_args=None):
        if aux_args is not None:
            for k in aux_args: self.args[k] = aux_args[k]
        nrobeamaverage(**self.args)
        self._get_data()

    def _get_data(self):
        self.i_tm, self.i_a1, self.i_a2, self.i_dd, self.i_sc, self.i_st = self._do_get_data(self.i_ms)
        self.o_tm, self.o_a1, self.o_a2, self.o_dd, self.o_sc, self.o_st = self._do_get_data(self.o_ms)

    def _do_get_data(self, msname):
        with tbmanager(msname) as tb:
            tm = tb.getcol('TIME')
            a1 = tb.getcol('ANTENNA1')
            a2 = tb.getcol('ANTENNA2')
            dd = tb.getcol('DATA_DESC_ID')
            sc = tb.getcol('SCAN_NUMBER')
            st = tb.getcol('STATE_ID')
        return tm, a1, a2, dd, sc, st

    def get_timebin(self, num_average):
        return str(num_average * self.interval) + 's'

    def check_num_data(self, num_ave=1):
        num_i_onsrc, num_i_others, num_o_onsrc, num_o_others = self._get_num_data()
        check_eq(num_i_onsrc, num_o_onsrc * num_ave)
        check_eq(num_i_others, num_o_others)

    def _get_num_data(self, stcol=None):
        if stcol is None:
            i_onsrc, i_others = self._get_num_data(self.i_st)
            o_onsrc, o_others = self._get_num_data(self.o_st)
            return i_onsrc, i_others, o_onsrc, o_others
        else:
            num_onsource = 0
            num_others = 0
            for i in range(len(stcol)):
                if (stcol[i] == self.st_onsrc):
                    num_onsource += 1
                else:
                    num_others += 1
            return num_onsource, num_others

    def check_values(self, num_ave=None, beam=None):
        if num_ave is None:
            for iidx in range(len(self.i_tm)):
                with tbmanager(self.i_ms) as tb:
                    self.i_dat = tb.getcell('FLOAT_DATA', iidx)
                oidx = self._get_index_outdata(iidx)
                with tbmanager(self.o_ms) as tb:
                    self.o_dat = tb.getcell('FLOAT_DATA', oidx)
                self._do_check_values(iidx, oidx, beam)
            return

        self.assertTrue(num_ave == 2)

        ival, oval = self._get_first_values(state=self.st_onsrc, spw=0)

        # time
        ref_tm = (ival['tm1'] + ival['tm2']) / float(num_ave)
        check_eq(oval['tm'], ref_tm, self.tol)

        # antenna ID
        ref_an = self.min_antid
        check_eq(oval['a1'], ref_an)
        check_eq(oval['a2'], ref_an)

        # spectrum
        for i in range(len(ival['dat'][0])):
            for j in range(len(ival['dat'][0][i])):
                ref_dat = (ival['dat'][0][i][j] + ival['dat'][1][i][j]) / float(num_ave)
                check_eq(oval['dat'][i][j], ref_dat, self.tol)

        # weight and sigma
        ref_wgt = float(num_ave)
        ref_sig = 1.0/math.sqrt(ref_wgt)
        for i in range(len(oval['wgt'])):
            check_eq(oval['wgt'][i], ref_wgt, self.tol)
            check_eq(oval['sig'][i], ref_sig, self.tol)

    def _get_index_outdata(self, iidx):
        res = None
        for oidx in range(len(self.o_tm)):
            if (self.o_dd[oidx] == self.i_dd[iidx]) and (self.o_sc[oidx] == self.i_sc[iidx]) and (self.o_st[oidx] == self.i_st[iidx]):
                if (self.o_st[oidx] != self.st_onsrc) and (self.o_a1[oidx] != self.i_a1[iidx]): continue
                res = oidx
                break
        if res is None: raise Exception('Output data not found.')
        return res

    def _do_check_values(self, iidx, oidx, beam=None):
        # spectrum shape
        o_npol = len(self.o_dat)
        check_eq(o_npol, len(self.i_dat))
        o_nchn = len(self.o_dat[0])
        check_eq(o_nchn, len(self.i_dat[0]))
        # spectrum value
        for ipol in range(o_npol):
            for ichan in range(o_nchn):
                check_eq(self.o_dat[ipol][ichan], self.i_dat[ipol][ichan])
        check_eq(self.o_tm[oidx], self.i_tm[iidx])
        # antenna ID
        if beam is None:
            lst_beam = self.antid
        else:
            lst_beam = beam.strip().split(',')
            for i in range(len(lst_beam)): lst_beam[i] = int(lst_beam[i])
            min_beam = lst_beam[0]
            for i in range(len(lst_beam)):
                if lst_beam[i] < min_beam: min_beam = lst_beam[i]
            self.min_antid = min_beam

        if (self.o_st[oidx] == self.st_onsrc) and (self.o_a1[oidx] in lst_beam):
            check_eq(self.o_a1[oidx], self.min_antid)
            check_eq(self.o_a2[oidx], self.min_antid)

    def _get_first_values(self, state, spw):
        ival = {}
        oval = {}
        # input data (first two data to be averaged into the first output data)
        ival['tm1'], ival['tm2'] = self._get_first_two_timestamps(spw)
        in_dat = []
        for i in range(len(self.i_tm)):
            if (self.i_st[i] == state) and (self.i_dd[i] == spw):
                if (self.i_tm[i] == ival['tm1']) or (self.i_tm[i] == ival['tm2']):
                    with tbmanager(self.i_ms) as tb:
                        in_dat.append(tb.getcell('FLOAT_DATA', i))
        ival['dat'] = in_dat

        # output data (only the first one)
        for i in range(len(self.o_tm)):
            if (self.o_st[i] == state) and (self.o_dd[i] == spw):
                oval['tm'] = self.o_tm[i]
                oval['a1'] = self.o_a1[i]
                oval['a2'] = self.o_a2[i]
                with tbmanager(self.o_ms) as tb:
                    oval['dat'] = tb.getcell('FLOAT_DATA', i)
                    oval['wgt'] = tb.getcell('WEIGHT', i)
                    oval['sig'] = tb.getcell('SIGMA', i)
                break

        return ival, oval

    def _get_first_two_timestamps(self, data_desc_id):
        time1 = None
        time2 = None
        for i in range(len(self.i_tm)):
            if (self.i_st[i] == self.st_onsrc) and (self.i_dd[i] == data_desc_id):
                if time1 is None:
                    time1 = self.i_tm[i]
                else:
                    if time2 is None:
                        if (self.i_tm[i] < time1):
                            time2 = time1
                            time1 = self.i_tm[i]
                        else:
                            time2 = self.i_tm[i]
                    else:
                        if (self.i_tm[i] < time1):
                            time2 = time1
                            time1 = self.i_tm[i]
                        elif (self.i_tm[i] < time2):
                            time2 = self.i_tm[i]
        return time1, time2

    def test_default(self): # no time averaging(timebin='0s'), rewriting beam IDs only
        self.run_task()
        self.check_num_data()
        self.check_values()

    def test_beam01(self): # beam='0,1': same as the default case
        beam = '0,1'
        self.run_task({'beam': beam})
        self.check_num_data()
        self.check_values(beam=beam)

    def test_beam0(self): # beam='0': no time averaging, no rewriting beam IDs
        beam = '0'
        self.run_task({'beam': beam})
        self.check_num_data()
        self.check_values(beam=beam)

    def test_beam1(self): # beam='1': no time averaging, no rewriting beam IDs
        beam = '1'
        self.run_task({'beam': beam})
        self.check_num_data()
        self.check_values(beam=beam)

    def test_time_averaging(self): # every two on-spectra are averaged into one specrum
        num_ave = 2
        self.run_task({'timebin': self.get_timebin(num_ave)})
        self.check_num_data(num_ave)
        self.check_values(num_ave=num_ave) # for the first data with state=on-source, spw=0


def suite():
    return [test_nrobeamaverage]
