import os
import sys
import casac
import unittest

import numpy


class calanalysis_tsys1_base(unittest.TestCase):

    """A Tsys calibration table is used in these tests."""

    ca = casac.casac.calanalysis()

    calName = 'uid___A002_X30a93d_X43e.ms.tsys.s3.tbl'
    msName = 'uid___A002_X30a93d_X43e.ms'
    parType = 'Float'
    polBasis = 'U'
    visCal = 'B TSYS'

    fieldName = ['J2253+161; 3c454.3', 'Callisto',
                 'B0007+106; J0010+109', 'GRB021004']
    fieldNumber = ['0', '1', '2', '3']
    numField = len(fieldNumber)

    antennaName = ['DA41', 'DA42', 'DA43', 'DV02', 'DV03', 'DV05',
                   'DV07', 'DV10', 'DV11', 'DV12', 'DV13', 'DV14', 'PM02', 'PM03']
    antennaNumber = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                     '10', '11', '12', '13']
    numAntenna = len(antennaNumber)

    antenna1Name = antennaName
    antenna1Number = antennaNumber
    numAntenna1 = numAntenna

    antenna2Name = ['NoName']
    antenna2Number = ['-1']
    numAntenna2 = 1

    feed = ['1', '2']
    numFeed = len(feed)

    spwName = ['', '', '', '', '', '', '', '', '', '', '', '', '', '',
               '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    spwNumber = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                 '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                 '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']
    numSPW = len(spwNumber)

    numChannel = numpy.array([4, 128, 1, 128, 1, 128, 1, 128, 1, 128, 1,
                              128, 1, 128, 1, 128, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                              4, 4, 4, 4, 4], dtype=numpy.int32)

    time = numpy.array([4827167647.3920002, 4827167780.5439997,
                        4827168227.2320004, 4827168353.7600002, 4827168887.7600002,
                        4827169009.632, 4827169543.6800003, 4827169665.6960001,
                        4827170204.3520002, 4827170332.4160004, 4827170866.5600004,
                        4827170989.5360003, 4827171523.776, 4827171647.2799997])
    numTime = len(time)

    def setUp(self):
        datapath = os.environ.get('CASAPATH').split()[0]
        datapath += '/data/regression/unittest/calanalysis/'

        os.system('cp -RL {0} {1}'.format(os.path.join(datapath, self.calName),
                                          self.calName))

        return self.ca.open(self.calName)

    def tearDown(self):
            os.system('rm -rf {0}'.format(self.calName))
            return self.ca.close()


class calanalysis_tsys1_introspective(calanalysis_tsys1_base):
    """ This is a very simple unit test for introspective methods. """

    def test_introspective(self):

        """Test of introspective member functions"""

        self.assertEqual(os.path.split(self.ca.calname())[1],
                         self.calName)
        self.assertEqual(os.path.split(self.ca.msname())[1],
                         self.msName)
        self.assertEqual(self.ca.partype(), self.parType)
        self.assertEqual(self.ca.polbasis(), self.polBasis)
        self.assertEqual(self.ca.viscal(), self.visCal)

        self.assertEqual(self.ca.numfield(), self.numField)
        self.assertEqual(self.ca.field(name=True), self.fieldName)
        self.assertEqual(self.ca.field(name=False), self.fieldNumber)

        self.assertEqual(self.ca.numantenna(), self.numAntenna)
        self.assertEqual(self.ca.antenna(name=True), self.antennaName)
        self.assertEqual(self.ca.antenna(name=False),
                         self.antennaNumber)

        self.assertEqual(self.ca.numantenna1(), self.numAntenna1)
        self.assertEqual(self.ca.antenna1(name=True),
                         self.antenna1Name)
        self.assertEqual(self.ca.antenna1(name=False),
                         self.antenna1Number)

        self.assertEqual(self.ca.numantenna2(), self.numAntenna2)
        self.assertEqual(self.ca.antenna2(name=True),
                         self.antenna2Name)
        self.assertEqual(self.ca.antenna2(name=False),
                         self.antenna2Number)

        self.assertEqual(self.ca.numfeed(), self.numFeed)
        self.assertEqual(self.ca.feed(), self.feed)

        self.assertEqual(self.ca.numspw(), self.numSPW)
        self.assertEqual(self.ca.spw(name=True), self.spwName)
        self.assertEqual(self.ca.spw(name=False), self.spwNumber)

        self.assertTrue(numpy.array_equal(self.ca.numchannel(), self.numChannel))
        self.assertEqual(self.ca.numtime(), self.numTime)
        self.assertTrue(numpy.allclose(self.ca.time(), self.time))


class calanalysis_tsys1_get(calanalysis_tsys1_base):

    def _check_ca_get_out(self, out):
        """ Checks one output item from calanalysis.get(), making sure expected entries are
            there, their types are correct, and some simple values specific to the next
            test cases. """

        for entry in ['valueErr', 'value', 'flag', 'frequency']:
            self.assertTrue(entry in out)
            self.assertEqual(type(out[entry]), numpy.ndarray)

        for entry in ['feed', 'rap', 'antenna1', 'antenna2', 'field', 'abscissa']:
            self.assertTrue(entry in out)
            self.assertEqual(type(out[entry]), str)

        self.assertEqual(out['rap'], 'REAL')

    def _check_stats_items_values(self, stats):
        """ Sanity checks on the stats output from a calanalysis.get """
        stats_len = 392
        self.assertEqual(type(stats), dict)
        self.assertEqual(len(stats), stats_len)
        for idx in range(0, len(stats)):
            self.assertEqual(len(stats[str(idx)]), 11)
            self._check_ca_get_out(stats[str(idx)])
        self.assertEqual(stats['0']['feed'], '1')
        self.assertEqual(stats['0']['field'], '0')
        self.assertEqual(stats['1']['feed'], '2')
        self.assertEqual(stats[str(stats_len-1)]['feed'], '2')

    def test_get_empty(self):
        """ Test tool get function with wrong selections """

        # This uses parameters in similar wasy as the pipeline does in tsyscal/renderer
        # SPW 10 is missing
        stats10 = self.ca.get(spw='10', antenna=self.antennaName, axis='TIME',
                              ap='AMPLITUDE')
        self.assertEqual(stats10, {})

        # SPW 12 also missing
        stats12 = self.ca.get(spw='12', antenna=self.antennaName, axis='TIME',
                              ap='AMPLITUDE')
        self.assertEqual(stats12, {})

    def test_get_one_spw(self):
        """ Test tool get function. Uses the main stuff in CalAnalysys/CalStats::stats """

        # SPW 13 should be there
        # This uses parameters in similar wasy as the pipeline does in tsyscal/renderer
        stats13 = self.ca.get(spw='13', antenna=self.antennaName, axis='TIME',
                              ap='AMPLITUDE')
        self._check_stats_items_values(stats13)

    def test_get_noparams(self):
        """ Test tool get function, no selection, no other params.
             Uses stuff in CalAnalysys/CalStats::stats """

        stats_all = self.ca.get()
        self._check_stats_items_values(stats_all)


class calanalysis_tsys1_fit(calanalysis_tsys1_base):
    """ Tests on the calanalysis.fit function. """

    def _check_ca_fit(self, fit):

        for entry in ['vars', 'frequency', 'res', 'valueErr', 'flag', 'covars',
                      'pars', 'value', 'model']:
            self.assertTrue(entry in fit['1'])
            self.assertEqual(type(fit['1'][entry]), numpy.ndarray)

        for entry in ['feed', 'rap', 'antenna1', 'antenna2', 'weight', 'field',
                      'abscissa', 'order']:
            self.assertTrue(entry in fit['1'])
            self.assertEqual(type(fit['1'][entry]), str)

        for entry in ['resMean', 'redChi2', 'time', 'resVar']:
            self.assertTrue(entry in fit['1'])
            self.assertEqual(type(fit['1'][entry]), float)

        self.assertTrue(fit['1']['abscissa'], 'frequency')
        self.assertTrue(fit['1']['order'], 'LINEAR')
        self.assertTrue(fit['1']['validFit'])
        self.assertTrue(numpy.all(fit['1']['flag'] == False))
        self.assertGreater(fit['1']['resVar'], 0)
        self.assertLess(fit['1']['resVar'], 200)

    def test_fit_amp(self):
        """ Test tool fit function (amp). Exercises stuff in CalAnalysys/CalStatsFitter """

        # An amp fit inspired by pipeline/qa/bpcal.py
        fit_amp = self.ca.fit(spw='13', axis='TIME', ap='AMPLITUDE', norm=True,
                              order='LINEAR', type='LSQ', weight=False)
        fit_len = 392
        self.assertEqual(len(fit_amp), fit_len)
        self._check_ca_fit(fit_amp)

    def test_fit_phase(self):
        """ Test tool fit function (phase). Exercises stuff in CalAnalysys/CalStatsFitter """

        # A phase fit inspired by pipeline/qa/bpcal.py
        fit_phase = self.ca.fit(spw='13', axis='TIME', ap='PHASE', unwrap=True,
                                jumpmax=0.1, order='LINEAR', type='LSQ', weight=False)
        fit_len = 392
        self.assertEqual(len(fit_phase), fit_len)
        self._check_ca_fit(fit_phase)

    def test_fit_amp_sel(self):
        """ Test tool fit function (amp + selection).
            Exercises stuff in CalAnalysys/CalStatsFitter """

        # A fit with additional field selection / less outputs
        fit_amp_field = self.ca.fit(field='Callisto', spw='13', axis='TIME', ap='AMPLITUDE',
                                    norm=True, order='LINEAR', type='LSQ', weight=False)
        fit_len_field = 28
        self.assertEqual(len(fit_amp_field), fit_len_field)
        self._check_ca_fit(fit_amp_field)


def suite():
    return [calanalysis_tsys1_introspective, calanalysis_tsys1_get, calanalysis_tsys1_fit]
