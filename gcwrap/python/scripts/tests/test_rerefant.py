import os
import shutil
from tasks import *
from taskinit import *
from __main__ import default
import unittest
import numpy

datadir = os.environ.get('CASAPATH').split()[0] + '/data/regression/evn/'
src = datadir + 'n08c1.ms'


class rerefant_test(unittest.TestCase):

    def test_fringefit(self):
        dst = "n08c1_reref.ms"
        calfile = "fringe.cal"
        recalfile = "reref.cal"

        shutil.copytree(src, dst)

        # Perform fringe fit with EF as the reference station.
        fringefit(vis=dst, caltable=calfile, refant="EF")
        self.assertTrue(os.path.exists(calfile))

        # Rereference the results using ON as the new reference station.
        rerefant(vis=dst, tablein=calfile, caltable=recalfile, refant="ON")
        self.assertTrue(os.path.exists(recalfile))

        # Check original calibration table.
        tb.open(calfile)
        ant2=tb.getcol('ANTENNA2')
        taql="ANTENNA1 == 3 && ANTENNA2 == 0"
        tsel=tb.query(taql)
        fparam=tsel.getcol('FPARAM')
        tb.close()
        # Reference anttena is EF, aka antenna number 0.
        self.assertTrue(len(set(ant2)) == 1)
        self.assertTrue(0 in set(ant2))

        # Check rereferenced calibration table.
        tb.open(recalfile)
        ant2=tb.getcol('ANTENNA2')
        taql="ANTENNA1 == 0 && ANTENNA2 == 3"
        tsel=tb.query(taql)
        refparam=tsel.getcol('FPARAM')
        tb.close()
        # Reference anttena is ON, aka antenna number 3.
        self.assertTrue(len(set(ant2)) == 1)
        self.assertTrue(3 in set(ant2))

        # Parameters on EF-ON baseline should be opposite.
        self.assertTrue(numpy.isclose(refparam, -fparam, 1e-15).all())

        # Rereferenced parameters for ON should all be zero.
        tb.open(recalfile)
        taql="ANTENNA1 == 3 && ANTENNA2 == 3"
        tsel=tb.query(taql)
        refparam=tsel.getcol('FPARAM')
        tb.close()
        self.assertTrue(not refparam.any())

        # Rereferenced parameters for EF should not be zero.
        tb.open(recalfile)
        taql="ANTENNA1 == 0 && ANTENNA2 == 3"
        tsel=tb.query(taql)
        refparam=tsel.getcol('FPARAM')
        tb.close()
        self.assertTrue(refparam.any())

        shutil.rmtree(recalfile)
        shutil.rmtree(calfile)
        shutil.rmtree(dst)


def suite():
    return [rerefant_test]
