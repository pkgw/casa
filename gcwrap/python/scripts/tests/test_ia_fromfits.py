##########################################################################
# imfit_test.py
#
# Copyright (C) 2008, 2009
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
# <author>
# Dave Mehringer
# </author>
#
# <summary>
# Test suite for the CASA tool method ia.fromfits()
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
# </ul>
# </prerequisite>
#
# <etymology>
# Test for the ia.fromfits() tool method
# </etymology>
#
# <synopsis>
# Test for the ia.fromfits() tool method
# </synopsis>
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
#
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_fromfits[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.fromfits() tool method to ensure
# coding changes do not break the associated bits
# </motivation>
#

###########################################################################
import shutil
import casac
from tasks import *
from taskinit import *
from __main__ import *
import unittest
import numpy

class ia_fromfits_test(unittest.TestCase):

    def setUp(self):
        self._myia = iatool()

    def tearDown(self):
        self._myia.done()
        self.assertTrue(len(tb.showcache()) == 0)

    def test_history(self):
        """test writing of history"""
        myia = self._myia
        ar1 = numpy.zeros([2, 3], numpy.float64)
        myia.fromarray("", ar1)
        fits = "myim.fits"
        myia.tofits(fits)
        myia.done()
        myia.fromfits("", fits)
        msgs = myia.history()
        self.assertTrue("ia.fromfits" in msgs[-2])
        self.assertTrue("ia.fromfits" in msgs[-1])

def suite():
    return [ia_fromfits_test]
