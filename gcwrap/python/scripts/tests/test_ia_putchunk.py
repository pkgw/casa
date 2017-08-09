##########################################################################
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
# Test suite for the CASA tool method ia.putchunk()
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
# Test for the ia.putchunk() tool method
# </etymology>
#
# <synopsis>
# Test for the ia.putchunk tool method
# </synopsis>
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
#
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_putchunk[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.putchunk() tool method to ensure
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

class ia_putchunk_test(unittest.TestCase):

    def setUp(self):
        self._myia = iatool()

    def tearDown(self):
        self._myia.done()
        self.assertTrue(len(tb.showcache()) == 0)

    def test_fromshape(self):
        """Test general functionality"""
        myia = self._myia
        shape = [2,3,4]
        fval = 2.7
        cval = 8.6-5.4j

        # complex valued image
        myia.fromshape("", shape, type='c')
        bb = myia.getchunk()
        bb[:] = cval
        myia.putchunk(bb)
        self.assertTrue((abs(abs(myia.getchunk()) - abs(cval)) < 1e-6).all())
        bb[:] = fval
        myia.putchunk(bb)
        self.assertTrue((abs(abs(myia.getchunk()) - abs(fval)) < 1e-6).all())

        # float valued image
        myia.fromshape("", shape, type='f')
        cc = myia.getchunk()
        cc[:] = fval
        myia.putchunk(cc)
        self.assertTrue((abs(myia.getchunk() - fval) < 1e-6).all())
        #can't put a complex valued array in a float valued image
        self.assertRaises(Exception, myia.putchunk, bb)

    def test_history(self):
        """Verify history is written"""
        myia = self._myia
        myia.fromshape("", [20,20])
        bb = myia.getchunk()
        bb[:] = 5
        myia.putchunk(bb)
        msgs = myia.history()
        myia.done()
        self.assertTrue("ia.putchunk" in msgs[-2])
        self.assertTrue("ia.putchunk" in msgs[-1])

def suite():
    return [ia_putchunk_test]
