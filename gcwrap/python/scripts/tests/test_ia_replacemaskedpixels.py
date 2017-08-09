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
# Test suite for the CASA tool method ia.replacemaskedpixels()
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
# Test for the ia.replacemaskedpixels() tool method
# </etymology>
#
# <synopsis>
# Test the ia.replacemaskedpixels() tool method
# </synopsis>
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
#
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_replacemaskedpixels[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.replacemaksedpixels() tool method to ensure
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

class ia_replacemaskedpixels_test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_stretch(self):
        """ ia.replacemaskedpixels(): Test stretch parameter"""
        yy = iatool()
        mymask = "maskim"
        yy.fromshape(mymask, [200, 200, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [200,200,1,20]
        yy.fromshape("", shape)
        yy.addnoise()
        self.assertRaises(
            Exception,
            yy.replacemaskedpixels, pixels=-255,
            mask=mymask + ">0", stretch=False
        )
        zz = yy.replacemaskedpixels(
            pixels=-255,
            mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz)
        yy.done()

    def test_history(self):
        """Verify history writing"""
        yy = iatool()
        mymask = "history.im"
        yy.fromshape(mymask, [20, 20])
        yy.addnoise()
        yy.replacemaskedpixels(pixels=-255, mask=mymask + ">0")
        msgs = yy.history()
        yy.done()
        self.assertTrue("ia.replacemaskedpixels" in msgs[-2])
        self.assertTrue("ia.replacemaskedpixels" in msgs[-1])

def suite():
    return [ia_replacemaskedpixels_test]
