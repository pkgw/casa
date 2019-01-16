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
# Test suite for the CASA tool method ia.putregion()
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
# Test for the ia.putregion() tool method
# </etymology>
#
# <synopsis>
# Test for the ia.putregion tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_putregion[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.putregion() tool method to ensure
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

class ia_putregion_test(unittest.TestCase):
    
    def setUp(self):
        self._myia = iatool()
    
    def tearDown(self):
        self._myia.done()
        self.assertTrue(len(tb.showcache()) == 0)
        
    def test_history(self):
        """Verify history is written"""
        myia = self._myia
        myia.fromshape("", [20,20])
        bb = myia.getchunk()
        bb[:] = 5
        myia.putregion(bb)
        msgs = myia.history()
        myia.done()
        self.assertTrue("ia.putregion" in msgs[-2])
        self.assertTrue("ia.putregion" in msgs[-1])
        
    def test_precision(self):
        """Test images with pixel values of various precisions"""
        myia = self._myia
        jj = 1.234567890123456789
        kk = jj*(1 + 2j)
        jj = numpy.array([[jj, jj], [jj, jj]], dtype=float)
        kk = numpy.array([[kk, kk], [kk, kk]], dtype=complex)
        for mytype in ('f', 'd', 'c', 'cd'):
            myia.fromshape("", [2,2], type=mytype)
            self.assertTrue(myia)
            for v in (0, 1):
                if v == 0:
                    gg = jj
                else:
                    gg = kk
                if (
                    (v == 1 and (mytype == 'f' or mytype == 'd'))
                    or (v == 0 and (mytype == 'c' or mytype == 'cd'))
                ):
                    self.assertRaises(Exception, myia.putregion, gg)
                else:
                    myia.putregion(gg)
                    bb = myia.getchunk()
                    if mytype == 'f' or mytype == 'c':
                        self.assertTrue(
                            numpy.isclose(bb, gg, 1e-8, 1e-8).all()
                        )
                        self.assertFalse(
                            numpy.isclose(bb, gg, 1e-9, 1e-9).all()
                        )
                    else:
                        self.assertTrue((bb == gg).all())
            myia.done()
                
def suite():
    return [ia_putregion_test]
