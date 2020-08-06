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
# Test suite for the CASA tool method ia.fromarray()
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
# Test for the ia.fromarray() tool method
# </etymology>
#
# <synopsis>
# Test for the ia.fromarray tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_fromarray[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.fromarray() tool method to ensure
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

class ia_fromarray_test(unittest.TestCase):
    
    def setUp(self):
        self._myia = iatool()
    
    def tearDown(self):
        self._myia.done()
        self.assertTrue(len(tb.showcache()) == 0)
        
    def test_fromarray(self):
        """Test general functionality"""
        myia = self._myia
        ar1 = numpy.zeros([2, 3], numpy.float64)
        fval = 2.2
        ar1[:] = fval
        ar2 = numpy.zeros([4, 4], numpy.complex)
        cval = 2 - 6j
        ar2[:] = cval
        i = 0
        for a in [ar1, ar2]:
            myia.fromarray("", a)
            self.assertTrue((myia.shape() == a.shape).all())
            bb = myia.getchunk()
            if (i == 0):
                self.assertTrue(abs(bb[0,0] - fval) < 1e-6)
            else:
                self.assertTrue(bb[0,0] == cval)
            i += 1
        myia.done()

    def test_history(self):
        """test writing of history"""
        myia = self._myia
        ar1 = numpy.zeros([2, 3], numpy.float64)
        myia.fromarray("", ar1)
        msgs = myia.history()
        self.assertTrue("ia.fromarray" in msgs[-2])
        self.assertTrue("ia.fromarray" in msgs[-1])
        
    def test_precision(self):
        """Test type parameter"""
        jj = 1.2345678901234567890123456789
        zz = numpy.array([jj, jj])
        myia = self._myia
        self.assertRaises(Exception, myia.fromarray, "", zz, type="x")
        ia2 = iatool()
        for i in [0, 1]:
            if i == 0:
                myia.fromarray("", zz, type="f")
            else:
                myia = ia2.newimagefromarray("", zz, type="f")
            kk = myia.getchunk()[0]
            myia.done()
            self.assertTrue(numpy.isclose(kk, jj, 1e-8, 1e-8))
            self.assertFalse(numpy.isclose(kk, jj, 1e-9, 1e-9))
            if i == 0:
                myia.fromarray("", zz, type="d")
            else:
                myia = ia2.newimagefromarray("", zz, type="d")
            kk = myia.getchunk()[0]
            myia.done()
            self.assertTrue(numpy.isclose(kk, jj, 1e-18, 1e-18))
            ia2.done()
        oc = 1 + 1j
        jj = jj * oc
        zz = numpy.array([jj, jj])
        for i in [0, 1]:
            if i == 0:
                myia.fromarray("", zz, type="f")
            else:
                myia = ia2.newimagefromarray("", zz, type="f")
            kk = myia.getchunk()[0]
            myia.done()
            print("diff", (1 - jj/kk))
            self.assertTrue(numpy.isclose(kk, jj, 1e-8*oc, 1e-8*oc))
            self.assertFalse(numpy.isclose(kk, jj, 1e-9*oc, 1e-9*oc))
            if i == 0:
                myia.fromarray("", zz, type="d")
            else:
                myia = ia2.newimagefromarray("", zz, type="d")
            kk = myia.getchunk()[0]
            myia.done()
            self.assertTrue(numpy.isclose(kk, jj, 1e-18*oc, 1e-18*oc))
            ia2.done()
        
def suite():
    return [ia_fromarray_test]
