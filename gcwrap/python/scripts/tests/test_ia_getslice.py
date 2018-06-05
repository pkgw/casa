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
# Test suite for the CASA tool method ia.getslice()
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
# Test for the ia.getslice() tool method
# </etymology>
#
# <synopsis>
# Test the ia.getslice() tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_getslice[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.getslice() tool method to ensure
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

class ia_getslice_test(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def test_basic(self):
        """Test basic functionality"""
        shape = [10, 10]
        dv = numpy.zeros(shape, dtype=numpy.float64)
        cv = numpy.zeros(shape, dtype=numpy.complex128)
        for i in range(10):
            for k in range(10):
                x = i+k
                dv[i, k] = x
                cv[i, k] = complex(x, 2*x)
        myia = iatool()
        expdist = [
            0.        ,  0.31426975,  0.62853932,  0.9428091 ,  1.25707865,
            1.57134843,  1.88561797,  2.19988775,  2.5141573 ,  2.82842708
        ]
        expxpos = [
            1.        ,  1.22222221,  1.44444442,  1.66666663,  1.88888884,
            2.11111116,  2.33333325,  2.55555558,  2.77777767,  3.
        ]
        expypos = [
            2.        ,  2.22222233,  2.44444442,  2.66666675,  2.88888884,
            3.11111116,  3.33333325,  3.55555558,  3.77777767,  4.
        ]
        exprpix = []
        for i in range(10):
            exprpix.append(3 + 4.0*i/9.0)
        expcpix = (1 + 2j) * numpy.array(exprpix)
        for mytype in ['f', 'd', 'c', 'cd']:
            myia.fromshape("", shape, type=mytype)
            if mytype == 'f' or mytype == 'd':
                myia.putchunk(dv)
            else:
                myia.putchunk(cv)
            bb = myia.getslice(x=[1,3], y=[2,4], npts=10, method='linear')
            self.assertTrue(
                numpy.isclose(bb['distance'], expdist).all(),
                "distance is wrong"
            )
            self.assertTrue(
                numpy.isclose(bb['xpos'], expxpos).all(),
                "xpos is wrong"
            )
            self.assertTrue(
                numpy.isclose(bb['ypos'], expypos).all(),
                "ypos is wrong"
            )
            self.assertTrue(
                bb['mask'].all(), "mask is wrong"
            )
            if mytype == 'f' or mytype == 'd':
                self.assertTrue(
                    numpy.isclose(bb['pixel'], exprpix).all(),
                    "data values are wrong"
                )
            else:
                self.assertTrue(
                    numpy.isclose(bb['pixel'], expcpix).all(),
                    "data values are wrong"
                )
            myia.done()

def suite():
    return [ia_getslice_test]
