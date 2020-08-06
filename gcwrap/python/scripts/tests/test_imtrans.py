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
# Test suite for the CASA task imtrans
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
#   <li> <linkto class="task_imtrans.py:description">imtrans</linkto> 
# </ul>
# </prerequisite>
#
# <etymology>
# Test for the imtrans task
# </etymology>
#
# <synopsis>
# Test the imtrans task and the ia.reorder() method upon which it is built.
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_imtrans[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the imtrans task to ensure
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

good_image = "reorder_in.fits"
cas_2364im = "CAS-2364.im"

def run_transpose(imagename, outfile, order):
    myia = iatool()
    myia.open(imagename)
    print("*** order " + str(order))
    res = myia.transpose(outfile=outfile, order=order)
    myia.done()
    return res

def run_imtrans(imagename, outfile, order):
    return imtrans(imagename=imagename, outfile=outfile, order=order)


class imtrans_test(unittest.TestCase):
    
    def setUp(self):
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/imtrans/'
        shutil.copy(datapath + good_image, good_image)
    
    def tearDown(self):
        os.remove(good_image)
        self.assertTrue(len(tb.showcache()) == 0)

    def test_exceptions(self):
        """imtrans: Test various exception cases"""
        
        def testit(imagename, outfile, order):
            for i in [0,1]:
                if (i==0):
                    self.assertRaises(Exception, run_transpose, imagename, outfile, order)
                else:
                    self.assertFalse(run_imtrans(imagename, outfile, order))

        # blank imagename
        testit("", "blah", "012")
        
        # not enough specified axes
        testit(good_image, "blah", "01")
        testit(good_image, "blah", 10)
        testit(good_image, "blah", ["r", "d"])
        
        # too many specified axes
        testit(good_image, "blah", "0123")
        testit(good_image, "blah", 1230)
        testit(good_image, "blah", ["r", "d", "f", "s"])

        # Bogus axes specification
        testit(good_image, "blah", "123")
        testit(good_image, "blah", ["r", "d", "s"])
        testit(good_image, "blah", ["r", "d", "r"])
        testit(good_image, "blah", 103)
        
    def test_straight_copy(self):
        """No actual transposing"""
        imagename = good_image
        myia = iatool()
        myia.open(imagename)
        expecteddata = myia.getchunk()
        expectednames = myia.coordsys().names()
        myia.close()
        count = 0
        for order in ["012", 12, ['r', 'd', 'f'], ["righ", "declin", "freq"]]:
            for code in [run_transpose, run_imtrans]:
                outfile = "straight_copy_" + str(count)
                newim = code(imagename, outfile, order)
                if (type(newim) == bool):
                    self.assertTrue(newim)
                    newim = iatool()
                    newim.open(outfile)
                gotdata = newim.getchunk()
                gotnames = newim.coordsys().names()
                newim.done()
                self.assertTrue((expecteddata == gotdata).all())
                self.assertTrue(expectednames == gotnames)
                count += 1

    def test_transpose(self):
        """Test transposing"""
        imagename = good_image
        myia = iatool()
        myia.open(imagename)
        expecteddata = myia.getchunk()
        expectednames = myia.coordsys().names()
        myia.done()
        count = 0
        for order in ["120", 120, ['d', 'f', 'r'], ["declin", "freq", "righ"]]:
            for code in [run_transpose, run_imtrans]:
                for outname in ["transpose_" + str(count), ""]:
                    newim = code(imagename, outname, order)
                    if code == run_imtrans and len(outname) == 0:
                        self.assertFalse(newim)
                    else:
                        if (type(newim) == bool):
                            self.assertTrue(newim)
                            newim = iatool()
                            newim.open(outname)
                        gotdata = newim.getchunk()
                        inshape = expecteddata.shape
                        for i in range(inshape[0]):
                            for j in range(inshape[1]):
                                for k in range(inshape[2]):
                                    self.assertTrue(expecteddata[i][j][k] == gotdata[j][k][i])
                        gotnames = newim.coordsys().names()
                        newim.done()
                        self.assertTrue(expectednames[0] == gotnames[2])
                        self.assertTrue(expectednames[1] == gotnames[0])
                        self.assertTrue(expectednames[2] == gotnames[1])
                    count += 1

    def test_cas_2364(self):
        "test CAS-2364 fix"
        datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/imtrans/'
        shutil.copytree(datapath + cas_2364im, cas_2364im)
        order="0132"
        out1 = "blahxx.im"
        myia = iatool()
        myia.open(cas_2364im)
        trans = myia.transpose(out1, order)
        myia.done()
        trans.done()
        self.assertTrue(len(tb.showcache()) == 0)

        # to verify fix, just open the image. bug was that exception was thrown when opening output from reorder
        myia.open(out1)
        self.assertTrue(myia)
        myia.done()
        out1 = "blah2.im"
        self.assertTrue(imtrans(imagename=cas_2364im, outfile=out1, order=order))
        myia.open(out1)
        self.assertTrue(myia)
        myia.done()
        shutil.rmtree(cas_2364im)

    def test_history(self):
        """Test history records are written"""
        myia = iatool()
        imagename = "zz.im"
        myia.fromshape(imagename, [10,10,4,10])
        order = "3210"
        kk = myia.transpose("",order=order)
        myia.done()
        msgs = kk.history()
        kk.done()
        teststr = "ia.transpose"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")    
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        
        outfile = "zz_out.im"
        self.assertTrue(
            imtrans(imagename=imagename, outfile=outfile, order=order)
        )
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "imtrans"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")

    def test_imageinfo(self):
        """Verify image info is copied"""
        myia = iatool()
        myia.fromshape("",[10,20,30])
        myia.setbrightnessunit("Jy/beam")
        myia.setrestoringbeam("4arcmin", "3arcmin", "0deg")
        kk = myia.transpose("", "201")
        self.assertEqual(
            kk.brightnessunit(), "Jy/beam",
            "brightness unit not copied"
        )
        self.assertTrue(
            kk.restoringbeam(), "restoring beam not copied"
        )

def suite():
    return [imtrans_test]
