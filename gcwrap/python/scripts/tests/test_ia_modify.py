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
# Test suite for the CASA tool method ia.modify()
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
# Test for the ia.modify() tool method
# </etymology>
#
# <synopsis>
# Test the ia.modify() tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_modify[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.modify() tool method to ensure
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

datapath=os.environ.get('CASAPATH').split()[0]+'/data/regression/unittest/ia_modify/'

class ia_modify_test(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def test_stretch(self):
        """ ia.histogram(): Test stretch parameter"""
        mycl = cltool()
        mycl.addcomponent(flux=1, dir=['J2000', '00:00:00.00', '00.00.00.0'])
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
            yy.modify, model=mycl.torecord(),
            mask=mymask + ">0", stretch=False
        )
        zz = yy.modify(
            model=mycl.torecord(), mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz and type(zz) == type(True))
        yy.done()
        mycl.done()
        
    def test_CAS5688(self):
        """verify output is the same after this performance fix"""
        myia = iatool()
        for i in [0, 1]:
            if i == 0:
                imagename = "CAS5688_1.im"
            elif i == 1:
                imagename = "CAS5688_2.im"
            myia.fromshape(imagename, [20,20,1,10])
            if i == 0:
                world = myia.toworld([4.4, 4.4])['numeric']
            elif i == 1:
                world = myia.toworld([4.51, 4.51])['numeric']
            myia.setbrightnessunit("Jy/pixel")
            v0 = qa.quantity(world[0], "arcmin")
            v1 = qa.quantity(world[1], "arcmin")
            dir = me.direction("J2000", v0, v1)
            mycl = cltool()
            mycl.addcomponent(
                [1,0, 0, 0], "Jy", dir=dir, shape="point",
                polarization="Stokes", spectrumtype="spectral index", index=2.5
            )
            myia.modify(model=mycl.torecord(), subtract=False)
            bb = myia.getchunk()
            myia.done()
            mycl.done()
            myia.open(datapath + imagename)
            cc = myia.getchunk()
            myia.done()
            self.assertTrue((bb == cc).all())

    def test_history(self):
        """Test history is added"""
        mycl = cltool()
        mycl.addcomponent(flux=1, dir=['J2000', '00:00:00.00', '00.00.00.0'])
        myia = iatool()
        myia.fromshape("", [200, 200, 1, 1])
        self.assertTrue(
            myia.modify(model=mycl.torecord()), "Failed to run ia.modify"
        )
        msgs = myia.history()
        mycl.done()
        myia.done()
        self.assertTrue("ia.modify" in msgs[-2], "History not written")
        self.assertTrue("ia.modify" in msgs[-1], "History not written")
 
    def test_oldtests(self):
        """tests moved from imagetest_regression test_25"""
        
        def gaussian(flux, major, minor, pa, dir=None):
            newcl = cltool()
            newcl.simulate(1,log=False);
            newcl.setshape(which=0, type='Gaussian',
                majoraxis=major, minoraxis=minor, positionangle=pa,
                majoraxiserror = '0arcsec', minoraxiserror = '0arcsec',
                positionangleerror = '0deg', log=False
            )
            flux2 = [flux, 0, 0, 0];
            newcl.setflux(
                which=0, value=flux2, unit='Jy',
                polarization='Stokes', log=False
            )
            if dir==None:
                dir = me.direction('J2000', '0rad', '0rad')
            values = me.getvalue(dir);
            newcl.setrefdir(which=0, ra=values['m0'], dec=values['m1'], log=False);
            return newcl;

        def compareComponentList(cl0, cl1, tol=0.005, dotype=True):
            n0 = cl0.length()
            n1 = cl1.length()
            errmsg = 'compareComponentList: '
            if (n0 != n1):
                errmsg += 'Number of components differ'
                print(cl0.torecord())
                print(cl1.torecord())
                info(errmsg)
                return False
            #
            for i in range(0,n0):
                f0 = cl0.getfluxvalue(i)
                f1 = cl1.getfluxvalue(i)
                t = tol * f0[0]
                self.assertTrue(
                    numpy.isclose(f1, f0, t).all(), 'Component fluxes differ'
                )
                shp0 = cl0.getshape(i)
                shp1 = cl1.getshape(i)
                type0 = cl0.shapetype(i)
                type1 = cl1.shapetype(i)
                if (dotype and type0!=type1):
                    errmsg+='Component types differ'
                    info(errmsg)
                    return False
                #
                dir0 = cl0.getrefdir(i)
                dir1 = cl1.getrefdir(i)
                #
                v0 = me.getvalue(dir0)
                v1 = me.getvalue(dir1)
                #
                d = abs(qa.convert(v1['m0'],v0['m0']['unit'])['value'] - v0['m0']['value'])
                t = tol * abs(v0['m0']['value'])
                if (d > t):
                    errmsg+='Longitudes differ'
                    info(errmsg)
                    return False
                #
                d = abs(qa.convert(v1['m1'],v0['m1']['unit'])['value'] - v0['m1']['value'])
                t = tol * abs(v0['m1']['value'])
                if (d > t):
                    errmsg+='Latitudes differ'
                    info(errmsg)
                    return False
                #
                if dotype and (type0=='Gaussian' or type1=='Disk'):
                    q0 = shp0['majoraxis']
                    q1 = shp1['majoraxis']
                    d = abs(qa.convert(q1,q0['unit'])['value']  - q0['value'])
                    t = tol * q0['value']
                    if (d > t):
                        errmsg+='Major axes differ'
                        info(errmsg)
                        return False
                    #
                    q0 = shp0['minoraxis']
                    q1 = shp1['minoraxis']
                    d = abs(qa.convert(q1,q0['unit'])['value']  - q0['value'])
                    t = tol * q0['value'];
                    if (d > t):
                        errmsg+='Minor axes differ'
                        info(errmsg)
                        return False
                    #
                    q0 = shp0['positionangle']
                    q1 = shp1['positionangle']
                    d = abs(qa.convert(q1,q0['unit'])['value']  - q0['value'])
                    t = tol * q0['value']
                    if (d > t):
                        errmsg+='Position angles differ'
                        info(errmsg)
                        return False
                    #
            return True

        imname = 'ia.fromshape.image'
        imshape = [128,128,1]
        myim = ia.newimagefromshape(imname, imshape)
        self.assertTrue(myim, 'ia.fromshape constructor 1 failed')
        self.assertTrue(myim.setbrightnessunit('Jy/beam'), 'failed in setbrightnessunit')
        self.assertTrue(
            myim.setrestoringbeam(
                major='5arcmin', minor='2.5arcmin', pa='60deg', log=False
            ), 'failed in setrestoringbeam'
        )
        #
        # Pretty hard to test properly.  Add model
        #
        qmaj = '10arcmin'  
        qmin = '5arcmin'   
        qpa = '45.0deg'   
        flux = 100.0
        cl0 = gaussian(flux, qmaj, qmin, qpa) 
        self.assertTrue(
            myim.modify(cl0.torecord(), subtract=False), 'failed in modify'
        )
        stats = myim.statistics(list=False)
        self.assertTrue(stats, 'failed to get statistics')
        diff = abs(stats['flux']-flux)/flux
        self.assertTrue(
            numpy.isclose(stats['flux'], flux, 0.001),
            'model image 1 has wrong values'
        )
        # Subtract it again
        # debug
        self.assertTrue(myim.modify(cl0.torecord(), subtract=True))
        stats = myim.statistics(list=False)
        self.assertTrue(stats)
        p = myim.getchunk()
        self.assertTrue(
            numpy.all(numpy.isclose(p, 0, atol=1e-6)),
            'model image 2 has wrong values'
        )
        #
        # Now add the model for fitting
        #
        self.assertTrue(myim.modify(cl0.torecord(), subtract=False))
        #

        cl1 = myim.fitcomponents(
            region=rg.box(blc=[32, 32, 0],
            trc=[96, 96, 0])
        )
        myim.done()
        self.assertTrue(cl1, 'fitcomponents 1 failed')
        self.assertTrue(cl1['converged'], 'fitcomponents 1 did not converge')
        cl1tool=cltool()
        cl1tool.fromrecord(cl1['results'])
        self.assertTrue(
            compareComponentList(cl0,cl1tool), 'failed fitcomponents 1'
        )
        
    def test_point_source_fix(self):
        """Test fix of x/y swap bug CAS-11502"""
        mycl = cltool()
        # pixel [10, 30]
        mycl.addcomponent(flux=1, dir=['J2000', '00:02:40.00', '-00.20.00.0'])
        myia = iatool()
        myia.fromshape("", [100, 100, 1, 1])
        self.assertTrue(
            myia.modify(model=mycl.torecord()), "Failed to run ia.modify"
        )
        mycl.done()
        stats = myia.statistics()
        myia.done()
        self.assertTrue(
            (stats['minpos'] == [10, 30, 0, 0]).all(),
            "Incorrect point source pixel posiiton"
        )

    def test_disk(self):
        """test disk gives the right flux, CAS-10887"""
        mycl = cltool()
        mycl.addcomponent(
            dir="J2000 0:00:00 0.00.00", flux=1.0, shape="disk", majoraxis="100arcmin",
            minoraxis="100arcmin", positionangle="0deg"
        )
        myia = iatool()
        myia.fromshape("", [101,101])
        myia.modify(mycl.torecord(), subtract=False)
        mycl.done()
        self.assertTrue(numpy.isclose(myia.statistics()['sum'][0], 1, 1e-2))
        myia.done()

def suite():
    return [ia_modify_test]
