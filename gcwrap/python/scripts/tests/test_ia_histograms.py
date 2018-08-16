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
# Test suite for the CASA tool method ia.histograms()
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
# Test for the ia.histograms() tool method
# </etymology>
#
# <synopsis>
# Test the ia.histograms() tool method
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_ia_histograms[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the ia.histograms() tool method to ensure
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

def alleqnum(x,num,tolerance=0):
    if len(x.shape)==1:
        for i in range(x.shape[0]):
            if not (abs(x[i]-num) < tolerance):
                print "x[",i,"]=", x[i]
                return False
    if len(x.shape)==2:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                if not (abs(x[i][j]-num) < tolerance):
                    print "x[",i,"][",j,"]=", x[i][j]
                    return False
    if len(x.shape)==3:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    if not (abs(x[i][j][k]-num) < tolerance):
                        print "x[",i,"][",j,"][",k,"]=", x[i][j][k]
                        return False
    if len(x.shape)==4:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    for l in range(x.shape[3]):
                        if not (abs(x[i][j][k][l]-num) < tolerance):
                            print "x[",i,"][",j,"][",k,"][",l,"]=", x[i][j][k]
                            return False
    if len(x.shape)>4:
        stop('unhandled array shape in alleq')
    return True

datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/ia_histograms/'

class ia_histograms_test(unittest.TestCase):
    
    def setUp(self):
        self._myia = iatool()
    
    def tearDown(self):
        self._myia.done()
    
    def test_general(self):
        """general tests"""
        # moved from imagetest_regression.py
        myia = self._myia
        imshape = [5,10]
        for mytype in ['f', 'd']:
            pixels = myia.makearray(0.0, imshape)
            pixels[0,0] = -100
            pixels[imshape[0]-1,imshape[1]-1] = 100
            imname = 'ia.fromarray_' + mytype + '.image'
            myim = myia.newimagefromarray(
                outfile=imname, pixels=pixels, type=mytype
            )
            self.assertTrue(myim)
            try:
                ok = myim.histograms(axes=[9,19])
            except Exception, e:
                print 'Caught expected Exception' + str(e)
                ok = False
            self.assertFalse(ok, 'Histograms unexpectedly did not fail (1)')
            
            nbins = 25
            idx = nbins/2+1
            hists = myim.histograms(nbins=nbins)
            self.assertTrue(hists, 'Histograms failed (1)')
            self.assertTrue(
                hists.has_key('values') and hists.has_key('counts'),
                'Histograms record does not have the correct fields'
            )
            self.assertTrue(
                len(hists['values'])==nbins and len(hists['counts'])==nbins
                , 'Histograms value arrays have the wrong shape (1)'
            )
            ok = hists['counts'][0]==1 and hists['counts'][nbins-1]==1
            ok = ok and (hists['counts'][idx-1]==(imshape[0]*imshape[1]-2))
            self.assertTrue(ok, 'histogram counts wrong')
            
            blc = [0,0]; trc = [4,4]
            r1 = rg.box(blc=blc, trc=trc)
            hists = myim.histograms(nbins=nbins, region=r1)
            self.assertTrue(hists, 'Histograms failed (2)')
            ok = (hists['counts'][0]==1) and (hists['counts'][nbins-1]==((trc[0]-blc[0]+1)*(trc[1]-blc[1]+1)-1))
            self.assertTrue(ok, 'Histograms values are wrong (2)')
    
            for j in range(imshape[1]):
                pixels[0,j] = -100*(j+1)
                pixels[imshape[0]-1,j] = 100*(j+1)
            ok = myim.putchunk(pixels)
            self.assertTrue(ok, 'putchunk failed (1)')
            hists = myim.histograms(nbins=nbins, axes=[0])
            self.assertTrue(hists, 'Histograms failed (3)')
            ok = list(hists['values'].shape)==[nbins,imshape[1]]
            ok = ok and list(hists['counts'].shape)==[nbins,imshape[1]]
            self.assertTrue(ok, 'Histograms value arrays have the wrong shape (2)')
            for j in range(imshape[1]):
                ok = hists['counts'][0,j]==1 and hists['counts'][nbins-1,j]==1
                ok = ok and alleqnum(hists['counts'][idx-1],(imshape[0]-2),tolerance=0.0001)
            self.assertTrue(ok, 'Histograms values are wrong (3)')
            
            hists = myim.histograms(includepix=[-5,5], nbins=25)
            self.assertTrue(hists, 'Histograms failed (4)')
            ok = hists['counts'][idx-1]==(imshape[0]*imshape[1]-(imshape[1]+imshape[1]))
            self.assertTrue(ok, 'Histograms values are wrong (4)')
            ok = ok and alleqnum(hists['counts'][0:(idx-2)],0,tolerance=0.0001)
            self.assertTrue(ok, 'Histograms values are wrong (4)')
            ok = ok and alleqnum(hists['counts'][idx:nbins],0,tolerance=0.0001)
            self.assertTrue(ok, 'Histograms values are wrong (4)')
            
            hists = myim.histograms()
            self.assertTrue(hists, 'histograms failed (4)')
            hists = myim.histograms()
            self.assertTrue(hists, 'histograms failed (5)')
            hists = myim.histograms(cumu=True, log=True)
            self.assertTrue(hists, 'histograms failed (6)')
    
            ok = myim.done()
            self.assertTrue(ok, 'Done failed (1)')
        
    def test_stretch(self):
        """ ia.histogram(): Test stretch parameter"""
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
            yy.histograms, mask=mymask + ">0", stretch=False
        )
        zz = yy.histograms(
            mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz and type(zz) == type({}))
        yy.done()
    
    def test_stats(self):
        """ test returned stats"""
        myia = iatool()
        myia.open(datapath + "myim.im")
        res = myia.histograms()
        self.assertTrue(abs(res['mean'][0] - -0.01586026) < 1e-7)
        self.assertTrue(abs(res['sigma'][0] - 0.99244785) < 1e-7)
        res = myia.histograms(axes=[0,1])
        bb = myia.getchunk()
        for i in range(20):
            self.assertTrue(abs(res['mean'][i] - numpy.mean(bb[:,:,i])) < 1e-7)
            self.assertTrue(abs(res['sigma'][i] - numpy.std(bb[:,:,i], ddof=1)) < 1e-7)
        myia.done()
    
def suite():
    return [ia_histograms_test]
