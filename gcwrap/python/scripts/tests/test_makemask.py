import shutil
import unittest
import os
import numpy
from recipes.pixelmask2cleanmask import *
from tasks import *
from taskinit import *
from __main__ import default
"""
Unit tests for task makemask 

"""

datapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/makemask/'
#datapath = '../localtest/'

#debug=True
debug=False
_ia = iatool( )

class makemaskTestBase(unittest.TestCase):
    """ base class for makemask unit test"""
    def compareimpix(self,refimage,inimage,evalmask=False):
	"""
	do pixel by pixel comparison of the cube image
	returns true if identical in terms of pixel values
	and shape (does not check coordinates)
	"""
	_ia.open(refimage)
	refdata=_ia.getchunk(getmask=evalmask)
	refshp=list(_ia.shape())
	_ia.close()
	_ia.open(inimage)
	indata=_ia.getchunk(getmask=evalmask)
	inshp=list(_ia.shape())
	_ia.close()
	if inshp==refshp:
	    for i in range(inshp[3]):
		for j in range(inshp[2]):
		    for k in range(inshp[1]):
			for l in range(inshp[0]):
			    if not indata[l][k][j][i]==refdata[l][k][j][i]:
				return False
	return True

 
class test_copy(makemaskTestBase):
    """test copy mode"""

    #data in repository
    inimage='ngc5921.cube1.mask'
    inimage2='ngc5921.cube2.mask'
    inimage3='ngc5921.cube1.bmask'
    inimage4='ngc5921.cube1.image' # actual cube image
    inimage5='ngc5921.cube1.UNKNOWNTEL.mask'# unknown telesocpe
    inimage6='3x3.image'

    outimage1='ngc5921.cube1.copy.mask'
    outimage2='ngc5921.cube1.copyinimage.mask'
    outimage3='ngc5921.cube2.copyinimage.mask'
    outimage4='ngc5921.cube2.copy.mask'
    outimage5='3x3b.image'
   
    refimage4=datapath+'reference/ngc5921.copytest4.ref.mask'
    refimage6=datapath+'reference/ngc5921.copytest6.ref.mask'

    def setUp(self):
        #for img in [self.inimage,self.outimage1,self.outimage2, self.outimage3]:
        #    if os.path.isdir(img):
        #        shutil.rmtree(img)
        for img in [self.inimage,self.inimage2,self.inimage3, self.inimage4,
            self.inimage5, self.inimage6]:
            if not os.path.isdir(img):
                shutil.copytree(datapath+img,img)

    def tearDown(self):
        if not debug:
            for img in [self.inimage,self.inimage2,self.inimage3,self.inimage4,self.inimage5, self.outimage1,self.outimage2,self.outimage3, self.outimage4]:
                #pass
                if os.path.isdir(img):
                    shutil.rmtree(img)
        else:
            print "debugging mode: clean-up did not performed" 
        
    def test1_copyimagemask(self):
        """ (copy mode) testcopy1: copying an image mask (1/0 mask) to a new image mask"""
        try:
            makemask(mode='copy',inpimage=self.inimage,inpmask=self.inimage,output=self.outimage1)
        except Exception, e:
            print "\nError running makemask"
            raise e
        
        self.assertTrue(os.path.exists(self.outimage1))           
        self.assertTrue(self.compareimpix(self.inimage,self.outimage1))           
         
    def test2_copyimagemask(self):
        """ (copy mode) testcopy2: copying an image mask (1/0 mask) to an existing image as a T/F mask"""
        shutil.copytree(self.inimage,self.outimage1)
        # or one can write to the image defined in inpimage if output=self.inimage+":masknew"
        # overwrite=True is necessary if the internal mask, masknew exist already...
        try:
            makemask(mode='copy',inpimage=self.inimage,inpmask=self.inimage,output=self.outimage1+":masknew")
        except Exception, e:
            print "\nError running makemask"
            raise e
       
        self.assertTrue(os.path.exists(self.outimage1))           
        self.assertTrue(self.compareimpix(self.inimage,self.outimage1))           
        _ia.open(self.outimage1)
        masknames=_ia.maskhandler('get')
        if masknames.count('masknew')==1:
            pixelmask2cleanmask(self.outimage1,'masknew','_tmp_im',False)
            self.assertTrue(self.compareimpix(self.inimage,'_tmp_im'))
            shutil.rmtree('_tmp_im')
            _ia.done()

    def test3_copyimagemask(self):
        """ (copy mode) testcopy3: copying an image mask (1/0 mask) to a new image with internal (T/F) mask"""
        # input mask is a 1/0 mask and 0s are interpreted as masked region in masknew...
        try:
            makemask(mode='copy',inpimage=self.inimage,inpmask=self.inimage,output=self.outimage2+":masknew")
        except Exception, e:
            print "\nError running makemask"
            raise e
         
        self.assertTrue(os.path.exists(self.outimage2))
        # check by converting the output back to a 1/0 mask to see if it 
        # is identical to the input 1/0 mask.
        if os.path.exists(self.outimage2):
            _ia.open(self.outimage2)
            masknames=_ia.maskhandler('get')
            if masknames.count('masknew')==1:
                pixelmask2cleanmask(self.outimage2,'masknew','_tmp_im',False)
                self.assertTrue(self.compareimpix(self.inimage,'_tmp_im')) 
                shutil.rmtree('_tmp_im')
            _ia.done()

    def test4_copyimagemask(self):
        """ (copy mode) testcopy4: copying an image mask (1/0 amsk) to a new image with different coordinates(regrid)""" 
        try:
            makemask(mode='copy',inpimage=self.inimage2,inpmask=self.inimage, output=self.outimage4)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage4))
        self.assertTrue(self.compareimpix(self.outimage4,self.refimage4)) 

    def test5_copyimagemask(self):
        """ (copy mode) testcopy5: copying a region txt to a new image with a copy of inpimage and store the mask as an in-image(T/F) mask"""
        # expected behavior: create an new image with name defined in self.outimage3, the input region will be
        # translated as a good/valid region and so unmasked and contains values of inpimage of the corresponding region.
        # All outside the ellipse will be masked as defined in newmask. The pixel values of inpimage is also copied over
        # output image. In this case, inpimage is an 1/0 mask image so the resulting ouput image contains 0's and masked pixels.

        try:
            makemask(mode='copy',inpimage=self.inimage2,
                     inpmask='ellipse [[13:30:15.79110, +030.13.51.8986], [340.4877arcsec, 299.4327arcsec], 0.00000000deg]', 
                     output=self.outimage3+':newmask')
        except Exception, e:
            print "\nError running makemask"
            raise e
           
        self.assertTrue(os.path.exists(self.outimage3))
        if os.path.exists(self.outimage3):
          _ia.open(self.outimage3)
          outmask = _ia.maskhandler('get')
          self.assertTrue(outmask[0]=='newmask')
          stats=_ia.statistics()
          # the inpimage is 1/0 image mask in this case and the ellipse region does not overlap with the '1' region in inpimage
          # so the ellipse region inside will be 0 in output and the rest is masked.
          self.assertTrue(stats['max'][0]==0. and stats['min'][0]==0.)
          _ia.close()

    def test6_copyimagemask(self):
        """ (copy mode) testcopy6: copying an internal mask to a new image mask with different coordinates(regrid)""" 
        # use cube2's coord. and shape and cube1's intrenal mask
        #shutil.copytree(self.inimage2,self.outimage1)
        try:
            makemask(mode='copy',inpimage=self.inimage2,inpmask=self.inimage3+":maskoo", output=self.outimage4)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage4))
        self.assertTrue(self.compareimpix(self.outimage4,self.refimage6))
        


    def test7_copyimagemask(self):
        """ (copy mode) testcopy7: copying an internal mask of one image to a new internal mask in a new image """ 
        # Note: The pixel values of inpimage is copied to the new image
        #shutil.copytree(self.inimage,self.outimage1)
        try:
            makemask(mode='copy',inpimage=self.inimage4,inpmask=self.inimage3+":maskoo", output=self.outimage1+":newmask")
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        if os.path.exists(self.outimage1):
            self.assertTrue(self.compareimpix(self.inimage4,self.outimage1)) 
            _ia.open(self.outimage1)
            masknames=_ia.maskhandler('get')
            _ia.close()
            if masknames.count('newmask')==1:
                self.assertTrue(self.compareimpix(self.inimage3,self.outimage1,True)) 


    def test8_copyimagemask(self):
        """ (copy mode) testcopy8: copying an image mask with UNKOWN telescope to a new image mask"""
        try:
            makemask(mode='copy',inpimage=self.inimage4,inpmask=self.inimage5,output=self.outimage1)
        except Exception, e:
            print "\nError running makemask"
            raise e
        
        self.assertTrue(os.path.exists(self.outimage1))           
        self.assertTrue(self.compareimpix(self.inimage5,self.outimage1))           


    def test9_copyimagemask(self):
        """ (copy mode) testcopy9: copying an image mask with file names containing '+' """
        inpimagename="copytest+1.0-0.5.image"
        inpmaskname="copytest+1.0-0.5.mask"
        try:
            shutil.copytree(self.inimage4,inpimagename)
            shutil.copytree(self.inimage,inpmaskname)

            makemask(mode='copy',inpimage=inpimagename,inpmask=inpmaskname,output=self.outimage1)
        except Exception, e:
            print "\nError running makemask"
            raise e
             
        self.assertTrue(os.path.exists(self.outimage1))           
        self.assertTrue(self.compareimpix(inpmaskname,self.outimage1))           
        if not debug:
            if os.path.isdir(inpimagename):
                shutil.rmtree(inpimagename) 
            if os.path.isdir(inpmaskname):
                shutil.rmtree(inpmaskname) 

    def test10_copyimagemask(self):
        """ (copy mode) testcopy10: copying a small image with an internal mask
        to another small image with the same shape"""
        try:
            makemask(mode='copy',inpimage=self.inimage6,inpmask=self.inimage6+':mask0',output=self.outimage5+':mask0')
        except Exception, e:
            print "\nError running makemask"
            raise e
        
        self.assertTrue(os.path.exists(self.outimage5))           
        self.assertTrue(self.compareimpix(self.inimage5,self.outimage5))           


class test_merge(makemaskTestBase):
    """test merging of multiple masks in copy mode"""

    #data in repository
    inimage='ngc5921.cube1.mask'
    #inimage2='ngc5921.cube1.bmask'
    inimage2='ngc5921.cube1.bmask3'
    inimage3='ngc5921.cube2.mask'
    infile1='ellipse_rg.txt'

    outimage1='ngc5921.cube1.merge.mask'
    outimage2='ngc5921.cube1.merge.copyinmage.mask'

    refimage1=datapath+'reference/ngc5921.mergetest1.ref.mask'
    refimage2=datapath+'reference/ngc5921.mergetest2.ref.mask'
    refimage3=datapath+'reference/ngc5921.mergetest3.ref.mask'

    def setUp(self):
        #for img in [self.inimage,self.outimage1,self.outimage2, self.outimage3]:
        #    if os.path.isdir(img):
        #        shutil.rmtree(img)
        for img in [self.inimage,self.inimage2,self.inimage3]:
            if not os.path.isdir(img):
                shutil.copytree(datapath+img,img)

    def tearDown(self):
        if not debug:
            for img in [self.inimage,self.inimage2,self.inimage3,self.outimage1,self.outimage2, self.infile1, '_tmp_im']:
                #pass
                if os.path.isdir(img):
                    shutil.rmtree(img)
                elif os.path.isfile(img):
                    os.system('rm '+img)
        else:
            print "debugging mode: clean-up did not performed"

    def test1_mergemasks(self):
        """ (copy mode) mergetest1: merging image mask (1/0 mask) and T/F mask and  overwrite to an existing image(1/0) mask"""

        #Only overwrapped regions is valid (whenever T/F mask involved)
        try:
            shutil.copytree(self.inimage,self.outimage1)
        #    makemask(mode='copy',inpimage=self.inimage,inpmask=[self.inimage,self.inimage2+':maskoo'], output=self.outimage1, overwrite=True)
            makemask(mode='copy',inpimage=self.inimage,inpmask=[self.inimage,self.inimage2+':maskformergetest'], output=self.outimage1, overwrite=True)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        self.assertTrue(self.compareimpix(self.refimage1,self.outimage1))
        #shutil.rmtree(self.outimage1)

    def test2_mergemasks(self):
        """ (copy mode) mergetest2 :merging two image mask (1/0 mask) with different chan width and a  T/F mask to create(overwrite) an image(1/0) mask"""
        # 1/0 masks - added ('OR')  and apply ('AND') with the T/F mask
        try:
            #shutil.copytree(self.inimage,self.outimage1)
            #makemask(mode='copy',inpimage=self.inimage,inpmask=[self.inimage, self.inimage3, self.inimage2+':maskoo'], output=self.outimage1)
            makemask(mode='copy',inpimage=self.inimage,inpmask=[self.inimage, self.inimage3, self.inimage2+':maskformergetest'], output=self.outimage1)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        self.assertTrue(self.compareimpix(self.refimage2,self.outimage1))
        #shutil.rmtree(self.outimage1)

    def test3_mergemasks(self):
        """ (copy mode) mergetest3: merging multiple masks (image mask, boolean mask, regions)"""
        try:
            #To test for existing outfile
            # Note: if make a copy of outfile from inimage, comparison with the current ref image will fail ....
            #shutil.copytree(self.inimage,self.outimage1)
            if not os.path.exists(self.infile1):
                shutil.copy(datapath+self.infile1, self.infile1)
            makemask(mode='copy',inpimage=self.inimage,\
            #        inpmask=[self.inimage3, self.inimage2+':maskoo','ellipse_rg.txt','box[[130pix,135pix],[160pix,165pix]]'],\
                    inpmask=[self.inimage3, self.inimage2+':maskformergetest','ellipse_rg.txt','box[[130pix,135pix],[160pix,165pix]]'],\
                    output=self.outimage1, overwrite=True)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        self.assertTrue(self.compareimpix(self.refimage3,self.outimage1))
        #shutil.rmtree(self.outimage1)


    def test4_mergemasks(self):
        """ (copy mode) mergetest4: merging multiple masks (image mask, boolean mask, regions) to a new internal mask"""
        # T/F mask converted to 1/0 image would be the same as the test3 output
        # but the output image itself different as the input image is copied to output and T/F mask is applied on it.
        try:
            #shutil.copytree(self.inimage,self.outimage1)
            if not os.path.exists(self.infile1):
                shutil.copy(datapath+self.infile1, self.infile1)
            makemask(mode='copy',inpimage=self.inimage,\
                    #inpmask=[self.inimage3, self.inimage2+':maskoo','ellipse_rg.txt','box[[130pix,135pix],[160pix,165pix]]'],\
                    inpmask=[self.inimage3, self.inimage2+':maskformergetest','ellipse_rg.txt','box[[130pix,135pix],[160pix,165pix]]'],\
                    output=self.outimage1+":newmask")
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        if os.path.exists(self.outimage1):
          _ia.open(self.outimage1)
          masknames=_ia.maskhandler('get')
          if masknames.count('newmask')==1:
              pixelmask2cleanmask(self.outimage1,'newmask','_tmp_im',False)
              self.assertTrue(self.compareimpix(self.refimage3,'_tmp_im')) 
              shutil.rmtree('_tmp_im')
          _ia.done()
        #self.assertTrue(self.compareimpix(self.refimage3,self.outimage1))

    def test5_mergemasks(self):
        """ (copy mode) mergetest5: same as mergetest1 but uses full paths
        -merging image mask (1/0 mask) and T/F mask and  overwrite to an
        existing image(1/0) mask (verification of CAS-8865 fix)"""     

        #Only overwrapped regions is valid (whenever T/F mask involved)
        try:
            shutil.copytree(self.inimage,self.outimage1)
        #    makemask(mode='copy',inpimage=self.inimage,inpmask=[self.inimage,self.inimage2+':maskoo'], output=self.outimage1, overwrite=True)
            makemask(mode='copy',inpimage=os.path.abspath(self.inimage),inpmask=[os.path.abspath(self.inimage),os.path.abspath(self.inimage2)+':maskformergetest'],
            output=os.path.abspath(self.outimage1), overwrite=True)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage1))
        self.assertTrue(self.compareimpix(self.refimage1,self.outimage1))
        #shutil.rmtree(self.outimage1)


class test_expand(makemaskTestBase):
    """test expand mode"""

    #data in repository
    inimage='ngc5921.cont.mask'
    inimage2='ngc5921.cube1.mask'
    inimage3='ngc5921.cube2.mask'
    inimage4='ngc5921.cube1.bmask' #T/F mask

    outimage='ngc5921.cube1.expand.mask'
    outimage2='ngc5921.cube1.copyinmage.mask'
    outimage3='ngc5921.cube2.expand.mask'

    refimage1=datapath+'reference/ngc5921.expandtest1.ref.mask'
    refimage2=datapath+'reference/ngc5921.expandtest2.ref.mask'
    refimage3=datapath+'reference/ngc5921.expandtest5.ref.mask'
    refimage4=datapath+'reference/ngc5921.expandtest6.ref.mask'
    refimage5=datapath+'reference/ngc5921.expandtest2b.ref.mask'
    refimage6=datapath+'reference/ngc5921.expandtest7.ref.mask'

    def setUp(self):
        #for img in [self.inimage,self.outimage1,self.outimage2, self.outimage3]:
        #    if os.path.isdir(img):
        #        shutil.rmtree(img)
        for img in [self.inimage,self.inimage2,self.inimage3,self.inimage4]:
            if not os.path.isdir(img):
                shutil.copytree(datapath+img,img)

    def tearDown(self):
        if not debug:
            for img in [self.inimage,self.inimage2,self.inimage3,self.inimage4,
                        self.outimage,self.outimage2,self.outimage3, '_tmp_im']:
                if os.path.isdir(img):
                    shutil.rmtree(img)
                    #pass
        else:
            print "debugging mode: clean-up did not performed"

    def test1_expandmask(self):
        """ (expand mode) test1: an image mask from continuum clean to a cube mask"""
        try:
            shutil.copytree(self.inimage2,self.outimage)
            makemask(mode='expand',inpimage=self.inimage,inpmask=self.inimage, output=self.outimage, overwrite=True)
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage))
        self.assertTrue(self.compareimpix(self.refimage1,self.outimage))

    #def test2_expandmask(self):
    #    """ (expand mode) test2: an image mask from continuum clean to a cube mask with outfreqs by channel indices with existing ouput image(over-written)"""
    #    # this is probably not recommended the resultant mask is added to the existing outimage which already have mask in some channels.
    #    # this works while inpimage here is a continuum mask image and not a correct template, since the output image exist it will use
    #    # as is when inpmask is a continuum image. 
    #    try:
    #        shutil.copytree(self.inimage2,self.outimage)
    #        makemask(mode='expand',inpimage=self.inimage,inpmask=self.inimage, output=self.outimage, outfreqs=[4,5,6,7],overwrite=True)
    #        shutil.copytree(self.outimage,'test2result.im')
    #    except Exception, e:
    #        print "\nError running makemask"
    #        raise e
    #
    #    self.assertTrue(os.path.exists(self.outimage))
    #    self.assertTrue(self.compareimpix(self.refimage2,self.outimage))

    # replacing with this
    def test2_expandmask(self):
        """ (expand mode) test2: an image mask from continuum clean to a cube mask with outfreqs by channel indices"""
        try:
            makemask(mode='expand',inpimage=self.inimage2,inpmask=self.inimage, output=self.outimage, outfreqs=[4,5,6,7],overwrite=False)
            #shutil.copytree(self.outimage,'test2bresult.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage))
        self.assertTrue(self.compareimpix(self.refimage5,self.outimage))


    def test3_expandmask(self):
        """ (expand mode) test3: an image mask from continuum clean to a cube mask with outfreqs by a frequency range"""
        # will be the same range as test2 will be masked (ch4,5,6,7) and one in original output mask (ch9)
        try:
            shutil.copytree(self.inimage2,self.outimage)
            makemask(mode='expand',inpimage=self.inimage,inpmask=self.inimage, output=self.outimage, 
              outfreqs='1413.007MHz~1413.08MHz',overwrite=True)
            #shutil.copytree(self.outimage,'test3result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage))
        self.assertTrue(self.compareimpix(self.refimage2,self.outimage))


    def test4_expandmask(self):
        """ (expand mode) test4: an image mask from continuum clean to a cube mask with outfreqs by a velocity range"""
        try:
            makemask(mode='expand',inpimage=self.inimage,inpmask=self.inimage, output=self.outimage, 
              outfreqs='1561.62km/s~1546.16km/s',overwrite=True)
            #shutil.copytree(self.outimage,'test4result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage))
        self.assertTrue(self.compareimpix(self.refimage2,self.outimage))

    def test5_expandmask(self):
        """ (expand mode) test5: an image mask from a cube mask to another cube that sepecified by a template"""
        # pick up a channel with a mask in inpimagei(20ch-cube), and copy and expand the mask to 
        # another cube (chan width=2, 10ch-cube) - corresponds to ch 4,5,6,7 of outimage
        try:
            makemask(mode='expand',inpimage=self.inimage3, inpmask=self.inimage2, inpfreqs='1413.029MHz~1413.229MHz', 
              output=self.outimage3, outfreqs='1413.117MHz~1413.263MHz')
            #shutil.copytree(self.outimage3,'test5result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage3))
        self.assertTrue(self.compareimpix(self.refimage3,self.outimage3))

    def test6_expandmask(self):
        """ (expand mode) test6: T/F mask from a cube mask to another cube that sepecified by a template"""
        # mask at chan 4, 1561.62km/s lsrk 1413.01MHz (will drop the empty channel planes)
        # expand to chan 2 - 6
        try:
            makemask(mode='expand',inpimage=self.inimage3,inpmask=self.inimage4+':maskoo', inpfreqs='1561km/s~1556km/s', 
              output=self.outimage3, outfreqs='1559.04km/s~1517.82km/s')
            #shutil.copytree(self.outimage3,'test6result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage3))
        self.assertTrue(self.compareimpix(self.refimage4,self.outimage3))

    def test7_expandmask(self):
        """ (expand mode) test7: T/F mask from a cube mask to overwrite to antoher existing cube with different coord."""
        # template image = "cube1" is specified but will be ignored and use the existing outfile image.
        # mask at chan 4, 1561.62km/s lsrk 1413.01MHz (will drop the empty channel planes)
        # expand to chan 2 - 6
        # The output file inimage3 contains mask at chan 1 since only the channel range 2-6 is choosen by outfreqs
        # chan1 will be unmodified and the original mask at chan 1 will remaine in the output file.
        
        # make a copy of cube2
        shutil.copytree(self.inimage3,self.outimage3)
        try:
            makemask(mode='expand',inpimage=self.inimage2,inpmask=self.inimage4+':maskoo', inpfreqs='1561km/s~1556km/s', 
              output=self.outimage3, outfreqs='1559.04km/s~1517.82km/s',overwrite=True)
            #shutil.copytree(self.outimage3,'test6result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        #self.assertTrue(os.path.exists(self.outimage3))
        self.assertTrue(self.compareimpix(self.refimage6,self.outimage3))

    def test8_expandmask(self):
        """ (expand mode) test8: inputs are the same as test7: T/F mask, but write to an internal mask(in a new image)"""
        # mask at chan 4, 1561.62km/s lsrk 1413.01MHz (will drop the empty channel planes)
        # expand to chan 2 - 6
        #shutil.copytree(self.inimage3,self.outimage3)
        try:
            makemask(mode='expand',inpimage=self.inimage3,inpmask=self.inimage4+':maskoo', inpfreqs='1561km/s~1556km/s', 
              output=self.outimage3+":newmask", outfreqs='1559.04km/s~1517.82km/s',overwrite=True)
            #shutil.copytree(self.outimage3,'test6result.im')
        except Exception, e:
            print "\nError running makemask"
            raise e

        self.assertTrue(os.path.exists(self.outimage3))
        if os.path.exists(self.outimage3):
          _ia.open(self.outimage3)
          masknames=_ia.maskhandler('get')
          if masknames.count('newmask')==1:
              pixelmask2cleanmask(self.outimage3,'newmask','_tmp_im',False)
              self.assertTrue(self.compareimpix(self.refimage4,'_tmp_im')) 
              shutil.rmtree('_tmp_im')
          _ia.done()


class test_inmask(makemaskTestBase):
    """internal mask manupilations"""

    #data in repository
    inimage='ngc5921.cube1.bmask2' #T/F mask

    def setUp(self):
        for img in [self.inimage]:
            if not os.path.isdir(img):
                shutil.copytree(datapath+img,img)

    def tearDown(self):
        if not debug:
            for img in [self.inimage]:
                if os.path.isdir(img):
                    shutil.rmtree(img)
        else:
            print "debugging mode: clean-up did not performed"

    def test_deletemask(self):
        """ (delete mode) delete an internal mask from the image"""
        try:
            makemask(mode='delete',inpmask=self.inimage+':'+'mask2')
        except Exception, e:
            print "\nError running makemask"
            raise e

        _ia.open(self.inimage)
        mlist=_ia.maskhandler('get')
        _ia.close()
        self.assertTrue(mlist.count('mask2')==0)

    def test_setdefault(self):
        """ (setdefaultmask mode) set an internal mask as a default mask"""
        try:
            makemask(mode='setdefaultmask',inpmask=self.inimage+':'+'mask2')
        except Exception, e:
            print "\nError running makemask"
            raise e

        _ia.open(self.inimage)
        defaultmask=_ia.maskhandler('default')[0]
        _ia.close()
        self.assertTrue(defaultmask=='mask2')

def suite():
    #return [test_inmask]
    return [test_merge,test_expand,test_copy,test_inmask]
