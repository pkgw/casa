#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for importfits/exportfits                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The proper handling of the FITS standard needs to be verified          #
#    using standard input files constructed based on the FITS               #
#    standard version 3.0 by Mark Calabretta                                #
#                                                                           # 
# Features tested:                                                          #
#    1) Is the import performed without raising exceptions?                 #
#    2) Is the imported image properly handling the WCS?                    #
#       (i.e. is the RA and DEC of the max. flux properly determined?)      #
#    3) Is the export performed without raising exceptions?                 #
#    4) Is the exported file compatible with the original?                  #
#    5) Does the bitpix=16 feature work                                     #
#    6) Is spectral WCS correctly read and written                          #
#    7) Does the defaultaxes feature work                                   #
#                                                                           #
# Input data:                                                               #
#    28 sample fits files constructed by Mark Calabretta                    #
#    for all standard projections                                           #
#                                                                           #
#############################################################################


casapath = os.environ['CASAPATH'].split()[0]
datapath = ""

# get the dataset name from the wrapper if possible

if 'datasets' not in (locals()):
    myname = 'fits-import-export_regression :'
    # The default dictionary of test datasets
    #  Each entry is a tuple of filename (w/o extension), expected pixel position of the maximum flux,
    #  expected absolute position of the maximum flux (RA, DEC), name.
    mydict = { 1: ('1904-66_AIR', [109, 167], '19:39:23.885, -63.45.36.905', 'Airy Projection (AIR)'),
               2: ('1904-66_AIT', [109, 168], '19:39:41.653, -63.43.54.147', 'Hammer-Aitoff Projection (AIT)'),  
               3: ('1904-66_ARC', [110, 170], '19:39:28.622, -63.41.53.659', 'Zenithal Equidistant Projection (ARC)'), 
               4: ('1904-66_AZP', [116, 186], '19:39:21.120, -63.44.26.642', 'Zenithal Perspective Projection (AZP)'), 
               5: ('1904-66_BON', [108, 173], '19:39:28.718, -63.41.12.383', 'Bonne\'s Equal Area Projection (BON)'), 
               6: ('1904-66_CAR', [113, 168], '19:39:42.371, -63.41.36.035', 'Plate Caree Projection (CAR)'), 
               7: ('1904-66_CEA', [113, 167], '19:39:35.136, -63.41.56.055', 'Cylindrical Equal Area Projection (CEA)'),
               8: ('1904-66_COD', [109, 166], '19:39:39.760, -63.42.02.640', 'Conic Equidistant Projection (COD)'), 
               9: ('1904-66_COE', [112, 172], '19:39:34.041, -63.44.23.296', 'Conic Equal-Area Projection (COE)'),
               10: ('1904-66_COO', [109, 161], '19:39:31.237, -63.44.09.556', 'Conic Orthomorphic Projection (COO)'), 
               11: ('1904-66_COP', [110, 161], '19:39:28.345, -63.44.40.626', 'Conic Perspective Projection (COP)'), 
               12: ('1904-66_CSC', [113, 180], '19:39:41.073, -63.43.25.624', 'COBE Quadrilateralized Spherical Cube Projection (CSC)'), 
               13: ('1904-66_CYP', [108, 157], '19:39:12.028, -63.43.07.315', 'Cylindrical Perspective Projection (CYP)'),  
               14: ('1904-66_HPX', [113, 179], '19:39:16.552, -63.42.47.347', 'HEALPix Grid Projection (HPX)'), 
               15: ('1904-66_MER', [113, 168], '19:39:16.276, -63.42.48.107', 'Mercator Projection (MER)'), 
               16: ('1904-66_MOL', [109, 175], '19:39:20.341, -63.41.44.201', 'Mollweide Projection (MOL)'), 
               17: ('1904-66_NCP', [107, 167], '19:39:38.614, -63.42.51.577', 'North Celetial Pole (SIN spcial case) Projection (NCP)'), 
               18: ('1904-66_PAR', [109, 171], '19:39:32.698, -63.42.04.737', 'Parabolic Projection (PAR)'), 
               19: ('1904-66_PCO', [108, 174], '19:39:21.403, -63.43.49.358', 'Polyconic Projection (PCO)'), 
               20: ('1904-66_QSC', [120, 182], '19:39:23.808, -63.41.22.666', 'Quadrilateralized Spherical Cube Projection (QSC)'), 
               21: ('1904-66_SFL', [108, 167], '19:39:16.950, -63.45.15.188', 'Samson-Flamsteed Projection (SFL)'), 
               22: ('1904-66_SIN', [107, 167], '19:39:38.614, -63.42.51.577', 'Slant Orthographic Projection (SIN)'), 
               23: ('1904-66_STG', [111, 171], '19:39:14.752, -63.44.20.882', 'Stereographic Projection (STG)'), 
               24: ('1904-66_SZP', [110, 177], '19:39:42.475, -63.42.13.751', 'Slant Zenithal Perspective Projection (SZP)'),
               25: ('1904-66_TAN', [116, 177], '19:39:30.753, -63.42.59.218', 'Gnomonic Projection (TAN)'), 
               26: ('1904-66_TSC', [112, 160], '19:39:39.997, -63.41.14.586', 'Tangential Spherical Cube Projection (TSC)'), 
               27: ('1904-66_ZEA', [109, 169], '19:39:26.872, -63.43.26.060', 'Zenithal Equal Area Projection (ZEA)'), 
               28: ('1904-66_ZPN', [94, 150], '19:39:24.948, -63.46.43.636', 'Zenithal Polynomial Projection (ZPN)'), 
               29: ('1904-66_AIT-obsgeo', [109, 168], '19:39:41.653, -63.43.54.147', 'Hammer-Aitoff Projection (AIT)')
           }
else: # the script has been called from the wrapper, don't need to output name
    myname = ' '
    mydict = (locals())['datasets']


def checkimage(myfitsimage_name, maxpos_expect, maxposf_expect):
    global myname
    global datapath
    subtest_passed = True

    if not os.path.exists(datapath+myfitsimage_name+'.fits'):
        datapath = casapath + "/data/regression/fits-import-export/input/"

    # import the image
    default('importfits')
    try:
        print(myname, ' Importing ', datapath+myfitsimage_name+'.fits', ' ...')
        importfits(fitsimage = datapath+myfitsimage_name+'.fits',
                   imagename = myfitsimage_name,
                   overwrite = True)
    except:
        print(myname, ' Error ', sys.exc_info()[0])
        raise    
    else:
        print(myname, ' No exceptions raised! Now checking image ...')
        # perform a basic check of the coordinate system
        ia.open(myfitsimage_name)
        mystat = imstat(imagename = myfitsimage_name)
        ia.close()
        if not (myname == ' '):
            print(mystat)
        if not ((mystat['maxpos']) == maxpos_expect).all():
            print(myname, ' Error in imported image ', myfitsimage_name, ':')
            print(myname, '    expected pixel position of maximum is ', maxpos_expect)
            print(myname, '                               but found ', mystat['maxpos']) 
            subtest_passed = False    
        if not (mystat['maxposf'] == maxposf_expect):
            print(myname, ' Error in imported image ', myfitsimage_name, ':')
            print(myname, '   expected absolute position of maximum is ', maxposf_expect)
            print(myname, '                                  but found ', mystat['maxposf']) 
            subtest_passed = False 
        if subtest_passed:
            print(myname, ' imported image ', myfitsimage_name, ' as expected.')

            # export the image
            default('exportfits')
            try:
                print('Exporting ', myfitsimage_name, '...')
                exportfits(fitsimage = myfitsimage_name + 'exp.fits',
                           imagename = myfitsimage_name,
                           overwrite = True)
            except:
                print(myname, ' Error ', sys.exc_info()[0])
                raise    
            else:
                print(myname, ' No exceptions raised! Now comparing exported image with original by re-importing it ...')
                # re-import the image    
                default('importfits')
                try:
                    print(myname, ' Re-importing ', myfitsimage_name+'exp.fits', ' ...')
                    importfits(fitsimage = myfitsimage_name+'exp.fits',
                               imagename = myfitsimage_name+'exp',
                               overwrite = True)
                except:
                    print(myname, ' Error ', sys.exc_info()[0])
                    raise    
                else:
                    print(myname, ' No exceptions raised! Now checking image ...')
                    ia.open(myfitsimage_name+'exp')
                    csm = ia.coordsys()
                    ia.close()
                    csm.summary() 
                    mystat = imstat(imagename = myfitsimage_name+'exp')
                    if not ((mystat['maxpos']) == maxpos_expect).all():
                        print(myname, ' Error in re-imported image ', myfitsimage_name+'exp', ':')
                        print(myname, '   expected pixel position of maximum is ', maxpos_expect)
                        print(myname, '                               but found ', mystat['maxpos']) 
                        subtest_passed = False    
                    if not (mystat['maxposf'] == maxposf_expect):
                        print(myname, ' Error in re-imported image ', myfitsimage_name+'exp', ':')
                        print(myname, '   expected absolute position of maximum is ', maxposf_expect)
                        print(myname, '                                  but found ', mystat['maxposf']) 
                        subtest_passed = False 
                    if subtest_passed:
                        print(myname, ' re-imported image ', myfitsimage_name+'exp', ' as expected.')
    return subtest_passed
# end def checkimage()

def checkimageb(myfitsimage_name):
    global myname
    subtest_passed = True

    # import the image
    default('importfits')
    try:
        print(myname, ' Importing ', myfitsimage_name+'.fits', ' ...')
        importfits(fitsimage = datapath+myfitsimage_name+'.fits',
                   imagename = myfitsimage_name,
                   overwrite = True)
    except:
        print(myname, ' Error ', sys.exc_info()[0])
        raise    
    else:
        print(myname, ' No exceptions raised!')
        mystat = imstat(imagename = myfitsimage_name)
        # export the image
        default('exportfits')
        try:
            print('Exporting ', myfitsimage_name, ' with bitpix=16 ...')
            exportfits(fitsimage = myfitsimage_name + 'exp.fits',
                       imagename = myfitsimage_name,
                       overwrite = True, bitpix=16) # with BITPIX=16 !
        except:
            print(myname, ' Error ', sys.exc_info()[0])
            raise    
        else:
            print(myname, ' No exceptions raised! Now testing minimum and maximum values ...')
            # re-import the image    
            default('importfits')
            try:
                print(myname, ' Re-importing ', myfitsimage_name+'exp.fits', ' ...')
                importfits(fitsimage = myfitsimage_name+'exp.fits',
                           imagename = myfitsimage_name+'exp',
                           overwrite = True)
            except:
                print(myname, ' Error ', sys.exc_info()[0])
                raise    
            else:
                print(myname, ' No exceptions raised! Now checking image ...')
                myotherstat = imstat(imagename = myfitsimage_name+'exp')
                if not (abs((mystat['min'] - myotherstat['min'])/mystat['min']) < 1E-6):
                    print(myname, ' Error in re-imported image ', myfitsimage_name+'exp', ':')
                    print(myname, '   expected  minimum is ', mystat['min'])
                    print(myname, '                               but found ', myotherstat['min']) 
                    subtest_passed = False    
                if not (abs((mystat['max'] - myotherstat['max'])/mystat['max']) < 1E-6):
                    print(myname, ' Error in re-imported image ', myfitsimage_name+'exp', ':')
                    print(myname, '   expected  maximum is ', mystat['max'])
                    print(myname, '                               but found ', myotherstat['max']) 
                    subtest_passed = False    
                if subtest_passed:
                    print(myname, ' re-imported image ', myfitsimage_name+'exp', ' as expected.')
    return subtest_passed
# end def checkimage()

failed_tests = []
passed_tests = []

for i in list(mydict.keys()):

    thefitsimage_name = mydict[i][0]
    themaxpos_expect = mydict[i][1]
    themaxposf_expect = mydict[i][2]
    print(myname, ' ***********************************************************')
    print(myname, ' Subtest ', i, ': ', mydict[i][3])
    passed = checkimage(thefitsimage_name, themaxpos_expect, themaxposf_expect)
    if passed:
        print(myname, ' Subtest ', i, ' passed.')
        passed_tests.append(mydict[i][3])
    else:
        print(myname, ' Subtest ', i, ' failed.')
        failed_tests.append(mydict[i][3])
    print(myname, ' ***********************************************************')

thefitsimage_name = mydict[1][0]
print(myname, ' ***********************************************************')
print(myname, ' Test of the BITPIX parameter:')
passed = checkimageb(thefitsimage_name)
if passed:
    print(myname, ' bitpix test passed.')
    passed_tests.append('bitpix')
else:
    print(myname, ' bitpix test failed.')
    failed_tests.append('bitpix')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of the stokeslast parameter and the SPECSYS keyword:')
exportfits(imagename=datapath+'stokeslast-test.image', fitsimage='stokeslast-test.fits', stokeslast=True, overwrite=True)
myresult0 = os.system('grep SPECSYS stokeslast-test.fits > /dev/null')
specsyspresent = (myresult0 == 0)
importfits(imagename='stokeslast-test2.image', fitsimage='stokeslast-test.fits', overwrite=True)
myrgn1 = rg.box([0,0,1,0],[64,64,1,0])
myrgn2 = rg.box([0,0,0,1],[64,64,0,1])
ia.open(datapath+'stokeslast-test.image')
ia.subimage(outfile='sub1.im', region = myrgn1, overwrite=True)
ia.close()
ia.open('stokeslast-test2.image')
ia.subimage(outfile='sub2.im', region = myrgn2, overwrite=True)
ia.close()
myresult1 = imstat('sub1.im', verbose=False)
myresult2 = imstat('sub2.im', verbose=False)
# imagecalc is on strike here because the formal coordinates of the slices disagree because the order is different
# so use min, max, and sum
passed = (myresult1['min']==myresult2['min']) and (myresult1['max']==myresult2['max']) \
         and (myresult1['sum']==myresult2['sum']) and specsyspresent
if passed:
    print(myname, ' stokeslast test passed.')
    passed_tests.append('stokeslast')
else:
    print(myname, ' stokeslast test failed.')
    failed_tests.append('stokeslast')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of the wavelength parameter:')
os.system('rm -f wavelength-test.fits')
ia.open(datapath+'stokeslast-test.image')
ia.tofits(outfile='wavelength-test.fits', wavelength=True)
ia.close()
passed = ia.open('wavelength-test.fits')
ia.close()

os.system('rm -f wavelength-test2.fits')
ia.open(datapath+'stokeslast-test.image')
ia.tofits(outfile='wavelength-test2.fits', wavelength=True, airwavelength=True)
ia.close()
passed = passed and ia.open('wavelength-test2.fits')
ia.close()

if passed:
    print(myname, ' wavelength parameter test passed.')
    passed_tests.append('wavelength')
else:
    print(myname, ' wavelength parameter test failed.')
    failed_tests.append('wavelength')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of export of a standard single-channel image from clean:')
os.system('rm -f xxx-test.fits xxx-test-w.fits')
ia.open(datapath+'xxx-clean.image')
ia.tofits(outfile='xxx-test.fits')
ia.tofits(outfile='xxx-test-w.fits', wavelength=True)
ia.close()
passed1 = ia.open('xxx-test.fits')
ia.close()
passed2 = ia.open('xxx-test-w.fits')
ia.close()
passed = passed1==True and passed2==True
if passed:
    print(myname, ' trivial export tests passed.')
    passed_tests.append('trivial')
else:
    print(myname, ' trivial export test failed.')
    failed_tests.append('trivial')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of import and export of FITs images with spectral WCS:')
expecta = {'freq': 2.190956850600E+11,
           'vrad': 2.183247880000E+11,
           'vopt': 2.183247880000E+11,
           'wave': 2.204356308820E+11,
           'awav': 2.203722484080E+11}
expectb = {'freq': 'Hz',
           'vrad': 'Hz',
           'vopt': 'Hz',
           'wave': 'Hz',
           'awav': 'Hz'}
expectc = {'freq': 'LSRK',
           'vrad': 'BARY',
           'vopt': 'LSRD',
           'wave': 'GALACTO',
           'awav': 'CMB'}
passedx = {'freq': False,
           'vrad': False,
           'vopt': False,
           'wave': False,
           'awav': False}
passed = True
#for ctype in ['freq','vrad','wave']:
for myctype in ['freq','vrad','vopt','wave','awav']:
    try:
        print('Testing CTYPE ', myctype, ' ...')
        importfits(imagename=myctype+'.im', fitsimage=datapath+'spec-test-'+myctype+'.fits', overwrite=True)
        cond0 = ia.open(myctype+'.im')
        coordm = ia.coordmeasures()
        cond1 = (abs(coordm['measure']['spectral']['frequency']['m0']['value']-expecta[myctype])<1.) # avoid Python precision problems
        print('value, expectation, diff:', coordm['measure']['spectral']['frequency']['m0']['value'],',',\
              expecta[myctype], ",", coordm['measure']['spectral']['frequency']['m0']['value']-expecta[myctype])
        cond2 = (coordm['measure']['spectral']['frequency']['m0']['unit']==expectb[myctype])
        cond3 = (coordm['measure']['spectral']['frequency']['refer']==expectc[myctype])
        print(cond0, cond1, cond2, cond3)
        passed1 = cond0 and cond1 and cond2 and cond3
        ia.close()
        if(myctype=='freq'):
            exportfits(imagename=myctype+'.im', fitsimage='spec-test-'+myctype+'-ex.fits', overwrite=True)
        elif(myctype=='vrad'):
            exportfits(imagename=myctype+'.im', fitsimage='spec-test-'+myctype+'-ex.fits', velocity=True, optical=False, overwrite=True)
        elif(myctype=='vopt'):
            exportfits(imagename=myctype+'.im', fitsimage='spec-test-'+myctype+'-ex.fits', velocity=True, optical=True, overwrite=True)
        elif(myctype=='wave'):
            os.system('rm -f spec-test-'+myctype+'-ex.fits')
            ia.open(myctype+'.im')
            ia.tofits(outfile='spec-test-'+myctype+'-ex.fits', wavelength=True)
            ia.close()
        elif(myctype=='awav'):
            os.system('rm -f spec-test-'+myctype+'-ex.fits')
            ia.open(myctype+'.im')
            ia.tofits(outfile='spec-test-'+myctype+'-ex.fits', wavelength=True, airwavelength=True)
            ia.close()
        else:
            print("Skipping export test for ", myctype)

        passed2 = True

        importfits(imagename=myctype+'2.im', fitsimage='spec-test-'+myctype+'-ex.fits', overwrite=True)
        cond0 = ia.open(myctype+'2.im')
        coordm = ia.coordmeasures()
        cond1 = (abs(coordm['measure']['spectral']['frequency']['m0']['value']-expecta[myctype])<1.) # avoid Python precision problems
        print('value, expectation, diff:', coordm['measure']['spectral']['frequency']['m0']['value'],',',\
              expecta[myctype], ",", coordm['measure']['spectral']['frequency']['m0']['value']-expecta[myctype])
        cond2 = (coordm['measure']['spectral']['frequency']['m0']['unit']==expectb[myctype])
        cond3 = (coordm['measure']['spectral']['frequency']['refer']==expectc[myctype])
        print(cond0, cond1, cond2, cond3)
        passed2 = cond0 and cond1 and cond2 and cond3
        ia.close()
            
        passedx[myctype] = passed1 and passed2
        passed = passed and passedx[myctype]
    except:
        print(myname, ' Error ', sys.exc_info()[0])
        passedx[myctype] = False
        passed = False
        
if passed:
    print(myname, ' Spectral WCS im/export tests passed.')
    passed_tests.append('spectralwcs')
else:
    print(myname, ' Spectral WCS im/export tests failed.')
    print(passedx)
    failed_tests.append('spectralwcs')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of import with use of parameter defaultaxes:')
passed = True
for name in ['dir_and_stokes.fits', 'dir_stokes_and_freq.fits', 'onlydir.fits', 'dir_and_freq.fits', 'dir_freq_and_stokes.fits', 'freq.fits']:
    print(name)
    importfits(fitsimage=datapath+name, imagename=name+'.im', overwrite=True, defaultaxes=True, defaultaxesvalues=['19h30m00',2.,'6GHz', 'Q'])
    ia.open(name+'.im')
    mycs = ia.coordsys()
    ia.close()
    if not mycs.axiscoordinatetypes() == ['Direction', 'Direction', 'Stokes', 'Spectral']:
        print("ERROR: unexpected value of axiscoordinatetypes: ", mycs.axiscoordinatetypes())
        passed = False
    else:
        mycsr = mycs.torecord()
        if name in ['dir_and_stokes.fits', 'onlydir.fits']:
            if not (('spectral1' in mycsr and mycsr['spectral1']['wcs']['crval'] == 6E9) \
                   or ('spectral2' in mycsr and mycsr['spectral2']['wcs']['crval'] == 6E9)):
                print("Error: Spectral axis CRVAL should be 6E9")
                passed = False
        if name in ['onlydir.fits', 'dir_and_freq.fits', 'freq.fits']:
            if not (('stokes1' in mycsr and mycsr['stokes1']['stokes'][0]=='Q') \
                    or ('stokes2' in mycsr and mycsr['stokes2']['stokes'][0]=='Q')):
                print("Error: Stokes axis CRVAL should be Q")
                passed = False
        if name in ['freq.fits']:
            if not ('direction1' in mycsr and (mycsr['direction1']['crval'][0]==292.5)):
                print("Error: RA axis CRVAL should be 292.5")
                passed = False

if passed:
    print(myname, ' tests of import with defaultaxes passed.')
    passed_tests.append('defaultaxes')
else:
    print(myname, ' tests of import with defaultaxes failed.')
    failed_tests.append('defaultaxes')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of export of FITs images with spectral conversion layer:')
expecta = 2.191033189757E+11 # BARY
expectb = 2.190956850603E+11 # LSRK
passed = True
try:
    os.system('rm -rf freq.im')
    importfits(imagename='freq.im', fitsimage=datapath+'spec-test-freq.fits', overwrite=True)
    os.system('rm -rf freqx.im')
    imreframe(imagename='freq.im', output='freqx.im', outframe='BARY')
    cond0 = ia.open('freqx.im')
    mycs1 = ia.coordsys()
    ia.close()
    myworld = mycs1.toworld([0,0,0,0])['numeric'][2]
    cond1 = (abs(myworld-expecta)<1.) # avoid Python precision problems
    print('value, expectation, diff:', myworld,',', expecta, ",", myworld-expecta)
    cond2 = (mycs1.referencecode()[2]=='LSRK')
    print(cond0, cond1, cond2)
    passed1 = cond0 and cond1 and cond2
    os.system('rm -rf spec-test-convbary-ex.fits')
    exportfits(imagename='freqx.im', fitsimage='spec-test-convbary-ex.fits')
    
    passed2 = True
    
    os.system('rm -rf freq2.im')
    importfits(imagename='freq2.im', fitsimage='spec-test-convbary-ex.fits', overwrite=True)
    cond0 = ia.open('freq2.im')
    mycs2 = ia.coordsys()
    ia.close()
    myworld = mycs2.toworld([0,0,0,0])['numeric'][2]
    cond1 = (abs(myworld-expecta)<1.) # avoid Python precision problems
    print('value, expectation, diff:', myworld,',', expecta, ",", myworld-expecta)
    cond2 = (mycs2.referencecode()[2]=='BARY')
    print(cond0, cond1, cond2)
    passed2 = cond0 and cond1 and cond2
    
    passed = passed1 and passed2

except:
    print(myname, ' Error ', sys.exc_info()[0])
    passed = False
        
if passed:
    print(myname, ' Export with conversion layer tests passed.')
    passed_tests.append('Export with conversion layer')
else:
    print(myname, ' Export with conversion layer tests failed.')
    failed_tests.append('Export with conversion layer')
print(myname, ' ***********************************************************')

print(myname, ' ***********************************************************')
print(myname, ' Test of the beam parameter:')
importfits(fitsimage=datapath+'1904-66_AIT.fits', imagename='1904-66_AIT_2.im', beam=['2arcsec','1arcsec','3.0deg'], overwrite=True)
passed = False
ia.open('1904-66_AIT_2.im')
x = ia.restoringbeam()
ia.close()
if(x != {} 
   and x['major']['value']==2. and  x['major']['unit']=='arcsec'
   and x['minor']['value']==1. and  x['minor']['unit']=='arcsec'
   and x['positionangle']['value']==3. and  x['positionangle']['unit']=='deg'):
    passed = True

if passed:
    print(myname, ' beam parameter test passed.')
    passed_tests.append('beam')
else:
    print(myname, ' beam parameter test failed.')
    failed_tests.append('beam')
print(myname, ' ***********************************************************')


if len(failed_tests)>0:
    if len(passed_tests)>0:
        print(myname, ' Summary of passed subtests:')
        for i in range(0, len(passed_tests)):
            print(myname, '          ', passed_tests[i])            
    print(myname, ' Summary of failed subtests:')
    for i in range(0, len(failed_tests)):
        print(myname, '          ', failed_tests[i])
    raise Exception(repr(len(failed_tests))+' out of '+repr(len(list(mydict.keys())))+' subtests failed.')
else:
    print('All subtests passed.')
    print('')
    print('Regression PASSED')
    print('')
print('Done.')


#End of Script
