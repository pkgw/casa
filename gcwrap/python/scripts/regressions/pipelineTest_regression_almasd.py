import os
import time
import datetime
import regression_utility as regutl
import sys, traceback
import numpy as np
import pipeline
import pipeline.infrastructure.casatools as casatools

pathname = os.environ.get('CASAPATH').split()[0]
rootdatapath = pathname+'/data/regression/pipeline/almasd/'

# ALMA Singledish pipeline regression
# modifled from VLA pipeline regression
# 2019-Mar-07 H. Ezawa
# 2019-Mar-08 H. Ezawa : Updated the goal values for CASA 5.5
# 2019-Mar-18 H. Ezawa : minor bug(s) fixed

'''Initial VLA pipeline regression
   B. Kent, May 2015
   Update September 1, 2015
   Update April 20, 2018
   Update June  01, 2018
'''

THISHOME  = "working/"
startTime=0.0
endTime=0.0
startProc=0.0
endProc=0.0
regstate = True
# standard_file = rootdatapath + 'VLApipeline44-standard'
# MIN_CASA_REVISION = 36095


def print_log( logfile, text ):
    print( text )
    print(text, file=logfile)
    ## for Python3 replace this with the following
    ## print( text, file=logfile )

def load_context(filename):
    with open(filename, 'rb') as picklefile:
        return pipeline.infrastructure.utils.pickle_load(picklefile)

#
#--------------------------------------------------------------
#
def pipeline_regression():
    global regstate
    # global MIN_CASA_REVISION
        
    
    #revision = int(casadef.subversion_revision)
    #if MIN_CASA_REVISION > revision:
    #    msg = ('Minimum CASA revision for the pipeline is r%s, '
    #           'CASA %s (r%s).' % (MIN_CASA_REVISION,
    #            cu.version_info( ),
    #            casadef.subversion_revision))
    #    print msg
    #    regstate = False
    #    raise EnvironmentError(msg)

    data_location = rootdatapath
    ASDM = data_location + "uid___A002_X85c183_X36f/"
    jyperk_file = 'jyperk.csv'
    jyperk = data_location + jyperk_file

    try:
        import pipeline.recipes.hsd as hsd
    except ImportError as e:
        print(e)
    
    # Check to see if the ASDM exists
    if not os.path.exists(ASDM):
        print(( "Unable to open ASDM {}".format(ASDM) ))
	regstate=False
        raise IOError
    else:
        print(( "Using {}".format(ASDM) ))
    
# Run the Pipeline standard recipe

    # copy jyperk file 
    try:
        shutil.copy2( jyperk, jyperk_file )
    except:
        print(( "Could not copy {0} file from {1}".format(jyperk_file, jyperk) ))
        regstate=False
        raise IOError

    # verify if jyperk.csv exists
    if not os.path.exists( jyperk ):
        print(( "Jy/K file {} not found.".format(jyperk) ))
        regstate=False
        raise IOError

    # run the receipe
    try:
        hsd.hsdms([ASDM], importonly=False)
    except Exception as e:
        print(e)
        regstage=False
        raise IOError

def run():
    '''Timing and running of main regression function
    '''
    global startTime, endTime, startProc, endProc, regstate;
    startTime = time.time()
    startProc = time.clock()
    regstate = False
    pipeline_regression()
    endTime = time.time()
    endProc = time.clock()
    print(("Run Time = {}".format(endTime-startTime,endProc-startProc)))

def stats():
    global startTime, endTime, startProc, endProc, regstate, standard_context_file
    
    datestring = datetime.datetime.isoformat(datetime.datetime.today())
    outfile = 'sd_pipeline.'+datestring+'.log'
    try:
        logfile = open(outfile,'w')
    except:
        print(( "Cannot open logfile \"{}\"".format(outfile) ))
        raise IOError

    casa_version = casalog.version()

    # import pipeline 
    try:
        import pipeline
    except ImportError as e:
        print_log( logfile, e )
        print_log( logfile, "Unable to import the CASA pipeline")
        regstate = False
    
    # Open context
    try:
        context = pipeline.Pipeline(context='last').context
###        context = h_resume()
###        print "!!!!! WARNING !!!!! TEST MODE !!!!!"

        regstate = True
        print_log( logfile, "Single Dish pipeline context created" )
        print_log( logfile, "Context verification - Single Dish pipeline regression PASSED." )
    except Exception as e:
        regstate = False
        print_log( logfile, "Single Dish pipeline context NOT created." )
        print_log( logfile, "Context verification - Single Dish pipeline regression FAILED." )
        print_log( logfile, e )

    # read context
    try:
        # Open context of regression pipeline run
        imageitem = context.results[11].read()[-1].outcome['image']
    except Exception as e:
        print_log( logfile, "Failed to read contet.results[11]" )
        print_log( logfile, e )

    # actual test for single-dish
    try:
        # parameters for analysis
        integ_radius = 5     # radius to average (in pixels)
                             # equivalent to beam sizae (9 pix)

        # define goal values and tolerances #
        im_text={}
        im_goal={} 
        im_rtol={}
        im_atol={}
        im_act={}

        im_text['pixel'] = "Peak position (pixel)"
        im_goal['pixel'] = (17,27)
        im_atol['pixel'] = ( 1, 1)
        im_rtol['pixel'] = ( 0, 0)   # n/a

        im_text['direction'] = "Peak direction (deg)"
        # update reference due to PIPE-313
        # up to r42407
        # im_goal['direction'] = me.direction( 'J2000', '-174.271746303deg', '15.8228893587deg' )
        # r42408 or later
        im_goal['direction'] = me.direction( 'J2000', '-174.27166297deg', '15.8229093587deg' )
        im_atol['direction'] = 0.01/3600.0   # separation angle in deg 
        im_rtol['direction'] = 0            # n/a

        im_text['flux'] = "Integrated Flux (Jy/beam)"       
        # update reference due to PIPE-313
        # up to r42407
        # im_goal['flux'] = 0.359638377757
        # r42408 or later
        im_goal['flux'] = 0.359533426907
        im_atol['flux'] = 1E-8
        im_rtol['flux'] = 1E-5

        sp_text={}
        sp_goal={} 
        sp_rtol={}
        sp_atol={}
        sp_act={}

        sp_text['channel'] = "Peak channel (ch)"
        sp_goal['channel'] = 2033
        sp_atol['channel'] = 1
        sp_rtol['channel'] = 0    # n/a

        sp_text['frequency'] = "Peak frequency (GHz)"
        # update reference due to PIPE-313
        # up to r42407
        # sp_goal['frequency'] = 114.687933618  # in GHz
        # r42408 or later
        sp_goal['frequency'] = 114.687933607  # in GHz
        sp_atol['frequency'] = 1E-8
        sp_rtol['frequency'] = 0   # n/a

        sp_text['flux'] = "Peak flux (Jy/beam)"
        # update reference due to PIPE-313
        # up to r42407
        # sp_goal['flux'] = 6.74114826797
        # r42408 or later
        sp_goal['flux'] = 6.73362001474
        sp_atol['flux'] = 1E-8
        sp_rtol['flux'] = 1E-5

        # open ia
        try:
            ia.open( imageitem.imagename )
            print_log( logfile, "Using image \"{}\"".format(imageitem.imagename) )
        except:
            print_log( logfile, "Unable to open image \"{}\"".format(imageitem.imagename) )
            regstate=False
            raise IOError

        # get image
        try:
            imdata  = ia.getchunk()
            imdata2 = imdata[:,:,:,(imdata.shape[3]/4):(imdata.shape[3]*3/4)]
            imdata2_2d = np.sum( np.sum( imdata2, axis=3 ), axis=2 )
        except:
            print_log( logfile, "Unable to get image data" )
            regstate=False
            raise IOError

        # get mask
        try:
            immask = ia.getchunk(getmask=True)
            immask2 = immask[:,:,:,(immask.shape[3]/4):(immask.shape[3]*3/4)]
            immask2_2d = 1.0*(np.sum( np.sum( immask2, axis=3 ), axis=2 )).astype(float)
        except:
            print_log( logfile, "Unable to get image data" )
            regstate=False
            raise IOError

        # find peak pixel and get its direction
        (cx0, cy0) = np.where( imdata2_2d == np.max(imdata2_2d) )
        (cx, cy) = (cx0[0], cy0[0])
        im_act['pixel'] = (cx, cy)

        meas = ia.coordmeasures([cx,cy,0,0])
        peak_direction = meas['measure']['direction']
        im_act['direction'] = peak_direction

        # generate integ_region mask
        integ_region = imdata2_2d.copy()
        for ix in range(imdata.shape[0]):
            for iy in range(imdata.shape[1]):
                integ_region[ix][iy] = 0
                if ( (cx-ix)**2 + (cy-iy)**2 ) < integ_radius**2:
                    integ_region[ix][iy] = 1
        
        # calc integrated map
        imdata2_2d_region = imdata2_2d * integ_region
        immask2_2d_region = immask2_2d * integ_region
        peak_average = np.sum(imdata2_2d_region) / np.sum(immask2_2d_region)
        im_act['flux'] = peak_average

        # integrate peak over beam-size for each channel
        spec_average = np.zeros( imdata.shape[3], dtype=imdata.dtype )
        npixel = 0
        for ix in range(imdata.shape[0]):
            for iy in range(imdata.shape[1]):
                if immask2_2d_region[ix][iy] > 0:
                    npixel += 1
                    spec = imdata[ix,iy,0,:]
                    spec_average += spec

        spec_average /= npixel
        peak_ch = np.where( spec_average ==  max(spec_average) )[0][0]
        sp_act['channel'] = peak_ch

        qa = pipeline.infrastructure.casatools.quanta
        peak_freq = ia.coordmeasures([cx,cy,0,peak_ch])['measure']['spectral']['frequency']['m0']
        sp_act['frequency'] = qa.convert( peak_freq, 'GHz')['value']

        peak_flux = spec_average[peak_ch]
        sp_act['flux'] = peak_flux

        # close ia
        ia.close()

        # evaluate the values for im_*
        im_result = {}
        for key in list(im_act.keys()):
            if key == 'direction':
                separation = me.separation( im_act['direction'], im_goal['direction'] )
                separation_deg = qa.convert( separation, 'deg' )['value']
                im_result['direction'] = np.isclose( separation_deg, 0.0, atol=im_atol['direction'], rtol=im_rtol['direction'] )
            else:
                im_result[key] = np.isclose( im_act[key], im_goal[key], atol=im_atol[key], rtol=im_rtol[key] )
                if type(im_result[key]) is not bool:
                    im_result[key] = all(im_result[key])

        # evaluate the values for sp_*
        sp_result = {}
        for key in list(sp_act.keys()):
            sp_result[key] = np.isclose( sp_act[key], sp_goal[key], atol=sp_atol[key], rtol=sp_rtol[key] )

        # summary
        print_log( logfile, "#### Results ####" )
        print_log( logfile, " integ_radius < {} pixels".format(integ_radius) )
        print_log( logfile, "  with integrated half-spw image:" )
        for key in list(im_result.keys()):
            print_log( logfile, "    {}".format(im_text[key]) )
            if ( key == "pixel" ):
                print_log( logfile, "      Actual value = ({0}, {1})".format(*im_act[key]) )
                print_log( logfile, "      Goal value   = ({0}, {1})".format(*im_goal[key]) )
            elif ( key == "direction" ):
                ra  = qa.convert( im_act[key]['m0'], 'deg')['value']
                dec = qa.convert( im_act[key]['m1'], 'deg')['value']
                print_log( logfile, "      Actual value = ({0}, {1}) ({2})".format(ra, dec, im_act[key]['refer']) )
                ra  = qa.convert( im_goal[key]['m0'], 'deg')['value']
                dec = qa.convert( im_goal[key]['m1'], 'deg')['value']
                print_log( logfile, "      Goal value   = ({0}, {1}) ({2})".format(ra, dec, im_goal[key]['refer']) )
                print_log( logfile, "      Separation   = {}".format(separation_deg) )
            else:
                print_log( logfile, "      Actual value = {}".format(im_act[key]) )
                print_log( logfile, "      Goal value   = {}".format(im_goal[key]) )
            print_log( logfile, "      Tolerances   = absolute {0}, relative {1}".format( im_atol[key], im_rtol[key] ) )
            print_log( logfile, "      Result       = {}".format(im_result[key]) )

        print_log( logfile, "  with integrated spectrum:" )
        for key in list(sp_result.keys()):
            print_log( logfile, "    {}:".format(sp_text[key]) )
            if ( key == "channel" ):
                print_log( logfile, "      Actual value = {}".format(sp_act[key]) )
                print_log( logfile, "      Goal value   = {}".format(sp_goal[key]) )
                print_log( logfile, "      Tolerances   = absolute {0}, relative {1}".format( sp_atol[key], sp_rtol[key] ) )
            else:
                print_log( logfile, "      Actual value = {}".format(sp_act[key]) )
                print_log( logfile, "      Goal value   = {}".format(sp_goal[key]) )
                print_log( logfile, "      Tolerances   = absolute {0}, relative {1}".format( sp_atol[key], sp_rtol[key] ) )
            print_log( logfile, "      Result       = {}".format(sp_result[key]) )

        regstate_bool = all(im_result.values()) and all(sp_result.values())

        print_log( logfile, "#### Summary ####" )
        uname = os.uname()
        print_log( logfile, " Platform : {0} {2} {3} {4}".format(*uname) )
        print_log( logfile, " CASA :     {}".format(casa_version) )
        if ( regstate_bool ):
            print_log( logfile, " Values are within tolerances." )
            print_log( logfile, " Sigle dish regression passed." )
            regstate = True
        else:
            print_log( logfile, " Values are **NOT** within tolerances." )
            print_log( logfile, " Sigle dish regression **NOT** passed." )
            regstate = False
            
    except Exception as e:
        regstate=False
        print_log( logfile, e )

    logfile.close()


def main():
    try:
        run()
        stats()

    except KeyboardInterrupt:
        print("Interrupt requested...exiting")
    except Exception:
        traceback.print_exc(file=sys.stdout)

if __name__ == "__main__":
    main()
    print("Regstate:" , regstate)
    if regstate:
        print("Regression PASSED")
    else:
        print("Regression FAILED")
