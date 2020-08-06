import os
import time
import regression_utility as regutl
import shutil
import sys
import pdb

THISHOME  = "ic2233_regression_data/";
IMSIZE    = 2048;
#IMSIZE    = 4096;
CELLSIZE  = "4arcsec";
STOKES    = "IV";
CHSTART   = 69;
NCHAN     = 11;
#CHSTART = 64;
#NCHAN =  1;
SPW   ="*:69~80";
FIELD  = "0";
PASTEP    = 360.0;
NITER     = 6000;
#IMAGE     = "imIC2233.pbm0";
DECONVOLVER = "hogbom";
FTMACHINE = "awproject"
PSTERM = False;
IMAGE     = "imIC2233."+FTMACHINE;
REUSEREPOS = False;
INTERACTIVE = False;
TIMERANGE = "";"*+0:30:0";
EPS       = 1e-5;  # Logical "zero"


#--------------------------------------------------------------
#
def ic2233_reg():
    MSFILE    = THISHOME  + "ic2233.lband.ms";
    MYIMAGE   = THISHOME  + IMAGE;

    #os.system("rm -rf "+MYIMAGE+"*");
    if (REUSEREPOS==False):
    	os.system("rm -rf "+THISHOME);
    	os.mkdir(THISHOME);
    	REPOSNAME = os.environ.get('CASAPATH').split()[0]
    	PREFIX    = REPOSNAME + "/data/regression/ic2233/";
    	FITSFILE  = PREFIX    + "ic2233.lband.fits";
        ms.fromfits(fitsfile=FITSFILE,msfile=MSFILE);
        ms.done();

    tclean(vis=MSFILE,selectdata=True,field=FIELD,spw=SPW,timerange=TIMERANGE,uvrange="",antenna="",scan="",observation="",
           intent="",datacolumn="corrected",imagename=MYIMAGE,imsize=IMSIZE,cell=CELLSIZE,phasecenter="",stokes=STOKES,projection="SIN",
           startmodel="",specmode="mfs",reffreq="",nchan=-1,start="",width="",outframe="LSRK",veltype="radio",restfreq=[],interpolation="linear",
           gridder=FTMACHINE,facets=1,chanchunks=1,wprojplanes=1,vptable="",aterm=True,psterm=PSTERM,wbawp=True,conjbeams=False,
           cfcache=MYIMAGE+".cf",computepastep=360.0,rotatepastep=5.0,pblimit=0.05,normtype="flatnoise",deconvolver=DECONVOLVER,scales=[],
           nterms=2,smallscalebias=0.6,restoration=True,restoringbeam=[],pbcor=False,outlierfile="",weighting="natural",robust=0.5,npixels=0,
           uvtaper=[],niter=NITER,gain=0.1,threshold=0.0,cycleniter=-1,cyclefactor=1.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=INTERACTIVE,
           usemask="user",mask="",pbmask=0.0,restart=True,savemodel="none",
           calcres=True,calcpsf=True,parallel=False);

#    tclean(vis=MSFILE,selectdata=True,field="0",spw="*:69~80",timerange="",uvrange="",antenna="",scan="",observation="",intent="",datacolumn="corrected",imagename=MYIMAGE,imsize=2048,cell="4arcsec",phasecenter="",stokes="IV",projection="SIN",startmodel="",specmode="mfs",reffreq="",nchan=-1,start="",width="",outframe="LSRK",veltype="radio",restfreq=[],interpolation="linear",gridder="awproject",facets=1,chanchunks=1,wprojplanes=1,vptable="",aterm=True,psterm=False,wbawp=True,conjbeams=True,cfcache="junk_newVR.cf",computepastep=360.0,rotatepastep=5.0,pblimit=0.2,normtype="flatnoise",deconvolver="hogbom",scales=[],nterms=2,smallscalebias=0.6,restoration=True,restoringbeam=[],pbcor=False,outlierfile="",weighting="natural",robust=0.5,npixels=0,uvtaper=[],niter=6000,gain=0.1,threshold=0.0,cycleniter=-1,cyclefactor=1.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=True,usemask="user",mask="",pbmask=0.0,maskthreshold="",maskresolution="",nmask=0,autoadjust=False,restart=True,savemodel="none",calcres=True,calcpsf=True,parallel=False);

        


    # im.open(MSFILE,usescratch=True);
    
    # im.selectvis(nchan=NCHAN,start=CHSTART,step=1,spw=SPW,time=TIMERANGE);
    # im.defineimage(nx=IMSIZE,ny=IMSIZE,cellx=CELLSIZE,celly=CELLSIZE,stokes=STOKES,
    #                nchan=1,start=CHSTART,step=NCHAN-1,phasecenter=0);
    # im.setoptions(cache=IMSIZE*IMSIZE*4,ftmachine=FTMACHINE,
    #               applypointingoffsets=False, dopbgriddingcorrections=True,
    #               cfcachedirname=MYIMAGE+".cf", pastep=360.0);
    # im.clean(algorithm="cs",niter=NITER,interactive=INTERACTIVE,
    #          model=MYIMAGE+".mod",
    #          image=MYIMAGE+".image",
    #          residual=MYIMAGE+".res");
    # im.done();
#
#--------------------------------------------------------------
# The objective truth!
#
#StokesIPeak         = 0.86263674;  # Jy/beam
#StokesIRMS          = 0.00138958;  # Jy/beam
#
# The following was changed to new values on July 14, 2009
# StokesIPeak         = 0.86230296;
# StokesIRMS          = 0.00138973;

# The following was changed to new values on Sept 17, 2009
# StokesIPeak          = 0.86240315;
# StokesIRMS           = 0.00138986;

# The following was changed to new values on Feb 23, 2010
# StokesIPeak          = 0.86264914  
# StokesIRMS           = 0.00138987

# The following was changed to new values on Mar. 23, 2011
#StokesIPeak          = 0.8555392
#StokesIRMS           = 0.00137654

# StokesIPeak          = 0.85550195
# StokesIRMS           = 0.00137647

# Change to this from the above after the change in 3rd party pkgs
# around Jan. 1st week, 2017.  Also using tclean now for the imaging
# instead of the old imager (im-tool).
#StokesIPeak          =  0.85711205
#StokesIRMS           =  0.00394761

# Changes to ensure merge of ARD-14 with casa-master

StokesIPeak          =  0.87293589
StokesIRMS           =  0.00401568


StokesIPeakPosWorld = '08:20:22.869, +44.40.38.993'; #J2000
StokesIPeakPos      = [942,1130,0,0]; #Pixels

#StokesVPeak         = 0.00063984;  # Jy/beam
#StokesVRMS          = 5.07250807e-05;  # Jy/beam
# The following was changed to new values on July 14, 2009
# StokesVPeak         = 0.00060623;  # Jy/beam
# StokesVRMS          = 5.06413598e-05;  # Jy/beam

# The following was changed to new values on Sept 17, 2009
# StokesVPeak         = 0.00061063;  # Jy/beam
# StokesVRMS          = 5.08003759e-05;  # Jy/beam

# The following was changed to new values on Feb 23, 2010
# StokesVPeak         = 0.00061679   # Jy/beam
# StokesVRMS          = 5.06523975e-05 # Jy/beam

# The following was changed to new values on Mar. 23, 2011
# StokesVPeak         = 0.00066234   # Jy/beam
# StokesVRMS          = 5.26905496e-05 # Jy/beam

# StokesVPeak         = 0.00065949   # Jy/beam
# StokesVRMS          = 5.26994445e-05 # Jy/beam

# Change to this from the above after the change in 3rd party pkgs
# around Jan. 1st week, 2017.  Also using tclean now for the imaging
# instead of the old imager (im-tool).
# StokesVPeak         =  0.00061497
# StokesVRMS          =  0.00012994

# Changes after fixing the w-term in prediction.  Stokes-V looks OK,
# but the peak and therefore the RMS has changed.
#StokesVPeak         = 0.00145028
#StokesVRMS          =  0.0001952

#Changes to ensure that ARD-14 merges with casa-master
StokesVPeak         =  0.00217099
StokesVRMS          =  0.00019424


StokesVPeakPosWorld = '08:11:29.219, +45.48.26.199'; #J2000
StokesVPeakPos      = [1184,1346,1,0]; #Pixels
#
#--------------------------------------------------------------
#
startTime=0.0;
endTime=0.0;
startProc=0.0;
endProc=0.0;
def run():
    global startTime, endTime, startProc, endProc;
    startTime = time.time();
    startProc = time.clock();
    regstate = False;
    ic2233_reg();
    endTime = time.time();
    endProc = time.clock();
    print("Run Time = ",endTime-startTime,endProc-startProc);

def stats():
    global startTime, endTime, startProc, endProc;
    try:
        # startTime = time.time();
        # startProc = time.clock();
        # regstate = False;
        # ic2233_reg();
        # endTime = time.time();
        # endProc = time.clock();
        #
        #--------------------------------------------------------------
        # Get the image statistics
        #
        ia.open(THISHOME+IMAGE+".image");
#        ibox=ia.setboxregion(blc=[0,0,0,0],trc=[2048,2048,0,0]);
#        vbox=ia.setboxregion(blc=[0,0,1,0],trc=[2048,2048,1,0]);
        ibox = rg.box(blc=[0,0,0,0],trc=[2048,2048,0,0]);
        vbox = rg.box(blc=[0,0,1,0],trc=[2048,2048,1,0]);
    
        istats=ia.statistics(region=ibox, list=True, verbose=True);
        vstats=ia.statistics(region=vbox, list=True, verbose=True);

        ia.done();
        #
        #--------------------------------------------------------------
        # The subjective truth!
        #
        import datetime
        datestring=datetime.datetime.isoformat(datetime.datetime.today())
        outfile='ic2233.'+datestring+'.log'
        logfile=open(outfile,'w')

        dIMax    = istats["max"]-StokesIPeak;
        dIRMS    = istats["rms"]-StokesIRMS;
        dIMaxPos = istats["maxpos"]-StokesIPeakPos;
        dVMax    = vstats["max"]-StokesVPeak;
        dVRMS    = vstats["rms"]-StokesVRMS;
        dVMaxPos = vstats["maxpos"]-StokesVPeakPos;

        print("Stokes-I Statistics:", file=logfile);
        print("-------------------------------------------------", file=logfile);
        print("Max = ",istats["max"], "  RMS = ",istats["rms"], file=logfile);
        print("MaxPos = ",istats["maxpos"], file=logfile);
        print("", file=logfile);
        print("Stokes-V Statistics:", file=logfile);
        print("-------------------------------------------------", file=logfile);
        print("Max = ",vstats["max"], "  RMS = ",vstats["rms"], file=logfile);
        print("MaxPos = ",vstats["maxpos"], file=logfile);
        print("", file=logfile);
        print("Stokes-I Delta-regression Statistics:", file=logfile)
        print("-------------------------------------------------", file=logfile);
        print("dMax = ",dIMax, "  dRMS = ",dIRMS, file=logfile);
        print("dMaxPos = ",dIMaxPos, file=logfile);
        print("", file=logfile);
        print("Stokes-V Delta-regression Statistics:", file=logfile)
        print("-------------------------------------------------", file=logfile);
        print("dMax = ",dVMax, "  dRMS = ",dVRMS, file=logfile);
        print("dMaxPos = ",dVMaxPos, file=logfile);
        print("", file=logfile);

        if ((abs(dIMax) < EPS) &
            (abs(dIRMS) < EPS) &
            (abs(dVMax) < EPS) &
            (abs(dVRMS) < EPS) &
            (dIMaxPos==0).all() 
            #        (dVMaxPos==0).all()
            ):
            regstate=True;
            print("IC2233 Regression passed.", file=logfile);
            print('')
            print('Regression PASSED')
            print('')
        else:
            regstate=False;
            print('')
            print('Regression FAILED')
            print('')
            print("IC2233 Regression failed.", file=logfile);
        
            print('', file=logfile)
            print('', file=logfile)
            print('********* Benchmarking *****************', file=logfile)
            print('*                                      *', file=logfile)
            print('Total wall clock time was: ', endTime - startTime, file=logfile)
            print('Total CPU        time was: ', endProc - startProc, file=logfile)

        logfile.flush()
        logfile.close();
        f=open(outfile, 'r')
        a=f.read()
        f.close() 
        casalog.post(a, origin="ic2233")
    except Exception as instance:
        print("###Error in ic2233 regression: ",instance);
        raise instance

for i in range(1):
     run();
     stats();
