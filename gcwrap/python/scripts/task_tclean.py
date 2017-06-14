################################################
# Refactored Clean task
#
# v1.0: 2012.10.05, U.R.V.
#
################################################

from taskinit import *

import os
import shutil
import numpy
from taskinit import *
import copy
import time;

#from refimagerhelper import PySynthesisImager
#from refimagerhelper import PyParallelContSynthesisImager,PyParallelCubeSynthesisImager
#from refimagerhelper import ImagerParameters

from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
from imagerhelpers.input_parameters import ImagerParameters



def tclean(
    ####### Data Selection
    vis,#='', 
    selectdata,
    field,#='', 
    spw,#='',
    timerange,#='',
    uvrange,#='',
    antenna,#='',
    scan,#='',
    observation,#='',
    intent,#='',
    datacolumn,#='corrected',

    ####### Image definition
    imagename,#='',
    imsize,#=[100,100],
    cell,#=['1.0arcsec','1.0arcsec'],
    phasecenter,#='J2000 19:59:28.500 +40.44.01.50',
    stokes,#='I',
    projection,#='SIN',
    startmodel,#='',

    ## Spectral parameters
    specmode,#='mfs',
    reffreq,#='',
    nchan,#=1,
    start,#='',
    width,#='',
    outframe,#='LSRK',
    veltype,#='',
    restfreq,#=[''],
#    sysvel,#='',
#    sysvelframe,#='',
    interpolation,#='',
    ## 
    ####### Gridding parameters
    gridder,#='ft', 
    facets,#=1,
    chanchunks,#=1,

    wprojplanes,#=1,

    ### PB
    vptable,

    aterm,#=True,
    psterm,#=True,
    wbawp ,#= True,
    conjbeams ,#= True,
    cfcache ,#= "",
    computepastep ,#=360.0,
    rotatepastep ,#=360.0,

    pblimit,#=0.01,
    normtype,#='flatnoise',

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,
    smallscalebias,#=0.6

    ### restoration options
    restoration,
    restoringbeam,#=[],
    pbcor,

    ##### Outliers
    outlierfile,#='',

    ##### Weighting
    weighting,#='natural',
    robust,#=0.5,
    npixels,#=0,
#    uvtaper,#=False,
    uvtaper,#=[],


    ##### Iteration control
    niter,#=0, 
    gain,#=0.1,
    threshold,#=0.0, 
    cycleniter,#=0, 
    cyclefactor,#=1.0,
    minpsffraction,#=0.1,
    maxpsffraction,#=0.8,
    interactive,#=False, 

    ##### (new) Mask parameters
    usemask,#='user',
    mask,#='',
    pbmask,#='',
    maskthreshold,#='',
    maskresolution,#='',
    nmask,#=0,

    ##### automask by multithresh
    sidelobethreshold,#=5.0,
    noisethreshold,#=3.0,
    lownoisethreshold,#=3.0,
    smoothfactor,#=1.0,
    minbeamfrac,#=0.3, 
    cutthreshold,#=0.01,

    ## Misc

    restart,#=True,

    savemodel,#="none",

#    makeimages,#="auto"
    calcres,#=True,
    calcpsf,#=True,

    ####### State parameters
    parallel):#=False):

    #####################################################
    #### Sanity checks and controls
    #####################################################
    
    ### Move these checks elsewhere ? 

    if specmode=='cont':
        specmode='mfs'

    if specmode=='mfs' and nterms==1 and deconvolver == "mtmfs":
        casalog.post( "The MTMFS deconvolution algorithm (deconvolver='mtmfs') needs nterms>1.Please set nterms=2 (or more). ", "WARN", "task_tclean" )
        return

    if specmode!='mfs' and deconvolver=="mtmfs":
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') applies only to specmode='mfs'.", "WARN", "task_tclean" )
        return

    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None

    # Put all parameters into dictionaries and check them. 
    paramList = ImagerParameters(
        msname =vis,
        field=field,
        spw=spw,
        timestr=timerange,
        uvdist=uvrange,
        antenna=antenna,
        scan=scan,
        obs=observation,
        state=intent,
        datacolumn=datacolumn,

        ### Image....
        imagename=imagename,
        #### Direction Image Coords
        imsize=imsize, 
        cell=cell, 
        phasecenter=phasecenter,
        stokes=stokes,
        projection=projection,
        startmodel=startmodel,

        ### Spectral Image Coords
        specmode=specmode,
        reffreq=reffreq,
        nchan=nchan,
        start=start,
        width=width,
        outframe=outframe,
        veltype=veltype,
        restfreq=restfreq,
        sysvel='', #sysvel,
        sysvelframe='', #sysvelframe,
        interpolation=interpolation,

        gridder=gridder,
#        ftmachine=ftmachine,
        facets=facets,
        chanchunks=chanchunks,

        wprojplanes=wprojplanes,

        vptable=vptable,

        ### Gridding....

        aterm=aterm,
        psterm=psterm,
        wbawp = wbawp,
        cfcache = cfcache,
        conjbeams = conjbeams,
        computepastep =computepastep,
        rotatepastep = rotatepastep,

        pblimit=pblimit,
        normtype=normtype,

        outlierfile=outlierfile,
        restart=restart,

        weighting=weighting,
        robust=robust,
        npixels=npixels,
        uvtaper=uvtaper,

        ### Deconvolution
        niter=niter,
        cycleniter=cycleniter,
        loopgain=gain,
        threshold=threshold,
        cyclefactor=cyclefactor,
        minpsffraction=minpsffraction, 
        maxpsffraction=maxpsffraction,
        interactive=interactive,

        deconvolver=deconvolver,
        scales=scales,
        nterms=nterms,
        scalebias=smallscalebias,
        restoringbeam=restoringbeam,
        
        ### new mask params
        usemask=usemask,
        mask=mask,
        pbmask=pbmask,
        maskthreshold=maskthreshold,
        maskresolution=maskresolution,
        nmask=nmask,

        ### automask multithresh params
        sidelobethreshold=sidelobethreshold,
        noisethreshold=noisethreshold,
        lownoisethreshold=lownoisethreshold,
        smoothfactor=smoothfactor,
        minbeamfrac=minbeamfrac,
        cutthreshold=cutthreshold,
 
        savemodel=savemodel
        )
    
    #paramList.printParameters()

    pcube=False
    concattype=''
    if parallel==True and specmode!='mfs':
        pcube=True
        parallel=False

    # catch non operational case (parallel cube tclean with interative=T)
    if pcube and interactive:
        casalog.post( "Interactive mode is not currently supported with parallel cube CLEANing, please restart by setting interactive=F", "WARN", "task_tclean" )
        return False
   

    ## Setup Imager objects, for different parallelization schemes.
    if parallel==False and pcube==False:
   
         imager = PySynthesisImager(params=paramList)
    elif parallel==True:
         imager = PyParallelContSynthesisImager(params=paramList)
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params=paramList)
         # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
         concattype='virtualcopy'
    else:
         print 'Invalid parallel combination in doClean.'
         return False
    
    retrec={}

    try: 
    #if (1):
        ## Init major cycle elements
        t0=time.time();
        imager.initializeImagers()
    
        # Construct the CFCache for AWProject-class of FTMs.  For
        # other choices the following three calls become NoOps.
        # imager.dryGridding();
        # imager.fillCFCache();
        # imager.reloadCFCache();

        imager.initializeNormalizers()
        imager.setWeighting()
        t1=time.time();
        casalog.post("***Time for initializing imager and normalizers: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ## Init minor cycle elements
        if niter>0 or restoration==True:
            t0=time.time();
            imager.initializeDeconvolvers()
            t1=time.time();
            casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        if niter>0:
            t0=time.time();
            imager.initializeIterationControl()
            t1=time.time();
            casalog.post("***Time for initializing iteration controller: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
            
        ## Make PSF
        if calcpsf==True:
            t0=time.time();

            imager.makePSF()

            t1=time.time();
            casalog.post("***Time for making PSF: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

            imager.makePB()

            t2=time.time();
            casalog.post("***Time for making PB: "+"%.2f"%(t2-t1)+" sec", "INFO3", "task_tclean");

        if niter >=0 : 

            ## Make dirty image
            if calcres==True:
                t0=time.time();
                imager.runMajorCycle()
                t1=time.time();
                casalog.post("***Time for major cycle (calcres=T): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

            ## In case of no deconvolution iterations....
            if niter==0 and calcres==False:
                if savemodel != "none":
                    imager.predictModel()

            ## Do deconvolution and iterations
            if niter>0 :

                isit = imager.hasConverged()
                imager.updateMask()

                while ( not imager.hasConverged() ):

#                    maskchanged = imager.updateMask()
#                    if maskchanged and imager.hasConverged() :
#                        break;

                    t0=time.time();
                    imager.runMinorCycle()
                    t1=time.time();
                    casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

                    t0=time.time();
                    imager.runMajorCycle()
                    t1=time.time();
                    casalog.post("***Time for major cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

                    imager.updateMask()

                ## Get summary from iterbot
                if type(interactive) != bool:
                    retrec=imager.getSummary();

            ## Restore images.
            if restoration==True:  
                t0=time.time();
                imager.restoreImages()
                t1=time.time();
                casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                if pbcor==True:
                    t0=time.time();
                    imager.pbcorImages()
                    t1=time.time();
                    casalog.post("***Time for pb-correcting images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                    


        if (pcube):
            print "running concatImages ..."
            casalog.post("Running virtualconcat (type=%s) of sub-cubes" % concattype,"INFO2", "task_tclean")
            # fixed to move subcubes
            imager.concatImages(type=concattype)

        ## Close tools.
        imager.deleteTools()

    except Exception as e:
        #print 'Exception : ' + str(e)
        casalog.post('Exception from task_tclean : ' + str(e), "SEVERE", "task_tclean")
        if imager != None:
            imager.deleteTools() 

        larg = list(e.args)
        larg[0] = 'Exception from task_tclean : ' + str(larg[0])
        e.args = tuple(larg)
        raise

    return retrec

##################################################
