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
import pdb
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
    perchanweightdensity, #=''
    ## 
    ####### Gridding parameters
    gridder,#='ft',
    facets,#=1,
    psfphasecenter, #=''
    chanchunks,#=1,

    wprojplanes,#=1,

    ### PB
    vptable,
    mosweight, #=True
    aterm,#=True,
    psterm,#=True,
    wbawp ,#= True,
    conjbeams ,#= True,
    cfcache ,#= "",
    usepointing, #=false
    computepastep ,#=360.0,
    rotatepastep ,#=360.0,
    pointingoffsetsigdev ,#=0.0,

    pblimit,#=0.01,
    normtype,#='flatnoise',

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,
    smallscalebias,#=0.0

    ### restoration options
    restoration,
    restoringbeam,#=[],
    pbcor,

    ##### Outliers
    outlierfile,#='',

    ##### Weighting
    weighting,#='natural',
    robust,#=0.5,
    noise,#0.0Jy
    npixels,#=0,
#    uvtaper,#=False,
    uvtaper,#=[],


    ##### Iteration control
    niter,#=0, 
    gain,#=0.1,
    threshold,#=0.0, 
    nsigma,#=0.0
    cycleniter,#=0, 
    cyclefactor,#=1.0,
    minpsffraction,#=0.1,
    maxpsffraction,#=0.8,
    interactive,#=False, 

    ##### (new) Mask parameters
    usemask,#='user',
    mask,#='',
    pbmask,#='',
    # maskthreshold,#='',
    # maskresolution,#='',
    # nmask,#=0,

    ##### automask by multithresh
    sidelobethreshold,#=5.0,
    noisethreshold,#=3.0,
    lownoisethreshold,#=3.0,
    negativethreshold,#=0.0,
    smoothfactor,#=1.0,
    minbeamfrac,#=0.3, 
    cutthreshold,#=0.01,
    growiterations,#=100
    dogrowprune,#=True
    minpercentchange,#=0.0
    verbose, #=False
    fastnoise, #=False

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
    inpparams=locals().copy()
    ###now deal with parameters which are not the same name 
    inpparams['msname']= inpparams.pop('vis')
    inpparams['timestr']= inpparams.pop('timerange')
    inpparams['uvdist']= inpparams.pop('uvrange')
    inpparams['obs']= inpparams.pop('observation')
    inpparams['state']= inpparams.pop('intent')
    inpparams['loopgain']=inpparams.pop('gain')
    inpparams['scalebias']=inpparams.pop('smallscalebias')
    if specmode=='cont':
        specmode='mfs'
        inpparams['specmode']='mfs'
#    if specmode=='mfs' and nterms==1 and deconvolver == "mtmfs":
#        casalog.post( "The MTMFS deconvolution algorithm (deconvolver='mtmfs') needs nterms>1.Please set nterms=2 (or more). ", "WARN", "task_tclean" )
#        return

    if (deconvolver=="mtmfs") and (specmode!='mfs') and (specmode!='cube' or nterms!=1) and (specmode!='cubedata' or nterms!=1):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') applies only to specmode='mfs' or specmode='cube' with nterms=1 or specmode='cubedata' with nterms=1.", "WARN", "task_tclean" )
        return
      
    if(deconvolver=="mtmfs" and (specmode=='cube' or specmode=='cubedata') and nterms==1 and parallel==True):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') with specmode='cube', nterms=1 currently only works in serial.", "WARN", "task_tclean" )
        return

    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None
    # Put all parameters into dictionaries and check them.
    ##make a dictionary of parameters that ImagerParameters take
    defparm=dict(list(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
    ###assign values to the ones passed to tclean and if not defined yet in tclean...
    ###assign them the default value of the constructor
    bparm={k:  inpparams[k] if k in inpparams else defparm[k]  for k in list(defparm.keys())}
    ###default mosweight=True is tripping other gridders as they are not
    ###expecting it to be true
    if(bparm['mosweight']==True and bparm['gridder'].find("mosaic") == -1):
        bparm['mosweight']=False
    paramList=ImagerParameters(**bparm)

    # deprecation message
    if usemask=='auto-thresh' or usemask=='auto-thresh2':
        casalog.post(usemask+" is deprecated, will be removed in CASA 5.4.  It is recommended to use auto-multithresh instead", "WARN") 

    #paramList.printParameters()
    
    if pointingoffsetsigdev!=0.0 and usepointing==False:
        casalog.post("pointingoffsetsigdev is only revelent when usepointing is True", "WARN") 

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
    imagerInst=PySynthesisImager
    if parallel==False and pcube==False:
         imager = PySynthesisImager(params=paramList)
         imagerInst=PySynthesisImager
    elif parallel==True:
         imager = PyParallelContSynthesisImager(params=paramList)
         imagerInst=PyParallelContSynthesisImager
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params=paramList)
         imagerInst=PyParallelCubeSynthesisImager
         # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
         concattype='virtualcopy'
    else:
         print('Invalid parallel combination in doClean.')
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

        ####now is the time to check estimated memory
        imager.estimatememory()
            
        if niter>0:
            t0=time.time();
            imager.initializeIterationControl()
            t1=time.time();
            casalog.post("***Time for initializing iteration controller: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
            
        ## Make PSF
        if calcpsf==True:
            t0=time.time();
             
            imager.makePSF()
            if((psfphasecenter != '') and (gridder=='mosaic')):
                print("doing with different phasecenter psf")
                imager.unlockimages(0)
                psfParameters=paramList.getAllPars()
                psfParameters['phasecenter']=psfphasecenter
                psfParamList=ImagerParameters(**psfParameters)
                psfimager=imagerInst(params=psfParamList)
                psfimager.initializeImagers()
                psfimager.setWeighting()
                psfimager.makeImage('psf', psfParameters['imagename']+'.psf')
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
            print("running concatImages ...")
            casalog.post("Running virtualconcat (type=%s) of sub-cubes" % concattype,"INFO2", "task_tclean")
            # fixed to move subcubes
            imager.concatImages(type=concattype)

        ## Close tools.
        imager.deleteTools()

        # CAS-10721 
        if niter>0 and savemodel != "none":
            casalog.post("Please check the casa log file for a message confirming that the model was saved after the last major cycle. If it doesn't exist, please re-run tclean with niter=0,calcres=False,calcpsf=False in order to trigger a 'predict model' step that obeys the savemodel parameter.","WARN","task_tclean")


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
