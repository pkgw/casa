import os
import sys
import string
import inspect
from casa_stack_manip import stack_find, find_casa

####---------------- return path to XML files ----------------
def static_var(varname, value):
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate

@static_var("path", None)
def xmlpath( ):
    if xmlpath.path is None:
        __casapath__ = os.environ['CASAPATH'].split(' ')[0]
        __casaarch__ = os.environ['CASAPATH'].split(' ')[1]
        if os.path.exists(__casapath__ + "/" + __casaarch__ + "/xml"):
            xmlpath.path = __casapath__ + "/" + __casaarch__ + "/xml"
        elif os.path.exists(__casapath__ + "/xml"):
            xmlpath.path = __casapath__ + "/xml"
        else:
            raise RuntimeError("Unable to find the XML constraints directory in your CASAPATH")

    return xmlpath.path
####----------------------------------------------------------

casa = find_casa( )

if 'state' in casa and 'init_version' in casa['state'] and casa['state']['init_version'] > 0:

    #
    ##allow globals for taskby default
    casaglobals=True

    casac = stack_find("casac")
    casalog = stack_find("casalog")
    gentools = stack_find("gentools")
    qatool = stack_find("qatool")
    qa = stack_find("qa")
    utilstool = stack_find("utilstool")
    cu = stack_find("cu")
    tbtool =  stack_find("tbtool")
    tb =  stack_find("tb")
    ms =  stack_find("ms")
    mstool =  stack_find("mstool")
    aftool =  stack_find("aftool")
    cbtool =  stack_find("cbtool")
    cltool =  stack_find("cltool")
    write_history = stack_find("write_history")
    me  = stack_find("me")
    metool = stack_find("metool")
    mttool = stack_find("mttool")
    imtool = stack_find("imtool")
    smtool = stack_find("smtool")
    at = stack_find("at")
    attool = stack_find("attool")
    msmdtool = stack_find("msmdtool")

    coordsystool = stack_find("coordsystool")
    cptool = stack_find("cptool")
    cstool = stack_find("cstool")
    dctool = stack_find("dctool")
    fitool = stack_find("fitool")
    fntool = stack_find("fntool")
    iatool = stack_find("iatool")
    imdtool = stack_find("imdtool")
    lmtool = stack_find("lmtool")
    mptool = stack_find("mptool")
    pmtool = stack_find("pmtool")
    pm = stack_find("pm")
    potool = stack_find("potool")
    rgtool = stack_find("rgtool")
    sdmstool = stack_find("sdmstool")
    sdms = stack_find("sdms")
    sltool = stack_find("sltool")
    tptool = stack_find("tptool")
    viewertool = stack_find("viewertool")
    vptool = stack_find("vptool")
    sbstool = stack_find("sbstool")

else:
    from casac import *
    import viewertool
    import os

    def __taskinit_setlogfile( logger ) :
        ####
        #### needed to allow pushing of the global 'casa' state dictionary
        ####
        a=inspect.stack()
        stacklevel=0
        for k in range(len(a)):
            if a[k][1] == "<string>":
                stacklevel=k

        myf=sys._getframe(stacklevel).f_globals

        if 'casa' in myf and 'files' in myf['casa'] and 'logfile' in myf['casa']['files'] :
            logger.setlogfile(myf['casa']['files']['logfile'])


    #
    ##allow globals for taskby default
    casaglobals=True
    # setup available tools
    imager = casac.imager
    imtool=imager
    calibrater = casac.calibrater
    cbtool=calibrater
    mstool = casac.ms
    tptool = casac.tableplot
    tp = tptool()
    mptool = casac.msplot
    mp = mptool()
    pmtool = casac.plotms
    pm = pmtool()
    cptool = casac.calplot
    cp = cptool()
    qatool = casac.quanta
    qa = casac.qa =  qatool()
    tbtool = casac.table
    #fgtool = casac.flagger
    aftool = casac.agentflagger
    af = aftool()
    metool = casac.measures
    iatool = casac.image
    potool = casac.imagepol
    lmtool= casac.linearmosaic
    smtool = casac.simulator
    cltool = casac.componentlist
    coordsystool = casac.coordsys
    cstool = casac.coordsys
    rgtool = casac.regionmanager
    sltool = casac.spectralline
    dctool = casac.deconvolver
    vptool = casac.vpmanager
    msmdtool = casac.msmetadata
    fitool = casac.fitter
    fntool = casac.functional
    imdtool = casac.imagemetadata

    utilstool = casac.utils
    cu = casac.cu = utilstool()
    vftask = casac.vlafillertask()
    vlafiller=vftask.fill
    attool = casac.atmosphere
    ca = casac.calanalysis()
    mttool = casac.mstransformer
    mt = mttool()
    sdmstool = casac.singledishms
    sdms = sdmstool()
    parallelimager = casac.parallelimager()
    sbstool = casac.sidebandseparator
    sbs = sbstool()

    # Log initialization ###################################################################################################

    # IMPORTANT: The following steps must be follow the described order,
    #            otherwise a seg fault occurs when setting the log file.
    # 1st Create casalog object, it will be used by tasks when importing taskinit
    casalog = casac.logsink()
    # 2nd Set log file accessing CASA dictionary of calling context via stack inspection
    __taskinit_setlogfile(casalog)
    # 3rd Set logger as global
    casalog.setglobal(True)

    # Set processor origin (normally "casa" but in the MPI case we use the hostname and rank involved)
    from mpi4casa.MPIEnvironment import MPIEnvironment
    processor_origin = MPIEnvironment.processor_origin
    casalog.processorOrigin(processor_origin)

    # Set showconsole to false for MPIServers
    casalog.showconsole(MPIEnvironment.log_to_console)


    # Log initialization ###################################################################################################

    def gentools(tools=None):
        """
        Generate a fresh set of tools; only the ones who have
        states..so globally sharing the same one can be unpredicatable
        im,cb,ms,tb,me,ia,po,sm,cl,cs,rg,sl,dc,vp,msmd,fi,fn,imd,sdms,lm,at=gentools()
        or if you want specific set of tools
        im, ia, cb=gentools(['im', 'ia', 'cb'])

        """
        tooldic={'im':'imager()', 'cb' :'calibrater()', 'ms':'mstool()',
                 'tb':'tbtool()',  'me' :'metool()',
                 'ia': 'iatool()', 'po':'potool()', 'sm' :'smtool()',
                 'cl': 'cltool()', 'cs' :'cstool()', 'rg':'rgtool()',
                 'sl':'sltool()', 'dc':'dctool()', 'vp':'vptool()',
                 'msmd':'msmdtool()','fi':'fitool()','fn':'fntool()',
                 'imd':'imdtool()','sdms':'sdmstool()', 'lm':'lmtool()',
                 'at':'attool()'}
        reqtools=[]
        if (not tools) or not hasattr(tools, '__iter__'):
            reqtools=['im', 'cb', 'ms','tb', 'me', 'ia', 'po',
                      'sm', 'cl', 'cs', 'rg','sl', 'dc', 'vp',
                      'msmd', 'fi', 'fn', 'imd', 'sdms', 'lm',
                      'at']
        else:
            reqtools=tools
        return tuple([eval(tooldic[reqtool]) for reqtool in reqtools])

    im,cb,ms,tb,me,ia,po,sm,cl,cs,rg,sl,dc,vp,msmd,fi,fn,imd,sdms,lm,at=gentools()

    from mstools import write_history

    ###done with common tools

    # setup viewer tool
    # jagonzal (CAS-4322): Don't load viewer at the engine level
    if 'CASA_ENGINE' not in os.environ:
        try:
            ving = viewertool.viewertool( False )
            if '--nogui' in casa['flags'] :
                vi = ving
            else:
                vi = viewertool.viewertool( True )
        except :
            print("Unable to start viewer, maybe no dbus available?")

    defaultsdir = {}
    defaultsdir['alma'] = casa['dirs']['xml'] + '/almadefaults.xml'
    defaultsdir['evla'] = casa['dirs']['xml'] + '/evladefaults.xml'


    def selectfield(vis,minstring):
        """Derive the fieldid from  minimum matched string(s): """

        tb.open(vis+'/FIELD')
        fields=list(tb.getcol('NAME'))#get fieldname list
        tb.close()              #close table
        indexlist=list()        #initialize list
        stringlist=list()

        fldlist=minstring.split()#split string into elements
        print('fldlist is ',fldlist)
        for fld in fldlist:     #loop over fields
            _iter=fields.__iter__() #create iterator for fieldnames
            while 1:
                try:
                    x=next(_iter) # has first value of field name
                except StopIteration:
                    break
                #
                if (x.find(fld)!=-1):
                    indexlist.append(fields.index(x))
                    stringlist.append(x)

        print('Selected fields are: ',stringlist)
        return indexlist

    def selectantenna(vis,minstring):
        """Derive the antennaid from matched string(s): """

        tb.open(vis+'/ANTENNA')
        ants=list(tb.getcol('NAME'))#get fieldname list
        tb.close()              #close table
        indexlist=list()        #initialize list
        stringlist=list()

        antlist=minstring.split()#split string into elements
        for ant in antlist:     #loop over fields
            try:
                ind=ants.index(ant)
                indexlist.append(ind)
                stringlist.append(ant)
            except ValueError:
                pass

        print('Selected reference antenna: ',stringlist)
        print('indexlist: ',indexlist)
        return indexlist[0]

    def readboxfile(boxfile):
        """ Read a file containing clean boxes (compliant with AIPS BOXFILE)

        Format is:
        #FIELDID BLC-X BLC-Y TRC-X TRC-Y
        0       110   110   150   150
        or
        0       hh:mm:ss.s dd.mm.ss.s hh:mm:ss.s dd.mm.ss.s

        Note all lines beginning with '#' are ignored.

        """
        union=[]
        f=open(boxfile)
        while 1:
            try:
                line=f.readline()
                if (line.find('#')!=0):
                    splitline=line.split('\n')
                    splitline2=splitline[0].split()
                    if (len(splitline2[1])<6):
                        boxlist=[int(splitline2[1]),int(splitline2[2]),
                                 int(splitline2[3]),int(splitline2[4])]
                    else:
                        boxlist=[splitline2[1],splitline2[2],splitline2[3],
                                 splitline2[4]]

                    union.append(boxlist)

            except:
                break

        f.close()
        return union


    # AUTHOR: S. Jaeger
    #
    # NAME: getimaxes
    #
    # DESCRIPTION:
    # This function uses the coordinate information associated
    # with an image to find where the directional (sky) axes are,
    # the spectral axes, and the stokes axes.
    #
    # INPUT:
    #    imagename   string   path to a file on disk.
    #
    # RETURN
    #    list of four lists, [list1, list2, list3, list4 ], as follows :
    #       list1: ['axis num of 1st sky axis', 'Name of axis' ]
    #       list2: ['axis num of 2nd sky axis', 'Name of axis' ]
    #       list3: ['axis num of spectral axis', 'Spectral' ]
    #       list4: ['axis num of stokes axis', 'Stokes' ]

    def getimaxes(imagename):
        """
        Open an image file, looking at its coordinate system information
        to determine which axes are directional, linear, spectral, and
        the stokes axies.

        The return list or lists contains the axis numbers and names in
        the following order:
             1. Directional or Linear
             2. Directional or Linear
             3. Spectral
             4. Stokes

        Note that if an axis type is not found an empty list is returned
        for that axis.
        """

        # Get the images coord. sys.
        csys=None
        ia.open( imagename )
        csys=ia.coordsys()

        # Find where the directional and channel axies are
        # Save the internal placement of the axies in a list
        # (which will be in the following order:
        #    direction1: RA, Longitude, Linear, el, ..
        #    direction2: DEC, Lattitude, Linear, az, ..
        #    spectral:
        #    stokes: I or V
        axisTypes=csys.axiscoordinatetypes()
        axisNames=csys.names()

        # Note that we make a potentially dangerous assumption here
        # that the first directional access is always RA, but it
        # may not be.  The names given to the axies are completely
        # arbitrary and can not be used to determine one axis from
        # another.
        # TODO throw exception??? if we find extra axies or
        #      unrecognized axies.
        retValue=[['',''],['',''],['',''],['','']]
        foundFirstDir=False
        for i in range( len( axisTypes ) ):
            if ( axisTypes[i]=='Direction' or axisTypes[i]=='Linear' ):
                if ( not foundFirstDir ) :
                    retValue[0]=[i,axisNames[i]]
                    foundFirstDir=True
                else:
                    retValue[1]=[i,axisNames[i]]
            elif ( axisTypes[i]=='Spectral' ) :
                retValue[2]=[i,axisNames[i]]
            elif ( axisTypes[i]=='Stokes' ) :
                retValue[3]=[i,axisNames[i]]

        if ( csys != None ):
            del csys
        if ( ia.isopen() ):
            ia.close()
        return retValue


    # from RR
    def announce_async_task(taskname):
        """Use the logger to herald the beginning of an asynchronous task."""
        casalog.origin(taskname)
        casalog.post('')
        casalog.post('###############################################')
        casalog.post('###  Begin Task: ' + taskname + ' ' * (27 - len(taskname)) + '###')
        casalog.post("")
        casalog.post("Use: ")
        casalog.post("      tm.retrieve(return_value) # to retrieve the status")
        casalog.post("")

    def write_task_obit(taskname):
        """Eulogize the task in the logger."""
        casalog.post('###  End Task: ' + taskname + ' ' * (29 - len(taskname)) + '###')
        casalog.post('###############################################')
        casalog.post('')


    def array2string( array ):
        returnValue=""
        for i in range( len(array) ):
            if ( i > 1 ):
                returnValue+=","
            if ( isinstance( array[i], str ) ):
                returnValue+=array[i]
            else:
                returnValue+=str(array[i])
        return returnValue

    def recursivermdir( top='' ):
        # Delete everything from the directory named in 'top',
        # assuming there are no symbolic links.
        for root, dirs, files in os.walk( top, topdown=False ):
            for name in files:
                os.remove( os.path.join( root, name ) )
            for name in dirs:
                os.rmdir( os.path.join( root, name ) )
        os.rmdir(top)
