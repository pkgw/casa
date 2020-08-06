try:
    from casac import casac
except ImportError as e:
    print("failed to load casa:\n", e)
    sys.exit(1)

from casa_system import casa

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
			  'msmd', 'fi', 'fn', 'imd', 'sdms', 'lm', 'at']
	else:
		reqtools=tools
	return tuple([eval(tooldic[reqtool]) for reqtool in reqtools])

from mstools import write_history

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
cutool = utilstool
cu = casac.cu = cutool( )
vftask = casac.vlafillertask()
vlafiller=vftask.fill
attool = casac.atmosphere
catool = casac.calanalysis
ca = catool( )
mttool = casac.mstransformer
mt = mttool()
sdmstool = casac.singledishms
sdms = sdmstool()
parallelimager = casac.parallelimager()
sbstool = casac.sidebandseparator
sbs = sbstool()


##
## viewer tool
##
from viewertool import viewertool
try:
    ving = viewertool( False )
    if casa['flags'].nogui :
        vi = ving
    else:
        vi = viewertool( True )
except:
    print("Unable to start viewer, maybe no dbus available?")


im,cb,ms,tb,me,ia,po,sm,cl,cs,rg,sl,dc,vp,msmd,fi,fn,imd,sdms,lm,at=gentools()
