import os as _os
import re as _re
import sys as _sys
#import time
from casac import casac as _casac

_homedir = _os.getenv('HOME')
if _homedir == None :
   print("Environment variable HOME is not set, please set it")
   _sys.exit(1)

#import casadef

##
## first set up CASAPATH
##
if 'CASAPATH' in _os.environ :
    __casapath__ = _os.environ['CASAPATH'].split(' ')[0]
    if not _os.path.exists(__casapath__ + "/data") :
        raise RuntimeError("CASAPATH environment variable is improperly set")
else :
    __casapath__ = _casac.__file__
    while __casapath__ and __casapath__ != "/" :
        if _os.path.exists( __casapath__ + "/data") :
            break
        __casapath__ = _os.path.dirname(__casapath__)
    if __casapath__ and __casapath__ != "/" :
        _os.environ['CASAPATH']=__casapath__ + " linux local host"
    else :
        raise RuntimeError("CASAPATH environment variable must be set")

#import __casac__
#cu = __casac__.utils.utils()
#xcu = __casac__.utils.utils

# this dict should transition to namespace (as with argparser)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#from casa_system import casa as config

#casa = { 'build': {
#             'time': casadef.build_time,
#             'version': cu.version_info( ),
#             'number': casadef.subversion_revision
#         },
#         'source': {
#             'url': casadef.subversion_url,
#             'revision': casadef.subversion_revision
#         },
#         'helpers': {
#             'logger': 'casalogger',
#             'viewer': 'casaviewer',
#             'info': None,
#             'dbus': None,
#             'ipcontroller': None,
#             'ipengine': None
#         },
#         'dirs': {
#             'rc': _homedir + '/.casa',
#             'data': __casapath__ + "/data",
#             'recipes': __casapath__ + "/lib/python2.7/recipes",
#             'root': __casapath__,
#             'python':  __casapath__ + "/lib/python2.7",
#             'pipeline': None,
#             'xml': __casapath__ + "/xml"
#         },
#         'flags': { },
#         'files': { 
#             'logfile': _os.getcwd( ) + '/casa-'+time.strftime("%Y%m%d-%H%M%S", time.gmtime())+'.log'
#         },
#         'state' : {
#             'startup': True,
#             'unwritable': set( )
#         }
#       }

##
## next adjust the PYTHONPATH
##
def _adapt_pythonpath(searchroot):
    # tarball location
    guess = _os.path.join(searchroot, 'lib/python2.7/site-packages/numpy')
    if _os.path.isdir(guess):
        _sys.path.append(_os.path.dirname(guess))
    else:
        for root, dirs, files in _os.walk(searchroot):
            # skip data folder which might be a network mount
            if root == searchroot and 'data' in dirs:
                del dirs[dirs.index('data')]
            if root.endswith("/numpy"):
                _sys.path.append(_os.path.dirname(root))
                break

if _re.match( r'.*/\d+\.\d+\.\d+\w*-\d+$', __casapath__ ) :
    _adapt_pythonpath(_os.path.dirname(__casapath__))
else:
    _adapt_pythonpath(__casapath__)

##
## next adjust PATH and LD_LIBRARY_PATH
##
def _setup_path():
    global __ipcontroller__, __ld_library_path__
    _rootdir = None
    if _os.path.exists(_os.path.join(__casapath__, 'bin', 'casapyinfo')):
        _rootdir = _os.path.join(__casapath__, 'bin')
    else:
        for root, dirs, files in _os.walk(__casapath__):
            # skip data folder which might be a network mount
            if root == __casapath__ and 'data' in dirs:
                del dirs[dirs.index('data')]
            if root.endswith("/bin") and "casapyinfo" in files :
                _rootdir = root
                break
    if _rootdir is None:
        return

    __ipcontroller__ = (lambda fd: fd.readline().strip('\n'))(_os.popen(_rootdir + "/casapyinfo --exec 'which ipcontroller'"))
    if _os.path.exists(__ipcontroller__) :
        _os.environ['PATH'] = _os.path.dirname(__ipcontroller__) + ":" + _os.environ['PATH']
    else :
        raise RuntimeError("cannot configure CASA tasking system")
    __ld_library_path__ = (lambda fd: fd.readline().strip('\n').split(':'))(_os.popen(_rootdir + "/casapyinfo --exec 'echo $LD_LIBRARY_PATH'"))
    for x in __ld_library_path__: _sys.path.append(x)

_setup_path()

##
## finally load tools
##
imager = _casac.imager
imtool=imager
calibrater = _casac.calibrater
cbtool=calibrater
mstool = _casac.ms
tptool = _casac.tableplot
mptool = _casac.msplot
pmtool = _casac.plotms
cptool = _casac.calplot
qatool = _casac.quanta
tbtool = _casac.table
aftool = _casac.agentflagger
metool = _casac.measures
iatool = _casac.image
potool = _casac.imagepol
lmtool= _casac.linearmosaic
sbstool = _casac.sidebandseparator
smtool = _casac.simulator
cltool = _casac.componentlist
coordsystool = _casac.coordsys
cstool = _casac.coordsys
rgtool = _casac.regionmanager
sltool = _casac.spectralline
dctool = _casac.deconvolver
vptool = _casac.vpmanager
msmdtool = _casac.msmetadata
fitool = _casac.fitter
fntool = _casac.functional
imdtool = _casac.imagemetadata

cutool = _casac.utils
mttool = _casac.mstransformer
sdmstool = _casac.singledishms
catool = _casac.calanalysis
attool = _casac.atmosphere

from accum import accum
from applycal import applycal
from asdmsummary import asdmsummary
from bandpass import bandpass
from blcal import  blcal
from browsetable import  browsetable
from calstat import  calstat
from caltabconvert import  caltabconvert
from clean import  clean
from clearcal import  clearcal
from clearstat import  clearstat
from concat import  concat
from conjugatevis import  conjugatevis
from cvel import  cvel
from cvel2 import  cvel2
from deconvolve import  deconvolve
from delmod import  delmod
from exportasdm import  exportasdm
from exportfits import  exportfits
from exportuvfits import  exportuvfits
from feather import  feather
from fixplanets import  fixplanets
from fixvis import  fixvis
from flagcmd import  flagcmd
from flagdata import  flagdata
from flagmanager import  flagmanager
from fluxscale import  fluxscale
from fringefit import  fringefit
from ft import  ft
from gaincal import  gaincal
from gencal import  gencal
from hanningsmooth import  hanningsmooth
from imcollapse import  imcollapse
from imcontsub import  imcontsub
from imdev import  imdev
from imfit import  imfit
from imhead import  imhead
from imhistory import  imhistory
from immath import  immath
from immoments import  immoments
from impbcor import  impbcor
from importatca import  importatca
from importasap import  importasap
from importasdm import  importasdm
from importfits import  importfits
from importfitsidi import  importfitsidi
from importgmrt import  importgmrt
from importmiriad import  importmiriad
from importnro import  importnro
from importuvfits import  importuvfits
from importvla import  importvla
from imrebin import  imrebin
from imreframe import  imreframe
from imregrid import  imregrid
from imsmooth import  imsmooth
from imstat import  imstat
from imsubimage import  imsubimage
from imtrans import  imtrans
from imval import  imval
from initweights import  initweights
from listcal import  listcal
from listhistory import  listhistory
from listfits import  listfits
from listobs import  listobs
from listpartition import  listpartition
from listsdm import  listsdm
from listvis import  listvis
from makemask import  makemask
from mstransform import  mstransform
from msuvbin import  msuvbin
from oldsplit import  oldsplit
from plotants import  plotants
from plotbandpass import  plotbandpass
from plotcal import  plotcal
from plotms import  plotms
from plotweather import  plotweather
from plotprofilemap import  plotprofilemap
from partition import  partition
from polcal import  polcal
from predictcomp import  predictcomp
from impv import  impv
from rmfit import  rmfit
from rmtables import  rmtables
from sdbaseline import  sdbaseline
from sdcal import  sdcal
from sdfit import  sdfit
from sdfixscan import  sdfixscan
from sdgaincal import  sdgaincal
from sdimaging import  sdimaging
from sdsidebandsplit import  sdsidebandsplit
from sdsmooth import  sdsmooth
from setjy import  setjy
from simalma import  simalma
from simobserve import  simobserve
from simanalyze import  simanalyze
from slsearch import  slsearch
from smoothcal import  smoothcal
from specfit import  specfit
from specflux import  specflux
from specsmooth import  specsmooth
from splattotable import  splattotable
from split import  split
#from split2 import split2
from spxfit import  spxfit
from oldstatwt import  oldstatwt
from statwt import  statwt
from tclean import  tclean
from tclean2 import  tclean2
from testconcat import  testconcat
from uvcontsub import  uvcontsub
from uvcontsub3 import  uvcontsub3
from uvmodelfit import  uvmodelfit
from uvsub import  uvsub
from wvrgcal import  wvrgcal
from virtualconcat import  virtualconcat
from vishead import  vishead
from visstat import  visstat
from widebandpbcor import  widebandpbcor
