import os
import sys
import time
import casadef
from procmgr import procmgr as _procmgr

import __casac__
_cu = __casac__.utils.utils()

from casa_system_defaults import casa
casa['version'] = _cu.version_string( )
casa['variant'] = _cu.version_variant( )
casa['version_desc'] = _cu.version_desc( )

###########################################################################
######## these will be removed at the CASA 5.1 code freeze time... ########
casa['build'] = {  'time': casadef.build_time,
                   'version': _cu.version_info( ),
                   'number': casadef.subversion_revision }
casa['source'] = { 'url': casadef.subversion_url,
                   'revision': casadef.subversion_revision }
###########################################################################

procmgr = _procmgr( )
