########################################################################3
#  task_immoments.py
#
#
# Copyright (C) 2008, 2009
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
# <author>
# Shannon Jaeger (University of Calgary)
# </author>
#
# <summary>
# CASA task for finding moments along a specified axis of a
# multi-dimentional CASA image.
# contents
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed

# <etymology>
# immoments stands for image momemnts
# </etymology>
#
# <synopsis>
# task_immoments.py is a Python script providing an easy to use task
# for generating momements along a user specified axis on a CASA image.
# This is a time-honoured spectral-line analysis technique for discovering
# spectral line information.
#
# In this task, moment, refers to collapsing an axis of the image,
# the moment axis, to a single pixel.
#
# The various moments that can be calculated are described in detail
# at http://casa.nrao.edu/docs/casaref/image.moments.html#x59-590001.1.1
#
# </synopsis>
#
# <example>
# <srcblock>
# # The following code snippet find the 1-moments, intensity-weighted
# # coordinate, often used for finding velocity fields.
# immoments( 'axis='spec', imagename='myimage', moment=1, outfile='velocityfields' )
#
# # Finding the spectral mean, -1 moment, on a specified portion of the image.
# # The region used is defined by the box and stokes parameters
# immoments( imagename='myimage', axis='spec', stokes='I', box=[55,12,97,32], moment=-1 )
#
# # The following example uses a second file to use as a mask.
# # The 0-moments, integrated values, are created on clean.im, but
# # the mask is based on all the data in  calibrated.im, all values
# # about the threshold 0.5 are used to create the mask.
#
# immoments( 'clean.image', axis='spec', mask='calibrated.im>0.5', outfile='mom_withmask.im' )
# </srblock>
# </example>
#
# <motivation>
# To provide a user-friendly method of calculating image moments.
# </motivation>
#
# <todo>
# </todo>

import os
from taskinit import *
_rg = rgtool( )

def immoments(
    imagename, moments, axis, region, box, chans, stokes,
    mask, includepix, excludepix, outfile, stretch
):
    retValue=None
    casalog.origin('immoments')
    _myia = iatool()
    try:
        # First check to see if the output file exists.  If it
        # does then we abort.  CASA doesn't allow files to be
        # over-written, just a policy.
        if ( len( outfile ) > 0 and os.path.exists( outfile ) ):
            raise Exception('Output file, '+outfile+\
              ' exists. immoment can not proceed, please\n'\
              'remove it or change the output file name.')
        elif ( len( outfile ) ==  1 ):
            raise Exception('outfile is not specified but must be')
        _myia.open(imagename)
        reg = _rg.frombcs(csys=_myia.coordsys().torecord(),
            shape=_myia.shape(), box=box, chans=chans, stokes=stokes,
            stokescontrol="a", region=region
        )
        if isinstance( axis, str ):
             axis = _myia.coordsys().findaxisbyname(axis)
        retValue = _myia.moments(
            moments=moments, axis=int(axis), mask=mask,
            region=reg, includepix=includepix,
            excludepix=excludepix, outfile=outfile, drop=False,
            stretch=stretch
        )
        retValue.done()
        return True
    except Exception as x:
        _myia.done()
        print('*** Error ***: ' + str(x))
        raise
    finally:
        _myia.done()

