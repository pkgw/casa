import os
from taskinit import *

def ssoflux(vis=None,field=None,spw=None,modimage=None,fluxdensity=None,standard=None):
       """
       *This is an experimental clone of setjy while flux calibration with
       Solar System objects is being tested.  It will eventually be merged
       back into setjy.*

       Fills the model column for flux density calibrators:

       The task places the model visibility amp and phase associated
       with a specified clean components image into the model column
       of the data set.  The simplest way is to enter the flux density
       (I,Q,U,V) explicitly.

       Models are available for 3C48, 3C138, and 3C286 between
       1.4 and 43 GHz.  3C147 is available above 13 GHz.  These models
       are scaled to the precise frequency of the data.  Only I models are
       presently available.

       The location of the models are system dependent:  At the AOC, the
       models are in the directory::/usr/lib/casapy/data/nrao/VLA/CalModels

       ssoflux need only be run on the calibrator sources with a know flux
       density and/or model.

       Keyword arguments:
       vis -- Name of input visibility file
               default: none.  example: vis='ngc5921.ms'
       field -- Select field using field id(s) or field name(s).
                 [run listobs to obtain the list id's or names of calibrators]
              default: ''=all fields
              Only one field can be used if a source model or specific fluxdensity
                 is inserted.
              If field string is a non-negative integer, it is assumed a field index
                otherwise, it is assumed a field name
              field='0~2'; field ids 0,1,2
              field='0,4,5~7'; field ids 0,4,5,6,7
              field='3C286,3C295'; fields named 3C286 and 3C295
              field = '3,4C*'; field id 3, all names starting with 4C
       spw -- Select spectral window/channels
              default: ''=all spectral windows and channels
              spw='0~2,4'; spectral windows 0,1,2,4 (all channels)
              spw='<2';  spectral windows less than 2 (i.e. 0,1)
              spw='0:5~61'; spw 0, channels 5 to 61
              spw='0,10,3:3~45'; spw 0,10 all channels, spw 3, channels 3 to 45.
              spw='0~2:2~6'; spw 0,1,2 with channels 2 through 6 in each.
              spw='0:0~10;15~60'; spectral window 0 with channels 0-10,15-60
              spw='0:0~10,1:20~30,2:1;2;4'; spw 0, channels 0-10,
                       spw 1, channels 20-30, and spw 2, channels, 1,2 and 4
       modimage -- Optional model image (I only) from which the model visibilities
              are calculated and placed in the model column.  Each field must
              be done separately.  The image flux density will be scaled from the
              frequency in the model to that actually used.
              In CV and the AOC, the models are located in the
              /usr/lib/casapy/data/nrao/VLA/CalModels directory
       fluxdensity -- Specified flux density [I,Q,U,V] in Jy
               default => -1, puts total flux density for known standard source
               places [1,0,0,0] for any other source.
               Otherwise,  places the specified flux density for the source.  This
               is the only way to insert a polarized flux density model at the
               present time (March 2008).
               example:  fluxdensity=[2.63,0.21,-0.33,0.02]
       standard -- Flux density standard, used if fluxdensity<0
               default: 'Perley-Taylor 99'; example: standard='Baars'
               Options: 'Baars','Perley 90','Perley-Taylor 95',
                  'Perley-Taylor 99'

       """

       try:
         casalog.origin('ssoflux')
         casalog.post('ssoflux is no longer supported.  Use setjy.', 'SEVERE')

         ## if ((type(vis)==str) & (os.path.exists(vis))):
         ##             im.open(vis, usescratch=True)
         ## else:
         ##             raise Exception, 'Visibility data set not found - please verify the name'

         ## im.ssoflux(field=field,spw=spw,modimage=modimage,fluxdensity=fluxdensity,standard=standard)
         ## im.close()


         ##      #write history
         ## ms.open(vis,nomodify=False)
         ## ms.writehistory(message='taskname = ssoflux',origin='ssoflux')
         ## ms.writehistory(message='vis         = "'+str(vis)+'"',origin='ssoflux')
         ## ms.writehistory(message='field       = "'+str(field)+'"',origin='ssoflux')
         ## ms.writehistory(message='spw       = '+str(spw),origin='ssoflux')
         ## ms.writehistory(message='modimage = '+str(modimage),origin='ssoflux')
         ## ms.writehistory(message='fluxdensity = '+str(fluxdensity),origin='ssoflux')
         ## ms.writehistory(message='standard    = "'+str(standard)+'"',origin='ssoflux')
         ## ms.close()

       except Exception as instance:
              print('*** Error ***',instance)

