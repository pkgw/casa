import sys
import os
import string
from locatescript import locatescript
import inspect

a=inspect.stack()
stacklevel=0
for k in range(len(a)):
    if (string.find(a[k][1], 'ipython console') > 0):
        stacklevel=k
        break
gl=sys._getframe(stacklevel).f_globals

def description():
    return "NGC 5921, VLA D-array, import, export, flagging, calibration, imaging, image statistics"

def run():
    lepath=locatescript('ngc5921_uvcontsub2.py')
    print('Script used is ',lepath)
    gl['regstate']=True
    exec(compile(open(lepath).read(), lepath, 'exec'), gl)
    print('regstate =', gl['regstate'])
    if not gl['regstate']:
        raise Exception('regstate = False')

###return the images that will be templated and compared in future runs
    return ['ngc5921uvcontsub2/ngc5921.clean.image']

def data():
    ### return the data files that is needed by the regression script
    return []
