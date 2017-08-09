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
    return "PdB, flagging, calibration, imaging"

def data():
    ### return the data files that is needed by the regression script
    return []

def run(fetch=False):
    lepath=locatescript('l02d_regression.py')
    print('Script used is ',lepath)
    gl['regstate']=True
    exec(compile(open(lepath).read(), lepath, 'exec'), gl)
    print('regstate =', gl['regstate'])
    if not gl['regstate']:
        raise Exception('regstate = False')

###return the images that will be templated and compared in future runs
    return ['l2.3mmcont.image', 'l2.1mm.image',  'l02d.3mmcont.image', 'l02d.1mmcont.image', 'l02d.c34s.image', 'l02d.ch3oh.image' ]
