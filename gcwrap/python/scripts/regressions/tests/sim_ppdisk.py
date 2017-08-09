# simulation tasks regression
#
#  single pointing
#  single interferometric measurement
#  no noise

import os, time, string, sys, inspect
from locatescript import locatescript

a=inspect.stack()
stacklevel=0
for k in range(len(a)):
    if (string.find(a[k][1], 'ipython console') > 0):
        stacklevel=k
        break
gl=sys._getframe(stacklevel).f_globals


# Short description
def description():
    return "Simulates a single 12m ALMA pointing from a model image. noise, imaged."


def data():
    ### return the data files that is needed by the regression script
    return []


def run(fetch=False):
    #####locate the regression script
    lepath=locatescript('ppdisk2_regression.py')
    print('Script used is ',lepath)
    gl['regstate']=True
    exec(compile(open(lepath).read(), lepath, 'exec'), gl)
    print('regstate =', gl['regstate'])
    if not gl['regstate']:
        raise Exception('regstate = False')
###return the images that will be templated and compared in future runs
    return ['psim2/psim2.alma.out20.noisy.ms','psim2/psim2.alma.out20.noisy.image','psim2/psim2.alma.out20.noisy.diff']

