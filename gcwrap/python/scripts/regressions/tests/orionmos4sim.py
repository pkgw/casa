import sys
import os
import string
from locatescript import copydata
from locatescript import locatescript
import inspect

a=inspect.stack()
stacklevel=0
for k in range(len(a)):
    if (string.find(a[k][1], 'ipython console') > 0):
        stacklevel=k
        break
gl=sys._getframe(stacklevel).f_globals

# Short description
def description():
    return "mosaicing"

def data():
    ### return the data files that is needed by the regression script
    return ['run_orionmos4sim.pickle', 'orion.ms']

def run(fetch=False):

    #####fetch data
    if fetch:
        for f in data( ):
            copydata( f, os.getcwd( ) )
    
    #####locate the regression script
    lepath=locatescript('run_orionmos4sim.py')
    print('Script used is ',lepath)
    gl['regstate']=True
    exec(compile(open(lepath).read(), lepath, 'exec'), gl)
    print('regstate =', gl['regstate'])
    if not gl['regstate']:
        raise Exception('regstate = False')
###return the images that will be templated and compared in future runs
    return []
