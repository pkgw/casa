import sys
import os
import string
import inspect
import pylab as pl

gl=sys._getframe(len(inspect.stack())-1).f_globals
def run(fetch=False):
    x=list(range(10000000))
    pl.ion()
    for k in range(10):
        pl.plot(x)
        pl.clf()
        pl.cla()
    print('')
    print('Regression PASSED')
    print('')
    return []

def data():
    return []
