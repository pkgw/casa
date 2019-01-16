import pdb
import shutil
import os
import numpy as np
import pylab as pl
import dateutil
from datetime import datetime
from casac import casac
me=casac.measures()
tb=casac.table()
qa=casac.quanta()

def fixpointing_vlass(vis, col='DIRECTION'):
    """
    this function will filter the pointing for jumps in pointing values that is seen in 
    VLASS pointing OTF pointing table and fix the interval interpretation to be MS V2.0 
    compliant
    """
    if(not os.path.exists(vis+'/POINTING_ORIG')):
        shutil.copytree(vis+'/POINTING', vis+'/POINTING_ORIG')
    tb.open(vis+'/POINTING', nomodify=False)
    ants=tb.getcol('ANTENNA_ID')
    ants=np.unique(ants)
    dref='J2000'
    eref=tb.getcolkeyword('TIME','MEASINFO')['Ref']
    colkeyw=tb.getcolkeyword(col,'MEASINFO')
    pref=tb.getcolkeyword(col,'MEASINFO')['Ref']
#    colkeyw['Ref']=dref
#    tb.putcolkeyword(col, 'MEASINFO', colkeyw)
    me.doframe(me.observatory('VLA'))
    for ant in ants:
        st=tb.query('ANTENNA_ID=='+str(ant))
        print st.nrows()
        ptimes=st.getcol('TIME')
        dirv_orig=st.getcol(col)
        sh=dirv_orig.shape
        dirv=dirv_orig.reshape((sh[0],sh[len(sh)-1]))
        rap=np.zeros(len(ptimes))
        t_rap=np.zeros(len(ptimes))
        decp=np.zeros(len(ptimes))
        epm=me.epoch(eref,str(ptimes[0])+'s')
        rdp=me.direction(pref,str(dirv[0,0])+'rad',str(dirv[1,0])+'rad')
        for i in range(len(ptimes)):
            epm=me.epoch('utc', qa.quantity(ptimes[i], 's'))
            me.doframe(epm)
            rdp['m0']['value']=dirv[0,i]
            rdp['m1']['value']=dirv[1,i]
            aep=me.measure(rdp,dref)  
            rap[i]=qa.convert(aep['m0'], 'rad')['value']
            decp[i]=qa.convert(aep['m1'], 'rad')['value']
            if(ant==0):
                t_rap[i]=pl.date2num(dateutil.parser.parse(qa.time(qa.quantity(rap[i], 'rad'), form='hms', prec=10)[0]))
            
        if(ant==0):
            pl.figure(1)
            pl.clf()
            pl.plot_date(t_rap, decp,'o')
            pl.title('BEFORE RA v/s DEC IN J2000')
            pl.ylabel('DEC in rad') 
        filterjump(rap,decp)
        filterjump2(rap,decp, ptimes)
        filterjump2(rap,decp, ptimes)
        if(ant==0):
            pl.figure(2)
            pl.clf()
            pl.plot_date(t_rap, decp,'o')
            pl.title('AFTER RA v/s DEC IN J2000')
            pl.ylabel('DEC in rad') 
        ###lets reconvert rap,decp back to original frame
        aep=me.direction(dref, str(rap[0])+'rad', str(decp[0])+'rad')
        for i in range(len(ptimes)):
            epm=me.epoch('utc', qa.quantity(ptimes[i], 's'))
            me.doframe(epm)
            aep['m0']['value']=rap[i]
            aep['m1']['value']=decp[i]
            rdp=me.measure(aep,pref)             
            rap[i]=qa.convert(rdp['m0'], 'rad')['value']
            decp[i]=qa.convert(rdp['m1'], 'rad')['value']
        ###
        if(len(dirv_orig.shape)==3):
            dirv_orig[0,0,:]=rap
            dirv_orig[1,0,:]=decp
        else:
            dirv_orig[0,:]=rap
            dirv_orig[1,:]=decp
        inter=st.getcol('INTERVAL')
        inter[:]=-1.0
        st.putcol('INTERVAL', inter)
        st.putcol(col, dirv_orig)
    #pl.plot(rap, decp, 'o')
    tb.done()
def filterjump(rap, decp,  nsigma=5):
    dy=np.diff(decp,1)
    dx=np.diff(rap,1)
    #pdb.set_trace()
    rajmp=np.where((np.abs(dx-np.mean(dx))) > np.abs(nsigma*np.std(dx)))[0]
    decjmp=np.where((np.abs(dy-np.mean(dy))) > np.abs(nsigma*np.std(dy)))[0]
    jmp=rajmp if(len(rajmp) > len(decjmp)) else decjmp
    jmprange=locatepairs(jmp)
    print 'jumps at', jmprange
    for k in range(len(jmprange)):
        for j in range(jmprange[k][0], jmprange[k][1]):
            #print 'before', j , rap[j+1], decp[j+1]
            rap[j+1]=rap[j]+dx[jmprange[k][0]-1]
            decp[j+1]=decp[j]+dy[jmprange[k][0]-1]
            #print 'after', rap[j+1], decp[j+1]
        #print 'Before', jmprange[k][1] , rap[jmprange[k][1]], decp[jmprange[k][1]]
        rap[jmprange[k][1]]= rap[jmprange[k][1]+1]-dx[jmprange[k][1]+1]
        decp[jmprange[k][1]]= decp[jmprange[k][1]+1]-dy[jmprange[k][1]+1]
        #print 'After', jmprange[k][1] , rap[jmprange[k][1]], decp[jmprange[k][1]]
    

def filterjump2(rap, decp, ptimes, nsigma=2):
    dy=np.diff(decp,1)
    dx=np.diff(rap,1)
    dt=np.diff(ptimes,1)
    #pdb.set_trace()
    dydt=dy/dt
    dxdt=dx/dt
    dxdt=np.sqrt(dydt*dydt+dxdt*dxdt)
    #jmp=np.where((np.abs(dxdt-np.mean(dxdt))) > np.abs(nsigma*np.std(dxdt)))[0]
    jmp=np.where(dxdt > nsigma*np.mean(dxdt))[0]
    print 'jumps ' , jmp
    for k in jmp:
            #print 'before', k , rap[k+1], decp[k+1]
            rap[k+1]=rap[k]+dx[k-1]
            decp[k+1]=decp[k]+dy[k-1]
            #print 'after', rap[k+1], decp[k+1]
    
def locatepairs(jmp):
    pairs=[]
    for k in range(len(jmp)-1):
        if( np.abs(jmp[k+1]-jmp[k]) < 4):
            pairs.append([jmp[k], jmp[k+1]])
    return pairs
     
        
