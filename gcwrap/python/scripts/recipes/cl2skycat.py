from casac import casac
import numpy as np
def cl2skycat(componentlist='', skycat=''):
    """
    Converts a componentlist dictionary or componnent list file on disk
    to  a skycatalog to overlay on image in the viewer
    
    """
    qa=casac.quanta()
    cl=casac.componentlist()
    tb=casac.table()
    if(type(componentlist)==str):
        cl.open(componentlist)
    elif(type(componentlist)==dict):
        cl.purge()
        cl.fromrecord(componentlist)
    if(cl.length()==0):
            print("no components found")
            return

    des={}
    des['Type']={'valueType':'string'}
    des['Long']={'valueType':'double'}
    des['Lat']={'valueType':'double'}
    des['COMP_ID']={'valueType':'string'}
    des['RA']={'valueType':'string'}
    des['DEC']={'valueType':'string'}
    des['FluxValue']={'valueType':'double'}
    tb.create(tablename=skycat, tabledesc=des, nrow=cl.length())
    eltype=[]
    nam=[]
    RA=[]
    DEC=[]
    lati=np.zeros((cl.length(),))
    longi=np.zeros((cl.length(),))
    fluxval=np.zeros((cl.length(),))
    for k in range(cl.length()):
        longi[k]=qa.convert(cl.getrefdir(k)['m0'],'deg')['value']
        lati[k]=qa.convert(cl.getrefdir(k)['m1'],'deg')['value']
        fluxval[k]=cl.getfluxvalue(k)[0]
        RA.append(qa.time(cl.getrefdir(k)['m0'], prec=10))
        DEC.append(qa.angle(cl.getrefdir(k)['m1'], prec=10))
        eltype.append(cl.getrefdir(k)['refer'])
        nam.append(str(k))
    tb.putcol('Type', eltype)
    tb.putcol('RA', RA)
    tb.putcol('DEC', DEC)
    tb.putcol('COMP_ID', nam)
    tb.putcol('Long', longi)
    tb.putcol('Lat', lati)
    tb.putcol('FluxValue', fluxval)
    tb.putcolkeyword(columnname='Long', keyword='UNIT', value='deg')
    tb.putcolkeyword(columnname='Lat', keyword='UNIT', value='deg')                     
    tb.putinfo({'type':'Skycatalog'})
    tb.done()
