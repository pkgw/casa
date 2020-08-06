from casac import casac
import numpy as np
import re;
def field2skycat(msname='', skycat='',fieldpattern=''):
    """
    Converts field table in ms to skycatalog
    allow overlay on image in the viewer
    """
    enumTypes =['J2000','JMEAN','JTRUE','APP','B1950',
                'B1950_VLA','BMEAN','BTRUE','GALACTIC','HADEC','AZEL','AZELSW',
                'AZELGEO','AZELSWGEO','JNAT','ECLIPTIC','MECLIPTIC','TECLIPTIC','SUPERGAL',
                'ITRF','TOPO','ICRS']

    qa=casac.quanta()
    tb=casac.table()
    tb.open(msname+'/FIELD')
    dir=tb.getcol('PHASE_DIR')
    nam=tb.getcol('NAME')
    nfield=tb.nrows()
    n=0;
    for k in range(nfield):
        if (re.match(fieldpattern,nam[k]) != None):
            n=n+1

    eltype=[]
    if('VarRefCol' in tb.getcolkeyword('PHASE_DIR', 'MEASINFO')):
        typeid=tb.getcol(tb.getcolkeyword('PHASE_DIR', 'MEASINFO')['VarRefCol'])
        for k in range(nfield):
            if (re.match(fieldpattern,nam[k]) != None):
                eltype.append(enumTypes[typeid[k]])
    else:
        eltype=[tb.getcolkeyword('PHASE_DIR', 'MEASINFO')['Ref']]*nfield
    unitra=tb.getcolkeyword('PHASE_DIR', 'QuantumUnits')[0]
    unitdec=tb.getcolkeyword('PHASE_DIR', 'QuantumUnits')[1]
    tb.done()
    des={}
    des['Type']={'valueType':'string'}
    des['Long']={'valueType':'double'}
    des['Lat']={'valueType':'double'}
    des['FIELD_NAME']={'valueType':'string'}
    des['FIELD_ID']={'valueType':'string'}
    des['RA']={'valueType':'string'}
    des['DEC']={'valueType':'string'}


    tb.create(tablename=skycat, tabledesc=des, nrow=n);
    tb.putcol('Type', eltype)
    lati=np.zeros((n,))
    longi=np.zeros((n,))
    RA=[]
    DEC=[]
    fid=[]
    fieldname=[]
    n=0;
    for k in range(nfield):
        if (re.match(fieldpattern,nam[k]) != None):
            longi[n]=qa.convert(qa.quantity(dir[0,0,k], unitra),'deg')['value']
            lati[n]=qa.convert(qa.quantity(dir[1,0,k], unitdec), 'deg')['value']
            RA.append(qa.time(qa.quantity(dir[0,0,k], unitra), prec=10))
            DEC.append(qa.angle(qa.quantity(dir[1,0,0], unitdec), prec=10))
            fid.append(str(k))
            fieldname.append(nam[k]);
            n=n+1;
    tb.putcol('RA', RA)
    tb.putcol('DEC', DEC)
#    tb.putcol('FIELD_NAME', nam)
    tb.putcol('FIELD_NAME', fieldname);
    tb.putcol('FIELD_ID', fid)
    tb.putcol('Long', longi)
    tb.putcol('Lat', lati)
    tb.putcolkeyword(columnname='Long', keyword='UNIT', value='deg')
    tb.putcolkeyword(columnname='Lat', keyword='UNIT', value='deg')                     
    tb.putinfo({'type':'Skycatalog'})
    tb.done()
