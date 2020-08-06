import os
import sys
import math
import shutil
import string
import time
import numpy as np
import math

try:
    # CASA 6
    from casatools import table
    # Most of this helper file is about table operations. The ms and image tools are used
    # only for a couple of functions (for which there might be a better place)
    from casatools import ms, image

    tb_local = table()
    tb_local2 = table()
    ms_local = ms()
    image_local = image()
except ImportError:
    # CASA 5
    from taskinit import tbtool
    from taskinit import mstool, iatool

    tb_local = tbtool()
    tb_local2 = tbtool()
    ms_local = mstool()
    image_local = iatool()


'''
A set of common helper functions for unit tests:
   compTables - compare two CASA tables
   compVarColTables - Compare a variable column of two tables
   DictDiffer - a class with methods to take a difference of two 
                Python dictionaries
   verify_ms - Function to verify spw and channels information in an MS   
   create_input - Save the string in a text file with the given name           
'''

def phasediffabsdeg(c1, c2):
    try:
        a = c1.imag
        a = c2.imag
    except:
        print("Phase difference of real numbers is always zero.")
        return 0.

    a = math.atan2(c1.imag, c1.real)
    b = math.atan2(c2.imag, c2.real)
    diff = abs(a-b)
    if diff>np.pi:
        diff = 2*np.py - diff
    return diff/np.pi*180. # (degrees)

def compTables(referencetab, testtab, excludecols, tolerance=0.001, mode="percentage", startrow = 0, nrow = -1, rowincr = 1):

    """
    compTables - compare two CASA tables
    
       referencetab - the table which is assumed to be correct

       testtab - the table which is to be compared to referencetab

       excludecols - list of column names which are to be ignored

       tolerance - permitted fractional difference (default 0.001 = 0.1 percent)

       mode - comparison is made as "percentage", "absolute", "phaseabsdeg" (for complex numbers = difference of the phases in degrees)  
    """

    rval = True

    tb_local.open(referencetab)
    cnames = tb_local.colnames()

    tb_local2.open(testtab)

    try:
        for c in cnames:
            if c in excludecols:
                continue
            
            print(("\nTesting column " + c))
            
            a = 0
            try:
                a = tb_local.getcol(c,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                rval = False
                print(('Error accessing column ', c, ' in table ', referencetab))
                print((sys.exc_info()[0]))
                break

            b = 0
            try:
                b = tb_local2.getcol(c,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                rval = False
                print(('Error accessing column ', c, ' in table ', testtab))
                print((sys.exc_info()[0]))
                break

            if not (len(a)==len(b)):
                print(('Column ',c,' has different length in tables ', referencetab, ' and ', testtab))
                print(a)
                print(b)
                rval = False
                break
            else:
                differs = False
                if not (a==b).all():
                    for i in range(0,len(a)):
                        if (isinstance(a[i],float)):
                            if ((mode=="percentage") and (abs(a[i]-b[i]) > tolerance*abs(a[i]))) or ((mode=="absolute") and (abs(a[i]-b[i]) > tolerance)):
                                print(("Column " + c + " differs"))
                                print(("Row=" + str(i)))
                                print(("Reference file value: " + str(a[i])))
                                print(("Input file value: " + str(b[i])))
                                if (mode=="percentage"):
                                    print(("Tolerance is {0}%; observed difference was {1} %".format (tolerance * 100, 100*abs(a[i]-b[i])/abs(a[i]))))
                                else:
                                    print(("Absolute tolerance is {0}; observed difference: {1}".format (tolerance, (abs(a[i]-b[i])))))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],int) or isinstance(a[i],np.int32)):
                            if (abs(a[i]-b[i]) > 0):
                                print(("Column " + c + " differs"))
                                print(("Row=" + str(i)))
                                print(("Reference file value: " + str(a[i])))
                                print(("Input file value: " + str(b[i])))
                                if (mode=="percentage"):
                                    print(("tolerance in % should be " + str(100*abs(a[i]-b[i])/abs(a[i]))))
                                else:
                                    print(("absolute tolerance should be " + str(abs(a[i]-b[i]))))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],str) or isinstance(a[i],np.bool_)):
                            if not (a[i]==b[i]):
                                print(("Column " + c + " differs"))
                                print(("Row=" + str(i)))
                                print(("Reference file value: " + str(a[i])))
                                print(("Input file value: " + str(b[i])))
                                if (mode=="percentage"):   
                                    print(("tolerance in % should be " + str(100*abs(a[i]-b[i])/abs(a[i]))))
                                else:
                                    print(("absolute tolerance should be " + str(abs(a[i]-b[i]))))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],list)) or (isinstance(a[i],np.ndarray)):
                            for j in range(0,len(a[i])):
                                if differs: break
                                if ((isinstance(a[i][j],float)) or (isinstance(a[i][j],int))):
                                    if ((mode=="percentage") and (abs(a[i][j]-b[i][j]) > tolerance*abs(a[i][j]))) or ((mode=="absolute") and (abs(a[i][j]-b[i][j]) > tolerance)):
                                        print(("Column " + c + " differs"))
                                        print(("(Row,Element)=(" + str(j) + "," + str(i) + ")"))
                                        print(("Reference file value: " + str(a[i][j])))
                                        print(("Input file value: " + str(b[i][j])))
                                        if (mode=="percentage"):
                                            print(("Tolerance in % should be " + str(100*abs(a[i][j]-b[i][j])/abs(a[i][j]))))
                                        else:
                                            print(("Absolute tolerance should be " + str(abs(a[i][j]-b[i][j]))))
                                        differs = True
                                        rval = False
                                        break
                                elif (isinstance(a[i][j],list)) or (isinstance(a[i][j],np.ndarray)):
                                    it = list(range(0,len(a[i][j])))
                                    if mode=="percentage":
                                        diff = np.abs(np.subtract(a[i][j], b[i][j])) > tolerance * np.abs(a[i][j])
                                        it = np.where(diff)[0]
                                    elif (mode=="absolute"):
                                        diff = np.abs(np.subtract(a[i][j], b[i][j])) > tolerance
                                        it = np.where(diff)[0]
                                    for k in it:
                                        if differs: break
                                        if ( ((mode=="percentage") and (abs(a[i][j][k]-b[i][j][k]) > tolerance*abs(a[i][j][k]))) \
                                                 or ((mode=="absolute") and (abs(a[i][j][k]-b[i][j][k]) > tolerance)) \
                                                 or ((mode=="phaseabsdeg") and (phasediffabsdeg(a[i][j][k],b[i][j][k])>tolerance)) \
                                                 ):
                                            print(("Column " + c + " differs"))
                                            print(("(Row,Channel,Corr)=(" + str(k) + "," + str(j) + "," + str(i) + ")"))
                                            print(("Reference file value: " + str(a[i][j][k])))
                                            print(("Input file value: " + str(b[i][j][k])))
                                            if (mode=="percentage"):
                                                print(("Tolerance in % should be " + str(100*abs(a[i][j][k]-b[i][j][k])/abs(a[i][j][k]))))
                                            elif (mode=="absolute"):
                                                print(("Absolute tolerance should be " + str(abs(a[i][j][k]-b[i][j][k]))))
                                            elif (mode=="phaseabsdeg"):
                                                print(("Phase tolerance in degrees should be " + str(phasediffabsdeg(a[i][j][k],b[i][j][k]))))
                                            else:
                                                print(("Unknown comparison mode: ",mode))
                                            differs = True
                                            rval = False
                                            break                                          
                                            
                        else:
                            print(("Unknown data type: ",type(a[i])))
                            differs = True
                            rval = False
                            break
                
                if not differs: print(("Column " + c + " PASSED"))
    finally:
        tb_local.close()
        tb_local2.close()

    return rval

    
def compVarColTables(referencetab, testtab, varcol, tolerance=0.):
    '''Compare a variable column of two tables.
       referencetab  --> a reference table
       testtab       --> a table to verify
       varcol        --> the name of a variable column (str)
       Returns True or False.
    '''
    
    retval = True

    tb_local.open(referencetab)
    cnames = tb_local.colnames()

    tb_local2.open(testtab)
    col = varcol
    if tb_local.isvarcol(col) and tb_local2.isvarcol(col):
        try:
            # First check
            if tb_local.nrows() != tb_local2.nrows():
                print(('Length of %s differ from %s, %s!=%s'%(referencetab,testtab,len(rk),len(tk))))
                retval = False
            else:
                for therow in range(tb_local.nrows()):
            
                    rdata = tb_local.getcell(col,therow)
                    tdata = tb_local2.getcell(col,therow)

#                    if not (rdata==tdata).all():
                    if not rdata.all()==tdata.all():
                        if (tolerance>0.):
                            differs=False
                            for j in range(0,len(rdata)):
###                                if (type(rdata[j])==float or type(rdata[j])==int):
                                if ((isinstance(rdata[j],float)) or (isinstance(rdata[j],int))):
                                    if (abs(rdata[j]-tdata[j]) > tolerance*abs(rdata[j]+tdata[j])):
#                                        print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
#                                        print(therow, j)
#                                        print(rdata[j])
#                                        print(tdata[j])
                                        differs = True
###                                elif (type(rdata[j])==list or type(rdata[j])==np.ndarray):
                                elif (isinstance(rdata[j],list)) or (isinstance(rdata[j],np.ndarray)):
                                    for k in range(0,len(rdata[j])):
                                        if (abs(rdata[j][k]-tdata[j][k]) > tolerance*abs(rdata[j][k]+tdata[j][k])):
#                                            print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
#                                            print(therow, j, k)
#                                            print(rdata[j][k])
#                                            print(tdata[j][k])
                                            differs = True
                                if differs:
                                    print(('ERROR: Column %s of %s and %s do not agree within tolerance %s'%(col,referencetab, testtab, tolerance)))
                                    retval = False
                                    break
                        else:
                            print(('ERROR: Column %s of %s and %s do not agree.'%(col,referencetab, testtab)))
                            print(('ERROR: First row to differ is row=%s'%therow))
                            retval = False
                            break
        finally:
            tb_local.close()
            tb_local2.close()
    
    else:
        print('Columns are not varcolumns.')
        retval = False

    if retval:
        print(('Column %s of %s and %s agree'%(col,referencetab, testtab)))
        
    return retval

    
        
class DictDiffer(object):
    """
    Calculate the difference between two dictionaries as:
    (1) items added
    (2) items removed
    (3) keys same in both but changed values
    (4) keys same in both and unchanged values
    Example:
            mydiff = DictDiffer(dict1, dict2)
            mydiff.changed()  # to show what has changed
    """
    def __init__(self, current_dict, past_dict):
        self.current_dict, self.past_dict = current_dict, past_dict
        self.set_current, self.set_past = set(current_dict.keys()), set(past_dict.keys())
        self.intersect = self.set_current.intersection(self.set_past)
    def added(self):
        return self.set_current - self.intersect 
    def removed(self):
        return self.set_past - self.intersect 
    def changed(self):
        return set(o for o in self.intersect if self.past_dict[o] != self.current_dict[o])            
    def unchanged(self):
        return set(o for o in self.intersect if self.past_dict[o] == self.current_dict[o])


def verifyMS(msname, expnumspws, expnumchan, inspw, expchanfreqs=[], ignoreflags=False):
    '''Function to verify spw and channels information in an MS
       msname        --> name of MS to verify
       expnumspws    --> expected number of SPWs in the MS
       expnumchan    --> expected number of channels in spw
       inspw         --> SPW ID
       expchanfreqs  --> numpy array with expected channel frequencies
       ignoreflags   --> do not check the FLAG column
           Returns a list with True or False and a state message'''
    
    msg = ''
    tb_local.open(msname+'/SPECTRAL_WINDOW')
    nc = tb_local.getcell("NUM_CHAN", inspw)
    nr = tb_local.nrows()
    cf = tb_local.getcell("CHAN_FREQ", inspw)
    tb_local.close()
    # After channel selection/average, need to know the exact row number to check,
    # ignore this check in these cases.
    if not ignoreflags:
        tb_local.open(msname)
        dimdata = tb_local.getcell("FLAG", 0)[0].size
        tb_local.close()
        
    if not (nr==expnumspws):
        msg =  "Found "+str(nr)+", expected "+str(expnumspws)+" spectral windows in "+msname
        return [False,msg]
    if not (nc == expnumchan):
        msg = "Found "+ str(nc) +", expected "+str(expnumchan)+" channels in spw "+str(inspw)+" in "+msname
        return [False,msg]
    if not ignoreflags and (dimdata != expnumchan):
        msg = "Found "+ str(dimdata) +", expected "+str(expnumchan)+" channels in FLAG column in "+msname
        return [False,msg]

    if not (expchanfreqs==[]):
        print("Testing channel frequencies ...")
#        print(cf)
#        print(expchanfreqs)
        if not (expchanfreqs.size == expnumchan):
            msg =  "Internal error: array of expected channel freqs should have dimension ", expnumchan
            return [False,msg]
        df = (cf - expchanfreqs)/expchanfreqs
        if not (abs(df) < 1E-8).all:
            msg = "channel frequencies in spw "+str(inspw)+" differ from expected values by (relative error) "+str(df)
            return [False,msg]

    return [True,msg]


def getChannels(msname, spwid, chanlist):
    '''From a list of channel indices, return their frequencies
       msname       --> name of MS
       spwid        --> spw ID
       chanlist     --> list of channel indices
    Return a numpy array, the same size of chanlist, with the frequencies'''
    
    try:
        try:
            tb_local.open(msname+'/SPECTRAL_WINDOW')
        except:
            print(('Cannot open table '+msname+'SPECTRAL_WINDOW'))
            
        cf = tb_local.getcell("CHAN_FREQ", spwid)
        
        # Get only the requested channels
        b = [cf[i] for i in chanlist]
        selchans = np.array(b)
    
    finally:
        tb_local.close()
        
    return selchans

def get_channel_freqs_widths(msname, spwid):
    '''
    Get frequencies and widths of all the channels for an spw ID
       msname       --> name of MS
       spwid        --> spw ID

    Return two numpy arrays (frequencies, widths), each of the same length as the number of
    channels'''

    try:
        spw_table = os.path.join(msname, 'SPECTRAL_WINDOW')
        try:
            tb_local.open(spw_table)
        except RuntimeError:
            print(('Cannot open table: {0}').format(spw_table))

        freqs = tb_local.getcell("CHAN_FREQ", spwid)
        widths = tb_local.getcell("CHAN_WIDTH", spwid)

    finally:
        tb_local.close()

    return freqs, widths

def getColDesc(table, colname):
    '''Get the description of a column in a table
       table    --> name of table or MS
       colname  --> column name
    Return a dictionary with the column description'''
    
    coldesc = {}
    try:
        try:
            tb_local.open(table)
            tcols = tb_local.colnames()
            if tcols.__contains__(colname):
                coldesc = tb_local.getcoldesc(colname)
        except:
            pass                        
    finally:
        tb_local.close()
        
    return coldesc

def getVarCol(table, colname):
    '''Return the requested variable column
       table    --> name of table or MS
       colname  --> column name
    Return the column as a dictionary'''
    
    col = {}
    try:
        try:
            tb_local.open(table)
            col = tb_local.getvarcol(colname)
        except:
            print(('Cannot open table '+table))

    finally:
        tb_local.close()
        
    return col
   
def createInput(str_text, filename):
    '''Save the string in a text file with the given name
    str_text    --> string to save
    filename    --> name of the file to save
            It will remove the filename if it exist!'''
    
    inp = filename    
    cmd = str_text
    
    # remove file first
    if os.path.exists(inp):
        os.system('rm -f '+ inp)
    
    try:
        # save to a file    
        with open(inp, 'w') as f:
            f.write(cmd)
            
    finally:  
        f.close()
    
    return

def calculateHanning(dataB,data,dataA):
    '''Calculate the Hanning smoothing of each element'''
    const0 = 0.25
    const1 = 0.5
    const2 = 0.25
    S = const0*dataB + const1*data + const2*dataA
    return S


def getTileShape(mydict, column='DATA'):
    '''Return the value of TileShape for a given column
       in the dictionary from data managers (tb.getdminfo).
       mydict --> dictionary from tb.getdminfo()
       column --> column where to look for TileShape'''
    
    tsh = {}
    for key, value in list(mydict.items()):
        if mydict[key]['COLUMNS'][0] == column:
             # Dictionary for requested column
            hyp = mydict[key]['SPEC']['HYPERCUBES']
                    
             # This is the HYPERCUBES dictionary
            for hk in list(hyp.keys()):
                tsh = hyp[hk]['TileShape']
                break
                    
            break
    
    return tsh

def checkwithtaql(taqlstring):
    os.system('rm -rf mynewtable.tab')
    tb_local.create('mynewtable.tab')
    tb_local.open('mynewtable.tab',nomodify=False)
    rval = tb_local.taql(taqlstring)
    tb_local.close()
    therval = rval.nrows()
    tmpname = rval.name()
    rval.close()
    os.system('rm -rf mynewtable.tab')
    os.system('rm -rf '+tmpname)
    print(("Found ", therval, " rows in selection."))
    return therval


def compcaltabnumcol(cal1, cal2, tolerance, colname1='CPARAM', colname2="CPARAM", testspw=None):
    print(("Comparing column "+colname1+" of caltable "+cal1))
    print(("     with column "+colname2+" of caltable "+cal2))
    if testspw!=None:
        print(("for SPW "+str(testspw)+" only."))
    print("Discrepant row search ...")
    rval = False
    try:
        discrepantrows = -1
        if(testspw==None):
            discrepantrows = checkwithtaql("select from [select from "+cal1+" orderby TIME, FIELD_ID, SPECTRAL_WINDOW_ID, ANTENNA1, ANTENNA2 ] t1, [select from "+cal2+" orderby TIME, FIELD_ID, SPECTRAL_WINDOW_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1."+colname1+",t2."+colname2+", "+str(tolerance)+")))")
        else:
            discrepantrows = checkwithtaql("select from [select from "+cal1+" where SPECTRAL_WINDOW_ID=="+str(testspw)+" orderby TIME, FIELD_ID, ANTENNA1, ANTENNA2 ] t1, [select from "+cal2+" where SPECTRAL_WINDOW_ID=="+str(testspw)+" orderby TIME, FIELD_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1."+colname1+",t2."+colname2+", "+str(tolerance)+")))")
        if discrepantrows==0:
            print("The two columns agree.")
            rval = True
    except Exception as instance:
        print(("Error: "+str(instance)))

    return rval

                    
def compmsmainnumcol(vis1, vis2, tolerance, colname1='DATA', colname2="DATA"):
    print(("Comparing column "+colname1+" of MS "+vis1))
    print(("     with column "+colname2+" of MS "+vis2))
    print("Discrepant row search ...")
    rval = False
    try:
        discrepantrows = checkwithtaql("select from [select from "+vis1+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "+vis2+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(near(t1."+colname1+",t2."+colname2+", "+str(tolerance)+")))")
        if discrepantrows==0:
            print("The two columns agree.")
            rval = True
    except Exception as instance:
        print(("Error: "+str(instance)))

    return rval

def compmsmainboolcol(vis1, vis2, colname1='FLAG', colname2='FLAG'):
    print(("Comparing column "+colname1+" of MS "+vis1))
    print(("     with column "+colname2+" of MS "+vis2))
    print("Discrepant row search ...")
    rval = False
    try:
        discrepantrows = checkwithtaql("select from [select from "+vis1+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t1, [select from "+vis2+" orderby TIME, DATA_DESC_ID, ANTENNA1, ANTENNA2 ] t2 where (not all(t1."+colname1+"==t2."+colname2+"))")
        if discrepantrows==0:
            print("The two columns agree.")
            rval = True
    except Exception as instance:
        print(("Error: "+str(instance)))

    return rval

def compareSubTables(input,reference,order=None,excluded_cols=[]):
    
    tbinput = tb_local
    tbinput.open(input)
    if order is not None:
        tbinput_sorted = tbinput.taql("SELECT * from " + input + " order by " + order)
    else:
        tbinput_sorted = tbinput
    
    tbreference = tb_local2
    tbreference.open(reference)
    if order is not None:
        tbreference_sorted = tbreference.taql("SELECT * from " + reference + " order by " + order)
    else:
        tbreference_sorted = tbreference
    
    columns = tbinput.colnames()
    for col in columns:
        if not col in excluded_cols:
            col_input = tbinput_sorted.getcol(col)
            col_reference = tbreference_sorted.getcol(col)
            if not (col_input == col_reference).all():
                tbinput.close()
                tbreference.close()
                return (False,col)

    tbinput.close()
    tbreference.close()
    if order is not None:
        tbinput_sorted.close()
        tbreference_sorted.close()

    return (True,"OK")

def getColShape(tab,col,start_row=0,nrow=1,row_inc=1):
    """ Get the shape of the given column.
    Keyword arguments:
        tab        --    input table or MS
        col        --    column to get the shape
        start_row  --    start row (default 0)
        nrow       --    number of rows to read (default 1)
        row_inc    --    increment of rows to read (default 1)
        
        Return a list of strings with the shape of each row in the column.
    
    """

    col_shape = []
    try:
        try:
            tb_local.open(tab)
            col_shape = tb_local.getcolshapestring(col,start_row,nrow,row_inc)
        except:
            print(('Cannot get shape of col %s from table %s '%(col,tab)))

    finally:
        tb_local.close()
            
    return col_shape

def findTemplate(testname,refimage,copy=False):
    """
    find a template image (or MS - it does assume its a directory)
    look in order in:
    REGRESSION_DATA/regression/testname/refimage
    CASAPATH/data/regression/testname/refimage
    REGRESSION_DATA/regression/testname/reference/refimage/
    CASAPATH/data/regression/testname/reference/refimage
    if copy=True, copy what's found to cwd
    """
    from os import F_OK
    try:
        datapaths=REGRESSION_DATA
    except:
        datapaths=[]
    datapaths.append(os.environ.get('CASAPATH').split()[0]+"/data/")
    possibilities=[x+'/regression/'+testname+'/'+refimage for x in datapaths]+[x+'/regression/'+testname+'/reference/'+refimage for x in datapaths] 

    #print(possibilities)
    from itertools import dropwhile
    try:
        found = next(dropwhile( lambda x: not os.access(x,F_OK),possibilities))
    except:
        raise IOError(" ERROR: "+refimage+" not found")
    if copy:
        from shutil import copytree
        print(("Copying "+found))
        copytree(found,msname)
    return found


# As opposed to most other functions in this file, this doesn't use the table tool but the
# image tool
def compImages(im0,im1,keys=['flux','min','max','maxpos','rms'],tol=1e-4,verbose=False):
    """
    compare two images using imstat and the specified keys, 
    to a tolerance tol, and printing the comparison if verbose==True
    note that the string keys like 'blcf' will fail
    """
    from os import F_OK
    if isinstance(tol,float):
        tol=tol+np.zeros(len(keys))
    ims=[im0,im1]
    s=[]
    for i in range(2):
        if not os.access(ims[i],F_OK): 
            print((ims[i]+" not found"))
            return False
        image_local.open(ims[1])
        s.append(image_local.statistics())
        image_local.done()
    status=True
    for ik in range(len(keys)):
        k=keys[ik]
        s0=s[0][k][0]
        s1=s[1][k][0]
        if abs(s0-s1)*2/(s0+s1)>tol[ik]: status=False
        if verbose:
            print((("%7s: "%k),s0,s1))
    return status


# As opposed to most other functions in this file, this doesn't use the table tool but the
# ms tool
def compMS(ms0,ms1,keys=['mean','min','max','rms'],ap="amp",tol=1e-4,verbose=False):
    """
    compare two MS using ms.statistics on amp or phase as specified, 
    and the specified keys, 
    to a tolerance tol, and printing the comparison if verbose==True
    """
    from os import F_OK
    if isinstance(tol,float):
        tol=tol+np.zeros(len(keys))
    mss=[ms0,ms1]
    s=[]
    for i in range(2):
        if not os.access(mss[i],F_OK): 
            print((mss[i]+" not found"))
            return False
        ms_local.open(mss[1])
        stats = ms_local.statistics("DATA",ap)
        s.append(stats[list(stats.keys())[0]])
        ms_local.done()
    status=True
    for ik in range(len(keys)):
        k=keys[ik]
        s0=s[0][k]
        s1=s[1][k]
        if abs(s0-s1)*2/(s0+s1)>tol[ik]: status=False
        if verbose:
            print((("%7s: "%k),s0,s1))
    return status

# A function object that can be passed to ignore parameter
# of shutil.copytree. It will ignore subversion directory
# when data are copied to working directory.
ignore_subversion = shutil.ignore_patterns('.svn')

def copytree_ignore_subversion(datadir, name, outname=None):
    if outname is None:
        outname = name
    if not os.path.exists(name):
        shutil.copytree(os.path.join(datadir, name), outname,
                        ignore=ignore_subversion)


def get_table_cache():
    cache = tb_local.showcache()
    # print('cache = {}'.format(cache))
    return cache


class TableCacheValidator(object):
    def __init__(self):
        self.original_cache = get_table_cache()

    def validate(self):
        cache = get_table_cache()
        #print 'original {} current {}'.format(self.original_cache, cache)
        return len(cache) == 0 or cache == self.original_cache
