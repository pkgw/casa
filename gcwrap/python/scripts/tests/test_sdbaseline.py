import os
import glob
import sys
import shutil
import numpy
from __main__ import default
from tasks import *
from taskinit import *
import unittest
#
import listing
from numpy import array

from sdbaseline import sdbaseline
from sdutil import tbmanager


try:
    from . import selection_syntax
except:
    import tests.selection_syntax as selection_syntax


### Utilities for reading blparam file
class FileReader(object):
    def __init__(self, filename):
        self.__filename = filename
        self.__data = None
        self.__nline = None

    def read(self):
        if self.__data is None:
            f = open(self.__filename, 'r')
            self.__data = f.readlines()
            f.close()
            self.__nline = len(self.__data)
        return

    def nline(self):
        self.read()
        return self.__nline

    def index(self, txt, start):
        return self.__data[start:].index(txt) + 1 + start

    def getline(self, idx):
        return self.__data[idx]

class BlparamFileParser(FileReader):
    def __init__(self, blfile):
        FileReader.__init__(self, blfile)
        self.__nrow = None
        self.__coeff = None
        self.__rms = None
        self.__ctxt = 'Baseline parameters\n'
        self.__rtxt = 'Results of baseline fit\n'

    def nrow(self):
        self.read()
        if self.__nrow is None:
            return self._nrow()
        else:
            return self.__nrow

    def coeff(self):
        self.read()
        if self.__coeff is None:
            self.parseCoeff()
        return self.__coeff

    def rms(self):
        self.read()
        if self.__rms is None:
            self.parseRms()
        return self.__rms

    def _nrow(self):
        self.__nrow = 0
        for i in range(self.nline()):
            if self.getline(i) == self.__ctxt:
                self.__nrow += 1
        return self.__nrow

    def parse(self):
        self.read()
        self.parseCoeff()
        self.parseRms()
        return
        
    def parseCoeff(self):
        self.__coeff = []
        nrow = self.nrow()
        idx = 0
        while (len(self.__coeff) < nrow):
            try:
                idx = self.index(self.__ctxt, idx)
                coeffs = []
                while(self.getline(idx) != self.__rtxt):
                    coeff = self.__parseCoeff(idx)
                    coeffs += coeff
                    idx += 1
                self.__coeff.append(coeffs)
            except:
                break
        return

    def parseRms(self):
        self.__rms = []
        nrow = self.nrow()
        idx = 0
        while (len(self.__rms) < nrow):
            try:
                idx = self.index(self.__rtxt, idx)
                self.__rms.append(self.__parseRms(idx))
            except:
                break   
        return

    def __parseCoeff(self, idx):
        return parseCoeff(self.getline(idx))

    def __parseRms(self, idx):
        return parseRms(self.getline(idx))

def parseCoeff(txt):
    clist = txt.rstrip('\n').split(',')
    ret = []
    for c in clist:
        ret.append(float(c.split('=')[1]))
    return ret
    
def parseRms(txt):
    t = txt.lstrip().rstrip('\n')[6:]
    return float(t)

class sdbaseline_unittest_base(unittest.TestCase):
    """
    Base class for sdbaseline unit test
    """
    # Data path of input/output
    datapath = os.environ.get('CASAPATH').split()[0] + \
              '/data/regression/unittest/tsdbaseline/'
    taskname = "sdbaseline"
    verboselog = False

    #complist = ['max','min','rms','median','stddev']

    blparam_order = ['row', 'pol', 'mask', 'nclip', 'cthre',
                     'uself', 'lthre', 'ledge', 'redge', 'chavg',
                     'btype', 'order', 'npiec', 'nwave']
    blparam_dic = {}
    blparam_dic['row']   = [0, 0, 1, 1, 2, 2, 3, 3]
    blparam_dic['pol']   = [0, 1, 0, 1, 0, 1, 0, 1]
    #blparam_dic['mask']  = ['0~4000;6000~8000']*3 + ['']*5
    blparam_dic['mask']  = ['500~2500;5000~7500']*8
    blparam_dic['nclip'] = [0]*8
    blparam_dic['cthre'] = ['3.']*8
    blparam_dic['uself'] = ['false']*4 + ['true'] + ['false']*3
    blparam_dic['lthre'] = ['0.']*4 + ['3.', '', '', '0.']
    blparam_dic['ledge'] = [0]*4 + [10, 50, '', 0]
    blparam_dic['redge'] = [0]*4 + [10, 50, '', 0]
    blparam_dic['chavg'] = [0]*4 + [4, '', '', 0]
    blparam_dic['btype'] = ['poly'] + ['chebyshev']*2 + ['poly', 'chebyshev', 'poly'] + ['cspline']*2
    blparam_dic['order'] = [0, 0, 1, 1, 2, 2, '', '']
    blparam_dic['npiec'] = [0]*6 + [1]*2
    blparam_dic['nwave'] = [[]]*3 + ['']*2 + [[]]*3

    ### helper functions for tests ###
    def _createBlparamFile(self, file, param_order, val, option=''):
        nspec = 8
        f = open(file, 'w')
        assert(len(param_order) == len(list(val.keys())))
        for key in list(val.keys()):
            assert(len(val[key]) == nspec)
        for i in range(nspec):
            do_write = True
            s = ''
            for key in param_order:
                v = val[key][i]
                if key == 'nwave':
                    if v != '':
                        s += ','
                        s += str(v)
                else:
                    s += str(v)
                    if key != 'npiec': s += ','
            s += '\n'
            if (option == 'r2p1less') and (val['row'][i] == 2) and (val['pol'][i] == 1):
                do_write = False
            if (option == 'r2p1cout') and (val['row'][i] == 2) and (val['pol'][i] == 1):
                s = '#' + s
            if do_write:
                f.write(s)
        f.close()

    def _checkfile(self, name, fail=True):
        """
        Check if the file exists.
        name : the path and file name to test
        fail : if True, Error if the file does not exists.
               if False, return if the file exists
        """
        isthere=os.path.exists(name)
        if fail:
            self.assertTrue(isthere,
                            msg='Could not find, %s'%(name))
        else: return isthere

    def _remove(self, names):
        """
        Remove a list of files and directories from disk
        """
        for name in names:
            if os.path.exists(name):
                if os.path.isdir(name):
                    shutil.rmtree(name)
                else:
                    os.remove(name)

    def _copy(self, names, from_dir=None, dest_dir=None):
        """
        Copy a list of files and directories from a directory (from_dir) to
        another (dest_dir) in the same name.
        
        names : a list of files and directories to copy
        from_dir : a path to directory from which search and copy files
                   and directories (the default is the current path)
        to_dir   : a path to directory to which copy files and directories
                   (the default is the current path)
        NOTE: it is not allowed to specify 
        """
        # Check for paths
        if from_dir==None and dest_dir==None:
            raise ValueError("Can not copy files to exactly the same path.")
        from_path = os.path.abspath("." if from_dir==None else from_dir.rstrip("/"))
        to_path = os.path.abspath("." if dest_dir==None else dest_dir.rstrip("/"))
        if from_path == to_path:
            raise ValueError("Can not copy files to exactly the same path.")
        # Copy a list of files and directories
        for name in names:
            from_name = from_path + "/" + name
            to_name = to_path + "/" + name
            if os.path.exists(from_name):
                if os.path.isdir(from_name):
                    shutil.copytree(from_name, to_name)
                else:
                    shutil.copyfile(from_name, to_name)
                if self.verboselog:
                    casalog.post("Copying '%s' FROM %s TO %s" % (name, from_path, to_path))
            else:
                casalog.post("Could not find '%s'...skipping copy" % from_name, 'WARN')
    
    def _getUniqList(self, val):
        """Accepts a python list and returns a list of unique values"""
        if not isinstance(val, list):
            raise Exception('_getUniqList: input value must be a list.')
        return list(set(val))

    def _getListSelection(self, val):
        """
        Converts input to a list of unique integers
        Input: Either comma separated string of IDs, an integer, or a list of values.
        Output: a list of unique integers in input arguments for string and integer input.
                In case the input is a list of values, output will be a list of unique values.
        """
        if isinstance(val, str):
            val_split = val.split(',')
            val_sel = []
            for j in range(len(val_split)):
                val_sel.append(int(val_split[j]))
        elif isinstance(val, int):
            val_sel = [val]
        elif isinstance(val, list) or isinstance(val, tuple):
            val_sel = val.copy()
        else:
            raise Exception('_getListSelection: wrong value ' + str(val) + ' for selection.')
        return self._getUniqList(val_sel)
    
    def _getListSelectedRowID(self, data_list, sel_list):
        """
        Returns IDs of data_list that contains values equal to one in
        sel_list.
        The function is used to get row IDs that corresponds to a
        selected IDs. In that use case, data_list is typically a list
        of values in a column of an MS (e.g., SCAN_NUMBER) and sel_list is
        a list of selected (scan) IDs.

        data_list : a list to test and get IDs from
        sel_list  : a list of values to look for existance in data_list
        """
        res = []
        for i in range(len(data_list)):
            if data_list[i] in sel_list:
                #idx = sel_list.index(data_list[i])
                res.append(i)
        return self._getUniqList(res)
    
    def _getEffective(self, spec, mask):
        """
        Returns an array made by selected elements in spec array.
        Only the elements in the ID range in mask are returned.

        spec : a data array
        mask : a mask list of the channel ranges to use. The format is
               [[start_idx0, end_idx0], [start_idx1, end_idx1], ...]
        """
        res = []
        for i in range(len(mask)):
            for j in range(mask[i][0], mask[i][1]):
                res.append(spec[j])
        return numpy.array(res)

    def _getStats(self, filename=None, spw=None, pol=None, colname=None, mask=None):
        """
        Returns a list of statistics dictionary of selected rows in an MS.

        filename : the name of MS
        spw      : spw ID selection (default: all spws in MS)
        pol      : pol ID selection (default: all pols in MS)
        colname  : the name of data column (default: 'FLOAT_DATA')
        mask     : a mask list of the channel ranges to use. The format is
                   [[start_idx0, end_idx0], [start_idx1, end_idx1], ...]
        
        The order of output list is in the ascending order of selected row IDs.
        The dictionary in output list has keys:
        'row' (row ID in MS), 'pol' (pol ID), 'rms', 'min', 'max', 'median',
        and 'stddev'
        """
        # Get selected row and pol IDs in MS. Also get spectrumn in the MS
        if not spw: spw = ''
        select_spw = (spw not in ['', '*'])
        if select_spw: spw_sel = self._getListSelection(spw)
        if not pol: pol = ''
        select_pol = (pol not in ['', '*'])
        if select_pol: pol_sel = self._getListSelection(pol)
        if not colname: colname='FLOAT_DATA'
        self._checkfile(filename)
        with tbmanager(filename) as tb:
            data = tb.getcol(colname)
            ddid = tb.getcol('DATA_DESC_ID')
        with tbmanager(filename+'/DATA_DESCRIPTION') as tb:
            spwid = tb.getcol('SPECTRAL_WINDOW_ID').tolist()
        if not select_spw: spw_sel = spwid
        # get the selected DD IDs from selected SPW IDs.
        dd_sel = self._getListSelectedRowID(spwid, spw_sel)
        # get the selected row IDs from selected DD IDs
        row_sel = self._getListSelectedRowID(ddid, dd_sel)
        if not select_spw: row_sel = list(range(len(ddid)))
        if not select_pol: pol_sel = list(range(len(data)))

        res = []
        for irow in row_sel:
            for ipol in pol_sel:
                spec = data[ipol,:,irow]
                res_elem = self._calc_stats_of_array(spec, mask=mask)
                res_elem['row'] = irow
                res_elem['pol'] = ipol
                
                res.append(res_elem)

        return res

    def _calc_stats_of_array(self, data, mask=None):
        """
        """
        if mask is not None:
            spec = self._getEffective(data, mask)
        else:
            spec = numpy.array(data)
        res_elem = {}
        res_elem['rms'] = numpy.sqrt(numpy.var(spec))
        res_elem['min'] = numpy.min(spec)
        res_elem['max'] = numpy.max(spec)
        spec_mea = numpy.mean(spec)
        res_elem['median'] = numpy.median(spec)
        res_elem['stddev'] = numpy.std(spec)
        return res_elem
        

    def _convert_statslist_to_dict(self, stat_list):
        """
        Returns a disctionary of statistics of selected rows in an MS.

        stat_list: a list of stats dictionary (e.g., return value of _getStats)

        The output dictionary is in form:
        {'max': [max0, max1, max2, ...], 'min': [min0, min1,...], ...}
        The order of elements are in ascending order of row and pol IDs pair, i.e.,
        (row0, pol0), (row0, pol1), (row1, pol0), ....
        """
        #if len(stat_list)==0: raise Exception, "No row selected in MS"
        keys=list(stat_list[0].keys())
        stat_dict={}
        for key in keys:
            stat_dict[key] = []
        for stat in stat_list:
            for key in keys:
                stat_dict[key].append(stat[key])
        return stat_dict

    def _compareStats(self, currstat, refstat, rtol=1.0e-2, atol=1.0e-5, complist=None):
        """
        Compare statistics results (dictionaries) and test if the values are within
        an allowed tolerance.

        currstat : the statistic values to test (either an MS name or
                   a dictionary)
        refstat  : the reference statistics values (a dictionary)
        rtol   : tolerance of relative difference
        atol   : tolerance of absolute difference
        complist : statistics to compare (default: keys in refstat)
        """
        # test if the statistics of baselined spectra are equal to
        # the reference values
        printstat = False #True
        # In case currstat is filename
        if isinstance(currstat, str) and os.path.exists(currstat):
            #print "calculating statistics from '%s'" % currstat
            currstat = self._getStats(currstat)

        self.assertTrue(isinstance(currstat,dict) and \
                        isinstance(refstat, dict),\
                        "Need to specify two dictionaries to compare")
        if complist:
            keylist = complist
        else:
            keylist = list(refstat.keys())
            #keylist = self.complist
        
        for key in keylist:
            self.assertTrue(key in currstat,\
                            msg="%s is not defined in the current results."\
                            % key)
            self.assertTrue(key in refstat,\
                            msg="%s is not defined in the reference data."\
                            % key)
            refval = refstat[key]
            currval = currstat[key]
            # Quantum values
            if isinstance(refval,dict):
                if 'unit' in refval and 'unit' in currval:
                    if printstat:
                        print("Comparing unit of '%s': %s (current run), %s (reference)" %\
                              (key,currval['unit'],refval['unit']))
                    self.assertEqual(refval['unit'],currval['unit'],\
                                     "The units of '%s' differs: %s (expected: %s)" % \
                                     (key, currval['unit'], refval['unit']))
                    refval = refval['value']
                    currval = currval['value']
                else:
                    raise Exception("Invalid quantum values. %s (current run) %s (reference)" %\
                                    (str(currval),str(refval)))
            currval = self._to_list(currval)
            refval = self._to_list(refval)
            if printstat:
                print("Comparing '%s': %s (current run), %s (reference)" %\
                      (key,str(currval),str(refval)))
            self.assertTrue(len(currval)==len(refval),"Number of elemnets in '%s' differs." % key)
            if isinstance(refval[0],str):
                for i in range(len(currval)):
                    if isinstance(refval[i],str):
                        self.assertTrue(currval[i]==refval[i],\
                                        msg="%s[%d] differs: %s (expected: %s) " % \
                                        (key, i, str(currval[i]), str(refval[i])))
            else:
                # numpy.allclose handles almost zero case more properly.
                self.assertTrue(numpy.allclose(currval, refval, rtol=rtol, atol=atol),
                                msg="%s differs: %s" % (key, str(currval)))
            del currval, refval

            
#     def _isInAllowedRange(self, testval, refval, reltol=1.e-2):
#         """
#         Check if a test value is within permissive relative difference from refval.
#         Returns a boolean.
#         testval & refval : two numerical values to compare
#         reltol           : allowed relative difference to consider the two
#                            values to be equal. (default 0.01)
#         """
#         denom = refval
#         if refval == 0:
#             if testval == 0:
#                 return True
#             else:
#                 denom = testval
#         rdiff = (testval-refval)/denom
#         del denom,testval,refval
#         return (abs(rdiff) <= reltol)

    def _to_list(self, input):
        """
        Convert input to a list
        If input is None, this method simply returns None.
        """
        import numpy
        listtypes = (list, tuple, numpy.ndarray)
        if input == None:
            return None
        elif type(input) in listtypes:
            return list(input)
        else:
            return [input]


    def _compareBLparam(self, out, reference):
        # test if baseline parameters are equal to the reference values
        # currently comparing every lines in the files
        # TO DO: compare only "Fitter range" and "Baseline parameters"
        self._checkfile(out)
        self._checkfile(reference)
        
        blparse_out = BlparamFileParser(out)
        blparse_out.parse()
        coeffs_out = blparse_out.coeff()
        rms_out = blparse_out.rms()
        blparse_ref = BlparamFileParser(reference)
        blparse_ref.parse()
        coeffs_ref = blparse_ref.coeff()
        rms_ref = blparse_ref.rms()
        allowdiff = 0.01
        print('Check baseline parameters:')
        for irow in range(len(rms_out)):
            print('Row %s:'%(irow))
            print('   Reference rms  = %s'%(rms_ref[irow]))
            print('   Calculated rms = %s'%(rms_out[irow]))
            print('   Reference coeffs  = %s'%(coeffs_ref[irow]))
            print('   Calculated coeffs = %s'%(coeffs_out[irow]))
            r0 = rms_ref[irow]
            r1 = rms_out[irow]
            rdiff = (r1 - r0) / r0
            self.assertTrue((abs(rdiff)<allowdiff),
                            msg='row %s: rms is different'%(irow))
            c0 = coeffs_ref[irow]
            c1 = coeffs_out[irow]
            for ic in range(len(c1)):
                rdiff = (c1[ic] - c0[ic]) / c0[ic]
                self.assertTrue((abs(rdiff)<allowdiff),
                                msg='row %s: coefficient for order %s is different'%(irow,ic))
        print('')
#         self.assertTrue(listing.compare(out,reference),
#                         'New and reference files are different. %s != %s. '
#                         %(out,reference))



class sdbaseline_basicTest(sdbaseline_unittest_base):
    """
    Basic unit tests for task sdbaseline. No interactive testing.

    List of tests:
    test000 --- default values for all parameters
    test001 --- polynominal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test002 --- Chebyshev polynominal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test003 --- cubic spline baselining with no mask (maskmode = 'list'). spw and pol specified.
    test004 --- sinusoidal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test050 --- existing file as outfile with overwrite=False (raises an exception)
    test051 --- no data after selection (raises an exception)

    Note: The input data 'OrionS_rawACSmod_calave.ms' is generated
          from a single dish regression data 'OrionS_rawACSmod' as follows:
          
          default(sdcal)
          sdcal(infile='OrionS_rawACSmod',scanlist=[20,21,22,23],
                calmode='ps',tau=0.09,outfile='temp.asap')
          default(sdcal)
          sdcal(infile='temp.asap',timeaverage=True,
                tweight='tintsys',outfile='temp2.asap')
          sdsave(infile='temp2.asap',outformat='MS2',
                 outfile='OrionS_rawACSmod_calave.ms')
    """
    # Input and output names
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_basictest'
    blrefroot = sdbaseline_unittest_base.datapath+'refblparam'
    tid = None

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)

        default(sdbaseline)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def test000(self):
        """Basic Test 000: default values for all parameters"""
        tid = '000'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        result = sdbaseline(infile=infile, datacolumn=datacolumn,
                             outfile=outfile)
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        row = 3
        pol = 1
        results = self._getStats(outfile, '')
        theresult = None
        for i in range(len(results)):
            if ((results[i]['row'] == int(row)) and (results[i]['pol'] == int(pol))):
                theresult = results[i]
        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }
        self._compareStats(theresult, reference)

    def test001(self):
        """Basic Test 001: simple successful case: blfunc = 'poly', maskmode = 'list' and masklist=[] (no mask)"""
        tid = '001'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        maskmode = 'list'
        blfunc = 'poly'
        spw = '3'
        pol = 'LL'
        overwrite = True
        result = sdbaseline(infile=infile, datacolumn=datacolumn,
                             maskmode=maskmode, blfunc=blfunc, 
                             spw=spw, pol=pol, outfile=outfile,
                             overwrite=overwrite)
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }

        self._compareStats(theresult, reference)
    
    def test002(self):
        """Basic Test 002: simple successful case: blfunc = 'chebyshev', maskmode = 'list' and masklist=[] (no mask)"""
        tid = '002'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        maskmode = 'list'
        blfunc = 'chebyshev'
        spw = '3'
        pol = 'LL'
        overwrite = True
        result = sdbaseline(infile=infile, datacolumn=datacolumn,
                             maskmode=maskmode, blfunc=blfunc, 
                             spw=spw, pol=pol, outfile=outfile,
                             overwrite=overwrite)
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }

        self._compareStats(theresult, reference)
    
    
    def test003(self):
        """Basic Test 003: simple successful case: blfunc = 'cspline', maskmode = 'list' and masklist=[] (no mask)"""
        print("")

        tid = '003'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'  
        maskmode = 'list'
        blfunc = 'cspline'
        overwrite = True
        npiece = 3
        spw='3'
        pol='LL'
        result = sdbaseline(infile=infile, datacolumn=datacolumn,
                             maskmode=maskmode, blfunc=blfunc, 
                             npiece=npiece,spw=spw, 
                             pol=pol,
                             outfile=outfile,overwrite=overwrite)
        
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,msg="The task returned '"+str(result)+"' instead of None")
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid) 
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16685959517745799,
                     'min': -2.5928177833557129,
                     'max': 1.3953156471252441,
                     'median': -0.00089824199676513672,
                     'stddev': 0.16685959517745766,
                    }

        self._compareStats(theresult, reference)

        #***
        #*** check if baseline is subtracted ***
        #***
        # Output MS only has the selected pol, LL
        in_pol=1
        out_pol=0
        sum_pol1=0.0
        sum_square_pol1 = 0.0

        # open the original MS
        tb.open(infile)
        orig_pol1_value = numpy.array(tb.getcell('FLOAT_DATA', int(spw))[in_pol,:])
        tb.close()
        variance_orig_pol1 = numpy.var(orig_pol1_value)
        
        #open the MS after sdbaseline
        tb.open(outfile)
        pol1_value = numpy.array(tb.getcell('FLOAT_DATA', 0)[out_pol,:])
        tb.close()
        variance_pol1 = numpy.var(pol1_value)

        #assert pol1_value < orig_pol1_value
        self.assertTrue((pol1_value<orig_pol1_value).all())
        
        #assert variance of pol1_value < variance of orig_pol1_value
        self.assertLess(variance_pol1**0.5, variance_orig_pol1**0.5)

        #print '1sigma before cspline (pol1)', variance_orig_pol1**0.5 
        #print '1sigma after cspline (pol1)',  variance_pol1**0.5 
        
    
    def test050(self):
        """Basic Test 050: failure case: existing file as outfile with overwrite=False"""
        infile = self.infile
        outfile = 'Dummy_Empty.ms'
        mode = 'list'
        os.mkdir(outfile)
        try:
            result = sdbaseline(infile=infile, outfile=outfile, overwrite=False, maskmode=mode)
        except Exception as e:
            #pos = str(e).find(outfile+' exists.')
            pos = str(e).find("outfile='" + outfile + "' exists, and cannot overwrite it.")
            self.assertNotEqual(pos, -1, msg='Unexpected exception was thrown: %s'%(str(e)))
        finally:
            shutil.rmtree(outfile)

    def test051(self):
        """Basic Test 051: failure case: no data after selection"""
        tid = '051'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        spw = '10' # non-existent IF value
        mode = 'list'
        try:
            sdbaseline(infile=infile, outfile=outfile, spw=spw, maskmode=mode)
        except Exception as e:
            self.assertIn('Spw Expression: No match found for 10,', str(e))


class sdbaseline_maskTest(sdbaseline_unittest_base):
    """
    Tests for various mask selections. No interactive testing.

    List of tests:
    test100 --- with masked ranges at the edges of spectrum. blfunc is cspline.
    test101 --- with masked ranges not touching spectrum edge

    Note: input data is generated from a single dish regression data,
    'OrionS_rawACSmod', as follows:
      default(sdcal)
      sdcal(infile='OrionS_rawACSmod',scanlist=[20,21,22,23],
                calmode='ps',tau=0.09,outfile='temp.asap')
      default(sdcal)
      sdcal(infile='temp.asap',timeaverage=True,
                tweight='tintsys',outfile='temp2.asap')
      sdsave(infile='temp2.asap',outformat='MS2',
                outfile='OrionS_rawACSmod_calave.ms')
    """
    # Input and output names
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_masktest'
    blrefroot = sdbaseline_unittest_base.datapath+'refblparam_mask'
    tid = None

    # Channel range excluding bad edge
    search = [[200,7599]]
    # Baseline channels. should be identical to one selected by 'auto' mode
    blchan0 = [[200,3979],[4152,7599]]
    blchan2 = [[200,2959],[3120,7599]]

    # reference values
    ref_pol0if0 =  {'linemaxpos': 4102.0, 'linesum': 103.81604766845703,
                    'linemax': 1.6280698776245117,
                    'baserms': 0.15021507441997528,
                    'basestd': 0.15022546052932739}
    ref_pol0if2 = {#'linemaxpos': 3045.0, 'linesum': 127.79755401611328,
                   #'linemax': 2.0193681716918945,
                   #'baserms': 0.13134850561618805,
                   #'basestd': 0.1313575953245163}
                   'rms': 0.13134850561618805,
                   'stddev': 0.1313575953245163}
     
    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')




    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def test100(self):
        """Mask Test 100: with masked ranges at the edges of spectrum. blfunc must be cspline."""
        self.tid='100'
        infile = self.infile
        outfile = self.outroot+self.tid+'.ms'
        datacolumn='float_data'
        mode = 'list'
        spw = '2:%s'%(';'.join(map(self._get_range_in_string,self.search)))
        pol = 'RR'
        blfunc = 'cspline'
        npiece = 4

        result = sdbaseline(infile=infile,datacolumn=datacolumn,maskmode=mode,
                             spw=spw,pol=pol,blfunc=blfunc,npiece=npiece,
                             outfile=outfile)
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        # Compare IF2
        testval = self._getStats(filename=outfile, spw='', pol=0, mask=self.search)
        ref100 = {'rms': 0.18957555661537034,
                  'min': -0.48668813705444336,
                  'max': 1.9516196250915527,
                  'median': -0.013428688049316406,
                  'stddev': 0.18957555661537034,
                  'row': 0,
                  'pol': 0}
        self._compareStats(testval[0], ref100)

    def test101(self):
        """Mask Test 101: with masked ranges not touching spectrum edge"""
        self.tid='101'
        infile = self.infile
        outfile = self.outroot+self.tid+'.ms'
        datacolumn='float_data'
        mode = 'list'
        spw = '2:%s'%(';'.join(map(self._get_range_in_string,self.blchan2)))
        pol = 'RR'

        print('spw =', spw)

        result = sdbaseline(infile=infile,datacolumn=datacolumn,maskmode=mode,
                            outfile=outfile,spw=spw,pol=pol)
        # sdbaseline returns None if it runs successfully
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        # Compare IF2
        testval = self._getStats(filename=outfile, spw='', pol=0, mask=self.blchan2)
        self._compareStats(testval[0],self.ref_pol0if2)
        #self._compareBLparam(self.outroot+self.tid+'.asap_blparam.txt',\
        #                     self.blrefroot+self.tid)

    def _get_range_in_string(self, valrange):
        if isinstance(valrange, list) or isinstance(valrange, tuple):
            return str(valrange[0])+'~'+str(valrange[1])
        else:
            return False


class sdbaseline_sinusoidTest(sdbaseline_unittest_base):
    """
    Tests for sinusoidal baseline fitting. No interactive testing.

    List of tests:
    test000 --- addwn as integer
    test001 --- addwn as list of an integer
    test002 --- addwn as list of integers
    test003 --- addwn as tuple of an integer
    test004 --- addwn as tuple of integers
    test005 --- addwn as string (single wave number)
    test006 --- addwn as string (comma-separated wave numbers)
    test007 --- addwn as string (wave number range specified with '-')
    test008 --- addwn as string (wave number range specified with '~')
    test009 --- addwn as string (less or equal pattern 1)
    test010 --- addwn as string (less or equal pattern 2)
    test011 --- addwn as string (less or equal pattern 3)
    test012 --- addwn as string (less or equal pattern 4)
    test013 --- addwn as string (less pattern 1)
    test014 --- addwn as string (less pattern 2)
    test015 --- addwn as string (greater or equal pattern 1)
    test016 --- addwn as string (greater or equal pattern 2)
    test017 --- addwn as string (greater or equal pattern 3)
    test018 --- addwn as string (greater or equal pattern 4)
    test019 --- addwn as string (greater pattern 1)
    test020 --- addwn as string (greater pattern 2)
    test021 --- specify fftthresh by 'sigma' + checking residual rms
    test022 --- specify fftthresh by 'top' + checking residual rms
    test023 --- sinusoid-related parameters with default values
    test024 --- addwn has too large value but rejwn removes it
    
    test100 --- no effective wave number set (addwn empty list, applyfft=False)
    test101 --- no effective wave number set (addwn empty list, applyfft=True)
    test102 --- no effective wave number set (addwn empty tuple, applyfft=False)
    test103 --- no effective wave number set (addwn empty tuple, applyfft=True)
    test104 --- no effective wave number set (addwn empty string, applyfft=False)
    test105 --- no effective wave number set (addwn empty string, applyfft=True)
    test106 --- no effective wave number set (addwn and rejwn identical, applyfft=False)
    test107 --- no effective wave number set (addwn and rejwn identical, applyfft=True)
    test108 --- no effective wave number set (rejwn covers wider range than that of addwn, applyfft=False)
    test109 --- no effective wave number set (rejwn covers wider range than that of addwn, applyfft=True)
    test110 --- wn range greater than upper limit
    test111 --- explicitly specify wn value (greater than upper limit)
    test112 --- explicitly specify wn value (negative)
    test113 --- explicitly specify wn value (addwn has negative and greater than upper limit)
    test114 --- explicitly specify wn value (both addwn/rejwn have negative and greater than upper limit)
    test115 --- wrong fftthresh (as list)
    test116 --- wrong fftthresh (as string 'asigma')
    test117 --- wrong fftthresh (as string 'topa')
    test118 --- wrong fftthresh (as string 'top3sigma')
    test119 --- wrong fftthresh (as string 'a123')
    test120 --- wrong fftthresh (as string '')
    test121 --- wrong fftthresh (as string '-3.0')
    test122 --- wrong fftthresh (as string '0.0')
    test123 --- wrong fftthresh (as string '-3')
    test124 --- wrong fftthresh (as string '0')
    test125 --- wrong fftthresh (as string '-3.0sigma')
    test126 --- wrong fftthresh (as string '0.0sigma')
    test127 --- wrong fftthresh (as string '-3sigma')
    test128 --- wrong fftthresh (as string '0sigma')
    test129 --- wrong fftthresh (as string 'top-3')
    test130 --- wrong fftthresh (as string 'top0')
    test131 --- wrong fftthresh (as string 'top1.5')
    test132 --- wrong fftthresh (as float -3.0)
    test133 --- wrong fftthresh (as float 0.0)
    test134 --- wrong fftthresh (as int -3)
    test135 --- wrong fftthresh (as int 0)

    Note: The input data 'sinusoidal.ms' has just two spectral data,
          which are actually identical and described as 
          spec[i] = sin(i*2*PI/8191) + 4 * sin(i*2*PI/8191*3)
                    + 8 * sin(i*2*PI/8191*5) + 2 * sin(i*2*PI/8191*12).
          addwn='1,3,5,12' will be enough to perfectly fit this spectrum, but
          applyfft=True and fftthresh='top4' will also do.
    """
    # Input and output names
    infile = 'sinusoidal.ms'
    outroot = sdbaseline_unittest_base.taskname + '_sinusoidtest'
    tid = None

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)

        default(sdbaseline)

        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def tearDown(self):
        if (os.path.exists(self.infile)):
            pass
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def test000(self):
        """Sinusoid Test 000: addwn as integer"""
        tid = '000'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = 0
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test001(self):
        """Sinusoid Test 001: addwn as list of an integer"""
        tid = '001'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [0]
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test002(self):
        """Sinusoid Test 002: addwn as list of integers"""
        tid = '002'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [0,1]
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test003(self):
        """Sinusoid Test 003: addwn as tuple of an integer"""
        tid = '003'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = (0)
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test004(self):
        """Sinusoid Test 004: addwn as tuple of integers"""
        tid = '004'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = (0,1)
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test005(self):
        """Sinusoid Test 005: addwn as string (single wave number)"""
        tid = '005'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test006(self):
        """Sinusoid Test 006: addwn as string (comma-separated wave numbers)"""
        tid = '006'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0,1'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test007(self):
        """Sinusoid Test 007: addwn as string (wave number range specified with '-')"""
        tid = '007'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0-2'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test008(self):
        """Sinusoid Test 008: addwn as string (wave number range specified with '~')"""
        tid = '008'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0~2'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test009(self):
        """Sinusoid Test 009: addwn as string (less or equal pattern 1)"""
        tid = '009'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '<=2'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test010(self):
        """Sinusoid Test 010: addwn as string (less or equal pattern 2)"""
        tid = '010'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '=<2'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test011(self):
        """Sinusoid Test 011: addwn as string (less or equal pattern 3)"""
        tid = '011'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '2>='
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test012(self):
        """Sinusoid Test 012: addwn as string (less or equal pattern 4)"""
        tid = '012'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '2=>'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test013(self):
        """Sinusoid Test 013: addwn as string (less pattern 1)"""
        tid = '013'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '<2'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test014(self):
        """Sinusoid Test 014: addwn as string (less pattern 2)"""
        tid = '014'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '2>'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test015(self):
        """Sinusoid Test 015: addwn as string (greater or equal pattern 1)"""
        tid = '015'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '4090<='
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test016(self):
        """Sinusoid Test 016: addwn as string (greater or equal pattern 2)"""
        tid = '016'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '4090=<'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test017(self):
        """Sinusoid Test 017: addwn as string (greater or equal pattern 3)"""
        tid = '017'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '>=4090'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test018(self):
        """Sinusoid Test 018: addwn as string (greater or equal pattern 4)"""
        tid = '018'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '=>4090'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test019(self):
        """Sinusoid Test 019: addwn as string (greater pattern 1)"""
        tid = '019'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '4090<'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test020(self):
        """Sinusoid Test 020: addwn as string (greater pattern 2)"""
        tid = '020'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '>4090'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=False)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

    def test021(self):
        """Sinusoid Test 021: specify fftthresh by 'sigma' + checking residual rms"""
        tid = '021'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0'
        fftthresh = '3.0sigma'
        torr = 1.0e-6
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=True,fftthresh=fftthresh)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        stat = self._getStats(filename=outfile, pol='0')
        self.assertTrue(stat[0]['rms'] < torr)

    def test022(self):
        """Sinusoid Test 022: specify fftthresh by 'top' + checking residual rms"""
        tid = '022'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '0'
        fftthresh = 'top4'
        torr = 1.0e-6
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',addwn=addwn,applyfft=True,fftthresh=fftthresh)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        stat = self._getStats(filename=outfile, pol='0')
        self.assertTrue(stat[0]['rms'] < torr)

    def test023(self):
        """Sinusoid Test 023: sinusoid-related parameters with default values"""
        tid = '023'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid')
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
    def test024(self):
        """Sinusoid Test 024: addwn has too large value but rejwn removes it"""
        tid = '024'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        applyfft = False
        addwn = [0,10000]
        rejwn = '4000<'
        result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                             blfunc='sinusoid',applyfft=applyfft,addwn=addwn,rejwn=rejwn)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        
    def test100(self):
        """Sinusoid Test 100: no effective wave number set (addwn empty list, applyfft=False)"""
        tid = '100'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = []
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=False)
        except Exception as e:
            self.assertEqual(e.message, 'addwn must contain at least one element.')
        
    def test101(self):
        """Sinusoid Test 101: no effective wave number set (addwn empty list, applyfft=True)"""
        tid = '101'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = []
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'addwn must contain at least one element.')
        
    def test102(self):
        """Sinusoid Test 102: no effective wave number set (addwn empty tuple, applyfft=False)"""
        tid = '102'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = ()
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=False)
        except Exception as e:
            self.assertEqual(e.message, 'addwn must contain at least one element.')
        
    def test103(self):
        """Sinusoid Test 103: no effective wave number set (addwn empty tuple, applyfft=True)"""
        tid = '103'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = ()
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'addwn must contain at least one element.')
        
    def test104(self):
        """Sinusoid Test 104: no effective wave number set (addwn empty string, applyfft=False)"""
        tid = '104'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = ''
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=False)
        except Exception as e:
            self.assertEqual(e.message, 'string index out of range')
        
    def test105(self):
        """Sinusoid Test 105: no effective wave number set (addwn empty string, applyfft=True)"""
        tid = '105'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = ''
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'string index out of range')
        
    def test106(self):
        """Sinusoid Test 106: no effective wave number set (addwn and rejwn identical, applyfft=False)"""
        tid = '106'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [0,1,2]
        rejwn = [0,1,2]
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=False)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')
        
    def test107(self):
        """Sinusoid Test 107: no effective wave number set (addwn and rejwn identical, applyfft=True)"""
        tid = '107'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [0,1,2]
        rejwn = [0,1,2]
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')
        
    def test108(self):
        """Sinusoid Test 108: no effective wave number set (rejwn covers wider range than that of addwn, applyfft=False)"""
        tid = '108'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '<5'
        rejwn = '<10'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=False)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')
        
    def test109(self):
        """Sinusoid Test 109: no effective wave number set (rejwn covers wider range than that of addwn, applyfft=True)"""
        tid = '109'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '<5'
        rejwn = '<10'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')

    def test110(self):
        """Sinusoid Test 110: wn range greater than upper limit"""
        tid = '110'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = '5000<'
        rejwn = '<5100'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')

    def test111(self):
        """Sinusoid Test 111: explicitly specify wn value (greater than upper limit)"""
        tid = '111'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [5000,5500]
        rejwn = []
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'No effective wave number given for sinusoidal fitting.')

    def test112(self):
        """Sinusoid Test 112: explicitly specify wn value (negative)"""
        tid = '112'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [-10,5]
        rejwn = []
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'wrong value given for addwn/rejwn')

    def test113(self):
        """Sinusoid Test 113: explicitly specify wn value (addwn has negative and greater than upper limit)"""
        tid = '113'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [-10,5000]
        rejwn = []
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'wrong value given for addwn/rejwn')

    def test114(self):
        """Sinusoid Test 114: explicitly specify wn value (both addwn/rejwn have negative and greater than upper limit)"""
        tid = '114'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        addwn = [-10,5000]
        rejwn = [-10,5500]
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',addwn=addwn,rejwn=rejwn,applyfft=True)
        except Exception as e:
            self.assertEqual(e.message, 'wrong value given for addwn/rejwn')

    def test115(self):
        """Sinusoid Test 115: wrong fftthresh (as list)"""
        tid = '115'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = [3.0]
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh must be float or integer or string.')

    def test116(self):
        """Sinusoid Test 116: wrong fftthresh (as string 'asigma')"""
        tid = '116'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'asigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test117(self):
        """Sinusoid Test 117: wrong fftthresh (as string 'topa')"""
        tid = '117'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'topa'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test118(self):
        """Sinusoid Test 118: wrong fftthresh (as string 'top3sigma')"""
        tid = '118'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'top3sigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test119(self):
        """Sinusoid Test 119: wrong fftthresh (as string 'a123')"""
        tid = '119'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'a123'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test120(self):
        """Sinusoid Test 120: wrong fftthresh (as string '')"""
        tid = '120'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = ''
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test121(self):
        """Sinusoid Test 121: wrong fftthresh (as string '-3.0')"""
        tid = '121'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '-3.0'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test122(self):
        """Sinusoid Test 122: wrong fftthresh (as string '0.0')"""
        tid = '122'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '0.0'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test123(self):
        """Sinusoid Test 123: wrong fftthresh (as string '-3')"""
        tid = '123'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '-3'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test124(self):
        """Sinusoid Test 124: wrong fftthresh (as string '0')"""
        tid = '124'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '0'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test125(self):
        """Sinusoid Test 125: wrong fftthresh (as string '-3.0sigma')"""
        tid = '125'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '-3.0sigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test126(self):
        """Sinusoid Test 126: wrong fftthresh (as string '0.0sigma')"""
        tid = '126'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '0.0sigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test127(self):
        """Sinusoid Test 127: wrong fftthresh (as string '-3sigma')"""
        tid = '127'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '-3sigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test128(self):
        """Sinusoid Test 128: wrong fftthresh (as string '0sigma')"""
        tid = '128'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = '0sigma'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test129(self):
        """Sinusoid Test 129: wrong fftthresh (as string 'top-3')"""
        tid = '129'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'top-3'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test130(self):
        """Sinusoid Test 130: wrong fftthresh (as string 'top0')"""
        tid = '130'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'top0'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test131(self):
        """Sinusoid Test 131: wrong fftthresh (as string 'top1.5')"""
        tid = '131'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 'top1.5'
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'fftthresh has a wrong format.')

    def test132(self):
        """Sinusoid Test 132: wrong fftthresh (as float -3.0)"""
        tid = '132'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = -3.0
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test133(self):
        """Sinusoid Test 133: wrong fftthresh (as float 0.0)"""
        tid = '133'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 0.0
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test134(self):
        """Sinusoid Test 134: wrong fftthresh (as int -3)"""
        tid = '134'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = -3
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')

    def test135(self):
        """Sinusoid Test 135: wrong fftthresh (as int 0)"""
        tid = '135'
        infile = self.infile
        outfile = self.outroot + tid + '.ms'
        datacolumn = 'float_data'
        fftthresh = 0
        try:
            result = sdbaseline(infile=infile,datacolumn=datacolumn,outfile=outfile,
                                 blfunc='sinusoid',applyfft=True,fftthresh=fftthresh)
        except Exception as e:
            self.assertEqual(e.message, 'threshold given to fftthresh must be positive.')


class sdbaseline_multi_IF_test(sdbaseline_unittest_base):
    """
    Unit tests for task sdbaseline. No interactive testing.

    This test intends to check whether sdbaseline task works fine
    for data that has multiple IFs whose nchan differ each other. 

    List of tests:
    test200 --- test multi IF data input
    """
    # Input and output names
    infile = 'testMultiIF.asap'
    blparamfile_suffix = '_blparam.txt'
    outroot = sdbaseline_unittest_base.taskname+'_multi'
    refblparamfile = 'refblparam_multiIF'

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')


    def tearDown(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def test200(self):
        """test200: Test the task works with multi IF data"""
        infile = self.infile
        mode = "list"
        blfunc = "poly"
        order = 1
        outfile = self.outroot+".asap"
        blparamfile = outfile+self.blparamfile_suffix
        
        result = sdbaseline(infile=infile,maskmode=mode,outfile=outfile,blfunc=blfunc,order=order)
        self.assertEqual(result, None, msg="The task returned '"+str(result)+"' instead of None")
        self._compareBLparam(blparamfile,self.datapath+self.refblparamfile)
        reference = {5: {'rms': 1.4250789880752563,
                         'min': -4.2702846527099609,
                         'max': 5.5566844940185547,
                         'max_abscissa': {'value': 823.0,
                                          'unit': 'channel'},
                         'median': 0.017315864562988281,
                         'min_abscissa': {'value': 520.0,
                                          'unit': 'channel'},
                         'stddev': 1.425775408744812},
                     7: {'rms': 1.4971292018890381,
                         'min': -4.7103700637817383,
                         'max': 5.4820127487182617,
                         'max_abscissa': {'value': 1335.0,
                                          'unit': 'channel'},
                         'median': 0.027227401733398438,
                         'min_abscissa': {'value': 1490.0,
                                          'unit': 'channel'},
                         'stddev': 1.4974949359893799}}
        for ifno in [5,7]:
            currstat = self._getStats(outfile,ifno)
            self._compareStats(currstat,reference[ifno])


class sdbaseline_outbltableTest(sdbaseline_unittest_base):
    """
    Tests for outputting baseline table

    List of tests
    test300 --- blmode='fit', bloutput='', dosubtract=False (no baselining, no bltable output)
    test301 --- blmode='fit', bloutput!='', dosubtract=True, blfunc='poly'/'chebyshev'/'cspline'
                (poly/chebyshev/cspline fit in MS, bltable is written)
    test302 --- blmode='fit', bloutput!='', dosubtract=True, blfunc='variable'
                (variable fit in MS, bltable is written)
                testing 3 cases:
                    (1) blparam contains values for all spectra
                    (2) no values for a spectrum (row=2,pol=1), which is to be skipped
                    (3) values commented out for a spectrum (row=2,pol=1), which is to be skipped
    test303 --- blmode='fit', bloutput!='', dosubtract=True, blfunc='poly','chebyshev','cspline'
                testing if bltable is shortened
                testing 3 cases:
                    (1) all spectra in row 2 are flagged entirely
                    (2) in row 2, entirely flagged for pol 0, also pol 1 is unselected
                    (3) in row 2, entirely flagged for pol 1, also pol 0 is unselected
    test304 --- same as test303, but for blfunc='variable'

    Note: input data is generated from a single dish regression data,
    'OrionS_rawACSmod', as follows:
      default(sdcal)
      sdcal(infile='OrionS_rawACSmod',scanlist=[20,21,22,23],
                calmode='ps',tau=0.09,outfile='temp.asap')
      default(sdcal)
      sdcal(infile='temp.asap',timeaverage=True,
                tweight='tintsys',outfile='temp2.asap')
      sdsave(infile='temp2.asap',outformat='MS2',
                outfile='OrionS_rawACSmod_calave.ms')
    """
    # Input and output names
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_bltabletest'
    tid = None
    ftype = {'poly': 0, 'chebyshev': 1, 'cspline': 2, 'sinusoid': 3}

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')


    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def _checkBltableVar(self, outms, bltable, blparam, option):
        npol = 2
        results = [[4.280704], [3.912475],
                   [4.323003, 0.00196013], [3.839441, -8.761247e-06],
                   [4.280719, 0.00255683, 0.00619966], [4.140454, -7.516477e-05, 6.538814e-09],
                   [4.221929, -8.751897e-06, -6.81303991e-09, 3.36383428e-13],
                   [3.983634, -6.322114e-06, -1.11215614e-08, 7.00922610e-13]
                   ]
        rms = [0.162739, 0.182507, 0.140955, 0.159999, 0.132135, 0.381708, 0.128761, 0.146849]
        tb.open(bltable)
        try:
            for i in range(npol*tb.nrows()):
                irow = i / npol
                ipol = i % npol
                is_skipped = (option != '') and (irow == 2) and (ipol == 1)

                self.assertEqual(not is_skipped, tb.getcell('APPLY', irow)[ipol][0]);
                if is_skipped: continue
            
                self.assertEqual(self.ftype[blparam['btype'][i]], tb.getcell('FUNC_TYPE', irow)[ipol][0]);
                fparam_key = 'order' if (blparam['btype'][i] != 'cspline') else 'npiec'
                self.assertEqual(blparam[fparam_key][i], tb.getcell('FUNC_PARAM', irow)[ipol][0])
            
                if (blparam['btype'][i] == 'cspline'):
                    for j in range(blparam['npiec'][i]):
                        self.assertEqual(0.0, tb.getcell('FUNC_FPARAM', irow)[ipol][j])
                else:
                    self.assertEqual(0, len(tb.getcell('FUNC_FPARAM', irow)[ipol]))
                for j in range(len(results[i])):
                    self._checkValue(results[i][j], tb.getcell('RESULT', irow)[ipol][j], 1.0e-5)
                self._checkValue(rms[i], tb.getcell('RMS', irow)[ipol][0], 1.0e-1)
                self._checkValue(float(blparam['cthre'][i]), tb.getcell('CLIP_THRESHOLD', irow)[ipol][0], 1.0e-6)
                self.assertEqual(blparam['nclip'][i], tb.getcell('CLIP_ITERATION', irow)[ipol][0])
                uself = (blparam['uself'][i] == 'true')
                self.assertEqual(uself, tb.getcell('USE_LF', irow)[ipol][0])
                lthre = 5.0 if ((blparam['lthre'][i] == '') or not uself) else float(blparam['lthre'][i])
                self._checkValue(lthre, tb.getcell('LF_THRESHOLD', irow)[ipol][0], 1.0e-6)
                chavg = 0 if (blparam['chavg'][i] == '') else int(blparam['chavg'][i])
                self.assertEqual(chavg, tb.getcell('LF_AVERAGE', irow)[ipol][0])
                ledge = 0 if ((blparam['ledge'][i] == '') or not uself) else int(blparam['ledge'][i])
                self.assertEqual(ledge, tb.getcell('LF_EDGE', irow)[ipol][0])
                redge = 0 if ((blparam['redge'][i] == '') or not uself) else int(blparam['redge'][i])
                self.assertEqual(redge, tb.getcell('LF_EDGE', irow)[ipol][1])
        finally:
            tb.close()
    
    def _checkBltable(self, outms, bltable, blfunc, order, mask):
        tb.open(bltable)
        for irow in range(tb.nrows()):
            for ipol in range(len(tb.getcell('RMS', irow))):
                self.assertEqual(tb.getcell('FUNC_TYPE', irow)[ipol], self.ftype[blfunc])
                self.assertEqual(tb.getcell('FUNC_PARAM', irow)[ipol], order)
                ref = self._getStats(filename=outms, spw=str(irow), pol=str(ipol), mask=mask[irow])
                #tolerance value in the next line is temporarily set a bit large 
                #since rms in bltable is smaller than expected because it is
                #calculated based on masklist currently stored in bltable, which 
                #is after an extra clipping.
                #this bug is already fixed in trunk of Sakura, so once libsakura
                #is updated we can set smaller tolerance value. (2015/4/22 WK)
                self._checkValue(ref[0]['rms'], tb.getcell('RMS', irow)[ipol][0], 2.0e-2)
        tb.close()

    def _checkValue(self, ref, out, tol=1.0e-02):
        #print '###################################'
        #print 'ref = ' + str(ref) + ', out = ' + str(out)
        if (abs(ref) > tol) or (abs(out) > tol):
            if ref != 0.0:
                rel = abs((out - ref)/ref)
            elif out != 0.0:
                rel = abs((out - ref)/out)
            else:
                rel = abs(out - ref)
            self.assertTrue((rel < tol), msg='the output ('+str(out)+') differs from reference ('+str(ref)+')')


    def test300(self):
        """test300: no baselining, no bltable output"""
        self.tid='300'
        infile = self.infile
        outfile = self.outroot+self.tid+'.ms'
        datacolumn='float_data'
        blmode='fit'
        bloutput=''
        dosubtract=False

        result = sdbaseline(infile=infile,datacolumn=datacolumn,
                             blmode=blmode,bloutput=bloutput,dosubtract=dosubtract,
                             outfile=outfile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        spec_in = []
        tb.open(infile)
        for i in range(tb.nrows()):
            spec_in.append(tb.getcell('FLOAT_DATA', i))
        tb.close()
        spec_out = []
        tb.open(outfile)
        for i in range(tb.nrows()):
            spec_out.append(tb.getcell('FLOAT_DATA', i))
        tb.close()
        for irow in range(len(spec_in)):
            for ipol in range(len(spec_in[0])):
                for ichan in range(len(spec_in[0][0])):
                    self.assertEqual(spec_in[irow][ipol][ichan], spec_out[irow][ipol][ichan],
                                     msg="output spectrum modified at row="+str(irow)+
                                     ",pol="+str(ipol)+",chan="+str(ichan))

    def test301(self):
        """test301: poly/chebyshev/cspline baselining, output bltable"""
        self.tid='301'
        infile = self.infile
        datacolumn='float_data'
        spw='0:1000~3500;5000~7500,1:500~7500,2:500~2500;3500~7500'
        mask=[ [[1000,3500],[5000,7500]],
               [[500,7500]],
               [[500,2500],[3500,7500]]
               ]
        blmode='fit'
        blformat='table'
        dosubtract=True
        blfunc=['poly','chebyshev','cspline']
        order=5
        npiece=4
        rms_s0p0_ms = [0.150905484071, 0.150905484071, 0.149185846787]

        for i in range(len(blfunc)):
            print('testing blfunc='+blfunc[i]+'...')
            outfile = self.outroot+self.tid+blfunc[i]+'.ms'
            bloutput= self.outroot+self.tid+blfunc[i]+'.bltable'
            result = sdbaseline(infile=infile,datacolumn=datacolumn,
                                 blmode=blmode,blformat=blformat,bloutput=bloutput,
                                 spw=spw,blfunc=blfunc[i],order=order,npiece=npiece,
                                 dosubtract=dosubtract,outfile=outfile)
            self.assertEqual(result,None,
                             msg="The task returned '"+str(result)+"' instead of None")
            msresult = self._getStats(filename=outfile, spw='0', pol='0', mask=mask[0])
            self._checkValue(rms_s0p0_ms[i], msresult[0]['stddev'], 1.0e-6)

            fparam = npiece if blfunc[i] == 'cspline' else order
            self._checkBltable(outfile, bloutput, blfunc[i], fparam, mask)

    def test302(self):
        """test302: per-spectrum baselining, output bltable"""
        self.tid='302'
        infile = self.infile
        datacolumn='float_data'
        blmode='fit'
        blformat='table'
        blfunc='variable'
        dosubtract=True

        for option in ['', 'r2p1less', 'r2p1cout']:
            bloutput= self.outroot+self.tid+option+'.bltable'
            outfile = self.outroot+self.tid+option+'.ms'
            blparam = self.outroot+self.tid+option+'.blparam'
            self._createBlparamFile(blparam, self.blparam_order, self.blparam_dic, option)
            result = sdbaseline(infile=infile,datacolumn=datacolumn,
                                 blmode=blmode,blformat=blformat,bloutput=bloutput,
                                 blfunc=blfunc,blparam=blparam,
                                 dosubtract=dosubtract,outfile=outfile)
            self.assertEqual(result,None,
                             msg="The task returned '"+str(result)+"' instead of None")
            self._checkBltableVar(outfile, bloutput, self.blparam_dic, option)

    def test303(self):
        """test303: testing shortening baseline table for poly,chebyshev,cspline"""
        self.tid = '303'
        infile = self.infile
        datacolumn='float_data'
        spw=''
        blmode='fit'
        blformat='table'
        dosubtract=True
        blfunc=['poly','chebyshev','cspline']
        order=5
        npiece=4
        with tbmanager(infile) as tb:
            nrow_data = tb.nrows()
        testmode = ['masked_masked', 'masked_unselect', 'unselect_masked']
        prange = [[0,1], [0], [1]]
        polval = ['', 'RR', 'LL']
        for i in range(len(blfunc)):
            for j in range(len(testmode)):
                print('testing blfunc='+blfunc[i]+', testmode='+testmode[j]+'...')
                #prepare input data
                if os.path.exists(infile):
                    shutil.rmtree(infile)
                shutil.copytree(self.datapath+self.infile, infile)

                tb.open(tablename=infile, nomodify=False)
                r2msk = tb.getcell('FLAG', 2)
                for ipol in prange[j]:
                    for ichan in range(len(r2msk[0])):
                        r2msk[ipol][ichan] = True
                tb.putcell('FLAG', 2, r2msk)
                tb.close()
                pol = polval[j]

                outfile = self.outroot+self.tid+blfunc[i]+testmode[j]+'.ms'
                bloutput= self.outroot+self.tid+blfunc[i]+testmode[j]+'.bltable'
                result = sdbaseline(infile=infile,datacolumn=datacolumn,
                                     blmode=blmode,blformat=blformat,bloutput=bloutput,
                                     spw=spw,pol=pol,blfunc=blfunc[i],order=order,npiece=npiece,
                                     dosubtract=dosubtract,outfile=outfile)
                self.assertEqual(result,None,
                                 msg="The task returned '"+str(result)+"' instead of None")
                with tbmanager(bloutput) as tb:
                    nrow_bltable = tb.nrows()
                self.assertTrue((nrow_bltable == nrow_data - 1), 
                                msg="The baseline table is not shortened...")
                #delete used data
                if (os.path.exists(self.infile)):
                    shutil.rmtree(self.infile)
                os.system('rm -rf '+self.outroot+'*')
    
    def test304(self):
        """test304: testing shortening baseline table for blfunc=variable"""
        self.tid = '304'
        infile = self.infile
        datacolumn='float_data'
        spw=''
        blmode='fit'
        blformat='table'
        blfunc='variable'
        dosubtract=True
        with tbmanager(infile) as tb:
            nrow_data = tb.nrows()
        testmode = ['masked_masked', 'masked_unselect', 'unselect_masked']
        prange = [[0,1], [0], [1]]
        polval = ['', 'RR', 'LL']
        for j in range(len(testmode)):
            print('testing blfunc='+blfunc+', testmode='+testmode[j]+'...')
            #prepare input data
            if os.path.exists(self.infile):
                shutil.rmtree(self.infile)
            shutil.copytree(self.datapath+self.infile, self.infile)

            blparam = self.outroot+'.blparam'
            self._createBlparamFile(blparam, self.blparam_order, self.blparam_dic, '')

            tb.open(tablename=infile, nomodify=False)
            r2msk = tb.getcell('FLAG', 2)
            for ipol in prange[j]:
                for ichan in range(len(r2msk[0])):
                    r2msk[ipol][ichan] = True
            tb.putcell('FLAG', 2, r2msk)
            tb.close()
            pol = polval[j]

            outfile = self.outroot+self.tid+blfunc+'.ms'
            bloutput= self.outroot+self.tid+blfunc+'.bltable'
            result = sdbaseline(infile=infile,datacolumn=datacolumn,
                                 blmode=blmode,blformat=blformat,bloutput=bloutput,
                                 spw=spw,pol=pol,blfunc=blfunc,blparam=blparam,
                                 dosubtract=dosubtract,outfile=outfile)
            with tbmanager(bloutput) as tb:
                nrow_bltable = tb.nrows()
            self.assertTrue((nrow_bltable == nrow_data - 1), 
                            msg="The baseline table is not shortened...")
            #delete used data
            if (os.path.exists(self.infile)):
                shutil.rmtree(self.infile)
            os.system('rm -rf '+self.outroot+'*')
    
class sdbaseline_applybltableTest(sdbaseline_unittest_base):
    """
    Tests for applying baseline table
    (blmode='apply' mode)

    List of tests
    test400 --- MS with no all-channel-flagged, bltable with apply=True for all spectra
    test401 --- MS with one spectrum (irow=2,ipol=1) with all channels flagged, while apply=True throughout bltable
    test402 --- MS with no all-channel-flagged, while apply=False for one spectrum (irow=2,ipol=1) in bltable
    test403 --- MS with no all-channel-flagger, while bltable lacks one row (irow=2)

    Note: for tests401-403, the spectrum with all channels flagged, or the corresponding
    data in baseline table has apply=False or is inexist, should not be subtracted baseline.
    """
    # Input and output names
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_bltabletest'
    reffile = outroot+'.ms'
    blmode = 'apply'
    bltable = outroot+'.bltable'
    tid = None
    
    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)
        
       
        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

        
        #create baseline table
        blparam = self.outroot+'.blparam'
        self._createBlparamFile(blparam, self.blparam_order, self.blparam_dic, '')
        result = sdbaseline(infile=self.infile,datacolumn='float_data',
                             blmode='fit',blformat='table',bloutput=self.bltable,
                             blfunc='variable',blparam=blparam,
                             dosubtract=True,outfile=self.reffile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        default(sdbaseline)

    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')

    def _checkResult(self, outfile, option):
        npol = 2
        with tbmanager(outfile) as tb:
            out_spec = tb.getcol('FLOAT_DATA')
            out_flag = tb.getcol('FLAG')
        with tbmanager(self.reffile) as tb:
            ref_spec = tb.getcol('FLOAT_DATA')
            ref_flag = tb.getcol('FLAG')
        with tbmanager(self.infile) as tb:
            in_spec = tb.getcol('FLOAT_DATA')
            in_flag = tb.getcol('FLAG')
            nrow     = tb.nrows()
            nchan    = len(in_spec[0][0])

        for ipol in range(npol):
            for ichan in range(nchan):
                for irow in range(nrow):
                    outspec = out_spec[ipol][ichan][irow]
                    outflag = out_flag[ipol][ichan][irow]
                    if ((option == 'r2p1msflagged') and (irow == 2) and (ipol == 1)) or \
                       ((option == 'r2p1bltnotapply') and (irow == 2) and (ipol == 1)) or \
                       ((option == 'r2p1bltinexist') and (irow == 2)):
                        ansspec = in_spec[ipol][ichan][irow]
                        ansflag = True
                    else:
                        ansspec = ref_spec[ipol][ichan][irow]
                        ansflag = ref_flag[ipol][ichan][irow]

                    self.assertTrue(abs(outspec-ansspec)<1e-6, msg='spec: result != answer')
                    self.assertEqual(outflag, ansflag, msg='flag: result != answer')


    def test400(self):
        """test400: apply baseline table. all bltable entries applied to all MS data."""
        self.tid = '400'
        outfile = self.outroot+self.tid+'.ms'
        result = sdbaseline(infile=self.infile,datacolumn='float_data',
                             blmode=self.blmode,bltable=self.bltable,
                             outfile=outfile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        self._checkResult(outfile, '')

    def test401(self):
        """test401: apply baseline table to MS with a spectrum totally flagged."""
        self.tid = '401'
        outfile = self.outroot+self.tid+'.ms'
        try:
            tb.open(tablename=self.infile, nomodify=False)
            tmpflag = tb.getcell('FLAG', 2)
            for ichan in range(len(tmpflag[0])):
                tmpflag[1][ichan] = True
            tb.putcell('FLAG', 2, tmpflag)
        finally:
            tb.close()
        
        result = sdbaseline(infile=self.infile,datacolumn='float_data',
                             blmode=self.blmode,bltable=self.bltable,
                             outfile=outfile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        self._checkResult(outfile, 'r2p1msflagged')

    def test402(self):
        """test402: apply baseline table containing apply=False data."""
        self.tid = '402'
        outfile = self.outroot+self.tid+'.ms'

        try:
            tb.open(tablename=self.bltable, nomodify=False)
            tmpapply = tb.getcell('APPLY', 2)
            tmpapply[1] = False
            tb.putcell('APPLY', 2, tmpapply)
        finally:
            tb.close()
        
        result = sdbaseline(infile=self.infile,datacolumn='float_data',
                             blmode=self.blmode,bltable=self.bltable,
                             outfile=outfile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        self._checkResult(outfile, 'r2p1bltnotapply')

    def test403(self):
        """test403: apply baseline table lacking data for a spectrum in MS."""
        self.tid = '403'
        outfile = self.outroot+self.tid+'.ms'

        try:
            tb.open(tablename=self.bltable, nomodify=False)
            tb.removerows([2])
            self.assertEqual(tb.nrows(), 3, msg='failed to remove a row in bltable.')
        finally:
            tb.close()
        
        result = sdbaseline(infile=self.infile,datacolumn='float_data',
                             blmode=self.blmode,bltable=self.bltable,
                             outfile=outfile)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        self._checkResult(outfile, 'r2p1bltinexist')

class sdbaseline_variableTest(sdbaseline_unittest_base):
    """
    Tests for blfunc='variable'

    List of tests necessary
    00: test baseline subtraction with variable baseline functions and orders
    01: test skipping rows by comment, i.e., lines start with '#' (rows should be flagged)
    02: test skipping rows by non-existent lines in blparam file (rows should be flagged)
    03: test mask selection
    04: test data selection
    05: test clipping
    06: duplicated fitting parameter in blparam file (the last one is adopted)
    NOT IMPLEMENTED YET
    * test dosubtract = False
    * line finder
    * edge flagging
    """
    outfile='variable_bl.ms'
    column='float_data'
    nspec = 4
    refstat0 = {'max': [0.0]*nspec, 'min': [0.0]*nspec,
                'rms': [0.0]*nspec, 'stddev': [0.0]*nspec}
    
    def setUp(self):
        if hasattr(self, 'infile'):
            self.__refetch_files(self.infile)

        default(sdbaseline)


    def tearDown(self):
        self._remove([self.infile, self.outfile])

    def _refetch_files(self, files, from_dir=None):
        if type(files)==str: files = [files]
        self._remove(files)
        self._copy(files, from_dir)

    def __select_stats(self, stats, idx_list):
        """
        Returns a dictionary with selected elements of statistics
        stats    : a dictionary of statistics
        idx_list : a list of indices to select in stats
        """
        ret_dict = {}
        for key in list(stats.keys()):
            ret_dict[key] = [stats[key][idx] for idx in idx_list]
        return ret_dict

    def _run_test(self, infile, reference, mask=None, rtol=1.e-5, atol=1.e-6, flag_spec=(), **task_param):
        """
        Run sdbaseline with mode='variable' and test output MS.

        infile    : input ms name
        reference : reference statistic values in form {'key': [value0, value1, ...], ...}
        mask      : list of masklist to calculate statistics of output MS (None=use all)
        rtol, atol: relative and absolute tolerance of comparison.
        flag_spec : a list of rowid and polid pair whose spectrum should be flagged in output MS
        **task_param : additional parameters to invoke task. blfunc and outfile are predefined.
        """
        self.infile = infile
        sdbaseline(infile=self.infile,blfunc='variable',outfile=self.outfile,**task_param)
        colname = (task_param['datacolumn'] if 'datacolumn' in task_param else 'data').upper()

        # calculate statistics of valid spectrum. Test flagged spectrum.
        ivalid_spec = 0
        ispec = 0
        stats_list = []
        valid_idx = []
        with tbmanager(self.outfile) as tb:
            for rowid in range(tb.nrows()):
                data = tb.getcell(colname, rowid)
                flag = tb.getcell('FLAG', rowid)
                npol = len(data)
                for polid in range(npol):
                    if (rowid, polid) in flag_spec:
                        # for flagged rows
                        self.assertTrue(flag[polid].all(),
                                        "row=%d, pol=%d should be flagged" % (rowid, polid))
                    else:
                        spec = data[polid,:]
                        masklist = mask[ivalid_spec] if mask is not None else None
                        stats_list.append(self._calc_stats_of_array(spec, masklist))
                        ivalid_spec += 1
                        valid_idx.append(ispec)
                    ispec += 1
        # shrink reference list if # of processed spectra is smaller than reference (selection)
        if len(stats_list) < len(reference[list(reference.keys())[0]]):
            self.assertEqual(len(valid_idx), len(stats_list),
                             "Internal error: len(valid_idx)!=len(stats_list)")
            reference = self.__select_stats(reference, valid_idx)

        currstat = self._convert_statslist_to_dict(stats_list)
        #print("cruustat=%s" % str(currstat))
        self._compareStats(currstat, reference, rtol=1.0e-6, atol=1.0e-6)

    def testVariable00(self):
        """Test blfunc='variable' with variable baseline functions and orders"""
        infile='analytic_variable.ms'
        paramfile='analytic_variable_blparam.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        self._run_test(infile,self.refstat0,blparam=paramfile,datacolumn=self.column)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def testVariable01(self):
        """Test blfunc='variable' with skipping rows by comment ('#') (rows should be flagged)"""
        infile='analytic_variable.ms'
        paramfile='analytic_variable_blparam_comment.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        self._run_test(infile,self.refstat0,flag_spec=[(0,0)],blparam=paramfile,datacolumn=self.column)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def testVariable02(self):
        """Test blfunc='variable' with non-existent lines in blparam file (rows should be flagged)"""
        infile='analytic_variable.ms'
        paramfile='analytic_variable_blparam_2lines.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        self._run_test(infile,self.refstat0,flag_spec=[(0,0),(1,1)],blparam=paramfile,datacolumn=self.column)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def testVariable03(self):
        """Test blfunc='variable' with mask selection"""
        infile='analytic_order3_withoffset.ms'
        paramfile='analytic_variable_blparam_mask.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        mask = [[[0,4000],[6000,8000]], [[0,5000],[6000,8000]], [[0,3000],[5000,8000]], None]
        self._run_test(infile,self.refstat0,mask=mask,blparam=paramfile,datacolumn=self.column)

    
        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')


    def testVariable04(self):
        """Test blfunc='variable' with data selection (spw='1')"""
        infile='analytic_variable.ms'
        paramfile='analytic_variable_blparam_spw1.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        self._run_test(infile,self.refstat0,spw='1',blparam=paramfile,datacolumn=self.column)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')


    def testVariable05(self):
        """Test blfunc='variable' with clipping"""
        infile='analytic_order3_withoffset.ms'
        paramfile='analytic_variable_blparam_clip.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        mask = [[[0,4000],[6000,8000]], [[0,5000],[6000,8000]], [[0,3000],[5000,8000]], None]
        self._run_test(infile,self.refstat0,atol=1.e-5,
                       mask=mask,blparam=paramfile,datacolumn=self.column)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def testVariable06(self):
        """Test blfunc='variable' with duplicated fitting parameters (the last one is adopted)"""
        infile='analytic_variable.ms'
        paramfile='analytic_variable_blparam_duplicate.txt'
        self._refetch_files([infile, paramfile], self.datapath)
        self._run_test(infile,self.refstat0,blparam=paramfile,datacolumn=self.column)



        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')



class sdbaseline_bloutputTest(sdbaseline_unittest_base):
    """
Basic unit tests for task sdbaseline. No interactive testing.

    List of tests:
    #'poly'
    test000 --- blformat=['csv','text','table'], bloutput=['test.csv','test.txt','test.table']
    test001 --- blformat=['text','csv','table'], bloutput=['test.txt','test.csv','test.table'] 
    test002 --- blformat=['table','text','csv'], bloutput=['test.table','test.txt','test.csv']
    test003 --- blformat=['table','text','csv'], bloutput=['','test.txt','test.csv'] 
    test004 --- blformat=['table','text'],       bloutput=['','','']
    test005 --- blformat=['table','text'],       bloutput=['','']
    test006 --- blformat=['table'],              bloutput=[''] 
    test007 --- blformat=['csv'],                bloutput=['']
    test008 --- blformat=['text'],               bloutput=[''] 
    test009 --- blformat=[''],                   bloutput=['']
    test010 --- blformat=['','csv'],             bloutput=['','test.csv'] 
    test010a --- blformat=['','csv'],             bloutput=['test.text',''] 
    test011 --- blformat='',                     bloutput='' 
    test012 --- blformat='',                     bloutput='test.csv'


    #'cspline'
    test016 --- blformat=['csv','text','table'], bloutput=['test.csv','test.txt','test.table']
    test017 --- blformat=['text','csv','table'], bloutput=['test.txt','test.csv','test.table'] 
    test018 --- blformat=['table','text','csv'], bloutput=['test.table','test.txt','test.csv']
    test019 --- blformat=['table','text','csv'], bloutput=['','test.txt','test.csv'] 
    test020 --- blformat=['table','text'],       bloutput=['','','']
    test021 --- blformat=['table','text'],       bloutput=['','']
    test022 --- blformat=['table'],              bloutput=[''] 
    test023 --- blformat=['csv'],                bloutput=['']
    test024 --- blformat=['text'],               bloutput=[''] 
    test025 --- blformat=[''],                   bloutput=['']
    test026 --- blformat=['','csv'],             bloutput=['','test.csv'] 
    test027 --- blformat='',                     bloutput='' 
    test028 --- blformat='',                     bloutput='test.csv'


    #'variable'
    test013 --- blformat=['csv','text','table'], bloutput=['test.csv','test.txt','test.table'] 
    test014 --- blformat=['table','text','csv'], bloutput=['test.table','','test.csv'] 
    test015 --- blformat=['table','text','csv'], bloutput=['test.table','test.txt','']

    # 'variable'
    test029 --- blformat=['csv','text','table'], bloutput=['test.csv','test.txt','test.table']
    test030 --- blformat=['text','csv','table'], bloutput=['test.txt','test.csv','test.table'] 
    test031 --- blformat=['table','text','csv'], bloutput=['test.table','test.txt','test.csv']
    test032 --- blformat=['table','text','csv'], bloutput=['','test.txt','test.csv'] 
    test033 --- blformat=['table','text'],       bloutput=['','','']
    test034 --- blformat=['table','text'],       bloutput=['','']
    test035 --- blformat=['table'],              bloutput=[''] 
    test036 --- blformat=['csv'],                bloutput=['']
    test037 --- blformat=['text'],               bloutput=[''] 
    test038 --- blformat=[''],                   bloutput=['']
    test039 --- blformat=['','csv'],             bloutput=['','test.csv'] 
    test040 --- blformat='',                     bloutput='' 
    test041 --- blformat='',                     bloutput='test.csv'



    # 'sinusoid'
    test042 --- blformat=['csv','text','table'], bloutput=['test.csv','test.txt','test.table']
    test043 --- blformat=['text','csv','table'], bloutput=['test.txt','test.csv','test.table'] 
    test044 --- blformat=['table','text','csv'], bloutput=['test.table','test.txt','test.csv']
    test045 --- blformat=['table','text','csv'], bloutput=['','test.txt','test.csv'] 
    test046 --- blformat=['table','text'],       bloutput=['','','']
    test047 --- blformat=['table','text'],       bloutput=['','']
    test048 --- blformat=['table'],              bloutput=[''] 
    test049 --- blformat=['csv'],                bloutput=['']
    test050 --- blformat=['text'],               bloutput=[''] 
    test051 --- blformat=[''],                   bloutput=['']
    test052 --- blformat=['','csv'],             bloutput=['','test.csv'] 
    test053 --- blformat='',                     bloutput='' 
    test054 --- blformat='',                     bloutput='test.csv'



    """

    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_bloutputtest'
    outfile = "test.ms"
    bloutput = "test.txt"
    blparam = 'analytic_variable_blparam.txt'
    bloutput_poly_txt ='bloutput_poly.txt'
    bloutput_poly_csv ='bloutput_poly.csv'
    bloutput_cspline_txt ='bloutput_cspline.txt'
    bloutput_cspline_csv ='bloutput_cspline.csv'
    bloutput_variable_txt ='bloutput_variable.txt'
    bloutput_variable_csv ='bloutput_variable.csv'
    blfunc ='poly'
    bloutput_sinusoid_txt ='bloutput_sinusoid.txt'
    bloutput_sinusoid_addwn012_rejwn0_txt = 'bloutput_sinusoid_addwn012_rejwn0.txt'
    bloutput_sinusoid_addwn012_rejwn02_txt = 'bloutput_sinusoid_addwn012_rejwn02.txt'
    bloutput_sinusoid_addwn012_rejwn1_txt = 'bloutput_sinusoid_addwn012_rejwn1.txt'

    bloutput_sinusoid_csv ='bloutput_sinusoid.csv'
    bloutput_sinusoid_addwn012_rejwn0_csv ='bloutput_sinusoid_addwn012_rejwn0.csv'
    bloutput_sinusoid_addwn012_rejwn02_csv ='bloutput_sinusoid_addwn012_rejwn02.csv'
    bloutput_sinusoid_addwn012_rejwn1_csv ='bloutput_sinusoid_addwn012_rejwn1.csv'

    bloutput_sinusoid_addwnGt4000_rejwn4005_txt ='bloutput_sinusoid_addwnGt4000_rejwn4005.txt'


    base_param = dict(infile=infile,
                      blfunc=blfunc,
                      datacolumn='float_data',
                      maskmode = 'list',
                      outfile=outfile,
                      blformat='text',
                      blparam=blparam,
                      bloutput=bloutput)


    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        
        if os.path.exists(self.blparam):
            #shutil.rmtree(self.blparam)
            os.system('rm '+ self.blparam)
        shutil.copyfile(self.datapath+self.blparam, self.blparam)
        shutil.copyfile(self.datapath+self.bloutput_poly_txt, self.bloutput_poly_txt)
        shutil.copyfile(self.datapath+self.bloutput_poly_csv, self.bloutput_poly_csv)        
        shutil.copyfile(self.datapath+self.bloutput_cspline_txt, self.bloutput_cspline_txt)
        shutil.copyfile(self.datapath+self.bloutput_cspline_csv, self.bloutput_cspline_csv)
        shutil.copyfile(self.datapath+self.bloutput_variable_txt, self.bloutput_variable_txt)
        shutil.copyfile(self.datapath+self.bloutput_variable_csv, self.bloutput_variable_csv)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_csv, self.bloutput_sinusoid_csv)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_txt, self.bloutput_sinusoid_txt)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn0_txt, self.bloutput_sinusoid_addwn012_rejwn0_txt)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn02_txt, self.bloutput_sinusoid_addwn012_rejwn02_txt)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn1_txt, self.bloutput_sinusoid_addwn012_rejwn1_txt)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn0_csv, self.bloutput_sinusoid_addwn012_rejwn0_csv)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn02_csv, self.bloutput_sinusoid_addwn012_rejwn02_csv)
        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwn012_rejwn1_csv, self.bloutput_sinusoid_addwn012_rejwn1_csv)

        shutil.copyfile(self.datapath+self.bloutput_sinusoid_addwnGt4000_rejwn4005_txt, self.bloutput_sinusoid_addwnGt4000_rejwn4005_txt)

        default(sdbaseline)


        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.bltable'):
            shutil.rmtree(self.infile+ '_blparam.bltable')

        if os.path.exists('test.txt'):
            os.remove('test.txt')
        if os.path.exists('test.csv'):
            os.remove('test.csv')
        if os.path.exists('test.table'):
            shutil.rmtree('test.table')


    def tearDown(self):
        if (os.path.exists(self.infile)):
            shutil.rmtree(self.infile)
        os.system('rm -rf '+self.outroot+'*')
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile)
        #print 'test'
        

    def run_test(self, **kwargs):
        task_param=self.base_param.copy()
        for key, value in list(kwargs.items()):
            task_param[key] = value
        result = sdbaseline(**task_param)


    def check_bloutput(self,bloutput):
        for fname in bloutput:
            if fname !='':
                 result_exist = os.path.exists(fname)
                 self.assertEqual(result_exist, True, msg=fname + 'does not exist!')


    def check_bloutputparam_csv(self,bloutputfile, ref_all):
        with open(bloutputfile,'r') as file:
            dataReader=csv.reader(file)     
            list_all=[]
            for row in dataReader:
                list_all.append(row)
 
            self.assertEqual(ref_all, list_all, msg='parameter values of the output csv file are not equivalent to referece values!')



    def test000(self):
        """Bloutput Test 000:blfunc='poly',blformat=['csv','text','table'],bloutput=['test.csv','test.txt','test.table']""" 
        blfunc='poly'
        blformat=['csv','text','table']
        bloutput=['test.csv','test.txt','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_poly_csv) 

        diff_value=os.system('diff test.txt ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_poly_txt) 




    def test001(self):
        """Bloutput Test 001: blfunc='poly', blformat=['text','csv','table'], bloutput=['test.txt','test.csv','test.table']"""
        blfunc='poly'
        blformat=['text','csv','table']
        bloutput=['test.txt','test.csv','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
       
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_poly_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_poly_txt)   


    def test002(self):
        """Bloutput Test 002: blfunc='poly', blformat=['table','text','csv'], bloutput=['test.table','test.txt','test.csv']"""
        blfunc='poly'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + ' is not equivalent to ' + self.bloutput_poly_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_poly_txt)   


    def test003(self):
        """Bloutput Test 003: blfunc='poly', blformat=['table','text','csv'], bloutput=['','test.txt','test.csv']"""
        blfunc='poly'
        blformat=['table','text','csv']
        bloutput=['','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_poly_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_poly_txt)   


    def test004(self):
        """Bloutput Test 004: blfunc='poly', blformat=['table','text','csv'], bloutput=['','','']"""
        blfunc='poly'
        blformat=['table','text','csv']
        bloutput=['','','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        #self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_poly_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + ' is not equivalent to ' + self.bloutput_poly_txt)   

        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + ' is not equivalent to ' + self.bloutput_poly_csv)   


    def test005(self):
        """Bloutput Test 005: blfunc='poly', blformat=['table','text'], bloutput=['','']"""
        blfunc='poly'    
        blformat=['table','text']
        bloutput=['','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=self.infile+'_blparam.txt' + ' is not equivalent to ' + self.bloutput_poly_txt)   


    def test006(self):
        """Bloutput Test 006: blfunc='poly', blformat=['table'], bloutput=['']"""
        blfunc='poly'
        blformat=['table']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')


    def test007(self):
        """Bloutput Test 007: blfunc='poly', blformat=['csv'], bloutput=['']"""
        blfunc='poly'
        blformat=['csv']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + ' is not equivalent to ' + self.bloutput_poly_csv)   


    def test008(self):
        """Bloutput Test 008: blfunc='poly', blformat=['text'], bloutput=['']"""
        blfunc='poly'
        blformat=['text']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + ' is not equivalent to ' + self.bloutput_poly_txt)   


    def test009(self):
        """Bloutput Test 009: blfunc='poly', blformat=[''], bloutput=['']"""
        blfunc='poly'
        blformat=['']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')


    def test010(self):
        """Bloutput Test 010: default values for all parameters except blformat=['','csv'] and bloutput=['','test.csv']"""
        blfunc='poly'
        blformat=['','csv']
        bloutput=['','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
    
        self.assertEqual(os.path.exists('test.csv'), True, msg='test.csv exists!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_poly_csv)   


    def test010a(self):
        """Bloutput Test 010a: default values for all parameters except blformat=['','csv'] and bloutput=['test.text','']"""
        blfunc='poly'
        blformat=['','csv']
        bloutput=['test.text','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
    
        self.assertEqual(os.path.exists('test.text'), False, msg='test.text exists!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + ' is not equivalent to ' + self.bloutput_poly_csv)   


    def test011(self):
        """Bloutput Test 011: default values for all parameters except blformat='' and  bloutput=''"""

        blfunc='poly'
        blformat=''
        bloutput=''

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does exist!')


    def test012(self):
        """Bloutput Test 012: default values for all parameters except blformat='' and  bloutput='test.csv'"""

        blfunc='poly'
        blformat=''
        bloutput='test.csv'

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
       
        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does exist!')
        self.assertEqual(os.path.exists(self.bloutput), False, msg= bloutput + ' exist!')


    def test013(self):
        """Bloutput Test 013: default values for all parameters except blfunc=variable, blparam='analytic_variable_blparam.txt', blformat=['csv','text','table'], and bloutput=['test.csv','test.txt','test.table']"""
        blfunc='variable'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')
        
        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_variable_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_variable_csv)   

        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_variable_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test014(self):
        """Bloutput Test 014: default values for all parameters except  blformat=['table','text','csv'] and  bloutput=['test.table','','test.csv']"""
        blfunc='poly'
        blformat=['table','text','csv']
        bloutput=['test.table','','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')
        
        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_poly_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + ' is not equivalent to ' + self.bloutput_poly_txt)   
        
        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_poly_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_poly_csv)   


    def test015(self):
        """Bloutput Test 015: default values for all parameters except blformat=['table','text','csv'] and bloutput=['test.table','test.txt','']"""
        
        blfunc='poly'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')
        
        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_poly_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + ' is not equivalent to ' + self.bloutput_poly_csv)   
        
        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_poly_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_poly_txt)   
        

    def test016(self):
        """Bloutput Test 016: default values for all parameters except blfunc='cspline'"""
       
        blfunc='cspline'
        blformat=['csv','text','table']
        bloutput=['test.csv','test.txt','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_cspline_csv) 

        diff_value=os.system('diff test.txt ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_cspline_txt) 
       
       
    def test017(self):
        """Bloutput Test 017: default values for all parameters except blfunc='cspline',blformat=['text','csv','table'] and bloutput=['test.txt','test.csv','test.table']"""
        blfunc='cspline'
        blformat=['text','csv','table']
        bloutput=['test.txt','test.csv','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_cspline_csv) 

        diff_value=os.system('diff test.txt ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_cspline_txt) 


    def test018(self):
        """Bloutput Test 018: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['test.table','test.txt','test.csv']"""
        blfunc='cspline'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + ' is not equivalent to ' + self.bloutput_cspline_csv) 

        diff_value=os.system('diff test.txt ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_cspline_txt) 


    def test019(self):
        """Bloutput Test 019: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['','test.txt','test.csv']"""
        blfunc='cspline'
        blformat=['table','text','csv']
        bloutput=['','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff test.csv ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_cspline_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_cspline_txt)   


    def test020(self):
        """Bloutput Test 020: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['','','']"""

        blfunc='cspline'
        blformat=['table','text','csv']
        bloutput=['','','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[0] + ' is not equivalent to ' + self.bloutput_cspline_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + ' is not equivalent to ' + self.bloutput_cspline_txt)   


    def test021(self):
        """Bloutput Test 021: default values for all parameters except blfunc='cspline',blformat=['table','text'] and bloutput=['','']"""
        blfunc='cspline'    
        blformat=['table','text']
        bloutput=['','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=self.infile+'_blparam.txt' + ' is not equivalent to ' + self.bloutput_cspline_txt)   


    def test022(self):
        """Bloutput Test 022: default values for all parameters except blfunc='cspline',blformat=['table'] and bloutput=['']"""
        blfunc='cspline'
        blformat=['table']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')


    def test023(self):
        """Bloutput Test 023: default values for all parameters except blfunc='cspline',blformat=['csv'] and bloutput=['']"""

        blfunc='cspline'
        blformat=['csv']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv' + ' ' + self.bloutput_cspline_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + ' is not equivalent to ' + self.bloutput_cspline_csv)   


    def test024(self):
        """Bloutput Test 024: default values for all parameters except blfunc='cspline',blformat=['text'] and bloutput=['']"""

        blfunc='cspline'
        blformat=['text']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt' + ' ' + self.bloutput_cspline_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + ' is not equivalent to ' + self.bloutput_cspline_txt)   


    def test025(self):
        """Bloutput Test 025: default values for all parameters except blfunc='cspline',blformat=[''] and bloutput=['']"""

        blfunc='cspline'
        blformat=['']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')


    def test026(self):
        """Bloutput Test 026: default values for all parameters except blfunc='cspline',blformat=['','csv'] and bloutput=['','test.csv']"""

        blfunc='cspline'
        blformat=['','csv']
        bloutput=['','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
    
        self.assertEqual(os.path.exists('test.csv'), True, msg='test.csv exists!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_cspline_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_cspline_csv)   


    def test027(self):
        """Bloutput Test 027: default values for all parameters except blfunc='cspline',blformat='' and  bloutput=''"""
        blfunc='cspline'
        blformat=''
        bloutput=''

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does exist!')


    def test028(self):
        """Bloutput Test 028: default values for all parameters except blfunc='cspline', blformat='' and  bloutput='test.csv'"""

        blfunc='cspline'
        blformat=''
        bloutput='test.csv'

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")

        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
       
        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.text')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does exist!')
        self.assertEqual(os.path.exists(self.bloutput), False, msg= bloutput + ' exist!')


    def test029(self):
        """Bloutput Test 029: default values for all parameters except blfunc='cspline'"""

        blfunc='variable'
        blformat=['csv','text','table']
        bloutput=['test.csv','test.txt','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_variable_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_variable_csv)   

        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_variable_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test030(self):
        """Bloutput Test 030: default values for all parameters except blfunc='cspline',blformat=['text','csv','table'] and bloutput=['test.txt','test.csv','test.table']"""

        blfunc='variable'
        blformat=['text','csv','table']
        bloutput=['test.txt','test.csv','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_variable_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_csv)   


        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_variable_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test031(self):
        """Bloutput Test 031: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['test.table','test.txt','test.csv']"""

        blfunc='variable'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_variable_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_csv)   


        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_variable_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test032(self):
        """Bloutput Test 032: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['','test.txt','test.csv']"""

        blfunc='variable'
        blformat=['table','text','csv']
        bloutput=['','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        diff_value=os.system('diff test.csv ' + self.bloutput_variable_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_variable_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_variable_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test033(self):
        """Bloutput Test 033: default values for all parameters except blfunc='cspline',blformat=['table','text','csv'] and bloutput=['','','']"""

        blfunc='variable'
        blformat=['table','text','csv']
        bloutput=['','','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test034(self):
        """Bloutput Test 034: default values for all parameters except blfunc='cspline',blformat=['table','text'] and bloutput=['','']"""

        blfunc='variable'
        blformat=['table','text']
        bloutput=['','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        #result_exist = os.path.exists(self.infile + '_blparam.csv')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test035(self):
        """Bloutput Test 035: default values for all parameters except blfunc='cspline',blformat=['table'] and bloutput=['']"""

        blfunc='variable'
        blformat=['table']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test036(self):
        """Bloutput Test 036: default values for all parameters except blfunc='cspline',blformat=['csv'] and bloutput=['']"""

        blfunc='variable'
        blformat=['csv']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test037(self):
        """Bloutput Test 037: default values for all parameters except blfunc='cspline',blformat=['text'] and bloutput=['']"""

        blfunc='variable'
        blformat=['text']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   

    def test038(self):
        """Bloutput Test 038: default values for all parameters except blfunc='cspline',blformat=[''] and bloutput=['']"""

        blfunc='variable'
        blformat=['']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   



    def test039(self):
        """Bloutput Test 039: default values for all parameters except blfunc='cspline',blformat=['','csv'] and bloutput=['','test.csv']"""

        blfunc='variable'
        blformat=['','csv']
        bloutput=['','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        #result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        #result_exist = not os.path.exists(self.infile + '_blparam.txt')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


        result_exist = os.path.exists(bloutput[1])
        self.assertEqual(result_exist, True, msg= bloutput[1] + ' does not exist!')
        diff_value=os.system('diff ' + bloutput[1] + ' ' + self.bloutput_variable_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_variable_csv)   



    def test040(self):
        """Bloutput Test 040: default values for all parameters except blfunc='cspline',blformat='' and  bloutput=''"""

        blfunc='variable'
        blformat=''
        bloutput=''

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


        #result_exist = os.path.exists(bloutput[1])
        #self.assertEqual(result_exist, True, msg= bloutput[1] + ' does not exist!')
        #diff_value=os.system('diff ' + bloutput[1] + ' ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_variable_csv)   

    def test041(self):
        """Bloutput Test 041: default values for all parameters except blfunc='cspline', blformat='' and  bloutput='test.csv'"""

        blfunc='variable'
        blformat=''
        bloutput='test.csv'

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')






    
    def test042(self):
        """Basic Test 042: default values for all parameters except blfunc='sinusoid'"""

        blfunc='sinusoid'
        blformat=['csv','text','table']
        bloutput=['test.csv','test.txt','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_sinusoid_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_sinusoid_csv)   

        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_sinusoid_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test043(self):
        """Basic Test 043: default values for all parameters except blfunc='sinusoid',blformat=['text','csv','table'] and bloutput=['test.txt','test.csv','test.table']"""

        blfunc='sinusoid'
        blformat=['text','csv','table']
        bloutput=['test.txt','test.csv','test.table']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_sinusoid_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_sinusoid_csv)   


        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_sinusoid_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test044(self):
        """Basic Test 044: default values for all parameters except blfunc='sinusoid',blformat=['table','text','csv'] and bloutput=['test.table','test.txt','test.csv']"""

        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['test.table','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        if(os.path.exists('test.csv') == True):
            diff_value=os.system('diff test.csv ' + self.bloutput_sinusoid_csv)
            self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_variable_csv)   


        if(os.path.exists('test.txt') == True):
            diff_value=os.system('diff test.txt ' + self.bloutput_sinusoid_txt)
            self.assertEqual(diff_value, 0, msg=bloutput[0] + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test045(self):
        """Basic Test 045: default values for all parameters except blfunc='sinusoid',blformat=['table','text','csv'] and bloutput=['','test.txt','test.csv']"""

        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['','test.txt','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')


        diff_value=os.system('diff test.csv ' + self.bloutput_sinusoid_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[2] + 'is not equivalent to ' + self.bloutput_sinusoid_csv)   

        diff_value=os.system('diff test.txt ' + self.bloutput_sinusoid_txt)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test046(self):
        """Basic Test 046: default values for all parameters except blfunc='sinusoid',blformat=['table','text','csv'] and bloutput=['','','']"""

        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['','','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_sinusoid_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_sinusoid_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test047(self):
        """Basic Test 047: default values for all parameters except blfunc='sinusoid',blformat=['table','text'] and bloutput=['','']"""

        blfunc='sinusoid'
        blformat=['table','text']
        bloutput=['','']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        #result_exist = os.path.exists(self.infile + '_blparam.csv')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   


    def test048(self):
        """Basic Test 048: default values for all parameters except blfunc='sinusoid',blformat=['table'] and bloutput=['']"""

        blfunc='sinusoid'
        blformat=['table']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_variable_csv)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test049(self):
        """Basic Test 049: default values for all parameters except blfunc='sinusoid',blformat=['csv'] and bloutput=['']"""

        blfunc='sinusoid'
        blformat=['csv']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' does not exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_sinusoid_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_sinusoid_csv)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


    def test050(self):
        """Basic Test 050: default values for all parameters except blfunc='sinusoid',blformat=['text'] and bloutput=['']"""

        blfunc='sinusoid'
        blformat=['text']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' does not exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_txt)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   

    def test051(self):
        """Basic Test 051: default values for all parameters except blfunc='sinusoid',blformat=[''] and bloutput=['']"""

        blfunc='sinusoid'
        blformat=['']
        bloutput=['']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   



    def test052(self):
        """Basic Test 052: default values for all parameters except blfunc='sinusoid',blformat=['','csv'] and bloutput=['','test.csv']"""

        blfunc='sinusoid'
        blformat=['','csv']
        bloutput=['','test.csv']

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)
        
        #result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        #result_exist = not os.path.exists(self.infile + '_blparam.txt')
        #self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


        result_exist = os.path.exists(bloutput[1])
        self.assertEqual(result_exist, True, msg= bloutput[1] + ' does not exist!')
        diff_value=os.system('diff ' + bloutput[1] + ' ' + self.bloutput_sinusoid_csv)
        self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_sinusoid_csv)   



    def test053(self):
        """Basic Test 053: default values for all parameters except blfunc='sinusoid',blformat='' and  bloutput=''"""

        blfunc='sinusoid'
        blformat=''
        bloutput=''

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

        #diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_variable_txt)
        #self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_variable_txt)   


        #result_exist = os.path.exists(bloutput[1])
        #self.assertEqual(result_exist, True, msg= bloutput[1] + ' does not exist!')
        #diff_value=os.system('diff ' + bloutput[1] + ' ' + self.bloutput_variable_csv)
        #self.assertEqual(diff_value, 0, msg=bloutput[1] + ' is not equivalent to ' + self.bloutput_variable_csv)   

    def test054(self):
        """Basic Test 054: default values for all parameters except blfunc='sinusoid', blformat='' and  bloutput='test.csv'"""

        blfunc='sinusoid'
        blformat=''
        bloutput='test.csv'

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
        self.assertEqual(result,None,
                         msg="The task returned '"+str(result)+"' instead of None")
        
        #if len(blformat)==len(bloutput):
        #    self.check_bloutput(bloutput)
        
        result_exist = not os.path.exists(self.infile + '_blparam.bltable')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.bltable'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.txt')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.txt'+' exist!')

        result_exist = not os.path.exists(self.infile + '_blparam.csv')
        self.assertEqual(result_exist, True, msg=self.infile + '_blparam.csv'+' exist!')

    

    def test0123(self):
        """Basic Test 0123: failure test"""
        blfunc='sinusoid'
        blformat=''
        bloutput='test.csv'
        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput)
                
    def test0124(self):
        """Basic Test 0124: addwn012, rejwn0 test"""
        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['','','']
        addwn=[0,1,2]
        rejwn=[0]

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput, addwn=addwn, rejwn=rejwn)
        self.assertEqual(result,None, msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_sinusoid_addwn012_rejwn0_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn0_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_addwn012_rejwn0_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn0_txt)

    def test0125(self):
        """Basic Test 0125: addwn012, rejwn02 test"""
        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['','','']
        addwn=[0,1,2]
        rejwn=[0,2]

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput, addwn=addwn, rejwn=rejwn)
        self.assertEqual(result,None, msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_sinusoid_addwn012_rejwn02_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn02_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_addwn012_rejwn02_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn02_txt)

    def test0126(self):
        """Basic Test 0126: addwn012, rejwn1 test"""
        blfunc='sinusoid'
        blformat=['table','text','csv']
        bloutput=['','','']
        addwn=[0,1,2]
        rejwn=[1]

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput, addwn=addwn, rejwn=rejwn)
        self.assertEqual(result,None, msg="The task returned '"+str(result)+"' instead of None")
        
        if len(blformat)==len(bloutput):
            self.check_bloutput(bloutput)

        diff_value=os.system('diff ' + self.infile + '_blparam.csv ' + self.bloutput_sinusoid_addwn012_rejwn1_csv)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.csv' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn1_csv)   

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_addwn012_rejwn1_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_addwn012_rejwn1_txt)


    def test0127(self):
        """Basic Test 0127: addwn>4000, rejwn4005 test"""
        
        blfunc='sinusoid'
        blformat=['text','csv']
        bloutput=['','']
        addwn='>4000'
        rejwn=[4005]
        spw='0'
        applyfft = False

        result = self.run_test(blfunc=blfunc, blformat=blformat, bloutput=bloutput, addwn=addwn, rejwn=rejwn, spw=spw, applyfft=applyfft)
        self.assertEqual(result,None, msg="The task returned '"+str(result)+"' instead of None")

        diff_value=os.system('diff ' + self.infile + '_blparam.txt ' + self.bloutput_sinusoid_addwnGt4000_rejwn4005_txt)
        self.assertEqual(diff_value, 0, msg=self.infile + '_blparam.txt' + 'is not equivalent to ' + self.bloutput_sinusoid_addwnGt4000_rejwn4005_txt)







class sdbaseline_autoTest(sdbaseline_unittest_base):
    """
    A class that tests maskmode='auto'.
    
    testAutoPolyNoMask : polynomial fitting using all channels but edge=(500, 500)
    testAutoChebNoMask : Chebyshev polynomial fitting using all channels but edge=(500, 500)
    testAutoCsplNoMask : cspline fitting using all channels but edge=(500, 500)
    testAutoSinuNoMask : sinusoidal fitting using all channels but edge=(500, 500)
    testAutoPolyMaskChan : polynomial fitting using 500~7691 channels (no edge mask)
    testAutoChebMaskChan : Chebyshev polynomial fitting using 500~7691 channels (no edge mask)
    testAutoCsplMaskChan : cspline fitting using 500~7691 channels (no edge mask)
    testAutoSinuMaskChan : sinusoidal fitting using 500~7691 channels (no edge mask)
    testAutoPolyMaskFreq : polynomial fitting using 500~7691 (no edge mask)
    testAutoChebMaskFreq : Chebyshev polynomial fitting using 500~7691 (no edge mask)
    testAutoCsplMaskFreq : cspline fitting using 500~7691 (no edge mask)
    testAutoSinuMaskFreq : sinusoidal fitting using 500~7691 (no edge mask)
    testAutoPolyChanFlag : polynomial fitting of all channels with channel flag in both edge
    testAutoChebChanFlag : Chebyshev polynomial fitting of all channels with channel flag in both edge
    testAutoCsplChanFlag : cspline fitting of all channels with channel flag in both edge
    testAutoSinuChanFlag : sinusoidal fitting of all channels with channel flag in both edge
    """
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_lftest'
    outfile = outroot+".ms"
    bloutput = outroot+"_blout"
    base_param = dict(infile=infile,
                      datacolumn='float_data',
                      pol='RR',
                      maskmode = 'auto',
                      thresh=5.0,
                      avg_limit=16,
                      minwidth=16,
                      outfile=outfile,
                      blformat='csv',
                      bloutput=bloutput)
    edge = [500,500]
    spw = '2'
    spwchan = '2:500~7691'
    spwfreq = '2:44052975469.940445~44096877113.524124Hz'#44052978522~44096874062Hz'
    # in either tests,
    statrange = [[1000, 7191]]
    polystat = {'rms': 0.20170082215673005, 'min': -0.42453908920288086,
                'max': 2.0263485908508301, 'median': 0.0034337043762207031,
                'stddev': 0.20170082215673005}
    chebstat = {'rms': 0.20170082215673005, 'min': -0.42453908920288086,
                'max': 2.0263485908508301, 'median': 0.0034337043762207031,
                'stddev': 0.20170082215673005}
    csplstat = {'rms': 0.20181625130943376, 'min': -0.42370939254760742,
                'max': 2.0274257659912109, 'median': 0.0038695335388183594,
                'stddev': 0.20181625130943376}
#     sinustat = {'max': , 'min': , 'median': , 'rms': , 'stddev': }

    def setUp(self):
        for prevout in glob.glob(self.outroot+'*'):
            if os.path.isdir(prevout): shutil.rmtree(prevout)
            else: os.remove(prevout)
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)

    def tearDown(self):
        if (os.path.exists(self.infile)):
           shutil.rmtree(self.infile)
        for outname in glob.glob(self.outroot+'*'):
           if os.path.isdir(outname): shutil.rmtree(outname)
           else: os.remove(outname)

    def flag(self, infile, edge=None, rowidx=None):
        rowflag = True if edge is None else False
        if type(rowidx)==int: rowidx = [rowidx]
        tb.open(infile, nomodify=False)
        if rowidx is None: rowidx = list(range(tb.nrows()))
        try:
            for idx in rowidx:
                specs = tb.getcell("FLAG", idx)
                if rowflag: 
                    specs = True
                else:
                    for ipol in range(len(specs)):
                        specs[ipol][0:edge[0]] = True
                        specs[ipol][-edge[1]:] = True
                tb.putcell('FLAG', idx, specs)
        finally:
            tb.close()

    def run_test(self, refstat, **kwargs):
        task_param = self.base_param.copy()
        for key, val in list(kwargs.items()):
            task_param[key] = val
        sdbaseline(**task_param)
        outfile = task_param['outfile']
        polid = 0 if task_param['pol'] in ['RR', 'LL'] else None
        currstat = self._getStats(outfile, spw='0', pol=polid,
                                   colname=task_param['datacolumn'].upper(),
                                   mask=self.statrange)
        self._compareStats(currstat[0],refstat)

        

    def testAutoPolyNoMask(self):
        """polynomial fitting using all channels but edge=[500, 500]"""
        self.run_test(self.polystat, spw=self.spw, edge=self.edge, blfunc='poly')
        
    def testAutoChebNoMask(self):
        """Chebyshev polynomial fitting using all channels but edge=[500, 500]"""
        self.run_test(self.chebstat, spw=self.spw, edge=self.edge, blfunc='chebyshev')

    def testAutoCsplNoMask(self):
        """cspline fitting using all channels but edge=[500, 500]"""
        self.run_test(self.csplstat, spw=self.spw, edge=self.edge, blfunc='cspline')

#     def testAutoSinuNoMask(self):
#         """sinusoidal fitting using all channels but edge=[500, 500]"""
#         self.run_test(self.sinustat, spw=self.spw, edge=self.edge, blfunc='sinusoid')

    def testAutoPolyMaskChan(self):
        """polynomial fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.polystat, spw=self.spwchan, edge=[0,0], blfunc='poly')
        
    def testAutoChebMaskChan(self):
        """Chebyshev polynomial fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.chebstat, spw=self.spwchan, edge=[0,0], blfunc='chebyshev')

    def testAutoCsplMaskChan(self):
        """cspline fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.csplstat, spw=self.spwchan, edge=[0,0], blfunc='cspline')

#     def testAutoSinuMaskChan(self):
#         """sinusoidal fitting using 500~7691 channels (no edge mask)"""
#         self.run_test(self.sinustat, spw=self.spwchan, edge=self.noedge, blfunc='sinusoid')

    def testAutoPolyMaskFreq(self):
        """polynomial fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.polystat, spw=self.spwfreq, edge=[0,0], blfunc='poly')
        
    def testAutoChebMaskFreq(self):
        """Chebyshev polynomial fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.chebstat, spw=self.spwfreq, edge=[0,0], blfunc='chebyshev')

    def testAutoCsplMaskFreq(self):
        """cspline fitting using 500~7691 channels (no edge mask)"""
        self.run_test(self.csplstat, spw=self.spwfreq, edge=[0,0], blfunc='cspline')

#     def testAutoSinuMaskFreq(self):
#         """sinusoidal fitting using 500~7691 channels (no edge mask)"""
#         self.run_test(self.sinustat, spw=self.spwfreq, edge=self.noedge, blfunc='sinusoid')

    def testAutoPolyChanFlag(self):
        """polynomial fitting of all channels with channel flag in both edge"""
        self.flag(self.infile,edge=self.edge)
        self.run_test(self.polystat, spw=self.spw, edge=[0,0], blfunc='poly')
        
    def testAutoChebChanFlag(self):
        """Chebyshev polynomial of all channels with channel flag in both edge"""
        self.flag(self.infile,edge=self.edge)
        self.run_test(self.chebstat, spw=self.spw, edge=[0,0], blfunc='chebyshev')

    def testAutoCsplChanFlag(self):
        """cspline fitting of all channels with channel flag in both edge"""
        self.flag(self.infile,edge=self.edge)
        self.run_test(self.csplstat, spw=self.spw, edge=[0,0], blfunc='cspline')

#     def testAutoSinuChanFlag(self):
#         """sinusoidal fitting of all channels with channel flag in both edge"""
#         self.flag(self.infile,edge=self.edge)
#         self.run_test(self.sinustat, spw=self.spw, edge=self.noedge, blfunc='sinusoid')

class sdbaseline_selection(unittest.TestCase):
    datapath = os.environ.get('CASAPATH').split()[0] + \
              '/data/regression/unittest/tsdbaseline/'
    infile = "analytic_type1.bl.ms"
    outfile = "baselined.ms"
    bloutfile = infile + "_blparam.txt"
    common_param = dict(infile=infile, outfile=outfile,
                        maskmode='list', blmode='fit', dosubtract=True,
                        blfunc='poly', order=1)
    selections=dict(intent=("CALIBRATE_ATMOSPHERE#OFF*", [1]),
                    antenna=("DA99", [1]),
                    field=("M1*", [0]),
                    spw=(">6", [1]),
                    timerange=("2013/4/28/4:13:21",[1]),
                    scan=("0~8", [0]),
                    pol=("YY", [1]))
    # baseline mask for each row of MS
    chan_mask = {'float_data': ("0~19;21~127", "0~39;41~127"),
                 'corrected': ("0~59;61~127", "0~79;81~127")}
    # data of line (chan, amp) for each pol and row of MS
    line_data = {'float_data': {'r0': ((20, 50.0), (20, 100.0)),
                                'r1': ((40, 150.0), (40, 200.0))},
                 'corrected': {'r0': ((60, 75.0), (60, 125.0)),
                               'r1': ((80, 175.0), (80, 225.0))} }
    templist = [infile, outfile, bloutfile]
    verbose = False
 
    def _clearup(self):
        for name in self.templist:
            if os.path.isdir(name):
                shutil.rmtree(name)
            elif os.path.exists(name):
                os.remove(name)

    def setUp(self):
        self._clearup()
        shutil.copytree(self.datapath+self.infile, self.infile)
        default(sdbaseline)

    def tearDown(self):
        self._clearup()

    def _get_selection_string(self, key):
        if key not in list(self.selections.keys()):
            raise ValueError("Invalid selection parameter %s" % key)
        return {key: self.selections[key][0]}

    def _get_selected_row_and_pol(self, key):
        if key not in list(self.selections.keys()):
            raise ValueError("Invalid selection parameter %s" % key)
        pols = [0,1]
        rows = [0,1]
        if key == 'pol':  #self.selection stores pol ids
            pols = self.selections[key][1]
        else: #self.selection stores row ids
            rows = self.selections[key][1]
        return (rows, pols)

    def _get_reference(self, nchan, irow, ipol, datacol):
        line_chan, line_amp = self.line_data[datacol][('r%d' % irow)][ipol]
        reference = numpy.zeros(nchan)
        reference[line_chan] = line_amp
        if self.verbose: print(("reference=%s" % str(reference)))
        return reference

    def _format_spw_mask(self, datacolumn, sel_param):
        (rowids, polids) = self._get_selected_row_and_pol(sel_param)
        spwstr = "*"
        if sel_param=="spw":
            spwstr = self._get_selection_string(sel_param)['spw']
        if len(rowids) == 1:
            return ("%s:%s" % (spwstr, self.chan_mask[datacolumn][rowids[0]]))
        else:
            spwids = ['6', '7']
            spwstr = ""
            for irow in rowids:
                if len(spwstr) > 0: spwstr = spwstr + ","
                spwstr = spwstr + \
                    ("%s:%s" % (spwids[irow], self.chan_mask[datacolumn][irow]))
            return spwstr
    
    def run_test(self, sel_param, datacolumn, reindex=True):
        inparams = self._get_selection_string(sel_param)
        inparams['spw'] = self._format_spw_mask(datacolumn, sel_param)
        inparams.update(self.common_param)
        print(("task param: %s" % str(inparams)))
        sdbaseline(datacolumn=datacolumn, reindex=reindex, **inparams)
        self._test_result(inparams["outfile"], sel_param, datacolumn)
        
    def _test_result(self, msname, sel_param, dcol, atol=1.e-5, rtol=1.e-5, applymode=False):
        # Make sure output MS exists
        self.assertTrue(os.path.exists(msname), "Could not find output MS")
        # Compare output MS with reference (nrow, npol, and spectral values)
        (rowids, polids) = self._get_selected_row_and_pol(sel_param)
        poltest = (sel_param == "pol")
        if dcol.startswith("float"):
            testcolumn = "FLOAT_DATA"
        else: #output is in DATA column
            testcolumn = "DATA"
        tb.open(msname)
        try:
            if not applymode: # normal fit
                self.assertEqual(tb.nrows(), len(rowids), "Row number is wrong %d (expected: %d)" % (tb.nrows(), len(rowids)))
            else: # in case of apply, rownumber does not change from input MS
                self.assertGreaterEqual(tb.nrows(), numpy.array(rowids).max(),
                                        'Reference row number is larger than table size.')
            for out_row in range(len(rowids)):
                in_row = rowids[out_row]
                if applymode: out_row = in_row
                sp = tb.getcell(testcolumn, out_row)
                if not poltest:
                    self.assertEqual(sp.shape[0], len(polids), "Number of pol is wrong in row=%d:  %d (expected: %d)" % (out_row,len(polids),sp.shape[0]))
                nchan = sp.shape[1]
                for out_pol in range(len(polids)):
                    in_pol = polids[out_pol]
                    reference = self._get_reference(nchan, in_row, in_pol, dcol)
                    if self.verbose: print(("data=%s" % str(sp[out_pol])))
                    self.assertTrue(numpy.allclose(sp[out_pol], reference,
                                                   atol=atol, rtol=rtol),
                                    "Baselined spectrum differs in row=%d, pol=%d" % (out_row, out_pol))
        finally:
            tb.close()
        
    def run_test_apply(self, sel_param, datacolumn, reindex=True):
        """BL table generation + application"""
        inparams = self._get_selection_string(sel_param)
        inparams['spw'] = self._format_spw_mask(datacolumn, sel_param)
        inparams.update(self.common_param)
        outms = inparams['outfile']
        bltable = outms+'.bl.cal'
        print('generate BL table')
        inparams.update(dict(dosubtract=False, blformat='table', bloutput=bltable, outfile=''))
        sdbaseline(datacolumn=datacolumn, reindex=reindex, **inparams)
        self.assertTrue(os.path.exists(bltable), 'Failed to generate BL caltable')
        self.assertFalse(os.path.exists(outms), 'Output MS should not be generated yet.')
        print('apply BL table')
        sdbaseline(datacolumn=datacolumn, reindex=reindex, infile=inparams['infile'],
                    outfile=outms, blmode='apply', bltable=bltable)
        self._test_result(outms, sel_param, datacolumn, applymode=True)

    def testIntentF(self):
        """Test selection by intent (float_data)"""
        self.run_test("intent", "float_data")

    def testIntentC(self):
        """Test selection by intent (corrected)"""
        self.run_test("intent", "corrected")

    def testAntennaF(self):
        """Test selection by antenna (float_data)"""
        self.run_test("antenna", "float_data")

    def testAntennaC(self):
        """Test selection by antenna (corrected)"""
        self.run_test("antenna", "corrected")

    def testFieldF(self):
        """Test selection by field (float_data)"""
        self.run_test("field", "float_data")

    def testFieldC(self):
        """Test selection by field (corrected)"""
        self.run_test("field", "corrected")

    def testSpwF(self):
        """Test selection by spw (float_data)"""
        self.run_test("spw", "float_data")

    def testSpwC(self):
        """Test selection by spw (corrected)"""
        self.run_test("spw", "corrected")

    def testTimerangeF(self):
        """Test selection by timerange (float_data)"""
        self.run_test("timerange", "float_data")

    def testTimerangeC(self):
        """Test selection by timerange (corrected)"""
        self.run_test("timerange", "corrected")

    def testScanF(self):
        """Test selection by scan (float_data)"""
        self.run_test("scan", "float_data")

    def testScanC(self):
        """Test selection by scan (corrected)"""
        self.run_test("scan", "corrected")

    def testPolF(self):
        """Test selection by pol (float_data)"""
        self.run_test("pol", "float_data")

    def testPolC(self):
        """Test selection by pol (corrected)"""
        self.run_test("pol", "corrected")

    def testReindexSpw(self):
        """Test reindex =T/F in spw selection"""
        outfile = self.common_param['outfile']
        for datacol in ['float_data', 'corrected']:
            print(("Test: %s" % datacol.upper()))
            for (reindex, ddid, spid) in zip([True, False], [0, 1], [0,7]):
                print(("- reindex=%s" % str(reindex)))
                self.run_test("spw", datacol, reindex=reindex)
                tb.open(outfile)
                try:
                    self.assertEqual(ddid, tb.getcell('DATA_DESC_ID', 0),
                                     "comparison of DATA_DESCRIPTION_ID failed.")
                finally: tb.close()
                tb.open(outfile+'/DATA_DESCRIPTION')
                try:
                    self.assertEqual(spid, tb.getcell('SPECTRAL_WINDOW_ID', ddid),
                                     "comparison of SPW_ID failed.")
                finally: tb.close()
                shutil.rmtree(outfile)
                os.remove('%s_blparam.txt' % self.common_param['infile'])

    def testReindexIntent(self):
        """Test reindex =T/F in intent selection"""
        outfile = self.common_param['outfile']
        for datacol in ['float_data', 'corrected']:
            print(("Test: %s" % datacol.upper()))
            for (reindex, idx) in zip([True, False], [0, 4]):
                print(("- reindex=%s" % str(reindex)))
                self.run_test("intent", datacol, reindex=reindex)
                tb.open(outfile)
                try:
                    self.assertEqual(idx, tb.getcell('STATE_ID', 0),
                                     "comparison of state_id failed.")
                finally: tb.close()
                shutil.rmtree(outfile)
                os.remove('%s_blparam.txt' % self.common_param['infile'])

def suite():
    return [sdbaseline_basicTest, 
            sdbaseline_maskTest,
            sdbaseline_sinusoidTest,
            sdbaseline_outbltableTest,
            sdbaseline_applybltableTest,
            sdbaseline_variableTest,
            sdbaseline_bloutputTest,
            sdbaseline_autoTest,
            sdbaseline_selection
            ]
