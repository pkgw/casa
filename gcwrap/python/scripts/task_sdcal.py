import sys
import os
import numpy
import numpy.random as random
import shutil

from taskinit import *
from applycal import applycal
from mstools import write_history
import types
import sdutil

# Calibrator tool
(cb,myms) = gentools(['cb','ms'])

def sdcal(infile=None, calmode='tsys', fraction='10%', noff=-1,
           width=0.5, elongated=False, applytable='',interp='', spwmap={},
           outfile='', overwrite=False, field='', spw='', scan='',intent=''): 
       
    """ Externally specify calibration solutions of various types
    """
    # Composite mode: compute calibration table and calibrate
    if ',' in calmode:
        handle_composite_mode(locals())
        return
            
    # Single mode: either calibrate or compute calibration table
    try:
        # Parameters check
        if calmode == 'tsys':
            if scan != '':
                raise UserWarning("Scan input must be ''(=all) in calmode='tsys'.")
            if spw != '':
                raise UserWarning("Spw input must be ''(=all) in calmode='tsys'.")

        if isinstance(infile,str) and os.path.exists(infile):
            # check if CORRECTED_DATA is necessary
            addcorr = calmode == 'apply'
            cb.setvi(old=True)
            cb.open(filename=infile,compress=False,addcorr=addcorr,addmodel=False)
            cb.selectvis(spw=spw, scan=scan, field=field)
        else:
            raise Exception('Infile data set not found - please verify the name')

        if not isinstance(calmode,str):
            raise Exception("Calmode must be a string")
        
        if calmode.lower() not in ['tsys', 'ps', 'otfraster', 'otf', 'apply']: 
            raise Exception("Calmode must be either 'ps' or 'otfraster' or  'otf' or 'tsys' or 'apply'.")

        if (not overwrite) and os.path.exists(outfile):
            raise RuntimeError("overwrite is False and output file exists: {}".format(outfile))

        # Calibration
        if calmode == 'apply': # Calibrate using existing tables
            # single calibration table
            if isinstance(applytable, str):
                _table_list = [applytable]

            # multiple calibration tables
            if isinstance(applytable, list) or isinstance(applytable, numpy.ndarray):
                _table_list = applytable
                
            # no calibration table
            if len(_table_list) == 0:
                raise Exception('Applytable name must be specified.')
            
            # check calibration table files
            for table in _table_list:
                # empty string
                if len(table) == 0:
                    raise Exception('Applytable name must be specified.')
                # unexisting table
                if not os.path.exists(table):
                    raise Exception("Table doesn't exist: {}".format(table))
            
            # warning on usage difference with asap.sdcal2
            if (outfile != ''):
                warning_msg = '\n'.join([
                    'The outfile you specified will NOT be created.',
                    "Calibrated data will be stored in a new 'CORRECTED_DATA' column",
                    'inserted in the main table of the input MS file.'
                ])
                casalog.post(warning_msg,priority="WARN")

            if(type(spwmap)!=list and (type(spwmap)!=dict)):
                raise Exception('Spwmap type must be list or dictionary.')

            if (type(spwmap)==dict):
                MS = infile
                tb.open(MS+'/SPECTRAL_WINDOW')
                total_spwID=tb.nrows()
                tb.close()
                
                spwmap_dict = spwmap
                spwmap_list = list(range(total_spwID))

                for key, value in list(spwmap_dict.items()):
                    for v in value:
                        if v in spwmap_list:
                            index = spwmap_list.index(v)
                            spwmap_list[index]=int(key)

                spwmap = spwmap_list
                
            # Setup calibrator
            for _table in _table_list:
                caltype = inspect_caltype(_table)
                if caltype == 'UNKNOWN':
                    raise RuntimeError('Applytable \'%s\' is not a caltable format'%(_table))
                elif caltype == 'B TSYS':
                    cb.setapply(table=_table, spwmap=spwmap, interp=interp, calwt=True)
                else:
                    # no spw mapping is needed for sky calibration
                    cb.setapply(table=_table, interp=interp, calwt=True)
                    
            # Calibrate
            cb.correct(applymode='calflag')
            
            # Write to HISTORY table of MS
            param_names = sdcal.__code__.co_varnames[:sdcal.__code__.co_argcount] 
            param_vals = [eval(p) for p in param_names]
            write_history(myms, infile, 'sdcal', param_names, 
                              param_vals, casalog) 
            

        else: # Compute calibration table
            # Reconciliating 'Python world' calmode with 'C++ world' calmode
            # cpp_calmode[python_calmode] = c++_calmode
            cpp_calmode = { 'tsys': 'tsys',
                            'ps': 'sdsky_ps',
                            'otfraster': 'sdsky_raster',
                            'otf': 'sdsky_otf'
                          }
            
            if len(outfile) == 0:
                raise RuntimeError('Output file name must be specified.')
            
            if calmode == 'tsys':
                cb.specifycal(caltable=outfile,time="",spw=spw,caltype=cpp_calmode[calmode])
            else:
                fraction_numeric = to_numeric_fraction(fraction)
                if noff <= 0 and fraction_numeric >= 0.5:
                    raise ValueError('Too many edge points. fraction must be < 0.5.')
                # Setup calibrator
                cb.selectvis(spw=spw, scan=scan, field=field, intent=intent)
                cb.setsolve(type=cpp_calmode[calmode], table=outfile, fraction=fraction_numeric, numedge=noff)
                # Compute calibration table
                cb.solve()

    except UserWarning as instance:
        print('*** Warning ***',instance)

    except Exception as instance:
        print('*** Error ***',instance)
        raise Exception(instance)

    finally:
        cb.close()

def inspect_caltype(table):
    caltype = 'UNKNOWN'
    with sdutil.tbmanager(table) as tb:
        if 'VisCal' in tb.keywordnames():
            caltype = tb.getkeyword('VisCal')
    return caltype

def to_numeric_fraction(fraction):
    try:
        if isinstance(fraction, str):
            if len(fraction.strip()) == 0:
                # invalid, use default
                fraction_numeric = 0.1
            else:
                pos = fraction.strip().find('%')
                if pos != -1:
                    # percentage 
                    fraction_numeric = float(fraction[:pos]) * 0.01
                else:
                    # direct value
                    fraction_numeric = float(fraction)
        else:
            fraction_numeric = float(fraction)
    except Exception as e:
        casalog.post(str(e), priority='SEVERE', origin='sdcal')
        raise RuntimeError('Invalid fraction value (original error message: "%s")'%(str(e)))

    return fraction_numeric
    
def temporary_name(calmode):
    num_trial = 100
    for i in range(num_trial):
        number = random.random_integers(num_trial)
        name = ('__sdcal_composite_mode_%s_%3s.tab'%(calmode,number)).replace(' ','0')
        if not os.path.exists(name):
            return name
    raise RuntimeError('Failed to configure temporary caltable name.')

def temporary_calibration(calmode, arg_template, **kwargs):
    caltable = temporary_name(calmode)
    myargs = arg_template.copy()
    myargs['calmode'] = calmode
    myargs['outfile'] = caltable
    # try to keep the existing file although
    # outfile should never point to existing file
    myargs['overwrite'] = False
    # optional argument for sdcal
    for (k,v) in list(kwargs.items()):
        if k in myargs:
            myargs[k] = v
    sdcal(**myargs)
    if not os.path.exists(caltable):
        raise RuntimeError('Failed to create temporary caltable.')
    return caltable

def fix_for_intent(calmodes, input_args):
    if 'tsys' in calmodes and ('otfraster' in calmodes or 'otf' in calmodes):
        casalog.post("Intent selection for 'otf' or 'otfraster' should be 'OBSERVE_TARGET#ON_SOURCE'. \n" 
                     "However, the task is not allowed to set global intent selection since calmode contains 'tsys'. \n" 
                     "As a workaround, set intent selection locally when 'otf' or 'otfraster' calibration is performed.",
                     priority='WARN')
        output_args = input_args.copy()
        output_args['intent'] = 'OBSERVE_TARGET#ON_SOURCE'
    else:
        output_args = input_args
    return output_args
        

def handle_composite_mode(args):
    kwargs = args.copy()
    calmodes = kwargs['calmode'].split(',')
    _applytable = kwargs['applytable']
    if isinstance(_applytable, str):
        if len(_applytable) > 0:
            precalibrations = [_applytable]
        else:
            precalibrations = []
    else:
        precalibrations = list(_applytable[:])
    
    applytable_list = []
    try:
        # sky calibration
        if 'ps' in calmodes:
            # ps calibration
            applytable_list.append(
                temporary_calibration('ps', kwargs, spwmap={})
                )
        elif 'otfraster' in calmodes:
            # otfraster calibration
            kwargs_local = fix_for_intent(calmodes, kwargs)
            applytable_list.append(
                temporary_calibration('otfraster', kwargs_local, spwmap={})
                )
        elif 'otf' in calmodes:
            # otf calibration
            kwargs_local = fix_for_intent(calmodes, kwargs)
            applytable_list.append(
                temporary_calibration('otf', kwargs_local, spwmap={})
                )

        # Tsys calibration
        if 'tsys' in calmodes:
            applytable_list.append(
                temporary_calibration('tsys', kwargs, field='', spw='', scan='')
                )

        # apply temporary caltables
        if 'apply' in calmodes:
            if len(applytable_list) == 0:
                raise RuntimeError("No applytable has been created/registered.")
            myargs = kwargs.copy()
            myargs['calmode'] = 'apply'
            myargs['applytable'] = precalibrations + applytable_list
            sdcal(**myargs)

    finally:
        # clean up temporary tables
        for _table in applytable_list:
            if os.path.exists(_table):
                casalog.post('removing \'%s\''%(_table), priority='DEBUG')
                shutil.rmtree(_table)
