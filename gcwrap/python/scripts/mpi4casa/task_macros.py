#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
from taskinit import * # casalog, casautils, etc
import sys, inspect, string # To navigate and access stack frames
import os # To handle files and environmental variables
import odict # To handle parameters dictionary

def update_params(func, printtext=True):

    a = inspect.stack()
    stacklevel = 0
    for k in range(len(a)):
        if a[k][1] == "<string>" or (string.find(a[k][1], 'ipython console') > 0 or string.find(a[k][1],"casapy.py") > 0):
            stacklevel = k
    myf = sys._getframe(stacklevel).f_globals

    # Set task to the one being called
    myf['taskname'] = func
    obj = myf[func]

    if (str(type(obj)) == "<type 'instance'>" and
         hasattr(obj, "check_params")):
        hascheck = True
    else:
        hascheck = False

    noerror = True

    # Check if task has defined a task_check_params function
    if (hascheck):
        has_othertasks = 'task_location' in myf
        if(has_othertasks) :
           has_task = myf['taskname'] in myf['task_location']
           if (has_task) :
                pathname = myf['task_location'][myf['taskname']]
           else :
                pathname = xmlpath( )
        else :
           pathname = xmlpath( )
        xmlfile = pathname + '/' + myf['taskname'] + '.xml'
        if(os.path.exists(xmlfile)) :
            cu.setconstraints('file://' + xmlfile);

    a = myf[myf['taskname']].defaults("paramkeys")

    params = a
    itsparams = {}
    for k in range(len(params)):
        paramval = obj.defaults(params[k])

        notdict = True
        ###if a dictionary with key 0, 1 etc then need to peel-open
        ###parameters
        if(type(paramval) == dict):
            if(0 in paramval):
                notdict = False

        if (notdict):
            if(params[k] not in myf):
                myf.update({params[k]:paramval})
            if(printtext):
                if(hascheck):
                    noerror = obj.check_params(params[k], myf[params[k]])
                if(myf[params[k]] == paramval):
                    print_params_col(params[k], myf[params[k]], obj.description(params[k]), 'ndpdef', 'black', noerror)
                else:
                    print_params_col(params[k], myf[params[k]], obj.description(params[k]), 'ndpnondef', 'black', noerror)
                itsparams[params[k]] = myf[params[k]]
        else:
            subdict = odict.odict(paramval)
            ##printtext is False....called most probably to set
            ##undefined params..no harm in doing it anyways
            if(not printtext):
                ##locate which dictionary is user selected
                userdict = {}
                subkeyupdated = {}
                for somekey in paramval:
                    somedict = dict(paramval[somekey])
                    subkeyupdated.update(dict.fromkeys(somedict, False))
                    if('value' in somedict and params[k] in myf):
                        if(somedict['value'] == myf[params[k]]):
                            userdict = somedict
                    elif('notvalue' in somedict and params[k] in myf):
                        if(somedict['notvalue'] != myf[params[k]]):
                            userdict = somedict
                ###The behaviour is to set to the first default
                ### all non set parameters and parameters that
                ### have no meaning for this selection
                for j in range(len(subdict)):
                    subkey = list(subdict[j].keys())
                    for kk in range(len(subkey)):

                        if((subkey[kk] != 'value') & (subkey[kk] != 'notvalue')):
                            #if user selecteddict
                            #does not have the key
                            ##put default
                            if((subkey[kk] not in userdict) and (not subkeyupdated[subkey[kk]])):
                                myf.update({subkey[kk]:subdict[j][subkey[kk]]})
                                subkeyupdated[subkey[kk]] = True

                    ###put default if not there
                            if(subkey[kk] not in myf):
                                myf.update({subkey[kk]:subdict[j][subkey[kk]]})

            ### need to do default when user has not set val
            if(params[k] not in myf):
                if('notvalue' in paramval[0]):
                    myf.update({params[k]:paramval[0]['notvalue']})
                else:
                    myf.update({params[k]:paramval[0]['value']})
            userval = myf[params[k]]
            choice = 0
            notchoice = -1
            valuekey = 'value'
            for j in range(len(subdict)):
                if('notvalue' in subdict[j]):
                    valuekey = 'notvalue'
                    if(subdict[j]['notvalue'] != userval):
                        notchoice = j;
                        break
                else:
                    if(subdict[j]['value'] == userval):
                        choice = j
                        notchoice = j
                        break
            subkey = list(subdict[choice].keys())
            if(hascheck):
                noerror = obj.check_params(params[k], userval)
            if(printtext):
                if(myf[params[k]] == paramval[0][valuekey]):
                    print_params_col(params[k], myf[params[k]], obj.description(params[k]), 'dpdef', 'black', noerror)
                else:
                    print_params_col(params[k], myf[params[k]], obj.description(params[k]), 'dpnondef', 'black', noerror)
                itsparams[params[k]] = myf[params[k]]
            for j in range(len(subkey)):
                if((subkey[j] != valuekey) & (notchoice > -1)):
                    ###put default if not there
                    if(subkey[j] not in myf):
                        myf.update({subkey[j]:subdict[choice][subkey[j]]})
                    paramval = subdict[choice][subkey[j]]
                    if (j == (len(subkey) - 1)):
                        # last subparameter - need to add an extra line to allow cut/pasting
                        comment = 'last'
                    else:
                        comment = 'blue'
                    if(hascheck):
                        noerror = obj.check_params(subkey[j], myf[subkey[j]])
                    if(printtext):
                        if(myf[subkey[j]] == paramval):
                            print_params_col(subkey[j], myf[subkey[j]], obj.description(subkey[j], userval), 'spdef', comment, noerror)
                        else:
                            print_params_col(subkey[j], myf[subkey[j]], obj.description(subkey[j], userval), 'spnondef', comment, noerror)
                        itsparams[params[k]] = myf[params[k]]

# EOF
