from taskinit import find_casa

def write_image_history(myia, tname, param_names, param_vals, myclog=None):
    """
    Update image attached to image tool with the parameters that task tname was called with.

    myia - attached image tool
    tname - name of the calling task.
    param_names - list of parameter names.
    param_vals - list of parameter values (in the same order as param_names).
    myclog - a casalog instance (optional)
    """
    if not hasattr(myia, 'sethistory'):
        return False
    try:
        if not myclog and hasattr(casalog, 'post'):
            myclog = casalog
    except Exception, instance:
        # There's no logger to complain to, and I don't want to exit
        # just because of that.
        pass
    try:
        vestr = 'version: '
        try:
            casa = find_casa()
            # Don't use myclog.version(); it also prints to the
            # logger, which is confusing.
            vestr += casa['build']['version'] + ' '
            #vestr += casa['source']['url']
            #vestr += ' rev. ' + casa['source']['revision']
            vestr += ' ' + casa['build']['time']
        except Exception, instance:
            if hasattr(myclog, 'version'):
                # Now give it a try.
                vestr += myclog.version()
            else:
                vestr += ' could not be determined' # We tried.

        myia.sethistory(tname, vestr)
        # Write the arguments.
        s = tname + "("
        n = len(param_names)
        for argnum in xrange(n):
            s += param_names[argnum] + "="
            val = param_vals[argnum]
            if type(val) == str:
                s += '"'
            s += str(val)
            if type(val) == str:
                s += '"'
            if argnum < n-1:
                s += ", "
        s += ")" 
        myia.sethistory(tname, s)
    except Exception, instance:
        if hasattr(myclog, 'post'):
            myclog.post("*** Error \"%s\" updating HISTORY of " % (instance),
                        'SEVERE')
        return False
    return True

