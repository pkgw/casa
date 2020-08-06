from taskinit import find_casa
from init_tools import iatool
import os
from stat import S_ISDIR, ST_MTIME, ST_MODE

def write_image_history(myia, tname, param_names, param_vals, myclog=None):
    """
    Update image attached to image tool with the parameters that task tname was called with.

    myia - attached image tool or image name (string)
    tname - name of the calling task.
    param_names - list of parameter names.
    param_vals - list of parameter values (in the same order as param_names).
    myclog - a casalog instance (optional)
    """
    
    myia_is_string = type(myia) == str
    if myia_is_string:
        if not myia:
            # empty string
            return
        _ia = iatool()
        _ia.open(myia)
    elif not hasattr(myia, 'sethistory'):
        return False
    else:
        _ia = myia
    try:
        if not myclog and hasattr(casalog, 'post'):
            myclog = casalog
    except Exception as instance:
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
        except Exception as instance:
            if hasattr(myclog, 'version'):
                # Now give it a try.
                vestr += myclog.version()
            else:
                vestr += ' could not be determined' # We tried.

        _ia.sethistory(tname, vestr)
        # Write the arguments.
        s = tname + "("
        n = len(param_names)
        for argnum in range(n):
            s += param_names[argnum] + "="
            val = param_vals[argnum]
            sval = str(val)
            if len(sval) > 300:
                s += "..."
            else:
                if type(val) == str:
                    s += '"'
                s += sval
                if type(val) == str:
                    s += '"'
            if argnum < n-1:
                s += ", "
        s += ")" 
        _ia.sethistory(tname, s)
    except Exception as instance:
        if hasattr(myclog, 'post'):
            myclog.post("*** Error \"%s\" updating HISTORY of " % (instance),
                        'SEVERE')
        return False
    finally:
        if myia_is_string:
            _ia.done()
    return True

def get_created_images(outfile, target_time):
    dirpath = os.path.dirname(outfile)
    if not dirpath:
        dirpath = "."
    base = os.path.basename(outfile)
    # get all entries in the directory w/ stats
    entries = []
    for fn in os.listdir(dirpath):
        if os.path.basename(fn).startswith(base):
            entries.append((os.stat(fn), fn))
    # leave only directories, insert creation date
    entries = ((stat.st_mtime, path)
        for stat, path in entries if S_ISDIR(stat[ST_MODE]))
    # reverse sort by time
    zz = sorted(entries)
    zz.reverse()
    created_images = []
    for mdate, path in zz:
        # kludge because all of a sudden, some mdates are less than time.time() value
        # that was gotten before these files were created on OSX. Weird.
        if mdate < target_time - 1:
            break
        if os.path.basename(path).startswith(base):
            created_images.append(path)
    return created_images

