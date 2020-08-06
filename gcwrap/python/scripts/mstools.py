
def write_history(myms, vis, tname, param_names, param_vals, myclog=None, debug=False):
        """
        Update vis with the parameters that task tname was called with.

        myms - an ms tool instance
        vis  - the MS to write to.
        tname - name of the calling task.
        param_names - list of parameter names.
        param_vals - list of parameter values (in the same order as param_names).
        myclog - a casalog instance (optional)
        debug - Turns on debugging print statements on errors if True.

        Example:
        The end of split does
        param_names = split.func_code.co_varnames[:split.func_code.co_argcount]
        param_vals = [eval(p) for p in param_names]  # Must be done in the task.
        write_history(myms, outputvis, 'split', param_names, param_vals,
                      casalog),
        which appends, e.g.,

        vis = 'TWHydra_CO3_2.ms'
        outputvis   = 'scan9.ms'
        datacolumn  = 'data'
        field       = ''
        spw         = ''
        width       = 1
        antenna     = ''
        timebin     = '0s'
        timerange   = ''
        scan        = '9'
        intent      = ''
        array       = ''
        uvrange     = ''
        correlation = ''
        keepflags   = True
        async       = False

        to the HISTORY of outputvis.
        """
        if not hasattr(myms, 'writehistory'):
                if debug:
                        print("write_history(myms, %s, %s): myms is not an ms tool" % (vis, tname))
                return False
        retval = True
        isopen = False
        try:
                if not myclog and hasattr(casalog, 'post'):
                        myclog = casalog
        except Exception as instance:
                # There's no logger to complain to, and I don't want to exit
                # just because of that.
                pass
        try:
                myms.open(vis, nomodify=False)
                isopen = True
                messages = ['taskname={0}'.format(tname)]
                vestr = 'version: '
                try:
                    # Don't use myclog.version(); it also prints to the
                    # logger, which is confusing.
                    vestr += casa['build']['version'] + ' '
                    vestr += casa['source']['url'].split('/')[-2]
                    vestr += ' rev. ' + casa['source']['revision']
                    vestr += ' ' + casa['build']['time']
                except Exception as instance:
                    if hasattr(myclog, 'version'):
                        # Now give it a try.
                        vestr += myclog.version()
                    else:
                        vestr += ' could not be determined' # We tried.

                messages.append(vestr)

                # Add the task arguments.
                for argnum in range(len(param_names)):
                    msg = "%-11s = " % param_names[argnum]
                    val = param_vals[argnum]
                    if type(val) == str:
                        msg += '"'
                    msg += str(val)
                    if type(val) == str:
                        msg += '"'
                    messages.append(msg)

                myms.writehistory_batch(messages=messages, origin=tname)

        except Exception as instance:
                if hasattr(myclog, 'post'):
                        myclog.post("*** Error \"%s\" updating HISTORY of %s" % (instance, vis),
                                    'SEVERE')
                retval = False
        finally:
                if isopen:
                        myms.close()
        return retval

