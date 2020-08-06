from casac import casac
def print_test(output, testname):
    outlist = output.split("\n")
    if len(outlist) > 1:
        for outline in outlist:
            print_test(outline, testname)
    elif len(output) == 0:
        # just a white line
        print("")
    elif len(testname) == 0:
        print(output)
    else:
        print(("%s: %s" % (testname, output)))

def get_band_center_width(numChan, refChan, refFreq, chanSep):
    """
    Returns center frequency of a band.
    numChan : the number of channels in a band
    refChan : the reference channel ID (0-based)
    refFreq : the frequency at refChan, quantum value
    chanSep : the channel width of a band
    """
    myqa = casac.quanta()
    centerChan = (numChan - 1)/2.
    centerFreq = myqa.add(myqa.mul(centerChan-refChan, chanSep), refFreq)
    bandwidth = myqa.mul(numChan, chanSep)
    return centerFreq, bandwidth
    