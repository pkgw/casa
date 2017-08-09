#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for the cvel task                               #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    The task cvel needs to be exercised and compared to clean              #
#                                                                           #
# Features tested:                                                          #
#    1) does cvel run without raising exceptions for input frame TOPO       #
#       and all possible output frames?                                     #
#    2) can clean process the cvel output?                                  #
#    3) does cvel+clean produce compatible results to clean-only?           #
#       (channel flux values, channel world coordinates)                    #
#                                                                           #
# Input data:                                                               #
#     one dataset, one scan of a VLA observation provided by Crystal Brogan #
#                                                                           #
#############################################################################

# copy module needed to create dictionaries to store the results
import copy

myname = 'cvel_regression'

# default dataset name
pathname=os.environ.get('CASAPATH').split()[0]
dataset_name_orig = pathname + '/data/regression/cvel/input/ANTEN_sort_hann_for_cvel_reg-thinned.ms'

# get the dataset name from the wrapper if possible
mydict = locals()
if "dataset_name" in mydict:
    dataset_name_orig = mydict["dataset_name"]

def isnear(a,b,p):
    print("  ", a, b)
    if(a==b):
        print("  exactly equal")
        return True
    dev = abs(a-b)
    print("  deviation = ", dev)
    if(dev<=p):
        return True
    return False

def isrnear(a,b,p):
    print("  ", a, b)
    if(a==b):
        print("  exactly equal")
        return True
    rdev = abs(a-b)/abs(max(a,b))
    print("  relative deviation = ", rdev)
    if(rdev<=p):
        return True
    return False


# start frequencies of the grids for the different output reference frames

freqmodestart = { 'TOPO': '2.21078e+10Hz',
                  'LSRK': '2.21067e+10Hz',
                  'LSRD': '2.21066e+10Hz',
                  'BARY': '2.21065e+10Hz',
                  'GALACTO': '2.21181e+10Hz',
                  'LGROUP': '2.2125e+10Hz',
                  'CMB': '2.20804e+10Hz'
                  }

# hanning smooth switches for the different output frames

dohanning = { 'TOPO': False,
              'LSRK': True,
              'LSRD': False,
              'BARY': False,
              'GALACTO': False,
              'LGROUP': False,
              'CMB': False
              }

#channel widths

freqmodewidth = { 'TOPO': '1.953125e+05Hz',
                  'LSRK': '1.953027e+05Hz',
                  'LSRD': '1.953017e+05Hz',
                  'BARY': '1.953005e+05Hz',
                  'GALACTO': '1.954028e+05Hz',
                  'LGROUP': '1.954639e+05Hz',
                  'CMB': '1.950702e+05Hz'
                  }

# channel to test with imstats

peakchan = '28'
otherchan1 = '5'
otherchan2 = '25'
otherchan3 = '30'
otherchan4 = '45'

testregion = '135,127,135,127'

# storage for results
imstats = { 'TOPO': 0,'LSRK': 0, 'LSRD': 0, 'BARY': 0, 'GALACTO': 0, 'LGROUP': 0, 'CMB': 0 }
mode_imstats = { peakchan: copy.deepcopy(imstats),
                 otherchan1: copy.deepcopy(imstats),
                 otherchan2: copy.deepcopy(imstats),
                 otherchan3: copy.deepcopy(imstats),
                 otherchan4: copy.deepcopy(imstats) }
cvel_imstats = { 'frequency': copy.deepcopy(mode_imstats),
                 'radio velocity': copy.deepcopy(mode_imstats),
                 'optical velocity': copy.deepcopy(mode_imstats),
                 'channel': copy.deepcopy(mode_imstats) }
cleanonly_imstats = copy.deepcopy(cvel_imstats)


# First we want to test the regridding and spectral frame transformation only.
#  So we combine the spws in the input dataset beforehand.
# Also: clean needs write access to the input MS, so we need a local copy anyway.

dataset_name = dataset_name_orig
os.system('rm -rf input.ms input2.ms')
os.system('cp -RL '+dataset_name_orig+' input.ms')
os.system('chmod -R u+w input.ms')
#os.system('cp -RL '+dataset_name_orig+' input2.ms')
os.system('chmod -R u+w input2.ms')
clean_inputvis_local_copy = 'input.ms'

clean_inputvis_local_copy2 = 'input2.ms' # we need a second copy for the hanning smoothed cases
hanningsmooth2(vis=dataset_name_orig, outputvis=clean_inputvis_local_copy2)


# loop over all possible output reference frames

# these are all possible frames:
#frames_to_do = ['TOPO','LSRK', 'LSRD', 'BARY', 'GALACTO', 'LGROUP', 'CMB']

# the most critical one is CMB (largest freq shift)
#frames_to_do = ['CMB']

# in order to shorten the test, leave out LSRD, GALACTO, and TOPO
frames_to_do = ['LGROUP', 'LSRK', 'BARY', 'CMB']
#frames_to_do = ['LSRK']

for frame in frames_to_do:

    restfrq = 22.235080E9
    restfreqstr = str(restfrq)+'Hz'

    ### frequency mode

    outvis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_freq.ms'
    os.system('rm -rf '+outvis)

    casalog.post(outvis, 'INFO')

    cvel(vis=dataset_name, outputvis=outvis,
         mode='frequency',nchan=50,
         start=freqmodestart[frame],
         width=freqmodewidth[frame],
         interpolation='linear',
         phasecenter='J2000 12h01m53.13s -18d53m09.8s',
         outframe=frame,
         hanning = dohanning[frame])

    invis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_freq.ms'
    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_freq_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    clean(vis=invis,
          imagename=iname,
          field='',spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='frequency',nchan=50,
          start=freqmodestart[frame],
          width=freqmodewidth[frame],
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False)

    cvel_imstats['frequency'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cvel_imstats['frequency'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cvel_imstats['frequency'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cvel_imstats['frequency'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cvel_imstats['frequency'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)

    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_freq_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    cvis = clean_inputvis_local_copy
    if(dohanning[frame]):
        casalog.post('Will Hanning smooth before cleaning ...', 'INFO')
        cvis = clean_inputvis_local_copy2

    clean(vis=cvis,
          imagename=iname,
          field='', spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='frequency',nchan=50,
          start=freqmodestart[frame],
          width=freqmodewidth[frame],
          outframe=frame,
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False)

    cleanonly_imstats['frequency'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cleanonly_imstats['frequency'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cleanonly_imstats['frequency'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cleanonly_imstats['frequency'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cleanonly_imstats['frequency'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)


    #### velocity mode (radio)

    f1 = qa.quantity(freqmodestart[frame])['value']
    f2 = f1+qa.quantity(freqmodewidth[frame])['value']

    vrads = (restfrq-f1)/restfrq *  2.99792E8
    vradstart = str(vrads)+'m/s'
    vradw = (restfrq-f2)/restfrq *  2.99792E8 - vrads
    vradwidth = str(vradw)+'m/s'

    outvis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vrad.ms'
    os.system('rm -rf '+outvis)

    casalog.post(outvis, 'INFO')

    cvel(vis=dataset_name,outputvis=outvis,
         mode='velocity',nchan=49,
         start=vradstart,
         width=vradwidth,
         interpolation='fftshift',
         phasecenter='J2000 12h01m53.13s -18d53m09.8s',
         restfreq=restfreqstr,
         outframe=frame,
         hanning=dohanning[frame])

    invis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vrad.ms'
    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vrad_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    clean(vis=invis,
          imagename=iname,
          field='',spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='velocity',nchan=49,
          start=vradstart,
          width=vradwidth,
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False)

    cvel_imstats['radio velocity'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cvel_imstats['radio velocity'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cvel_imstats['radio velocity'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cvel_imstats['radio velocity'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cvel_imstats['radio velocity'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)

    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_vrad_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    cvis = clean_inputvis_local_copy
    if(dohanning[frame]):
        casalog.post('Will Hanning smooth before cleaning ...', 'INFO')
        cvis = clean_inputvis_local_copy2

    clean(vis=cvis,
          imagename=iname,
          field='', spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='velocity',nchan=49,
          start=vradstart,
          width=vradwidth,
          outframe=frame,
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False)

    cleanonly_imstats['radio velocity'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cleanonly_imstats['radio velocity'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cleanonly_imstats['radio velocity'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cleanonly_imstats['radio velocity'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cleanonly_imstats['radio velocity'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)

    #### velocity mode (optical)

    lambda0 = 2.99792E8/restfrq
    lambda1 = 2.99792E8/f1
    lambda2 = 2.99792E8/f2
    vopts = (lambda1-lambda0)/lambda0 * 2.99792E8
    voptw = (lambda2-lambda0)/lambda0 * 2.99792E8 - vopts
    voptstart = str(vopts)+'m/s'
    voptwidth = str(voptw)+'m/s'

    outvis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vopt.ms'
    os.system('rm -rf '+outvis)

    casalog.post(outvis, 'INFO')

    cvel(vis=dataset_name, outputvis=outvis,
         mode='velocity',nchan=49,
         start=voptstart,
         width=voptwidth,
         interpolation='linear',
         phasecenter='J2000 12h01m53.13s -18d53m09.8s',
         restfreq=restfreqstr,
         outframe=frame,
         veltype='optical',
         hanning=dohanning[frame])

    invis = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vopt.ms'
    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_cvel_vopt_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    clean(vis=invis,
          imagename=iname,
          field='',spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='velocity',nchan=49,
          start=voptstart,
          width=voptwidth,
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False,
          veltype='optical')

    cvel_imstats['optical velocity'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cvel_imstats['optical velocity'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cvel_imstats['optical velocity'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cvel_imstats['optical velocity'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cvel_imstats['optical velocity'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)


    iname = 'ANTEN_ALLRE_sort_hann_'+frame+'_vopt_clean'
    os.system('rm -rf '+iname+'.*')

    casalog.post(iname, 'INFO')

    cvis = clean_inputvis_local_copy
    if(dohanning[frame]):
        casalog.post('Will Hanning smooth before cleaning ...', 'INFO')
        cvis = clean_inputvis_local_copy2

    clean(vis=cvis,
          imagename=iname,
          field='', spw='',
          cell=[0.01,0.01],imsize=[256,256],
          stokes='I',
          mode='velocity',nchan=49,
          start=voptstart,
          width=voptwidth,
          outframe=frame,
          interpolation='linear',
          psfmode='clark',imagermode='csclean',
          scaletype='SAULT',
          niter=0,threshold='1.5mJy',
          restfreq=restfreqstr,
          phasecenter='J2000 12h01m53.13s -18d53m09.8s',
          mask='',
          weighting='briggs',
          interactive=False,
          minpb=0.3,pbcor=False,
          veltype='optical')

    cleanonly_imstats['optical velocity'][peakchan][frame] = imstat(iname+'.image', box=testregion, chans=peakchan)
    cleanonly_imstats['optical velocity'][otherchan1][frame] = imstat(iname+'.image', box=testregion, chans=otherchan1)
    cleanonly_imstats['optical velocity'][otherchan2][frame] = imstat(iname+'.image', box=testregion, chans=otherchan2)
    cleanonly_imstats['optical velocity'][otherchan3][frame] = imstat(iname+'.image', box=testregion, chans=otherchan3)
    cleanonly_imstats['optical velocity'][otherchan4][frame] = imstat(iname+'.image', box=testregion, chans=otherchan4)

# end loop over frames

# Analysis

passed = True
# normal values
#tolerance = 0.0013
#avtolerance = 0.00035
# temporarily relaxed values, should be reduced to the normal ones after March 2011
tolerance = 0.012
avtolerance = 0.0039
numpoints = 0.
avdev = 0.
maxdev = 0.
maxdevat = " "
problems = 0
for frame in frames_to_do:
    for mode in ['frequency', 'radio velocity', 'optical velocity']:
        for chan in list(mode_imstats.keys()):
            print("Testing ", frame, ", ",  mode, ", Hanning ", dohanning[frame], ", box ", testregion, ", channel ", chan, " ...")

            isok = True
            c1 = 0.
            c2 = 0.
            s1 = 0.
            s2 = 0.
            try:
                c1 = cleanonly_imstats[mode][chan][frame]['max']
                c2 = cvel_imstats[mode][chan][frame]['max']
                s1 = cleanonly_imstats[mode][chan][frame]['maxposf']
                s2 = cvel_imstats[mode][chan][frame]['maxposf']
            except:
                isok = False
                passed = False
                problems +=1
                print("Error: this subtest failed")

            if isok:
                if(abs(c1-c2) > maxdev):
                    maxdev = abs(c1-c2)
                    maxdevat = mode+" mode for output frame "+frame\
                               +":\n    cvel+clean finds max flux in channel "+str(chan)+" to be "+str(c2)\
                               +"\n    clean-only finds max flux in channel "+str(chan)+" to be "+str(c1)
                if(not isnear(c1,c2, tolerance)):
                    print(" ** Problem in ", mode, " mode for output frame ", frame, ":")
                    print("     cvel+clean finds max flux in channel ", chan, " to be ", c2)
                    print("     clean-only finds max flux in channel ", chan, " to be ", c1)
                    passed = False
                    isok = False
                    problems +=1

                avdev += abs(c1-c2)
                numpoints += 1.

                if(not s1 == s2):
                    print(" ** Problem in ", mode, " mode for output frame ", frame, ":")
                    print("     cvel+clean finds world coordinates for channel ", chan, " to be ", s2)
                    print("     clean-only finds world coordinates for channel ", chan, " to be ", s1)
                    passed = False
                    isok = False
                    problems +=1
                else:
                    print("  World coordinates identical == ", s2)

                if isok:
                    print("... OK")

if(numpoints > 0.):
    avdev = avdev/numpoints
    print(numpoints, " spectral points compared, average deviation = ", avdev, " Jy")
    if(avdev>avtolerance):
        passed = False
        print(" ** Problem: average deviation too large. Expected is value < ", avtolerance)
        problems += 1
    print("   maximum deviation = ", maxdev, " in ", maxdevat)

if passed:
    print('')
    print('Regression PASSED')
    print('')
else:
    print("Execution successful but found ", problems, " issues in analysis of results.")
    print('')
    print('Regression FAILED')
    print('')
    raise
