#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for ASDM version 1 import to MS                 #
#    and the "inverse filler" task exportasdm, and subsequent               #
#    data analysis.                                                         #
#                                                                           # 
# Rationale for Inclusion:                                                  #
#    The conversion of ASDM to MS and back needs to be verified.            #
#                                                                           # 
# Features tested:                                                          #
#    1) Is the import performed without raising exceptions?                 #
#    2) Do all expected tables exist?                                       #
#    3) Can the MS be opened?                                               #
#    4) Do the tables contain expected values?                              #
#    5) Is exportasdm performed without raising exceptions?                 #
#    6) Is the created ASDM well-formed (XML) and complete?                 #
#    7) Can the resulting ASDM be reimported without raising exceptions?    #
#    8) Does it have the same number of integrations as the original?       #
#    9) Does an imported ASDM pass a serious analysis?                      #
#   10) After exporting the MS generated in (9) to an ASDM, does the        #
#       re-imported MS pass the same analysis with the same results?        #
#                                                                           #
# Input data:                                                               #
#     two datasets for the filler of ASDM 1                                 #
#     one WVR correction cal table                                          #
#     one simulated MS dataset                                              #
#                                                                           #
#############################################################################

myname = 'alma_asdm+reimported-asdm_sf'

# default ASDM dataset name
myasdm_dataset_name = 'uid___X5f_X18951_X1'
myms_dataset_name = 'M51.ms'
myasdm_dataset2_name = 'uid___X02_X56142_X1'
mywvr_correction_file = 'N3256_B6_0.WVR'

# get the dataset name from the wrapper if possible
mydict = locals()
if "asdm_dataset_name" in mydict:
    myasdm_dataset_name = mydict["asdm_dataset_name"]
if "ms_dataset_name" in mydict:
    myms_dataset_name = mydict["ms_dataset_name"]
if "asdm_dataset2_name" in mydict:
    myasdm_dataset_name = mydict["asdm_dataset2_name"]
if "wvr_dataset_name" in mydict:
    mywvr_correction_file = mydict["wvr_correction_file"]

# name of the resulting MS
msname = myasdm_dataset_name+'.ms'

# name of the exported ASDM
asdmname = myms_dataset_name+'.asdm'

# name of the reimported MS
reimp_msname = 'reimported-'+myms_dataset_name

def checktable(thename, theexpectation):
    global msname, myname
    tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print(myname, ": comparing ", mycell)
        value = tb.getcell(mycell[0], mycell[1])
        # see if value is array
        try:
            isarray = value.__len__
        except:
            # it's not an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement = (value == mycell[2])
            else:
                in_agreement = ( abs(value - mycell[2]) < mycell[3]) 
        else:
            # it's an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement =  (value == mycell[2]).all() 
            else:
                in_agreement = (abs(value - mycell[2]) < mycell[3]).all() 
        if not in_agreement:
            print(myname, ":  Error in MS subtable", thename, ":")
            print("     column ", mycell[0], " row ", mycell[1], " contains ", value)
            print("     expected value is ", mycell[2])
            tb.close()
            raise
    tb.close()
    print(myname, ": table ", thename, " as expected.")
    return

#########################

def verify_asdm(asdmname, withPointing):
    print("Verifying asdm ", asdmname)
    if(not os.path.exists(asdmname)):
        print("asdm ", asdmname, " doesn't exist.")
        raise Exception
    # test for the existence of all obligatory tables
    allTables = [ "Antenna.xml",
                  "ASDM.xml",
                 # "CalData.xml",
                 # "CalDelay.xml",
                 # "CalReduction.xml",
                  "ConfigDescription.xml",
                  "CorrelatorMode.xml",
                  "DataDescription.xml",
                  "ExecBlock.xml",
                  "Feed.xml",
                  "Field.xml",
                 #"FocusModel.xml",
                 #"Focus.xml",
                  "Main.xml",
                  "PointingModel.xml",
                  "Polarization.xml",
                  "Processor.xml",
                  "Receiver.xml",
                  "SBSummary.xml",
                  "Scan.xml",
                  "Source.xml",
                  "SpectralWindow.xml",
                  "State.xml",
                  "Station.xml",
                  "Subscan.xml",
                  "SwitchCycle.xml"
                  ]
    isOK = True
    # test if xmllint is available and can be used in the following
    xmllint_ok = (os.system('xmllint --version') == 0)
    for fileName in allTables:
        filePath = asdmname+'/'+fileName
        if(not os.path.exists(filePath)):
            print("ASDM table file ", filePath, " doesn't exist.")
            isOK = False
        elif(xmllint_ok):
            # test if well formed
            rval = os.system('xmllint --noout '+filePath)
            if(rval !=0):
                print("Table ", filePath, " is not a well formed XML document.")
                isOK = False
    if(isOK and not xmllint_ok):
        print("Note: Test of XML well-formedness not possible since xmllint not available.")
    else:
        print("Note: xml validation not possible since ASDM DTDs (schemas) not yet online.")
        
    if(not os.path.exists(asdmname+"/ASDMBinary")):
        print("ASDM binary directory "+asdmname+"/ASDMBinary doesn't exist.")
        isOK = False

    if(withPointing and not os.path.exists(asdmname+"/Pointing.bin")):
        print("ASDM binary file "+asdmname+"/Pointing.bin doesn't exist.")
        isOK = False

    if (not isOK):
        raise Exception


def analyseASDM(basename, caltablename0, genwvr=True):
    # Reduction of NGC3256 Band 6
    # M. Zwaan, May 2010
    # D. Petry, May 2010
    
    # We ignore flux calibration for now.
    # The script does bandpass, gain calibration (gaincal), WVR correction and a delay corrections.
    # Calibration tables are applied with applycal and images of the calibrator and the
    # galaxy are made.

    isOK = True
    
    msname=basename+'.ms'
    contimage=basename+'_cont'
    lineimage=basename+'_line'
    bname=basename+'_B'
    caltablename=bname+'_cal'
    msn=bname+'.ms'

    calfield="J1037*" # which field is the phase calibrator
    sciencefield="ngc3256" # which field is the target

    avspw=["0:17~118","1:17~118"] # which channels to average for the gain calibration
    delay=[0,-8] # there is an 8ns delay problem in the second BB

        
    # Find the asdm
    asdm=basename
    os.system('rm -rf '+bname+'_* '+msn)
    print(">> Importing the asdm: ", asdm, " as measurement set: ",msn) 
    importasdm(
        asdm=asdm,
        vis=msn,
        corr_mode="all",
        srt="all",
        time_sampling="all",ocorr_mode="ca",compression=False,asis="",
        wvr_corrected_data="no",
        verbose=True,
        showversion=False,
        useversion='v3'
        )

    if(cu.compare_version('>=',[3,4,0]) and genwvr): # generate the WVR correction table
        print(">> Generating the WVR caltable ", caltablename0)
        os.system('rm -rf '+caltablename0)
        wvrgcal(vis=msn, caltable=caltablename0, segsource=False, toffset=-2, reversespw='0~7')
    else:
        print(">> Reusing existing WVR caltable ", caltablename0)
    
    # Delay correction
    
    # Could have done this for both spw simultanesouly
    print("\n>> 8ns delay corrections: Calculate K tables")
    for i in range(2):
        os.system('rm -rf '+caltablename+"_spw"+str(i)+'.K')
        gencal(
            vis=msn,
            caltable=caltablename+"_spw"+str(i)+".K",
            caltype="sbd",
            spw=str(i),
            antenna="0",
            pol="",
            parameter=delay[i]
            )

    # Bandpass and Gain calibration
    
    # Bandpass calibration with 'bandpass' 
	
    print(">> Find B solutions")
    for i in range(2):
        print(">> SPW: ",i)
        os.system('rm -rf '+caltablename+'_spw'+str(i)+'.B')
        
        # Use the bright phase calibrator for the bandpass
        bandpass(
            vis=msn,
            spw=avspw[i],  #str(i),
            field=calfield,
            caltable=caltablename+"_spw"+str(i)+".B",
            bandtype="B",
            solint="inf",
            combine="scan,obs",
            minblperant=2,
            solnorm=False,
            minsnr=-2, # i.e. use all solutions, new parameter in 3.0.2
            gaintable=caltablename+"_spw"+str(i)+".K",
            spwmap=[i]
            )

    print(">> Plot the bandpass solutions")
    pl.clf()
    for i in range(2):
        plotcal(
            caltable=caltablename+"_spw"+str(i)+".B",xaxis="chan",yaxis="phase",
            spw="",subplot=221+i,overplot=False,plotsymbol='.',
            timerange="", showgui=False
            )
        tget(plotcal)
        subplot=223+i
        yaxis="amp"
        plotrange=[0,0,0,0.02]
        figfile=caltablename+"_bandpass.png"
        plotcal()
        
    # For GSPLINE solutions, have to do it for different spws separately
    
    print(">> Find G solutions")
    for i in range(2):
        print(">> SPW: ",i)
        os.system('rm -rf '+caltablename+'_spw'+str(i)+'.G')
        gaincal(
            vis=msn,
            caltable=caltablename+"_spw"+str(i)+".G",
            field=calfield,
            spw=avspw[i],
            selectdata=True,
            solint="60s",
            gaintable=[caltablename+'_spw'+str(i)+'.K',caltablename+'_spw'+str(i)+'.B'],
            spwmap=[[i],[i]],
            combine="",refant="0",minblperant=2,minsnr=-1,solnorm=False,
            gaintype="G",calmode="ap",
            )

    print(">> Find G solutions, using WVR corrections")
    for i in range(2):
        print(">> SPW: ",i)
        os.system('rm -rf '+caltablename+'_spw'+str(i)+'.G_WVR')
        wvrspw = 0
        if(cu.compare_version('>=',[3,4,0])):
            wvrspw = i
        gaincal(
            vis=msn,
            caltable=caltablename+"_spw"+str(i)+".G_WVR",
            field=calfield,
            spw=avspw[i],
            selectdata=True,
            solint="60s",
            # Pre-apply the WVR correction, the bandpass and the delay correction
            gaintable=[caltablename0,caltablename+'_spw'+str(i)+'.K',caltablename+'_spw'+str(i)+'.B'],
            spwmap=[[wvrspw],[i],[i]],
            combine="",refant="0",minblperant=2,minsnr=-1,solnorm=False,
            gaintype="G",calmode="ap",
            )	


    print(">> Find GSPLINE solutions, using WVR corrections")
    for i in range(2):
        print(">> SPW: ",i)
        os.system('rm -rf '+caltablename+'.GSPLINE_WVR_'+str(i))
        wvrspw = 0
        if(cu.compare_version('>=',[3,4,0])):
            wvrspw = i
        gaincal(
            vis=msn,
            caltable=caltablename+".GSPLINE_WVR_"+str(i),
            field=calfield,
            spw=avspw[i],
            selectdata=True,
            solint="300s",
            splinetime=10000,
            # Pre-apply the WVR correction, the bandpass and the delay correction
            gaintable=[caltablename0,caltablename+'_spw'+str(i)+'.K',caltablename+'_spw'+str(i)+'.B'],
            spwmap=[[wvrspw],[i],[i]],
            combine="",refant="0",minblperant=2,minsnr=-1,solnorm=False,
            gaintype="GSPLINE",
            npointaver=2,
            preavg=300,
            calmode="ap",
            )	


    print("\n>> Plot the solutions: small points: uncorrected")
    print(">>                     large points: WVR corrected")
    print(">>                     lines: spline fits, WVR corrected")
    pl.clf()
    for i in range(2):
        plotcal(
            caltable=caltablename+"_spw"+str(i)+".G_WVR",xaxis="time",yaxis="phase",
            spw="",subplot=211+i,overplot=False,plotsymbol='o',
            timerange="", showgui=False
            )
        tget(plotcal)
        overplot=True
        plotsymbol='.'; caltable=caltablename+"_spw"+str(i)+".G"; plotcal()
        plotsymbol=':'; caltable=caltablename+".GSPLINE_WVR_"+str(i);
        figfile=caltablename+"_gaincals.png"
        plotcal()


    # Apply the gain calibrations
    
    # Have to apply the G_WVR table and the WVR table simultaneously
    # Or just the G table
    # (Do not apply just the G_WVR table or G and WVR tables)
    
    # Use spwmap to apply to all spws
    # Use spwmap=[[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]
    # for the unsplit ms.
    # There are eight numbers because wvr has four spws, although only one is shown in listobs
    
    print(">> Apply the G and WVR solutions to",msn)
    for i in range (2):
        print(">> SPW: ",i)
        wvrspw = 0
        if(cu.compare_version('>=',[3,4,0])):
            wvrspw = i
        applycal(
            vis=msn,
            field="",
            spw=avspw[i],
            selectdata=False,timerange="",uvrange="",
            antenna="",scan="",msselect="",
            parang=False,calwt=True,
            # We use G after all, not Gspline
            gaintable=[caltablename+"_spw"+str(i)+".G_WVR",caltablename0,caltablename+"_spw"+str(i)+".K",caltablename+"_spw"+str(i)+".B"],
            gainfield=['','','',''],
            interp=['linear','','',''],
            spwmap=[[i],[wvrspw],[i],[i]]
            )

    # Clean the calibrator

    rms = [0.,0.]
    peak = [0.,0.]
    
    for i in range (2):
        
        # tb.open(msn0+'/FIELD')
        # coords= tb.getcol('PHASE_DIR')
        # calCoords=coords[:,:,int(calfield)]
        # calRa=qa.formxxx(str(float(calCoords[0]))+'rad','hms')
        # calDec=qa.formxxx(str(float(calCoords[1]))+'rad','dms')
        
        # print ">> The coordinates of the phase calibrator are:", calRa,",",calDec
        
        calimage=caltablename+"cal_spw"+str(i)
        
        os.system('rm -rf '+calimage+'.*')
	
        print("\n>> Clean the phase calibrator, spw",str(i))
        print(">> The imagename is",calimage)
        
        default(clean)
        clean(
            vis=msn,
            imagename=calimage,
            field=calfield,
            spw=avspw[i],
            niter=500,
            gain=0.15,
            threshold="0.0mJy",
            psfmode="hogbom",
            imagermode='',
            interactive=False,
            mask=[298,298,302,302],
            nchan=-1, start=0, width=1,
            imsize=[600, 600],
            cell=['0.2arcsec','0.2arcsec']
            )
        
        calstat=imstat(imagename=calimage+".image",region="",box="100,100,240,500")
        rms[i]=(calstat['rms'][0])
        print(">> rms in calibrator image: "+str(rms[i]))
        calstat=imstat(imagename=calimage+".image",region="")
        peak[i]=(calstat['max'][0])
        print(">> Peak in calibrator image: "+str(peak[i]))
        print(">> Dynamic range in calibrator image: "+str(peak[i]/rms[i]))

    #reference_rms = [0.0016516, 0.0009388]
    #reference_peak = [1.00014091, 1.00002384]
    # the follow test values introduced w/ data weight updates
    reference_rms = [0.0019042, 0.0008856]
    reference_peak = [0.999954, 1.0001681]

    for i in range(2):
        print(">> image rms ", i, " is ", rms[i], " expected value is ", reference_rms[i]) 
        if(abs(rms[i] - reference_rms[i])/reference_rms[i] > 0.01):
            print(">> ERROR.") 
            isOK = False
        print(">> image peak ", i, " is ", peak[i], " expected value is ", reference_peak[i]) 
        if(abs(peak[i] - reference_peak[i])/reference_peak[i] > 0.01):
            print(">> ERROR.") 
            isOK = False
            
    if (isOK):
        print('')
        print('Regression PASSED')
        print('')
    else:
        print('')
        print('Regression FAILED')
        print('')

    return isOK


###########################
# beginning of actual test 

part1 = True

# part 1
try:
    importasdm(myasdm_dataset_name, useversion='v3')
except:
    print(myname, ": Error ", sys.exc_info()[0])
    part1 = False
else:
    print(myname, ": Success! Now checking output ...")
    # when SubMS::setupMS() is used, the main table
    # somehow starts from table.f1 instead of table.f0
    mscomponents = set(["table.dat",
                        #"table.f0",
                        "table.f1",
                        "table.f2",
                        "table.f3",
                        "table.f4",
                        "table.f5",
                        "table.f6",
                        "table.f7",
                        "table.f8",
                        "ANTENNA/table.dat",
                        "DATA_DESCRIPTION/table.dat",
                        "FEED/table.dat",
                        "FIELD/table.dat",
                        "FLAG_CMD/table.dat",
                        "HISTORY/table.dat",
                        "OBSERVATION/table.dat",
                        "POINTING/table.dat",
                        "POLARIZATION/table.dat",
                        "PROCESSOR/table.dat",
                        "SOURCE/table.dat",
                        "SPECTRAL_WINDOW/table.dat",
                        "STATE/table.dat",
                        "SYSCAL/table.dat",
                        "ANTENNA/table.f0",
                        "DATA_DESCRIPTION/table.f0",
                        "FEED/table.f0",
                        "FIELD/table.f0",
                        "FLAG_CMD/table.f0",
                        "HISTORY/table.f0",
                        "OBSERVATION/table.f0",
                        "POINTING/table.f0",
                        "POLARIZATION/table.f0",
                        "PROCESSOR/table.f0",
                        "SOURCE/table.f0",
                        "SPECTRAL_WINDOW/table.f0",
                        "STATE/table.f0",
                        "SYSCAL/table.f0"
                        ])
    for name in mscomponents:
        if not os.access(msname+"/"+name, os.F_OK):
            print(myname, ": Error  ", msname+"/"+name, "doesn't exist ...")
            part1 = False
        else:
            print(myname, ": ", name, "present.")
    print(myname, ": MS exists. All tables present. Try opening as MS ...")
    try:
        ms.open(msname)
    except:
        print(myname, ": Error  Cannot open MS table", tablename)
        part1 = False
    else:
        ms.close()
        print(myname, ": OK. Checking tables in detail ...")

        # check main table first
        name = ""
        #             col name, row number, expected value, tolerance
        expected = [
                     ['UVW',       42, [ 0., 0., 0. ], 1E-7],
                     ['EXPOSURE',  42, 2.016, 0],
                     ['DATA',      42, [ [ 0.11164435+0.j],[ 1.01130176+0.j] ], 1E-7]
                     ]
        checktable(name, expected)

        expected = [
                     ['UVW',       557, [-172.89958856,  -80.72977367,   63.39288203], 1E-6], # based on geodetic data r3846 
                     ['EXPOSURE',  557, 6.048, 0],
                     ['DATA',      557,
                      [[ -4.17647697e-03 +3.08606686e-05j,  -1.18642126e-03 +7.54371868e-05j,
                         9.60109683e-05 -1.37158531e-05j,  -5.48634125e-05 -7.20082244e-05j,
                         7.54371868e-05 -1.37158531e-05j,   5.48634125e-05 +2.74317063e-05j,
                         2.74317063e-05 +3.77185934e-05j,  -1.02868898e-05 +3.08606686e-05j,
                         1.37158531e-05 -2.40027421e-05j,   5.82923749e-05 -6.17213373e-05j,
                         9.25820059e-05 -5.82923749e-05j,   8.22951188e-05 -3.42896328e-06j,
                         8.57240811e-05 +3.42896310e-05j,   1.13155787e-04 -3.42896328e-06j,
                         9.94399306e-05 -6.85792656e-06j,   6.17213373e-05 +4.45765218e-05j,
                         2.40027421e-05 +2.40027421e-05j,  -2.74317063e-05 -5.82923749e-05j,
                         -3.08606686e-05 -7.20082244e-05j,   3.42896328e-06 -1.71448155e-05j,
                         -2.40027421e-05 -6.85792656e-06j,  -7.54371868e-05 -2.40027421e-05j,
                         -4.80054841e-05 -1.02868898e-05j,  -1.37158531e-05 +3.42896310e-05j,
                         0.00000000e+00 +3.77185934e-05j,   6.85792656e-06 +2.74317063e-05j,
                         -1.37158531e-05 +7.54371868e-05j,  -1.02868898e-05 +8.22951188e-05j,
                         5.82923749e-05 +1.37158531e-05j,   5.48634125e-05 +3.42896328e-06j,
                         -6.85792656e-06 +2.40027421e-05j,  -3.77185934e-05 +3.42896310e-05j,
                         -4.11475594e-05 +1.37158531e-05j,  -3.42896328e-06 +2.05737797e-05j,
                         1.71448155e-05 +5.48634125e-05j,   1.37158531e-05 +2.05737797e-05j,
                         2.40027421e-05 -4.45765218e-05j,   7.88661564e-05 -5.48634125e-05j,
                         8.91530435e-05 -3.77185934e-05j,   2.05737797e-05 -3.77185934e-05j,
                         1.02868898e-05 +3.42896328e-06j,   1.02868898e-05 +5.48634125e-05j,
                         2.05737797e-05 +7.54371868e-05j,   7.20082244e-05 +5.82923749e-05j,
                         4.11475594e-05 +6.51502996e-05j,  -5.48634125e-05 +1.16584750e-04j,
                         -9.25820059e-05 +8.57240811e-05j,  -7.54371868e-05 +3.42896310e-05j,
                         -4.80054841e-05 +2.74317063e-05j,  -1.02868898e-05 +0.00000000e+00j,
                         5.48634125e-05 +0.00000000e+00j,   1.13155787e-04 +1.71448155e-05j,
                         1.09726825e-04 +2.05737797e-05j,   8.57240811e-05 +3.77185934e-05j,
                         3.08606686e-05 +4.45765218e-05j,  -4.45765218e-05 +4.45765218e-05j,
                         -2.05737797e-05 +7.88661564e-05j,  -6.85792656e-06 +8.91530435e-05j,
                         -4.80054841e-05 +6.85792620e-05j,  -3.42896328e-06 +5.14344465e-05j,
                         4.45765218e-05 +1.02868898e-05j,   5.14344465e-05 -5.14344465e-05j,
                         6.17213373e-05 +1.37158531e-05j,   5.82923749e-05 +1.44016449e-04j,
                         4.11475594e-05 +7.88661564e-05j,   3.42896310e-05 -1.02868898e-05j,
                         2.74317063e-05 +1.37158531e-05j,   2.05737797e-05 +6.85792656e-06j,
                         1.37158531e-05 -3.42896328e-06j,   2.05737797e-05 +0.00000000e+00j,
                         4.45765218e-05 +0.00000000e+00j,   9.60109683e-05 -3.42896328e-06j,
                         8.57240811e-05 +6.85792656e-06j,  -3.42896310e-05 +2.40027421e-05j,
                         -6.85792620e-05 -2.05737797e-05j,  -1.02868898e-05 -7.20082244e-05j,
                         0.00000000e+00 -8.22951188e-05j,   1.02868898e-05 -2.40027421e-05j,
                         3.77185934e-05 +6.85792656e-06j,   1.71448155e-05 -2.40027421e-05j,
                         -6.17213373e-05 +2.05737797e-05j,  -7.20082244e-05 +2.05737797e-05j,
                         1.37158531e-05 -2.05737797e-05j,   3.77185934e-05 -2.40027421e-05j,
                         4.11475594e-05 -3.08606686e-05j,   9.60109683e-05 -5.14344465e-05j,
                         1.09726825e-04 -7.88661564e-05j,   3.08606686e-05 -7.88661564e-05j,
                         -4.11475594e-05 -2.40027421e-05j,  -2.05737797e-05 +4.45765218e-05j,
                         1.37158531e-05 +2.74317063e-05j,  -2.74317063e-05 -1.71448155e-05j,
                         -5.48634125e-05 -1.02868898e-05j,  -2.74317063e-05 +4.80054841e-05j,
                         0.00000000e+00 +1.16584750e-04j,  -6.85792656e-06 +9.94399306e-05j,
                         -4.45765218e-05 +1.71448155e-05j,  -4.11475594e-05 +3.42896328e-06j,
                         0.00000000e+00 +3.77185934e-05j,   3.08606686e-05 +2.40027421e-05j,
                         3.77185934e-05 +1.37158531e-05j,   5.14344465e-05 +7.88661564e-05j,
                         6.85792620e-05 +1.26871644e-04j,   6.51502996e-05 +8.91530435e-05j,
                         6.51502996e-05 +4.45765218e-05j,   3.77185934e-05 +5.48634125e-05j,
                         1.71448155e-05 +7.20082244e-05j,   1.02868898e-05 +3.77185934e-05j,
                         2.40027421e-05 +2.05737797e-05j,   5.82923749e-05 +1.71448155e-05j,
                         3.42896328e-06 +3.42896328e-06j,  -6.17213373e-05 -1.02868898e-05j,
                         -4.80054841e-05 -3.42896310e-05j,  -6.17213373e-05 -4.80054841e-05j,
                         -7.54371868e-05 -4.45765218e-05j,   0.00000000e+00 -2.40027421e-05j,
                         2.05737797e-05 -2.40027421e-05j,  -8.22951188e-05 -3.42896328e-06j,
                         -8.22951188e-05 +2.05737797e-05j,   3.42896328e-06 +2.40027421e-05j,
                         -2.40027421e-05 +7.54371868e-05j,  -5.48634125e-05 +9.60109683e-05j,
                         2.40027421e-05 +2.05737797e-05j,   7.88661564e-05 -3.42896328e-06j,
                         4.45765218e-05 +7.54371868e-05j,  -3.08606686e-05 +8.22951188e-05j,
                         -5.48634125e-05 +1.02868898e-05j,   1.02868898e-05 -1.02868898e-05j
                         ],
                       [ -2.24771962e-01 +1.02868898e-05j,  -6.19887970e-02 +6.17213373e-05j,
                         8.64784513e-03 +8.22951188e-05j,  -2.27340264e-03 +0.00000000e+00j,
                         9.36106953e-04 -4.45765218e-05j,  -5.21202397e-04 -1.71448155e-05j,
                         2.81174987e-04 -3.77185934e-05j,  -2.16024680e-04 -4.80054841e-05j,
                         3.08606686e-05 -1.37158531e-05j,  -2.09166756e-04 +0.00000000e+00j,
                         2.05737797e-05 -3.42896328e-06j,  -4.80054841e-05 -6.17213373e-05j,
                         5.14344465e-05 -1.02868893e-04j,  -7.20082244e-05 -7.20082244e-05j,
                         -2.05737797e-05 -8.57240811e-05j,  -1.02868898e-05 -1.16584750e-04j,
                         5.48634125e-05 -9.60109683e-05j,   1.02868898e-05 -8.57240811e-05j,
                         6.85792620e-05 -2.05737797e-05j,  -1.02868898e-05 +1.06297855e-04j,
                         -6.85792656e-06 +7.20082244e-05j,  -3.42896328e-06 -2.05737797e-05j,
                         -3.77185934e-05 -2.05737797e-05j,  -6.51502996e-05 +4.11475594e-05j,
                         -6.51502996e-05 +1.09726825e-04j,  -1.50874374e-04 +8.91530435e-05j,
                         -8.91530435e-05 +2.74317063e-05j,   3.42896328e-06 +3.42896310e-05j,
                         3.42896328e-06 +3.77185934e-05j,  -5.48634125e-05 +2.40027421e-05j,
                         -6.17213373e-05 +5.82923749e-05j,  -6.17213373e-05 +7.88661564e-05j,
                         -1.71448155e-05 +2.05737797e-05j,  -2.40027421e-05 -5.82923749e-05j,
                         -4.80054841e-05 -7.88661564e-05j,  -1.02868898e-05 -1.37158531e-05j,
                         -2.74317063e-05 -2.74317063e-05j,  -9.94399306e-05 -1.23442675e-04j,
                         -7.20082244e-05 -9.25820059e-05j,  -4.80054841e-05 -3.42896328e-06j,
                         -2.74317063e-05 -3.77185934e-05j,  -5.82923749e-05 -6.85792620e-05j,
                         -1.20013712e-04 -4.80054841e-05j,  -9.94399306e-05 -6.17213373e-05j,
                         -4.45765218e-05 -5.14344465e-05j,  -8.22951188e-05 -1.71448155e-05j,
                         -7.20082244e-05 -4.45765218e-05j,   1.02868898e-05 -6.17213373e-05j,
                         3.42896310e-05 +1.71448155e-05j,   1.02868898e-05 +7.54371868e-05j,
                         0.00000000e+00 +1.71448155e-05j,  -1.02868898e-05 -4.80054841e-05j,
                         3.08606686e-05 -4.80054841e-05j,   4.11475594e-05 -1.02868898e-05j,
                         4.45765218e-05 +2.05737797e-05j,   3.08606686e-05 +6.17213373e-05j,
                         -4.45765218e-05 +8.91530435e-05j,  -7.54371868e-05 +5.82923749e-05j,
                         -5.48634125e-05 -4.80054841e-05j,  -2.05737797e-05 -1.30300599e-04j,
                         4.45765218e-05 -9.25820059e-05j,   0.00000000e+00 -2.40027421e-05j,
                         -7.54371868e-05 +1.02868898e-05j,  -7.88661564e-05 +3.42896310e-05j,
                         -2.05737797e-05 +2.05737797e-05j,   1.37158531e-05 -5.48634125e-05j,
                         1.37158531e-05 -8.22951188e-05j,   2.74317063e-05 -4.80054841e-05j,
                         1.71448155e-05 -5.82923749e-05j,  -1.71448155e-05 -6.17213373e-05j,
                         2.74317063e-05 -5.48634125e-05j,   6.17213373e-05 -9.60109683e-05j,
                         -1.02868898e-05 -1.06297855e-04j,  -7.20082244e-05 -8.57240811e-05j,
                         1.02868898e-05 -8.22951188e-05j,   7.88661564e-05 -5.48634125e-05j,
                         1.02868893e-04 -2.74317063e-05j,   1.23442675e-04 -3.42896328e-06j,
                         9.94399306e-05 +3.77185934e-05j,   3.42896310e-05 +5.14344465e-05j,
                         2.40027421e-05 +4.11475594e-05j,   3.77185934e-05 +0.00000000e+00j,
                         -3.42896328e-06 -5.48634125e-05j,  -1.37158531e-05 -8.91530435e-05j,
                         3.77185934e-05 -8.57240811e-05j,   2.40027421e-05 -5.82923749e-05j,
                         -1.02868898e-05 -6.17213373e-05j,   3.77185934e-05 -3.77185934e-05j,
                         4.11475594e-05 +3.42896310e-05j,  -2.40027421e-05 -2.05737797e-05j,
                         -2.74317063e-05 -1.09726825e-04j,  -2.40027421e-05 -8.57240811e-05j,
                         3.42896310e-05 -2.74317063e-05j,   1.13155787e-04 -1.37158531e-05j,
                         1.30300599e-04 -3.77185934e-05j,   1.13155787e-04 -6.17213373e-05j,
                         8.57240811e-05 -3.42896310e-05j,   6.51502996e-05 +2.05737797e-05j,
                         4.45765218e-05 +6.85792656e-06j,  -1.37158531e-05 -8.57240811e-05j,
                         -7.54371868e-05 -8.91530435e-05j,  -3.42896310e-05 -3.08606686e-05j,
                         3.42896328e-06 -5.48634125e-05j,  -3.77185934e-05 -1.16584750e-04j,
                         -2.40027421e-05 -9.60109683e-05j,  -2.40027421e-05 -3.08606686e-05j,
                         -4.45765218e-05 +6.51502996e-05j,   2.40027421e-05 +1.06297855e-04j,
                         7.20082244e-05 +2.05737797e-05j,  -1.71448155e-05 -6.51502996e-05j,
                         -1.02868893e-04 -5.82923749e-05j,  -1.16584750e-04 -4.11475594e-05j,
                         -4.45765218e-05 -1.71448155e-05j,   4.45765218e-05 +3.42896328e-06j,
                         1.71448155e-05 +1.37158531e-05j,  -2.74317063e-05 +3.77185934e-05j,
                         -1.02868898e-05 +4.45765218e-05j,   3.42896328e-06 +1.71448155e-05j,
                         3.42896328e-06 +2.74317063e-05j,  -6.17213373e-05 +5.48634125e-05j,
                         -9.60109683e-05 +3.42896328e-06j,  -3.42896310e-05 -2.40027421e-05j,
                         -3.42896310e-05 +1.71448155e-05j,  -8.91530435e-05 -2.40027421e-05j,
                         0.00000000e+00 -5.82923749e-05j,   1.57732313e-04 -3.42896328e-06j,
                         1.47445418e-04 +1.71448155e-05j,   1.02868893e-04 -1.02868898e-05j
                         ]
                       ],
                      1E-7
                      ]
                     ]
        checktable(name, expected)
        
        name = "ANTENNA"
        expected = [ ['OFFSET',       1, [ 0.,  0.,  0.], 0],
                     #['POSITION',     1, [2224600.1130,  -5440364.6494, -2481417.4709], 0.001],
                     ['POSITION',     1,  [2224602.5538,  -5440370.6185, -2481420.1935], 0.001],
                     ['DISH_DIAMETER',1, 12.0, 0]
                     ]
        checktable(name, expected)
        
        name = "POINTING"
        expected = [ ['DIRECTION',       10, [[-0.97039579],[ 0.88554736]], 1E-7],
                     ['INTERVAL',        10, 0.048, 0],
                     ['TARGET',          10, [[-0.97039579],[ 0.88554736]], 1E-7],
                     ['TIME',            10, 4775780726.88, 0.01],
                     ['TIME_ORIGIN',     10, 0., 0],
                     ['POINTING_OFFSET', 10, [[ 0.],[ 0.]], 0],
                     ['ENCODER',         10, [-0.96686338, 0.88392100], 1E-7 ]
                     ]
        checktable(name, expected)

if (not part1):
    print("Part 1 failed.")

part2 = True
        
myvis = myms_dataset_name
os.system('rm -rf exportasdm-output.asdm myinput.ms')
os.system('cp -R ' + myvis + ' myinput.ms')
default('exportasdm')
try:
    print("\n>>>> Test of exportasdm: input MS  is ", myvis)
    print("(a simulated input MS with pointing table)")
    rval = exportasdm(
        vis = 'myinput.ms',
        asdm = 'exportasdm-output.asdm',
        archiveid="S002",
        apcorrected=False,
        useversion='v3'
        )
    print("rval is ", rval)
    if not rval:
        raise Exception
    os.system('rm -rf '+asdmname+'; mv exportasdm-output.asdm '+asdmname)
    verify_asdm(asdmname, True)
except:
    print(myname, ': *** Unexpected error exporting MS to ASDM, regression failed ***')   
    raise
    
try:
    print("Reimporting the created ASDM ....")
    importasdm(asdm=asdmname, vis=reimp_msname, wvr_corrected_data='no', useversion='v3')
    print("Testing existence of reimported MS ....")
    if(not os.path.exists(reimp_msname)):
        print("MS ", reimp_msname, " doesn't exist.")
        raise Exception
    print("Testing equivalence of the original and the reimported MS.")
    tb.open(myms_dataset_name)
    nrowsorig = tb.nrows()
    print("Original MS contains ", nrowsorig, "integrations.")
    tb.close()
    tb.open(reimp_msname)
    nrowsreimp = tb.nrows()
    print("Reimported MS contains ", nrowsreimp, "integrations.")
    if(not nrowsreimp==nrowsorig):
        print("Numbers of integrations disagree.")
        part2 = False
except:
    print(myname, ': *** Unexpected error reimporting the exported ASDM, regression failed ***')   
    part2 = False
    
#############
# Now import an ASDM and do a serious analysis

print()
print('===================================================================')
print("Serious analysis of an ASDM ...")

rval = True
part3 = True

try:
    rval = analyseASDM(myasdm_dataset2_name, mywvr_correction_file)
except:
    print(myname, ': *** Unexpected error analysing ASDM, regression failed ***')   
    part3 = False

if(not rval):
    print(myname, ': *** Unexpected error analysing ASDM, regression failed ***')   
    part3 = False

#############
#

dopart4 = True

part4 = True

if dopart4:

    print()
    print('===================================================================')
    print("Export the previously imported ASDM ...")
    default('exportasdm')
    try:
        # os.system('rm -rf '+myasdm_dataset2_name+'-re-exported*')
        os.system('rm -rf '+myasdm_dataset2_name+'-re-exported* '+myasdm_dataset2_name+'-split*')
        
        # use only the actual visibility data, not the WVR data
        split(vis=myasdm_dataset2_name+'.ms',
              outputvis=myasdm_dataset2_name+'-split.ms',
              datacolumn='data',
              spw='0~3'
              )
        
        rval = exportasdm(
#            vis = myasdm_dataset2_name+'.ms',
            vis = myasdm_dataset2_name+'-split.ms',
            asdm = myasdm_dataset2_name+'-re-exported',
            archiveid="X001",
            apcorrected=False,
            datacolumn='DATA', # the default
            useversion='v3'
            )
        print("rval is ", rval)
        if not rval:
            raise Exception
        verify_asdm(myasdm_dataset2_name+'-re-exported', True)
    except:
        print(myname, ': *** Unexpected error re-exporting MS to ASDM, regression failed ***')   
        raise

#############
    # repeat analysis on re-exported ASDM

    print()
    print('===================================================================')
    print("Serious analysis of an exported and re-imported ASDM ...")
    
    rval = True
    try:
        rval = analyseASDM(myasdm_dataset2_name+'-re-exported', mywvr_correction_file,
                           False # do not regenerate the WVR table
                           )
    except:
        print(myname, ': *** Unexpected error analysing re-exported ASDM, regression failed ***')   
        part4 = False
        
    if(not rval):
        print(myname, ': *** Unexpected error analysing re-exported ASDM, regression failed ***')   
        part4 = False

# end if dopart4

print()
print('===================================================================')
print("Summary")

if(not part1):
    print("Part 1: ASDM import failed.")
else:
    print("Part 1: ASDM import passed.")
if(not part2):
    print("Part 2: ASDM export failed.")
else:
    print("Part 2: ASDM export passed.")
if(not part3):
    print("Part 3: serious analysis of an ASDM failed.")
else:
    print("Part 3: serious analysis of an ASDM passed.")
if dopart4:
    if not part4:
        print("Part 4: serious analysis of an exported and re-imported ASDM failed.")
    else:
        print("Part 4: serious analysis of an exported and re-imported ASDM passed.")
else:
        print("Part 4: serious analysis of an exported and re-imported ASDM not executed.")    

if(not (part1 and part2 and part3 and part4)):
    print("Regression failed.")
    raise
else:
    print("Regression passed.")
