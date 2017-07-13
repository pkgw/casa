from TestUtil import print_test, get_band_center_width
TestName = "RefractiveIndexProfileTest"

myat = casac.atmosphere()
myqa = casac.quanta()
tmp = myqa.quantity(270.0,'K')
pre = myqa.quantity(560.0,'mbar')
hum = 20.0
alt = myqa.quantity(5000,'m')
h0  = myqa.quantity(2.0,'km')
wvl = myqa.quantity(-5.6,'K/km')
mxA = myqa.quantity(48,'km')
dpr = myqa.quantity(5.0,'mbar')
dpm = 1.1
att = 1
myatm = myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
print_test(myatm, TestName)
print_test("  ", TestName)
print_test("  ", TestName)
print_test("Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries", TestName)
print_test("  ", TestName)
p = myat.getProfile()
print_test(p[0], TestName)
w = myat.getGroundWH2O()
print_test("First guess precipitable water vapor content: %f %s" % \
		( w['value'][0], w['unit']), TestName)
print_test("(This value is estimated from the relative humidity at ground level and the water vapor scale height)", TestName)
print_test("  ", TestName)
print_test("  ", TestName)

print_test("Example 1: Absorption profile for a single frequency: 850 GHz", TestName)
nb = 1
fC = myqa.quantity([850.0],'GHz')
fW = myqa.quantity([1.0],'Hz')
fR = myqa.quantity([0.0],'Hz')  # force single frequency
myat.initSpectralWindow(nb,fC,fW,fR)
do=myat.getDryOpacitySpec()
rchan=myat.getRefChan()
print("")
print_test("Absorption Profile built from RefractiveIndexProfile CONSTRUCTOR. Summary of results:", TestName)
print_test("Total Dry Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getDryOpacity()), TestName)
print("")
print_test("Total Dry Cont Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getDryContOpacity()), TestName)
print_test("Total O2 lines Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getO2LinesOpacity()), TestName)
print_test("Total O3 lines Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getO3LinesOpacity()), TestName)
print_test("Total CO lines Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getCOLinesOpacity()), TestName)
print_test("Total N2O lines Opacity at %f %s for 1.0 air mass: %f" % \
		(fC['value'][0], fC['unit'], myat.getN2OLinesOpacity()), TestName)
print("")
wo=myat.getWetOpacitySpec()
rchan=myat.getRefChan()
print_test("Total Wet Opacity at %f %s for 1.0 air mass: %f %s" % (fC['value'][0], fC['unit'], myat.getWetOpacity()['value'][0], myat.getWetOpacity()['unit']), TestName)
print("")
print_test("Total H2O lines Opacity at %f %s for 1.0 air mass: %f" % (fC['value'][0], fC['unit'], myat.getH2OLinesOpacity()), TestName)
print_test("Total H2O Cont Opacity at %f %s for 1.0 air mass: %f" % (fC['value'][0], fC['unit'], myat.getH2OContOpacity()), TestName)
print("")
print("")
print_test("Total Dispersive PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getDispersiveWetPathLength()['value'][0], myat.getDispersiveWetPathLength()['unit']), TestName)
print_test("Total Non-Dispersive PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getNonDispersiveWetPathLength()['value'][0], myat.getNonDispersiveWetPathLength()['unit']), TestName)
print_test("Ratio Dispersive/Non-Dispersive PathLength at %f %s for 1.0 air mass: %f" % \
     (fC['value'][0], fC['unit'], myat.getDispersiveWetPathLength()['value'][0]/myat.getNonDispersiveWetPathLength()['value'][0]), TestName)

print_test("Total Dry PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getNonDispersiveDryPathLength()['value'][0], myat.getNonDispersiveDryPathLength()['unit']), TestName)
print_test("Total O2 lines PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getO2LinesPathLength()['value'][0], myat.getO2LinesPathLength()['unit']), TestName)
print_test("Total O3 lines PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getO3LinesPathLength()['value'][0], myat.getO3LinesPathLength()['unit']), TestName)
print_test("Total CO lines PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getCOLinesPathLength()['value'][0], myat.getCOLinesPathLength()['unit']), TestName)
print_test("Total N2O lines PathLength at %f %s for 1.0 air mass: %f %s" % \
		(fC['value'][0], fC['unit'], myat.getN2OLinesPathLength()['value'][0], myat.getN2OLinesPathLength()['unit']), TestName)
print("")
print("")
print_test("(your actual water vapor column is %f %s)" % (myat.getGroundWH2O()['value'][0], myat.getGroundWH2O()['unit']), TestName)
print("")

# print_test("change a basic parameter", TestName)
# print_test("========================", TestName)
# print_test("Old ground temperature: %f %s" % (myat.getBasicAtmParms()[2]['value'][0], myat.getBasicAtmParms()[2]['unit']), TestName)
# new_tmp = myqa.quantity(275.0,'K')
# print_test("New ground temperature: %f %s" % (new_tmp['value'], new_tmp['unit']), TestName)
# print_test(myat.updateAtmProfile(alt, new_tmp, pre, hum, wvl, h0), TestName)
# print_test("Absorption Profile with this new temperature.  Summary of results:", TestName)
# do=myat.getDryOpacitySpec()
# print_test("Total Dry Opacity at %f %s for 1.0 air mass: %f" % (fC['value'][0], fC['unit'], myat.getDryOpacity()), TestName)
# wo=myat.getWetOpacitySpec()
# print_test("Total Wet Opacity at %f %s for 1.0 air mass: %f %s" % (fC['value'][0], fC['unit'], myat.getWetOpacity()['value'][0] / w['value'][0], myat.getWetOpacity()['unit']), TestName)
# print_test("Total Dispersive Delay at %f %s for 1.0 air mass: %f meters per mm of water vapor" % (fC['value'][0], fC['unit'], myat.getDispersivePathLength()['value'][0] /  w['value'][0]), TestName)
# print_test("(%f %% of the Non-dispersive one )" % (100*(myat.getDispersivePathLength()['value'][0] /  w['value'][0])/(myat.getNonDispersivePathLength()['value'][0] / w['value'][0])), TestName)
# print_test("(your actual water vapor column is %f %s)" % (myat.getGroundWH2O()['value'][0], myat.getGroundWH2O()['unit']), TestName)
# print("")
# 
# print_test("Add a spectral window", TestName)
# print_test("=====================", TestName)
# numChan = 4;
# refChan = 2;
# refFreq = myqa.quantity(284.97346,"GHz"); # 350
# chanSep = myqa.quantity(2.0,"MHz");
# fC2, fW2 = get_band_center_width(numChan, refChan, refFreq, chanSep)
# fR2 = chanSep
# 
# nc = myat.addSpectralWindow(fC2,fW2,fR2)
# #print_test("New spectral window has %d channels" % nc, TestName)
# #w = myat.getStartupWaterContent()
# w = myat.getGroundWH2O()
# numSpw = myat.getNumSpectralWindows()
# print_test("There are now %d spectral windows" % numSpw, TestName)
# print_test("Absorption profiles including this new spectral window.  Summary of results:", TestName)
# print_test("Total Dry Opacity at %f %s for 1.0 air mass: %f" % (fC['value'][0], fC['unit'], myat.getDryOpacity()), TestName)
# print("")
# for spwid in range(numSpw):
# 	numCh = myat.getNumChan(spwid)
# 	print_test("Spectral window %d has %d frequency channels" % (spwid, numCh), TestName)
# 	for n in range(numCh):
# 		freq = myat.getChanFreq(n, spwid)
# 		print_test("Total Wet Opacity at %f %s for 1.0 air mass: %f %s" % \
# 				(freq['value'][0], freq['unit'], myat.getWetOpacity(n,spwid)['value'][0] / w['value'][0], myat.getWetOpacity()['unit']), TestName)
# 		print_test("Total Non-Dispersive Delay at %f %s for 1.0 air mass: %f meters per %s of water vapor" % \
# 				(freq['value'][0], freq['unit'], myat.getNonDispersivePathLength(n,spwid)['value'][0]/w['value'][0], w['unit']), TestName)
# 		print_test("Total Dispersive Delay at at %f %s for 1.0 air mass: %f meters per %s of water vapor %f %% of the Non-dispersive one )" % \
# 				(freq['value'][0], freq['unit'], (myat.getDispersivePathLength(n,spwid)['value'][0])/(w['value'][0]), w['unit'],(100*myat.getDispersivePathLength(n,spwid)['value'][0]/myat.getNonDispersivePathLength(n,spwid)['value'][0])), TestName)
# 		#print_test("Total Dry Delay at %f %s for 1.0 air mass: %f %s" % (freq['value'][0], freq['unit'], myat.getDryPathLength(n,spwid)['value'][0] / w['value'][0], myat.getDryPathLength()['unit']), TestName)
# 		print_test("(your actual water vapor column is %f %s of water vapor)." % (w['value'][0], w['unit']), TestName)
# 		print("")
# print_test("=====================", TestName)
# for spwid in range(numSpw):
# 	numCh = myat.getNumChan(spwid)
# 	print_test("Spectral window %d has %d frequency channels" % (spwid, numCh), TestName)
# 	for n in range(numCh):
# 		freq = myat.getChanFreq(n, spwid)
# 		#print_test("Total Dispersive Phase Delay at ", freq['value'][0], freq['unit'], " for 1.0 air mass: ", (myat.getDispersivePhaseDelay(n,spwid)['value'][0])/(w['value'][0]), " degrees per mm of water vapor (", ((100*myat.getDispersivePhaseDelay(n,spwid)['value'][0])/(w['value'][0]))/(myat.getNonDispersivePhaseDelay(n,spwid)['value'][0]/w['value'][0]), "% of the Non-dispersive one )"
# 		print("")
#exit()
