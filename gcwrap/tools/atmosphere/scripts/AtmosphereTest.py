myat = casac.atmosphere()
myqa = casac.quanta()

alt = myqa.quantity(2550.,'m')
tmp = myqa.quantity(270.32,'K')
pre = myqa.quantity(73585.,'Pa')
mxA = myqa.quantity(45000,'m')
hum = 20.0
wvl = myqa.quantity(-0.0056,'K/m')
dpr = myqa.quantity(500.0,'Pa')
dpm = 1.25
h0  = myqa.quantity(2000.0,'m')
att = 2
print("Test of constructor Atmosphere(altitude,temperature,pressure,maxAltitude,humidity,dTempdH,dP,dPm,h0,MIDLAT_SUMMER)")
myatm = myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
print(myatm)

print()
print("Test: getGroundWH2O()")
w = myat.getGroundWH2O()
print("Guessed water content: ", w['value'][0], w['unit'])

print()
print("Test: getProfile()")
p = myat.getProfile()
print("altitude\tthickness\ttemperature\twatermassdensity\tpressure")
alt = 0
for i in range(myat.getNumLayers()):
	alt += p[1]['value'][i]
	print(alt, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])

print()
print("Test of initSpectralWindow()")
nbands=1
fC=myqa.quantity(88.,'GHz')
fW=myqa.quantity(.5,'GHz')
fR=myqa.quantity(.5,'GHz')
myat.initSpectralWindow(nbands,fC,fW,fR)
n=myat.getNumChan(0)
print(n, "channel(s) in band 0")

print()
print("Test: getOpacity()")
print(" - dryOpacity ", myat.getDryOpacity(), "wetOpacity ", myat.getWetOpacity()['value'][0],myat.getWetOpacity()['unit'])

print()
print("AbsCoeff getAbsCoeff() (3 first layers):", myat.getAbsH2OLines(0,0)['value'][0],myat.getAbsH2OLines(1,0)['value'][0],myat.getAbsH2OLines(2,0)['value'][0])

print()
print("Test of SkyBrightness calculations")
myat.setAirMass(1.51)
myat.setSkyBackgroundTemperature(myqa.quantity(2.73,'K'))
myat.setUserWH2O(myqa.quantity(4.05,'mm'))
print("SkyBrightness =", myat.getTebbSky()['value'][0], myat.getTebbSky()['unit'], " TEBB")
print("")
print("==========================================================")
print("Test with spectral data")
nb=2
fC=myqa.quantity([88.,90.],'GHz')
fW=myqa.quantity([.5,.5],'GHz')
fR=myqa.quantity([.125,.125],'GHz')
print("Test of initSpectralWindow")
myat.initSpectralWindow(nb,fC,fW,fR)
n=myat.getNumChan(0)
print(nb, " bands ", n, "channels(s)")
print("Test: Opacity getOpacitySpec")
for s in range(myat.getNumSpectralWindows()):
	print("band", s)
	for i in range(n):
		print(" - dryOpacity ", myat.getDryOpacitySpec(s)[1][i], " wet Opacity ", myat.getWetOpacitySpec(s)[1]['value'][i],myat.getWetOpacitySpec(s)[1]['unit'])
print("")
print("Test of SkyBrightness calculations")
myat.setAirMass(1.51)
myat.setSkyBackgroundTemperature(myqa.quantity(2.73,'K'))
myat.setUserWH2O(myqa.quantity(4.05,'mm'))
for s in range(myat.getNumSpectralWindows()):
	for i in range(n):
		print("Band", s, " channel ", i, "TebbSky = ", myat.getTebbSky(i,s)['value'][0], myat.getTebbSky()['unit'])
print("")
#exit()
