from TestUtil import print_test, get_band_center_width
TestName = "SkyStatusTest"

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
#print_test(" SkyStatusTest: size of Temperature Profile: %d" % TBD, TestName)
w = myat.getGroundWH2O()
print_test("First guess precipitable water vapor content: %f %s" % (w['value'][0], w['unit']), TestName)
print_test("(This value is estimated from the relative humidity at ground level and the water vapor scale height)", TestName)
print_test("  ", TestName)
print_test("  ", TestName)
print_test("TEST FOR JUST 1 FREQUENCY", TestName)
print_test("=========================", TestName)
print_test("  ", TestName)
nb = 1
fC = myqa.quantity([850.0],'GHz')
fW = myqa.quantity([0.005],'GHz')
fR = myqa.quantity([0.000],'GHz')  # force single frequency
myat.initSpectralWindow(nb,fC,fW,fR)
myat.setUserWH2O(myqa.quantity(1.0,'mm'))
for i in range(myat.getNumSpectralWindows()):
	for j in range(myat.getNumChan(i)):
		print_test("Frequency: %f %s" % \
				(myat.getChanFreq(j,i)['value'][0], myat.getChanFreq(j,i)['unit']), TestName)
		print_test("Wet opacity: %f %s for %f %s H2O" % \
				(myat.getWetOpacity(j,i)['value'][0], myat.getWetOpacity(j,i)['unit'],
				myat.getUserWH2O()['value'][0], myat.getUserWH2O()['unit']), TestName)
		print_test("Dry opacity: %f nepers" % (myat.getDryOpacity(j,i)), TestName)
		print_test("Sky brightness: %f %s" % \
				(myat.getAverageTebbSky(j)['value'][0], myat.getAverageTebbSky(j)['unit']), TestName)

wh2o=myqa.quantity(0.45,'mm')
print_test("(EXTERNAL CHANGE) water vapor column: %f %s" % (wh2o['value'], wh2o['unit']), TestName)
print_test("(NEW OUTPUT) T_EBB = %f %s" % (myat.getAverageTebbSky(0,wh2o)['value'][0], myat.getAverageTebbSky(0,wh2o)['unit']), TestName)
print_test("Current water vapor column in Radiance Class: %f %s" % (myat.getUserWH2O()['value'][0],myat.getUserWH2O()['unit']), TestName)
print_test("  ", TestName)

myat.setAirMass(2.0)
print_test("(INTERNAL CHANGE) Air mass: %f" % (myat.getAirMass()), TestName)
print_test("(NEW OUTPUT) T_EBB = %f %s" % (myat.getAverageTebbSky()['value'][0], myat.getAverageTebbSky()['unit']), TestName)
print("")
wh2o=myqa.quantity(0.8,"mm")
myat.setUserWH2O(wh2o)
print_test("(INTERNAL CHANGE) water vapor column: %f %s" % (wh2o['value'], wh2o['unit']), TestName)
print_test("(NEW OUTPUT) T_EBB = %f %s" % (myat.getAverageTebbSky()['value'][0], myat.getAverageTebbSky(0,wh2o)['unit']), TestName)
print("")
print_test("TEST FOR ONE SPECTRAL WINDOW WITH SEVERAL CHANNELS", TestName)
print_test("=====================================================", TestName)
print_test("  ", TestName)
nb2 = 1
numchan=5
refchan=3
chansep = myqa.quantity(0.01,'GHz')
reffreq = myqa.quantity(616.50,"GHz")
fC2, fW2 = get_band_center_width(numchan, refchan, reffreq, chansep)
fR2 = chansep
myat.initSpectralWindow(nb2,fC2,fW2,fR2)
print_test("water vapor column: %f %s" % (myat.getUserWH2O()['value'][0], myat.getUserWH2O()['unit']), TestName)
print_test("Air mass          : %f" % (myat.getAirMass()), TestName)
print("")
for i in range(myat.getNumChan(0)):
	print_test("Freq: %f %s / T_EBB= %f %s" % \
			(myat.getChanFreq(i)['value'][0], myat.getChanFreq(i)['unit'], myat.getTebbSky(i)['value'][0], myat.getTebbSky(i)['unit']), TestName)
print_test("  ", TestName)

myat.setAirMass(1.0)
myat.setUserWH2O(myqa.quantity(0.45,"mm"))
print_test("water vapor column: %f %s" % (myat.getUserWH2O()['value'][0], myat.getUserWH2O()['unit']), TestName)
print_test("Air mass          : %f" % (myat.getAirMass()), TestName)
print("")
for i in range(myat.getNumChan(0)):
	print_test("Freq: %f %s / T_EBB= %f %s  Dry opacity: %f np  Wet opacity: %f %s  Total opacity: %f %s" % \
			(myat.getChanFreq(i, 0)['value'][0], myat.getChanFreq(i, 0)['unit'],
			myat.getTebbSky(i, 0)['value'][0], myat.getTebbSky(i, 0)['unit'],
			myat.getDryOpacity(i,0),
			myat.getWetOpacity(i, 0)['value'][0], myat.getWetOpacity(i, 0)['unit'],
			myat.getDryOpacity(i,0)+myat.getWetOpacity(i, 0)['value'][0],
			myat.getWetOpacity(i, 0)['unit']), TestName)

#exit()
