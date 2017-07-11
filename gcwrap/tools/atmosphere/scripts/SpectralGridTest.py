from TestUtil import print_test, get_band_center_width
TestName = "SpectralGridTest"

myat = casac.atmosphere()
myqa = casac.quanta()

myat.initAtmProfile()
#
numChan = 64
refChan = 32
myRefFreq = myqa.quantity(90.0,"GHz")
myChanSep = myqa.quantity(0.01,"GHz")
nb = 1
fC, fW = get_band_center_width(numChan, refChan, myRefFreq, myChanSep)
fR = myChanSep

print_test("Test 1:", TestName)
myat.initSpectralWindow(nb,fC,fW,fR)
print_test("Number of channels retrieved: %d Input: %d" % (myat.getNumChan(),numChan), TestName)
rf = myat.getRefFreq()
print_test("Reference frequency retrieved: %f %s Input: %f %s" % \
		(rf['value'][0], rf['unit'],
		myat.getChanFreq((numChan-1)/2, 0)['value'][0],
		myat.getChanFreq((numChan-1)/2, 0)['unit']), TestName)
rc = myat.getRefChan()
print_test("Reference channel retrieved: %d Input: %d" % \
		(rc, (numChan-1)/2), TestName)
cs = myat.getChanSep()
print_test("Channel separation retrieved: %f %s Input: %f %s" % (cs['value'][0], cs['unit'], fR['value'], fR['unit']), TestName)
print_test("Number of spectral windows: %d Input: %d" % (myat.getNumSpectralWindows(),nb), TestName)
temp_expected_bw = myqa.sub(fW,fR) # ATM lib returns maxFreq - minFreq
print_test("Total bandwidth retrieved: %f %s Input: %f %s" % \
		(myat.getBandwidth()['value'][0], myat.getBandwidth()['unit'],
		temp_expected_bw['value'], temp_expected_bw['unit']), TestName)
print_test("Frequency range: from %f to %f %s" % \
		(myat.getMinFreq()['value'][0], myat.getMaxFreq()['value'][0], myat.getMinFreq()['unit']), TestName)
print_test("Number of channels for spectral window identifier 0: %d" % myat.getNumChan(0), TestName)
print_test("Expect an error here", TestName)
try:
	myat.getNumChan(1)  # spwId is 0-based.
except Exception as e:
	print_test(e.message, TestName)
print_test("End of Test 1", TestName)
print("")

print_test("Test 2", TestName)
print_test("Adding a new spectral window", TestName)
numChan1  = 128;
refChan1 = 32;
myNewRefFreq = myqa.quantity(215.0,"GHz");
myNewChanSep = myqa.quantity(0.02,"GHz");
fC1, fW1 = get_band_center_width(numChan1, refChan1, myNewRefFreq, myNewChanSep)
fR1 = myNewChanSep
# #
# trc1 = 32 # reference channel
# trf1 = 215.0 #GHz frequency of trc1 channel.
# tfc1 = trf1 + (numChan1/2 - trc1)*chansep1
# #
# fC1 = myqa.quantity(tfc1, 'GHz')
# fR1 = myqa.quantity(chansep1, 'GHz')
# fW1 = myqa.quantity(numChan1*chansep1, 'GHz')
nc = myat.addSpectralWindow(fC1, fW1, fR1)
print_test("New spectral window has %d channels. Input: %d" % (nc, numChan1), TestName)
print_test("Number of spectral windows: %d Expected: %d" % (myat.getNumSpectralWindows(), nb+1), TestName)

spwId=0
print_test("Number of channels retrieved for spwID %d: %d Input: %d" % \
		(spwId, myat.getNumChan(spwId), numChan), TestName)
print_test("Reference frequency retrieved: %f %s Input: %f %s" % \
		(myat.getRefFreq(spwId)['value'][0],myat.getRefFreq(spwId)['unit'],
		myat.getChanFreq((numChan-1)/2, spwId)['value'][0],
		myat.getChanFreq((numChan-1)/2, spwId)['unit']), TestName)
print_test("Reference channel retrieved: %d Expected: %d" % \
		(myat.getRefChan(spwId), (numChan-1)/2), TestName)
print_test("Channel separation retrieved: %f %s Input: %f %s" % \
		(myat.getChanSep(spwId)['value'][0], myat.getChanSep(spwId)['unit'],
		fR['value'], fR['unit']), TestName)
spwId=1
print_test("Number of channels retrieved for spwID %d: %d Input: %d" % \
		(spwId, myat.getNumChan(spwId), numChan1), TestName)
print_test("Reference frequency retrieved: %f %s Input: %f %s" % \
		(myat.getRefFreq(spwId)['value'][0],myat.getRefFreq(spwId)['unit'],
		myat.getChanFreq((numChan1-1)/2, spwId)['value'][0],
		myat.getChanFreq((numChan1-1)/2, spwId)['unit']), TestName)
print_test("Reference channel retrieved: %d Expected: %d" % \
		(myat.getRefChan(spwId), (numChan1-1)/2), TestName)
print_test("Channel separation retrieved: %f %s Input: %f %s" % \
		(myat.getChanSep(spwId)['value'][0], myat.getChanSep(spwId)['unit'],
		fR1['value'], fR1['unit']), TestName)

print_test("Number of spectral windows: %d Expected: 2"% myat.getNumSpectralWindows(), TestName)
print_test("End of Test 2", TestName)


print_test("Channel frequency and number for the first spectral window: ", TestName)
n = myat.getNumChan()
rc = myat.getRefChan()
for i in range(n):
	cf = myat.getChanFreq(i)
	print_test("%d channel: %d freq: %f %s" % (i, i-rc, cf['value'][0], cf['unit']), TestName)

print("")
print_test("Test 2:", TestName)
print_test("Initializing SpectralGrid", TestName)
spwId=0
myat.initSpectralWindow(nb,fC,fW,fR)
print_test("Number of channels retrieved for spwID %d: %d Input: %d" % \
		(spwId, myat.getNumChan(spwId), numChan), TestName)
print_test("Reference frequency retrieved: %f %s Input: %f %s" % \
		(myat.getRefFreq(spwId)['value'][0],myat.getRefFreq(spwId)['unit'],
		myat.getChanFreq((numChan-1)/2, spwId)['value'][0],
		myat.getChanFreq((numChan-1)/2, spwId)['unit']), TestName)
print_test("Reference channel retrieved: %d Expected: %d" % \
		(myat.getRefChan(spwId), (numChan-1)/2), TestName)
print_test("Channel separation retrieved: %f %s Input: %f %s" % \
		(myat.getChanSep(spwId)['value'][0], myat.getChanSep(spwId)['unit'],
		fR['value'], fR['unit']), TestName)
print_test("Number of spectral windows: %d Expected: %d" % (myat.getNumSpectralWindows(),nb), TestName)
print("")
print_test("Skipping irregular spectral grid test...", TestName)
print_test("CASA atmosphere tool currently does not have ability to handle irregular spectral grid", TestName)


print("")
print("")
print_test("Test 3:", TestName)
print_test("Initializing SpectralGrid", TestName)
spwId=0
myat.initSpectralWindow(nb,fC,fW,fR)
print_test("Number of channels retrieved for spwID %d: %d Input: %d" % \
		(spwId, myat.getNumChan(spwId), numChan), TestName)
print_test("Reference frequency retrieved: %f %s Input: %f %s" % \
		(myat.getRefFreq(spwId)['value'][0],myat.getRefFreq(spwId)['unit'],
		myat.getChanFreq((numChan-1)/2, spwId)['value'][0],
		myat.getChanFreq((numChan-1)/2, spwId)['unit']), TestName)
print_test("Reference channel retrieved: %d Expected: %d" % \
		(myat.getRefChan(spwId), (numChan-1)/2), TestName)
print_test("Channel separation retrieved: %f %s Input: %f %s" % \
		(myat.getChanSep(spwId)['value'][0], myat.getChanSep(spwId)['unit'],
		fR['value'], fR['unit']), TestName)
print_test("Number of spectral windows: %d Expected: %d" % (myat.getNumSpectralWindows(),nb), TestName)

chan=16.123456
chanfreq=myqa.add(myqa.quantity(rf['value'][0],rf['unit']),myqa.mul(myqa.quantity(cs['value'][0],cs['unit']),chan))
print_test("Position (GU) retrieved: %f Exact: %f" % \
		(myat.getChanNum(myqa.quantity(chanfreq['value'],chanfreq['unit']),spwId), chan), TestName)
temp_expected_bw = myqa.sub(fW,fR) # ATM lib returns maxFreq - minFreq
print_test("Total bandwidth retrieved: %f %s Input: %f %s" % \
		(myat.getBandwidth(spwId)['value'][0], myat.getBandwidth(spwId)['unit'],
		temp_expected_bw['value'], temp_expected_bw['unit']), TestName)
print_test("Frequency range: from %f to %f %s" % \
		(myat.getMinFreq(spwId)['value'][0], myat.getMaxFreq(spwId)['value'][0],
		myat.getMinFreq(spwId)['unit']), TestName)

print("")
print("")
print_test("Skipping double sideband spectral grid test...", TestName)
print_test("CASA atmosphere tool currently does not have ability to handle double sideband spectra", TestName)

#exit()
