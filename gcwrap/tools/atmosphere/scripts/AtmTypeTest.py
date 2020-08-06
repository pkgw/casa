myat = casac.atmosphere()
myqa = casac.quanta()

print("Test AtmType class")

alt = myqa.quantity(2550.,'m')
tmp = myqa.quantity(270.32,'K')
pre = myqa.quantity(73585.,'Pa')
mxA = myqa.quantity(45000.,'m')
hum = 20.0
wvl = myqa.quantity(-0.0056,'K/m')
dpr = myqa.quantity(500.,'Pa')
#dpm = 1.25
dpm = 1.5
h0  = myqa.quantity(2000.,'m')
att = 1
myatm=myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
l = myat.getBasicAtmParms()
print("Test for type ", l[10])
print(myatm)
p = myat.getProfile()
print("layer\taltitude\tthickness\ttemperature\twatermassdensity\tpressure")
height = 0
for i in range(myat.getNumLayers()):
	height += p[1]['value'][i]
	print(i, '\t', height, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])
print()

att = 2
myatm=myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
l = myat.getBasicAtmParms()
print("Test for type ", l[10])
print(myatm)
p = myat.getProfile()
print("layer\taltitude\tthickness\ttemperature\twatermassdensity\tpressure")
height = 0
for i in range(myat.getNumLayers()):
	height += p[1]['value'][i]
	print(i, '\t', height, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])
print()

att = 3
myatm=myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
l = myat.getBasicAtmParms()
print("Test for type ", l[10])
print(myatm)
p = myat.getProfile()
print("layer\taltitude\tthickness\ttemperature\twatermassdensity\tpressure")
height = 0
for i in range(myat.getNumLayers()):
	height += p[1]['value'][i]
	print(i, '\t', height, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])
print() 

att = 4
myatm=myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
l = myat.getBasicAtmParms()
print("Test for type ", l[10])
print(myatm)
p = myat.getProfile()
print("layer\taltitude\tthickness\ttemperature\twatermassdensity\tpressure")
height = 0
for i in range(myat.getNumLayers()):
	height += p[1]['value'][i]
	print(i, '\t', height, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])
print()

att = 5
myatm=myat.initAtmProfile(alt, tmp, pre, mxA, hum, wvl, dpr, dpm, h0, att)
l = myat.getBasicAtmParms()
print("Test for type ", l[10])
print(myatm)
p = myat.getProfile()
print("layer\taltitude\tthickness\ttemperature\twatermassdensity\tpressure")
height = 0
for i in range(myat.getNumLayers()):
	height += p[1]['value'][i]
	print(i, '\t', height, '\t', p[1]['value'][i],'\t', p[2]['value'][i],'\t', p[3]['value'][i],'\t', p[5]['value'][i])
#exit()
