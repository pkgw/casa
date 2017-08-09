"""
<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" ?>
<casaxml xmlns="http://casa.nrao.edu/schema/psetTypes.html"
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
xsi:schemaLocation="http://casa.nrao.edu/schema/casa.xsd
file:///opt/casa/code/xmlcasa/xml/casa.xsd">




        <tool name="measures" module="measures">
        <shortdescription>measures tool</shortdescription>


<keyword>measures</keyword>

<code>
        <include>xmlcasa/measures/measures_forward.h</include>
<private>
        #include <xmlcasa/measures/measures_private.h>
</private>
</code>



<!--
   <method type="constructor" name="measures">
   <shortdescription>Construct measures tool</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="host">
     <description>host on which to run tool</description>
     <value></value>
     </param>

     <param xsi:type="bool" direction="in" name="forceneweserver">
     <description>force the use of a new
server</description>
     <value>false</value>
     </param>

</input>
<returns xsi:type="measures"/>
<description>
Create a measures tool on the specified host (or by default the
host you are running on).
</description>
</method>
-->



   <method type="function" name="dirshow">
   <shortdescription>
<!-- Format a direction using globally set formats -->
Show direction measure as a string.
   </shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>a direction measure value to be converted to string
     </description>
     </param>
</input>
<returns xsi:type="string"/>
<description>
dirshow will convert a direction measure to a string
</description>
<!--
%%, using the formats set
%%globally for longitude and latitude (see \casa\
%%quanta:quanta.setformat.function).
-->
<example>
\begin{verbatim}
"""
#
print("\t----\t dirshow Ex 1 \t----")
print(me.dirshow(me.direction('venus')))
#[0, 90] deg  VENUS
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="show">
   <shortdescription> Show a measure as a string
   </shortdescription>
<!--
%%Format a measure using globally set formats
-->


<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure value to be converted to string</description>
     </param>

     <param xsi:type="bool" direction="in" name="refcode">
     <description>add the reference code to output</description>
     <value>true</value>
     </param>
</input>
<returns xsi:type="string"/>
<description>
show will convert a measure to a string.
<!--
%% , using the formats set
%%globally for various types of variables (see
%%<link anchor="quanta:quanta.setformat.function">setformat</link>
%%).
-->
All measures are catered for (at this moment {\em direction, position, epoch,
radialvelocity, frequency, doppler, baseline, uvw, earthmagnetic} ).
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t show Ex 1 \t----")
print(me.show(me.frequency('lsrk', qa.constants('HI'))))
#1.42041e+09 Hz LSRK
print(me.show(me.frequency('lsrk', qa.constants('HI')), refcode=false))
#1.42041e+09 Hz
#
"""
\end{verbatim}
</example>
<!--
%%#: qa.setformat('freq', 'keV')
%%#T
%%#: me.show(me.frequency('lsrk', qa.constants('HI')))
%%#5.87432838e-09 keV LSRK
-->
<!--
%%#5.87432838e-09 keV
-->
</method>




   <method type="function" name="epoch">
   <shortdescription>define an epoch measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>UTC</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>epoch value</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset epoch measure</description>
     <value></value>
     </param>

</input>
<returns xsi:type="any">
<shortdescription>epoch measure</shortdescription>
<any type="record"/>
</returns>
<description>
epoch defines an epoch measure from the CLI. It has to specify a
reference code, an epoch quantity value (see introduction for the
action on a scalar quantity with either a vector or scalar value),
<!-- , and when a vector of
quantities is given)
-->
and optionally it can specify an offset, which in itself has to be an
epoch. Allowable reference codes are:
{\em UTC TAI LAST LMST GMST1 GAST UT1 UT2 TDT TCG TDB TCB}.\\
Note that additional ones may become available. Check in \casa\ with:

\begin{verbatim}
"""
#
print("\t----\t epoch Ex 1 \t----")
print(me.listcodes(me.epoch()))
#{'normal': ['LAST', 'LMST', 'GMST1', 'GAST', 'UT1', 'UT2', 'UTC', 'TAI',
# 'TDT', 'TCG', 'TDB', 'TCB', 'IAT', 'GMST', 'TT', 'ET', 'UT'], 'extra': []}
#
"""
\end{verbatim}
See <link anchor="quanta">quantity</link> for possible time formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t epoch Ex 2 \t----")
print(me.epoch('utc','today'))
#{'m0': {'value': 54048.861237743055, 'unit': 'd'},
# 'refer': 'UTC',
# 'type': 'epoch'}
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="direction">
   <shortdescription>define a direction measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>J2000</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>longitude</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v1">
     <any type="variant"/>
     <description>latitude</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset direction measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<any type="record"/>
<shortdescription>direction measure</shortdescription>
</returns>
<description>
direction defines a direction measure from the CLI. It has to specify a
reference code, direction quantity values (see introduction for the action on a
scalar quantity with either a vector or scalar value),
<!-- , and when a vector of
quantities is given), -->
and optionally it can specify an
offset, which in itself has to be a direction. Allowable reference codes are:
{\em J2000 JMEAN JTRUE APP B1950 BMEAN BTRUE GALACTIC HADEC AZEL
SUPERGAL ECLIPTIC MECLIPTIC TECLIPTIC MERCURY
VENUS MARS JUPITER SATURN URANUS NEPTUNE PLUTO MOON SUN COMET}.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t direction Ex 1 \t----")
print(me.listcodes(me.direction()))
#{'normal': ['J2000', 'JMEAN', 'JTRUE', 'APP', 'B1950', 'BMEAN',
#'BTRUE', 'GALACTIC', 'HADEC', 'AZEL', 'AZELSW', 'AZELNE', 'AZELGEO',
#'AZELSWGEO', 'AZELNEGEO', 'JNAT', 'ECLIPTIC', 'MECLIPTIC',
#'TECLIPTIC', 'SUPERGAL', 'ITRF', 'TOPO', 'ICRS'], 'extra': ['MERCURY',
#'VENUS', 'MARS', 'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE', 'PLUTO',
#'SUN', 'MOON', 'COMET']}
#
"""
\end{verbatim}
The direction quantity values should be longitude(angle) and
latitude(angle) (none needed for planets: the frame epoch defines coordinates).
See <link anchor="quanta">quantity</link> for possible angle formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t direction Ex 2 \t----")
print(me.direction('j2000','30deg','40deg'))
#{'m0': {'value': 0.52359877559829882, 'unit': 'rad'},
# 'm1': {'value': 0.69813170079773168, 'unit': 'rad'},
# 'refer': 'J2000',
# 'type': 'direction'}
#
print(me.direction('mars'))
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': 1.5707963267948966, 'unit': 'rad'},
# 'refer': 'MARS',
# 'type': 'direction'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="getvalue">
   <shortdescription>get the value of a measure</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure (array of measures)</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/>
<shortdescription>array of quantities</shortdescription>
</returns>
<description>
getvalue gets the actual implementation value of the measure.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t getvalue Ex 1 \t----")
b=me.direction('j2000','0deg','80deg')
print(me.getvalue(b))
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': 1.3962634015954634, 'unit': 'rad'}}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="gettype">
   <shortdescription>get the type of a measure</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure (array of measures)</description>
     </param>
</input>
<returns xsi:type="string"/>
<description>
gettype gets the actual type of the measure.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t gettype Ex 1 \t----")
b=me.direction('j2000','0deg','80deg')
print(me.getvalue(b))
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': 1.3962634015954634, 'unit': 'rad'}}
print(me.gettype(b))
#'Direction'
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="getref">
   <shortdescription>get the reference code of a measure</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure (array of measures)</description>
     </param>
</input>
<returns xsi:type="string"/>
<description>
gettype gets the actual reference code of the measure.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t getref Ex 1 \t----")
b=me.direction('j2000','0deg','80deg')
print(me.getvalue(b))
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': 1.3962634015954634, 'unit': 'rad'}}
print(me.gettype(b))
#'Direction'
print(me.getref(b))
#'J2000'
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="getoffset">
   <shortdescription>get the offset of a measure</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure (array of measures)</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>measure or boolean</shortdescription>
<any type="record"/>
</returns>
<description>
getoff gets the actual offset of the measure (as a measure) or F if no offset
given.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t getoffset Ex 1 \t----")
b=me.direction('j2000','0deg','80deg')
print(me.getvalue(b))
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': 1.3962634015954634, 'unit': 'rad'}}
print(me.gettype(b))
#'Direction'
print(me.getref(b))
#'J2000'
print(me.getoffset(b))
#{}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="cometname">
   <shortdescription>get the current comet name</shortdescription>

<input>
</input>
<returns xsi:type="string"/>
<description>
cometname gets the name of the current comet (if any).
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t cometname Ex 1 \t----")
print(me.cometname())
#Thu Nov 9 21:27:25 2006      WARN :
#Method cometname fails! No Comet table present
#''
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="comettype">
   <shortdescription>get the current comet table type</shortdescription>

<input>
</input>
<returns xsi:type="string"/>
<description>
comettype gets the comet table type (apparent or topocentric)
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t comettype Ex 1 \t----")
print(me.comettype())
# 'none'
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="comettopo">
   <shortdescription>get the current comet table coordinates</shortdescription>

<input>
</input>
<returns xsi:type="any"><any type="record"/>
<shortdescription>position measure or fail</shortdescription>
</returns>
<description>
comettopo gets the comet table's topographic coordinates used.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t comettopo Ex 1 \t----")
print(me.comettopo())
#Thu Nov 9 21:45:40 2006      WARN :
#Method comettopo fails!  No Topocentric Comet table present
#{'value': [0.0], 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="framecomet">
   <shortdescription>set the current comet table</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="v">
     <description>name of a table</description>
     <value></value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
framecomet will put the specified comet table in the frame.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t framecomet Ex 1 \t----")
print(me.framecomet('VGEO'))
#True
print(me.showframe())
#'Frame: VENUS comet between MJD 50802.7 and 50803.1'
print(me.cometname())
#'VENUS'
print(me.comettype())
#'APP'
print(me.doframe(me.epoch('et',qa.quantity('1997/12/20/17:30:0'))))
#True
print(me.measure(me.direction('comet'),'app'))
#{'m0': {'value': -0.94936485919663083, 'unit': 'rad'},
# 'm1': {'value': -0.34710256485894436, 'unit': 'rad'},
# 'refer': 'APP',
# 'type': 'direction'}
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="position">
   <shortdescription>define a position measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>WGS84</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>longitude or x</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v1">
     <any type="variant"/>
     <description>latitude or y</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v2">
     <any type="variant"/>
     <description>height or z</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset position measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>position measure</shortdescription>
<any type="record"/>
</returns>
<description>
position defines a position measure from the CLI. It has to specify a
reference code, position quantity values (see introduction for the action on a
scalar quantity with either a vector or scalar value),
<!-- and when a vector of
quantities is given), -->
and optionally it can specify an
offset, which in itself has to be a position. Allowable reference codes are:
{\em WGS84 ITRF} (World Geodetic System and International Terrestrial
Reference Frame).\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t position Ex 1 \t----")
print(me.listcodes(me.position()))
#{'normal': ['ITRF', 'WGS84'], 'extra': []}
#
"""
\end{verbatim}
 The position quantity values should be either longitude
(angle), latitude(angle) and height(length); or x,y,z (length).
See <link anchor="quanta">quantity</link> for possible angle formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t position Ex 2 \t----")
print(me.position('wgs84','30deg','40deg','10m'))
#{'m0': {'value': 0.52359877559829882, 'unit': 'rad'},
# 'm1': {'value': 0.6981317007977319, 'unit': 'rad'},
# 'm2': {'value': 9.9999999999999982, 'unit': 'm'},
# 'refer': 'WGS84',
# 'type': 'position'}
print(me.observatory('ATCA'))
#{'m0': {'value': 2.6101423190348916, 'unit': 'rad'},
# 'm1': {'value': -0.5261379196128062, 'unit': 'rad'},
# 'm2': {'value': 6372960.2577234386, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'position'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="observatory">
   <shortdescription>get position of an observatory</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="name">
     <description>observatory name - case insensitive</description>
     <value>ATCA</value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>position measure</shortdescription>
<any type="record"/>
</returns>
<description>
observatory will give you the position of an observatory as given in the
system. At the time of writing the following observatories are recognised
(but check e.g. the position GUI for currently known ones, or the
me.obslist() tool function):
{\em ALMA ARECIBO ATCA BIMA CLRO DRAO DWL GB GBT GMRT IRAM PDB IRAM_PDB
 JCMT MOPRA MOST NRAO12M NRAO_GBT PKS SAO SMA VLA VLBA WSRT}.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t observatory Ex 1 \t----")
print(me.observatory('ATCA'))
#{'m0': {'value': 2.6101423190348916, 'unit': 'rad'},
# 'm1': {'value': -0.5261379196128062, 'unit': 'rad'},
# 'm2': {'value': 6372960.2577234386, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'position'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="obslist">
   <shortdescription>get a list of known observatories</shortdescription>

<returns xsi:type="string">
<shortdescription>position measure</shortdescription>
</returns>
<description>
obslist will give you a string with the space separated list of
observatories known in the Observatories table.
</description>
<!-- (see <link
anchor="measuresdata">measuresdata</link> module). -->
<example>
\begin{verbatim}
"""
#
print("\t----\t obslist Ex 1 \t----")
print(me.obslist())
#'ALMA ARECIBO ATCA BIMA CLRO DRAO DWL GB GBT GMRT IRAM
# PDB IRAM_PDB JCMT MOPRA MOST NRAO12M NRAO_GBT PKS SAO
# SMA VLA VLBA WSRT'
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="linelist">
   <shortdescription>get a list of known spectral lines</shortdescription>

<returns xsi:type="string"/>
<description>
linelist will give you a string with a space separated list of spectral lines
known in the Lines table.
<!-- (see <link anchor="measuresdata">measuresdata</link>
module).-->
A number of lines are available now, but tables with many lines are
already online, and will be interfaced once a nomenclature can be defined for
the tens of thousands of lines.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t linelist Ex 1 \t----")
print(me.linelist())
#'C109A CI CII166A DI H107A H110A H138B H166A H240A H272A
# H2CO HE110A HE138B HI OH1612 OH1665 OH1667 OH1720'
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="spectralline">
   <shortdescription>get frequency of a spectral line</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="name">
     <description>name</description>
     <value>HI</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
spectralline will give you the frequency of a spectral line. The known list
can be obtained by <link anchor="measures:measures.linelist.function">me.linelist()</link>.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t spectralline Ex 1 \t----")
print(me.spectralline('HI'))
#{'m0': {'value': 1420405751.786, 'unit': 'Hz'},
# 'refer': 'REST',
# 'type': 'frequency'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="sourcelist">
   <shortdescription>get a list of known sources</shortdescription>

<returns xsi:type="string"/>
<description>
sourcelist will give you a string with the space separated list of sources
known in the Sources table.
</description>
<!-- (see <link anchor="measuresdata">measuresdata</link> module). -->
<example>
\begin{verbatim}
"""
#
print("\t----\t sourcelist Ex 1 \t----")
print(me.sourcelist()[0:62])
#'0002-478 0003+380 0003-066 0007+106 0007+171 0008-264 0008-421'
#......
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="source">
   <shortdescription>get direction of a source</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="name">
     <description>name</description>
     <any type="variant"/>
     <value>1934-638</value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>direction measure</shortdescription>
<any type="record"/>
</returns>
<description>
source will give you the direction of a source. The known list
can be obtained by <link anchor="measures:measures.sourcelist.function">me.sourcelist()</link>.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t source Ex 1 \t----")
print(me.source())
print(me.source('1934-638'))
#  Out[19]:
#{'m0': {'value': -1.1370073467795063, 'unit': 'rad'},
# 'm1': {'value': -1.1119959323803881, 'unit': 'rad'},
# 'refer': 'ICRS',
# 'type': 'direction'}
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="frequency">
   <shortdescription>define a frequency measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>LSRK</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>frequency/wavelength/\ldots</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset frequency measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>frequency measure</shortdescription>
<any type="record"/>
</returns>
<description>
frequency defines a frequency measure from the CLI. It has to specify a
reference code, frequency quantity value (see introduction for the action on a
scalar quantity with either a vector or scalar value),
<!-- and when a vector of
quantities is given), -->
and optionally it can specify an
offset, which in itself has to be a frequency. Allowable reference codes are:
{\em REST LSRK LSRD BARY GEO TOPO GALACTO LGROUP CMB}.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t frequency Ex 1 \t----")
print(me.listcodes(me.frequency()))
#{'normal': ['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO',
# 'GALACTO', 'LGROUP', 'CMB'], 'extra': []}
#
"""
\end{verbatim}
The frequency quantity values should be in one of the recognised units
(examples all give same frequency):
\begin{itemize}
\item value with time units: a period (0.5s)
\item value as frequency: 2Hz
\item value in angular frequency: 720deg/s
\item value as length: 149896km
\item value as wave number: 4.19169e-8m-1
\item value as energy (h.nu): 8.27134e-9ueV
\item value as momentum: 4.42044e-42kg.m
\end{itemize}
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t frequency Ex 2 \t----")
print(me.frequency('lsrk','5GHz'))
#{'m0': {'value': 5000000000.0, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
print(me.frequency('lsrk','21cm'))
#{'m0': {'value': 1427583133.3333333, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="doppler">
   <shortdescription>define a doppler measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>RADIO</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <value></value>
     <description>doppler ratio/velocity</description>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset doppler measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>doppler measure</shortdescription>
<any type="record"/>
</returns>
<description>
doppler defines a doppler measure from the CLI. It has to specify a
reference code, doppler quantity value (see introduction for the action on a
scalar quantity with either a vector or scalar value),
<!-- and when a vector of
quantities is given),-->
and optionally it can specify an offset,
which in itself has to be a doppler. Allowable reference codes are:
{\em RADIO Z RATIO BETA GAMMA OPTICAL TRUE RELATIVISTIC}.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t doppler Ex 1 \t----")
print(me.listcodes(me.doppler()))
#{'normal': ['RADIO', 'Z', 'RATIO', 'BETA', 'GAMMA', 'OPTICAL',
# 'TRUE', 'RELATIVISTIC'], 'extra': []}
#
"""
\end{verbatim}
The doppler quantity values should be either non-dimensioned to specify a
ratio of the light velocity, or in velocity.
</description>
<example>
Examples both give same doppler:
\begin{verbatim}
"""
#
print("\t----\t doppler Ex 2 \t----")
print(me.doppler('radio','0.4'))
#{'m0': {'value': 119916983.2, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
print(me.doppler('radio',qa.mul(qa.quantity('0.4'),qa.constants('c'))))
#{'m0': {'value': 119916983.2, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="radialvelocity">
   <shortdescription>define a radialvelocity measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>LSRK</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>radial velocity</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset radialvelocity measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>radialvelocity measure</shortdescription>
<any type="record"/>
</returns>
<description>
radialvelocity defines a radialvelocity measure from the CLI. It has to
specify a reference code, radialvelocity quantity value (see introduction for
the action on a
scalar quantity with either a vector or scalar value),
<!-- and when a vector of
quantities is given), -->
and optionally it
can specify an offset, which in itself has to be a radialvelocity.
Allowable reference codes are:
{\em LSRK LSRD BARY GEO TOPO GALACTO LGROUP CMB}.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t radialvelocity Ex 1 \t----")
print(me.listcodes(me.radialvelocity()))
#  Out[17]:
#{'extra': [],
# 'normal': ['LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO',
# 'LGROUP', 'CMB']}
#
"""
\end{verbatim}
The radialvelocity quantity values should be given as velocity.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t radialvelocity Ex 2 \t----")
print(me.radialvelocity('lsrk','20km/s'))
#  Out[18]:
#{'m0': {'value': 20000.0, 'unit': 'm/s'},
# 'refer': 'LSRK',
# 'type': 'radialvelocity'}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="uvw">
   <shortdescription>define a uvw measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>ITRF</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>longitude or x</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v1">
     <any type="variant"/>
     <description>latitude or y</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v2">
     <any type="variant"/>
     <description>height or z</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record" />
     <description>optional offset uvw measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>uvw measure</shortdescription>
<any type="record" />
</returns>
<description>
uvw defines a uvw measure from the CLI. It has to specify a
reference code, uvw quantity values (see introduction for the action on a
scalar quantity with either a vector or scalar value), <!-- and when
a vector of quantities is given), --> and optionally it can specify an
offset, which in itself has to be a uvw. Allowable reference codes are
ITRF and the direction ones.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t uvw Ex 1 \t----")
print(me.listcodes(me.uvw()))
#{'normal': ['J2000', 'JMEAN', 'JTRUE', 'APP', 'B1950', 'BMEAN',
# 'BTRUE', 'GALACTIC', 'HADEC', 'AZEL', 'AZELSW', 'AZELNE',
# 'AZELGEO', 'AZELSWGEO', 'AZELNEGEO', 'JNAT', 'ECLIPTIC',
# 'MECLIPTIC', 'TECLIPTIC', 'SUPERGAL', 'ITRF', 'TOPO',
# 'ICRS'], 'extra': []}
#
"""
\end{verbatim}
The uvw quantity values should be either longitude
(angle), latitude(angle) and height(length); or x,y,z (length).
See <link anchor="quanta">quantity</link> for possible angle formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t uvw Ex 2 \t----")
print(me.uvw('itrf','30deg','40deg','10m'))
#{'m0': {'value': 0.52359877559829882, 'unit': 'rad'},
# 'm1': {'value': 0.6981317007977319, 'unit': 'rad'},
# 'm2': {'value': 9.9999999999999982, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'uvw'}
print(me.doframe(me.epoch('utc','today')))
#True
print(me.doframe(me.observatory('ALMA')))
#True
print(me.doframe(me.direction('mars')))
#True
print(me.measure(me.uvw('itrf','30deg','40deg','10m'), 'j2000'))
#{'m0': {'value': 0.52321924738347259, 'unit': 'rad'},
# 'm1': {'value': 0.69813169995801672, 'unit': 'rad'},
# 'm2': {'value': 10.0, 'unit': 'm'},
# 'refer': 'J2000',
# 'type': 'uvw'}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="touvw">
   <shortdescription>calculate a uvw measure from a baseline</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>baseline measure</description>
     </param>

</input>
<output>
     <param xsi:type="any" direction="out" name="dot">
     <any type="record"/>
     <description>uvw-dot (quantity array)</description>
     </param>

     <param xsi:type="any" direction="out" name="xyz">
     <any type="record"/>
     <description>uvw (quantity array)</description>
     </param>

</output>
<returns xsi:type="any">
<shortdescription>uvw measure</shortdescription>
<any type="record"/>
</returns>
<description>
touvw calculates a uvw measure from a baseline. <!-- The baseline can
consist of a vector of actual baseline positions. -->  Note that the
baseline does not have to be a proper {\em baseline}, but can be a
series of positions (to call positions baselines see
<link anchor="measures:measures.asbaseline.function">asbaseline</link> ) for speed reasons:
operations are linear and can be done on positions, which are
converted to baseline values at the end (with
<link anchor="measures:measures.expand.function">expand</link> ).

Whatever the reference code of the baseline, the returned {\em uvw} will be
given in J2000. If the {\em dot} argument is given, that variable
will be filled with a quantity array consisting of the time
derivative of the uvw (note that only the sidereal rate is taken
into account; not precession, earth tides and similar variations,
which are much smaller). If the {\em xyz} variable is given, it will
be filled with the quantity values of the uvw measure.

The values of the input baselines can be given as a quantity
vector per x, y or z value. <!--; or as an array as explained in the introduction.-->

uvw coordinates are calculated for a certain direction in the sky;
hence the frame has to contain the direction for the calculation to
work. Since the baseline and the sky rotate with respect of each
other, the time should be specified as well.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t touvw Ex 1 \t----")
print(me.doframe(me.observatory('atca')))
#True
print(me.doframe(me.source('1934-638')))
#True
print(me.doframe(me.epoch('utc',qa.unit('today'))))
#True
b=me.baseline('itrf','10m','20m','30m')
print(me.touvw(b))
#{'dot': {'unit': 'm/s',
#         'value': [-0.0011912452908351659,
#                   -0.00098731747136827593,
#                   -0.00048769097314181744]},
# 'return': {'m0': {'value': -0.094777304811312649, 'unit': 'rad'},
#            'm1': {'value': -1.1509286139398101, 'unit': 'rad'},
#            'm2': {'value': 37.416573867739416, 'unit': 'm'},
#            'refer': 'J2000',
#            'type': 'uvw'},
# 'xyz': {'unit': 'm',
#         'value': [15.184026188402472,
#                   -1.4434256399579168,
#                   -34.166677788919138]}}
print(me.getvalue(me.touvw(b)))
#{'m0': {'value': -0.094777304811312649, 'unit': 'rad'},
# 'm1': {'value': -1.1509286139398101, 'unit': 'rad'},
# 'm2': {'value': 37.416573867739416, 'unit': 'm'}}
print(me.getvalue(me.touvw(b))['m0'])
#{'value': -0.094777304811312649, 'unit': 'rad'}
#
"""
\end{verbatim}
<!--
An example with more than one value:
\begin{verbatim}
#
sb = me.baseline('itrf',qa.unit([10,50],'m'),qa.unit([20,100],'m'),
                qa.unit([30,150],'m'))
print me.touvw(sb,d,x); print d; print x;
#[type=uvw, refer=J2000, m2=[value=[37.4165739 187.082869] , unit=m],
#           m1=[unit=rad, value=[-0.743811234 -0.743811234] ],
#           m0=[unit=rad, value=[2.50094148 2.50094148] ]]
#[value=[[1:3,]
#    0.00025643414  0.0012821707
#    0.00143537137  0.00717685684
#    0.000709009712 0.00354504856], unit=m/s]
#[value=[[1:3,]
#    -22.0746793 -110.373397
#    16.45792    82.2895998
#    -25.334668  -126.67334], unit=m]
print me.getvalue(me.touvw(sb))
#[*7=[unit=rad, value=[2.50094148 2.50094148] ],
#               *8=[unit=rad, value=[-0.743811234 -0.743811234] ],
#                *9=[value=[37.4165739 187.082869] , unit=m]]
print me.getvalue(me.touvw(sb))[1]
#[unit=rad, value=[2.50094148 2.50094148] ]
print qa.getvalue(me.getvalue(me.touvw(sb))[1])[2]
#2.50094148
print me.doframe(me.epoch('utc','today'))
#T
print me.expand(sb)
#[type=baseline, refer=ITRF, m2=[value=149.666295, unit=m],
#                m1=[unit=rad, value=0.930274014],
#                m0=[unit=rad, value=1.10714872]]
print me.expand(sb,x)
#[type=baseline, refer=ITRF, m2=[value=149.666295, unit=m],
#               m1=[unit=rad, value=0.930274014],
#               m0=[unit=rad, value=1.10714872]]
print x
#[value=[[1:3,]
#    40
#    80
#    120], unit=m]
print me.expand(me.touvw(sb),x); x
#[type=uvw, refer=J2000, m2=[value=149.666295, unit=m],
#          m1=[unit=rad, value=-0.654614537],
#          m0=[unit=rad, value=2.32532487]]
#[value=[[1:3,]
#    -81.3219596
#     86.5043397
#    -91.124849], unit=m]
print me.touvw(me.expand(sb),xyz=x); x
#[type=uvw, refer=J2000, m2=[value=149.666295, unit=m],
#          m1=[unit=rad, value=-0.654614537],
#          m0=[unit=rad, value=2.32532487]]
#[value=[[1:3,]
#    -81.3219596
#    86.5043397
#    -91.124849], unit=m]
a = me.touvw(sb, xyz=x)#; x
print a
#[value=[[1:3,]
#    -20.3304899 -101.652449
#    21.6260849  108.130425
#    -22.7812122 -113.906061], unit=m]
#
\end{verbatim}
-->
</example>
</method>




   <method type="function" name="expand">
   <shortdescription>expand n positions to n*(n-1)/2 baselines</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure (baseline, position or uvw measure)</description>
     </param>
</input>
<output>
     <param xsi:type="any" direction="out" name="xyz">
     <any type="record"/>
     <description>uvw (quantity array)</description>
     </param>
</output>

<returns xsi:type="any"><shortdescription>measure</shortdescription>
<any type="record"/>
</returns>
<description>
expand calculates the differences between a series of given measure
values: it calculates baseline values from position values. The
returned value is a measure, but the value of the optional output
variable {\em xyz} will be set to an array of values.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t expand Ex 1 \t----")
b=me.baseline('itrf','10m','20m','30m')
print(me.expand(b))
#{'return': {'m0': {'value': 1.1071487177940904, 'unit': 'rad'},
#            'm1': {'value': 0.93027401411547195, 'unit': 'rad'},
#            'm2': {'value': 37.416573867739416, 'unit': 'm'},
#            'refer': 'J2000',
#            'type': 'uvw'},
# 'xyz': {'value': [10.000000000000004, 20.000000000000004, 30.0], 'unit': 'm'}}
#
"""
\end{verbatim}
<!--
sb = me.baseline('itrf',qa.unit([10,50],'m'),qa.unit([20,100],'m'),
+    qa.unit([30,150],'m'))
print me.expand(sb,x); x
[type=baseline, refer=ITRF, m2=[value=149.666295, unit=m],
                m1=[unit=rad, value=0.930274014],
                m0=[unit=rad, value=1.10714872]]
[value=[[1:3,]
    40
    80
    120], unit=m]
-
-->
</example>
</method>




   <method type="function" name="earthmagnetic">
   <shortdescription>define an earthmagnetic measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>IGRF</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <value></value>
     <description>Field strength</description>
     </param>

     <param xsi:type="any" direction="in" name="v1">
     <any type="variant"/>
     <description>longitude</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v2">
     <any type="variant"/>
     <description>latitude</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset earthmagnetic measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>earthmagnetic measure</shortdescription>
<any type="record"/>
</returns>
<description>
earthmagnetic defines an earthmagnetic measure from the CLI. It needs
a reference code, earthmagnetic quantity values
(see introduction for the action on a
scalar quantity with either a vector or scalar value) <!-- and when a vector
of quantities is given) --> if the reference code is not
for a model, and optionally it
can specify an offset, which in itself has to be a earthmagnetic. In general
you specify a model (IGRF is the default and the only one known) and convert
it to an explicit field.  (See
\begin{verbatim}
 http://fdd.gsfc.nasa.gov/IGRF.html
\end{verbatim}
for information on the International Geomagnetic Reference Field). The
earthmagnetic quantity values should be either longitude (angle),
latitude(angle) and length(field strength); or x,y,z (field).
See <link anchor="quanta">quantity</link> for possible angle formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t earthmagnetic Ex 1 \t----")
print(me.earthmagnetic('igrf'))
#{'type': 'earthmagnetic', 'refer': 'IGRF', 'm1': {'value': 0.0, 'unit': 'nT'},
# 'm0': {'value': 6.1230317691118855e-23, 'unit': 'nT'},
# 'm2': {'value': 9.9999999999999995e-07, 'unit': 'nT'}}
print(me.doframe(me.observatory('atca')))
print(me.doframe(me.source('1934-638')))
print(me.doframe(me.epoch('utc',qa.unit('today'))))
print(me.measure(me.earthmagnetic('igrf'), 'j2000'))
#{'type': 'earthmagnetic', 'refer': 'J2000',
# 'm1': {'value': -8664.8767628222304, 'unit': 'nT'},
# 'm0': {'value': 50544.054410564473, 'unit': 'nT'},
# 'm2': {'value': 1799.5131920958615, 'unit': 'nT'}}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="baseline">
   <shortdescription>define a baseline measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>reference code</description>
     <value>ITRF</value>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="variant"/>
     <description>longitude or x</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v1">
     <any type="variant"/>
     <description>latitude or y</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="v2">
     <any type="variant"/>
     <description>height or z</description>
     <value></value>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional offset baseline measure</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>baseline measure</shortdescription>
<any type="record"/>
</returns>
<description>
baseline defines a baseline measure from the CLI. It has to specify a
reference code, baseline quantity values (see introduction for the action on a
scalar quantity with either a vector or scalar value, and when a vector of
quantities is given), and optionally it can specify an
offset, which in itself has to be a baseline. Allowable reference codes are
ITRF and the direction ones.\\
Note that additional ones may become available. Check in \casa\ with:
\begin{verbatim}
"""
#
print("\t----\t baseline Ex 1 \t----")
print(me.listcodes(me.baseline()))
#{'normal': ['J2000', 'JMEAN', 'JTRUE', 'APP', 'B1950', 'BMEAN', 'BTRUE',
# 'GALACTIC', 'HADEC', 'AZEL', 'AZELSW', 'AZELNE', 'AZELGEO', 'AZELSWGEO',
# 'AZELNEGEO', 'JNAT', 'ECLIPTIC', 'MECLIPTIC', 'TECLIPTIC', 'SUPERGAL',
# 'ITRF', 'TOPO', 'ICRS'], 'extra': []}
#
"""
\end{verbatim}
The baseline quantity values should be either longitude
(angle), latitude(angle) and height(length); or x,y,z (length).
See <link anchor="quanta">quantity</link> for possible angle formats.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t Ex 2 \t----")
print(me.baseline('itrf','30deg','40deg','10m'))
#{'m0': {'value': 0.52359877559829882, 'unit': 'rad'},
# 'm1': {'value': 0.6981317007977319, 'unit': 'rad'},
# 'm2': {'value': 9.9999999999999982, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'baseline'}
print(me.doframe(me.observatory('atca')))
print(me.doframe(me.source('1934-638')))
print(me.doframe(me.epoch('utc',qa.unit('today'))))
print(me.measure(me.baseline('itrf','30deg','40deg','10m'), 'J2000'))
#{'m0': {'value': 0.58375325605991979, 'unit': 'rad'},
# 'm1': {'value': 0.69758519780286155, 'unit': 'rad'},
# 'm2': {'value': 9.9999999999999964, 'unit': 'm'},
# 'refer': 'J2000',
# 'type': 'baseline'}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="asbaseline">
   <shortdescription>define a baseline from a position measure</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="pos">
     <any type="record"/>
     <description>position measure</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>baseline measure</shortdescription>
<any type="record"/>
</returns>
<description>
asbaseline converts a position measure into a baseline measure. No
actual baseline is calculated, since operations can be done on
positions, with subtractions to obtain baselines at a later stage.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t asbaseline Ex 1 \t----")
b = me.position('itrf', '10m', '20m', '30m');
print(b)
#{'m0': {'value': 1.1071487177940904, 'unit': 'rad'},
# 'm1': {'value': 0.93027401411547195, 'unit': 'rad'},
# 'm2': {'value': 37.416573867739416, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'position'}
print(me.asbaseline(b))
#{'m0': {'value': 1.1071487177940904, 'unit': 'rad'},
# 'm1': {'value': 0.93027401411547195, 'unit': 'rad'},
# 'm2': {'value': 37.416573867739416, 'unit': 'm'},
# 'refer': 'J2000',
# 'type': 'baseline'}
#
"""
\end{verbatim}
<!--
sb = me.position('itrf',qa.unit([10,50],'m'),qa.unit([20,100],'m'),
+ qa.unit([30,150],'m'));
print b; print sb;
[type=position, refer=ITRF, m2=[value=37.4165739, unit=m],
                 m1=[unit=rad, value=0.930274014],
                 m0=[unit=rad, value=1.10714872]]
[type=position, refer=ITRF, m2=[value=[37.4165739 187.082869] , unit=m],
         m1=[unit=rad, value=[0.930274014 0.930274014] ],
         m0=[unit=rad, value=[1.10714872 1.10714872] ]]
print b; print me.asbaseline(b); print sb; print me.asbaseline(sb)
[type=position, refer=ITRF, m2=[value=37.4165739, unit=m],
                m1=[unit=rad, value=0.930274014],
                m0=[unit=rad, value=1.10714872]]
[type=baseline, refer=ITRF, m2=[value=37.4165739, unit=m],
                m1=[unit=rad, value=0.930274014],
                m0=[unit=rad, value=1.10714872]]
[type=position, refer=ITRF, m2=[value=[37.4165739 187.082869] , unit=m],
         m1=[unit=rad, value=[0.930274014 0.930274014] ],
         m0=[unit=rad, value=[1.10714872 1.10714872] ]]
[type=baseline, refer=ITRF, m2=[value=[37.4165739 187.082869] , unit=m],
         m1=[unit=rad, value=[0.930274014 0.930274014] ],
         m0=[unit=rad, value=[1.10714872 1.10714872] ]]
-->
</example>
</method>




   <method type="function" name="listcodes">
   <shortdescription>get known reference codes</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="ms">
     <any type="record"/>
     <description>the measure type for which to list</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>record with two string vectors</shortdescription>
<any type="record"/>
</returns>
<description>
listcodes will produce the known reference codes for a specified measure
type. It will return a record with two entries. The first is a string vector
of all normal codes; the second a string vector (maybe empty) with all extra
codes (like planets).
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t listcodes Ex 1 \t----")
# Generate some direction
# Note that an empty or non-specified reference code will produce the
# measure with the default code for that measure type
a=me.direction()
print(me.getref(a))
#'J2000'
print(me.ismeasure(a))
#True
# Get the known reference codes for direction
print(me.listcodes(a))
#{'normal': ['J2000', 'JMEAN', 'JTRUE', 'APP', 'B1950', 'BMEAN',
#'BTRUE', 'GALACTIC', 'HADEC', 'AZEL', 'AZELSW', 'AZELNE', 'AZELGEO',
#'AZELSWGEO', 'AZELNEGEO', 'JNAT', 'ECLIPTIC', 'MECLIPTIC',
#'TECLIPTIC', 'SUPERGAL', 'ITRF', 'TOPO', 'ICRS'],
#'extra': ['MERCURY', 'VENUS', 'MARS', 'JUPITER', 'SATURN', 'URANUS',
#'NEPTUNE', 'PLUTO', 'SUN', 'MOON', 'COMET']}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="measure">
   <shortdescription>convert a measure to another reference</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>measure to be converted</description>
     </param>

     <param xsi:type="string" direction="in" name="rf">
     <description>output reference code</description>
     </param>

     <param xsi:type="any" direction="in" name="off">
     <any type="record"/>
     <description>optional output offset measure</description>
     <value></value>
     </param>

</input>
<returns xsi:type="any">
<shortdescription>measure</shortdescription>
<any type="record"/>
</returns>
<description>
measure converts measures (epoch, direction etc.) from one reference to
another. It will, for instance, convert a direction from J2000 to AZEL
representation. \\
Its arguments are a measure, an output reference code (see the individual
measures for the allowable codes (<link anchor="measures:measures.direction.function">direction</link>,
<link anchor="measures:measures.position.function">position</link>,
<link anchor="measures:measures.epoch.function">epoch</link>,
<link anchor="measures:measures.frequency.function">frequency</link>,
<link anchor="measures:measures.doppler.function">doppler</link>,
<link anchor="measures:measures.radialvelocity.function">radialvelocity</link>,
<link anchor="measures:measures.baseline.function">baseline</link>,
<link anchor="measures:measures.uvw.function">uvw</link>,
<link anchor="measures:measures.earthmagnetic.function">earthmagnetic</link>)), and an optional offset of
the same type as the main measure. The offset will be subtracted from the
result before it is returned.\\
In some cases (see the individual measures for when), more information than
just a reference code is necessary. E.g. the above example of a conversion to
AZEL, needs to know for when, and where on Earth we want it. This information
is stored in a reference frame. Measures are set in the reference frame with
the <link anchor="measures:measures.doframe.function">doframe</link> function. The frame is tool
 wide.\\
<!-- An optional argument {em qv} can contain a quantity with a vector of values
and a unit. When present, the vector will be converted (and returned) in the
same way as the main reference. Not all Measures accept this argument yet. If
needed please ask an enhancement. -->
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t measure Ex 1 \t----")
a = me.epoch('utc','today')                             # a time
print(a)
#{'m0': {'value': 54054.872957673608, 'unit': 'd'},
# 'refer': 'UTC',
# 'type': 'epoch'}
print(me.doframe(me.source('1934-638')))
print(me.measure(a, 'tai'))                                     # convert to IAT
#{'m0': {'value': 54054.873339618054, 'unit': 'd'},
# 'refer': 'TAI',
# 'type': 'epoch'}
print(me.doframe(a))                                            # set time in frame
#True
print(me.doframe(me.observatory('ALMA')))                       # set position in frame
#True
b=me.direction('j2000', qa.toangle('0h'), '-30deg')  # a direction
print(b)
#{'m0': {'value': 0.0, 'unit': 'rad'},
# 'm1': {'value': -0.52359877559829882, 'unit': 'rad'},
# 'refer': 'J2000',
# 'type': 'direction'}
print(me.measure(b, 'azel'))                                    # convert to AZEL
#{'m0': {'value': 1.9244096810822324, 'unit': 'rad'},
# 'm1': {'value': 0.76465385681363052, 'unit': 'rad'},
# 'refer': 'AZEL',
# 'type': 'direction'}
print(qa.angle(me.getvalue(me.measure(b,'azel'))['m0']))     # show as angles
#['+110.15.38']
print(qa.angle(me.getvalue(me.measure(b,'azel'))['m1']))
#['+043.48.41']
#
"""
\end{verbatim}
<!-- In the following the {\em qv} argument is used.-->
Another example:
\begin{verbatim}
"""
#
print("\t----\t measure Ex 2 \t----")
# Fill the frame with necessary information
print(me.doframe(me.epoch('utc','today')))
#True
print(me.doframe(me.observatory('ALMA')))
#True
print(me.doframe(me.direction('mars')))
#True
a=qa.unit('1GHz')
print(a)
#{'value': 1.0, 'unit': 'GHz'}
m=me.frequency('lsrk',qa.quantity(qa.getvalue(a),qa.getunit(a)))
print(m)
#{'m0': {'value': 1000000000.0, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
print(me.measure(m,'lsrd'))
#{'m0': {'value': 1000001766.3928765, 'unit': 'Hz'},
# 'refer': 'LSRD',
# 'type': 'frequency'}
#
"""
\end{verbatim}
<!--
# Make a list of frequencies to be converted
a=qa.unit([1,1.1,1.2,1.3],'GHz')
print a
[value=[1 1.1 1.2 1.3] , unit=GHz]
# Make a frequency measure. Although any can be used, it is advisable
# to use an element of the list, to make sure all units are correct
print m=me.frequency('lsrk',qa.unit(qa.getvalue(a)[1], qa.getunit(a)))
# Convert all
print me.measure(m,'lsrd',qv=a)
[type=frequency, refer=LSRD, m0=[value=999995118, unit=Hz]]
# And check
print a
[value=[0.999995118 1.09999463 1.19999414 1.29999365] , unit=GHz]
-->
</example>
</method>

   <method type="function" name="doframe">
   <shortdescription>save a measure as frame reference</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v">
     <any type="record" />
     <description>measure to be set in frame</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
doframe will set the measure specified as part of a frame. <!-- The
 {\em Frame it} button on the different measure GUI's has the same effect. -->

If conversion from one type to another is necessary, <!-- (either in the GUI
with the {\em Convert$->$} button or --> with the
<link anchor="measures:measures.measure.function">measure</link> function,
<!-- ), --> the following frames
should be set if one of the reference types involved in the conversion is as
in the following lists.\\
{\em Epoch}
\begin{verbatim}
 UTC
 TAI
 LAST   position
 LMST   position
 GMST1
 GAST
 UT1
 UT2
 TDT
 TCG
 TDB
 TCD
\end{verbatim}
{\em Direction}
\begin{verbatim}
 J2000
 JMEAN          epoch
 JTRUE          epoch
 APP            epoch
 B1950
 BMEAN          epoch
 BTRUE          epoch
 GALACTIC
 HADEC          epoch   position
 AZEL           epoch   position
 SUPERGALACTIC
 ECLIPTIC
 MECLIPTIC      epoch
 TECLIPTIC      epoch
 PLANET         epoch   [position]
\end{verbatim}
{\em Position}
\begin{verbatim}
 WGS84
 ITRF
\end{verbatim}
{\em Radial Velocity}
\begin{verbatim}
 LSRK           direction
 LSRD           direction
 BARY           direction
 GEO            direction       epoch
 TOPO           direction       epoch   position
 GALACTO        direction
\end{verbatim}
{\em Doppler}
\begin{verbatim}
 RADIO
 OPTICAL
 Z
 RATIO
 RELATIVISTIC
 BETA
 GAMMA
\end{verbatim}
{\em Frequency}
\begin{verbatim}
 REST           direction                       radialvelocity
 LSRK           direction
 LSRD           direction
 BARY           direction
 GEO            direction       epoch
 TOPO           direction       epoch   position
 GALACTO
\end{verbatim}
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t doframe Ex 1 \t----")
a = me.epoch('utc', 'today')                    # a time
print(a)
#{'m0': {'value': 54054.91671484954, 'unit': 'd'},
# 'refer': 'UTC',
# 'type': 'epoch'}
print(me.doframe(a))                                    # set time in frame
#True
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="framenow">
   <shortdescription>set the active frame time at now</shortdescription>

<input>
</input>
<returns xsi:type="bool"/>
<description>
framenow will fill the active frame time with the current date and time.
The different frame values necessary are described in the
<link anchor="measures:measures.doframe.function">doframe</link> function
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t framenow Ex 1 \t----")
print(me.framenow())                    # specify now as frame reference
#True
print(me.showframe())           # and show the current frame
#'Frame: Epoch: 54054::22:01:42.2880'
#
"""
\end{verbatim}
<!-- epoch: 2002/10/18/00:15:49 UTC
     T
-->
</example>
</method>

   <method type="function" name="showframe">
   <shortdescription>show the currently active frame reference</shortdescription>

<input>
</input>
<returns xsi:type="string"/>
<description>
showframe will display the currently active reference frame values <!-- in
either a separate screen (if working in GUI environment), or --> on the
terminal. <!-- if no GUI active, or a forcing argument {\em F} is given.--> The
different frame values necessary are described in the
<link anchor="measures:measures.doframe.function">doframe</link> function.
The frame is
displayed on the terminal using the formatting as done for the
<link anchor="measures:measures.show.function">show</link> function.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t showframe Ex 1 \t----")
print(me.doframe(me.epoch('utc','today')))              # specify now as frame reference
#T
print(me.showframe())                           # and show the current frame
#'Frame: Epoch: 54054::22:01:42.2880'
#
"""
\end{verbatim}
<!-- epoch: 2002/10/18/00:15:49 UTC
T -->
</example>
</method>


   <method type="function" name="toradialvelocity">
   <shortdescription>convert a doppler type value to a real
radial velocity</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>radial velocity reference type</description>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="record"/>
     <description>doppler value measure</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>radialvelocity measure</shortdescription>
<any type="record"/>
</returns>
<description>
toradialvelocity will convert a Doppler type value (e.g. in radio mode) to a
real radialvelocity. The type of velocity (e.g. LSRK) should be specified
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t toradialvelocity Ex 1 \t----")
a = me.doppler('radio','0.4')
print(a)
#  Out[4]:
#{'m0': {'value': 119916983.2, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
print(me.toradialvelocity('topo',a))
#{'m0': {'value': 141078803.7647059, 'unit': 'm/s'},
# 'refer': 'TOPO',
# 'type': 'radialvelocity'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="tofrequency">
   <shortdescription>convert a doppler type value to a frequency</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>frequency reference type</description>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="record"/>
     <description>doppler measure value</description>
     </param>

     <param xsi:type="any" direction="in" name="rfq">
     <any type="record"/>
     <description>rest frequency (frequency measure or freuency quantity)</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>frequency measure</shortdescription>
<any type="record"/>
</returns>
<description>
tofrequency will convert a Doppler type value (e.g. in radio mode) to a
frequency. The type of frequency (e.g. LSRK) and a rest frequency (either as a
frequency quantity (e.g. qa.constants('HI')) or a frequency measure (e.g.
me.frequency('rest','5100MHz')) should be specified
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t tofrequency Ex 1 \t----")
a=me.doppler('radio','0.4')
print(a)
#{'m0': {'value': 119916983.2, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
print(me.tofrequency('lsrk',a,qa.constants('HI')))
#{'m0': {'value': 852243451.07159996, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="todoppler">
   <shortdescription>convert a frequency or radialvelocity measure
to a doppler measure</shortdescription>

<input>

     <param xsi:type="string" direction="in" name="rf">
     <description>doppler reference type</description>
     </param>

     <param xsi:type="any" direction="in" name="v0">
     <any type="record"/>
     <description>radial velocity or frequency measure</description>
     </param>

     <param xsi:type="any" direction="in" name="rfq">
     <any type="record"/>
     <description>rest frequency (frequency measure or frequency quantity)</description>
     <value></value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>frequency measure</shortdescription>
<any type="record"/>
</returns>
<description>
todoppler will convert a radialvelocity measure or a frequency measure to a
doppler measure. In the case of a frequency, a rest frequency has to be
specified. The type of doppler wanted (e.g. RADIO) has to be specified.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t todoppler Ex 1 \t----")
f = me.frequency('lsrk','1410MHz')     # specify a frequency
print(f)
#{'m0': {'value': 1410000000.0, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
print(me.todoppler('radio', f, qa.constants('HI'))) # give doppler, using HI rest
#{'m0': {'value': 2196249.8401180855, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="torestfrequency">
   <shortdescription>convert a frequency and doppler measure
to a rest frequency</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="v0">
     <any type="record"/>
     <description>frequency reference type</description>
     </param>

     <param xsi:type="any" direction="in" name="d0">
     <any type="record"/>
     <description>doppler measure value</description>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>frequency measure</shortdescription>
<any type="record"/>
</returns>
<description>
torestfrequency will convert a frequency measure and a doppler measure
(e.g. obtained from another spectral line with a known rest frequency) to a
rest frequency.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t torestfrequency Ex 1 \t----")
dp = me.doppler('radio', '2196.24984km/s')  # a measured doppler speed
print(dp)
#{'m0': {'value': 2196249.8399999999, 'unit': 'm/s'},
# 'refer': 'RADIO',
# 'type': 'doppler'}
f = me.frequency('lsrk','1410MHz')    # a measured frequency
print(f)
#{'m0': {'value': 1410000000.0, 'unit': 'Hz'},
# 'refer': 'LSRK',
# 'type': 'frequency'}
print(me.torestfrequency(f, dp))                   # the corresponding rest frequency
#{'m0': {'value': 1420405751.7854364, 'unit': 'Hz'},
# 'refer': 'REST',
# 'type': 'frequency'}
#
"""
\end{verbatim}
</example>
</method>

   <method type="function" name="rise">
   <shortdescription>get rise and set sidereal time</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="crd">
     <any type="variant"/>
     <description>direction of source (direction measure)</description>
     <value></value>
     </param>
     <param xsi:type="any" direction="in" name="ev">
     <description>elevation angle limit</description>
     <any type="variant"/>
     <value>0.0deg</value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>record with rise and set siodereal time quantities; or a 2 element
string with below or above</shortdescription>
<any type="record"/>
</returns>
<description>
rise will give the rise/set hour-angles of a source. It needs the position
in the frame, and a time. If the latter is not set, the current time will be
used.
</description>
<example>
\begin{verbatim}
# NOT IMPLEMENTED
print "\t----\t rise Ex 1 \t----"
print me.rise(me.direction('sun'))
#[rise=[value=267.12445, unit=deg], set=[value=439.029964, unit=deg]]
print qa.form.long(me.rise(me.direction('sun')).rise)
#17:48:29.868
#
\end{verbatim}
</example>
</method>

   <method type="function" name="riseset">
   <shortdescription>get rise and set times </shortdescription>

<input>

     <param xsi:type="any" direction="in" name="crd">
     <any type="variant"/>
     <description>direction of source (direction measure)</description>
     </param>

     <param xsi:type="any" direction="in" name="ev">
     <description>elevation limit</description>
     <any type="variant"/>
     <value>0.0deg</value>
     </param>
</input>
<returns xsi:type="any">
<shortdescription>record with rise and set epoch or strings</shortdescription>
<any type="record"/>
</returns>
<description>
rise will give the rise/set times of a source. It needs the position
in the frame, and a time. If the latter is not set, the current time will be
used. The returned value is a record with a 'solved' field, which is F if the
source is always below or above the horizon. In that case the rise and set
fields will all have a string value. The record also returns a rise and set
record, with 'last' and 'utc' fields showing the rise and set times as epochs.
</description>
<example>
\begin{verbatim}
# NOT IMPLEMENTED
print "\t----\t riseset Ex 1 \t----"
print me.riseset(me.direction('sun'))
#[solved=T,
# rise=[last=[type=epoch, refer=LAST, m0=[value=0.0731388605, unit=d]],
#       utc=[type=epoch, refer=UTC, m0=[value=52085.8964, unit=d]]],
# set=[last=[type=epoch, refer=LAST, m0=[value=0.455732593, unit=d]],
#       utc=[type=epoch, refer=UTC, m0=[value=52086.2779, unit=d]]]]
print me.riseset(me.direction('sun'), qa.unit('80deg'))
#[solved=F,
# rise=[last=below, utc=below],
# set=[last=below, utc=below]]
print qa.form.long(me.riseset(me.direction('sun')).rise.utc.m0)
#21:30:47.439
#
\end{verbatim}
</example>
</method>

   <method type="function" name="posangle">
   <shortdescription>get position angle of two directions</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="m1">
     <any type="record"/>
     <description>direction of source (direction measure)</description>
     </param>

     <param xsi:type="any" direction="in" name="m2">
     <any type="record"/>
     <description>direction of other source (direction measure)</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/>
<shortdescription>position angle as angle quantity</shortdescription>
</returns>
<description>
posangle will give the position angle from a direction to another. I.e. the
angle in a direction between the direction to the North pole and the other
direction.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t posangle Ex 1 \t----")
a=me.direction('j2000','0deg','70deg')
b=me.direction('j2000','0deg','80deg')
print(me.posangle(a,b))
#{'value': -0.0, 'unit': 'deg'}
print(me.separation(a,b))
#{'value': 9.9999999999999893, 'unit': 'deg'}
tim=me.epoch('utc','today')
print(me.doframe(tim))
#True
pos=me.observatory('ATCA')
print(me.doframe(pos))
#True
print(me.posangle(a,b))
#{'value': -0.0, 'unit': 'deg'}
#
"""
\end{verbatim}
<!-- : me.posangle(a,b)
[value=7.99647705, unit=deg] ??? -->
</example>
</method>

   <method type="function" name="separation">
   <shortdescription>get separation angle between two directions</shortdescription>

<input>

     <param xsi:type="any" direction="in" name="m1">
     <any type="record"/>
     <description>direction of source (direction measure)</description>
     </param>

     <param xsi:type="any" direction="in" name="m2">
     <any type="record"/>
     <description>direction of other source (direction measure)</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/>
<shortdescription>position angle as angle quantity</shortdescription>
</returns>
<description>
separation will give the separation of a direction from another as an angle.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t separation Ex 1 \t----")
a=me.direction('j2000','0deg','70deg')
b=me.direction('j2000','0deg','80deg')
print(me.separation(a,b))
#{'value': 9.9999999999999893, 'unit': 'deg'}
tim = me.epoch('utc','today')              # set the time
print(me.doframe(tim))
#True
pos = me.observatory('ATCA')               # set where
print(me.doframe(pos))
#True
c=me.measure(b,'azel')                     # try with different type
print(me.separation(a,c))
#{'value': 10.000000000062277, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>




   <method type="function" name="addxvalue">
   <shortdescription>get some additional measure information</shortdescription>

<input>
     <param xsi:type="any" direction="in" name="a">
     <any type="record"/>
     <description>measures for which extra information is to be gotten</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
addxvalue will give some additional information about some measures as a vector
of quantities. It is used internally to get the rectangular coordinates of
measures that are normally given in angles. The casual user will probably in
general not interested in this function.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t addxvalue Ex 1 \t----")
a=me.observatory('atca')
print(a)
#{'m0': {'value': 2.6101423190348916, 'unit': 'rad'},
# 'm1': {'value': -0.5261379196128062, 'unit': 'rad'},
# 'm2': {'value': 6372960.2577234386, 'unit': 'm'},
# 'refer': 'ITRF',
# 'type': 'position'}
print(me.addxvalue(a))
#{'value': [-4750915.8370000012, 2792906.1819999996, -3200483.747], 'unit': 'm'}
print(me.addxvalue(me.epoch('utc','today')))
#{}
#
"""
\end{verbatim}
</example>
</method>



   <method type="function" name="type">
   <shortdescription>type of tool</shortdescription>

<input>
</input>
<returns xsi:type="string"/>
<description>
type will return the tool name.
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t type Ex 1 \t----")
print(me.type())
#'measures'
#
"""
\end{verbatim}
</example>
</method>





   <method type="function" name="done">
   <shortdescription>free resources used by tool.</shortdescription>

<input>
</input>
<returns xsi:type="bool"/>
<description>
In general you will not want to call this method.  It removes and then
recreates the default measures tool.
<!-- done will free the resources used by the tool. If the tool is the
default tool ({\em dm}) the done will only be executed if the kill
argument is true. -->
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t done Ex 1 \t----")
print(me.done())
#True
#
"""
\end{verbatim}
</example>
</method>



   <method type="function" name="ismeasure">
   <shortdescription>Check if measure</shortdescription>

<input>
     <param xsi:type="any" direction="in" name="v">
     <any type="record"/>
     <description>value to be tested</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
Checks if the operand is a correct measure
</description>
<example>
\begin{verbatim}
"""
#
print("\t----\t ismeasure Ex 1 \t----")
x=me.epoch('utc','today')
print(x)
#{'m0': {'value': 54056.043754386577, 'unit': 'd'},
# 'refer': 'UTC',
# 'type': 'epoch'}
print(me.ismeasure(x))
#True
y=me.getvalue(x)
print(y)
#{'m0': {'value': 54056.043754386577, 'unit': 'd'}}
print(me.ismeasure(y))
#False
print("Last example, exiting!")
exit()
#
"""
\end{verbatim}
</example>
</method>

</tool>



</casaxml>
"""
