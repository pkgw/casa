"""
<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" ?>
<casaxml xmlns="http://casa.nrao.edu/schema/psetTypes.html"
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
xsi:schemaLocation="http://casa.nrao.edu/schema/casa.xsd
file:///opt/casa/code/xmlcasa/xml/casa.xsd">

 
   <tool name="quanta" module="quanta">
      <shortdescription>quanta tool handles units and quantities
      </shortdescription>
       

<keyword>quanta</keyword>


<code>
<include>xmlcasa/casa/quanta_forward.h</include>
<private>
#include <xmlcasa/casa/quanta_private.h>
</private>
</code>


<!-- 
   <method type="constructor" name="quanta">
   <shortdescription>Construct quanta tool</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="host">
     <description>host on which to run tool</description>
     <value>valid host name string</value>
     <value></value>
     </param>

     <param xsi:type="bool" direction="in" name="forceneweserver">
     <description>force the use of a new server</description>
     <value>false</value>
     </param>
</input>
<returns xsi:type="quanta"/>
<description>
Create a quanta tool on the specified host (or by default the
host you are running on). 
</description>
</method>
-->


 
   <method type="function" name="convertfreq">
   <shortdescription>convert a frequency quantity to another unit</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
       <any type="variant"/>
     <description>quantity to convert</description>
     <value>1.0</value>
     </param>

     <param xsi:type="string" direction="in" name="outunit">
     <description>unit to convert to</description>
     <value>Hz</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
convertfreq converts a frequency quantity to another unit.
<!-- All units allowed in 
<link anchor="quanta:quanta.setformat.function">setformat</link> are allowed.
-->
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t convertfreq Ex 1 \t----"
print qa.convertfreq('5GHz','cm')
#{'value': 5.9958491599999997, 'unit': 'cm'}
print qa.convertfreq('5cm','GHz')
#{'value': 5.9958491599999997, 'unit': 'GHz'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="convertdop">
   <shortdescription>convert a doppler velocity quantity to another unit</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
         <any type="variant"/>
     <description>quantity to convert</description>
     <value>0.0</value>
     </param>

     <param xsi:type="string" direction="in" name="outunit">
     <description>unit to convert to</description>
     <value>km/s</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
convertfreq converts a velocity quantity to another unit. Units are either
 velocity or dimensionless.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t convertdop Ex 1 \t----"
print qa.convertdop('1','km/s')
#{'value': 299792.45799999998, 'unit': 'km/s'}
print qa.convertdop('10km/s','1') 
#{'value': 3.3356409519815205e-05, 'unit': '1'}
#
"""
\end{verbatim}
</example>
</method>


 
   <method type="function" name="quantity">
   <shortdescription>make a quantity from a string or from a numeric value
    and a unit string</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
	<any type="variant"/>
     <description>quantity or numeric or string to convert to quantity</description>
     </param>

     <param xsi:type="string" direction="in" name="unitname">
     <description>unit string if v numeric</description>
     <value></value>
     </param>

</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
quantity makes a quantity from a string, or from a value and a
string. Note that a function unit exists which is a synonym for
quantity. If only a string is given, it can be a scalar string.
The result will be a scalar quantity.
<!--, a vector of strings, or an array of strings. -->
<!--, or an array of quantities.-->
If a numeric value and a unit string
are given, the numeric value can be any numeric type, and can also be
a vector of numeric values.
<!--an array or a vector of them. See the introduction for an example of
all the possibilities.-->
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t quantity Ex 1 \t----"
tu = qa.quantity('1Jy')			# make quantity
print tu
#{'value': 1.0, 'unit': 'Jy'}
print qa.quantity(tu)			# also accepts a quantity
#{'value': 1.0, 'unit': 'Jy'}
tu = qa.unit('1Jy')			# make quantity with synonym
print tu
#{'value': 1.0, 'unit': 'Jy'}
print qa.quantity(-1.3, 'Jy')		# make quantity with separate value
#{'value': -1.3, 'unit': 'Jy'}
q1 = qa.quantity([8.57132661e+09, 1.71426532e+10], 'km/s')
print q1
#{'value': array([  8.57132661e+09,   1.71426532e+10]), 'unit': 'km/s'}
#
"""
\end{verbatim}
</example>
</method>
 

   <method type="function" name="getvalue">
   <shortdescription>get the internal value of a quantity</shortdescription>
   
<input>
     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>quantity</description>
     </param>
</input>
<returns xsi:type="doubleArray"/>
<description>
getvalue returns the internal value of a quantity. It also can handle an array
of quantities.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t getvalue Ex 1 \t----"
tu = qa.quantity(-1.3, 'Jy')         # make quantity
print tu
#{'value': -1.3, 'unit': 'Jy'}
print qa.getvalue(tu)
#-1.3 
print qa.getunit(tu)
#Jy 
a = qa.quantity([3,5],'cm')
print a
#{'value': array([ 3.,  5.]), 'unit': 'cm'}
print qa.getvalue(a)
#[3.0, 5.0]
#
"""
\end{verbatim}
</example>
<!--
%%a = qa.quantity("5m/s 7A")        # NOT IMPLEMENTED YET-->
%%print a
%%#[id=quant, shape=2] 
%%print qa.getvalue(a)
%%#[5 7]  
-->

</method>
 
   <method type="function" name="getunit">
   <shortdescription>get the internal unit of a quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>quantity</description>
     </param>
</input>
<returns xsi:type="string"/>
<description>
getunit returns the internal unit string of a quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t getunit Ex 1 \t----"
tu = qa.quantity(-1.3, 'Jy')         # make quantity
print tu
#{'value': -1.3, 'unit': 'Jy'}
print qa.getvalue(tu)
#-1.3 
print qa.getunit(tu)
#Jy 
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="canonical">
   <shortdescription>get canonical value of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value to convert</description>
     <value>1.0</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
canonical (with alias canon) gets the canonical value of a quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t canonical Ex 1 \t----"
print qa.canonical('1Jy')			# canonical value of a string
#{'value': 1e-26, 'unit': 'kg.s-2'}
print qa.canon(qa.quantity('1Jy'))		# canonical value of a unit
#{'value': 1e-26, 'unit': 'kg.s-2'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="canon">
   <shortdescription>get canonical value of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value to convert</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
canon gets the canonical value of a quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t canon Ex 1 \t----"
print qa.canon('1Jy')			        # canonical value of a string
#{'value': 1e-26, 'unit': 'kg.s-2'}
print qa.canonical(qa.quantity('1Jy'))		# canonical value of a unit
#{'value': 1e-26, 'unit': 'kg.s-2'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="convert">
   <shortdescription>convert a quantity to another unit</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
	<any type="variant"/>
     <description>quantity to convert</description>
     </param>

     <param xsi:type="any" direction="in" name="outunit">
	<any type="variant"/>
     <description>unit to convert to</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
convert converts a quantity to another unit. If no output unit given,
conversion is to canonical units
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t convert Ex 1 \t----"
tu = qa.quantity('5Mm/s')		# specify a quantity
print tu
#{'value': 5.0, 'unit': 'Mm/s'}
print qa.convert(tu, 'pc/a')		# convert it to parsec per year
#{'value': 0.0051135608266237404, 'unit': 'pc/a'}
print qa.convert(tu)			# convert to canonical units
#{'value': 5000000.0, 'unit': 'm.s-1'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="define">
   <shortdescription>define a new unit name</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="name">
     <description>name of unit to define</description>
     </param>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>quantity value of new unit</description>
     <value>1</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
define defines the name and value of a user defined unit
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t define Ex 1 \t----"
print qa.define('JY','1Jy')			# your misspelling
#True
print qa.define('VLAunit', '0.898 JY')	# a special unit using it
#True
print qa.quantity('5 VLAunit') 			# check its use
#{'value': 5.0, 'unit': 'VLAunit'}
print qa.convert('5 VLAunit','Jy')
#{'value': 4.4900000000000002, 'unit': 'Jy'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="map">
   <shortdescription>list known unit names and constants</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="v">
     <description>type of information to list - coded string</description>
     <value>all</value>
     </param>
</input>
<returns xsi:type="string">
</returns>
<description>
map lists the known mapping of units and constants. It has a single argument,
which can be a coded string (no-case, minimax match):
\begin{description}
\item[all] all of the following units (not constants): also the default 
\item[Prefix] known decimal prefixes
\item[SI] known SI units
\item[Customary] a set of customary units known to programs
\item[User] units defined by the user
\item[Constants] known constants (note: only 'const', 'Const', 'constants'
and 'Constants' recognised).
\end{description}
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t map Ex 1 \t----"
print qa.map('pre')			# list decimal prefixes
#        == Prefix ==== 20 ====
#        E         (exa)                        1e+18
#        G         (giga)                       1000000000
#        M         (mega)                       1000000
#        P         (peta)                       1e+15
#        T         (tera)                       1e+12
#        Y         (yotta)                      1e+24
#        Z         (zetta)                      1e+21
#        a         (atto)                       1e-18
#        c         (centi)                      0.01
#        d         (deci)                       0.1
#        da        (deka)                       10
#        f         (femto)                      1e-15
#        h         (hecto)                      100
#        k         (kilo)                       1000
#        m         (milli)                      0.001
#        n         (nano)                       1e-09
#        p         (pico)                       1e-12
#        u         (micro)                      1e-06
#        y         (yocto)                      1e-24
#        z         (zepto)                      1e-21
print qa.map('Constants')			# list known constants
#        == Constants ====
#        pi    3.14..                    3.14159 
#        ee    2.71..                    2.71828 
#        c     light vel.                2.99792e+08 m/s
#        G     grav. const               6.67259e-11 N.m2/kg2
#        h     Planck const              6.62608e-34 J.s
#        HI    HI line                   1420.41 MHz
#        R     gas const                 8.31451 J/K/mol
#        NA    Avogadro #                6.02214e+23 mol-1
#        e     electron charge           1.60218e-19 C
#        mp    proton mass               1.67262e-27 kg
#        mp_me mp/me                     1836.15 
#        mu0   permeability vac.         1.25664e-06 H/m
#        eps0  permittivity vac.         1.60218e-19 C
#        k     Boltzmann const           1.38066e-23 J/K
#        F     Faraday const             96485.3 C/mol
#        me    electron mass             9.10939e-31 kg
#        re    electron radius           2.8179e-15 m
#        a0    Bohr's radius             5.2918e-11 m
#        R0    solar radius              6.9599e+08 m
#        k2    IAU grav. const^2         0.000295912 AU3/d2/S0
#
"""
\end{verbatim}
</example>
</method>


   <method type="function" name="maprec">
   <shortdescription>create record containing list of known unit names and
   constants</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="v">
     <description>type of information to list - coded string</description>
     <value>all</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/>returns a list in a record
</returns>
<description>
maprec returns a record with the known mapping of units and constants. It has a single argument,
which can be a coded string (no-case, minimax match):
\begin{description}
\item[all] all of the following units (not constants): also the default 
\item[Prefix] known decimal prefixes
\item[SI] known SI units
\item[Customary] a set of customary units known to programs
\item[User] units defined by the user
\end{description}
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t maprec Ex 1 \t----"
p = qa.maprec('pre')			# list decimal prefixes
print p['Prefix_G']
#        G         (giga)               1000000000
s = qa.maprec('SI')		        # list SI units
print s['SI_Jy']
#Jy        (jansky)                     1e-26 kg.s-2
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="fits">
   <shortdescription>define some FITS units</shortdescription>
   
<input>
</input>
<returns xsi:type="bool"/>
<description>
fits defines some unit names used in reading and writing FITS files.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t fits Ex 1 \t----"
print qa.fits()
#True
print qa.map('user')
#        == User ====
#        BEAM      (dimensionless beam)         1 _
#        DAYS      (day)                        86400 s
#        DEG       (degree)                     0.0174532925199 rad
#        DEGREES   (degree)                     0.0174532925199 rad
#        HZ        (hertz)                      1 s-1
#        JY        (jansky)                     1e-26 kg.s-2
#        KELVIN    (kelvin)                     1 K
#        KELVINS   (kelvin)                     1 K
#        KM        (km)                         1000 m
#        M         (meter)                      1 m
#        METERS    (meter)                      1 m
#        PASCAL    (pascal)                     1 m-1.kg.s-2
#        PIXEL     (dimensionless pixel)        1 _
#        S         (second)                     1 s
#        SEC       (second)                     1 s
#        SECONDS   (second)                     1 s
#        VOLTS     (volt)                       1 m2.kg.s-3.A-1
#        YEAR      (year)                       31557600 s
#        YEARS     (year)                       31557600 s
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="angle">
   <shortdescription>show an angle as a formatted string</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle quantity value to output</description>
     </param>

     <param xsi:type="int" direction="in" name="prec">
     <description>number of digits shown</description>
     <value>0</value>
     </param>

     <param xsi:type="stringArray" direction="in" name="form">
     <description>formatting information in coded string array</description>
     <value></value>
     </param>

     <param xsi:type="bool" direction="in" name="showform">
     <description>show square brackets and separating ,</description>
     <value>false</value>
     </param>
</input>
<returns xsi:type="stringArray"/>
<description>
angle converts an angle quantity to a formatted string. The formatting
information is a precision (0 is default, 6 includes +-ddd.mm.ss) and a
string array of codes (no-case, minimax match):
Codes include:
\begin{description}
\item[clean] delete leading/trailing superfluous separators
\item[no\_d] do not show degrees part
\item[no\_dm] do not show degrees and minutes part
\item[dig2] show only 2 digits of degrees in angle format
\item[time] show as time (hh:mm:ss.ttt) rather than as angle
\end{description}
If a multi-dimensional value is given for the value $v$, the returned value
is a string vector of a length equal to last dimension. Each string has a
number of fields equal to the number of elements in all earlier
dimensions. If the {\em showform} is $T$, each vector element is surrounded
by a pair of square brackets if there is more than one entry, and fields are
separated by a ','.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t angle Ex 1 \t----"
tu = qa.quantity('5.7.12.345678')		# define an angle
print tu
#{'value': 5.1200960216666669, 'unit': 'deg'}
print qa.angle(tu)    				# default output
#+005.07.12 
print qa.angle(tu, prec=7)			# 7 digits
#+005.07.12.3 
print qa.angle(tu, prec=4)			# 4 digits
#+005.07. 
print qa.angle(tu, form=["tim","no_d"])		# as time, no hours shown
#:20:29 
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="time">
   <shortdescription>show a time (or date) as a formatted string</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>time quantity value to output</description>
     </param>

     <param xsi:type="int" direction="in" name="prec">
     <description>number of digits shown</description>
     <value>0</value>
     </param>

     <param xsi:type="stringArray" direction="in" name="form">
     <description>formatting information in coded string array</description>
     <value></value>
     </param>

     <param xsi:type="bool" direction="in" name="showform">
     <description>show square brackets and separating ,</description>
     <value>false</value>
     </param>
</input>
<returns xsi:type="stringArray"/>
<description>
time converts a time quantity to a formatted string. The formatting
information is a precision (0 is default, 6 includes hh.mm.ss) and a
string array of codes (no-case, minimax match):
Codes include:
\begin{description}
\item[clean] delete leading/trailing superfluous separators
\item[no\_d] do not show hours part
\item[no\_dm] do not show hours and minutes part
\item[ymd] include a date as yyyy/mm/dd (date is by default not shown)
\item[dmy] include a date as ddMMMyyyy (date is by default not shown)
\item[mjd] include a date as Modified Julian Day (date is by default not shown)
\item[fits] include a date and show time in FITS format: le from OS
\item[angle] show in angle (dd.mm.ss.ttt) rather than time format
\item[day] prefix day-of-week to output
\item[local] show local time rather than UTC (add timezone offset)
\item[no\_time] suppress printing of time part
\end{description}
If a multi-dimensional value is given for the value $v$, the returned value
is a string vector of a length equal to last dimension. Each string has a
number of fields equal to the number of elements in all earlier
dimensions. If the {\em showform} is $T$, each vector element is surrounded
by a pair of square brackets if there is more than one entry, and fields are
separated by a ','.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t time Ex 1 \t----"
tu = qa.quantity('today')		# a time
print tu
#{'value': 54175.708981504627, 'unit': 'd'}
print qa.time(tu)			# default format
#17:00:56
print qa.time(tu,form="dmy")  		# show date
#16-Mar-2007/17:00:56
print qa.time(tu,form=["ymd","day"])	# and day
#Fri-2007/03/16/17:00:56
print qa.time(tu,form="fits")           # FITS format    
#2007-03-16T17:00:56
print qa.time(tu,form=["fits","local"]) # local FITS format
#2007-03-16T10:00:56-07:00
print qa.time(tu,form=["ymd","local"])  # local time         
#2007/03/16/10:00:56
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="add">
   <shortdescription>add quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
add adds two quantities
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t add Ex 1 \t----"
print qa.add('5m', '2yd')   
#{'value': 6.8288000000000002, 'unit': 'm'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="sub">
   <shortdescription>subtract quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
sub subtracts two quantities
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t sub Ex 1 \t----"
print qa.sub('5m', '2yd')   
#{'value': 3.1712000000000002, 'unit': 'm'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="mul">
   <shortdescription>multiply quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>1</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
mul multiplies two quantities
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t mul Ex 1 \t----"
print qa.mul('5m', '3s')
#{'value': 15.0, 'unit': 'm.s'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="div">
   <shortdescription>divides quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>1</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
div divides two quantities
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t div Ex 1 \t----"
print qa.div('5m', '3s') 
#{'value': 1.6666666666666667, 'unit': 'm/(s)'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="neg">
   <shortdescription>negate quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     <value>1</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
neg negates a quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t neg Ex 1 \t----"
print qa.neg('5m')   
#{'value': -5.0, 'unit': 'm'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="norm">
   <shortdescription>normalise angle</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle quantity</description>
     </param>

     <param xsi:type="double" direction="in" name="a">
     <description>lower interval boundary</description>
     <value>-0.5</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
norm normalise angles in interval of $2\pi$ radians. The default interval is
from -0.5 to +0.5 of a full interval (i.e. from -180 to +180 degrees). The
lower end of the interval can be set as a fraction of $2\pi$
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t norm Ex 1 \t----"
print qa.norm('713deg')			#default normalisation
#{'value': -6.9999999999999716, 'unit': 'deg'}
print qa.norm('713deg', -2.5) 		# normalise to interval -900 - -540 deg
#{'value': -727.0, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="le">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
le compares two quantities for less than or equal.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t le Ex 1 \t----"
print qa.le('5m', '2yd')   
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="lt">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
lt compares two quantities for less than.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t lt Ex 1 \t----"
print qa.lt('5m', '2yd') 
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="eq">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
eq compares two quantities for equality.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t eq Ex 1 \t----"
print qa.eq('5m', '2yd')  
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="ne">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
ne compares two quantities for non equality.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t ne Ex 1 \t----"
print qa.ne('5m', '2yd')   
#True
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="gt">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
gt compares two quantities for greater than.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t gt Ex 1 \t----"
print qa.gt('5m', '2yd')   
#True
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="ge">
   <shortdescription>compare quantities</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     <value>0</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
ge  compares two quantities for greater than or equal.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t ge Ex 1 \t----"
print qa.ge('5m', '2yd') 
#True
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="sin">
   <shortdescription>sine of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
sin gives sine of angle quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t sin Ex 1 \t----"
print qa.sin('7deg')
#{'value': 0.12186934340514748, 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="cos">
   <shortdescription>cosine of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
cos gives cosine of angle quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t cos Ex 1 \t----"
print qa.cos('7deg')
#{'value': 0.99254615164132198, 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="tan">
   <shortdescription>tangent of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
tan gives tangent of angle quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t tan Ex 1 \t----"
print qa.tan('7deg')
#{'value': 0.1227845609029046, 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="asin">
   <shortdescription>arcsine of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>non-dimensioned quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
asin gives arcsine of non-dimensioned quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t asin Ex 1 \t----"
print qa.convert(qa.asin(qa.sin('7deg')), 'deg')
#{'value': 7.0, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="acos">
   <shortdescription>arccosine of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>non-dimensioned quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
acos gives arccosine of non-dimensioned quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t acos Ex 1 \t----"
print qa.convert(qa.acos(qa.cos('7deg')), 'deg')
#{'value': 7.0000000000000249, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="atan">
   <shortdescription>arctangent of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>non-dimensioned quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
atan gives arctangent of non-dimensioned quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t atan Ex 1 \t----"
print qa.convert(qa.atan(qa.tan('7deg')), 'deg')
#{'value': 7.0, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="atan2">
   <shortdescription>arctangent of two quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>non-dimensioned quantity</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>non-dimensioned quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
atan gives arctangent of two non-dimensioned quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t atan2 Ex 1 \t----"
print qa.convert(qa.atan2(qa.sin('7deg'), qa.cos('7deg')), 'deg')
#{'value': 7.0, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="abs">
   <shortdescription>absolute value of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
abs gives absolute value of quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t abs Ex 1 \t----"
print qa.abs('-5km/s')
#{'value': 5.0, 'unit': 'km/s'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="ceil">
   <shortdescription>ceil value of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
ceil gives ceiling value of quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t ceil Ex 1 \t----"
print qa.ceil('5.1AU')
#{'value': 6.0, 'unit': 'AU'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="floor">
   <shortdescription>floor value of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
floor gives flooring value of quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t floor Ex 1 \t----"
print qa.floor('-5.1AU')
#{'value': -6.0, 'unit': 'AU'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="log">
   <shortdescription>logarithm of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>dimensionless quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
log gives natural logarithm of dimensionless quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t log Ex 1 \t----"
print qa.log('2')
#{'value': 0.69314718055994529, 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="log10">
   <shortdescription>logarithm of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>dimensionless quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
log10 gives logarithm of dimensionless quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t log10 Ex 1 \t----"
print qa.log10('2')
#{'value': 0.3010299956639812, 'unit': ''}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="exp">
   <shortdescription>exponential of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>dimensionless quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
exp gives exponential value of dimensionless quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t exp Ex 1 \t----"
print qa.exp('2')
#{'value': 7.3890560989306504, 'unit': ''}
try:
  print qa.exp('2m')
except Exception, e:
  print "Caught an expected exception", e
#Caught an expected exception Quantum::exp illegal unit type 'm'
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="sqrt">
   <shortdescription>square root of quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>dimensionless quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
sqrt gives square root of quantity with only even powered dimensions
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t sqrt Ex 1 \t----"
print qa.sqrt('2m2')
#{'value': 1.4142135623730951, 'unit': 'm'}
try:
  print qa.sqrt('2s')
except Exception, e:
  print "Caught an expected exception", e
#Caught an expected exception UnitVal::UnitVal Illegal unit dimensions for root
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="compare">
   <shortdescription>compare dimensionality of units</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="any" direction="in" name="a">
        <any type="variant"/>
     <description>value</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
compare compares the dimensionality of units of two qauntities
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t compare Ex 1 \t----"
print qa.compare('5yd/a', '6m/s')  		# equal dimensions
#True
print qa.compare('5yd', '5s')		# unequal dimensions
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="check">
   <shortdescription>check for proper unit string</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="v">
     <description>value</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
check checks if the argument has a properly defined unit string
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t check Ex 1 \t----"
print qa.check('5AE/Jy.pc5/s')
#True
print qa.check('7MYs')
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="checkfreq">
   <shortdescription>check for proper frequency unit</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="cm">
        <any type="variant"/>
     <description>value</description>
     </param>
</input>
<returns xsi:type="bool">bool</returns>
<description>
checkfreq checks if the argument has a properly defined frequency interpretable
unit string
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t checkfreq Ex 1 \t----"
print qa.checkfreq('5GHz')
#True
print qa.checkfreq('5cm')  
#True
print qa.checkfreq('5cm/s2')
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="pow">
   <shortdescription>raise quantity to power</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="int" direction="in" name="a">
     <description>power</description>
     <value>1</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
pow raises a quantity to an integer power
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t pow Ex 1 \t----"
print qa.pow('7.2km/s', -3)
#{'value': 0.0026791838134430724, 'unit': '(km/s)-3'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="constants">
   <shortdescription>get a constant</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="v">
     <description>name</description>
     <value>pi</value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
constants gets a named constant quantity. Names (no-case, minimax) are:
\begin{verbatim}
        pi    3.14..                    3.14159 
        ee    2.71..                    2.71828 
        c     light vel.                2.99792e+08 m/s
        G     grav. const               6.67259e-11 N.m2/kg2
        h     Planck const              6.62608e-34 J.s
        HI    HI line                   1420.41 MHz
        R     gas const                 8.31451 J/K/mol
        NA    Avogadro #                6.02214e+23 mol-1
        e     electron charge           1.60218e-19 C
        mp    proton mass               1.67262e-27 kg
        mp_me mp/me                     1836.15 
        mu0   permeability vac.         1.25664e-06 H/m
        eps0  permittivity vac.         1.60218e-19 C
        k     Boltzmann const           1.38066e-23 J/K
        F     Faraday const             96485.3 C/mol
        me    electron mass             9.10939e-31 kg
        re    electron radius           2.8179e-15 m
        a0    Bohr's radius             5.2918e-11 m
        R0    solar radius              6.9599e+08 m
        k2    IAU grav. const^2         0.000295912 AU3/d2/S0
\end{verbatim}
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t constants Ex 1 \t----"
print qa.constants()
#{'unit': '', 'value': 3.1415926535897931}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="isangle">
   <shortdescription>check if valid angle or time quantity</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle/time quantity</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
isangle checks if the argument is a valid angle/time quantity.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t isangle Ex 1 \t----"
print qa.isangle(qa.constants('pi'))
#False
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="totime">
   <shortdescription>convert an angle (or a time) to a time</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle/time quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
totime converts an angle quantity (or a time) to a time quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t totime Ex 1 \t----"
print qa.totime('2d5m')
#{'value': 0.0057870370370370376, 'unit': 'd'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="toangle">
   <shortdescription>convert a time (or an angle) to an angle</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle/time quantity</description>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
toangle converts a time quantity (or an angle) to an angle quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t toangle Ex 1 \t----"
print qa.toangle('5h30m12.6')
#{'value': 82.552499999999995, 'unit': 'deg'}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="splitdate">
   <shortdescription>split a date/time into a record</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>angle/time quantity</description>
     </param>
</input>
<returns xsi:type="any"> <any type="record"/></returns>
<description>
splitdate splits a date/time quantity into a record with constituent fields
like year, yearday, month etc. All fields will be integer (to enable use as
index and easy personal formatting), with the exception of the {\em s} field
which is a double float. See the example for the fields returned.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t splitdate Ex 1 \t----"
print qa.splitdate('today')

#{'mjd': 54175.752367291658, 'week': 11, 'usec': 533999, 'hour': 18,
# 'min': 3, 'yearday': 75, 'msec': 533, 'month': 3, 's':
# 24.533999226987362, 'sec': 24, 'weekday': 5, 'year': 2007, 'monthday':
# 16} print qa.splitdate('183.33333333deg')
#{'mjd': 0.50925925925000004, 'week': 46, 'usec': 999999, 'hour': 12,
# 'min': 13, 'yearday': 321, 'msec': 999, 'month': 11, 's':
# 19.999999200003487, 'sec': 19, 'weekday': 3, 'year': 1858,
# 'monthday': 17}
#
"""
\end{verbatim}
</example>
</method>
 
   <method type="function" name="tos">
   <shortdescription>convert quantity to string</shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value</description>
     </param>

     <param xsi:type="int" direction="in" name="prec">
     <description>convert precision of value</description>
     <value>9</value>
     </param>
</input>
<returns xsi:type="string"/>
<description>
tos converts a quantity to a string with the precision defined with
the {\em setformat('prec')} (which defaults to 9). If the optional
{\em prec} argument is set to an integer value greater than 1, that
precision is used in the conversion
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t tos Ex 1 \t----"
a = qa.quantity('2.56 yd/s')
print a
#{'value': 2.5600000000000001, 'unit': 'yd/s'}
print qa.tos(a)
#2.560000000yd/s
a=qa.quantity(1./7, 'km/s')
print qa.tos(a)
#0.142857143km/s
print qa.tos(a,2)
#0.14km/s
print qa.tos(a,20)
#0.14285714285714284921km/s
print qa.tos(a)   
#0.142857143km/s
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
print "\t----\t type Ex 1 \t----"
print qa.type()
#
"""
\end{verbatim}
</example>
</method>


 
   <method type="function" name="done">
   <shortdescription>Free resources used by tool.  Current implementation
    ignores input parameter, does nothing and returns true</shortdescription>
   
<input>

     <param xsi:type="bool" direction="in" name="kill">
     <description>force kill of the default tool (ignored)</description>
     <value>false</value>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
Currently, this method is an NOP.
<!--
done will free the resources used by the tool. If the tool is the
default tool ({\em qa}) the done will only be executed if the kill
argument is true.
-->
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t done Ex 1 \t----"
print qa.done()
#True
print qa.done()
#True
print qa.done(kill=T)
#True
#
"""
\end{verbatim}
</example>
</method>

<method type="function" name="unit">
<shortdescription>quantity from value v and unit string</shortdescription>

<input>
     <param xsi:type="any" direction="in" name="v">
         <any type="variant"/>
     <description></description>
     </param>

     <param xsi:type="string" direction="in" name="unitname">
     <description></description>
     <value></value>
     </param>
</input>
<returns xsi:type="any"><any type="record"/></returns>
<description>
unit makes a quantity from a string, or from a value and a string.
Note that unit is a synonym for quantity.
</description>
</method>

<method type="function" name="isquantity">
<shortdescription>Check if quantity</shortdescription>

<input>
     <param xsi:type="any" direction="in" name="v">
         <any type="variant"/>
     <description>value to be tested</description>
     </param>
</input>
<returns xsi:type="bool"/>
<description>
Checks if the operand is a correct quantity
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t isQuantity Ex 1 \t----"
a = qa.quantity("5Jy")               # make a quantity
print a
#{'value': 5.0, 'unit': 'Jy'}
print qa.isquantity(a)                  # is it one?
#True
print qa.isquantity("5Jy")              # and this string?
#True
#
"""
\end{verbatim}
</example>
</method>

 
   <method type="function" name="setformat">
   <shortdescription>set format for output of numbers. 
   (NOT IMPLEMENTED YET!)</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="t">
     <description>type -coded string indicating which format
parameter to set</description>
     <value></value>
     </param>

     <param xsi:type="string" direction="in" name="v">
     <description>format parameter value - numeric or coded string,
depending on format type to be set</description>
     <value>F</value>
     </param>
</input>
<returns xsi:type="bool"/>

<!--
<description>
setformat allows to define output format parameters for a variety of
value types. The {\em type} is indicated by the first parameter, and
coded (see following list). The actual format parameter set depends
on the value parameter {\em v} given. If this value is specified as
the '...' string, the value will be asked interactively from the
user. The formats specified (or the defaults) are used in the
<link anchor="quanta:quanta.form.xxx.function">form</link> functions.
The formatting that can be specified:
\begin{description}
\item[prec] set the output precision for general numbers. If the
value given is an integer less than one, the precision will be set
to 6.
\item[aprec] set the output prescion for angles. If the
value given is an integer less than one, the precision will be set
to 6. A value of 6 will output up to the seconds field, with decimal
places if larger than 6.
\item[tprec] set the output prescion for times. If the
value given is an integer less than one, the precision will be set
to 6. A value of 6 will output up to the seconds field, with decimal
places if larger than 6.
\item[long] longitude display. Valid formats: 'hms', 'dms', 'deg',
'rad', 'd', '+deg'.
\item[lat] latitude display. Valid formats: 'hms', 'dms', 'deg',
'rad', 'd'.
\item[len] length display. Valid format: any length dimension
\item[dtime] time type display: 'local', 'utc', 'last' or 'solar'
\item[elev] elevation limit rise/set: angle quantity
\item[auto] auto frame update interval: time quantity
\item[vel] velocity: velocity unit
\item[freq] wave characteristics: frequency, length, time, angle/time,
/length, energy, impulse units
\item[dop] doppler type display: 'radio', 'opt', 'true'
\item[unit] default units: any valid unit 
\end{description}
</description>
<example>
\begin{verbatim}
- 1.234567890123456		# note that quanta set precision to 9
1.23456789			# (as compared to default Glish 6)
- a=qa.quantity(1.2345678901234567, 'km/s')
- qa.form.vel(a)                # display as velocity in default units
1.23456789 km/s 
- qa.setformat('vel', 'yd/a')   # try another unit
T 
- qa.form.vel(a)
4.26071737e+10 yd/a 
- qa.setformat('prec',12)       # or another precision
T 
- qa.form.vel(a)          
42607173719.8 yd/a 
- qa.setformat('unit', 'm/s')   # set generic units
T 
- qa.form.unit(a)               # and see what it looks like
1234.56789012 m/s 
- b=qa.quantity('5d2m3')       # an angle
- qa.form.long(b)               # as longitude gives it in hms
00:20:08.200 
- qa.setformat('long', 'deg')   # display longitude in degrees
                                # (rather than default hms)
T
- qa.form.long(b)               # and show it   
5.03416667 deg 
- qa.form.lat(b)                # as a latitude it is default dms
+005.02.03.000 
-
\end{verbatim}
</example>
-->
</method>



 
   <method type="function" name="getformat">
   <shortdescription>get current output format
   (NOT IMPLEMENTED YET!)</shortdescription>
   
<input>

     <param xsi:type="string" direction="in" name="t">
     <description>type - coded string</description>
     <value></value>
     </param>
</input>
<returns xsi:type="string"/>
<description>
getformat returns the current format value set for the different
format possibilities. See the
<link anchor="quanta:quanta.setformat.function">setformat</link> function for the
different format type descriptions. The known types are: \\
prec, aprec, tprec, long, lat, len, dtime, elev, auto, vel, freq,
dop, unit.
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t getformat Ex 1 \t----"
print qa.getformat('prec')
#6 
#setformat is NOT IMPLEMENTED YET!
#qa.setformat('prec', 12)	# set precision to 12 significant digits
#T 
#print qa.getformat('prec')                                             
#12 
print qa.getformat('long')
#hms
#
"""
\end{verbatim}
</example>
</method>



 
   <method type="function" name="formxxx">
   <shortdescription>Format a quantity using given format, allowed are hms, dms, deg, rad, +deg.
   </shortdescription>
   
<input>

     <param xsi:type="any" direction="in" name="v">
        <any type="variant"/>
     <description>value to be converted</description>
     </param>


     <param xsi:type="string" direction="in" name="format">
     <description>xxx can be hms, dms, deg, rad or +deg 
        </description>
     <value>dms</value>
     </param>

     
     <param type="int" direction="in" name="prec">
     <description># digits in fractional part of output string for dms,hms</description>
     <value>2</value>
     </param>

     
<!--Do we need this parameter
     <param xsi:type="bool" direction="in" name="showform">
     <description>show square brackets and separating , (not for
all)</description>
     <value>false</value>
     </param>
-->
</input>
<returns xsi:type="string"/>
If a multi-dimensional value is given for the value $v$ in the case of {\em
dtime}, {\em long} or {\em lat}, the returned value
is a string vector of a length equal to last dimension. Each string has a
number of fields equal to the number of elements in all earlier
dimensions. If the {\em showform} is $T$, each vector element is surrounded
by a pair of square brackets if there is more than one entry, and fields are
separated by a ','.
<description>
form.xxx (xxx can be lat, long, len, vel, freq, dtime, unit) will format the
input into a string using the global format information set by setformat().
</description>
<example>
\begin{verbatim}
"""
#
print "\t----\t formxxx Ex 1 \t----"
#qa.setformat('freq','cm')
#T 
qa.formxxx('freq',qa.quantity('5GHz'))
#form_xxx NOT IMPLEMENTED YET!
#5.99584916 cm 
print "Last example, exiting! ..."
exit()
#
"""
\end{verbatim}
</example>
</method>



</tool>


</casaxml>
"""
