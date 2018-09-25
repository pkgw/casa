<?xml version="1.0"?>
<xsl:stylesheet version="2.0"
          xmlns:aps="http://casa.nrao.edu/schema/psetTypes.html"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output omit-xml-declaration="yes"></xsl:output>
<xsl:param name="needscomma"/>
<xsl:param name="taskname"/>
<xsl:param name="paramname"/>
<xsl:param name="setme"/>
<xsl:param name="taskdescription"/>
<xsl:param name="unitsare"/>
<xsl:param name="arraytype"/>
<xsl:template match="*">
<xsl:apply-templates select="aps:task"/>
</xsl:template>
<xsl:template match="aps:task">
<xsl:param name="taskname"><xsl:value-of select="@name"/></xsl:param>
<xsl:param name="taskdescription"><xsl:value-of select="aps:shortdescription"/></xsl:param>
<xsl:text disable-output-escaping="yes">#
# This file was generated using xslt from its XML file
#
# Copyright 2007, Associated Universities Inc., Washington DC
#
import sys
import os
import string
import inspect
from odict import odict
import task_</xsl:text><xsl:value-of select="$taskname"/>
<xsl:text>
def </xsl:text><xsl:value-of select="@name"/><xsl:text>(</xsl:text><xsl:apply-templates select="aps:input"/>
<xsl:text>
        """</xsl:text>
<xsl:apply-templates select="aps:shortdescription"></xsl:apply-templates>
<xsl:apply-templates select="aps:example"/>
<xsl:text>
        """</xsl:text>
<xsl:text disable-output-escaping="yes">
        a=inspect.stack()
        stacklevel=0
        for k in range(len(a)):
          if (string.find(a[k][1], &apos;ipython console&apos;) &gt; 0) or (string.find(a[k][1], &apos;&lt;string&gt;&apos;) &gt;= 0):
                stacklevel=k
                break
        myf=sys._getframe(stacklevel).f_globals
        myf['__last_task'] = '</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">'
        myf['taskname'] = '</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">'
        ###
        myf['update_params'](func=myf['taskname'],printtext=False)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        function_signature_defaults=dict(zip(</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">.func_code.co_varnames,</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">.func_defaults))
        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
                        myf[key]=keyVal

</xsl:text>
<xsl:for-each select="aps:input">
<xsl:apply-templates select="aps:param"/>
<xsl:for-each select="aps:param">
<xsl:choose>
<xsl:when test="lower-case(@type)='boolarray'">
        if type(<xsl:value-of select="@name"/>)==bool: <xsl:value-of select="@name"/>=[<xsl:value-of select="@name"/>]<xsl:text/>
</xsl:when>
<xsl:when test="lower-case(@type)='intarray'">
        if type(<xsl:value-of select="@name"/>)==int: <xsl:value-of select="@name"/>=[<xsl:value-of select="@name"/>]<xsl:text/>
</xsl:when>
<xsl:when test="lower-case(@type)='stringarray'">
        if type(<xsl:value-of select="@name"/>)==str: <xsl:value-of select="@name"/>=[<xsl:value-of select="@name"/>]<xsl:text/>
</xsl:when>
<xsl:when test="lower-case(@type)='doublearray'">
        if type(<xsl:value-of select="@name"/>)==float: <xsl:value-of select="@name"/>=[<xsl:value-of select="@name"/>]<xsl:text/>
</xsl:when>
</xsl:choose>
</xsl:for-each>

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}
<xsl:for-each select="aps:param">
<xsl:choose>

<xsl:when test="@units">
        if type(<xsl:value-of select="@name"/>) == str :
           mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = myf[&apos;qa&apos;].quantity(<xsl:value-of select="@name"/>)
        else :
           mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = <xsl:value-of select="@name"/>
</xsl:when>
<xsl:otherwise>
        mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = <xsl:value-of select="@name"/>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>

<xsl:text disable-output-escaping="yes">
        pathname='file:///'+os.environ.get('CASAPATH').split()[0]+'/'+os.environ.get('CASAPATH").split()[1]+'/xml/'
        trec = myf['cu'].torecord(pathname+</xsl:text>&apos;<xsl:value-of select="$taskname"></xsl:value-of><xsl:text disable-output-escaping="yes">.xml&apos;)
</xsl:text>
<xsl:text disable-output-escaping="yes">
        if trec.has_key(&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;) and myf['cu'].verify(mytmp, trec[&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;]) :
          result = task_<xsl:value-of select="$taskname"/>.<xsl:value-of select="$taskname"/>(<xsl:call-template name="doargs2"/>)
</xsl:for-each>
<xsl:text disable-output-escaping="yes">
          saveinputs = myf['saveinputs']
          saveinputs(</xsl:text>&apos;<xsl:value-of select="$taskname"/>&apos;, &apos;<xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">.last&apos;)
        else :
          result = False
        return result
#
#
#
def </xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">_defaults(param=None):
        a=inspect.stack()
        stacklevel=0
        for k in range(len(a)):
          if (string.find(a[k][1], &apos;ipython console&apos;) &gt; 0) or (string.find(a[k][1], &apos;&lt;string&gt;&apos;) &gt;= 0):
                stacklevel=k
                break
        myf=sys._getframe(stacklevel).f_globals
        a = odict()
</xsl:text>
<xsl:for-each select="aps:input">
<xsl:call-template name="setdefaults"/>

<xsl:for-each select="aps:constraints">
<xsl:call-template name="setdefaults2"/>
</xsl:for-each>
</xsl:for-each>
<xsl:text disable-output-escaping="yes">

### This function sets the default values but also will return the list of
### parameters or the default value of a given parameter
        if(param == None):
                myf['__set_default_parameters'](a)
        elif(param == 'paramkeys'):
                return a.keys()
        else:
                if(a.has_key(param)):
                        return a[param]

#
#
def </xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">_check_params(param=None, value=None):
    a=inspect.stack()
    stacklevel=0
    for k in range(len(a)):
        if (string.find(a[k][1], &apos;ipython console&apos;) &gt; 0) or (string.find(a[k][1], &apos;&lt;string&gt;&apos;) &gt;= 0):
            stacklevel=k
            break
    myf=sys._getframe(stacklevel).f_globals

    return myf['cu'].verifyparam({param:value})
#
#
def </xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">_description(key=&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;, subkey=None):
        desc={&apos;</xsl:text><xsl:value-of select="$taskname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;: &apos;</xsl:text><xsl:value-of select="$taskdescription"/><xsl:text disable-output-escaping="yes">&apos;,
</xsl:text>
<xsl:for-each select="aps:input">
<xsl:call-template name="oneliners"/>

<xsl:text disable-output-escaping="yes">
              }
</xsl:text>
<xsl:for-each select="aps:constraints">
<xsl:text disable-output-escaping="yes">
#
# Set subfields defaults if needed
#
</xsl:text>
<xsl:call-template name="oneliners2"></xsl:call-template>
</xsl:for-each>
<xsl:text disable-output-escaping="yes">
        if(desc.has_key(key)) :
           return desc[key]
</xsl:text>
</xsl:for-each>
</xsl:template>

<xsl:template match="aps:input"> <xsl:call-template name="doargs"></xsl:call-template>
</xsl:template>
<xsl:template match="aps:shortdescription"><xsl:value-of select="."/></xsl:template>
<xsl:template match="aps:example"><xsl:value-of select="replace(., '\\.*\{verbatim\}', '')" disable-output-escaping="yes"/></xsl:template>

<xsl:template name="doargs">
<xsl:for-each select="aps:param"><xsl:value-of select="@name"/>=None</xsl:for-each>):
</xsl:template>

<xsl:template name="doargs2">
<xsl:for-each select="aps:param"><xsl:value-of select="@name"/><xsl:if test="position()&lt;last()">, </xsl:if></xsl:for-each>
</xsl:template>

<xsl:template match="aps:param">
<xsl:text>        </xsl:text><xsl:value-of select="@name"/> = myf[&apos;<xsl:value-of select="@name"/>&apos;]
</xsl:template>

<xsl:template name="oneliners">
<xsl:for-each select="aps:param">
<xsl:text disable-output-escaping="yes">               &apos;</xsl:text><xsl:value-of select="@name"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;: &apos;</xsl:text><xsl:value-of select="aps:description" disable-output-escaping="yes"/><xsl:text disable-output-escaping="yes">&apos;</xsl:text><xsl:if test="position()&lt;last()">,&#10;</xsl:if>
</xsl:for-each>
</xsl:template>

<xsl:template name="oneliners2">
<xsl:for-each select="aps:when">
<xsl:for-each select="aps:equals">
<xsl:call-template name="contextdesc">
<xsl:with-param name="paramname"><xsl:value-of select="@value"/></xsl:with-param>
</xsl:call-template>
</xsl:for-each>
</xsl:for-each>
</xsl:template>

<xsl:template name="contextdesc">
<xsl:param name="paramname"></xsl:param>
<xsl:for-each select="aps:default">
<xsl:if test="aps:description">
<xsl:text disable-output-escaping="yes">        if(subkey == &apos;</xsl:text><xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;):
          desc[&apos;</xsl:text><xsl:value-of select="@param"/><xsl:text disable-output-escaping="yes">&apos;] = &apos;</xsl:text><xsl:value-of select="aps:description" disable-output-escaping="yes"></xsl:value-of>&apos;
</xsl:if>
</xsl:for-each>
</xsl:template>

<xsl:template name="setdefaults">
<xsl:for-each select="aps:param">
<xsl:choose>
<xsl:when test="lower-case(@subparam)='yes'"></xsl:when>
<xsl:when test="lower-case(@subparam)='true'"></xsl:when>
<xsl:otherwise>
<xsl:choose>
<xsl:when test="lower-case(@type)='any'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = <xsl:call-template name="handlevalue"/><xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='variant'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = <xsl:call-template name="handlevalue"/><xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='record'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = <xsl:call-template name="handlevalue"/><xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='string'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = &apos;<xsl:value-of select="aps:value"/>&apos;<xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='boolarray'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = [<xsl:value-of select="aps:value"/>]<xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='stringarray'">
        <xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = [<xsl:apply-templates select="aps:value"><xsl:with-param name="arraytype"><xsl:value-of>string</xsl:value-of></xsl:with-param></xsl:apply-templates>]<xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='intarray'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = [<xsl:apply-templates select="aps:value"><xsl:with-param name="unitsare"><xsl:if test="@units"><xsl:value-of select="@units"/></xsl:if></xsl:with-param></xsl:apply-templates>]<xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='doublearray'">
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = [<xsl:apply-templates select="aps:value"><xsl:with-param name="unitsare"><xsl:if test="@units"><xsl:value-of select="@units"/></xsl:if></xsl:with-param></xsl:apply-templates>]<xsl:text>&#10;</xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:text>        a[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;]  = <xsl:if test="@units!=''">&apos;</xsl:if><xsl:value-of select="aps:value"/><xsl:if test="@units!=''"><xsl:value-of select="@units"/>&apos;</xsl:if><xsl:text>&#10;</xsl:text>
</xsl:otherwise>
</xsl:choose>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>
</xsl:template>

<xsl:template name="setdefaults2">

<xsl:for-each select="aps:when">
<xsl:call-template name="contextdefs">
<xsl:with-param name="paramname"><xsl:value-of select="@param"/></xsl:with-param>
</xsl:call-template>
</xsl:for-each>

</xsl:template>
<xsl:template name="contextdefs">
<xsl:param name="paramname"></xsl:param>
<xsl:text disable-output-escaping="yes">
        a[&apos;</xsl:text><xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] = {</xsl:text>
<xsl:for-each select="aps:equals">
<xsl:choose>
<xsl:when test="aps:default">
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;value&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">&apos;}, </xsl:text>
</xsl:when>
<xsl:when test="@type">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;value&apos;:</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">}, </xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;value&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">&apos;}, </xsl:text>
</xsl:otherwise>
</xsl:choose>

<xsl:for-each select="aps:default">
<xsl:text disable-output-escaping="yes">{&apos;</xsl:text> <xsl:value-of select="@param"/>&apos;:<xsl:call-template name="handlevalue"></xsl:call-template>}<xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:for-each>
<xsl:text>])</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:when>
<xsl:otherwise>
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;value&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">&apos;}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:when>
<xsl:when test="@type">
<xsl:if test="lower-case(@type)!='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;value&apos;:</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:if>
</xsl:when>
<xsl:otherwise>
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;value&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">&apos;}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:otherwise>
</xsl:choose>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>

<xsl:for-each select="aps:notequals">
<xsl:choose>
<xsl:when test="aps:default">
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;notvalue&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">&apos;}, </xsl:text>
</xsl:when>
<xsl:when test="@type">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;notvalue&apos;:</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">}, </xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"/><xsl:text disable-output-escaping="yes">:odict([{&apos;notvalue&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:call-template name="handlevalue"></xsl:call-template><xsl:text disable-output-escaping="yes">&apos;}, </xsl:text>
</xsl:otherwise>
</xsl:choose>

<xsl:for-each select="aps:default">
<xsl:text disable-output-escaping="yes">{&apos;</xsl:text> <xsl:value-of select="@param"/>&apos;:<xsl:call-template name="handlevalue"></xsl:call-template>}<xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:for-each>
<xsl:text>])</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:when>
<xsl:otherwise>
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;notvalue&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">&apos;}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:when>
<xsl:when test="@type">
<xsl:if test="lower-case(@type)!='string'">
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;notvalue&apos;:</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:if>
</xsl:when>
<xsl:otherwise>
<xsl:text disable-output-escaping="yes">
                    </xsl:text><xsl:value-of select="position()-1"></xsl:value-of><xsl:text disable-output-escaping="yes">:{&apos;notvalue&apos;:&apos;</xsl:text><xsl:value-of select="@value"/><xsl:text disable-output-escaping="yes">&apos;}</xsl:text><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:otherwise>
</xsl:choose>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>

<xsl:text disable-output-escaping="yes">}</xsl:text>
</xsl:template>

<xsl:template name="handlevalue">
<xsl:for-each select="aps:value">
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
<xsl:text>&apos;</xsl:text><xsl:value-of select="." disable-output-escaping="yes"></xsl:value-of><xsl:text>&apos;</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='vector' or ends-with(lower-case(@type), 'array')">
<xsl:text>[</xsl:text><xsl:choose>
<xsl:when test="count(aps:value)">
<xsl:for-each select="aps:value"><xsl:value-of select="."/><xsl:if test="position()&lt;last()">, </xsl:if></xsl:for-each>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="."/>
</xsl:otherwise>
</xsl:choose><xsl:text>]</xsl:text>
</xsl:when>
<xsl:when test="lower-case(@type)='record'">
<xsl:text>[</xsl:text><xsl:choose>
<xsl:when test="count(aps:value)">
<xsl:call-template name="handlevalue"/>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="."/>
</xsl:otherwise>
</xsl:choose><xsl:text>]</xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="." disable-output-escaping="yes"></xsl:value-of>
</xsl:otherwise>
</xsl:choose>
<xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:for-each>
</xsl:template>

<xsl:template match="aps:value">
<xsl:param name="unitsare"/>
<xsl:param name="arraytype"/>
<xsl:choose>
<xsl:when test="count(aps:value)">
        <xsl:for-each select="aps:value">
        <xsl:choose>
                <xsl:when test="$unitsare!=''">&apos;<xsl:value-of select="."/><xsl:value-of select="$unitsare"></xsl:value-of>&apos;<xsl:if test="position()&lt;last()">, </xsl:if></xsl:when>
                <xsl:when test="lower-case(@type)='string'">&apos;<xsl:value-of select="."/>&apos;<xsl:if test="position()&lt;last()">, </xsl:if></xsl:when>
                <xsl:otherwise><xsl:if test="$arraytype='string'">&apos;</xsl:if><xsl:value-of select="."/><xsl:if test="$arraytype='string'">&apos;</xsl:if><xsl:if test="position()&lt;last()">, </xsl:if></xsl:otherwise>
        </xsl:choose>
</xsl:for-each>
</xsl:when>
<xsl:otherwise><xsl:if test="$arraytype='string'">&apos;</xsl:if><xsl:value-of select="."/><xsl:if test="$arraytype='string'">&apos;</xsl:if><xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:otherwise>
</xsl:choose>
</xsl:template>

  <!-- templates go here -->
</xsl:stylesheet>
