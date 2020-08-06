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
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
import datetime
#from casac import *
import casac
import string
import time
import inspect
import numpy
from casa_stack_manip import stack_frame_find
from odict import odict
from types import *
from task_</xsl:text><xsl:value-of select="$taskname"/> import <xsl:value-of select="$taskname"/>
<xsl:text>
class </xsl:text><xsl:value-of select="@name"/><xsl:text>_cli_:</xsl:text>
<xsl:text>
    __name__ = "</xsl:text><xsl:value-of select="$taskname"/><xsl:text>"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (</xsl:text><xsl:value-of select="@name"/><xsl:text>_cli_,)
       self.__doc__ = self.__call__.__doc__</xsl:text>

       self.parameters=<xsl:text>{</xsl:text><xsl:apply-templates select="aps:input" mode="quotes"/>}
<xsl:text disable-output-escaping="yes">

    def result(self, key=None):
            #### and add any that have completed...
            return None

</xsl:text>
    def __call__<xsl:text>(self, </xsl:text><xsl:apply-templates select="aps:input" mode="noquotes"/>):
<xsl:text>
        """</xsl:text>
<xsl:apply-templates select="aps:shortdescription"></xsl:apply-templates>
<xsl:apply-templates select="aps:description"/>
<xsl:if test="aps:input">
<xsl:text>
        Arguments :
</xsl:text>
</xsl:if>
<xsl:for-each select="aps:input">
        <xsl:for-each select="aps:param">
                <xsl:if test="not(@visibility) or @visibility!='hidden'">
                <xsl:choose>
                <xsl:when test="aps:description">
                <xsl:text>              </xsl:text><xsl:value-of select="@name"/><xsl:text>:    </xsl:text><xsl:value-of select="aps:description" disable-output-escaping="yes"/>
                </xsl:when>
                <xsl:otherwise>
                <xsl:text>              </xsl:text><xsl:value-of select="@name"/><xsl:text>:    </xsl:text><xsl:value-of select="aps:shortdescription" disable-output-escaping="yes"/>
                </xsl:otherwise>
                </xsl:choose>
                <xsl:text>
</xsl:text>
<xsl:text>                 Default Value: </xsl:text><xsl:value-of select="aps:value"/>
<xsl:if test="aps:allowed">
        <xsl:text>
                   Allowed Values:</xsl:text>
        <xsl:for-each select="aps:allowed/aps:value">
        <xsl:text>
                                </xsl:text><xsl:value-of select="."/>
        </xsl:for-each>
</xsl:if>
<xsl:text>

</xsl:text>
</xsl:if>
        </xsl:for-each>
</xsl:for-each>
<xsl:if test="aps:returns">
        <xsl:for-each select="aps:returns">
                <xsl:text>      Returns: </xsl:text><xsl:value-of select="@type"/><xsl:text>
</xsl:text>
</xsl:for-each>
</xsl:if>
<xsl:text>
        Example :
</xsl:text>
<xsl:apply-templates select="aps:example"/>
<xsl:text>
        """</xsl:text>
<xsl:text disable-output-escaping="yes">
        if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        #casac = self.__globals__['casac']
        casalog = self.__globals__['casalog']
        casa = self.__globals__['casa']
        #casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = '</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">'
        self.__globals__['taskname'] = '</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">'
        ###
        self.__globals__['update_params'](func=self.__globals__['taskname'],printtext=False,ipython_globals=self.__globals__)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults={}
        else:
            function_signature_defaults=dict(zip(self.__call__.func_code.co_varnames[1:],self.__call__.func_defaults))
        useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
                        if (key != 'self') :
                           useLocalDefaults = True

        myparams = {}
        if useLocalDefaults :
           for item in function_signature_defaults.iteritems():
               key,val = item
               keyVal = eval(key)
               exec('myparams[key] = keyVal')
               self.parameters[key] = keyVal
               if (keyVal == None):
                   exec('myparams[key] = '+ key + ' = self.itsdefault(key)')
                   keyVal = eval(key)
                   if(type(keyVal) == dict) :
                      if len(keyVal) > 0 :
                         exec('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
                      else :
                         exec('myparams[key] = ' + key + ' = {}')

        else :
            print('')

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

        result = None

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}
<xsl:for-each select="aps:param">
<xsl:choose>

<xsl:when test="@units">
        if type(<xsl:value-of select="@name"/>) == str :
           mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = casac.casac.qa.quantity(<xsl:value-of select="@name"/>)
        else :
           mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = <xsl:value-of select="@name"/>
</xsl:when>
<xsl:otherwise>
        mytmp[&apos;<xsl:value-of select="@name"/>&apos;] = <xsl:value-of select="@name"/>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>
<xsl:text disable-output-escaping="yes">
        pathname='file://' + casa['dirs']['xml'] + '/'
        trec = casac.casac.utils().torecord(pathname+</xsl:text>&apos;<xsl:value-of select="$taskname"></xsl:value-of><xsl:text disable-output-escaping="yes">.xml&apos;)
</xsl:text>

<xsl:text disable-output-escaping="yes">
        casalog.origin(&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;)
        try :
          #if not trec.has_key(&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;) or not casac.casac.utils().verify(mytmp, trec[&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;]) :
            #return False

          casac.casac.utils().verify(mytmp, trec[&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']

          # Save .last file for this task execution. MPI servers don't write it (CASR-329).
          from mpi4casa.MPIEnvironment import MPIEnvironment
          do_full_logging = MPIEnvironment.is_mpi_disabled_or_client()
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs(</xsl:text>&apos;<xsl:value-of select="$taskname"/>&apos;, &apos;<xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">.last&apos;, myparams, self.__globals__,scriptstr=scriptstr, do_save_inputs=do_full_logging)</xsl:text>

          tname = '<xsl:value-of select="$taskname"/>'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          # Don't do telemetry from MPI servers (CASR-329)
          if do_full_logging and casa['state']['telemetry-enabled']:
              #casalog.poststat('Begin Task: ' + tname)
              task_starttime = str(datetime.datetime.now())
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else:
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')

          # Effective call to the task as defined in gcwrap/python/scripts/task_*
          result = <xsl:value-of select="$taskname"/>(<xsl:call-template name="doargs2"/>)

          if do_full_logging and casa['state']['telemetry-enabled']:
              task_endtime = str(datetime.datetime.now())
              casalog.poststat( 'Task ' + tname + ' complete. Start time: ' + task_starttime + ' End time: ' + task_endtime )
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')
</xsl:for-each>
<xsl:text disable-output-escaping="yes">
        except Exception as instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print('**** Error **** ',instance)
             tname = </xsl:text>'<xsl:value-of select="$taskname"/>'<xsl:text disable-output-escaping="yes">
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass
        casalog.origin('')
</xsl:text>
<xsl:for-each select="aps:output">
   <xsl:call-template name="checkoutput"/>
</xsl:for-each>
<xsl:text disable-output-escaping="yes">
        return result
#
#
#
#    def paramgui(self, useGlobals=True, ipython_globals=None):
#        """
#        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
#        """
#        import paramgui
#        if not hasattr(self, "__globals__") or self.__globals__ == None :
#           self.__globals__=stack_frame_find( )
#
#        if useGlobals:
#            if ipython_globals == None:
#                myf=self.__globals__
#            else:
#                myf=ipython_globals
#
#            paramgui.setGlobals(myf)
#        else:
#            paramgui.setGlobals({})
#
#        paramgui.runTask(&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;, myf['_ip'])
#        paramgui.setGlobals({})
#
#
#
#
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
        if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        if ipython_globals == None:
            myf=self.__globals__
        else:
            myf=ipython_globals

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
            if(paramvalue==None and subparam==None):
               if(a.has_key(param)):
                  return a[param]
               else:
                  return self.itsdefault(param)
            else:
               retval=a[param]
               if(type(a[param])==dict):
                  for k in range(len(a[param])):
                     valornotval='value'
                     if(a[param][k].has_key('notvalue')):
                        valornotval='notvalue'
                     if((a[param][k][valornotval])==paramvalue):
                        retval=a[param][k].copy()
                        retval.pop(valornotval)
                        if(subparam != None):
                           if(retval.has_key(subparam)):
                              retval=retval[subparam]
                           else:
                              retval=self.itsdefault(subparam)
                     else:
                        retval=self.itsdefault(subparam)
               return retval


#
#
    def check_params(self, param=None, value=None, ipython_globals=None):
      if ipython_globals == None:
          myf=self.__globals__
      else:
          myf=ipython_globals
#      print('param:', param, 'value:', value)
      try :
         if str(type(value)) != "&lt;type &apos;instance&apos;&gt;" :
            value0 = value
            value = myf['cu'].expandparam(param, value)
            matchtype = False
            if(type(value) == numpy.ndarray):
               if(type(value) == type(value0)):
                  myf[param] = value.tolist()
               else:
                  #print('value:', value, 'value0:', value0)
                  #print('type(value):', type(value), 'type(value0):', type(value0))
                  myf[param] = value0
                  if type(value0) != list :
                     matchtype = True
            else :
               myf[param] = value
            value = myf['cu'].verifyparam({param:value})
            if matchtype:
               value = False
      except Exception as instance:
         #ignore the exception and just return it unchecked
         myf[param] = value
      return value
#
#
    def description(self, key=&apos;</xsl:text><xsl:value-of select="$taskname"/><xsl:text disable-output-escaping="yes">&apos;, subkey=None):
        desc={&apos;</xsl:text><xsl:value-of select="$taskname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;: &apos;</xsl:text><xsl:value-of select="$taskdescription"/><xsl:text disable-output-escaping="yes">&apos;,
</xsl:text>
<xsl:for-each select="aps:input">
        <xsl:call-template name="oneliners"/><xsl:text disable-output-escaping="yes">
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
<xsl:text>
    def itsdefault(self, paramname) :
        a = {}
</xsl:text>
<xsl:for-each select="aps:input">
<xsl:call-template name="setdefaults3"/>
        #a = sys._getframe(len(inspect.stack())-1).f_globals
<xsl:for-each select="aps:constraints">
<xsl:call-template name="setdefaults4"/>
</xsl:for-each>
<xsl:text>
        if a.has_key(paramname) :
              return a[paramname]
</xsl:text>
</xsl:for-each>


<xsl:value-of select="$taskname"/>_cli = <xsl:value-of select="$taskname"/>_cli_()
</xsl:template>

<xsl:template match="aps:input" mode="noquotes"> <xsl:call-template name="doargs"></xsl:call-template>
</xsl:template>
<xsl:template match="aps:input" mode="quotes"> <xsl:call-template name="doqargs"></xsl:call-template>
</xsl:template>
<xsl:template match="aps:shortdescription"><xsl:value-of select="."/></xsl:template>
<xsl:template match="aps:description"><xsl:text>

        Detailed Description:
</xsl:text><xsl:value-of select="."/></xsl:template>
<xsl:template match="aps:example"><xsl:value-of select="replace(., '\\.*\{verbatim\}', '')" disable-output-escaping="yes"/></xsl:template>

<xsl:template name="checkoutput">
<xsl:choose>
        <xsl:when test="count(aps:param) &gt; 1">
        for arg in result :
           if not result.has_key(arg) :
                 throw('Missing output value '+arg)
</xsl:when>
</xsl:choose>
</xsl:template>

<xsl:template name="doargs">
        <xsl:for-each select="aps:param"><xsl:value-of select="@name"/>=None, </xsl:for-each>
</xsl:template>

<xsl:template name="doqargs">
        <xsl:for-each select="aps:param">&apos;<xsl:value-of select="@name"/>&apos;:None, </xsl:for-each>
</xsl:template>

<xsl:template name="doargs2">
<xsl:for-each select="aps:param"><xsl:value-of select="@name"/><xsl:if test="position()&lt;last()">, </xsl:if></xsl:for-each>
</xsl:template>

<xsl:template match="aps:param">
        <xsl:text>            myparams[&apos;</xsl:text><xsl:value-of select="@name"/>&apos;] = <xsl:value-of select="@name"/> = self.parameters[&apos;<xsl:value-of select="@name"/>&apos;]
</xsl:template>

<xsl:template name="oneliners">
<xsl:for-each select="aps:param">
        <xsl:text disable-output-escaping="yes">               &apos;</xsl:text><xsl:value-of select="@name"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;: &apos;</xsl:text><xsl:choose><xsl:when test="aps:shortdescription"><xsl:value-of select="aps:shortdescription" disable-output-escaping="yes"/></xsl:when><xsl:otherwise><xsl:value-of select="aps:description" disable-output-escaping="yes"/></xsl:otherwise></xsl:choose><xsl:text disable-output-escaping="yes">&apos;,&#10;</xsl:text></xsl:for-each>
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
<xsl:text>{</xsl:text><xsl:choose>
<xsl:when test="count(aps:value)">
<xsl:call-template name="handlevalue"/>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="."/>
</xsl:otherwise>
</xsl:choose><xsl:text>}</xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="." disable-output-escaping="yes"></xsl:value-of>
</xsl:otherwise>
</xsl:choose>
<xsl:if test="position()&lt;last()">, </xsl:if>
</xsl:for-each>
</xsl:template>

<xsl:template name="setdefaults3">
<xsl:for-each select="aps:param">
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
</xsl:for-each>
</xsl:template>

<xsl:template name="setdefaults4">
<xsl:for-each select="aps:when">
<xsl:call-template name="contextdefs2">
<xsl:with-param name="paramname"><xsl:value-of select="@param"/></xsl:with-param>
</xsl:call-template>
</xsl:for-each>
</xsl:template>


<xsl:template name="contextdefs2">
<xsl:param name="paramname"></xsl:param>
<xsl:for-each select="aps:equals">
<xsl:choose>
<xsl:when test="count(aps:default)">
<xsl:choose>
<xsl:when test="@type">
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> == &apos;<xsl:value-of select="@value"/>&apos;:
</xsl:when>
<xsl:otherwise>
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> == <xsl:value-of select="@value"/>:
</xsl:otherwise>
</xsl:choose>
</xsl:when>
<xsl:otherwise>
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> == &apos;<xsl:value-of select="@value"/>&apos;:
</xsl:otherwise>
</xsl:choose>
</xsl:when>
</xsl:choose>
<xsl:for-each select="aps:default">
<xsl:text disable-output-escaping="yes">            a[&apos;</xsl:text> <xsl:value-of select="@param"/>&apos;] = <xsl:call-template name="handlevalue"/>
<xsl:text>
</xsl:text>
</xsl:for-each>


</xsl:for-each>

<xsl:for-each select="aps:notequals">

<xsl:choose>
<xsl:when test="count(aps:default)">
<xsl:choose>
<xsl:when test="@type">
<xsl:choose>
<xsl:when test="lower-case(@type)='string'">
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> != &apos;<xsl:value-of select="@value"/>&apos;:
</xsl:when>
<xsl:otherwise>
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> != <xsl:value-of select="@value"/>:
</xsl:otherwise>
</xsl:choose>
</xsl:when>
<xsl:otherwise>
        if self.parameters[&apos;<xsl:value-of select="$paramname"></xsl:value-of><xsl:text disable-output-escaping="yes">&apos;] </xsl:text> != &apos;<xsl:value-of select="@value"/>&apos;:
</xsl:otherwise>
</xsl:choose>
</xsl:when>
</xsl:choose>

<xsl:for-each select="aps:default">
<xsl:text disable-output-escaping="yes">            a[&apos;</xsl:text> <xsl:value-of select="@param"/>&apos;] = <xsl:call-template name="handlevalue"/>
<xsl:text>
</xsl:text>
</xsl:for-each>
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
