import urllib2
import re
import ssl
import webbrowser

class __doc(object):
    "command-line Plone help"

    def __init__( self ):

        version = "%d.%d.%d" % tuple(cu.version( )[:3])

        self.__task_url = "https://casa.nrao.edu/PloneResource/casa-" + version + "/taskXml/"
        self.__tool_url = "https://casa.nrao.edu/PloneResource/casa-" + version + "/toolXml/"
        self.__toc_url = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-task-list"
        self.__start_url = "https://casa.nrao.edu/casadocs/casa-" + version

        self.__task_prefix = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-task-list/task_"
        self.__tool_prefix = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-tool-list/tool_"

        self.__task_list = [ ]
        self.__tool_list = [ ]

    def __call__( self, topic=None ):
        "open browser with documentation, try \"doc('toc')\""

        if len(self.__task_list) == 0:
            try:
                ### osx rejects NRAO's CERT
                self.__task_list = re.findall("\w+.xml", urllib2.urlopen(self.__task_url).read().decode(),context=ssl.SSLContext())
            except:
                self.__task_list = [ ]

        if len(self.__tool_list) == 0:
            try:
                ### osx rejects NRAO's CERT
                self.__tool_list = re.findall("\w+.xml", urllib2.urlopen(self.__tool_url).read().decode(),context=ssl.SSLContext())
            except:
                self.__tool_list = [ ]

        if type(topic) != str or topic == "toc":
            webbrowser.open_new_tab(self.__toc_url)
        elif topic == "start":
            webbrowser.open_new_tab(self.__start_url)
        elif topic+'.xml' in self.__task_list:
            webbrowser.open_new_tab(self.__task_prefix+topic+"/parameters")
        elif topic+'.xml' in self.__tool_list:
            webbrowser.open_new_tab(self.__tool_prefix+topic+"/methods")
        else:
            webbrowser.open_new_tab(self.__toc_url)

doc = __doc( )
