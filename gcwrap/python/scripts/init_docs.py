import urllib2
import re
import webbrowser

class __doc(object):
    "command-line Plone help"

    def __init__( self ):

        version = "%d.%d.%d" % tuple(cu.version( )[:3])

        self.task_url = "https://casa.nrao.edu/PloneResource/casa-" + version + "/taskXml/"
        self.tool_url = "https://casa.nrao.edu/PloneResource/casa-" + version + "/toolXml/"
        self.toc_url = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-task-list"
        self.start_url = "https://casa.nrao.edu/casadocs/casa-" + version

        self.task_prefix = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-task-list/task_"
        self.tool_prefix = "https://casa.nrao.edu/casadocs/casa-" + version + "/global-tool-list/tool_"

        self.tasklist = [ ]
        self.toollist = [ ]

    def __call__( self, topic=None ):
        "open browser with documentation, try \"doc('toc')\""

        if len(self.tasklist) == 0:
            try:
                self.tasklist = re.findall("\w+.xml", urllib2.urlopen(self.task_url).read().decode())
            except:
                self.tasklist = [ ]

        if len(self.toollist) == 0:
            try:
                self.toollist = re.findall("\w+.xml", urllib2.urlopen(self.tool_url).read().decode())
            except:
                self.toollist = [ ]

        if type(topic) != str or topic == "toc":
            webbrowser.open_new_tab(self.toc_url)
        elif topic == "start":
            webbrowser.open_new_tab(self.start_url)
        elif topic+'.xml' in self.tasklist:
            webbrowser.open_new_tab(self.task_prefix+topic+"/parameters")
        elif topic+'.xml' in (toollist):
            webbrowser.open_new_tab(self.tool_prefix+topic+"/methods")
        else:
            webbrowser.open_new_tab(self.toc_url)

doc = __doc( )
