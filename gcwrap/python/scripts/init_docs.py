import os
import datetime
import webbrowser
import subprocess
import sys
import xml.etree.cElementTree as ET
import urllib2
from urlparse import urlparse

class __doc(object):
    "command-line Plone help"

    def __init__( self ):
        self.local_toc = None
        self.remote_toc = None
        version = "casa-%d.%d.%d" % tuple(cu.version( )[:3])
        self.remote_source_url = "https://casa.nrao.edu/casadocs/%s" % version
        self.remote_source_url_components = urlparse(self.remote_source_url)
        self.remote_toc_url = 'https://%s/PloneResource/%s/toc.xml' % (self.remote_source_url_components[1],version)

        self.local_toc_url = None if casa['dirs']['doc'] is None else casa['dirs']['doc'] + '/casa.nrao.edu/casadocs/toc.xml'
        self.local_start_path = "usingcasa/starting-casa.html"
        self.local_base = "/casadocs/%s/" % version

    def __welcome( self, welcome="\nOpening packaged documentation.\n" ):
        if welcome is not None:
            print welcome
        print "The most recent version of all CASA documentation is available online from:"
        print "\thttps://casa.nrao.edu/casadocs/\n"

    def __call__( self, sec=None ):
        "open browser with documentation, try \"doc('toc')\""

        ## for now, access to the plone site is turned off
        remote=False
        def show_toc( toc_dict ):
            width = max(len(key) for key in toc_dict.keys( ))+3
            for i in sorted(toc_dict.iterkeys( )):
                if "external" in toc_dict[i]['visibility']:
                    print "".join([i.ljust(width),toc_dict[i]['desc'].replace('\n','')])

        def entry_to_dict(acc,e):
            if e.tag == 'entry':
                acc[e.find('key').text] = {
                    'desc': e.find('desc').text,
                    'type': e.find('type').text,
                    'visibility': e.find('visibility').text,
                    'path': e.find('path').text }
            return acc

        if remote:
            if sec is None:
                return webbrowser.open("https://casa.nrao.edu/casadocs/")
            else:
                if self.remote_toc is None:
                    self.remote_toc =  reduce( entry_to_dict, ET.ElementTree(file=urllib2.urlopen(self.remote_toc_url)).getroot( ).getchildren( ), { } )
                if sec == 'toc':
                    show_toc(self.remote_toc)
                elif self.remote_toc.has_key(sec):
                    return webbrowser.open("https://casa.nrao.edu/casadocs/stable/" + self.remote_toc[sec]['path'])
                else:
                    print "Sorry '%s' is not a recognized section..." % sec
                    print "------------------------------------------------------------------------------"
                    show_toc(self.remote_toc)
        else:
            path = casa['dirs']['doc'] + "/casa.nrao.edu"
            if sec is None:
                homepage = "%s%s.html" % (path,self.remote_source_url_components[2])
                if os.path.exists(path):
                    self.__welcome( )
                    return webbrowser.open("file://" + homepage)
                else:
                    print "local documentation tree not found..."
                    self.__welcome(None)
                    return False
            else:
                if self.local_toc is None:
                    if self.local_toc_url is not None:
                        self.local_toc = reduce( entry_to_dict, ET.ElementTree(file=urllib2.urlopen("file://" + self.local_toc_url)).getroot( ).getchildren( ), { } )
                    else:
                        print "local documentation tree not found..."
                        self.__welcome(None)
                        return False
                if sec == 'toc':
                    show_toc(self.local_toc)
                elif sec == 'start':
                    self.__welcome( )
                    return webbrowser.open("file://" + path + self.local_base + self.local_start_path)
                elif self.local_toc.has_key(sec):
                    self.__welcome( )
                    return webbrowser.open("file://" + path + self.local_base + self.local_toc[sec]['path'])
                else:
                    self.__welcome(None)
                    print "Sorry '%s' is not a recognized section..." % sec
                    print "------------------------------------------------------------------------------"
                    show_toc(self.local_toc)


    def fetch( self ):
        if casa['dirs']['doc'] is None:
            print "casa['dirs']['doc'] has not been set..."
            return False
        if not os.path.exists(casa['dirs']['doc']):
            print ("directory %s does not exist..." % casa['dirs']['doc'])
            return False

        ## rename existing directory
        path = casa['dirs']['doc'] + "/casa.nrao.edu"
        if os.path.exists(casa['dirs']['doc'] + "/casa.nrao.edu"):
            now = datetime.datetime.now( ).isoformat('-')
            os.rename(path, path + "." + now)

        print "               source:  %s" % self.remote_source_url
        print "    table of contents:  %s" % self.remote_toc_url
        print "       download point:  %s" % path
        print "--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---"
        print "this will take some time..."
        print "relax..."
        print "do not hit ^C ..."
        print "do not expect output..."

        wgetcmd = "wget"

        if sys.platform == "darwin":
            wgetcmd = casa['dirs']['root'] + "/Resources/wget"

        tree = subprocess.call( [ wgetcmd, "--no-parent", "--no-check-certificate", "--html-extension", "--convert-links", "--recursive",
                                  "--level=inf", "--page-requisites", "-e", "robots=off", "--wait=0", "--quota=inf", "--reject",
                                  '*_form,RSS,*login*,logged_in,*logout*,logged_out,createObject*,select_default_page,selectViewTemplate*,object_cut,object_copy,object_rename,delete_confirmation,content_status_*,addtoFavorites,pdf.html,print.html',
                                  "--exclude-directories='search,*com_mailto*'", "--directory-prefix=" + casa['dirs']['doc'],
                                  "--convert-links", self.remote_source_url], stderr=subprocess.STDOUT, stdout=open(os.devnull,"w") )
        toc = subprocess.call( [ wgetcmd, self.remote_toc_url, "-O", self.local_toc_url], stderr=subprocess.STDOUT, stdout=open(os.devnull,"w") )
        if self.remote_source_url_components[1] != 'casa.nrao.edu':
            orig = os.getcwd()
            os.chdir(casa['dirs']['doc'])
            if os.path.exists(self.remote_source_url_components[1]):
                if os.path.exists('casa.nrao.edu'):
                    os.remove('casa.nrao.edu')
                os.symlink(self.remote_source_url_components[1],'casa.nrao.edu')
            else:
                print "warning, could not find mirror (%s/%s)" % (casa['dirs']['doc'],self.remote_source_url_components[1])
            os.chdir(orig)
        return (tree, toc)


doc = __doc( )
