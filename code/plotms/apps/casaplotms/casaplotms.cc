//# casaplotms.cc: Application for running plotms standalone.
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$
#include <display/QtViewer/QtApp.h>
#include <plotms/Gui/PlotMSPlotter.qo.h>
#include <plotms/PlotMS/PlotMS.h>
#if ! defined(WITHOUT_DBUS)
#include <plotms/PlotMS/PlotMSDBusApp.h>
#include <casadbus/utilities/Diagnostic.h>
#else
#include <plotms/PlotMS/grpcPlotMSAdaptor.qo.h>
#endif
#include <casacore/casa/System/AppState.h>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Plots/PlotMSPlot.h>
#ifndef NO_CRASH_REPORTER
#include <stdcasa/StdCasa/CrashReporter.h>
#endif

#include <signal.h>
#include <QDebug>
#include <casa/namespace.h>

#if defined(__APPLE__)
// for executable_path( )
#include <mach-o/dyld.h>
#else
// for executable_path( )
std::string read_link( const std::string &path ) {
    int buffer_size = 128;
    char *buffer = new char[buffer_size+1];
    int nchars = readlink( path.c_str( ), buffer, buffer_size );
    while ( nchars == buffer_size ) {
        buffer_size *= 2;
        delete [] buffer;
        buffer = new char[buffer_size+1];
        nchars = readlink( path.c_str( ), buffer, buffer_size );
    }
    std::string result;
    if ( nchars > 0 ) {
        buffer[nchars] = '\0';
        char *exe = realpath(buffer,NULL);
        result = exe;
        free(exe);
    }
    delete [] buffer;
    return result;
}
#endif

static std::string executable_path( ) {
#if defined(__APPLE__)
    uint32_t size = PATH_MAX;
    char *buffer = (char *) malloc(sizeof(char)*size);
    if ( _NSGetExecutablePath(buffer, &size) == -1 ) {
        ++size;
        buffer = (char *) realloc(buffer,sizeof(char)*size);
        if ( _NSGetExecutablePath(buffer, &size) != 0 ) {
            free(buffer);
            fprintf( stderr, "cannot discover path to executable...\n" );
            return "";
        }
    }
    char *exepath = realpath(buffer,NULL);
    std::string result(exepath);
    free(buffer);
    free(exepath);
    return result;
#else
    char buffer[256];
    sprintf( buffer, "/proc/%d/exe", getpid( ) );
    struct stat statbuf;
    if ( lstat( buffer, &statbuf ) == 0 && S_ISLNK(statbuf.st_mode) ) {
        return read_link(buffer);
    }
    return "";
#endif
}

using namespace casa;

class PlotmsDataState: public casacore::AppState {
public:

    PlotmsDataState(const std::list<std::string> &path ) : data_path(path) { }
    virtual bool initialized( ) const { return true; }
    virtual std::list<std::string> dataPath( ) const { return data_path; }
private:
    std::list<std::string> data_path;
};


int main(int argc, char* argv[]) {    
    std::string exepath(executable_path( ));

#if (0)
setbuf(stdout, NULL); /* for debugging - forces all printf() to flush immediately always */
#endif

#ifndef NO_CRASH_REPORTER
    String s = CrashReporter::initializeFromApplication(argv[0]);
    if (! s.empty()){
        cerr << "*** Warning during crash reporter initialization: " << s << endl;
    }
#endif
#if ! defined(WITHOUT_DBUS)
    // log arguments if CASA_DIAGNOSTIC_FILE environment variable is set...
    casa::dbus::diagnostic.argv( argc, argv );
#endif
    // Parameter defaults.
    String ms    = "",
           xaxis = PMS::axis(PMS::DEFAULT_XAXIS),
           yaxis = PMS::axis(PMS::DEFAULT_YAXIS),
           logfile = "", logfilter = "";
    PlotMSSelection select;
    PlotMSAveraging averaging;
    bool cachedImageSizeToScreenResolution = false, usePixels = false,
         casapy = false, debug = false,
         nopopups = false;
    bool showGui = true;
  
    // Parse arguments.
    String arg, arg2, arg3;
    size_t index;
    String ARG_HELP1 = "-h", 
           ARG_HELP2 = "--help", 
           ARG_VIS = "vis",
           ARG_XAXIS = "xaxis", 
           ARG_YAXIS = "yaxis", 
           ARG_CISTSR = "-c",
           ARG_CISTSR2 = "--cachedimagesizetoscreenresolution",
           ARG_PIXELS1 = "-p", 
           ARG_PIXELS2 = "--pixels",
#if ! defined(WITHOUT_DBUS)
           ARG_CASAPY = PlotMSDBusApp::APP_CASAPY_SWITCH,
           ARG_LOGFILE = PlotMSDBusApp::APP_LOGFILENAME_SWITCH,
           ARG_LOGFILTER = PlotMSDBusApp::APP_LOGFILTER_SWITCH,
#else
           casapy_address = "",
           ARG_SERVER = grpcPlotMS::APP_SERVER_SWITCH,
#endif
           ARG_DEBUG1 = "-d",
           ARG_DEBUG2 = "--debug",
           ARG_NOPOPUPS = "--nopopups",
    	   ARG_SHOWGUI = "--nogui";

	//
	// configure datapath for casacore and colormaps...
	//
	auto ends_with = []( const std::string& str, const std::string& ending ) {
		return ( str.size( ) >= ending.size( ) ) && equal( ending.rbegin( ), ending.rend( ), str.rbegin( ) );
	};

	if ( ends_with(exepath, "Contents/MacOS/CASAplotms") ||
         ends_with(exepath, "Contents/MacOS/casaplotms") ) {
		// initialize CASAviewer app data...
		if ( ! casacore::AppStateSource::fetch( ).initialized( ) ) {
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// Mac OSX  --  path is specific to package format
			// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
			// initialize CASAviewer app data...
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// generate path to data...
			std::string datapath(exepath);
			datapath.erase( datapath.end( ) -  16, datapath.end( ) );
			std::string pluginpath = datapath;             // save for later...
			datapath += "Resources/casa-data";
			// initialize casacore...
			std::list<std::string> datadirs;
			datadirs.push_back(datapath);
			casacore::AppStateSource::initialize(new PlotmsDataState(datadirs));
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// set up Qt Plugin Path
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			pluginpath += "Plugins";
			QCoreApplication::addLibraryPath(QString(pluginpath.c_str( )));
		}

	} else if ( ends_with(exepath, "/AppRun") ||
                ends_with(exepath, "/CASAplotms.app/usr/bin/CASAplotms") ||
                ends_with(exepath, "/casaplotms.app/usr/bin/casaplotms") ) {

		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// linux  --  path is specific to package format
		//
		//    .../AppRun implies AppImage bash script startup, e.g. from an unpacked AppImage
		//    .../CASAplotms.app/usr/bin/CASAplotms implies debugging or running from the
		//                                          build tree before it has been packaged
		//
		// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
		// initialize CASAplotms app data...
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// generate path to data...
		bool packed_app = ends_with(exepath, "/AppRun");
		std::string datapath(exepath);
		//     packed_app -> .../AppRun
		// not packed_app -> .../CASAplotms.app/usr/bin/CASAplotms
		datapath.erase( datapath.end( ) -  (packed_app ? 6 : 18), datapath.end( ) );
		std::string pluginpath = datapath;			   // save for later...
		datapath += "data";
		// initialize casacore...
		std::list<std::string> datadirs;
		datadirs.push_back(datapath);
		casacore::AppStateSource::initialize(new PlotmsDataState(datadirs));
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// set up Qt Plugin Path
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		pluginpath += "usr/lib/plugins";
		QCoreApplication::addLibraryPath(QString(pluginpath.c_str( )));
	}

    const vector<String>& selectFields = PlotMSSelection::fieldStrings(),
                          averagingFields = PlotMSAveraging::fieldStrings();

    
    for(int i = 1; i < argc; i++) {
        arg = arg2 = argv[i];
        arg2.downcase();
        
        if(arg.empty()) continue;
        
        if(arg2 == ARG_HELP1 || arg2 == ARG_HELP2) {
            cout << argv[0] << ": Stand-alone executable for CASA PlotMS."
                 << "\nAvailable arguments:\n"
                 
                 << "* " << ARG_HELP1 << " or " << ARG_HELP2 << "\n     "
                 << "Prints this message then exits."
                 
                 << "\n* " << ARG_VIS << "=[ms]\n     "
                 << "Path to MS used for initial plot."
                 
                 << "\n* " << ARG_XAXIS << "=[axis], " << ARG_YAXIS
                 << "=[axis]\n     "
                 << "Axes for initial plot (see documentation).";
            
            if(selectFields.size() > 0) cout << "\n*";
            for(unsigned int i = 0; i < selectFields.size(); i++) {
                if(i > 0) cout << ',';
                cout << " " << selectFields[i] << "=[val]";           
            }
            if(selectFields.size() > 0)
                cout << "\n     MS Selection parameters for initial plot.";
            
            if(averagingFields.size() > 0) cout << "\n*";
            for(unsigned int i = 0; i < averagingFields.size(); i++) {
                if(i > 0) cout << ',';
                cout << " avg" << averagingFields[i];
                if(PlotMSAveraging::fieldHasValue(PlotMSAveraging::field(
                        averagingFields[i]))) {
                    cout << "[=val]";
                }      
            }
            if(averagingFields.size() > 0)
                cout << "\n     MS Averaging parameters for initial plot.";
            

            
		 cout << "\n* "  << ARG_NOPOPUPS  << "\n     "
		 << "Use logger and status bar instead of showing error messages in a popup."
			     
                 << "\n* " << ARG_PIXELS1 << " or " << ARG_PIXELS2 << "\n     "
                 << "Use pixels instead of symbols for initial plot."
            
                 << "\n* " << ARG_CISTSR << " or " << ARG_CISTSR2 << "\n     "
                 << "Toggles setting the cached image size to screen "
                 << "resolution."
#if ! defined(WITHOUT_DBUS)
                 /*
                 // Don't advertise this switch...
                 << "\n* " << ARG_CASAPY << "\n     "
                 << "Sets up application to be controlled by casapy."
                 */
                 
                 << "\n* " << ARG_LOGFILE << "=[filename]\n     "
                 << "Sets the log file location (blank to use global)."
                 
                 << "\n* " << ARG_LOGFILTER << "=[priority]\n     "
                 << "Sets the log minimum priority filter."
#endif
                 << "\n* " << ARG_DEBUG1 << " or " << ARG_DEBUG2 << "\n     "
                 << "Turn on debugging log messages."

                 << "\n* " << ARG_SHOWGUI << "\n           "
                 << "Run without showing anything on the screen"

                 << endl;
            return 0;
        }
        
        if((index = arg.find("=")) < arg.size() - 1) {
            arg2 = arg.before(index); arg2.downcase();
            arg3 = arg.after(index);
#if defined(WITHOUT_DBUS)
            if (arg2 == ARG_SERVER) {
                casapy = true;
                casapy_address = arg3;
                continue;
            }
#endif
        } else if(arg2 == ARG_CISTSR || arg2 == ARG_CISTSR2) {
            cachedImageSizeToScreenResolution = true;
            continue;
        } else if(arg2 == ARG_PIXELS1 || arg2 == ARG_PIXELS2) {
            usePixels = true;
            continue;
        } else if(arg2 == ARG_NOPOPUPS) {
            nopopups = true;
            continue;
#if ! defined(WITHOUT_DBUS)
        } else if(arg2 == ARG_CASAPY) {
            casapy = true;
            continue;
#endif
        } else if (arg2 == ARG_SHOWGUI){
        	showGui = false;
            continue;
        }
        else if(arg2 == ARG_DEBUG1 || arg2 == ARG_DEBUG2) {
            debug = true;
            continue;
        } else if(i < argc - 1) {
            arg3 = argv[++i];
            if(arg3.size() == 0) continue;
            if(arg3[0] == '=' && arg3.size() > 1) arg3 = arg3.after(0);
            else if(arg3 == "=" && i < argc - 1) arg3 = argv[++i];
            else if(index != arg.size() - 1) continue;
        } else {
            arg3 = "";
            bool found = false;
            for(unsigned int i = 0; !found && i < averagingFields.size(); i++)
                if(PMS::strEq(arg2, "avg" + averagingFields[i], true))
                    found = true;
            if(!found) {
                cout << "Unrecognized argument: " << arg2 << endl;
                cout << "Please type 'casaplotms --help' for more info" << endl;
                return 0;
            }
        }
        
        if(arg2[arg2.size() - 1] == '=') arg2.erase(arg2.size() - 1);
        
        if(arg2 == ARG_VIS)         ms = arg3;
        else if(arg2 == ARG_XAXIS) xaxis = arg3;
        else if(arg2 == ARG_YAXIS) yaxis = arg3;
#if ! defined(WITHOUT_DBUS)
        else if(arg2 == ARG_LOGFILE) logfile = arg3;
        else if(arg2 == ARG_LOGFILTER) logfilter = arg3;
#endif
        else {
            bool found = false;
            for(unsigned int i = 0; !found && i < selectFields.size(); i++) {
                if(PMS::strEq(arg2, selectFields[i], true)) {
                    select.setValue(PlotMSSelection::field(selectFields[i]),
                                    arg3);
                    found = true;
                }
            }
            PlotMSAveraging::Field f;
            double d;
            for(unsigned int i = 0; !found && i < averagingFields.size(); i++){
                if(PMS::strEq(arg2, "avg" + averagingFields[i], true)) {
                    f = PlotMSAveraging::field(averagingFields[i]);
                    averaging.setFlag(f, true);
                    if(PlotMSAveraging::fieldHasValue(f) && !arg3.empty()) {
                        if(sscanf(arg3.c_str(), "%lf", &d) >= 1)
                            averaging.setValue(f, d);
                    }
                    found = true;
                }
            }

        }
    }

    // If run from casapy, don't let Ctrl-C kill the application.
    if(casapy) signal(SIGINT,SIG_IGN);
    
    // Although the qwt plotter backend can manage the QApplication, since
    // the PlotMSPlotter is a QMainWindow we need the QApplication regardless
    // of the plotting backend.
    QtApp::init();
    
    // Get priority filter.
    LogMessage::Priority p = debug? LogMessage::DEBUGGING: LogMessage::NORMAL;
    if(!logfilter.empty()) {
        for(int i = LogMessage::DEBUGGING; i < LogMessage::SEVERE; i++) {
            if(LogMessage::toString(LogMessage::Priority(i)) == logfilter) {
                p = LogMessage::Priority(i);
                break;
            }
        }
    }

    // Set up parameters for plotms.
    PlotMSParameters params(logfile, PlotLogger::FLAG_FROM_PRIORITY(p), p);
    
    if(cachedImageSizeToScreenResolution)
        params.setCachedImageSizeToResolution();
    
    // Set up plotms object.
    PlotMSApp plotmsapp(params, casapy, showGui
#if defined(WITHOUT_DBUS)
                        , casapy_address
#endif
                        );
   /* if(!casapy){
    	plotmsapp.showGUI(showGui); // don't automatically show for casapy
    }*/
    if (nopopups)
      plotmsapp.allowPopups(false);
    
    
    // check for hackjob env var CASAPLOTMS_NOPOPUPS set to "yes" (or
    // anything starting with 'Y')
    const char *ev = getenv("CASAPLOTMS_NOPOPUPS");
    if (ev && (ev[0]=='Y' || ev[0]=='y'))
      plotmsapp.allowPopups(false);
    
    // Set up parameters for plot.
    PlotMSPlotParameters plotparams = PlotMSPlot::makeParameters(&plotmsapp);

    PMS_PP_CALL(plotparams, PMS_PP_MSData, setFilename, ms)
    PMS_PP_CALL(plotparams, PMS_PP_MSData, setSelection, select)
    PMS_PP_CALL(plotparams, PMS_PP_MSData, setAveraging, averaging)
    PMS_PP_CALL(plotparams, PMS_PP_Cache, setXAxis, PMS::axis(xaxis))
    PMS_PP_CALL(plotparams, PMS_PP_Cache, setYAxis, PMS::axis(yaxis))
    
    if(usePixels) {
        //PlotSymbolPtr sym = PMS_PP_RETCALL(plotparams, PMS_PP_Display,
        //        unflaggedSymbol, NULL);
        PlotSymbolPtr sym = NULL;
        if(plotparams.typedGroup< PMS_PP_Display >() != NULL){
           sym = plotparams.typedGroup< PMS_PP_Display >()-> unflaggedSymbol ( ); 
        }
        sym->setSymbol(PlotSymbol::PIXEL);
        PMS_PP_CALL(plotparams, PMS_PP_Display, setUnflaggedSymbol, sym)
    }
    
    // Add the plot to plotms.
    plotmsapp.addOverPlot(&plotparams);
    
    // If we're connected to DBus, don't quite the application when the window
    // is closed.  This is somewhat risky in that if the remote applications
    // forget to tell this application to quit, it never will.
    if(casapy) QApplication::setQuitOnLastWindowClosed(false);
    
    //return plotms.execLoop();
    return QtApp::exec();
}
