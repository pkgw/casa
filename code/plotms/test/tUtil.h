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
#include <casa/BasicSL/String.h>
#include <casa/BasicMath/Math.h>
#include <casa/OS/EnvVar.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlot.h>

#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <QApplication>
#include <QDebug>

namespace casa {
class tUtil {
public:
	static casacore::String getFullPath( casacore::String fileName, casacore::String directory="plotms" ){
		//casacore::Path for data
		casacore::String aipsPath = casacore::EnvironmentVariable::get("CASAPATH");
		if (aipsPath.empty()) {
			aipsPath = casacore::EnvironmentVariable::get("AIPSPATH");
		}
		casacore::String dataPath;
		qDebug() << "Aips path="<<aipsPath.c_str();
		if ( !aipsPath.empty() ){
			//Last part of path returned appears to be a space rather
			//than a slash.  "/casa/trunk linux_64b"
			int spaceIndex = aipsPath.find( " ");
			if ( spaceIndex > 0 ){
				dataPath = aipsPath.before( spaceIndex );
			}
			else {
				dataPath = aipsPath;
			}
			dataPath.append( "/data/regression/unittest/");
			dataPath.append( directory );
			dataPath.append( "/" );
			qDebug() << "Datapath="<<dataPath.c_str();
			dataPath.append( fileName );
		}
		return dataPath;
	}

    static casacore::String getExportPath() {
        // create unique path in /tmp using pid
        casacore::String path = "/tmp/" + casacore::String::toString(getpid()) + "/";
        clearFile(path);
        mkdir(path.c_str(), 0775);
        return path;
    }

	static void updatePlot( PlotMSApp* app ){
		PlotMSPlotManager& manager = app->getPlotManager();
		int plotCount = manager.numPlots();
		if ( plotCount > 0 ){
		    PlotMSPlot* currentPlot = manager.plot( 0 );
		    PlotMSPlotParameters& params = currentPlot->parameters();
		    currentPlot->parametersHaveChanged(params,PlotMSWatchedParameters::ALL_UPDATE_FLAGS());
		}
	}

	static int clearFile( const casacore::String& fileName ){
        // remove works on files and directories
		int result = remove( fileName.c_str());
		if ( result == 0 ){
			qDebug() << fileName.c_str()<<" was deleted.";
		}
		return result;
	}

	static int exitMain( bool showGui ){
		int exitCode = 1;
		if ( showGui ){
		    QApplication::setQuitOnLastWindowClosed(false);
		    exitCode = QApplication::exec();
		}
		return exitCode;
	}

	static bool checkFile( casacore::String fileName, int minBytes, int maxBytes, int digest ){
		bool fileOK = true;
		ifstream ifile( fileName.c_str() );
		if ( ! ifile ){
			cerr<< "FAIL output file did not exist!"<<endl;
			fileOK = false;
		}
		else {
			struct stat filestatus;
			stat( fileName.c_str(), &filestatus );
			qDebug() << "Output file size is "<<filestatus.st_size;
			if ( filestatus.st_size < minBytes ){
				qDebug() << "FAIL output file size was too small min="<<minBytes;
				fileOK = false;
			}
			else if ( filestatus.st_size > maxBytes ){
				qDebug() << "FAIL output file size was too large max="<<maxBytes;
				fileOK = false;
			}
			else if ( digest > 0 ){
			/*
				        if(self.plotfile_hash):
				            self.assertEqual(
				                sha.new(open(self.plotfile_jpg, 'r').read()).hexdigest(),
				                self.plotfile_hash
				            )
				        else:
				            # store to check against following test results
				            self.plotfile_hash = sha.new(open(self.plotfile_jpg, 'r').read()).hexdigest()
			*/
			}
		}
		return fileOK;
	}

	static bool allEQDiv(const casacore::Array<casacore::Float> arr1,
		const casacore::Array<casacore::Float> arr2) {
		// compare elements of a divided array, in which some elements may be 
		// nan or inf (after divide by 0)
		casacore::IPosition arr1shape(arr1.shape()), arr2shape(arr2.shape());
		bool equal(arr1shape==arr2shape);
		if (equal) { // so far so good
			casacore::uInt ndata = arr1.nelements();
			casacore::Vector<casacore::Float> vec1 =
				arr1.reform(casacore::IPosition(1, ndata));
			casacore::Vector<casacore::Float> vec2 =
				arr2.reform(casacore::IPosition(1, ndata));
			// element-by-element comparison
			for (casacore::uInt i=0; i<ndata; ++i) {
				Float f1(vec1(i)), f2(vec2(i));
				if (casacore::isNaN(f1) || casacore::isNaN(f2) || 
					casacore::isInf(f1) || casacore::isInf(f2)) {
					// if one is nan/inf, both must be nan/inf
					equal &= ((casacore::isNaN(f1) || casacore::isInf(f1)) && 
							  (casacore::isNaN(f2) || casacore::isInf(f2)));
				} else {
					equal &= (f1==f2);
				}
			}
		}
		return equal;
	  }

	static casacore::Array<casacore::Float> getWtAmp(
		const casacore::Array<casacore::Float>& wt, 
		const casacore::Array<casacore::Float>& amp) {
		// apply per-corr wt to each row
		casacore::Array<casacore::Float> wtamp(amp.shape());
		casacore::uInt nrow(wtamp.shape()(2)), ncorr(2);
		for (casacore::uInt row=0; row<nrow; ++row) {
			for (casacore::uInt corr=0; corr<ncorr; ++corr) {
				casacore::Float corrWt = wt(IPosition(2,corr,row));
				casacore::Slicer rowslicer = Slicer(Slice(corr), Slice(), Slice(row));
				wtamp(rowslicer) = amp(rowslicer) * corrWt;
			}
		}
		return wtamp;
	}

	static void makeSigmaSpFromSigma(Array<Float>& sigmasp,
			Array<Float>& sigma, Int nchan) {
		// add channel axis to [ncorr,nrow] sigma array
		casacore::IPosition sigmashape = sigma.shape();
		casacore::uInt ncorr(sigmashape(0)), nrow(sigmashape(1));
		sigmasp.resize(IPosition(3, ncorr, nchan, nrow));
		for (uInt row=0; row<nrow; ++row) {
			for (Int chan=0; chan<nchan; ++chan) {
				sigmasp[row][chan] = sigma[row];
			}
		}
	}

private:
	tUtil(){};
	~tUtil();
};

}
