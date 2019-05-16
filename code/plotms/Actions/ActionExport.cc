//# Copyright (C) 2009
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

#include "ActionExport.h"
#include <plotms/Client/Client.h>
#include <plotms/Actions/ActionFactory.h>
#include <plotms/Threads/ThreadController.h>
#include <plotms/Threads/ExportThread.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <casa/OS/Path.h>
#include <iomanip>
#include <fstream>

#include <QDebug>

using namespace casacore;
namespace casa {

ActionExport::ActionExport( Client* client )
: PlotMSAction( client ), format( PlotExportFormat::JPG, ""){
	itsType_= PLOT_EXPORT;
	interactive = false;
}

void ActionExport::setExportFormat( const PlotExportFormat& exportFormat ){
	format = exportFormat;
}

void ActionExport::setInteractive( bool interactive ){
	this->interactive = interactive;
}

bool ActionExport::loadParameters(){
	bool parametersLoaded = false;
	if ( client != NULL ){
		plots = client->getCurrentPlots();
		interactive = client->isInteractive();
		if ( plots.size() > 0  ){
			if ( !format.location.empty() ){
				parametersLoaded = true;
			}
		}
	}
    return parametersLoaded;
}

bool ActionExport::exportText( PlotMSApp* plotms ){
	Record rec;
	CountedPtr<PlotMSAction> action = ActionFactory::getAction( SEL_INFO, client );
	bool ok = action->doActionWithResponse(plotms, rec);
	if (rec.nfields() < 1) return ok;

	ofstream csv_file;
	csv_file.open(format.location.c_str());
	if (rec.isDefined("vis")) {
		casacore::Path fullvis(rec.asString("vis"));
		csv_file << "# vis: " << fullvis.baseName() << endl;
		rec.removeField("vis");
	}
	if (rec.isDefined("selection")) {
		Record selRec = rec.asRecord("selection");
	   	String field(selRec.asString("field")),
			spw(selRec.asString("spw")),
			timerange(selRec.asString("timerange")),
			uvrange(selRec.asString("uvrange")),
			antenna(selRec.asString("antenna")),
			scan(selRec.asString("scan")),
			corr(selRec.asString("corr")),
			array(selRec.asString("array")),
			observation(selRec.asString("observation")),
			intent(selRec.asString("intent")),
			feed(selRec.asString("feed")),
			msselect(selRec.asString("msselect"));
		if (!field.empty()) csv_file << "# field: " << field << endl;
		if (!spw.empty()) csv_file << "# spw: " << spw << endl;
		if (!timerange.empty()) csv_file << "# timerange: " << timerange << endl;
		if (!uvrange.empty()) csv_file << "# uvrange: " << uvrange << endl;
		if (!antenna.empty()) csv_file << "# antenna: " << antenna << endl;
		if (!scan.empty()) csv_file << "# scan: " << scan << endl;
		if (!corr.empty()) csv_file << "# correlation: " << corr << endl;
		if (!array.empty()) csv_file << "# array: " << array << endl;
		if (!observation.empty()) csv_file << "# observation: " << observation << endl;
		if (!intent.empty()) csv_file << "# intent: " << intent << endl;
		if (!feed.empty()) csv_file << "# feed: " << feed << endl;
		if (!msselect.empty()) csv_file << "# msselect: " << msselect << endl;
		rec.removeField("selection");
	}

	if (rec.isDefined("averaging")) {
		Record avgRec = rec.asRecord("averaging");
		if (avgRec.asBool("channel")) {
			csv_file << "# channel average: " << avgRec.asDouble("channelValue") << endl;
		} else {
			csv_file << "# channel average: None" << endl;
		}
		if (avgRec.asBool("time")) {
			csv_file << "# time average: " << avgRec.asDouble("timeValue") << "s";
			if (avgRec.asBool("scan"))
				csv_file << " scan: True";
			if (avgRec.asBool("field"))
				csv_file << " field: True";
			csv_file << endl;
		} else {
			csv_file << "# time average: None" << endl;
		}
		if (avgRec.asBool("baseline"))
			csv_file << "# baseline average: True" << endl;
		if (avgRec.asBool("antenna"))
			csv_file << "# antenna average: True" << endl;
		if (avgRec.asBool("spw"))
			csv_file << "# spw average: True" << endl; 
		if (avgRec.asBool("scalar"))
			csv_file << "# scalar average: True" << endl;
		rec.removeField("averaging");
	}

	if (rec.isDefined("transformations")) {
		Record transRec = rec.asRecord("transformations");
		String frame(transRec.asString("Frame")),
			   veldef(transRec.asString("Veldef"));
		Double restfreq(transRec.asDouble("RestFreq")),
			   dx(transRec.asDouble("XpcOffset")),
			   dy(transRec.asDouble("YpcOffset"));
		if (!frame.empty())
			csv_file << "# frequency frame: " << frame << endl;
		csv_file << "# veldef: " << veldef << endl;
		if (restfreq > 0.0)
			csv_file << "# restfreq: " << restfreq << endl;
		if (dx > 0.0 || dy > 0.0)
			csv_file << "# shift: [" << dx << "," << dy << "]" << endl;
		rec.removeField("transformations");
	}
	if (rec.isDefined("calibration")) {
		csv_file << "# cal library: " << rec.asString("calibration") << endl;
		rec.removeField("calibration");
	}
	if (rec.isDefined("iteraxis")) {
		csv_file << "# iteraxis: " << rec.asString("iteraxis") << endl;
		rec.removeField("iteraxis");
	}

	bool verbose = format.verbose;
	const String X_AXIS("xaxis");
	const String Y_AXIS("yaxis");
	// for iterated plots on a grid:
	for (uInt iplot=0; iplot<rec.nfields(); ++iplot) {
		Record plotRecord = rec.subRecord(iplot);
		// print plot number and iteration
		csv_file << "# From plot " << iplot;
		if (plotRecord.isDefined("iteration")) {
			csv_file << " iteration " << plotRecord.asInt("iteration");
			plotRecord.removeField("iteration");
		}
		csv_file << endl;
		for(uInt datacount = 0; datacount < plotRecord.nfields(); ++datacount) {
			// one record per datacount (x/y data, in case of multiple y-axes)
			Record dataRecord = plotRecord.subRecord(datacount);
			// when no points are plotted (e.g. all flagged) no data is returned!
			// print header
			csv_file << "# x y ";
			bool hasData(false), hasChan(false), hasAnt2(false), hasFreq(false);
			if (dataRecord.isDefined("0")) {
			   	hasData = true;
				Record firstPoint = dataRecord.subRecord( "0" );
				hasChan = firstPoint.isDefined("chan");
				hasAnt2 = firstPoint.isDefined("ant2");
				hasFreq = firstPoint.isDefined("freq");
			}

			if (verbose) {
				// print more header
				if (hasChan) csv_file << "chan ";
				csv_file << "scan field ant1 ";
				if (hasAnt2) csv_file << "ant2 ";
				csv_file << "ant1name ";
				if (hasAnt2) csv_file << "ant2name ";
				csv_file << "time ";
				if (hasFreq) csv_file << "freq ";
				csv_file << "spw corr obs";
			}
			csv_file << endl;

			// print x/y axes and units for other headers
			String xaxisstr(""), yaxisstr("");
			if ( dataRecord.isDefined( X_AXIS )){
				xaxisstr = dataRecord.asString(X_AXIS);
				dataRecord.removeField( X_AXIS );
			}
			if ( dataRecord.isDefined( Y_AXIS)){
				yaxisstr = dataRecord.asString(Y_AXIS);
				dataRecord.removeField( Y_AXIS );
			}
			if ( xaxisstr.length() > 0 || yaxisstr.length() > 0 )
				csv_file << "# " << xaxisstr << " " << yaxisstr;
			if (verbose) {
				if (hasChan) csv_file << " None"; // chan
				csv_file << " None None None"; // scan field ant1
				if (hasAnt2) csv_file << " None"; // ant2
				csv_file << " None"; // ant1name
				if (hasAnt2) csv_file << " None"; // ant2name
				csv_file << " MJD(seconds)"; // time
				if (hasFreq) csv_file << " GHz"; // freq
				csv_file << " None None None";  // spw corr obs
			}
			csv_file << endl;
			if (!hasData)
				csv_file << "# No data to export" << endl;

			// one field per point (if any)
			for(uInt _field = 0; _field < dataRecord.nfields(); ++_field) {
				String field_str = String::toString(_field);
				Record fieldRecord = dataRecord.subRecord( field_str );
				Double x = fieldRecord.asDouble("x");
				Double y = fieldRecord.asDouble("y");
				int precision = csv_file.precision();

				if(xaxisstr == "Time") {
					csv_file << std::setprecision(3) << std::fixed
							<< x << " ";
					csv_file.unsetf(ios_base::fixed);
					csv_file.precision(precision);
				} else if(xaxisstr == "Frequency") {
					csv_file << std::setprecision(9) << std::fixed
							<< x << " ";
					csv_file.unsetf(ios_base::fixed);
					csv_file.precision(precision);
				} else {
					csv_file << x << " ";
				}

				if(yaxisstr == "Time") {
					csv_file << std::setprecision(3) << std::fixed
							<< y << " ";
					csv_file.unsetf(ios_base::fixed);
					csv_file.precision(precision);
				} else if(yaxisstr == "Frequency") {
					csv_file << std::setprecision(9) << std::fixed
							<< y << " ";
					csv_file.unsetf(ios_base::fixed);
					csv_file.precision(precision);
				} else {
					csv_file << y;
				}

				if (verbose) {
					// collect values from Record...
					Int chan = -1, ant2 = -1;
					String ant2name;
					Double freq;
					if (hasChan) chan = fieldRecord.asInt("chan");
					Int scan = fieldRecord.asInt("scan");
					Int field = fieldRecord.asInt("field");
					Int ant1 = fieldRecord.asInt("ant1");
					if (hasAnt2) ant2 = fieldRecord.asInt("ant2");
					String ant1name = fieldRecord.asString("ant1name");
					if (hasAnt2) ant2name = fieldRecord.asString("ant2name");
					Double time = fieldRecord.asDouble("time");
					Int spw = fieldRecord.asInt("spw");
					if (hasFreq) freq = fieldRecord.asDouble("freq");
					String corr = fieldRecord.asString("corr");
					Int obsId = fieldRecord.asInt("obsid");
					// ... and print them
					csv_file << " ";
					if (hasChan) csv_file << chan << " ";
					csv_file << scan << " " << field << " " << ant1 << " ";
					if (hasAnt2) csv_file  << ant2 << " "; 
					csv_file << ant1name << " ";
					if (hasAnt2) csv_file << ant2name << " ";
					csv_file << std::setprecision(3) << std::fixed
							<< time << " ";
					if (hasFreq) csv_file << std::setprecision(9) << std::fixed
							<< freq << " ";
					csv_file.unsetf(ios_base::fixed);
					csv_file.precision(precision);
					csv_file << spw << " " << corr << " " << obsId; 
				}
				csv_file << endl;
			}
		}
	}
	csv_file.close();
	return ok;
}

PlotExportFormat ActionExport::adjustFormat( PlotExportFormat::Type t){
	PlotExportFormat exportFormat(t, format.location);
	exportFormat.resolution = format.resolution;

	exportFormat.dpi = format.dpi;
	if(exportFormat.dpi <= 0){
		exportFormat.dpi = -1;
	}
	exportFormat.width = format.width;
	if(exportFormat.width <= 0){
		exportFormat.width = -1;
	}
	exportFormat.height = format.height;
	if(exportFormat.height <= 0){
		exportFormat.height = -1;
	}
	return exportFormat;
}

bool ActionExport::doActionSpecific(PlotMSApp* plotms){
	bool ok = true;

	String form = PlotExportFormat::exportFormat( format.type );
	PlotExportFormat::Type t = PlotExportFormat::exportFormat(form, &ok);
	if(!ok) {
		t = PlotExportFormat::typeForExtension(format.location, &ok);
		if(!ok) {
			itsDoActionResult_ = "Invalid format extension for filename '"+
						format.location + "'!";
			return ok;
		}
	}

	if ( t == PlotExportFormat::TEXT ){
		ok = exportText(plotms);
	}
	else {

		PlotExportFormat exportFormat = adjustFormat( t );
		// TODO export fix screen resolution
		// Quick hack for screen resolution images.  Taking a screenshot without
		// drawing the items is basically impossible in the non-main (non-GUI) thread,
		// so for now just turn on high resolution so that it has to draw each
		// items.  This isn't ideal because it is slow, but for now it's better to
		// have something that works and is slow than something that doesn't work.

		if((exportFormat.type == PlotExportFormat::JPG ||
			exportFormat.type == PlotExportFormat::PNG) &&
			exportFormat.resolution == PlotExportFormat::SCREEN) {
			cout << "NOTICE: Exporting to images in screen resolution is currently"
					<< " not working.  Switching to high resolution (which is slower,"
					<< " but works)." << endl;
			exportFormat.resolution = PlotExportFormat::HIGH;
		}

		// Tell parent the format so it can be used in msg to user
		plotms->setExportFormat(exportFormat);

		ExportThread* exportThread = new ExportThread();
		exportThread->setExportFormat( exportFormat );
		exportThread->setPlots( plots );
		setUpClientCommunication( exportThread, -1 );

		setUseThreading( false );

		ok = initiateWork( exportThread );
		if ( threadController == NULL ){
			delete exportThread;
		}
	}
	return ok;
}



ActionExport::~ActionExport() {
}

using namespace casacore;
} /* namespace casa */
