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

#include "ActionSummary.h"
#include <plotms/Client/Client.h>
#include <plotms/PlotMS/PlotMSLabelFormat.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlotParameters.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTSummary.h>
#include <synthesis/CalTables/CalSummary.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MSOper/MSSummary.h>
#include <casa/Logging/LogFilter.h>
#include <QDebug>

using namespace casacore;
namespace casa {

ActionSummary::ActionSummary( Client* client )
	:PlotMSAction( client ){

	itsType_=MS_SUMMARY;
	verbose = false;
	summaryType = PMS::S_ALL;
}

bool ActionSummary::loadParameters(){
	bool parametersLoaded = false;
	if ( !filename.empty() ){
		parametersLoaded = true;
	}
	return parametersLoaded;
}

void ActionSummary::setSummaryType( PMS::SummaryType type ){
	summaryType = type;
}

void ActionSummary::setCTSummaryType( PMS::CTSummaryType type ){
	CTsummaryType = type;
}

void ActionSummary::setFile( String file ){
	filename = file;
}

void ActionSummary::setVerbose( bool verbose ){
	this->verbose = verbose;
}


bool ActionSummary::doActionSpecific(PlotMSApp* plotms) {
	bool success = false;
	bool reenableGlobal = false;
	try {
		// If not, exit.
		if(filename.empty()) {
			itsDoActionResult_ = "MS has not been opened/set yet!";
			return false;
		}

		// Set up log objects.
		LogSink sink(LogFilter(plotms->getLogger()->filterMinPriority()));
		if(!plotms->getLogger()->usingGlobalSink()) {
			LogSinkInterface* ic = plotms->getLogger()->localSinkCopy();
			sink.localSink(ic);
			// Temporarily disable global log sink if we're not using it,
			// since MSSummary posts to both (how annoying).
			PlotLogger::disableGlobalSink();
			reenableGlobal = true;
		}
		LogIO log(LogOrigin(PMS::LOG_ORIGIN,PMS::LOG_ORIGIN_SUMMARY),sink);

		bool vb = verbose;
		Table tab(filename);
		if (tab.keywordSet().isDefined("ParType")) {
			NewCalTable ct = NewCalTable(filename, Table::Old, Table::Plain);
			CTSummary cts(ct);
			// Log summary of the appropriate type and verbosity.
			switch( CTsummaryType ) {
				case PMS::S_ALL_CT:          cts.list(log, vb); break;
				case PMS::S_WHERE_CT:        cts.listWhere(log, vb); break;
				case PMS::S_WHAT_CT:         cts.listWhat(log, vb); break;
				case PMS::S_HOW_CT:          cts.listHow(log, vb); break;
				case PMS::S_MAIN_CT:         cts.listMain(log, vb); break;
				case PMS::S_TABLES_CT:       cts.listTables(log, vb); break;
				case PMS::S_ANTENNA_CT:      cts.listAntenna(log, vb); break;
				case PMS::S_FIELD_CT:        cts.listField(log, vb); break;
				case PMS::S_OBSERVATION_CT:  cts.listObservation(log, vb); break;
				case PMS::S_HISTORY_CT:      cts.listHistory(log); break;
				case PMS::S_SPW_CT:          cts.listSpectralWindow(log, vb); break;
				default: 					 cts.listMain(log, vb); break;
			}
			success = true;
		} else if (tab.keywordSet().isDefined("CAL_DESC")) {
			CalTable ct = CalTable(filename, Table::Old);
			CalSummary calsum(ct);
			// Log summary of the appropriate type and verbosity.
			switch( CTsummaryType ) {
				case PMS::S_ALL_CT:          calsum.list(log, vb); break;
				case PMS::S_WHERE_CT:        calsum.listWhere(log, vb); break;
				case PMS::S_WHAT_CT:         calsum.listWhat(log, vb); break;
				case PMS::S_HOW_CT:          calsum.listHow(log, vb); break;
				case PMS::S_MAIN_CT:         calsum.listMain(log, vb); break;
				case PMS::S_TABLES_CT:       calsum.listTables(log, vb); break;
				case PMS::S_ANTENNA_CT:      calsum.listAntenna(log, vb); break;
				case PMS::S_FIELD_CT:        calsum.listField(log, vb); break;
				case PMS::S_OBSERVATION_CT:  calsum.listObservation(log, vb); break;
				case PMS::S_HISTORY_CT:      calsum.listHistory(log); break;
				case PMS::S_SPW_CT:          calsum.listSpectralWindow(log, vb); break;
				default: 					 calsum.listMain(log, vb); break;
			}
			success = true;
		} else {
			MeasurementSet ms = MeasurementSet(filename, 
				TableLock(TableLock::AutoLocking), Table::Old);
			MSSummary mss(ms);
			// Log summary of the appropriate type and verbosity.
			switch( summaryType ) {
				case PMS::S_ALL:          mss.list(log, vb); break;
				case PMS::S_WHERE:        mss.listWhere(log, vb); break;
				case PMS::S_WHAT:         mss.listWhat(log, vb); break;
				case PMS::S_HOW:          mss.listHow(log, vb); break;
				case PMS::S_MAIN:         mss.listMain(log, vb); break;
				case PMS::S_TABLES:       mss.listTables(log, vb); break;
				case PMS::S_ANTENNA:      mss.listAntenna(log, vb); break;
				case PMS::S_FEED:         mss.listFeed(log, vb); break;
				case PMS::S_FIELD:        mss.listField(log, vb); break;
				case PMS::S_OBSERVATION:  mss.listObservation(log, vb); break;
				case PMS::S_HISTORY:      mss.listHistory(log); break;
				case PMS::S_POLARIZATION: mss.listPolarization(log, vb); break;
				case PMS::S_SOURCE:       mss.listSource(log, vb); break;
				case PMS::S_SPW:          mss.listSpectralWindow(log, vb); break;
				case PMS::S_SPW_POL:      mss.listSpectralAndPolInfo(log,vb);break;
				case PMS::S_SYSCAL:       mss.listSysCal(log, vb); break;
				case PMS::S_WEATHER:      mss.listWeather(log, vb); break;
			}
			success = true;
		}
	}
	catch(AipsError& x) {
		qDebug()<<"ActionSummary AipsError: "<<x.getMesg().c_str();
		itsDoActionResult_ = x.getMesg();
	}
	// Cleanup.
	if(reenableGlobal) PlotLogger::enableGlobalSink();
	return success;
}

ActionSummary::~ActionSummary() {
}

using namespace casacore;
} /* namespace casa */
