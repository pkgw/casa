//# MSCache.cc: Specialized PlotMSCache for filling MSs
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
//# $Id: $
#include <plotms/Data/MSCache.h>
#include <plotms/Data/PlotMSIndexer.h>

#include <casa/OS/Timer.h>
#include <casa/OS/Memory.h>
#include <casa/Quanta/MVTime.h>
#include <casa/System/Aipsrc.h>
#include <casa/Utilities/Sort.h>
#include <casa/Arrays/ArrayMath.h>
#include <tables/Tables/ScalarColumn.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>
#include <ms/MeasurementSets/MSColumns.h> 	 
#include <msvis/MSVis/VisSet.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <mstransform/MSTransform/MSTransformIteratorFactory.h>
#include <plotms/Data/PlotMSVBAverager.h>
#include <plotms/Data/MSCacheVolMeter.h>
#include <plotms/PlotMS/PlotMS.h>
#include <tables/Tables/Table.h>
#include <measures/Measures/Stokes.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MFrequency.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MSAntennaColumns.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogFilter.h>

#include <ctime>

using namespace casacore;
namespace casa {

MSCache::MSCache(PlotMSApp* parent):
		  PlotMSCacheBase(parent)
{
	ephemerisAvailable = false;
	vi_p = NULL;
	vm_ = NULL;
}

MSCache::~MSCache() {}

String MSCache::polname(Int ipol) {
	return Stokes::name(Stokes::type(ipol));
}


void MSCache::loadIt(vector<PMS::Axis>& loadAxes,
		vector<PMS::DataColumn>& loadData,
		ThreadCommunication* thread) {

	// process selected columns			
	dataColumn_ = getDataColumn(loadAxes, loadData);

    // Get strings stored in MS 
    Table::TableOption tabopt(Table::Old);
    MeasurementSet* inputMS = new MeasurementSet(filename_, TableLock(TableLock::AutoLocking), tabopt);
    getNamesFromMS(*inputMS);

    // Apply selections to MS to create selection MS and channel/correlation selections
    Vector<Vector<Slice> > chansel;
    Vector<Vector<Slice> > corrsel;
    MeasurementSet* selMS = new MeasurementSet();
    try {
        // get chansel, corrsel
        selection_.apply(*inputMS, *selMS, chansel, corrsel);
    } catch(AipsError& log) {
        // improper selection can cause exception
        delete inputMS;
        delete selMS;
        loadError(log.getMesg());
    }

    // Make volume meter for countChunks to estimate memory requirements
    vm_ = new MSCacheVolMeter(*inputMS, averaging_, chansel, corrsel);
    delete inputMS;
    delete selMS;
    Vector<Int> nIterPerAve;

    // make sure user set avgtime value, else get AveragingTVI error
    if (averaging_.time() && averaging_.timeValue()==0.0) {
        averaging_.setTime(false);
        logWarn(PMS::LOG_ORIGIN_LOAD_CACHE, "Time averaging disabled: value is zero."); 
    }
    // only use scalarAve if other averaging enabled
    bool useScalarAve = averaging_.scalarAve() && (averaging_.time() ||
        averaging_.baseline() || averaging_.antenna() ||  averaging_.spw());
    if (averaging_.scalarAve() && !useScalarAve)
        logWarn(PMS::LOG_ORIGIN_LOAD_CACHE, "Scalar averaging ignored: no other averaging is enabled.");
    averaging_.setScalarAve(useScalarAve);

	if ( averaging_.baseline() || averaging_.antenna() || useScalarAve) {
        // Averaging with PlotMSVBAverager
        // Create visibility iterator vi_p
        setUpVisIter(selection_, calibration_, dataColumn_, 
            loadAxes, loadData);
        // Set nIterPerAve (number of chunks per average)
		bool chunksCounted = countChunks(*vi_p, nIterPerAve, loadAxes, 
            loadData, thread);
        if (chunksCounted) {
            try {
                trapExcessVolume(pendingLoadAxes_);  // check mem req using VolMeter
                deleteVm();
                loadChunks(*vi_p, averaging_, nIterPerAve, loadAxes, loadData, thread);
            } catch(AipsError& log) {
                deleteVm();
                loadError(log.getMesg());	
            }
        }    
	} else {
        // Averaging with TransformingVI2 
		try {
			// setUpVisIter also gets the VB shapes and calls trapExcessVolume:
			setUpVisIter(selection_, calibration_, dataColumn_, 
                loadAxes, loadData, false, true, thread);
			loadChunks(*vi_p, loadAxes, loadData, thread);
		} catch(AipsError& log) {
			loadError(log.getMesg());
		}	
	}
	// Remember # of VBs per Average
	nVBPerAve_.resize();
	if (nIterPerAve.nelements()>0) {
		nVBPerAve_ = nIterPerAve;
    } else {
		nVBPerAve_.resize(nChunk_);
		nVBPerAve_.set(1);
	}
	deleteVi(); // close any open tables
}

void MSCache::loadError(String mesg) {
	// catch load error, clear the existing cache, and rethrow
	logLoad(mesg);
	clear();
	deleteVi();  // close any open tables
	stringstream ss;
	ss << mesg;
	throw(AipsError(ss.str()));
}

void MSCache::deleteVi() {
	if (vi_p) delete vi_p;
	vi_p = NULL;
}

void MSCache::deleteVm() {
	if (vm_) delete vm_;
	vm_ = NULL;
}

String MSCache::getDataColumn(vector<PMS::Axis>& loadAxes,
                              vector<PMS::DataColumn>& loadData)
{	// Check data column choice and determine which column to pass to VisIter
	String dataColumn = "NONE";  // default is none - CAS-7506
    std::set<String> dataCols;
    
    // Get datacolumn for visibility & weight axes only
	for (uInt i=0; i<loadAxes.size(); ++i) {
        PMS::Axis thisAxis = loadAxes[i];
        if (PMS::axisIsData(thisAxis) || PMS::axisIsWeight(thisAxis)) {
            // check if requested data column exists in MS
            PMS::DataColumn adjustedCol = checkReqDataColumn(loadData[i]);
            if (adjustedCol != loadData[i])
                adjustCurrentAxes(thisAxis, loadData[i], adjustedCol);
            loadData[i] = adjustedCol;
            
            // check if axis/datacol combo valid
            switch (thisAxis) {
                case PMS::AMP:
                case PMS::REAL:
                case PMS::WTxAMP:
                case PMS::SIGMA:
                case PMS::SIGMASP: {
                    dataColumn = PMS::dataColumn(loadData[i]);
                    break;
                }
                case PMS::PHASE:
                case PMS::IMAG: {
                    // These axes not valid with float data
                    if (loadData[i] == PMS::FLOAT_DATA) {
                        throw(AipsError("Chosen axis not valid for FLOAT_DATA, please use AMP or change Data Column"));
                    } else {
                        dataColumn = PMS::dataColumn(loadData[i]);
                    }
                    break;
                }
                case PMS::WTSP: {
                    // CAS-7517 wtsp col exists but is empty - plot weight instead
	                Table thisTable(filename_);
	                const ColumnDescSet cds = thisTable.tableDesc().columnDescSet();
	                if (cds.isDefined("WEIGHT_SPECTRUM")) {
                        ArrayColumn<Float> weightSpectrum;
                        weightSpectrum.attach(thisTable,
                            MS::columnName(MS::WEIGHT_SPECTRUM));
                        if (!weightSpectrum.hasContent()) {
                            logWarn("load_cache", "Plotting WEIGHT column, WEIGHT_SPECTRUM (WTSP) has not been initialized (this can be changed with initweights task)");
                            cout << "WARNING: Plotting WEIGHT column, WEIGHT_SPECTRUM (WTSP) has not been initialized (this can be changed with initweights task)" << endl;
                        }
                    }
                    dataColumn = PMS::dataColumn(loadData[i]);
                    break;
                }
                case PMS::WT: {
                    // CAS-8895 warn user requested mean WEIGHT not channelized WEIGHT_SPECTRUM
	                Table thisTable(filename_);
	                const ColumnDescSet cds = thisTable.tableDesc().columnDescSet();
	                if (cds.isDefined("WEIGHT_SPECTRUM")) {
                        ArrayColumn<Float> weightSpectrum;
                        weightSpectrum.attach(thisTable,
                            MS::columnName(MS::WEIGHT_SPECTRUM));
                        if (weightSpectrum.hasContent()) {
                            logLoad("Plotting mean WEIGHT column but WEIGHT_SPECTRUM is available. Request 'WtSp' axis to see channelized weights.");
                            cout << "INFO: Plotting mean WEIGHT column but WEIGHT_SPECTRUM is available. Request 'WtSp' axis to see channelized weights." << endl;
                        }
                    }
                    dataColumn = PMS::dataColumn(loadData[i]);
                }
                default:
                    break;
            }
            dataCols.insert(dataColumn);
        }
	} 

    // loadAxes is only new axes to load; if no datacolumn,
    // try datacolumn of already-loaded axes
    if (dataColumn == "NONE") {
        dataColumn = checkLoadedAxesDatacol();
        // might still be "NONE"!
    }
    // convert to mstransform datacolumn string
    if (dataColumn != "NONE") {
        if (dataCols.size() > 1)
            dataColumn = "ALL";
        else
            dataColumn = normalizeColumnName(dataColumn);
    }
	return dataColumn;
}

PMS::DataColumn MSCache::checkReqDataColumn(PMS::DataColumn reqDataCol) {
	// Check if requested data, scratch, or float cols exist
    PMS::DataColumn datacol = reqDataCol;
	Bool datacolOk(false), corcolOk(false), floatcolOk(false);
	Table thisTable(filename_);
	const ColumnDescSet cds = thisTable.tableDesc().columnDescSet();
	datacolOk  = cds.isDefined("DATA");
	corcolOk   = cds.isDefined("CORRECTED_DATA");
	floatcolOk = cds.isDefined("FLOAT_DATA");

    switch (reqDataCol) {
        case PMS::DATA: {
            // CAS-7482 - for singledish, use FLOAT if no DATA 
            if (!datacolOk && floatcolOk) {
                datacol = PMS::FLOAT_DATA;
                logWarn( "load_cache", "DATA column not present; will use FLOAT_DATA instead.");
            }
            break;
        }
        case PMS::CORRECTED:
        case PMS::CORRECTED_DIVIDE_MODEL:
        case PMS::CORRMODEL: {
            // requested corrected data but no (real or OTF) corrected column
            if (!corcolOk && !calibration_.useCallib()) {
                if (datacolOk) {
                    // CAS-5214 - use DATA if no CORRECTED_DATA with warning
                    datacol = PMS::DATA;
                    logWarn( "load_cache", "CORRECTED_DATA column not present and calibration library not set or enabled; will use DATA instead.");
                } else if (floatcolOk) {
                    // CAS-7761 - for singledish, use FLOAT if no CORRECTED or DATA
                    datacol = PMS::FLOAT_DATA;
                    logWarn( "load_cache", "CORRECTED_DATA column not present and calibration library not set or enabled; will use FLOAT_DATA instead.");
                }
            }
            break;
        }
        case PMS::FLOAT_DATA: {
            // requested float data but no FLOAT column 
            if (!floatcolOk) {
                throw(AipsError("FLOAT_DATA not present, please use DATA"));
            }
            break;
        }
        default:
            break;
    } // switch
    return datacol;
}

void MSCache::adjustCurrentAxes(PMS::Axis axis,
        PMS::DataColumn olddata, PMS::DataColumn newdata) {
    // adjust data column for currentX or currentY
    for (uInt a=0; a<currentX_.size(); ++a) {
        if (currentX_[a] == axis && currentXData_[a]==olddata) {
            currentXData_[a] = newdata;
        }
        if (currentY_[a] == axis && currentYData_[a]==olddata) {
            currentYData_[a] = newdata;
        }
    }
}

String MSCache::checkLoadedAxesDatacol() {
    // Check data column of plotted axes
    String loadedCol = "NONE";
    std::set<String> loadedColumns;
    int axesCount = currentX_.size();
    for (int i=0; i<axesCount; ++i) {
        if (PMS::axisIsData(currentX_[i]) || PMS::axisIsWeight(currentX_[i])) 
            loadedColumns.insert(PMS::dataColumn(currentXData_[i]));
        if (PMS::axisIsData(currentY_[i]) || PMS::axisIsWeight(currentY_[i]))
            loadedColumns.insert(PMS::dataColumn(currentYData_[i]));
    }
    if (loadedColumns.size() == 1)
        loadedCol = *loadedColumns.begin();
    else if (loadedColumns.size() > 1)
        loadedCol = "ALL";
    return loadedCol;
}

String MSCache::normalizeColumnName(String plotmscol)
{
	// Convert datacolumn as needed for MSTransformManager
    String colname = plotmscol;
	if ((plotmscol == "corrected-model") || 
	    (plotmscol == "data-model") ||
	    (plotmscol == "data/model") || 
	    (plotmscol == "corrected/model")) {
			colname = "ALL";
	} else if (plotmscol == "float") {
		colname = "FLOAT_DATA";
	} else {   // "data", "corrected", "model"	
		colname.upcase();
	}
	return colname;
}

void MSCache::getNamesFromMS(MeasurementSet& ms)
{
    ROMSColumns msCol(ms);
    antnames_.resize();
    stanames_.resize();
    antstanames_.resize();
    fldnames_.resize();
    intentnames_.resize();
    positions_.resize();

    antnames_    = msCol.antenna().name().getColumn();
    stanames_    = msCol.antenna().station().getColumn();
    antstanames_ = antnames_+String("@")+stanames_;
    positions_   = msCol.antenna().position().getColumn();

    fldnames_    = msCol.field().name().getColumn();

    intentnames_ = msCol.state().obsMode().getColumn();
    mapIntentNamesToIds();  // eliminate duplicate intent names
}

void MSCache::setUpVisIter(PlotMSSelection& selection,
		PlotMSCalibration& calibration,
		String dataColumn, 
        vector<PMS::Axis>& loadAxes,
		vector<PMS::DataColumn>& loadData, 
        Bool interactive,
        Bool estimateMemory,
        ThreadCommunication* thread) {
	/* Create plain or averaging (time or channel) VI with 
           configuration Record and MSTransformIterator factory */

	// Create configuration:
	// Start with data selection; rename fields with expected keywords
	Record configuration = selection.toRecord();
	configuration.renameField("correlation", configuration.fieldNumber("corr"));

	// Add needed fields
	configuration.define("inputms", filename_);
	configuration.define("datacolumn", dataColumn);
	configuration.define("buffermode", true);
	configuration.define("reindex", false);
    configuration.define("interactive", interactive);

	// Add transformation selection with expected keywords and string value
	configuration.merge(transformations_.toRecord());
	Double restfreq = configuration.asDouble(configuration.fieldNumber("RestFreq"));
	configuration.removeField(configuration.fieldNumber("RestFreq"));
	configuration.define("restfreq", String::toString(restfreq));
	configuration.renameField("outframe", configuration.fieldNumber("Frame"));
	configuration.renameField("veltype", configuration.fieldNumber("Veldef"));

	// Add calibration library if set
	if (calibration.useCallib()) {
		configuration.define("callib", calibration.calLibrary());
	}	

	// Apply averaging
	if (averaging_.time()){
        if (averaging_.scalarAve()) {
		    configuration.define("scalaraverage", true);
        } else {
		    configuration.define("timeaverage", true);
		    String timespanStr = "state";
		    if (averaging_.field())
			    timespanStr += ",scan,field";
		    else if (averaging_.scan())
			    timespanStr += ",scan";
		    configuration.define("timespan", timespanStr);
        }
		configuration.define("timebin", averaging_.timeStr());
	}
	if (averaging_.channel()) {
        int chanBin;
        double chanVal = averaging_.channelValue();
        if (chanVal > INT_MAX) {
            chanBin = INT_MAX;
            logWarn(PMS::LOG_ORIGIN_LOAD_CACHE, "avgchannel value exceeds maximum integer allowed (" + String::toString(INT_MAX) + ")");
        } else if (chanVal == 0.0) {
            chanBin = 1;
            logLoad("Cannot average 0 channels, using 1 instead (no averaging).");
        } else {
            chanBin = static_cast<int>(chanVal);
        }
		configuration.define("chanaverage", true);
		configuration.define("chanbin", chanBin);
	}
    if (averaging_.spw()) {
        configuration.define("spwaverage", true);
    }

    LogFilter oldFilter(plotms_->getParameters().logPriority());
	MSTransformIteratorFactory* factory = NULL;
	try {
        // Filter out MSTransformManager setup messages
        LogFilter filter(LogMessage::WARN);
        LogSink().globalSink().filter(filter);
		factory = new MSTransformIteratorFactory(configuration);
		if (estimateMemory) {
            if (thread != NULL)
                updateEstimateProgress(thread);
			visBufferShapes_ = factory->getVisBufferStructure();
			Int chunks = visBufferShapes_.size();
			setCache(chunks, loadAxes, loadData);
			trapExcessVolume(pendingLoadAxes_);
		} else {
            visBufferShapes_.clear();
        }
		vi_p = new vi::VisibilityIterator2(*factory);
	} catch(AipsError& log) {
        // now put filter back
        LogSink().globalSink().filter(oldFilter);
		try {
			if (factory) delete factory;
		} catch(AipsError ae) {}
		throw(AipsError(log.getMesg()));
	}
    // now put filter back
    LogSink().globalSink().filter(oldFilter);
	if (factory) delete factory;
}

vi::VisibilityIterator2* MSCache::setUpVisIter(MeasurementSet& selectedMS,
	Vector<Vector<Slice> > chansel, Vector<Vector<Slice> > corrsel) {
	// Plain VI2 for chunk counting with baseline/antenna/spw averaging
	// (need to set SortColumns manually)
	Bool combscan(averaging_.scan());
	Bool combfld(averaging_.field());
	Bool combspw(averaging_.spw());

	// Set iterInterval
	Double iterInterval(0.0);
	if (averaging_.time()){
		iterInterval = averaging_.timeValue();
	}
	if (combspw || combfld) iterInterval = DBL_MIN;  // force per-timestamp chunks

	// Create SortColumns
	Int nsortcol(4 + Int(!combscan));  // include room for scan
	Block<Int> columns(nsortcol);
	Int i(0);
	columns[i++]                = MS::ARRAY_ID;
	if (!combscan) columns[i++] = MS::SCAN_NUMBER;  // force scan boundaries
	if (!combfld) columns[i++]  = MS::FIELD_ID;      // force field boundaries
	if (!combspw) columns[i++]  = MS::DATA_DESC_ID;  // force spw boundaries
	columns[i++]                = MS::TIME;
	if (combfld) columns[i++]   = MS::FIELD_ID;      // effectively ignore field boundaries
	if (combspw) columns[i++]   = MS::DATA_DESC_ID;  // effectively ignore spw boundaries
	vi::SortColumns sortcol(columns, false);

	vi::VisibilityIterator2* vi2 = new vi::VisibilityIterator2(selectedMS, sortcol, false, 0, iterInterval);
	vi::FrequencySelectionUsingChannels fs;
	setUpFrequencySelectionChannels(fs, chansel);
	fs.addCorrelationSlices(corrsel);
	// Add FrequencySelection to VI
	vi2->setFrequencySelection(fs);
	return vi2;
}

void MSCache::setUpFrequencySelectionChannels(vi::FrequencySelectionUsingChannels fs,
						Vector<Vector<Slice> > chansel) {
	/* For the plain VI2 for chunk counting */
	int nSpws = static_cast<int>(chansel.nelements());
	int nChansels;
	for (int spw=0; spw<nSpws; spw++) {
		nChansels = static_cast<int>(chansel[spw].nelements());
		// Add channel selections to FrequencySelection
		for ( int sel=0; sel<nChansels; sel++) {
			fs.add(spw, 
			       chansel[spw][sel].start(), 
			       chansel[spw][sel].length(),
			       chansel[spw][sel].inc());
		}
	}
}

void MSCache::updateEstimateProgress(ThreadCommunication* thread) {
    thread->setStatus("Establishing cache size.  Please wait...");
    thread->setAllowedOperations(false,false,true);
    thread->setProgress(2);
}

bool MSCache::countChunks(vi::VisibilityIterator2& vi,
        Vector<Int>& nIterPerAve,
        vector<PMS::Axis>& loadAxes,
		vector<PMS::DataColumn>& loadData, 
        ThreadCommunication* thread) {
    // Let plotms count the chunks for memory estimation 
    //   when baseline/antenna/spw/scalar averaging
    if (thread != NULL)
        updateEstimateProgress(thread);

    Bool verby(False);
    stringstream ss;

    Bool combscan(averaging_.scan());
    Bool combfld(averaging_.field());
    Bool combspw(averaging_.spw());

    vi::VisBuffer2* vb = vi.getVisBuffer();
    vi.originChunks();
    vi.origin();

    // Keeping time
    Double time1(0.0), avetime1(-1.0);
    Double interval(0.0);
    if (averaging_.time())
        interval = averaging_.timeValue();
    // Keep track of other boundaries
    Int thisscan(-1),lastscan(-1);
    Int thisfld(-1), lastfld(-1);
    Int thisspw(-1),lastspw(-1);
    Int thisddid(-1),lastddid(-1);
    Int thisobsid(-1),lastobsid(-1);
    // Averaging stats
    Int chunk(0), subchunk(0);
    Int maxAveNRows(0);
    nIterPerAve.resize(100);
    nIterPerAve = 0;
    Int nAveInterval(-1);

    for (vi.originChunks(); vi.moreChunks(); vi.nextChunk(), chunk++) {
        subchunk = 0;
        for (vi.origin(); vi.more(); vi.next(), subchunk++) {
            // If a thread is given, check if the user canceled.
            if (thread != NULL) {
                if (thread->wasCanceled()) {
                    dataLoaded_ = false;
                    userCanceled_ = true;
                    return false;
                } else {
                    // else users think it's hung...
                    if ((chunk % 100) == 0)
                        thread->setProgress(chunk/100);
                }
            }

            time1 = vb->time()(0); // first timestamp in this vb
            thisscan = vb->scan()(0);
            thisfld = vb->fieldId()(0);
            thisspw = vb->spectralWindows()(0);
            thisddid = vb->dataDescriptionIds()(0);
            thisobsid = vb->observationId()(0);

            // New ave interval if:
            if ( ((time1-avetime1) > interval) ||          // exceeded time interval
                 ((time1-avetime1) < 0.0) ||               // negative timestep
                 (!combscan && (thisscan != lastscan)) ||  // not combing scans, and new scan encountered OR
                 (!combfld && (thisfld != lastfld)) ||     // not combing fields, and new field encountered OR
                 (!combspw && (thisspw != lastspw)) ||     // not combing spws, and new spw encountered  OR
                 (thisobsid != lastobsid) ||               // don't average over obs id
                 (nAveInterval == -1)) {                   // this is the first interval

                if (verby) {
                    ss << "--------------------------------\n";
                    ss << boolalpha << interval << " "
                       << ((time1 - avetime1)>interval) << " "
                       << ((time1 - avetime1)<0.0) << " "
                       << (!combscan && (thisscan!=lastscan)) << " "
                       << (!combspw && (thisspw!=lastspw)) << " "
                       << (!combfld && (thisfld!=lastfld)) << " "
                       << (!combspw && (thisspw!=lastspw)) << " "
                       << (thisobsid!=lastobsid) << " "
                       << (nAveInterval == -1) << "\n";
                }

                // If we have accumulated enough info, poke the volume meter,
                //  with the _previous_ info, and reset the ave'd row counter
                if (nAveInterval > -1) {
                    vm_->add(lastddid, maxAveNRows);
                    maxAveNRows = 0;
                }

                nAveInterval++;
                if (verby) ss << "ave = " << nAveInterval << "\n";

                // increase size of nIterPerAve array, if needed
                if (nIterPerAve.nelements() < uInt(nAveInterval+1))
                    nIterPerAve.resize(nIterPerAve.nelements()+100, true);
                // initialize next ave interval
                nIterPerAve(nAveInterval) = 0;
                avetime1 = time1;  // first timestamp in this averaging interval
            }

            // Keep track of the maximum # of rows that might get averaged
            maxAveNRows = max(maxAveNRows, vb->nRows());
            // Increment chunk-per-sol count for current solution
            nIterPerAve(nAveInterval)++;

            if (verby) {
                ss << "     chunk=" << chunk << " subchunk " << subchunk << "\n";
                ss << "         time=" << vb->time()(0) << " ";
                ss << "arrayId=" << vb->arrayId()(0) << " ";
                ss << "scan" << thisscan << " ";
                ss << "fieldId=" << thisfld << " ";
                ss << "spw=" << thisspw << " ";
                ss << "obsId=" << thisobsid << "\n";
            }

            lastscan = thisscan;
            lastfld  = thisfld;
            lastspw  = thisspw;
            lastddid = thisddid;
            lastobsid = thisobsid;
        }
    }
    // Add in the last iteration
    vm_->add(lastddid,maxAveNRows);

    Int nAve(nAveInterval+1);
    nIterPerAve.resize(nAve, True);
    setCache(nAve, loadAxes, loadData);  // sets nChunk_

    if (verby) {
        ss << "nIterPerAve = " << nIterPerAve << "\n";
        ss << "Found " << nChunk_ << " chunks." << endl;
        logInfo("count_chunks", ss.str());
    }
    return true;
}

void MSCache::trapExcessVolume(map<PMS::Axis,Bool> pendingLoadAxes) {
	try {
		String s;
		if (visBufferShapes_.size() > 0) {
			s = vm_->evalVolume(visBufferShapes_, pendingLoadAxes); }
		else {
			Vector<Bool> mask(4, false);
			int dataCount = getDataCount();

			for ( int i = 0; i < dataCount; i++ ){
				Vector<Bool> subMask = netAxesMask( currentX_[i], currentY_[i]);
				mask = mask || subMask;
			}
			s = vm_->evalVolume(pendingLoadAxes, mask);
		}
		logLoad(s);

	} catch(AipsError& log) {
		// catch detected volume excess, clear the existing cache, and rethrow
		logLoad(log.getMesg());
		deleteVm();
		stringstream ss;
		ss << "Please try selecting less data or averaging and/or" << endl
		   << " 'force reload' (to clear unneeded cache items) and/or" << endl
		   << " letting other memory-intensive processes finish.";
		throw(AipsError(ss.str()));
	}
}

void MSCache::updateProgress(ThreadCommunication* thread, Int chunk) {
	double progress;
	// Update progress meter
	if((nChunk_ <= (int)THREAD_SEGMENT || chunk % THREAD_SEGMENT == 0)) {
		thread->setStatus("Loading chunk " + String::toString(chunk) +
				" / " + String::toString(nChunk_) + ".");
		progress = ((double)chunk+1) / nChunk_;
		thread->setProgress((unsigned int)((progress * 100) + 0.5));
    }
}

void MSCache::loadChunks(vi::VisibilityIterator2& vi,
		const vector<PMS::Axis> loadAxes,
		const vector<PMS::DataColumn> loadData,
		ThreadCommunication* thread) {

	// permit cancel in progress meter:
	if(thread != NULL)
		thread->setAllowedOperations(false,false,true);
	logLoad("Loading chunks......");

	// Initialize VI and get VB
	vi::VisBuffer2* vb = vi.getVisBuffer();

	vi.originChunks();
    try {
	    vi.origin();
    } catch (ArraySlicerError & err) {
	    throw(AipsError("PlotMS averaging error: Data shapes do not conform."));
    } catch (ArrayConformanceError & err) {
	    throw(AipsError("PlotMS averaging error: Data shapes do not conform."));
    }

	nAnt_ = vb->nAntennas();  // needed to set up indexer
	// set frame; VB2 does not handle N_Types, just passes it along
	// and fails check in MFrequency so handle it here
	freqFrame_ = transformations_.frame();
	if (freqFrame_ == MFrequency::N_Types)
		freqFrame_ = static_cast<MFrequency::Types>(vi.getReportingFrameOfReference());

	Int chunk = 0;
	chshapes_.resize(4,nChunk_);
	goodChunk_.resize(nChunk_);
	goodChunk_.set(false);

	for(vi.originChunks(); vi.moreChunks(); vi.nextChunk()) {
		for(vi.origin(); vi.more(); vi.next()) {
            if (vb->nRows() > 0) {
                if (chunk >= nChunk_) {  // nChunk_ was just an estimate
                    setCache(chunk+1, loadAxes, loadData);
                    chshapes_.resize(4, nChunk_, true);
                    goodChunk_.resize(nChunk_, true);
                }

                // If a thread is given, update its chunk number and progress bar
                if(thread != NULL)
                    updateProgress(thread, chunk);

                // Cache the data shapes
                chshapes_(0,chunk) = vb->nCorrelations();
                chshapes_(1,chunk) = vb->nChannels();
                chshapes_(2,chunk) = vb->nRows();
                chshapes_(3,chunk) = vb->nAntennas();
                goodChunk_(chunk)  = true;
                for(unsigned int i = 0; i < loadAxes.size(); i++) {
                    // If a thread is given, check if the user canceled.
                    if(thread != NULL && thread->wasCanceled()) {
                        dataLoaded_ = false;
                        userCanceled_ = true;
                        goodChunk_(chunk) = false; //only partially loaded
                        return;
                    }
                    loadAxis(vb, chunk, loadAxes[i], loadData[i]);
                }

                chunk++;
            }
		}
	}
    // Report averaged channels per spw in log
    // (should match MSTransformManager output in console)
    map<Int,Int>::iterator it;
    for (it=chansPerSpw_.begin(); it!=chansPerSpw_.end(); ++it) {
        logLoad("SPW " + String::toString(it->first) + ": number of channels averaged = " + String::toString(it->second));
    }
}

void MSCache::loadChunks(vi::VisibilityIterator2& vi,
        const PlotMSAveraging& averaging,
        const Vector<Int>& nIterPerAve,
        const vector<PMS::Axis> loadAxes,
        const vector<PMS::DataColumn> loadData,
        ThreadCommunication* thread) {

    // permit cancel in progress meter:
    if(thread != NULL)
        thread->setAllowedOperations(false,false,true);
    logLoad("Loading chunks with averaging.....");

    Bool verby(false);

    vi::VisBuffer2* vb = vi.getVisBuffer();
    vi.originChunks();
    vi.origin();
    nAnt_ = vb->nAntennas();  // needed to set up indexer

    // set frame; VB2 does not handle N_Types, just passes it along
    // and fails check in MFrequency so handle it here
    freqFrame_ = transformations_.frame();
    if (freqFrame_ == MFrequency::N_Types)
        freqFrame_ = static_cast<MFrequency::Types>(vi.getReportingFrameOfReference());

    chshapes_.resize(4, nChunk_);
    goodChunk_.resize(nChunk_);
    goodChunk_.set(false);
    Int nAnts;
    vi::VisBuffer2* vbToUse = NULL;

    for (Int chunk=0; chunk<nChunk_; ++chunk) {
        if (chunk >= nChunk_) {  // nChunk_ was just an estimate!
            setCache(chunk, loadAxes, loadData);
            chshapes_.resize(4, nChunk_, true);
            goodChunk_.resize(nChunk_, true);
        }

        // If a thread is given, update it.
        if(thread != NULL)
            updateProgress(thread, chunk);

        // Set up VB averager
        nAnts = vb->nAntennas();
        PlotMSVBAverager pmsvba(nAnts);
        pmsvba.setBlnAveraging(averaging.baseline());
        pmsvba.setAntAveraging(averaging.antenna());
        pmsvba.setScalarAve(averaging.scalarAve());
        // Sort out which data to read
        discernData(loadAxes,loadData,pmsvba);

        stringstream ss;
        if (verby) ss << "Chunk " << chunk << " ----------------------------------\n";

        Int iter=0;
        while (iter < nIterPerAve(chunk)) {
            if (verby) {
                ss << "chunk=" << chunk << " iter=" << iter << " nIterPerAve =" << nIterPerAve(chunk) << " ";
                ss << "scan=" << vb->scan()(0) << " ";
                ss << "fieldId=" << vb->fieldId()(0) << " ";
                ss << "spw=" << vb->spectralWindows()(0) << " ";
                ss << "amp=" << vb->visCube().shape() << " ";
                ss << "freq=" << vb->getFrequencies(0, freqFrame_).shape() << " ";
            }
            // Accumulate into the averager
            pmsvba.accumulate(*vb);

            // Advance to next VB unless you are going to finalize
            if (iter+1 < nIterPerAve(chunk)) {
                if (verby) ss << " next VB "; 
                vi.next(); 
                if (!vi.more() && vi.moreChunks()) { 
                     // go to first vb in next chunk 
                     if (verby) ss << "  stepping VI"; 
                     vi.nextChunk(); 
                     vi.origin(); 
                }
            }
            if (verby) ss << "\n";
            ++iter;
        }

        if (verby) ss << "Finalize average\n";
        // Finalize the averaging
        pmsvba.finalizeAverage();
        // The averaged VisBuffer
        vi::VisBuffer2& avb(pmsvba.aveVisBuff());
        // Only if the average yielded some data:
        if (avb.nRows() > 0) {
            // Cache the data shapes
            chshapes_(0,chunk) = avb.nCorrelations();
            chshapes_(1,chunk) = avb.nChannels();
            chshapes_(2,chunk) = avb.nRows();
            chshapes_(3,chunk) = nAnts;
            goodChunk_(chunk)  = true;

            for(unsigned int i = 0; i < loadAxes.size(); i++) {
                // If a thread is given, check if the user canceled.
                if(thread != NULL && thread->wasCanceled()) {
                    dataLoaded_ = false;
                    userCanceled_ = true;
                    goodChunk_(chunk)  = false;
                    return;
                }
                if (useAveragedVisBuffer(loadAxes[i])) {
                    vbToUse = &avb;
                } else {
                    vbToUse = vb;
                }
                loadAxis(vbToUse, chunk, loadAxes[i], loadData[i]);
            }
        } else {
            // no points in this chunk
            goodChunk_(chunk) = false;
            chshapes_.column(chunk) = 0;
        }
        // Now advance to next chunk
        if (verby) ss << " next VB "; 
        vi.next();
        if (!vi.more() && vi.moreChunks()) { 
             // go to first vb in next chunk 
             if (verby) ss << "  stepping VI"; 
             vi.nextChunk(); 
             vi.origin(); 
        }
        if(verby) {
            ss << "\n";
            logLoad(ss.str());
        }
    }
    //cout << boolalpha << "goodChunk_ = " << goodChunk_ << endl;
}

bool MSCache::useAveragedVisBuffer(PMS::Axis axis) {
	// Some axes should be obtained from the VB2 provided by the VI2
	// rather than the averaged VB2, which is not attached
	bool useAvg(true);
	switch(axis) {
	case PMS::CHANNEL:
	case PMS::FREQUENCY:
	case PMS::VELOCITY:
	case PMS::UVDIST_L:
	case PMS::UWAVE:
	case PMS::VWAVE:
	case PMS::WWAVE:
	case PMS::AZ0:
	case PMS::EL0:
	case PMS::HA0:
	case PMS::PA0:
	case PMS::ANTENNA:
	case PMS::AZIMUTH:
	case PMS::ELEVATION:
	case PMS::PARANG:
	case PMS::ROW: {
		useAvg = false;
		break;
	}
	default:
		break;
	}
	return useAvg;
}

void MSCache::forceVBread(vi::VisBuffer2* vb,
		vector<PMS::Axis> loadAxes,
		vector<PMS::DataColumn> loadData) {

	// pre-load requisite pieces of VisBuffer for averaging
	for(unsigned int i = 0; i < loadAxes.size(); i++) {
		switch (loadAxes[i]) {
		case PMS::AMP:
		case PMS::PHASE:
		case PMS::REAL:
		case PMS::IMAG: {
			switch(loadData[i]) {
			case PMS::DATA: {
				vb->visCube();
				break;
			}
			case PMS::MODEL: {
				vb->visCubeModel();
				break;
			}
			case PMS::CORRECTED: {
				vb->visCubeCorrected();
				break;
			}
			case PMS::CORRECTED_DIVIDE_MODEL:
			case PMS::CORRMODEL: {
				vb->visCubeCorrected();
				vb->visCubeModel();
				break;
			}
			case PMS::DATAMODEL: {
				vb->visCube();
				vb->visCubeModel();
				break;
			}
			case PMS::DATA_DIVIDE_MODEL: {
				vb->visCube();
				vb->visCubeModel();
				break;
			}
			case PMS::FLOAT_DATA: {
				vb->visCubeFloat();
				break;
			}
			default:
				break;
			}
			break;
		}
		default:
			break;
		}
	}

	// Always need flags
	vb->flagRow();
	vb->flagCube();

}

void MSCache::discernData(vector<PMS::Axis> loadAxes,
		vector<PMS::DataColumn> loadData,
		PlotMSVBAverager& vba) {

	// Turn off
	vba.setNoData();

	// Tell the averager which data column to read
	for(unsigned int i = 0; i < loadAxes.size(); i++) {
		switch (loadAxes[i]) {
		case PMS::AMP:
		case PMS::PHASE:
		case PMS::REAL:
		case PMS::IMAG: {
			switch(loadData[i]) {
			case PMS::DATA: {
				vba.setDoVC();
				break;
			}
			case PMS::MODEL: {
				vba.setDoMVC();
				break;
			}
			case PMS::CORRECTED: {
				vba.setDoCVC();
				break;
			}
			case PMS::CORRECTED_DIVIDE_MODEL:
			case PMS::CORRMODEL: {
				vba.setDoCVC();
				vba.setDoMVC();
				break;
			}
			case PMS::DATAMODEL:
			case PMS::DATA_DIVIDE_MODEL: {
				vba.setDoVC();
				vba.setDoMVC();
                break;
			}
			case PMS::FLOAT_DATA:
				vba.setDoFC();
				break;
			default:
				break;
			}
			break;
		}
		case PMS::UVDIST:
		case PMS::UVDIST_L:
		case PMS::U:
		case PMS::V:
		case PMS::W:
		case PMS::UWAVE:
		case PMS::VWAVE:
		case PMS::WWAVE: {
			//  cout << "Arranging to load UVW
			vba.setDoUVW();
            break;
		}
		default:
			break;
		}
	}

}




void MSCache::loadAxis(vi::VisBuffer2* vb, Int vbnum, PMS::Axis axis,
		PMS::DataColumn data) {

	switch(axis) {

	case PMS::SCAN: { // assumes scan unique in VB
		scan_(vbnum) = vb->scan()(0);
		break;
    }

	case PMS::FIELD: { // assumes field unique in VB
		field_(vbnum) = vb->fieldId()(0);
		break;
    }

	case PMS::TIME: { // assumes time unique in VB
		time_(vbnum) = vb->time()(0);
		break;
    }

	case PMS::TIME_INTERVAL: { // assumes timeInterval unique in VB
		timeIntr_(vbnum) = vb->timeInterval()(0);
		break;
    }

	case PMS::SPW: {
		spw_(vbnum) = vb->spectralWindows()(0);
		break;
    }

	case PMS::CHANNEL: {
        Vector<Int> chans = vb->getChannelNumbers(0);
        *chan_[vbnum] = chans;
	    if (averaging_.channel()) {
            // Save which channels are being averaged, for Locate
            Int numBins = chans.size();
            Int numChans = vb->getChannelNumbersSelected(0).size();
            Array<Int> chansPerBin(IPosition(2, numChans, numBins));
            // when there is only one channel in the spw selected chans is empty array
            if (numChans==0 && numBins==1) {
                numChans = 1;
                chansPerBin.resize(IPosition(2, 1, 1));
                chansPerBin[0] = chans;
            } else {
                for (Int bin=0; bin < numBins; ++bin) {
                    chansPerBin[bin] = vb->getChannelNumbersSelected(bin);
                }
            }
            *chansPerBin_[vbnum] = chansPerBin;
            chansPerSpw_[vb->spectralWindows()(0)] = numChans;
        }
		break;
    }
	case PMS::FREQUENCY: {
		// Convert freq to desired frame
  		*freq_[vbnum] = vb->getFrequencies(0, freqFrame_);
		(*freq_[vbnum]) /= 1.0e9; // in GHz
		break;
	}
	case PMS::VELOCITY: {
		*vel_[vbnum] = calcVelocity(vb); 
		break;
	}

	case PMS::CORR: {
        Vector<Stokes::StokesTypes> corrTypes = vb->getCorrelationTypesSelected();
        Vector<Int> corrTypesInt;
        corrTypesInt.resize(corrTypes.size());
        for (uInt i=0; i<corrTypes.size(); ++i) {
            corrTypesInt[i] = static_cast<Int>(corrTypes[i]);
        }
		*corr_[vbnum] = corrTypesInt;
		break;
    }

	case PMS::ANTENNA1: {
		*antenna1_[vbnum] = vb->antenna1();
		break;
	}
	case PMS::ANTENNA2: {
		*antenna2_[vbnum] = vb->antenna2();
		break;
	}
	case PMS::BASELINE: {
		Vector<Int> a1(vb->antenna1());
		Vector<Int> a2(vb->antenna2());
		baseline_[vbnum]->resize(vb->nRows());
		Vector<Int> bl(*baseline_[vbnum]);
		for (Int irow = 0; irow < vb->nRows(); ++irow) {
			if (a1(irow)<0) a1(irow)=chshapes_(3,0);
			if (a2(irow)<0) a2(irow)=chshapes_(3,0);
			bl(irow)=(chshapes_(3,0)+1)*a1(irow) - (a1(irow) * (a1(irow) + 1)) / 2 + a2(irow);
		}
		break;
	}
	case PMS::UVDIST: {
		Array<Double> u(vb->uvw().row(0));
		Array<Double> v(vb->uvw().row(1));
		*uvdist_[vbnum] = sqrt(u*u+v*v);
		break;
	}
	case PMS::U: {
		*u_[vbnum] = vb->uvw().row(0);
		break;
	}
	case PMS::V: {
		*v_[vbnum] = vb->uvw().row(1);
		break;
	}
	case PMS::W: {
		*w_[vbnum] = vb->uvw().row(2);
		break;
	}
	case PMS::UVDIST_L: {
		Array<Double> u(vb->uvw().row(0));
		Array<Double> v(vb->uvw().row(1));
		Vector<Double> uvdistM = sqrt(u*u + v*v);
		uvdistM /=C::c;
		uvdistL_[vbnum]->resize(vb->nChannels(), vb->nRows());
		Vector<Double> uvrow;
		for (Int irow = 0; irow < vb->nRows(); ++irow) {
			uvrow.reference(uvdistL_[vbnum]->column(irow));
			uvrow.set(uvdistM(irow));
			uvrow *= vb->getFrequencies(irow, freqFrame_);
		}
		break;
	}

	case PMS::UWAVE: {
		Vector<Double> uM(vb->uvw().row(0));
		uM/=C::c;
		uwave_[vbnum]->resize(vb->nChannels(), vb->nRows());
		Vector<Double> urow;
		for (Int irow = 0; irow < vb->nRows(); ++irow) {
			urow.reference(uwave_[vbnum]->column(irow));
			urow.set(uM(irow));
			urow *= vb->getFrequencies(irow, freqFrame_);
		}
		break;
	}
	case PMS::VWAVE: {
		Vector<Double> vM(vb->uvw().row(1));
		vM/=C::c;
		vwave_[vbnum]->resize(vb->nChannels(), vb->nRows());
		Vector<Double> vrow;
		for (Int irow = 0; irow < vb->nRows(); ++irow) {
			vrow.reference(vwave_[vbnum]->column(irow));
			vrow.set(vM(irow));
			vrow *= vb->getFrequencies(irow, freqFrame_);
		}
		break;
	}
	case PMS::WWAVE: {
		Vector<Double> wM(vb->uvw().row(2));
		wM/=C::c;
		wwave_[vbnum]->resize(vb->nChannels(), vb->nRows());
		Vector<Double> wrow;
		for (Int irow = 0; irow < vb->nRows(); ++irow) {
			wrow.reference(wwave_[vbnum]->column(irow));
			wrow.set(wM(irow));
			wrow *= vb->getFrequencies(irow, freqFrame_);
		}
		break;
	}
	case PMS::AMP: {
		switch(data) {
		case PMS::DATA: {
			//CAS-5730.  For single dish data, absolute value of
			//points should not be plotted.
			MeasurementSet ms( filename_);
			if ( ms.isColumn( MS::FLOAT_DATA ) || averaging_.scalarAve()){
				*amp_[vbnum]=real(vb->visCube());
			}
			else {
				*amp_[vbnum] = amplitude(vb->visCube());
			}
			// TEST fft on freq axis to get delay
			if (false) {

				// Only transform frequency axis
				//   (Should avoid cross-hand data, too?)
				Vector<Bool> ax(3,false);
				ax(1) = true;

				// Support padding for higher delay resolution
				Int fact(4);
				IPosition ip = vb->visCube().shape();
				Int nch = ip(1);
				ip(1) *= fact;

				Slicer sl(Slice(),Slice(nch*(fact-1)/2,nch,1),Slice());

				Array<Complex> vpad(ip);
				vpad.set(Complex(0.0));
				vpad(sl) = vb->visCube();


				cout << "vpad.shape() = " << vpad.shape() << endl;
				cout << "vpad(sl).shape() = " << vpad(sl).shape() << endl;

				Vector<Complex> testf(64,Complex(1.0));
				FFTServer<Float,Complex> ffts;
				cout << "FFTServer..." << flush;
				ffts.fft(testf,true);
				cout << "done." << endl;

				ArrayLattice<Complex> tf(testf);
				cout << "tf.isWritable() = " << boolalpha << tf.isWritable() << endl;

				LatticeFFT::cfft(tf,false);
				cout << "testf = " << testf << endl;


				cout << "Starting ffts..." << flush;

				ArrayLattice<Complex> c(vpad);
				cout << "c.shape() = " << c.shape() << endl;
				//	  LatticeFFT::cfft(c,ax);
				LatticeFFT::cfft2d(c,false);

				cout << "done." << endl;

				*amp_[vbnum] = amplitude(vpad(sl));
			}
			break;
		}
		case PMS::MODEL: {
            if (averaging_.scalarAve()) 
			    *ampModel_[vbnum] = real(vb->visCubeModel());
            else
			    *ampModel_[vbnum] = amplitude(vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED: {
            if (averaging_.scalarAve()) 
			    *ampCorr_[vbnum] = real(vb->visCubeCorrected());
            else
			    *ampCorr_[vbnum] = amplitude(vb->visCubeCorrected());
			break;
		}
		case PMS::CORRMODEL: {
            if (averaging_.scalarAve()) 
			  *ampCorrModel_[vbnum] = 
                real(vb->visCubeCorrected() - vb->visCubeModel());
            else
			  *ampCorrModel_[vbnum] = 
                amplitude(vb->visCubeCorrected() - vb->visCubeModel());
			break;
		}
		case PMS::DATAMODEL: {
            if (averaging_.scalarAve()) 
			  *ampDataModel_[vbnum] = 
                real(vb->visCube() - vb->visCubeModel());
            else
			  *ampDataModel_[vbnum] = 
                amplitude(vb->visCube() - vb->visCubeModel());
			break;
		}
		case PMS::DATA_DIVIDE_MODEL: {
            if (averaging_.scalarAve()) 
			  *ampDataDivModel_[vbnum] = 
                real( vb->visCube() / vb->visCubeModel());
            else
			  *ampDataDivModel_[vbnum] = 
                amplitude( vb->visCube() / vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED_DIVIDE_MODEL: {
            if (averaging_.scalarAve()) 
			  *ampCorrDivModel_[vbnum] = 
                real( vb->visCubeCorrected() / vb->visCubeModel());
            else
			  *ampCorrDivModel_[vbnum] = 
                amplitude( vb->visCubeCorrected() / vb->visCubeModel());
			break;
		}
		case PMS::FLOAT_DATA: {
			*ampFloat_[vbnum] = vb->visCubeFloat();
			break;
		}
		}
		break;
	}
	case PMS::PHASE: {
		switch(data) {
		case PMS::DATA: {
            if (averaging_.scalarAve()) 
                *pha_[vbnum]=imag(vb->visCube()) * 180.0 / C::pi;
            else
			    *pha_[vbnum] = phase(vb->visCube()) * 180.0 / C::pi;
			break;
		}
		case PMS::MODEL: {
            if (averaging_.scalarAve()) 
			    *phaModel_[vbnum] = imag(vb->visCubeModel()) * 180.0 / C::pi;
            else
			    *phaModel_[vbnum] = phase(vb->visCubeModel()) * 180.0 / C::pi;
			break;
		}
		case PMS::CORRECTED: {
            if (averaging_.scalarAve()) 
                *phaCorr_[vbnum]=imag(vb->visCubeCorrected()) * 180.0 / C::pi;
            else
			    *phaCorr_[vbnum] = phase(vb->visCubeCorrected()) * 180.0 / C::pi;
			break;
		}
		case PMS::CORRMODEL: {
            if (averaging_.scalarAve()) 
			*phaCorrModel_[vbnum] = 
                imag(vb->visCubeCorrected() - vb->visCubeModel()) * 180.0 / C::pi;
            else
			*phaCorrModel_[vbnum] = 
                phase(vb->visCubeCorrected() - vb->visCubeModel()) * 180.0 / C::pi;
			break;
		}
		case PMS::DATAMODEL: {
            if (averaging_.scalarAve()) 
			    *phaDataModel_[vbnum] = 
                    imag(vb->visCube() - vb->visCubeModel()) * 180.0 / C::pi;
            else
			    *phaDataModel_[vbnum] = 
                    phase(vb->visCube() - vb->visCubeModel()) * 180.0 / C::pi;
			break;
		}
		case PMS::DATA_DIVIDE_MODEL: {
            if (averaging_.scalarAve()) 
			    *phaDataDivModel_[vbnum] = 
                    imag(vb->visCube() / vb->visCubeModel()) * 180.0 / C::pi;
            else
			    *phaDataDivModel_[vbnum] = 
                    phase(vb->visCube() / vb->visCubeModel()) * 180.0 / C::pi;
			break;
		}
		case PMS::CORRECTED_DIVIDE_MODEL: {
            if (averaging_.scalarAve()) 
			    *phaCorrDivModel_[vbnum] = 
                  imag(vb->visCubeCorrected() / vb->visCubeModel()) * 180.0 / C::pi;
            else
			    *phaCorrDivModel_[vbnum] = 
                  phase(vb->visCubeCorrected() / vb->visCubeModel()) * 180.0 / C::pi;
			break;
		}
		case PMS::FLOAT_DATA:  // should have caught this already
			break;
		}
		break;
	}

	case PMS::REAL: {
		switch(data) {
		case PMS::DATA: {
			*real_[vbnum] = real(vb->visCube());
			break;
		}
		case PMS::MODEL: {
			*realModel_[vbnum] = real(vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED: {
			*realCorr_[vbnum] = real(vb->visCubeCorrected());
			break;
		}
		case PMS::CORRMODEL: {
			*realCorrModel_[vbnum] = 
                real(vb->visCubeCorrected()) - real(vb->visCubeModel());
			break;
		}
		case PMS::DATAMODEL: {
			*realDataModel_[vbnum] = 
                real(vb->visCube()) - real(vb->visCubeModel());
			break;
		}
		case PMS::DATA_DIVIDE_MODEL: {
			*realDataDivModel_[vbnum] = 
                real(vb->visCube()) / real(vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED_DIVIDE_MODEL: {
			*realCorrDivModel_[vbnum] = 
                real(vb->visCubeCorrected()) / real(vb->visCubeModel());
			break;
		}
		case PMS::FLOAT_DATA: {
            *real_[vbnum] = vb->visCubeFloat();  // float data is real
			break;
		}
        }
		break;
	}
	case PMS::IMAG: {
		switch(data) {
		case PMS::DATA: {
			*imag_[vbnum] = imag(vb->visCube());
			break;
		}
		case PMS::MODEL: {
			*imagModel_[vbnum] = imag(vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED: {
			*imagCorr_[vbnum] = imag(vb->visCubeCorrected());
			break;
		}
		case PMS::CORRMODEL: {
			*imagCorrModel_[vbnum] = 
                imag(vb->visCubeCorrected()) - imag(vb->visCubeModel());
			break;
		}
		case PMS::DATAMODEL: {
			*imagDataModel_[vbnum] = 
                imag(vb->visCube()) - imag(vb->visCubeModel());
			break;
		}
		case PMS::DATA_DIVIDE_MODEL: {
			*imagDataDivModel_[vbnum] = 
                imag(vb->visCube()) / imag(vb->visCubeModel());
			break;
		}
		case PMS::CORRECTED_DIVIDE_MODEL: {
			*imagCorrDivModel_[vbnum] = 
                imag(vb->visCubeCorrected()) / imag(vb->visCubeModel());
			break;
		}
		case PMS::FLOAT_DATA:  // should have caught this already
			break;
		}
		break;
	}

	case PMS::FLAG: {
		*flag_[vbnum] = vb->flagCube();
		break;
    }
	case PMS::FLAG_ROW: {
		*flagrow_[vbnum] = vb->flagRow();
		break;
	}

	case PMS::WT: {
		*wt_[vbnum] = vb->weight();
		break;
	}
	case PMS::WTSP: {
	    if (vb->weightSpectrum().nelements()>0)
			*wtsp_[vbnum] = vb->weightSpectrum();
		else
	  		throw(AipsError("This MS does not have a valid WEIGHT_SPECTRUM column."));
        break;
	}

	case PMS::WTxAMP: {
		uInt nchannels = vb->nChannels();
		switch(data) {
		case PMS::DATA: {
			*wtxamp_[vbnum] = amplitude(vb->visCube());
		    Cube<Float> wtA(*wtxamp_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::MODEL: {
			*wtxampModel_[vbnum] = amplitude(vb->visCubeModel());
		    Cube<Float> wtA(*wtxampModel_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::CORRECTED: {
			*wtxampCorr_[vbnum] = amplitude(vb->visCubeCorrected());
		    Cube<Float> wtA(*wtxampCorr_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::CORRMODEL: {
			*wtxampCorrModel_[vbnum] = 
                amplitude(vb->visCubeCorrected() - vb->visCube());
		    Cube<Float> wtA(*wtxampCorrModel_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::DATAMODEL: {
			*wtxampDataModel_[vbnum] = 
                amplitude(vb->visCube() - vb->visCubeModel());
		    Cube<Float> wtA(*wtxampDataModel_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::DATA_DIVIDE_MODEL: {
			*wtxampDataDivModel_[vbnum] = 
                amplitude(vb->visCube() / vb->visCubeModel());
		    Cube<Float> wtA(*wtxampDataDivModel_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::CORRECTED_DIVIDE_MODEL: {
			*wtxampCorrDivModel_[vbnum] = 
                amplitude(vb->visCubeCorrected() / vb->visCubeModel());
		    Cube<Float> wtA(*wtxampCorrDivModel_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		case PMS::FLOAT_DATA: {
			*wtxampFloat_[vbnum] = vb->visCubeFloat();
		    Cube<Float> wtA(*wtxampFloat_[vbnum]);
		    for(uInt c = 0; c < nchannels; ++c) {
			    wtA.xzPlane(c) = wtA.xzPlane(c) * vb->weight();
		    }
			break;
		}
		}
        break;
	}

	case PMS::SIGMA: {
		*sigma_[vbnum] = vb->sigma();
		break;
	}
	case PMS::SIGMASP: {
		*sigmasp_[vbnum] = vb->sigmaSpectrum();
		break;
	}

	case PMS::AZ0:
	case PMS::EL0: {
		MDirection azelMDir = vb->azel0(vb->time()(0));
		Vector<Double> azelVec = azelMDir.getAngle("deg").getValue();
		az0_(vbnum) = azelVec(0);
		el0_(vbnum) = azelVec(1);
		break;
	}

	case PMS::HA0: {
		ha0_(vbnum) = vb->hourang(vb->time()(0))*12/C::pi;  // in hours
		break;
	}
	case PMS::PA0: {
		pa0_(vbnum) = vb->parang0(vb->time()(0))*180.0/C::pi; // in degrees
		if (pa0_(vbnum)<0.0) pa0_(vbnum)+=360.0;
		break;
	}
	case PMS::ANTENNA: {
		antenna_[vbnum]->resize(vb->nAntennas());
		indgen(*antenna_[vbnum]);
		break;
	}
	case PMS::AZIMUTH:
	case PMS::ELEVATION: {
		Vector<MDirection> azelVec = vb->azel(vb->time()(0));
		Matrix<Double> azelMat;
                azelMat.resize(2, azelVec.nelements());
		for (uInt iant = 0; iant < azelVec.nelements(); ++iant) {
			azelMat.column(iant) = (azelVec(iant).getAngle("deg").getValue());
		}
		*az_[vbnum] = azelMat.row(0);
		*el_[vbnum] = azelMat.row(1);
		break;
	}
	case PMS::RADIAL_VELOCITY: {
		Int fieldId = vb->fieldId()(0);
		const ROMSFieldColumns& fieldColumns = vi_p->subtableColumns().field();
		MRadialVelocity radVelocity = fieldColumns.radVelMeas(fieldId, vb->time()(0));
		radialVelocity_(vbnum) = radVelocity.get("AU/d").getValue( "km/s");
		break;
	}
	case PMS::RHO:{
		Int fieldId = vb->fieldId()(0);
		const ROMSFieldColumns& fieldColumns = vi_p->subtableColumns().field();
		Quantity rhoQuantity = fieldColumns.rho(fieldId, vb->time()(0));
		rho_(vbnum ) = rhoQuantity.getValue( "km");
		break;
	}
	case PMS::PARANG: {
		*parang_[vbnum] = vb->feedPa(vb->time()(0))*(180.0/C::pi);  // in degrees
		break;
	}

	case PMS::ROW: {
		*row_[vbnum] = vb->rowIds();
		break;
	}

	case PMS::OBSERVATION: {
		*obsid_[vbnum] = vb->observationId();
		break;
	}

	case PMS::INTENT:{
		Vector<Int> states = vb->stateId();
		*intent_[vbnum] = assignIntentIds(states);
		break;
	}

	case PMS::FEED1:{
		*feed1_[vbnum] = vb->feed1();
		break;
	}
	case PMS::FEED2:{
		*feed2_[vbnum] = vb->feed2();
		break;
	}

	default: {
		throw(AipsError("Axis choice not supported for MS"));
		break;
	}
	}
}

bool MSCache::isEphemeris(){
	if ( !ephemerisInitialized ){
	        Table::TableOption tabopt(Table::Old);
		MeasurementSet ms(filename_,TableLock(TableLock::AutoLocking), tabopt);
		ROMSColumns msc(ms);

                // Check the field subtable for ephemeris fields
		const ROMSFieldColumns& fieldColumns = msc.field();
                uInt nrow = fieldColumns.nrow();

		ephemerisAvailable = false;
                String ephemPath;
                for (uInt i=0; i<nrow; ++i) {
                        ephemPath = fieldColumns.ephemPath(i);
		        if (!ephemPath.empty()) ephemerisAvailable = true;
		}
		ephemerisInitialized = true;
	}
	return ephemerisAvailable;
}


void MSCache::flagToDisk(const PlotMSFlagging& flagging,
		Vector<Int>& flchunks, Vector<Int>& flrelids,
		Bool setFlag, PlotMSIndexer* indexer, int dataIndex) {

	// Sort the flags by chunk and relative index:
	Sort sorter;
	sorter.sortKey(flchunks.data(),TpInt);
	sorter.sortKey(flrelids.data(),TpInt);

	Vector<uInt> order;  // holds the sort order (indices) for flchunks
	uInt nflag = sorter.sort(order, flchunks.nelements());
    uInt iflag(0);  // index into 'order' array

    bool pmsavg = averaging_.baseline() || averaging_.antenna() || 
        averaging_.spw() || averaging_.scalarAve();
    // check if extending flags and subset selected
    bool extendcorr = flagging.corr() && !selection_.corr().empty();
    bool extendchan = flagging.channel() && !selection_.spw().empty();

    // get loadedAxes and loadedData for setting up vis iter
    vector<PMS::Axis> loadedAxes;
	vector<PMS::DataColumn> loadedData;
    for (std::map<PMS::Axis, bool>::iterator mapit=loadedAxes_.begin(); 
            mapit!=loadedAxes_.end(); ++mapit) {
        if (mapit->second) {
            PMS::Axis loadedAxis = mapit->first;
            if (loadedAxesData_.find(loadedAxis) != loadedAxesData_.end()) {
                std::set<PMS::DataColumn> datacols=loadedAxesData_[loadedAxis];
                for (auto it=datacols.begin(); it!=datacols.end(); ++it) {
                    loadedAxes.push_back(loadedAxis);
                    loadedData.push_back(*it);
                }
            } else {
                loadedAxes.push_back(loadedAxis);
                loadedData.push_back(PMS::DATA);
            }
        }
    }

    // If extending flags, select all correlations or channels in vis iter
    PlotMSSelection flagSel = selection_;
    if (extendcorr) {
        // select all corrs
        flagSel.setCorr("");
    }
    if (extendchan) {
        // select spws only (all channels per spw)
        MeasurementSet ms(filename_);
        MSSelection mssel(ms);
        mssel.setSpwExpr(selection_.spw());
        Vector<Int> spws = mssel.getSpwList();
        String spwExpr = casacore::String::toString(spws(0));
        for (uInt spw=1; spw<spws.size(); ++spw) 
            spwExpr += "," + casacore::String::toString(spws(spw));
        flagSel.setSpw(spwExpr);
    }

    setUpVisIter(flagSel, calibration_, dataColumn_, loadedAxes,
        loadedData, true, false);
    vi_p->originChunks();
    vi_p->origin();
    vi::VisBuffer2* vb = vi_p->getVisBuffer();

    // iterate through chunks and find the ones that need flags changed
    for (Int ichk=0; ichk<nChunk_; ++ichk) {
        if (ichk != flchunks(order[iflag])) {
            // Step over current chunk if not in flchunks
            for (Int i=0;i<nVBPerAve_(ichk);++i) {
                vi_p->next();
                if (!vi_p->more() && vi_p->moreChunks()) {
                    vi_p->nextChunk();
                    vi_p->origin();
                }
            }
        } else if ( pmsavg || extendcorr || extendchan) {
            // This chunk requires flag-setting but have to handle chunks
            // (shape has changed: VBs averaged together or extending flags)
            Cube<Bool> vbflag;
            Vector<Bool> vbflagrow;
            Vector<Int> a1, a2;

            for (Int i=0; i<nVBPerAve_(ichk); ++i) {

                // Refer to VB pieces we need
                vbflag = vb->flagCube();
                vbflagrow = vb->flagRow();
                a1 = vb->antenna1();
                a2 = vb->antenna2();
                Int ncorr = vb->nCorrelations();
                Int nchan = vb->nChannels();
                Int nrow  = vb->nRows();

                // Apply all flags in this chunk to this VB
                while (iflag<nflag && flchunks(order[iflag])==ichk) {

                    Int currChunk = flchunks(order[iflag]);
                    Int irel = flrelids(order[iflag]);
                    Slice corr,chan,bsln;

                    // Set flag range on correlation axis:
                    if (netAxesMask_[dataIndex](0) && !flagging.corrAll()) {
                        // A specific single correlation
                        Int icorr = indexer->getIndex1000(currChunk, irel);
                        corr = Slice(icorr, 1, 1);
                    } else {
                        // Extend to all correlations
                        corr=Slice(0, ncorr, 1);
                    }

                    // Set Flag range on channel axis:
                    if (netAxesMask_[dataIndex](1) && !flagging.channel()) {
                        // A single specific channel
                        Int ichan = indexer->getIndex0100(currChunk, irel);
                        chan = Slice(ichan, 1, 1);
                    } else {
                        // Extend to all channels
                        chan = Slice(0, nchan, 1);
                    }

                    // Set Flags on the baseline axis:
                    Int thisAnt1 = Int(getAnt1(currChunk,
                        indexer->getIndex0010(currChunk, irel)));
                    Int thisAnt2 = Int(getAnt2(currChunk,
                        indexer->getIndex0010(currChunk, irel)));
                    if (netAxesMask_[dataIndex](2) &&
                        !flagging.antennaBaselinesBased() && (thisAnt1 > -1) ) {
                        // i.e., if baseline is an explicit data axis,
                        //       full baseline extension is OFF
                        //       and the first antenna in the selected point is > -1
                        // Do some variety of detailed per-baseline flagging
                        for (Int irow=0; irow<nrow; ++irow) {
                            if (thisAnt2 > -1) {
                                // match a baseline exactly
                                if ((a1(irow)==thisAnt1) &&
                                    (a2(irow)==thisAnt2)) {
                                    vbflag(corr, chan, Slice(irow,1,1)) = setFlag;
                                    // unset flag_row when unflagging
                                    if (!setFlag) vbflagrow(irow) = false;
                                    break;  // found the one baseline, escape from for loop
                                }
                            } else {
                                // either antenna matches the one specified antenna
                                //  (don't break because there will be more than one)
                                //  TBD: this doesn't get cross-hands quite right when
                                //    averaging 'by antenna'...
                                if ((a1(irow) == thisAnt1) ||
                                    (a2(irow)==thisAnt1)) {
                                    vbflag(corr, chan, Slice(irow,1,1)) = setFlag;
                                    // unset flag_row when unflagging
                                    if (!setFlag) vbflagrow(irow) = false;
                                }
                            }
                        }
                    } else {
                        // Set flags for all baselines, because the plot
                        //  is ordinarily implicit in baseline, we've turned on baseline
                        //  extension, or we've averaged over all baselines
                        bsln = Slice(0, nrow, 1);
                        vbflag(corr, chan, bsln) = setFlag;
                        // unset flag_row when unflagging
                        if (!setFlag) vbflagrow(bsln) = false;
                    }
                    ++iflag;
                } // done with this flchunk

                // Put the flags back into the MS
                vi_p->writeFlag(vbflag);
                // Advance to the next vb
                vi_p->next();
                if (!vi_p->more() && vi_p->moreChunks()) {
                    vi_p->nextChunk();
                    vi_p->origin();
                }
            }  // VBs in this averaging chunk

            // Escape if we are already finished
            if (iflag >= nflag) break;

        } else {
            // same shape (possibly averaging handled by MSTransformIterator)
            // just write out flag chunk
            vi_p->writeFlag(flag(ichk));

            // Advance to the next vb
            vi_p->next();
            if (!vi_p->more() && vi_p->moreChunks()) {
                vi_p->nextChunk();
                vi_p->origin();
            }
            // Advance to the next flagged chunk
            ++iflag;
            while ((iflag < nflag) && 
                   (flchunks(order[iflag]) == flchunks(order[iflag-1]))) {
                ++iflag;
            }

            // Escape if we are finished
            if (iflag >= nflag) break;
        } 
    }
	// Delete the VisIter so lock is released
	deleteVi();
}

void MSCache::mapIntentNamesToIds() {
	intentIds_.clear();
	Int newId = 0;

	for (uInt i=0; i<intentnames_.size(); i++) {
		if (intentIds_.find(intentnames_[i]) == intentIds_.end()) {
			intentIds_[intentnames_[i]] = newId;
			newId++;
		}
	}
}

Vector<Int> MSCache::assignIntentIds(Vector<Int>& stateIds) {
	Vector<Int> intents;
	intents.resize(stateIds.size());

	for (uInt i=0; i<stateIds.size(); i++) {
		if ((intentnames_.size() > 0) && (stateIds[i] >= 0))
			intents[i] = intentIds_[intentnames_[stateIds[i]]];
		else
			intents[i] = stateIds[i];
	}

	return intents;
}

Vector<Double> MSCache::calcVelocity(vi::VisBuffer2* vb) {
	// Convert freq in the vb to velocity
	Vector<Double> outVel;
	VisBufferUtil vbu = VisBufferUtil();
	vbu.toVelocity(outVel,
			*vb,
			freqFrame_,
			MVFrequency(transformations_.restFreqHz()),
			transformations_.veldef());
	outVel /= 1.0e3;  // in km/s
	return outVel;
}

}
