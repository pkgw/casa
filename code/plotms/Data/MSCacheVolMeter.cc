//# MSCacheVolMeter.cc: Implementation of MSCache Volume meter
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
#include <plotms/Data/MSCacheVolMeter.h>

#include <casa/OS/Timer.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/Memory.h>
#include <casa/Quanta/MVTime.h>
#include <casa/System/Aipsrc.h>
#include <casa/Utilities/Sort.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>
#include <ms/MeasurementSets/MSColumns.h> 	 
#include <msvis/MSVis/VisSet.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisBufferUtil.h>
#include <plotms/Data/PlotMSVBAverager.h>
#include <plotms/PlotMS/PlotMS.h>
#include <tables/Tables/Table.h>
#include <synthesis/Utilities/AppRC.h>

using namespace casacore;
namespace casa {

MSCacheVolMeter::MSCacheVolMeter():
		  nDDID_(0),
		  nPerDDID_(),
		  nRowsPerDDID_(),
		  nChanPerDDID_(),
		  nCorrPerDDID_(),
		  nAnt_(0) {}


MSCacheVolMeter::MSCacheVolMeter(const MeasurementSet& ms, 
		const PlotMSAveraging ave,
		const Vector<Vector<Slice> >& chansel,
		const Vector<Vector<Slice> >& corrsel):
		nDDID_(0),
		nPerDDID_(),
		nRowsPerDDID_(),
		nChanPerDDID_(),
		nCorrPerDDID_(),
		nAnt_(0) {

	ROMSColumns msCol(ms);

	// Initialize chunks and rows counters
	nDDID_ = msCol.dataDescription().nrow();
	nPerDDID_.resize(nDDID_);
	nPerDDID_.set(0);
	nRowsPerDDID_.resize(nDDID_);
	nRowsPerDDID_.set(0);
	nChanPerDDID_.resize(nDDID_);
	nCorrPerDDID_.resize(nDDID_);

	// Fill Corr/Chan-per-DDID Vectors
	Vector<Int> nChanPerSpw;
	msCol.spectralWindow().numChan().getColumn(nChanPerSpw);
	Vector<Int> nCorrPerPol;
	msCol.polarization().numCorr().getColumn(nCorrPerPol);
	Vector<Int> polPerDDID;
	msCol.dataDescription().polarizationId().getColumn(polPerDDID);
	Vector<Int> spwPerDDID;
	msCol.dataDescription().spectralWindowId().getColumn(spwPerDDID);

	// nAnt:
	nAnt_ = msCol.antenna().nrow();

	Bool chave = (ave.channel() && ave.channelValue()>1.0);

	for (Int iddid=0; iddid<nDDID_; ++iddid) {
		// ncorr is simple (for now, maybe Stokes later?):
		Int ipol = polPerDDID(iddid);
		// Unselected corr counting
		nCorrPerDDID_(iddid) = nCorrPerPol(ipol);

		if (corrsel.nelements()>0 && corrsel(ipol).nelements()>0) {
			Vector<Slice> s(corrsel(ipol));
			Int ncorr0 = 0;
			for (uInt j=0; j<s.nelements(); ++j) ncorr0 += s(j).length();
			nCorrPerDDID_(iddid) = ncorr0;
		}

		Int ispw = spwPerDDID(iddid);
		// Unselected channel counting
		Int nchan0 = nChanPerSpw(spwPerDDID(iddid));
		Int nchanA = ( chave ? Int(ceil(Double(nchan0)/ave.channelValue())) : nchan0 );

		// Override channel counting if specific selection provided
		if (chansel.nelements()>0 && chansel(ispw).nelements()>0) {
			Vector<Slice> s(chansel(ispw));
			Int minchan(INT_MAX),maxchan(-1);
			nchan0 = 0;
			for (uInt j=0; j<s.nelements(); ++j) {
				nchan0 += s(j).length();  // unaveraged chan count
				minchan = min(minchan, Int(s(j).start()));
				maxchan = max(maxchan, Int(s(j).end()));
			}
			// Average sees full range of channels
			nchanA = Int(ceil(Double(maxchan-minchan+1)/ave.channelValue()));
		}
		// record channel counting result
		nChanPerDDID_(iddid) = (chave ? nchanA : nchan0);
	}
}


MSCacheVolMeter::~MSCacheVolMeter() {}


void MSCacheVolMeter::reset() {
	nDDID_=0;
	nPerDDID_.resize();
	nRowsPerDDID_.resize();
	nChanPerDDID_.resize();
	nCorrPerDDID_.resize();
	nAnt_=0;
}


void MSCacheVolMeter::add(Int DDID, Int nRows) {
	++nPerDDID_(DDID);
	nRowsPerDDID_(DDID) += nRows;
}

void MSCacheVolMeter::add(const vi::VisBuffer2* vb) {
	this->add(vb->dataDescriptionIds()(0), vb->nRows());
}

String MSCacheVolMeter::evalVolume(map<PMS::Axis,Bool> axes, Vector<Bool> axesmask) {
  /*
  cout << "nPerDDID_     = " << nPerDDID_ << endl;
  cout << "nRowsPerDDID_ = " << nRowsPerDDID_ << endl;
  cout << "nChanPerDDID_ = " << nChanPerDDID_ << endl;
  cout << "nCorrPerDDID_ = " << nCorrPerDDID_ << endl;

  cout << "Product = " << nRowsPerDDID_*nChanPerDDID_*nCorrPerDDID_ << endl;
  cout << "Sum     = " << sum(nRowsPerDDID_*nChanPerDDID_*nCorrPerDDID_) << endl;

  cout << "sizeof(Int)    = " << sizeof(Int) << endl;
  cout << "sizeof(Long)   = " << sizeof(Long) << endl;
  cout << "sizeof(uInt64)   = " << sizeof(uInt64) << endl;
  cout << "sizeof(Float)  = " << sizeof(Float) << endl;
  cout << "sizeof(Double) = " << sizeof(Double) << endl;
  */
	uInt64 totalVol(0);
	for (map<PMS::Axis,Bool>::iterator pAi=axes.begin();
			pAi!=axes.end(); ++pAi) {
		if (pAi->second) {
			uInt64 axisVol(0);
			switch(pAi->first) {
			case PMS::SCAN:
			case PMS::FIELD:
			case PMS::SPW:
			case PMS::OBSERVATION:
			case PMS::INTENT:
			case PMS::FEED1:
			case PMS::FEED2:
				axisVol = sizeof(Int) * sum(nPerDDID_);
				break;
			case PMS::TIME:
			case PMS::TIME_INTERVAL:
				axisVol = sizeof(Double) * sum(nPerDDID_);
				break;
			case PMS::CHANNEL:
				axisVol = sizeof(Int) * sum(nPerDDID_ * nChanPerDDID_);
				break;
			case PMS::FREQUENCY:
			case PMS::VELOCITY:
				axisVol = sizeof(Double) * sum(nPerDDID_ * nChanPerDDID_);
				break;
			case PMS::CORR:
				axisVol = sizeof(Int) * sum(nPerDDID_ * nCorrPerDDID_);
				break;
			case PMS::ANTENNA1:
			case PMS::ANTENNA2:
			case PMS::BASELINE:
				axisVol = sizeof(Int) * sum(nRowsPerDDID_);
				break;
			case PMS::UVDIST:
			case PMS::U:
			case PMS::V:
			case PMS::W:
				axisVol = sizeof(Double) * sum(nRowsPerDDID_);
				break;
			case PMS::UVDIST_L:
			case PMS::UWAVE:
			case PMS::VWAVE:
			case PMS::WWAVE:
				axisVol = sizeof(Double) * sum(nRowsPerDDID_ * nChanPerDDID_);
				break;
			case PMS::AMP:
			case PMS::PHASE:
			case PMS::REAL:
			case PMS::IMAG:
			case PMS::WTSP:
			case PMS::SIGMASP:
				axisVol = sizeof(Float) * sum(nRowsPerDDID_ * nChanPerDDID_ * nCorrPerDDID_);
				break;
			case PMS::FLAG:
				axisVol = sizeof(Bool) * sum(nRowsPerDDID_ * nChanPerDDID_ * nCorrPerDDID_);
				break;
			case PMS::FLAG_ROW:
				axisVol = sizeof(Bool) * sum(nRowsPerDDID_);
				break;
			case PMS::WT:
			case PMS::SIGMA:
				axisVol = sizeof(Float) * sum(nRowsPerDDID_ * nCorrPerDDID_);
				break;
			case PMS::AZ0:
			case PMS::EL0:
			case PMS::RADIAL_VELOCITY:
			case PMS::RHO:
			case PMS::HA0:
			case PMS::PA0:
				axisVol = sizeof(Double) * sum(nPerDDID_);
				break;
			case PMS::ANTENNA:
				axisVol = sizeof(Int) * nAnt_ * sum(nPerDDID_);
				break;
			case PMS::AZIMUTH:
			case PMS::ELEVATION:
				axisVol = sizeof(Double) * nAnt_ * sum(nPerDDID_);
				break;
			case PMS::PARANG:
				axisVol = sizeof(Float) * nAnt_ * sum(nPerDDID_);
				break;
			case PMS::ROW:
				axisVol = sizeof(uInt) * sum(nRowsPerDDID_);
				break;
			case PMS::GAMP:
			case PMS::GPHASE:
			case PMS::GREAL:
			case PMS::GIMAG:
			case PMS::DELAY:
			case PMS::SWP:
			case PMS::TSYS:
			case PMS::OPAC:
			case PMS::SNR:
			case PMS::TEC:
			case PMS::WTxAMP:
			case PMS::NONE:
				break;
			} // switch
			totalVol += axisVol;
			//cout << " " << PMS::axis(pAi->first) << " volume = " << axisVol << " bytes." << endl;
		}
	} // for

	// Add in the plotting mask
	//  (TBD: only if does not reference the flags)
	if (true) {  // ntrue(axesmask)<2) {
		Vector<uInt64> nplmaskPerDDID(nDDID_, 0);
		nplmaskPerDDID(nPerDDID_>uInt64(0)) = 1;
		if (axesmask(0)) nplmaskPerDDID *= nCorrPerDDID_;
		if (axesmask(1)) nplmaskPerDDID *= nChanPerDDID_;
		if (axesmask(2)) nplmaskPerDDID *= nRowsPerDDID_;
		if (axesmask(3)) nplmaskPerDDID *= uInt64(nAnt_);
		uInt64 plmaskVol = sizeof(Bool) * sum(nplmaskPerDDID);
		//    cout << " Collapsed flag (plot mask) volume = " << plmaskVol << " bytes." << endl;
		totalVol += plmaskVol;
	}

	// Finally, count the total points for the plot:
	Vector<uInt64> nPointsPerDDID(nDDID_,0);
	nPointsPerDDID(nPerDDID_>uInt64(0)) = 1;
	if (axesmask(0)) nPointsPerDDID *= nCorrPerDDID_;
	if (axesmask(1)) nPointsPerDDID *= nChanPerDDID_;
	if (axesmask(2)) nPointsPerDDID *= nRowsPerDDID_;
	if (axesmask(3)) nPointsPerDDID *= uInt64(nAnt_);

	uInt64 totalPoints = sum(nPointsPerDDID);

	Double totalVolGB = Double(totalVol)/1.0e9;  // in GB
	Double bytesPerPt = Double(totalVol)/Double(totalPoints);  // bytes/pt

	// Memory limit from CASA_MAX_MEMORY env variable, 
    // .casarc, or host MemTotal (in that order) in Bytes converted to GB
	Double memAvailGB = AppRC::getMemoryAvailable() / 1.0e9;
	Double fracMem = 100.0 * totalVolGB / memAvailGB;  // fraction require in %
 
	// Report number of points to be plotted.
	stringstream ss;
	ss << "Data selection will yield a total of " << totalPoints
			<< " plottable points (flagged and unflagged).";
	// Report required memory
	ss << endl<< "The plotms cache will require an estimated "
			<< totalVolGB << " GB of memory (" << bytesPerPt << " bytes/point)." << endl
			<< "This is " << fracMem << "% of the memory avail. to CASA ("
			<< memAvailGB << " GB).";

	// Trap too many points
	Bool toomany(false);
	if (totalPoints > UINT_MAX) {
		ss << endl << "Too many points!  CASA plotms cannot plot more than " << UINT_MAX << " points.";
		toomany=true;
	}
	// Trap insufficient memory
	Bool toomuch(false);
	if (totalVolGB > memAvailGB) {
		ss << endl << "Insufficient memory!";
		toomuch=true;
	}
	// Throw exception if toomuch or toomany
	if (toomuch || toomany)
		throw(AipsError(ss.str()));
	return ss.str();
}

// =======================================================================
// Same method with visbuffer structure from MSTransformIteratorFactory
String MSCacheVolMeter::evalVolume(std::vector<IPosition> vbShapes,
	map<PMS::Axis,Bool> axes) {

	Int nChunks = vbShapes.size();
	Vector<Int> nCorr(nChunks);
	Vector<Int> nChan(nChunks);
	Vector<Int> nRows(nChunks);
	uInt64 nElements = 0;

	for (uInt i=0; i < vbShapes.size(); i++)
	{
		nCorr(i) = vbShapes[i][0];
		nChan(i) = vbShapes[i][1];
		nRows(i) = vbShapes[i][2];
		nElements += nCorr(i) * nChan(i) * nRows(i);
	}
	//cout << "nElements = " << nElements << endl;

	uInt64 totalVol(0);
	for (map<PMS::Axis,Bool>::iterator pAi=axes.begin();
			pAi!=axes.end(); ++pAi) {
		if (pAi->second) {
			uInt64 axisVol(0);
			switch(pAi->first) {
			case PMS::SCAN:
			case PMS::FIELD:
			case PMS::SPW:
			case PMS::OBSERVATION:
			case PMS::INTENT:
			case PMS::FEED1:
			case PMS::FEED2:
				axisVol = sizeof(Int) * nChunks;
				break;
			case PMS::TIME:
			case PMS::TIME_INTERVAL:
				axisVol = sizeof(Double) * nChunks;
				break;
			case PMS::CHANNEL:
				axisVol = sizeof(Int) * sum(nChan);
				break;
			case PMS::FREQUENCY:
			case PMS::VELOCITY:
				axisVol = sizeof(Double) * sum(nChan);
				break;
			case PMS::CORR:
				axisVol = sizeof(Int) * sum(nCorr);
				break;
			case PMS::ANTENNA1:
			case PMS::ANTENNA2:
			case PMS::BASELINE:
				axisVol = sizeof(Int) * sum(nRows);
				break;
			case PMS::UVDIST:
			case PMS::U:
			case PMS::V:
			case PMS::W:
				axisVol = sizeof(Double) * sum(nRows);
				break;
			case PMS::UVDIST_L:
			case PMS::UWAVE:
			case PMS::VWAVE:
			case PMS::WWAVE:
				axisVol = sizeof(Double) * sum(nRows * nChan);
				break;
			case PMS::AMP:
			case PMS::PHASE:
			case PMS::REAL:
			case PMS::IMAG:
			case PMS::WTxAMP:
			case PMS::WTSP:
			case PMS::SIGMASP:
				axisVol = sizeof(Float) * nElements;
				break;
			case PMS::FLAG:
				axisVol = sizeof(Bool) * nElements;
				break;
			case PMS::FLAG_ROW:
				axisVol = sizeof(Bool) * sum(nRows);
				break;
			case PMS::WT:
			case PMS::SIGMA:
				axisVol = sizeof(Float) * sum(nRows * nCorr);
				break;
			case PMS::AZ0:
			case PMS::EL0:
			case PMS::RADIAL_VELOCITY:
			case PMS::RHO:
			case PMS::HA0:
			case PMS::PA0:
				axisVol = sizeof(Double) * nChunks;
				break;
			case PMS::ANTENNA:
				axisVol = sizeof(Int)*nAnt_ * nChunks;
				break;
			case PMS::AZIMUTH:
			case PMS::ELEVATION:
				axisVol = sizeof(Double)*nAnt_ * nChunks;
				break;
			case PMS::PARANG:
				axisVol = sizeof(Float)*nAnt_ * nChunks;
				break;
			case PMS::ROW:
				axisVol = sizeof(uInt) * sum(nRows);
				break;
			case PMS::GAMP:
			case PMS::GPHASE:
			case PMS::GREAL:
			case PMS::GIMAG:
			case PMS::DELAY:
			case PMS::SWP:
			case PMS::TSYS:
			case PMS::OPAC:
			case PMS::SNR:
			case PMS::TEC:
			case PMS::NONE:
				break;
			} // switch
			totalVol += axisVol;
			//cout << PMS::axis(pAi->first) << " volume = " << axisVol << " bytes." << endl;
		}
	} // for

	// Add in the plotting mask
	//  (TBD: only if does not reference the flags)
	if (true) {  // ntrue(axesmask)<2) {
		uInt64 plmaskVol = sizeof(Bool) * nElements;
		//cout << " Collapsed flag (plot mask) volume = " << plmaskVol << " bytes." << endl;
		totalVol += plmaskVol;
	}

	uInt64 totalPoints = nElements;

	Double totalVolGB = Double(totalVol) / 1.0e9;  // in GB
	Double bytesPerPt = Double(totalVol) / Double(totalPoints);  // bytes/pt

	// Memory limit from CASA_MAX_MEMORY env variable, 
    // .casarc, or host MemTotal (in that order) in Bytes converted to GB
	Double memAvailGB = AppRC::getMemoryAvailable() / 1.0e9;
	Double fracMem = 100.0 * totalVolGB / memAvailGB;  // fraction require in %

	// Report number of points to be plotted.
	stringstream ss;
	ss << "Data selection will yield a total of " << totalPoints
			<< " plottable points (flagged and unflagged).";
	// Report required memory
	ss << endl<< "The plotms cache will require an estimated "
			<< totalVolGB << " GB of memory (" << bytesPerPt << " bytes/point)." << endl
			<< "This is " << fracMem << "% of the memory avail. to CASA ("
			<< memAvailGB << " GB).";

	// Trap too many points
	Bool toomany(false);
	if (totalPoints > UINT_MAX) {
		ss << endl << "Too many points!  CASA plotms cannot plot more than " << UINT_MAX << " points.";
		toomany=true;
	}

	// Trap insufficient memory
	Bool toomuch(false);
	if (totalVolGB > memAvailGB) {
		ss << endl << "Insufficient memory!";
		toomuch=true;
	}

	// Throw exception if toomuch or toomany
	if (toomuch || toomany)
		throw(AipsError(ss.str()));

	return ss.str();
}

using namespace casacore;
}  // namespace casa
