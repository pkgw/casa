//# CalCache.cc: Specialized PlotMSCache for filling CalTables
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
#include <plotms/Data/CalCache.h>
#include <plotms/Data/PlotMSIndexer.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Threads/ThreadCommunication.h>

#include <casa/OS/Timer.h>
#include <casa/OS/HostInfo.h>
#include <casa/OS/Memory.h>
#include <casa/Quanta/MVTime.h>
#include <casa/System/Aipsrc.h>
#include <casa/Utilities/Sort.h>
#include <casa/Arrays/ArrayMath.h>
#include <tables/Tables/Table.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/MeasurementComponents/VisCalGlobals.h>
#include <synthesis/MeasurementComponents/BPoly.h>
#include <synthesis/MeasurementComponents/GSpline.h>

using namespace casacore;

namespace casa {

CalCache::CalCache(PlotMSApp* parent):
  PlotMSCacheBase(parent),
  divZero_(False),
  ci_p(NULL),
  wci_p(NULL),
  basis_("unknown"),
  parsAreComplex_(False),
  msname_("")
{
}

CalCache::~CalCache() {}


String CalCache::polname(Int ipol) {
  if (polnRatio_) return "/";
  if (basis_=="Linear")
    return ( (ipol%2==0) ? String("X") : String("Y") );
  else if (basis_=="Circular")
    return ( (ipol%2==0) ? String("R") : String("L") );
  else // "unknown"
    return ( (ipol%2==0) ? String("0") : String("1") );
}


//*********************************
// protected method implementations

void CalCache::loadIt(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& /*loadData*/,
    ThreadCommunication* thread) {

  // update cal type
  setFilename(filename_);

  // Trap unsupported modes: cal types, averaging, transforms, poln ratio
  if (calType_[0]=='A' || calType_[0]=='M' || calType_[0]=='X')
    throw AipsError("Cal table type " + calType_ + " is unsupported in plotms. Please continue to use plotcal.");
  if (calType_=="KAntPos Jones")
    throw AipsError("Cannot plot " + calType_ + " tables.");
  logLoad("Plotting a " + calType_ + " calibration table.");
  // Warn that averaging and transformations will be ignored
  if (averaging().anyAveraging())
    logWarn("CalCache::loadIt",
      "Averaging ignored: not supported for calibration tables");
  if (transformations().anyTransform())
    logWarn("CalCache::loadIt",
      "Transformations ignored: not supported for calibration tables");
  // poln ratio
  if (selection_.corr()=="/") {
    if (calType_=="BPOLY" || calType_[0] == 'T' || calType_[0] == 'F')
      throw(AipsError("Polarization ratio plots not supported for " + calType_ + " tables."));
    else
      polnRatio_ = true;
  }

  antnames_.resize();
  stanames_.resize();
  antstanames_.resize();
  fldnames_.resize();
  positions_.resize();
  polnRatio_ = false;

  vector<PMS::DataColumn> loadData(loadAxes.size());
  for (uInt i=0; i<loadData.size(); ++i) 
    loadData[i] = PMS::DEFAULT_DATACOLUMN;

  if (calType_=="BPOLY") {
    loadBPoly(loadAxes, loadData, thread);
  } else if (calType_=="GSPLINE") {
    checkAxes(loadAxes);  // check for invalid axis before proceeding
    loadGSpline(loadAxes, loadData, thread);
  } else {
    loadNewCalTable(loadAxes, loadData, thread);
  }
}

// ======================== NewCalTable ==========================

void CalCache::loadNewCalTable(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // Get various names, properties from cal table
  NewCalTable ct(NewCalTable::createCT(filename_, Table::Old, Table::Plain));
  basis_ = ct.polBasis();
  parsAreComplex_ = ct.isComplex();
  ROCTColumns ctCol(ct);
  antnames_ = ctCol.antenna().name().getColumn();
  stanames_ = ctCol.antenna().station().getColumn();
  antstanames_ = antnames_ + String("@") + stanames_;
  fldnames_ = ctCol.field().name().getColumn();
  positions_ = ctCol.antenna().position().getColumn();    
  nAnt_ = ctCol.antenna().nrow();

  setUpCalIter(filename_,selection_, True,True,True);
  countChunks(*ci_p, loadAxes, loadData, thread);
  loadCalChunks(*ci_p, loadAxes, thread);
  if (ci_p)
    delete ci_p;
  ci_p = NULL;
}

void CalCache::setUpCalIter(const String& ctname,
    PlotMSSelection& selection, Bool readonly,
    Bool /*chanselect*/, Bool /*corrselect*/) {
  // for NewCalTable
  Int nsortcol(4);
  Block<String> columns(nsortcol);
  columns[0]="SCAN_NUMBER";
  columns[1]="FIELD_ID";   
  columns[2]="SPECTRAL_WINDOW_ID";
  columns[3]="TIME";

  // Now open the MS, select on it, make the VisIter
  TableLock lock(TableLock::AutoNoReadLocking);
  Table::TableOption tabopt(Table::Update);
  if (readonly) tabopt = Table::Old;
  NewCalTable ct(ctname, lock, tabopt, Table::Plain);

  // Apply selection
  NewCalTable selct;
  Vector<Vector<Slice> > chansel;
  Vector<Vector<Slice> > corrsel;
  selection.apply(ct, selct, chansel, corrsel);

  // setup the volume meter
  //  vm_.reset();
  //  vm_= PMSCacheVolMeter(ms,averaging_,chansel,corrsel);

  if (readonly) {
    // Readonly version, for caching
    ci_p = new ROCTIter(selct,columns);
    wci_p =NULL;
  }
  else {
    // Writable, e.g. for flagging
    wci_p = new CTIter(selct,columns);
    ci_p = wci_p;  // const access
  }

  // Apply chan/corr selection
  //if (chanselect) rci_p->selectChannel(chansel);
  //if (corrselect) rci_p->selectCorrelation(corrsel);

}
      
void CalCache::countChunks(ROCTIter& ci,
    vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData,
    ThreadCommunication* thread) {
  // for NewCalTable

  if (thread!=NULL) {
    thread->setStatus("Establishing cache size.  Please wait...");
    thread->setAllowedOperations(false,false,false);
  }

  // Count number of chunks.
  int chunk(0);
  ci.reset();
  while (!ci.pastEnd()) {
    ++chunk;
    ci.next0();
  }

  setCache(chunk, loadAxes, loadData);
  //  cout << " Found " << nChunk_ << " chunks." << endl;
}

  /*
void MSCache::trapExcessVolume(map<PMS::Axis,Bool> pendingLoadAxes) {
  try {
    String s;
    s=vm_.evalVolume(pendingLoadAxes,netAxesMask(currentX_,currentY_));
    logLoad(s);
  } catch(AipsError& log) {
    // catch detected volume excess, clear the existing cache, and rethrow
    logLoad(log.getMesg());
    clear();
    stringstream ss;
    ss << "Please try selecting less data," << endl
       << "averaging, checking 'Reload' (to clear unneeded cache items)," << endl
       << "or letting other memory-intensive processes finish.";
    throw(AipsError(ss.str()));
  }
}
  */ 

void CalCache::loadCalChunks(ROCTIter& ci,
     const vector<PMS::Axis> loadAxes,
     ThreadCommunication* thread) {
  // for NewCalTable
  // permit cancel in progress meter:
  if(thread != NULL)
    thread->setAllowedOperations(false,false,true);
    
  logLoad("Loading chunks......");

  // Initialize the freq/vel calculator (in case we use it)
  //  vbu_=VisBufferUtil(vb);

  Int chunk(0);
  chshapes_.resize(4,nChunk_);
  goodChunk_.resize(nChunk_);
  goodChunk_.set(False);
  double progress;

  // Reset iterator
  ci.reset();

  while (!ci.pastEnd()) {

      // If a thread is given, check if the user canceled.
      if(thread != NULL && thread->wasCanceled()) {
        dataLoaded_ = false;
        return;
      }
      
      // If a thread is given, update it.
      if(thread != NULL && (nChunk_ <= (int)THREAD_SEGMENT ||
         chunk % THREAD_SEGMENT == 0)) {
          thread->setStatus("Loading chunk " + String::toString(chunk) +
              " / " + String::toString(nChunk_) + ".");
      }
      
      // Discern npar/nchan shape
      IPosition pshape(ci.flag().shape());
      size_t nPol;
      String pol = selection_.corr();
      if (pol=="" || pol=="RL" || pol=="XY") { // no selection
        nPol = pshape[0];
        // half the data for EVLASWP table is swp, half is tsys
        if (calType_.contains("EVLASWP")) nPol = pshape[0]/2;
        pol = "";
      } else { // poln selection using calParSlice
        String paramAxis = toVisCalAxis(PMS::AMP);
        if (polnRatio_)  // pick one!
            nPol = getParSlice(paramAxis, "R").length();
        else 
            nPol = getParSlice(paramAxis, pol).length();
      }

      // Cache the data shapes
      chshapes_(0,chunk) = nPol;
      chshapes_(1,chunk) = pshape[1];
      chshapes_(2,chunk) = ci.nrow();
      chshapes_(3,chunk) = nAnt_;
      goodChunk_(chunk) = True;

      for(unsigned int i = 0; i < loadAxes.size(); i++) {
        loadCalAxis(ci, chunk, loadAxes[i], pol);
      }
        chunk++;
        ci.next();
      
        // If a thread is given, update it.
        if(thread != NULL && (nChunk_ <= (int)THREAD_SEGMENT ||
            chunk % THREAD_SEGMENT == 0)) {
            progress = ((double)chunk+1) / nChunk_;
            thread->setProgress((unsigned int)((progress * 100) + 0.5));
        }
  }
  if (divZero_)
    logWarn("CalCache::loadIt", "Caught divide-by-zero exception in ratio plots; result(s) set to 1.0 and flagged");
}

void CalCache::loadCalAxis(ROCTIter& cti, Int chunk, PMS::Axis axis, String pol) {
    // for NewCalTable  
    Slice parSlice1 = Slice();
    Slice parSlice2 = Slice();
    if (PMS::axisNeedsCalSlice(axis)) {
        String calAxis = toVisCalAxis(axis);
        if (polnRatio_) {
            parSlice1 = getParSlice(calAxis, "R");
            parSlice2 = getParSlice(calAxis, "L");
        } else {
            parSlice1 = getParSlice(calAxis, pol);
        }
    }

  switch(axis) {
    
  case PMS::SCAN: // assumes scan unique
    scan_(chunk) = cti.thisScan();
    break;
    
  case PMS::FIELD:
    field_(chunk) = cti.thisField();
    break;
    
  case PMS::TIME: { // assumes time unique 
    time_(chunk) = cti.thisTime();
    break;
  }
  /*        
    case PMS::TIME_INTERVAL: // assumes timeInterval unique in VB
        timeIntr_(chunk) = cti.interval()(0); 
    break;
  */
  case PMS::SPW:
    spw_(chunk) = cti.thisSpw();
    break;
  case PMS::CHANNEL: {
    cti.chan(*chan_[chunk]);
    break;
  }
    case PMS::FREQUENCY: {
      // TBD: Convert freq to desired frame
      cti.freq(*freq_[chunk]);
      (*freq_[chunk])/=1.0e9; // in GHz
      break;
    }

  /*
    case PMS::VELOCITY: {
      // Convert freq in the vb to velocity
      vbu_.toVelocity(*vel_[chunk], vb, transformations_.frame(),
        MVFrequency(transformations_.restFreqHz()),
        transformations_.veldef());
      (*vel_[chunk]) /= 1.0e3;  // in km/s
      break;
    }
  */

  case PMS::CORR: {
    corr_[chunk]->resize(chshapes_(0,chunk));
    if (pol=="" || pol=="RL" || pol=="XY") {
        indgen(*corr_[chunk]);
    } else if (pol=="/") {
      corr_[chunk]->resize(1);
      corr_[chunk]->set(-1);
    } else {
      Int poln = ((pol=="R" || pol=="X") ? 0 : 1);
      corr_[chunk]->resize(1);
      corr_[chunk]->set(poln);
    }
    break;
  }
        
  case PMS::ANTENNA1:
    *antenna1_[chunk] = cti.antenna1(); 
    break;
  case PMS::ANTENNA2:
    *antenna2_[chunk] = cti.antenna2(); 
    break;
    
  case PMS::BASELINE: {
    Vector<Int> a1(cti.antenna1());
    Vector<Int> a2(cti.antenna2());
    baseline_[chunk]->resize(cti.nrow());
    Vector<Int> bl(*baseline_[chunk]);
    for (Int irow=0;irow<cti.nrow();++irow) {
      if (a1(irow)<0) a1(irow)=chshapes_(3,0);
      if (a2(irow)<0) a2(irow)=chshapes_(3,0);
      bl(irow)=(chshapes_(3,0)+1)*a1(irow) - (a1(irow) * (a1(irow) + 1)) / 2 + a2(irow);
    }
    break;
  }

  case PMS::GAMP:
  case PMS::AMP: {
    if (parsAreComplex()) {
        Cube<Complex> cArray = cti.cparam();
        if (polnRatio_) {
            Array<Float> ampRatio = amplitude(cArray(parSlice1, Slice(), Slice()) / cArray(parSlice2, Slice(), Slice()));
            checkRatioArray(ampRatio, chunk);
            *amp_[chunk] = ampRatio;
        } else {
            *amp_[chunk] = amplitude(cArray(parSlice1, Slice(), Slice()));
        }
    } else {
        Cube<Float> fArray = cti.fparam();
        if (polnRatio_) {
            Array<Float> ampRatio = fArray(parSlice1, Slice(), Slice()) / fArray(parSlice2, Slice(), Slice());
            checkRatioArray(ampRatio, chunk);
            *amp_[chunk] = ampRatio;
        } else {        
            *amp_[chunk] = fArray(parSlice1, Slice(), Slice());
        }
        if (calType_[0] == 'F') // TEC table
            (*amp_[chunk]) /= Float(1e+16);
    }
    break;
  }
    
  case PMS::GPHASE:
  case PMS::PHASE: {
    if (parsAreComplex()) {
        Cube<Complex> cArray = cti.cparam();
        if (polnRatio_) {
            Array<Float> phaseRatio = phase(cArray(parSlice1, Slice(), Slice()) / cArray(parSlice2, Slice(), Slice()));
            checkRatioArray(phaseRatio, chunk);
            *pha_[chunk] = phaseRatio;
        } else {
            *pha_[chunk] = phase(cArray(parSlice1, Slice(), Slice()));
        }
        (*pha_[chunk]) *= Float(180.0/C::pi);
    } else {
        throw(AipsError("phase has no meaning for this table"));
    }
    break;
  }

  case PMS::GREAL:
  case PMS::REAL: {
    if (parsAreComplex()) {
        Cube<Complex> cArray = cti.cparam();
        if (polnRatio_) {
            Array<Float> realRatio = real(cArray(parSlice1, Slice(), Slice()) / cArray(parSlice2, Slice(), Slice()));
            checkRatioArray(realRatio, chunk);
            *real_[chunk] = realRatio;
        } else {        
            *real_[chunk] = real(cArray(parSlice1, Slice(), Slice()));
        }
    } else {  // allow float for single dish cal tables
        Cube<Float> fArray = cti.fparam();
        if (polnRatio_) {
            Array<Float> ampRatio = fArray(parSlice1, Slice(), Slice()) / fArray(parSlice2, Slice(), Slice());
            checkRatioArray(ampRatio, chunk);
            *real_[chunk] = ampRatio;
        } else {        
            *real_[chunk] = fArray(parSlice1, Slice(), Slice());
        }
    }
    break;
  }

  case PMS::GIMAG:
  case PMS::IMAG: {
    if (parsAreComplex()) {
        Cube<Complex> cArray = cti.cparam();
        if (polnRatio_) {
            Array<Float> imagRatio = imag(cArray(parSlice1, Slice(), Slice()) / cArray(parSlice2, Slice(), Slice()));
            checkRatioArray(imagRatio, chunk);
            *imag_[chunk] = imagRatio;
        } else {        
            *imag_[chunk] = imag(cArray(parSlice1, Slice(), Slice()));
        }
    } else
        throw(AipsError("imag has no meaning for this table"));
    break;
  }

  case PMS::DELAY:{
    if (!parsAreComplex()) {
        Cube<Float> fArray = cti.fparam();
        if (polnRatio_) {
            Array<Float> delayRatio = fArray(parSlice1, Slice(), Slice()) - fArray(parSlice2, Slice(), Slice());
            checkRatioArray(delayRatio, chunk);
            *par_[chunk] = delayRatio;
        } else {
           *par_[chunk] = fArray(parSlice1, Slice(), Slice());
        }
    } else
        throw(AipsError( "delay has no meaning for this table"));
    break;
  }

  case PMS::OPAC: {
    if (!parsAreComplex() && calType_.contains("Opac")) {
        Cube<Float> fArray = cti.fparam();
        *par_[chunk] = fArray(parSlice1, Slice(), Slice());
    } else
        throw(AipsError( "opacity has no meaning for this table"));
    break;
  }

  case PMS::SWP: {   // "SPGAIN" in plotcal
    if ( !parsAreComplex() && calType_.contains("EVLASWPOW")) {
        Cube<Float> fArray = cti.fparam();
        if (polnRatio_) {
            Array<Float> swpRatio = fArray(parSlice1, Slice(), Slice()) / fArray(parSlice2, Slice(), Slice());
            checkRatioArray(swpRatio, chunk);
            *par_[chunk] = swpRatio;
        } else {
           *par_[chunk] = fArray(parSlice1, Slice(), Slice());
        }
    } else
        throw( AipsError( "SwPower has no meaning for this table"));
    break;
  }

  case PMS::TSYS: {
    if ( !parsAreComplex() &&
         (calType_.contains("EVLASWPOW") || calType_.contains("TSYS"))) {
        Cube<Float> fArray = cti.fparam();
        if (polnRatio_) {
            Array<Float> tsysRatio = fArray(parSlice1, Slice(), Slice()) / fArray(parSlice2, Slice(), Slice());
            checkRatioArray(tsysRatio, chunk);
            *par_[chunk] = tsysRatio;
        } else {
           *par_[chunk] = fArray(parSlice1, Slice(), Slice());
        }
    } else
        throw(AipsError( "Tsys has no meaning for this table"));
    break;
  }

  case PMS::SNR: {
    if (polnRatio_) {
        Array<Float> snrRatio = cti.snr()(parSlice1, Slice(), Slice()) / cti.snr()(parSlice2, Slice(), Slice());
        checkRatioArray(snrRatio, chunk);
        *snr_[chunk] = snrRatio;
    } else {
        *snr_[chunk] = cti.snr()(parSlice1, Slice(), Slice());
    }
    break;
  }

   case PMS::TEC: {
    if ( !parsAreComplex() && calType_[0]=='F') {
        Cube<Float> fArray = cti.fparam();
        *par_[chunk] = (fArray(parSlice1, Slice(), Slice()))/1e+16;
    } else
        throw(AipsError( "TEC has no meaning for this table"));
    break;
  }

  case PMS::FLAG:
    if (polnRatio_)
        *flag_[chunk] = cti.flag()(parSlice1, Slice(), Slice()) | cti.flag()(parSlice2, Slice(), Slice());
    else
        *flag_[chunk] = cti.flag()(parSlice1, Slice(), Slice());
    break;
 /*
    case PMS::WT: {
      *wt_[chunk] = cti.weightMat();
      break;
    }
    case PMS::AZ0:
    case PMS::EL0: {
      Vector<Double> azel;
      cti.azel0Vec(cti.time()(0),azel);
      az0_(chunk) = azel(0);
      el0_(chunk) = azel(1);
      break;
    }
    case PMS::HA0: 
      ha0_(chunk) = cti.hourang(cti.time()(0))*12/C::pi;  // in hours
      break;
    case PMS::PA0: {
      pa0_(chunk) = cti.parang0(cti.time()(0))*180.0/C::pi; // in degrees
      if (pa0_(chunk)<0.0) pa0_(chunk)+=360.0;
      break;
    }
 */
  case PMS::ANTENNA: {
    antenna_[chunk]->resize(nAnt_);
    indgen(*antenna_[chunk]);
    break;
  }

 /*
    case PMS::AZIMUTH:
    case PMS::ELEVATION: {
      Matrix<Double> azel;
      cti.azelMat(cti.time()(0),azel);
      *az_[chunk] = azel.row(0);
      *el_[chunk] = azel.row(1);
      break;
    }
    case PMS::PARANG:
      *parang_[chunk] = cti.feed_pa(cti.time()(0))*(180.0/C::pi);  // in degrees
      break;
    case PMS::ROW: {
      *row_[chunk] = cti.rowIds();
      break;
    }
 */
  case PMS::OBSERVATION: {
    (*obsid_[chunk]).resize(1);
    *obsid_[chunk] = cti.thisObs();
    break;
  }
  case PMS::INTENT: {
    throw(AipsError( "INTENT has no meaning for this table"));
    break;
  }
  default:
    throw(AipsError("Axis choice not supported for Cal Tables"));
    break;
  }
}
// ======================== end NewCalTable ==========================


// ======================== CalTable ==========================

// **** for BPOLY CalTable ****
void CalCache::loadBPoly(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // set up BPoly table from filename and load cache
  BJonesPolyTable ct = BJonesPolyTable(filename_);
  BJonesPolyTable selct;
  if (!selection_.isEmpty()) {
    Vector<Vector<Slice> > chansel;
    Vector<Vector<Slice> > corrsel;
    try {
      selection_.apply(ct, selct, chansel, corrsel);
    } catch (AipsError& err) {
      String msg(err.getMesg());
      if (msg == "No MS")
        throw(AipsError("Associated MS is not available, cannot plot solutions."));
      else
        throw(AipsError(msg));
    }
  }

  ROBJonesPolyMCol mainCol(ct);
  ROCalDescColumns calDescCol(ct);
  setMSname(calDescCol.msName()(0)); // throws error if not valid
  getNamesFromMS(); // field and antenna

  // count and load chunks
  Int nrow = ct.nRowMain(); // iterate per row to load cache
  countChunks(nrow, loadAxes, loadData, thread); // set up cache size
  loadCalChunks(mainCol, calDescCol, nrow, loadAxes, thread); 
}

// **** for GSPLINE CalTable ****
void CalCache::loadGSpline(vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  GJonesSplineTable ct = GJonesSplineTable(filename_);
  GJonesSplineTable selct;
  if (!selection_.isEmpty()) {
    Vector<Vector<Slice> > chansel;
    Vector<Vector<Slice> > corrsel;
    try {
      selection_.apply(ct, selct, chansel, corrsel);
    } catch (AipsError& err) {
      String msg(err.getMesg());
      if (msg == "No MS")
        throw(AipsError("Associated MS is not available, cannot plot solutions."));
      else
        throw(AipsError(msg));
    }
  }

  ROGJonesSplineMCol mainCol(ct);
  ROCalDescColumns calDescCol(ct);
  setMSname(calDescCol.msName()(0));  // throws error if not valid
  getNamesFromMS(); // field and antenna

  // count and load chunks
  Int nsample(1000); // make time samples to load cache
  countChunks(nsample, loadAxes, loadData, thread);
  loadCalChunks(mainCol, calDescCol, nsample, loadAxes, thread);
}

void CalCache::countChunks(Int nchunks, vector<PMS::Axis>& loadAxes,
    vector<PMS::DataColumn>& loadData, ThreadCommunication* thread) {
  // for CalTable
  if (thread!=NULL) {
    thread->setStatus("Establishing cache size.  Please wait...");
    thread->setAllowedOperations(false,false,false);
  }
  setCache(nchunks, loadAxes, loadData);
}

// **** for BPOLY CalTable ****
void CalCache::loadCalChunks(ROBJonesPolyMCol& mcol, ROCalDescColumns& dcol,
    Int nrow, const vector<PMS::Axis> loadAxes, ThreadCommunication* thread) {
  Slice parslice;
  setUpLoad(thread, parslice);

  // freq info from ms
  chanfreqs_.resize();
  try {
    getChanFreqsFromMS(); // set chanfreqs_; throws exception if no ms
  } catch (AipsError& err) {
  }

  MSMetaInfoForCal msmeta(msname_);
  BJonesPoly* bpoly = new BJonesPoly(msmeta);
  Record rec;  // for solving params
  rec.define("caltable", filename_);
  bpoly->setSpecify(rec);  // solves and makes data & flag cubes

  // each row of main table is a "chunk"
  for (Int row = 0; row < nrow; row++) {  // main loop to load axes
    // bpoly solutions are per spw
    Int calDescId = mcol.calDescId()(row);
    Vector<Int> spwIds = dcol.spwId()(calDescId);
    Int spw = spwIds(0);
    // nchan for chunk shapes
    Slicer chanslicer = Slicer(Slice(), Slice(spw));
    Int nChan = chanfreqs_(chanslicer).nelements();
    // use ant1 id for cube slicer (vis, flag, snr)
    Int ant1 = mcol.antenna1()(row);
    Slicer slicer(Slicer(parslice, Slice(), Slice(ant1)));

    // Cache the data shapes
    chshapes_(0,row) = parslice.length();
    chshapes_(1,row) = nChan;
    chshapes_(2,row) = 1;  // one row at a time
    chshapes_(3,row) = 1;  // one antenna per row
    goodChunk_(row) = True;

    // load axes for each row
    for(unsigned int i = 0; i < loadAxes.size(); i++) {
      PMS::Axis axis = loadAxes[i];
      if (PMS::axisIsData(axis)) {
        Cube<Complex> cpar, viscube;
        bpoly->solveAllCPar(spw, cpar);
        viscube = cpar(slicer);
        // get amp, phase, real, imag from viscube
        getCalDataAxis(axis, viscube, row);
      } else if (axis==PMS::FLAG) {
        Cube<Bool> parOK;
        bpoly->solveAllParOK(spw, parOK);
        // OK=true means flag=false
        *flag_[row] = !(parOK(slicer));
      } else if (axis==PMS::SNR) {
        Cube<Float> snrcube;
        bpoly->solveAllParSNR(spw, snrcube);
        *snr_[row] = snrcube(slicer);
      } else {
        loadCalAxis(mcol, dcol, row, axis);
      }
    }
    // If a thread is given, update it.
    if(thread != NULL) {
      double progress = ((double)row) / nrow;
      thread->setProgress((unsigned int)((progress * 100) + 0.5));
    }
  }
}
// **** for GSPLINE CalTable ****
void CalCache::loadCalChunks(ROGJonesSplineMCol& mcol, ROCalDescColumns& dcol,
    Int nsample, const vector<PMS::Axis> loadAxes, ThreadCommunication* thread) {
  // GSPLINE does not load per chunk or row as in other cal tables. 
  // Its "chunks" are per-time sample, generated from SPLINE_KNOTS_AMP/PHASE.
  // Therefore:
  //   Scalar columns are read once and all values plotted per time sample.
  //   Columns that normally have one value per chunk (field, spw, scan, time,
  //   and interval) will use the value in the column if all are the same; 
  //   else multiple values result in setting the value to -1.
  //   Otherwise locate does not work.

  Slice parslice, parslice2;
  setUpLoad(thread, parslice);
  Int nPol = parslice.length();  // for chunk shapes
  if (polnRatio_) 
    parslice2 = getParSlice(toVisCalAxis(PMS::AMP), "L");

  // load main table metadata once; use vector/value for every timestamp
  // field
  Vector<Int> fields(mcol.fieldId().getColumn());
  Int nItems = GenSort<Int>::sort(fields, Sort::Ascending, Sort::NoDuplicates);
  Int field = (nItems>1 ? -1 : fields(0));
  // spw
  Array<Int> spws(dcol.spwId().getColumn());
  nItems = GenSort<Int>::sort(spws, Sort::Ascending, Sort::NoDuplicates);
  Int spw = (nItems>1 ? -1 : spws(IPosition(2,0,0)));
  // scan
  Vector<Int> scans(mcol.scanNo().getColumn());
  nItems = GenSort<Int>::sort(scans, Sort::Ascending, Sort::NoDuplicates);
  Int scan = (nItems>1 ? -1 : scans(0));
  // antenna1
  Vector<Int> ant1(mcol.antenna1().getColumn());
  nItems = GenSort<Int>::sort(ant1, Sort::Ascending, Sort::NoDuplicates);
  ant1.resize(nItems, true);
  Int nAnt(nItems);
  // obsid
  Vector<Int> obsid(mcol.obsId().getColumn());
  nItems = GenSort<Int>::sort(obsid, Sort::Ascending, Sort::NoDuplicates);
  obsid.resize(nItems, true);
  // feed1
  Vector<Int> feed1(mcol.feed1().getColumn());
  nItems = GenSort<Int>::sort(feed1, Sort::Ascending, Sort::NoDuplicates);
  feed1.resize(nItems, true);
  // interval
  Vector<Double> intervals(mcol.interval().getColumn());
  nItems = GenSort<Double>::sort(intervals, Sort::Ascending, Sort::NoDuplicates);
  Double interval = (nItems>1 ? -1 : intervals(0));
  // poln
  Vector<Int> polns(parslice.length());
  if (polnRatio_) {
    polns = -1;
  } else {
    Int start(parslice.start()), inc(parslice.inc());
    indgen(polns, start, inc);
  }

  // set up gspline
  MSMetaInfoForCal msmeta(msname_);
  GJonesSpline* gspline = new GJonesSpline(msmeta);
  String mode(mcol.polyMode()(0));
  Record rec;  // for solving params
  rec.define("table", filename_);
  rec.define("apmode", mode);
  gspline->setApply(rec); // set up calbuffer for calcPar
  gspline->setSolve(rec); // set mode

  // values for setMeta: make sure ids are valid
  Int fieldId = (field>-1 ? field : 0);
  Int spwId = (spw>-1 ? spw : 0);
  Int obs(obsid(0));

  // get first row for chosen field
  uInt row;
  for (row=0; row<fields.size(); ++row)
    if (fields(row) == fieldId) break;

  MFrequency refFreq = mcol.refFreqMeas()(row)(IPosition(3,0,0,0));
  Double refFreqHz = refFreq.get("Hz").getValue();
  Vector<Double> freq(1, refFreqHz);    // for setMeta

  // Create 1000 time samples from spline knots: see PlotCal::virtualGSpline
  Vector<Double> splineKnots, times(nsample);
  // get splineKnots for that row
  if (mode.contains("AMP") || mode.contains("A&P")) {
    splineKnots = mcol.splineKnotsAmp()(row);
  } else if (mode.contains("PHAS") || mode.contains("A&P")) {
    splineKnots = mcol.splineKnotsPhase()(row);
  }
  // make samples based on spline knots
  Double dt((max(splineKnots)-min(splineKnots)) / Double(nsample));
  Double mintime(splineKnots(0) + (dt/2.0));
  for (Int sample=0; sample<nsample; sample++) {
    times(sample) = mintime + sample*dt;
  }

  // Now ready to load cache
  for (Int sample=0; sample<nsample; sample++) {
    // set time and field for calcPar
    Double time = times(sample);
    gspline->setMeta(obs, scan, time, spwId, freq, fieldId);
    gspline->doCalcPar();

    // Cache the data shapes
    chshapes_(0,sample) = nPol;
    chshapes_(1,sample) = 1; // nChan
    chshapes_(2,sample) = nAnt;
    chshapes_(3,sample) = 1; 
    goodChunk_(sample) = True;

    // load axes for each row
    for(unsigned int i = 0; i < loadAxes.size(); i++) {
      PMS::Axis axis = loadAxes[i];
      // slice viscube for poln selection, ant1
      Slicer slicer1(Slicer(parslice, Slice(), Slice()));
      if (PMS::axisIsData(axis)) { // amp, phase, real, imag 
        Cube<Complex> cpar, viscube;
        cpar = gspline->currCPar(); 
        viscube = cpar(slicer1);
        if (polnRatio_) {
          Slicer slicer2(Slicer(parslice2, Slice(), Slice()));
          viscube /= cpar(slicer2);
        }
        getCalDataAxis(axis, viscube, sample);
      } else {
        switch(axis) {
          case PMS::SCAN:
            scan_(sample) = scan;
            break;
          case PMS::FIELD:
            field_(sample) = field;
            break;
          case PMS::TIME: 
            time_(sample) = time;
            break;
          case PMS::TIME_INTERVAL:
            timeIntr_[sample] = interval;
            break;
          case PMS::SPW:
            spw_(sample) = spw;
            break;
          case PMS::CORR: {
            corr_[sample]->resize(nPol);
            String pol = selection_.corr();
            if (pol=="" || pol=="RL" || pol=="XY") {
              indgen(*corr_[sample]);
            } else if (pol=="/") {
              corr_[sample]->resize(1);
              corr_[sample]->set(-1);
            } else {  // R/X or L/Y
              Int poln = ((pol=="R" || pol=="X") ? 0 : 1);
              corr_[sample]->resize(1);
              corr_[sample]->set(poln);
            }
            break;
          }
          case PMS::ANTENNA1:
            *antenna1_[sample] = ant1;
            break;
          case PMS::FLAG: {
            Cube<Bool> parOK = gspline->currParOK();
            // OK=true means flag=false
            Cube<Bool> flagcube(!parOK(slicer1));
            if (polnRatio_) {
              Slicer slicer2(Slicer(parslice2, Slice(), Slice()));
              flagcube &= !parOK(slicer2);
            }
            *flag_[sample] = flagcube;
            break;
          }
          case PMS::ROW: {
            Vector<uInt> sampleRow(nAnt, sample);
            *row_[sample] = sampleRow;
            break;
          }
          case PMS::OBSERVATION:
            *obsid_[sample] = obsid;
            break;
          case PMS::FEED1:
            *feed1_[sample] = feed1;
            break;
          default:
            break;
          }
      }
    }

    // If a thread is given, update it.
    if(thread != NULL) {
      double progress = ((double)sample) / nsample;
      thread->setProgress((unsigned int)((progress * 100) + 0.5));
    }
  }
}

// **** for BPOLY CalTable ****
void CalCache::loadCalAxis(ROSolvableVisJonesMCol& mcol,
  ROCalDescColumns& dcol, Int chunk, PMS::Axis axis) {
    switch(axis) {
    case PMS::SCAN:
      scan_(chunk) = mcol.scanNo()(chunk);
      break;
    case PMS::FIELD:
      field_(chunk) = mcol.fieldId()(chunk);
      break;
    case PMS::TIME: 
      time_(chunk) = mcol.time()(chunk);
      break;
    case PMS::TIME_INTERVAL:
      timeIntr_(chunk) = mcol.interval()(chunk);
      break;
    case PMS::SPW: {
      Int calDescId = mcol.calDescId()(chunk);
      Vector<Int> spws = dcol.spwId()(calDescId);
      spw_(chunk) = spws(0);
      break;
    }
    case PMS::CHANNEL: {
      Int calDescId = mcol.calDescId()(chunk);
      Vector<Int> spwIds = dcol.spwId()(calDescId);
      Slicer slicer = Slicer(Slice(), Slice(spwIds(0)));
      Int nchanPerSpw = chanfreqs_(slicer).nelements();
      chan_[chunk]->resize(nchanPerSpw);
      indgen(*chan_[chunk]);
      break;
    }
    case PMS::FREQUENCY: {
      Int calDescId = mcol.calDescId()(chunk);
      Vector<Int> spwIds = dcol.spwId()(calDescId);
      Slicer slicer = Slicer(Slice(), Slice(spwIds(0)));
      *freq_[chunk] = chanfreqs_(slicer);
      // TBD: Convert freq to desired frame
      (*freq_[chunk]) /= 1.0e9; // in GHz
      break;
    }
    //case PMS::VELOCITY:
    case PMS::CORR: {
      corr_[chunk]->resize(chshapes_(0,chunk));
      String pol = selection_.corr();
      if (pol=="" || pol=="RL" || pol=="XY") {
        indgen(*corr_[chunk]);
      } else {
        Int poln = ((pol=="R" || pol=="X") ? 0 : 1);
        corr_[chunk]->resize(1);
        corr_[chunk]->set(poln);
      }
      break;
    }
    case PMS::ANTENNA1: { // holds a Vector of antenna ids
      Vector<Int> ant1(1, mcol.antenna1()(chunk));
      *antenna1_[chunk] = ant1;
      break;
    }
    case PMS::ANTENNA:
      throw(AipsError( "Antenna has no meaning for this table; use Antenna1"));
      break;
    case PMS::ANTENNA2:
      throw(AipsError( "Antenna2 has no meaning for this table"));
      break;
    case PMS::BASELINE:
      throw(AipsError( "Baseline has no meaning for this table"));
      break;
    case PMS::GAMP:
    case PMS::AMP: {
      // handled in getCalDataAxis
      break;
    }
    case PMS::GPHASE:
    case PMS::PHASE: {
      // handled in getCalDataAxis
      break;
    }
    case PMS::GREAL:
    case PMS::REAL: {
      // handled in getCalDataAxis
      break;
    }
    case PMS::GIMAG:
    case PMS::IMAG: {
      // handled in getCalDataAxis
      break;
    }
    case PMS::DELAY:
      throw(AipsError( "Delay has no meaning for this table"));
      break;
    case PMS::OPAC:
      throw(AipsError( "Opac has no meaning for this table"));
      break;
    case PMS::SWP:   // "SPGAIN" in plotcal
      throw(AipsError( "SwPower has no meaning for this table"));
      break;
    case PMS::TSYS:
      throw(AipsError( "Tsys has no meaning for this table"));
      break;
    case PMS::SNR: {
      // handled in loadCalChunks for BPOLY
      break;
    }
    case PMS::TEC:
      throw(AipsError( "TEC has no meaning for this table"));
      break;
    case PMS::FLAG: {
      // handled in loadCalChunks
      break;
    }
    //case PMS::WT: 
    //case PMS::AZ0:
    //case PMS::EL0:
    //case PMS::HA0: 
    //case PMS::PA0: 
    //case PMS::AZIMUTH:
    //case PMS::ELEVATION:
    //case PMS::PARANG:
    case PMS::ROW: {
      Vector<uInt> rows(1, chunk);
      *row_[chunk] = rows;
      break;
    }
    case PMS::OBSERVATION: {
      Vector<Int> obsIds(1, mcol.obsId()(chunk));
      *obsid_[chunk] = obsIds;
      break;
    }
    case PMS::INTENT: 
      throw(AipsError( "INTENT has no meaning for this table"));
      break;
    case PMS::FEED1: { 
      Vector<Int> feedIds(1, mcol.feed1()(chunk));
      *feed1_[chunk] = feedIds;
      break;
    }
    default:
      throw(AipsError("Axis choice not supported for Cal Tables"));
      break;
    } // switch
}

// **** for GSPLINE CalTable ****
void CalCache::checkAxes(const vector<PMS::Axis>& loadAxes) {
  // trap user-requested axes that are invalid for GSPLINE
  for(unsigned int i = 0; i < loadAxes.size(); i++) {
    PMS::Axis axis(loadAxes[i]);
    switch (axis) {
      case PMS::SCAN:
      case PMS::FIELD:
      case PMS::TIME:
      case PMS::SPW:
      case PMS::TIME_INTERVAL:
      case PMS::CORR:
      case PMS::ANTENNA1:
      case PMS::AMP:
      case PMS::GAMP:
      case PMS::PHASE:
      case PMS::GPHASE:
      case PMS::REAL:
      case PMS::GREAL:
      case PMS::IMAG:
      case PMS::GIMAG:
      case PMS::FLAG:
      case PMS::OBSERVATION:
      case PMS::FEED1: {
        // allowed
        break;
      }
      case PMS::ANTENNA:
      case PMS::CHANNEL:
      case PMS::FREQUENCY:
      case PMS::SNR: {
        String msg("GSPLINE plotting does not support " + PMS::axis(axis) + " axis");
        if (axis==PMS::ANTENNA) msg += "; use ANTENNA1";
        throw(AipsError(msg));
        break;
      }
      case PMS::ANTENNA2:
      case PMS::BASELINE:
      case PMS::ROW:
      case PMS::DELAY:
      case PMS::OPAC:
      case PMS::SWP:
      case PMS::TSYS:
      case PMS::TEC:
      case PMS::INTENT: {
        String msg(PMS::axis(axis) + " has no meaning for this table");
        throw(AipsError(msg));
        break;
      }
      default:
        throw(AipsError("Axis choice not supported for Cal Tables"));
        break;
    }
  }
}

void CalCache::setMSname(String msname) {
  // set msname_ (with path) if valid
  Path filepath(filename_);
  String path(filepath.dirName()), mspath("");
  if (!path.empty()) mspath += path + "/";
  mspath += msname;
  if (msname.empty() || !Table::isReadable(mspath))
    throw(AipsError("Associated MS is not available, cannot plot solutions."));
  msname_ = mspath;
}

void CalCache::getNamesFromMS() {
  // Set antenna and field names for Locate.
  MSMetaInfoForCal msmeta(msname_);
  msmeta.antennaNames(antnames_);
  msmeta.fieldNames(fldnames_);
  nAnt_ = msmeta.nAnt();
}

void CalCache::setUpLoad(ThreadCommunication* thread, Slice& parSlice) {
  // common setup for CalTables
  // permit cancel in progress meter:
  if(thread != NULL)
    thread->setAllowedOperations(false,false,true);
  logLoad("Loading chunks......");

  // initial chunk info
  chshapes_.resize(4,nChunk_);
  goodChunk_.resize(nChunk_);
  goodChunk_.set(False);

  // get selected npol
  String polSelection(selection_.corr());
  String paramAxis(toVisCalAxis(PMS::AMP));
  // getParSlice() checks for valid axis and pol sel
  if (polnRatio_) {  // just pick one for length
    parSlice = getParSlice(paramAxis, "R");
  } else { 
    parSlice = getParSlice(paramAxis, polSelection);
  }
}

void CalCache::getChanFreqsFromMS() {
  // shape is (nchan, nspw)
  MeasurementSet ms(msname_);
  ROMSColumns mscol(ms);
  chanfreqs_ = mscol.spectralWindow().chanFreq().getColumn();
}

void CalCache::getCalDataAxis(PMS::Axis axis, Cube<Complex>& viscube,
    Int chunk) {
  // Get axes derived from calculated data cube; 
  // poln selection (parSlice) already applied to viscube
  switch(axis) {
    case PMS::GAMP:
    case PMS::AMP: {
      *amp_[chunk] = amplitude(viscube);
      break;
    }
    case PMS::GPHASE:
    case PMS::PHASE: {
      *pha_[chunk] = phase(viscube);
      (*pha_[chunk]) *= Float(180.0/C::pi);
      break;
    }
    case PMS::GREAL:
    case PMS::REAL: {
      *real_[chunk] = real(viscube);
      break;
    }
    case PMS::GIMAG:
    case PMS::IMAG: {
      *imag_[chunk] = imag(viscube);
      break;
    }
    default:
      throw(AipsError("Axis choice not supported for Cal Tables"));
      break;
  }
}

// ======================== end CalTable ==========================

void CalCache::flagToDisk(const PlotMSFlagging& flagging,
    Vector<Int>& flchunks, Vector<Int>& flrelids,
    Bool flag, PlotMSIndexer* indexer, int dataIndex ) {
  
  // Sort the flags by chunk:
  Sort sorter;
  sorter.sortKey(flchunks.data(),TpInt);
  sorter.sortKey(flrelids.data(),TpInt);
  Vector<uInt> order;
  uInt nflag;
  nflag = sorter.sort(order,flchunks.nelements());

  stringstream ss;

  // Make the VisIterator writable, with selection revised as appropriate
  Bool selectchan(netAxesMask_[dataIndex](1) && !flagging.channel());
  Bool selectcorr(netAxesMask_[dataIndex](0) && !flagging.corrAll());

  // Set up a write-able CTIter (ci_p also points to it)
  setUpCalIter(filename_,selection_,False,selectchan,selectcorr);
  ci_p->reset();

  Int iflag(0);
  for (Int ichk=0;ichk<nChunk_;++ichk) {
    if (ichk!=flchunks(order[iflag]) && !ci_p->pastEnd())
      // nothing to flag this chunk, just advance
      ci_p->next();
    else {
      // This chunk requires flag-setting
      Int ifl(iflag);
      
      // Get bits we need from the table
      Cube<Bool> ctflag;
      Vector<Int> channel,a1,a2;
      ci_p->flag(ctflag);
      ci_p->chan(channel);
      ci_p->antenna1(a1);
      ci_p->antenna2(a2);

      // Apply poln selection
      Int npar;
      String pol = selection_.corr();
      if (pol=="" || pol=="RL" || pol=="XY" || pol=="/") { // both axes
        npar = ctflag.shape()(0);
      } else { // poln selection using calParSlice
        String paramAxis = toVisCalAxis(PMS::FLAG);
        npar = getParSlice(paramAxis, pol).length();
      }    
      Int nchan = channel.nelements();
      Int nrow = ci_p->nrow();

      if (True) {
        Int currChunk=flchunks(order[iflag]);
        Double time=getTime(currChunk,0);
        Double cttime=ci_p->time()(0);
        Int spw=Int(getSpw(currChunk,0));
        Int ctspw=ci_p->thisSpw();
        Int field=Int(getField(currChunk,0));
        Int ctfld=ci_p->thisField();
        ss << "Time diff:  " << time-cttime << " " << time  << " " << cttime << "\n";
        ss << "Spw diff:   " << spw-ctspw   << " " << spw   << " " << ctspw  << "\n";
        ss << "Field diff: " << field-ctfld << " " << field << " " << ctfld  << "\n";
      }

      // Apply all flags in this chunk to this VB
      ifl=iflag;
      while (ifl<Int(nflag) && flchunks(order[ifl])==ichk) {
        Int currChunk=flchunks(order[ifl]);
        Int irel=flrelids(order[ifl]);
        Slice par1,chan,bsln;
        Slice par2 = Slice();

        // Set flag range on par axis:
        if (netAxesMask_[dataIndex](0) && !flagging.corrAll()) {
          // A specific single par
          if (pol=="" || pol=="RL" || pol=="XY") {  // flag both axes
            Int ipar=indexer->getIndex1000(currChunk,irel);
            par1 = Slice(ipar,1,1);
          } else if (polnRatio_) {
            par1 = getParSlice(toVisCalAxis(PMS::AMP), "R");
            par2 = getParSlice(toVisCalAxis(PMS::AMP), "L");
          } else {
            par1 = getParSlice(toVisCalAxis(PMS::AMP), pol);
          }
        } else {
          // all on par axis
          par1 = Slice(0,npar,1);
        }

        // Set Flag range on channel axis:
        if (netAxesMask_[dataIndex](1) && !flagging.channel()) {
          // A single specific channel
          Int ichan=indexer->getIndex0100(currChunk,irel);
          chan=Slice(ichan,1,1);
        } else {
          // Extend to all channels
          chan=Slice(0,nchan,1);
        }

        // Set Flags on the baseline axis:
        Int thisA1=Int(getAnt1(currChunk,indexer->getIndex0010(currChunk,irel)));
        Int thisA2=Int(getAnt2(currChunk,indexer->getIndex0010(currChunk,irel)));
        if (netAxesMask_[dataIndex](2) &&
            !flagging.antennaBaselinesBased() &&
            thisA1>-1 ) {
          // i.e., if baseline is an explicit data axis,
          //       full baseline extension is OFF
          //       and the first antenna in the selected point is > -1
          // Do some variety of detailed per-baseline flagging
          for (Int irow=0;irow<nrow;++irow) {
            if (thisA2>-1) {
              // match a baseline exactly
              if (a1(irow)==thisA1 &&
                  a2(irow)==thisA2) {
                ctflag(par1,chan,Slice(irow,1,1)) = flag;
                if (par2.length() > 0)
                  ctflag(par2,chan,Slice(irow,1,1)) = flag;
                break;  // found the one baseline, escape from for loop
              }
            } else {
              // either antenna matches the one specified antenna
              //  (don't break because there will be more than one)
              if (a1(irow)==thisA1 ||
                  a2(irow)==thisA1) {
                ctflag(par1,chan,Slice(irow,1,1)) = flag;
                if (par2.length() > 0)
                  ctflag(par2,chan,Slice(irow,1,1)) = flag;
              }
            }
          }
        } else {
          // Set flags for all baselines, because the plot
          //  is ordinarily implicit in baseline, we've turned on baseline
          //  extension, or we've avaraged over all baselines
          bsln=Slice(0,nrow,1);
          ctflag(par1,chan,bsln) = flag;
          if (par2.length() > 0)
            ctflag(par2,chan,bsln) = flag;
        }

      ++ifl;
      } // while
      
      // Put the flags back into the MS
      wci_p->setflag(ctflag);
      
      // Advance to the next vb
      if (!ci_p->pastEnd())
        ci_p->next();
      else
        // we are done, so escape chunk loop
        break;

      // step over the flags we've just done
      iflag=ifl;
      
      // Escape if we are already finished
      if (uInt(iflag)>=nflag) break;
    } // flaggable chunk
  } // ichk

  // Delete the VisIter so lock is released
  if (wci_p)
    delete wci_p;
  wci_p=NULL;
  ci_p=NULL;
  logFlag(ss.str());
}

String CalCache::toVisCalAxis(PMS::Axis axis) {
    switch (axis) {
        // FLAG and SNR have same shape as AMP 
        // and should be sliced the same way
        case PMS::AMP:
        case PMS::GAMP:
        case PMS::FLAG:
        case PMS::SNR:
            if (calType_.contains("EVLASWP")) return "GAINAMP";
            if (calType_.contains("TSYS")) return "TSYS";
            if (calType_[0] == 'K') return "DELAY";
            if (calType_[0] == 'F') return "TEC";
            if (calType_ == "TOpac") return "OPAC";
            return "AMP";
            break;
        case PMS::PHASE:
        case PMS::GPHASE:
            return "PHASE";
            break;
        case PMS::REAL:
        case PMS::GREAL:
            return "REAL";
            break;
        case PMS::IMAG:
        case PMS::GIMAG:
            return "IMAG";
            break;
        default:
            return PMS::axis(axis);
            break;
    }
}

Slice CalCache::getParSlice(String axis, String polnSel) {
    Slice parSlice = Slice();
    try {
        parSlice = viscal::calParSlice(filename_, axis, polnSel);
    } catch(AipsError& err) {
        if (err.getMesg().contains("Unsupported value type")) {
            // Message a bit vague at top level, add some explanation
            String errMsg = err.getMesg() + ". Invalid axis or polarization selection for cal table type.";
            throw(AipsError(errMsg));
        } else { // unsupported cal type
            throw(AipsError(err));
        }
    }
    return parSlice;
}

void CalCache::setFilename(String filename) { 
    filename_ = filename;
    NewCalTable ct(NewCalTable::createCT(filename_,Table::Old,Table::Plain));
    calType_= ct.tableInfo().subType();
}

void CalCache::checkRatioArray(Array<Float>& array, Int chunk) {
    Cube<Float> ratioCube;
    ratioCube.reference(array);
    Cube<Bool> flags;
    flags.reference(*flag_[chunk]);

    IPosition cubeShape = ratioCube.shape();
    for (uInt i=0; i<cubeShape[0]; ++i) {
        for (uInt j=0; j<cubeShape[1]; ++j) {
            for (uInt k=0; k<cubeShape[2]; ++k) {
                if (isInf(ratioCube(i,j,k))) {
                    ratioCube(i,j,k) = 1.0;
                    flags(i,j,k) = True;
                    divZero_ = True;
                }
            }
        }
    }
}

}
