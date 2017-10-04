//# FringeJones.cc: Implementation of FringeJones
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003,2011
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

#include <synthesis/MeasurementComponents/FringeJones.h>

#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/VisBuffAccumulator.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <synthesis/MeasurementEquations/VisEquation.h>  // *
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
// FIXME: ?
#include <casa/Arrays/ArrayLogical.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/System/Aipsrc.h>

#include <casa/sstream.h>

#include <measures/Measures/MCBaseline.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasTable.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>

#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/MaskArrMath.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_multilarge_nlinear.h>

#define DEVDEBUG 0

using namespace casa::vi;
using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN
    
static void unitize(Array<Complex>& vC) 
{
    Array<Float> vCa(amplitude(vC));
    // Divide by non-zero amps
    vCa(vCa<FLT_EPSILON)=1.0;
    vC /= vCa;
}


SDBListGridManager::SDBListGridManager(SDBList& sdbs_) :
    sdbs (sdbs_)
{
    std::set<Double> fmaxes;
    std::set<Double> fmins;
    Float dfn;
    Int totalChans0( 0 ) ;
    Int nchan;

    
    for (Int i=0; i != sdbs.nSDB(); i++) {
        SolveDataBuffer& sdb = sdbs(i);
        Int spw = sdb.spectralWindow()(0);
        Double t = sdbs(i).time()(0);
        times.insert( t ); 
        if ( spwins.find( spw ) == spwins.end() ) {
            spwins.insert(spw);
            const Vector<Double>& fs = sdb.freqs();
            spwIdToFreqMap[spw] = &(sdb.freqs());
            nchan = sdb.nChannels();
            fmaxes.insert(fs(nchan-1));
            fmins.insert(fs(0));
            // We assume they're all at the same time.

            totalChans0 += nchan;
            Float df0 = fs(1) - fs(0);
            dfn = (fs(nchan-1) - fs(0))/(nchan-1);
            if (DEVDEBUG) {
                cerr << "Spectral window " << spw << " has " << nchan << " channels" << endl;
                cerr << "df0 "<< df0 << "; " << "dfn " << dfn << endl;
            }
        } else {
            continue;
        }
    }
    nt = sdbs.nSDB()/spwins.size();
    tmin = *(times.begin());
    tmax = *(times.rbegin());
    dt = (tmax - tmin) / (nt - 1);
    nSPWChan = nchan;
    fmin = *(fmins.begin());
    fmax = *(fmaxes.rbegin());
    totalChans = round((fmax - fmin)/dfn + 1);
    df = (fmax - fmin)/(totalChans-1);
    if (DEVDEBUG) {
        cerr << "Global fmin " << fmin << " global max " << fmax << endl;
        cerr << "nt " << nt << " dt " << dt << endl;
        cerr << "tmin " << tmin << " tmax " << tmax << endl;
        cerr << "Global df " << df << endl;
        cerr << "I guess we'll need " << totalChans << " freq points in total." << endl;
        cerr << "Compared to " << totalChans0 << " with simple-minded concatenation." << endl;
    }
}

void
SDBListGridManager::checkAllGridpoints() {
        map<Int , Vector<Double> const * >::iterator it;
        for (it = spwIdToFreqMap.begin(); it != spwIdToFreqMap.end(); it++) {
            Int spwid = it->first;
            Vector<Double> const* fs = it->second;
            Int length;
            fs->shape(length);
            for (Int i=0; i!=length; i++) {
                Double f = (*fs)(i);
                Int j = bigFreqGridIndex(f);
                cerr << "spwid, i = (" << spwid << ", " << i << ") => " << j << " (" << f << ")" << endl;
            }
        }
        cerr << "[1] spwins.size() " << nSPW() << endl;
        cerr << "[2] spwins.size() " << spwins.size() << endl;
    }

Int
SDBListGridManager::swStartIndex(Int spw) {
    Vector<Double> const* fs = spwIdToFreqMap[spw];
    Double f0 = (*fs)(0);
    return bigFreqGridIndex( f0 );
}
    


// **************************************************************
// DelayRateFFT modeled on DelayFFT(const VisBuffer&, Double padBW, Int refant) in KJones.{cc|h}

// Remark: Double, Int et al are in casacore namespace. Need to be
// qualified in headers, I guess, but here we are not in headers and we
// are using that namespace.



DelayRateFFT::DelayRateFFT(SDBList& sdbs, Int refant) :
    refant_( refant ),
    gm_ ( sdbs ),
    nPadFactor_ ( max(2, 8  / gm_.nSPW() )), 
    nt_( gm_.nt ),
    nPadT_( nPadFactor_ * nt_ ),
    nChan_ ( gm_.nChannels() ),
    nPadChan_ ( nPadFactor_*nChan_ ),
    dt_ ( gm_.dt ),
    f0_( gm_.fmin / 1.e9),      // GHz
    df_( gm_.df / 1.e9),
    Vpad_(),
    xcount_(),
    sumw_(),
    sumww_() {
    // Actual code!
    // gm_.checkAllGridpoints();
    Int nCorrOrig( sdbs(0).nCorrelations() );
    nCorr_ = (nCorrOrig> 1 ? 2 : 1); // number of p-hands

    for (Int i=0; i<nCorr_; i++) {
        activeAntennas_[i].insert(refant_);
    }
    
    // when we get the visCubecorrected it is already
    // reduced to parallel hands, but there isn't a
    // corresponding method for flags.
    Int corrStep = (nCorrOrig > 2 ? 3 : 1); // step for p-hands

    SolveDataBuffer& s0 ( sdbs(0) );
    nElem_ =  1 + max( max(s0.antenna1()), max(s0.antenna2())) ;
    //

    IPosition aggregateDim(2, nCorr_, nElem_);
    xcount_.resize(aggregateDim);
    sumw_.resize(aggregateDim);
    sumww_.resize(aggregateDim);

    xcount_ = 0;
    sumw_ = 0.0;
    sumww_ = 0.0;
        
    // Can't get timeInterval, fails with error
    // "Caught exception: Exception: Can't fill VisBuffer component TimeInterval:
    // Not attached to VisibilityIterator."

    if (DEVDEBUG) {
        cerr << "Filling FFT grid with " << sdbs.nSDB() << " data buffers." << endl;
    }
    if (sdbs.nSDB() < 2) {
        throw(AipsError("Not enough sdbs!"));
    }
        
    IPosition paddedDataSize(4, nCorr_, nElem_, nPadT_, nPadChan_);
    Vpad_.resize(paddedDataSize);

    // FIXME: There are now multiple SDBs per time step!
    for (Int ibuf=0; ibuf != sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s ( sdbs(ibuf) );
        if ( !s.Ok() )
            continue;

        Int nr = 0;
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if ( s.flagRow()(irow) )
                continue;
            Int iant;
            Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
            if (a1 == a2) {
                continue;
            }
            else if (a1 == refant_) {
                iant = a2;
            }
            else if (a2 == refant_) {
                iant = a1;
            }
            else {
                continue;
            }
            // OK, we're not skipping this one so we have to do something.

            // v has shape (nelems, ?, nrows, nchannels)
            Cube<Complex> v = s.visCubeCorrected();
            Cube<Float> w = s.weightSpectrum();
            Cube<Bool> fl = s.flagCube();
            Int spw = s.spectralWindow()( 0 );
            Int f_index = gm_.swStartIndex( spw );    // ditto!
            Int t_index = gm_.getTimeIndex( s.time()(0) );
            Int spwchans = gm_.nSPWChan;
            IPosition start( 4,      0,      iant, t_index, f_index);
            IPosition stop(  4,      nCorr_,    1,       1, spwchans);
            IPosition stride(4,      1,         1,       1, 1);
            Slicer sl1(start,     stop, stride, Slicer::endIsLength);
            Slicer sl2(IPosition(3, 0,         0, irow),
                       IPosition(3, nCorr_, spwchans, 1),
                       IPosition(3, corrStep,        1,  1), Slicer::endIsLength);
                
            Slicer flagSlice(IPosition(3, 0,         0, irow),
                             IPosition(3, nCorr_, spwchans, 1),
                             IPosition(3, corrStep,        1, 1), Slicer::endIsLength);
            nr++;
            if ( DEVDEBUG && 0) {
                cerr << "nr " << nr
                     << " irow " << endl
                     << "Vpad shape " << Vpad_.shape() << endl
                     << "v shape " << v.shape() << endl
                     << "sl2 " << sl2 << endl
                     << "sl1 " << sl1 << endl
                     << "flagSlice " << flagSlice << endl;
            }
            if ( DEVDEBUG ) {
                cerr << "SolveDataBuffer[" << ibuf << "] has spw=" << spw << "; freq, time indices ("
                     << f_index << ", " << t_index << ")" << endl;
            }
            Array<Complex> rhs = v(sl2).nonDegenerate();
            Array<Float> weights = w(sl2).nonDegenerate();
                
            unitize(rhs);
            Vpad_(sl1).nonDegenerate() = rhs * weights;
            // Zero flagged entries.

            Array<Bool> flagged( fl(flagSlice).nonDegenerate() );
                
            if ( allTrue(flagged) ) {
                ; // cerr << "irow " << irow << " Whoopsie!" << endl;
            } else {
                for (Int icorr=0; icorr<nCorr_; ++icorr) {
                    IPosition p(2, icorr, iant);
                    activeAntennas_[icorr].insert(iant);
                    for (size_t ichan=0; ichan != spwchans+1; ichan++) {
                        IPosition pchan(2, ichan, irow);
                        if (!flagged(pchan)) {
                            Float w = weights(pchan);
                            xcount_(p)++;
                            sumw_(p) += w;
                            sumww_(p) += w*w;
                        }
                    }
                }
            }                
            if (DEVDEBUG && 0) {
                cerr << "flagSlice " << flagSlice << endl
                     << "fl.shape() " << fl.shape() << endl
                     << "Vpad_.shape() " << Vpad_.shape() << endl
                     << "flagged.shape() " << flagged.shape() << endl
                     << "sl1 " << sl1 << endl;
            }
            Vpad_(sl1).nonDegenerate()(flagged) = Complex(0.0);
        }
    }
    if (DEVDEBUG) {
        cerr << "Constructed a DelayRateFFT object." << endl;
    }
}

void DelayRateFFT::printActive() {
    for (Int icorr=0; icorr != nCorr_; icorr++) {
        cerr << "Antennas found for correlation " << icorr << ": ";
        for (std::set<Int>::iterator it = activeAntennas_[icorr].begin(); it != activeAntennas_[icorr].end(); it++) {
            cerr << *it << ", ";
        }
        cerr << endl;
    }
    cerr << endl;
}    

void DelayRateFFT::FFT() {
    // Axes are 0: correlation (i.e., hand of polarization), 1: antenna, 2: time, 3: channel
    Vector<Bool> ax(4, false);
    ax(2) = true;
    ax(3) = true;
    // Also copied from DelayFFT in KJones.
    ArrayLattice<Complex> c(Vpad_);
    // But variable c is not returned and is not a member?
    // IT SEEMS that the FFT is in place, and that FFTing c actually changes Vpad_
    LatticeFFT::cfft0(c, ax, true);
}

std::pair<Bool, Float>
DelayRateFFT::xinterp(Float alo, Float amax, Float ahi) {
    Float denom (alo-2.0*amax+ahi);
    Bool cond = amax>0.0 && abs(denom)>0.0 ;
    Float fpk = cond ? 0.5-(ahi-amax)/denom : 0.0;
    return std::make_pair(cond, fpk);
}
    
void
DelayRateFFT::searchPeak() {
    // Recall param_ -> [phase, delay, rate] for each correlation
    param_.resize(3*nCorr_, nElem_); // Srsly, why we do this here?
    param_.set(0.0);
    // Note: It looks like we have to have one flag per parameter?
    flag_.resize(3*nCorr_, nElem_);
    flag_.set(true);  // all flagged initially
    if (DEVDEBUG) {
        cerr << "nt_ " << nt_ << " nPadChan_ " << nPadChan_ << endl;
    }
    for (Int icorr=0; icorr<nCorr_; ++icorr) {
        flag_(icorr*3 + 0, refant()) = false; 
        flag_(icorr*3 + 1, refant()) = false;
        flag_(icorr*3 + 2, refant()) = false;
        for (Int ielem=0; ielem<nElem_; ++ielem) {
            if (ielem==refant()) {
                continue;
            }
            // NB: Time, Channel
            // And once again we fail at slicing
            IPosition start(4, icorr, ielem,      0,         0);
            IPosition stop( 4,     1,     1, nPadT_, nPadChan_);
            IPosition step( 4,     1,     1,       1,        1);
            Slicer sl(start, stop, step, Slicer::endIsLength);
            Matrix<Complex> aS = Vpad_(sl).nonDegenerate();
            // cerr << "aS.shape()=" <<aS.shape() << endl;

            Matrix<Float> amp( amplitude(aS) );
            Int ipkch(0);
            Int ipkt(0);
            Float amax(-1.0);
            // Unlike KJones we have to iterate in time too
            for (Int itime=0; itime != nPadT_; itime++) {
                for (Int ich=0; ich != nPadChan_; ich++) {
                    if (amp(itime, ich) > amax) {
                        ipkch = ich;
                        ipkt  = itime;
                        amax=amp(itime, ich);
                    }
                }
            }
            // Finished grovelling. Now we have the location of the
            // maximum amplitude.
            Float alo_ch = amp(ipkt, (ipkch > 0) ? ipkch-1 : nPadChan_-1);
            Float ahi_ch = amp(ipkt, ipkch<(nPadChan_-1) ? ipkch+1 : 0);
            // cerr << "In channel dimension ipkch " << ipkch << " alo " << alo_ch  << " amax " << amax << " ahi " << ahi_ch << endl;
            std::pair<Bool, Float> maybeFpkch = xinterp(alo_ch, amax, ahi_ch);
            // We handle wrapping while looking for neighbours
            Float alo_t = amp(ipkt > 0 ? ipkt-1 : nPadT_ -1,     ipkch);
            Float ahi_t = amp(ipkt < (nPadT_ -1) ? ipkt+1 : 0,   ipkch);
            // cerr << "In time dimension ipkt " << ipkt << " alo " << alo_t  << " amax " << amax << " ahi " << ahi_t << endl;
            std::pair<Bool, Float> maybeFpkt = xinterp(alo_t, amax, ahi_t);

            Int sgn = (ielem < refant()) ? 1 : -1;
            if (maybeFpkch.first and maybeFpkt.first) {
                // Phase
                Complex c = aS(ipkt, ipkch);
                Float phase = arg(c);
                // FIXME: Here we go again.
                param_(icorr*3 + 0, ielem) = sgn*phase;
                // Delay
                // FIXME!:
                // Float delay = (ipkch + maybeFpkch.second)/Float(nPadChan_);
                    
                Float delay = (ipkch)/Float(nPadChan_);
                if (delay > 0.5) delay -= 1.0;           // fold
                delay /= df_;                           // nsec
                param_(icorr*3 + 1, ielem) = sgn*delay; //
                // FIXME!:
                // Double rate = (ipkt + maybeFpkt.second)/Float(nPadT_);
                Double rate = (ipkt)/Float(nPadT_);
                if (rate > 0.5) rate -= 1.0;
                Double rate0 = rate/dt_;
                Double rate1 = rate0/(1e9 * f0_); 

                param_(icorr*3 + 2, ielem) = Float(sgn*rate1); 
                if (DEVDEBUG) {
                    cerr << "maybeFpkch.second=" << maybeFpkch.second
                         << ", df_ " << df_ 
                         << " fpkch " << (ipkch + maybeFpkch.second) << endl;
                    cerr << " maybeFpkt.second=" << maybeFpkt.second
                         << " rate0 " << rate
                         << " 1e9 * f0_ " << 1e9 * f0_ 
                         << ", dt_ " << dt_
                         << " fpkt " << (ipkt + maybeFpkt.second) << endl;
                        
                }
                if (DEVDEBUG) {
                    cerr << "Found peak for element " << ielem << " correlation " << icorr
                         << " ipkt=" << ipkt << "/" << nPadT_ << ", ipkch=" << ipkch << "/" << nPadChan_
                         << " peak=" << amax 
                         << "; delay " << delay << ", rate " << rate
                         << ", phase " << arg(c) << " sign= " << sgn << endl;
                }
                // Set 3 flags.
                flag_(icorr*3 + 0, ielem)=false; 
                flag_(icorr*3 + 1, ielem)=false;
                flag_(icorr*3 + 2, ielem)=false;
            }
            else {
                if (DEVDEBUG) {
                    cerr << "No peak in 2D FFT for element " << ielem << " correlation " << icorr << endl;
                }
            }
        }
    }
}


Float
DelayRateFFT::snr(Int icorr, Int ielem, Float delay, Float rate) {
    // Have to convert delay and rate back into indices on the padded 2D grid.
    Int sgn = (ielem < refant()) ? 1 : -1;
    delay *= sgn*df_;
    if (delay < 0.0) delay += 1;
    Int ichan = Int(delay*nPadChan_ + 0.5); 
    if (ichan == nPadChan_) ichan = 0;
        
    rate *= sgn*1e9 * f0_;
    rate *= dt_;
    if (rate < 0.0) rate += 1;
    Int itime = Int(rate*nPadT_ + 0.5);
    if (itime == nPadT_) itime = 0;
    // FIXME!:
    // Double rate = (ipkt + maybeFpkt.second)/Float(nPadT_);
        
    // What about flags? If the datapoint closest to the computed
    // delay and rate values is flagged we probably shouldn't use
    // it, but what *should* we use?
    IPosition ipos(4, icorr, ielem, itime, ichan);
    IPosition p(2, icorr, ielem);
    Complex v = Vpad_(ipos);
    Float peak = abs(v);
    if (peak > 0.999*sumw_(p)) peak=0.999*sumw_(p);
    // xcount is number of data points for baseline to ielem
    // sumw is sum of weights,
    // sumww is sum of squares of weights

    Float x = C::pi/2*peak/sumw_(p);
    Float cwt = (pow(tan(x), 1.163) * sqrt(sumw_(p)/sqrt(sumww_(p)/xcount_(p))));

    if (DEVDEBUG) {
        cerr << "Correlation " << icorr << " antenna " << ielem << " ipos " << ipos
             << " peak=" << peak << "; xang=" << x << "; xcount=" << xcount_(p) << "; sumw=" << sumw_(p) << "; sumww=" << sumww_(p) 
             << " snr " << cwt << endl;
    }
    return cwt;
}
    

void DelayRateFFT::removeAntennasCorrelation(Int icor, std::set< Int > s) {
    std::set< Int > & as = activeAntennas_.find(icor)->second;
    for (std::set< Int >::iterator it=s.begin(); it!=s.end(); it++) {
        as.erase(*it);
    }
}


// Start of GSL compliant solver
// This function is supposed to evaluate the vector for xi-squared vector

class AuxParamBundle {
public:
    SDBList &sdbs;    
    size_t nCalls;
private:
    // We make sure there are no copy or default constructors to
    // preserve the integrity of our reference member.
    AuxParamBundle();
    AuxParamBundle(AuxParamBundle const&);
    AuxParamBundle const& operator=(AuxParamBundle const&);

    size_t refant;
    size_t nCorrelations;
    size_t corrStep;
    Double t0;
    Double reftime;
    std::map< Int, std::set< Int > > activeAntennas;
    std::map< Int, Int > antennaIndexMap;
    size_t activeCorr;
public:
    AuxParamBundle(SDBList& sdbs_, size_t refant, const std::map< Int, std::set<Int> >& activeAntennas_) :
        sdbs(sdbs_),
        nCalls(0),
        refant(refant),
        nCorrelations( sdbs.nCorrelations() > 1 ? 2 : 1 ),
        // nCorrelations( sdbs.nCorrelations() ),
        // nCorrelations( 2 ),
        corrStep( sdbs.nCorrelations() > 2 ? 3 : 1),
        activeAntennas( activeAntennas_ ),
        activeCorr( -1 )
        // corrStep( 3 )
        {
            Int last_index = sdbs.nSDB() - 1 ;
            t0 = sdbs(0).time()(0);
            Double tlast = sdbs(last_index).time()(0);
            reftime = 0.5*(t0 + tlast);
            // cerr << "AuxParamBundle reftime " << reftime << " t0 " << t0 <<" dt " << tlast - t0 << endl;
        }

    Double get_t0() {
        return t0;
    }
    Double
    get_ref_time() {
        return reftime;
    }
    size_t
    get_num_corrs() {
        //return sdbs.nCorrelations() > 1 ? 2 : 1;
        return nCorrelations;
    }
    size_t
    get_num_antennas() {
        if (activeCorr < 0) {
            throw(AipsError("Correlation out of range."));
        }
        std::set< Int > ants = activeAntennas.find(activeCorr)->second;
        return (size_t) ants.size();
    }
    size_t
    get_max_antenna_index() {
        return max( max(sdbs(0).antenna1()), max(sdbs(0).antenna2())) ;
    }
    // Sometimes there is Int, sometimes size_t; the following ones are casacore::Int.
    Int
    get_num_data_points() {
        Int nTotalRows = 0;
        for (Int i = 0; i != sdbs.nSDB(); i++) {
            nTotalRows += sdbs(i).nRows();
        }
        return 2 * nTotalRows * nCorrelations * sdbs.nChannels();
    }
    size_t 
    get_data_corr_index(size_t icorr) {
        if (icorr > nCorrelations) {
            throw(AipsError("Correlation out of range."));
        }
        size_t dcorr = icorr * corrStep;
        return dcorr;
    }
    bool
    isActive(size_t iant) {
        std::set<Int> ants = activeAntennas.find(activeCorr)->second;
        if (iant == refant) return true;
        else return (ants.find( iant ) != ants.end() );
    }
    Int
    get_param_index(size_t iant, size_t icor) {
        // here we use parallel correlation indices, because parameters
        // by definition only have one hand.
        if (iant == refant) return -1;
        int ipar = antennaIndexMap[iant];
        if (iant > refant) ipar -= 1;
        return 3*(ipar*nCorrelations + icor);
    }
    Int
    get_param_corr_index(size_t iant) {
        if (iant == refant) return -1;
        int ipar = antennaIndexMap[iant];
        if (iant > refant) ipar -= 1;
        return 3*ipar;        
    }
    size_t
    get_active_corr() {
        return activeCorr;
    }
    void
    set_active_corr(size_t icorr) {
        activeCorr = icorr;
        antennaIndexMap.clear();
        Int i = 0;
        std::set<Int>::iterator it;
        std::set<Int> ants = activeAntennas.find(activeCorr)->second;
        for (it = ants.begin(); it != ants.end(); it++) {
            antennaIndexMap[*it] = i++;
        }
    }
};
    
int
expb_f(const gsl_vector *param, void *d, gsl_vector *f)
{
    AuxParamBundle *bundle =  (AuxParamBundle *)d;
    SDBList& sdbs = bundle->sdbs;
    Double refTime = bundle->get_t0();

    gsl_vector_set_zero(f);
    Vector<Double> freqs = sdbs.freqs();
    
    size_t count = 0; // This is the master index.
    for (Int ibuf=0; ibuf < sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s ( sdbs(ibuf) );
        if ( !s.Ok() ) continue;

        Cube<Complex> v = s.visCubeCorrected();
        Cube<Bool> fl = s.flagCube();
        Cube<Float> weights = s.weightSpectrum();
           
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if ( s.flagRow()(irow) ) continue;

            Int ant1(s.antenna1()(irow));
            Int ant2(s.antenna2()(irow));
            if (!bundle->isActive(ant1) || !bundle->isActive(ant2))
                continue;            
            if (ant1==ant2) continue;

            // VisBuffer.h seems to suggest that a vb.visCube may have shape
            // (nCorr(), nChannel(), nRow())
            size_t icorr0 = bundle->get_active_corr();
            size_t dcorr = bundle->get_data_corr_index(icorr0);
            // We also need to get the right parameters for this,
            // polarization (icorr is an encoding of the
            // polarization of the correlation products).
            Int iparam1 = bundle->get_param_corr_index(ant1);
            Double phi0_1, tau1, r1;
            if (iparam1 >= 0) {
                phi0_1 = gsl_vector_get(param, iparam1+0);
                tau1 =   gsl_vector_get(param, iparam1+1);
                r1 =     gsl_vector_get(param, iparam1+2);
            } else {
                phi0_1 = 0.0;
                tau1 = 0.0;
                r1 = 0.0;
            }
            Int iparam2 = bundle->get_param_corr_index(ant2);
            Double phi0_2, tau2, r2;
            if (iparam2 >= 0) {
                phi0_2 = gsl_vector_get(param, iparam2+0);
                tau2 =   gsl_vector_get(param, iparam2+1);
                r2 =     gsl_vector_get(param, iparam2+2);
                // cerr << "phi0_2 " << phi0_2 << " tau2 " << tau2 << " r2 " << r2 << endl;
            } else {
                phi0_2 = 0.0;
                tau2 = 0.0;
                r2 = 0.0;
            }
            
            Float phi0 = phi0_2 - phi0_1;
            Float tau  = tau2 - tau1;
            Float r    = r2 - r1;
            for (size_t ichan = 0; ichan != v.ncolumn(); ichan++) {
                if ( fl(dcorr, ichan, irow) ) continue;
                Complex vis = v(dcorr, ichan, irow);
                Double w = weights(dcorr, ichan, irow);
                if (fabs(w) < FLT_EPSILON) continue;
                
                // We have to turn the delay back into seconds from nanoseconds.
                // Freq difference is in Hz, which comes out typically as 1e6 bands
                Double wDf = C::_2pi*(freqs(ichan) - freqs(0))*1e-9;
                //
                Double t1 = s.time()(0);
                Double ref_freq = freqs(0);
                Double wDt = C::_2pi*(t1 - refTime) * ref_freq; 

                Double mtheta = -(phi0 + tau*wDf + r*wDt); 
                Double vtheta = arg(vis);
                
                gsl_vector_set(f, count, w*(cos(mtheta) - cos(vtheta)));
                gsl_vector_set(f, count+1, w*(sin(mtheta)  - sin(vtheta)));
                
                count += 2;
            }
        }
    }
    return GSL_SUCCESS;
}


void
print_baselines(std::set<std::pair< Int, Int > > baselines) {
    cerr << "Baselines encountered ";
    std::set<std::pair< Int, Int > >::iterator it;
    for (it=baselines.begin(); it != baselines.end(); ++it) {
        cerr << "(" << it->first << ", " << it->second << ") ";
    }
    cerr << endl;
}

    
int
expb_df(CBLAS_TRANSPOSE_t TransJ, const gsl_vector* x, const gsl_vector *u, void *bundle_, gsl_vector *v, gsl_matrix *JTJ)
{

    // x is the current vector for which we're finding the jacobian.
    // if TransJ is true, evaluate J^T u and store in v.
    // Also store J^T . J in lower half of JTJ.
    std::set <std::pair < Int, Int> > baselines;
    AuxParamBundle *bundle =  (AuxParamBundle *)bundle_;

    SDBList& sdbs = bundle->sdbs;
    Vector<Double> freqs = sdbs.freqs();

    size_t count = 0; // This is the master index.

    gsl_vector_set_zero(v);
    gsl_matrix_set_zero(JTJ);
    
    Double refTime = bundle->get_t0();
    std::set< Int > params;
    for (Int ibuf=0; ibuf < sdbs.nSDB(); ibuf++) {
        // cerr << "OK so count = " << count << endl;
        SolveDataBuffer& s ( sdbs(ibuf) );
        if ( !s.Ok() ) continue;

        Cube<Complex> vis = s.visCubeCorrected();
        Cube<Bool> fl = s.flagCube();
        Cube<Float> weights = s.weightSpectrum();

        Double t1 = s.time()(0);
        // cerr << "ibuf " << ibuf << " t1 - t0 = " << t1 - t0 << endl;
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if ( s.flagRow()(irow) ) continue;

            Int ant1(s.antenna1()(irow));
            Int ant2(s.antenna2()(irow));

            if (ant1==ant2) continue;
            if (!bundle->isActive(ant1) || !bundle->isActive(ant2)) {
                // cerr << "Skipping " << ant1 << ", " << ant2 << endl;                   
                continue;
            }

            // VisBuffer.h seems to suggest that a vb.visCube may have shape
            // (nCorr(), nChannel(), nRow()) 

            size_t icorr0 = bundle->get_active_corr();
            size_t dcorr = bundle->get_data_corr_index(icorr0);
            // We also need to get the right parameters for this
            // polarization (icorr is an encoding of the
            // polarization of the correlation products).

            Int iparam1 = bundle->get_param_corr_index(ant1);
            Double phi0_1, tau1, r1;
            if (iparam1 >= 0) {
                phi0_1 = gsl_vector_get(x, iparam1+0);
                tau1 =   gsl_vector_get(x, iparam1+1);
                r1 =     gsl_vector_get(x, iparam1+2);
            } else {
                phi0_1 = 0.0;
                tau1 = 0.0;
                r1 = 0.0;
            }
            Int iparam2 = bundle->get_param_corr_index(ant2);
            Double phi0_2, tau2, r2;
            if (iparam2 >= 0) {
                phi0_2 = gsl_vector_get(x, iparam2+0);
                tau2 =   gsl_vector_get(x, iparam2+1);
                r2 =     gsl_vector_get(x, iparam2+2);
            } else {
                phi0_2 = 0.0;
                tau2 = 0.0;
                r2 = 0.0;
            }
            Double phi0 = phi0_2 - phi0_1;
            Double tau = tau2 - tau1;
            Double r = r2-r1;

            Double ref_freq = freqs(0); 
            Double wDt = C::_2pi*(t1 - refTime) * ref_freq; 
            // cerr << "Dt " << t1 - refTime << " ref_freq " << ref_freq << " wDt " << wDt << endl;
            bool found_data = false;
            
            for (size_t ichan = 0; ichan != vis.ncolumn(); ichan++) {
                if ( fl(dcorr, ichan, irow) ) continue;
                Double w = weights(dcorr, ichan, irow);
                if (fabs(w) < FLT_EPSILON) continue;
                found_data = true;
                // Add a 1e-9 factor because tau parameter is in nanoseconds.
                Double wDf = C::_2pi*(freqs(ichan) - freqs(0))*1e-9;
                //
                Double mtheta = -(phi0 + tau*wDf + r*wDt);
                Double ws = sin(mtheta);
                Double wc = cos(mtheta);

                if (iparam2 >= 0) {
                    // cerr << "insert param2! " << iparam2
                    // << " into array of size (" << JTJ->size1 << ", " <<  JTJ->size2 << ")"
                    //      << endl;
                    params.insert(iparam2);
                    /* 
                       What we want to express is just:
                       J[count + 0, iparam2 + 0] = w*-ws*-1.0; 
                       J[count + 1, iparam2 + 0] = w*+wc*-1.0;
                       J[count + 0, iparam2 + 1] = w*-ws*-wDf;
                       J[count + 1, iparam2 + 1] = w*+wc*-wDf;
                       J[count + 0, iparam2 + 2] = w*-ws*-wDt;
                       J[count + 1, iparam2 + 2] = w*+wc*-wDt;

                       But in the GSL multilarge framework we have to
                       be ready to calculate either J*u for a given u
                       or J^T*u, depending on the flag TransJ, and we also have to fill in the 
                       
                       v[iparam + ...] = J[count + ..., iparam + ...] * u[iparam + ...]

                       or
                       
                       v[iparam + ...] = J^T[iparam + ..., count + ...] * u[count + ...]

                       <https://www.gnu.org/software/gsl/doc/html/nls.html#c.gsl_multifit_nlinear_default_parameters>

                       "Additionally, the normal equations matrix J^T J should be stored in the lower half of JTJ."

                       So we should also use
                       JTJ[iparam + ..., iparam + ...] += J^T[iparam + ..., count + ...] J[count + ..., iparam + ...] 

                    */
                    if (TransJ==CblasNoTrans) {
                        // v = J u expressed here as v[count + 0] += J[count + 0, iparam2+0] *u[iparam2 + 0] etc.
                        // where we have to use gsl syntax and 
                        // Jacobians listed here 
                        // J[count + 0, iparam2 + 0]:
                        (*gsl_vector_ptr(v, count + 0)) += (w*-ws*-1.0) * gsl_vector_get(u, iparam2 + 0);
                        // J[count + 0, iparam2 + 1]:
                        (*gsl_vector_ptr(v, count + 0)) += (w*-ws*-wDf) * gsl_vector_get(u, iparam2 + 1);
                        // J[count + 0, iparam2 + 2]:
                        (*gsl_vector_ptr(v, count + 0)) += (w*-ws*-wDt) * gsl_vector_get(u, iparam2 + 2);
                        
                        // J[count + 1, iparam2 + 0]:
                        (*gsl_vector_ptr(v, count + 1)) += (w*+wc*-1.0) * gsl_vector_get(u, iparam2 + 0);
                        // J[count + 1, iparam2 + 1]:
                        (*gsl_vector_ptr(v, count + 1)) += (w*+wc*-wDf) * gsl_vector_get(u, iparam2 + 1);
                        // J[count + 1, iparam2 + 2]:
                        (*gsl_vector_ptr(v, count + 1)) += (w*+wc*-wDt) * gsl_vector_get(u, iparam2 + 2);
                    } else {
                        (*gsl_vector_ptr(v, iparam2 + 0)) += (w*-ws*-1.0) * gsl_vector_get(u, count + 0);
                        (*gsl_vector_ptr(v, iparam2 + 0)) += (w*+wc*-1.0) * gsl_vector_get(u, count + 1);
                        (*gsl_vector_ptr(v, iparam2 + 1)) += (w*-ws*-wDf) * gsl_vector_get(u, count + 0);
                        (*gsl_vector_ptr(v, iparam2 + 1)) += (w*+wc*-wDf) * gsl_vector_get(u, count + 1);
                        (*gsl_vector_ptr(v, iparam2 + 2)) += (w*-ws*-wDt) * gsl_vector_get(u, count + 0);
                        (*gsl_vector_ptr(v, iparam2 + 2)) += (w*+wc*-wDt) * gsl_vector_get(u, count + 1);
                    }
                    // JTJ part. I'm calculating these from the CblasNoTrans version.
                    // JTJ[i, k] = sum_j JT[i, j] * J[j, k] , which is to say
                    // JTJ[i, k] = sum_j J[j, i] * J[j, k].
                    // We're summing over the count + i vectors.
                    //
                    // Note that terms like (-ws*-1.0) * (-ws*-1.0) + (+wc*-1.0) * (+wc*-1.0)
                    // can be simplified to (ws*ws) + (wc*wc) and that that reduces to 1.
                    // But I'd like to have regression tests for some of that I guess.
                    if (JTJ) {
                        (*gsl_matrix_ptr(JTJ, iparam2 + 0, iparam2 + 0)) += w*w*-1.0*-1.0*((-ws) * (-ws) +  // J[count + 0, iparam2 + 0] * J[count + 0, iparam2 + 0]
                                                                                           (+wc) * (+wc));  // J[count + 1, iparam2 + 0] * J[count + 1, iparam2 + 0]
                        (*gsl_matrix_ptr(JTJ, iparam2 + 1, iparam2 + 0)) += w*w*-wDf*-1.0*((-ws) * (-ws) +  // J[count + 0, iparam2 + 1] * J[count + 0, iparam2 + 0]
                                                                                           (+wc) * (+wc));  // J[count + 1, iparam2 + 1] * J[count + 1, iparam2 + 0]
                        (*gsl_matrix_ptr(JTJ, iparam2 + 1, iparam2 + 1)) += w*w*-wDf*-wDf*((-ws) * (-ws) +  // J[count + 0, iparam2 + 1]^2
                                                                                           (+wc) * (+wc));  // J[count + 1, iparam2 + 1]^2 
                        (*gsl_matrix_ptr(JTJ, iparam2 + 2, iparam2 + 0)) += w*w*-wDt*-1.0*((-ws) * (-ws) +  // J[count + 0, iparam2 + 2] * J[count + 0, iparam2 + 0]
                                                                                           (+wc) * (+wc) ); // J[count + 1, iparam2 + 2] * J[count + 1, iparam2 + 0] 
                        (*gsl_matrix_ptr(JTJ, iparam2 + 2, iparam2 + 1)) += w*w*-wDt*-wDf*((-ws) * (-ws) +  // J[count + 0, iparam2 + 2] * J[count + 0, iparam2 + 1]
                                                                                           (+wc) * (+wc));  // J[count + 1, iparam2 + 2] * J[count + 1, iparam2 + 1]
                        (*gsl_matrix_ptr(JTJ, iparam2 + 2, iparam2 + 2)) += w*w*-wDt*-wDt*((-ws) * (-ws) +  // J[count + 0, iparam2 + 2]^2 
                                                                                           (+wc) * (+wc));  // J[count + 1, iparam2 + 2]^2
                    }

                }
                if (iparam1 >= 0) {
                    params.insert(iparam1);
                    if (TransJ==CblasNoTrans) {
                        (*gsl_vector_ptr(v, count + 0)) += gsl_vector_get(u, iparam1 + 0) * (w*-ws*+1.0);
                        (*gsl_vector_ptr(v, count + 0)) += gsl_vector_get(u, iparam1 + 1) * (w*-ws*+wDf); 
                        (*gsl_vector_ptr(v, count + 0)) += gsl_vector_get(u, iparam1 + 2) * (w*-ws*+wDt);
                        // 
                        (*gsl_vector_ptr(v, count + 1)) += gsl_vector_get(u, iparam1 + 0) * (w*+wc*+1.0);
                        (*gsl_vector_ptr(v, count + 1)) += gsl_vector_get(u, iparam1 + 1) * (w*+wc*+wDf);
                        (*gsl_vector_ptr(v, count + 1)) += gsl_vector_get(u, iparam1 + 2) * (w*+wc*+wDt);
                    } else {
                        // Transpose
                        (*gsl_vector_ptr(v, iparam1 + 0)) += gsl_vector_get(u, count + 0) * (w*-ws*+1.0);
                        (*gsl_vector_ptr(v, iparam1 + 0)) += gsl_vector_get(u, count + 1) * (w*+wc*+1.0);
                        (*gsl_vector_ptr(v, iparam1 + 1)) += gsl_vector_get(u, count + 0) * (w*-ws*+wDf);
                        (*gsl_vector_ptr(v, iparam1 + 1)) += gsl_vector_get(u, count + 1) * (w*+wc*+wDf);
                        (*gsl_vector_ptr(v, iparam1 + 2)) += gsl_vector_get(u, count + 0) * (w*-ws*+wDt);
                        (*gsl_vector_ptr(v, iparam1 + 2)) += gsl_vector_get(u, count + 1) * (w*+wc*+wDt);
                    } 
                    if (JTJ) {
                        (*gsl_matrix_ptr(JTJ, iparam1 + 0, iparam1 + 0)) += w*w*1.0*1.0*((-ws) * (-ws) +  // J[count + 0, iparam1 + 0]^2
                                                                                         (+wc) * (+wc));  // J[count + 1, iparam1 + 0]^2
                        (*gsl_matrix_ptr(JTJ, iparam1 + 1, iparam1 + 0)) += w*w*wDf*1.0*((-ws) * (-ws) +  // J[count + 0, iparam1 + 1] * J[count + 0, iparam1 + 0]
                                                                                         (+wc) * (+wc));  // J[count + 1, iparam1 + 1] * J[count + 1, iparam1 + 0]
                        (*gsl_matrix_ptr(JTJ, iparam1 + 1, iparam1 + 1)) += w*w*wDf*wDf*((-ws) * (-ws) +  // J[count + 0, iparam1 + 1]^2
                                                                                         (+wc) * (+wc));  // J[count + 1, iparam1 + 1]^2 
                        (*gsl_matrix_ptr(JTJ, iparam1 + 2, iparam1 + 0)) += w*w*wDt*1.0*((-ws) * (-ws) +  // J[count + 0, iparam1 + 2] * J[count + 0, iparam1 + 0]
                                                                                         (+wc) * (+wc) ); // J[count + 1, iparam1 + 2] * J[count + 1, iparam1 + 0] 
                        (*gsl_matrix_ptr(JTJ, iparam1 + 2, iparam1 + 1)) += w*w*wDt*wDf*((-ws) * (-ws) +  // J[count + 0, iparam1 + 2] * J[count + 0, iparam1 + 1]
                                                                                         (+wc) * (+wc));  // J[count + 1, iparam1 + 2] * J[count + 1, iparam1 + 1]
                        (*gsl_matrix_ptr(JTJ, iparam1 + 2, iparam1 + 2)) += w*w*wDt*wDt*((-ws) * (-ws) +  // J[count + 0, iparam1 + 2]^2 
                                                                                         (+wc) * (+wc));  // J[count + 1, iparam1 + 2]^2
                    }
                }
                count += 2;
            } // loop over rows
            if (found_data) {
                std::pair<Int, Int> antpair = std::make_pair(ant1, ant2);
                bool newBaseline = (baselines.find( antpair ) == baselines.end());
                if ( newBaseline ) {
                    // print_baselines(baselines);
                    // cerr << "Adding (" << ant1 << ", " << ant2 << ")" << endl;
                    baselines.insert( antpair );
                }
                // only print weights to ref ant.
                if (0 && newBaseline && ((iparam1 == -1) || (iparam2 == -1))) {
                    cerr << "baseline (" << ant1 << ", " << ant2 << ") "
                         << "weight " << weights(dcorr, vis.ncolumn()/2, irow) << endl;
                }
            }
        }
    }
    if (DEVDEBUG && 0) {
        cerr << "Param indices ";
        std::copy(
            params.begin(),
            params.end(),
            std::ostream_iterator<Int>(std::cerr, " ")
            );
        cerr << endl;
        print_baselines(baselines);
        cerr << "count " << count << endl;
        cerr <<"JTJ " << std::scientific << endl;
        for (int i=0; i!=JTJ->size1; i++) {
            for (int j=0; j!=JTJ->size2; j++) {
                cerr << gsl_matrix_get(JTJ, i, j) << " ";
            }
            cerr << endl;
        }
        cerr << endl;
    }
    return GSL_SUCCESS;
}
    

void
least_squares_driver(SDBList& sdbs, Matrix<Float>& param, Int refant, const std::map< Int, std::set<Int> >& activeAntennas, LogIO& logSink) {
    // n below is number of variables,
    // p is number of parameters

    // Int midTimeIndex = sdbs.nSDB() / 2;
    // Double refTime = sdbs(midTimeIndex).time()(0);

    AuxParamBundle bundle( sdbs, refant, activeAntennas );
    for (size_t icor=0; icor != bundle.get_num_corrs(); icor++ ) {
        bundle.set_active_corr(icor);
        if (bundle.get_num_antennas() == 0) {
            logSink << "No antennas for correlation " << icor << "; not running least-squares solver." << LogIO::POST;
            continue;
        }
        if (bundle.get_num_antennas() == 1) {
            logSink << "No baselines for correlation " << icor << "; not running least-squares solver." << LogIO::POST;
            continue;
        }
        // Three parameters for every antenna.
        size_t p = 3 * (bundle.get_num_antennas() - 1);
        // We need to store complex visibilities in a real matrix so we
        // just store real and imaginary components separately.
        size_t n = 2 * bundle.get_num_data_points();

        if (DEVDEBUG) {
            cerr << "p " << p << " n " << n << endl;
        }
        // Parameters for the least-squares solver.
        const double param_tol = 1.0e-50;
        //const double gtol = 1e-50; // 1e-30; 
        const double gtol = 1.0e-50; // 1e-30; 
        const double ftol = 1.0e-50;   // eps rel
        const size_t max_iter = 200;

        const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
        gsl_multilarge_nlinear_parameters params = gsl_multilarge_nlinear_default_parameters();
        // params.trs = gsl_multilarge_nlinear_trs_lm;
        params.scale = gsl_multilarge_nlinear_scale_more;
        params.solver = gsl_multilarge_nlinear_solver_cholesky;

        gsl_multilarge_nlinear_workspace *w = gsl_multilarge_nlinear_alloc(T, &params, n, p);
        gsl_multilarge_nlinear_fdf f;

        f.f = &expb_f;
        /* Can't set to NULL for finite-difference Jacobian in multilarge case. */
        f.df =  &expb_df;   
        f.n = n;    /* number of data points */
        f.p = p;    /* number of parameters */
        f.params = &bundle;

        // Our original param is a matrix of (3*nCorr, nElem).
        // We have to transcribe it to a vector.

        gsl_vector *gp = gsl_vector_alloc(p);
        gsl_vector_set_zero(gp);

        for (size_t iant=0; iant != bundle.get_max_antenna_index()+1; iant++) {
            if ( !bundle.isActive(iant) ) {
                logSink << "Skipping antenna " << iant << " for correlation " << icor << "." << LogIO::POST;
                continue;
            }
            Int ind = bundle.get_param_corr_index(iant);
            if (ind < 0) continue;
            gsl_vector_set( gp, ind+0, param(3*icor + 0, iant) );
            gsl_vector_set( gp, ind+1, param(3*icor + 1, iant) );
            gsl_vector_set( gp, ind+2, param(3*icor + 2, iant) );
        }
        gsl_vector *gp_orig = gsl_vector_alloc( p );
        // Keep a copy of original parameters
        gsl_vector_memcpy (gp_orig, gp);
        // initialise workspace
        gsl_multilarge_nlinear_init(gp, &f, w);
    
        // compute initial residual norm */
        gsl_vector *res_f = gsl_multilarge_nlinear_residual(w);
        double chi0 = gsl_blas_dnrm2(res_f);

        int info;
        int status = gsl_multilarge_nlinear_driver(max_iter, param_tol, gtol, ftol,
                                                   NULL, NULL, &info, w);
    
        double chi1 = gsl_blas_dnrm2(res_f);

        gsl_vector_sub( gp_orig, w->x);
        gsl_vector *diff = gp_orig;
        double diffsize = gsl_blas_dnrm2(diff);
    
        gsl_vector *res = gsl_multilarge_nlinear_position(w);
        
        for (size_t iant=0; iant != bundle.get_max_antenna_index()+1; iant++) {
            if ( !bundle.isActive(iant) ) continue;
            Int iparam = bundle.get_param_corr_index(iant);
            if (iparam<0) continue;
            if (1) {
                bool flag = false;
                if ( fabs(gsl_vector_get(diff, iparam + 0) > FLT_EPSILON) ) {
                    flag = true;
                }
                if ( fabs(gsl_vector_get(diff, iparam + 1) > FLT_EPSILON) ) {
                    flag = true;
                }
                if ( fabs(gsl_vector_get(diff, iparam + 2) > 1e-30) ) {
                    flag = true;
                }
                logSink << "Old values for ant " << iant << " correlation " << icor 
                // if (flag) cerr << "Changed! "<< endl;
                        << ": Angle " << param(3*icor + 0, iant)
                        << " delay " << param(3*icor + 1, iant) << " ns "
                        << " rate " << param(3*icor + 2, iant) << "."
                        << endl
                        << "New values for ant " << iant << " correlation " << icor 
                        << ": Angle " << gsl_vector_get(res, iparam+0)
                        << " delay " << gsl_vector_get(res, iparam+1) << " ns "
                        << " rate " << gsl_vector_get(res, iparam+2) << "."
                        << LogIO::POST;
            }
            param(3*icor + 0, iant) = gsl_vector_get( res, iparam+0 );
            param(3*icor + 1, iant) = gsl_vector_get( res, iparam+1 );
            param(3*icor + 2, iant) = gsl_vector_get( res, iparam+2 );
        }

        logSink <<  "Least squares complete for correlation " << icor << "," << endl
                << "number of iterations: " <<  gsl_multilarge_nlinear_niter(w) << endl
            // << "reason for stopping: " << ( (info == 1) ? "small step size" : "small gradient" ) << endl
                << "initial |f(x)| = " << chi0 << endl
                << "final   |f(x)| = " << chi1 << endl
                << "final step taken = " << diffsize 
                << "." << LogIO::POST;
        gsl_vector_free( gp );
        gsl_multilarge_nlinear_free( w );
    }    
}

    


// **********************************************************
//  CTRateAwareTimeInterp1 Implementations
//

CTRateAwareTimeInterp1::CTRateAwareTimeInterp1(NewCalTable& ct,
					       const casacore::String& timetype,
					       casacore::Array<Float>& result,
					       casacore::Array<Bool>& rflag) :
  CTTimeInterp1(ct,timetype,result,rflag)
{}

// Destructor (nothing to do locally)
CTRateAwareTimeInterp1::~CTRateAwareTimeInterp1() {}

Bool CTRateAwareTimeInterp1::interpolate(Double newtime) {
  
  // Call generic first
  if (CTTimeInterp1::interpolate(newtime)) {
    // Only if generic yields new calibration
    // NB: lastWasExact_=exact in generic
    applyPhaseRate(timeType().contains("nearest") || lastWasExact_);
    return true;
  }
  else
    // No change
    return false;

}

// Do the phase rate math
void CTRateAwareTimeInterp1::applyPhaseRate(Bool single)
{
  Int ispw=mcols_p->spwId()(0);
  MSSpectralWindow msSpw(ct_.spectralWindow());
  ROMSSpWindowColumns msCol(msSpw);
  Vector<Double> refFreqs;
  msCol.refFrequency().getColumn(refFreqs,True);

  //  cout << "time = " << (currTime_ - timeRef_) << endl;

  if (single) {
    for (Int ipol=0;ipol<2;ipol++) {
      Double dtime=(currTime_-timeRef_)-timelist_(currIdx_);
      Double phase=result_(IPosition(2,ipol*3,0));
      Double rate=result_(IPosition(2,ipol*3+2,0));
      phase+=2.0*C::pi*rate*refFreqs(ispw)*dtime;
      result_(IPosition(2,ipol*3,0))=phase;
    }
  } else {
    Vector<uInt> rows(2); indgen(rows); rows+=uInt(currIdx_);
    Cube<Float> r(mcols_p->fparamArray("",rows));

    Vector<Double> dtime(2);
    dtime(0)=(currTime_-timeRef_)-timelist_(currIdx_);
    dtime(1)=(currTime_-timeRef_)-timelist_(currIdx_+1);
    Double wt=dtime(1) / (dtime(1)-dtime(0));


    for (Int ipol=0;ipol<2;ipol++) {
      Vector<Double> phase(2), rate(2);
      phase(0)=r.xyPlane(0)(IPosition(2,ipol*3,0));
      phase(1)=r.xyPlane(1)(IPosition(2,ipol*3,0));
      rate(0)=r.xyPlane(0)(IPosition(2,ipol*3+2,0));
      rate(1)=r.xyPlane(1)(IPosition(2,ipol*3+2,0));

      phase(0)+=2.0*C::pi*rate(0)*refFreqs(ispw)*dtime(0);
      phase(1)+=2.0*C::pi*rate(1)*refFreqs(ispw)*dtime(1);

      Vector<Complex> ph(2);
      ph(0)=Complex(cos(phase(0)),sin(phase(0)));
      ph(1)=Complex(cos(phase(1)),sin(phase(1)));
      ph(0)=Float(wt)*ph(0) + Float(1.0-wt)*ph(1);
      result_(IPosition(2,ipol*3,0))=arg(ph(0));
    }
  }
}




// **********************************************************
//  FringeJones Implementations
//
FringeJones::FringeJones(VisSet& vs) :
    VisCal(vs),             // virtual base
    VisMueller(vs),         // virtual base
    GJones(vs)       // immediate parent
{
    if (prtlev()>2) cout << "FringeJones::FringeJones(vs)" << endl;
}

FringeJones::FringeJones(String msname,Int MSnAnt,Int MSnSpw) :
    VisCal(msname,MSnAnt,MSnSpw),             // virtual base
    VisMueller(msname,MSnAnt,MSnSpw),         // virtual base
    GJones(msname,MSnAnt,MSnSpw)    // immediate parent
{
    if (prtlev()>2) cout << "FringeJones::FringeJones(msname,MSnAnt,MSnSpw)" << endl;
}

FringeJones::FringeJones(const MSMetaInfoForCal& msmc) :
    VisCal(msmc),             // virtual base
    VisMueller(msmc),         // virtual base
    GJones(msmc)    // immediate parent
{
    if (prtlev()>2) cout << "FringeJones::FringeJones(msmc)" << endl;
}

FringeJones::FringeJones(Int nAnt) :
    VisCal(nAnt), 
    VisMueller(nAnt),
    GJones(nAnt)
{
    if (prtlev()>2) cout << "FringeJones::FringeJones(nAnt)" << endl;
}

FringeJones::~FringeJones() {
    if (prtlev()>2) cout << "FringeJones::~FringeJones()" << endl;
}

void FringeJones::setApply(const Record& apply) {
    // Call parent to do conventional things
    GJones::setApply(apply);

    if (calWt()) 
        logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;

    // Enforce calWt() = false for delays
    calWt()=false;

    // Extract per-spw ref Freq for phase(delay) calculation
    //  from the CalTable
    // TBD:  revise as per refFreq decisions
    MSSpectralWindow msSpw(ct_->spectralWindow());
    ROMSSpWindowColumns msCol(msSpw);
    msCol.refFrequency().getColumn(KrefFreqs_,true);
    KrefFreqs_/=1.0e9;  // in GHz

    /// Re-assign KrefFreq_ according spwmap (if any)
    if (spwMap().nelements()>0) {
        Vector<Double> tmpfreqs;
        tmpfreqs.assign(KrefFreqs_);
        for (uInt ispw=0;ispw<spwMap().nelements();++ispw)
            if (spwMap()(ispw)>-1)
                KrefFreqs_(ispw)=tmpfreqs(spwMap()(ispw));
    }
}

void FringeJones::setCallib(const Record& callib,
                            const MeasurementSet& selms) {

    // Call parent to do conventional things
    SolvableVisCal::setCallib(callib,selms);

    /*
    if (calWt()) 
        logSink() << " (" << this->typeName() << ": Enforcing calWt()=false for phase/delay-like terms)" << LogIO::POST;
    */
    // Enforce calWt() = false for delays
    calWt()=false;

    // Extract per-spw ref Freq for phase(delay) calculation
    //  from the CalTable 
   KrefFreqs_.assign(cpp_->refFreqIn());
    KrefFreqs_/=1.0e9;  // in GHz

    // Re-assign KrefFreq_ according spwmap (if any)
    if (spwMap().nelements()>0) {
        Vector<Double> tmpfreqs;
        tmpfreqs.assign(KrefFreqs_);
        for (uInt ispw=0;ispw<spwMap().nelements();++ispw)
            if (spwMap()(ispw)>-1)
                KrefFreqs_(ispw)=tmpfreqs(spwMap()(ispw));
    }
}

void FringeJones::setSolve(const Record& solve) {

    // Call parent to do conventional things
    GJones::setSolve(solve);
    // Trap unspecified refant:
    if (refant()<0)
        throw(AipsError("Please specify a good reference antenna (refant) explicitly."));
    if (solve.isDefined("zerorates")) {
        zeroRates() = solve.asBool("zerorates");
    }
}

void FringeJones::calcAllJones() {

  if (prtlev()>6) cout << "       FringeJones::calcAllJones()" << endl;

  // Should handle OK flags in this method, and only
  //  do Jones calc if OK

  Vector<Complex> oneJones;
  Vector<Bool> oneJOK;
  Vector<Float> onePar;
  Vector<Bool> onePOK;

  ArrayIterator<Complex> Jiter(currJElem(),1);
  ArrayIterator<Bool>    JOKiter(currJElemOK(),1);
  ArrayIterator<Float>   Piter(currRPar(),1);
  ArrayIterator<Bool>    POKiter(currParOK(),1);

  Double phase;
  for (Int iant=0; iant<nAnt(); iant++) {

    for (Int ich=0; ich<nChanMat(); ich++) {
      
      oneJones.reference(Jiter.array());
      oneJOK.reference(JOKiter.array());
      onePar.reference(Piter.array());
      onePOK.reference(POKiter.array());

      for (Int ipar=0;ipar<nPar();ipar+=3) {
        if (onePOK(ipar)) {
          phase=onePar(ipar);
          phase+=2.0*C::pi*onePar(ipar+1)*
            (currFreq()(ich)-KrefFreqs_(currSpw()));
          phase+=2.0*C::pi*onePar(ipar+2)*KrefFreqs_(currSpw())*1e9*
            (currTime() - refTime());
          oneJones(ipar/3)=Complex(cos(phase),sin(phase));
          oneJOK(ipar/3)=True;
        }
      }
      
      // Advance iterators
      Jiter.next();
      JOKiter.next();
    }
    // Step to next antenns's pars
    Piter.next();
    POKiter.next();
  }
}


void
FringeJones::calculateSNR(Int nCorr, DelayRateFFT drf, Float ref_freq, Float df0, Float dt0 ) {
    Matrix<Float> sRP(solveRPar().nonDegenerate(1));
    Matrix<Bool> sPok(solveParOK().nonDegenerate(1));
    Matrix<Float> sSNR(solveParSNR().nonDegenerate(1));
    
    for (size_t icor=0; icor != nCorr; icor++ ) {
        const set<Int>& activeAntennas = drf.getActiveAntennasCorrelation(icor);
        for (Int iant=0; iant != nAnt(); iant++) {
            if (iant == refant()) {
                Double maxsnr = 1.0e6;
                sSNR(3*icor + 0, iant) = maxsnr;
                sSNR(3*icor + 1, iant) = maxsnr;
                sSNR(3*icor + 2, iant) = maxsnr;
            }
            else if ( activeAntennas.find( iant ) != activeAntennas.end() ) {
                Double delay = sRP(3*icor + 1, iant);
                Double rate = sRP(3*icor + 2, iant);
                // Note that DelayRateFFT::snr is also used to calculate SNRs for the least square values!
                Float snrval = drf.snr(icor, iant, delay, rate);
                sSNR(3*icor + 0, iant) = snrval;
                sSNR(3*icor + 1, iant) = snrval;
                sSNR(3*icor + 2, iant) = snrval;
            } else {
                sPok(3*icor + 0, iant) = false;
                sPok(3*icor + 1, iant) = false;
                sPok(3*icor + 2, iant) = false;
            }
        }
    }
}



// void FringeJones::solveLotsOfSDBs(SDBList& sdbs)

void
FringeJones::selfSolveOne(SDBList& sdbs) {
    solveRPar()=0.0;
    solveParOK()=false; 
    solveParErr()=1.0; // Does nothing?
    // Maybe we put refFreq, refTime stuff in here?
    Vector<Double> myRefFreqs;
    // FIXME: Update for multiple SWs!
    // FIMXE: No, really!
    MSSpectralWindow msSpw(ct_->spectralWindow());
    ROMSSpWindowColumns spwCol(msSpw);
    spwCol.refFrequency().getColumn(myRefFreqs, true);
    Double ref_freq = myRefFreqs(currSpw());
    Double ref_time = refTime();
    Double dt0 = (ref_time - sdbs(0).time()(0));
    Double df0 = ref_freq - sdbs.freqs()(0);

    // ROMSSpWindowColumns are a casacore type!
    // http://casacore.github.io/casacore/classcasacore_1_1ROMSSpWindowColumns.html


    /* We'd like to learn about the spectral windows in the actual
     * dataset but I haven't succeeded yet. */
    // cerr << "nSpw() " << nSpw() << " nelements() " << spwMap().nelements() << endl;
    // cerr << "spwMap() " << endl;
    // cerr << "combine() " << combine() << endl;
    // cerr << "combspw() " << combspw() << endl;
    // throw(AipsError("That's quite enough of that."));
    
    if (0) {
        /* The MSSpectralWindow above is constructed from the
         * calibration table (ct_). In the case of the fringe jones
         * table this reports only one frequency channel for each
         * spectral window. This may be a bug in the definition of the
         * spectral window in the calibration table; until it is
         * addressed I need to work around it.
         */

        MeasurementSet ms( msName() );
        ROMSColumns mscol(ms);
        const ROMSSpWindowColumns& spwcol(mscol.spectralWindow());


        cerr << "combine() " << combine() << endl;
        cerr << "combspw() " << combspw() << endl;
        
        // spwMap seems to report 1, even when you are combining spectral windows;
        // nSPW() reports 30 even when you aren't.

        // Can we get the spectral windows without grovelling over the entire set of SolveDataBuffers?
        for (uInt ispw=0; ispw != spwMap().nelements(); ++ispw) {
            uInt jspw = spwMap()(ispw);
            const Vector<Double>& chanfreqs = spwcol.chanFreq()(ispw);
            // const ScalarColumn<int>& nchan = spwCol.numChan();
            int nchan = spwcol.numChan()(ispw);
        
            cerr << "ispw " << ispw << " freq(0) " << chanfreqs(0) << endl;
            cerr << "shape " << chanfreqs.shape() << endl;
            cerr << "nelements " << chanfreqs.nelements() << endl;
            cerr << "nchan " << nchan << endl;
            // ncc.spectralWindow().chan
            // spwins.insert(spw);
        }
        throw(AipsError("That's quite enough of that."));
    }
    
    // Pausing here:
    // throw(AipsError("Just checking ref values."));
    // the values seemed reasonable
    DelayRateFFT drf( sdbs, refant());
    drf.FFT(); 
    logSink() << "FFT completed." << LogIO::POST;
    drf.searchPeak();
    logSink() << "Searched for peaks in FFT." << LogIO::POST;
    Matrix<Float> sRP(solveRPar().nonDegenerate(1));
    Matrix<Bool> sPok( solveParOK().nonDegenerate(1) );
    Matrix<Float> sSNR(solveParSNR().nonDegenerate(1));

    
    // Map from MS antenna number to index 
    Int ncol = drf.param().ncolumn();
    for (Int i=0; i!=ncol; i++) {
        IPosition start(2, 0,                  i);
        IPosition stop( 2, drf.param().nrow(), 1);
        IPosition step( 2, 1,                  1);
        Slicer sl(start, stop, step, Slicer::endIsLength);
        sRP(sl) = drf.param()(sl);
        sPok(sl) = !(drf.flag()(sl));
    }


    size_t nCorrOrig( sdbs(0).nCorrelations() );
    size_t nCorr = (nCorrOrig> 1 ? 2 : 1); // number of p-hands

    // cerr << "(Reminder:) Current spectral window =" << currSpw() << endl;
    calculateSNR(nCorr, drf, ref_freq, df0, dt0);

    set<Int> belowThreshold;
    // In the Python interface there is a parameter minsnr.
    Float threshold = minSNR();
    
    for (size_t icor=0; icor != nCorr; icor++ ) {
        const set<Int>& activeAntennas = drf.getActiveAntennasCorrelation(icor);
        for (Int iant=0; iant != nAnt(); iant++) {
            if (iant != refant() && ( activeAntennas.find( iant ) != activeAntennas.end() )) {
                Float s = sSNR(3*icor + 0, iant);
		// Start the log message; finished below
		logSink() << "Antenna " << iant << " correlation has (FFT) SNR of " << s;
                if (s < threshold) {
                    belowThreshold.insert( iant );
                    logSink() << " below threshold (" << threshold << ")";
                    // Don't assume these will be flagged later; do it right away.
                    // (The least squares routine will eventually become optional.)
                    sPok(3*icor + 0, iant) = false;
                    sPok(3*icor + 1, iant) = false;
                    sPok(3*icor + 2, iant) = false;
                }
		// Finish the log message
		logSink() << "." << LogIO::POST;
            }
        }
        // We currently remove the antennas below SNR threshold from
        // the object used to handle the FFT fringe search.
        drf.removeAntennasCorrelation(icor, belowThreshold);
        if (DEVDEBUG) {
            drf.printActive();
        }
    }
    if (1) {
        logSink() << "Starting least squares optimization." << LogIO::POST;
        // Note that least_squares_driver is *not* a method of
        // FringeJones so we pass everything in, including the logSink
        // reference.  Note also that sRP is passed by reference and
        // altered in place.
        least_squares_driver(sdbs, sRP, refant(), drf.getActiveAntennas(), logSink());
    }


    // We can zero the rates here (if needed) whether we did least squares or not.
    
    if ( zeroRates() ) {
        logSink() << "Zeroing delay rates in calibration table." << LogIO::POST;
        
        for (size_t icor=0; icor != nCorr; icor++ ) {
            for (Int iant=0; iant != nAnt(); iant++) {
                sRP(3*icor + 2, iant) = 0.0;
            }
        }
    }

    if (DEVDEBUG) {
        cerr << "Just sayin': nAnt()*nCorr=" << nAnt() * nCorr << " while sRP has shape " << sRP.shape() << endl;
        cerr << "Full shape of solveRPar(): " << solveRPar().shape() << endl;
        cerr << "Shape of solveAllRPar(): " << solveAllRPar().shape() << endl;
    }
    
    if (DEVDEBUG) {
        cerr << "Ref time " << MVTime(refTime()/C::day).string(MVTime::YMD,7) << endl;
        cerr << "df0 " << df0 << " dt0 " << dt0 << " ref_freq*dt0 " << ref_freq*dt0 << LogIO::POST;
        cerr << "ref_freq " << ref_freq << endl;
        cerr << "df0 " << df0 << " dt0 " << dt0 << " ref_freq*dt0 " << ref_freq*dt0 << endl;
    }
    // calculateSNR is a method of FringeJones so it has access to the sSNR (et al) arrays.
    calculateSNR(nCorr, drf, ref_freq, df0, dt0);
}


void
FringeJones::solveOneVB(const VisBuffer&) {
    throw(AipsError("VisBuffer interface not supported!"));
}


} //# NAMESPACE CASA - END


/*

fringefit(vis="n14c2.ms", caltable="fail.fj", field="",spw="1",intent="",
          selectdata=True, timerange="", antenna="", scan="5", observation="",
          msselect="", solint="inf", refant="EF", minsnr=3.0, append=False,
          gaintable=['n14c2.gcal'], parang=False)
*/
