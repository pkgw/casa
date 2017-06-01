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
#include <gsl/gsl_multifit_nlin.h>


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

// **************************************************************
// DelayRateFFT modeled on DelayFFT(const VisBuffer&, Double padBW, Int refant) in KJones.{cc|h}

// Remark: Double, Int et al are in casacore namespace. Need to be
// qualified in headers, I guess, but here we are not in headers and we
// are using that namespace.
class DelayRateFFT {
private:
    Int refant_;
    Int nPadFactor_;
    Int nt_;
    Int nPadT_;
    Int nChan_;
    Int nPadChan_;
    Int nElem_;
    Double f0_, df_;
    Double t0_, t1_, dt_;
    Double padBW_;
    Array<Complex> Vpad_; 
    Int nCorr_;
   // 
    Matrix<Float> param_;
    Matrix<Bool> flag_; //?

public:
    DelayRateFFT(SDBList& sdbs, Int refant) :
        refant_( refant ),
        nPadFactor_ ( 8 ), // , // FIXME: 8 is surely better, just debugging
        nt_( sdbs.nSDB() ),
        nPadT_( nPadFactor_ * nt_ ),
        nChan_ ( sdbs.nChannels() ),
        nPadChan_ ( nPadFactor_*nChan_ ),
        f0_(sdbs.freqs()(0)/1.e9),      // GHz
        df_(sdbs.freqs()(1)/1.e9-f0_),
        Vpad_()
        {
        Int veryDebug ( 1 );
        Int nCorrOrig( sdbs(0).nCorrelations() );
        nCorr_ = (nCorrOrig> 1 ? 2 : 1); // number of p-hands
        // when we get the visCubecorrected it is already
        // reduced to parallel hands, but there isn't a
        // corresponding method for flags.
        Int sC = (nCorrOrig > 2 ? 3 : 1); // step for p-hands

        SolveDataBuffer& s0 ( sdbs(0) );
        nElem_ =  1 + max( max(s0.antenna1()), max(s0.antenna2())) ;
        // Can't get timeInterval, fails with error
        // "Caught exception: Exception: Can't fill VisBuffer component TimeInterval:
        // Not attached to VisibilityIterator."
        // dt_ = s0.timeInterval()(0);
        cerr << "There are " << sdbs.nSDB() << " SDBs." << endl;
        if (sdbs.nSDB() < 2) {
            const Vector<Double>& time = sdbs(0).time();
            int nt = time.shape()[0];
            for (int i = 0; i < min(20, nt); i++) {
                cerr << "DTime(" << i << ") = " << time(i) -time(0) << endl;
            }
            cerr << "Antenna1 list: " << sdbs(0).antenna1() << endl;
            cerr << "Antenna2 list: " << sdbs(0).antenna2() << endl;
            cerr << "Times for first SolveDataBuffer: " << time - time(0) << endl;
            cerr << "# entries " << nt << endl;
            cerr << "Time step: " << time(1) - time(0) << endl;
            throw(AipsError("Not enough sdbs!"));
        }
        dt_ = sdbs(1).time()(0) - s0.time()(0);
        cerr << "Time step=" << dt_ << "s" << endl;
            
        
        IPosition paddedDataSize(4, nCorr_, nElem_, nPadT_, nPadChan_);
        Vpad_.resize(paddedDataSize);
                
        for (Int ibuf=0; ibuf!=nt_; ibuf++) {
            SolveDataBuffer& s ( sdbs(ibuf) );
            if ( !s.Ok() )
                continue;

            Int nr = 0;
            for (Int irow=0; irow!=s.nRows(); irow++) {
                if ( s.flagRow()(irow) )
                    continue;
                Int iant;
                Int sign;
                Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
                if (a1 == a2) {
                    continue;
                }
                else if (a1 == refant_) {
                    iant = a2;
                    sign = -1;
                }
                else if (a2 == refant_) {
                    iant = a1;
                    sign = 1;
                }
                else {
                    continue;
                }
                // v has shape (nelems, ?, nrows, nchannels)
                Cube<Complex> v = s.visCubeCorrected();
                Cube<Bool> fl = s.flagCube();


                IPosition start( 4,      0,      iant, ibuf, 0);
                IPosition stop(  4,      nCorr_,    1,    1, nChan_);
                IPosition stride(4,      1,         1,    1, 1);
                Slicer sl1(start,     stop, stride, Slicer::endIsLength);


                // cerr << "nCorrOrig " << nCorrOrig << endl;
                Slicer sl2(IPosition(3, 0,         0, irow),
                           IPosition(3, nCorr_,    nChan_, 1),
                           IPosition(3, sC,        1,  1), Slicer::endIsLength);
                
                Slicer flagSlice(IPosition(3, 0,         0, irow),
                                 IPosition(3, nCorr_, nChan_, 1),
                                 IPosition(3, sC,        1, 1), Slicer::endIsLength);
                nr++;
                if (veryDebug) {
                    cerr << "nr " << nr
                         << " irow " << endl
                         << "Vpad shape " << Vpad_.shape() << endl
                         << "v shape " << v.shape() << endl
                         << "sl2 " << sl2 << endl
                         << "sl1 " << sl1 << endl
                         << "flagSlice " << flagSlice << endl;
                }
                Array<Complex> rhs = v(sl2).nonDegenerate();
                unitize(rhs);
                Vpad_(sl1).nonDegenerate() = sign ? rhs : -rhs;
                // Zero flagged entries.
                // FIXME: Stolen from the VisBuf version
                // and not yet ported

                Array<Bool> flagged( fl(flagSlice).nonDegenerate() );
                if (veryDebug) {
                    cerr << "flagSlice " << flagSlice << endl
                         << "fl.shape() " << fl.shape() << endl
                         << "Vpad_.shape() " << Vpad_.shape() << endl
                         << "flagged.shape() " << flagged.shape() << endl
                         << "sl1 " << sl1 << endl;
                    cerr << "halfflagged" << endl;
                }
                Vpad_(sl1).nonDegenerate()(flagged) = Complex(0.0);
            }
        }
        cerr << "Constructed a DelayRateFFT object." << endl;
    }

        
    DelayRateFFT(const VisBuffer& vb, Int nPadFactor, Int refant) :
        refant_(refant),
        nPadFactor_(nPadFactor),

        nPadT_(0),    // set in body
        nPadChan_(nPadFactor_*vb.nChannel()),
        nElem_(), // antenna-based
        f0_(vb.frequency()(0)/1.e9),      // GHz
        df_(vb.frequency()(1)/1.e9-f0_),
        Vpad_(),
        nCorr_(0),    // set in body
        param_(), //? really
        flag_() {
        
        // VisBuffer facts
        Int nCorrOrig = vb.nCorr();
        Int nChan = vb.nChannel();
        Int nRow = vb.nRow();

        // Discern effective shapes
        nCorr_=(nCorrOrig> 1 ? 2 : 1); // number of p-hands
        Int sC=(nCorrOrig > 2 ? 3 : 1); // step for p-hands
        // double interval_length = e3.getTime("s").getValue();

        t0_ = min(vb.time());
        t1_ = max(vb.time());
        dt_ = vb.timeInterval()(0); // We will check this is a constant as we go.
        nt_ = round((t1_ - t0_) / dt_ + 0.5); // Round up to be sure.

        MEpoch e1;
        MVEpoch e2, e3;
        vb.timeRange(e1, e2, e3);

        cerr << "Times " << t0_ << ", " << t1_ << ", " << t1_ - t0_ << ", " << dt_ << ", => " << nt_ << endl;
        cerr << "Epochery " << e1 << ", " << e2 <<  ", " << e3 << endl;
            
        nPadT_ = nPadFactor_*nt_;
        cerr << "*****************************************************************************"
             << endl << "resizing..." << endl
             << "*****************************************************************************"
             << endl;

        // Shape the data

        // Note: the first argument to IPosition constructor is the number of
        // dimensions; optional further arguments are their sizes in each
        // dimension.
        // IPosition ip1(4, nCorr_, nt_, nPadChan_, nElem_);
        // FIXME: Shuffling dimensions
        
        IPosition ip1(4, nCorr_, nElem_, nPadT_, nPadChan_);
        Vpad_.resize(ip1);
        Vpad_.set(0.0);

        cerr << "*****************************************************************************"
             << endl << "resized..." << endl
             << "*****************************************************************************"
             << endl;
        // We grovel over the rows
        for (Int irow=0; irow != nRow; irow++) {
            Int sgn, iant;
            Int a1(vb.antenna1()(irow));
            Int a2(vb.antenna2()(irow));
            if (!vb.flagRow()(irow) && a1!=a2) {
                if (a1==refant) {
                    iant = a2;
                    sgn = -1;
                }
                else if (a2==refant) {
                    iant = a2;
                    sgn = +1;
                }
                else
                    continue; // Ignore baselines not to refant
                Float t( vb.time()(irow) );
                // Note: we do NOT multiply itime up by nPadFactor_! We
                // fill the bottom left corner of the original Vpad_
                // matrix and leave it to the FFT to smear the response
                // over the full width of the matrix.
                Int itime( round(0.5 + (t - t0_)/dt_) ); 
                // Still following the KJones DelayFFT constructor.
                // source slice:
                Slicer sl0( Slice(0,nCorr_,sC), Slice(), Slice(irow,1,1) );
                // target slice:
                // FIXME: Can't have four(4) slice indices!
                //Slicer sl1( Slice(), Slice(iant,1,1), Slice(itime,1,1), Slice(0,nChan,1) ); 
                IPosition start(4, 0,      iant, itime, 0);
                IPosition stop( 4, nCorr_, 1, 1, nChan);
                IPosition stride(4, sC, 1, 1, 1);
                Slicer sl1(start, stop, stride, Slicer::endIsLength);

                Cube<Complex> vC(vb.visCube()(sl0));
                unitize(vC);

                // Zero flagged channels.
                // (<small@jive.eu> Note to self: Array::nonDegenerate removes degenerate axes;
                // if passed an IPosition this denotes axes to ignore.

                cerr << "*****************************************************************************"
                     << endl
                     << "vC shape " << vC.shape() << endl
                     << "Vpad_ shape " << Vpad_.shape() << endl
                     << "*****************************************************************************"
                     << endl;

                
                Slicer fsl0(Slice(),Slice(irow,1,1));                
                Array<Bool> fl(vb.flag()(fsl0).nonDegenerate(IPosition(1,0)));
                for (Int icor=0; icor != nCorr_; ++icor) {
                    // Used to be IPosition(1,1) which was the chan axis; that's moved to 2 now.
                    Slicer sl3 = Slicer(Slice(icor,1,1), Slice(), Slice());
                    vC(sl3).nonDegenerate(IPosition(1, 1))(fl) = Complex(0.0);
                }
                cerr << "*****************************************************************************"
                     << endl << "zeroed flags..." << endl
                     << "*****************************************************************************"
                     << endl;

                // Normalise the sign. Can't seem to multiply by sgn; can do the below.
                // Vpad_(sl1) = Float(sgn)*vC;
                Vpad_(sl1) = (sgn == 1) ? vC : conj(vC);
                cerr << "*****************************************************************************"
                     << endl << "assigned slice." << endl
                     << "*****************************************************************************"
                     << endl;
                
            }
        }
    } // DelayRateFFT
    
    // The following are copied from KJones.h definition of DelayFFT.
    // I'm putting them here because I haven't yet split out the header version.

    const Array<Bool>& flag() const { return flag_; }
    const Array<Complex>& Vpad() const { return Vpad_; }
    const Matrix<Float>& param() const { return param_; }
    Int refant() const { return refant_; }
    
    void FFT() {
        // Axes are 0: correlation (i.e., hand of polarization), 1: antenna, 2: time, 3: channel
        Vector<Bool> ax(4, false);
        ax(2) = true;
        ax(3) = true;
        // Also copied from DelayFFT in KJones.
        ArrayLattice<Complex> c(Vpad_);
        // But variable c is not returned and is not a member?
        // IT SEEMS that the FFT is in place, and that FFTing c actually changes Vpad_
        LatticeFFT::cfft0(c, ax, true);
        //LatticeFFT::cfft(c,ax,true);
        // 
    }

    // Are we allowed to use std::pair in Casa?
    std::pair<Bool, Float>  xinterp(Float alo, Float amax, Float ahi) {
        Float denom (alo-2.0*amax+ahi);
        Bool cond = amax>0.0 && abs(denom)>0.0 ;
        Float fpk = cond ? 0.5-(ahi-amax)/denom : 0.0;
        return std::make_pair(cond, fpk);
    }
    
    void searchPeak() {
        int veryDebug( 1 );
        // Recall param_ -> [phase, delay, rate] for each correlation
        param_.resize(3*nCorr_, nElem_); // Srsly, why we do this here?
        param_.set(0.0);
        // Note: It looks like we have to have one flag per parameter?
        flag_.resize(3*nCorr_, nElem_);
        flag_.set(true);  // all flagged initially 
        cerr << "nt_ " << nt_ << " nPadChan_ " << nPadChan_ << endl;
        for (Int icorr=0;icorr<nCorr_;++icorr) {
            flag_(icorr*3 + 0, refant())=false; 
            flag_(icorr*3 + 1, refant())=false;
            flag_(icorr*3 + 2, refant())=false;
            for (Int ielem=0;ielem<nElem_;++ielem) {
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
                cerr << "aS.shape()=" <<aS.shape() << endl;

                Matrix<Float> amp( amplitude(aS) );
                Int ipkch(0);
                Int ipkt(0);
                Float amax(-1.0);
                // Unlike KJones we have to iterate in time too
                for (Int itime=0; itime != nt_; itime++) {
                    for (Int ich=0; ich != nPadChan_; ++ich) {
                        if (amp(itime, ich) > amax) {
                            ipkch=ich;
                            ipkt = itime;
                            amax=amp(itime, ich);
                        }
                    }
                }
                // We used to print out a slice once, but that's not needed now.
                // cerr << "Slice: "
                //      << amplitude(
                //          Vpad_(Slicer(
                //                    IPosition( 4, icorr, ielem,      ipkt,     0),
                //                    IPosition( 4,     1,     1,      1, nPadChan_),
                //                    IPosition( 4,     1,     1,      1,        1),
                //                    Slicer::endIsLength)).nonDegenerate())
                //      << endl;

                // Finished grovelling. Now we have the location of the
                // maximum amplitude.
                Float alo_ch = amp(ipkt, (ipkch > 0) ? ipkch-1 : nPadChan_-1);
                Float ahi_ch = amp(ipkt, ipkch<(nPadChan_-1) ? ipkch+1 : 0);
                std::pair<Bool, Float> maybeFpkch = xinterp(alo_ch, amax, ahi_ch);
                // We handle wrapping while looking for neighbours
                Float alo_t = amp(ipkt > 0 ? ipkt-1 : nPadT_ -1,     ipkch);
                Float ahi_t = amp(ipkt < (nPadT_ -1) ? ipkt+1 : 0,   ipkch);
                std::pair<Bool, Float> maybeFpkt = xinterp(alo_t, amax, ahi_t);

                // Int sgn = (icorr==0) ? -1 : +1;
                Int sgn = -1;
                if (maybeFpkch.first and maybeFpkt.first) {
                    // Phase
                    Complex c = aS(ipkt, ipkch);
                    param_(icorr*3 + 0, ielem) = sgn*arg(c);
                    // Delay
                    Float delay = (ipkch + maybeFpkch.second)/Float(nPadChan_);
                    if (delay > 0.5) delay -= 1.0;           // fold
                    delay /= (df_);                            // nsec
                    param_(icorr*3 + 1, ielem) = sgn*delay; // FIXME: Looks like -delay is better?
                    // Rate
                    Float rate = (ipkt + maybeFpkt.second)/Float(nPadT_);
                    if (rate > 0.5) rate -= 1.0;
                    rate /= dt_;
                    rate /= (1e9 * f0_);
                    param_(icorr*3 + 2, ielem) = sgn*rate; // FIXME: so probably also -rate?
                    cerr << "maybeFpkch.second=" << maybeFpkch.second << ", "
                         << "df_ " << df_ << ", nt_ " << nt_
                         << " fpkch " << (ipkch + maybeFpkch.second) << endl;
                    cerr << "Found peak for element " << ielem << " correlation " << icorr
                         << " ipkt=" << ipkt << "/" << nPadT_ << ", ipkch=" << ipkch << "/" << nPadChan_
                         << "; delay " << delay << ", rate " << rate
                         << ", phase " << arg(c) << endl;
                    // Set 3 flags.
                    flag_(icorr*3 + 0, ielem)=false; 
                    flag_(icorr*3 + 1, ielem)=false;
                    flag_(icorr*3 + 2, ielem)=false;
                   
                }
                else {
                    cerr << "No peak for element " << ielem << " correlation " << icorr << endl;
                }
            }
        }
    }
    // The idiom used in KJones solvers is:
    // DelayFFT delfft1(vbga(ibuf), ptbw, refant());
    // delfft1.FFT();
    // delfft1.shift(f0[0]);
    // delfft1.searchPeak();
}; // End of class DelayRateFFT.



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
public:
    AuxParamBundle(SDBList& sdbs_, size_t refant) :
        sdbs(sdbs_),
        nCalls(0),
        refant(refant),
        nCorrelations( sdbs.nCorrelations() > 1 ? 2 : 1 ),
        // nCorrelations( sdbs.nCorrelations() ),
        // nCorrelations( 2 ),
        corrStep( sdbs.nCorrelations() > 2 ? 3 : 1)
        // corrStep( 3 )
        {
            Int last_index = sdbs.nSDB() - 1 ;
            t0 = sdbs(0).time()(0);
            Double tlast = sdbs(last_index).time()(0);
            reftime = 0.5*(t0 + tlast);
            cerr << "AuxParamBundle reftime " << reftime << " dt " << tlast - t0 << endl;
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
        SolveDataBuffer& s0 = sdbs(0);
        return 1 + max( max(s0.antenna1()), max(s0.antenna2()) );
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
    int
    get_param_index(size_t iant, size_t icor) {
        // here we use parallel correlation indices, because parameters
        // by definition only have one hand.
        if (iant == refant) return -1;
        if (iant > refant) iant -= 1;
        return 3*(iant*nCorrelations + icor);
    }
};
    
int
expb_f(const gsl_vector *param, void *d, gsl_vector *f)
{
    
    AuxParamBundle *bundle =  (AuxParamBundle *)d;
    SDBList& sdbs = bundle->sdbs;
    // Double refTime = bundle->get_ref_time();
    Double refTime = bundle->get_t0();

    // if ((bundle->nCalls % 10) == 0) {
    if (0) {
        cerr << "Entering expb_f (call " << bundle->nCalls++ << ")." << endl;
    }

    for (size_t i=0; i != f->size; i++) {
        gsl_vector_set(f, i, 0.0);
    }
    Vector<Double> freqs = sdbs.freqs();
    cerr << "freqs(0) = " << freqs(0) << endl;
    
    size_t count = 0; // This is the master index.
    for (Int ibuf=0; ibuf < sdbs.nSDB(); ibuf++) {
        // cerr << "OK so count = " << count << endl;
        SolveDataBuffer& s ( sdbs(ibuf) );
        if ( !s.Ok() ) continue;

        Cube<Complex> v = s.visCubeCorrected();
        Cube<Bool> fl = s.flagCube();
        Cube<Float> weights = s.weightSpectrum();
           
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if ( s.flagRow()(irow) ) continue;

            Int ant1(s.antenna1()(irow));
            Int ant2(s.antenna2()(irow));
            if (ant1==ant2) continue;

            // VisBuffer.h seems to suggest that a vb.visCube may have shape
            // (nCorr(), nChannel(), nRow()) 
            for (size_t icorr0 = 0; icorr0 != bundle->get_num_corrs(); icorr0++) {

                size_t dcorr = bundle->get_data_corr_index(icorr0);
                // We also need to get the right parameters for this,
                // polarization (icorr is an encoding of the
                // polarization of the correlation products).

                Int iparam1 = bundle->get_param_index(ant1, icorr0);
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
                Int iparam2 = bundle->get_param_index(ant2, icorr0);
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
                Float tau = (tau2 - tau1);
                Float r = r2-r1;
                for (size_t ichan = 0; ichan != v.ncolumn(); ichan++) {
                    // if (1) {
                    if (count >= f->size-4) {
                        cerr << "ibuf= " << ibuf << "/" << sdbs.nSDB()
                             << " irow= " << irow << "/" << s.nRows()
                             << " (corr) dcorr= " << dcorr << "/" << v.nrow()
                             << " (column) ichan= " << ichan << "/" << v.ncolumn()
                             << " count " << count 
                             << "." << endl;
                        // throw(AipsError("Index too big."));
                    }
                    if ( fl(dcorr, ichan, irow) ) continue;
                    Complex vis = v(dcorr, ichan, irow);
                    Float w = weights(dcorr, ichan, irow);
                    // We have to turn the delay back into seconds from nanoseconds.
                    // Freq difference is in Hz, which comes out typically as 1e6 bands
                    Double wDf = 2.0*C::pi*(freqs(ichan) - freqs(0))*1e-9;
                    //
                    Double t1 = s.time()(irow);
                    Double ref_freq = freqs(0);
                    
                    Double wDt = 2.0*C::pi*(t1 - refTime) * ref_freq; 
                    //
                    //Float mtheta = -(phi0 + tau*wDf); 
                    Float mtheta = -(phi0 + tau*wDf + r*wDt); 
                    Float vtheta = arg(vis);
                    // FIXME! Weights!
                    // gsl_vector_set(f, count, w*(cos(mtheta) - cos(vtheta)));
                    // gsl_vector_set(f, count+1, w*(sin(mtheta)  - sin(vtheta)));
                    gsl_vector_set(f, count, w*(cos(mtheta) - cos(vtheta)));
                    gsl_vector_set(f, count+1, w*(sin(mtheta)  - sin(vtheta)));

                    count += 2;
                }
            }
        }
    }
    if (1) {
        cerr << "Count " << count << endl;
        cerr << "Leaving expb_f." << endl;
    }
    return GSL_SUCCESS;
}


int
expb_df(const gsl_vector *param, void *d, gsl_matrix *J)
{
    
    cerr << "Entering expb_df." << endl;
    AuxParamBundle *bundle =  (AuxParamBundle *)d;

    SDBList& sdbs = bundle->sdbs;
    Vector<Double> freqs = sdbs.freqs();

    size_t count = 0; // This is the master index.

    for (size_t i=0; i!=J->size1; i++) {
        for (size_t j=0; j!=J->size2; j++) {
            gsl_matrix_set(J, i, j, 0);
        }
    }

    //Double refTime = bundle->get_ref_time();
    Double refTime = bundle->get_t0();
    
    std::set< Int > params;
    std::set <std::pair < Int, Int> > baselines;
    for (Int ibuf=0; ibuf < sdbs.nSDB(); ibuf++) {
        // cerr << "OK so count = " << count << endl;
        SolveDataBuffer& s ( sdbs(ibuf) );
        if ( !s.Ok() ) continue;

        Cube<Complex> v = s.visCubeCorrected();
        Cube<Bool> fl = s.flagCube();
        Cube<Float> weights = s.weightSpectrum();


        Double t1 = s.time()(ibuf);
        // cerr << "ibuf " << ibuf << " t1 - t0 = " << t1 - t0 << endl;
            
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if ( s.flagRow()(irow) ) continue;

            Int ant1(s.antenna1()(irow));
            Int ant2(s.antenna2()(irow));
            if (ant1==ant2) continue;

            std::pair<Int, Int> antpair = std::make_pair(ant1, ant2);
            bool newBaseline = (baselines.find( antpair ) == baselines.end());
            
            baselines.insert( antpair );

            // VisBuffer.h seems to suggest that a vb.visCube may have shape
            // (nCorr(), nChannel(), nRow()) 

            for (size_t icorr0 = 0; icorr0 != bundle->get_num_corrs(); icorr0++) {
                // 
                size_t dcorr = bundle->get_data_corr_index(icorr0);
                // We also need to get the right parameters for this
                // polarization (icorr is an encoding of the
                // polarization of the correlation products).

                Int iparam1 = bundle->get_param_index(ant1, icorr0);
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
                
                Int iparam2 = bundle->get_param_index(ant2, icorr0);
                Double phi0_2, tau2, r2;
                if (iparam2 >= 0) {
                    phi0_2 = gsl_vector_get(param, iparam2+0);
                    tau2 =   gsl_vector_get(param, iparam2+1);
                    r2 =     gsl_vector_get(param, iparam2+2);
                } else {
                    phi0_2 = 0.0;
                    tau2 = 0.0;
                    r2 = 0.0;
                }
                Double phi0 = phi0_2 - phi0_1;
                Double tau = tau2 - tau1;
                Double r = r2-r1;


                
                Double ref_freq = freqs(0); 
                Double wDt = 2.0*C::pi*(t1 - refTime) * ref_freq; 

                for (size_t ichan = 0; ichan != v.ncolumn(); ichan++) {
                    if ( fl(dcorr, ichan, irow) ) continue;
                    Double w = 1.0; // VERY FIXME: This is the weights.
                    // Double w = weights(dcorr, ichan, irow);
                    // Add a 1e-9 factor because tau parameter is in milliseconds.
                    Double wDf = 2.0*C::pi*(freqs(ichan) - freqs(0))*1e-9;
                    //
                    Double mtheta = -(phi0 + tau*wDf + r*wDt);
                    Double ws = sin(mtheta);
                    Double wc = cos(mtheta);


                    if (0 && newBaseline) {
                        Double eps = sqrt(FLT_EPSILON)*abs(r);
                        Double mtheta2 = -(phi0 + tau*wDf + (r+eps)*wDt);
                        Double mtheta0 = -(phi0 + tau*wDf + (r-eps)*wDt);
                        Double rderivest = (cos(mtheta2) - cos(mtheta0))/(2*eps);
                        Double rderivofficial = -ws*-wDt;
                        cerr << "refTime " << refTime << " dt " << t1 - refTime << " wDt " << wDt << " eps " << eps << endl
                             << "rderiv est " << rderivest << " calc " << rderivofficial << endl;
                    }
                    
                    if (iparam2 >= 0) {
                        params.insert(iparam2);
                        gsl_matrix_set (J, count + 0, iparam2 + 0, w*-ws*-1.0);
                        gsl_matrix_set (J, count + 1, iparam2 + 0, w*+wc*-1.0);
                        gsl_matrix_set (J, count + 0, iparam2 + 1, w*-ws*-wDf);
                        gsl_matrix_set (J, count + 1, iparam2 + 1, w*+wc*-wDf);
                        gsl_matrix_set (J, count + 0, iparam2 + 2, w*-ws*-wDt);
                        gsl_matrix_set (J, count + 1, iparam2 + 2, w*+wc*-wDt);
                    }
                    if (iparam1 >= 0) {
                        params.insert(iparam1);
                        gsl_matrix_set (J, count + 0, iparam1 + 0, w*-ws*+1.0);
                        gsl_matrix_set (J, count + 1, iparam1 + 0, w*+wc*+1.0);
                        gsl_matrix_set (J, count + 0, iparam1 + 1, w*-ws*+wDf);
                        gsl_matrix_set (J, count + 1, iparam1 + 1, w*+wc*+wDf);
                        gsl_matrix_set (J, count + 0, iparam1 + 2, w*-ws*+wDt);
                        gsl_matrix_set (J, count + 1, iparam1 + 2, w*+wc*+wDt);
                    }
                    count += 2;
                }
            }
        }
    }
    cerr << "Leaving expb_df." << endl;
    cerr << "Count " << count << endl;
    cerr << "Param indices ";
    std::copy(
        params.begin(),
        params.end(),
        std::ostream_iterator<Int>(std::cerr, " ")
        );
    cerr << endl;
    
    cerr << "Jacobian slice" << endl;
    for (size_t i=10; i!=20; i++) {
        for (size_t j=0; j!=J->size2; j++) {
            Float f = gsl_matrix_get(J, i, j);
            if (fabs(f) > FLT_EPSILON) {
                cerr << "J(" <<i << ", " << j <<") = " << f << endl;
            }
        }
    }

    for (size_t i=0; i!=J->size1; i++) {
        size_t j=10;
        Float f = gsl_matrix_get(J, i, j);
        if (fabs(f) > FLT_EPSILON) {
            cerr << "J(" <<i << ", " << j <<") = " << f << endl;
        }
    }

    
    cerr << "Baselines encountered ";
    std::set<std::pair< Int, Int > >::iterator it;
    for (it=baselines.begin(); it != baselines.end(); ++it) {
        cerr << "(" << it->first << ", " << it->second << ") ";
    }
    cerr << endl;
    return GSL_SUCCESS;
}

int
expb_fdf(const gsl_vector *param, void *data, gsl_vector *f, gsl_matrix *J)
{
     expb_f(param, data, f);
     expb_df(param, data, J);

     return GSL_SUCCESS;
}

void
least_squares_driver(SDBList& sdbs, Matrix<Float>& param, Int refant) {
    cerr << "Enter the least_squares_driver" << endl;
    // n below is number of variables,
    // p is number of parameters


    // Int midTimeIndex = sdbs.nSDB() / 2;
    // Double refTime = sdbs(midTimeIndex).time()(0);

    AuxParamBundle bundle( sdbs, refant );
    cerr << "Constructed a bundle" << endl;
    cerr << "num corrs: " << bundle.get_num_corrs() << " num antennas: " << bundle.get_num_antennas() << endl;

    cerr << "Reftime again " << bundle.get_ref_time() << endl;
    // throw(AipsError("Yeah."));

    size_t p = 3 * bundle.get_num_corrs() * (bundle.get_num_antennas() - 1);
    size_t n = 2 * bundle.get_num_data_points();

    cerr << "n: " << n << " p: " << p << endl;
    
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, n, p);

    gsl_multifit_function_fdf f;
    f.f = &expb_f;
    f.df = &expb_df;   /* set to NULL for finite-difference Jacobian */
    //f.fdf = &expb_fdf;
    // f.df = NULL;
    
    f.n = n;    /* number of data points */
    f.p = p;    /* number of parameters */
    f.params = &bundle;

    gsl_vector *gp = gsl_vector_alloc( p );
    cerr << "Allocated vector gp of size " << p << "." <<endl;
    for (size_t icor=0; icor != bundle.get_num_corrs(); icor++ ) {
        for (size_t iant=0; iant != bundle.get_num_antennas(); iant++) {
            Int ind = bundle.get_param_index(iant, icor);
            if (ind < 0) continue;
            // cerr << "icor " << icor << " iant " << iant << " ind " << ind << "." << endl;
            gsl_vector_set( gp, ind+0, param(3*icor + 0, iant) );
            gsl_vector_set( gp, ind+1, param(3*icor + 1, iant) );
            gsl_vector_set( gp, ind+2, param(3*icor + 2, iant) );
        }
    }
    // It is said that param is a matrix of (3*nCorr, nElem)
    cerr << "Initialized parameter vector." << endl;
    gsl_multifit_fdfsolver_set(s, &f, gp);
    
    // compute initial residual norm */
    gsl_vector *res_f = gsl_multifit_fdfsolver_residual(s);
    double chi0 = gsl_blas_dnrm2(res_f);

    // int status = gsl_multifit_fdfsolver_iterate (s);

    // Calculate the jacobian manually
    // gsl_vector *y2 = gsl_vector_alloc( n );
    // gsl_matrix *J2 = gsl_matrix_alloc( n, p);
    // f.df(s->x, f.params, J2);


    // I keep trying to do this. It keeps not working.
    // gsl_vector *f3 = gsl_vector_alloc( n );
    // gsl_matrix *J3 = gsl_matrix_alloc( n, p);
    // // Calculate the jacobian and the function
    // gsl_multifit_fdfsolver_dif_fdf(gp, &f, f3, J3);
    
    // solve the system with a maximum of max_iter iterations */
    int info;
    const double param_tol = 1e-30;
    const double gtol = 1e-30; 
    const double ftol = 0.0;   // eps rel
    const size_t max_iter = 5;

    int status = gsl_multifit_fdfsolver_driver(s, max_iter, param_tol, gtol, ftol, &info);
    
    gsl_vector *g = s->g;
    for (size_t iparam=0; iparam!=p; iparam++) {
        cerr << "g(" << iparam << ") = " << gsl_vector_get(g, iparam) << endl;
    }

    // gsl_multifit_covar(J, 0.0, covar);
    // compute final residual norm
    // double chi1 = gsl_blas_dnrm2(f3);
    double chi1 = gsl_blas_dnrm2(res_f);

    gsl_vector *diff = s->dx;
    double diffsize = gsl_blas_dnrm2(diff);
    
    // Was I really getting the output?
    // A day of poking results in this:
    // For a test case, only one parameter (delay for antenna 11) changes,
    // and that change is not applied to actual parameter vector.
    gsl_vector *res = gsl_multifit_fdfsolver_position(s);
    for (size_t icor=0; icor != bundle.get_num_corrs(); icor++) {
        for (size_t iant=0; iant!=bundle.get_num_antennas(); iant++) {
            Int iparam = bundle.get_param_index(iant, icor);
            if (iparam<0) continue;
            cerr.precision(20);
            if (1) {
                bool flag = false;
                if ( fabs(gsl_vector_get(diff, iparam + 0) > FLT_EPSILON) ) {
                    flag = true;
                    Float psi0 = param(3*icor + 0, iant);
                    Float psi1 = gsl_vector_get(res, iparam+0);

                    cerr << "Psi changed by " << gsl_vector_get(diff, iparam + 0)
                         << " for " << "Antenna "  << iant << " correlation " << icor << " " <<  endl;
                    cerr << "psi0 = "    << psi0 << " psi1 = " << psi1 << " diff= " << psi1 - psi0  << endl;
                }
                if ( fabs(gsl_vector_get(diff, iparam + 1) > FLT_EPSILON) ) {
                    flag = true;
                    Float tau0 = param(3*icor + 1, iant);
                    Float tau1 = gsl_vector_get(res, iparam+1);
                    cerr << "Delay changed by " << gsl_vector_get(diff, iparam + 1)
                         << " for " << "Antenna "  << iant << " correlation " << icor << " " <<  endl;
                    cerr << "delay0 = " << tau0 << " delay1 = " << tau1 << " diff= " << tau1 - tau0 << " " << endl;
                    cerr << "diff2 = " << gsl_vector_get(diff, iparam + 1)  << endl;
                }
                if ( fabs(gsl_vector_get(diff, iparam + 2) > FLT_EPSILON) ) {
                    flag = true;
                    Float r0 = param(3*icor + 2, iant);
                    Float r1 = gsl_vector_get(res, iparam+2);
                    cerr << "Rate changed by " << gsl_vector_get(diff, iparam + 2)
                         << " for " << "Antenna "  << iant << " correlation " << icor << " " <<  endl;
                    cerr << "rate0 = " << r0 << " rate1 = " << r1 << " diff = " << r1 - r0 << " " << endl;
                }
                if (flag) cerr << endl;
            }
            param(3*icor + 0, iant) = gsl_vector_get( res, iparam+0);
            param(3*icor + 1, iant) = gsl_vector_get( res, iparam+1);
            param(3*icor + 2, iant) = gsl_vector_get( res, iparam+2);
        }
    }

    cerr << "Summary from method '" << gsl_multifit_fdfsolver_name(s) << "'" << endl
         << "number of iterations: " <<  gsl_multifit_fdfsolver_niter(s) << endl
         << "function evaluations: " << f.nevalf << endl
         << "Jacobian evaluations: " << f.nevaldf << endl
         << "reason for stopping: " << ( (info == 1) ? "small step size" : "small gradient" ) << endl
         << "initial |f(x)| = " << chi0 << endl
         << "final   |f(x)| = " << chi1 << endl
         << "final step taken = " << diffsize 
         << "." << endl;
    //<< "status = " << gsl_strerror(status) << endl;
    
    // // gsl_matrix_free(covar);
    // // gsl_matrix_free(J);
    // // gsl_rng_free(r);
    gsl_multifit_fdfsolver_free( s );
    gsl_vector_free( gp );
    
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

FringeJones::FringeJones(const Int& nAnt) :
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

    void FringeJones::selfSolveOne(SDBList& sdbs) {
        solveLotsOfSDBs(sdbs);
        
        // Implement actual solve here!!!
        //  E.g., gather data from visCubeCorrected in the SolveDataBuffers 
        //   in the SDBList and feed to solving code
        //  Typically, each SolveDataBuffer will contain a single timestamp 
        //    in a single spw
        //  E.g., see KJones::solveOneSDBmbd(SDBList&)
}

void FringeJones::solveLotsOfSDBs(SDBList& sdbs) {
    solveRPar()=0.0;
    solveParOK()=false; 
    solveParErr()=1.0; // Does nothing?
    cerr << "In solveLotsOfSDBs with " 
         << sdbs.nSDB() << " SolveDataBuffers " << endl;

    // Maybe we put refFreq, refTime stuff in here?
    Vector<Double> myRefFreqs;
    MSSpectralWindow msSpw(ct_->spectralWindow());
    ROMSSpWindowColumns msCol(msSpw);
    msCol.refFrequency().getColumn(myRefFreqs,true);
    Double ref_freq = myRefFreqs(currSpw());
    Double ref_time = refTime();
    Double dt0 = (ref_time - sdbs(0).time()(0));
    Double df0 = ref_freq - sdbs.freqs()(0);

    cerr << "ref_freq " << ref_freq << ", ref_time " << ref_time << endl;
    cerr << "dt0 " << dt0 << " df0 " << df0 << endl;
    // Pausing here:
    // throw(AipsError("Just checking ref values."));
    // the values seemed reasonable
    DelayRateFFT drf( sdbs, refant());
    drf.FFT(); 
    cerr << "Actually FFTed!" << endl;
    drf.searchPeak();
    cerr << "Actually searchPeaked!" << endl;
    Matrix<Float> sRP(solveRPar().nonDegenerate(1));
    Matrix<Bool> sPok(solveParOK().nonDegenerate(1));
    
    cerr << "Apparently"
         << " sRP has size " << sRP.size() << endl 
         << " parameters have shape " << drf.param().ncolumn() << endl
         << " flag matrix has shape " << sPok.shape() 
         << " whereas I have " << drf.flag().shape() << endl;


    Int ncol = drf.param().ncolumn();
    for (Int i=0; i!=ncol; i++) {
        IPosition start(2, 0,                  i);
        IPosition stop( 2, drf.param().nrow(), 1);
        IPosition step( 2, 1,                  1);
        Slicer sl(start, stop, step, Slicer::endIsLength);
        sRP(sl) = drf.param()(sl);
        sPok(sl) = !(drf.flag()(sl));
    }


    //sRP = drf.param();
    // cerr << "Flags " << drf.flag();
    //sPok = (!drf.flag());

    if (1) {
        cerr << "(Finally) having a go at least squares." << endl;
        least_squares_driver(sdbs, sRP, refant());
    }
    
    size_t nCorrOrig( sdbs(0).nCorrelations() );
    size_t nCorr = (nCorrOrig> 1 ? 2 : 1); // number of p-hands

    cerr << "df0 " << df0 << " dt0 " << dt0 << " ref_freq*dt0 " << ref_freq*dt0 << endl;
    for (size_t iant=0; iant != nAnt(); iant++) {
        for (size_t icor=0; icor != nCorr; icor++ ) {
            Double phi0 = sRP(3*icor + 0, iant);
            Double delay = sRP(3*icor + 1, iant);
            Double rate = sRP(3*icor + 2, iant);
            Double delta1 = df0*delay;
            Double delta2 = ref_freq*dt0*rate;
            // FIXME: or should we be subtracting?
            sRP(3*icor + 0, iant) -= 2*C::pi*(delta1+delta2);
            cerr << "For " << iant << " correlation " << icor << "." << endl;
            cerr << "phi0 " << phi0 << " delay " << delay << " rate " << rate << endl;
            // cerr << "Adding " << 360*delta1 << " and " << 360*delta2 << " degrees." << endl << endl;
            cerr << "Subtracting " << 360*delta1 << " and " << 360*delta2 << " degrees." << endl << endl;
        }
    }
    
       
    
}

void FringeJones::solveOneVB(const VisBuffer& vb) {
    // This seems to be the place to do things.
    solveRPar()=0.0;
    solveParOK()=false;
    
    // throw(AipsError("FringeJones: FIXME by implementing
    // solveOneVB")); My theory is that it's easiest to start from
    // solveOneVBmbd in KJones.h|cc which uses the auxilliary class
    // DelayFFT to simplify the local logic.  The alternative is
    // solveOneVB which looks like a version unfactored

    Int nch( vb.nChannel() );
    Vector<Double> chf( vb.frequency() ); 

    Double f0( chf(0)/1.0e9 );           // GHz
    Double df( (chf(1)-chf(0))/1.0e9 );  // GHz
    Double flo( f0 ), fhi( f0+nch*df );
    Double tbw( fhi-flo );
    cerr << "tbw = " << tbw << "  (" << flo << "-" << fhi << ")" << endl;
    Int padFactor( 8 );
    // Double ptbw( tbw*padFactor );  // pad total bw by 8X

    casacore::MEpoch e1;
    casacore::MVEpoch e2, e3;
    vb.timeRange(e1, e2, e3);
    cerr << "Epochs" << e1 << ", " << e2 << ", " << e3 << endl;

    cerr << vb.flagRow() << endl;
    LogicalArray mask(!vb.flagRow());

    MaskedArray<Double> maskTime(vb.time(), mask);

    Double minTime = min(maskTime);
    Double maxTime = max(maskTime);
    cerr << "Min time: "<< minTime << " max time: " << maxTime << endl;
    cerr << "Time matrix shape: " << vb.time().shape() << endl;
    
    double interval_length = e3.getTime("s").getValue();
    if (interval_length <= FLT_EPSILON) {
        cerr << "VisBuffer starting at " << e1 << " has length " << interval_length << " s; skipping."<< endl;
        return;
    }

    
    cerr << "*****************************************************************************"
         << endl << "solving..." << endl
         << "*****************************************************************************"
         << endl;

    DelayRateFFT drfft1(vb, padFactor, refant());

    cerr << "*****************************************************************************"
         << endl << "constructed..." << endl
         << "*****************************************************************************"
         << endl;


    drfft1.FFT();
    // Don't bother shifting.
    drfft1.searchPeak(); 
    // cout << delfft1.delay()(Slice(0,1,1),Slice()) << endl;
    

    // sRP is a matrix of corr, 0, nantenna size
    Matrix<Float> sRP( solveRPar().nonDegenerate(1) );
    sRP = drfft1.param();  
    Matrix<Bool> sPok( solveParOK().nonDegenerate(1) );
    sPok = (!drfft1.flag());

}


} //# NAMESPACE CASA - END


/*

fringefit(vis="n14c2.ms", caltable="fail.fj", field="",spw="1",intent="",
          selectdata=True, timerange="", antenna="", scan="5", observation="",
          msselect="", solint="inf", refant="EF", minsnr=3.0, append=False,
          gaintable=['n14c2.gcal'], parang=False)
*/
