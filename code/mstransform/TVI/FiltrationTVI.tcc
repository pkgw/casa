//# FiltrationTVI.tcc: Template class for data filtering TVI
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the Implied warranty of MERCHANTABILITY or
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
#ifndef _MSVIS_FILTRATIONTVI_TCC_
#define _MSVIS_FILTRATIONTVI_TCC_

//#include <mstransform/TVI/FiltrationTVI.h>

#include <climits>

#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Arrays/Array.h>

#include <msvis/MSVis/VisibilityIterator2.h>

using namespace casacore;

namespace {
template<class T>
inline void FiltrateVector(Vector<T> const &feed,
    Vector<bool> const &is_filtrate, Vector<T> &filtrate) {
  AlwaysAssert(feed.conform(is_filtrate), AipsError);
  // filter_flag: true --> filtrate, false --> residue
  filtrate.resize(ntrue(is_filtrate));
  Int k = 0;
  for (size_t i = 0; i < feed.nelements(); ++i) {
    if (is_filtrate[i]) {
      filtrate[k] = feed[i];
      ++k;
    }
  }
}

template<class T, class Func>
inline void FiltrateVector2(Func feeder, Vector<bool> const &is_filtrate,
    Vector<T> &filtrate) {
  Vector<T> feed;
  feeder(feed);
  FiltrateVector(feed, is_filtrate, filtrate);
}

template<class T>
inline void FiltrateMatrix(Matrix<T> const &feed,
    Vector<bool> const &is_filtrate, Matrix<T> &filtrate) {
  AlwaysAssert(feed.conform(is_filtrate), AipsError);
  // filter_flag: true --> filtrate, false --> residue
  filtrate.resize(ntrue(is_filtrate));
  Int k = 0;
  for (size_t i = 0; i < feed.nelements(); ++i) {
    if (is_filtrate[i]) {
      filtrate.column(k) = feed.column(i);
      ++k;
    }
  }
}

template<class T, class Func>
inline void FiltrateMatrix2(Func feeder, Vector<bool> const &is_filtrate,
    Matrix<T> &filtrate) {
  Matrix<T> feed;
  feeder(feed);
  FiltrateMatrix(feed, is_filtrate, filtrate);
}

template<class T>
inline void FiltrateCube(Cube<T> const &feed, Vector<bool> const &is_filtrate,
    Cube<T> &filtrate) {
  AlwaysAssert(feed.conform(is_filtrate), AipsError);
  // filter_flag: true --> filtrate, false --> residue
  filtrate.resize(ntrue(is_filtrate));
  Int k = 0;
  for (size_t i = 0; i < feed.nelements(); ++i) {
    if (is_filtrate[i]) {
      filtrate.xyPlane(k) = feed.xyPlane(i);
      ++k;
    }
  }
}

template<class T, class Func>
inline void FiltrateCube2(Func feeder, Vector<bool> const &is_filtrate,
    Cube<T> &filtrate) {
  Cube<T> feed;
  feeder(feed);
  FiltrateCube(feed, is_filtrate, filtrate);
}

template<class T>
inline void FiltrateArray(Array<T> const &feed, Vector<bool> const &is_filtrate,
    Array<T> &filtrate) {
  // filter_flag: true --> filtrate, false --> residue
  filtrate.resize(ntrue(is_filtrate));
  Int k = 0;
  uInt const ndim = feed.ndim();
  AlwaysAssert(feed.shape()[ndim - 1] == is_filtrate.nelements(), AipsError);
  IPosition iter_axis(1, ndim - 1);
  ArrayIterator<T> from_iter(feed, iter_axis, False);
  ArrayIterator<T> to_iter(filtrate, iter_axis, False);
  for (size_t i = 0; i < is_filtrate.nelements(); ++i) {
    if (is_filtrate[i]) {
      to_iter.array() = from_iter.array();
      to_iter.next();
    }
    from_iter.next();
  }
}

template<class T, class Func>
inline void FiltrateArray2(Func feeder, Vector<bool> const &is_filtrate,
    Array<T> &filtrate) {
  Array<T> feed;
  feeder(feed);
  FiltrateArray(feed, is_filtrate, filtrate);
}
}

namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE vi - BEGIN

// forward declaration
template<class Filter>
class FiltrationTVI;

// FiltrationTVI implementation
template<class Filter>
FiltrationTVI<Filter>::FiltrationTVI(ViImplementation2 * inputVi,
    Filter *filter) :
    TransformingVi2(inputVi), filter_p(filter), num_filtrates_p(0), is_filtrate_p() {
  // Initialize attached VisBuffer
  setVisBuffer(createAttachedVisBuffer(VbPlain, VbRekeyable));
}

template<class Filter>
FiltrationTVI<Filter>::~FiltrationTVI() {
}

template<class Filter>
void FiltrationTVI<Filter>::origin() {
  auto const vii = getVii();
  vii->origin();

  // Synchronize own VisBuffer -- is it required?
  //configureNewSubchunk();

  // filtration
  filter();
}

template<class Filter>
void FiltrationTVI<Filter>::next() {
  // filtration
  filter();
}

template<class Filter>
Int FiltrationTVI<Filter>::nRows() const {
  return num_filtrates_p;
}

template<class Filter>
void FiltrationTVI<Filter>::getRowIds(Vector<uInt> &rowids) const {
  ::FiltrateVector2(getVii()->getRowIds, is_filtrate_p, rowids);
//  Vector<uInt> rowids_org;
//  getVii()->getRowIds(rowids_org);
//  ::FiltrateVector(rowids_org, is_filtrate_p, rowids);
}

template<class Filter>
void FiltrationTVI<Filter>::antenna1(Vector<Int> &ant1) const {
  ::FiltrateVector2(getVii()->antenna1, is_filtrate_p, ant1);
//  Vector<Int> ant1_org;
//  getVii()->antenna1(ant1_org);
//  ::FiltrateVector(ant1_org, is_filtrate_p, ant1);
}

template<class Filter>
void FiltrationTVI<Filter>::antenna2(Vector<Int> &ant2) const {
  ::FiltrateVector2(getVii()->antenna2, is_filtrate_p, ant2);
//  Vector<Int> ant2_org;
//  getVii()->antenna2(ant2_org);
//  ::FiltrateVector(ant2_org, is_filtrate_p, ant2);
}

template<class Filter>
void FiltrationTVI<Filter>::corrType(Vector<Int> &corrTypes) const {
  ::FiltrateVector2(getVii()->corrType, is_filtrate_p, corrTypes);
}

template<class Filter>
void FiltrationTVI<Filter>::exposure(Vector<double> &expo) const {
  ::FiltrateVector2(getVii()->exposure, is_filtrate_p, expo);
}

template<class Filter>
void FiltrationTVI<Filter>::feed1(Vector<Int> &fd1) const {
  ::FiltrateVector2(getVii()->feed1, is_filtrate_p, fd1);
}

template<class Filter>
void FiltrationTVI<Filter>::feed2(Vector<Int> &fd2) const {
  ::FiltrateVector2(getVii()->feed2, is_filtrate_p, fd2);
}

template<class Filter>
void FiltrationTVI<Filter>::fieldIds(Vector<Int> &fld) const {
  ::FiltrateVector2(getVii()->fieldIds, is_filtrate_p, fld);
}

template<class Filter>
void FiltrationTVI<Filter>::arrayIds(Vector<Int> &arr) const {
  ::FiltrateVector2(getVii()->arrayIds, is_filtrate_p, arr);
}

template<class Filter>
void FiltrationTVI<Filter>::flag(Cube<Bool> &flags) const {
  ::FiltrateCube2(getVii()->flag, is_filtrate_p, flags);
}

template<class Filter>
void FiltrationTVI<Filter>::flag(Matrix<Bool> &flags) const {
  ::FiltrateMatrix2(getVii()->flag, is_filtrate_p, flags);
}

template<class Filter>
void FiltrationTVI<Filter>::flagCategory(Array<Bool> & flagCategories) const {
  ::FiltrateArray2(getVii()->flagCategory, is_filtrate_p, flagCategories);
}

template<class Filter>
void FiltrationTVI<Filter>::flagRow(Vector<Bool> &rowflags) const {
  ::FiltrateVector2(getVii()->falgRow, is_filtrate_p, rowflags);
}

template<class Filter>
void FiltrationTVI<Filter>::observationId(Vector<Int> &obsids) const {
  ::FiltrateVector2(getVii()->observationId, is_filtrate_p, obsids);
}

template<class Filter>
void FiltrationTVI<Filter>::processorId(Vector<Int> &procids) const {
  ::FiltrateVector2(getVii()->processorId, is_filtrate_p, procids);
}

template<class Filter>
void FiltrationTVI<Filter>::scan(Vector<Int> & scans) const {
  ::FiltrateVector2(getVii()->scan, is_filtrate_p, scans);
}

template<class Filter>
void FiltrationTVI<Filter>::stateId(Vector<Int> & stateids) const {
  ::FiltrateVector2(getVii()->stateId, is_filtrate_p, stateids);
}

template<class Filter>
void FiltrationTVI<Filter>::jonesC(
      Vector<SquareMatrix<Complex, 2> > & cjones) const {
  ::FiltrateVector2(getVii()->jonesC, is_filtrate_p, cjones);
}

template<class Filter>
void FiltrationTVI<Filter>::sigma(Matrix<Float> & sigmat) const {
  ::FiltrateMatrix2(getVii()->sigma, is_filtrate_p, sigmat);
}

template<class Filter>
void FiltrationTVI<Filter>::spectralWindows(Vector<Int> & spws) const {
  ::FiltrateVector2(getVii()->spectralWindows, is_filtrate_p, spws);
}

template<class Filter>
void FiltrationTVI<Filter>::time(Vector<double> & t) const {
  ::FiltrateVector2(getVii()->time, is_filtrate_p, t);
}

template<class Filter>
void FiltrationTVI<Filter>::timeCentroid(Vector<double> & t) const {
  ::FiltrateVector2(getVii()->timeCentroid, is_filtrate_p, t);
}

template<class Filter>
void FiltrationTVI<Filter>::timeInterval(Vector<double> & ti) const {
  ::FiltrateVector2(getVii()->timeInterval, is_filtrate_p, ti);
}

template<class Filter>
void FiltrationTVI<Filter>::uvw(Matrix<double> & uvwmat) const {
  ::FiltrateMatrix2(getVii()->uvw, is_filtrate_p, uvwmat);
}

template<class Filter>
void FiltrationTVI<Filter>::visibilityCorrected(
  Cube<Complex> & vis) const {
  ::FiltrateCube2(getVii()->visibilityCorrected, is_filtrate_p, vis);
}

template<class Filter>
void FiltrationTVI<Filter>::visibilityModel(Cube<Complex> & vis) const {
  ::FiltrateCube2(getVii()->visibilityModel, is_filtrate_p, vis);
}

template<class Filter>
void FiltrationTVI<Filter>::visibilityObserved(
  Cube<Complex> & vis) const {
  ::FiltrateCube2(getVii()->visibilityObserved, is_filtrate_p, vis);
}

template<class Filter>
void FiltrationTVI<Filter>::floatData(Cube<Float> & fcube) const {
  ::FiltrateCube2(getVii()->floatData, is_filtrate_p, fcube);
}

template<class Filter>
IPosition FiltrationTVI<Filter>::visibilityShape() const {
  IPosition shape = getVii()->visibilityShape();
  shape[shape.nelements() - 1] = num_filtrates_p;
  return shape;
}

template<class Filter>
void FiltrationTVI<Filter>::weight(Matrix<Float> & wtmat) const {
  ::FiltrateMatrix2(getVii()->weight, is_filtrate_p, wtmat);
}

template<class Filter>
void FiltrationTVI<Filter>::weightSpectrum(Cube<Float> & wtsp) const {
  ::FiltrateCube2(getVii()->weightSpectrum, is_filtrate_p, wtsp);
}

template<class Filter>
void FiltrationTVI<Filter>::sigmaSpectrum(Cube<Float> & sigmasp) const {
  ::FiltrateCube2(getVii()->sigmaSpectrum, is_filtrate_p, sigmasp);
}

template<class Filter>
void FiltrationTVI<Filter>::dataDescriptionIds(Vector<Int> &ddids) const {
  ::FiltrateVector2(getVii()->dataDescriptionIds, is_filtrate_p, ddids);
}



template<class Filter>
void FiltrationTVI<Filter>::filter() {
  auto const vii = getVii();
  auto const vb = vii->getVisBuffer();
  for (; vii->more() && filter_p->isResidue(vb); vii->next()) {
    // Synchronize own VisBuffer -- is it required to do inside the loop?
    //configureNewSubchunk();
  }

  // Synchronize own VisBuffer
  configureNewSubchunk();

  // update filter information
  num_filtrates_p = filter_p->isFiltratePerRow(vb, is_filtrate_p);
}

} //# NAMESPACE vi - END

} //# NAMESPACE CASA - END

#endif // _MSVIS_FILTRATIONTVI_TCC_
