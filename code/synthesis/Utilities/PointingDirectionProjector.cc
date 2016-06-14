
#include <synthesis/Utilities/PointingDirectionProjector.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

using namespace casacore;

namespace casa {
OrthographicProjector::~OrthographicProjector()
{
  // Do nothing
}

OrthographicProjector::OrthographicProjector(Float pixel_scale)
  : Projector(),
  pixel_scale_(pixel_scale),
  p_center_(2,0.0),
  p_size_(2,0.0)

{
}


void Projector::setDirection( const Matrix<Double> &dir )
{
  dir_.reference(dir.copy());
  Vector<Double> ra( dir_.row(0) ) ;
  rotateRA( ra ) ;
}


void Projector::rotateRA( Vector<Double> &v )
{
  uInt len = v.nelements() ;
  Vector<Double> work( len ) ;

  for ( uInt i = 0 ; i < len ; i++ ) {
  work[i] = fmod( v[i], C::_2pi ) ;
  if ( work[i] < 0.0 ) {
    work[i] += C::_2pi ;
  }
  }

  Vector<uInt> quad( len ) ;
  Vector<uInt> nquad( 4, 0 ) ;
  for ( uInt i = 0 ; i < len ; i++ ) {
  uInt q = uInt( work[i] / C::pi_2 ) ;
  nquad[q]++ ;
  quad[i] = q ;
  }

  Vector<Bool> rot( 4, False ) ;
  if ( nquad[0] > 0 && nquad[3] > 0
     && ( nquad[1] == 0 || nquad[2] == 0 ) ) {
  //cout << "need rotation" << endl ;
  rot[3] = True ;
  rot[2] = (nquad[1]==0 && nquad[2]>0) ;
  }

  for ( uInt i = 0 ; i < len ; i++ ) {
  if ( rot[quad[i]] ) {
    v[i] = work[i] - C::_2pi ;
  }
  else {
    v[i] = work[i] ;
  }
  }
}


const Matrix<Double>& OrthographicProjector::project()
{
  scale_and_center();
  // using DirectionCoordinate
  Matrix<Double> identity(2,2,Double(0.0)) ;
  identity.diagonal() = 1.0 ;
  DirectionCoordinate coord( MDirection::J2000,
               Projection( Projection::SIN ),
               cenx_, ceny_,
               dx_, dy_,
               identity,
               pcenx_, pceny_);

  Double *pdir_p = new Double[dir_.nelements()] ;
  pdir_.takeStorage( dir_.shape(), pdir_p, TAKE_OVER ) ;
  uInt len = dir_.ncolumn() ;
  Bool b ;
  Double *dir_p = dir_.getStorage( b ) ;
  Double *wdir_p = dir_p ;
  Vector<Double> world ;
  Vector<Double> pixel ;
  IPosition vshape( 1, 2 ) ;
  for ( uInt i = 0 ; i < len ; i++ ) {
    world.takeStorage( vshape, wdir_p, SHARE ) ;
    pixel.takeStorage( vshape, pdir_p, SHARE ) ;
    coord.toPixel( pixel, world ) ;
    pdir_p += 2 ;
    wdir_p += 2 ;
  }
  dir_.putStorage( dir_p, b ) ;
  return pdir_;
}

void OrthographicProjector::scale_and_center()
{
    os_.origin(LogOrigin( "OrthographicProjector", "scale_and_center", WHERE )) ;

    Double xmax, xmin, ymax, ymin ;
    minMax( xmin, xmax, dir_.row( 0 ) ) ;
    minMax( ymin, ymax, dir_.row( 1 ) ) ;
    Double wx = ( xmax - xmin ) * 1.1 ;
    Double wy = ( ymax - ymin ) * 1.1 ;

    cenx_ = 0.5 * ( xmin + xmax ) ;
    ceny_ = 0.5 * ( ymin + ymax ) ;
    Double decCorr = cos( ceny_ ) ;

    // Renaud: uInt len = time_.nelements() ;
    uInt len = dir_.ncolumn();
    Matrix<Double> dd = dir_.copy() ;
    for ( uInt i = len-1 ; i > 0 ; i-- ) {
      //dd(0,i) = ( dd(0,i) - dd(0,i-1) ) * decCorr ;
      dd(0,i) = ( dd(0,i) - dd(0,i-1) ) * cos( 0.5*(dd(1,i-1)+dd(1,i)) ) ;
      dd(1,i) = dd(1,i) - dd(1,i-1) ;
    }
    Vector<Double> dr( len-1 ) ;
    Bool b ;
    const Double *dir_p = dd.getStorage( b ) ;
    const Double *x_p = dir_p + 2 ;
    const Double *y_p = dir_p + 3 ;
    for ( uInt i = 0 ; i < len-1 ; i++ ) {
      dr[i] = sqrt( (*x_p) * (*x_p) + (*y_p) * (*y_p) ) ;
      x_p += 2 ;
      y_p += 2 ;
    }
    dir_.freeStorage( dir_p, b ) ;
    Double med = median( dr, False, True, True ) ;
    dy_ = med * pixel_scale_ ;
    dx_ = dy_ / decCorr ;

    Double nxTemp = ceil(wx / dx_);
    Double nyTemp = ceil(wy / dy_);

    os_ << LogIO::DEBUGGING
      << "len = " << len
        << "range x = (" << xmin << "," << xmax << ")" << endl
        << "range y = (" << ymin << "," << ymax << ")" << endl
      << "direction center = (" << cenx_ << "," << ceny_ << ")" << endl
      << "declination correction: cos(dir_center.y)=" << decCorr << endl
        << "median separation between pointings: " << med << endl
        << "dx=" << dx_ << ", dy=" << dy_ << endl
        << "wx=" << wx  << ", wy=" << wy  << endl
      << "nxTemp=" << nxTemp  << ", nyTemp=" << nyTemp  << LogIO::POST ;

    if (nxTemp > (Double)UINT_MAX || nyTemp > (Double)UINT_MAX) {
      throw AipsError("Error in setup: Too large number of pixels.");
    }
    nx_ = uInt( nxTemp ) ;
    ny_ = uInt( nyTemp ) ;

    // Renaud debug
    p_size_[0] = nxTemp;
    p_size_[1] = nyTemp;

    pcenx_ = 0.5 * Double( nx_ - 1 ) ;
    pceny_ = 0.5 * Double( ny_ - 1 ) ;

    // Renaud debug
    p_center_[0] = pcenx_;
    p_center_[1] = pceny_;

    os_ << LogIO::DEBUGGING
      << "pixel center = (" << pcenx_ << "," << pceny_ << ")" << endl
        << "nx=" << nx_ << ", ny=" << ny_
      << "n_pointings=" << len << " must be < n_pixels=" << nx_ * ny_ << LogIO::POST ;
}

}
