/*
 * SideBandSeparator.h
 *
 *  Created on: 2017/07/19
 *      Author: kana
 */

#ifndef SIDEBANDSEPARATOR_H_
#define SIDEBANDSEPARATOR_H_
// STL
#include <iostream>
#include <string>
#include <vector>
// casacore
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>

#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/scimath/Mathematics/FFTServer.h>

#include <imageanalysis/ImageTypedefs.h>

namespace casa { //# NAMESPACE CASA - BEGIN

class SimpleSideBandSeparator {
public:

	  /**
	   * constructors and a destructor
	   **/
	  SimpleSideBandSeparator(const std::vector<std::string>& imagename);
	  virtual ~SimpleSideBandSeparator();

	  /**
	   * Set the number of channels shifted in image side band
	   * of each of scantable.
	   **/
	  void setShift(const std::vector<double> &shift, const bool signal = true);

	  /**
	   * Set rejection limit of solution.
	   **/
	  void setThreshold(const double limit);

	  /**
	   * Resolve both image and signal sideband when true is set.
	   **/
	  void solveBoth(const bool flag) { doboth_ = flag; };

	  /**
	   * Obtain spectra by subtracting the solution of the other sideband.
	   **/
	  void solvefromOther(const bool flag) { otherside_ = flag; };

	  /**
	  * invoke sideband separation
	  **/
	  void separate(const std::string& outfile, const bool overwrite);

protected:
	  /** Initialize member variables **/
	  void init();
	  void initshift();
	  void setImage(const std::vector<std::string>& imagename);

	  /** Return if the path exists (optionally, check file type) **/
	  casacore::Bool checkFile(const std::string name, std::string type="");

	  std::vector<float> solve(const casacore::Matrix<float> &specMat,
			      const std::vector<casacore::uInt> &tabIdvec,
			      const bool signal = true);

	  casacore::Vector<bool> collapseMask(const casacore::Matrix<bool> &flagMat,
				    const std::vector<casacore::uInt> &tabIdvec,
				    const bool signal = true);
	  void shiftSpectrum(const casacore::Vector<float> &invec, double shift,
			     casacore::Vector<float> &outvec);

	  void shiftFlag(const casacore::Vector<bool> &invec, double shift,
			     casacore::Vector<bool> &outvec);

	  void deconvolve(casacore::Matrix<float> &specmat, const std::vector<double> shiftvec,
			  const double threshold, casacore::Matrix<float> &outmat);

	  void aggregateMat(casacore::Matrix<float> &inmat, std::vector<float> &outvec);

	  void subtractFromOther(const casacore::Matrix<float> &shiftmat,
				 const std::vector<float> &invec,
				 const std::vector<double> &shift,
				 std::vector<float> &outvec);
	  ////
	  size_t setupShift();
	  bool getImageCoordinate(const string& imagename, CoordinateSystem &csys, IPosition &npix);
	  bool compareImageAxes(const string& imagename, const CoordinateSystem &refcsys, const IPosition &refnpix);
	  bool getSpectraToSolve(const vector<SPIIF> &images, const Slicer &slicer,
			  Matrix<float>& specMat, Matrix<bool>& maskMat, vector<uInt>& imgIdvec);

	  /** Member variables **/
	  // name of images
	  std::vector<std::string> imageNames_;
	  // frequency and direction setup to select data.
	  std::vector<double> sigShift_, imgShift_;
	  unsigned int nshift_, nchan_;
	  // solution parameters
	  bool otherside_, doboth_;
	  double rejlimit_;

	  casacore::FFTServer<casacore::Float, casacore::Complex> fftsf, fftsi;


};



} //# NAMESPACE CASA - END
#endif /* SIDEBANDSEPARATOR_H_ */
