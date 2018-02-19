#ifndef POINTINGINTERPOLATIONTVI_H_
#define POINTINGINTERPOLATIONTVI_H_



#include <casacore/casa/aips.h>

#include <casa/Containers/Record.h>

#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/ViImplementation2.h>
#include <msvis/MSVis/TransformingVi2.h>

#include <measures/Measures/MDirection.h>


namespace casa {

namespace vi {

using casacore::Vector;
using casacore::Double;

class PointingInterpolationTVI: public TransformingVi2 {

public:
	PointingInterpolationTVI(ViImplementation2 * inputVi);
	~PointingInterpolationTVI();

	virtual casacore::String ViiType() const ;

	// Sub-chunk navigation methods
	virtual void origin();
	virtual void next();

	virtual std::pair<bool, casacore::MDirection> getPointingAngle (int antenna, double time) const;

private:
	casacore::Unit lonUnit_;
	casacore::Unit latUnit_;
	casacore::MDirection::Ref dirRef_;
	void setupInterpolator();
	// Utilities
	casacore::String taQLSet(const std::set<int> &);

	// Interpolator based on SDPosInterpolator
	// TODO: handle corner cases (nPointings < 4) properly
public:
	class Interpolator {
	public:
		using PointingTimes = Vector<Double>;      // pointingsCount[antId]
		using PointingDir   = Vector<Double>;      // 2
		using PointingDirs  = Vector<PointingDir>; // pointingsCount[antId]
		Interpolator(){}

		void setData(const Vector<PointingTimes> &antsTimes, // nAnt
					 const Vector<PointingDirs> &antsDirs,   // nAnt
					 const Vector<bool> &antSelected);       // nAnt

		// Interpolation
		Vector<Double> pointingDir(int ant,double time) const;
	private:
		void clear();
		void init(const Vector<bool> &antSelected = Vector<bool>());
		using SplineSegmentCoeffs = Vector<Double>;           // 4
		using DirSegmentCoeffs = Vector<SplineSegmentCoeffs>; // 2
		using SplineCoeffs = Vector<DirSegmentCoeffs>;        // pointingsCount[antId] - 1
		void computeSplineCoeffs(const PointingTimes& timeStamps,
							 const PointingDirs& dirs,
							 SplineCoeffs& coeffs);
		Vector<PointingTimes> antsTimes_;        // nAnt
		Vector<bool> isSelected_;                // nAnt
		Vector<bool> isInterpolated_;            // nAnt
		Vector<SplineCoeffs> antsSplinesCoeffs_; // nAnt

	};

private:
	Interpolator interpolator_;

};

class PointingInterpolationVi2Factory: public ViFactory {

public:
	// Constructor
	PointingInterpolationVi2Factory(casacore::Record const &configuration,
			ViImplementation2 *inputVII);
	PointingInterpolationVi2Factory(casacore::Record const &configuration,
			const casacore::MeasurementSet *ms, const SortColumns &sortColumns,
			casacore::Double timeInterval, casacore::Bool isWritable);

	// Destructor
	~PointingInterpolationVi2Factory();

	ViImplementation2 * createVi() const;

private:

	ViImplementation2 *inputVII_p;
	casacore::Record configuration_;

};

class PointingInterpolationTVILayerFactory: public ViiLayerFactory {

public:
	PointingInterpolationTVILayerFactory(casacore::Record const &configuration);
	virtual ~PointingInterpolationTVILayerFactory();

protected:

	ViImplementation2 * createInstance(ViImplementation2* vii0) const;

	casacore::Record configuration_p;

};

} // namespace vi

} // namespace casa

#endif /* POINTINGINTERPOLATIONTVI_H_ */
