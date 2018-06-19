//# TestUtilsTVI.h This file contains the interface definition of the TestUtilsTVI class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#ifndef TestUtilsTVI_H_
#define TestUtilsTVI_H_

// Google test
#include <gtest/gtest.h>

// casacore containers
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Record.h>

// Measurement Set enumerations
#include <mstransform/MSTransform/MSTransformManager.h>

// VI/VB framework
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>


namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// FreqAxisTVITest class
//////////////////////////////////////////////////////////////////////////
class FreqAxisTVITest: public ::testing::Test {

public:

	FreqAxisTVITest();
	FreqAxisTVITest(casacore::Record configuration);
    virtual ~FreqAxisTVITest();

    void SetUp();
    void TearDown();
    casacore::Bool getTestResult() {return testResult_p;}

protected:

    void init(casacore::Record &configuration);
    virtual void generateTestFile() = 0;
    virtual void generateReferenceFile() = 0;
    virtual void initTestConfiguration(casacore::Record &configuration) = 0;
    virtual void initReferenceConfiguration(casacore::Record &configuration) = 0;

    casacore::Bool autoMode_p;
    casacore::Bool testResult_p;
    casacore::String inpFile_p;
    casacore::String testFile_p;
    casacore::String referenceFile_p;
    casacore::Record refConfiguration_p;
    casacore::Record testConfiguration_p;
};


//////////////////////////////////////////////////////////////////////////
// Convenience methods
//////////////////////////////////////////////////////////////////////////
template <class T> void compareVector(const casacore::Char* column,
										const casacore::Vector<T> &inp,
										const casacore::Vector<T> &ref,
										const casacore::Vector<casacore::uInt> &rowIds,
										casacore::Float tolerance = FLT_EPSILON);

template <class T> void compareMatrix(const casacore::Char* column,
										const casacore::Matrix<T> &inp,
										const casacore::Matrix<T> &ref,
										const casacore::Vector<casacore::uInt> &rowIds,
										casacore::Float tolerance = FLT_EPSILON);

template <class T> void compareCube(const casacore::Char* column,
									const casacore::Cube<T> &inp,
									const casacore::Cube<T> &ref,
									const casacore::Vector<casacore::uInt> &rowIds,
									casacore::Float tolerance = FLT_EPSILON);

void compareVisibilityIterators(VisibilityIterator2 &testTVI,
								VisibilityIterator2 &refTVI,
								VisBufferComponents2 &columns,
								casacore::Float tolerance = FLT_EPSILON,
								dataColMap *datacolmap = NULL);

void copyTestFile(casacore::String &path,casacore::String &filename,casacore::String &outfilename);

const casacore::Cube<casacore::Complex> & getViscube(VisBuffer2 *vb,
									casacore::MS::PredefinedColumns datacol,
									dataColMap *datacolmap);

void flagEachOtherChannel(VisibilityIterator2 &vi);

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif /* TestUtilsTVI_H_ */
