//# FiltrationTVI_GTest:   test of Filtration TVIs
//# Copyright (C) 1995,1999,2000,2001,2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
#include <gtest/gtest.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <cmath>

#include <mstransform/TVI/FiltrationTVI.h>
#include <mstransform/TVI/FiltrationTVI.tcc>

#include <casacore/casa/aips.h>
#include <casacore/casa/OS/EnvVar.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/RegularFile.h>
#include <casacore/casa/OS/SymLink.h>
#include <casacore/casa/OS/Directory.h>
#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/iostream.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/iomanip.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/GenSort.h>

#include <casacore/tables/Tables/ArrColData.h>
#include <casacore/tables/DataMan/TiledShapeStMan.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <msvis/MSVis/VisibilityIteratorImpl2.h>
#include <msvis/MSVis/LayeredVi2Factory.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;

namespace {
string GetCasaDataPath() {
  if (casacore::EnvironmentVariable::isDefined("CASAPATH")) {
    string casapath = casacore::EnvironmentVariable::get("CASAPATH");
    size_t endindex = casapath.find(" ");
    if (endindex != string::npos) {
      string casaroot = casapath.substr(0, endindex);
      cout << "casaroot = " << casaroot << endl;
      return (casaroot + "/data/");
    } else {
      cout << "hit npos" << endl;
      return "/data/";
    }
  } else {
    cout << "CASAPATH is not defined" << endl;
    return "";
  }
}

template<class T>
struct VerboseDeleterForNew {
  void operator()(T *p) {
    cout << "Destructing " << typeid(p).name() << endl;
    delete p;
  }
};

// dummy
class FiltrationTestTVIFactory;
template<class Filter>
class FiltrationTVIWrapper: public FiltrationTVI<Filter> {
  FiltrationTVIWrapper(ViImplementation2 * inputVi, Filter *filter) :
      FiltrationTVI<Filter>(inputVi, filter) {
  }

  virtual ~FiltrationTVIWrapper() {
  }

private:
  friend class FiltrationTestTVIFactory;
};

// Filter class for testing
class FilterTypeLocal {
public:
  enum {
    Porous = -1, Nonporous = -2
  };
};

/**
 * PorousFilter
 *
 * PorousFilter is an implementaiton of the filter that pass through
 * everything. It is equivalent to the case when no filter is
 * inserted to the TVI layer.
 */
class PorousFilter {
public:
  // constructor
  PorousFilter(MeasurementSet const &/*ms*/, Record const &/*configuration*/) {
  }

  // destructor
  ~PorousFilter() {
  }

  // return string representation of the filter type
  String filterType() const {
    return String("Porous");
  }

  // filter query
  // isResidue returns true if given vb doesn't pass through the filter
  bool isResidue(VisBuffer2 const *vb) {
    return !isFiltrate(vb);
  }

  // isFiltrate returns true if given vb does pass through the filter
  // (either fully and partly)
  bool isFiltrate(VisBuffer2 const *vb) {
    return true;
  }

  // row-wise filtration information
  // it fills in is_filtrate vector (resize if necessary)
  // and returns number of rows that pass through the filter
  int isFiltratePerRow(VisBuffer2 const *vb, Vector<bool> &is_filtrate) {
    int nrows = vb->nRows();
    is_filtrate.resize(nrows);
    is_filtrate = true;
    return nrows;
  }

private:
  void initFilter() {
  }
};

/**
 * NonporousFilter
 *
 * NonporousFilter is an implementaiton of the filter that filter out
 * everything. No data will be emerged if it is inserted to the TVI
 * layer.
 */
class NonporousFilter {
public:
  // constructor
  NonporousFilter(MeasurementSet const &/*ms*/,
      Record const &/*configuration*/) {
  }

  // destructor
  ~NonporousFilter() {
  }

  // return string representation of the filter type
  String filterType() const {
    return String("Nonporous");
  }

  // filter query
  // isResidue returns true if given vb doesn't pass through the filter
  bool isResidue(VisBuffer2 const *vb) {
    return !isFiltrate(vb);
  }

  // isFiltrate returns true if given vb does pass through the filter
  // (either fully and partly)
  bool isFiltrate(VisBuffer2 const *vb) {
    return false;
  }

  // row-wise filtration information
  // it fills in is_filtrate vector (resize if necessary)
  // and returns number of rows that pass through the filter
  int isFiltratePerRow(VisBuffer2 const *vb, Vector<bool> &is_filtrate) {
    int nrows = vb->nRows();
    is_filtrate.resize(nrows);
    is_filtrate = false;
    return nrows;
  }

private:
  void initFilter() {
  }
};

class FiltrationTestTVIFactory: public ViFactory {

public:
  // Constructor
  FiltrationTestTVIFactory(Record const &configuration,
      ViImplementation2 *inputVII) :
      inputVII_p(inputVII), configuration_p(configuration) {
  }

  // Destructor
  ~FiltrationTestTVIFactory() {
  }

  virtual ViImplementation2 * createVi() const {
    ViImplementation2 *vii = nullptr;

    Bool is_porous = true;
    if (configuration_p.isDefined("type")
        && configuration_p.asInt("type") == (Int) FilterTypeLocal::Nonporous) {
      is_porous = false;
    }
    cout << "type_enum = " << configuration_p.asInt("type") << endl;
    cout << "is_porous = " << is_porous << endl;

    MeasurementSet const &ms = inputVII_p->ms();

    if (is_porous) {
      PorousFilter *filter = new PorousFilter(ms, configuration_p);
      vii = new FiltrationTVIWrapper<PorousFilter>(inputVII_p, filter);
    } else {
      NonporousFilter *filter = new NonporousFilter(ms, configuration_p);
      vii = new FiltrationTVIWrapper<NonporousFilter>(inputVII_p, filter);
    }

    return vii;
  }

private:
  ViImplementation2 *inputVII_p;
  Record configuration_p;
};

//template<class T>
//struct Filler {
//  static void FillArrayReference(MeasurementSet const &ms,
//      String const &columnName, Vector<uInt> const &rowIds, Array<T> &data) {
//    //cout << "Start " << __func__ << endl;
//    ArrayColumn<T> col(ms, columnName);
//    RefRows rows(rowIds);
//    col.getColumnCells(rows, data);
//  }
//  static void FillScalarReference(MeasurementSet const &ms,
//      String const &columnName, Vector<uInt> const &rowIds, Vector<T> &data) {
//    ScalarColumn<T> col(ms, columnName);
//    RefRows rows(rowIds);
////    Vector<T> vref(data);
//    col.getColumnCells(rows, data);
//  }
//  static void FillWeightSp(MeasurementSet const &ms, Vector<uInt> const &rowIds,
//      Array<Float> &weightSp) {
//    // weightSp must be resized to appropriate shape
//    if (ms.tableDesc().isColumn("WEIGHT_SPECTRUM")) {
//      Filler<Float>::FillArrayReference(ms, "WEIGHT_SPECTRUM", rowIds,
//          weightSp);
//    } else {
//      Array<Float> weight;
//      Filler<Float>::FillArrayReference(ms, "WEIGHT", rowIds, weight);
//      auto const wshape(weight.shape());
//      ASSERT_EQ(wshape.size(), (uInt )2);
//      ASSERT_EQ(wshape[0], weightSp.shape()[0]);
//      ASSERT_EQ(wshape[1], weightSp.shape()[2]);
//      Matrix<Float> wtMat(weight);
//      Cube<Float> wtCube(weightSp);
//      cout << "wtCube.shape() = " << wtCube.shape() << " weightSp.shape() = "
//          << weightSp.shape() << endl;
//      size_t nPol = wshape[0];
//      size_t nRow = wshape[1];
//      cout << "wtMat.shape() = " << wtMat.shape() << " nrow = " << wtMat.nrow()
//          << " ncol = " << wtMat.ncolumn() << endl;
//      for (size_t i = 0; i < nPol; ++i) {
//        for (size_t j = 0; j < nRow; ++j) {
//          cout << "wtCube.yzPlane(" << i << ").shape() = "
//              << wtCube.yzPlane(i).shape() << endl;
//          cout << "wtCube.yzPlane(i).column(j).shape() = "
//              << wtCube.yzPlane(i).column(j).shape() << endl;
//          cout << "wtMat(i, j) = " << wtMat(i, j) << endl;
//          wtCube.yzPlane(i).column(j) = wtMat(i, j);
//        }
//      }
//    }
//  }
//};

struct ValidatorUtil {
  template<class T>
  static void ValidateArray(Array<T> const &ref, Array<T> const &result) {
    ASSERT_TRUE(ref.conform(result));
    Bool b1, b2;
    auto const p_ref = ref.getStorage(b1);
    auto const p_result = result.getStorage(b2);
    size_t arraySize = ref.size();
    for (size_t i = 0; i < arraySize; ++i) {
      ValidateScalar(p_ref[i], p_result[i]);
    }
  }
private:
  template<class T>
  static void ValidateScalar(T const ref, T const result) {
    EXPECT_EQ(ref, result);
  }

};

template<>
void ValidatorUtil::ValidateScalar<Float>(Float const ref, Float const result) {
  EXPECT_FLOAT_EQ(ref, result);
}

template<>
void ValidatorUtil::ValidateScalar<Complex>(Complex const ref,
    Complex const result) {
  EXPECT_FLOAT_EQ(ref.real(), result.real());
  EXPECT_FLOAT_EQ(ref.imag(), result.imag());
}

// Base class for validating polarization average
//struct ValidatorBase {
//public:
//  static void ValidateData(Array<Complex> const &data, MeasurementSet const &ms,
//      Vector<uInt> const &rowIds) {
//    ValidateDataColumn(data, ms, "DATA", rowIds);
//  }
//
//  static void ValidateCorrected(Cube<Complex> const &data,
//      MeasurementSet const &ms, Vector<uInt> const &rowIds) {
//    ValidateDataColumn(data, ms, "CORRECTED_DATA", rowIds);
//  }
//
//  static void ValidateModel(Cube<Complex> const &data, MeasurementSet const &ms,
//      Vector<uInt> const &rowIds) {
//    ValidateDataColumn(data, ms, "MODEL_DATA", rowIds);
//  }
//
//  static void ValidateFloat(Cube<Float> const &data, MeasurementSet const &ms,
//      Vector<uInt> const &rowIds) {
//    ValidateDataColumn(data, ms, "FLOAT_DATA", rowIds);
//  }
//
//  static void ValidateFlag(Cube<Bool> const &flag, MeasurementSet const &ms,
//      Vector<uInt> const &rowIds) {
//    Array<Bool> baseFlag;
//    Filler<Bool>::FillArrayReference(ms, "FLAG", rowIds, baseFlag);
//    Cube<Bool> ref(flag.shape(), True);
//    IPosition const baseShape = baseFlag.shape();
//    ASSERT_EQ(baseShape.size(), (uInt )3);
//    ASSERT_EQ(baseShape[1], flag.shape()[1]);
//    ASSERT_EQ(baseShape[2], flag.shape()[2]);
//    size_t nPol = baseShape[0];
//    size_t nChan = baseShape[1];
//    size_t nRow = baseShape[2];
//    for (size_t i = 0; i < nPol; ++i) {
//      IPosition start(3, i, 0, 0);
//      IPosition end(3, i, nChan - 1, nRow - 1);
//      auto fslice = baseFlag(start, end);
//      ref &= fslice;
//    }
//    ValidatorUtil::ValidateArray(ref, flag);
//  }
//
//  static void ValidateFlagRow(Vector<Bool> const &flag,
//      MeasurementSet const &ms, Vector<uInt> const &rowIds) {
//    ASSERT_EQ(flag.size(), rowIds.size());
//    Vector<Bool> ref;
//    Filler<Bool>::FillScalarReference(ms, "FLAG_ROW", rowIds, ref);
//    ASSERT_EQ(flag.shape(), ref.shape());
//    EXPECT_TRUE(allEQ(flag, ref));
//  }
//
//  static void ValidateWeight(Matrix<Float> const &weight,
//      MeasurementSet const &ms, Vector<uInt> const &rowIds) {
//    ValidateWeightColumn(weight, ms, "WEIGHT", rowIds);
//  }
//
//  static void ValidateWeightSp(Cube<Float> const &weight,
//      MeasurementSet const &ms, Vector<uInt> const &rowIds) {
//    if (ms.tableDesc().isColumn("WEIGHT_SPECTRUM")) {
//      // WEIGHT_SPECTRUM exists
//      ValidateWeightColumn(weight, ms, "WEIGHT_SPECTRUM", rowIds);
//    } else {
//      // only WEIGHT exists
//      Matrix<Float> scalarWeight = weight.xzPlane(0);
//      ValidateWeightColumn(scalarWeight, ms, "WEIGHT", rowIds);
//      IPosition const weightShape = weight.shape();
//      size_t nPol = weightShape[0];
//      size_t nChan = weightShape[1];
//      size_t nRow = weightShape[2];
//      for (size_t i = 0; i < nPol; ++i) {
//        for (size_t j = 0; j < nRow; ++j) {
//          IPosition start(3, i, 0, j);
//          IPosition end(3, i, nChan - 1, j);
//          auto wslice = weight(start, end);
//          auto wref = scalarWeight(i, j);
//          EXPECT_TRUE(allEQ(wref, wslice));
//        }
//      }
//    }
//  }
//};

struct PorousValidator {
  static Int GetMode() {
    return FilterTypeLocal::Porous;
  }

  static String GetTypePrefix() {
    return "FiltrationTVI<Porous>(";
  }

  static bool IsResidue(VisBuffer2 const *vb) {
    return false;
  }

  static bool IsFiltrate(VisBuffer2 const *vb) {
    return true;
  }
};

struct NonporousValidator {
  static Int GetMode() {
    return FilterTypeLocal::Nonporous;
  }

  static String GetTypePrefix() {
    return "FiltrationTVI<Nonporous>(";
  }

  static bool IsResidue(VisBuffer2 const *vb) {
    return true;
  }

  static bool IsFiltrate(VisBuffer2 const *vb) {
    return false;
  }
};

template<class Impl>
class Manufacturer {
public:
  struct Product {
    ViFactory *factory;
    ViImplementation2 *vii;
  };
  static VisibilityIterator2 *ManufactureVI(MeasurementSet *ms,
      Int const &type_enum) {
    cout << "### Manufacturer: " << endl << "###   " << Impl::GetTestPurpose()
        << endl;

    Record type_rec;
    type_rec.define("type", type_enum);

    // build factory object
    Product p = Impl::BuildFactory(ms, type_rec);
    std::unique_ptr<ViFactory> factory(p.factory);

    std::unique_ptr<VisibilityIterator2> vi;
    try {
      vi.reset(new VisibilityIterator2(*factory.get()));
    } catch (...) {
      cout << "Failed to create VI at factory" << endl;
      // vii must be deleted since it is never managed by vi
      if (p.vii) {
        delete p.vii;
      }
      throw;
    }

    cout << "Created VI type \"" << vi->ViiType() << "\"" << endl;

    return vi.release();
  }
};

class TestManufacturer: public Manufacturer<TestManufacturer> {
public:
  static Product BuildFactory(MeasurementSet *ms, Record const &mode) {
    // create read-only VI impl
    Block<MeasurementSet const *> const mss(1, ms);
    SortColumns defaultSortColumns;

    std::unique_ptr<ViImplementation2> inputVii(
        new VisibilityIteratorImpl2(mss, defaultSortColumns, 0.0, VbPlain,
            False));
    std::unique_ptr<ViFactory> factory(
        new FiltrationTestTVIFactory(mode, inputVii.get()));

    Product p;

    // vi will be responsible for releasing inputVii so unique_ptr
    // should release the ownership here
    p.vii = inputVii.release();
    p.factory = factory.release();

    return p;
  }

  static String GetTestPurpose() {
    return "Test FiltrationTestTVIFactory(Record const &, ViImplementation2 *)";
  }
};

class BasicManufacturer : public Manufacturer<BasicManufacturer> {
public:
  static Product BuildFactory(MeasurementSet *ms, Record const &mode) {
    // create read-only VI impl
    Block<MeasurementSet const *> const mss(1, ms);
    SortColumns defaultSortColumns;

    std::unique_ptr<ViImplementation2> inputVii(
        new VisibilityIteratorImpl2(mss, defaultSortColumns, 0.0, VbPlain,
            False));
    std::unique_ptr<ViFactory> factory(
        new FiltrationTVIFactory(mode, inputVii.get()));

    Product p;

    // vi will be responsible for releasing inputVii so unique_ptr
    // should release the ownership here
    p.vii = inputVii.release();
    p.factory = factory.release();

    return p;

  }

  static String GetTestPurpose() {
    return "Test FiltrationTVIFactory(Record const &, ViImplementation2 *)";
  }
};

class LayerManufacturer: public Manufacturer<LayerManufacturer> {
public:
  class LayerFactoryWrapper: public ViFactory {
  public:
    LayerFactoryWrapper(MeasurementSet *ms, Record const &mode) :
        ms_(ms), mode_(mode) {
    }

    ViImplementation2 *createVi() const {
      Vector<ViiLayerFactory *> v(1);
      auto layer0 = VisIterImpl2LayerFactory(ms_, IteratingParameters(0.0),
          false);
      auto layer1 = FiltrationTVILayerFactory(mode_);
      v[0] = &layer0;
      return layer1.createViImpl2(v);
    }
  private:
    MeasurementSet *ms_;
    Record const mode_;
  };

  static Product BuildFactory(MeasurementSet *ms, Record const &mode) {
    // create read-only VI impl
    std::unique_ptr<ViFactory> factory(new LayerFactoryWrapper(ms, mode));

    Product p;
    p.vii = nullptr;
    p.factory = factory.release();
    return p;
  }

  static String GetTestPurpose() {
    return "Test FiltrationTVILayerFactory";
  }
};

// copy & paste from Calibrater::initWeights
//void initWeights(MeasurementSet *ms) {
//  // add columns
//  TableDesc mstd = ms->actualTableDesc();
//  String colWtSp = MS::columnName(MS::WEIGHT_SPECTRUM);
//  Bool wtspexists = mstd.isColumn(colWtSp);
//  String colSigSp = MS::columnName(MS::SIGMA_SPECTRUM);
//  Bool sigspexists = mstd.isColumn(colSigSp);
//
//  if (!wtspexists) {
//    // Nominal defaulttileshape
//    IPosition dts(3, 4, 32, 1024);
//
//    // Discern DATA's default tile shape and use it
//    const Record dminfo = ms->dataManagerInfo();
//    for (uInt i = 0; i < dminfo.nfields(); ++i) {
//      Record col = dminfo.asRecord(i);
//      //if (upcase(col.asString("NAME"))=="TILEDDATA") {
//      if (anyEQ(col.asArrayString("COLUMNS"), String("DATA"))) {
//        dts = IPosition(col.asRecord("SPEC").asArrayInt("DEFAULTTILESHAPE"));
//        //cout << "Found DATA's default tile: " << dts << endl;
//        break;
//      }
//    }
//
//    // Add the column
//    String colWtSp = MS::columnName(MS::WEIGHT_SPECTRUM);
//    TableDesc tdWtSp;
//    tdWtSp.addColumn(ArrayColumnDesc<Float>(colWtSp, "weight spectrum", 2));
//    TiledShapeStMan wtSpStMan("TiledWgtSpectrum", dts);
//    ms->addColumn(tdWtSp, wtSpStMan);
//  }
//
//  if (!sigspexists) {
//    // Nominal defaulttileshape
//    IPosition dts(3, 4, 32, 1024);
//
//    // Discern DATA's default tile shape and use it
//    const Record dminfo = ms->dataManagerInfo();
//    for (uInt i = 0; i < dminfo.nfields(); ++i) {
//      Record col = dminfo.asRecord(i);
//      //if (upcase(col.asString("NAME"))=="TILEDDATA") {
//      if (anyEQ(col.asArrayString("COLUMNS"), String("DATA"))) {
//        dts = IPosition(col.asRecord("SPEC").asArrayInt("DEFAULTTILESHAPE"));
//        //cout << "Found DATA's default tile: " << dts << endl;
//        break;
//      }
//    }
//
//    // Add the column
//    String colSigSp = MS::columnName(MS::SIGMA_SPECTRUM);
//    TableDesc tdSigSp;
//    tdSigSp.addColumn(ArrayColumnDesc<Float>(colSigSp, "sigma spectrum", 2));
//    TiledShapeStMan sigSpStMan("TiledSigtSpectrum", dts);
//    ms->addColumn(tdSigSp, sigSpStMan);
//    {
//      TableDesc loctd = ms->actualTableDesc();
//      String loccolSigSp = MS::columnName(MS::SIGMA_SPECTRUM);
//      AlwaysAssert(loctd.isColumn(loccolSigSp), AipsError);
//    }
//
//    ArrayColumn<Float> weightColumn(*ms, "WEIGHT");
//    ArrayColumn<Float> sigmaColumn(*ms, "SIGMA");
//    ArrayColumn<Float> weightSpColumn(*ms, "WEIGHT_SPECTRUM");
//    ArrayColumn<Float> sigmaSpColumn(*ms, "SIGMA_SPECTRUM");
//    ArrayColumn<Bool> const flagColumn(*ms, "FLAG");
//    ROScalarColumn<Double> exposureColumn(*ms, "EXPOSURE");
//    for (size_t i = 0; i < ms->nrow(); ++i) {
//      IPosition const cellShape = flagColumn.shape(i);
//      Double const exposure = exposureColumn(i);
//      Matrix<Float> weightSp(cellShape, exposure);
//      Vector<Float> weight(cellShape[0], exposure);
//      Matrix<Float> sigmaSp = 1.0f / sqrt(weightSp);
//      Vector<Float> sigma = 1.0f / sqrt(weight);
//      weightColumn.put(i, weight);
//      sigmaColumn.put(i, sigma);
//      weightSpColumn.put(i, weightSp);
//      sigmaSpColumn.put(i, sigmaSp);
//    }
//
//  }
//}

} // anonymous namespace

// explicitly instantiate test TVI
namespace casa {
namespace vi {
template class FiltrationTVI<PorousFilter> ;
template class FiltrationTVI<NonporousFilter> ;
}
}

class FiltrationTVITestBase: public ::testing::Test {
public:
  FiltrationTVITestBase() :
      my_ms_name_("filtration_test.ms"), my_data_name_(), ms_(nullptr) {
  }

  virtual void SetUp() {
//    my_data_name_ = "analytic_spectra.ms";
    my_data_name_ = GetDataName(); //"analytic_type1.bl.ms";
    std::string const data_path = ::GetCasaDataPath() + "/regression/unittest/"
        + GetRelativeDataPath() + "/";
//        + "/regression/unittest/tsdbaseline/";
//    + "/regression/unittest/singledish/";

    ASSERT_TRUE(Directory(data_path).exists());
    cout << "data_path = " << data_path << endl;
    copyDataFromRepository(data_path);
    ASSERT_TRUE(File(my_data_name_).exists());
    deleteTable(my_ms_name_);

    // create MS
    ms_ = new MeasurementSet(my_data_name_, Table::Update);
  }

  virtual void TearDown() {
    // delete MS explicitly to detach from MS on disk
    delete ms_;

    // just to make sure all locks are effectively released
    Table::relinquishAutoLocks();

    cleanup();
  }

protected:
  std::string const my_ms_name_;
  std::string my_data_name_;
  MeasurementSet *ms_;

  virtual std::string GetDataName() {
    return "";
  }

  virtual std::string GetRelativeDataPath() {
    return "";
  }

  template<class Validator, class Manufacturer = TestManufacturer>
  void TestTVI() {
    cout << "TestTVI" << endl;

    // Create VI
    // VI with filter
    std::unique_ptr<VisibilityIterator2> vi(
        Manufacturer::ManufactureVI(ms_, Validator::GetMode()));
    ASSERT_TRUE(vi->ViiType().startsWith(Validator::GetTypePrefix()));

    // reference VI
    std::unique_ptr<VisibilityIterator2> refvi(new VisibilityIterator2(*ms_));

    // MS property
    auto ms = vi->ms();
    auto refms = refvi->ms();
    uInt const nRowMs = ms.nrow();
    uInt const nRowRefMs = refms.nrow();
    EXPECT_EQ(nRowMs, nRowRefMs);
////    uInt const nRowPolarizationTable = ms.polarization().nrow();
    auto const desc = ms.tableDesc();
    auto const correctedExists = desc.isColumn("CORRECTED_DATA");
    auto const modelExists = desc.isColumn("MODEL_DATA");
    auto const dataExists = desc.isColumn("DATA");
    auto const floatExists = desc.isColumn("FLOAT_DATA");
    //auto const weightSpExists = desc.isColumn("WEIGHT_SPECTRUM");
    cout << "MS Property" << endl;
    cout << "\tMS Name: \"" << ms.tableName() << "\"" << endl;
    cout << "\tNumber of Rows: " << nRowMs << endl;
    cout << "\tNumber of Spws: " << vi->nSpectralWindows() << endl;
    cout << "\tNumber of Polarizations: " << vi->nPolarizationIds() << endl;
    cout << "\tNumber of DataDescs: " << vi->nDataDescriptionIds() << endl;
    cout << "\tChannelized Weight Exists? "
        << (vi->weightSpectrumExists() ? "True" : "False") << endl;
    //cout << "\tChannelized Sigma Exists? " << (vi->sigmaSpectrumExists() ? "True" : "False") << endl;

//    // mv-VI consistency check
////    EXPECT_EQ(nRowPolarizationTable + 1, (uInt )vi->nPolarizationIds());
//
//    // VI iteration
//    Vector<uInt> swept(nRowMs, 0);
//    uInt nRowChunkSum = 0;
    VisBuffer2 *vb = vi->getVisBuffer();
    VisBuffer2 *vb_ref = refvi->getVisBuffer();
    vi->originChunks();
    refvi->originChunks();
    // iteration loop is based on refvi
    while (refvi->moreChunks()) {
      // make sure there is a chunk in vi
      EXPECT_TRUE(vi->moreChunks());

      // initialize subchunk iterator
      refvi->origin();
      vi->origin();

      // increment chunk until refvi iteration hits filtrate subchunk
      while (refvi->more() && Validator::IsResidue(vb_ref)) {
        refvi->next();
      }

      // nRowsInChunk returns number of rows in chunk regardless of
      // whether they are filtered out or not
      Int const nrow_chunk = vi->nRowsInChunk();
      Int const nrow_chunk_ref = refvi->nRowsInChunk();
      EXPECT_EQ(nrow_chunk_ref, nrow_chunk);

      // chunk id should be the same
      auto const chunk_id = vi->getSubchunkId().chunk();
      auto const chunk_id_ref = refvi->getSubchunkId().chunk();
      EXPECT_EQ(chunk_id_ref, chunk_id);
//      nRowChunkSum += nRowChunk;
      cout << "*************************" << endl;
      cout << "*** Start loop on chunk " << chunk_id << endl;
      cout << "*** Number of Rows: " << nrow_chunk << endl;
      cout << "*************************" << endl;
//
//      Int nRowSubchunkSum = 0;
//
      // no valid subchunk exists
      if (!refvi->more()) {
        cout << "No valid chunk exists." << endl;
        EXPECT_FALSE(vi->more());

        refvi->nextChunk();
        vi->nextChunk();
        continue;
      }

      // again, iteration loop is based on refvi
      while (refvi->more()) {
        EXPECT_TRUE(vi->more());

        auto const subchunk = vi->getSubchunkId();
        auto const subchunk_ref = refvi->getSubchunkId();
        auto const subchunk_id = subchunk.subchunk();
        auto const subchunk_id_ref = subchunk_ref.subchunk();
        cout << "=== Start loop on subchunk " << subchunk_id_ref << " ==="
            << endl;

        EXPECT_EQ(subchunk_id_ref, subchunk_id);
//
//        // cannot use getInterval due to the error
//        // "undefined reference to VisibilityIterator2::getInterval"
//        // even if the code is liked to libmsvis.so.
//        //cout << "Interval: " << vi->getInterval() << endl;
//
//        cout << "Antenna1: " << vb->antenna1() << endl;
//        cout << "Antenna2: " << vb->antenna2() << endl;
//        cout << "Array Id: " << vb->arrayId() << endl;
//        cout << "Data Desc Ids: " << vb->dataDescriptionIds() << endl;
//        cout << "Polarization Id: " << vb->polarizationId() << endl;
//        cout << "Exposure: " << vb->exposure() << endl;
//        cout << "Feed1: " << vb->feed1() << endl;
//        cout << "Feed2: " << vb->feed2() << endl;
//        cout << "Field Id: " << vb->fieldId() << endl;
//        cout << "Flag Row: " << vb->flagRow() << endl;
//        cout << "Observation Id: " << vb->observationId() << endl;
//        cout << "Processor Id: " << vb->processorId() << endl;
//        cout << "Scan: " << vb->scan() << endl;
//        cout << "State Id: " << vb->stateId() << endl;
//        cout << "Time: " << vb->time() << endl;
//        cout << "Time Centroid: " << vb->timeCentroid() << endl;
//        cout << "Time Interval: " << vb->timeInterval() << endl;
//        auto const corrTypes = vb->correlationTypes();
//        auto toStokes = [](Vector<Int> const &corrTypes) {
//          Vector<String> typeNames(corrTypes.size());
//          for (size_t i = 0; i < corrTypes.size(); ++i) {
//            typeNames[i] = Stokes::name((Stokes::StokesTypes)corrTypes[i]);
//          }
//          return typeNames;
//        };
//        cout << "Correlation Types: " << toStokes(corrTypes) << endl;
//        //cout << "UVW: " << vb->uvw() << endl;
//
//        cout << "---" << endl;
//        Int nRowSubchunk = vb->nRows();
//        Vector<uInt> rowIds = vb->rowIds();
//        for (auto iter = rowIds.begin(); iter != rowIds.end(); ++iter) {
//          swept[*iter] += 1;
//        }
//        nRowSubchunkSum += nRowSubchunk;
//        Int nAnt = vb->nAntennas();
//        Int nChan = vb->nChannels();
//        Int nCorr = vb->nCorrelations();
//        IPosition visShape = vb->getShape();
//        cout << "Number of Subchunk Rows: " << nRowSubchunk << endl;
//        cout << "Number of Antennas: " << nAnt << endl;
//        cout << "Number of Channels: " << nChan << endl;
//        cout << "Number of Correlations: " << nCorr << endl;
//        cout << "Row Ids: " << rowIds << endl;
//        cout << "Spectral Windows: " << vb->spectralWindows() << endl;
//        cout << "Visibility Shape: " << visShape << endl;
//        cout << "---" << endl;
//        Cube<Complex> visCube = vb->visCube();
//        cout << "DATA Shape: " << visCube.shape() << endl;
//        Cube<Complex> visCubeCorrected = vb->visCubeCorrected();
//        cout << "CORRECTED_DATA Shape: " << visCubeCorrected.shape() << endl;
//        Cube<Complex> visCubeModel = vb->visCubeModel();
//        cout << "MODEL_DATA Shape: " << visCubeModel.shape() << endl;
//        Cube<Float> visCubeFloat = vb->visCubeFloat();
//        cout << "FLOAT_DATA Shape: " << visCubeFloat.shape() << endl;
//        Cube<Bool> flagCube = vb->flagCube();
//        cout << "FLAG Shape: " << flagCube.shape() << endl;
//        Vector<Bool> flagRow = vb->flagRow();
//        cout << "FLAG_ROW Shape: " << flagRow.shape() << endl;
//        Matrix<Float> weight = vb->weight();
//        cout << "WEIGHT Shape: " << weight.shape() << endl;
//        Cube<Float> weightSp = vb->weightSpectrum();
//        cout << "WEIGHT_SPECTRUM Shape: " << weightSp.shape() << endl;
//        cout << "===" << endl;
//
//        // internal consistency check
//        EXPECT_EQ(nRowSubchunk, visShape[2]);
//        EXPECT_EQ(nChan, visShape[1]);
//        EXPECT_EQ(nCorr, visShape[0]);
//        EXPECT_EQ(!dataExists, visCube.empty());
//        if (!visCube.empty()) {
//          EXPECT_EQ(visShape, visCube.shape());
//        }
//        EXPECT_EQ(!correctedExists, visCubeCorrected.empty());
//        if (!visCubeCorrected.empty()) {
//          EXPECT_EQ(visShape, visCubeCorrected.shape());
//        }
//        EXPECT_EQ(!modelExists, visCubeModel.empty());
//        if (!visCubeModel.empty()) {
//          EXPECT_EQ(visShape, visCubeModel.shape());
//        }
//        EXPECT_EQ(!floatExists, visCubeFloat.empty());
//        if (!visCubeFloat.empty()) {
//          EXPECT_EQ(visShape, visCubeFloat.shape());
//        }
//        EXPECT_EQ((ssize_t )nRowSubchunk, weight.shape()[1]);
//        // NB: weight spectrum is created on-the-fly based on WEIGHT
//        //     so that weightSp is always non-empty.
//        //     see VisBufferImpl2::fillWeightSpectrum.
//        //EXPECT_EQ(!weightSpExists, weightSp.empty());
//        EXPECT_FALSE(weightSp.empty());
//        if (!weightSp.empty()) {
//          EXPECT_EQ((ssize_t )nChan, weightSp.shape()[1]);
//          EXPECT_EQ((ssize_t )nRowSubchunk, weightSp.shape()[2]);
//        }
//
//        // polarization averaging specific check
//        // polarization id always points to the row to be appended
//        ASSERT_EQ(nRowPolarizationTable, (uInt )vb->polarizationId());
//        Validator::ValidatePolarization(corrTypes);
//
//        // validation of polarization average
//        if (!visCube.empty()) {
//          cout << "validate DATA" << endl;
//          Validator::ValidateData(visCube, ms, rowIds);
//        }
//        if (!visCubeCorrected.empty()) {
//          cout << "validate CORRECTED_DATA" << endl;
//          Validator::ValidateCorrected(visCubeCorrected, ms, rowIds);
//        }
//        if (!visCubeModel.empty()) {
//          cout << "validate MODEL_DATA" << endl;
//          Validator::ValidateModel(visCubeModel, ms, rowIds);
//        }
//        if (!visCubeFloat.empty()) {
//          cout << "validate FLOAT_DATA" << endl;
//          Validator::ValidateFloat(visCubeFloat, ms, rowIds);
//        }
//        Validator::ValidateFlag(flagCube, ms, rowIds);
//        Validator::ValidateFlagRow(flagRow, ms, rowIds);
//        Validator::ValidateWeight(weight, ms, rowIds);
//        Validator::ValidateWeightSp(weightSp, ms, rowIds);
//
        // next round of iteration
        refvi->next();
        vi->next();
      }
//
//      // chunk-subchunk consistency check
//      EXPECT_EQ(nRowChunk, nRowSubchunkSum);
//
      refvi->nextChunk();
      vi->nextChunk();
    }

    // make sure there is no chunk remaining
    EXPECT_FALSE(vi->moreChunks());
//
//    // chunk-ms consistency check
//    EXPECT_EQ(nRowMs, nRowChunkSum);
//
//    // iteration check
//    EXPECT_TRUE(allEQ(swept, (uInt )1));

  }

  template<class Manufacturer>
  void TestFactory(Int const &type_enum, String const &expectedClassName) {

    cout << "Type \"" << type_enum << "\" expected class name \""
        << expectedClassName << "\"" << endl;

    if (expectedClassName.size() > 0) {
      std::unique_ptr<VisibilityIterator2> vi(Manufacturer::ManufactureVI(ms_, type_enum));

      // Verify type string
      String viiType = vi->ViiType();
      EXPECT_TRUE(viiType.startsWith(expectedClassName));
    } else {
      cout << "Creation of VI via factory will fail" << endl;
      // exception must be thrown
      EXPECT_THROW( {
            std::unique_ptr<VisibilityIterator2> vi(Manufacturer::ManufactureVI(ms_, type_enum)); //new VisibilityIterator2(factory));
          },
          AipsError)<< "The process must throw AipsError";
    }
  }

private:
  void copyRegular(String const &src, String const &dst) {
    RegularFile r(src);
    r.copy(dst);
  }
  void copySymLink(String const &src, String const &dst) {
    Path p = SymLink(src).followSymLink();
    String actual_src = p.absoluteName();
    File f(actual_src);
    if (f.isRegular()) {
      copyRegular(actual_src, dst);
    } else if (f.isDirectory()) {
      copyDirectory(actual_src, dst);
    }
  }
  void copyDirectory(String const &src, String const &dst) {
    Directory dsrc(src);
    Directory ddst(dst);
    ddst.create();
    DirectoryIterator iter(dsrc);
    while (!iter.pastEnd()) {
      String name = iter.name();
      if (name.contains(".svn")) {
        iter++;
        continue;
      }
      File f = iter.file();
      Path psrc(src);
      Path pdst(dst);
      psrc.append(name);
      String sub_src = psrc.absoluteName();
      pdst.append(name);
      String sub_dst = pdst.absoluteName();
      if (f.isSymLink()) {
        copySymLink(sub_src, sub_dst);
      } else if (f.isRegular()) {
        copyRegular(sub_src, sub_dst);
      } else if (f.isDirectory()) {
        copyDirectory(sub_src, sub_dst);
      }
      iter++;
    }
  }
  void copyDataFromRepository(std::string const &data_dir) {
    if (my_data_name_.size() > 0) {
      std::string full_path = data_dir + my_data_name_;
      std::string work_path = my_data_name_;
      File f(full_path);
      ASSERT_TRUE(f.exists());
      if (f.isSymLink()) {
        copySymLink(full_path, work_path);
      } else if (f.isRegular()) {
        copyRegular(full_path, work_path);
      } else if (f.isDirectory()) {
        copyDirectory(full_path, work_path);
      }
    }
  }
  void cleanup() {
    if (my_data_name_.size() > 0) {
      deleteTable(my_data_name_);
    }
    deleteTable(my_ms_name_);
  }
  void deleteTable(std::string const &name) {
    File file(name);
    if (file.exists()) {
      std::cout << "Removing " << name << std::endl;
      Table::deleteTable(name, true);
    }
  }
};

class FiltrationTVITest: public FiltrationTVITestBase {
protected:
  virtual std::string GetDataName() {
    return "analytic_type1.bl.ms";
  }

  virtual std::string GetRelativeDataPath() {
    return "tsdbaseline";
  }
};

TEST_F(FiltrationTVITest, PorousTest) {
  TestTVI<PorousValidator, TestManufacturer>();
}

TEST_F(FiltrationTVITest, NonporousTest) {
  TestTVI<NonporousValidator, TestManufacturer>();
}

TEST_F(FiltrationTVITest, FactoryTest) {
  Int const type_enum = (Int)FilteringType::SDDoubleCircleFilter;
  String const expected_name("FiltrationTVI<SDDoubleCircle>");
  TestFactory<BasicManufacturer>(type_enum, expected_name);
  TestFactory<LayerManufacturer>(type_enum, expected_name);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "FiltrationTVI test " << std::endl;
  return RUN_ALL_TESTS();
}
