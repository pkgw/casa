//# tViiLayerFactory.cc: Tests Recursive factory for layered VI2s
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


#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casacore/casa/OS/EnvVar.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casa/iostream.h>
#include <msvis/MSVis/IteratingParameters.h>
#include <msvis/MSVis/ViiLayerFactory.h>
#include <msvis/MSVis/LayeredVi2Factory.h>
#include <msvis/MSVis/TransformingVi2.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/test/TestUtilsTVI.h>
#include <casa/iomanip.h>
#include <gtest/gtest.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;
using namespace casa::vi::test;


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST( ViiLayerFactoryTest , ViiLayerFactoryBasicTest ) {
 
  // A very rudimentary test of a single layer

  SimpleSimVi2Parameters s0;
  SimpleSimVi2LayerFactory fac(s0);

  Vector<ViiLayerFactory*> facts(1);
  facts[0]=&fac;

  std::unique_ptr<VisibilityIterator2> vi(new VisibilityIterator2(facts));
  VisBuffer2 *vb = vi->getImpl()->getVisBuffer();

  // Has a viable VI2 been generated?
  Int niter(0);
  for (vi->originChunks();vi->moreChunks();vi->nextChunk()) {
    for (vi->origin();vi->more();vi->next()) {
      ASSERT_EQ(4,vb->nAntennas());
      ++niter;
    }
  }
  ASSERT_EQ(1,niter);
}


TEST( ViiLayerFactoryTest , ViiLayerFactoryRealDataBasicTest ) {
 
  // A very rudimentary test of a single layer, using a real MS


  String *casapath = new String[2];
  split(EnvironmentVariable::get("CASAPATH"), casapath, 2, String(" "));
  // Use of Path().absoluteName() absorbs relative stuff in casapath
  String mspath(Path(casapath[0]+"/data/regression/unittest/flagdata/Four_ants_3C286.ms").absoluteName());

  MeasurementSet ms(mspath);

  Double interval(60000.0);
  IteratingParameters ipar(interval);   // seems to include SCAN_NUMBER automatically?
  VisIterImpl2LayerFactory fac(&ms,ipar,False); 

  Vector<ViiLayerFactory*> facts(1);
  facts[0]=&fac;

  std::unique_ptr<VisibilityIterator2> vi(new VisibilityIterator2(facts));
  VisBuffer2 *vb = vi->getImpl()->getVisBuffer();

  // Has a viable VI2 been generated?
  Int chunk(0),niter(0);
  for (vi->originChunks();vi->moreChunks();vi->nextChunk(),++chunk) {
    vi->origin();
    /*
    cout << "ch="<< chunk
	 << " scan="<<vb->scan()(0)
	 << " field="<< vb->fieldId()(0)
	 << " spw="<< vb->spectralWindows()(0)
	 << endl;
    */

    for (vi->origin();vi->more();vi->next(),++niter) {

      /*
      cout << "*************************************" << endl;
      cout << "chunk="<< chunk << ";  niter=" << niter << endl;

      cout << " scan=" << vb->scan()(0) << endl;
      cout << " fieldId=" << vb->fieldId()(0) << endl;

      cout << "nAntennas=" << vb->nAntennas() << endl;
      cout << "nRows=" << vb->nRows() << endl;
      cout << "nChannels=" << vb->nChannels() << endl;
      cout << "nCorrelations=" << vb->nCorrelations() << endl;
      cout << "spectralWindows=" << vb->spectralWindows() << endl;
      */

      ASSERT_EQ(4,vb->nAntennas());
      ASSERT_EQ(6,vb->nRows());
      ASSERT_EQ(64,vb->nChannels());  // all spws
      ASSERT_EQ(4,vb->nCorrelations());

    }
  }

  //cout << "chunk=" << chunk << endl;
  //  cout << "niter=" << niter << endl;

  ASSERT_EQ(32,chunk);
  ASSERT_EQ(2864,niter);
  delete[] casapath;
}

/* 
 * This class is a very simple TVI that passes all the requests for the
 * underlying VI (delegating to TransformingVi2) except that it swaps
 * the access to corrected data and data, i.e., when accessing the 
 * data corrected it returns the data column of the underlying VI.
 */
class DataSwappingTVI : public TransformingVi2
{
public:
  
  //Constructor
  DataSwappingTVI(ViImplementation2 * inputVii) :
    TransformingVi2 (inputVii)
  {
    setVisBuffer(createAttachedVisBuffer (VbRekeyable));
  }

  //Access to DATA CORRECTED returns underlying DATA column
  virtual void visibilityCorrected(casacore::Cube<casacore::Complex>& vis) const
  {
    VisBuffer2* vb = getVii()->getVisBuffer();
    vis = vb->visCube();
  }
  
  //Access to DATA returns underlying DATA CORRECTED column
  virtual void visibilityObserved(casacore::Cube<casacore::Complex>& vis) const
  {
    VisBuffer2* vb = getVii()->getVisBuffer();
    vis = vb->visCubeCorrected();
  }
  
  void origin()
  {
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

  void next()
  {
    // Drive underlying ViImplementation2
    getVii()->next();

    // Synchronize own VisBuffer
    configureNewSubchunk();
  }

};

/*
 * Factory that allows the creation of DataSwappingTVI classes.
 * This factory doesn't have any parameter to configure
 */
class DataSwappingTVILayerFactory : public ViiLayerFactory
{

public:

  DataSwappingTVILayerFactory()
  {
  }

  virtual ~DataSwappingTVILayerFactory() {};

protected:

  virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const
  {
    ViImplementation2 *vii = new DataSwappingTVI(vii0);
    return vii;
  }
};

/*
 * Gtest fixture used to test the access to the different data columns.
 * This class will create a synthetic MS in a temporary directory with or 
 * without DATA CORRECTED columns. 
 * It will also create a stack of TVIs with two TVIs: 
 * a disk access layer and a swapping data TVI (DataSwappingTVI). This TVI
 * is optional.
 */
class DataAccessTest : public MsFactoryTVITester
{
public:
  
  /*
   * Constructor: create the temporary dir and the MsFactory used later on
   * to create the MS.
   */
  DataAccessTest() :
    MsFactoryTVITester("test_tViiLayerFactory","DataAccessTest")
  {
  }

  /*
   * Do not create DATA CORRECTED column. 
   * This function must be called before createTVIs()
   */
  void removeCorrectedData()
  {
    msf_p->removeColumn(MS::CORRECTED_DATA);
  }

  /*
   * Do not create DATA column. 
   * This function must be called before createTVIs()
   */
  void removeData()
  {
    msf_p->removeColumn(MS::DATA);
  }

  /*
   * Switch the creation of the DataSwappingTVI on top the disk access layer.
   * This function must be called before createTVIs(). 
   */
  void addSwappingDataTVI()
  {
    withSwappingDataTVI_p = true;
  }
  
  /*
   * Create the synthetic MS and the TVI stack to access it. 
   */
  void createTVIs()
  {
    createMS();
    
    //Create a disk layer type VI Factory
    IteratingParameters ipar;
    VisIterImpl2LayerFactory diskItFac(ms_p.get(),ipar,false);

    //Create a SwappingDataTVI Factory if requested
    std::unique_ptr<DataSwappingTVILayerFactory> swapFac;
    if(withSwappingDataTVI_p)
      swapFac.reset(new DataSwappingTVILayerFactory());

    //Create a layered factory with all the layers of factories
    size_t nFac = 1;
    if(withSwappingDataTVI_p) 
      nFac++;
    std::vector<ViiLayerFactory*> factories(nFac);
    factories[0]=&diskItFac;
    if(withSwappingDataTVI_p)
      factories[1]= swapFac.get();

    instantiateVI(factories);
  }

  //Destructor
  ~DataAccessTest()
  {
  }

  //Wether to use the DataSwappingTVI
  bool withSwappingDataTVI_p = false;
};
 

/*
 * This test will simply access the corrected data column of a
 * synthetic created MS
 */
TEST_F(DataAccessTest, AccessCorrectedData)
{
  createTVIs();

  //Traverse the iterator accessing the corrected data cube
  visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();});
}

/*
 * This test will check that an exception is thrown if the
 * CORRECTED DATA column is missing in the MS.
 */
TEST_F(DataAccessTest, AccessCorrectedDataWhenMissing)
{
  removeCorrectedData();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This should
  //throw, since it has been removed from the MS.
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();}),
               AipsError);
}

/*
 * This test will access the corrected data column of a
 * synthetic MS that doesn't contain that column. The upper TVI, however
 * will swap the access to DATA column instead so the test should succeed
 */
TEST_F(DataAccessTest, AccessCorrectedDataInSwappingDataTVI)
{
  removeCorrectedData();

  addSwappingDataTVI();
  
  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This works
  //despite removing the corrected data from disk because the upper TVI layer
  //(DataSwappingTVI) delivers DATA when CORRECTED DATA is requested.
  visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();});
}

/*
 * This test will access the corrected data column of a
 * synthetic that doesn't contain column DATA. The upper TVI, however,
 * will deliver the underlying DATA column when asking for CORRECTED DATA,
 * and therefore the test will fail.
 */
TEST_F(DataAccessTest, AccessCorrectedDataInSwappingDataTVIWhenMissingData)
{
  removeData();

  addSwappingDataTVI();
  
  createTVIs();

  //Traverse the iterator accessing the corrected data cube. 
  //The swapping TVI will access the underlying DATA column when accessing
  //the CORRECTED DATA. Since DATA has been removed from the MS this should
  //throw
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();}),
               AipsError);
}

/* 
 * This test will simply access the data column of a
 * synthetic created MS
 */
TEST_F(DataAccessTest, AccessData)
{
  createTVIs();

  //Traverse the iterator accessing the corrected data cube
  visitIterator([&]() -> void {vb_p->visCube().shape();});
}

/*
 * This test will check that an exception is thrown if the
 * CORRECTED DATA column is missing in the MS.
 */
TEST_F(DataAccessTest, AccessDataWhenMissing)
{
  removeData();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This should
  //throw, since it has been removed from the MS.
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCube().shape();}),
               AipsError);
}

/*
 * This test will access the corrected data column of a
 * synthetic that doesn't contain that column. The upper TVI, however
 * will swap the access to DATA column instead so the test should succeed
 */
TEST_F(DataAccessTest, AccessDataInSwappingDataTVI)
{
  removeData();

  addSwappingDataTVI();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This works
  //despite removing the corrected data from disk because the upper TVI layer
  //(DataSwappingTVI) delivers DATA when CORRECTED DATA is requested.
  visitIterator([&]() -> void {vb_p->visCube().shape();});
}

TEST_F(DataAccessTest, AccessDataInSwappingDataTVIWhenMissingCorrectedData)
{
  removeCorrectedData();

  addSwappingDataTVI();

  createTVIs();

  //Traverse the iterator accessing the corrected data cube.
  //The swapping TVI will access the underlying DATA column when accessing
  //the CORRECTED DATA. Since DATA has been removed from the MS this should
  //throw
  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCube().shape();}),
               AipsError);
}

/*
 * This class is a simplistic TVI that modifies the subtables
 * antenna, spw and dd
 */
class SubtableChangerTVI : public TransformingVi2
{
public:

    //Constructor
    SubtableChangerTVI(ViImplementation2 * inputVii) :
        TransformingVi2 (inputVii)
    {
        setVisBuffer(createAttachedVisBuffer (VbRekeyable));
        resetSubtables();
    }

    void origin()
    {
        // Drive underlying ViImplementation2
        getVii()->origin();

        // Synchronize own VisBuffer
        configureNewSubchunk();
    }

    void next()
    {
        // Drive underlying ViImplementation2
        getVii()->next();

        // Synchronize own VisBuffer
        configureNewSubchunk();
    }

    void
    originChunks(Bool forceRewind) override
    {
        // Drive underlying ViImplementation2
        getVii()->originChunks(forceRewind);

        // Potentially the new chunk can be from a different MS
        resetSubtables();
    }

    void
    nextChunk() override
    {
        // Drive underlying ViImplementation2
        getVii()->nextChunk();

        // Potentially the new chunk can be from a different MS
        resetSubtables();
    }

    void resetSubtables()
    {
        // Note that the creation of a new subtables is done using a
        // copy of the original subtables. However, to access these we
        // need to use the method table() of a given column (e. g. name() )
        // It would be better if the ROMSAntennaColumns object had
        // an getter to the MSAntenna object.
        // The same applies to the other subtables

        // Create antenna subtable
        auto& underlyingAntennaSubtablecols = getVii()->antennaSubtablecols();
        auto underlyingAntennaSubtable = underlyingAntennaSubtablecols.name().table();
        newAntennaSubtable_p = underlyingAntennaSubtable.copyToMemoryTable("SubtableChangerAntennaSubtable");
        newAntennaSubtablecols_p.reset(new MSAntennaColumns(newAntennaSubtable_p));
        // Add one antenna
        newAntennaSubtable_p.addRow();

        // Create DD subtable
        auto& underlyingDDSubtablecols = getVii()->dataDescriptionSubtablecols();
        auto underlyingDDSubtable = underlyingDDSubtablecols.spectralWindowId().table();
        newDDSubtable_p = underlyingDDSubtable.copyToMemoryTable("SubtableChangerDDSubtable");
        newDDSubtablecols_p.reset(new MSDataDescColumns(newDDSubtable_p));
        // Double the rows
        auto nrowDD = newDDSubtable_p.nrow();
        for(size_t irow = 0 ; irow < nrowDD; irow++)
            newDDSubtable_p.addRow();

        // Create spw subtable
        auto& underlyingSPWSubtablecols = getVii()->spectralWindowSubtablecols();
        auto underlyingSPWSubtable = underlyingSPWSubtablecols.name().table();
        newSPWSubtable_p = underlyingSPWSubtable.copyToMemoryTable("SubtableChangerSPWSubtable");
        newSPWSubtablecols_p.reset(new MSSpWindowColumns(newSPWSubtable_p));
        // Double the rows
        auto nrowSPW = newSPWSubtable_p.nrow();
        for(size_t irow = 0 ; irow < nrowSPW; irow++)
            newSPWSubtable_p.addRow();
    }


    const casacore::ROMSAntennaColumns& antennaSubtablecols() const override
    {
        return *newAntennaSubtablecols_p;
    }

    const casacore::ROMSDataDescColumns& dataDescriptionSubtablecols() const override
    {
        return *newDDSubtablecols_p;
    }

    const casacore::ROMSSpWindowColumns& spectralWindowSubtablecols() const override
    {
        return *newSPWSubtablecols_p;
    }

private:

    casacore::MSAntenna newAntennaSubtable_p;
    std::unique_ptr<casacore::ROMSAntennaColumns> newAntennaSubtablecols_p;

    casacore::MSSpectralWindow newSPWSubtable_p;
    std::unique_ptr<casacore::ROMSSpWindowColumns> newSPWSubtablecols_p;

    casacore::MSDataDescription newDDSubtable_p;
    std::unique_ptr<casacore::ROMSDataDescColumns> newDDSubtablecols_p;

};

/*
 * Factory that allows the creation of SubtableChangerTVI classes.
 * This factory doesn't have any parameter to configure
 */
class SubtableChangerTVILayerFactory : public ViiLayerFactory
{

public:

    SubtableChangerTVILayerFactory()
  {
  }

  virtual ~SubtableChangerTVILayerFactory() {};

protected:

  virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const
  {
    ViImplementation2 *vii = new SubtableChangerTVI(vii0);
    return vii;
  }
};

/*
 * Gtest fixture used to test the creation of subtables by a TVI.
 * This class will create a synthetic MS in a temporary directory.
 * It will also create a stack of TVIs with two TVIs:
 * a disk access layer and a TVI that modifies subtables (SubtableChangerTVI).
 * The function checkSubtables can actually check that the subtables
 * have been changed as expected.
 */
class SubtableChangerTest : public MsFactoryTVITester
{
public:

    /*
     * Constructor: create the temporary dir and the MsFactory used later on
     * to create the MS.
     */
    SubtableChangerTest() : 
        MsFactoryTVITester("tViiLayerFactory","SubtableChangerTest")
    {
    }

    /*
     * Create the synthetic MS and the TVI stack to access it.
     */
    void createTVIs()
    {
        // Set the number of antennas and SPWs for the MS generated on disk
        nAntennas = 10;
        nSPWs = 10;

        msf_p->addAntennas(nAntennas);
        msf_p->addSpectralWindows(nSPWs);

        // Create synthethic MS using the msf_p factory
        createMS();

        // Create a disk layer type VI Factory
        IteratingParameters ipar;
        VisIterImpl2LayerFactory diskItFac(ms_p.get(),ipar,false);

        // Create a SubtableChangerTVI Factory
        SubtableChangerTVILayerFactory subtableChangerFac;

        // Create a layered factory with all the layers of factories
        size_t nFac = 2;
        std::vector<ViiLayerFactory*> facts(nFac);
        facts[0]=&diskItFac;
        facts[1]= &subtableChangerFac;

        // Finally create the top VI
        instantiateVI(facts);
    }

    void checkSubtables()
    {
        // Check the antenna tables size (the TVI has added one antenna)
        EXPECT_EQ(nAntennas + 1, vi_p->antennaSubtablecols().nrow());

        // Check the SPW tables size (the TVI has doubled the number)
        EXPECT_EQ(nSPWs * 2, vi_p->spectralWindowSubtablecols().nrow());

        // Check the DD tables size (the TVI has doubled the number)
        EXPECT_EQ(nSPWs * 2, vi_p->dataDescriptionSubtablecols().nrow());
    }

    // The number of antennas originally created
    size_t nAntennas;
    // The number of SPWs originally created
    size_t nSPWs;
};


/*
 * This test will simply access the corrected data column of a
 * synthetic created MS
 */
TEST_F(SubtableChangerTest, CheckSubtables)
{

  createTVIs();

  checkSubtables();
}

