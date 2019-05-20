//# tHanningSmoothTVI: This file contains the unit tests of the HanningSmoothTVI class.
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


#include <mstransform/TVI/test/tHanningSmoothTVI.h>
#include <msvis/MSVis/PassThroughTVI.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;


//////////////////////////////////////////////////////////////////////////
// HanningSmoothTVITest class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::generateTestFile()
{
    String path("");
    if (autoMode_p) path = String("/data/regression/unittest/flagdata/");
    copyTestFile(path,inpFile_p,testFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::generateReferenceFile()
{
    String path("");
    if (autoMode_p) path = String("/data/regression/unittest/flagdata/");
    copyTestFile(path,inpFile_p,referenceFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::initTestConfiguration(Record &configuration)
{
    testConfiguration_p = configuration;
    testConfiguration_p.define ("inputms", testFile_p);

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::initReferenceConfiguration(Record &configuration)
{
    refConfiguration_p = configuration;
    refConfiguration_p.define ("inputms", referenceFile_p);
    refConfiguration_p.define ("reindex", false);
    refConfiguration_p.define ("hanning", true);
    refConfiguration_p.define ("datacolumn", String("ALL"));

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
HanningSmoothTVITest::HanningSmoothTVITest(): FreqAxisTVITest ()
{
    inpFile_p = String("Four_ants_3C286.ms");
    testFile_p = String("Four_ants_3C286.ms.test");
    referenceFile_p = String("Four_ants_3C286.ms.ref");

    Record configuration;
    configuration.define ("spw", "1");

    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
HanningSmoothTVITest::HanningSmoothTVITest(Record configuration): FreqAxisTVITest(configuration)
{
    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::TestBody()
{
    SetUp();
    testCompareTransformedData();
    TearDown();

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void HanningSmoothTVITest::testCompareTransformedData()
{
    // Declare working variables
    Float tolerance = 1E-5; // FLT_EPSILON is 1.19209290e-7F

    // Create MSTransformIterator pointing to reference file
    refConfiguration_p.define("factory",False);
    MSTransformIteratorFactory refFactory(refConfiguration_p);
    VisibilityIterator2 refTVI(refFactory);

    // Use MSTransformFactory to create a plain input VII
    testConfiguration_p.define("factory",True);
    MSTransformIteratorFactory plainVIFactory(testConfiguration_p);
    ViImplementation2 *inputVI = plainVIFactory.getInputVI()->getImpl();

    // Generate TVI to test
    HanningSmoothTVIFactory testFactory(inputVI);
    VisibilityIterator2 testTVI(testFactory);

    // Determine columns to check
    VisBufferComponents2 columns;
    columns += VisBufferComponent2::NRows;
    columns += VisBufferComponent2::NChannels;
    columns += VisBufferComponent2::NCorrelations;
    columns += VisBufferComponent2::FlagRow;
    columns += VisBufferComponent2::FlagCube;
    columns += VisBufferComponent2::VisibilityCubeObserved;
    columns += VisBufferComponent2::VisibilityCubeCorrected;
    columns += VisBufferComponent2::VisibilityCubeModel;
    columns += VisBufferComponent2::WeightSpectrum;
    columns += VisBufferComponent2::SigmaSpectrum;
    columns += VisBufferComponent2::Weight;
    columns += VisBufferComponent2::Sigma;

    // Compare
    SCOPED_TRACE("Comparing transformed data");
    compareVisibilityIterators(testTVI,refTVI,columns,tolerance);
}

//////////////////////////////////////////////////////////////////////////
// Googletest macros
//////////////////////////////////////////////////////////////////////////
TEST_F(HanningSmoothTVITest, testCompareTransformedData)
{
    testCompareTransformedData();
}

HanningSmoothTVISpwChannTest::HanningSmoothTVISpwChannTest() :
    MsFactoryTVITester("HanningSmoothTVIMSSelTest","MSFactoryCreated"),
    useMSSelection_p(false), addPassThroughTVI_p(false), addExtraHanningTVI_p(false)
{
}

void HanningSmoothTVISpwChannTest::useMSSelection(bool use)
{
    useMSSelection_p = use;
}

void HanningSmoothTVISpwChannTest::addPassThroughTVI(bool use)
{
    addPassThroughTVI_p = use;
}

void HanningSmoothTVISpwChannTest::addExtraHanningTVI(bool use)
{
    addExtraHanningTVI_p = use;
}

void HanningSmoothTVISpwChannTest::createTVIs()
{
    //Setting the parameters to generate a synthetic MS 
    double timeInterval = 10;
    double timeSpan = 100;
    msf_p->setTimeInfo (0, timeSpan, timeInterval);
    msf_p->addAntennas(6);
    //Adding two spectral windows with different number of channels
    int nChannels = 100;
    double frequency = 1e11;
    double frequencyDelta = 1e9;
    std::string stokes("XX YY");
    msf_p->addSpectralWindow("SPW0", nChannels, 
                             frequency, frequencyDelta, stokes);
    nChannels = 50;
    frequency = 2e11;
    msf_p->addSpectralWindow("SPW1", nChannels, 
                             frequency, frequencyDelta, stokes);
    msf_p->addFeeds (10); //Need antenna and spw to be set up first

    //Creates the synthetic MS 
    createMS();

    //If there is selection, create a new MS as a selection on the original one
    std::unique_ptr<casacore::MeasurementSet> msSelected;
    std::shared_ptr<vi::FrequencySelectionUsingChannels> freqSel;
    if(useMSSelection_p)
    {
        MSSelection thisSelection;
        std::string spwSelection("0:1~40");
        thisSelection.setSpwExpr(spwSelection);
        TableExprNode exprNode = thisSelection.toTableExprNode(ms_p.get());
        msSelected.reset(new MeasurementSet((*ms_p.get())(exprNode)));
        freqSel = std::make_shared<vi::FrequencySelectionUsingChannels>();
        freqSel->add(thisSelection, ms_p.get());
    }

    //Create a VI Factory that access directly the MS (or the selected MS)
    IteratingParameters ipar;
    std::unique_ptr<VisIterImpl2LayerFactory> diskItFac;
    if(useMSSelection_p)
    {
        diskItFac.reset(new VisIterImpl2LayerFactory(msSelected.get(),ipar, false));
        auto selections = std::make_shared<vi::FrequencySelections>();
        selections->add(*freqSel);
        diskItFac->setFrequencySelections(selections);
    }
    else
        diskItFac.reset(new VisIterImpl2LayerFactory(ms_p.get(),ipar, false));

    //Create a HanningSmoothTVI Factory
    std::unique_ptr<HanningSmoothTVILayerFactory> hannningFac;
    hannningFac.reset(new HanningSmoothTVILayerFactory());

    //Create a layered factory with all the layers of factories
    //If requested a PassThroughTVI is added to test that HanningSmooth
    //forwards properly all the information to layers above
    std::vector<ViiLayerFactory*> factories;
    factories.push_back(diskItFac.get());
    factories.push_back(hannningFac.get());
    std::unique_ptr<PassThroughTVILayerFactory> passThroughFactory;
    if(addPassThroughTVI_p)
    {
        passThroughFactory.reset(new PassThroughTVILayerFactory());
        factories.push_back(passThroughFactory.get());
    }
    //Add a second averaging layer if requested
    std::unique_ptr<HanningSmoothTVILayerFactory> hanningFac2;
    if(addExtraHanningTVI_p)
    {
      hanningFac2.reset(new HanningSmoothTVILayerFactory());
        factories.push_back(hanningFac2.get());
    }
    
    //Finally create the VI resulting from all the layered TVIs
    instantiateVI(factories);
}

TEST_F(HanningSmoothTVISpwChannTest, CheckOutputSpwChannels)
{
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after smoothing the first spectral window has 100
    //channels and the second spectral window has 50 channels (the original ones)
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 100);
                                 }
                                 else if(allEQ(vb_p->spectralWindows(), 1)) //SPW1
                                 {
                                     nRowsSpw1+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 50);
                                 }
                                 ASSERT_EQ(vi_p->nSpectralWindows(), 2);
                                 nRows+=vb_p->visCube().shape()[2];});

    //All the original rows
    size_t expectedRows = 300;
    //The synthetic MS has half of the rows in each SPW
    size_t expectedRowsSpw0 = expectedRows / 2; 
    size_t expectedRowsSpw1 = expectedRows / 2;

    //Check that we get the number of expected rows for each spw
    ASSERT_EQ(nRows, expectedRows);
    ASSERT_EQ(nRowsSpw0, expectedRowsSpw0);
    ASSERT_EQ(nRowsSpw1, expectedRowsSpw1);
}

TEST_F(HanningSmoothTVISpwChannTest, CheckMSSelOutputSpwChannels)
{
    useMSSelection(true);
    addPassThroughTVI(true);
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after smoothing and spw selection "0:1~40",
    //the first spectral window has 40 channels 
    //and there are no rows with the second spectral window 
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
    
    std::cout<<vb_p->spectralWindows()<<vb_p->nChannels()<<std::endl;
    std::cout<<vb_p->antenna1()<<std::endl;
    std::cout<<vb_p->dataDescriptionIds()<<std::endl;
    std::cout<<shape<<std::endl;
    
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 40);
                                 }
                                 //Note that as per CAS-10294 this should be 3
                                 ASSERT_EQ(vi_p->nSpectralWindows(), 2);
                                 nRows+=vb_p->visCube().shape()[2];});

    //After selection we have only half of the original rows
    //and all of them belong to SPW0.
    size_t expectedRows = 150;
    size_t expectedRowsSpw0 = expectedRows;
    size_t expectedRowsSpw1 = 0;

    //Check that we get the number of expected rows for each spw
    ASSERT_EQ(nRows, expectedRows);
    ASSERT_EQ(nRowsSpw0, expectedRowsSpw0);
    ASSERT_EQ(nRowsSpw1, expectedRowsSpw1);
}

TEST_F(HanningSmoothTVISpwChannTest, CheckMSSelOutputTwoAvgSpwChannels)
{
    //useMSSelection(true);
    addPassThroughTVI(true);
    addExtraHanningTVI(true);
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after smoothing twice the first spectral window has 100
    //channels and the second spectral window has 50 channels (the original ones)
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 100);
                                 }
                                 else if(allEQ(vb_p->spectralWindows(), 1)) //SPW1
                                 {
                                     nRowsSpw1+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 50);
                                 }
                                 //Note that as per CAS-10294 this should be 3
                                 ASSERT_EQ(vi_p->nSpectralWindows(), 2);
                                 nRows+=vb_p->visCube().shape()[2];});

    //All the original rows
    size_t expectedRows = 300;
    //The synthetic MS has half of the rows in each SPW
    size_t expectedRowsSpw0 = expectedRows / 2;
    size_t expectedRowsSpw1 = expectedRows / 2;

    //Check that we get the number of expected rows for each spw
    ASSERT_EQ(nRows, expectedRows);
    ASSERT_EQ(nRowsSpw0, expectedRowsSpw0);
    ASSERT_EQ(nRowsSpw1, expectedRowsSpw1);
}

//////////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    int ret;
    string parameter,value;
    Record configuration;
    Bool autoMode = true;

    for (unsigned short i=0;i<argc-1;i++)
    {
        parameter = string(argv[i]);
        value = string(argv[i+1]);

        if (parameter == string("-vis"))
        {
            configuration.define ("inputms", value);
            autoMode = false;
        }
        else if (parameter == string("-spw"))
        {
            configuration.define ("spw", value);
        }
    }

    if (autoMode)
    {
        ::testing::InitGoogleTest(&argc, argv);
        ret = RUN_ALL_TESTS();
    }
    else
    {
        HanningSmoothTVITest test(configuration);
        test.TestBody();
        if (test.getTestResult()) ret = 0;
        else ret = 1;
    }

    return ret;
}
