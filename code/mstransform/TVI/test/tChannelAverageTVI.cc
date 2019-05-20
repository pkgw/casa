//# tChannelAverageTVI: This file contains the unit tests of the ChannelAverageTVI class.
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


#include <limits>
#include <mstransform/TVI/test/tChannelAverageTVI.h>
#include <msvis/MSVis/SimpleSimVi2.h>
#include <msvis/MSVis/PassThroughTVI.h>
#include <msvis/MSVis/test/MsFactory.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;


//////////////////////////////////////////////////////////////////////////
// ChannelAverageTVITest class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::generateTestFile()
{
	String path("");
	if (autoMode_p) path = String("/data/regression/unittest/flagdata/");
	copyTestFile(path,inpFile_p,testFile_p);

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::generateReferenceFile()
{
	String path("");
	if (autoMode_p) path = String("/data/regression/unittest/flagdata/");
	copyTestFile(path,inpFile_p,referenceFile_p);

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::initTestConfiguration(Record &configuration)
{
	testConfiguration_p = configuration;
	testConfiguration_p.define ("inputms", testFile_p);

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::initReferenceConfiguration(Record &configuration)
{
	refConfiguration_p = configuration;
	refConfiguration_p.define ("inputms", referenceFile_p);
	refConfiguration_p.define ("reindex", false);
	refConfiguration_p.define ("chanaverage", true);
	refConfiguration_p.define ("datacolumn", String("ALL"));

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
ChannelAverageTVICompareTest::ChannelAverageTVICompareTest(): FreqAxisTVITest ()
{
	  inpFile_p = String("Four_ants_3C286.ms");
    testFile_p = String("Four_ants_3C286.ms.test");
    referenceFile_p = String("Four_ants_3C286.ms.ref");

    Record configuration;
    configuration.define ("spw", "1:8~63,4:16~63");
    configuration.define ("chanbin", 8);

    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
ChannelAverageTVICompareTest::ChannelAverageTVICompareTest(Record configuration): FreqAxisTVITest(configuration)
{
    init(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::TestBody()
{
	SetUp();
	testCompareMSTransformTransformedData();
	TearDown();

	SetUp();
	testCompareMSTransformPropagatedFlags();
	TearDown();

	return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::testCompareMSTransformTransformedData()
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
	ChannelAverageTVIFactory testFactory(testConfiguration_p,inputVI);
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
	columns += VisBufferComponent2::Frequencies;

	// Compare
    SCOPED_TRACE("Comparing transformed data");
	compareVisibilityIterators(testTVI,refTVI,columns,tolerance);

}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::testCompareMSTransformPropagatedFlags()
{
    
    // Declare working variables
    Float tolerance = 1E-5; // FLT_EPSILON is 1.19209290e-7F

    // Propagate flags
    propagateFlags();

    // Use MSTransformIteratorFactory to create a plain input VI pointing to the test file
    testConfiguration_p.define("factory",False);
    MSTransformIteratorFactory testFactory(testConfiguration_p);
    VisibilityIterator2 *testTVI = testFactory.getInputVI();

    // Use MSTransformIteratorFactory to create a plain input VI pointing to the reference file
    refConfiguration_p.define("factory",False);
    MSTransformIteratorFactory refFactory(refConfiguration_p);
    VisibilityIterator2 *refTVI = refFactory.getInputVI();

    // Determine columns to check
    VisBufferComponents2 columns;
    columns += VisBufferComponent2::FlagCube;

    // Compare
    SCOPED_TRACE("Comparing propagated flags");
    compareVisibilityIterators(*testTVI,*refTVI,columns,tolerance);

}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::testWriteFlags()
{
    //Tolerance for booleans is just 0
    Float tolerance = 0;

    // Create VisibilityIterator2 pointing to reference file
    MeasurementSet msRef(referenceFile_p, Table::Update);
    vi::VisibilityIterator2 refVI(msRef, SortColumns(), true);

    // Propagate flags in raw MS (takes into account the chanbin)
    flagEachOtherChannel(refVI, true, refConfiguration_p.asInt("chanbin"));

    //Create channel average TVI pointing to the testing file
    {
        MeasurementSet msTest(testFile_p, Table::Update);
        vi::VisibilityIterator2*  testVI =
                new vi::VisibilityIterator2(msTest, SortColumns(), true);
        Record configuration;
        configuration.define ("chanbin", 8);
        ChannelAverageTVIFactory testFactory(configuration,
                                             testVI->getImpl());
        VisibilityIterator2 testTVI(testFactory);

        // Propagate flags using the TVI
        flagEachOtherChannel(testTVI, false);

        //testVI is deleted by testTVI destructor (!)
    }

    // Create VisibilityIterator2 pointing to test file
    //(after flags have been written)
    MeasurementSet msTestAfter(testFile_p);
    vi::VisibilityIterator2 testVIAfter(msTestAfter);

    // Determine columns to check
    VisBufferComponents2 columns;
    columns += VisBufferComponent2::FlagCube;

    SCOPED_TRACE("Comparing propagated flags");
    compareVisibilityIterators(testVIAfter, refVI, columns, tolerance);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void ChannelAverageTVICompareTest::propagateFlags()
{
	// Create MSTransformIterator pointing to reference file
	refConfiguration_p.define("factory",False);
	MSTransformIteratorFactory refFactory(refConfiguration_p);
	VisibilityIterator2 refTVI(refFactory);

	// Use MSTransformFactory to create a plain input VII
	testConfiguration_p.define("factory",True);
	MSTransformIteratorFactory plainVIFactory(testConfiguration_p);
	ViImplementation2 *inputVI = plainVIFactory.getInputVI()->getImpl();

	// Generate TVI to test
	ChannelAverageTVIFactory testFactory(testConfiguration_p,inputVI);
	VisibilityIterator2 testTVI(testFactory);

	// Propagate flags with MSTransformIterator
	flagEachOtherChannel(refTVI, false);

	// Propagate flags with TVI to test
	flagEachOtherChannel(testTVI, false);

	return;
}

//////////////////////////////////////////////////////////////////////////
// Googletest macros
//////////////////////////////////////////////////////////////////////////
TEST_F(ChannelAverageTVICompareTest, CompareMSTransformTransformedData)
{
	testCompareMSTransformTransformedData();
}

TEST_F(ChannelAverageTVICompareTest, CompareMSTransformPropagatedFlags)
{
    testCompareMSTransformPropagatedFlags();
}

TEST_F(ChannelAverageTVICompareTest, TestWriteFlags)
{
    testWriteFlags();
}

TEST(ChannelAverageTVIConfTest, NoChanbinParam)
{
    Record configuration;
    ChannelAverageTVIFactory testFactory(configuration, nullptr);
    ASSERT_THROW(VisibilityIterator2 testTVI(testFactory), AipsError);
}

TEST(ChannelAverageTVIConfTest, WrongChanbinType)
{
    //Checks that an exception is thrown if chanbin parameter 
    //has an invalid type. Only Int and Array<Int> are allowed
    Record configuration;
    configuration.define ("chanbin", true); //Checking boolean
    ChannelAverageTVIFactory testFactory(configuration, nullptr);
    ASSERT_THROW(VisibilityIterator2 testTVI(testFactory), AipsError);
    Record configuration2;
    configuration2.define ("chanbin", 4.5); //Checking double
    ChannelAverageTVIFactory testFactory2(configuration2, nullptr);
    ASSERT_THROW(VisibilityIterator2 testTVI2(testFactory2), AipsError);
}

TEST(ChannelAverageTVIConfTest, WrongChanbinValue)
{
    //Checks that an exception is thrown if chanbin parameter 
    //has an invalid value (the string cannot be converted to Array<Int>)
    Record configuration;
    configuration.define ("chanbin", "invalid");
    ChannelAverageTVIFactory testFactory(configuration, nullptr);
    ASSERT_THROW(VisibilityIterator2 testTVI(testFactory), AipsError);
}

TEST(ChannelAverageTVIConfTest, NullInput)
{
    //Check that an exception is thrown if the input ViImplementation is null
    Record configuration;
    configuration.define ("chanbin", 2);
    ChannelAverageTVIFactory testFactory(configuration, nullptr);
    ASSERT_THROW(VisibilityIterator2 testTVI(testFactory), AipsError);
}

TEST(ChannelAverageTVIConfTest, WrongChanbinForMultipleSpw)
{
    Record configuration;
    configuration.define ("chanbin", "2,2,2");
    //Generates a simulated Vi with  nField = 1 , nScan = 1, nSpw = 2, nAnt = 4,
    //nCorr = 4, nTimePerField = (1), nChan = (10, 10)
    SimpleSimVi2Parameters simParam(1, 1, 2, 4, 4 , 
                                    Vector<Int>(1,1), Vector<Int>(2,10));
    SimpleSimVi2Factory simFactory(simParam);
    VisibilityIterator2 *simVi = new VisibilityIterator2(simFactory);
    ChannelAverageTVIFactory testFactory(configuration, simVi->getImpl());
    ASSERT_THROW(VisibilityIterator2 testTVI(testFactory), AipsError);
}

TEST(ChannelAverageTVIExecuteSimulatedTest, UniformMS)
{
    for(int nField = 1; nField < 4; nField+=2)
    {
        for(int nScan = 1; nScan < 3; nScan++)
        {
            for(int nSpw = 1; nSpw < 4; nSpw+=2)
            {
                for(int nAnt = 4; nAnt <= 8; nAnt*=2)
                {
                    for(int nTimePerField = 1; nTimePerField <= 10; nTimePerField*=10)
                    {
                        for(int nChan = 4; nChan <= 64; nChan*=2)
                        {
                            for(int chanbin = 0; chanbin <= nChan; 
                                    chanbin = (chanbin == 0 ? 1 : chanbin*4))
                            {
                                int nCorr = 4;
                                casacore::Complex visValue(1.0, 2.0);
                                SCOPED_TRACE(string("Channel averaging data with nField=") + 
                                             to_string(nField) + " nScan=" + to_string(nScan) +
                                             " nAnt " + to_string(nAnt) + " nCorr=" + to_string(nCorr) +
                                             " nTimePerField " + to_string(nTimePerField) + 
                                             " nChan=" + to_string(nChan));
                                //Generating a uniform simulated Vi2 
                                SimpleSimVi2Parameters simParam(nField, nScan, 
                                    nSpw, nAnt, nCorr,
                                    Vector<Int>(nField, nTimePerField),
                                    Vector<Int>(nSpw, nChan), visValue);
                                SimpleSimVi2Factory simFactory(simParam);
                                VisibilityIterator2 *simVi = 
                                        new VisibilityIterator2(simFactory);
                                
                                //Chaining a ChannelAverageTVI 
                                //after the simulated Vi2
                                Record configuration;
                                configuration.define ("chanbin", chanbin);
                                ChannelAverageTVIFactory testFactory(configuration, 
                                                                     simVi->getImpl());
                                VisibilityIterator2 testTVI(testFactory);
                                
                                //Generating a simulated Vi2 with the  
                                //expected result (which is also uniform)
                                int nAveragedChannels = 
                                    chanbin == 0 ? nChan : (chanbin == 1 ? 1 : nChan / chanbin);
                                //chanbin == 0 means no averaging
                                //chanbin == 1 means do full averaging across the whole SPW
                                SimpleSimVi2Parameters simResultParam(nField, nScan, 
                                    nSpw, nAnt, nCorr, 
                                    Vector<Int>(nField, nTimePerField),
                                    Vector<Int>(nSpw, nAveragedChannels),
                                    visValue);
                                SimpleSimVi2Factory simResultFactory(simResultParam);
                                VisibilityIterator2 *simResultVi = 
                                        new VisibilityIterator2(simResultFactory);

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

                                // Compare the channel average 
                                SCOPED_TRACE("Comparing transformed data for simulated uniform ms");
                                double tolerance = std::numeric_limits<double>::epsilon();
                                compareVisibilityIterators(testTVI,*simResultVi,
                                                           columns, tolerance);
                            }
                        }
                    }
                }
            }
        }
    }
}

ChannelAverageTVISpwChannTest::ChannelAverageTVISpwChannTest() :
    MsFactoryTVITester("ChannelAverageTVIMSSelTest","MSFactoryCreated"),
    useMSSelection_p(false), addPassThroughTVI_p(false), addExtraAvgTVI_p(false)
{
}

void ChannelAverageTVISpwChannTest::useMSSelection(bool use)
{
    useMSSelection_p = use;
}

void ChannelAverageTVISpwChannTest::addPassThroughTVI(bool use)
{
    addPassThroughTVI_p = use;
}

void ChannelAverageTVISpwChannTest::addExtraAvgTVI(bool use)
{
    addExtraAvgTVI_p = use;
}

void ChannelAverageTVISpwChannTest::createTVIs()
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

    //Create a ChannelAverageTVI Factory
    std::unique_ptr<ChannelAverageTVILayerFactory> chanAvgFac;
    casacore::Record configuration;
    configuration.define ("chanbin", 5);
    chanAvgFac.reset(new ChannelAverageTVILayerFactory(configuration));

    //Create a layered factory with all the layers of factories
    //If requested a PassThroughTVI is added to test that ChannelAverage
    //forwards properly all the information to layers above
    std::vector<ViiLayerFactory*> factories;
    factories.push_back(diskItFac.get());
    factories.push_back(chanAvgFac.get());
    std::unique_ptr<PassThroughTVILayerFactory> passThroughFactory;
    if(addPassThroughTVI_p)
    {
        passThroughFactory.reset(new PassThroughTVILayerFactory());
        factories.push_back(passThroughFactory.get());
    }
    //Add a second averaging layer if requested
    std::unique_ptr<ChannelAverageTVILayerFactory> chanAvgFac2;
    if(addExtraAvgTVI_p)
    {
        casacore::Record configuration2;
        Array<int> chanArr(IPosition(1,2));
        chanArr[0] = 2;
        chanArr[1] = 5;
        configuration2.define ("chanbin", chanArr);
        chanAvgFac2.reset(new ChannelAverageTVILayerFactory(configuration2));
        factories.push_back(chanAvgFac2.get());
    }
    
    //Finally create the VI resulting from all the layered TVIs
    instantiateVI(factories);
}

TEST_F(ChannelAverageTVISpwChannTest, CheckOutputSpwChannels)
{
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after averaging with chanbin=5 the first spectral window has 20
    //channels and the second spectral window has 10 channel
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 20);
                                 }
                                 else if(allEQ(vb_p->spectralWindows(), 1)) //SPW1
                                 {
                                     nRowsSpw1+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 10);
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

TEST_F(ChannelAverageTVISpwChannTest, CheckMSSelOutputSpwChannels)
{
    useMSSelection(true);
    addPassThroughTVI(true);
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after averaging with chanbin=5 and spw selection "0:1~40",
    //the first spectral window has 8 channels 
    //and there are no rows with the second spectral window 
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
    
    std::cout<<vb_p->spectralWindows()<<vb_p->nChannels()<<std::endl;
    std::cout<<vb_p->antenna1()<<std::endl;
    std::cout<<vb_p->dataDescriptionIds()<<std::endl;
    std::cout<<shape<<std::endl;
    
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 8);
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

TEST_F(ChannelAverageTVISpwChannTest, CheckMSSelOutputTwoAvgSpwChannels)
{
    //useMSSelection(true);
    addPassThroughTVI(true);
    addExtraAvgTVI(true);
    createTVIs();
    size_t nRows = 0;
    size_t nRowsSpw0 = 0;
    size_t nRowsSpw1 = 0;
    //Check that after averaging with chanbin=5 and further averaging with
    //chanbin=(2,5), the first spectral window has 10
    //channels and the second spectral window has 2 channels
    visitIterator([&]() -> void {auto shape = vb_p->visCube().shape(); 
                                 if(allEQ(vb_p->spectralWindows(), 0)) //SPW0 
                                 {
                                     nRowsSpw0+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 10);
                                 }
                                 else if(allEQ(vb_p->spectralWindows(), 1)) //SPW1
                                 {
                                     nRowsSpw1+=shape[2];
                                     ASSERT_EQ(vb_p->nChannels(), 2);
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
        else if (parameter == string("-chanbin"))
        {
            Int tmp = Int(atoi(value.c_str()));
            configuration.define ("chanbin", tmp);
        }
    }

    if (autoMode)
    {
        ::testing::InitGoogleTest(&argc, argv);
        ret = RUN_ALL_TESTS();
    }
    else
    {
        ChannelAverageTVICompareTest test(configuration);
        test.TestBody();
        if (test.getTestResult()) ret = 0;
        else ret = 1;
    }

  return ret;
}
