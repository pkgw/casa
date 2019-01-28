//# TestUtilsTVI.cc:  This file contains the implementation of the ChannelAverageTVI class.
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

#define _XOPEN_SOURCE 700 //For nftw(), stpcpy(), mkdtemp()
#define _DARWIN_C_SOURCE //im macOS mkdtemp() is not available if _POSIX_C_SOURCE=200809L (Apple bug report #35851865)

#include <ftw.h>
#include <string.h>
#include <limits.h>
#include <unistd.h> //im macOS mkdtemp() is not defined in stdlib.h as POSIX dictates (Apple bug report #35830645)
#include <msvis/MSVis/test/TestUtilsTVI.h>
#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/TransformingVi2.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

int removeFile(const char *fpath, const struct stat *sb, int typeflag, 
               struct FTW* ftwbuf);

int removeFile(const char *fpath, const struct stat *sb, int typeflag, 
               struct FTW* ftwbuf)
{
  (void)sb;  //Unused vars
  (void)typeflag;
  (void)ftwbuf;
    
  int rv = remove(fpath);
  if(rv)
    perror(fpath);
  return rv;
}


//////////////////////////////////////////////////////////////////////////
// FreqAxisTVITest class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
FreqAxisTVITest::FreqAxisTVITest():
		autoMode_p(true), testResult_p(true)
{

}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
FreqAxisTVITest::FreqAxisTVITest(Record configuration):
		autoMode_p(false), testResult_p(true)
{
    configuration.get (configuration.fieldNumber ("inputms"), inpFile_p);
    testFile_p = inpFile_p + String(".test");
    referenceFile_p = inpFile_p + String(".ref");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
FreqAxisTVITest::~FreqAxisTVITest()
{
    TearDown();

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVITest::init(Record &configuration)
{
    initTestConfiguration(configuration);
    initReferenceConfiguration(configuration);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVITest::SetUp()
{
    // Generate test file
    generateTestFile();

    // Generate reference file
    generateReferenceFile();

    return;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void FreqAxisTVITest::TearDown()
{
    //Recursively remove the test and reference files
    nftw(testFile_p.c_str(), removeFile, 64, FTW_DEPTH | FTW_PHYS);
    nftw(referenceFile_p.c_str(), removeFile, 64, FTW_DEPTH | FTW_PHYS);

    if (autoMode_p)
        nftw(inpFile_p.c_str(), removeFile, 64, FTW_DEPTH | FTW_PHYS);

    return;
}

MsFactoryTVITester::MsFactoryTVITester(const std::string& testSubdir, 
                                       const std::string& msName) :
    msName_p(msName)
{
    //Use the system temp dir, if not defined or too long resort to /tmp
    char * sys_tmpdir = getenv("TMPDIR");
    if(sys_tmpdir != NULL &&
       strlen(sys_tmpdir) < _POSIX_PATH_MAX - 1 - testSubdir.size() + 8)
        strncpy(tmpdir_p, sys_tmpdir, strlen(sys_tmpdir)+1);
    else
        strncpy(tmpdir_p, "/tmp", 5);
    stpcpy (tmpdir_p+strlen(tmpdir_p), (std::string("/")+testSubdir+"_XXXXXX").c_str());
    
}

void MsFactoryTVITester::SetUp()
{
    mkdtemp(tmpdir_p);
    msf_p.reset(new casa::vi::test::MsFactory(String::format
        ("%s/%s.ms", tmpdir_p,msName_p.c_str())));
}
  
void MsFactoryTVITester::createMS()
{
    //Create MS using the simulator MsFactory
    std::pair<MeasurementSet *, Int> p = msf_p->createMs();
    ms_p.reset(p.first); //MsFactory has given up ownership
}

void MsFactoryTVITester::instantiateVI(std::vector<ViiLayerFactory*>& factories)
{
    //Create the top VI using the factories provided
    vi_p.reset(new VisibilityIterator2(factories));

    vb_p = vi_p->getVisBuffer();
}

void MsFactoryTVITester::visitIterator(std::function<void(void)> visitor)
{
    for (vi_p->originChunks (); vi_p->moreChunks(); vi_p->nextChunk())
    {
        for (vi_p->origin(); vi_p->more (); vi_p->next())
        {
            visitor();
        }
    }
}

casa::vi::test::MsFactory& MsFactoryTVITester::getMsFactory()
{
    return *msf_p;
}
  
void MsFactoryTVITester::TearDown()
{
    //The MS destructor will update the file system, so deleting it before removing the directory 
    msf_p.reset();
    ms_p.reset();
    vi_p.reset();
    //This will recursively remove everything in the directory
    nftw(tmpdir_p, removeFile, 64, FTW_DEPTH | FTW_PHYS);
}

MsFactoryTVITester::~MsFactoryTVITester()
{
}


//////////////////////////////////////////////////////////////////////////
// Convenience methods
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template <class T> void compareVector(	const Char* column,
										const Vector<T> &inp,
										const Vector<T> &ref,
										Float tolerance)
{
	// Check matching shape
    ASSERT_EQ(inp.size(), ref.size()) 
        << " test and reference vectors don't have the same size";

	// Compare values
	for (uInt index=0;index < inp.size(); index++)
	{
	    ASSERT_NEAR(abs(inp(index) - ref(index)), 0, tolerance)
            << column << " does not match in position ="
            << index
            << " test=" << inp(index)
            << " reference=" << ref(index);
	}
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template <class T> void compareMatrix(	const Char* column,
										const Matrix<T> &inp,
										const Matrix<T> &ref,
										Float tolerance)
{
	// Check matching shape
    ASSERT_EQ(inp.shape(), ref.shape()) 
        << " test and reference matrices don't have the same shape";

	// Compare values
	const IPosition &shape = inp.shape();
	for (uInt row=0;row < shape(1); row++)
	{
		for (uInt col=0;col < shape(0); col++)
		{
		    ASSERT_NEAR(abs(inp(col,row) - ref(col,row)), 0, tolerance)
                << column << " does not match in position (row,col)="
                << "("<< row << "," << col << ")"
                << " test=" << inp(col,row)
                << " reference=" << ref(col,row);
		}
	}
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template <class T> void compareCube(const Char* column,
                                    const Cube<T> &inp,
                                    const Cube<T> &ref,
                                    Float tolerance)
{
    // Check matching shape
    ASSERT_EQ(inp.shape(), ref.shape()) 
         << " test and reference cubes don't have the same shape";

    // Compare values
    const IPosition &shape = inp.shape();
    for (uInt row=0;row < shape(2); row++)
    {
        for (uInt chan=0;chan < shape(1); chan++)
        {
            for (uInt corr=0;corr < shape(0); corr++)
            {
                ASSERT_NEAR(abs(inp(corr,chan,row) - ref(corr,chan,row)), 0, tolerance)
	                << column << " does not match in position (corr,chan,row)="
			        << "("<< corr << "," << chan << "," << row << ")"
                    << " test=" << inp(corr,chan,row)
                    << " reference=" << ref(corr,chan,row);
            }
        }
    }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void compareVisibilityIterators(VisibilityIterator2 &testTVI,
                                VisibilityIterator2 &refTVI,
                                VisBufferComponents2 &columns,
                                Float tolerance,
                                std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap)
{
    // Declare working variables
    String columnName;
    Int chunk = 0,buffer = 0;

    // Get VisBuffers
    VisBuffer2 *refVb = refTVI.getImpl()->getVisBuffer();
    VisBuffer2 *testVb = testTVI.getImpl()->getVisBuffer();

    // Compare selected columns
    refTVI.originChunks();
    testTVI.originChunks();
    while (refTVI.moreChunks() and testTVI.moreChunks())
    {
        chunk += 1;
        buffer = 0;

        refTVI.origin();
        testTVI.origin();

        while (refTVI.more() and testTVI.more())
        {
            buffer += 1;
            SCOPED_TRACE(string("Comparing chunk ") + std::to_string(chunk) + 
                         " buffer " + std::to_string(buffer) +
                         " Spw " + std::to_string(refVb->spectralWindows()[0]) + 
                         " scan " + std::to_string(refVb->scan()[0]));

            if (columns.contains(VisBufferComponent2::NRows))
            {
                SCOPED_TRACE("Comparing NRows component ");
                ASSERT_EQ(testVb->nRows() , refVb->nRows());
            }

            if (columns.contains(VisBufferComponent2::NChannels))
            {
                SCOPED_TRACE("Comparing NChannels component ");
                ASSERT_EQ(testVb->nChannels(), refVb->nChannels());
            }

            if (columns.contains(VisBufferComponent2::NCorrelations))
            {
                SCOPED_TRACE("Comparing NCorrelations component ");
                ASSERT_EQ(testVb->nCorrelations(), refVb->nCorrelations());
            }

            if (columns.contains(VisBufferComponent2::NAntennas))
            {
                SCOPED_TRACE("Comparing NAntennas component ");
                ASSERT_EQ(testVb->nAntennas(), refVb->nAntennas());
            }

            if (columns.contains(VisBufferComponent2::Time))
            {
                SCOPED_TRACE("Comparing Time component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Time);
                compareVector(columnName.c_str(),testVb->time(),refVb->time(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::TimeCentroid))
            {
                SCOPED_TRACE("Comparing TimeCentroid component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::TimeCentroid);
                compareVector(columnName.c_str(),testVb->timeCentroid(),refVb->timeCentroid(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::TimeInterval))
            {
                SCOPED_TRACE("Comparing TimeInterval component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::TimeInterval);
                compareVector(columnName.c_str(),testVb->timeInterval(),refVb->timeInterval(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::Exposure))
            {
                SCOPED_TRACE("Comparing Exposure component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Exposure);
                compareVector(columnName.c_str(),testVb->exposure(),refVb->exposure(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::SpectralWindows))
            {
                SCOPED_TRACE("Comparing SpectralWindows component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::SpectralWindows);
                compareVector(columnName.c_str(),testVb->spectralWindows(),
                              refVb->spectralWindows(), 0);
            }

            if (columns.contains(VisBufferComponent2::Antenna1))
            {
                SCOPED_TRACE("Comparing Antenna1 component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Antenna1);
                compareVector(columnName.c_str(),testVb->antenna1(),
                              refVb->antenna1(), 0);
            }

            if (columns.contains(VisBufferComponent2::Antenna2))
            {
                SCOPED_TRACE("Comparing Antenna2 component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Antenna2);
                compareVector(columnName.c_str(),testVb->antenna2(),
                              refVb->antenna2(), 0);
            }

            if (columns.contains(VisBufferComponent2::DataDescriptionIds))
            {
                SCOPED_TRACE("Comparing DataDescriptionIds component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::DataDescriptionIds);
                compareVector(columnName.c_str(),testVb->dataDescriptionIds(),
                              refVb->dataDescriptionIds(), 0);
            }

            if (columns.contains(VisBufferComponent2::PolarizationId))
            {
                SCOPED_TRACE("Comparing PolarizationId component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::PolarizationId);
                ASSERT_EQ(testVb->polarizationId(), refVb->polarizationId());
            }

            if (columns.contains(VisBufferComponent2::RowIds))
            {
                SCOPED_TRACE("Comparing RowIds component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::RowIds);
                compareVector(columnName.c_str(),testVb->rowIds(), 
                              refVb->rowIds(), 0);
            }

            if (columns.contains(VisBufferComponent2::Uvw))
            {
                SCOPED_TRACE("Comparing Uvw component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Uvw);
                compareMatrix(columnName.c_str(),testVb->uvw(),refVb->uvw(), 0);
            }

            if (columns.contains(VisBufferComponent2::FlagRow))
            {
                SCOPED_TRACE("Comparing FlagRow component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::FlagRow);
                compareVector(columnName.c_str(),testVb->flagRow(),refVb->flagRow(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::FlagCube))
            {
                SCOPED_TRACE("Comparing FlagCube component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::FlagCube);
                compareCube(columnName.c_str(),testVb->flagCube(),refVb->flagCube(),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::VisibilityCubeObserved))
            {
                SCOPED_TRACE("Comparing VisibilityCubeObserved component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeObserved);
                compareCube(columnName.c_str(),testVb->visCube(),getViscube(refVb,MS::DATA,datacolmap),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::VisibilityCubeCorrected))
            {
                SCOPED_TRACE("Comparing VisibilityCubeCorrected component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeCorrected);
                compareCube(columnName.c_str(),testVb->visCubeCorrected(),getViscube(refVb,MS::CORRECTED_DATA,datacolmap),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::VisibilityCubeModel))
            {
                SCOPED_TRACE("Comparing VisibilityCubeModel component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeModel);
                compareCube(columnName.c_str(),testVb->visCubeModel(),getViscube(refVb,MS::MODEL_DATA,datacolmap),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::VisibilityCubeFloat))
            {
                SCOPED_TRACE("Comparing VisibilityCubeFloat component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeFloat);
                compareCube(columnName.c_str(),testVb->visCubeFloat(),refVb->visCubeFloat(),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::WeightSpectrum))
            {
                SCOPED_TRACE("Comparing WeightSpectrum component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::WeightSpectrum);
                compareCube(columnName.c_str(),testVb->weightSpectrum(),refVb->weightSpectrum(),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::SigmaSpectrum))
            {
                SCOPED_TRACE("Comparing SigmaSpectrum component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::SigmaSpectrum);
                compareCube(columnName.c_str(),testVb->sigmaSpectrum(),refVb->sigmaSpectrum(),
                            tolerance);
            }

            if (columns.contains(VisBufferComponent2::Weight))
            {
                SCOPED_TRACE("Comparing Weight component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Weight);
                compareMatrix(columnName.c_str(),testVb->weight(),refVb->weight(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::Sigma))
            {
                SCOPED_TRACE("Comparing Sigma component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Sigma);
                compareMatrix(columnName.c_str(),testVb->sigma(),refVb->sigma(),
                              tolerance);
            }

            if (columns.contains(VisBufferComponent2::Frequencies))
            {
                SCOPED_TRACE("Comparing Frequencies component ");
                columnName = VisBufferComponents2::name(VisBufferComponent2::Frequencies);
                compareVector(columnName.c_str(),testVb->getFrequencies(0),refVb->getFrequencies(0),
                              tolerance);
            }

            refTVI.next();
            testTVI.next();
        }

        refTVI.nextChunk();
        testTVI.nextChunk();
    }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void copyTestFile(String &path,String &filename,String &outfilename)
{
    if (path.size() > 0)
    {
        char* pathChar = getenv ("CASAPATH");
        if (pathChar != NULL)
        {
            // Get base path
            String pathStr(pathChar);
            String res[2];
            casacore::split(pathChar,res,2,String(" "));

            // Generate full qualified filename
            String fullfilename(res[0]);
            fullfilename += path + "/" + filename;

            // Remove any previously existing copy
            String rm_command = String ("rm -rf ") + outfilename;
            ASSERT_TRUE(system(rm_command.c_str()) == 0)
                << "Failed to remove " <<outfilename;

            // Make a copy of the file in the working directory
            String cp_command = String ("cp -r ") + fullfilename + String(" ") + outfilename;
            ASSERT_TRUE(system(cp_command.c_str()) == 0)
                << "Failed to copy " <<fullfilename<<" to "<<outfilename;
        }
        else
        {
            FAIL() << "CASAPATH environmental variable not defined ";
        }
    }
    else
    {
        // Remove any previously existing copy
        String rm_command = String ("rm -rf ") + outfilename;
        ASSERT_TRUE(system(rm_command.c_str()) == 0)
            << "Failed to remove " <<outfilename;

        // Make a copy of the file in the working directory
        String cp_command = String ("cp -r ") + filename + String(" ") + outfilename;
        ASSERT_TRUE(system(cp_command.c_str()) != 0)
            << "Failed to copy " << filename << " to "<<outfilename;
    }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
const Cube<Complex> & getViscube(	VisBuffer2 *vb,
									MS::PredefinedColumns datacol,
									std::map<casacore::MS::PredefinedColumns,casacore::MS::PredefinedColumns> *datacolmap)
{
    MS::PredefinedColumns mappeddatacol;
    if (datacolmap == NULL)
    {
        mappeddatacol = datacol;
    }
    else
    {
        if (datacolmap->find(datacol) != datacolmap->end())
        {
            mappeddatacol = datacolmap->at(datacol);
        }
        else
        {
            mappeddatacol = datacol;
        }
    }


	switch (mappeddatacol)
	{
		case MS::DATA:
		{
			return vb->visCube();
			break;
		}
		case MS::CORRECTED_DATA:
		{
			return vb->visCubeCorrected();
			break;
		}
		case MS::MODEL_DATA:
		{
			return vb->visCubeModel();
			break;
		}
		default:
		{
			return vb->visCube();
			break;
		}
	}
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void flagEachOtherChannel(VisibilityIterator2 &vi, bool unfoldChanbin, int chanbin)
{
    // Declare working variables
    Int chunk = 0,buffer = 0;

    // Get VisBuffer
    VisBuffer2 *vb = vi.getVisBuffer();

    // Propagate flags
    vi.originChunks();
    while (vi.moreChunks())
    {
        chunk += 1;
        buffer = 0;

        vi.origin();
        vi.origin();

        while (vi.more())
        {
            buffer += 1;

            // Initialize flag cube
            IPosition shape = vb->getShape();
            Cube<Bool> flagCube(shape,false);

            // Switch each other buffer the sign of the flag of the first block of channels
            Bool firstChanBlockFlag = buffer % 2? true:false;

            // Fill flag cube alternating flags per blocks channels
            size_t nCorr = shape(0);
            size_t nChan = shape(1);
            size_t nRows = shape(2);
            for (size_t row_i =0;row_i<nRows;row_i++)
            {
                // Row completely flagged
                if (row_i % 2)
                {
                    flagCube.xyPlane(row_i) = true;
                }
                else
                {
                    for (size_t chan_i =0;chan_i<nChan;chan_i++)
                    {
                        // Set the flags in each other block of channels
                        Bool chanBlockFlag;
                        if(unfoldChanbin)
                            chanBlockFlag = ((chan_i / chanbin) % 2 )? firstChanBlockFlag:!firstChanBlockFlag;
                        else
                            chanBlockFlag = chan_i % 2? firstChanBlockFlag:!firstChanBlockFlag;

                        for (size_t corr_i =0;corr_i<nCorr;corr_i++)
                        {
                            flagCube(corr_i,chan_i,row_i) = chanBlockFlag;
                        }
                    }
                }
            }

            vi.writeFlag(flagCube);
            vi.next();
        }

        vi.nextChunk();
    }

    return;
}

} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END
