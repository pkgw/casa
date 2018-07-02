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

#include <msvis/MSVis/test/TestUtilsTVI.h>

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
                                dataColMap *datacolmap)
{
    // Declare working variables
    String columnName;
    Int chunk = 0,buffer = 0;

    // Get VisBuffers
    VisBuffer2 *refVb = refTVI.getImpl()->getVisBuffer();
    VisBuffer2 *testVb = testTVI.getImpl()->getVisBuffer();

    // Compare selected columns
    try
    {
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
                SCOPED_TRACE(string("Comparing chunk ") + to_string(chunk) + 
                             " buffer " + to_string(buffer) +
                             " Spw " + to_string(refVb->spectralWindows()[0]) + 
                             " scan " + to_string(refVb->scan()[0]));

                if (columns.contains(VisBufferComponent2::NRows))
                    ASSERT_EQ(testVb->nRows() , refVb->nRows());

                if (columns.contains(VisBufferComponent2::NChannels))
                    ASSERT_EQ(testVb->nChannels(), refVb->nChannels());

                if (columns.contains(VisBufferComponent2::NCorrelations))
                    ASSERT_EQ(testVb->nCorrelations(), testVb->nCorrelations());

                if (columns.contains(VisBufferComponent2::FlagRow))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::FlagRow);
                    compareVector(columnName.c_str(),testVb->flagRow(),refVb->flagRow(),
                                  tolerance);
                }

                if (columns.contains(VisBufferComponent2::FlagCube))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::FlagCube);
                    compareCube(columnName.c_str(),testVb->flagCube(),refVb->flagCube(),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::VisibilityCubeObserved))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeObserved);
                    compareCube(columnName.c_str(),testVb->visCube(),getViscube(refVb,MS::DATA,datacolmap),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::VisibilityCubeCorrected))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeCorrected);
                    compareCube(columnName.c_str(),testVb->visCubeCorrected(),getViscube(refVb,MS::CORRECTED_DATA,datacolmap),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::VisibilityCubeModel))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeModel);
                    compareCube(columnName.c_str(),testVb->visCubeModel(),getViscube(refVb,MS::MODEL_DATA,datacolmap),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::VisibilityCubeFloat))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::VisibilityCubeFloat);
                    compareCube(columnName.c_str(),testVb->visCubeFloat(),refVb->visCubeFloat(),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::WeightSpectrum))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::WeightSpectrum);
                    compareCube(columnName.c_str(),testVb->weightSpectrum(),refVb->weightSpectrum(),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::SigmaSpectrum))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::SigmaSpectrum);
                    compareCube(columnName.c_str(),testVb->sigmaSpectrum(),refVb->sigmaSpectrum(),
                                tolerance);
                }

                if (columns.contains(VisBufferComponent2::Weight))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::Weight);
                    compareMatrix(columnName.c_str(),testVb->weight(),refVb->weight(),
                                  tolerance);
                }

                if (columns.contains(VisBufferComponent2::Sigma))
                {
                    columnName = VisBufferComponents2::name(VisBufferComponent2::Sigma);
                    compareMatrix(columnName.c_str(),testVb->sigma(),refVb->sigma(),
                                  tolerance);
                }

                if (columns.contains(VisBufferComponent2::Frequencies))
                {
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
    catch (AipsError &ex)
    {
        FAIL()<< "Exception comparing visibility iterators: " << ex.getMesg() 
                      << endl << "Stack Trace: " << ex.getStackTrace();
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
									dataColMap *datacolmap)
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
