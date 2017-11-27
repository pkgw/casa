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

#define _POSIX_C_SOURCE 200809L //For mkdtemp(), stpcpy(), nftw()

#include <ftw.h>
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
#include <msvis/MSVis/test/MsFactory.h>
#include <casa/iomanip.h>
#include <gtest/gtest.h>

using namespace std;
using namespace casa;
using namespace casacore;
using namespace casa::vi;
using namespace casa::vi::test;

int removeFile(const char *fpath, const struct stat *sb, int typeflag, 
               struct FTW* ftwbuf);

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

int removeFile(const char *fpath, const struct stat *sb, int typeflag, 
               struct FTW* ftwbuf)
{
  int rv = remove(fpath);
  if(rv)
    perror(fpath);
  return rv;
}
 
TEST( ViiLayerFactoryTest , ViiLayerFactoryBasicTest ) {
 
  // A very rudimentary test of a single layer

  SimpleSimVi2Parameters s0;
  SimpleSimVi2LayerFactory fac(s0);

  Vector<ViiLayerFactory*> facts(1);
  facts[0]=&fac;

  VisibilityIterator2 *vi = new VisibilityIterator2(facts);
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

  VisibilityIterator2 *vi = new VisibilityIterator2(facts);
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
}

class DataSwappingTVI : public TransformingVi2
{
public:
  
  DataSwappingTVI(ViImplementation2 * inputVii) :
    TransformingVi2 (inputVii)
  {
    setVisBuffer(createAttachedVisBuffer (VbRekeyable));
  }

  virtual void visibilityCorrected (casacore::Cube<casacore::Complex> & vis) const
  {
    ViImplementation2* vii = getVii();
    VisBuffer2* vb = vii->getVisBuffer();
    vis = vb->visCube();
    return;
  }
  
  void origin()
  {
    // Drive underlying ViImplementation2
    getVii()->origin();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
  }

  void next()
  {
    // Drive underlying ViImplementation2
    getVii()->next();

    // Synchronize own VisBuffer
    configureNewSubchunk();

    return;
  }

};

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

class DataAccessTest : public ::testing::Test
{
public:
  DataAccessTest()
  {
    //Use the system temp dir, if not defined resort to /tmp
    char * sys_tmpdir = getenv("TMPDIR");
    if(sys_tmpdir != NULL)
      strcpy(tmpdir_p, sys_tmpdir);
    else
      strcpy(tmpdir_p, "/tmp");
    stpcpy (tmpdir_p+strlen(tmpdir_p), "/test_tViiLayerFactory_XXXXXX");
    mkdtemp(tmpdir_p);

    msf_p = std::move(
        std::unique_ptr<MsFactory>(
            new MsFactory(String::format("%s/DataAccessTest.ms", tmpdir_p))));
  }

  void removeCorrectedData()
  {
    msf_p->removeColumn(MS::CORRECTED_DATA);
  }

  void visitIterator(std::function<void(void)> visitor)
  {
    for (vi_p->originChunks (); vi_p->moreChunks(); vi_p->nextChunk()){
      for (vi_p->origin(); vi_p->more (); vi_p->next()){
        visitor();
      }
    }
  }

  void addSwappingDataTVI()
  {
    withSwappingDataTVI_p = true;
  }
  
  void createTVIs()
  {
    //Create MS using the simulator MsFactory
    pair<MeasurementSet *, Int> p = msf_p->createMs();
    ms_p.reset(p.first);

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
    Vector<ViiLayerFactory*> facts(nFac);
    facts[0]=&diskItFac;
    if(withSwappingDataTVI_p)
      facts[1]= swapFac.get();

    //Finally create the top VI
    vi_p.reset(new VisibilityIterator2(facts));

    vb_p = vi_p->getVisBuffer();
  }

  ~DataAccessTest()
  {
    //This will recursively remove everything in the directory
    nftw(tmpdir_p, removeFile, 64, FTW_DEPTH | FTW_PHYS);
  }

  char tmpdir_p[100];
  bool withSwappingDataTVI_p = false;
  std::unique_ptr<casa::vi::test::MsFactory> msf_p;
  std::unique_ptr<VisibilityIterator2> vi_p;
  VisBuffer2 * vb_p;
  std::unique_ptr<MeasurementSet> ms_p;
};
 
TEST_F(DataAccessTest, AccessCorrectedData)
{
  createTVIs();

  //Traverse the iterator accessing the corrected data cube
  visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();std::cout<< vb_p->getSubchunk()<<std::endl;});
}

TEST_F(DataAccessTest, AccessCorrectedDataWhenMissing)
{
//  removeCorrectedData();

//  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This should
  //throw, since it has been removed from the MS.
//  ASSERT_THROW(visitIterator([&]() -> void {vb_p->visCubeCorrected().shape();}),
//               AipsError);
}

TEST_F(DataAccessTest, AccessCorrectedDataInSwappingDataTVI)
{
  removeCorrectedData();

  addSwappingDataTVI();
  
  createTVIs();

  //Traverse the iterator accessing the corrected data cube. This works
  //despite removing the corrected data from disk because the upper TVI layer
  //delivers data when corrected data is requested.
  SCOPED_TRACE("LL");
  visitIterator([&]() -> void {SCOPED_TRACE("II"); vb_p->visCubeCorrected().shape();});
}
