
#include <gtest/gtest.h>

#include <singledish/Filler/test/SampleReader.h>
#include <singledish/Filler/SingleDishMSFiller.h>

#include <casacore/casa/Logging/LogIO.h>

// This is kind of demonstration for beginner of filler framework
TEST(SingleDishMSFillerTest, DemonstrateSampleReader) {
	casacore::LogIO os(LogOrigin("GTEST", "SampleReaderTest", WHERE));
	os << "Create filler object" << casacore::LogIO::POST;
	casa::SingleDishMSFiller<SampleReader> filler("mysampledata.nro", false);
	filler.fill();
	os << "Write MS as \"mysampledata.ms\"" << casacore::LogIO::POST;
	filler.save("mysampledata.ms");
}

int main(int nArgs, char * args[]) {
  ::testing::InitGoogleTest(&nArgs, args);
  return RUN_ALL_TESTS();
}
