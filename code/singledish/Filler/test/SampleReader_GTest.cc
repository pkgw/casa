#include <gtest/gtest.h>

#include <singledish/Filler/test/SampleReader.h>
#include <singledish/Filler/SingleDishMSFiller.h>

#include <casacore/casa/Logging/LogIO.h>

// This is kind of demonstration for beginner of filler framework
TEST(SingleDishMSFillerTest, DemonstrateSampleReader) {
  casacore::LogIO os(LogOrigin("GTEST", "SampleReaderTest", WHERE));
  os << "This is a test program to demonstrate how SingleDishMSFiller<SampleReader>\n"
     << "works. Usage of the filler is only three steps below:\n"
     << "\n"
     << "    int main(int argc, char *argv[]) {\n"
     << "(1)   casa::SingleDishMSFiller<SampleReader> filler(\"mysampledata.nro\", false);\n"
     << "(2)   filler.fill();\n"
     << "(3)   filler.save(\"mysampledata.ms\");\n"
     << "      return 0;\n"
     << "    }\n"
     << "\n" << casacore::LogIO::POST;
  os << "Here I will show you how Filler/Reader works:\n"
     << "\n" << casacore::LogIO::POST;
  os << "=== Step 1. Create filler object with SampleReader\n"
     << "=== (1)   casa::SingleDishMSFiller<SampleReader> filler(\"mysampledata.nro\", false);\n"
     << "  (NB: input data name \"mysampledata.nro\" is dummy. SampleReader generates\n"
     << "       the data on-the-fly)\n"
     << "\n" << casacore::LogIO::POST;
  casa::SingleDishMSFiller<SampleReader> filler("mysampledata.nro", false);
  os << "=== Step 2. Fill MS\n"
     << "=== (2)   filler.fill();\n"
     << "\n" << casacore::LogIO::POST;
  filler.fill();
  os << "=== Step 3. Write MS as \"mysampledata.ms\"\n"
     << "=== (3)   filler.save(\"mysampledata.ms\");\n"
     << "\n" << casacore::LogIO::POST;
  filler.save("mysampledata.ms");
  os << "Now you should have \"mysampledata.m\" on your current working directory.\n"
     << "Please open it with, e.g., casabrowser to see the contents!\n"
     << casacore::LogIO::POST;
}

int main(int nArgs, char * args[]) {
  ::testing::InitGoogleTest(&nArgs, args);
  return RUN_ALL_TESTS();
}
