// How to compile:
//    * requirement is that produced executable is portable, i.e.
//      it can run without casacore libraries at runtime
//    1. build casacore to produce static libraries (-DENABLE_SHARED=NO)
//    2. compile PThreadUtil.cc to produce PThreadUtil.o
//    3. compile SingleDishMSFiller.cc to produce SingleDishFiller.o
//    4. compile this file, statically link casacore libraries and object files
//       created at the steps 2 and 3
//
// Makefile which implements the above steps will be provided externally.

#include <singledishfiller/Filler/SingleDishMSFiller.h>
#include <casacore/casa/Logging/LogIO.h>
#include "SampleReader.h"

#include <string>
#include <iostream>

void usage(char const *command) {
  std::string const command_string(command);
  auto const pos = command_string.find_last_of('/');
  std::string basename(command);
  if (pos != std::string::npos) {
    basename = command_string.substr(pos+1);
  }
  std::cout << "Usage: " << std::endl;
  std::cout << "    " << basename << " [parallel|serial]" << std::endl;
}

int main(int argc, char *argv[]) {
  casacore::LogIO os(LogOrigin("", "SampleReaderDemo", WHERE));

  // parsing command line option
  // usage is
  //     SampleReaderDemo [parallel|serial]
  bool parallel = false;
  if (argc == 2) {
    std::string const parallel_key("parallel");
    std::string const serial_key("serial");
    std::string const arg(argv[1]);
    if (arg == parallel_key) {
      parallel = true;
    } else if (arg == serial_key) {
      parallel = false;
    } else {
      usage(argv[0]);
      return 1;
    }
  } else if (argc > 2) {
    usage(argv[0]);
    return 1;
  }

  // start processing
  os <<"This is a test program to demonstrate how SingleDishMSFiller<SampleReader>\n"
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
  casa::SingleDishMSFiller<SampleReader> filler("mysampledata.nro", parallel);
  os << "=== Step 2. Fill MS\n"
     << "=== (2)   filler.fill();\n"
     << "\n" << casacore::LogIO::POST;
  filler.fill();
  os << "=== Step 3. Write MS as \"mysampledata.ms\"\n"
     << "=== (3)   filler.save(\"mysampledata.ms\");\n"
     << "\n" << casacore::LogIO::POST;
  filler.save("mysampledata.ms");
  os << "Now you should have \"mysampledata.ms\" on your current working directory.\n"
     << "Please open it with, e.g., casabrowser to see the contents!\n"
     << casacore::LogIO::POST;

  return 0;
}

