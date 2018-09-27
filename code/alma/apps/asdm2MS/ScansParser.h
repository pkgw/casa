// Boost grammer code to handle the scan selection argument to asdm2MS and bdflags2MS

//
// A collection of declarations and functions used for the parsing of the 'scans' option.
//

#if !defined(ALMA_SCANSPARSER_H)
#define ALMA_SCANSPARSER_H

#include <map>
#include <set>
#include <string>

// Parse the scans options string used by asdm2MS and bdflags2MS.
// sets a map of scans set associated with Execution Block numbers using the scansOptionsValue
// returns 1 on success and 0 on failure.
// ebScansMap is always empty on failure
int scansParser(const std::string &scansOptionValue, std::map<int, std::set<int> > &ebScansMap);

#endif
