
#include <alma/apps/asdm2MS/ScansParser.h>

using namespace std;

#include <vector>
#include <cstdlib>

#include <alma/ASDM/Misc.h>  // asdm::trim_copy comes from here

int strToInt(const string& valStr, int &status) {
  // strip any leading and trailing whitespace.
  string numStr = asdm::trim_copy(valStr);
  int result = -1;

  // and empty string is a problem
  if (numStr.size() == 0) {
    status = 0;
  } else {
    // make sure it really is just digits
    if (numStr.find_first_not_of("0123456789") == string::npos) {
      // and just convert it
      result = atoi(numStr.c_str());
    } else {
      // not allowed
      status = 0;
    }
  }

  return result;
}

// returns 1 on success, 0 on failure. Adds the scan range to scanSet
int parseScanRange(const string& scanRange, set<int> &scanSet) {
  int status = 1;
  // does it have a tilde
  size_t tilPos = scanRange.find("~");
  if (tilPos == 0) {
    // that's a problem
    status = 0;
  } else {
    int scanNumber0, scanNumber1;
    scanNumber0 = scanNumber1 = -1;
    if (tilPos != string::npos) {
      scanNumber0 = strToInt(scanRange.substr(0,tilPos),status);
      if (status > 0) scanNumber1 = strToInt(scanRange.substr(tilPos+1),status);
    } else {
      scanNumber0 = scanNumber1 = strToInt(scanRange,status);
    }
    if (status > 0) {
      for (int i=scanNumber0;i<(scanNumber1+1);i++) {
	scanSet.insert(i);
      }
    }
  }
  return status;
}

// returns 1 on success, 0 on failure. Always clears scanSet before setting in from the list.
int parseScanList(const string& scanList, set<int>& scanSet) {
  int status=1;
  scanSet.clear();

  string strippedScanList = asdm::trim_copy(scanList);
  // an empty list is OK here, empty scan set will be associated with the current EB number
  // so only do something if it has size
  if (strippedScanList.size() > 0) {
  
    size_t last=0;
    size_t next=0;
    while ((status>0) && (next = scanList.find(",",last)) != string::npos) {
      // this also inserts the range into scanSet
      status = parseScanRange(scanList.substr(last,(next-last)),scanSet);
      last = next+1;
    }
    if (status>0) status = parseScanRange(scanList.substr(last),scanSet);
  }
  return status;
}

// find the optional EB number and associated scan list from ebScan
// Add the resulting scan set to ebScanMap, keyed on the EB number (-1 if no EB number)
// returns 1 on success, 0 on failure
int parseEBScan(const string &ebScan, map<int, set<int>> &ebScanMap) {
  set<int> scanSet;
  int status = 1;

  // start from the generic
  int ebNumber = -1;
  string scanList = "";

  // is there a colon ?
  size_t colPos = ebScan.find(":");

  if (colPos == 0) {
    // that's a problem
    status = 0;
  } else {
    if (colPos != string::npos) {
      string ebString = ebScan.substr(0,colPos);
      ebNumber = strToInt(ebString,status);
      scanList = ebScan.substr(colPos+1);
    } else {
      scanList = ebScan;
    }
    if (status >0 ) {
      // parseScanList always clears scanSet before adding to it
      status = parseScanList(scanList, scanSet);
      if (status > 0) {
	ebScanMap[ebNumber].insert(scanSet.begin(),scanSet.end());
      }
    }
  }
  return status;
}

// associates EB number with a set of scans using the scansOptionsValue
// ebScanMap is always cleared before being set
// returns 1 if OK and 0 if there was a problem parsing scansOptionsValue
// ebScanMap will always be empty if the returned status is 0
int scansParser(const string &scansOptionValue, map<int, set<int>> &ebScansMap) {
  int status = 1; // set to 0 on failure
  ebScansMap.clear();

  // break up the string by semi-colon, parsing each piece in turn
  size_t last=0;
  size_t next=0;
  while ((status>0) && (next = scansOptionValue.find(";",last)) != string::npos) {
    status = parseEBScan(scansOptionValue.substr(last,(next-last)),ebScansMap);
    last = next+1;
  }
  if (status>0) {
    status = parseEBScan(scansOptionValue.substr(last),ebScansMap);
  }
  
  if (status<=0) {
    ebScansMap.clear();
  }
  return status;
}
