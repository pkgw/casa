#include <iostream>
#include <alma/ASDM/TableStreamReader.h>
#include <alma/ASDM/PointingTable.h>
#include <alma/ASDM/PointingRow.h>
#include <alma/ASDM/SysPowerRow.h>
#include <alma/ASDM/SysPowerTable.h>

using namespace std;
using namespace asdm;

int main(int argC, char *argV[]) {
  if (argC < 2) {
    cout << "Usage : TestTableStreamReader <asdm-directory>" << endl;
    exit (1);
  }

  TableStreamReader<PointingTable, PointingRow> tsrPointing;
  try {
    tsrPointing.open(string(argV[1]));

    int numRows = 0;
    while (tsrPointing.hasRows()) {
      const vector<PointingRow *>& rows = tsrPointing.nextNRows(10);
      cout << "I have read " << rows.size() << " rows." << endl;
      numRows += rows.size();
    }
    cout << numRows << " read in total." << endl;
    tsrPointing.close();
  }
  catch (ConversionException e) {
    cout << e.getMessage() << endl;
  }

  TableStreamReader<SysPowerTable, SysPowerRow> tsrSysPower;
  try {
    tsrSysPower.open(string(argV[1]));
    
    unsigned int numRows = 0;
    unsigned int nBytes = 10 * 1024 * 1024;
    while (tsrSysPower.hasRows()) {
      const vector<SysPowerRow*>& rows = tsrSysPower.untilNBytes(nBytes);
      cout << "nBytes = " << nBytes << ", I have read " << rows.size() << " rows." << endl;
      numRows += rows.size();
    }
    cout << numRows << " read in total." << endl;
    tsrSysPower.close();
  }
  catch (ConversionException e) {
    cout << e.getMessage() << endl;
  }

}
