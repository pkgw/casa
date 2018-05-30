/*
 * OptionalTables.h
 *
 * All OptionalTables classes should have a public static method, Generate.
 * It takes references to casacore::Table object and reader object as an argument.
 *
 *  Created on: May 30, 2018
 *      Author: nakazato
 */

#ifndef _SINGLEDISHFILLER_FILLER_OPTIONALTABLES_H_
#define _SINGLEDISHFILLER_FILLER_OPTIONALTABLES_H_

#include <iostream>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ScaColDesc.h>
//#include <casacore/ms/MeasurementSets/MeasurementSet.h>

// empty OptionalTables class
template<class Reader>
class NullOptionalTables {
public:
  static void Generate(casacore::Table &/*table*/, Reader const &/*reader*/) {
    //std::cout << "This is default. NullOptionalTables::Generate" << std::endl;
  }
};


// OptionalTables class for NRO data
template<class Reader>
class NROOptionalTables {
public:
  static void Generate(casacore::Table &table, Reader const &reader) {
    // generate NRO_ARRAY table
    Generate_NRO_ARRAY(table, reader);
  }

private:
  static void Generate_NRO_ARRAY(casacore::Table &table, Reader const &reader) {
    casacore::String const nro_tablename = "NRO_ARRAY";

    casacore::TableDesc td(nro_tablename, casacore::TableDesc::Scratch);
    td.addColumn(casacore::ScalarColumnDesc<casacore::Int>("ARRAY"));
    td.addColumn(casacore::ScalarColumnDesc<casacore::Int>("BEAM"));
    td.addColumn(casacore::ScalarColumnDesc<casacore::Int>("POLARIZATION"));
    td.addColumn(casacore::ScalarColumnDesc<casacore::Int>("SPECTRAL_WINDOW"));
    casacore::String tabname = table.tableName() + "/" + nro_tablename;
    casacore::SetupNewTable newtab(tabname, td, casacore::Table::Scratch);
    table.rwKeywordSet().defineTable(nro_tablename,
        Table(newtab, reader.getNROArraySize()));

    casacore::Table nro_table = table.rwKeywordSet().asTable(nro_tablename);
    casacore::ScalarColumn<int> arr(nro_table, "ARRAY");
    casacore::ScalarColumn<int> bea(nro_table, "BEAM");
    casacore::ScalarColumn<int> pol(nro_table, "POLARIZATION");
    casacore::ScalarColumn<int> spw(nro_table, "SPECTRAL_WINDOW");
    for (int iarr = 0; iarr < reader.getNROArraySize(); ++iarr) {
      arr.put(iarr, iarr);
      if (reader.isNROArrayUsed(iarr)) {
        bea.put(iarr, reader.getNROArrayBeamId(iarr));
        pol.put(iarr, reader.getNROArrayPol(iarr));
        spw.put(iarr, reader.getNROArraySpwId(iarr));
      } else {
        // array is not used, fill with -1
        bea.put(iarr, -1);
        pol.put(iarr, -1);
        spw.put(iarr, -1);
      }
    }
  }
};



#endif /* CODE_SINGLEDISHFILLER_FILLER_OPTIONALTABLES_H_ */
