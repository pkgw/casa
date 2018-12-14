
#ifndef _NAME2TABLE_H_
#define _NAME2TABLE_H_
/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File Name2Table.h
 */
#include "ASDMTableBase.h"

#include <map>
#include <set>

using namespace std;

class Name2Table {
 private:
  static map<string, ASDM_TABLE_BASE*> name2Table_;
  static bool init_;
  static bool init();

  static set<ASDM_TABLE_BASE*> table_;

 public:
  static const set<ASDM_TABLE_BASE*>& find(const std::vector<string>& name,bool verbose=false);
};
#endif // _NAME2TABLE_H_
