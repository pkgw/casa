/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.cxx

*/


#include "bnmin_main.hxx"
#include "../config.h"

namespace Minim {

  const char * version(void)
  {
    return PACKAGE_VERSION;
  }

  BaseErr::BaseErr(const std::string &s):
    std::runtime_error(s)
  {
  }

  NParsErr::NParsErr(const std::string &fname, size_t expected, size_t received) :
      BaseErr( std::string("In function ") + fname + std::string(" expected ") +
               std::to_string(expected) + std::string(" but received ") +
               std::to_string(received) + std::string(" pars "))
  {
  }
    



}


