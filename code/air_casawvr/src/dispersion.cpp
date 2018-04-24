/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version August 2010.
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2

   \file dispersion.cpp
*/
#include <iostream>
#include <fstream>
#include <algorithm>

#include <string>
#include <vector>
#include <cctype>
#include <cstring>

#include "dispersion.hpp"

namespace LibAIR2 {

  double DispersionTab::operator() (double fnu)
  {
    std::pair<double, double> b= *lower_bound(fnu);
    if (b.first==fnu)
      return b.second;

    std::pair<double, double> l= *(--lower_bound(fnu));
    std::pair<double, double> u= *upper_bound(fnu);
    
    const double f=(fnu-l.first)/(u.first-l.first);
    return l.second+ f*(u.second-l.second);
    
  }

  inline std::string trim(const std::string &s)
  {
    auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){return std::isspace(c);});
    auto wsback=std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
    return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
  }

  std::vector<std::string> tokLine(const std::string &input) {
    char *buf = strdup(input.c_str( ));
    std::vector<std::string> result;
    char *elem = strtok(buf,",;");
    while ( elem ) {
      result.push_back(std::move(std::string(elem)));
      elem = strtok(0,",;");
    }
    free(buf);
    return result;
  }

  void loadCSV(const char *fname,
	       DispersionTab &dt)
  {
    std::ifstream ifs(fname);
    if (not ifs.good())
    {
      throw std::runtime_error(std::string("Could not open dispersion table ")+fname);
    }
    std::string   scratch;

    while(ifs.good())
    {
      std::getline(ifs, scratch);
      if (scratch.size() < 5) continue;

      auto tok = tokLine(scratch);
      if (tok.size() < 2) continue;

      std::string first  = trim(tok[0]);
      std::string second = trim(tok[1]);

      try {
          dt.insert( dt.end(), std::pair<double, double>(std::stod(first), std::stod(second)) );
      } catch (const std::bad_cast &bc) {
          std::cerr << "Could not interpret " << first << " and " << second << std::endl;
      }
    }
  }


}



