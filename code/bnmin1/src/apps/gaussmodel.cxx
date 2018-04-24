/**
g   \file gaussmodel.cpp
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>

*/

#include "gaussmodel.hxx"

namespace Minim {
  
  void GaussObs::AddParams(std::vector< Minim::DParamCtr > &pars)
  {
    using namespace Minim;
    for (size_t i=0; i<p.size(); ++i)
    {
      ParamCtr<double> pa( &p[i], 
                           std::string("p") + std::to_string(i), 
                           true, 
                           "n-th parameter");
      pars.push_back(pa);
  }
    
  }
  
}
