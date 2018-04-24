/**
   \file gaussmodel.hpp
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>

   
*/

#ifndef _BNMIN1_TEST_GAUSSMODEL_HPP__
#define _BNMIN1_TEST_GAUSSMODEL_HPP__

#include <numeric>
#include <vector>
#include <cmath>

#include "../minimmodel.hxx"

namespace Minim {
  
  /**
     Very simple likelihood model in which the negative log-likelihood
     is proportional to the distance from the origin of the parameters
     of the model. Useful for thesting of likelihood characterisation
     algorithms.
   */
  class GaussObs:
    public Minim::MLikelihood
  {
    
  public:
    
    /// Parameters of the model
    std::vector<double> p;    
    
    /// Scale parameter
    double sigma;
    
    /**
       \param n Number of parameters/dimensions
    */
    GaussObs(size_t n):
      p(n),
      sigma(1.0)
    {
    }
    
    double lLikely(void) const
    {
      double p_p = 0;

      std::accumulate( p.begin( ), p.end( ), p.begin( ),
                       [&](std::vector<double>::const_iterator &p1, double p2) -> std::vector<double>::const_iterator &{
                           p_p += *p1 * p2;
                           return ++p1;
                       }
                     );

      return p.size()* 0.5*std::log(2*M_PI*std::pow(sigma,2))+ p_p/(2*std::pow(sigma,2));
    }
    
    void AddParams(std::vector< Minim::DParamCtr > &pars);
    
    
  };
  
}

#endif

