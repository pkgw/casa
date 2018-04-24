/** 
    \file twoerrline_imp.cxx

    Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>
    Initial version 2009    


*/

#include "twoerrline_ml.hxx"
#include <algorithm>
#include <numeric>
#include <cmath>


namespace Minim {

  LineTwoErrML::LineTwoErrML(const std::vector<double> &xvals,
			     const std::vector<double> &yvals,
			     double sigmax,
			     double sigmay):
    xobs(xvals.size()),
    yobs(yvals.size()),
    sigmax(sigmax),
    sigmay(sigmay),
    nobs(xvals.size())
  {
    for (size_t i=0; i<xobs.size(); ++i)
    {
      xobs[i]=xvals[i];
      yobs[i]=yvals[i];
    }
  }

  void LineTwoErrML::residuals(std::vector<double> &res) const
  {
    res.resize(0);
    std::accumulate( yobs.begin( ), yobs.end( ), xobs.begin( ),
                     [&](std::vector<double>::const_iterator &x, double y) -> std::vector<double>::const_iterator &{
                         res.push_back(y-(*x)*a-b);
                         return ++x;
                     }
                   );
  }

  double LineTwoErrML::lLikely(void) const
  {
    std::vector<double> ub(xobs.size(),b);
    double ressq = 0;
    accumulate( yobs.begin( ), yobs.end( ), xobs.begin( ),
                [&](std::vector<double>::const_iterator &x, double y) -> std::vector<double>::const_iterator &{
                    ressq += std::pow( y - (*x)*a-b, 2 );
                    return ++x;
                }
              );
    const double r=0.5*ressq/(pow(sigmay,2)+ pow(sigmax*a,2));
    return r;
  }

  void LineTwoErrML::lGrd(std::vector< double > &res) const
  {
    res.resize(2);
    std::vector<double> rr;

    std::accumulate( yobs.begin( ), yobs.end( ), xobs.begin( ),
                     [&](std::vector<double>::const_iterator &x, double y) -> std::vector<double>::const_iterator &{
                         rr.push_back(y-(*x)*a-b);
                         return ++x;
                     }
                   );
    
    const double st=(std::pow(sigmay,2)+ std::pow(sigmax*a,2));
    double xobs_rr = 0;
    std::accumulate( rr.begin( ), rr.end( ), xobs.begin( ),
                     [&](std::vector<double>::const_iterator &x, double r) -> std::vector<double>::const_iterator &{
                         xobs_rr += *x * r;
                         return ++x;
                     }
                   );

    double rr_rr = 0;
    std::accumulate( rr.begin( ), rr.end( ), rr.begin( ),
                     [&](std::vector<double>::iterator &r1, double r2) -> std::vector<double>::iterator &{
                         rr_rr += *r1 * r2;
                         return ++r1;
                     }
                   );

    res[0]= -1.0* xobs_rr / st - rr_rr / std::pow(st,2)* a* std::pow(sigmax,2);

    double sum_rr = 0;
    for (auto &v : rr) sum_rr += v;

    res[1]= -1.0* sum_rr / st;
  }

  LineTwoErr_LavMarq::LineTwoErr_LavMarq(const std::vector<double> &xvals,
					 const std::vector<double> &yvals,
					 double sigmax,
					 double sigmay):
    m(xvals,
      yvals,
      sigmax,
      sigmay)
  {
    m.a=1.0;
    m.b=1.0;
  }

  unsigned LineTwoErr_LavMarq::nres (void) const
  {
    return m.nobs;
  }


  void  LineTwoErr_LavMarq::residuals (std::vector<double> &res) const
  {
    res.resize(m.nobs);
    std::vector<double> rr(m.nobs);
    m.residuals(rr);
    const double st= std::pow((std::pow(m.sigmay,2)+ std::pow(m.sigmax*m.a,2)), 0.5);

    std::transform( rr.begin(), rr.end(), rr.begin(), [&](double v){ return v / st; } );

    std::copy(rr.begin(),
	      rr.end(),
	      res.begin());
  }

  void LineTwoErr_LavMarq::AddParams(std::vector< Minim::DParamCtr > &pars)
  {
    m.AddParams(pars);
  }



}




