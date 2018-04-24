/**
   \file nestedgauss.hpp

   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>
   
   Nested sampling exploration of the simple gaussian model, for
   testing purposes
*/

#include "../minimmodel.hxx"
#include "../apps/gaussmodel.hxx"
#include "../nestedsampler.hxx"
#include "../mcmonitor.hxx"
#include "../priors.hxx"
#include <memory>

struct pdesc {
  std::shared_ptr<Minim::IndependentFlatPriors> obs;
  std::shared_ptr<Minim::NestedS> s;
  std::shared_ptr<Minim::SOutMCMon> mon;
};

pdesc mkDesc(double l_sigma,
	     bool monitor,
	     double /*sigma*/=0.1,
	     size_t ndim=3)
{
  using namespace Minim;
  GaussObs *gp = new GaussObs(ndim);
  gp->sigma=l_sigma;

  pdesc res;
  res.obs=std::shared_ptr<IndependentFlatPriors>(new  IndependentFlatPriors(gp));

  for (size_t i=0; i<ndim; ++i)
  {
      res.obs->AddPrior(std::string("p")+std::to_string(i), -1.01,1.01);
  }

  
  std::list<Minim::MCPoint> startset;

  startSetDirect(*res.obs,
		 20,
		 startset);

  res.s=std::shared_ptr<NestedS>(new NestedS(*res.obs,
					       startset));

  if (monitor)
  {
    SOutMCMon *pp=new SOutMCMon();
    res.s->mon=pp;
    res.mon=std::shared_ptr<SOutMCMon>(pp);
  }
  return res;
};

double getEvidence(double l_sigma,
		   size_t nsample,
		   bool monitor)
{
  pdesc d=mkDesc(l_sigma,
		 monitor);

  const double res=d.s->sample(nsample);

  return res;
};



