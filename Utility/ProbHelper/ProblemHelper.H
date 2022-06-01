#ifndef _PROBLEMHELPER_H_
#define _PROBLEMHELPER_H_

#include <AMReX_ParmParse.H>
#include "PelePhysics.H"
#include "Factory.H"
#include "mechanism.H"
#include "turbinflow.H"
#include "PMFData.H"

namespace pele {
namespace physics {

class ProblemHelper : public Factory<ProblemHelper>
{
public:
  static const std::string identifier() { return "ProblemHelper"; }

  virtual void initdata() = 0;

  virtual void bcnormal() = 0;

  virtual ~ProblemHelper() = default;

#ifdef AMREX_USE_EB
  // Tag cut cells type: dirichlet, Neumann & what not 
  virtual void tagEBtype() = 0;

  virtual void getEBState() = 0;
#endif

protected:
  PMF::PmfData pmf_data; 
  turbinflow::TurbInflow tp;
};

} // namespace physics
} // namespace pele
#endif
