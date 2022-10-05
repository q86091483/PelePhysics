#ifndef PROBLEMHELPER_H
#define PROBLEMHELPER_H

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
  static std::string base_identifier() { return "ProblemHelper"; }

  void init();

  void init(amrex::Real &a_pInit);

  void initData(amrex::MultiFab &a_state,
                amrex::GeometryData const& a_geomdata);

  void initData(amrex::MultiFab &a_state,
                amrex::GeometryData const& a_geomdata,
                amrex::Real &a_pInit);

  pele::physics::PMF::PmfData::DataContainer const* getPmfData() {return pmf_data.getDeviceData();};

  ~ProblemHelper() override = default;

  virtual void initParams(void **a_params, void **a_params_d);
  virtual void initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit);

  virtual void bcnormal(const amrex::Real* x,
                        const int m_nAux,
                        amrex::Vector<amrex::Real> &s_ext,
                        const int idir,
                        const int sgn,
                        const amrex::Real time,
                        amrex::GeometryData const& geomdata,
                        void* a_params,
                        pele::physics::PMF::PmfData::DataContainer const * pmf_data);

  virtual void initdata_k(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& state,
                          amrex::GeometryData const& geomdata,
                          void* a_params,
                          pele::physics::PMF::PmfData::DataContainer const * pmf_data) = 0;

  virtual void zero_visc(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& beta,
                         amrex::GeometryData const& geomdata,
                         amrex::Box const& domainBox,
                         const int  dir,
                         const int  beta_comp,
                         const int  nComp);

#ifdef AMREX_USE_EB
  // Tag cut cells type: dirichlet, Neumann & what not 
  virtual void tagEBtype() = 0;

  virtual void getEBState() = 0;
#endif

  void *m_baseParams = nullptr;
  void *m_baseParams_d = nullptr;
  pele::physics::PMF::PmfData pmf_data;
};

} // namespace physics
} // namespace pele
#endif
