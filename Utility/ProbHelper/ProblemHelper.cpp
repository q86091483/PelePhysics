#include "ProblemHelper.H"
#include "AMReX_MultiFab.H"

using namespace amrex;

namespace pele {
namespace physics {

void ProblemHelper::init()
{
    initParams(&m_baseParams,&m_baseParams_d);
}

void ProblemHelper::init(amrex::Real &a_pInit)
{
    initParams(&m_baseParams,&m_baseParams_d,a_pInit);
}

void ProblemHelper::initParams(void **a_params, void **a_params_d)
{
    amrex::ignore_unused(a_params,a_params_d);
    amrex::Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::initParams(void **a_params, void **a_params_d, amrex::Real &a_pInit)
{
    amrex::ignore_unused(a_params,a_params_d,a_pInit);
    amrex::Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::initData(amrex::MultiFab &a_state,
                             amrex::GeometryData const& a_geomdata)
{
   pele::physics::PMF::PmfData::DataContainer const* lpmfdata   = pmf_data.getDeviceData();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto  const &state_arr   = a_state.array(mfi);
      amrex::ParallelFor(bx, [=]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         initdata_k(i, j, k, state_arr, a_geomdata, m_baseParams_d, lpmfdata);
      });
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ProblemHelper::bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int /*m_nAux*/,
  amrex::Vector<amrex::Real> &/*s_ext*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  void* /*a_params*/,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ProblemHelper::zero_visc(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& /*beta*/,
                         amrex::GeometryData const& /*geomdata*/,
                         amrex::Box const& /*domainBox*/,
                         const int  /*dir*/,
                         const int  /*beta_comp*/,
                         const int  /*nComp*/)
{
  amrex::ignore_unused(i,j,k);
}

#ifdef AMREX_USE_EB
void ProblemHelper::tagEBtype()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::getEBState()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}
#endif

} // namespace physics
} // namespace pele
