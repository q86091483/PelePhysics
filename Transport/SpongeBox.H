#ifndef SPONGEBOX_H
#define SPONGEBOX_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

namespace pele {
namespace physics {
namespace transport {

// template <typename ParmType>
// class Params
// {
// public:
//   Params() {}

//   ~Params() {}

//   void allocate()
//   {
//     m_h_parm.initialize();
//     if (!m_device_allocated) {
//       m_d_parm = (ParmType)amrex::The_Device_Arena()->alloc(sizeof(m_h_parm));
//       m_device_allocated = true;
//       sync_to_device();
//     }
//   }

//   void deallocate()
//   {
//     if (m_device_allocated)
//       amrex::The_Device_Arena()->free(m_d_parm);
//   }

//   void sync_to_device()
//   {
//     if (!m_device_allocated) {
//       amrex::Abort("Device params not allocated yet");
//     } else {

//       // would have liked to use this but apparently only works for trivial
//       // types
//       // amrex::Gpu::copy(
//       //   amrex::Gpu::hostToDevice, &m_h_parm, &m_h_parm + 1,
//       //   m_d_parm);
// #ifdef AMREX_USE_GPU
//       amrex::Gpu::htod_memcpy(m_d_parm, &m_h_parm, sizeof(m_h_parm));
// #else
//       std::memcpy(m_d_parm, &m_h_parm, sizeof(m_h_parm));
// #endif
//     }
//   }

//   ParmType& host_trans_parm() { return m_h_parm; }
//   const ParmType* device_trans_parm() { return m_d_parm; }

// private:
//   ParmType m_h_parm;
//   ParmType* m_d_parm;
//   bool m_device_allocated{false};
// };

struct SpongeBoxParm
{
  amrex::Box box;
  amrex::Real tanh_offset_x{0.0};
  amrex::Real tanh_offset_y{0.0};
  amrex::Real tanh_offset_z{0.0};
  amrex::Real tanh_slope_x{0.0};
  amrex::Real tanh_slope_y{0.0};
  amrex::Real tanh_slope_z{0.0};
  amrex::Real multiplier{1.0};
  SpongeBoxParm() {}
  ~SpongeBoxParm() {}
  void initialize()
  {
    amrex::ParmParse pp("sponge");
    pp.query("tanh_offset_x", tanh_offset_x);
    pp.query("tanh_offset_y", tanh_offset_y);
    pp.query("tanh_offset_z", tanh_offset_z);
    pp.query("tanh_slope_x", tanh_slope_x);
    pp.query("tanh_slope_y", tanh_slope_y);
    pp.query("tanh_slope_z", tanh_slope_z);
    pp.query("multiplier", multiplier);
  }
};

struct NoOpSpongeStruct
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void sponge() noexcept {}
};

struct SpongeBoxStruct
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void sponge(
    const amrex::Box bx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
    SpongeBoxParm const& parms,
    amrex::Array4<amrex::Real> const& mu) noexcept
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      const amrex::Real z = plo[2] + (k + 0.5) * dx[2];
      for (int j = lo.y; j <= hi.y; ++j) {
        const amrex::Real y = plo[1] + (j + 0.5) * dx[1];
        for (int i = lo.x; i <= hi.x; ++i) {
          const amrex::Real x = plo[0] + (i + 0.5) * dx[0];
          const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
          if (parms.box.contains(iv)) {
            AMREX_D_TERM(
              const amrex::Real etax =
                0.5 *
                (1.0 + tanh(-(x - parms.tanh_offset_x) / parms.tanh_slope_x));
              , const amrex::Real etay =
                  0.5 *
                  (1.0 + tanh(-(y - parms.tanh_offset_y) / parms.tanh_slope_y));
              ,
              const amrex::Real 0.5 *
                (1.0 + tanh(-(z - parms.tanh_offset_z) / parms.tanh_slope_z)););

            const amrex::Real weight = AMREX_D_TERM(etax, *etay, *etaz);
            mu(iv) =
              (1.0 - weight) * mu(iv) + weight * parms.multiplier * mu(iv);
          }
        }
      }
    }
  }
};

template <typename SpongeStruct, typename... Args>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
sponge(Args&&... args)
{
  SpongeStruct::sponge(std::forward<Args>(args)...);
}

} // namespace transport
} // namespace physics
} // namespace pele
#endif
