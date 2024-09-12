#ifndef REACTORUTILS_H
#define REACTORUTILS_H

#include <utility>
#include <AMReX_Gpu.H>
#include "mechanism.H"
#include "PelePhysics.H"
#include "ReactorTypes.H"
#include "AMReX_Sundials.H"
#include "PeleLMeX_Index.H"

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <cvode/cvode.h>

#include <nvector/nvector_serial.h>

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_SYCL
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele::physics::reactions::utils {
int check_flag(void* flagvalue, const char* funcname, int opt);

struct YCOrder
{
  // default: [ [Y_0,...Y_n]_{cell 0}, ..., [Y_0,...Y_n]_{cell m}]
  using ordering_type = YCOrder;
};

struct CYOrder
{
  // [ [Y_{cell 0},...Y_m]_{species 0}, ..., [Y_{cell 0},...Y_m]_{species n}]
  using ordering_type = CYOrder;
};

template <typename T>
struct is_ordering_type
{
  static constexpr bool value =
    std::is_base_of<typename T::ordering_type, T>::value;
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * (NUM_SPECIES + 1 + NUMAUX) + specie_idx;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * NUM_SPECIES + specie_idx;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
aux_index(const int aux_idx, const int icell, const int /*ncells*/)
{
  return icell * NUMAUX + aux_idx;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return specie_idx * ncells + icell;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return vec_index<CYOrder>(specie_idx, icell, ncells);
}

template <typename OrderType>
void
set_sundials_solver_tols(
  // cppcheck-suppress constParameter
  sundials::Context& sunctx,
  void* sundials_mem,
  const int ncells,
  const amrex::Real relTol,
  const amrex::Real absTol,
  const amrex::GpuArray<amrex::Real, NUM_SPECIES + 1 + NUMAUX>& typ_vals,
  const std::string& solvername,
  const int verbose)
{
  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  const int neq_tot = (NUM_SPECIES + 1 + NUMAUX) * ncells;

#if defined(AMREX_USE_CUDA)
  N_Vector atol = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_Vector atol = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Hip(atol);
#elif defined(AMREX_USE_SYCL)
  N_Vector atol = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Sycl(atol);
#else
  N_Vector atol = N_VNew_Serial(neq_tot, sunctx);
  amrex::Real* ratol = N_VGetArrayPointer(atol);
#endif

  if (typ_vals[0] > 0.0) {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances with TypVals rtol = " << relTol
                     << " atolfact = " << absTol << " in PelePhysics \n";
    }

    amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) {
      for (int k = 0; k < NUM_SPECIES + 1 + NUMAUX; k++) {
        ratol[vec_index<OrderType>(k, icell, ncells)] = typ_vals[k] * absTol;
      }
    });
  } else {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances rtol = " << relTol << " atol = " << absTol
                     << " in PelePhysics \n";
    }
    amrex::ParallelFor(
      neq_tot, [=] AMREX_GPU_DEVICE(int i) { ratol[i] = absTol; });
  }

  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  int flag;
  if (solvername == "cvode") {
    flag = CVodeSVtolerances(sundials_mem, relTol, atol);
  } else if (solvername == "arkstep") {
    flag = ARKStepSVtolerances(sundials_mem, relTol, atol);
  } else if (solvername == "erkstep") {
    flag = ERKStepSVtolerances(sundials_mem, relTol, atol);
  } else {
    amrex::Abort("setSundialsSolverTols not implemented for this solver type");
  }
  if (utils::check_flag(&flag, "SVtolerances", 1) != 0) {
    amrex::Abort("Problem in setSundialsSolverTols");
  }

  N_VDestroy(atol);
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_flatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  const int reactor_type,
  const bool clean_init_massfrac,
  amrex::Array4<const amrex::Real> const& rhoY,
  amrex::Array4<const amrex::Real> const& frcExt,
  amrex::Array4<const amrex::Real> const& temperature,
  amrex::Array4<const amrex::Real> const& rhoE,
  amrex::Array4<const amrex::Real> const& frcEExt,
#if (NUMAUX > 0)
  amrex::Array4<const amrex::Real> const& rhoAux,
  amrex::Array4<const amrex::Real> const& frcAuxExt,
#endif
  amrex::Real* y_vect,
  amrex::Real* src_vect,
  amrex::Real* vect_energy,
  amrex::Real* src_vect_energy
#if (NUMAUX > 0)
  , amrex::Real* src_vec_aux
  , amrex::Real* vec_aux_init
#endif
  )
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  amrex::Real rho_old = 0.0;
  if (clean_init_massfrac) {
    for (int n = 0; n < NUM_SPECIES; n++) {
      rho_old += rhoY(i, j, k, n);
    }
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    const int idx = vec_index<OrderType>(n, icell, ncells);
    const int idx2 = spec_index<OrderType>(n, icell, ncells);
    const amrex::Real rhoYn =
      (clean_init_massfrac)
        ? amrex::min<amrex::Real>(
            amrex::max<amrex::Real>(0.0, rhoY(i, j, k, n)), rho_old)
        : rhoY(i, j, k, n);
    y_vect[idx] = rhoYn;
    src_vect[idx2] = frcExt(i, j, k, n);
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (amrex::Real& Yn : mass_frac) {
    Yn *= rho_inv;
  }

  amrex::Real temp = temperature(i, j, k, 0);
  vect_energy[icell] = rhoE(i, j, k, 0);
  src_vect_energy[icell] = frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = vect_energy[icell] * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }
  y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] = temp;

#if (NUMAUX > 0)
  for (int n = 0; n < NUMAUX; n++) {
    const int idx = vec_index<OrderType>(NUM_SPECIES + 1 + n, icell, ncells);
    const int idx_src = aux_index<OrderType>(n, icell, ncells);
    y_vect[idx] = rhoAux(i, j, k, idx_src);
    src_vec_aux[idx_src] = frcAuxExt(i, j, k, idx_src);
    vec_aux_init[idx_src] = rhoAux(i, j, k, idx_src);
  }

  if (false) {
    amrex::Print() << "In box_flatten: ----------------" << std::endl;
    amrex::Print() << " rho = " << rho << std::endl;
    for (int n = 0; n < NUM_SPECIES; n++) {
      const int idx = vec_index<OrderType>(n, icell, ncells);
      const int idx_src = spec_index<OrderType>(n, icell, ncells);
      amrex::Print() << " y[" << n << "] = " << y_vect[idx] << ", idx=" << idx;
      amrex::Print() << ", ydot[" << n << "] = " << src_vect[idx_src] << ", idx_src=" << idx_src << std::endl;
    }
    amrex::Print() << " y[NUM_SPECIES] = " << y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)]
                    << ", idx=" << vec_index<OrderType>(NUM_SPECIES, icell, ncells);
    amrex::Print() << ", ydot[NUM_SPECIES] = " << src_vect_energy[icell] << ", idx_src=" << icell << std::endl;
    for (int n = 0; n < NUMAUX; n++) {
      const int idx = vec_index<OrderType>(NUM_SPECIES + 1 + n, icell, ncells);
      const int idx_src = aux_index<OrderType>(n, icell, ncells);
      amrex::Print() << " y[NUM_SPECIES + " << n << "] = " << y_vect[idx] << ", idx=" << idx;
      amrex::Print() << ", ydot[NUM_SPECIES + " << n << "] = " << src_vec_aux[idx_src] << ", idx_src=" << idx_src << std::endl;
    }
    amrex::Print() << " " << std::endl;
  } // if(true or false) for output
#endif // #if (NUMAUX > 0)
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_unflatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  const int reactor_type,
  const bool clean_init_massfrac,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
#if (NUMAUX > 0)
  amrex::Array4<amrex::Real> const& rhoAux,
#endif
  amrex::Array4<amrex::Real> const& FC_in,
  const amrex::Real* y_vect,
  const amrex::Real* vect_energy,
  long int const& FCunt,
  amrex::Real dt)
{
  // Get the old Ys & rho as well as the clipped version
  // that went into CVODE
  amrex::Real rhoY_oldclipped[NUM_SPECIES];
  amrex::Real rho_old = 0.0;
  amrex::Real rho_oldclipped = 0.0;
  if (clean_init_massfrac) {
    for (int n = 0; n < NUM_SPECIES; n++) {
      rho_old += rhoY(i, j, k, n);
    }
    for (int n = 0; n < NUM_SPECIES; n++) {
      rhoY_oldclipped[n] = amrex::min<amrex::Real>(
        amrex::max<amrex::Real>(0.0, rhoY(i, j, k, n)), rho_old);
      rho_oldclipped += rhoY_oldclipped[n];
    }
  }

  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    // With clean_init_massfrac, new rhoY equal old one + CVODE increment
    const amrex::Real rhoYn =
      (clean_init_massfrac)
        ? rhoY(i, j, k, n) + (y_vect[vec_index<OrderType>(n, icell, ncells)] -
                              rhoY_oldclipped[n]) *
                               rho_old / rho_oldclipped
        : y_vect[vec_index<OrderType>(n, icell, ncells)];
    rhoY(i, j, k, n) = rhoYn;
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (amrex::Real& Yn : mass_frac) {
    Yn *= rho_inv;
  }

  amrex::Real temp = y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];
  rhoE(i, j, k, 0) = vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = rhoE(i, j, k, 0) * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }
  temperature(i, j, k, 0) = temp;
#if (NUMAUX > 0)
  for (int n = 0; n < NUMAUX; n++) {
    rhoAux(i, j, k, n) = y_vect[vec_index<OrderType>(NUM_SPECIES + 1 + n, icell, ncells)];
  }
#endif
  FC_in(i, j, k, 0) = static_cast<amrex::Real>(FCunt);
}

template <typename OrderType>
class FlattenOps
{
public:
  static_assert(
    is_ordering_type<OrderType>::value,
    "FlattenOps must have ordering type as its template argument");

  void flatten(
    const amrex::Box& box,
    const int ncells,
    const int reactor_type,
    const bool clean_init_massfrac,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
#if (NUMAUX > 0)
    amrex::Array4<const amrex::Real> const& rhoAux,
    amrex::Array4<const amrex::Real> const& frcAuxExt,
#endif
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy
#if (NUMAUX > 0)
    , amrex::Real* src_vect_aux
    , amrex::Real* vec_aux_init
#endif
    )
  {
    BL_PROFILE("Pele::flatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_flatten<OrderType>(
        icell, i, j, k, ncells, reactor_type, clean_init_massfrac,
        rhoY, frcExt, temperature, rhoE, frcEExt,
#if (NUMAUX > 0)
        rhoAux, frcAuxExt,
#endif
        y_vect, src_vect, vect_energy, src_vect_energy
#if (NUMAUX > 0)
        , src_vect_aux
        , vec_aux_init
#endif
        );
    });
  }

  void unflatten(
    const amrex::Box& box,
    const int ncells,
    const int reactor_type,
    const bool clean_init_massfrac,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
#if (NUMAUX > 0)
    amrex::Array4<amrex::Real> const& rhoAux,
#endif
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
    long int* FCunt,
    amrex::Real dt)
  {
    BL_PROFILE("Pele::unflatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_unflatten<OrderType>(
        icell, i, j, k, ncells, reactor_type, clean_init_massfrac, rhoY,
        temperature, rhoE, frcEExt,
#if (NUMAUX > 0)
        rhoAux,
#endif
        FC_in, y_vect, vect_energy, FCunt[icell],
        dt);
    });
  }
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
fKernelSpec(
  const int icell,
  const int ncells,
  const amrex::Real dt_save,
  const int reactor_type,
  const amrex::Real* yvec_d,
  amrex::Real* ydot_d, // NOLINT(readability-non-const-parameter)
  const amrex::Real* rhoe_init,
  const amrex::Real* rhoesrc_ext,
  const amrex::Real* rYs
#if (NUMAUX > 0)
  , const amrex::Real* rhoAuxsrc_ext
  , const amrex::Real* rhoAux_init
#endif
  )
{
  amrex::Real rho_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = yvec_d[vec_index<OrderType>(n, icell, ncells)];
    rho_pt += massfrac[n];
  }
  const amrex::Real rho_pt_inv = 1.0 / rho_pt;

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] *= rho_pt_inv;
  }

  const amrex::Real nrg_pt =
    (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) * rho_pt_inv;

  amrex::Real temp_pt =
    yvec_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];

  amrex::Real Cv_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Ei(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cv(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Hi(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cp(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> cdots_pt = {0.0};
  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  amrex::Real rhoesrc = rhoesrc_ext[icell];
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real cdot_rYs =
      cdots_pt[n] + rYs[spec_index<OrderType>(n, icell, ncells)];
    ydot_d[vec_index<OrderType>(n, icell, ncells)] = cdot_rYs;
    rhoesrc -= cdot_rYs * ei_pt[n];
  }
  ydot_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] =
    rhoesrc * (rho_pt_inv / Cv_pt);

#if (NUMAUX > 0) // Advance ODE dU/dt = f(U) + g

  // Add g to r.h.s.
  for (int n = 0; n < NUMAUX; n++) {
    const int idx = vec_index<OrderType>(NUM_SPECIES + 1 + MIXF_IN_AUX + n, icell, ncells);
    const int idx_src = aux_index<OrderType>(n, icell, ncells);
    ydot_d[idx] = rhoAuxsrc_ext[idx_src];
  }

  // Ad f(U) to r.h.s.
  for (int n = 0; n < NUMAGE; n++) {
    const int idx_rhoAge = vec_index<OrderType>(NUM_SPECIES + 1 + AGE_IN_AUX + n, icell, ncells);
    const int idx_rhoZ = vec_index<OrderType>(NUM_SPECIES + 1 + MIXF_IN_AUX + n, icell, ncells);
    amrex::Real alpha_1 = yvec_d[idx_rhoAge] / yvec_d[idx_rhoZ];
    amrex::Real alpha_2 = dt_save + rhoAux_init[idx_rhoAge] / (rhoAux_init[idx_rhoZ] + 1E-12);
    amrex::Real alpha_min = std::min(alpha_1, alpha_2);
    amrex::Print() << "alpha_1: " << alpha_1 << ", alpha_2: " << alpha_2 << std::endl;
    if ((yvec_d[idx_rhoZ] / rho_pt) > 1E-3) {
      ydot_d[idx_rhoAge] = yvec_d[idx_rhoZ] + alpha_1 * ydot_d[idx_rhoZ];
    } else {
      ydot_d[idx_rhoAge] = 0.0;
    }
  }

  if (false) {
    amrex::Print() << "In fKernelSpec: ----------------" << std::endl;
    //amrex::Print() << "MIXF_IN_AUX, " << MIXF_IN_AUX << ", " << NUMMIXF << std::endl;
    //amrex::Print() << "AGE_IN_AUX, " << AGE_IN_AUX << ", " << NUMAGE << std::endl;
    //amrex::Print() << "AGEPV_IN_AUX, " << AGEPV_IN_AUX << ", " << NUMAGEPV << std::endl;
    amrex::Print() << " rho = " << rho_pt << std::endl;
    for (int n = 0; n < NUM_SPECIES; n++) {
      const int idx = vec_index<OrderType>(n, icell, ncells);
      const int idx_src = spec_index<OrderType>(n, icell, ncells);
      amrex::Print() << " y[" << n << "] = " << yvec_d[idx] << ", idx=" << idx;
      amrex::Print() << ", ydot[" << n << "] = " << ydot_d[idx_src] << ", idx_src=" << idx_src << std::endl;
    }
    amrex::Print() << " y[NUM_SPECIES] = " << yvec_d[NUM_SPECIES] << ", idx=" << NUM_SPECIES;
    amrex::Print() << ", ydot[NUM_SPECIES] = " << ydot_d[NUM_SPECIES] << ", idx_src=" << NUM_SPECIES << std::endl;
    for (int n = 0; n < NUMAUX; n++) {
      const int idx = vec_index<OrderType>(NUM_SPECIES + 1 + n, icell, ncells);
      const int idx_src = aux_index<OrderType>(n, icell, ncells);
      amrex::Print() << " y[FIRSTAUX + " << n << "] = " << yvec_d[idx] << ", idx=" << idx;
      amrex::Print() << ", ydot[FIRSTAUX + " << n << "] = " << ydot_d[idx_src] << ", idx_src=" << idx_src << std::endl;
    }
    amrex::Print() << " " << std::endl;
  } // if (true or not) for output
#endif
}

#ifdef AMREX_USE_GPU
N_Vector
setNVectorGPU(int nvsize, int atomic_reduction, amrex::gpuStream_t stream);
#endif
} // namespace pele::physics::reactions::utils
#endif
