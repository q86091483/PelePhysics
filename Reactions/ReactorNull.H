#ifndef REACTORNULL_H
#define REACTORNULL_H

#include "ReactorBase.H"

namespace pele::physics::reactions {

class ReactorNull : public ReactorBase::Register<ReactorNull>
{
public:
  static std::string identifier() { return "ReactorNull"; }

  using Ordering = utils::YCOrder;

  int init(int reactor_type, int ncells) override;

  void close() override {}

  void print_final_stats(void* /*mem*/) override {}

  int react(
    const amrex::Box& box,
    amrex::Array4<amrex::Real> const& rY_in,
    amrex::Array4<amrex::Real> const& rYsrc_in,
    amrex::Array4<amrex::Real> const& T_in,
    amrex::Array4<amrex::Real> const& rEner_in,
    amrex::Array4<amrex::Real> const& rEner_src_in,
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
    amrex::Array4<amrex::Real> const& rAux_in,
    amrex::Array4<amrex::Real> const& rAux_src_in,
#endif
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Array4<int> const& mask,
    amrex::Real& dt_react,
    amrex::Real& time
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  int react(
    amrex::Real* rY_in,
    amrex::Real* rYsrc_in,
    amrex::Real* rX_in,
    amrex::Real* rX_src_in,
    amrex::Real& dt_react,
    amrex::Real& time,
    int ncells
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  void flatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
    amrex::Array4<const amrex::Real> const& rhoAux,
    amrex::Array4<const amrex::Real> const& frcAuxExt,
#endif
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy
#if defined (PELE_USE_AUX) && (NUMNEW > 0)
    , amrex::Real* y_vect_aux
    , amrex::Real* src_vect_aux
    , amrex::Real* vec_aux_init
#endif
    ) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac,
      rhoY, frcExt, temperature, rhoE, frcEExt,
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
      rhoAux, frcAuxExt,
#endif
      y_vect, src_vect, vect_energy, src_vect_energy
#if defined (PELE_USE_AUX) && (NUMNEW > 0)
      , y_vect_aux
      , src_vect_aux
      , vec_aux_init
#endif
      );
  }

  void unflatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
    amrex::Array4<amrex::Real> const& rhoAux,
#endif
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
#if defined (PELE_USE_AUX) && (NUMNEW > 0)
    amrex::Real* y_vect_aux,
#endif
    long int* FCunt,
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac,
      rhoY, temperature, rhoE, frcEExt,
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
      rhoAux,
#endif
      FC_in, y_vect, vect_energy,
#if defined (PELE_USE_AUX) && (NUMNEW > 0)
      y_vect_aux,
#endif
      FCunt, dt);
  }

private:
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
  bool m_clean_init_massfrac{false};
};
} // namespace pele::physics::reactions
#endif
