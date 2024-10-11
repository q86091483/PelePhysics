#include "ReactorBase.H"

namespace pele::physics::reactions {

void
ReactorBase::set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals)
{
  const int size_ETV = static_cast<int>(ExtTypVals.size());
  AMREX_ALWAYS_ASSERT(size_ETV == static_cast<int>(m_typ_vals.size()));
  amrex::Vector<std::string> kname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV; i++) {
    m_typ_vals[i] = ExtTypVals[i];
  }

  // cppcheck-suppress knownConditionTrueFalse
  if ((omp_thread == 0) && (verbose > 0)) {
    amrex::Print() << "Set the typical values for PelePhysics ODE integration:"
                   << std::endl;
    for (int i = 0; i < size_ETV - 1; i++) {
      amrex::Print() << kname[i] << " : " << m_typ_vals[i] << std::endl;
    }
    amrex::Print() << "Temp : " << m_typ_vals[size_ETV - 1] << std::endl;
  }
}

#if defined (PELE_USE_AUX) && (NUMAUX > 0)
void
ReactorBase::set_typ_vals_ode_aux(const std::vector<amrex::Real>& ExtTypVals_aux)
{
  const int size_ETV_aux = static_cast<int>(ExtTypVals_aux.size());
  AMREX_ALWAYS_ASSERT(size_ETV_aux == static_cast<int>(m_typ_vals_aux.size()));
  //amrex::Vector<std::string> kname;
  //pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV_aux; i++) {
    m_typ_vals_aux[i] = ExtTypVals_aux[i];
  }

  // cppcheck-suppress knownConditionTrueFalse
  //if ((omp_thread == 0) && (verbose > 0)) {
  //  amrex::Print() << "Set the typical values for PelePhysics ODE integration:"
  //                 << std::endl;
  //  for (int i = 0; i < size_ETV - 1; i++) {
  //    amrex::Print() << kname[i] << " : " << m_typ_vals[i] << std::endl;
  //  }
  //  amrex::Print() << "Temp : " << m_typ_vals[size_ETV - 1] << std::endl;
  //}
}
#endif // #if (NUMAUX > 0) set_typ_vals_ode_aux

} // namespace pele::physics::reactions
