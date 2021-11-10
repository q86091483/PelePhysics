#ifndef _TABULAR_H_
#define _TABULAR_H_

#include "Table.H"

namespace pele {
namespace physics {
namespace eos {

struct Manifold
{
  using eos_type = Manifold;

  static std::string identifier() { return "Manifold"; }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void molecular_weight(amrex::Real* /*mw[]*/)
  {
    amrex::Error("molecular_weight for species does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void inv_molecular_weight(amrex::Real* /*imw[]*/)
  {
    amrex::Error("inv_molecular_weight for species does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& /*WBAR*/)
  {
    amrex::Error("Y2WBAR is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real* /*Ei[]*/)
  {
    amrex::Error("RTY2Ei does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real /*T*/, amrex::Real* /*Hi[]*/)
  {
    amrex::Error("T2Hi does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real* /*Hi[]*/)
  {
    amrex::Error("RTY2Hi does not have significance for Manifold EOS");
  }
   
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real /*Cv[]*/)
  {
    amrex::Error("TY2Cv is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(
	      amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real /*Cv[]*/)
  {
    amrex::Error("RTY2Cv is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real /*Cp[]*/)
  {
    amrex::Error("TY2Cp is not yet implemented for Manifold EOS");
  }  
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(
	      amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real /*Cp[]*/)
  {
    amrex::Error("RTY2Cp is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real /*R*/, amrex::Real /*E*/, amrex::Real* /*Y[]*/, amrex::Real& /*T*/)
  {
    amrex::Error("REY2T is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real /*H*/, amrex::Real Y[], amrex::Real& T)
  {
    tabfunc.get_from_table(tabfunc.idx_T, Y, T);
    // amrex::Error("HY2T is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real /*H*/, amrex::Real* /*Y[]*/, amrex::Real& /*T*/)
  {
    amrex::Error("RHY2T is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real& /*H*/)
  {
    amrex::Error("RHY2T is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(amrex::Real /*X[]*/, amrex::Real* /*Y[]*/)
  {
    amrex::Error("X2Y does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(
	      amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real* /*WDOT[]*/)
  {
    amrex::Error("RTY2WDOT is not yet implemented for Manifold EOS");
  }
    
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(
	       amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real* /*JAC[]*/, const int /*HP*/)
  {
    amrex::Error("RTY2JAC is not yet implemented for Manifold EOS");
  }
    
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(
	      amrex::Real /*R*/, amrex::Real /*T*/, amrex::Real* /*Y[]*/, amrex::Real* /*JAC[]*/)
  {
    amrex::Error("RTY2C does not have significance for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(
	      amrex::Real R, amrex::Real /*T*/, amrex::Real Y[], amrex::Real& P)
  {
    // Get value of Rho consistent with other state variables
    amrex::Real rho_from_table;
    tabfunc.get_from_table(tabfunc.idx_density, Y, rho_from_table);

    // Return pressure is nominal pressure scaled due to desnity mismatch
    P = tabfunc.Pnom * R/rho_from_table;
    
    // amrex::Error("RTY2P is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(amrex::Real /*P*/, amrex::Real* /*Y[]*/, amrex::Real /*T*/, amrex::Real& /*R*/)
  {
    amrex::Error("PYT2R is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2dRdY(amrex::Real Y[], amrex::Real dRdY[])
  {
    tabfunc.get_table_dir_derivs(tabfunc.idx_density, Y, dRdY);
    // amrex::Error("Y2dRdY is not yet implemented for Manifold EOS");
  }
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RY2R(amrex::Real* /*RY[]*/, amrex::Real& /*R*/)
  {
    amrex::Error("RY2R is not yet implemented for Manifold EOS");
  }
  
  AMREX_GPU_HOST_DEVICE
  Manifold(pele::physics::TabFuncParams::TabFuncData* tf_data_in)
  {
    tabfunc = pele::physics::TabFunc(tf_data_in);
  }
  
  template <class... Args>
  AMREX_GPU_HOST_DEVICE
  Manifold(Args...)
  {
  }

private:
  pele::physics::TabFunc tabfunc;
};
  

}
}
}
#endif
