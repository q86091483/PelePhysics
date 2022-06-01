#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

void
CKAWT(amrex::Real* /*awt*/)
{
  amrex::Abort("CKAWT not implemented for Manifold chemistry.");
}

void
CKNCF(int* /*ncf*/)
{
  amrex::Abort("CKNCF not implemented for Manifold chemistry.");
}

void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(1);
  ename[0] = "Null";
}
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(4);
  kname[0] = "X0";
  kname[1] = "X1";
  kname[2] = "X2";
  kname[3] = "XRHO";
}

// compute sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* /*nJdata*/, const int* /*consP*/, int /*NCELLS*/)
{
}

// compute sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* /*nJdata*/, const int* /*consP*/, int /*NCELLS*/)
{
}

// compute sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* /*nJdata*/, const int* /*consP*/)
{
}

// compute sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void
SPARSITY_PREPROC_CSC(
  int* /*rowVals*/
  ,
  int* /*colPtrs*/,
  const int* /*consP*/,
  int /*NCELLS*/)
{
}

// compute sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void
SPARSITY_PREPROC_CSR(
  int* /*colVals*/,
  int* /*rowPtrs*/,
  const int* /*consP*/,
  int /*NCELLS*/,
  int /*base*/)
{
}

// compute sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* /*colVals*/,
  int* /*rowPtr*/,
  const int* /*consP*/,
  int /*NCELLS*/,
  int /*base*/)
{
}

// compute sparsity pattern of the simplified (for precond) system Jacobian on
// CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* /*rowVals*/, int* /*colPtrs*/, int* /*indx*/, const int* /*consP*/)
{
}

// compute sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* /*colVals*/
  ,
  int* /*rowPtr*/,
  const int* /*consP*/,
  int /*base*/)
{
}
#endif
