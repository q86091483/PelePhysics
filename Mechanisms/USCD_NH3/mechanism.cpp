#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  8,  15, 45, 49, 52, 4,  5,  6,  7,  40, 57, 0,  1,  2,  3,  9,
  10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
  27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43,
  44, 46, 47, 48, 50, 51, 53, 54, 55, 56, 58, 59, 60, 61, 62, 63};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[NUM_GAS_REACTIONS] = {
    4, 4, 4, 3, 2, 3, 2, 3, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 3, 4, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4};
  const int kiv[NUM_GAS_REACTIONS * 4] = {
    1,  2,  3,  4,  0,  3,  1,  4,  0,  4,  1,  6,  6,  3,  4,  0,  1,  0,  0,
    0,  1,  4,  6,  0,  3,  2,  0,  0,  1,  3,  4,  0,  1,  2,  5,  0,  1,  5,
    4,  0,  1,  5,  0,  2,  1,  5,  6,  3,  5,  3,  2,  4,  5,  4,  6,  2,  5,
    4,  6,  2,  4,  7,  0,  0,  5,  7,  2,  0,  5,  7,  2,  0,  1,  7,  0,  5,
    1,  7,  6,  4,  7,  4,  6,  5,  7,  4,  6,  5,  7,  3,  5,  4,  18, 3,  8,
    9,  8,  2,  9,  3,  8,  4,  1,  9,  1,  12, 0,  8,  12, 3,  1,  9,  12, 4,
    1,  16, 12, 4,  6,  8,  12, 2,  16, 3,  12, 9,  1,  11, 12, 9,  18, 4,  1,
    13, 0,  12, 13, 3,  1,  16, 13, 2,  17, 3,  13, 4,  6,  12, 8,  13, 1,  18,
    13, 9,  6,  18, 13, 9,  14, 4,  15, 1,  13, 0,  1,  15, 0,  13, 15, 3,  13,
    4,  15, 4,  6,  13, 14, 2,  5,  18, 14, 1,  18, 0,  1,  14, 0,  18, 14, 3,
    1,  11, 14, 4,  6,  18, 1,  9,  16, 0,  1,  16, 0,  9,  16, 4,  6,  9,  11,
    18, 3,  0,  1,  11, 18, 4,  1,  11, 18, 4,  11, 3,  9,  0,  11, 4,  5,  18,
    10, 9,  3,  0,  5,  9,  10, 4,  1,  10, 9,  4,  10, 3,  9,  2,  17, 3,  16,
    4,  17, 2,  16, 5,  17, 5,  7,  16};
  const int nuv[NUM_GAS_REACTIONS * 4] = {
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 2, 0, -2, 1,  0, 0,
    -1, -1, 1, 0, -2, 1,  0, 0, -1, -1, 1, 0, -1, -1, 1, 0, -1, -1, 2, 0,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -2, 1,  0, 0, -2, 1,  1, 0, -2, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 2, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 0,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 2, 0, -1, -1, 1, 1, -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 4;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 4 + j] + 1;
        nu[j] = nuv[(i - 1) * 4 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[21]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 21; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 64; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[21];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.007000; // N
  awt[1] = 39.950000; // Ar
  awt[2] = 4.002602;  // He
  awt[3] = 1.008000;  // H
  awt[4] = 15.999000; // O
  awt[5] = 12.011000; // C
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 6;
  // Zero ncf
  for (int id = 0; id < kd * 21; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 3] = 2; // H

  // H
  ncf[1 * kd + 3] = 1; // H

  // O2
  ncf[2 * kd + 4] = 2; // O

  // O
  ncf[3 * kd + 4] = 1; // O

  // OH
  ncf[4 * kd + 3] = 1; // H
  ncf[4 * kd + 4] = 1; // O

  // HO2
  ncf[5 * kd + 3] = 1; // H
  ncf[5 * kd + 4] = 2; // O

  // H2O
  ncf[6 * kd + 3] = 2; // H
  ncf[6 * kd + 4] = 1; // O

  // H2O2
  ncf[7 * kd + 3] = 2; // H
  ncf[7 * kd + 4] = 2; // O

  // N
  ncf[8 * kd + 0] = 1; // N

  // NO
  ncf[9 * kd + 0] = 1; // N
  ncf[9 * kd + 4] = 1; // O

  // NO2
  ncf[10 * kd + 0] = 1; // N
  ncf[10 * kd + 4] = 2; // O

  // N2O
  ncf[11 * kd + 0] = 2; // N
  ncf[11 * kd + 4] = 1; // O

  // NH
  ncf[12 * kd + 3] = 1; // H
  ncf[12 * kd + 0] = 1; // N

  // NH2
  ncf[13 * kd + 3] = 2; // H
  ncf[13 * kd + 0] = 1; // N

  // N2H
  ncf[14 * kd + 3] = 1; // H
  ncf[14 * kd + 0] = 2; // N

  // NH3
  ncf[15 * kd + 3] = 3; // H
  ncf[15 * kd + 0] = 1; // N

  // HNO
  ncf[16 * kd + 3] = 1; // H
  ncf[16 * kd + 0] = 1; // N
  ncf[16 * kd + 4] = 1; // O

  // H2NO
  ncf[17 * kd + 3] = 2; // H
  ncf[17 * kd + 0] = 1; // N
  ncf[17 * kd + 4] = 1; // O

  // N2
  ncf[18 * kd + 0] = 2; // N

  // AR
  ncf[19 * kd + 1] = 1; // Ar

  // HE
  ncf[20 * kd + 2] = 1; // He
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(6);
  ename[0] = "N";
  ename[1] = "Ar";
  ename[2] = "He";
  ename[3] = "H";
  ename[4] = "O";
  ename[5] = "C";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(21);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O2";
  kname[3] = "O";
  kname[4] = "OH";
  kname[5] = "HO2";
  kname[6] = "H2O";
  kname[7] = "H2O2";
  kname[8] = "N";
  kname[9] = "NO";
  kname[10] = "NO2";
  kname[11] = "N2O";
  kname[12] = "NH";
  kname[13] = "NH2";
  kname[14] = "N2H";
  kname[15] = "NH3";
  kname[16] = "HNO";
  kname[17] = "H2NO";
  kname[18] = "N2";
  kname[19] = "AR";
  kname[20] = "HE";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (Jac[22 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 22;
    int offset_col = nc * 22;
    for (int k = 0; k < 22; k++) {
      for (int l = 0; l < 22; l++) {
        if (Jac[22 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (Jac[22 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (Jac[22 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[22 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 22;
      for (int l = 0; l < 22; l++) {
        for (int k = 0; k < 22; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[22 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 22; k++) {
    for (int l = 0; l < 22; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 22 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[22 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 22 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 484> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 21> conc = {0.0};
  for (int n = 0; n < 21; n++) {
    conc[n] = 1.0 / 21.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 22; l++) {
      for (int k = 0; k < 22; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[22 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 22; l++) {
      for (int k = 0; k < 22; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[22 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
