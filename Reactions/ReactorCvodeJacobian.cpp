#include "ReactorCvodeJacobian.H"
#include "PeleLMeX_Index.H"
#include "mechanism.H"

namespace pele::physics::reactions::cvode {
#ifdef AMREX_USE_GPU
int
cJac(
  amrex::Real /*t*/,
  N_Vector y_in,
  N_Vector /*fy*/,
  SUNMatrix J,
  void* user_data,
  N_Vector /*tmp1*/,
  N_Vector /*tmp2*/,
  N_Vector /*tmp3*/)
{
  BL_PROFILE("Pele::ReactorCvode::cJac()");
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  auto solveType = udata->solve_type;
  auto ncells = udata->ncells;
  auto NNZ = udata->NNZ;
  auto stream = udata->stream;
  auto nbThreads = udata->nbThreads;
  auto nbBlocks = udata->nbBlocks;
  auto react_type = udata->reactor_type;

  if (solveType == sparseDirect) {
#ifdef AMREX_USE_CUDA
    amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y_in);
    amrex::Real* Jdata = SUNMatrix_cuSparse_Data(J);
    int* csr_row_count_d = SUNMatrix_cuSparse_IndexPointers(J);
    int* csr_col_index_d = SUNMatrix_cuSparse_IndexValues(J);

    // Checks
    AMREX_ASSERT(
      (SUNMatrix_cuSparse_Rows(J) == (NUM_SPECIES + 1) * ncells) &&
      (SUNMatrix_cuSparse_Columns(J) == (NUM_SPECIES + 1) * ncells) &&
      (SUNMatrix_cuSparse_NNZ(J) == ncells * NNZ));

    const auto ec = amrex::Gpu::ExecutionConfig(ncells);

    AMREX_ALWAYS_ASSERT(nbThreads == CVODE_NB_THREADS);
    amrex::launch_global<CVODE_NB_THREADS>
      <<<nbBlocks, CVODE_NB_THREADS, ec.sharedMem, stream>>>(
        [=] AMREX_GPU_DEVICE() noexcept {
          for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
                   stride = blockDim.x * gridDim.x;
               icell < ncells; icell += stride) {
            fKernelComputeAJchem(
              icell, NNZ, react_type, csr_row_count_d, csr_col_index_d, yvec_d,
              Jdata);
          }
        });
    amrex::Gpu::Device::streamSynchronize();
#else
    amrex::Abort(
      "Calling cJac with solve_type = sparse_direct only works with CUDA !");
#endif
  } else if (solveType == magmaDirect) {
#ifdef PELE_USE_MAGMA
    amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y_in);
    amrex::Real* Jdata = SUNMatrix_MagmaDense_Data(J);
    const auto ec = amrex::Gpu::ExecutionConfig(ncells);
    AMREX_ALWAYS_ASSERT(nbThreads == CVODE_NB_THREADS);
    amrex::launch_global<CVODE_NB_THREADS>
      <<<nbBlocks, CVODE_NB_THREADS, ec.sharedMem, stream>>>(
        [=] AMREX_GPU_DEVICE() noexcept {
          for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
                   stride = blockDim.x * gridDim.x;
               icell < ncells; icell += stride) {
            fKernelDenseAJchem(icell, react_type, yvec_d, Jdata);
          }
        });
    amrex::Gpu::Device::streamSynchronize();
#else
    amrex::Abort(
      "Calling cJac with solve_type = magma_direct requires PELE_USE_MAGMA = "
      "TRUE !");
#endif
  }

  return (0);
}

#else

int
cJac(
  amrex::Real /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{
  BL_PROFILE("Pele::ReactorCvode::cJacDense()");

  // Make local copies of pointers to input data
  amrex::Real* ydata = N_VGetArrayPointer(u);

  // Make local copies of pointers in user_data
  auto* udata = static_cast<CVODEUserData*>(user_data);
  auto ncells = udata->ncells;
  auto reactor_type = udata->reactor_type;
#if defined (PELE_USE_AUX) && (NUMUDA > 0)
  auto* rhoAuxsrc_ext = udata->rhoAuxsrc_ext;
  auto* rhoAux_init = udata->rhoAux_init;
#endif

  for (int tid = 0; tid < ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1 + NUMODE);

    // rho MKS
    amrex::Real rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }

    amrex::Real temp = ydata[offset + NUM_SPECIES];

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    // Yks
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] / rho;
    }

    // Jac
    amrex::Real Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
    const int consP =
      static_cast<int>(reactor_type == ReactorTypes::h_reactor_type);
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);

    // fill the sunMat and scale
    for (int i = 0; i < NUM_SPECIES; i++) {
      // cppcheck-suppress cstyleCast
      amrex::Real* J_col = SM_COLUMN_D(J, offset + i);
      for (int k = 0; k < NUM_SPECIES; k++) {
        J_col[offset + k] =
          Jmat_tmp[i * (NUM_SPECIES + 1) + k] * mw(k) * imw(i);
      }
      J_col[offset + NUM_SPECIES] =
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] * imw(i);
    }
    // cppcheck-suppress cstyleCast
    amrex::Real* J_col = SM_COLUMN_D(J, offset + NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++) {
      J_col[offset + i] = Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] * mw(i);
    }
    // J_col = SM_COLUMN_D(J, offset); // Never read

#if defined (PELE_USE_AUX) && (NUMAUX > 0)
    /*STAR
    for (int ii = 0; ii < NUM_SPECIES + 1 + NUMAUX ; ii++) {
      amrex::Print() << "Line ii = " << ii << ", offset = "  << offset << std::endl;
      amrex::Print() << "  ";
      for (int kk = offset; kk < offset + NUM_SPECIES + 1 + NUMAUX; kk++) {
        amrex::Real* J_col_temp = SM_COLUMN_D(J, kk);
        const int ROW_J = offset + ii;
        amrex::Print() << J_col_temp[ROW_J] << ", ";
      }
      amrex::Print() << std::endl;
    }
    */
#if (NUMAGE > 0)
    /*
    for (int i = 0; i < NUMAGE; i++) {
      const int MIXF_IN_J = offset + NUM_SPECIES + 1 + MIXF_IN_AUX + i;
      const int AGE_IN_J  = offset + NUM_SPECIES + 1 + AGE_IN_AUX  + i;
      amrex::Real* J_col_mixf = SM_COLUMN_D(J, MIXF_IN_J);
      amrex::Real* J_col_age  = SM_COLUMN_D(J, AGE_IN_J);
      J_col_mixf[AGE_IN_J] = 1 - (ydata[AGE_IN_J]/ydata[MIXF_IN_J]/ydata[MIXF_IN_J]) * rhoAuxsrc_ext[MIXF_IN_AUX+i];
      J_col_age[AGE_IN_J] = rhoAuxsrc_ext[MIXF_IN_AUX+i] / ydata[MIXF_IN_J];
    }
    */
#endif // #if (NUMAGE > 0)
#if (NUMAGEPV > 0)
    const amrex::Real Tc      = 1750.;
    const amrex::Real dT      = 10.;
    const amrex::Real tanh_T  = std::tanh((temp - Tc) / dT);
    const amrex::Real f_T     = 0.5 * (1 + tanh_T);
    const amrex::Real dfdT    = 0.5 * (1 - tanh_T * tanh_T) / dT;
    for (int i = 0; i < NUMAGEPV; i++) {

      const int TEMP_IN_J   = offset + NUM_SPECIES;
      const int MIXF_IN_J   = offset + NUM_SPECIES + 1 + MIXF_IN_AUX  + i;
      const int AGE_IN_J    = offset + NUM_SPECIES + 1 + AGE_IN_AUX   + i;
      const int AGEPV_IN_J  = offset + NUM_SPECIES + 1 + AGEPV_IN_AUX + i;

      amrex::Real* J_col_temp   = SM_COLUMN_D(J,  TEMP_IN_J);
      amrex::Real* J_col_mixf   = SM_COLUMN_D(J,  MIXF_IN_J);
      amrex::Real* J_col_age    = SM_COLUMN_D(J,   AGE_IN_J);
      amrex::Real* J_col_agepv  = SM_COLUMN_D(J, AGEPV_IN_J);

      J_col_temp[AGEPV_IN_J]  = ydata[MIXF_IN_J] * dfdT;
      J_col_mixf[AGEPV_IN_J]  = f_T - (ydata[AGEPV_IN_J]/ydata[MIXF_IN_J]/ydata[MIXF_IN_J]) * rhoAuxsrc_ext[MIXF_IN_AUX+i];
      J_col_age[AGEPV_IN_J]   = 0.0;
      J_col_agepv[AGEPV_IN_J] = rhoAuxsrc_ext[MIXF_IN_AUX+i] / ydata[MIXF_IN_J];
    }
#endif // #if (NUMAGEPV > 0)
#endif // #if (NUMAUX > 0)
  }

  return (0);
}

#if defined (PELE_USE_AUX) && (NUMAUX > 0)
int
cJac_aux(
  amrex::Real /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{
  BL_PROFILE("Pele::ReactorCvode::cJacDense()");

  // Make local copies of pointers to input data
  amrex::Real* ydata = N_VGetArrayPointer(u);
  // Make local copies of pointers in user_data
  auto* udata = static_cast<CVODEUserData*>(user_data);
  auto ncells = udata->ncells;
  auto reactor_type = udata->reactor_type;
#if defined (PELE_USE_AUX) && (NUMAUX > 0)
  auto* rhoAuxsrc_ext = udata->rhoAuxsrc_ext;
  auto* rhoAux_init = udata->rhoAux_init;
  auto* rhoY_T_init = udata->rhoY_T_init;
#endif

  for (int tid = 0; tid < ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUMAUX);

    /*
    // rho MKS
    amrex::Real rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }

    amrex::Real temp = ydata[offset + NUM_SPECIES];

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    // Yks
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] / rho;
    }

    // Jac
    amrex::Real Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
    const int consP =
      static_cast<int>(reactor_type == ReactorTypes::h_reactor_type);
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);

    // fill the sunMat and scale
    for (int i = 0; i < NUM_SPECIES; i++) {
      // cppcheck-suppress cstyleCast
      amrex::Real* J_col = SM_COLUMN_D(J, offset + i);
      for (int k = 0; k < NUM_SPECIES; k++) {
        J_col[offset + k] =
          Jmat_tmp[i * (NUM_SPECIES + 1) + k] * mw(k) * imw(i);
      }
      J_col[offset + NUM_SPECIES] =
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] * imw(i);
    }
    // cppcheck-suppress cstyleCast
    amrex::Real* J_col = SM_COLUMN_D(J, offset + NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++) {
      J_col[offset + i] = Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] * mw(i);
    }
    // J_col = SM_COLUMN_D(J, offset); // Never read
    */

#if defined (PELE_USE_AUX) && (NUMAUX > 0)
    //for (int ii = 0; ii < NUM_SPECIES + 1 + NUMAUX ; ii++) {
    //  amrex::Print() << "Line ii = " << ii << ", offset = "  << offset << std::endl;
    //  amrex::Print() << "  ";
    //  for (int kk = offset; kk < offset + NUM_SPECIES + 1 + NUMAUX; kk++) {
    //    amrex::Real* J_col_temp = SM_COLUMN_D(J, kk);
    //    const int ROW_J = offset + ii;
    //    amrex::Print() << J_col_temp[ROW_J] << ", ";
    //  }
    //  amrex::Print() << std::endl;
    //}

#if (NUMAGE > 0)
    for (int i = 0; i < NUMAGE; i++) {
      //const int MIXF_IN_J = offset + NUM_SPECIES + 1 + MIXF_IN_AUX + i;
      //const int AGE_IN_J  = offset + NUM_SPECIES + 1 + AGE_IN_AUX  + i;
      const int MIXF_IN_J = offset + MIXF_IN_AUX + i;
      const int AGE_IN_J  = offset + AGE_IN_AUX  + i;
      amrex::Real* J_col_mixf = SM_COLUMN_D(J, MIXF_IN_J);
      amrex::Real* J_col_age  = SM_COLUMN_D(J, AGE_IN_J);
      J_col_mixf[AGE_IN_J] = 1 - (ydata[AGE_IN_J]/ydata[MIXF_IN_J]/ydata[MIXF_IN_J]) * rhoAuxsrc_ext[MIXF_IN_AUX+i];
      J_col_age[AGE_IN_J] = rhoAuxsrc_ext[MIXF_IN_AUX+i] / ydata[MIXF_IN_J];
    }
#endif // #if (NUMAGE > 0)
#if (NUMAGEPV > 0)
    const amrex::Real Tc      = 1750.;
    const amrex::Real dT      = 10.;
    const amrex::Real tanh_T  = std::tanh((temp - Tc) / dT);
    const amrex::Real f_T     = 0.5 * (1 + tanh_T);
    const amrex::Real dfdT    = 0.5 * (1 - tanh_T * tanh_T) / dT;
    for (int i = 0; i < NUMAGEPV; i++) {

      const int TEMP_IN_J   = offset + NUM_SPECIES;
      const int MIXF_IN_J   = offset + NUM_SPECIES + 1 + MIXF_IN_AUX  + i;
      const int AGE_IN_J    = offset + NUM_SPECIES + 1 + AGE_IN_AUX   + i;
      const int AGEPV_IN_J  = offset + NUM_SPECIES + 1 + AGEPV_IN_AUX + i;

      amrex::Real* J_col_temp   = SM_COLUMN_D(J,  TEMP_IN_J);
      amrex::Real* J_col_mixf   = SM_COLUMN_D(J,  MIXF_IN_J);
      amrex::Real* J_col_age    = SM_COLUMN_D(J,   AGE_IN_J);
      amrex::Real* J_col_agepv  = SM_COLUMN_D(J, AGEPV_IN_J);

      J_col_temp[AGEPV_IN_J]  = ydata[MIXF_IN_J] * dfdT;
      J_col_mixf[AGEPV_IN_J]  = f_T - (ydata[AGEPV_IN_J]/ydata[MIXF_IN_J]/ydata[MIXF_IN_J]) * rhoAuxsrc_ext[MIXF_IN_AUX+i];
      J_col_age[AGEPV_IN_J]   = 0.0;
      J_col_agepv[AGEPV_IN_J] = rhoAuxsrc_ext[MIXF_IN_AUX+i] / ydata[MIXF_IN_J];
    }
#endif // #if (NUMAGEPV > 0)
#endif // #if (NUMAUX > 0)
  }
  return (0);
}
#endif // #if (NUMNEW > 0) cJac_aux

// Analytical SPARSE CSR Jacobian evaluation
int
cJac_sps(
  amrex::Real /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{
  BL_PROFILE("Pele::ReactorCvode::cJacSparse()");
  // Make local copies of pointers to input data
  amrex::Real* ydata = N_VGetArrayPointer(u);

  // Make local copies of pointers in user_data (cell M)*/
  auto* udata = static_cast<CVODEUserData*>(user_data);
  auto NNZ = udata->NNZ;
  auto reactor_type = udata->reactor_type;
  auto ncells = udata->ncells;
  auto* colVals_c = udata->colVals_c;
  auto* rowPtrs_c = udata->rowPtrs_c;

  sunindextype* rowPtrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype* colIndx_tmp = SUNSparseMatrix_IndexValues(J);
  amrex::Real* Jdata = SUNSparseMatrix_Data(J);
  // Fixed colVal
  for (int i = 0; i < NNZ * ncells; i++) {
    colIndx_tmp[i] = (sunindextype)colVals_c[i];
  }
  rowPtrs_tmp[0] = (sunindextype)rowPtrs_c[0];
  // Fixed rowPtrs
  for (int i = 0; i < ncells * (NUM_SPECIES + 1); i++) {
    rowPtrs_tmp[i + 1] = (sunindextype)rowPtrs_c[i + 1];
  }

  // Temp vectors
  // Save Jac from cell to cell if more than one
  amrex::Real temp_save_lcl = 0.0;
  for (int tid = 0; tid < ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);
    int offset_J = tid * NNZ;
    // rho MKS
    amrex::Real rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }
    // Yks
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real rhoinv = 1.0 / rho;
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] * rhoinv;
    }
    amrex::Real temp = ydata[offset + NUM_SPECIES];

    // Do we recompute Jac ?
    amrex::Real Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
    if (fabs(temp - temp_save_lcl) > 1.0) {
      const int consP =
        static_cast<int>(reactor_type == ReactorTypes::h_reactor_type);
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
      temp_save_lcl = temp;
      // rescale
      for (int i = 0; i < NUM_SPECIES; i++) {
        for (int k = 0; k < NUM_SPECIES; k++) {
          Jmat_tmp[k * (NUM_SPECIES + 1) + i] *= mw(i) * imw(k);
        }
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] *= imw(i);
      }
      for (int i = 0; i < NUM_SPECIES; i++) {
        Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] *= mw(i);
      }
    }
    // Go from Dense to Sparse
    for (int i = 1; i < NUM_SPECIES + 2; i++) {
      int nbVals = rowPtrs_c[i] - rowPtrs_c[i - 1];
      for (int j = 0; j < nbVals; j++) {
        int idx = colVals_c[rowPtrs_c[i - 1] + j];
        Jdata[offset_J + rowPtrs_c[i - 1] + j] =
          Jmat_tmp[(i - 1) + (NUM_SPECIES + 1) * idx];
      }
    }
  }

  return (0);
}

#ifdef PELE_USE_KLU
// Analytical SPARSE KLU CSC Jacobian evaluation
int
cJac_KLU(
  amrex::Real /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{
  BL_PROFILE("Pele::ReactorCvode::cJacSparseKLU()");

  // Make local copies of pointers to input data
  amrex::Real* ydata = N_VGetArrayPointer(u);

  // Make local copies of pointers in user_data (cell M)
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  auto NNZ = udata->NNZ;
  auto reactor_type = udata->reactor_type;
  auto ncells = udata->ncells;
  auto colPtrs = udata->colPtrs;
  auto rowVals = udata->rowVals;

  // Fixed RowVals
  sunindextype* colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype* rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  amrex::Real* Jdata = SUNSparseMatrix_Data(J);
  for (int i = 0; i < NNZ; i++) {
    rowvals_tmp[i] = rowVals[0][i];
  }
  // Fixed colPtrs
  colptrs_tmp[0] = colPtrs[0][0];
  for (int i = 0; i < ncells * (NUM_SPECIES + 1); i++) {
    colptrs_tmp[i + 1] = colPtrs[0][i + 1];
  }

  // Save Jac from cell to cell if more than one
  amrex::Real temp_save_lcl = 0.0;
  for (int tid = 0; tid < ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);
    // rho
    amrex::Real rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }
    // Yks
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real rhoinv = 1.0 / rho;
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] * rhoinv;
    }
    amrex::Real temp = ydata[offset + NUM_SPECIES];

    // Do we recompute Jac ?
    amrex::Real Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
    if (fabs(temp - temp_save_lcl) > 1.0) {
      const int consP = reactor_type == ReactorTypes::h_reactor_type;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
      temp_save_lcl = temp;
      // rescale
      for (int i = 0; i < NUM_SPECIES; i++) {
        for (int k = 0; k < NUM_SPECIES; k++) {
          Jmat_tmp[k * (NUM_SPECIES + 1) + i] *= mw(i) * imw(k);
        }
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] *= imw(i);
      }
      for (int i = 0; i < NUM_SPECIES; i++) {
        Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] *= mw(i);
      }
    }
    // Go from Dense to Sparse
    BL_PROFILE_VAR("DensetoSps", DtoS);
    for (int i = 1; i < NUM_SPECIES + 2; i++) {
      int nbVals = colPtrs[0][i] - colPtrs[0][i - 1];
      for (int j = 0; j < nbVals; j++) {
        int idx = rowVals[0][colPtrs[0][i - 1] + j];
        Jdata[colPtrs[0][offset + i - 1] + j] =
          Jmat_tmp[(i - 1) * (NUM_SPECIES + 1) + idx];
      }
    }
    BL_PROFILE_VAR_STOP(DtoS);
  }

  return (0);
}
#endif
#endif
} // namespace pele::physics::reactions::cvode
