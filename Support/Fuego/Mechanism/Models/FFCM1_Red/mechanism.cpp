#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

const int rmap[124] =
     {12,33,40,58,64,65,72,75,85,86,96,97,19,5,7,9,10,24,0,1,2,3,4,
      6,8,11,13,14,15,16,17,18,20,21,22,23,25,26,27,28,29,30,31,32,
      34,35,36,37,38,39,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
      55,56,57,59,60,61,62,63,66,67,68,69,70,71,73,74,76,77,78,79,
      80,81,82,83,84,87,88,89,90,91,92,93,94,95,98,99,100,101,102,
      103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,
      118,119,120,121,122,123};

/*Returns 0-based map of reaction order */
void GET_RMAP(int * _rmap)
{
    for (int j=0; j<124; ++j) {
        _rmap[j] = rmap[j];
    }
}

/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    const int ns[124] =
     {3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,4,4,4,4,3,4,4,4,4,3,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4,4,5,4,4,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4,4,4,4,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
      4,4,4,4};
    const int kiv[620] =
     {2,4,7,0,0,10,1,13,0,0,11,2,13,0,0,12,6,19,0,0,15,2,16,0,0,16,
      1,8,0,0,13,2,14,0,0,13,5,19,0,0,18,2,16,0,0,18,2,19,0,0,17,2,
      16,0,0,17,2,19,0,0,8,3,9,0,0,1,2,0,0,0,3,4,0,0,0,3,2,5,0,0,6,
      2,5,0,0,15,2,8,0,0,2,4,3,5,0,3,1,2,5,0,3,1,2,5,0,5,1,2,6,0,5,
      3,6,0,0,1,20,2,20,0,3,20,4,20,0,6,2,5,6,0,7,2,1,4,0,7,2,5,0,
      0,7,2,3,6,0,7,3,5,4,0,7,5,6,4,0,7,5,6,4,0,8,4,3,9,0,8,5,2,9,
      0,8,5,2,9,0,8,7,5,9,0,15,2,1,8,0,15,3,5,8,0,15,3,2,9,0,15,5,
      6,8,0,15,4,7,8,0,10,3,2,8,0,10,5,2,15,0,10,1,2,11,0,10,6,2,
      16,0,10,4,3,15,0,10,4,9,2,0,10,4,8,5,0,10,4,3,2,8,10,9,15,8,
      0,11,3,2,8,0,11,5,2,16,0,11,5,10,6,0,11,7,5,16,0,11,1,2,13,0,
      11,4,5,2,8,11,4,2,9,0,11,4,3,16,0,11,4,1,9,0,11,4,6,8,0,12,0,
      11,0,0,12,20,11,20,0,12,2,10,1,0,12,3,2,8,0,12,5,2,16,0,12,1,
      13,2,0,12,4,11,4,0,12,6,11,6,0,12,6,1,16,0,12,8,11,8,0,12,9,
      11,9,0,12,9,8,16,0,16,2,15,1,0,16,3,5,15,0,16,5,15,6,0,16,4,
      7,15,0,16,11,13,15,0,16,12,13,15,0,13,3,2,16,0,13,3,2,1,8,13,
      5,11,6,0,13,5,12,6,0,13,5,1,16,0,13,7,4,14,0,13,7,5,18,0,13,
      4,3,18,0,13,4,5,16,0,13,15,14,8,0,13,16,15,14,0,18,2,2,17,0,
      18,2,1,16,0,18,2,5,13,0,18,2,12,6,0,18,3,5,16,0,18,5,6,16,0,
      18,4,7,16,0,18,13,14,16,0,18,8,13,9,0,17,2,1,16,0,17,2,5,13,
      0,17,2,12,6,0,17,3,5,16,0,17,5,6,16,0,17,4,7,16,0,17,13,14,
      16,0,14,2,13,1,0,14,3,5,13,0,14,5,13,6,0,14,11,13,0,0,14,12,
      13,0,0,19,2,17,1,0,19,2,18,1,0,19,3,5,17,0,19,3,5,18,0,19,5,
      17,6,0,19,5,18,6,0,19,4,17,7,0,19,10,13,16,0,19,11,13,17,0,
      19,11,13,18,0,19,12,13,18,0,19,12,13,17,0,19,13,17,14,0,19,
      13,18,14,0};
    const int nuv[620] =
     {-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,
      -1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,
      1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,2,0,0,0,-2,1,0,0,0,-1,-1,
      1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,2,1,0,-2,-1,1,1,0,-2,1,1,1,
      0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,2,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      1,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,2,0,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0};
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 5;
    } else {
        if (*i > 124) {
            *nspec = -1;
        } else {
            *nspec = ns[*i-1];
            for (int j=0; j<*nspec; ++j) {
                ki[j] = kiv[(*i-1)*5 + j] + 1;
                nu[j] = nuv[(*i-1)*5 + j];
            }
        }
    }
}


/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 4.002600; /*HE */
}



/*get atomic weight for all elements */
void CKAWT( amrex::Real *  awt)
{
    atomicWeight(awt);
}



/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 5; 
    /*Zero ncf */
    for (id = 0; id < kd * 21; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 3 ] = 2; /*N */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 3 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 5 * kd + 0 ] = 1; /*O */
    ncf[ 5 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 7 * kd + 1 ] = 1; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*CO */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 2; /*O */

    /*CH */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 1; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 4; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*HE */
    ncf[ 20 * kd + 4 ] = 1; /*HE */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "HE";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(21);
    kname[0] = "N2";
    kname[1] = "H2";
    kname[2] = "H";
    kname[3] = "O";
    kname[4] = "O2";
    kname[5] = "OH";
    kname[6] = "H2O";
    kname[7] = "HO2";
    kname[8] = "CO";
    kname[9] = "CO2";
    kname[10] = "CH";
    kname[11] = "CH2";
    kname[12] = "CH2(S)";
    kname[13] = "CH3";
    kname[14] = "CH4";
    kname[15] = "HCO";
    kname[16] = "CH2O";
    kname[17] = "CH2OH";
    kname[18] = "CH3O";
    kname[19] = "CH3OH";
    kname[20] = "HE";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(Jac[ 22 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 22 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 22 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 22;
        int offset_col = nc * 22;
        for (int k=0; k<22; k++) {
            for (int l=0; l<22; l++) {
                if(Jac[22*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(Jac[22*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtrs[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(Jac[22*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[22*k + l] != 0.0) {
                            colVals[nJdata_tmp-1] = k+1 + offset; 
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
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[22*k + l] != 0.0) {
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

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 22*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[22*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 22*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
    amrex::GpuArray<amrex::Real,484> Jac = {0.0};
    amrex::GpuArray<amrex::Real,21> conc = {0.0};
    for (int n=0; n<21; n++) {
        conc[n] = 1.0/ 21.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[22*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[22*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }
}

#endif
