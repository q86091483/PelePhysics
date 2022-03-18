#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

const int rmap[177] =
     {49,51,53,55,56,58,62,69,70,71,73,75,82,84,94,130,139,146,157,
      173,0,1,11,32,38,42,166,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,
      18,19,20,21,22,23,24,25,26,27,28,29,30,31,33,34,35,36,37,39,
      40,41,43,44,45,46,47,48,50,52,54,57,59,60,61,63,64,65,66,67,
      68,72,74,76,77,78,79,80,81,83,85,86,87,88,89,90,91,92,93,95,
      96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,
      112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,
      127,128,129,131,132,133,134,135,136,137,138,140,141,142,143,
      144,145,147,148,149,150,151,152,153,154,155,156,158,159,160,
      161,162,163,164,165,167,168,169,170,171,172,174,175,176};

/*Returns 0-based map of reaction order */
void GET_RMAP(int * _rmap)
{
    for (int j=0; j<177; ++j) {
        _rmap[j] = rmap[j];
    }
}

/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    const int ns[177] =
     {3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,3,3,3,2,3,2,3,3,3,2,3,3,4,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,
      4,3,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
      4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
      3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,3,4,4,4,5,4,4,
      4,4,3,4,4,4,4,4,4,4,3,4,4,4,4,4,4,5,4,4,4,4,4,4,4,4,3};
    const int kiv[885] =
     {1,10,12,0,0,1,12,13,0,0,1,16,17,0,0,1,17,18,0,0,1,17,19,0,0,
      1,18,20,0,0,1,19,20,0,0,1,21,22,0,0,1,22,23,0,0,1,23,24,0,0,
      1,24,25,0,0,1,25,26,0,0,0,14,17,0,0,4,7,0,0,0,4,12,20,0,0,9,
      14,27,0,0,10,14,28,0,0,11,5,20,0,0,12,26,0,0,0,24,0,22,0,0,2,
      3,0,0,0,2,1,4,0,0,2,14,15,0,0,1,3,6,0,0,1,0,0,0,0,1,4,5,0,0,
      16,1,14,0,0,2,0,1,4,0,2,6,4,3,0,2,7,4,6,0,2,9,1,14,0,2,10,1,
      16,0,2,11,0,14,0,2,11,1,16,0,2,12,1,17,0,2,13,4,12,0,2,16,4,
      14,0,2,16,1,15,0,2,17,4,16,0,2,18,4,17,0,2,19,4,17,0,2,20,4,
      18,0,2,20,4,19,0,2,21,9,14,0,2,22,1,27,0,2,22,4,21,0,2,22,14,
      10,0,2,23,1,28,0,2,24,12,16,0,2,25,12,17,0,2,26,4,25,0,2,27,
      1,14,0,2,28,4,27,0,2,28,10,15,0,3,14,2,15,0,3,17,6,16,0,1,3,
      6,3,0,1,3,5,6,5,1,3,30,6,30,1,3,31,6,31,1,3,2,4,0,1,0,0,0,0,
      1,5,0,5,0,1,15,0,15,0,1,6,2,5,0,1,6,3,0,0,1,6,4,0,0,1,7,6,0,
      0,1,7,4,5,0,1,9,8,0,0,1,11,9,0,0,1,13,12,0,0,1,16,0,14,0,1,
      17,16,0,0,1,18,0,17,0,1,18,4,12,0,1,18,11,5,0,1,19,1,18,0,1,
      19,0,17,0,1,19,4,12,0,1,19,11,5,0,1,20,18,0,0,1,20,19,0,0,1,
      23,0,22,0,1,24,23,0,0,1,25,0,24,0,1,26,25,0,0,1,27,11,14,0,1,
      28,27,0,0,1,28,12,14,0,1,29,1,28,0,4,0,1,5,0,4,2,5,0,0,4,6,3,
      5,0,4,7,6,5,0,4,7,6,5,0,4,8,1,14,0,4,9,1,16,0,4,10,1,17,0,4,
      10,9,5,0,4,11,1,17,0,4,12,10,5,0,4,12,11,5,0,4,13,12,5,0,4,
      14,1,15,0,4,16,5,14,0,4,17,16,5,0,4,18,5,17,0,4,19,5,17,0,4,
      20,18,5,0,4,20,19,5,0,4,21,1,27,0,4,22,1,28,0,4,22,1,29,0,4,
      22,21,5,0,4,22,12,14,0,4,23,5,22,0,4,24,23,5,0,4,26,25,5,0,4,
      28,27,5,0,6,3,7,0,0,6,3,7,0,0,6,10,4,17,0,6,12,3,13,0,6,12,4,
      19,0,6,14,4,15,0,6,17,16,7,0,8,3,2,14,0,8,10,1,21,0,8,12,1,
      22,0,9,3,2,16,0,9,0,1,10,0,9,5,1,17,0,9,10,1,22,0,9,12,1,23,
      0,9,13,1,24,0,9,15,16,14,0,9,17,1,28,0,9,27,14,22,0,10,3,4,
      16,0,10,0,1,12,0,10,0,22,0,0,10,12,1,24,0,10,13,12,0,0,10,27,
      23,14,0,11,30,10,30,0,11,31,10,31,0,11,3,1,4,14,11,3,14,5,0,
      11,0,12,1,0,11,5,10,5,0,11,12,1,24,0,11,13,12,0,0,11,14,10,
      14,0,11,15,10,15,0,11,15,14,17,0,11,26,12,25,0,12,3,2,19,0,
      12,3,4,17,0,12,7,6,13,0,12,1,25,0,0,12,16,13,14,0,12,17,16,
      13,0,12,20,18,13,0,12,20,19,13,0,12,24,23,13,0,12,26,25,13,0,
      16,5,1,14,5,16,3,6,14,0,18,3,6,17,0,19,3,6,17,0,21,3,16,14,0,
      21,0,1,22,0,23,3,16,17,0,25,3,6,24,0,27,3,4,14,0,27,14,22,0,
      0};
    const int nuv[885] =
     {-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,
      -1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,
      -1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,1,0,0,
      -1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,1,1,0,0,-2,
      1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,
      1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-2,1,1,0,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,-1,
      1,1,0,-2,-1,2,0,0,-2,-1,1,1,0,-2,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,
      0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -2,1,1,0,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,
      -2,2,1,0,0};
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 5;
    } else {
        if (*i > 177) {
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
    awt[4] = 39.948000; /*AR */
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
    for (id = 0; id < kd * 32; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*C */
    ncf[ 8 * kd + 2 ] = 1; /*C */

    /*CH */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 16 * kd + 1 ] = 1; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 3; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 20 * kd + 2 ] = 1; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */
    ncf[ 20 * kd + 0 ] = 1; /*O */

    /*C2H */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 1; /*H */

    /*C2H2 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 2 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 24 * kd + 2 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 25 * kd + 2 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 26 * kd + 2 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */

    /*HCCO */
    ncf[ 27 * kd + 1 ] = 1; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*C */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*CH2CO */
    ncf[ 28 * kd + 2 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 2; /*H */
    ncf[ 28 * kd + 0 ] = 1; /*O */

    /*HCCOH */
    ncf[ 29 * kd + 2 ] = 2; /*C */
    ncf[ 29 * kd + 0 ] = 1; /*O */
    ncf[ 29 * kd + 1 ] = 2; /*H */

    /*N2 */
    ncf[ 30 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 31 * kd + 4 ] = 1; /*AR */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(32);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "C";
    kname[9] = "CH";
    kname[10] = "CH2";
    kname[11] = "CH2(S)";
    kname[12] = "CH3";
    kname[13] = "CH4";
    kname[14] = "CO";
    kname[15] = "CO2";
    kname[16] = "HCO";
    kname[17] = "CH2O";
    kname[18] = "CH2OH";
    kname[19] = "CH3O";
    kname[20] = "CH3OH";
    kname[21] = "C2H";
    kname[22] = "C2H2";
    kname[23] = "C2H3";
    kname[24] = "C2H4";
    kname[25] = "C2H5";
    kname[26] = "C2H6";
    kname[27] = "HCCO";
    kname[28] = "CH2CO";
    kname[29] = "HCCOH";
    kname[30] = "N2";
    kname[31] = "AR";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(Jac[ 33 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 33 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 33 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 33;
        int offset_col = nc * 33;
        for (int k=0; k<33; k++) {
            for (int l=0; l<33; l++) {
                if(Jac[33*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if(Jac[33*k + l] != 0.0) {
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
            int offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if(Jac[33*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[33*k + l] != 0.0) {
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
            int offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[33*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 33*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[33*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 33*k + l;
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
    amrex::GpuArray<amrex::Real,1089> Jac = {0.0};
    amrex::GpuArray<amrex::Real,32> conc = {0.0};
    for (int n=0; n<32; n++) {
        conc[n] = 1.0/ 32.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<33; l++) {
            for (int k=0; k<33; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[33*k + l] != 0.0) {
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
        for (int l=0; l<33; l++) {
            for (int k=0; k<33; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[33*k + l] != 0.0) {
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
