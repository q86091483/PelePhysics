#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

const int rmap[387] =
     {0,14,21,23,48,55,57,24,37,38,81,1,2,3,4,5,6,7,8,9,10,11,12,
      13,15,16,17,18,19,20,22,25,26,27,28,29,30,31,32,33,34,35,36,
      39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,56,58,59,60,61,
      62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,82,
      83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,
      102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
      117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,
      132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,
      147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,
      162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,
      177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
      192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,
      207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,
      222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,
      237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,
      252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,
      267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,
      282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,
      297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,
      312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,
      327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,
      342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,
      357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,
      372,373,374,375,376,377,378,379,380,381,382,383,384,385,386};

/*Returns 0-based map of reaction order */
void GET_RMAP(int * _rmap)
{
    for (int j=0; j<387; ++j) {
        _rmap[j] = rmap[j];
    }
}

/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    const int ns[387] =
     {3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,
      3,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4,4,4,3,4,3,4,4,4,4,4,4,3,
      3,4,3,5,4,4,4,3,4,4,4,4,4,4,3,3,3,4,4,4,4,3,3,3,3,4,4,4,4,4,
      4,4,4,3,3,3,4,4,4,4,4,4,3,3,4,3,3,4,4,4,4,4,4,3,4,4,4,3,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,3,3,5,4,4,4,4,4,4,4,4,5,
      4,3,4,3,4,4,4,4,4,4,3,3,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,
      4,3,3,4,2,4,4,5,4,4,4,4,5,4,4,4,4,4,4,4,3,4,4,4,4,4,3,4,4,4,
      4,4,3,3,5,4,4,4,5,4,4,4,4,4,4,4,4,4,4,3,3,3,4,4,4,4,4,4,4,3,
      4,4,5,4,5,4,4,4,4,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,
      3,3,3,4,4,4,4,2,2,2,2,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,2,2,3,
      3,3,3,3,3,3,4,4,4,3,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5};
    const int kiv[1935] =
     {9,0,10,0,0,0,16,17,0,0,9,9,15,0,0,0,3,11,0,0,4,4,12,0,0,19,0,
      16,0,0,18,0,19,0,0,6,2,8,0,0,14,13,0,0,0,16,18,1,0,0,22,9,6,
      0,0,10,0,9,1,0,10,4,9,5,0,10,2,9,4,0,15,9,17,10,0,7,4,6,5,0,
      6,4,8,0,0,0,3,2,4,0,2,1,0,4,0,2,5,4,4,0,4,1,0,5,0,7,0,6,0,0,
      12,4,5,11,0,16,2,9,7,0,15,0,17,1,0,17,3,16,11,0,15,4,17,5,0,
      15,2,17,4,0,9,11,14,4,0,6,11,8,4,0,5,0,4,0,0,6,3,8,2,0,7,0,6,
      1,0,7,2,6,4,0,13,7,0,0,0,13,4,7,5,0,13,0,7,1,0,13,2,7,4,0,9,
      4,13,1,0,9,2,13,0,0,9,3,14,2,0,13,9,7,10,0,7,9,10,6,0,11,2,4,
      3,0,7,11,13,3,0,14,3,13,11,0,9,11,10,3,0,7,3,6,11,0,11,0,4,4,
      0,11,0,1,3,0,11,4,5,3,0,12,3,11,11,0,12,0,5,4,0,10,11,9,12,0,
      13,11,7,12,0,4,2,0,0,0,3,2,2,0,0,1,0,0,0,0,17,19,16,16,0,16,
      0,19,1,0,16,4,19,5,0,18,3,21,4,0,19,3,18,11,0,12,2,4,11,0,18,
      2,21,0,0,18,4,20,0,0,20,0,9,6,0,20,2,21,4,0,20,4,21,5,0,20,0,
      21,1,0,21,4,7,7,0,21,0,50,6,0,21,2,0,6,6,15,3,17,11,0,15,11,
      17,12,0,9,19,10,18,0,9,17,10,16,0,19,0,18,1,0,17,0,9,9,0,19,
      3,13,7,0,15,17,0,0,0,16,9,19,10,0,24,9,7,0,0,24,3,22,11,0,24,
      4,22,5,0,24,0,22,1,0,24,2,22,4,0,24,11,22,12,0,24,9,22,10,0,
      47,18,9,0,0,26,19,9,0,0,18,9,25,0,0,26,47,0,0,0,26,2,20,9,0,
      26,2,17,7,0,26,11,47,12,0,26,4,47,5,0,27,19,19,0,0,27,4,17,
      20,0,27,4,13,47,0,27,4,19,24,0,27,2,16,20,0,27,2,13,25,0,16,
      3,19,11,0,13,6,1,0,0,28,9,16,0,0,28,0,26,0,0,28,3,26,11,0,26,
      2,47,4,0,26,0,47,1,0,26,0,16,9,0,35,47,16,0,0,35,27,9,0,0,29,
      27,0,0,0,29,16,19,0,0,30,27,29,29,0,29,9,27,10,0,47,29,26,27,
      0,29,3,27,11,0,0,29,27,1,0,17,29,27,15,0,17,29,16,30,0,19,29,
      16,27,0,30,47,9,0,0,30,19,17,0,0,30,0,29,0,0,30,9,29,10,0,30,
      0,29,1,0,30,4,29,5,0,30,47,29,26,0,30,11,29,12,0,30,3,29,11,
      0,31,17,16,0,0,31,30,0,0,0,31,3,30,11,0,32,20,9,0,0,34,17,6,
      0,0,33,0,34,1,0,33,2,34,4,0,33,4,34,5,0,33,9,34,10,0,33,11,
      34,12,0,33,17,34,15,0,33,17,7,0,0,33,3,34,11,0,33,19,34,16,0,
      33,47,34,26,0,36,17,47,0,0,36,0,35,1,0,36,2,35,4,0,36,2,31,7,
      0,36,2,28,22,0,36,4,35,5,0,36,4,31,13,0,36,4,28,24,0,36,9,35,
      10,0,12,0,1,11,0,7,2,8,0,0,22,0,20,1,0,22,2,20,4,0,22,9,20,
      10,0,16,2,23,0,0,17,2,24,0,0,19,16,27,0,0,37,9,13,0,0,37,3,
      24,11,0,12,3,11,11,0,19,3,23,2,0,38,9,3,0,0,39,14,4,0,0,49,3,
      21,6,0,38,13,39,7,0,16,38,19,39,0,10,38,9,39,0,17,11,37,4,0,
      38,9,14,14,0,38,17,14,37,0,38,11,39,3,0,12,4,5,11,0,38,38,3,
      14,14,15,38,17,39,0,37,24,0,0,0,38,24,39,22,0,40,19,6,0,0,41,
      4,40,5,0,41,0,40,1,0,41,2,40,4,0,41,11,40,12,0,41,9,40,10,0,
      41,38,40,39,0,42,41,0,0,0,42,19,13,0,0,42,3,41,11,0,47,11,42,
      4,0,47,38,42,14,0,26,38,47,39,0,9,4,50,5,0,43,24,19,0,0,43,
      41,9,0,0,29,11,43,4,0,29,38,43,14,0,30,4,28,13,0,30,2,26,13,
      0,30,2,24,16,0,30,2,22,17,0,30,4,24,17,0,30,4,22,15,0,30,2,
      34,9,0,30,4,33,9,0,30,4,34,10,0,30,38,29,39,0,44,45,3,0,0,44,
      51,4,0,0,45,4,13,26,0,46,45,0,0,0,33,38,34,39,0,33,29,34,30,
      0,25,11,16,6,4,25,11,48,12,0,26,3,47,11,0,26,9,47,10,0,26,17,
      47,15,0,47,11,19,13,4,47,0,25,1,0,47,9,25,10,0,47,17,15,25,0,
      47,17,16,26,0,47,19,16,25,0,25,26,47,47,0,47,3,41,4,0,25,48,
      0,0,0,25,3,48,11,0,48,0,49,1,0,25,4,48,5,0,25,2,16,6,0,49,4,
      18,7,0,47,25,0,0,0,25,0,48,1,0,25,9,48,10,0,25,47,48,26,0,48,
      4,49,5,0,48,3,20,7,0,46,31,3,0,0,23,20,0,0,0,23,3,13,6,4,51,
      24,23,4,0,47,3,25,11,0,47,3,23,13,0,47,3,18,13,4,21,3,8,7,0,
      9,3,13,4,0,16,1,9,9,0,52,3,53,11,0,52,4,53,5,0,52,0,53,1,0,
      52,2,53,4,0,52,11,53,12,0,52,9,53,10,0,52,38,53,39,0,53,28,6,
      0,0,54,20,17,0,0,55,28,20,0,0,56,3,57,11,0,56,4,57,5,0,56,0,
      57,1,0,56,2,57,4,0,56,11,57,12,0,56,9,57,10,0,56,38,57,39,0,
      57,31,6,0,0,50,10,9,9,0,50,15,9,17,0,50,3,6,4,0,50,1,9,0,0,
      50,2,6,0,0,50,4,13,0,0,50,8,13,6,0,50,9,16,0,0,50,20,16,6,0,
      58,0,59,0,0,58,0,60,0,0,58,0,61,0,0,58,0,62,0,0,58,31,28,0,0,
      58,0,59,1,0,58,0,60,1,0,58,0,61,1,0,58,0,62,1,0,58,2,59,4,0,
      58,2,60,4,0,58,2,61,4,0,58,2,62,4,0,58,4,59,5,0,58,4,60,5,0,
      58,4,61,5,0,58,4,62,5,0,58,11,59,12,0,58,11,60,12,0,58,11,61,
      12,0,58,11,62,12,0,58,9,59,10,0,58,9,60,10,0,58,9,61,10,0,58,
      9,62,10,0,58,3,59,11,0,58,3,60,11,0,58,3,61,11,0,58,3,62,11,
      0,58,17,59,15,0,58,17,60,15,0,58,17,61,15,0,58,17,62,15,0,58,
      19,59,16,0,58,19,60,16,0,58,19,61,16,0,58,19,62,16,0,58,38,
      59,39,0,58,38,60,39,0,58,38,61,39,0,58,38,62,39,0,58,59,60,
      58,0,58,59,61,58,0,58,59,62,58,0,58,60,61,58,0,58,60,62,58,0,
      58,61,62,58,0,60,31,26,0,0,60,63,0,0,0,61,30,28,0,0,61,63,0,
      0,0,61,64,0,0,0,62,17,36,0,0,62,64,0,0,0,60,3,63,11,0,61,3,
      63,11,0,61,3,64,11,0,62,3,64,11,0,59,61,0,0,0,59,62,0,0,0,60,
      61,0,0,0,59,60,0,0,0,63,4,33,31,0,64,4,33,31,0,63,2,24,36,0,
      64,2,24,36,0,63,29,28,0,0,64,29,28,0,0,65,59,3,0,0,66,60,3,0,
      0,67,61,3,0,0,68,62,3,0,0,65,69,0,0,0,66,70,0,0,0,66,71,0,0,
      0,67,72,0,0,0,67,73,0,0,0,67,74,0,0,0,68,75,0,0,0,68,76,0,0,
      0,70,63,11,0,0,72,63,11,0,0,73,64,11,0,0,76,64,11,0,0,69,81,
      4,0,0,71,82,4,0,0,75,82,4,0,0,71,4,24,36,0,74,4,33,30,0,75,4,
      52,26,0,77,69,3,0,0,78,71,3,0,0,79,74,3,0,0,80,75,3,0,0,77,
      83,4,0,0,78,84,4,0,0,79,85,4,0,0,80,86,4,0,0,83,56,23,4,0,84,
      52,32,4,0,85,33,54,4,0,86,24,55,4,0,82,4,22,36,5,81,4,16,57,
      5,82,4,26,53,5,82,11,22,36,12,81,11,16,57,12,82,11,26,53,12};
    const int nuv[1935] =
     {-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,
      -1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,
      1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,
      0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,
      -1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,
      -1,1,1,0,-1,1,1,0,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,
      1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,
      0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,1,
      1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,
      1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,1,0,
      -1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,
      -1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,
      -1,1,1,1,-1,1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,
      1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,
      0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,
      1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,
      0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,
      -1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,
      0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,
      0,-1,-1,1,1,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,
      -1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,
      1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,0,0,0,
      -1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,
      0,0,0,-1,1,0,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,
      -1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,
      1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,
      -1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,
      1,1,0,-1,1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,
      1,1,-1,-1,1,1,1,-1,-1,1,1,1};
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 5;
    } else {
        if (*i > 387) {
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
    awt[0] = 1.007970; /*H */
    awt[1] = 12.011150; /*C */
    awt[2] = 15.999400; /*O */
    awt[3] = 14.006700; /*N */
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
    int kd = 4; 
    /*Zero ncf */
    for (id = 0; id < kd * 88; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 0 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 2 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*H */
    ncf[ 4 * kd + 2 ] = 1; /*O */

    /*H2O */
    ncf[ 5 * kd + 0 ] = 2; /*H */
    ncf[ 5 * kd + 2 ] = 1; /*O */

    /*CO */
    ncf[ 6 * kd + 1 ] = 1; /*C */
    ncf[ 6 * kd + 2 ] = 1; /*O */

    /*HCO */
    ncf[ 7 * kd + 0 ] = 1; /*H */
    ncf[ 7 * kd + 1 ] = 1; /*C */
    ncf[ 7 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 8 * kd + 1 ] = 1; /*C */
    ncf[ 8 * kd + 2 ] = 2; /*O */

    /*CH3 */
    ncf[ 9 * kd + 1 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 3; /*H */

    /*CH4 */
    ncf[ 10 * kd + 1 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 4; /*H */

    /*HO2 */
    ncf[ 11 * kd + 0 ] = 1; /*H */
    ncf[ 11 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 12 * kd + 0 ] = 2; /*H */
    ncf[ 12 * kd + 2 ] = 2; /*O */

    /*CH2O */
    ncf[ 13 * kd + 1 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*CH3O */
    ncf[ 14 * kd + 1 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 3; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 15 * kd + 1 ] = 2; /*C */
    ncf[ 15 * kd + 0 ] = 6; /*H */

    /*C2H4 */
    ncf[ 16 * kd + 1 ] = 2; /*C */
    ncf[ 16 * kd + 0 ] = 4; /*H */

    /*C2H5 */
    ncf[ 17 * kd + 1 ] = 2; /*C */
    ncf[ 17 * kd + 0 ] = 5; /*H */

    /*C2H2 */
    ncf[ 18 * kd + 1 ] = 2; /*C */
    ncf[ 18 * kd + 0 ] = 2; /*H */

    /*C2H3 */
    ncf[ 19 * kd + 1 ] = 2; /*C */
    ncf[ 19 * kd + 0 ] = 3; /*H */

    /*CH2CO */
    ncf[ 20 * kd + 1 ] = 2; /*C */
    ncf[ 20 * kd + 0 ] = 2; /*H */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*HCCO */
    ncf[ 21 * kd + 0 ] = 1; /*H */
    ncf[ 21 * kd + 1 ] = 2; /*C */
    ncf[ 21 * kd + 2 ] = 1; /*O */

    /*CH3CO */
    ncf[ 22 * kd + 1 ] = 2; /*C */
    ncf[ 22 * kd + 0 ] = 3; /*H */
    ncf[ 22 * kd + 2 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 23 * kd + 2 ] = 1; /*O */
    ncf[ 23 * kd + 0 ] = 3; /*H */
    ncf[ 23 * kd + 1 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 24 * kd + 1 ] = 2; /*C */
    ncf[ 24 * kd + 2 ] = 1; /*O */
    ncf[ 24 * kd + 0 ] = 4; /*H */

    /*C3H4-A */
    ncf[ 25 * kd + 0 ] = 4; /*H */
    ncf[ 25 * kd + 1 ] = 3; /*C */

    /*C3H6 */
    ncf[ 26 * kd + 1 ] = 3; /*C */
    ncf[ 26 * kd + 0 ] = 6; /*H */

    /*C4H6 */
    ncf[ 27 * kd + 1 ] = 4; /*C */
    ncf[ 27 * kd + 0 ] = 6; /*H */

    /*NC3H7 */
    ncf[ 28 * kd + 1 ] = 3; /*C */
    ncf[ 28 * kd + 0 ] = 7; /*H */

    /*C4H7 */
    ncf[ 29 * kd + 1 ] = 4; /*C */
    ncf[ 29 * kd + 0 ] = 7; /*H */

    /*C4H8-1 */
    ncf[ 30 * kd + 1 ] = 4; /*C */
    ncf[ 30 * kd + 0 ] = 8; /*H */

    /*PC4H9 */
    ncf[ 31 * kd + 1 ] = 4; /*C */
    ncf[ 31 * kd + 0 ] = 9; /*H */

    /*CH3COCH2 */
    ncf[ 32 * kd + 1 ] = 3; /*C */
    ncf[ 32 * kd + 0 ] = 5; /*H */
    ncf[ 32 * kd + 2 ] = 1; /*O */

    /*C2H5CHO */
    ncf[ 33 * kd + 1 ] = 3; /*C */
    ncf[ 33 * kd + 0 ] = 6; /*H */
    ncf[ 33 * kd + 2 ] = 1; /*O */

    /*C2H5CO */
    ncf[ 34 * kd + 1 ] = 3; /*C */
    ncf[ 34 * kd + 0 ] = 5; /*H */
    ncf[ 34 * kd + 2 ] = 1; /*O */

    /*C5H9 */
    ncf[ 35 * kd + 1 ] = 5; /*C */
    ncf[ 35 * kd + 0 ] = 9; /*H */

    /*C5H10-1 */
    ncf[ 36 * kd + 1 ] = 5; /*C */
    ncf[ 36 * kd + 0 ] = 10; /*H */

    /*C2H5O */
    ncf[ 37 * kd + 1 ] = 2; /*C */
    ncf[ 37 * kd + 0 ] = 5; /*H */
    ncf[ 37 * kd + 2 ] = 1; /*O */

    /*CH3O2 */
    ncf[ 38 * kd + 1 ] = 1; /*C */
    ncf[ 38 * kd + 0 ] = 3; /*H */
    ncf[ 38 * kd + 2 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 39 * kd + 1 ] = 1; /*C */
    ncf[ 39 * kd + 0 ] = 4; /*H */
    ncf[ 39 * kd + 2 ] = 2; /*O */

    /*C2H3CO */
    ncf[ 40 * kd + 1 ] = 3; /*C */
    ncf[ 40 * kd + 0 ] = 3; /*H */
    ncf[ 40 * kd + 2 ] = 1; /*O */

    /*C2H3CHO */
    ncf[ 41 * kd + 1 ] = 3; /*C */
    ncf[ 41 * kd + 0 ] = 4; /*H */
    ncf[ 41 * kd + 2 ] = 1; /*O */

    /*C3H5O */
    ncf[ 42 * kd + 1 ] = 3; /*C */
    ncf[ 42 * kd + 0 ] = 5; /*H */
    ncf[ 42 * kd + 2 ] = 1; /*O */

    /*C4H7O */
    ncf[ 43 * kd + 1 ] = 4; /*C */
    ncf[ 43 * kd + 0 ] = 7; /*H */
    ncf[ 43 * kd + 2 ] = 1; /*O */

    /*C4H8OOH1-3O2 */
    ncf[ 44 * kd + 1 ] = 4; /*C */
    ncf[ 44 * kd + 0 ] = 9; /*H */
    ncf[ 44 * kd + 2 ] = 4; /*O */

    /*C4H8OOH1-3 */
    ncf[ 45 * kd + 1 ] = 4; /*C */
    ncf[ 45 * kd + 0 ] = 9; /*H */
    ncf[ 45 * kd + 2 ] = 2; /*O */

    /*PC4H9O2 */
    ncf[ 46 * kd + 1 ] = 4; /*C */
    ncf[ 46 * kd + 0 ] = 9; /*H */
    ncf[ 46 * kd + 2 ] = 2; /*O */

    /*C3H5-A */
    ncf[ 47 * kd + 1 ] = 3; /*C */
    ncf[ 47 * kd + 0 ] = 5; /*H */

    /*C3H3 */
    ncf[ 48 * kd + 1 ] = 3; /*C */
    ncf[ 48 * kd + 0 ] = 3; /*H */

    /*C3H2 */
    ncf[ 49 * kd + 0 ] = 2; /*H */
    ncf[ 49 * kd + 1 ] = 3; /*C */

    /*CH2(S) */
    ncf[ 50 * kd + 1 ] = 1; /*C */
    ncf[ 50 * kd + 0 ] = 2; /*H */

    /*NC4KET13 */
    ncf[ 51 * kd + 1 ] = 4; /*C */
    ncf[ 51 * kd + 0 ] = 8; /*H */
    ncf[ 51 * kd + 2 ] = 3; /*O */

    /*NC3H7CHO */
    ncf[ 52 * kd + 1 ] = 4; /*C */
    ncf[ 52 * kd + 0 ] = 8; /*H */
    ncf[ 52 * kd + 2 ] = 1; /*O */

    /*NC3H7CO */
    ncf[ 53 * kd + 1 ] = 4; /*C */
    ncf[ 53 * kd + 0 ] = 7; /*H */
    ncf[ 53 * kd + 2 ] = 1; /*O */

    /*C2H5COCH2 */
    ncf[ 54 * kd + 1 ] = 4; /*C */
    ncf[ 54 * kd + 0 ] = 7; /*H */
    ncf[ 54 * kd + 2 ] = 1; /*O */

    /*NC3H7COCH2 */
    ncf[ 55 * kd + 1 ] = 5; /*C */
    ncf[ 55 * kd + 0 ] = 9; /*H */
    ncf[ 55 * kd + 2 ] = 1; /*O */

    /*NC4H9CHO */
    ncf[ 56 * kd + 1 ] = 5; /*C */
    ncf[ 56 * kd + 0 ] = 10; /*H */
    ncf[ 56 * kd + 2 ] = 1; /*O */

    /*NC4H9CO */
    ncf[ 57 * kd + 1 ] = 5; /*C */
    ncf[ 57 * kd + 0 ] = 9; /*H */
    ncf[ 57 * kd + 2 ] = 1; /*O */

    /*NC7H16 */
    ncf[ 58 * kd + 1 ] = 7; /*C */
    ncf[ 58 * kd + 0 ] = 16; /*H */

    /*C7H15-1 */
    ncf[ 59 * kd + 1 ] = 7; /*C */
    ncf[ 59 * kd + 0 ] = 15; /*H */

    /*C7H15-2 */
    ncf[ 60 * kd + 1 ] = 7; /*C */
    ncf[ 60 * kd + 0 ] = 15; /*H */

    /*C7H15-3 */
    ncf[ 61 * kd + 1 ] = 7; /*C */
    ncf[ 61 * kd + 0 ] = 15; /*H */

    /*C7H15-4 */
    ncf[ 62 * kd + 1 ] = 7; /*C */
    ncf[ 62 * kd + 0 ] = 15; /*H */

    /*C7H14-2 */
    ncf[ 63 * kd + 1 ] = 7; /*C */
    ncf[ 63 * kd + 0 ] = 14; /*H */

    /*C7H14-3 */
    ncf[ 64 * kd + 1 ] = 7; /*C */
    ncf[ 64 * kd + 0 ] = 14; /*H */

    /*C7H15O2-1 */
    ncf[ 65 * kd + 1 ] = 7; /*C */
    ncf[ 65 * kd + 0 ] = 15; /*H */
    ncf[ 65 * kd + 2 ] = 2; /*O */

    /*C7H15O2-2 */
    ncf[ 66 * kd + 1 ] = 7; /*C */
    ncf[ 66 * kd + 0 ] = 15; /*H */
    ncf[ 66 * kd + 2 ] = 2; /*O */

    /*C7H15O2-3 */
    ncf[ 67 * kd + 1 ] = 7; /*C */
    ncf[ 67 * kd + 0 ] = 15; /*H */
    ncf[ 67 * kd + 2 ] = 2; /*O */

    /*C7H15O2-4 */
    ncf[ 68 * kd + 1 ] = 7; /*C */
    ncf[ 68 * kd + 0 ] = 15; /*H */
    ncf[ 68 * kd + 2 ] = 2; /*O */

    /*C7H14OOH1-3 */
    ncf[ 69 * kd + 1 ] = 7; /*C */
    ncf[ 69 * kd + 0 ] = 15; /*H */
    ncf[ 69 * kd + 2 ] = 2; /*O */

    /*C7H14OOH2-3 */
    ncf[ 70 * kd + 1 ] = 7; /*C */
    ncf[ 70 * kd + 0 ] = 15; /*H */
    ncf[ 70 * kd + 2 ] = 2; /*O */

    /*C7H14OOH2-4 */
    ncf[ 71 * kd + 1 ] = 7; /*C */
    ncf[ 71 * kd + 0 ] = 15; /*H */
    ncf[ 71 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-2 */
    ncf[ 72 * kd + 1 ] = 7; /*C */
    ncf[ 72 * kd + 0 ] = 15; /*H */
    ncf[ 72 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-4 */
    ncf[ 73 * kd + 1 ] = 7; /*C */
    ncf[ 73 * kd + 0 ] = 15; /*H */
    ncf[ 73 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-5 */
    ncf[ 74 * kd + 1 ] = 7; /*C */
    ncf[ 74 * kd + 0 ] = 15; /*H */
    ncf[ 74 * kd + 2 ] = 2; /*O */

    /*C7H14OOH4-2 */
    ncf[ 75 * kd + 1 ] = 7; /*C */
    ncf[ 75 * kd + 0 ] = 15; /*H */
    ncf[ 75 * kd + 2 ] = 2; /*O */

    /*C7H14OOH4-3 */
    ncf[ 76 * kd + 1 ] = 7; /*C */
    ncf[ 76 * kd + 0 ] = 15; /*H */
    ncf[ 76 * kd + 2 ] = 2; /*O */

    /*C7H14OOH1-3O2 */
    ncf[ 77 * kd + 1 ] = 7; /*C */
    ncf[ 77 * kd + 0 ] = 15; /*H */
    ncf[ 77 * kd + 2 ] = 4; /*O */

    /*C7H14OOH2-4O2 */
    ncf[ 78 * kd + 1 ] = 7; /*C */
    ncf[ 78 * kd + 0 ] = 15; /*H */
    ncf[ 78 * kd + 2 ] = 4; /*O */

    /*C7H14OOH3-5O2 */
    ncf[ 79 * kd + 1 ] = 7; /*C */
    ncf[ 79 * kd + 0 ] = 15; /*H */
    ncf[ 79 * kd + 2 ] = 4; /*O */

    /*C7H14OOH4-2O2 */
    ncf[ 80 * kd + 1 ] = 7; /*C */
    ncf[ 80 * kd + 0 ] = 15; /*H */
    ncf[ 80 * kd + 2 ] = 4; /*O */

    /*C7H14O1-3 */
    ncf[ 81 * kd + 1 ] = 7; /*C */
    ncf[ 81 * kd + 0 ] = 14; /*H */
    ncf[ 81 * kd + 2 ] = 1; /*O */

    /*C7H14O2-4 */
    ncf[ 82 * kd + 1 ] = 7; /*C */
    ncf[ 82 * kd + 0 ] = 14; /*H */
    ncf[ 82 * kd + 2 ] = 1; /*O */

    /*NC7KET13 */
    ncf[ 83 * kd + 1 ] = 7; /*C */
    ncf[ 83 * kd + 0 ] = 14; /*H */
    ncf[ 83 * kd + 2 ] = 3; /*O */

    /*NC7KET24 */
    ncf[ 84 * kd + 1 ] = 7; /*C */
    ncf[ 84 * kd + 0 ] = 14; /*H */
    ncf[ 84 * kd + 2 ] = 3; /*O */

    /*NC7KET35 */
    ncf[ 85 * kd + 1 ] = 7; /*C */
    ncf[ 85 * kd + 0 ] = 14; /*H */
    ncf[ 85 * kd + 2 ] = 3; /*O */

    /*NC7KET42 */
    ncf[ 86 * kd + 1 ] = 7; /*C */
    ncf[ 86 * kd + 0 ] = 14; /*H */
    ncf[ 86 * kd + 2 ] = 3; /*O */

    /*N2 */
    ncf[ 87 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "H";
    ename[1] = "C";
    ename[2] = "O";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(88);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "CO";
    kname[7] = "HCO";
    kname[8] = "CO2";
    kname[9] = "CH3";
    kname[10] = "CH4";
    kname[11] = "HO2";
    kname[12] = "H2O2";
    kname[13] = "CH2O";
    kname[14] = "CH3O";
    kname[15] = "C2H6";
    kname[16] = "C2H4";
    kname[17] = "C2H5";
    kname[18] = "C2H2";
    kname[19] = "C2H3";
    kname[20] = "CH2CO";
    kname[21] = "HCCO";
    kname[22] = "CH3CO";
    kname[23] = "CH2CHO";
    kname[24] = "CH3CHO";
    kname[25] = "C3H4-A";
    kname[26] = "C3H6";
    kname[27] = "C4H6";
    kname[28] = "NC3H7";
    kname[29] = "C4H7";
    kname[30] = "C4H8-1";
    kname[31] = "PC4H9";
    kname[32] = "CH3COCH2";
    kname[33] = "C2H5CHO";
    kname[34] = "C2H5CO";
    kname[35] = "C5H9";
    kname[36] = "C5H10-1";
    kname[37] = "C2H5O";
    kname[38] = "CH3O2";
    kname[39] = "CH3O2H";
    kname[40] = "C2H3CO";
    kname[41] = "C2H3CHO";
    kname[42] = "C3H5O";
    kname[43] = "C4H7O";
    kname[44] = "C4H8OOH1-3O2";
    kname[45] = "C4H8OOH1-3";
    kname[46] = "PC4H9O2";
    kname[47] = "C3H5-A";
    kname[48] = "C3H3";
    kname[49] = "C3H2";
    kname[50] = "CH2(S)";
    kname[51] = "NC4KET13";
    kname[52] = "NC3H7CHO";
    kname[53] = "NC3H7CO";
    kname[54] = "C2H5COCH2";
    kname[55] = "NC3H7COCH2";
    kname[56] = "NC4H9CHO";
    kname[57] = "NC4H9CO";
    kname[58] = "NC7H16";
    kname[59] = "C7H15-1";
    kname[60] = "C7H15-2";
    kname[61] = "C7H15-3";
    kname[62] = "C7H15-4";
    kname[63] = "C7H14-2";
    kname[64] = "C7H14-3";
    kname[65] = "C7H15O2-1";
    kname[66] = "C7H15O2-2";
    kname[67] = "C7H15O2-3";
    kname[68] = "C7H15O2-4";
    kname[69] = "C7H14OOH1-3";
    kname[70] = "C7H14OOH2-3";
    kname[71] = "C7H14OOH2-4";
    kname[72] = "C7H14OOH3-2";
    kname[73] = "C7H14OOH3-4";
    kname[74] = "C7H14OOH3-5";
    kname[75] = "C7H14OOH4-2";
    kname[76] = "C7H14OOH4-3";
    kname[77] = "C7H14OOH1-3O2";
    kname[78] = "C7H14OOH2-4O2";
    kname[79] = "C7H14OOH3-5O2";
    kname[80] = "C7H14OOH4-2O2";
    kname[81] = "C7H14O1-3";
    kname[82] = "C7H14O2-4";
    kname[83] = "NC7KET13";
    kname[84] = "NC7KET24";
    kname[85] = "NC7KET35";
    kname[86] = "NC7KET42";
    kname[87] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(Jac[ 89 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 89 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 89 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 89;
        int offset_col = nc * 89;
        for (int k=0; k<89; k++) {
            for (int l=0; l<89; l++) {
                if(Jac[89*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if(Jac[89*k + l] != 0.0) {
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
            int offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if(Jac[89*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[89*k + l] != 0.0) {
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
            int offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[89*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 89*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[89*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 89*k + l;
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
    amrex::GpuArray<amrex::Real,7921> Jac = {0.0};
    amrex::GpuArray<amrex::Real,88> conc = {0.0};
    for (int n=0; n<88; n++) {
        conc[n] = 1.0/ 88.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<89; l++) {
            for (int k=0; k<89; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[89*k + l] != 0.0) {
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
        for (int l=0; l<89; l++) {
            for (int k=0; k<89; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[89*k + l] != 0.0) {
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
