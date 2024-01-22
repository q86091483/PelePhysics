#include "mechanism.H"
const int rmap[210] = {
  11,  20,  26,  78,  113, 122, 129, 138, 142, 148, 153, 165, 6,   7,   8,
  9,   10,  91,  92,  102, 147, 189, 0,   1,   2,   3,   4,   5,   12,  13,
  14,  15,  16,  17,  18,  19,  21,  22,  23,  24,  25,  27,  28,  29,  30,
  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,
  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
  76,  77,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  93,
  94,  95,  96,  97,  98,  99,  100, 101, 103, 104, 105, 106, 107, 108, 109,
  110, 111, 112, 114, 115, 116, 117, 118, 119, 120, 121, 123, 124, 125, 126,
  127, 128, 130, 131, 132, 133, 134, 135, 136, 137, 139, 140, 141, 143, 144,
  145, 146, 149, 150, 151, 152, 154, 155, 156, 157, 158, 159, 160, 161, 162,
  163, 164, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
  179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 190, 191, 192, 193, 194,
  195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 210; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[210] = {
    4, 4, 4, 4, 3, 3, 2, 3, 2, 3, 3, 3, 4, 3, 4, 4, 4, 4, 3, 3, 2, 4, 4, 4,
    4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3,
    4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 3, 4,
    4, 4, 5, 3, 3, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 2, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 2, 4, 4, 4, 4, 4, 5, 5, 4, 4, 4};
  const int kiv[1050] = {
    2,  1,  3,  4,  0, 0,  3,  2,  4,  0, 0,  3,  2,  4,  0, 0,  4,  2,  6,  0,
    4,  6,  3,  0,  0, 4,  6,  3,  0,  0, 0,  2,  0,  0,  0, 2,  3,  4,  0,  0,
    3,  1,  0,  0,  0, 6,  2,  4,  0,  0, 6,  2,  4,  0,  0, 2,  1,  5,  0,  0,
    2,  5,  0,  1,  0, 2,  5,  4,  0,  0, 2,  5,  6,  3,  0, 5,  3,  1,  4,  0,
    5,  4,  6,  1,  0, 5,  4,  6,  1,  0, 5,  7,  1,  0,  0, 5,  7,  1,  0,  0,
    7,  4,  0,  0,  0, 2,  7,  6,  4,  0, 2,  7,  0,  5,  0, 7,  3,  5,  4,  0,
    7,  4,  6,  5,  0, 7,  4,  6,  5,  0, 2,  10, 9,  0,  0, 2,  9,  0,  10, 0,
    9,  3,  10, 4,  0, 9,  4,  6,  10, 0, 5,  9,  7,  10, 0, 0,  11, 2,  10, 0,
    10, 3,  2,  21, 0, 10, 3,  11, 4,  0, 10, 3,  11, 4,  0, 10, 4,  6,  11, 0,
    5,  10, 9,  1,  0, 5,  10, 19, 4,  0, 5,  10, 6,  21, 0, 5,  10, 6,  21, 0,
    5,  10, 6,  22, 0, 10, 1,  19, 3,  0, 10, 1,  21, 4,  0, 10, 11, 9,  0,  0,
    11, 10, 12, 9,  0, 12, 10, 2,  29, 0, 21, 10, 9,  8,  0, 10, 8,  6,  29, 0,
    10, 8,  13, 4,  0, 24, 10, 9,  23, 0, 10, 23, 19, 8,  0, 10, 23, 6,  28, 0,
    2,  11, 0,  12, 0, 11, 3,  2,  8,  0, 11, 4,  2,  21, 0, 11, 4,  6,  12, 0,
    11, 1,  21, 3,  0, 11, 1,  8,  4,  0, 11, 12, 10, 0,  0, 12, 11, 2,  29, 0,
    11, 8,  2,  28, 0, 11, 8,  29, 4,  0, 24, 11, 10, 23, 0, 11, 23, 28, 4,  0,
    11, 23, 21, 8,  0, 12, 4,  2,  8,  0, 12, 1,  8,  3,  0, 12, 8,  29, 3,  0,
    13, 2,  29, 0,  0, 2,  13, 0,  29, 0, 13, 3,  2,  28, 0, 13, 3,  29, 4,  0,
    13, 3,  11, 8,  0, 13, 4,  6,  29, 0, 13, 1,  5,  29, 0, 11, 13, 29, 10, 0,
    10, 13, 29, 9,  0, 13, 8,  21, 29, 0, 18, 10, 4,  0,  0, 2,  18, 0,  20, 0,
    2,  18, 0,  19, 0, 18, 3,  20, 4,  0, 18, 3,  19, 4,  0, 18, 4,  6,  20, 0,
    18, 4,  19, 6,  0, 10, 18, 20, 9,  0, 10, 18, 19, 9,  0, 11, 18, 20, 10, 0,
    11, 18, 19, 10, 0, 5,  18, 7,  20, 0, 5,  18, 19, 7,  0, 19, 2,  21, 0,  0,
    19, 20, 0,  0,  0, 2,  19, 0,  21, 0, 2,  19, 10, 4,  0, 19, 3,  21, 4,  0,
    19, 4,  6,  21, 0, 19, 5,  7,  21, 0, 19, 1,  21, 5,  0, 19, 10, 21, 9,  0,
    19, 8,  21, 0,  0, 19, 23, 21, 24, 0, 20, 2,  21, 0,  0, 2,  20, 10, 4,  0,
    2,  20, 0,  21, 0, 20, 3,  21, 4,  0, 20, 3,  21, 4,  0, 20, 4,  6,  21, 0,
    20, 5,  7,  21, 0, 20, 5,  18, 1,  0, 20, 1,  21, 5,  0, 20, 10, 21, 9,  0,
    20, 23, 21, 24, 0, 2,  8,  21, 0,  0, 2,  21, 0,  8,  0, 21, 3,  8,  4,  0,
    21, 4,  6,  8,  0, 21, 5,  25, 4,  0, 21, 1,  5,  8,  0, 21, 6,  28, 0,  0,
    21, 23, 24, 8,  0, 5,  8,  23, 4,  0, 8,  3,  23, 0,  0, 2,  23, 8,  4,  0,
    23, 3,  8,  1,  0, 5,  23, 24, 1,  0, 5,  23, 25, 1,  0, 23, 8,  1,  0,  0,
    23, 8,  26, 0,  0, 8,  4,  24, 0,  0, 0,  23, 2,  24, 0, 0,  23, 2,  24, 0,
    2,  24, 21, 4,  0, 2,  24, 6,  8,  0, 24, 3,  23, 4,  0, 24, 4,  6,  23, 0,
    24, 23, 27, 8,  0, 24, 6,  8,  23, 0, 25, 24, 0,  0,  0, 0,  23, 2,  25, 0,
    25, 3,  23, 4,  0, 25, 4,  6,  23, 0, 23, 3,  26, 0,  0, 2,  26, 23, 4,  0,
    26, 3,  23, 1,  0, 26, 4,  5,  23, 0, 5,  26, 23, 1,  4, 26, 8,  1,  0,  0,
    23, 4,  27, 0,  0, 2,  27, 0,  26, 0, 2,  27, 6,  23, 0, 2,  27, 24, 4,  0,
    27, 4,  6,  26, 0, 28, 29, 3,  0,  0, 2,  28, 29, 4,  0, 28, 3,  8,  0,  0,
    28, 3,  29, 1,  0, 28, 4,  5,  29, 0, 28, 4,  21, 8,  0, 28, 8,  29, 23, 0,
    10, 0,  16, 0,  0, 10, 0,  17, 0,  0, 11, 10, 2,  16, 0, 20, 10, 15, 4,  0,
    20, 10, 17, 6,  0, 10, 14, 0,  0,  0, 14, 0,  17, 0,  0, 2,  14, 0,  15, 0,
    2,  14, 10, 9,  0, 14, 3,  15, 4,  0, 14, 3,  6,  16, 0, 14, 4,  6,  15, 0,
    14, 10, 15, 9,  0, 14, 8,  21, 15, 0, 14, 23, 24, 15, 0, 14, 23, 25, 15, 0,
    2,  16, 15, 0,  0, 2,  15, 0,  16, 0, 15, 3,  16, 4,  0, 15, 3,  21, 10, 0,
    15, 3,  2,  10, 8, 15, 4,  6,  16, 0, 15, 4,  17, 6,  0, 15, 4,  21, 9,  0,
    5,  15, 7,  16, 0, 5,  15, 14, 1,  0, 15, 10, 16, 9,  0, 15, 10, 17, 9,  0,
    15, 11, 16, 10, 0, 16, 2,  13, 0,  0, 2,  16, 0,  13, 0, 16, 3,  13, 4,  0,
    16, 3,  10, 8,  0, 16, 4,  6,  13, 0, 16, 10, 9,  13, 0, 16, 11, 10, 13, 0,
    16, 8,  28, 10, 0, 17, 2,  13, 0,  0, 17, 0,  29, 0,  0, 17, 16, 0,  0,  0,
    2,  17, 0,  13, 0, 2,  17, 2,  16, 0, 17, 3,  13, 4,  0, 17, 3,  10, 8,  0,
    17, 4,  6,  13, 0, 17, 4,  2,  10, 8, 17, 5,  10, 8,  4, 17, 5,  7,  13, 0,
    17, 1,  10, 23, 0, 17, 10, 9,  13, 0};
  const int nuv[1050] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, 2,  0, 0, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 2,  1, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 1, 0, -1, 1,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 210) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
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
  amrex::Real c[30]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 30; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 210; ++id) {
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
  amrex::Real g_RT[30];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 1.008000;  // H
  awt[1] = 15.999000; // O
  awt[2] = 14.007000; // N
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
  int kd = 3;
  // Zero ncf
  for (int id = 0; id < kd * 30; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 0] = 2; // H

  // O2
  ncf[1 * kd + 1] = 2; // O

  // H
  ncf[2 * kd + 0] = 1; // H

  // O
  ncf[3 * kd + 1] = 1; // O

  // OH
  ncf[4 * kd + 0] = 1; // H
  ncf[4 * kd + 1] = 1; // O

  // HO2
  ncf[5 * kd + 0] = 1; // H
  ncf[5 * kd + 1] = 2; // O

  // H2O
  ncf[6 * kd + 0] = 2; // H
  ncf[6 * kd + 1] = 1; // O

  // H2O2
  ncf[7 * kd + 0] = 2; // H
  ncf[7 * kd + 1] = 2; // O

  // NO
  ncf[8 * kd + 2] = 1; // N
  ncf[8 * kd + 1] = 1; // O

  // NH3
  ncf[9 * kd + 0] = 3; // H
  ncf[9 * kd + 2] = 1; // N

  // NH2
  ncf[10 * kd + 0] = 2; // H
  ncf[10 * kd + 2] = 1; // N

  // NH
  ncf[11 * kd + 0] = 1; // H
  ncf[11 * kd + 2] = 1; // N

  // N
  ncf[12 * kd + 2] = 1; // N

  // NNH
  ncf[13 * kd + 0] = 1; // H
  ncf[13 * kd + 2] = 2; // N

  // N2H4
  ncf[14 * kd + 0] = 4; // H
  ncf[14 * kd + 2] = 2; // N

  // N2H3
  ncf[15 * kd + 0] = 3; // H
  ncf[15 * kd + 2] = 2; // N

  // N2H2
  ncf[16 * kd + 0] = 2; // H
  ncf[16 * kd + 2] = 2; // N

  // H2NN
  ncf[17 * kd + 0] = 2; // H
  ncf[17 * kd + 2] = 2; // N

  // NH2OH
  ncf[18 * kd + 0] = 3; // H
  ncf[18 * kd + 2] = 1; // N
  ncf[18 * kd + 1] = 1; // O

  // H2NO
  ncf[19 * kd + 0] = 2; // H
  ncf[19 * kd + 2] = 1; // N
  ncf[19 * kd + 1] = 1; // O

  // HNOH
  ncf[20 * kd + 0] = 2; // H
  ncf[20 * kd + 2] = 1; // N
  ncf[20 * kd + 1] = 1; // O

  // HNO
  ncf[21 * kd + 0] = 1; // H
  ncf[21 * kd + 2] = 1; // N
  ncf[21 * kd + 1] = 1; // O

  // HON
  ncf[22 * kd + 0] = 1; // H
  ncf[22 * kd + 2] = 1; // N
  ncf[22 * kd + 1] = 1; // O

  // NO2
  ncf[23 * kd + 2] = 1; // N
  ncf[23 * kd + 1] = 2; // O

  // HONO
  ncf[24 * kd + 0] = 1; // H
  ncf[24 * kd + 2] = 1; // N
  ncf[24 * kd + 1] = 2; // O

  // HNO2
  ncf[25 * kd + 0] = 1; // H
  ncf[25 * kd + 2] = 1; // N
  ncf[25 * kd + 1] = 2; // O

  // NO3
  ncf[26 * kd + 2] = 1; // N
  ncf[26 * kd + 1] = 3; // O

  // HONO2
  ncf[27 * kd + 0] = 1; // H
  ncf[27 * kd + 2] = 1; // N
  ncf[27 * kd + 1] = 3; // O

  // N2O
  ncf[28 * kd + 2] = 2; // N
  ncf[28 * kd + 1] = 1; // O

  // N2
  ncf[29 * kd + 2] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(3);
  ename[0] = "H";
  ename[1] = "O";
  ename[2] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(30);
  kname[0] = "H2";
  kname[1] = "O2";
  kname[2] = "H";
  kname[3] = "O";
  kname[4] = "OH";
  kname[5] = "HO2";
  kname[6] = "H2O";
  kname[7] = "H2O2";
  kname[8] = "NO";
  kname[9] = "NH3";
  kname[10] = "NH2";
  kname[11] = "NH";
  kname[12] = "N";
  kname[13] = "NNH";
  kname[14] = "N2H4";
  kname[15] = "N2H3";
  kname[16] = "N2H2";
  kname[17] = "H2NN";
  kname[18] = "NH2OH";
  kname[19] = "H2NO";
  kname[20] = "HNOH";
  kname[21] = "HNO";
  kname[22] = "HON";
  kname[23] = "NO2";
  kname[24] = "HONO";
  kname[25] = "HNO2";
  kname[26] = "NO3";
  kname[27] = "HONO2";
  kname[28] = "N2O";
  kname[29] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 31;
    int offset_col = nc * 31;
    for (int k = 0; k < 31; k++) {
      for (int l = 0; l < 31; l++) {
        if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (Jac[31 * k + l] != 0.0) {
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
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[31 * k + l] != 0.0) {
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
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[31 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 31 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 31 * k + l;
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
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 31; l++) {
      for (int k = 0; k < 31; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[31 * k + l] != 0.0) {
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
    for (int l = 0; l < 31; l++) {
      for (int k = 0; k < 31; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[31 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
