//
// Created by logan on 8/10/20.
//

#ifndef DARKSUN_DC_DECSOL_HPP
#define DARKSUN_DC_DECSOL_HPP

#include "stiff/common.hpp"
#include "stiff/decsol.hpp"

namespace stiff {

static int c__1 = 1;

int decomr(int *n, double *fjac, int *ldjac, const double *fmas,
           const int *ldmas, const int *mlmas, const int *mumas, const int *m1,
           const int *m2, int *nm1, const double *fac1, double *e1, int *lde1,
           int *ip1, int *ier, const int *ijob, bool *calhes, int *iphes,
           StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset, i_1,
      i_2, i_3, i_4, i_5, i_6;

  /* Local variables */
  int i_, j, k, j1, ib, mm, jm1;
  double sum;

  /* Parameter adjustments */
  --iphes;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip1;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L14;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *n;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e1[i_ + j * e1_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    e1[j + j * e1_dim1] += *fac1;
  }
  dec(n, lde1, &e1[e1_offset], &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = *nm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e1[i_ + j * e1_dim1] = -fjac[i_ + jm1 * fjac_dim1];
    }
    e1[j + j * e1_dim1] += *fac1;
  }
L45:
  mm = *m1 / *m2;
  i_1 = *m2;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *nm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      sum = 0.;
      i_3 = mm - 1;
      for (k = 0; k <= i_3; ++k) {
        sum = (sum + fjac[i_ + (j + k * *m2) * fjac_dim1]) / *fac1;
      }
      e1[i_ + j * e1_dim1] -= sum;
    }
  }
  dec(nm1, lde1, &e1[e1_offset], &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e1[i_ + linal_1.mle + j * e1_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    e1[linal_1.mdiag + j * e1_dim1] += *fac1;
  }
  decb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e1[i_ + linal_1.mle + j * e1_dim1] = -fjac[i_ + jm1 * fjac_dim1];
    }
    e1[linal_1.mdiag + j * e1_dim1] += *fac1;
  }
L46:
  mm = *m1 / *m2;
  i_1 = *m2;
  for (j = 1; j <= i_1; ++j) {
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      sum = 0.;
      i_3 = mm - 1;
      for (k = 0; k <= i_3; ++k) {
        sum = (sum + fjac[i_ + (j + k * *m2) * fjac_dim1]) / *fac1;
      }
      e1[i_ + linal_1.mle + j * e1_dim1] -= sum;
    }
  }
  decb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *n;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e1[i_ + j * e1_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    /* Computing std::max */
    i_2 = 1, i_3 = j - *mumas;
    /* Computing std::min */
    i_5 = *n, i_6 = j + *mlmas;
    i_4 = std::min(i_5, i_6);
    for (i_ = std::max(i_2, i_3); i_ <= i_4; ++i_) {
      e1[i_ + j * e1_dim1] +=
          *fac1 * fmas[i_ - j + linal_1.mbdiag + j * fmas_dim1];
    }
  }
  dec(n, lde1, &e1[e1_offset], &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_4 = *nm1;
    for (i_ = 1; i_ <= i_4; ++i_) {
      e1[i_ + j * e1_dim1] = -fjac[i_ + jm1 * fjac_dim1];
    }
    /* Computing std::max */
    i_4 = 1, i_2 = j - *mumas;
    /* Computing std::min */
    i_5 = *nm1, i_6 = j + *mlmas;
    i_3 = std::min(i_5, i_6);
    for (i_ = std::max(i_4, i_2); i_ <= i_3; ++i_) {
      e1[i_ + j * e1_dim1] +=
          *fac1 * fmas[i_ - j + linal_1.mbdiag + j * fmas_dim1];
    }
  }
  goto L45;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_3 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_3; ++i_) {
      e1[i_ + linal_1.mle + j * e1_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    i_3 = linal_1.mbb;
    for (i_ = 1; i_ <= i_3; ++i_) {
      ib = i_ + linal_1.mdiff;
      e1[ib + j * e1_dim1] += *fac1 * fmas[i_ + j * fmas_dim1];
    }
  }
  decb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L14:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_3 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_3; ++i_) {
      e1[i_ + linal_1.mle + j * e1_dim1] = -fjac[i_ + jm1 * fjac_dim1];
    }
    i_3 = linal_1.mbb;
    for (i_ = 1; i_ <= i_3; ++i_) {
      ib = i_ + linal_1.mdiff;
      e1[ib + j * e1_dim1] += *fac1 * fmas[i_ + j * fmas_dim1];
    }
  }
  goto L46;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_3 = *n;
    for (i_ = 1; i_ <= i_3; ++i_) {
      e1[i_ + j * e1_dim1] =
          fmas[i_ + j * fmas_dim1] * *fac1 - fjac[i_ + j * fjac_dim1];
    }
  }
  dec(n, lde1, &e1[e1_offset], &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_3 = *nm1;
    for (i_ = 1; i_ <= i_3; ++i_) {
      e1[i_ + j * e1_dim1] =
          fmas[i_ + j * fmas_dim1] * *fac1 - fjac[i_ + jm1 * fjac_dim1];
    }
  }
  goto L45;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  if (*calhes) {
    elmhes(ldjac, n, &c__1, n, &fjac[fjac_offset], &iphes[1]);
  }
  *calhes = false;
  i_1 = *n - 1;
  for (j = 1; j <= i_1; ++j) {
    j1 = j + 1;
    e1[j1 + j * e1_dim1] = -fjac[j1 + j * fjac_dim1];
  }
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_3 = j;
    for (i_ = 1; i_ <= i_3; ++i_) {
      e1[i_ + j * e1_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    e1[j + j * e1_dim1] += *fac1;
  }
  dech(n, lde1, &e1[e1_offset], &c__1, &ip1[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
}

int decomc(int *n, const double *fjac, const int *ldjac, const double *fmas,
           const int *ldmas, const int *mlmas, const int *mumas, const int *m1,
           const int *m2, int *nm1, const double *alphn, const double *betan,
           double *e2r, double *e2i, int *lde1, int *ip2, int *ier,
           const int *ijob, StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, e2r_offset,
      e2i_dim1, e2i_offset, i_1, i_2, i_3, i_4, i_5, i_6;
  double d_1, d_2;

  /* Local variables */
  int i_, j, k, j1;
  double bb;
  int ib, mm, jm1;
  double bet, alp;
  double ffma, abno;
  int imle;
  double sumi, sumr, sums;

  /* Parameter adjustments */
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip2;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e2i_dim1 = *lde1;
  e2i_offset = 1 + e2i_dim1;
  e2i -= e2i_offset;
  e2r_dim1 = *lde1;
  e2r_offset = 1 + e2r_dim1;
  e2r -= e2r_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L14;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *n;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + j * e2r_dim1] = -fjac[i_ + j * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = 0.;
    }
    e2r[j + j * e2r_dim1] += *alphn;
    e2i[j + j * e2i_dim1] = *betan;
  }
  decc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = *nm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + j * e2r_dim1] = -fjac[i_ + jm1 * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = 0.;
    }
    e2r[j + j * e2r_dim1] += *alphn;
    e2i[j + j * e2i_dim1] = *betan;
  }
L45:
  mm = *m1 / *m2;
  /* Computing 2nd power */
  d_1 = *alphn;
  /* Computing 2nd power */
  d_2 = *betan;
  abno = d_1 * d_1 + d_2 * d_2;
  alp = *alphn / abno;
  bet = *betan / abno;
  i_1 = *m2;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *nm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      sumr = 0.;
      sumi = 0.;
      i_3 = mm - 1;
      for (k = 0; k <= i_3; ++k) {
        sums = sumr + fjac[i_ + (j + k * *m2) * fjac_dim1];
        sumr = sums * alp + sumi * bet;
        sumi = sumi * alp - sums * bet;
      }
      e2r[i_ + j * e2r_dim1] -= sumr;
      e2i[i_ + j * e2i_dim1] -= sumi;
    }
  }
  decc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      imle = i_ + linal_1.mle;
      e2r[imle + j * e2r_dim1] = -fjac[i_ + j * fjac_dim1];
      e2i[imle + j * e2i_dim1] = 0.;
    }
    e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
    e2i[linal_1.mdiag + j * e2i_dim1] = *betan;
  }
  decbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + linal_1.mle + j * e2r_dim1] = -fjac[i_ + jm1 * fjac_dim1];
      e2i[i_ + linal_1.mle + j * e2i_dim1] = 0.;
    }
    e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
    e2i[linal_1.mdiag + j * e2i_dim1] += *betan;
  }
L46:
  mm = *m1 / *m2;
  /* Computing 2nd power */
  d_1 = *alphn;
  /* Computing 2nd power */
  d_2 = *betan;
  abno = d_1 * d_1 + d_2 * d_2;
  alp = *alphn / abno;
  bet = *betan / abno;
  i_1 = *m2;
  for (j = 1; j <= i_1; ++j) {
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      sumr = 0.;
      sumi = 0.;
      i_3 = mm - 1;
      for (k = 0; k <= i_3; ++k) {
        sums = sumr + fjac[i_ + (j + k * *m2) * fjac_dim1];
        sumr = sums * alp + sumi * bet;
        sumi = sumi * alp - sums * bet;
      }
      imle = i_ + linal_1.mle;
      e2r[imle + j * e2r_dim1] -= sumr;
      e2i[imle + j * e2i_dim1] -= sumi;
    }
  }
  decbc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle,
        &linal_1.mue, &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *n;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + j * e2r_dim1] = -fjac[i_ + j * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = 0.;
    }
  }
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    /* Computing std::max */
    i_2 = 1, i_3 = j - *mumas;
    /* Computing std::min */
    i_5 = *n, i_6 = j + *mlmas;
    i_4 = std::min(i_5, i_6);
    for (i_ = std::max(i_2, i_3); i_ <= i_4; ++i_) {
      bb = fmas[i_ - j + linal_1.mbdiag + j * fmas_dim1];
      e2r[i_ + j * e2r_dim1] += *alphn * bb;
      e2i[i_ + j * e2i_dim1] = *betan * bb;
    }
  }
  decc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_4 = *nm1;
    for (i_ = 1; i_ <= i_4; ++i_) {
      e2r[i_ + j * e2r_dim1] = -fjac[i_ + jm1 * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = 0.;
    }
    /* Computing std::max */
    i_4 = 1, i_2 = j - *mumas;
    /* Computing std::min */
    i_5 = *nm1, i_6 = j + *mlmas;
    i_3 = std::min(i_5, i_6);
    for (i_ = std::max(i_4, i_2); i_ <= i_3; ++i_) {
      ffma = fmas[i_ - j + linal_1.mbdiag + j * fmas_dim1];
      e2r[i_ + j * e2r_dim1] += *alphn * ffma;
      e2i[i_ + j * e2i_dim1] += *betan * ffma;
    }
  }
  goto L45;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_3 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_3; ++i_) {
      imle = i_ + linal_1.mle;
      e2r[imle + j * e2r_dim1] = -fjac[i_ + j * fjac_dim1];
      e2i[imle + j * e2i_dim1] = 0.;
    }
    /* Computing std::max */
    i_3 = 1, i_4 = *mumas + 2 - j;
    /* Computing std::min */
    i_5 = linal_1.mbb, i_6 = *mumas + 1 - j + *n;
    i_2 = std::min(i_5, i_6);
    for (i_ = std::max(i_3, i_4); i_ <= i_2; ++i_) {
      ib = i_ + linal_1.mdiff;
      bb = fmas[i_ + j * fmas_dim1];
      e2r[ib + j * e2r_dim1] += *alphn * bb;
      e2i[ib + j * e2i_dim1] = *betan * bb;
    }
  }
  decbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L14:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = linal_1.mbjac;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + linal_1.mle + j * e2r_dim1] = -fjac[i_ + jm1 * fjac_dim1];
      e2i[i_ + linal_1.mle + j * e2i_dim1] = 0.;
    }
    i_2 = linal_1.mbb;
    for (i_ = 1; i_ <= i_2; ++i_) {
      ib = i_ + linal_1.mdiff;
      ffma = fmas[i_ + j * fmas_dim1];
      e2r[ib + j * e2r_dim1] += *alphn * ffma;
      e2i[ib + j * e2i_dim1] += *betan * ffma;
    }
  }
  goto L46;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = *n;
    for (i_ = 1; i_ <= i_2; ++i_) {
      bb = fmas[i_ + j * fmas_dim1];
      e2r[i_ + j * e2r_dim1] = bb * *alphn - fjac[i_ + j * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = bb * *betan;
    }
  }
  decc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i_1 = *nm1;
  for (j = 1; j <= i_1; ++j) {
    jm1 = j + *m1;
    i_2 = *nm1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2r[i_ + j * e2r_dim1] =
          *alphn * fmas[i_ + j * fmas_dim1] - fjac[i_ + jm1 * fjac_dim1];
      e2i[i_ + j * e2i_dim1] = *betan * fmas[i_ + j * fmas_dim1];
    }
  }
  goto L45;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i_1 = *n - 1;
  for (j = 1; j <= i_1; ++j) {
    j1 = j + 1;
    e2r[j1 + j * e2r_dim1] = -fjac[j1 + j * fjac_dim1];
    e2i[j1 + j * e2i_dim1] = 0.;
  }
  i_1 = *n;
  for (j = 1; j <= i_1; ++j) {
    i_2 = j;
    for (i_ = 1; i_ <= i_2; ++i_) {
      e2i[i_ + j * e2i_dim1] = 0.;
      e2r[i_ + j * e2r_dim1] = -fjac[i_ + j * fjac_dim1];
    }
    e2r[j + j * e2r_dim1] += *alphn;
    e2i[j + j * e2i_dim1] = *betan;
  }
  dechc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &ip2[1], ier);
  return 0;

L55:
  return 0;
}

int slvrar(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, int *m1, int *m2,
           int *nm1, double *fac1, double *e1, int *lde1, double *z1,
           double *f1, int *ip1, int *iphes, int *ier, int *ijob,
           StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset, i__1,
      i__2, i__3, i__4, i__5, i__6;

  /* Local variables */
  int i__, j, k;
  double s1;
  int mm, mp, im1, mp1, jkm;
  double sum1;
  double zsafe;

  /* Parameter adjustments */
  --iphes;
  --f1;
  --z1;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip1;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L13;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
L48:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum1 = (z1[jkm] + sum1) / *fac1;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
      }
    }
  }
  sol(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
L49:
  for (i__ = *m1; i__ >= 1; --i__) {
    z1[i__] = (z1[i__] + z1[*m2 + i__]) / *fac1;
  }
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
L45:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum1 = (z1[jkm] + sum1) / *fac1;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        z1[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum1;
      }
    }
  }
  solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
       &ip1[1]);
  goto L49;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
    }
    z1[i__] += s1 * *fac1;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s1 = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__2 = std::min(i__5, i__6);
    for (j = std::max(i__3, i__4); j <= i__2; ++j) {
      s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1];
    }
    z1[im1] += s1 * *fac1;
  }
  if (*ijob == 14) {
    goto L45;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    /* Computing std::max */
    i__2 = 1, i__3 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__4 = std::min(i__5, i__6);
    for (j = std::max(i__2, i__3); j <= i__4; ++j) {
      s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
    }
    z1[i__] += s1 * *fac1;
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
  return 0;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    i__4 = *n;
    for (j = 1; j <= i__4; ++j) {
      s1 -= fmas[i__ + j * fmas_dim1] * f1[j];
    }
    z1[i__] += s1 * *fac1;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s1 = 0.;
    i__4 = *nm1;
    for (j = 1; j <= i__4; ++j) {
      s1 -= fmas[i__ + j * fmas_dim1] * f1[j + *m1];
    }
    z1[im1] += s1 * *fac1;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    z1[i__] -= f1[i__] * *fac1;
  }
  for (mm = *n - 2; mm >= 1; --mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L746;
    }
    zsafe = z1[mp];
    z1[mp] = z1[i__];
    z1[i__] = zsafe;
  L746:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      z1[i__] -= fjac[i__ + mp1 * fjac_dim1] * z1[mp];
    }
  }
  solh(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
  i__1 = *n - 2;
  for (mm = 1; mm <= i__1; ++mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__4 = *n;
    for (i__ = mp + 1; i__ <= i__4; ++i__) {
      z1[i__] += fjac[i__ + mp1 * fjac_dim1] * z1[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L750;
    }
    zsafe = z1[mp];
    z1[mp] = z1[i__];
    z1[i__] = zsafe;
  L750:;
  }
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
} /* slvrar_ */

int slvrai(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, int *m1, int *m2,
           int *nm1, double *alphn, double *betan, double *e2r, double *e2i,
           int *lde1, double *z2, double *z3, double *f2, double *f3,
           double *cont, int *ip2, int *iphes, int *ier, int *ijob,
           StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, e2r_offset,
      e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5, i__6;
  double d__1, d__2;

  /* Local variables */
  int i__, j, k;
  double s2, s3, bb;
  int mm, mp, im1, jm1, mp1;
  double z2i, z3i;
  int jkm, mpi;
  double sum2, sum3, abno;
  int iimu;
  double sumh, e1imp;
  double zsafe;

  /* Parameter adjustments */
  --iphes;
  --f3;
  --f2;
  --z3;
  --z2;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip2;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e2i_dim1 = *lde1;
  e2i_offset = 1 + e2i_dim1;
  e2i -= e2i_offset;
  e2r_dim1 = *lde1;
  e2r_offset = 1 + e2r_dim1;
  e2r -= e2r_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L13;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
L48:
  /* Computing 2nd power */
  d__1 = *alphn;
  /* Computing 2nd power */
  d__2 = *betan;
  abno = d__1 * d__1 + d__2 * d__2;
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum2 = 0.;
    sum3 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sumh = (z2[jkm] + sum2) / abno;
      sum3 = (z3[jkm] + sum3) / abno;
      sum2 = sumh * *alphn + sum3 * *betan;
      sum3 = sum3 * *alphn - sumh * *betan;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
        z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
      }
    }
  }
  solc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1],
       &z3[*m1 + 1], &ip2[1]);
L49:
  for (i__ = *m1; i__ >= 1; --i__) {
    mpi = *m2 + i__;
    z2i = z2[i__] + z2[mpi];
    z3i = z3[i__] + z3[mpi];
    z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
    z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
  }
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
L45:
  /* Computing 2nd power */
  d__1 = *alphn;
  /* Computing 2nd power */
  d__2 = *betan;
  abno = d__1 * d__1 + d__2 * d__2;
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum2 = 0.;
    sum3 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sumh = (z2[jkm] + sum2) / abno;
      sum3 = (z3[jkm] + sum3) / abno;
      sum2 = sumh * *alphn + sum3 * *betan;
      sum3 = sum3 * *alphn - sumh * *betan;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        iimu = i__ + *mujac + 1 - j;
        z2[im1] += fjac[iimu + jkm * fjac_dim1] * sum2;
        z3[im1] += fjac[iimu + jkm * fjac_dim1] * sum3;
      }
    }
  }
  solbc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle,
        &linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
  goto L49;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__2 = std::min(i__5, i__6);
    for (j = std::max(i__3, i__4); j <= i__2; ++j) {
      jm1 = j + *m1;
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s2 -= bb * f2[jm1];
      s3 -= bb * f3[jm1];
    }
    z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
    z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
  }
  if (*ijob == 14) {
    goto L45;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__2 = 1, i__3 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__4 = std::min(i__5, i__6);
    for (j = std::max(i__2, i__3); j <= i__4; ++j) {
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = 0.;
    s3 = 0.;
    i__4 = *n;
    for (j = 1; j <= i__4; ++j) {
      bb = fmas[i__ + j * fmas_dim1];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s2 = 0.;
    s3 = 0.;
    i__4 = *nm1;
    for (j = 1; j <= i__4; ++j) {
      jm1 = j + *m1;
      bb = fmas[i__ + j * fmas_dim1];
      s2 -= bb * f2[jm1];
      s3 -= bb * f3[jm1];
    }
    z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
    z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  for (mm = *n - 2; mm >= 1; --mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L746;
    }
    zsafe = z2[mp];
    z2[mp] = z2[i__];
    z2[i__] = zsafe;
    zsafe = z3[mp];
    z3[mp] = z3[i__];
    z3[i__] = zsafe;
  L746:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      e1imp = fjac[i__ + mp1 * fjac_dim1];
      z2[i__] -= e1imp * z2[mp];
      z3[i__] -= e1imp * z3[mp];
    }
  }
  solhc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
        &ip2[1]);
  i__1 = *n - 2;
  for (mm = 1; mm <= i__1; ++mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__4 = *n;
    for (i__ = mp + 1; i__ <= i__4; ++i__) {
      e1imp = fjac[i__ + mp1 * fjac_dim1];
      z2[i__] += e1imp * z2[mp];
      z3[i__] += e1imp * z3[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L750;
    }
    zsafe = z2[mp];
    z2[mp] = z2[i__];
    z2[i__] = zsafe;
    zsafe = z3[mp];
    z3[mp] = z3[i__];
    z3[i__] = zsafe;
  L750:;
  }
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
} /* slvrai_ */

int slvrad(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, int *m1, int *m2,
           int *nm1, double *fac1, double *alphn, double *betan, double *e1,
           double *e2r, double *e2i, int *lde1, double *z1, double *z2,
           double *z3, double *f1, double *f2, double *f3, double *cont,
           int *ip1, int *ip2, int *iphes, int *ier, int *ijob,
           StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset,
      e2r_dim1, e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5,
      i__6;
  double d__1, d__2;

  /* Local variables */
  int i__, j, k;
  double s1, s2, s3, bb;
  int mm, mp, j1b, j2b, im1, jm1, mp1;
  double z2i, z3i;
  int jkm, mpi;
  double sum1, sum2, sum3, ffja, abno;
  double sumh, e1imp;
  double zsafe;

  /* Parameter adjustments */
  --iphes;
  --f3;
  --f2;
  --f1;
  --z3;
  --z2;
  --z1;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip2;
  --ip1;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e2i_dim1 = *lde1;
  e2i_offset = 1 + e2i_dim1;
  e2i -= e2i_offset;
  e2r_dim1 = *lde1;
  e2r_offset = 1 + e2r_dim1;
  e2r -= e2r_offset;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L13;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
L48:
  /* Computing 2nd power */
  d__1 = *alphn;
  /* Computing 2nd power */
  d__2 = *betan;
  abno = d__1 * d__1 + d__2 * d__2;
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    sum2 = 0.;
    sum3 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum1 = (z1[jkm] + sum1) / *fac1;
      sumh = (z2[jkm] + sum2) / abno;
      sum3 = (z3[jkm] + sum3) / abno;
      sum2 = sumh * *alphn + sum3 * *betan;
      sum3 = sum3 * *alphn - sumh * *betan;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
        z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
        z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
      }
    }
  }
  sol(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
  solc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1],
       &z3[*m1 + 1], &ip2[1]);
L49:
  for (i__ = *m1; i__ >= 1; --i__) {
    mpi = *m2 + i__;
    z1[i__] = (z1[i__] + z1[mpi]) / *fac1;
    z2i = z2[i__] + z2[mpi];
    z3i = z3[i__] + z3[mpi];
    z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
    z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
  }
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
  solbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
L45:
  /* Computing 2nd power */
  d__1 = *alphn;
  /* Computing 2nd power */
  d__2 = *betan;
  abno = d__1 * d__1 + d__2 * d__2;
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    sum2 = 0.;
    sum3 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum1 = (z1[jkm] + sum1) / *fac1;
      sumh = (z2[jkm] + sum2) / abno;
      sum3 = (z3[jkm] + sum3) / abno;
      sum2 = sumh * *alphn + sum3 * *betan;
      sum3 = sum3 * *alphn - sumh * *betan;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        ffja = fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1];
        z1[im1] += ffja * sum1;
        z2[im1] += ffja * sum2;
        z3[im1] += ffja * sum3;
      }
    }
  }
  solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
       &ip1[1]);
  solbc(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle,
        &linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
  goto L49;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s1 -= bb * f1[j];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z1[i__] += s1 * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    j1b = std::max(i__3, i__4);
    /* Computing std::min */
    i__3 = *nm1, i__4 = i__ + *mumas;
    j2b = std::min(i__3, i__4);
    i__3 = j2b;
    for (j = j1b; j <= i__3; ++j) {
      jm1 = j + *m1;
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s1 -= bb * f1[jm1];
      s2 -= bb * f2[jm1];
      s3 -= bb * f3[jm1];
    }
    z1[im1] += s1 * *fac1;
    z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
    z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
  }
  if (*ijob == 14) {
    goto L45;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__2 = std::min(i__5, i__6);
    for (j = std::max(i__3, i__4); j <= i__2; ++j) {
      bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
      s1 -= bb * f1[j];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z1[i__] += s1 * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
  solbc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue,
        &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
      bb = fmas[i__ + j * fmas_dim1];
      s1 -= bb * f1[j];
      s2 -= bb * f2[j];
      s3 -= bb * f3[j];
    }
    z1[i__] += s1 * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  sol(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
  solc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    im1 = i__ + *m1;
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    i__2 = *nm1;
    for (j = 1; j <= i__2; ++j) {
      jm1 = j + *m1;
      bb = fmas[i__ + j * fmas_dim1];
      s1 -= bb * f1[jm1];
      s2 -= bb * f2[jm1];
      s3 -= bb * f3[jm1];
    }
    z1[im1] += s1 * *fac1;
    z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
    z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    s2 = -f2[i__];
    s3 = -f3[i__];
    z1[i__] -= f1[i__] * *fac1;
    z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
    z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
  }
  for (mm = *n - 2; mm >= 1; --mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L746;
    }
    zsafe = z1[mp];
    z1[mp] = z1[i__];
    z1[i__] = zsafe;
    zsafe = z2[mp];
    z2[mp] = z2[i__];
    z2[i__] = zsafe;
    zsafe = z3[mp];
    z3[mp] = z3[i__];
    z3[i__] = zsafe;
  L746:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      e1imp = fjac[i__ + mp1 * fjac_dim1];
      z1[i__] -= e1imp * z1[mp];
      z2[i__] -= e1imp * z2[mp];
      z3[i__] -= e1imp * z3[mp];
    }
  }
  solh(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
  solhc(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
        &ip2[1]);
  i__1 = *n - 2;
  for (mm = 1; mm <= i__1; ++mm) {
    mp = *n - mm;
    mp1 = mp - 1;
    i__2 = *n;
    for (i__ = mp + 1; i__ <= i__2; ++i__) {
      e1imp = fjac[i__ + mp1 * fjac_dim1];
      z1[i__] += e1imp * z1[mp];
      z2[i__] += e1imp * z2[mp];
      z3[i__] += e1imp * z3[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L750;
    }
    zsafe = z1[mp];
    z1[mp] = z1[i__];
    z1[i__] = zsafe;
    zsafe = z2[mp];
    z2[mp] = z2[i__];
    z2[i__] = zsafe;
    zsafe = z3[mp];
    z3[mp] = z3[i__];
    z3[i__] = zsafe;
  L750:;
  }
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
} /* slvrad_ */

int estrad(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, double *h__,
           double *dd1, double *dd2, double *dd3, F_fcn fcn, int *nfcn,
           double *y0, double *y, int *ijob, double *x, int *m1, int *m2,
           int *nm1, double *e1, int *lde1, double *z1, double *z2, double *z3,
           double *cont, double *f1, double *f2, int *ip1, int *iphes,
           double *scal, double *err, bool *first, bool *reject, double *fac1,
           StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset, i__1,
      i__2, i__3, i__4, i__5, i__6;
  double d__1;

  /* Local variables */
  int i__, j, k, mm, mp, im1;
  double sum, hee1, hee2, hee3, sum1;
  double zsafe;

  /* Parameter adjustments */
  --scal;
  --iphes;
  --f2;
  --f1;
  --cont;
  --z3;
  --z2;
  --z1;
  --y;
  --y0;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip1;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;

  /* Function Body */
  hee1 = *dd1 / *h__;
  hee2 = *dd2 / *h__;
  hee3 = *dd3 / *h__;
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L14;
  case 15:
    goto L15;
  }

L1:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L11:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
L48:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      sum1 = (cont[j + k * *m2] + sum1) / *fac1;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
      }
    }
  }
  sol(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
  }
  goto L77;

L2:
  /* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
  goto L77;

L12:
  /* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
L45:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      sum1 = (cont[j + k * *m2] + sum1) / *fac1;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        cont[im1] +=
            fjac[i__ + *mujac + 1 - j + (j + k * *m2) * fjac_dim1] * sum1;
      }
    }
  }
  solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1],
       &ip1[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
  }
  goto L77;

L3:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
    }
    f2[i__] = sum;
    cont[i__] = sum + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L13:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__2 = std::min(i__5, i__6);
    for (j = std::max(i__3, i__4); j <= i__2; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1];
    }
    im1 = i__ + *m1;
    f2[im1] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L48;

L4:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__2 = 1, i__3 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__4 = std::min(i__5, i__6);
    for (j = std::max(i__2, i__3); j <= i__4; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
    }
    f2[i__] = sum;
    cont[i__] = sum + y0[i__];
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
  goto L77;

L14:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1];
    }
    im1 = i__ + *m1;
    f2[im1] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L45;

L5:
  /* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *n;
    for (j = 1; j <= i__3; ++j) {
      sum += fmas[i__ + j * fmas_dim1] * f1[j];
    }
    f2[i__] = sum;
    cont[i__] = sum + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L15:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *nm1;
    for (j = 1; j <= i__3; ++j) {
      sum += fmas[i__ + j * fmas_dim1] * f1[j + *m1];
    }
    im1 = i__ + *m1;
    f2[im1] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L48;

L6:
  /* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ------  THIS OPTION IS NOT PROVIDED */
  return 0;

L7:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    cont[i__] = f2[i__] + y0[i__];
  }
  for (mm = *n - 2; mm >= 1; --mm) {
    mp = *n - mm;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L310;
    }
    zsafe = cont[mp];
    cont[mp] = cont[i__];
    cont[i__] = zsafe;
  L310:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
    }
  }
  solh(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
  i__1 = *n - 2;
  for (mm = 1; mm <= i__1; ++mm) {
    mp = *n - mm;
    i__3 = *n;
    for (i__ = mp + 1; i__ <= i__3; ++i__) {
      cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L440;
    }
    zsafe = cont[mp];
    cont[mp] = cont[i__];
    cont[i__] = zsafe;
  L440:;
  }

  /* -------------------------------------- */

L77:
  *err = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing 2nd power */
    d__1 = cont[i__] / scal[i__];
    *err += d__1 * d__1;
  }
  /* Computing std::max */
  d__1 = sqrt(*err / *n);
  *err = std::max(d__1, 1e-10);

  if (*err < 1.) {
    return 0;
  }
  if (*first || *reject) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + cont[i__];
    }
    fcn(n, x, &cont[1], &f1[1]);
    ++(*nfcn);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = f1[i__] + f2[i__];
    }
    switch (*ijob) {
    case 1:
      goto L31;
    case 2:
      goto L32;
    case 3:
      goto L31;
    case 4:
      goto L32;
    case 5:
      goto L31;
    case 6:
      goto L32;
    case 7:
      goto L33;
    case 8:
      goto L55;
    case 9:
      goto L55;
    case 10:
      goto L55;
    case 11:
      goto L41;
    case 12:
      goto L42;
    case 13:
      goto L41;
    case 14:
      goto L42;
    case 15:
      goto L41;
    }
  /* ------ FULL MATRIX OPTION */
  L31:
    sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L88;
  /* ------ FULL MATRIX OPTION, SECOND ORDER */
  L41:
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
      sum1 = 0.;
      for (k = mm - 1; k >= 0; --k) {
        sum1 = (cont[j + k * *m2] + sum1) / *fac1;
        i__3 = *nm1;
        for (i__ = 1; i__ <= i__3; ++i__) {
          im1 = i__ + *m1;
          cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
        }
      }
    }
    sol(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
      cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L88;
  /* ------ BANDED MATRIX OPTION */
  L32:
    solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1],
         &ip1[1]);
    goto L88;
  /* ------ BANDED MATRIX OPTION, SECOND ORDER */
  L42:
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
      sum1 = 0.;
      for (k = mm - 1; k >= 0; --k) {
        sum1 = (cont[j + k * *m2] + sum1) / *fac1;
        /* Computing std::max */
        i__3 = 1, i__4 = j - *mujac;
        /* Computing std::min */
        i__5 = *nm1, i__6 = j + *mljac;
        i__2 = std::min(i__5, i__6);
        for (i__ = std::max(i__3, i__4); i__ <= i__2; ++i__) {
          im1 = i__ + *m1;
          cont[im1] +=
              fjac[i__ + *mujac + 1 - j + (j + k * *m2) * fjac_dim1] * sum1;
        }
      }
    }
    solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1],
         &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
      cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L88;
  /* ------ HESSENBERG MATRIX OPTION */
  L33:
    for (mm = *n - 2; mm >= 1; --mm) {
      mp = *n - mm;
      i__ = iphes[mp];
      if (i__ == mp) {
        goto L510;
      }
      zsafe = cont[mp];
      cont[mp] = cont[i__];
      cont[i__] = zsafe;
    L510:
      i__1 = *n;
      for (i__ = mp + 1; i__ <= i__1; ++i__) {
        cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
      }
    }
    solh(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
      mp = *n - mm;
      i__2 = *n;
      for (i__ = mp + 1; i__ <= i__2; ++i__) {
        cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
      }
      i__ = iphes[mp];
      if (i__ == mp) {
        goto L640;
      }
      zsafe = cont[mp];
      cont[mp] = cont[i__];
      cont[i__] = zsafe;
    L640:;
    }
  /* ----------------------------------- */
  L88:
    *err = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      /* Computing 2nd power */
      d__1 = cont[i__] / scal[i__];
      *err += d__1 * d__1;
    }
    /* Computing std::max */
    d__1 = sqrt(*err / *n);
    *err = std::max(d__1, 1e-10);
  }
  return 0;
  /* ----------------------------------------------------------- */
L55:
  return 0;
}

int estrav(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, double *h__,
           double *dd, F_fcn fcn, int *nfcn, double *y0, double *y, int *ijob,
           double *x, int *m1, int *m2, int *nm1, int *ns, int *nns, double *e1,
           int *lde1, double *zz, double *cont, double *ff, int *ip1,
           int *iphes, double *scal, double *err, bool *first, bool *reject,
           double *fac1, StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset, i__1,
      i__2, i__3, i__4, i__5, i__6;
  double d__1;

  /* Local variables */
  int i__, j, k, mm, mp, im1;
  double sum, sum1;
  double zsafe;

  /* Parameter adjustments */
  --scal;
  --iphes;
  --cont;
  --y;
  --y0;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip1;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  --dd;
  --ff;
  --zz;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L14;
  case 15:
    goto L15;
  }

L1:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__2 = *ns;
    for (k = 1; k <= i__2; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L11:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__2 = *ns;
    for (k = 1; k <= i__2; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
L48:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      sum1 = (cont[j + k * *m2] + sum1) / *fac1;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
      }
    }
  }
  sol(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
  }
  goto L77;

L2:
  /* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__2 = *ns;
    for (k = 1; k <= i__2; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
  goto L77;

L12:
  /* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__2 = *ns;
    for (k = 1; k <= i__2; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
L45:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum1 = 0.;
    for (k = mm - 1; k >= 0; --k) {
      sum1 = (cont[j + k * *m2] + sum1) / *fac1;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        cont[im1] +=
            fjac[i__ + *mujac + 1 - j + (j + k * *m2) * fjac_dim1] * sum1;
      }
    }
  }
  solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1],
       &ip1[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
  }
  goto L77;

L3:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__4 = *ns;
    for (k = 1; k <= i__4; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
    }
    ff[i__ + *n] = sum;
    cont[i__] = sum + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L13:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__3 = 1, i__4 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__2 = std::min(i__5, i__6);
    for (j = std::max(i__3, i__4); j <= i__2; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *m1];
    }
    im1 = i__ + *m1;
    ff[im1 + *n] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L48;

L4:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__2 = *ns;
    for (k = 1; k <= i__2; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__2 = 1, i__3 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *n, i__6 = i__ + *mumas;
    i__4 = std::min(i__5, i__6);
    for (j = std::max(i__2, i__3); j <= i__4; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
    }
    ff[i__ + *n] = sum;
    cont[i__] = sum + y0[i__];
  }
  solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
  goto L77;

L14:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__4 = *ns;
    for (k = 1; k <= i__4; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__4 = *ns;
    for (k = 1; k <= i__4; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    /* Computing std::max */
    i__4 = 1, i__2 = i__ - *mlmas;
    /* Computing std::min */
    i__5 = *nm1, i__6 = i__ + *mumas;
    i__3 = std::min(i__5, i__6);
    for (j = std::max(i__4, i__2); j <= i__3; ++j) {
      sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *m1];
    }
    im1 = i__ + *m1;
    ff[im1 + *n] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L45;

L5:
  /* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *n;
    for (j = 1; j <= i__3; ++j) {
      sum += fmas[i__ + j * fmas_dim1] * ff[j];
    }
    ff[i__ + *n] = sum;
    cont[i__] = sum + y0[i__];
  }
  sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
  goto L77;

L15:
  /* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  i__1 = *m1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  i__1 = *n;
  for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__] = sum / *h__;
  }
  i__1 = *nm1;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *nm1;
    for (j = 1; j <= i__3; ++j) {
      sum += fmas[i__ + j * fmas_dim1] * ff[j + *m1];
    }
    im1 = i__ + *m1;
    ff[im1 + *n] = sum;
    cont[im1] = sum + y0[im1];
  }
  goto L48;

L6:
  /* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ------  THIS OPTION IS NOT PROVIDED */
  return 0;

L7:
  /* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    sum = 0.;
    i__3 = *ns;
    for (k = 1; k <= i__3; ++k) {
      sum += dd[k] * zz[i__ + (k - 1) * *n];
    }
    ff[i__ + *n] = sum / *h__;
    cont[i__] = ff[i__ + *n] + y0[i__];
  }
  for (mm = *n - 2; mm >= 1; --mm) {
    mp = *n - mm;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L310;
    }
    zsafe = cont[mp];
    cont[mp] = cont[i__];
    cont[i__] = zsafe;
  L310:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
    }
  }
  solh(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
  i__1 = *n - 2;
  for (mm = 1; mm <= i__1; ++mm) {
    mp = *n - mm;
    i__3 = *n;
    for (i__ = mp + 1; i__ <= i__3; ++i__) {
      cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L440;
    }
    zsafe = cont[mp];
    cont[mp] = cont[i__];
    cont[i__] = zsafe;
  L440:;
  }

  /* -------------------------------------- */

L77:
  *err = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing 2nd power */
    d__1 = cont[i__] / scal[i__];
    *err += d__1 * d__1;
  }
  /* Computing std::max */
  d__1 = sqrt(*err / *n);
  *err = std::max(d__1, 1e-10);

  if (*err < 1.) {
    return 0;
  }
  if (*first || *reject) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + cont[i__];
    }
    fcn(n, x, &cont[1], &ff[1]);
    ++(*nfcn);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = ff[i__] + ff[i__ + *n];
    }
    switch (*ijob) {
    case 1:
      goto L31;
    case 2:
      goto L32;
    case 3:
      goto L31;
    case 4:
      goto L32;
    case 5:
      goto L31;
    case 6:
      goto L32;
    case 7:
      goto L33;
    case 8:
      goto L55;
    case 9:
      goto L55;
    case 10:
      goto L55;
    case 11:
      goto L41;
    case 12:
      goto L42;
    case 13:
      goto L41;
    case 14:
      goto L42;
    case 15:
      goto L41;
    }
  /* ------ FULL MATRIX OPTION */
  L31:
    sol(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L88;
  /* ------ FULL MATRIX OPTION, SECOND ORDER */
  L41:
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
      sum1 = 0.;
      for (k = mm - 1; k >= 0; --k) {
        sum1 = (cont[j + k * *m2] + sum1) / *fac1;
        i__3 = *nm1;
        for (i__ = 1; i__ <= i__3; ++i__) {
          im1 = i__ + *m1;
          cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
        }
      }
    }
    sol(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
      cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L88;
  /* ------ BANDED MATRIX OPTION */
  L32:
    solb(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1],
         &ip1[1]);
    goto L88;
  /* ------ BANDED MATRIX OPTION, SECOND ORDER */
  L42:
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
      sum1 = 0.;
      for (k = mm - 1; k >= 0; --k) {
        sum1 = (cont[j + k * *m2] + sum1) / *fac1;
        /* Computing std::max */
        i__3 = 1, i__4 = j - *mujac;
        /* Computing std::min */
        i__5 = *nm1, i__6 = j + *mljac;
        i__2 = std::min(i__5, i__6);
        for (i__ = std::max(i__3, i__4); i__ <= i__2; ++i__) {
          im1 = i__ + *m1;
          cont[im1] +=
              fjac[i__ + *mujac + 1 - j + (j + k * *m2) * fjac_dim1] * sum1;
        }
      }
    }
    solb(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1],
         &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
      cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L88;
  /* ------ HESSENBERG MATRIX OPTION */
  L33:
    for (mm = *n - 2; mm >= 1; --mm) {
      mp = *n - mm;
      i__ = iphes[mp];
      if (i__ == mp) {
        goto L510;
      }
      zsafe = cont[mp];
      cont[mp] = cont[i__];
      cont[i__] = zsafe;
    L510:
      i__1 = *n;
      for (i__ = mp + 1; i__ <= i__1; ++i__) {
        cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
      }
    }
    solh(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
      mp = *n - mm;
      i__2 = *n;
      for (i__ = mp + 1; i__ <= i__2; ++i__) {
        cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
      }
      i__ = iphes[mp];
      if (i__ == mp) {
        goto L640;
      }
      zsafe = cont[mp];
      cont[mp] = cont[i__];
      cont[i__] = zsafe;
    L640:;
    }
  /* ----------------------------------- */
  L88:
    *err = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      /* Computing 2nd power */
      d__1 = cont[i__] / scal[i__];
      *err += d__1 * d__1;
    }
    /* Computing std::max */
    d__1 = sqrt(*err / *n);
    *err = std::max(d__1, 1e-10);
  }
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
}

int slvrod(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, int *m1, int *m2,
           int *nm1, double *fac1, double *e, int *lde, int *ip, double *dy,
           double *ak, double *fx, double *ynew, double *hd, int *ijob,
           bool *stage1, StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, i__1,
      i__2, i__3, i__4, i__5, i__6;

  /* Local variables */
  int i__, j, k, mm, im1, jkm;
  double sum;

  /* Parameter adjustments */
  --ynew;
  --fx;
  --ak;
  --dy;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e_dim1 = *lde;
  e_offset = 1 + e_dim1;
  e -= e_offset;

  /* Function Body */
  if (*hd == 0.) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] = dy[i__];
    }
  } else {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] = dy[i__] + *hd * fx[i__];
    }
  }

  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L3;
  case 4:
    goto L4;
  case 5:
    goto L5;
  case 6:
    goto L6;
  case 7:
    goto L55;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L13;
  case 14:
    goto L13;
  case 15:
    goto L15;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
  }
  sol(n, lde, &e[e_offset], &ak[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
  }
L48:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum = (ak[jkm] + sum) / *fac1;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        ak[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
      }
    }
  }
  sol(nm1, lde, &e[e_offset], &ak[*m1 + 1], &ip[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
  }
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
  }
  solb(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
  }
L45:
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum = (ak[jkm] + sum) / *fac1;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        ak[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum;
      }
    }
  }
  solb(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[*m1 + 1],
       &ip[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
  }
  return 0;

  /* ----------------------------------------------------------- */

L3:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      /* Computing std::max */
      i__4 = 1, i__2 = i__ - *mlmas;
      /* Computing std::min */
      i__5 = *n, i__6 = i__ + *mumas;
      i__3 = std::min(i__5, i__6);
      for (j = std::max(i__4, i__2); j <= i__3; ++j) {
        sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[j];
      }
      ak[i__] += sum;
    }
  }
  sol(n, lde, &e[e_offset], &ak[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L13:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  if (*stage1) {
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      /* Computing std::max */
      i__3 = 1, i__4 = i__ - *mlmas;
      /* Computing std::min */
      i__5 = *nm1, i__6 = i__ + *mumas;
      i__2 = std::min(i__5, i__6);
      for (j = std::max(i__3, i__4); j <= i__2; ++j) {
        sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[j + *m1];
      }
      im1 = i__ + *m1;
      ak[im1] += sum;
    }
  }
  if (*ijob == 14) {
    goto L45;
  }
  goto L48;

  /* ----------------------------------------------------------- */

L4:
  /* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      /* Computing std::max */
      i__2 = 1, i__3 = i__ - *mlmas;
      /* Computing std::min */
      i__5 = *n, i__6 = i__ + *mumas;
      i__4 = std::min(i__5, i__6);
      for (j = std::max(i__2, i__3); j <= i__4; ++j) {
        sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[j];
      }
      ak[i__] += sum;
    }
  }
  solb(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L5:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      i__4 = *n;
      for (j = 1; j <= i__4; ++j) {
        sum += fmas[i__ + j * fmas_dim1] * ynew[j];
      }
      ak[i__] += sum;
    }
  }
  sol(n, lde, &e[e_offset], &ak[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L15:
  /* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
  if (*stage1) {
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
      ak[i__] += ynew[i__];
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      i__4 = *nm1;
      for (j = 1; j <= i__4; ++j) {
        sum += fmas[i__ + j * fmas_dim1] * ynew[j + *m1];
      }
      im1 = i__ + *m1;
      ak[im1] += sum;
    }
  }
  goto L48;

  /* ----------------------------------------------------------- */

L6:
  /* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
  /* ---  THIS OPTION IS NOT PROVIDED */
  if (*stage1) {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = 0.;
      i__4 = *n;
      for (j = 1; j <= i__4; ++j) {
        /* L623: */
        sum += fmas[i__ + j * fmas_dim1] * ynew[j];
      }
      /* L624: */
      ak[i__] += sum;
    }
    solb(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
  }
  return 0;

L55:
  return 0;
}

int slvseu(int *n, double *fjac, int *ldjac, int *mljac, int *mujac,
           double *fmas, int *ldmas, int *mlmas, int *mumas, int *m1, int *m2,
           int *nm1, double *fac1, double *e, int *lde, int *ip, int *iphes,
           double *del, int *ijob, StiffLinAlg &linal_1) {
  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, i__1,
      i__2, i__3, i__4, i__5, i__6;

  /* Local variables */
  int i__, j, k, mm, mp, im1, mp1, jkm, mmm;
  double sum;
  double zsafe;

  /* Parameter adjustments */
  --del;
  --iphes;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  --ip;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;
  e_dim1 = *lde;
  e_offset = 1 + e_dim1;
  e -= e_offset;

  /* Function Body */
  switch (*ijob) {
  case 1:
    goto L1;
  case 2:
    goto L2;
  case 3:
    goto L1;
  case 4:
    goto L2;
  case 5:
    goto L1;
  case 6:
    goto L55;
  case 7:
    goto L7;
  case 8:
    goto L55;
  case 9:
    goto L55;
  case 10:
    goto L55;
  case 11:
    goto L11;
  case 12:
    goto L12;
  case 13:
    goto L11;
  case 14:
    goto L12;
  case 15:
    goto L11;
  }

  /* ----------------------------------------------------------- */

L1:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
  sol(n, lde, &e[e_offset], &del[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L11:
  /* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum = (del[jkm] + sum) / *fac1;
      i__2 = *nm1;
      for (i__ = 1; i__ <= i__2; ++i__) {
        im1 = i__ + *m1;
        del[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
      }
    }
  }
  sol(nm1, lde, &e[e_offset], &del[*m1 + 1], &ip[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
  }
  return 0;

  /* ----------------------------------------------------------- */

L2:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
  solb(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[1], &ip[1]);
  return 0;

  /* ----------------------------------------------------------- */

L12:
  /* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
  mm = *m1 / *m2;
  i__1 = *m2;
  for (j = 1; j <= i__1; ++j) {
    sum = 0.;
    for (k = mm - 1; k >= 0; --k) {
      jkm = j + k * *m2;
      sum = (del[jkm] + sum) / *fac1;
      /* Computing std::max */
      i__2 = 1, i__3 = j - *mujac;
      /* Computing std::min */
      i__5 = *nm1, i__6 = j + *mljac;
      i__4 = std::min(i__5, i__6);
      for (i__ = std::max(i__2, i__3); i__ <= i__4; ++i__) {
        im1 = i__ + *m1;
        del[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum;
      }
    }
  }
  solb(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[*m1 + 1],
       &ip[1]);
  for (i__ = *m1; i__ >= 1; --i__) {
    del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
  }
  return 0;

  /* ----------------------------------------------------------- */

L7:
  /* ---  HESSENBERG OPTION */
  for (mmm = *n - 2; mmm >= 1; --mmm) {
    mp = *n - mmm;
    mp1 = mp - 1;
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L110;
    }
    zsafe = del[mp];
    del[mp] = del[i__];
    del[i__] = zsafe;
  L110:
    i__1 = *n;
    for (i__ = mp + 1; i__ <= i__1; ++i__) {
      del[i__] -= fjac[i__ + mp1 * fjac_dim1] * del[mp];
    }
  }
  solh(n, lde, &e[e_offset], &c__1, &del[1], &ip[1]);
  i__1 = *n - 2;
  for (mmm = 1; mmm <= i__1; ++mmm) {
    mp = *n - mmm;
    mp1 = mp - 1;
    i__4 = *n;
    for (i__ = mp + 1; i__ <= i__4; ++i__) {
      del[i__] += fjac[i__ + mp1 * fjac_dim1] * del[mp];
    }
    i__ = iphes[mp];
    if (i__ == mp) {
      goto L240;
    }
    zsafe = del[mp];
    del[mp] = del[i__];
    del[i__] = zsafe;
  L240:;
  }
  return 0;

  /* ----------------------------------------------------------- */

L55:
  return 0;
} /* slvseu_ */

} // namespace stiff

#endif // DARKSUN_DC_DECSOL_HPP
