//
// Created by logan on 8/10/20.
//

#ifndef DARKSUN_STIFF_DECSOL_HPP
#define DARKSUN_STIFF_DECSOL_HPP

#include <algorithm>
#include <cmath>

namespace stiff {

int dec(const int *n, const int *ndim, double *a, int *ip, int *ier) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3;
  double d_1, d_2;

  /* Local variables */
  int i_, j, k, m;
  double t;
  int nm1, kp1;

  --ip;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  if (*n == 1) {
    goto L70;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = k;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      if ((d_1 = a[i_ + k * a_dim1], std::abs(d_1)) >
          (d_2 = a[m + k * a_dim1], std::abs(d_2))) {
        m = i_;
      }
    }
    ip[k] = m;
    t = a[m + k * a_dim1];
    if (m == k) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    a[m + k * a_dim1] = a[k + k * a_dim1];
    a[k + k * a_dim1] = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      /* L30: */
      a[i_ + k * a_dim1] = -a[i_ + k * a_dim1] * t;
    }
    i_2 = *n;
    for (j = kp1; j <= i_2; ++j) {
      t = a[m + j * a_dim1];
      a[m + j * a_dim1] = a[k + j * a_dim1];
      a[k + j * a_dim1] = t;
      if (t == 0.) {
        goto L45;
      }
      i_3 = *n;
      for (i_ = kp1; i_ <= i_3; ++i_) {
        /* L40: */
        a[i_ + j * a_dim1] += a[i_ + k * a_dim1] * t;
      }
    L45:;
    }
  }
L70:
  k = *n;
  if (a[*n + *n * a_dim1] == 0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
}

int sol(const int *n, const int *ndim, const double *a, double *b,
        const int *ip) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2;

  /* Local variables */
  int i_, k, m;
  double t;
  int kb, km1, nm1, kp1;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    B = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --b;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  if (*n == 1) {
    goto L50;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = ip[k];
    t = b[m];
    b[m] = b[k];
    b[k] = t;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      /* L10: */
      b[i_] += a[i_ + k * a_dim1] * t;
    }
    /* L20: */
  }
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    km1 = *n - kb;
    k = km1 + 1;
    b[k] /= a[k + k * a_dim1];
    t = -b[k];
    i_2 = km1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      /* L30: */
      b[i_] += a[i_ + k * a_dim1] * t;
    }
    /* L40: */
  }
L50:
  b[1] /= a[a_dim1 + 1];
  return 0;
  /* ----------------------- END OF SUBROUTINE SOL ------------------------- */
}

int dech(const int *n, const int *ndim, double *a, const int *lb, int *ip,
         int *ier) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3;
  double d_1, d_2;

  /* Local variables */
  int i_, j, k, m;
  double t;
  int na, nm1, kp1;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIstd::minATION OF A HESSENBERG */
  /*  MATRIX WITH LOWER BANDWIDTH LB */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX A. */
  /*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*     A = MATRIX TO BE TRIANGULARIZED. */
  /*     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1). */
  /*  OUTPUT.. */
  /*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
  /*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A SLIGHT MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  if (*n == 1) {
    goto L70;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = k;
    /* Computing std::min */
    i_2 = *n, i_3 = *lb + k;
    na = std::min(i_2, i_3);
    i_2 = na;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      if ((d_1 = a[i_ + k * a_dim1], std::abs(d_1)) >
          (d_2 = a[m + k * a_dim1], std::abs(d_2))) {
        m = i_;
      }
    }
    ip[k] = m;
    t = a[m + k * a_dim1];
    if (m == k) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    a[m + k * a_dim1] = a[k + k * a_dim1];
    a[k + k * a_dim1] = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i_2 = na;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      /* L30: */
      a[i_ + k * a_dim1] = -a[i_ + k * a_dim1] * t;
    }
    i_2 = *n;
    for (j = kp1; j <= i_2; ++j) {
      t = a[m + j * a_dim1];
      a[m + j * a_dim1] = a[k + j * a_dim1];
      a[k + j * a_dim1] = t;
      if (t == 0.) {
        goto L45;
      }
      i_3 = na;
      for (i_ = kp1; i_ <= i_3; ++i_) {
        /* L40: */
        a[i_ + j * a_dim1] += a[i_ + k * a_dim1] * t;
      }
    L45:;
    }
  }
L70:
  k = *n;
  if (a[*n + *n * a_dim1] == 0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
}

int solh(const int *n, const int *ndim, const double *a, const int *lb,
         double *b, const int *ip) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3;

  /* Local variables */
  int i_, k, m;
  double t;
  int kb, na, km1, nm1, kp1;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX A. */
  /*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH. */
  /*    LB = LOWER BANDWIDTH OF A. */
  /*    B = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DECH HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --b;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  if (*n == 1) {
    goto L50;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = ip[k];
    t = b[m];
    b[m] = b[k];
    b[k] = t;
    /* Computing std::min */
    i_2 = *n, i_3 = *lb + k;
    na = std::min(i_2, i_3);
    i_2 = na;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      /* L10: */
      b[i_] += a[i_ + k * a_dim1] * t;
    }
  }
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    km1 = *n - kb;
    k = km1 + 1;
    b[k] /= a[k + k * a_dim1];
    t = -b[k];
    i_2 = km1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      /* L30: */
      b[i_] += a[i_ + k * a_dim1] * t;
    }
    /* L40: */
  }
L50:
  b[1] /= a[a_dim1 + 1];
  return 0;
}

int decc(const int *n, const int *ndim, double *ar, double *ai, int *ip,
         int *ier) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2, i_3;
  double d_1, d_2, d_3, d_4;

  /* Local variables */
  int i_, j, k, m;
  double ti, tr;
  int nm1, kp1;
  double den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIstd::minATION */
  /*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
  /*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
  /*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
  /*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    REAL PART. */
  /*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    IMAGINARY PART. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  if (*n == 1) {
    goto L70;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = k;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      if ((d_1 = ar[i_ + k * ar_dim1], std::abs(d_1)) +
              (d_2 = ai[i_ + k * ai_dim1], std::abs(d_2)) >
          (d_3 = ar[m + k * ar_dim1], std::abs(d_3)) +
              (d_4 = ai[m + k * ai_dim1], std::abs(d_4))) {
        m = i_;
      }
    }
    ip[k] = m;
    tr = ar[m + k * ar_dim1];
    ti = ai[m + k * ai_dim1];
    if (m == k) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
    ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
    ar[k + k * ar_dim1] = tr;
    ai[k + k * ai_dim1] = ti;
  L20:
    if (std::abs(tr) + std::abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      ar[i_ + k * ar_dim1] = -prodr;
      ai[i_ + k * ai_dim1] = -prodi;
    }
    i_2 = *n;
    for (j = kp1; j <= i_2; ++j) {
      tr = ar[m + j * ar_dim1];
      ti = ai[m + j * ai_dim1];
      ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
      ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
      ar[k + j * ar_dim1] = tr;
      ai[k + j * ai_dim1] = ti;
      if (std::abs(tr) + std::abs(ti) == 0.) {
        goto L48;
      }
      if (ti == 0.) {
        i_3 = *n;
        for (i_ = kp1; i_ <= i_3; ++i_) {
          prodr = ar[i_ + k * ar_dim1] * tr;
          prodi = ai[i_ + k * ai_dim1] * tr;
          ar[i_ + j * ar_dim1] += prodr;
          ai[i_ + j * ai_dim1] += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i_3 = *n;
        for (i_ = kp1; i_ <= i_3; ++i_) {
          prodr = -ai[i_ + k * ai_dim1] * ti;
          prodi = ar[i_ + k * ar_dim1] * ti;
          ar[i_ + j * ar_dim1] += prodr;
          ai[i_ + j * ai_dim1] += prodi;
          /* L45: */
        }
        goto L48;
      }
      i_3 = *n;
      for (i_ = kp1; i_ <= i_3; ++i_) {
        prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
        prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
        ar[i_ + j * ar_dim1] += prodr;
        ai[i_ + j * ai_dim1] += prodi;
      }
    L48:;
    }
  }
L70:
  k = *n;
  if ((d_1 = ar[*n + *n * ar_dim1], std::abs(d_1)) +
          (d_2 = ai[*n + *n * ai_dim1], std::abs(d_2)) ==
      0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
}

int solc(const int *n, const int *ndim, const double *ar, const double *ai,
         double *br, double *bi, const int *ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2;

  /* Local variables */
  int i_, k, m, kb;
  double ti, tr;
  int km1, nm1, kp1;
  double den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
  /*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    (BR,BI) = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --bi;
  --br;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  if (*n == 1) {
    goto L50;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = ip[k];
    tr = br[m];
    ti = bi[m];
    br[m] = br[k];
    bi[m] = bi[k];
    br[k] = tr;
    bi[k] = ti;
    i_2 = *n;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[i_] += prodr;
      bi[i_] += prodi;
    }
  }
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    km1 = *n - kb;
    k = km1 + 1;
    den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] +
          ai[k + k * ai_dim1] * ai[k + k * ai_dim1];
    prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
    prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
    br[k] = prodr / den;
    bi[k] = prodi / den;
    tr = -br[k];
    ti = -bi[k];
    i_2 = km1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[i_] += prodr;
      bi[i_] += prodi;
    }
  }
L50:
  den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 1];
  prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
  prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
  br[1] = prodr / den;
  bi[1] = prodi / den;
  return 0;
}

int dechc(const int *n, const int *ndim, double *ar, double *ai, const int *lb,
          int *ip, int *ier) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2, i_3;
  double d_1, d_2, d_3, d_4;

  /* Local variables */
  int i_, j, k, m, na;
  double ti, tr;
  int nm1, kp1;
  double den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIstd::minATION */
  /*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
  /*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
  /*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
  /*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    REAL PART. */
  /*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    IMAGINARY PART. */
  /*     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  if (*lb == 0) {
    goto L70;
  }
  if (*n == 1) {
    goto L70;
  }
  nm1 = *n - 1;
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = k;
    /* Computing std::min */
    i_2 = *n, i_3 = *lb + k;
    na = std::min(i_2, i_3);
    i_2 = na;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      if ((d_1 = ar[i_ + k * ar_dim1], std::abs(d_1)) +
              (d_2 = ai[i_ + k * ai_dim1], std::abs(d_2)) >
          (d_3 = ar[m + k * ar_dim1], std::abs(d_3)) +
              (d_4 = ai[m + k * ai_dim1], std::abs(d_4))) {
        m = i_;
      }
      /* L10: */
    }
    ip[k] = m;
    tr = ar[m + k * ar_dim1];
    ti = ai[m + k * ai_dim1];
    if (m == k) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
    ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
    ar[k + k * ar_dim1] = tr;
    ai[k + k * ai_dim1] = ti;
  L20:
    if (std::abs(tr) + std::abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i_2 = na;
    for (i_ = kp1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      ar[i_ + k * ar_dim1] = -prodr;
      ai[i_ + k * ai_dim1] = -prodi;
      /* L30: */
    }
    i_2 = *n;
    for (j = kp1; j <= i_2; ++j) {
      tr = ar[m + j * ar_dim1];
      ti = ai[m + j * ai_dim1];
      ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
      ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
      ar[k + j * ar_dim1] = tr;
      ai[k + j * ai_dim1] = ti;
      if (std::abs(tr) + std::abs(ti) == 0.) {
        goto L48;
      }
      if (ti == 0.) {
        i_3 = na;
        for (i_ = kp1; i_ <= i_3; ++i_) {
          prodr = ar[i_ + k * ar_dim1] * tr;
          prodi = ai[i_ + k * ai_dim1] * tr;
          ar[i_ + j * ar_dim1] += prodr;
          ai[i_ + j * ai_dim1] += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i_3 = na;
        for (i_ = kp1; i_ <= i_3; ++i_) {
          prodr = -ai[i_ + k * ai_dim1] * ti;
          prodi = ar[i_ + k * ar_dim1] * ti;
          ar[i_ + j * ar_dim1] += prodr;
          ai[i_ + j * ai_dim1] += prodi;
          /* L45: */
        }
        goto L48;
      }
      i_3 = na;
      for (i_ = kp1; i_ <= i_3; ++i_) {
        prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
        prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
        ar[i_ + j * ar_dim1] += prodr;
        ai[i_ + j * ai_dim1] += prodi;
        /* L47: */
      }
    L48:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = *n;
  if ((d_1 = ar[*n + *n * ar_dim1], std::abs(d_1)) +
          (d_2 = ai[*n + *n * ai_dim1], std::abs(d_2)) ==
      0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
  /* ----------------------- END OF SUBROUTINE DECHC ----------------------- */
}

int solhc(const int *n, const int *ndim, const double *ar, const double *ai,
          const int *lb, double *br, double *bi, const int *ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2, i_3, i_4;

  /* Local variables */
  int i_, k, m, kb;
  double ti, tr;
  int km1, nm1, kp1;
  double den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
  /*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
  /*    LB = LOWER BANDWIDTH OF A. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    (BR,BI) = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --bi;
  --br;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  if (*n == 1) {
    goto L50;
  }
  nm1 = *n - 1;
  if (*lb == 0) {
    goto L25;
  }
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    kp1 = k + 1;
    m = ip[k];
    tr = br[m];
    ti = bi[m];
    br[m] = br[k];
    bi[m] = bi[k];
    br[k] = tr;
    bi[k] = ti;
    /* Computing std::min */
    i_3 = *n, i_4 = *lb + k;
    i_2 = std::min(i_3, i_4);
    for (i_ = kp1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[i_] += prodr;
      bi[i_] += prodi;
      /* L10: */
    }
    /* L20: */
  }
L25:
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    km1 = *n - kb;
    k = km1 + 1;
    den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] +
          ai[k + k * ai_dim1] * ai[k + k * ai_dim1];
    prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
    prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
    br[k] = prodr / den;
    bi[k] = prodi / den;
    tr = -br[k];
    ti = -bi[k];
    i_2 = km1;
    for (i_ = 1; i_ <= i_2; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[i_] += prodr;
      bi[i_] += prodi;
      /* L30: */
    }
    /* L40: */
  }
L50:
  den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 1];
  prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
  prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
  br[1] = prodr / den;
  bi[1] = prodi / den;
  return 0;
}

int decb(const int *n, const int *ndim, double *a, const int *ml, const int *mu,
         int *ip, int *ier) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3, i_4;
  double d_1, d_2;

  /* Local variables */
  int i_, j, k, m;
  double t;
  int md, jk, mm, ju, md1, nm1, kp1, mdl, ijk;

  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIstd::minATION OF A BANDED */
  /*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
  /*  INPUT.. */
  /*     N       ORDER OF THE ORIGINAL MATRIX A. */
  /*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
  /*     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
  /*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND */
  /*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
  /*                ML+1 THROUGH 2*ML+MU+1 OF  A. */
  /*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*  OUTPUT.. */
  /*     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
  /*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
  /*     IP      INDEX VECTOR OF PIVOT INDICES. */
  /*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
  /*                SINGULAR AT STAGE K. */
  /*  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
  /*  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  md = *ml + *mu + 1;
  md1 = md + 1;
  ju = 0;
  if (*ml == 0) {
    goto L70;
  }
  if (*n == 1) {
    goto L70;
  }
  if (*n < *mu + 2) {
    goto L7;
  }
  i_1 = *n;
  for (j = *mu + 2; j <= i_1; ++j) {
    i_2 = *ml;
    for (i_ = 1; i_ <= i_2; ++i_) {
      /* L5: */
      a[i_ + j * a_dim1] = 0.;
    }
  }
L7:
  nm1 = *n - 1;
  i_2 = nm1;
  for (k = 1; k <= i_2; ++k) {
    kp1 = k + 1;
    m = md;
    /* Computing std::min */
    i_1 = *ml, i_3 = *n - k;
    mdl = std::min(i_1, i_3) + md;
    i_1 = mdl;
    for (i_ = md1; i_ <= i_1; ++i_) {
      if ((d_1 = a[i_ + k * a_dim1], std::abs(d_1)) >
          (d_2 = a[m + k * a_dim1], std::abs(d_2))) {
        m = i_;
      }
      /* L10: */
    }
    ip[k] = m + k - md;
    t = a[m + k * a_dim1];
    if (m == md) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    a[m + k * a_dim1] = a[md + k * a_dim1];
    a[md + k * a_dim1] = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i_1 = mdl;
    for (i_ = md1; i_ <= i_1; ++i_) {
      /* L30: */
      a[i_ + k * a_dim1] = -a[i_ + k * a_dim1] * t;
    }
    /* Computing std::min */
    /* Computing std::max */
    i_3 = ju, i_4 = *mu + ip[k];
    i_1 = std::max(i_3, i_4);
    ju = std::min(i_1, *n);
    mm = md;
    if (ju < kp1) {
      goto L55;
    }
    i_1 = ju;
    for (j = kp1; j <= i_1; ++j) {
      --m;
      --mm;
      t = a[m + j * a_dim1];
      if (m == mm) {
        goto L35;
      }
      a[m + j * a_dim1] = a[mm + j * a_dim1];
      a[mm + j * a_dim1] = t;
    L35:
      if (t == 0.) {
        goto L45;
      }
      jk = j - k;
      i_3 = mdl;
      for (i_ = md1; i_ <= i_3; ++i_) {
        ijk = i_ - jk;
        /* L40: */
        a[ijk + j * a_dim1] += a[i_ + k * a_dim1] * t;
      }
    L45:
        /* L50: */
        ;
    }
  L55:
      /* L60: */
      ;
  }
L70:
  k = *n;
  if (a[md + *n * a_dim1] == 0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
  /* ----------------------- END OF SUBROUTINE DECB ------------------------ */
} /* decb_ */

int solb(const int *n, const int *ndim, const double *a, const int *ml,
         const int *mu, double *b, const int *ip) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3;

  /* Local variables */
  int i_, k, m;
  double t;
  int kb, md, lm, md1, nm1, imd, kmd, mdl, mdm;

  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N      ORDER OF MATRIX A. */
  /*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
  /*    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB. */
  /*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    B      RIGHT HAND SIDE VECTOR. */
  /*    IP     PIVOT VECTOR OBTAINED FROM DECB. */
  /*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B      SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --b;
  a_dim1 = *ndim;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  md = *ml + *mu + 1;
  md1 = md + 1;
  mdm = md - 1;
  nm1 = *n - 1;
  if (*ml == 0) {
    goto L25;
  }
  if (*n == 1) {
    goto L50;
  }
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    m = ip[k];
    t = b[m];
    b[m] = b[k];
    b[k] = t;
    /* Computing std::min */
    i_2 = *ml, i_3 = *n - k;
    mdl = std::min(i_2, i_3) + md;
    i_2 = mdl;
    for (i_ = md1; i_ <= i_2; ++i_) {
      imd = i_ + k - md;
      /* L10: */
      b[imd] += a[i_ + k * a_dim1] * t;
    }
    /* L20: */
  }
L25:
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    k = *n + 1 - kb;
    b[k] /= a[md + k * a_dim1];
    t = -b[k];
    kmd = md - k;
    /* Computing std::max */
    i_2 = 1, i_3 = kmd + 1;
    lm = std::max(i_2, i_3);
    i_2 = mdm;
    for (i_ = lm; i_ <= i_2; ++i_) {
      imd = i_ - kmd;
      /* L30: */
      b[imd] += a[i_ + k * a_dim1] * t;
    }
    /* L40: */
  }
L50:
  b[1] /= a[md + a_dim1];
  return 0;
  /* ----------------------- END OF SUBROUTINE SOLB ------------------------ */
} /* solb_ */

int decbc(const int *n, const int *ndim, double *ar, double *ai, const int *ml,
          const int *mu, int *ip, int *ier) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2, i_3, i_4;
  double d_1, d_2, d_3, d_4;

  /* Local variables */
  int i_, j, k, m, md, jk, mm;
  double ti;
  int ju;
  double tr;
  int md1, nm1, kp1;
  double den;
  int mdl, ijk;
  double prodi, prodr;

  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIstd::minATION OF A BANDED COMPLEX
   */
  /*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
  /*  INPUT.. */
  /*     N       ORDER OF THE ORIGINAL MATRIX A. */
  /*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
  /*     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
  /*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL */
  /*                PART) AND AI (IMAGINARY PART)  AND */
  /*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
  /*                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI. */
  /*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*  OUTPUT.. */
  /*     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
  /*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
  /*     IP      INDEX VECTOR OF PIVOT INDICES. */
  /*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
  /*                SINGULAR AT STAGE K. */
  /*  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
  /*  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  *ier = 0;
  ip[*n] = 1;
  md = *ml + *mu + 1;
  md1 = md + 1;
  ju = 0;
  if (*ml == 0) {
    goto L70;
  }
  if (*n == 1) {
    goto L70;
  }
  if (*n < *mu + 2) {
    goto L7;
  }
  i_1 = *n;
  for (j = *mu + 2; j <= i_1; ++j) {
    i_2 = *ml;
    for (i_ = 1; i_ <= i_2; ++i_) {
      ar[i_ + j * ar_dim1] = 0.;
      ai[i_ + j * ai_dim1] = 0.;
      /* L5: */
    }
  }
L7:
  nm1 = *n - 1;
  i_2 = nm1;
  for (k = 1; k <= i_2; ++k) {
    kp1 = k + 1;
    m = md;
    /* Computing std::min */
    i_1 = *ml, i_3 = *n - k;
    mdl = std::min(i_1, i_3) + md;
    i_1 = mdl;
    for (i_ = md1; i_ <= i_1; ++i_) {
      if ((d_1 = ar[i_ + k * ar_dim1], std::abs(d_1)) +
              (d_2 = ai[i_ + k * ai_dim1], std::abs(d_2)) >
          (d_3 = ar[m + k * ar_dim1], std::abs(d_3)) +
              (d_4 = ai[m + k * ai_dim1], std::abs(d_4))) {
        m = i_;
      }
      /* L10: */
    }
    ip[k] = m + k - md;
    tr = ar[m + k * ar_dim1];
    ti = ai[m + k * ai_dim1];
    if (m == md) {
      goto L20;
    }
    ip[*n] = -ip[*n];
    ar[m + k * ar_dim1] = ar[md + k * ar_dim1];
    ai[m + k * ai_dim1] = ai[md + k * ai_dim1];
    ar[md + k * ar_dim1] = tr;
    ai[md + k * ai_dim1] = ti;
  L20:
    if (std::abs(tr) + std::abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i_1 = mdl;
    for (i_ = md1; i_ <= i_1; ++i_) {
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      ar[i_ + k * ar_dim1] = -prodr;
      ai[i_ + k * ai_dim1] = -prodi;
      /* L30: */
    }
    /* Computing std::min */
    /* Computing std::max */
    i_3 = ju, i_4 = *mu + ip[k];
    i_1 = std::max(i_3, i_4);
    ju = std::min(i_1, *n);
    mm = md;
    if (ju < kp1) {
      goto L55;
    }
    i_1 = ju;
    for (j = kp1; j <= i_1; ++j) {
      --m;
      --mm;
      tr = ar[m + j * ar_dim1];
      ti = ai[m + j * ai_dim1];
      if (m == mm) {
        goto L35;
      }
      ar[m + j * ar_dim1] = ar[mm + j * ar_dim1];
      ai[m + j * ai_dim1] = ai[mm + j * ai_dim1];
      ar[mm + j * ar_dim1] = tr;
      ai[mm + j * ai_dim1] = ti;
    L35:
      if (std::abs(tr) + std::abs(ti) == 0.) {
        goto L48;
      }
      jk = j - k;
      if (ti == 0.) {
        i_3 = mdl;
        for (i_ = md1; i_ <= i_3; ++i_) {
          ijk = i_ - jk;
          prodr = ar[i_ + k * ar_dim1] * tr;
          prodi = ai[i_ + k * ai_dim1] * tr;
          ar[ijk + j * ar_dim1] += prodr;
          ai[ijk + j * ai_dim1] += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i_3 = mdl;
        for (i_ = md1; i_ <= i_3; ++i_) {
          ijk = i_ - jk;
          prodr = -ai[i_ + k * ai_dim1] * ti;
          prodi = ar[i_ + k * ar_dim1] * ti;
          ar[ijk + j * ar_dim1] += prodr;
          ai[ijk + j * ai_dim1] += prodi;
          /* L45: */
        }
        goto L48;
      }
      i_3 = mdl;
      for (i_ = md1; i_ <= i_3; ++i_) {
        ijk = i_ - jk;
        prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
        prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
        ar[ijk + j * ar_dim1] += prodr;
        ai[ijk + j * ai_dim1] += prodi;
      }
    L48:;
    }
  L55:;
  }
L70:
  k = *n;
  if ((d_1 = ar[md + *n * ar_dim1], std::abs(d_1)) +
          (d_2 = ai[md + *n * ai_dim1], std::abs(d_2)) ==
      0.) {
    goto L80;
  }
  return 0;
L80:
  *ier = k;
  ip[*n] = 0;
  return 0;
}

int solbc(const int *n, const int *ndim, const double *ar, const double *ai,
          const int *ml, const int *mu, double *br, double *bi, const int *ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i_1, i_2, i_3;

  /* Local variables */
  int i_, k, m, kb, md, lm;
  double ti, tr;
  int md1, nm1;
  double den;
  int imd, kmd, mdl, mdm;
  double prodi, prodr;

  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B , */
  /*                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION. */
  /*  INPUT.. */
  /*    N      ORDER OF MATRIX A. */
  /*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
  /*    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
   */
  /*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART). */
  /*    IP     PIVOT VECTOR OBTAINED FROM DECBC. */
  /*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART). */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --ip;
  --bi;
  --br;
  ai_dim1 = *ndim;
  ai_offset = 1 + ai_dim1;
  ai -= ai_offset;
  ar_dim1 = *ndim;
  ar_offset = 1 + ar_dim1;
  ar -= ar_offset;

  /* Function Body */
  md = *ml + *mu + 1;
  md1 = md + 1;
  mdm = md - 1;
  nm1 = *n - 1;
  if (*ml == 0) {
    goto L25;
  }
  if (*n == 1) {
    goto L50;
  }
  i_1 = nm1;
  for (k = 1; k <= i_1; ++k) {
    m = ip[k];
    tr = br[m];
    ti = bi[m];
    br[m] = br[k];
    bi[m] = bi[k];
    br[k] = tr;
    bi[k] = ti;
    /* Computing std::min */
    i_2 = *ml, i_3 = *n - k;
    mdl = std::min(i_2, i_3) + md;
    i_2 = mdl;
    for (i_ = md1; i_ <= i_2; ++i_) {
      imd = i_ + k - md;
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[imd] += prodr;
      bi[imd] += prodi;
      /* L10: */
    }
    /* L20: */
  }
L25:
  i_1 = nm1;
  for (kb = 1; kb <= i_1; ++kb) {
    k = *n + 1 - kb;
    den = ar[md + k * ar_dim1] * ar[md + k * ar_dim1] +
          ai[md + k * ai_dim1] * ai[md + k * ai_dim1];
    prodr = br[k] * ar[md + k * ar_dim1] + bi[k] * ai[md + k * ai_dim1];
    prodi = bi[k] * ar[md + k * ar_dim1] - br[k] * ai[md + k * ai_dim1];
    br[k] = prodr / den;
    bi[k] = prodi / den;
    tr = -br[k];
    ti = -bi[k];
    kmd = md - k;
    /* Computing std::max */
    i_2 = 1, i_3 = kmd + 1;
    lm = std::max(i_2, i_3);
    i_2 = mdm;
    for (i_ = lm; i_ <= i_2; ++i_) {
      imd = i_ - kmd;
      prodr = ar[i_ + k * ar_dim1] * tr - ai[i_ + k * ai_dim1] * ti;
      prodi = ai[i_ + k * ai_dim1] * tr + ar[i_ + k * ar_dim1] * ti;
      br[imd] += prodr;
      bi[imd] += prodi;
      /* L30: */
    }
    /* L40: */
  }
  den =
      ar[md + ar_dim1] * ar[md + ar_dim1] + ai[md + ai_dim1] * ai[md + ai_dim1];
  prodr = br[1] * ar[md + ar_dim1] + bi[1] * ai[md + ai_dim1];
  prodi = bi[1] * ar[md + ar_dim1] - br[1] * ai[md + ai_dim1];
  br[1] = prodr / den;
  bi[1] = prodi / den;
L50:
  return 0;
  /* ----------------------- END OF SUBROUTINE SOLBC ------------------------ */
} /* solbc_ */

int elmhes(const int *nm, const int *n, const int *low, const int *igh,
           double *a, int *int_) {
  /* System generated locals */
  int a_dim1, a_offset, i_1, i_2, i_3;
  double d_1;

  /* Local variables */
  int i_, j, m;
  double x, y;
  int la, mm1, kp1, mp1;

  a_dim1 = *nm;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --int_;

  /* Function Body */
  la = *igh - 1;
  kp1 = *low + 1;
  if (la < kp1) {
    goto L200;
  }

  i_1 = la;
  for (m = kp1; m <= i_1; ++m) {
    mm1 = m - 1;
    x = 0.;
    i_ = m;

    i_2 = *igh;
    for (j = m; j <= i_2; ++j) {
      if ((d_1 = a[j + mm1 * a_dim1], std::abs(d_1)) <= std::abs(x)) {
        goto L100;
      }
      x = a[j + mm1 * a_dim1];
      i_ = j;
    L100:;
    }

    int_[m] = i_;
    if (i_ == m) {
      goto L130;
    }
    /*    :::::::::: interchange rows and columns of a :::::::::: */
    i_2 = *n;
    for (j = mm1; j <= i_2; ++j) {
      y = a[i_ + j * a_dim1];
      a[i_ + j * a_dim1] = a[m + j * a_dim1];
      a[m + j * a_dim1] = y;
      /* L110: */
    }

    i_2 = *igh;
    for (j = 1; j <= i_2; ++j) {
      y = a[j + i_ * a_dim1];
      a[j + i_ * a_dim1] = a[j + m * a_dim1];
      a[j + m * a_dim1] = y;
      /* L120: */
    }
  /*    :::::::::: end interchange :::::::::: */
  L130:
    if (x == 0.) {
      goto L180;
    }
    mp1 = m + 1;

    i_2 = *igh;
    for (i_ = mp1; i_ <= i_2; ++i_) {
      y = a[i_ + mm1 * a_dim1];
      if (y == 0.) {
        goto L160;
      }
      y /= x;
      a[i_ + mm1 * a_dim1] = y;

      i_3 = *n;
      for (j = m; j <= i_3; ++j) {
        /* L140: */
        a[i_ + j * a_dim1] -= y * a[m + j * a_dim1];
      }

      i_3 = *igh;
      for (j = 1; j <= i_3; ++j) {
        /* L150: */
        a[j + m * a_dim1] += y * a[j + i_ * a_dim1];
      }

    L160:;
    }

  L180:;
  }

L200:
  return 0;
  /*    :::::::::: last card of elmhes :::::::::: */
} /* elmhes_ */
} // namespace stiff

#endif // DARKSUN_DC_DECSOL_HPP
